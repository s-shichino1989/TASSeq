#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
AnnotateReads is based off of the functions realting to reads,
originally located in files GetAnnotation and ResolveAnalysis.
"""

from __future__ import unicode_literals
import argparse
import csv
import time
import sys, os, shutil
import cellkeys
import Levenshtein
import _version as _v
from difflib import SequenceMatcher
import utils
import gzip

def package_main():
    main(sys.argv[1:])


def main(argv):

    des = 'Annotate R1, ' + _v.desc
    parser = argparse.ArgumentParser(description=des, add_help=True)

    parser.add_argument('--R1', action='store', dest='R1', required=True,
                        help='Read 1 FASTQ file from the split_fastqs node.')
    parser.add_argument('--label-version', action='store', type=int, dest='label_version', default=2, choices=[1, 2, 3,4],
                        help="Specify which version of the cell label you are using: '1' for 8mer, "
                             "'2' for 9mer (default), '3' for Precise targeted, '4' for Precise WTA.")

    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)
    R1 = args.R1
    sample = os.path.basename(R1).split('_R1_')[0]
    label = args.label_version

    # output
    annotation_R1 = "./{}_Annotation_R1.csv.gz".format(sample)

    #start = time.strftime("%Y-%m-%d %H:%M:%S")
    #print ("Start read 1 annotation: {}".format(start))

    # Precise
    if label in (3,4):
        annotateR1_precise(annotation_R1, unzipped_R1)

    # Rhapsody
    else:
        # define the default sections, such as startpos, reference seqs depending on label version
        default_section_startpos, appended_startpos, refs, appended_refs = getDefaultSections(label)

        with gzip.open(annotation_R1, mode='wt', compresslevel=1) as f, gzip.open(R1, mode='rt', encoding='utf-8') as r:
            ra1 = csv.writer(f, lineterminator='\n')
            inx = 0
            for read in parse_fastq(r):
                inx += 1
                #if inx % 10000000 == 0:
                #    print ("Annotated", inx, "reads of R1")

                # get cell labels and indels to determine MI start position
                cells, mismatch, mol_polyT_idx = checkMatches(read, refs, appended_refs, default_section_startpos,
                                                              appended_startpos)

                # define MI sequence and polyT
                mol = read[mol_polyT_idx[0]:mol_polyT_idx[1]]
                polyT = 'T'

                ra1.writerow([cells, mismatch, mol, polyT])

        #print ("Annotated all", inx, "reads of R1")

    time.sleep(1)
    os.remove(R1)
    
    #end = time.strftime("%Y-%m-%d %H:%M:%S")
    #print ("Finished read 1 annotation: {}".format(end))
    
    return 0

def parse_fastq(f):

    while True:
       try:
        next(f)
       except:
        return
       else:
        read = str(next(f).strip())
        next(f)
        next(f)
        yield read

def errorcorrect_96(barcode):
    """
    function that corrects sample barcodes with single substitution errors to correct barcode
    using modulo 4 properties/checksums of Hamming barcodes [Bystrykh 2012]
    input: sample barcode (wrong or correct)
    output: correct sample barcode
    """

    dict = {"A": 0, "C": 1, "G": 2, "T": 3}
    reverse = {dict[x]: x for x in dict}

    bc = [dict[letter] for letter in barcode]

    p_1 = (bc[0] + bc[2] + bc[4] + bc[6]) % 4
    p_2 = (bc[1] + bc[2] + bc[5] + bc[6]) % 4
    p_3 = (bc[3] + bc[4] + bc[5] + bc[6]) % 4
    p_4 = sum(bc) % 4

    p_1b = int(p_1 > 0)
    p_2b = int(p_2 > 0)
    p_3b = int(p_3 > 0)
    p_4b = int(p_4 > 0)

    binary = [p_4b, p_3b, p_2b, p_1b]
    binary = [str(x) for x in binary]
    binary = ''.join(binary)

    error_type = max([p_1, p_2, p_3, p_4])
    if not error_type:
        return barcode

    else:
        error_pos = int(binary, 2) % 8 - 1
        error_val = bc[error_pos]
        true_val = (error_val - error_type) % 4

        bc[error_pos] = true_val
        bc = [reverse[x] for x in bc]
        bc = "".join(bc)

        return bc


def getDefaultSections(label):
    """ return the default starting position of each section, including the CL, linker, UMI and polyT.
        also return the reference sequences from cellkeys, as well as the polyT cutoff """

    if label == 1:
        start_pos = [0, 8, 20, 28, 40, 48, 56, 64]
        appended_startpos = [0, 20, 40, 48, 56, 64]

        refs = [[str(ref) for ref in cellkeys.cell_key1[:96]], [str(ref) for ref in cellkeys.linker1],
                [str(ref) for ref in cellkeys.cell_key2[:96]], [str(ref) for ref in cellkeys.linker2],
                [str(ref) for ref in cellkeys.cell_key3[:96]]]

    elif label == 2:
        start_pos = [0, 9, 21, 30, 43, 52, 60, 68]
        appended_startpos = [0, 21, 43, 52, 60, 68]

        refs = [[str(ref) for ref in cellkeys.v2_cell_key1[:96]], [str(ref) for ref in cellkeys.linker1],
                [str(ref) for ref in cellkeys.v2_cell_key2[:96]], [str(ref) for ref in cellkeys.linker2],
                [str(ref) for ref in cellkeys.v2_cell_key3[:96]]]

    # CL1 + L1
    cl1_l1 = [ref + refs[1][0] for ref in refs[0]]

    # CL2 + L2
    if label == 1:
        cl2_l2 = [ref + refs[3][0] for ref in refs[2]]
    else:
        cl2_l2 = [str(ref + refs[3][0] + 'A') for ref in refs[2]]

    # CL3 alone
    cl3 = refs[4]

    appended_refs = [cl1_l1, cl2_l2, cl3]

    return start_pos, appended_startpos, refs, appended_refs


def checkMatches(seq, refs, appended_refs, section_startpos, appended_startpos):
    """Returns the assigned cell label section index, or 'x' if no match found.
    Returns the number of mismatches for each cell label section: '0' is a perfect match.
    Returns the start position of each cell label section.

    Rhapsody:
     1. Check for perfect match (PM) for Cell Label sections only (CL1, Cl2, CL3). Return if no mismatches found.

     2. Find the section that is not PM, and check one section prior to it and all later sections for 2nd round of
     matching by taking into account of indels in previous sections and shift later sections as appropriate.

     """

    # check for PMs in CL sections only
    CL_sections = []
    for i in [0, 2, 4]:
        sections = seq[section_startpos[i]:section_startpos[i+1]]
        try:
            CL_sections += str(refs[i].index(sections) + 1),
        except ValueError:
            CL_sections += 'x',

    # if all CL sections are perfect, return
    if 'x' not in CL_sections:
        return '-'.join(CL_sections), '0-0-0', section_startpos[5:]

    # if has mismatch in any of CL sections, account for indels to find the best match
    # combine CL1+L1 as well as CL2+L2, and use appended_startpos for combined sections
    section_startpos = list(appended_startpos)
    section_dists = ['0', '0', '0']
    indels = ['0', '0', '0']

    # smaller threshold is used for cell label section 3 as only cell label is considered (no linker appended)
    allowed_mismatch = [4, 4, 2]

    # find first section with x and start one section before. shift all following sections if indel found
    first_section = max(CL_sections.index('x') - 1, 0)
    for i in range(first_section, 3):

        append_seq = seq[section_startpos[i]: section_startpos[i+1]]

        try:
            CL_sections[i] = str(appended_refs[i].index(append_seq) + 1)

        except ValueError:
            CL_sections[i], section_dists[i], indels[i] = map_cell_sections(append_seq, appended_refs[i],
                                                                            allowed_mismatch[i])

            if indels[i] != 0:
                # if have indels in previous section, shift the start pos for all later sections
                section_startpos[(i + 1):] = [x - indels[i] for x in section_startpos[(i + 1):]]

            # if still no match, exit early
            if CL_sections[i] == 'x':
                mol_polyT_idx = section_startpos[3:]
                return '-'.join(CL_sections), '-'.join(section_dists), mol_polyT_idx

    mol_polyT_idx = section_startpos[3:]
    return '-'.join(CL_sections), '-'.join(section_dists), mol_polyT_idx


def find_shift(seq1, seq2):
    """Given two seqs, get the operations needed to convert one seq to another; get the shift btw the two seqs.
       Example: seq1 = GGTAGCGGTGAC, seq2 = GTAGCGGCTGAC
                operations to turn seq1 to seq2 by get_opcodes() are:
                [('delete', 0, 1, 0, 0), ('equal', 1, 8, 0, 7), ('insert', 8, 8, 7, 8), ('equal', 8, 12, 8, 12)]
                The usage and output of get_opcodes() can be found at
                https://docs.python.org/2/library/difflib.html#difflib.SequenceMatcher.get_opcodes
                check the last operation in the output list of get_opcodes() and determine the num of shift btw two
                seqs: in the above example, num of shift is 8 - 8 = 0."""

    match = SequenceMatcher(None, seq1, seq2).get_opcodes()
    num_shift = match[-1][1] - match[-1][3]
    return num_shift


def map_cell_sections(seq, references, allowed_mismatches):
    """Given a seq, search against the candidate match in targets by calculating the edit distance"""

    scores = []
    for target in references:
        score = Levenshtein.distance(seq, target)
        if score == 1 and 1 in scores:  # if multiple editDist of 1, return  early
            return 'x', 'x', 0
        scores += score,

    mismatch = min(scores)

    if mismatch <= allowed_mismatches and scores.count(mismatch) == 1:
        cl_section = scores.index(mismatch) + 1
        indels = find_shift(references[cl_section - 1], seq)
        return str(cl_section), str(mismatch), indels

    else:
        return 'x', 'x', 0


if __name__ == '__main__':
    sys.exit(package_main())
