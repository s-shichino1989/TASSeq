#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import unicode_literals

'''
Module for Fast annotating and parse BAM file.

@author: Lenore Pafford
mlpafford@cellular-research.com
@author: Shigeyuki Shichino
s_shichino@rs.tus.ac.jp modified at 2020.05.01

'''
import argparse
import pandas as pd
import time
import sys
import os
import utils
import glob
import regex as re
import _version as _v
import gzip
import pysam
import csv
import gzip

def main():

    des = 'Add to Sam, ' + _v.desc
    parser = argparse.ArgumentParser(description=des, add_help=True)
    
    #only require annotated R1 and R2 files
    parser.add_argument('--annot-R1', action='store', dest='annotR1', required=True,
                        help='Read 1 annotation files generated by AnnotateR1 node')
    parser.add_argument('--bam', action='store', dest='R2_bam', required=True,
                        help='Bowtie2 generated BAM file')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()

    #start = time.strftime("%Y-%m-%d %H:%M:%S")
    #print ("Start adding annotations to SAM: {}".format(start))
    annotR1 = args.annotR1
    matching_R1 = annotR1

    # output and intermediate files
    R2_BAM = args.R2_bam
    final_txt = "{}_count_R2.txt".format(os.path.basename(R2_BAM).split('_mapping_R2')[0])

    addR1toSAM(matching_R1, R2_BAM, final_txt)
    
    #os.remove(matching_R1)
    #os.remove(R2_BAM)
    #end = time.strftime("%Y-%m-%d %H:%M:%S")
    #print ("Finished adding annotations to sam and reshaping: {}".format(end))
    return 0

def addR1toSAM(annotR1, R2_bam, final_txt):
    """Add Annotations from experiment to the SAM read line it is associated with"""

    #print ("Adding R1 information to parsed BAM file...")
    fnew = open(final_txt, mode='wt', encoding='utf-8')
    fread_R1 = gzip.open(annotR1, mode='rt', encoding='utf-8')
    fbam = pysam.AlignmentFile(R2_bam, 'rb', threads=4)
    ra1 = csv.writer(fnew, delimiter='\t', lineterminator='\n')
    for rname in fbam:
        try:
            readAnnot_R1 = next(fread_R1)
            readAnnot_R1 = readAnnot_R1.rstrip().split(',')
        except StopIteration:
            #print ("End of annotR1, now exit loop")
            break

        # cell label to index
        if 'x' in readAnnot_R1[0]:
            cell_index = '*'
        else:
            cell_index = utils.label2index(readAnnot_R1[0])

        # mol_label
        mol_label = readAnnot_R1[2]
        if mol_label:
            mol_label = mol_label
        else:
            mol_label = '0'

        #parse reference_name (gene symbol)
        ref_name = rname.get_tag("XF")
        if ref_name:
            ref_name = ref_name
            if ref_name.startswith("__"):
            # if annotation (gene symbol) is empty or ambiguous
               ref_name = '*'
            else:
               ref_name = ref_name
        else:
            ref_name = '*'
        ra1.writerow([ref_name, cell_index])
    
    fnew.close()
    fread_R1.close()

    return

if __name__ == '__main__':
    sys.exit(main())
