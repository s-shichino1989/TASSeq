import colorsys
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pandas as pd
import regex as re
import csv
import subprocess32 as subprocess
import os
import shutil
import gzip
import sys

def clean_up_decimals(metric_list):
    """ returns list with ints if whole number, otherwise float rounded to 2 decimal places
    All list members must be floats"""

    new_list = [int(x) if x.is_integer() else round(x, 2) for x in metric_list]
    return new_list


# Load csv file (skipcol = 4 for rhapsody)
def load_csv(filename,skipcol=4):
    X = np.genfromtxt(filename,delimiter=',',skip_header=1)
    X = X[:,skipcol:]
    print ('We have %d samples and %d features per sample.'%np.shape(X))
    print ('%.3f%% of the entries are 0.'%(100*np.sum(X == 0)/float(np.prod(np.shape(X)))))
    return X

# Check if all entries in a numpy matrix are whole numbers
def all_entries_are_whole_nums(X):
    a = np.product(np.shape(X))
    b = sum([1 for i in np.ndarray.flatten(X) if float.is_integer(i)])
    return a == b

# When printing, decide whether to use scientific notation (if value gets too big)
def sn(i,j=2):
    if i > 1000 or i < 0.01: return '%.1E'%(i)
    s = '%.'+str(j)+'f'
    return s%(i)

# Simple way to plot colors and labels with legend                                        
def plot_labels_legend(x1,x2,y,labels=None,title=None,save_name=None,label_singletons=True):
    if label_singletons:
        y=map_singleton_to_label(y)
    Ncolors = len(np.unique(y))
    if labels is None: labels = np.unique(y)
    HSVs = [(x*1.0/Ncolors, 0.8, 0.9) for x in range(Ncolors)]
    RGBs = map(lambda x: colorsys.hsv_to_rgb(*x), HSVs)
    for j,i in enumerate(labels):
        if i != -1:
            lab = 'Cluster #{}\nCells: {}'.format(str(i), str(np.sum(y==i)))
            plt.plot(x1[y==i], x2[y==i], '.', c=RGBs[j], label=lab)
        else:
            plt.plot(x1[y==i], x2[y==i], '.', c=RGBs[j], label='Singletons:'+str(np.sum(y==i)), markeredgecolor='k')
    _ = plt.axis('off')
    plt.subplots_adjust(left=0.05, right=0.7, top=0.95, bottom=0.05)
    plt.legend(bbox_to_anchor=(1.4, 1.0), prop={'size': 10}, frameon=False, columnspacing=2, numpoints=1,
               markerscale=3, labelspacing=1)
    if title:
        plt.title(title)
    if save_name is not None:
        plt.savefig(save_name+'.png', format='png', bbox_inches='tight', dpi=300)

# For each feature, determine the index of the cluster with the greatest expression
def compare_feature_means(X,Y):
    Nc = len(np.unique(Y))
    N,M = np.shape(X)
    m = np.zeros((Nc,M))
    for i,c in enumerate(np.unique(Y)):
        m[i,:] = np.mean(X[Y == c,:],0)
    return np.argmax(m,0)

# Map labels to integers
def str_labels_to_ints(y_str):
    y_int = np.zeros(len(y_str))
    for i,label in enumerate(np.unique(y_str)):
        y_int[y_str == label] = i
    return y_int.astype(int)

# Get all off-diagonal entries in a distance matrix
def flatten_distance_matrix(D,inds=None):
    if inds is not None: D2 = cut_matrix_along_both_axes(D,inds)
    else: D2 = D
    d = D2.reshape(-1,1)
    return np.delete(d,np.where(d==0)[0])

# index a 2D matrix along both axes
def cut_matrix_along_both_axes(X,inds):
    Z = X[inds,:]
    return Z[:,inds]

# map singleton string labels to -1
def map_singleton_to_label(y):
    for i,c in enumerate(np.unique(y)):
        if np.sum(y == c) == 1:
            y[y == c] = -1
    return y


# return median of a correlation distribution
def median_cdist_corr(D, i, j, z):
    dM = flatten_distance_matrix(D, np.logical_or(z == i, z == j))
    return np.median(dM)


def label2index(label):
    """return the appropriately formatted cell_index based on input cell_label"""
    if '-' in label:
        cell_label = [int(n) for n in label.split('-')]
        return str((cell_label[0] - 1) * 96 * 96 + (cell_label[1] - 1) * 96 + cell_label[2])
    else:
        return str(label)


def get_assay(fname):
    if 'WTA' in fname:
        assay = 'WTA'

    elif 'Targeted' in fname:
        assay = 'Targeted'

    return assay


def findExpressedGenes(dt, len_header):
    """Drops the genes not expressed in data tables."""
    df = pd.read_csv(dt, header=len_header)
    not_expressed = []
    df.set_index('Cell_Index', inplace=True)
    for gene in df:
        if gene.endswith('pAbO') or np.sum(df[gene]) == 0:  # only include genes, exclude pAbOs
            not_expressed.append(gene)
    df = df.drop(not_expressed, axis=1)
    return df


def findExpressedAbOs(dt, len_header):
    """Drops the genes not expressed in data tables."""
    df = pd.read_csv(dt, header=len_header)
    not_expressed = []
    df.set_index('Cell_Index', inplace=True)
    for gene in df:
        if not gene.endswith('pAbO') or np.sum(df[gene]) == 0:
            not_expressed.append(gene)
    df = df.drop(not_expressed, axis=1)
    return df


def openMetricsFile(metrics_file):
    """ Fetches metrics from intermediate metrics files. Skip over main header, Grabs as many metrics rows as it finds.
    Files must follow the format with a section header followed by metrics + \n
    """
    def is_number_or_na(putative_number):
        if putative_number == 'NA':
            return True
        else:
            try:
                float(putative_number)
            except ValueError:
                return False
            else:
                return True

    if metrics_file:  # algo file for non-precise pipeline
        with open(metrics_file, 'r') as f:
            reader = csv.reader(f)
            non_header_rows = [row for row in reader if is_number_or_na(row[0])]
            return non_header_rows


def grab_main_header(file):
    # get main output header from file
    with open(file, 'r') as f:
        output_header = []
        reader = csv.reader(f)
        for line in reader:
            if line[0].startswith("##"):
                output_header.append([str(line[0]).strip()])
            else:
                break

    if 'WTA' in output_header[1][0]:
        assay = 'WTA'
    else:
        assay = 'Targeted'

    if 'Rhapsody' in output_header[1][0]:
        label_version = 2
    else:
        label_version = 4 if assay == 'WTA' else 3

    if 'Multiplex' in output_header[1][0]:
        trueno = True
    else:
        trueno = False

    sample = output_header[3][0].split(': ')[1]
    reference = output_header[4][0].split(': ')[1]
    bam_input = None

    for row in output_header[5:]:
        if 'Bam Input' in row[0]:
            bam_input = row[0].split(': ')[1]
        elif 'Cell Label' in row[0]:
            cell_label = row[0].split(': ')[1]
            if cell_label == '8-mer':
                label_version = 1

    run_info = [assay, label_version, sample, reference, bam_input, trueno]
    return output_header, run_info


def len_match(CIGAR):
    '''Take input of CIGAR string from SAM file and return length of match'''
    re_match = re.findall("([0-9]*)M", CIGAR)
    matched = [int(match) for match in re_match]
    len_match = sum(matched)
    return len_match


def execute_shell_commands(cmds):
    # Spawn shell processes
    processes = [subprocess.Popen(cmd, shell=True) for cmd in cmds]
    # Wait for their completion
    for process in processes:
        process.wait()
    return 0


def get_control(cfile):

    if cfile == 'phix':
        control_file = 'phix_genome.fasta'

    elif cfile == 'precise':
        control_file = 'precise_internal_control.fasta'

    elif cfile.lower() in ['hs', 'human']:
        control_file = 'SampleTagSequences_HomoSapiens_ver1.fasta'

    elif cfile.lower() in ['mm', 'mouse']:
        control_file = 'SampleTagSequences_MusMusculus_ver1.fasta'

    else:
        print ('Error: Unknown Sample Tag version selection. \n' \
              'For human specify any of: hs or human .\n' \
              'For mouse specify any of: mm or mouse.\n' \
              'Not case-sensitive.')
        sys.exit()

    return '/mist/control_files/' + control_file


def plot_signal_noise_histogram(noise_list, signal_list, plot_name, plot_title, xlab='RSEC reads', ylab='Histogram', colors=['blue', 'green'], labels=['Noise', 'Signal'], n_bins=15):
    '''given the signal and noise data stored in two separate array/list, plot the histogram for both signal and noise in one plot '''
    x_multi = [x for x in [noise_list, signal_list]]
    plt.hist(x_multi, n_bins, normed=1, histtype='bar', color=colors, label=labels)
    plt.legend(loc='upper right')
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(plot_title)
    plt.savefig(plot_name, format='png', dpi=300)
    plt.close()


def csv_input(fps):
    for fp in fps:
        with open(fp, mode='rU') as f:
            csv_reader = csv.reader(f)
            for row in csv_reader:
                yield row


def extract_name(fp, regex):
    """template function to extract names and fail gracefully"""
    base_name = os.path.basename(fp)
    try:
        lib_name_matches = re.findall(regex, base_name)
        lib_name = lib_name_matches[0]
    except IndexError:
        if len(fp) == 0 or len(base_name):
            raise NameError('Invalid filepath: {}'.format(fp))
        else:
            raise NameError('Unrecognized naming convention for: {}'.format(base_name))
    else:
        return lib_name


def extract_library_name_r1(fp):
    return extract_name(fp, r'(?<=-).*?(?=_[0-9]*_Annotation_R1\.csv|_subsampled_Annotation_R1\.csv)')


def remove_illumina_tags(fp):
    return extract_name(fp, r'^.*?(?=_S[0-9]*_|_L[0-9]{3}_|_R[12]_)')


def extract_library_name_quality_filter(fp):
    return extract_name(fp, r'^.*(?=_.*-.*_read_quality\.csv)')


def collect_quality_yield_metrics(r2_bam, fp_out=None):
    """wrapper for picard tools' CollectQualityYieldMetrics program"""
    qual_file = 'quality_yield_metrics.txt'

    subprocess.call(['java', '-jar', '/opt/sequencing/bin/picard.jar', 'CollectQualityYieldMetrics',
                     'INPUT={}'.format(r2_bam), 'OUTPUT={}'.format(qual_file)])

    with open(qual_file) as qmet:
        qmet_reader = csv.reader(qmet, delimiter=str(u'\t').encode('utf-8'))  # replace with '\t' after upgrade to py 3
        for i in xrange(6):
            next(qmet_reader)
        header_row = next(qmet_reader)
        data_row = next(qmet_reader)
        quality_metrics = {header_entry: float(data_entry) for header_entry, data_entry in zip(header_row, data_row)}

        ratio_q30_bases = float(quality_metrics['PF_Q30_BASES']) / float(quality_metrics['PF_BASES'])
        pct_q30_bases = ratio_q30_bases * 100.0
        quality_metrics['PCT_Q30_BASES'] = pct_q30_bases

    if fp_out is not None:
        with open(fp_out, mode='w+') as q30_stats_f_out:
            qmet_writer = csv.DictWriter(q30_stats_f_out, fieldnames=quality_metrics.keys())
            qmet_writer.writeheader()
            qmet_writer.writerow(quality_metrics)

    return quality_metrics


def plot_signal_noise_histogram(noise_list, signal_list, plot_name, plot_title, xlab='RSEC reads', ylab='Histogram', colors=['blue', 'green'], labels=['Noise', 'Signal'], n_bins=15):
    '''given the signal and noise data stored in two separate array/list, plot the histogram for both signal and noise in one plot '''
    x_multi = [x for x in [noise_list, signal_list]]
    plt.hist(x_multi, n_bins, normed=1, histtype='bar', color=colors, label=labels)
    plt.legend(loc='upper right')
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(plot_title)
    plt.savefig(plot_name, format='png', dpi=300)
    plt.close()


def compress_file(_file):
    """

    Args:
        file: file to be compressed

    Returns: path of the compressed file

    """
    out_file = '{}.gz'.format(_file)
    with open(_file, 'rb') as f_in, gzip.open(out_file, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(_file)
    return out_file

def concate_files(file_list):
    target_file = file_list[0]
    with open(target_file, 'wb') as tf:
        for i in xrange(1, len(file_list)):
            with open(file_list[i], 'rb') as fh:
                shutil.copyfileobj(fh, tf)

    return target_file
