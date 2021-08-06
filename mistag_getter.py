import glob, os
import argparse

from mistag_utils import *


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', nargs=1, required=True, help='Input fasta \
                        file name [or library folder containing the \
                        `_fwd.fastq` and `_rev.fastq` files and the \
                        `_mistag_R1.fastq` and `_mistag_R2.fastq` files] \
                        (required))')
    parser.add_argument('-d', nargs=1, required=True, help='Multiplexing \
                        design file name (required)')
    parser.add_argument('-o', nargs='?', default=argparse.SUPPRESS, help =
                        "Output fasta file (default=input appended with \
                        'mistagFiltered.fasta')")
    parser.add_argument('-a', nargs='?', default=0.05, type=float,
                        metavar='float between 0 and 1, max. 3 decimals',\
                        choices=[float(x)/1000 for x in range(1,1001)], help =
                        "Alpha level for finding the Student's T critical value \
                        for the modified Thompson Tau rejection region \
                        calulation (default=0.05)")
    parser.add_argument('--out', action='store_true', default=False, help =
                        "Leave expected sample sequences out of non-critical \
                        mistags distribution for calculations of the rejection \
                        region (default=not active)")
    parser.add_argument('-m', nargs='?', choices=['pandaseq', 'vsearch'],
                        default='pandaseq', help='Reads merging software - \
                        which must be installed and running (default=pandaseq)')
    parser.add_argument('--fastq', action='store_true', default=False,
                        help="Output per-sample fastq files")
    parse = parser.parse_args()
    args = vars(parse)
    return args


def get_fastqs(fastin, design):
    """Return a dict with per-sample _fwd.fastq and _rev.fastq files as a list
    value under the generic sample name as key
    """
    samples_fastqs = {}
    # get the fwd/rev fasta files per sample
    files = glob.glob('%s/*_fwd.fastq' % fastin)
    for f in sorted(files):
        fsplit = f.split('_')
        samples_fastqs['_'.join(fsplit[:-1])] = [f, f.replace('_fwd.fastq', '_rev.fastq')]
    # invert the combito-sample to samples-to-combi dict
    design_rev = rev_dict(design)
    print()
    print()
    print(design_rev)
    # add the sample name for each fwd/rev fastqs list
    add_to_dict_list(samples_fastqs, design_rev)
    print()
    print()
    print(samples_fastqs)

    combi_fastqs = {}
    for _, (fwd, rev, cmb) in samples_fastqs.items():
        print(fwd, rev, cmb)
        if cmb in combi_fastqs:
            print('%s affected twice...' % ' + '.join(list(cmb)))
            sys.exit()
        combi_fastqs[cmb] = [fwd, rev]

    return samples_fastqs, combi_fastqs


def get_mistags_files(fastin):
    # get unexpected samples information
    fs = glob.glob('%s/*fastq' % fastin)
    mistag_fastqs = [f for f in fs if re.search(r'_mistag_R[1-2]\.fastq$', f)]
    if len(mistag_fastqs) == 2:
        return mistag_fastqs
    print('No pair of mistag files detected.\nExiting...')
    sys.exit()


def get_output_fp(args):
    f = args['i'][0]
    if 'o' in args:
        o_fas = args['o']
        o_rad = '%s_mistagFilt' % os.path.splitext(o_fas)[0]
    else:
        curdir = os.path.split(os.path.abspath(f))
        if os.path.isdir(f):
            o_rad = '%s/%s_mistagFilt' % (os.path.abspath(f), curdir[1])
        else:
            o_rad = '%s_mistagFilt' % os.path.splitext(os.path.abspath(f))[0]
        o_fas='%s.fasta' % o_rad
    return o_rad, o_fas


def get_saturation(design, primers):
    n_f = len(primers['F'])
    n_r = len(primers['R'])
    fr = float(n_f * n_r)
    n_s = float(len(design))
    sat_fraction = (n_s / fr) * 100
    sat_percent = round(sat_fraction, 2)
    return fr, sat_percent