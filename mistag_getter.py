import glob, os
import argparse

from mistag_utils import *

def get_args():
    parser=argparse.ArgumentParser()
    parser.add_argument('-i', nargs = 1, required = True, help = 'Input fasta \
                        file name [or folder if pre-formatting needed] \
                        (required)')
    parser.add_argument('-d', nargs = 1, required = True, help = 'Multiplexing \
                        design file name (required)')
    parser.add_argument('-o', nargs = '?', default = argparse.SUPPRESS, help =
                        "Output fasta file (default = input appended with \
                        'mistagFiltered.fasta')")
    parser.add_argument('-a', nargs = '?', default = 0.05, type = float,
                        metavar='float between 0 and 1, max. 3 decimals',\
                        choices = [float(x)/1000 for x in range(1,1001)], help =
                        "Alpha level for finding the Student's T critical value \
                        for the modified Thompson Tau rejection region \
                        calulation (default = 0.05)")
    parser.add_argument('--out', action = 'store_true', default = False, help =
                        "Leave expected sample sequences out of non-critical \
                        mistags distribution for calculations of the rejection \
                        region (default = not active)")
    parser.add_argument('-m', nargs = '?', choices = ['pandaseq', 'vsearch'],
                        default = 'pandaseq', help = 'Reads merging software - \
                        which must be installed and running (default = pandaseq)')
    parse=parser.parse_args()
    args=vars(parse)
    return args


def get_fastqs(fastin, design):
    """Return a dict with per-sample _fwd.fastq and _rev.fastq files as a list
    value under the generic sample name as key
    """
    samples_fastqs = {}
    # get the fwd/rev fasta files per sample
    files =  glob.glob('%s/*fastq' % fastin)
    for f in sorted(files):
        fsplit = f.split('_')
        if fsplit[-1] in ['rev.fastq', 'fwd.fastq']:
            samples_fastqs.setdefault('_'.join(fsplit[:-1]), []).append(f)
    # invert the combito-sample to samples-to-combi dict
    design_rev = rev_dict(design)
    # add the sample name for each fwd/rev fastqs list
    add_to_dict_list(samples_fastqs, design_rev)
    return samples_fastqs


def get_mistags_files(fastin, primers_rad):
    fs =  glob.glob('%s/*fastq' % fastin)
    mistag_fastqs = [f for f in fs if re.search(r'_mistag_R[1-2]\.fastq$', f)]
    if len(mistag_fastqs) == 2:
        return mistag_fastqs
        # get unexpected samples information
    print('No pair of mistag files detected.\nExiting...')
    sys.exit()


def get_output_fp(args):
    f = args['i'][0]
    if 'o' in args:
        oFas = args['o']
        oRad = '%s_mistagFilt' % os.path.splitext(oFasta)[0]
    else:
        curdir = os.path.split(os.path.abspath(f))
        if os.path.isdir(f):
            oRad = '%s/%s_mistagFilt' % (os.path.abspath(f), curdir[1])
        else:
            oRad = '%s_mistagFilt' % os.path.splitext(os.path.abspath(f))[0]
        oFas = '%s.fasta' % oRad
    return oRad, oFas


def get_saturation(design, primers):
    n_F = len(primers['F'])
    n_R = len(primers['R'])
    FR = float(n_F * n_R)
    n_S = float(len(design))
    sat_frac = (n_S / (FR))
    sat_perc = round((sat_frac * 100), 2)
    return FR, sat_perc

