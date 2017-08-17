#!/usr/bin/env python

import re, time
import argparse
import subprocess
from math import sqrt
from scipy import stats
from numpy import array, std, mean

__author__ = "Franck Lejzerowicz"
__copyright__ = "Copyright 2017, The Deep-Sea Microbiome Project"
__credits__ = ["Philippe Esling"]
__license__ = "GPL V3"
__version__ = "1.0"
__maintainer__ = "Franck Lejzerowicz"
__email__ = "franck.lejzerowicz@unige.ch"

def mistag_filter():
    """
    Filter the critical mistags as in Esling et al. 2015 (NAR, https://doi.org/10.1093/nar/gkv107)
    """
    parser=argparse.ArgumentParser()
    parser.add_argument('-i', nargs = 1, required = True, help = 'Input fasta file name (required)')
    parser.add_argument('-d', nargs = 1, required = True, help = 'Multiplexing design file name (required)')
    parser.add_argument('-o', nargs = '?', default = argparse.SUPPRESS, help = "Output fasta file (default = input appended with 'mistagFiltered.txt')")
    parser.add_argument('-a', nargs = '?', default = 0.05, type = float, metavar='float between 0 and 1, max. 3 decimals', choices = [float(x)/1000 for x in range(1,1001)], help = "Alpha level for finding the Student's T critical value for the modified Thompson Tau rejection region calulation (default = 0.05)")
    parser.add_argument('-s', nargs = '?', default = ',', help = "Field separator in multiplexing design file (default = ',')")
    parser.add_argument('--out', action = 'store_true', default = False, help = "Leave expected sample sequences out of non-critical mistags distribution for calculations of the rejection region (default = not active)")
    parse=parser.parse_args()
    args=vars(parse)

    # parse design and collect all primers
    designFilin = args['d'][0]
    sep = args['s']
    design = get_design(designFilin, sep)
    primers = get_primers(design)

    # parse input fasta and collect the abundance distributions
    fastaFilin = args['i'][0]
    #   - for sequences in expected samples
    #   - for sequences in unexpected samples (non-critical mistags)
    expected, nonCritic = parse_input(fastaFilin, design, primers)
    # unexpected samples sharing one tagged primer with each expected sample
    ortho_samples = sample_orthogonals(design, expected, nonCritic)

    # mistag filter main
    alpha_level = args['a']
    notIn = args['out']
    filtered = filter_mistags(expected, nonCritic, ortho_samples, alpha_level, notIn)

    # outputs
    if args.has_key('o'):
        outputFasta = args['o']
        outputStats = args['o'] + '.stats'
    else:
        outputRad = get_output_filename(fastaFilin)
        outputFasta = outputRad + '.txt'
        outputStats = outputRad + '.stats'
    make_outputs(outputFasta, outputStats, filtered, design)

#def make_outputs(outFasta, outStats, filt, design):


def filter_mistags(expected, nonCritic, ortho_samples, alpha_level, notIn):
    """
    Mistag filter
    For each sample;
        get the othogonal samples (sharing the forward / the reverse primer)
        For each ISU;
            compute the distribution moments of this ISU abundance in these samples group
            compare the deviation of the current ISU with the rejection region of the Thompson Tau test
            (see e.g. http://www.statisticshowto.com/modified-thompson-tau-test/)
    Returns the dict with kept (key "1") or removed (key "0") ISU lists for each sample (nested key)
    """
    filt = {0: {}, 1: {}}
    for sample in expected:
        curOrtho_samples = ortho_samples[sample]
        for ISU in expected[sample]:
            curAbund = int(expected[sample][ISU][0])
            orthoAbunds = [int(nonCritic[x][ISU][0]) for x in curOrtho_samples if nonCritic[x].has_key(ISU)]
            if notIn:
                totalAbund = orthoAbunds
            else:
                totalAbund = [curAbund] + orthoAbunds
            # get the properties of the distribution
            meanAbund = mean(totalAbund)
            devAbund = std(totalAbund)
            n = len(totalAbund)
            if n == 1:
                filt[1].setdefault(sample,[]).append(ISU)
                continue
            df = (n-1)
            curDev = abs(meanAbund - curAbund)
            # get the critical student's t-test value
            tAlpha = stats.t.isf(alpha_level, df)
            # rejection region
            r = ((tAlpha*(n-1))/(sqrt(n)*sqrt(n-1+(tAlpha**2))))*devAbund
            if curDev <= r:
                filt[0].setdefault(sample,[]).append(ISU)
            else:
                filt[1].setdefault(sample,[]).append(ISU)
    return filt

def sample_orthogonals(design, expected, nonCritic):
    """
    For each expected primer combination of the original design (dict key)
    Return as value the list of the primer combinations that are not in the original design and that share the same forward or reverse primers
    """
    ortho = {}
    for combi in expected:
        curF = combi[0]
        curR = combi[1]
        same_fwd_idx = [x for x in nonCritic.keys() if x[0] == curF]
        same_rev_idx = [x for x in nonCritic.keys() if x[1] == curR]
        ortho[combi] = same_fwd_idx + same_rev_idx
    return ortho

def get_output_filename(filin):
    """
    Return automatic output file name based on input file name
    """
    filinStrip = filin.rstrip('txt')
    if filin == filinStrip:
        return '%smistagFiltered' % filin
    else:
        return '%s.mistagFiltered' % filinStrip

def get_primers(design):
    """
    Get all the forward and all the reverse primer name in a dict
    """
    primers = {}
    for frx, FR in enumerate(['F','R']):
        primers[FR] = {}
        for primer in list(set([x[frx] for x in design.keys()])):
            primers[FR][primer] = 1
    return primers

def add_key(d, pair, seq, n, seqID):
    """
    Fill a nested dict with given info
    Here for registering the sequences that are associated with an expected primer combinations, and with non-critical combinations
    """
    if d.has_key(pair):
        if d[pair].has_key(seq):
            print 'Warning: sequence "%s" already encountered for the current primer combination "%s"' % (seq, ' + '.join(pair))
            d[pair][seq] = [n, seqID]
        else:
            d[pair][seq] = [n, seqID]
    else:
        d[pair] = {seq: [n, seqID]}

def parse_input(fastaFilin, design, primers):
    """
    Read the fasta input file provided by the pipeline developed by Yoann Dufresne
    Format:
            >SEQUENCE_ID;size=NUMBER_OF_READ_AFTER_DEREPLICATION;for=FORWARD_PRIMER_NAME;rev=REVERSE_PRIMER_NAME
            TACTAGTCGTGTAGCTGATCGTAGCTAGCTAGCTAGCTAGCTGTAGATCGTAGCTGACTGAC
    Returns two dicts with [abundance, sequence ID] values for each sequence key, nested under each sample key
    One dict called "expected" for the sequences belonging to an expected samples of the original design
    One dict called "nonCritic" for the sequences associated with to a unexpected primer combination
    """
    expected = {}
    nonCritic = {}
    start = 0
    curSeq = {}
    seq = ''
    for line in open(fastaFilin, 'rU'):
        if line[0] == '>':
            if len(seq):
                # update the dicts for expected tage combis and non-critical mistag combis
                if design.has_key(pair):
                    add_key(expected, pair, seq, n, seqID)
                elif primers['F'].has_key(f) and primers['R'].has_key(r):
                    add_key(nonCritic, pair, seq, n, seqID)
                curSeq = {}
                seq = ''
                start = 1
            splitID = line.strip()[1:].split(';')
            seqID = splitID[0]
            for field in splitID[1:]:
                if field.startswith('size='):
                    n = re.search('\d+', field).group(0)
                elif field.startswith('fwd=') or field.startswith('for='):
                    f = field.split('=')[-1].strip()
                elif field.startswith('rev=') or field.startswith('rv'):
                    r = field.split('=')[-1].strip()
            pair = tuple([f,r])
        else:
            seq += line.strip()
    return expected, nonCritic

def get_design(designFile, sep):
    """
    Rarse the deign file and collect the info about which primer combinations are expected
    Returns a dict with each expected primer combination as key (value set to sample name for output)
    """
    design = {}
    for ldx, line in enumerate(open(designFile)):
        splitLine = line.strip().split(sep)
        if ldx:
            f = splitLine[fields['f']]
            r = splitLine[fields['r']]
            s = splitLine[fields['s']]
            design[tuple([f,r])] = s
        else:
            fields = {}
            ranks = splitLine
            for rdx, rank in enumerate(ranks):
                if rank.startswith('sam'):
                    fields['s'] = rdx
                elif rank.startswith('for'):
                    fields['f'] = rdx
                elif rank.startswith('rev'):
                    fields['r'] = rdx
    return design

mistag_filter()
