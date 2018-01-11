#!/usr/bin/env python

import glob
import re, time
import os, sys
import argparse
import subprocess
from math import sqrt
from scipy import stats
from numpy import array, std, mean
from multiprocessing import Process, Manager, cpu_count, current_process
import cProfile

__author__ = "Franck Lejzerowicz"
__copyright__ = "Copyright 2017, The Deep-Sea Microbiome Project"
__credits__ = ["Philippe Esling"]
__license__ = "GPL V3"
__version__ = "1.0"
__maintainer__ = "Franck Lejzerowicz"
__email__ = "franck.lejzerowicz@unige.ch"


def get_args():
    """
    Filter the critical mistags as in Esling et al. 2015 (NAR, https://doi.org/10.1093/nar/gkv107)
    """
    parser=argparse.ArgumentParser()
    parser.add_argument('-i', nargs = 1, required = True, help = 'Input fasta file name [or folder if pre-formatting needed] (required)')
    parser.add_argument('-d', nargs = 1, required = True, help = 'Multiplexing design file name (required)')
    parser.add_argument('-o', nargs = '?', default = argparse.SUPPRESS, help = "Output fasta file (default = input appended with 'mistagFiltered.fasta')")
    parser.add_argument('-a', nargs = '?', default = 0.05, type = float, metavar='float between 0 and 1, max. 3 decimals', choices = [float(x)/1000 for x in range(1,1001)], help = "Alpha level for finding the Student's T critical value for the modified Thompson Tau rejection region calulation (default = 0.05)")
    parser.add_argument('--out', action = 'store_true', default = False, help = "Leave expected sample sequences out of non-critical mistags distribution for calculations of the rejection region (default = not active)")
    parser.add_argument('-m', nargs = '?', choices = ['pandaseq', 'vsearch'], default = 'pandaseq', help = 'Reads merging software - which must be installed and running (default = pandaseq)')
    parse=parser.parse_args()
    args=vars(parse)
    return args


def rev_dict(d):
    D = {}
    for k,v in d.items():
        D[v] = k
    return D


def get_inputs_to_format(fastaFolder, design, primers_rad):
    design_rev = rev_dict(design)
    files = {}
    folderFiles =  glob.glob('%s/*fastq' % fastaFolder)
    for f in folderFiles:
        fsplit = f.split('_')
        if fsplit[-1] in ['rev.fastq', 'fwd.fastq']:
            files.setdefault('_'.join(fsplit[:-1]), []).append(f)
    for f in files:
        for s in design_rev:
            if s in f:
                files[f].append(design_rev[s])
                break
    mistag_files = [x for x in folderFiles if re.search(r'_mistag_R[1-2]\.fastq$', x)]
    if len(mistag_files) == 2:
        mistag_IDs, mistag_all = get_mistags_unexpected(mistag_files, primers_rad)
    return files, [mistag_files, mistag_IDs, mistag_all]


def get_mistags_unexpected(mistags, primers_rad):
    d = {}
    D = {'weird': {}, 'unexpected': {}}
    with open(mistags[0]) as f1, open(mistags[1]) as f2:
        idx = 0
        for id1, id2 in zip(f1, f2):
            idx += 1
            if idx == 1:
                t1 = id1.strip().split(';tag:')[-1]
                t2 = id2.strip().split(';tag:')[-1]
                combi = tuple([t1,t2])
                if 'unknown' in combi:
                    continue
                typ = 'weird'
                if t1.split('-')[0] != t2.split('-')[0]:
                    combi = tuple([[x for x in combi if x.split('-')[0] == primers_rad['F']][0], [x for x in combi if x.split('-')[0] == primers_rad['R']][0]])
                    d[id1.split()[0][1:]] = combi
                    typ = 'unexpected'
                if D[typ].has_key(combi):
                    D[typ][combi] += 1
                else:
                    D[typ][combi] = 1
            if idx == 4:
                idx = 0
    return d, D


def make_outputs(expected, nonCritic, ortho_samples, designFilin, fastaFilin, outFasta, outStats, filtered, design, primers, stats_merging):
    ost = open(outStats, 'w')
    ost.write('# script: %s\n' % os.path.abspath(sys.argv[0]))
    ost.write('# github: https://github.com/FranckLejzerowicz/mistag_filter\n')
    curTime = time.strftime("%D %H:%M:%S", time.localtime())
    ost.write('# date: %s\n' % curTime)
    ost.write('# input fasta: %s\n' % fastaFilin)
    ost.write('# input design: %s\n' % designFilin)
    ost.write('# output fasta: %s\n' % outFasta)
    ost.write('# Design\n')
    ost.write('>Saturation:\t%s %%\t(%s samples, %s possible combinations)\n' % (round((float(len(design))/(len(primers['F'])*len(primers['R'])))*100, 2), len(design), len(primers['F'])*len(primers['R'])))
    # filtering data
    reads1_per_ISU = [[y[-1] for y in x] for x in filtered[1].values()]
    reads0_per_ISU = [[y[-1] for y in x] for x in filtered[0].values()]
    reads1 = sum(sum(reads1_per_ISU, []))
    reads0 = sum(sum(reads0_per_ISU, []))
    # samples data
    total_expect_reads = float(reads0 + reads1)
    # non-critical data
    total_unexpect_reads = sum([sum([float(val[0]) for val in nonCritic[unexp].values()]) for unexp in nonCritic])
    kept_removed = output_fasta_and_matrices(filtered, design, outFasta, ost, primers)
    output_filtering_results(filtered, kept_removed, design, ost, total_expect_reads, total_unexpect_reads, reads0, reads1, stats_merging)
    output_nonCritical_data(expected, nonCritic, ortho_samples, primers, design, ost, total_expect_reads, total_unexpect_reads)
    ost.close()

def output_nonCritical_data(expected, nonCritic, ortho_samples, primers, design, ost, total_expect_reads, total_unexpect_reads):
    ost.write('# Non-critical mistags / Unexpected samples\n')
    ost.write('>Total reads in unexpected samples:\t%s\t%s %% (of all %s reads)\n' % (total_unexpect_reads, round((total_unexpect_reads/(total_expect_reads+total_unexpect_reads))*100, 2), total_expect_reads+total_unexpect_reads))
    ost.write('>Unexpected samples:\t%s\t%s %% of possible unexpected\n' % (len(nonCritic), round(((len(nonCritic)/float((len(primers['F'])*len(primers['R']))-len(design)))*100), 2)))
    n_unexpected_per_expected = [len(val) for key,val in sorted(ortho_samples.items())]
    ost.write('>Unexpected per expected:\t%s\n' % ','.join(map(str, n_unexpected_per_expected)))
    ost.write('>Mean(Unexpected per expected):\t%s\n' % round(mean(n_unexpected_per_expected),2))
    ost.write('>Standard_deviation(Unexpected per expected):\t%s\n' % round(std(n_unexpected_per_expected), 2))
    ost.write('>Range(Unexpected per expected):\t%s - %s\n' % (min(n_unexpected_per_expected), max(n_unexpected_per_expected)))
    ost.write('>Unexpected samples per expected sample\n')
    ost.write('Forward\tReverse\tExpected sample\tNumber of unexpected samples\tUnexpected samples (same forward)\tUnexpected samples (same reverse)\n')
    for combi in sorted(expected):
        curOrthos = ortho_samples[combi]
        ost.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (combi[0], combi[1], design[combi], len(curOrthos),
                                                ','.join(sorted(['+'.join(x) for x in sorted(curOrthos) if x[0]==combi[0]])),
                                                ','.join(sorted(['+'.join(x) for x in sorted(curOrthos) if x[1]==combi[1]]))))
    ost.write('>Non-critical mistags per unexpected sample\n')
    ost.write('Forward\tReverse\tISUs\tReads\tMean(Reads per ISU)\tStd(Reads per ISU)\tMin(Reads per ISU)\tMax(reads per ISU)\n')
    for combi, isus in sorted(nonCritic.items()):
        ost.write('%s\t%s' % (combi[0], combi[1]))
        non_critics_reads_per_ISU = [float(x[0]) for x in isus.values() if float(x[0])]
        non_critics_reads = sum(non_critics_reads_per_ISU)
        ost.write('\t%s\t%s\t%s\t%s\t%s\t%s\n' % (len(non_critics_reads_per_ISU),
                                  non_critics_reads,
                                  mean(non_critics_reads_per_ISU),
                                  std(non_critics_reads_per_ISU),
                                  min(non_critics_reads_per_ISU),
                                  max(non_critics_reads_per_ISU)))

def output_filtering_results(filtered, kept_removed, design, ost, total_expect_reads, total_unexpect_reads, reads0, reads1, stats_merging):
    ost.write('# Filtering\n')
    ost.write('>Total reads in expected samples:\t%s\t%s %% (of all %s reads)\n' % (total_expect_reads, round((total_expect_reads/(total_expect_reads+total_unexpect_reads))*100, 2), total_expect_reads+total_unexpect_reads))
    ost.write('>Kept reads in expected samples:\t%s\t%s %%\n' % (reads1, round((reads1/total_expect_reads)*100, 2)))
    ost.write('>Removed reads from expected samples:\t%s\t%s %%\n' % (reads0, round((reads0/total_expect_reads)*100, 2)))
    ost.write('>Column output for expected samples\n')
    if len(stats_merging):
        ost.write('Forward\tReverse\tSample\tdtd_demultiplexed_reads\tmerged_reads\tRemoved reads\tRemoved ISUs\tPass-filter reads\tPass-filter ISUs\tPass-filter (%)\n')
    else:
        ost.write('Forward\tReverse\tSample\tRemoved reads\tRemoved ISUs\tPass-filter reads\tPass-filter ISUs\tPass-filter (%)\n')
    for combi in sorted(kept_removed):
        pass_in = [x for x in kept_removed[combi][1] if x]
        pass_out = [x for x in kept_removed[combi][0] if x]
        sample = design[combi]
        fileRad = [x for x in stats_merging.keys() if sample in x][0]
        if len(stats_merging):
            ost.write('%s\n' % '\t'.join(map(str, [combi[0], combi[1], sample,
                                                stats_merging[fileRad][0],
                                                stats_merging[fileRad][1],
                                                sum(pass_out), len(pass_out),
                                                sum(pass_in), len(pass_in),
                                                (sum(pass_in)/(sum(pass_in)+sum(pass_out)))*100])))
        else:
            ost.write('%s\n' % '\t'.join(map(str, [combi[0], combi[1], sample,
                                                sum(pass_out), len(pass_out),
                                                sum(pass_in), len(pass_in),
                                                (sum(pass_in)/(sum(pass_in)+sum(pass_out)))*100])))

def output_fasta_and_matrices(filtered, design, outFasta, ost, primers):
    ofas = open(outFasta, 'w')
    kept_removed = {}
    for typdx, typ in enumerate(['Design', 'Critical mistags (%)']):
        for combi in filtered[1]:
            sample = design[combi]
            F,R = combi
            check_duplicates(filtered[1][combi], combi, sample)
            for sdx, sequence in enumerate(filtered[1][combi]):
                seq_1,seqID_1,reads_1 = tuple(sequence)
                seq_0,seqID_0,reads_0 = tuple(filtered[0][combi][sdx])
                if typdx == 0 and float(reads_1):
                    ofas.write('>%s_%s;size=%s;for=%s;rev=%s\n%s\n' % (seqID_1, sample, reads_1, F, R, seq_1))
        ost.write('>%s\n' % typ)
        ost.write('\t%s\n' % '\t'.join(sorted(primers['F'])))
        for r in sorted(primers['R']):
            if typdx == 0:
                matData = ['1' if design.has_key((f,r)) else '' for f in sorted(primers['F'])]
            else:
                matData = []
                for f in sorted(primers['F']):
                    if design.has_key((f,r)):
                        removed, kept = tuple([[float(x[-1]) for x in filtered[kerm][(f,r)]] for kerm in [0, 1]])
                        percentKept = round(sum(removed)/(sum(kept)+float(sum(removed)))*100, 2)
                        kept_removed[(f,r)] = [removed, kept]
                        matData.append(percentKept)
                    else:
                        matData.append('')
            ost.write('%s\t%s\n' % (r, '\t'.join(map(str, matData))))
    ofas.close()
    return kept_removed

def check_duplicates(liste, combi, sample):
    """
    Check if one of the sample unique sequences has been found twice --- should be useless
    """
    if len(set([x[0] for x in liste])) != len(liste):
        print 'Error: duplicate sequences found for the sample %s (%s)' % (sample, ' + '.join(combi))
        print '  - %s' % '\n'.join([x for x in sorted([x for x in set(l) if liste.count(x)>1])])
        print 'Exiting'
        return -1

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
    # for each sample
    for sample in expected:
        # get the orthogonal samples from which to compute the abundance distributions
        curOrtho_samples = ortho_samples[sample]
        # for each ISU in the expected samples
        for ISU in expected[sample]:
            isuName = str(expected[sample][ISU][1])
            curAbund = int(expected[sample][ISU][0])
            orthoAbunds = [int(nonCritic[x][ISU][0]) for x in curOrtho_samples if nonCritic[x].has_key(ISU)]
            if notIn:
                totalAbund = orthoAbunds
            else:
                totalAbund = [curAbund] + orthoAbunds
            # get the number of abundances in the distribution
            n = len(totalAbund)
            # keep the ISU if only one abundance value (i.e. may be rare but not a cross-conta)
            if n <= 1:
                filt[1].setdefault(sample,[]).append((ISU, isuName, curAbund))
                filt[0].setdefault(sample,[]).append((ISU, isuName, 0))
                continue
            # get distribution paramters
            meanAbund = mean(totalAbund)
            devAbund = std(totalAbund)
            df = (n-1)
            curDev = abs(meanAbund - curAbund)
            # get the critical student's t-test value
            tAlpha = stats.t.isf(alpha_level, df)
            # rejection region
            r = ((tAlpha*(n-1))/(sqrt(n)*sqrt(n-1+(tAlpha**2))))*devAbund
            # remove ISU if its abundance is in the rejection region
            if curDev <= r:
                filt[1].setdefault(sample,[]).append((ISU, isuName, 0))         # keep 0 reads of the current ISU
                filt[0].setdefault(sample,[]).append((ISU, isuName, curAbund))  # remove the current abundance of the current ISU
            # keep ISU if its abundance is in the rejection region
            else:
                filt[1].setdefault(sample,[]).append((ISU, isuName, curAbund))  # keep the current abundance of the current ISU
                filt[0].setdefault(sample,[]).append((ISU, isuName, 0))         # remove 0 reads of the current ISU
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


def get_primers(design):
    """
    Get all the forward and all the reverse primer name in a dict
    """
    primers = {}
    for frx, FR in enumerate(['F', 'R']):
        primers[FR] = {}
        for primer in list(set([x[frx] for x in design.keys()])):
            primers[FR][primer] = 1
    return primers

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
    curSeq = {}
    seq = ''
    for lindx,line in enumerate(open(fastaFilin, 'rU')):
        if lindx and lindx % 10 == 0:
            print '\r', 'Number of parsed per-sample ISUs:', lindx,
            sys.stdout.flush()
        if line[0] == '>':
            if len(seq):
                # update the dicts for expected tage combis and non-critical mistag combis
                update_dict(pair, expected, nonCritic, seq, n, seqID, primers, f, r, design)
                curSeq = {}
                seq = ''
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
    update_dict(pair, expected, nonCritic, seq, n, seqID, primers, f, r, design)
    return expected, nonCritic

def update_dict(pair, expected, nonCritic, seq, n, seqID, primers, f, r, design):
    """
    Add the recorded sequence data into either the expected or nonCritic
    """
    if design.has_key(pair):
        add_key(expected, pair, seq, n, seqID)
    elif primers['F'].has_key(f) and primers['R'].has_key(r):
        add_key(nonCritic, pair, seq, n, seqID)

def add_key(d, pair, seq, n, seqID):
    """
    Fill a nested dict with given info
    Register sequences associated with:
        - expected primer combinations (first call in update_dict)
        - non-critical combinations (second call in update_dict)
    """
    if d.has_key(pair):
        if d[pair].has_key(seq):
            print 'Error: sequence "%s" already encountered (with size=%s) for the current primer combination "%s" (now with size=%s)\nExiting' % (seq, d[pair][seq][0], ' + '.join(pair), n)
            sys.exit()
        else:
            d[pair][seq] = [n, seqID]
    else:
        d[pair] = {seq: [n, seqID]}

def get_design(designFile):
    """
    Rarse the design file and collect the info about which primer combinations are expected
    Returns a dict with each expected primer combination as key (value set to sample name for output)
    """
    design = {}
    fields = ['run', 'forward', 'reverse', 'sample']
    for ldx, line in enumerate(open(designFile)):
        splitLine = line.strip().split(',')
        # get the design dict
        if ldx:
            if check_fields(fields, info, designFile):
                r,F,R,s = tuple([splitLine[info[x]] for x in fields])
                design[tuple([F,R])] = s
        # get the headers indices
        else:
            info = {}
            ranks = splitLine
            for field in fields:
                info[field] = ranks.index(field)
    return design

def check_fields(fields, info, filin):
    """
    Check if the minimum required info is present in the design file:
    """
    for field in fields:
        if info.has_key(field)==False:
            print 'No field "%s" in %s\nExiting' % (field, filin)
            sys.exit()
    else:
        return 1



def get_merging_config(soft, folder):
    merging_config_fp = '%s/merging_%s.conf' % (folder, soft)
    merging = {}
    if os.path.isfile(merging_config_fp) == False:
        o=open(merging_config_fp, 'w')
        if soft == 'vsearch':
            o.write('min_qual,--fastq_qmin,0\n')
            o.write('max_qual,--fastq_qmax,41\n')
            o.write('min_overlap,--fastq_minovlen,15\n')
            o.write('max_diffs,--fastq_maxdiffs,5\n')
            o.write('min_size,--fastq_minmergelen,1\n')
            o.write('max_size,--fastq_maxmergelen,1000\n')
            o.write('threads,--threads,%s\n' % cpu_count())
        elif soft == 'pandaseq':
            o.write('thresh,-t,0.6\n')
            o.write('algo,-A,pear\n')
            o.write('threads,-T,%s\n' % cpu_count())
            o.write('min_len,-l,\n')
            o.write('max_len,-L,\n')
            o.write('min_overap,-o,30\n')
            o.write('comment,-d,bfsrk\n')
        o.close()
    if soft == 'vsearch':
        merging['for'] = ['--fastq_mergepairs', None]
        merging['rev'] = ['--reverse', None]
        merging['out'] = ['--fastaout', None]
    elif soft == 'pandaseq':
        merging['for'] = ['-f', None]
        merging['rev'] = ['-r', None]
        merging['out'] = ['-w', None]
        merging['log'] = ['-G', None]
    with open(merging_config_fp) as f:
        for line in f:
            line_split = line.strip().split(',')
            merging[line_split[0]] = line_split[1:]
    return merging


def run_merging(fastqs, merging_soft, merging_config, return_dict):
    curProc_name = current_process().name
    cmd = [merging_soft]
    for config_opt, config_val in merging_config.items():
        if config_val[-1]:
            cmd.extend(config_val)
        else:
            if config_opt == 'for':
                cmd.extend([config_val[0], fastqs[0]])
            if config_opt == 'rev':
                cmd.extend([config_val[0], fastqs[1]])
            if config_opt == 'log':
                cmd.extend([config_val[0], '%s_log.bz2' % os.path.splitext(fastqs[0])[0]])
            if config_opt == 'out':
                out = '%s.fasta' % '_'.join(fastqs[0].split('_')[:-1])
                return_dict[out] = [curProc_name]
                cmd.extend([config_val[0], out])
    subprocess.call(cmd)


def add_seq_combi_keys(derep, seq, combi):
    if derep.has_key(seq):
        if derep[seq].has_key(combi):
            derep[seq][combi] += 1
        else:
            derep[seq][combi] = 1
    else:
        derep[seq] = {combi: 1}


def format_fastas(folder, return_dict, samples_fastqs, mistag_unexpected, stats_counts):
    ## perform counting HERE!!! no grep
    derep = {}
    for fasta in return_dict.keys():
        if fasta.endswith('mistag.fasta'):
            combi = 0
        else:
            combi = samples_fastqs[fasta.replace('.fasta', '')][-1]
        seq = ''
        c = 0
        with open(fasta) as f:
            for line in f:
                if line[0] == '>':
                    if len(seq):
                        if add:
                            c += 1
                            add_seq_combi_keys(derep, seq, combi)
                        seq = ''
                    add = 1
                    if fasta.endswith('_mistag.fasta'):
                        ID = ':'.join(line[1:].split(':')[:7])
                        if mistag_unexpected.has_key(ID):
                            combi = mistag_unexpected[ID]
                        else:
                            add = 0
                else:
                    seq += line.strip()
            if add:
                c += 1
                add_seq_combi_keys(derep, seq, combi)
            if fasta.endswith('_mistag.fasta') == False:
                stats_counts[fasta].append(str(c))
    fasta_out = '%s/%s_toMistagFilter.fasta' % (folder, folder)
    o = open(fasta_out, 'w')
    counter = {}
    for idx, seq_combis in enumerate(derep.items()):
        seq = seq_combis[0]
        combis = seq_combis[1]
        for combi, n in combis.items():
            o.write('>seq%s;size=%s;fwd=%s;rev=%s\n%s\n' % (idx, n, combi[0], combi[1], seq))
    o.close()
    return fasta_out


def get_merging_stats(return_dict, folder, stats_counts):
    with open(return_dict.keys()[0]) as f:
        for line in f:
            grep_key = '\:'.join(line[1:].split(':')[:3])
            break
    cmd = 'grep -cR "%s" %s/*_fwd.fastq' % (grep_key, folder)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    for fastq in out.strip().split('\n'):
        f = fastq.split(':')[0]
        n = fastq.split(':')[1]
        stats_counts['%s.fasta' % f.split('_fwd.fastq')[0]] = [n]


if __name__ == '__main__':
    args = get_args()
    # parse design and collect all primers
    designFilin = args['d'][0]
    design = get_design(designFilin)
    primers = get_primers(design)
    primers_rad = dict([x, y.keys()[0].split('-')[0]] for x,y in primers.items())

    stats_counts = {}
    fastaFilin = args['i'][0]
    # perform formatting - may include reads merging step
    if os.path.isdir(fastaFilin):
        fastaFilin = fastaFilin.rstrip('/')
        samples_fastqs, mistags = get_inputs_to_format(fastaFilin, design, primers_rad)
        mistag_files = mistags[0]
        mistag_unexpected = mistags[1]
        jobs_merging = []
        manager = Manager()
        return_dict = manager.dict()
        merging_soft = args['m']
        merging_config = get_merging_config(merging_soft, fastaFilin)
        for sample_fastq, sample_fastqs in samples_fastqs.items():
            p = Process(target = run_merging, args = (sample_fastqs, merging_soft, merging_config, return_dict,), name = 'merging_%s' % sample_fastq.split('/')[-1])
            jobs_merging.append(p)
            p.start()
        p = Process(target = run_merging, args = (mistag_files, merging_soft, merging_config, return_dict,), name = 'merging_mistags')
        jobs_merging.append(p)
        p.start()
        for job_merging in jobs_merging:
            job_merging.join()
        get_merging_stats(return_dict, fastaFilin, stats_counts)
        fastaFilin = format_fastas(fastaFilin, return_dict, samples_fastqs, mistag_unexpected, stats_counts)

    # parse input fasta and collect the abundance distributions
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
        outputRad = os.path.splitext(outputFasta)[0]
    else:
        outputRad = '%s_mistagFilt' % os.path.splitext(fastaFilin)[0]
        outputFasta = outputRad + '.fasta'
    outputStats = outputRad + '_stats.tsv'
    make_outputs(expected, nonCritic, ortho_samples, designFilin, fastaFilin, outputFasta, outputStats, filtered, design, primers, stats_counts)
    print 'Outputs:'
    print os.path.abspath(outputFasta)
    print os.path.abspath(outputStats)
