import time
import sys, os
import pandas as pd
import numpy as np

from mistag_getter import get_saturation


def make_outputs(filt, oFilt, ost, exp, nonCritic, orthos, design, primers, stats_merging):
    # kept_removed = {(f,r): [removed, kept]}
    output_filtering_results(filt, kept_removed, design, ost,
                             total_expect_reads, total_unexpect_reads, reads0,
                             reads1, stats_merging)
    output_nonCritical_data(expected, nonCritic, orthos, primers, design,
                            ost, total_expect_reads, total_unexpect_reads)


def get_pandas_tables(oFilt, expected, nonCritic):
    filt_table = pd.read_csv(oFilt, header=None, sep = '\t',
                                names=['sample', 'combi', 'isu',
                                       'reads', 'status'])
    expected_d = dict([k, pd.Series([int(v[0]) for v in Vs.values()])] for k,Vs in
                       expected.items())
    expected_df = pd.concat(expected_d, axis=1)
    nonCritic_d = dict([k, pd.Series([int(v[0]) for v in Vs.values()])] for k,Vs in
                        nonCritic.items())
    nonCritic_df = pd.concat(nonCritic_d, axis=1)
    return filt_table, expected_df, nonCritic_df


def write_output(args, oFilt, oFas, expected, nonCritic, orthos, design,
                  primers, stats_merging):
    # write the first lines of the output stats file
    ost = write_stats_header(args, oFas, design, primers)
    # get the pandas dataframe
    filt_table, expected_df, nonCritic_df = get_pandas_tables(oFilt, expected,
                                                              nonCritic)

    # calculate basic statistics (tot_exp, tot_nCr, total)
    basic_stats = calculate_mistag_stats(filt_table, nonCritic_df)
    exp_d = write_designs(ost, design, primers, filt_table, nonCritic_df,
                          basic_stats)
    write_filtering(ost, filt_table, stats_merging, basic_stats, exp_d)
    write_non_critical(ost, nonCritic_df, expected_df, basic_stats, primers,
                       design, orthos)
    ost.close()



def write_non_critical(ost, nCr_df, exp_df, basic_stats, primers, design,
                       orthos):
    Fs = sorted(primers['F'])
    Rs = sorted(primers['R'])
    n_combis = len(Fs) * len(Rs)
    n_nCr_combis = len(nCr_df.columns.values)
    tot_exp, tot_nCr, total = basic_stats
    ost.write('# Non-critical mistags / Unexpected samples\n')
    ost.write('>Total reads in unexpected samples:\t%s\t%s %% (of all %s \
reads)\n' % (tot_nCr, round((tot_nCr/float(total))*100, 2), total))
    ost.write('>tUnexpected samples:\t%s\t%s %% of possible unexpected\n' %
              (n_nCr_combis,
               round((n_nCr_combis/(float(n_combis)-len(design)))*100, 2)))
    unexp_per_exp = [len(v) for k, v in sorted(orthos.items())]
    unexp_per_exp = [len(v) for k, v in sorted(orthos.items())]
    ost.write('>Unexpected per expected:\t%s\n' % ','.join(map(str,
                                                               unexp_per_exp)))
    ost.write('>Mean(Unexpected per expected):\t%s\n' %
              round(np.mean(unexp_per_exp),2))
    ost.write('>Standard_deviation(Unexpected per expected):\t%s\n' %
              round(np.std(unexp_per_exp), 2))
    ost.write('>Range(Unexpected per expected):\t%s - %s\n' %
              (min(unexp_per_exp), max(unexp_per_exp)))
    ost.write('>Unexpected samples per expected sample\n')
    ost.write('Forward\tReverse\tExpected sample\tNumber of unexpected \
samples\tUnexpected samples (same forward)\tUnexpected samples \
(same reverse)\n')
    for combi in sorted(exp_df.columns.values):
        ortho = orthos[combi]
        F, R = combi
        sam = design[combi]
        sameF = ','.join(sorted(['+'.join(x) for x in sorted(ortho) if
                                 x[0]==combi[0]]))
        sameR = ','.join(sorted(['+'.join(x) for x in sorted(ortho) if
                                 x[1]==combi[1]]))

        ost.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (F, R, sam, len(ortho), sameF,
                                                sameR))
    ost.write('>Non-critical mistags per unexpected sample\n')
    ost.write('Forward\tReverse\tISUs\tReads\tMean(Reads per \
ISU)\tStd(Reads per ISU)\tMin(Reads per ISU)\tMax(reads per ISU)\n')
    nCr_df_stats = nCr_df.agg(['sum', 'count', 'mean', 'std', 'min', 'max']).T
    for combi in nCr_df.columns.values:
        F, R = combi
        row = map(str, nCr_df_stats.loc[combi,:])
        ost.write('%s\t%s\t%s\n' % (F,  R, '\t'.join(row)))


def write_stats_header(args, oFas, design, primers):
    oStats = oFas.replace('.fasta', '_stats.tsv')
    o = open(oStats, 'w')
    o.write('# script: %s\n' % os.path.abspath(sys.argv[0]))
    o.write('# github: https://github.com/FranckLejzerowicz/mistag_filter\n')
    curTime = time.strftime("%D %H:%M:%S", time.localtime())
    o.write('# date: %s\n' % curTime)
    o.write('# input fasta: %s\n' % args['i'][0])
    o.write('# input design: %s\n' % args['d'][0])
    o.write('# output fasta: %s\n' % oFas)
    FR, sat_perc = get_saturation(design, primers)
    o.write('# Design\n')
    o.write('>Saturation:\t%s %%\n' % sat_perc)
    o.write('>%s samples in design\n' % len(design))
    o.write('>%s possible primers combinations (samples)\n' % FR)
    return o


def write_fastas(folder, derep, stats_counts):
    curdir = os.path.abspath(folder).split('/')[-1]
    fastout = '%s/%s_mistagFormat.fasta' % (folder, curdir)
    with open(fastout, 'w') as o:
        for idx, seq_combis in enumerate(derep.items()):
            seq = seq_combis[0]
            combis = seq_combis[1]
            for combi, n in combis.items():
                o.write('>seq%s;size=%s;fwd=%s;rev=%s\n%s\n' % (idx, n,
                                                                combi[0],
                                                                combi[1],
                                                                seq))
    return fastout


def write_merging_cmd(out, cmd):
    if os.path.isfile(out):
        with open(out, 'a') as o:
            o.write('%s\n' % ' '.join(cmd))
    else:
        with open(out, 'w') as o:
            o.write('%s\n' % ' '.join(cmd))


def calculate_mistag_stats(tab, nCr):
    tot_exp = sum(tab.groupby(['status', 'combi']).agg('reads').sum())
    tot_nCr = sum(nCr.sum())
    total = float(tot_nCr + tot_exp)
    return tot_exp, tot_nCr, total


def write_designs(ost, design, primers, tab, nCr, basic_stats):
    tab_sums = tab.groupby(['status', 'combi']).agg('reads').sum()
    tab_counts= tab.groupby(['status', 'combi']).agg(np.count_nonzero)['reads']
    Fs = sorted(primers['F'])
    Rs = sorted(primers['R'])
    d = {}
    tot_exp, tot_nCr, total = basic_stats
    types = ['Design', 'Non-critical mistags (% of total)', 'Critical mistags (%)']
    for typdx, typ in enumerate(types):
        ost.write('>%s\n' % typ)
        ost.write('\t%s\n' % '\t'.join(Fs))
        for r in Rs:
            # matrix representation of the samples
            if typdx == 0:
                matData = ['X' if (f,r) in design else '' for f in Fs]
            # matrix representation of the detectable mistags
            elif typdx == 1:
                matData = []
                for f in Fs:
                    if (f,r) in nCr.columns:
                        n_nCr = nCr.loc[:,(f,r)].sum()
                        p_nCr = round(float(n_nCr/total)*100, 5)
                        matData.append(p_nCr)
                    else:
                        matData.append('')
            # matrix representation of the detected cross-contaminations
            else:
                matData = []
                for f in Fs:
                    if (f,r) in design:
                        k = ('kept', '%s__%s' % (f,r))
                        rm = ('removed', '%s__%s' % (f,r))
                        nKept = fecth_df_info(tab_sums, k)
                        NKept = fecth_df_info(tab_counts, k)
                        nRemoved = fecth_df_info(tab_sums, rm)
                        NRemoved = fecth_df_info(tab_counts, rm)
                        nTotal = float(nKept + nRemoved)
                        if nTotal:
                            pKept = round((nKept/nTotal)*100, 2)
                            PKept = round((NKept/float(NKept+NRemoved))*100, 2)
                        else:
                            pKept = 0.00
                            PKept = 0.00
                        matData.append(pKept)
                        d[(f,r)] = [design[(f,r)], nTotal,
                                    nRemoved, NRemoved,
                                    nKept, NKept, pKept, PKept]
                    else:
                        matData.append('')
            ost.write('%s\t%s\n' % (r, '\t'.join(map(str, matData))))
    return d


def fecth_df_info(tab, k):
    if k in list(tab.index.values):
        return tab.loc[k]
    else:
        return 0


def write_filtering(ost, tab, stats_merging, basic_stats, exp_d):
    tot_exp, tot_nCr, total = basic_stats
    ost.write('# Filtering\n')
    ost.write('>Total reads in expected samples:\t%s\t%s %% (of all %s reads)\n'
              % (tot_exp, round(float(tot_exp/total), 2), total))
    tot_kept = sum([x[4] for x in exp_d.values()])
    tot_rm = sum([x[2] for x in exp_d.values()])
    perc_tot_kept = round((tot_kept/float(tot_exp))*100, 2)
    perc_tot_rm = round((tot_kept/float(tot_exp))*100, 2)
    ost.write('>Kept reads in expected samples:\t%s\t%s %%\n' % (tot_kept,
                                                                 perc_tot_kept))
    ost.write('>Removed reads from expected samples:\t%s\t%s %%\n' % (tot_rm,
                                                                      perc_tot_rm))
    ost.write('>Column output for expected samples\n')
    if len(stats_merging):
        ost.write('Forward\tReverse\tSample\tdtd_demultiplexed_reads\t\
merged_reads\tRemoved reads\tRemoved ISUs\tPass-filter reads\t\
Pass-filter ISUs\tPass-filter ISUs (%)\tPass-filter reads (%)\n')
    else:
        ost.write('Forward\tReverse\tSample\tdtd_demultiplexed_reads\tRemoved \
reads\tRemoved ISUs\tPass-filter reads\tPass-filter \
ISUs\tPass-filter ISUs (%)\tPass-filter reads (%)\n')
    for combi in sorted(exp_d):
        F, R = combi
        sam = exp_d[combi][0]
        if len(stats_merging):
            fil = [x for x in stats_merging if sam in x][0]
            ost.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (F, R, sam,
                                        stats_merging[fil][0],
                                        stats_merging[fil][1],
                                        '\t'.join(map(str, exp_d[combi][2:]))))
        else:
            ost.write('%s\t%s\t%s\n' % (F, R, '\t'.join(map(str,
                                                            exp_d[combi]))))

