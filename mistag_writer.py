import time
import sys, os
import pandas as pd
import numpy as np

from mistag_getter import get_saturation


def get_pandas_tables(o_filt, expected, nonCritic):
    filt_table = pd.read_csv(
        o_filt,
        header=None,
        sep='\t',
        names=[
            'sample', 'combi', 'isu', 'reads', 'status'
        ]
    )

    expected_d = dict(
        (k, pd.Series([int(v[0]) for v in vs.values()]))
        for k, vs in expected.items()
    )
    expected_df = pd.concat(expected_d, axis=1)

    non_critic_d = dict(
        (k, pd.Series([int(v[0]) for v in vs.values()]))
        for k, vs in nonCritic.items()
    )
    non_critic_df = pd.concat(non_critic_d, axis=1)

    return filt_table, expected_df, non_critic_df


def write_output(args, o_filt, o_fas, expected, non_critic,
                 orthos, design, primers, stats_merging):
    # write the first lines of the output stats file
    ost = write_stats_header(args, o_fas, design, primers)
    # get the pandas dataframe
    filt_table, expected_df, non_critic_df = get_pandas_tables(
        o_filt, expected, non_critic
    )
    # calculate basic statistics (tot_exp, tot_nCr, total)
    basic_stats = calculate_mistag_stats(filt_table, non_critic_df)
    exp_design = write_designs(
        ost, design, primers, filt_table, non_critic_df, basic_stats
    )
    write_filtering(
        ost, stats_merging, basic_stats, exp_design
    )
    write_non_critical(
        ost, non_critic_df, expected_df, basic_stats, primers, design, orthos
    )
    ost.close()


def write_non_critical(
        ost, ncr_df, exp_df, basic_stats, primers, design, orthos):

    fs = sorted(primers['F'])
    rs = sorted(primers['R'])
    n_combis = len(fs) * len(rs)
    n_ncr_combis = len(ncr_df.columns.values)
    tot_exp, tot_ncr, total = basic_stats
    ost.write('# Non-critical mistags / Unexpected samples\n')
    ost.write('>Total reads in unexpected samples:\t%s\t%s %% (of all %s reads)\n' % (
        tot_ncr, round((tot_ncr/float(total))*100, 2), total
    ))
    ost.write('>tUnexpected samples:\t%s\t%s %% of possible unexpected\n' % (
        n_ncr_combis, round((n_ncr_combis/(float(n_combis)-len(design)))*100, 2)
    ))
    unexp_per_exp = [len(v) for k, v in sorted(orthos.items())]
    unexp_str = ','.join(map(str, unexp_per_exp))
    ost.write('>Unexpected per expected:\t%s\n' % unexp_str)
    unexp_mean = round(float(np.mean(unexp_per_exp)), 2)
    ost.write('>Mean(Unexpected per expected):\t%s\n' % unexp_mean)
    unexp_std = round(float(np.std(unexp_per_exp)), 2)
    ost.write('>Standard_deviation(Unexpected per expected):\t%s\n' % unexp_std)
    ost.write('>Range(Unexpected per expected):\t%s - %s\n' % (
        min(unexp_per_exp), max(unexp_per_exp)
    ))
    ost.write('>Unexpected samples per expected sample\n')
    ost.write('Forward\tReverse\tExpected sample\tNumber of unexpected \
samples\tUnexpected samples (same forward)\tUnexpected samples \
(same reverse)\n')
    for combi in sorted(exp_df.columns.values):
        ortho = orthos[combi]
        f, r = combi
        sam = design[combi]
        same_f = ','.join(
            sorted(['+'.join(x) for x in sorted(ortho) if x[0] == combi[0]])
        )
        same_r = ','.join(
            sorted(['+'.join(x) for x in sorted(ortho) if x[1] == combi[1]])
        )
        ost.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (
            f, r, sam, len(ortho), same_f, same_r
        ))

    ost.write('>Non-critical mistags per unexpected sample\n')
    ost.write('Forward\tReverse\tISUs\tReads\tMean(Reads per \
ISU)\tStd(Reads per ISU)\tMin(Reads per ISU)\tMax(reads per ISU)\n')
    ncr_df_stats = ncr_df.agg(['sum', 'count', 'mean', 'std', 'min', 'max']).T
    for combi in ncr_df.columns.values:
        f, r = combi
        row = map(str, ncr_df_stats.loc[combi, :])
        ost.write('%s\t%s\t%s\n' % (f, r, '\t'.join(row)))


def write_stats_header(args, o_fas, design, primers):
    o_stats = o_fas.replace('.fasta', '_stats.tsv')
    o = open(o_stats, 'w')
    o.write('# script: %s\n' % os.path.abspath(sys.argv[0]))
    o.write('# github: https://github.com/FranckLejzerowicz/mistag_filter\n')
    cur_time = time.strftime("%D %H:%M:%S", time.localtime())
    o.write('# date: %s\n' % cur_time)
    o.write('# input fasta: %s\n' % args['i'][0])
    o.write('# input design: %s\n' % args['d'][0])
    o.write('# output fasta: %s\n' % o_fas)
    fr, sat_percent = get_saturation(design, primers)
    o.write('# Design\n')
    o.write('>Saturation:\t%s %%\n' % sat_percent)
    o.write('>%s samples in design\n' % len(design))
    o.write('>%s possible primers combinations (samples)\n' % fr)
    return o


def write_fastas(folder, derep):
    curdir = os.path.abspath(folder).split('/')[-1]
    fastout = '%s/%s_mistagFormat.fasta' % (folder, curdir)
    with open(fastout, 'w') as o:
        for idx, seq_combis in enumerate(derep.items()):
            seq = seq_combis[0]
            combis = seq_combis[1]
            for combi, (n, headers) in combis.items():
                o.write('>seq%s;size=%s;fwd=%s;rev=%s\n%s\n' % (
                    idx, n, combi[0], combi[1], seq
                ))
    return fastout


def write_merging_cmd(out, cmd):
    if os.path.isfile(out):
        with open(out, 'a') as o:
            o.write('%s\n' % ' '.join(cmd))
    else:
        with open(out, 'w') as o:
            o.write('%s\n' % ' '.join(cmd))


def calculate_mistag_stats(tab, ncr):
    tot_exp = sum(tab.groupby(['status', 'combi']).agg('reads').sum())
    tot_ncr = sum(ncr.sum())
    total = float(tot_ncr + tot_exp)
    return tot_exp, tot_ncr, total


def write_designs(ost, design, primers, tab, ncr, basic_stats):
    tab_sums = tab.groupby(['status', 'combi']).agg('reads').sum()
    tab_counts = tab.groupby(['status', 'combi']).agg(np.count_nonzero)['reads']
    fs = sorted(primers['F'])
    rs = sorted(primers['R'])
    exp_design = {}
    tot_exp, tot_ncr, total = basic_stats
    types = ['Design', 'Non-critical mistags (% of total)', 'Critical mistags (%)']
    for tdx, typ in enumerate(types):
        ost.write('>%s\n' % typ)
        ost.write('\t%s\n' % '\t'.join(fs))
        for r in rs:
            # matrix representation of the samples
            if tdx == 0:
                mat_data = ['X' if (f, r) in design else '' for f in fs]
            # matrix representation of the detectable mistags
            elif tdx == 1:
                mat_data = []
                for f in fs:
                    if (f, r) in ncr.columns:
                        n_ncr = ncr.loc[:, (f, r)].sum()
                        p_ncr = round(float(n_ncr / total) * 100, 5)
                        mat_data.append(p_ncr)
                    else:
                        mat_data.append('')
            # matrix representation of the detected cross-contaminations
            else:
                mat_data = []
                for f in fs:
                    if (f, r) in design:
                        k = ('kept', '%s__%s' % (f, r))
                        rm = ('removed', '%s__%s' % (f, r))
                        n_kept = fecth_df_info(tab_sums, k)
                        m_kept = fecth_df_info(tab_counts, k)
                        n_removed = fecth_df_info(tab_sums, rm)
                        m_removed = fecth_df_info(tab_counts, rm)
                        n_total = float(n_kept + n_removed)
                        if n_total:
                            p_kept = round((n_kept / n_total) * 100, 2)
                            q_kept = round((m_kept / float(m_kept + m_removed)) * 100, 2)
                        else:
                            p_kept = 0.00
                            q_kept = 0.00
                        mat_data.append(p_kept)
                        exp_design[(f, r)] = [
                            design[(f, r)],
                            n_total, n_removed, m_removed,
                            n_kept, m_kept, p_kept, q_kept
                        ]
                    else:
                        mat_data.append('')
            ost.write('%s\t%s\n' % (r, '\t'.join(map(str, mat_data))))
    return exp_design


def fecth_df_info(tab, k):
    if k in list(tab.index.values):
        return tab.loc[k]
    else:
        return 0


def write_filtering(ost, stats_merging, basic_stats, exp_design):
    tot_exp, tot_ncr, total = basic_stats
    ost.write('# Filtering\n')
    ost.write('>Total reads in expected samples:\t%s\t%s %% (of all %s reads)\n'
              % (tot_exp, round(float(tot_exp/total), 2), total))
    tot_kept = sum([x[4] for x in exp_design.values()])
    tot_rm = sum([x[2] for x in exp_design.values()])
    perc_tot_kept = round((tot_kept/float(tot_exp))*100, 2)
    perc_tot_rm = round((tot_kept/float(tot_exp))*100, 2)
    ost.write('>Kept reads in expected samples:\t%s\t%s %%\n' % (
        tot_kept, perc_tot_kept))
    ost.write('>Removed reads from expected samples:\t%s\t%s %%\n' % (
        tot_rm, perc_tot_rm))
    ost.write('>Column output for expected samples\n')
    if len(stats_merging):
        ost.write('Forward\tReverse\tSample\tdtd_demultiplexed_reads\t\
merged_reads\tRemoved reads\tRemoved ISUs\tPass-filter reads\t\
Pass-filter ISUs\tPass-filter ISUs (%)\tPass-filter reads (%)\n')
    else:
        ost.write('Forward\tReverse\tSample\tdtd_demultiplexed_reads\tRemoved \
reads\tRemoved ISUs\tPass-filter reads\tPass-filter \
ISUs\tPass-filter ISUs (%)\tPass-filter reads (%)\n')
    for combi in sorted(exp_design):
        F, R = combi
        sam = exp_design[combi][0]
        if len(stats_merging):
            fil = [x for x in stats_merging if sam in x][0]
            ost.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (
                F, R, sam, stats_merging[fil][0], stats_merging[fil][1],
                '\t'.join(map(str, exp_design[combi][2:]))
            ))
        else:
            ost.write('%s\t%s\t%s\n' % (
                F, R, '\t'.join(map(str, exp_design[combi]))
            ))

