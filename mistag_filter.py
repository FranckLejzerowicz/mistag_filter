import os
from os.path import dirname, isdir
from math import sqrt
from scipy import stats
from numpy import std, mean
from skbio.io import read

from mistag_parser import parse_formatted_input
from mistag_getter import get_output_fp


def filter_mistags(args, expected, non_critic, design,
                   orthos, alpha, exclude, combi_fastqs, derep):
    """
    Mistag filter
    For each sample;
        get the othogonal samples (sharing the forward / the reverse primer)
        For each ISU;
            compute the distribution moments of this ISU abundance in these
            samples group
            compare the deviation of the current ISU with the rejection region
            of the Thompson Tau test
            (http://www.statisticshowto.com/modified-thompson-tau-test/)
    Returns the dict with kept (key "1") or removed (key "0") ISU lists for each
    sample (nested key)
    """
    o_rad, o_fas = get_output_fp(args)
    o_filt = '%s_tmpFilt.tsv' % o_rad
    o_exp = '%s_expected.tsv' % o_rad   # to be removed
    o_ncr = '%s_nonCritic.tsv' % o_rad   # to be removed
    # filt = {0: {}, 1: {}}
    with open(o_ncr, 'w') as ob:
        for combi in non_critic:
            com = '__'.join(combi)
            for isu, n_ID in non_critic[combi].items():
                cur_abund = int(n_ID[0])
                isu_name = str(n_ID[1])
                ob.write('%s\t%s\t%s\n' % (com, isu_name, cur_abund))
        # filt[1] => kept
        # filt[0] => removed
        # for each combi
    with open(o_filt, 'w') as ofilt, open(o_fas, 'w') as ofas, open(o_exp, 'w') as oa:
        for combi in sorted(expected):
            # get orthogonal combis to compute abundance distributions
            f_tag, r_tag = combi
            com = '__'.join(combi)
            sample = design[combi]
            ortho = orthos[combi]
            # for each ISU in the expected combis
            kept_reads = set()
            for isu, n_ID in expected[combi].items():
                cur_abund = int(n_ID[0])
                isu_name = str(n_ID[1])
                oa.write('%s\t%s\t%s\n' % (com, isu_name, cur_abund))
                # get the abundances of the sequence in the mistag copies
                ortho_abunds = [
                    int(non_critic[x][isu][0]) for x in ortho
                    if isu in non_critic[x]
                ]
                # decide whether to keep the current sequence abundance
                if exclude:
                    total_abund = ortho_abunds
                else:
                    total_abund = [cur_abund] + ortho_abunds
                # get the number of abundances in the distribution
                n = len(total_abund)
                # keep ISU if only one count (i.e. rare but not cross-conta?)
                if n <= 1:
                    ofilt.write('%s\t%s\t%s\t%s\tkept\n' % (
                        sample, com, isu_name, cur_abund))
                    ofilt.write('%s\t%s\t%s\t%s\tremoved\n' % (
                        sample, com, isu_name, 0))
                    ofas.write('>%s_%s;size=%s;for=%s;rev=%s\n%s\n' % (
                        isu_name, sample, cur_abund, f_tag, r_tag, isu))
                    if combi in derep[isu]:
                        kept_reads.update(derep[isu][combi][1])
                    continue
                # get distribution parameters
                mean_abund = mean(total_abund)
                dev_abund = std(total_abund)
                df = (n-1)
                cur_dev = abs(mean_abund - cur_abund)
                # get the critical student's t-test value
                a = stats.t.isf(alpha, df)
                # rejection region
                r = (( a * (n-1)) / (sqrt(n) * sqrt(n - 1 + (a ** 2)))) * dev_abund
                # fill the filtering results dict
                # ***************************************************************
                # REMOVE (abundance in rejection region)
                if cur_dev <= r:
                    ofilt.write('%s\t%s\t%s\t0\tkept\n' % (
                        sample, com, isu_name))
                    ofilt.write('%s\t%s\t%s\t%s\tremoved\n' % (
                        sample, com, isu_name, cur_abund))
                # KEEP (abundance in rejection region)
                else:
                    ofilt.write('%s\t%s\t%s\t%s\tkept\n' % (
                        sample, com, isu_name, cur_abund))
                    ofilt.write('%s\t%s\t%s\t0\tremoved\n' % (
                        sample, com, isu_name))
                    ofas.write('>%s_%s;size=%s;for=%s;rev=%s\n%s\n' % (
                        isu_name, sample, cur_abund, f_tag, r_tag, isu))
                    if combi in derep[isu]:
                        kept_reads.update(derep[isu][combi][1])
                # ***************************************************************

            if not args['fastq']:
                continue
            print(sample, len(kept_reads))
            fwd_fp, rev_fp = combi_fastqs[combi]
            up_dir = '%s/mistagFilt_fastqs' % dirname(fwd_fp)
            if not isdir(up_dir):
                os.makedirs(up_dir)
            fwd_fpo = '%s/%s_mistagFilt_R1.fastq' % (up_dir, sample)
            rev_fpo = '%s/%s_mistagFilt_R2.fastq' % (up_dir, sample)
            with open(fwd_fpo, 'w') as fwd_o, open(rev_fpo, 'w') as rev_o:
                rev = open(rev_fp)
                fdx = 0
                write_it = False
                for f in open(fwd_fp):
                    r = rev.readline()
                    if fdx == 0:
                        f_id = f.strip()[1:].split()[0]
                        if f_id in kept_reads:
                            write_it = True
                    if write_it:
                        fwd_o.write(f)
                        rev_o.write(r)
                    fdx += 1
                    if fdx == 4:
                        write_it = False
                        fdx = 0
    return o_filt, o_fas


def orthogonal_combis(expected, non_critic):
    """
    For each expected primer combination of the original design (expected),
    get the list of the primer combinations that are not in the original design
    and that share the same forward or reverse primers.
    Return dict: {(for-x, rev-x): []}
    """
    ortho = {}
    for combi in expected:
        cur_f = combi[0]
        cur_r = combi[1]
        same_fwd_idx = [x for x in non_critic.keys() if x[0] == cur_f]
        same_rev_idx = [x for x in non_critic.keys() if x[1] == cur_r]
        ortho[combi] = same_fwd_idx + same_rev_idx
    return ortho


def perform_filtering(
        args, fastin, design, primers, alpha,
        exclude, combi_fastqs, derep):

    # parse input fasta and collect the abundance distributions
    #   - for sequences in expected combis
    #   - for sequences in unexpected combis (non-critical mistags)
    expected, non_critic = parse_formatted_input(fastin, design, primers)
    # {(for-x, rev-x): {'seq1': [reads, ID], 'seq2': ... },
    #  (for-y, rev-y): {...}, ... }
    # expected: (for-x,rev-x) of expected combis sequences
    # nonCritic: (for-x,rev-x) of unexpected mistag sequences

    # unexpected combis sharing one tagged primer with each expected combi
    orthos = orthogonal_combis(expected, non_critic)

    o_filt, o_fas = filter_mistags(
        args, expected, non_critic, design, orthos,
        alpha, exclude, combi_fastqs, derep
    )
    return o_filt, o_fas, expected, non_critic, orthos
