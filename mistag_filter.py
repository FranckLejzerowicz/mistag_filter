from math import sqrt
from scipy import stats
from numpy import array, std, mean

from mistag_parser import parse_formatted_input
from mistag_getter import get_output_fp

def filter_mistags(args, expected, nonCritic, design, orthos, alpha, exclude):
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
    oRad, oFas = get_output_fp(args)
    oFilt = '%s_tmpFilt.tsv' % oRad
    oA = '%s_expected.tsv' % oRad   ## to be removed
    oB = '%s_nonCritic.tsv' % oRad   ## to be removed
   # filt = {0: {}, 1: {}}
    with open(oFilt, 'w') as ofilt, open(oFas, 'w') as ofas, open(oA, 'w') as oa, open(oB, 'w') as ob:
        for combi in nonCritic:
            com = '__'.join(combi)
            for ISU, n_ID in nonCritic[combi].items():
                curAbund = int(n_ID[0])
                isuName = str(n_ID[1])
                ob.write('%s\t%s\t%s\n' % (com, isuName, curAbund))
        # filt[1] => kept
        # filt[0] => removed
        # for each combi
        for combi in expected:
            # get orthogonal combis to compute abundance distributions
            F = combi[0]
            R = combi[1]
            com = '__'.join(combi)
            sample = design[combi]
            ortho = orthos[combi]
            # for each ISU in the expected combis
            for ISU, n_ID in expected[combi].items():
                curAbund = int(n_ID[0])
                isuName = str(n_ID[1])
                oa.write('%s\t%s\t%s\n' % (com, isuName, curAbund))
                # get the abundances of the sequence in the mistag copies
                orthoAbunds = [int(nonCritic[x][ISU][0]) for x in ortho if ISU
                               in nonCritic[x]]
                # decide whether to keep the current sequence abundance
                if exclude:
                    totalAbund = orthoAbunds
                else:
                    totalAbund = [curAbund] + orthoAbunds
                # get the number of abundances in the distribution
                n = len(totalAbund)
                # keep ISU if only one count (i.e. rare but not cross-conta?)
                if n <= 1:
                    ofilt.write('%s\t%s\t%s\t%s\tkept\n' % (sample, com,
                                                            isuName, curAbund))
                    ofilt.write('%s\t%s\t%s\t%s\tremoved\n' % (sample, com,
                                                            isuName, 0))
                    ofas.write('>%s_%s;size=%s;for=%s;rev=%s\n%s\n' % (isuName,
                                                                       sample,
                                                                       curAbund,
                                                                       F, R,
                                                                       ISU))
               #     filt[1].setdefault(combi,[]).append((ISU, isuName, curAbund))
                #    filt[0].setdefault(combi,[]).append((ISU, isuName, 0))
                    continue
                # get distribution parameters
                meanAbund = mean(totalAbund)
                devAbund = std(totalAbund)
                df = (n-1)
                curDev = abs(meanAbund - curAbund)
                # get the critical student's t-test value
                a = stats.t.isf(alpha, df)
                # rejection region
                r = (( a* (n-1)) / (sqrt(n) * sqrt(n - 1 + (a ** 2)))) * devAbund
                # fill the filtering results dict
                # ***************************************************************
                # REMOVE (abundance in rejection region)
                if curDev <= r:
                    ofilt.write('%s\t%s\t%s\t%s\tkept\n' % (sample, com,
                                                            isuName, 0))
                    ofilt.write('%s\t%s\t%s\t%s\tremoved\n' % (sample, com,
                                                               isuName,
                                                               curAbund))
             #       filt[1].setdefault(combi,[]).append((ISU, isuName, 0))
              #      filt[0].setdefault(combi,[]).append((ISU, isuName, curAbund))
                # KEEP (abundance in rejection region)
                else:
                    ofilt.write('%s\t%s\t%s\t%s\tkept\n' % (sample, com,
                                                            isuName, curAbund))
                    ofilt.write('%s\t%s\t%s\t%s\tremoved\n' % (sample, com,
                                                                   isuName, 0))
                    ofas.write('>%s_%s;size=%s;for=%s;rev=%s\n%s\n' % (isuName,
                                                                       sample,
                                                                       curAbund,
                                                                       F, R,
                                                                       ISU))
            #        filt[1].setdefault(combi,[]).append((ISU, isuName, curAbund))
             #       filt[0].setdefault(combi,[]).append((ISU, isuName, 0))
                # ***************************************************************
    return oFilt, oFas

def orthogonal_combis(design, expected, nonCritic):
    """
    For each expected primer combination of the original design (expected),
    get the list of the primer combinations that are not in the original design
    and that share the same forward or reverse primers.
    Return dict: {(for-x, rev-x): []}
    """
    ortho = {}
    for combi in expected:
        curF = combi[0]
        curR = combi[1]
        same_fwd_idx = [x for x in nonCritic.keys() if x[0] == curF]
        same_rev_idx = [x for x in nonCritic.keys() if x[1] == curR]
        ortho[combi] = same_fwd_idx + same_rev_idx
    return ortho


def perform_filtering(args, fastin, design, primers, alpha, exclude):
    # parse input fasta and collect the abundance distributions
    #   - for sequences in expected combis
    #   - for sequences in unexpected combis (non-critical mistags)
    expected, nonCritic = parse_formatted_input(fastin, design, primers)
    # {
    #   (for-x, rev-x):
    #     {'seq1':
    #       [reads, ID],
    #      'seq2': ...
    #     },
    #   (for-y, rev-y):
    #     {...},
    #   ...
    # }
    ## expected: (for-x,rev-x) of expected combis sequences
    ## nonCritic: (for-x,rev-x) of unexpected mistag sequences

    # unexpected combis sharing one tagged primer with each expected combi
    orthos = orthogonal_combis(design, expected, nonCritic)
    oFilt, oFas = filter_mistags(args, expected, nonCritic, design,
                                 orthos, alpha, exclude)
    return oFilt, oFas, expected, nonCritic, orthos
