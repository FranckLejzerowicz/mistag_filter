import sys
import re

def update_dict(pair, expected, nonCritic, seq,
                n, seqID, primers, f, r, design):
    """
    Add the demultiplexed sequence data into either the expected or nonCritic
    depending on whether the primer pair/sample is in the design or not
    """
    if pair in design:
        add_key(expected, pair, seq, n, seqID)
    elif f in primers['F'] and r in primers['R']:
        add_key(nonCritic, pair, seq, n, seqID)


def add_key(d, pair, seq, n, seqID):
    """
    Fill a nested dict with given info
    Register sequences associated with:
        - expected primer combinations (first call in update_dict)
        - non-critical combinations (second call in update_dict)
    """
    if pair in d:
        if seq in d[pair]:
            print('Error: sequence "%s" already encountered (with size=%s) for \
the current primer combination "%s" (now with size=%s)\nExiting' %
            (seq, d[pair][seq][0], ' + '.join(pair), n))
            sys.exit()
        else:
            d[pair][seq] = [int(n), seqID]
    else:
        d[pair] = {seq: [int(n), seqID]}


def check_primer_rad(primers_rad):
    for FR,list_rads in primers_rad.items():
        if len(set(list_rads)) != 1:
            print('More than one amplification primer in design\nExiting...')
            sys.exit()
    primers_rad = dict([x,y[0]] for x,y in primers_rad.items())
    return primers_rad


def check_fields(fields, info, filin):
    """
    Check if the minimum required info is present in the design file:
    """
    for field in fields:
        if field not in info:
            print('No field "%s" in %s\nExiting' % (field, filin))
            sys.exit()
    else:
        return 1


def rev_dict(d):
    """Simply reverse the key/value of a dict into another tree
    """
    D = {}
    for k,v in d.items():
        D[v] = k
    return D


def add_to_dict_list(S, d):
    """Append the list of fastq files with the tagged primers combination used to
    generate the PCR products labeled to demultiplex these fastq files
    """
    for s in S:
        for k, v in d.items():
            if k in s:
                S[s].append(v)
                break


def increment_nested(d, s, c):
    """Update the dict with the number of sequence copies per combination
    """
    if s in d:
        if c in d[s]:
            d[s][c] += 1
        else:
            d[s][c] = 1
    else:
        d[s] = {c: 1}


def get_mistag_combi_type(id1, id2, primers_rad):
    # combi in R1/R2 file sequences
    t1 = id1.strip().split(';tag:')[-1]
    t2 = id2.strip().split(';tag:')[-1]
    combi = tuple([t1,t2])
    # skip unassigned tags
    if 'unknown' in combi:
        return combi, 0
    # default to weird combination
    typ = 'weird'
    # unless it is a plausible combination
    if t1.split('-')[0] != t2.split('-')[0]:
        # sort the combination to have the forward primer first
        if t1.split('-')[0] == primers_rad['R']:
            combi = tuple([t2,t1])
        # collect the mistag combination for the current sequence
        typ = 'unexpected'
    return combi, typ

