import glob

from mistag_utils import *

def parse_primers(design):
    """
    Get all the forward and all the reverse primer name in a dict
    """
    primers = {}
    primers_rad = {}
    for frx, FR in enumerate(['F', 'R']):
        primers[FR] = {}
        # frx: 0 => 'F'
        # frx: 1 => 'R'
        for primer in list(set([x[frx] for x in design.keys()])):
            # primer will be either each of the pforward primers (frx = 0) of
            # the reverse primers (frx = 1)
            primers[FR][primer] = 1
            primers_rad.setdefault(FR, []).append(primer.split('-')[0])
    primers_rad = check_primer_rad(primers_rad)
    return primers, primers_rad


def parse_design(designFile):
    """
    Rarse the design file and collect the info about which primer combinations
    are expected Returns a dict with each expected primer combination as key
    (value set to sample name for output)
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


def parse_unexpected(mistag_fastqs, primers_rad):
    """Return the dicts with the unexpected mistag combinations per type of
    mistag and for each sequence of the fastq files
    """
    mistag_IDs = {}
    mistags_all = {'weird': {}, 'unexpected': {}}
    # parse both R1 and R2 mistag fastq files
    with open(mistag_fastqs[0]) as f1, open(mistag_fastqs[1]) as f2:
        idx = 0
        for id1, id2 in zip(f1, f2):
            idx += 1
            # fastq sequence header
            if idx == 1:
                combi, typ = get_mistag_combi_type(id1, id2, primers_rad)
                if typ:
                    mistag_IDs[id1.split()[0][1:]] = combi
                    # count the number of combinations per type of mistag
                    if combi in mistags_all[typ]:
                        mistags_all[typ][combi] += 1
                    else:
                        mistags_all[typ][combi] = 1
            if idx == 4:
                idx = 0
    return mistag_IDs, mistags_all


def parse_formatted_input(fastin, design, primers):
    """
    Read the fasta input file provided by the pipeline developed by Yoann
    Dufresne
    Sequence header format:
    >SEQ_ID;size=DEREPS_READS;for=TAGGED_PRIMER;rev=TAGGED_PRIMER
    Returns two dicts with [abundance, sequence ID] values for each sequence
    key, nested under each sample key
    One dict called "expected" for the sequences belonging to an expected
    samples of the original design
    One dict called "nonCritic" for the sequences associated with to a
    unexpected primer combination
    """
    seq = ''
    curSeq = {}
    expected = {}
    nonCritic = {}
    # parse each sequence of the appropriately formatted fasta file
    for lindx,line in enumerate(open(fastin, 'rU')):
        # progrssion bar
        if lindx and lindx % 10 == 0:
            print('\r', 'Number of parsed per-sample ISUs:', lindx, end='')
            sys.stdout.flush()
        # only check the header for size / fwd / rev
        if line[0] == '>':
            # fill 'expected' and 'nonCritic' dicts with a parsed sequence
            if len(seq):
                # update the dicts for expected tage combis and non-critical mistag combis
                update_dict(pair, expected, nonCritic, seq, n, seqID, primers, f, r, design)
                curSeq = {}
                seq = ''
            # current sequence info
            splitID = line.strip()[1:].split(';')
            seqID = splitID[0]
            for field in splitID[1:]:
                if field.startswith('size='):
                    n = re.search('\d+', field).group(0)
                elif field.startswith('fwd=') or field.startswith('for='):
                    f = field.split('=')[-1].strip()
                elif field.startswith('rev=') or field.startswith('rv'):
                    r = field.split('=')[-1].strip()
            # get the primer pair
            pair = tuple([f,r])
        else:
            seq += line.strip()
    # fill 'expected' and 'nonCritic' dicts with the last parsed sequence
    update_dict(pair, expected, nonCritic, seq, n, seqID, primers, f, r, design)
    print()
    return expected, nonCritic

