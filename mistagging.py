#!/usr/bin/env python3

import os

from mistag_parser import *
from mistag_getter import *
from mistag_writer import *
from mistag_merger import perform_merging
from mistag_filter import perform_filtering

__author__ = "Franck Lejzerowicz"
__copyright__ = "Copyright 2017, The Deep-Sea Microbiome Project"
__credits__ = ["Philippe Esling"]
__license__ = "GPL V3"
__version__ = "1.0"
__maintainer__ = "Franck Lejzerowicz"
__email__ = "franck.lejzerowicz@unige.ch"


if __name__ == '__main__':
    """
    Filter the critical mistags as in Esling et al. 2015 (NAR,
    https://doi.org/10.1093/nar/gkv107)
    To be used only for one amplification primer pair (for-/rev-)
    => Not possibly used for e.g. for1-:for2-/rev-
    """
    #
    # ==============
    # mistag parsing
    # ==============
    args = get_args()
    designFilin = args['d'][0]
    design = parse_design(designFilin)
    primers, primers_rad = parse_primers(design)
    # design: {('for-X, rev-X'): sample, (...): ...,}
    # primers: {'F': {'for|rev-x': 1, 'for|rev-y': 1', ...}, 'R': {...}}
    # primers_rad: {'F': 'for', 'R': 'rev'}

    stats_merging = {}
    fastin = args['i'][0]
    # =================
    # sequences merging
    # =================
    if os.path.isdir(fastin):
        fastin = fastin.rstrip('/')
        samples_fastqs = get_fastqs(fastin, design)
        mistag_fastqs = get_mistags_files(fastin, primers_rad)
        mistag_unexpected, mistag_all = parse_unexpected(mistag_fastqs,
                                                         primers_rad)
        # launch the samples merging as parallel multiprocesses
        multiproc = 0
        fastin = perform_merging(args, fastin, samples_fastqs,
                                    mistag_fastqs, stats_merging,
                                    mistag_unexpected, multiproc)
    # ================
    # mistag filtering
    # ================
    alpha = args['a']
    exclude = args['out']
    oFilt, oFas, expected, nonCritic, orthos = perform_filtering(args, fastin,
                                                                design,
                                                                primers,
                                                                alpha,
                                                                exclude)
    # ============
    # write output
    # ============
    print('Written:', oFas)
    write_output(args, oFilt, oFas, expected, nonCritic, orthos, design,
                  primers, stats_merging)
    print('Written:', oFilt)

