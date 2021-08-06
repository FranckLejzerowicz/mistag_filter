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

    # ==============
    # mistag parsing
    # ==============
    args = get_args()

    designFilin = args['d'][0]

    design = parse_design(designFilin)
    # design: {('for-X, rev-X'): sample, (...): ...,}
    # e.g. {('V4F-C', 'V4R-K'): 'V4.DNA.16.A', ('V4F-D', 'V4R-N'): 'V4.DNA.21.A', ...

    primers, primers_rad = parse_primers(design)
    # primers: {'F': {'for|rev-x': 1, 'for|rev-y': 1', ...}, 'R': {...}}
    # e.g. {'R': {'V4R-V': 1, 'V4R-G': 1, 'V4R-Y': 1, ...
    # primers_rad: {'F': 'for', 'R': 'rev'}
    # e.g. {'R': 'V4R', 'F': 'V4F'}

    stats_merging = {}
    fastin = args['i'][0]

    # =================
    # sequences merging
    # =================
    if os.path.isdir(fastin):

        fastin = fastin.rstrip('/')

        samples_fastqs, combi_fastqs = get_fastqs(fastin, design)
        # e.g. {'folder/sample': ['folder/sample_fwd.fastq',
        #                         'folder/sample_rev.fastq',
        #                         ('V4F-K', 'V4R-O')], ...

        mistag_fastqs = get_mistags_files(fastin)
        # e.g. ['folder/*_mistag_R1.fastq', 'folder/*_mistag_R2.fastq']

        mistag_unexpected, mistag_all = parse_unexpected(
            mistag_fastqs, primers_rad
        )
        # mistag_unexpected:
        # e.g. {'M02442:25:000000000-AR1YW:1:2113:23652:16643':
        #           ('V4F-D', 'V4R-C'),
        #       'M02442:25:000000000-AR1YW:1:1104:4368:11203':
        #           ('V4F-K', 'V4R-A'), ...
        # mistag_all:
        # e.g. {'unexpected': {('V4F-M', 'V4R-A'): 3,
        #                      ('V4F-L', 'V4R-G'): 10, ...,
        #       'weird': {('V4R-Q', 'V4R-E'): 1,
        #                 ('V4R-U', 'V4R-V'): 1, ...

        fastin, derep = perform_merging(
            args, fastin, samples_fastqs, mistag_fastqs,
            stats_merging, mistag_unexpected
        )

    # ================
    # mistag filtering
    # ================
    alpha = args['a']
    exclude = args['out']
    o_filt, o_fas, expected, non_critic, orthos = perform_filtering(
        args, fastin, design, primers, alpha, exclude, combi_fastqs, derep
    )
    print('Written:', o_fas)

    # ============
    # write output
    # ============
    write_output(
        args, o_filt, o_fas, expected, non_critic, orthos, design,
        primers, stats_merging
    )
    print('Written:', o_filt)


