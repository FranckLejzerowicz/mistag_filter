### mistag_filter
Metabarcoding sequence abundance filter for projects involving combinations of short tagged primers for the labeling of amplicons and samples multiplexing

## Mandatory reading
Esling, P., Lejzerowicz, F., & Pawlowski, J. (2015). Accurate multiplexing and filtering for high-throughput amplicon-sequencing. _Nucleic acids research_, *43*(5), 2513-2524.
Reference url: https://academic.oup.com/nar/article-lookup/doi/10.1093/nar/gkv107

# Abstract extract
Tagging amplicons with tag sequences appended to PCR primers allow the multiplexing of numerous samples for high-throughput sequencing (HTS). This approach is routinely used in HTS-based diversity analyses, especially in microbial ecology and biomedical diagnostics. However, amplicon library preparation is subject to pervasive sample sequence cross-contaminations as a result of tag switching events referred to as mistagging. Here, we sequenced seven amplicon libraries prepared using various multiplexing designs in order to measure the magnitude of this phenomenon and its impact on diversity analyses. Up to 28.2% of the unique sequences correspond to undetectable (critical) mistags in single- or saturated double-tagging libraries. We show the advantage of multiplexing samples following Latin Square Designs in order to optimize the detection of mistags and maximize the information on their distribution across samples. We use this information in designs incorporating PCR replicates to **filter the critical mistags and to recover the exact composition** of mock community samples. Being **parameter-free and data-driven, our approach can provide more accurate and reproducible HTS data sets**, improving the reliability of their interpretations.

## Usage

The script necessitates a particular input format that corresponds to a fasta layout of each sample unique sequence. The sequences must be labelled with both a forward tagged primer and a reverse tagged primer. 
It only works for short primer constructs composed of a N-nt long tag appended at the 5'-end of the primers are allowed (N can be any length).
two found 

```
```

### Optional arguments
