# mistag_filter

Metabarcoding sequence abundance filter for projects involving combinations of short tagged primers for the labeling of amplicons and samples multiplexing.

The script necessitates at least a specific **fasta file** and an experimental **design file** inputs.

## Fasta file
Should contain a copy of each unique sequence for each sample. The header sequence descriptions must include at least three semi-colon separated fields starting with:
* ```size=```: number of sequence read copies (integer)
* ```for=```: name of the *forward* tagged primer found in the sample sequence
* ```rev=```: name of the *reverse* tagged primer found in the sample sequence

#### example:
```
>unique_sequence_1;size=1000;for=F515-X;rev=R806-Y
ATTGCGGATTATTCGGGAGGGCGGCGAGAGCGTATATTCTAGGCGGATTCTGAC
>unique_sequence_1;rev=R806-B;size=1;for=F515-R
ATTGCGGATTATTCGGGAGGGCGGCGAGAGCGTATATTCTAGGCGGATTCTGAC
```
with ```F515``` being the forward amplification primer name and ```F515-X``` and ```F515-R``` being two different tagged versions of this primer, and ```R806``` being the reverse amplification primer name and ```R806-Y``` and ```R806-B``` being two different tagged versions of this primer.

The sequences in the fasta must be labelled with both a forward tagged primer and a reverse tagged primer. Hence, this script only works for short primer constructs composed of N-nt long tag sequences appended at the 5'-ends of the primers (N can be any length), as in e.g. Gloor, Gregory B., et al. "Microbiome profiling by illumina sequencing of combinatorial sequence-tagged PCR products." _PloS one_ **5.10** (2010): e15406.

**Note**: the primer sequences are not necessary, as a demultiplexing tool must have been used that is able to keep the count of the non-critical mistags (e.g. the sequences labeled with forward/reverse tagged primers combinations corresponding to unexpected samples absent from the experimental design).

## Design file
Should be a table containing the tags-to-samples information with samples as rows and at least three columns named in the header as follows:
* ```sample```: unique sample name
* ```for```: forward tagged primer used to PCR-amplify the sequences of the sample
* ```rev```: reverse tagged primer used to PCR-amplify the sequences of the sample

#### example

```
sample,for,rev
S1,F515-X;R806-Y
S2,F515-R,R806-B
...
```
A design provides the information about the expected samples and could be represented as a matrix, as follows for the above example:

n | F515-A | F515-R | F515-? | F515-X | F515-Y
:---:|:---:|:---:|:---:|:---:|:---:
R806-A | 0 | R | ... | X | 0 
R806-B | B | 1 | ... | BX | B 
R806-? | ... | ... | ... | ... | ...
R806-X | 0 | R | ... | X | 0 
R806-Y | Y | RY | ... | 1 | Y

with the possible expected samples being labeled by a "1" while the resulting unexpected samples are labeled by the letter of one or the other tagged primer they have in common with the expected sample. These unexpected sample contain the non-critical mistag sequences. It is based on the distribution of each sequence in these unexpected samples ("orthogonal samples") that the filter computes the modified Thompson Tau test rejection region to decide whether a sequence in an expected sample is also a mistag and should be removed.

## Reading and citing
Esling, P., Lejzerowicz, F., & Pawlowski, J. (2015). Accurate multiplexing and filtering for high-throughput amplicon-sequencing. _Nucleic acids research_, **43**(5), 2513-2524.
Reference url: https://academic.oup.com/nar/article-lookup/doi/10.1093/nar/gkv107

##### Abstract
Tagging amplicons with tag sequences appended to PCR primers allow the multiplexing of numerous samples for high-throughput sequencing (HTS). This approach is routinely used in HTS-based diversity analyses, especially in microbial ecology and biomedical diagnostics. However, amplicon library preparation is subject to pervasive sample sequence cross-contaminations as a result of tag switching events referred to as mistagging. Here, we sequenced seven amplicon libraries prepared using various multiplexing designs in order to measure the magnitude of this phenomenon and its impact on diversity analyses. Up to 28.2% of the unique sequences correspond to undetectable (critical) mistags in single- or saturated double-tagging libraries. We show the advantage of multiplexing samples following Latin Square Designs in order to optimize the detection of mistags and maximize the information on their distribution across samples. We use this information in designs incorporating PCR replicates to __**filter the critical mistags and to recover the exact composition**__ of mock community samples. Being __**parameter-free and data-driven, our approach can provide more accurate and reproducible HTS data sets**__, improving the reliability of their interpretations.

## Usage

```
mistag_filter.py [-h] -i I -d D [-o [O]] [-sep [SEP]] [-a [float between 0 and 1, max. 3 decimals]] [--out]
```

### Optional arguments

```
  -h, --help            show this help message and exit
  -i I                  Input fasta file name (required)
  -d D                  Multiplexing design file name (required)
  -o [O]                Output fasta file (default = input appended with
                        'mistagFiltered.txt')
  -a [float between 0 and 1, max. 3 decimals]
                        Alpha level for finding the Student's T critical value
                        for the Thomson Tau rejection region calulation
                        (default = 0.05)
  -s [S]                Field separator in multiplexing design file (default =
                        ',')
  --out                 Leave expected sample sequences out of non-critical
                        mistags distribution for calculations of the rejection
                        region (default = not active)
```

## Requirements
Python2.7<br />
numpy<br />
scipy<br />

