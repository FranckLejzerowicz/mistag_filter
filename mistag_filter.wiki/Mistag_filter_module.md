# Mistag filter module

## The mistag filter module

This module allows you to discard demultiplexed sequences representing cross-contaminations due to the mistagging phenomenon (see [this article]() and this [wiki](https://github.com/FranckLejzerowicz/mistag_filter/wiki/Mistagging)). It is necessary that the
* sequences were multiplexed using both forward _and_ reverse tagged primers \<n-nt long tag sequence\> + \<amplification primer sequence\> prior to the library preparation step,
* dereplicated sequenced reads found with primer combination not corresponding to an expected sample are counted and present in the input fasta file.

### Module interactions

_**Attention**: the filtering must be library-wise, i.e. it must be applied only for filtering the mistag sequences of samples multiplexed in one library at a time._

#### Inputs
* Fasta file: The demultiplexed, dereplicated sequences (i.e. Individual Sequence Units, or "ISUs"), including one copy of each ISU for each tagger primers combination. The sequence headers must include at least three semi-colon separated fields starting with: **size=**_<int>_, **for=**_<forward primer>_, **rev=**_<reverse primer>_.
###### example:
```
>unique_sequence_1;size=1000;for=F515-X;rev=R806-Y
ATTGCGGATTATTCGGGAGGGCGGCGAGAGCGTATATTCTAGGCGGATTCTGAC
>unique_sequence_1;rev=R806-B;size=1;for=F515-R
ATTGCGGATTATTCGGGAGGGCGGCGAGAGCGTATATTCTAGGCGGATTCTGAC
>unique_sequence_2;rev=R806-A;size=190;for=F515-C
ATTGCGGATTCGAGGGGGTTTTAGGGCTAG
>unique_sequence_2;size=10;for=F515-G;rev=R806-A
ATTGCGGATTCGAGGGGGTTTTAGGGCTAG
>unique_sequence_2;size=4;rev=R806-T;for=F515-G
ATTGCGGATTCGAGGGGGTTTTAGGGCTAG
```
* Design file: a tab-separated file providing the information about which combination og forward/reverse tagged primers correspond to a sample. A header must be present with as least three columns named **sample** (sample name), **for** (forward tagged primer), and **rev** (reverse tagged primer).

 #### Outputs

* Mistag-filtered fasta file: Same format and header information as the input fasta file, but only for the sequences associated with tagged primers combinations corresponding to expected samples of the design file.

* Statistics file: Report file containing descriptive statistics about the filtering process. This file is a text file containing several sections:
- _"#"-commented header_ lines: inputs/outputs and filtering run info.
- "_Design_": design summaries including saturation level (% of employed combinations) and matrix representation for summarizing the mutliplexing design an the 
- "_# Filtering_": filtering results about sequences in expected samples (i.e. referred in the design file)
- "_# Non-critical mistags / Unexpected samples_": filtering results about sequences in unexpected samples (i.e. not referred in the design file)

#### Options

* alpha: floating number for the \alpha level for finding the Student's T critical value for the modified Thompson Tau rejection region calulation. Default to 0.05.
* out: boolean to decide whether to leave expected sample sequences out of non-critical mistags distributions for the calculations of Thompson's rejection regions (not recommended)

### Reference and citing

* [Esling _et al._ 2015](https://academic.oup.com/nar/article-lookup/doi/10.1093/nar/gkv107)
* Formula in [Supplementay data](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/nar/43/5/10.1093_nar_gkv107/1/gkv107_Supplementary_Data.zip?Expires=1503414493&Signature=LUZ3ydiM0jRgWYx3i2kzS-TsbmwJaA6j6NGOevQRW5siNCMegoyhlp6WUj2yykbMP9JZRGE-EMrNsT-XAdejZydc1RYme1RGHt0Dm1qv5tCH2g1cn4W~NSSQ1EEwdT1-xSIe5~Vf78GZQ8geYziFv~WLyjq~Hkd1Drp6~3OeEnxuDY-Am65o0Lu8vFyp~EiaCKj3bYq6ALGCxLV8y0oCj~2BqUITYOUFwLI9D3aWT6MC3KxXUjFJqgB-6vlDKqChMaqE3WI5Wlw6fn9lAMi~BeW8PzUKjSfVWhWUrCT~53D2~eGWeXdSVpDik4XUq0g0qZWPaYge4Cp1E7OROO8CtA__&Key-Pair-Id=APKAIUCZBIA4LVPAVW3Q)
