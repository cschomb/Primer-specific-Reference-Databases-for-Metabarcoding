# Primer specific Reference Databases for Metabarcoding
## Introduction
This repository contains primer specific reference databases and score sheets for the analysis of Metabarcoding data. They are best used together with [Metabarcoding-Analysis-with-a-Custom-Reference](https://github.com/cschomb/Metabarcoding-Analysis-with-a-Custom-Reference) (see its README for details).
The zipped files for each primer pair contain the Reference fasta file, the scoring file, an overview of the taxonomies in the reference, as well as metadata for the reference sequences, and information about the different steps during database generation. Here is a detailed overview for the generation of the reference database for the ITS-3p62plF1 and ITS-4unR1 primers for Metabarcoding in plants:

## Database generation
### Sequence download from NCBI
Sequences were downloaded from [NCBI Nucleotide](https://www.ncbi.nlm.nih.gov/nuccore) with the following search terms:
```
5.8S ribosomal DNA gene
28S ribosomal DNA gene
18S ribosomal DNA gene
internal transcribed spacer 1
internal transcribed spacer 2
ITS
ITS1
ITS2
NOT whole genome
NOT complete genome
NOT partial genome
```

### Taxonomy
In order to generate a consistent taxonomy for the references, the downloaded sequences were assigned their NCBI taxonomy with the [taxonomizR](https://cran.r-project.org/web/packages/taxonomizr/index.html) R-package and subsequently matched against the GBIF taxonomic backbone with [rgbif](https://www.gbif.org/tool/81747/rgbif). Only sequences with species level information were kept and their corresponding NCBI and GBIF taxids were written into the fasta header along with theit Accession numbers and taxonomy in the SINTAX format.

```[Accession number];tax=p:[PHYLUM],c:[CLASS],o:[ORDER],f:[FAMILY],g:[GENUS],s:[SPECIES],gbiftax:[GBIF TAXONOMY],ncbitax:[NCBI TAXONOMY]```

Example:

```OR234967.1;tax=p:Tracheophyta,c:Magnoliopsida,o:Boraginales,f:Boraginaceae,g:Borago,s:Borago_officinalis,gbiftax:2926110,ncbitax:13363```



### Primer specificity
Since the focus of these reference sequences was not to be limited by a only a certain subset of clades, but to contain all available sequences, that can be amplified with a specific set of primers, the remaining sequences were cut with [cutadapt](https://cutadapt.readthedocs.io/en/stable/) at the primer binding sites (in this case ITS-3p62plF1 and ITS-4unR1 from [Kolter & Gemeinholzer, 2020](https://cdnsciencepub.com/doi/10.1139/gen-2019-0198?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed))

Further sequences were then excluded if they did not fall into an acceptable length range (in this case 200-500bp).

To remove redundant information, a custom script was used along with [SeqKit](https://bioinf.shenwei.me/seqkit/) to identify records that showed the exact same taxonomy and identical nucleotide sequence.

### Quality of the reference sequences
The biggest caveat in using publicly available sequences from big repositories is their variable quality. Therefore, two steps were taken to identify whether a certain sequence can be reliably used as a reference.

#### Selfblast
In this step the whole database was blasted against itself. Assuming that there are enough different sequences for any given species, the BLAST results should largely consist of the same species (or at least species from the same genus if the barcode gap is too small). We set thresholds to classify the sequences into distinct categories and analysed the results with a custom script:

```
BLAST results >= 75% same species as query → accepted
BLAST results <= 25% same species as query → rejected
BLAST results between 25% and 75% same species as query → check manually
Species with <= 3 sequences → accept as low_count

```
Sequences in the 25% to 75% bracket were checked manually if the difference to the other BLAST hits from the same query was from phlyum to family level. Sequences where the difference was in the genus or species level were provisionally accepted, but will receive a different score (see Scoring part). This whole step was repeated until there were no more sequences automatically rejected (In this case a total of two times). 

### Scoring
Scoring of the reference sequences was done by allocating points depending on the categories from the selfblast step.
```
Selfblast:
accepted:               4 pts --> gold
hits in the same genus: 2 pts --> silver
low_count / singleton:  2 pts --> silver
hits in other genera:   0 pts --> bronze

```


Species in the reference file were checked against GBIF for at least 10 occurrences since 1980 in Germany. Five additional points were added to the sequences of species with occurrences. The "Score" column can be used as an unbiased score, that includes the whole reference database, while the "Score_ger" column is for use as a local reference.


### Final clean-up
Sequences with a bad quality indicator (no points from either the selfblast or the mPTP analysis) were removed from the reference file.

## How to use
This reference database and corresponding are made for use with the [Metabarcoding-Analysis-with-a-Custom-Reference](https://github.com/cschomb/Metabarcoding-Analysis-with-a-Custom-Reference) repository. The reference fasta and scoring file can be put into the "input" folder of a project and set in the "settings.xlsx"
