# Primer specific Reference Databases for Metabarcoding
## Introduction

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
In this step the whole database was blasted against itself. Assuming that there are enough different sequences for any given species, the BLAST results should largely consist of the same species (or at least species from the same genus if the barcode gap is too small). We set thresholds to classify the sequences into distinct categories:

```
BLAST results >= 75% same species as query → accepted
BLAST results <= 25% same species as query → rejected
BLAST results between 25% and 75% same species as query → check manually
Species with <= 3 sequences → accept as low_count

```

#### Species delimitation

### Scoring

### Final clean-up


## How to use
