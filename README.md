# Primer specific Reference Databases for Metabarcoding
## Introduction
This repository contains Python and R-scripts for the generation of a custom reference database for the analysis of Metabarcoding datasets. The Method uses primarily the scripts provided, but several steps are most easily executed through a command line interface, such as WSL (if you're working on a windows machine). Primer specific reference databases and score sheets for the analysis of Metabarcoding data are best used together with [Metabarcoding-Analysis-with-a-Custom-Reference](https://github.com/cschomb/Metabarcoding-Analysis-with-a-Custom-Reference) (see its README for details).
The folders for each primer pair contain the scoring file, an overview of the taxonomies in the reference, as well as the link to the Reference fasta file, and information about the different steps during database generation. Here is a detailed overview for the generation of the reference database for the ITS-3p62plF1 and ITS-4unR1 primers for Metabarcoding in plants:

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

### Primer specificity
Since the focus of these reference sequences was not to be limited by a only a certain subset of taxa, but to contain all available sequences, that can be amplified with a specific set of primers, the downloaded sequences were cut with [cutadapt](https://cutadapt.readthedocs.io/en/stable/) at the primer binding sites (in this case ITS-3p62plF1 and ITS-4unR1 from [Kolter & Gemeinholzer, 2020](https://cdnsciencepub.com/doi/10.1139/gen-2019-0198?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed)).
Further sequences were then excluded if they did not fall into an acceptable length range (in this case 200-500bp).

### Taxonomy
In order to generate a consistent taxonomy for the references, the downloaded sequences were assigned their NCBI taxonomy with the [taxonomizR](https://cran.r-project.org/web/packages/taxonomizr/index.html) R-package and subsequently matched against the Catalogue of Life (COL) taxonomic backbone with [rgbif](https://www.gbif.org/tool/81747/rgbif). Only sequences with species level information were kept and their corresponding NCBI, COL and GBIF taxids were written into the fasta header along with theit Accession numbers and taxonomy in the SINTAX format.

```[Accession number];tax=p:[PHYLUM],c:[CLASS],o:[ORDER],f:[FAMILY],g:[GENUS],s:[SPECIES],ncbitax:[NCBI TAXONOMY],gbiftax:[GBIF TAXONOMY],coltax:[COL TAXONOMY]```

Example:

```OR234967.1;tax=p:Tracheophyta,c:Magnoliopsida,o:Boraginales,f:Boraginaceae,g:Borago,s:Borago_officinalis,ncbitax:13363,gbiftax:2926110,coltax:MJ92```

To remove redundant information, a custom script was used along with [SeqKit](https://bioinf.shenwei.me/seqkit/) to identify records that showed the exact same taxonomy and identical nucleotide sequence.

### Quality of the reference sequences
The biggest caveat in using publicly available sequences from big repositories is their variable quality. Therefore, two steps were taken to identify whether a certain sequence can be reliably used as a reference.

#### Selfblast
In this step the whole database was blasted against itself. Assuming that there are enough different sequences for any given species, the BLAST results should largely consist of the same species (or at least species from the same genus if the barcode gap is too small). We set thresholds to classify the sequences into distinct categories, depending on the frequency of BLAST results with identical taxonomies. Sequences with >75% identical taxonomies, were accepted, while sequences with <25% identical taxonomies were rejected. Sequences in the 25% to 75% bracket were checked manually if the difference to the other BLAST hits from the same query was from phlyum to family level. Sequences where the difference was in the genus or species level were provisionally accepted, but will receive a different score (see Scoring part). This whole step was repeated until there were no more sequences automatically rejected (In this case a total of two times). 

### Scoring
Scoring of the reference sequences was done by allocating points depending on the categories from the selfblast step. Species in the reference file were checked against GBIF for at least 50 occurrences since 1980 in Germany. Five additional points were added to the sequences of species with occurrences. The "Score" column can be used as an unbiased score, that includes the whole reference database, while the "Score_ger" column is for use as a local reference (= Score + local score).

## How to use
### Database generation
The R-script in the scripts folder is designed to guide you through the process.

Dependencies:
- R (including the libraries: Biostrings, taxonomizr, tidyr, dplyr, rgbif, tibble, purr, data.table)
- CLI (including Python 3.x (see individual scripts for further Python packages), cutadapt, Seqkit, ncbi-blast+)

### Workflow and scripts
- Starting with a FASTA file containing sequences from NCBI, the first step is to cut the sequences with the Primers used in the Metabarcoding study with cutadapt, and filter the resulting sequences for the desired length with SeqKit.

- In order to assign taxonomic information to the sequences, you need the Accessions to Taxid table provided by NCBI ("nucl_gb_accession2taxid"). (can be downloaded from the NCBI ftp server [here](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/). With this file you can search the accession numbers in your database for their taxid ("grep" the list of accession numbers created in the R script).

- With the list of taxids you can use taxonomizr to get the full taxonomic lineage of each sequence, from which you can chose the ranks of interest (i.e. Phylum, Class, Order, Family, Genus, Species).

- For a consistent taxonomic backbone you can match the resulting taxonomy table against the Catalogue of Life (or GBIF) backbone.

- From the results of this match you can contruct new FASTA headers for your sequences containing their Accession number, taxonomic information, and Taxon IDs from different backbones.

- The next steps will remove redundant information and ensure a high quality of reference sequences. First, duplicate sequences (same sequences with same taxonomic information) are removed with SeqKit and a custom Python script.
  ```
  seqkit rmdup sequences.fasta -s -i -D dupliucates.txt
  ```
  With the resulting file, you can run the 02_remove_duplicates.py script
  ```
  python3 remove_duplicates.py [--out-prefix TARGET_FOLDER] duplicates.txt
  ```
  This will results in two text files containing the set of accession numbers from the database, which are accepted (to keep) and one to remove (delete from database), which can also be done with Seqkit grep function.

- With a clean and non-redundant database, you can now identify sequences, which are likely to be wrongly annotated. To that end, you can blast the database against itself and analyse the results for consistency with the selfblast_analysis.py script
  ```
  makeblastdb -in sequence.fasta -dbtype nucl -out DB_NAME
  blastn -query sequences.fasta -db DB_NAME -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" -max_target_seqs 50 -num_threads 6 -out selfblast.out
  ```
  The blast output format is a custom one, because we want to include the information about the query coverage.

- The BLAST output can then further be analysed and put into the categories described above, by feeding it into the selfblast_analysis.py script.
  ```
  python3 selfblast_analysis.py [-o OUTPUT_DIR] [-t LENGTH_THRESHOLD] selfblast.out [--resume]
  ```
  The script will produce a number of output files:
    - filtered_blast_results.csv (BLAST hits, which fall within 3% of identity of the best BLAST hit)
      
    | q_acc | q_phylum | q_class | q_order | q_family | q_genus | q_species | s_acc | s_phylum | s_class | s_order | s_family | s_genus | s_species | pident | length | mismatches | gapopens | qstart | qend | sstart | send | evalue | bitscore | qcov | score | score_interval |
    |-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|
    | Accession number of the query | Phylum of the query | Class of the query | Order of the query | Family of the query | Genus of the query | Species of the query | Accession of the subject | Phylum of the subject | Class of the subject | Order of the subject | Family of the subject | Genus of the subject | Species of the subject | Percent identity of the alignment | Length of the alignment | Mismatches in the alignment | Opened gaps in the alignment | Alignment start in the query | Alignment end in the query | Alignment start in the subject | Alignment end in the subject | E-Value of the Alignment |Bitscore of the ALignment | Query coverage | Score of the alignment (= pident + (qcov /5)) | Interval of the Score |
  
    - species_count.txt (a text file containing the amount of sequences per species; used to check for low coverage and sigleton species)
   
    - summarised_blast_results.csv (summary of the filtered BLAST results, one query sequence per line)
      
      | q_acc | q_phylum | q_class | q_order | q_family | q_genus | q_species | diff_counts | highest_diff_level | group_size | score interval | difference_values |
      |-|-|-|-|-|-|-|-|-|-|-|-|
      | Accession number of the query | Phylum of the query | Class of the query | Order of the query | Family of the query | Genus of the query | Species of the query | Overview on how many sequences with different levels of taxonomy were found in the filteres BLAST hits | Highest level of taxonomic defferences in the BLAST hits | Amount of sequences retained from the BLAST hits after filtering | Interval of the score from the BLAST | Detailed view of different taxonomies within the filtered BLAST hits, incl. amounts |

    - Five .csv files with categorised BLAST results:
        - accepted - >= 75% of BLAST hits have the same taxonomy as the query sequence
        - low_coverage - <= 3 sequences with the same taxonomy as the query sequence are in the database
        - singleton - sequences with unique taxonomies
        - check - the BLAST hits have between 25% and 75% the same taxonomy as the query sequence
        - refused - <= 25% of the BLAST hits have the same taxonomy as the query sequence
       
        | q_acc | q_phylum | q_class | q_order | q_family | q_genus | q_species | diff_counts | highest_diff_level | group_size | score interval | difference_values | species_count | category | 
        |-|-|-|-|-|-|-|-|-|-|-|-|-|-|
        | Accession number of the query | Phylum of the query | Class of the query | Order of the query | Family of the query | Genus of the query | Species of the query | Overview on how many sequences with different levels of taxonomy were found in the filteres BLAST hits | Highest level of taxonomic defferences in the BLAST hits | Amount of sequences retained from the BLAST hits after filtering | Interval of the score from the BLAST | Detailed view of different taxonomies within the filtered BLAST hits, incl. amounts | Amount of sequences with the same taxonomy in the database | Category of the sequence (accepted, low_coverage, singleton, check, refused)

    - Query sequences which were refused can be removed from the reference database (similar as described above with the inverted grep function of SeqKit). Similarly, the user can decide whether to incorporate the low abundance and singleton sequences into the reference database (since the selfblast quality control mechanism does not apply here). For the entries categorised under 'check', you can manually inspect the sequences and decide to keep or remove them. Sequences where the taxonomic differences in the BLAST are on the species or genus level can be expected to be very similar, since they are closely related (they should be incorporated into the reference database). The BLAST hits with higher taxonomic differences (Phylum - Family) usually contain a handful of wrongly annotated sequences, which can be identified by BLAST against the nr database on the BLAST homepage and then removed from you r sequences.

- These database pruning steps of the database should be repeated until you are satisfied with results (usually two or three rounds of selfblast and analysis).

- With the reference database finished, you can go back to the R-script and create a scroring file with the analysis from the selfblast steps (the categories translate to quality information of the reference sequence)

  | Category | Score | Quality |
  |-|-|-|
  | accepted | 4 | gold |
  | low_coverage | 2 | silver |
  | singleton | 2 | silver |
  | check | 0 | bronze |
  | check with highest difference level at species | 2 | silver |

- For the incorporation of a local occurence score into the reference database, the R-script has an example to query GBIF for occurences with parameters (in this case a country, a year after which occurrences are considered and a threshold for occurrences). The sequences of species considered to have local occurrences according to these parameters are assigned additional score points in a new column (combined with the normal score) so they will have a higher overall score than non-local species.

- The reference database and scoring file generated with the steps above are made to be used with the [Metabarcoding-Analysis-with-a-Custom-Reference](https://github.com/cschomb/Metabarcoding-Analysis-with-a-Custom-Reference) repository. This script can takes the quality information from the scoring file into account, so that the resulting taxonomic assignment contains high quality sequences. With the additional information of local occurrences, the reference database will prefer these local species, while not sacrificing information from sequences, which may introduced or at least unexpected in the region of interest.
