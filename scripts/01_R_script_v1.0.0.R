library(Biostrings)
library(taxonomizr)
library(tidyr)
library(dplyr)
library(rgbif)
library(tibble)
library(purrr)
library(data.table)

#-- DOWNLOAD SEQUENCES FROM NCBI --#
#-- CUT SEQUENCES WITH PRIMERS THAT YOU ARE USING TO AMPLIFY THE MARKER (USE CUTADAPT OR SIMILAR) --#
#-- CHECK FOR LENGTH (ca. 200 bp UP TO EXPECTED LENGTH OF YOUR MARKER) --#


#-- LOAD SEQUENCES --#

sequences <- readDNAStringSet("CUT_SEQUENCES_OF_CORRECT_LENGTH.fasta")
sequences <- sequences[order(names(sequences), decreasing = F)]


#-- EXTRACT ONLY ACCESSION NUMBERS FROM THE HEADERS AND RENAME SEQUENCES AND SAVE ACCESSIONS AS TEXT --#

accessions <- as.vector(names(sequences))
accessions <- gsub(" .*","", accessions)
write(accessions, file="accessions.txt")
names(sequences) <- accessions


#-- DOWNLOAD THE ACCESSIONS TO TAXIDS TABLE FROM NCBI "nucl_gb.accession2taxid" FROM THE FTP SERVER HERE: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/ AND UNPACK THE FILE --#
#-- GREP THE TAXONOMY INFORMATION OF YOUR ACCESSION NUMBERS FROM THE FILE --#

taxids <- system("grep -f accessions.txt nucl_gb.accession2taxid", intern = T)
taxids <- fread(text = taxids, sep = "\t", header = F)


#-- ONLY RETAIN SEQUENCES THAT RETURNED A NCBI TAXID --#

sequences2 <- sequences[names(sequences) %in% taxids$V2]


#-- USE TAXONOMIZR TO GET FULL NCBI TAXONOMY --#
#-- FOR THE PREPARATION OF THE DATABASE THE FILES "names.dmp" AND "nodes.dmp" ARE NEEDED. TAXONOMIZR SHOULD DOWNLOAD THESE FILES AUTOMATICALLY --#
#-- IF NOT THEY CAN BE FOUND INSIDE THE TAXDUMP ARCHIVES FROM NCBI https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/ --#

prepareDatabase(getAccessions = F)
taxa <- getTaxonomy(taxids$V3, "nameNode.sqlite") #taxids$V3 are the NCBI taxids associated with your accession numbers
taxa <- as.data.frame(taxa)


#-- REPLACE SPECIAL CHARACTERS IN NAMES BY "_" --#

taxa$species <- gsub("[:;(),.]", "_", taxa$species)
taxa$species <- gsub("__", "_", taxa$species)


#-- MAKE A DATA FRAME WITH NCBI TAXONOMY FOR CHECKING AGAINST OTHER TAXONOMIES I.E. WITH RGBIF --#

ncbi_taxonomy <- data.frame(scientificName = gsub("_", " ", taxa$species),
                            rank = "SPECIES",
                            phylum = taxa$phylum,
                            class = taxa$class,
                            order= taxa$order,
                            family = taxa$family,
                            genus = taxa$genus)

ncbi_taxonomy_unique <- unique(ncbi_taxonomy)


#-- FUNCTION FOR MATCHING AGAINST GBIF BACKBONE WITH RGBIF SINCE THERE ARE PROBABLY MANY ENTRIES TO CHECK, NAME MATCHING IN CHUNKS WILL BE MORE EFFICIENT --#

chunk_backbone <- function(dat, chunk_size = 1000, pause = 0.2) {
  stopifnot(is.data.frame(dat), chunk_size > 0)
  
  if (nrow(dat) == 0) {
    return(data.frame())
  }
  
  chunk_id <- (seq_len(nrow(dat)) - 1L) %/% chunk_size + 1L
  chunks <- split(dat, chunk_id)
  
  results <- lapply(seq_along(chunks), function(i) {
    message(
      "Processing chunk ", i, " of ", length(chunks),
      " (", nrow(chunks[[i]]), " rows)"
    )
    
    result <- rgbif::name_backbone_checklist(chunks[[i]], verbose=F)
    
    # Small delay to avoid sending requests too quickly
    Sys.sleep(pause)
    
    result
  })
  
  bind_rows(results, .id = "chunk")
}


#-- MATCH NCBI TAXONOMY AGAINST GBIF TAXONOMY --#

gbif_taxonomy <- chunk_backbone(ncbi_taxonomy_unique, chunk_size = 1000)
taxonomy <- data.frame(ncbi_taxonomy_unique, gbif_taxonomy)


#-- IN CASE AN ENTRY CONTAINS AN ACCEPTED USAGE KEY CHANGE USAGE KEY TO ACCEPTED USAGE KEY --#

taxonomy <- taxonomy %>%
  mutate(usageKey = ifelse(!is.na(acceptedUsageKey), acceptedUsageKey, usageKey))

#-- SINCE THIS SCRIPT IS FROM BEFORE THE PUBLICATION OF RGBIF 3.9, THE MATCH AGAINST THE CATALOGUE OF LIFE TAXONOMY WILL BE DONE WITH THIS WORKAROUND --#
#-- AFTER THE PUBLICATION OF RGBIF 3.9 NAME MATCHING AGAINST THE CATALOGIE OF LIFE BACKBONE SHOULD WORK FROM THE PACKAGE --#
#-- CATALOGUE OF LIFE PROVIDES A MAPPING FILE BETWEEN GBIF TAXONOMY AND CATALOGUE OF LIFE HERE: https://download.checklistbank.org/col/gbif/README.html --#

#-- EXPORT GBIF USAGE KEYS TO GREP LINES FROM THE MAPPING FILE --#

write(gbif_taxonomy$usageKey, file="gbif_usageKeys.txt")


#-- NOW WE USE COMMANDLINE TOOLS TO GET THE MATCHES FROM THE FILE --#
##- awk -F'\t' 'NR==FNR{ids[$1]; next} $1 in ids' gbif_usageKeys.txt gbif-col265.tsv > gbif2col.tsv --#
#-- AND LOAD THE RESULTS BACK IN R --#

gbif2col <- read.delim("gbif2col.tsv", header = T)
gbif2col$gbif.ID <- as.character(gbif2col$gbif.ID)
gbif2col[gbif2col==""] <- NA

#-- AGAIN CHANGE ACCEPTED NAMES AND IDs WHERE APPLICABLE --#

gbif2col <- gbif2col %>%
  mutate(col.ID = ifelse(!is.na(col.acceptedID), col.acceptedID, col.ID)) %>%
  mutate(col.scientificName = ifelse(!is.na(col.acceptedScientificName), col.acceptedScientificName, col.scientificName))


#-- PUT INFORMATION FROM GBIF AND COL TOGETHER --#

taxonomy <- left_join(taxonomy, gbif2col, join_by(usageKey == gbif.ID))


#-- MAP THE RESULTS TO ALL SEQUENCES --#

all_taxonomy <- ncbi_taxonomy %>%
  left_join(taxonomy, by=(c("phylum", "class", "order", "family", "genus", "scientificName")))


#-- PREPARE DATA FRAMES WITH TAXONOMY YOU WANT TO USE (IN THIS CASE THE COL TAXONOMY), BUT KEEP THE TAX IDENTIFIERS FOR ALL THREE --#

col_headers_df <- data.frame(accession=names(sequences2),
                             phylum=all_taxonomy$col.phylum,
                             class=all_taxonomy$col.class,
                             order=all_taxonomy$col.order,
                             family=all_taxonomy$col.family,
                             genus=all_taxonomy$col.genus,
                             species=all_taxonomy$col.scientificName,
                             ncbitax=taxids$V3,
                             gbiftax=all_taxonomy$usageKey,
                             coltax=all_taxonomy$col.ID)


#-- IDENTIFY SEQUENCES THAT CONTAIN INFORMATION UP TO THE SPECIES LEVEL --#

col_headers_df_species <- subset(col_headers_df, !is.na(col_headers_df$species))


#-- OPTIONALLY WRITE REMAINING TAXONOMIES AS CSV FOR FURTHER MANUAL CHECK (ONLY A GENUS NAME IN THE SPECIES COLUMN, ETC.) --#

write.csv2(col_headers_df_species, file="col_headers_species.csv")
#-- MAKE CHANGES TO THE FILE AS NECCESSARY AND LOAD IT BACK INTO R --#
col_headers_df_species_2 <- read.csv2("col_headers_species.csv")

#-- col_headers_df_species_2 CONTAINS THE INFORMATION OF ALL THE SEQUENCES YOU WANT TO KEEP --#
#-- SUBSET YOUR SEQUENCES WITH THE INFORMATION FROM THIS DATA FRAME --#

col_sequences <- sequences2[names(sequences2) %in% col_headers_df_species_2$accession]


#-- MAKE NEW HEADERS FOR SEQUENCES IN YOU REFERENCE DATABASE AND RENAME THEM --#

col_headers <- paste0(col_headers_df_species_2$accession,
                       ";tax=p:", col_headers_df_species_2$phylum,
                       ",c:", col_headers_df_species_2$class,
                       ",o:", col_headers_df_species_2$order,
                       ",f:", col_headers_df_species_2$family,
                       ",g:", col_headers_df_species_2$genus,
                       ",s:", gsub(" ", "_",col_headers_df_species_2$species),
                       ",ncbitax:", col_headers_df_species_2$ncbitax,
                       ",gbiftax:", col_headers_df_species_2$gbiftax,
                       ",coltax:", col_headers_df_species_2$coltax)

names(col_sequences) <- col_headers

#-- WRITE THE SEQUENCES WITH THE NEW HEADERS TO A FASTA FILE --#

writeXStringSet(col_sequences, file="col_sequences.fasta")


#-- FROM HERE CHANGE TO THE COMMAND LINE --#

#-- REMOVE DUPLICATE SEQUENCES (SEQUENCES, WHICH HAVE THE SAME TAXONOMY AND THE EXACT SAME SEQUENCE, SEQUENCES WITH 100% SIMILARITY, BUT DIFFERENT TAXONOMIES ARE KEPT) --#
#-- FIND DUPLICATES WITH SEQKIT --#
#-- seqkit rmdup col_sequences.fasta -s -i -D duplicates.txt
#-- RUN THE "remove_duplicates.py" SCRIPT PROVIDED FROM THE GIT REPOSITORY --#
#-- THIS WILL GIVE A LIST OF ACCESSION NUMBERS (*.remove.txt), WHICH CAN BE REMOVED FROM THE FASTA FILE WITH "seqkit" --#
#-- (EASIEST WAY TO GREP THE ACCESSIONS FROM THE FASTA TO GET THE FULL HEADER AND THEN USE "seqkit" WITH INVERSE MATCH TO REMOVE THEM FROM THE FASTA FILE) --#
#-- (IF THE CATALOGUE OF LIFE TAXONOMY IS USED, IT MAY CONTAIN NON-ASCII CHARACTERS, WHICH HAVE TO BE CHANGED) --#

#-- USE THE NCBI COMMAND LINE TOOLS TO PERFORM A SELFBLAST OF YOU REFERENCE DATABASE --#
#-- USE "makblastdb -in col_sequences.fasta" -dbtype nucl -out NAME_OF_BLAST_DB" TO GENERATE A BLAST DATABASE FROM YOUR SEQUENCES --#
#-- BLAST YOUR SEQUENCES AGAINST THIS DATABASE: "blastn -query col_sequences.fasta -db NAME_OF_BLAST_DB -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" -max_target_seqs 50 -num_threads 6 -out selfblast1.out" --#
#-- THE CUSTOM OUTPUT FORMAT IS USED BECAUSE THE "qcovs" INFORMATION IS USED FOR SCORING --#

#-- ANALYSE THE OUTPUT FILE WITH THE "selfblast_analysis.py" SCRIPT PROVIDED IN THE GIT REPOSITORY --#
#-- THIS GENERATES CSV FILES "blast_refused.csv" CONTAINS SEQUENCES, WHICH CAN IMMEDIATELY BE REMOVED FROM THE FASTA FILE, "blast_check.csv" CONTAINS SEQUENCES WHICH NEED TO BE CHECKED MANUALLY (DEPENDING ON THE HIGHEST LEVEL OF TAXONOMIC DIFFERENCE: USUALLY GENUS AND SPECIES LEVEL ARE KEPT, SINCE THE PROBABILITY OF CLOSELY RELATED SPECIES HAVING VERY SIMILAR SEQUENCES IS HIGH) --#
#-- AFTER THE REMOVAL OF UNWANTED SEQUENCES (BEST: "seqkit" SIMILAR TO THE REMOVAL OF DUPLICATES) THE SELFBLAST STEPS ARE REPEATED UNTIL THE RESULTS ARE SATISFACTORY --#



#-- CREATING (LOCAL) SCORING FILES FOR USE WITH MEAN-CURE ANALYSIS PIPELINE --#
#-- LOAD THE TAXONOMIC INFORMATION FROM THE HEADERS OF THE FASTA AFTER THE LAST ROUND OF SELFBLAST AND PREPARE SCORING TABLE --#

selfblast <- readDNAStringSet("selfblast.fasta")
selfblast_headers <- data.frame(V1 = names(selfblast))
selfblast_headers <- selfblast_headers %>%
  separate_wider_delim(V1, ";", names_sep = "", too_few = "align_start") %>%
  separate_wider_delim(V12, ":", names_sep = "", too_few = "align_start") %>%
  separate_wider_delim(c(V122, V123, V124, V125, V126, V127, V128, V129), ",", names_sep = "", too_few = "align_start")

score <- data.frame(Accession_Number = selfblast_headers$V11, 
                    Phylum = selfblast_headers$V1221, 
                    Class=selfblast_headers$V1231, 
                    Order=selfblast_headers$V1241, 
                    Family=selfblast_headers$V1251, 
                    Genus=selfblast_headers$V1261,
                    Species=selfblast_headers$V1271,
                    ncbitax=selfblast_headers$V1281,
                    gbiftax=selfblast_headers$V1291,
                    coltax=selfblast_headers$V1210
                    )


#-- LOAD THE SELFBLAST OUTPUT --#

self_accepted <- read.csv2("col_selfblast/sb3_out/blast_accepted.csv", sep=",")
self_low <- read.csv2("col_selfblast/sb3_out/blast_low_coverage.csv", sep=",")
self_single <- read.csv2("col_selfblast/sb3_out/blast_singletons.csv", sep=",")
self_check <- read.csv2("col_selfblast/sb3_out/blast_check2.csv", sep=",")


#-- SET SCORES FOR THE DIFFERENT CATEGORIES AND COMBINE (IF THE HIGHEST LEVEL OF DIFFERENCE IN THE CATEGORY "CHECK" IS SPECIES POINTS WILL BE AWARDED SINCE CLOSELY RELATED SPECIES ARE LIKELY TO HAVE VERY SIMILAR SEQUENCES --#

self_accepted$score <- 4
self_low$score <- 2
self_single$score <- 2
self_check$score <- 0
self_check$score[self_check$highest_diff_level =="S"] <- 2  
self <- bind_rows(self_accepted, self_low, self_single, self_check)
self$quality <- 0
self$quality[self$score == 4] <- "gold"
self$quality[self$score == 2] <- "silver"
self$quality[self$score == 0] <- "bronze"
score <- score[,-c(11:22)]


#-- ADD LOCAL SCORE FOR SPECIES OCCURRING IN AREA OF INTEREST --#
#-- FOR NOW DEV VERSION OF RGBIF --#

install.packages("devtools")
devtools::install_github("ropensci/rgbif", ref = "v3.9.0")


#-- EXTRACT CATALOGUE OF LIFE IDENTIFIERS AND MAKE AN OCCURRENCE DOWNLOAD WITH RGBIF (SEE RGBIF DOCUMENTATION FOR PARAMETERS) --#
#-- IN THIS CASE WE'RE LOKKING FOR OCCURRENCES FROM GERMANY, LATER THAN 1980 --#

col_keys <- unique(score$coltax)

occ_down <- occ_download(
  pred_default(), 
  pred_in("taxonKey", col_keys),
  pred("country", "DE"), 
  pred_gte("year", 1980),
  format = "SPECIES_LIST",
  user = "",
  pwd="",
  email=""
)


#-- IMPORT THE SPECIES LIST WITH OCCURENCES --#

dl <- occ_download_get(occ_down, overwrite = T) %>%
  occ_download_import


#-- SET A THRESHOLD FOR THE MINIMUM NUMBER OF OCCURRENCES FOR WHICH YOU ACCEPT "LOCAL PRESENCE" --#

dl$score[dl$numberOfOccurrences >= 50] <- 5
dl2 <- subset(dl, dl$score == 5)


#-- JOIN THE LOCAL SCORES TO THE SCORE DATA FRAME AND CREATE COLUMNS FOR INTEGRATED LOCAL SCORE AND WRITE AS CSV FILE --#

score_dl <- left_join(score, dl2, join_by(coltax == taxonKey))
score_dl <- score_dl[,-c(14:34)]
score_dl$score.y[is.na(score_dl$score.y)] <- 0
score_dl$Score_ger <- score_dl$score.x + score_dl$score.y
score_dl <- score_dl %>%
  select(Accession_Number, Phylum, Class, Order, Family, Genus, Species, ncbitax, gbiftax, coltax, category, score.x, score.y, Score_ger, quality)
colnames(score_dl) <- c("Accession_Number", "Phylum", "Class", "Order", "Family", "Genus", "Species", "ncbitax", "gbiftax", "coltax", "selfblast_category", "Score", "occurrence_score_ger", "Score_ger", "quality")
write.csv2(score_dl, file="score.csv", quote = F, row.names = F)
