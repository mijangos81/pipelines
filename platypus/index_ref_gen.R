
library(tictoc)
library(tools)

path_ref_genome_zip <- "/scratch/mo73/Luis/GCF_004115215.2_mOrnAna1.pri.v4_genomic.fna.gz"
path_ref_genome_unzip <- "/scratch/mo73/Luis/GCF_004115215.2_mOrnAna1.pri.v4_genomic.fna"

tic("Indexing reference genome")
 system(paste(
   "bwa index -a bwtsw", path_ref_genome_zip
 ))
 toc()
 
 tic("Generate the sequence dictionary")
 system(paste(
   "java -jar picard.jar CreateSequenceDictionary",
   "R=", path_ref_genome_zip
 ))
 toc()
 
 tic("Creating the FASTA index file,using the unzipped file")
 system(paste(
   "samtools faidx", path_ref_genome_unzip
 ))
 toc()