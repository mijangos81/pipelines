
library(tictoc)
library(tools)

path_ref_genome_zip <- "/scratch/mo73/Luis/reference/oryx_raca_all.fa.gz"
path_ref_genome_unzip <- "/scratch/mo73/Luis/reference/oryx_raca_all.fa"

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