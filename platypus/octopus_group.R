library(stringr)
library(data.table)
library(tictoc)
library(tools)

#path_ref_genome_unzip <- "/Users/s441489/new_pipe/reference/GCF_004115215.2_mOrnAna1.pri.v4_genomic.fa"
path_ref_genome_unzip <- "/scratch/mo73/Luis/reference/GCF_004115215.2_mOrnAna1.pri.v4_genomic.fa"
forest <- "/scratch/mo73/octopus/aab_build/src/resources/forests/octopus_forests_germline.v0.7.4.forest"
#forest <-"/Users/s441489/octopus/resources/forests/octopus_forests_germline.v0.7.4.forest"
 
# path_output <- "/scratch/mo73/octopus_group/"
# path_folder_samples <- "/scratch/mo73/octopus_samples/"
# 
# fq_files <- list.files(path_folder_samples)
# fq_files <- fq_files[order(fq_files)]
# 
# fq_files_names <- str_split(fq_files, pattern = "_")
# fq_files_names <- rbindlist(lapply(lapply(fq_files_names, t), as.data.frame))
# colnames(fq_files_names) <- c("name", "pat1","pat_2","pattern", "extension")
# fq_files_names$complete_name <- fq_files
# split_names <- split(fq_files_names, fq_files_names$name)

#cores <- 16
#ram <- 80
 cores <- 8
 # ram <- 

tic("Step_7")
  output_step_7 <- "/scratch/mo73/fin_chr/chr_6.vcf"

 #group_vcf <- list.files("/Users/s441489/final_pipe/octopus_group",pattern ="^octopus",full.names = T)
 # group_vcf <- list.files("/scratch/mo73/octopus_group_test",pattern ="^octopus",full.names = T)
# group_vcf <- paste0(path_output,group_vcf)
# group_vcf <- do.call(paste,as.list(group_vcf))
#input_files <- list.files("/Users/s441489/final_pipe/map",pattern ="^mapped_sorted.*\\.bam$",full.names = T)

input_files <- list.files("/scratch/mo73/bam/_6",pattern =".*\\.bam$",full.names = T)

 input_files <- do.call(paste,as.list(input_files))

system(paste(
  # "/Users/s441489/octopus/bin/octopus",
   "/scratch/mo73/octopus/aab_build/src/octopus" ,
  "-R" ,  path_ref_genome_unzip ,
  "-I ", input_files,
  # "-T", "NC_041728.1",
  "-T", "NC_041733.1",
  # "--disable-denovo-variant-discovery" ,
  # "-c", group_vcf ,
  # "--min-variant-posterior 20",
  "-o" ,output_step_7,
 "-B 4Gb",
  #paste0("-X ",ram,"GB"),
  #paste0("-B ",ram,"GB"),
  "--threads ",cores,
  "--very-fast",
  "--sequence-error-model PCR-FREE.NOVASEQ",
  "--forest",forest
))

toc()
