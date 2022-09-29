library(stringr)
library(data.table)
library(tictoc)
library(tools)

chr_info <- read.csv(paste0(getwd(),"/chrom_platypus.csv"))

chr_info <- chr_info[17,]

path_folder_samples <- "/scratch/mo73/map"
#path_folder_samples <- "/Users/s441489/final_pipe/map"
path_out <- "/scratch/mo73/bam"
# path_out <- "/Users/s441489/final_pipe/bam"

bam_files <- list.files(path_folder_samples,pattern = "\\.bam$",full.names = T)
bam_files_names <- unlist(lapply(bam_files,function(x){
  tmp <- unlist(strsplit(x, "_|\\."))
  return(tmp[length(tmp)-1])
}))

bam_files_names <- bam_files_names[1]

chr_name <- paste0(chr_info$Type,"_",chr_info$Name)
chr_folder <- paste0(path_out,"/",chr_name)
lapply(chr_folder,dir.create)
tic("split bams")
for(i in 1:length(bam_files_names)){
  for(x in 1:length(chr_folder)){
    
    system(
      paste(
        "samtools view -bh -@8 -m16g",bam_files[i],chr_info[x,"RefSeq"],">", 
        paste0(chr_folder[x],"/",chr_name[x],"_",bam_files_names[i],".bam")
      )
    )
    
    system(paste(
      "samtools index -@8 -m16g", paste0(chr_folder[x],"/",chr_name[x],"_",bam_files_names[i],".bam")
    ))
    
  }
}

toc()


