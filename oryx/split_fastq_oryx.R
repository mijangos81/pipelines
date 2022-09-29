library(tictoc)

tic()
file_names_R1 <- list.files("/scratch/mo73/Qatar",pattern = "R1",full.names = TRUE)
file_names_R1 <- file_names_R1[order(file_names_R1)]
b_name_R1 <- basename(file_names_R1)

file_names_R2 <- list.files("/scratch/mo73/Qatar",pattern = "R2",full.names = TRUE)
file_names_R2 <- file_names_R2[order(file_names_R2)]
b_name_R2 <- basename(file_names_R2)


for(i in 1:length(file_names_R1)){
  
  system(paste(
    
    "/home/549/jm9807/fastp",
    "-i", file_names_R1[i],
    "-I", file_names_R2[i],
    "-AGQL", 
    "-w 4",
    "-S 2000000",
    "-d 0", 
    "--out1", paste0("/scratch/mo73/Fastq_split_oryx/",b_name_R1[i]),
    "--out2", paste0("/scratch/mo73/Fastq_split_oryx/",b_name_R2[i])
    
  ))
}

toc()

