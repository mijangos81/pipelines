library(stringr)
library(data.table)
library(tictoc)
library(tools)

sample_names <- c("T21","T29","T36","T8")

path_folder_samples <- "/scratch/mo73/Align_split_3"

for(i in 1:length(sample_names)){
sample_name <- sample_names[i]
print(sample_name)
merge_files <- list.files(path_folder_samples,pattern = paste0("(",sample_name,"_)(.*sort.bam)"),full.names = T)
merge_files <- paste("I=",merge_files,sep="")
merge_files <- do.call(paste,as.list(merge_files))

 output_step_3 <- paste0("/scratch/mo73/Align_merged_2/",sample_name,".merged.nameSorted.bam")
tic("merge")
system(paste(
  "java -Xmx48g -jar /home/549/jm9807/picard.jar MergeSamFiles",
  merge_files,
 "O=", output_step_3 ,
 "USE_THREADING=true",
 "MERGE_SEQUENCE_DICTIONARIES=true",
 "VALIDATION_STRINGENCY=SILENT"
))
toc()

 tic("Step_4")
output_step_4 <-
  paste0("/scratch/mo73/dup/", "output_duplicates_", sample_name, ".bam")
output_step_4b <-
  paste0("/scratch/mo73/dup/", "marked_dup_metrics", sample_name, ".txt")
system(
  paste(
    "java -Xmx48g -jar /home/549/jm9807/picard.jar MarkDuplicates REMOVE_DUPLICATES=true",
    "I=",
    output_step_3,
    "O=",
    output_step_4,
    "M=",
    output_step_4b
  )
)
toc()

 # file.remove(output_step_3)

system(paste("samtools index",
             output_step_4
)
)

tic("Step_5")
output_step_5 <-
  paste0("/scratch/mo73/map/", "mapped_sorted_", sample_name, ".bam")
system(
  paste(
    "java -Xmx48g -jar /home/549/jm9807/picard.jar SortSam",
    "I=",
    output_step_4,
    "O=",
    output_step_5,
    "SORT_ORDER=coordinate"
  )
)
toc()

system(paste("samtools index",
             output_step_5
))

# file.remove(output_step_4)

}