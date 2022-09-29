# loading packages
library(ape)
library(dartR)
library(stringr)
library(dplyr)
# installing the developing version of dartR
# gl.install.vanilla.dartR(flavour = "dev")
path_dir <- getwd()
# reading chromosome information
chrom_platypus <- read.csv(paste0(path_dir,"/final_platypus/info/chrom_platypus.csv"))
chrom_platypus <- chrom_platypus[1:31,]

# Creating folders
chr_names <- paste0(path_dir,"/final_platypus/chr_",chrom_platypus$Name)
lapply(chr_names,dir.create)

############################################################################
#################### Splitting FASTA ######################################
############################################################################
system(paste("gzip -dk ",paste0(path_dir,"/final_platypus/reference_genomes/mOrnAna1_v4/GCF_004115215.2_mOrnAna1.pri.v4_genomic.fna.gz")))
fasta <- read.dna("final_platypus/reference_genomes/mOrnAna1_v4/GCF_004115215.2_mOrnAna1.pri.v4_genomic.fna", format = "fasta")

for(i in 1:nrow(chrom_platypus)){
  
  cat("Spllitting Chromosome",chrom_platypus[i,"Name"],"\n")
  
  fasta_tmp <- fasta[ grep(chrom_platypus[i,"RefSeq_mOrnAna_v4"], names(fasta) ) ]

  file_fasta <- paste0(chr_names[i],"/chr_",chrom_platypus[i,"Name"],"_mOrnAna1.pri.v4.fa")

  write.FASTA(fasta_tmp,file = file_fasta)

  system(paste("gzip",file_fasta))
  
}

############################################################################
#################### Splitting GFF ######################################
############################################################################
system(paste("gzip -dk ",paste0(path_dir,"/final_platypus/original_files/GCF_004115215.2_mOrnAna1.pri.v4_genomic.gff.gz")))
gff <- read.table(paste0(path_dir,"/final_platypus/original_files/GCF_004115215.2_mOrnAna1.pri.v4_genomic.gff"), sep="\t", quote="")

for(i in 1:nrow(chrom_platypus)){
  
  cat("Spllitting Chromosome",chrom_platypus[i,"Name"],"\n")
  
  gff_tmp <- gff[grep(chrom_platypus[i,"RefSeq_mOrnAna_v4"], gff[,1]),]
  
  file_gff <- paste0(chr_names[i],"/chr_",chrom_platypus[i,"Name"],"_mOrnAna1.pri.v4.gff")
  
  write.table(gff_tmp,file = file_gff,row.names = FALSE,col.names = FALSE,quote = FALSE,sep="\t" )
  
  system(paste("gzip",file_gff))
  
}

############################################################################
#################### Splitting DArT ######################################
############################################################################
# loading data
x_dart <- gl.read.dart(paste0(path_dir,"/final_platypus/original_files/Report_DPla19-4171_SNP_mapping_2.csv"), ind.metafile = paste0(path_dir,"/final_platypus/info/ID_Pop_Platypus.csv"), nas = "-", probar = TRUE)
x_dart <- gl.recode.pop(x_dart, pop.recode = paste0(path_dir,"/final_platypus/info/new_pop_assignments.csv"))
x_dart$other$loc.metrics$Chrom_Platypus_Chrom_NCBIv1 <- 
  str_extract(x_dart$other$loc.metrics$Chrom_Platypus_Chrom_NCBIv1, "[^_]*_[^_]*")
#amending data
#Gilad said: "There was never T42 and I suspect it was a recapture of T28".
# samples T3 (reported as male) and T5 (reported as female) are the same sample.
# sample T3 has more missing data. 
# E32 and V34 have missing data in more than 50% of loci
x_dart <- gl.drop.ind(x_dart,ind.list=c("T42","T3","E32","V34") )
# these samples were mislabeled, they belong to the opposite population.
# just the population is changed:
V32 <- which(x_dart$ind.names=="V32")
V30 <- which(x_dart$ind.names=="V30")

gps_temp_1 <- x_dart$other$latlong[V32,]
gps_temp_2 <- x_dart$other$latlong[V30,]
x_dart$other$latlong[V32,] <- gps_temp_2
x_dart$other$latlong[V30,] <- gps_temp_1
x_dart$pop[V32] <- "OVENS"   
x_dart$pop[V30]  <- "MITTA_ABOVE" 

#filtering
x_dart <- gl.filter.callrate(x_dart)
# gl.report.rdepth(x_dart)
x_dart <- gl.filter.rdepth(x_dart,lower = 10,upper = 100)
x_dart <- gl.filter.reproducibility(x_dart)
x_dart <- gl.filter.allna(x_dart)
x_dart <- gl.filter.monomorphs(x_dart)

for(i in 1:nrow(chrom_platypus)){
  
  cat("Spllitting Chromosome",chrom_platypus[i,"Name"],"\n")
  
  chrom_dart <- which(x_dart$other$loc.metrics$Chrom_Platypus_Chrom_NCBIv1 %in% chrom_platypus[i,"RefSeq_mOrnAna_v1"] == TRUE)

  dart_tmp <- gl.keep.loc(x_dart,loc.list = locNames(x_dart)[chrom_dart])
  
  if(nLoc(dart_tmp)>1){
    
    dart_file <- paste0(chr_names[i], outfile = paste0("/chr_",chrom_platypus[i,"Name"],"_dart.rds"))
    
    gl.save(dart_tmp,file = dart_file)
  
  # # the executable of plink should be in the working directory
  # gl2vcf(dart_tmp,
  #        outpath = chr_names[i],
  #        outfile = paste0("/chr_",chrom_platypus[i,"Name"],"_dart"),
  #        snp_pos = "ChromPos_Platypus_Chrom_NCBIv1",
  #        snp_chr = "Chrom_Platypus_Chrom_NCBIv1",
  #        sex_code =x_dart$other$ind.metrics$Sex )
  }else{
    next()
  }
}

############################################################################
#################### Splitting Oxford ######################################
############################################################################

#indexing the vcf file
system(
  paste(
    paste0(path_dir,"/final_platypus/programs/bcftools index -f"),
  paste0(path_dir,"/final_platypus/original_files/popgen.vcf.gz")
  )
)
# loading the information of chromosomes mapped from Oxford's reference genome
# to the latest reference genome
oxf_v4 <- read.csv(paste0(path_dir,"/final_platypus/info/contigs_Oxford.csv"))
oxf_v4 <- oxf_v4[order(oxf_v4$GenBank_Accn),]

# the system used by the BED format is zero-based for the coordinate start and 
# one-based for the coordinate end. Thus, the nucleotide with the coordinate 1 in 
# a genome will have a value of 0 in column 2 and a value of 1 in column 3.
# A thousand-base BED interval with the following start and end:
# chr7    0    1000
# writing bed file to convert coordinates between genomes using liftOver
bed_ox_contigs <- cbind(oxf_v4$GenBank_Accn,0,oxf_v4$Sequence_Length)
bed_file <- paste0(path_dir,"/final_platypus/info/preLift_ox_contigs.bed")
write.table(bed_ox_contigs,file = bed_file,col.names = FALSE,row.names = FALSE,quote = FALSE,sep = "\t")
liftOver_exe <- paste0(path_dir,"/final_platypus/programs/liftOver")
chain_file_ox <-  paste0(path_dir,"/final_platypus/info/GCA_002966995.1ToGCF_004115215.2.over.chain.gz")
conversion_file <- paste0(path_dir,"/final_platypus/info/ox_contigs_mOrnAna1_v4.bed")
unMapped_file <- paste0(path_dir,"/final_platypus/info/ox_contigs_unMapped")

# Granting access to the executable 
# and making it executable 
# This should be done just once
# system(paste0("chmod +x ",liftOver_exe))
# system(paste0("chmod 755 ",liftOver_exe))

# running liftOver
system(
  paste(
    liftOver_exe,
    bed_file,
    chain_file_ox,
    conversion_file,
    unMapped_file
  )
)

#loading mapped contigs
map <- read.table(conversion_file)

# loading unmapped contigs
unmap <- read.table(unMapped_file)

contig_unmap <- unlist(lapply(unmap$V1,function(x){
  which(oxf_v4$GenBank_Accn == x)
}))

oxf_v4 <- oxf_v4[-contig_unmap,]
oxf_v4$mOrnAna1_v4 <- map$V1

oxf_4_tmp <- which(oxf_v4$mOrnAna1_v4 %in% chrom_platypus$RefSeq_mOrnAna_v4)

oxf_v4 <- oxf_v4[oxf_4_tmp,]
colnames(oxf_v4) <-  c("Sequence_Name","GenBank_Accn","Sequence_Length","RefSeq_mOrnAna_v4")

oxf_v4_2 <- merge(oxf_v4,chrom_platypus[,c("Type","Name","RefSeq_mOrnAna_v4")],by='RefSeq_mOrnAna_v4')

split_oxf_v4 <- split(oxf_v4_2,f=oxf_v4_2$RefSeq_mOrnAna_v4)

# writing the regions file for bcftools
lapply(split_oxf_v4,function(x){
  tmp <- x[,c("Sequence_Name","Sequence_Length")]
  tmp_2 <- as.data.frame(cbind(tmp$Sequence_Name,1,tmp$Sequence_Length))
  
  write.table(tmp_2,
              file = paste0(path_dir,
                            "/final_platypus/",
                            "chr_",
                            x$Name[1],
                            "/",
                            "chr_",
                            x$Name[1],
                            "_Oxford_to_mOrnAna1_v4_regions.bed"),
              col.names = FALSE,
              row.names = FALSE,
              quote = FALSE,
              sep = "\t")
})

#runninf bcftools
for(i in 1:nrow(chrom_platypus)){
  cat("Spllitting Chromosome",chrom_platypus[i,"Name"],"\n")
  
  system(paste(
   paste0(path_dir,"/final_platypus/programs/bcftools view -O z"),
    "-o", paste0(path_dir,
                 "/final_platypus/",
                  "chr_",chrom_platypus[i,"Name"],
                 "/",
                 "chr_",chrom_platypus[i,"Name"],
                 "_oxford.vcf.gz"),
    paste0(path_dir,"/final_platypus/original_files/popgen.vcf.gz"), 
    "-R",paste0(path_dir,
                "/final_platypus/",
                "chr_",chrom_platypus[i,"Name"],
                "/",
                "chr_",chrom_platypus[i,"Name"],
                "_Oxford_to_mOrnAna1_v4_regions.bed")
  ))
  
}
