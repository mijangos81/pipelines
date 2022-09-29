
load_platy <- function(chrom){
  
  tic("Loading process took")
  
cat("LOADING CHROMOSOME", chrom,"...\n\n")
# working directory
path_dir <- getwd()
# reading chromosome information
chrom_platypus <-
  read.csv(paste0(path_dir, "/final_platypus/info/chrom_platypus.csv"))
chrom_platypus <- chrom_platypus[1:31, ]
# Chromosome to analyse
chrom <- chrom
chrom_row <- which(chrom_platypus$Name == chrom)
chrom <- paste0("chr_", chrom)
chrom_dir <- paste0(path_dir, "/final_platypus/", chrom, "/")

########################################################################
################### loading GFF file ####################################
########################################################################
cat("LOADING GFF FILE...\n")

# decompressing GFF
system(paste0("gzip -dk ", chrom_dir, chrom, "_mOrnAna1.pri.v4.gff.gz"))
gff <- read.gff(paste0(chrom_dir, chrom, "_mOrnAna1.pri.v4.gff"))
genes <- gff[gff$type == "gene", ]
genes_meta_tmp <- genes$attributes
genes_meta <-
  as.data.frame(str_split(genes_meta_tmp, pattern = "=|;|-", simplify = TRUE))
genes_meta <- genes_meta[, c(3, 5)]
colnames(genes_meta) <- c("name", "ID")
genes$attributes <- genes_meta$name
genes$length <- genes$end - genes$start

cat("  Objects 'genes'and 'gff' created\n\n")
print(table(gff$type))

########################################################################
################### loading Jenna samples vcf ##########################
########################################################################
cat("\nLOADING 26 GENOMES DATASET...\n")

# samples information
samples_info <-
  read.csv(paste0(path_dir, "/final_platypus/info/", "samples_jenna.csv"))
# decompressing vcf
system(paste0("gzip -dk ", chrom_dir, chrom, ".vcf.gz"))
gl_vcf <- gl.read.vcf(paste0(chrom_dir, chrom, ".vcf"),verbose = 0)
# filtering
# INFO codes in vcf file
# RFGQ_ALL "Empirical quality score (phred scaled) for the call - the geometric
# mean of all sample RFGQ probabilities"
# AC "Allele count in genotypes for each ALT allele, in the same order as listed"
# AN "Total number of alleles in called genotypes"
# DP "Combined depth across samples"
# MQ "RMS mapping quality"
# MQ0 "Number of MAPQ == 0 reads covering this record"
# NS "Number of samples with data"
# END "End position on CHROM"
# MP "Model posterior"
# names(gl_vcf$other$loc.metrics)
gl_vcf <- gl.filter.monomorphs(gl_vcf,verbose = 0)
#gl.report.locmetric(gl_vcf,metric ="RFGQ_ALL" )
gl_vcf <-
  gl.filter.locmetric(
    gl_vcf,
    metric = "RFGQ_ALL",
    lower = 20,
    upper = max(gl_vcf$other$loc.metrics$RFGQ_ALL),
    verbose = 0
  )
# assigning information to individuals and populations
gl_vcf$other$ind.metrics <- samples_info
indNames(gl_vcf) <-
  c(
    "S11",
    "S13",
    "S14",
    "S17",
    "S20",
    "S22",
    "S23",
    "S27",
    "S3",
    "S31",
    "S35",
    "S36",
    "S5",
    "S8",
    "S9",
    "T1",
    "T16",
    "T18",
    "T21",
    "T23",
    "T27",
    "T29",
    "T36",
    "T38",
    "T5",
    "T8"
  )
pop(gl_vcf) <- gl_vcf$other$ind.metrics$Group

gl_26 <- gl_vcf

cat("  Object 'gl_26' created with", nLoc(gl_26),"loci\n\n")

########################################################################
################### loading Oxford samples #############################
########################################################################
cat("LOADING OXFORD DATASET...\n")

# samples information
samples_info_ox <-
  read.csv(paste0(path_dir, "/final_platypus/info/", "locations_oxford.csv"))
# decompressing vcf
system(paste0("gzip -dk ", chrom_dir, chrom, "_oxford.vcf.gz"))
file_info <- file.info(paste0(chrom_dir, chrom, "_oxford.vcf"))
if(file_info$size > 14000){

gl_vcf_ox <- gl.read.vcf(paste0(chrom_dir, chrom, "_oxford.vcf"),verbose = 0)
# filtering
# names(gl_vcf_ox$other$loc.metrics)
gl_vcf_ox$other$loc.metrics$QUAL <-
  as.numeric(gl_vcf_ox$other$loc.metrics$QUAL)
#gl.report.locmetric(gl_vcf_ox,metric ="QUAL")

gl_vcf_ox <- gl_vcf_ox[order(gl_vcf_ox$ind.names), ]

# assigning information to individuals and populations
gl_vcf_ox$other$ind.metrics <- samples_info_ox
pop(gl_vcf_ox) <- gl_vcf_ox$other$ind.metrics$Region

#changing name of chromosome names to the used in genbank
contigs_ox <-
  read.csv(paste0(path_dir, "/final_platypus/info/contigs_Oxford.csv"))
t1 <- as.data.frame(as.character(gl_vcf_ox$chromosome))
colnames(t1) <- "Sequence_Name"
t2 <- left_join(t1, contigs_ox, by = 'Sequence_Name')
gl_vcf_ox$chromosome <- as.factor(t2$GenBank_Accn)

cat("  Remapping Oxford dataset...\n")

# the system used by the BED format is zero-based for the coordinate start and
# one-based for the coordinate end. Thus, the nucleotide with the coordinate 1 in
# a genome will have a value of 0 in column 2 and a value of 1 in column 3.
# A thousand-base BED interval with the following start and end:
# chr7    0    1000
# writing bed file to convert coordinates between genomes using liftOver
bed_ox <-
  cbind(as.character(gl_vcf_ox$chromosome),
        (gl_vcf_ox$position - 1),
        gl_vcf_ox$position)
bed_file <- paste0(chrom_dir, chrom , "_preLift_ox.bed")
write.table(
  bed_ox,
  file = bed_file,
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)
liftOver_exe <- paste0(path_dir, "/final_platypus/programs/liftOver")
chain_file_ox <-
  paste0(path_dir,
         "/final_platypus/info/GCA_002966995.1ToGCF_004115215.2.over.chain.gz")
conversion_file <- paste0(chrom_dir, chrom, "_mOrnAna1_v4.bed")
unMapped_file <- paste0(chrom_dir, chrom, "_unMapped")

# Granting access to the executable
# and making it executable
# This should be done just once
# system(paste0("chmod +x ",liftOver_exe))
# system(paste0("chmod 755 ",liftOver_exe))

# running liftOver
system(paste(
  liftOver_exe,
  bed_file,
  chain_file_ox,
  conversion_file,
  unMapped_file
))

#loading mapped SNPs
map <- read.table(conversion_file)

# loading unmapped SNPs
unmap <- read.table(unMapped_file)

loc_unmap <- unlist(lapply(unmap$V2, function(x) {
  which(gl_vcf_ox$position == x)
}))

gl_vcf_ox <-
  gl.drop.loc(gl_vcf_ox, loc.list = locNames(gl_vcf_ox)[loc_unmap],verbose = 0)

gl_vcf_ox$chromosome <- as.factor(map$V1)
gl_vcf_ox$position <- map$V2

loc_chr <-
  which(as.character(gl_vcf_ox$chromosome) != chrom_platypus[chrom_row, "RefSeq_mOrnAna_v4"])

gl_vcf_ox <-
  gl.drop.loc(gl_vcf_ox, loc.list = locNames(gl_vcf_ox)[loc_chr],verbose = 0)

gl_ox <- gl_vcf_ox

cat("  Object 'gl_ox' created with", nLoc(gl_ox),"loci\n")

}else{
  cat(" No loci in the Oxford dataset for this chromosome\n\n")
  gl_ox <- NULL
}

########################################################################
################### loading DArT samples ###############################
########################################################################
cat("LOADING DArT DATASET...\n")

gl_vcf_dart <- gl.load(paste0(chrom_dir, chrom, "_dart.rds"),verbose = 0)
gl_vcf_dart$chromosome <-
  as.factor(gl_vcf_dart$other$loc.metrics$Chrom_Platypus_Chrom_NCBIv1)
gl_vcf_dart$position <-
  gl_vcf_dart$other$loc.metrics$ChromPos_Platypus_Chrom_NCBIv1

# decompressing fasta
ref_gen_file_zip <- paste0(chrom_dir, chrom, "_mOrnAna1.pri.v4.fa.gz")
system(paste0("gzip -dk ", ref_gen_file_zip))
ref_gen_file_unzip <- paste0(chrom_dir, chrom, "_mOrnAna1.pri.v4.fa")
cat("  Remapping DArT dataset...\n")

# remapping using BLAST
# Installing BLAST
#'  You can download the BLAST installs from:
#'  https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
#'  It is important to install BLAST in a path that does not contain spaces for
#'  this function to work.

gl_vcf_dart <- gl.blast(gl_vcf_dart, ref_genome = ref_gen_file_unzip,verbose = 0)
pos <-
  gl_vcf_dart$other$loc.metrics$SnpPosition + gl_vcf_dart$other$loc.metrics$sstart

gl_vcf_dart$position <- as.integer(pos)

gl_dart <- gl_vcf_dart

cat("  Object 'gl_dart' created with", nLoc(gl_dart),"loci\n\n")

########################################################################
################### separating samples from Jenna's rivers #############
########################################################################
cat("SUBSETTING gl_dart TO COINCIDE WITH gl_26 ...\n")

gl_vcf_dart_jenna <-
  gl.keep.pop(gl_vcf_dart,
              pop.list = c("TENTERFIELD",
                           "SEVERN_BELOW",
                           "SEVERN_ABOVE"),verbose = 0)
gl_vcf_dart_jenna <-
  gl.filter.monomorphs(gl_vcf_dart_jenna, verbose = 0)

#changing the names so they coincide with Jenna's names
indNames(gl_vcf_dart_jenna) <-
  c(
    "T27",
    "T35",
    "S4",
    "S12",
    "S20",
    "S28",
    "S36",
    "T11",
    "T19",
    "T28",
    "T36",
    "S5",
    "S13",
    "S21",
    "S29",
    "S37",
    "T4",
    "T12",
    "T20",
    "T29",
    "S6",
    "S14",
    "S22",
    "S30",
    "S38",
    "T5",
    "T13",
    "T21",
    "T30",
    "T38",
    "S7",
    "S15",
    "S23",
    "S31",
    "S39",
    "T6",
    "T14",
    "T22",
    "T31",
    "T39",
    "S8",
    "S16",
    "S24",
    "S32",
    "S40",
    "T7",
    "T15",
    "T23",
    "T32",
    "T40",
    "S9",
    "S17",
    "S25",
    "S33",
    "S41",
    "T8",
    "T16",
    "T24",
    "S2",
    "T33",
    "T41",
    "S10",
    "S18",
    "S26",
    "S34",
    "T1",
    "T9",
    "T17",
    "T25",
    "S3",
    "T34",
    "S11",
    "S19",
    "S27",
    "S35",
    "T2",
    "T10",
    "T18",
    "T26"
  )

gl_dart_26 <- gl_vcf_dart_jenna

cat("  Object 'gl_dart_26' created with", nLoc(gl_dart_26),"loci\n\n")

########################################################################
# testing percentage of correct genotypes between DArT and new dataset##
########################################################################
cat("Testing percentage of correct genotypes between DArT and 26 genomes\n")

x_dart <- gl.keep.ind(gl_vcf_dart_jenna, ind.list = indNames(gl_vcf),verbose = 0)
x_dart <- gl.filter.callrate(x_dart, threshold = 1,verbose=0)
x_dart <- gl.filter.monomorphs(x_dart,verbose = 0)

x_jenna <- gl.filter.callrate(gl_vcf, threshold = 1,verbose = 0)

t2 <- which(x_dart$position %in% x_jenna$position)
t3 <- which(x_jenna$position %in% x_dart$position)

x_dart <- gl.keep.loc(x_dart, loc.list = locNames(x_dart)[t2],verbose = 0)
x_dart <- x_dart[, order(x_dart$position)]
x_dart <- x_dart[order(indNames(x_dart)), ]

x_jenna <- gl.keep.loc(x_jenna, loc.list = locNames(x_jenna)[t3],verbose = 0)
x_jenna <- x_jenna[, order(x_jenna$position)]
x_jenna <- x_jenna[order(indNames(x_jenna)), ]

s1 <- as.matrix(x_dart)
s2 <- as.matrix(x_jenna)

s1[s1 == 2] <- 0
s2[s2 == 2] <- 0

x_dart$position == x_jenna$position

# percentage of incorrect called genotypes between DArT and new dataset
 genotypes <- nInd(x_jenna) * nLoc(x_jenna)
 correct_genotypes <- sum(colSums(s1 == s2))

incorrect_percentage <-(correct_genotypes / genotypes) * 100

incorrect_percentage <-  round(incorrect_percentage,2)
cat(incorrect_percentage,"percentage of", genotypes,"genotypes were the same in DArT dataset and the 26 genomes\n\n")

return(list(genes = genes,
            gff = gff,
            gl_26 = gl_26,
            gl_ox = gl_ox,
            gl_dart = gl_dart,
            gl_dart_26 = gl_dart_26
              ))

toc()

}