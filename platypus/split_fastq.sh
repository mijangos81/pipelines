#!/bin/bash 

#Split each fastq pair into 2,000,000 lines (500,000 reads)

fastp=/home/549/jm9807/fastp

fqpair=`echo $1 | cut -d ',' -f 1`
file=$(basename $fqpair)

$fastp -i /scratch/mo73/fastq/S8R_HTH7NDSXX_GGACTTGG-CGCAGACG_L001_R1.fastq.gz \
	-I /scratch/mo73/fastq/S8R_HTH7NDSXX_GGACTTGG-CGCAGACG_L001_R2.fastq.gz \
	-AGQL \
	-w 4 \
	-S 2000000 \
	-d 0 \
	--out1 /scratch/mo73/Fastq_split_2/S8R_HTH7NDSXX_GGACTTGG-CGCAGACG_L001_R1.fastq.gz \
	--out2 /scratch/mo73/Fastq_split_2/S8R_HTH7NDSXX_GGACTTGG-CGCAGACG_L001_R2.fastq.gz
