#!/bin/bash

module load samtools/1.9

ref=/scratch/mo73/Luis/reference/oryx_raca_all.fa

fqpair=`echo $1 | cut -d ',' -f 1` 
fq1=${fqpair}_R1.fq.gz
fq2=${fqpair}_R2.fq.gz
sampleID=`echo $1 | cut -d ',' -f 2`
centre=`echo $1 | cut -d ',' -f 3`
lib=`echo $1 | cut -d ',' -f 4`
platform=`echo $1 | cut -d ',' -f 5`
flowcell=`echo $1 | cut -d ',' -f 6`
lane=`echo $1 | cut -d ',' -f 7`
NCPUS=`echo $1 | cut -d ',' -f 8` 	

outPrefix=$(basename $fqpair)

/scratch/mo73/Luis/NextGenMap-0.5.2/bin/ngm-0.5.2/ngm \
-r $ref \
-1 $fq1  \
-2 $fq2 \
-o /scratch/mo73/Align_split_oryx/${outPrefix}.aln.sam  \
--min-identity 0.95 \
--min-mq 10 \
--strata \
--paired \
--no-unal \
--rg-id ${flowcell}.${lane}_${sampleID}_${lib}  \
--rg-sm  ${sampleID}\
--rg-lb  ${sampleID}_${lib}\
--rg-pl ${platform}\
--rg-cn ${centre}\
--rg-pu  ${flowcell}.${lane} \
-t $NCPUS 

samtools view -b \
-@ $NCPUS \
/scratch/mo73/Align_split_oryx/${outPrefix}.aln.sam \
> \
/scratch/mo73/Align_split_oryx/${outPrefix}.aln.bam 

java -jar picard.jar SortSam \
I= /scratch/mo73/Align_split_oryx/${outPrefix}.aln.bam \
O= /scratch/mo73/Align_split_oryx/${outPrefix}.sort.bam \
SORT_ORDER=coordinate








