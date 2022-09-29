#!/bin/bash

#Create input text file for parallel fast splitting

rm -f split_fastq_oryx.input

ls /scratch/mo73/align_redo/*.fq.gz | sed  's/_R1.*\|_R2.*//' | uniq > split_fastq_oryx.input

pairs=`ls /scratch/mo73/align_redo/*.fq.gz | sed  's/_R1.*\|_R2.*//' | uniq | wc -l`

printf "Number of fastq pairs to split: ${pairs}\n"
