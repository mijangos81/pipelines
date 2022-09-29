#!/bin/bash
#PBS -l wd
#PBS -P mo73 
#PBS -N pipe_1
#PBS -l walltime=06:00:00,mem=16Gb,ncpus=6,jobfs=10Gb
#PBS -M jose.mijangosaraujo@sydney.edu.au
#PBS -m abe
#PBS -q normal

module load R/3.6.1
module load python3/3.7.4
module load bwa/0.7.17
module load samtools/1.9
module load java/jdk-8.40 
module load gatk/4.2.1.0
export R_LIBS=/home/549/jm9807/.R

Rscript index_ref_gen.R