#!/bin/bash
#PBS -l wd
#PBS -P mo73 
#PBS -N chr6
#PBS -l walltime=10:00:00,mem=32Gb,ncpus=8,jobfs=1Gb
#PBS -M jose.mijangosaraujo@sydney.edu.au
#PBS -m abe
#PBS -q normal

module load R/3.6.1
module load python3/3.7.4
module load boost/1.71.0
module load samtools/1.9
module load java/jdk-8.40 
module load gatk/4.2.1.0

export R_LIBS=/home/549/jm9807/.R

Rscript octopus_group.R

