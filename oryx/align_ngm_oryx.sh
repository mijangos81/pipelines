#!/bin/bash
#PBS -l wd
#PBS -P mo73 
#PBS -N align_NGM_oryx1
#PBS -l walltime=08:00:00,mem=64Gb,ncpus=16,jobfs=4Gb
#PBS -W umask=022
#PBS -M jose.mijangosaraujo@sydney.edu.au
#PBS -m abe
#PBS -q normal

module load parallel/20191022

set -e

SCRIPT=/home/549/jm9807/align_oryx.sh  #Script to run
INPUTS=/home/549/jm9807/1.txt   # Each line in this file is used as arguments to ${SCRIPT}. Line includes comma delimited list of arguments

NCPUS=4 #cpus per task
sed -i "s/$/,$NCPUS/" $INPUTS #Update the inputs file to include the cpus per task (so the sh script does not have to be manually edited)

#mkdir -p /scratch/mo73/Align_split_oryx

parallel -j $((${PBS_NCPUS}/${NCPUS})) --rpl '{%} 1 $_=($job->slot()-1)'"*${NCPUS}" pbsdsh -n {%}  -- bash -l -c "'${SCRIPT} {}'" :::: ${INPUTS}
