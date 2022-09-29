#!/bin/bash

#Create input list for parallel alignment with BWAkit
#Read sample info from samples.config
#Unsplit fastq in Fastq directory, split fastq in Fastq_split
#Assumes flowcell ID is field 3 of ':' delim read ID
#Assumes flowcell lane number is field 4 of ':' delim read ID

rm -f /home/549/jm9807/align.input.txt

awk 'NR>1' /home/549/jm9807/samples.config_2.txt | while read LINE
do 
	sample=`echo $LINE | cut -d ',' -f 1`
	labSampleID=`echo $LINE | cut -d ',' -f 2`
	centre=`echo $LINE | cut -d ',' -f 3`
	lib=`echo $LINE | cut -d ',' -f 4`
	platform=illumina
	
	if [ ! "$lib" ]
	then
	    	lib=1
	fi
	
	fqpairs=$(ls /scratch/mo73/Fastq_split/*${sample}*fastq.gz | sed  's/_R1.*\|_R2.*//' | uniq)
	fqpairs=($fqpairs)
	
	#Get the flowcell and lane info from the original pairs:	
	for ((i=0; i<${#fqpairs[@]}; i++))
	do
		flowcell=$(zcat ${fqpairs[i]}_R1.fastq.gz  | head -1 | cut -d ':' -f 3)
		lane=$(zcat ${fqpairs[i]}_R1.fastq.gz | head -1 | cut -d ':' -f 4)
		
		#Print each of the split chunks with flowcell and lane info to inputs file:
		set=$(basename ${fqpairs[i]})
		splitpairs=$(ls /scratch/mo73/Fastq_split/*${set}* | sed  's/_R1.*\|_R2.*//' | uniq)
		splitpairs=($splitpairs)
		
		for ((c=0; c<${#splitpairs[@]}; c++))
		do
			printf "${splitpairs[c]},${sample},${centre},${lib},${platform},${flowcell},${lane}\n" >> /home/549/jm9807/align.input.txt
		done			
	done			
done	

tasks=`wc -l < /home/549/jm9807/align.input.txt`
printf "Number of alignment tasks to run: ${tasks}\n"
	
