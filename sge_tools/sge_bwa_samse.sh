#!/bin/bash
#$ -cwd -l h_rt=5:59:00
#$ -o .
#$ -e .
#$ -pe memory-2G 1
. /etc/profile
. /exports/work/biology_ieb_mblaxter/software/.softwarerc
if [ -e $1 -a -e $2.sai -a -e $2.sanfastq ]
then
	bwa samse -r "@RG\tID:$3\tSM:$3\tPL:Illumina" -f $3.sam $1 $2.sai $2.sanfastq
        samtools view -bS $3.sam -o $3.bam
fi
