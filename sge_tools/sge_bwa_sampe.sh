#!/bin/bash
#$ -cwd -l h_rt=5:59:00
#$ -o .
#$ -e .
#$ -pe memory-2G 1
. /etc/profile
. /exports/work/biology_ieb_mblaxter/software/.softwarerc
if [ -e $1 -a -e $2_1.sai -a -e $2_1.sanfastq -a -e $2_2.sai -a -e $2_2.sanfastq ]
then
	bwa sampe -r "@RG\tID:$3\tSM:$3\tPL:Illumina" -f $3.sam $1 $2_1.sai $2_2.sai $2_1.sanfastq $2_2.sanfastq
        samtools view -bS $3.sam -o $3.bam
fi
