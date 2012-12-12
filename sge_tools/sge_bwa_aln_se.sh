#!/bin/bash
#$ -cwd -l h_rt=5:59:00
#$ -o .
#$ -e .
#$ -pe memory-2G 1
. /etc/profile
. /exports/work/biology_ieb_mblaxter/software/.softwarerc
if [ -e $1 -a -e $2 ]
then
	bwa aln -f ${2%fastq.gz}sai $1 <(gunzip -c $2)
fi
