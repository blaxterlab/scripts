#!/bin/bash
#$ -cwd -l h_rt=5:59:00
#$ -o .
#$ -e .
#$ -pe memory-2G 1
. /etc/profile
. /exports/work/biology_ieb_mblaxter/software/.softwarerc
if [ -e $1 -a -e $2 ]
then
	java -Xmx1500m -jar $GATKPATH/GenomeAnalysisTK.jar -T UnifiedGenotyper -et NO_ET -out_mode EMIT_ALL_CONFIDENT_SITES -baq CALCULATE_AS_NECESSARY -L ${2%fasta}intervals/$SGE_TASK_ID.intervals -R $2 -o ${1%.bam}.${2%.fasta}.$SGE_TASK_ID.vcf -I $1 
fi
