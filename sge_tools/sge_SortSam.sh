#!/bin/bash
#$ -cwd -l h_rt=5:59:00
#$ -o .
#$ -e .
#$ -pe memory-2G 4
. /etc/profile
. /exports/work/biology_ieb_mblaxter/software/.softwarerc
if [ -e $1.bam ]
then
	java -Xmx7G -jar $PICARDPATH/SortSam.jar TMP_DIR=$TMPDIR MAX_RECORDS_IN_RAM=100000 SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT INPUT=$1.bam OUTPUT=$1_sorted.bam
fi
