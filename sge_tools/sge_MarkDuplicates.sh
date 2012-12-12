#!/bin/bash
#$ -cwd -l h_rt=5:59:00
#$ -o .
#$ -e .
#$ -pe memory-2G 2
. /etc/profile
. /exports/work/biology_ieb_mblaxter/software/.softwarerc
if [ -e $1 ]
then
	java -Xmx3500m -jar $PICARDPATH/MarkDuplicates.jar TMP_DIR=$TMPDIR REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 MAX_RECORDS_IN_RAM=50000 INPUT=$1 OUTPUT=${1%bam}rmdup.bam METRICS_FILE=${1%bam}MarkDuplicates.metrics
	samtools index ${1%bam}rmdup.bam
fi
