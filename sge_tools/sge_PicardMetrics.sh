#!/bin/bash
#$ -cwd -l h_rt=5:59:00
#$ -o .
#$ -e .
#$ -pe memory-2G 1
. /etc/profile
. /exports/work/biology_ieb_mblaxter/software/.softwarerc
if [ -e $1.bam -a -e $2 ]
then
	java -Xmx1500m -jar $PICARDPATH/CollectAlignmentSummaryMetrics.jar TMP_DIR=$TMPDIR VALIDATION_STRINGENCY=LENIENT MAX_INSERT_SIZE=5000 REFERENCE_SEQUENCE=$2 INPUT=$1.bam OUTPUT=$1.CollectAlignmentSummary.metrics
	java -Xmx1500m -jar $PICARDPATH/CollectGcBiasMetrics.jar TMP_DIR=$TMPDIR VALIDATION_STRINGENCY=LENIENT REFERENCE_SEQUENCE=$2 INPUT=$1.bam OUTPUT=$1.CollectGcBias.metrics CHART_OUTPUT=$1.CollectGcBias.pdf SUMMARY_OUTPUT=$1.CollectGcBias.summary
	java -Xmx1500m -jar $PICARDPATH/CollectInsertSizeMetrics.jar TMP_DIR=$TMPDIR VALIDATION_STRINGENCY=LENIENT REFERENCE_SEQUENCE=$2 INPUT=$1.bam OUTPUT=$1.CollectInsertSize.metrics HISTOGRAM_FILE=$1.CollectInsertSize.pdf MINIMUM_PCT=0
fi
