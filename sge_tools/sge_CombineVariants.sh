#!/bin/bash
#$ -cwd -l h_rt=5:59:00
#$ -o .
#$ -e .
#$ -pe memory-2G 1
. /etc/profile
. /exports/work/biology_ieb_mblaxter/software/.softwarerc
if [ -e $2 ]
then
	input=""
	for ((i=1; i<=$(ls $1\.*\.vcf | wc -l);i++))
	do
		ivcf="$1.$i.vcf"
		if [ -e $ivcf ]
		then
			input="$input -B:$i,VCF $ivcf "
		fi
	done
fi
eval "java -Xmx1500m -jar $GATKPATH/GenomeAnalysisTK.jar -T CombineVariants -et NO_ET -assumeIdenticalSamples -R $2 -o $1\.vcf $input" 
