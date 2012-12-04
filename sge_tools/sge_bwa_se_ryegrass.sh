#!/bin/bash
#$ -cwd -l h_rt=5:59:00
#$ -o .
#$ -e .
#$ -pe memory-2G 1
. /etc/profile
. /exports/work/biology_ieb_mblaxter/software/.softwarerc
if [ -e $1 -a -e $2 ]
then
        readstub=${2%.fastq.gz}
        filename=$(basename "$2")
        samplename=${filename%_1.fastq.gz}
	bwa aln -f $readstub.sai $1 <(gunzip -c $2 | cut -c1-90)
	bwa samse -r "@RG\tID:$samplename\tSM:$samplename\tPL:Illumina" -f $readstub.sam $1 $readstub.sai <(gunzip -c $2 | cut -c1-90)
        samtools view -F4 -bS $readstub.sam -o $readstub.bam
        samtools sort $readstub.bam ${readstub}_sorted
        rm $readstub.sai $readstub.sam $readstub.bam
fi
