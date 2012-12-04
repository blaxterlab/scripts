#!/bin/bash
#$ -cwd -l h_rt=5:59:00
#$ -o .
#$ -e .
#$ -pe memory-2G 1 -R y
. /etc/profile
. /exports/work/biology_ieb_mblaxter/software/.softwarerc
if [ -e "/exports/work/scratch/jdavey/WGS.$1.Hmel1-1_primaryScaffolds.$SGE_TASK_ID.sam" ]
then
        rm /exports/work/scratch/jdavey/WGS.$1.Hmel1-1_primaryScaffolds.$SGE_TASK_ID.sam
fi
stampy.py --substitutionrate=0.01 --gatkcigarworkaround --baq  --alignquals --processpart $SGE_TASK_ID/$2 --readgroup=ID:WGS.$1,SM:WGS.$1,PL:Illumina -g ./Hmel1-1_primaryScaffolds -h ./Hmel1-1_primaryScaffolds -o /exports/work/scratch/jdavey/WGS.$1.Hmel1-1_primaryScaffolds.$SGE_TASK_ID.sam -M /exports/work/scratch/jdavey/$1.R1.fastq.gz /exports/work/scratch/jdavey/$1.R2.fastq.gz
samtools view -Sb -o /exports/work/scratch/jdavey/WGS.$1.Hmel1-1_primaryScaffolds.$SGE_TASK_ID.bam /exports/work/scratch/jdavey/WGS.$1.Hmel1-1_primaryScaffolds.$SGE_TASK_ID.sam
