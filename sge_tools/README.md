# Sun Grid Engine scripts

These scripts are mostly designed to run various tools in parallel on a Sun Grid Engine cluster, namely Eddie, the Edinburgh compute cluster. Many of the shorter scripts are best thought of as templates that will probably need rewriting for actual use, depending on the file structure and filename conventions of a particular project. A couple of the scripts are more robust, or are pipelines to string multiple tools together.

## Crude shell scripts

These tools run one tool only in a very crude fashion.

### BWA/Samtools tools

**sge_bwa_aln_se.sh** - run bwa aln on single end data  
**sge_bwa_sampe.sh** - run bwa sampe  
**sge_bwa_samse.sh** - run bwa samse  
**sge_flagstat.sh** - run samtools flagstat

### GATK/Picard tools

The following tools run one component of GATK or Picard. All hardcode the number of slots requested and the amount of memory allocated. Both of these values should be turned into parameters as these values often vary.

**sge_CombineVariants.sh**  
**sge_MarkDuplicates.sh**  
**sge_PicardMetrics.sh**  
**sge_SortSam.sh**    
**sge_UnifiedGenotyper.sh**  

### Others

**sge_blast** - Splits FASTA files and runs blast jobs as an array  
**sge_iprscan** - Splits FASTA files and runs InterProScan 5 jobs as an array  
**sge_stampy** - run Stampy as a job array, splitting the input by feeding $SGE_TASK_ID to --processpart. Convert SAM output to BAM afterwards using **samtools view**.
**jobsum.sh** - run some crude accounting on a completed array job to peek at the cpu and memory usage

## Pipelines

### sge_cegma

**sge_cegma** generates BLAST output for CEGMA with **sge_blast** and then runs CEGMA on this BLAST output. The script **sge_cegma_blast** only runs the BLAST, if CEGMA itself is to be run on a local machine. (It would be more elegant to add an option to the *sge_cegma* script to end after the BLAST.)

### sge_bwa_se_ryegrass

This script runs BWA from raw reads to sorted BAM file, single end: bwa aln, bwa samse, samtools view to convert to BAM, samtools sort. It was written to process ryegrass RAD data and will crop reads to 90bp. Please rewrite before using on another project.

### sge_align

This script runs Stampy, GATK MarkDuplicates, GATK IndelRealigner and GATK UnifiedGenotyper on multiple genome references, samples and FASTQ files. See header for brief description and code for options. Contact [John Davey](johnomics@gmail.com) for more details.



