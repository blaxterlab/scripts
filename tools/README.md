bam2fastx.pl
=

Extracts FASTA or FASTQ files from a bam file, with options to get output as pairs, singles or both

bigmem_blat.pl
=
Runs blat in parallel on a multi-core machine.

interval_mask.pl
=
Takes fasta file and interval file as inputs.
Masks intervals in fasta file as upper or lower case
  * fasta_extract_lc.pl
  * fasta_extract_uc.pl

blast_separate_taxa.pl
=
Takes two sets of blast tabular results.
Separates contigs into two separate files, too close to call (min bit
score difference) into a third file

blast_taxonomy_report.pl
=
no description

blastm8_filter.pl
=
handy for checking "completeness" etc (eg tblastn protein seq against test genome, how many proteins are present
at 70% completeness...)

bowtie2_extract_reads_mapped_to_specific_contigs.pl
=
no description

clc_assembly_table_a_chrontig_st_en.pl
=
converts clc assembly_table -a full text alignments into a 3 column tabular format: reference
CONTIGNAME START END

donkey.pl / scaffold_stats.pl
=
use instead of contig_stats.pl

fastaqual_multiline_to_singleline.pl
=
no description

fastaqual_select.pl
=
no description

fastx_filterheader.pl
=
for pulling out reads from a fasta or fastq file based on regexp or include lists, or exclude lists

fastx_separate.pl
=
Usage: -l [2|4] -p "pattern1" -p "pattern2" -p "pattern3" <STDIN>

fgrep.pl
=
same as fgrep, but written in perl. Search for fgrep.pl in thesis github

freplace.pl
=
Uses a map file (old_value<tab>new_value) - replaces all instances of old_value with new_value in STDIN stream (delimited
by spaces/punctuation). To do possibly - add -d delimiter option like in fgrep.pl

gff_remove_id_and_parent.pl
=
Check exact description

gff_remove_orphans.pl
=
no description

gff3_select.pl
=
outputs gff3 with only those contigs that are specified in a file

join.pl
=
Test it. Check with tim. Compare to unix join

nclust.pl
=
single linkage clustering. but, consider using mcl or some other clever clusterer

pick_random.pl
= picks blocks of lines (4 for fastq, 8 for paired fastq, 2 for single line fasta, etc)

remove_overlapping_blast.pl
=
Removes blast tabular rows if they overlap given intervals

seq_st_en_merge_overlapping.pl
=
Merges sequence intervals (reorders so that st < en). Starts are 1-based not 0-based

taxid_parents_list.pl
=
no description

taxid_children.pl
=
no description

taxid2gid.pl
=
For each taxon id, get all the gid (genbank id) from the gi_taxid dmp files. Look up
https://www.wiki.ed.ac.uk/display/BlaxterLab/Local+blast+against+subsets+of+NCBI+databases+based+on+taxid
for how to use
