#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case);

#globals
my $bam;
my $type = 'fastq';
my $format = "i";
my $prefix = "";

GetOptions (
  "bam:s" => \$bam,
  "type:s" => \$type,
  "format:s" => \$format,
  "out:s" => \$prefix,
);
die "Usage:\n\tbam2fastx.pl -b bam file \nOptional:\n\t-t output type fasta or fastq (default = fastq)\n\t-o output file prefix (default is input file name)\n\t-f output format i = interleave, s = split, b = both (default = i)\n"
	unless ($bam && ($type eq 'fasta' || $type eq 'fastq') && ($format eq 'i' || $format eq 's' || $format eq 'b'));

#set the prefix
if ($prefix eq ""){
	$prefix = $bam;
}	
#sort the file by read name and then read
print "Sorting bam file...\n";
`samtools sort -m 2000000000 -n $bam $bam.readsorted`;
print "Parsing bam file...\n";
open(B, "samtools view $bam.readsorted.bam |") or die;

#get file type
my $suff;
if ($type eq "fastq"){
	$suff = "fastq";
}else{
	$suff = "fna";
} 
#open files and create array of pair files for zipping later
my @outfiles;
if ($format eq 'i' || $format eq 'b'){
	push (@outfiles,"$prefix.pairs.$suff");
	open P,">$prefix.pairs.$suff";
}if ($format eq 's' || $format eq 'b'){
	push (@outfiles,"$prefix.pair1.$suff");
	push (@outfiles,"$prefix.pair2.$suff");
	open P1,">$prefix.pair1.$suff";
	open P2,">$prefix.pair2.$suff";
}
push (@outfiles,"$prefix.singletons.$suff");
open S,">$prefix.singletons.$suff";

my $head1=""; my $seq1=""; my $qual1=""; 
my $head2=""; my $seq2=""; my $qual2="";
my %paired;
while (<B>){
	chomp;
	my @a = split('\t');
	$head2=$a[0]; $qual2 = $a[10]; $seq2 = $a[9];
	#check for bam files with no quality information and create fake one
	if ($qual2 eq "*"){$qual2 = "I" x length($seq2);} 
	if ($head2 eq $head1){
		%paired=();
		#print fastq pairs
		if ($type eq 'fastq'){
			if ($format eq 'i' || $format eq 'b'){
				print P "@".$head1."/1\n".$seq1."\n+".$head1."/1\n".$qual1."\n";
				print P "@".$head2."/2\n".$seq2."\n+".$head2."/2\n".$qual2."\n"
			}if ($format eq 's' || $format eq 'b'){
				print P1 "@".$head1."/1\n".$seq1."\n+".$head1."/1\n".$qual1."\n";
				print P2 "@".$head2."/2\n".$seq2."\n+".$head2."/2\n".$qual2."\n"
			}
		#print fasta pairs	
		}else{
			if ($format eq 'i' || $format eq 'b'){
				print P ">".$head1."/1\n".$seq1."\n";
				print P ">".$head2."/2\n".$seq2."\n";
			}if ($format eq 's' || $format eq 'b'){
				print P1 ">".$head1."/1\n".$seq1."\n";
				print P2 ">".$head2."/2\n".$seq2."\n";
			}
		}
		$paired{$head1}="";
		next;
	#print singletons	
	}elsif ($head1 ne "" && !exists $paired{$head1}){
		if ($type eq 'fastq'){
			print S "@".$head1."/1\n".$seq1."\n+".$head1."/1\n".$qual1."\n";
		}else{
			print S ">".$head1."/1\n".$seq1."\n";
		}
		%paired=();
	}
	$head1=$a[0]; $qual1 = $a[10]; $seq1 = $a[9];
	#check for bam files with no quality information and create fake one 
	if ($qual1 eq "*"){$qual1 = "I" x length($seq1);}
}	
#there might be a singleton right at the very end!
if (!exists $paired{$head1}){
	if ($type eq 'fastq'){
		print S "@".$head1."/1\n".$seq1."\n+".$head1."/1\n".$qual1."\n";
	}else{
		print S ">".$head1."/1\n".$seq1."\n";
	}
}

#`rm $bam.readsorted.*.bam`;
#print "array = @outfiles\n";
print "Zipping up...\n";
`parallel gzip -f ::: @outfiles`;