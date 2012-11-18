#!/usr/bin/perl -w

use strict;
use warnings;

die "Usage: pick_random.pl <lines per read> <proportion between 0 to 1> <file(s) with reads>\n" .
    "Examples:\n".
    "To get 20% of the reads from a fasta file:               pick_random.pl 2 0.2 reads.fasta\n".  
    "To get half the reads from a fastq file:                 pick_random.pl 4 0.5 reads.fastq\n".
    "To get half the reads from many interleaved fastq files: pick_random.pl 8 0.5 readsA.fastq readsB.fastq readsC.fastq\n"
unless @ARGV;

my $lines = shift @ARGV;
my $prop  = shift @ARGV;

my ($picked, $total) = (0,0);

while (<>)
{
	my $read = $_;
	$total++;
	for (2..$lines) { $_ = <>; $read .= $_ }
	if (rand() <= $prop) {
		print $read;
		$picked++;
	};
}

print STDERR "Picked $picked chunks of $lines lines each out of total $total chunks\n"
