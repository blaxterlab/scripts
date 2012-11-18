#!/usr/bin/perl -w

use strict;
use warnings;

my  $usage = "taxid_parents_list.pl <file_with_list_of_taxids>\n";
die $usage if @ARGV == 0;

open NODES, "</exports/work/blast/nodes.dmp" or die $usage;
my %parentof;
my %taxlevel;
while (my $line = <NODES>)
{
	next if $line !~ /^(\d+)\s*\|\s*(\d+)\s*\|\s*(\w+)/;
	next if $1 == $2;	# if any node has itself as a parent (eg taxid 1) then the recursion doesn't stop and the script falls over
	$parentof{$1} = $2;
	$taxlevel{$1} = $3;
}
close NODES;

open NAMES, "</exports/work/blast/names.dmp" or die $usage;
my %taxid2name;
while (my $line = <NAMES>)
{
	next unless $line =~ /^(\d+)\s*\|\s*(.+?)\s*\|\s*.+?\s*\|\s*scientific name/;
	$taxid2name {$1} = $2;
}
close NAMES;

my $taxid;
while ($taxid = <>)
{
	chomp($taxid);
	while (exists $taxid2name{$taxid}) {
		print "$taxid2name{$taxid} | ";
		if (exists $parentof{$taxid}) {
			$taxid = $parentof{$taxid}
		} else { last }
	}
        print "\n"
}
