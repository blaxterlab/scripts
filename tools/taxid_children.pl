#!/usr/bin/perl -w

=pod

=head1 NAME

taxid_children.pl - Perl script that takes an NCBI taxonomy ID and returns a list of all taxids that are children of that taxid

=head1 SYNOPSIS

Usage: taxid_children.pl <taxonomy_id> <nodes.dmp>

=head1 DESCRIPTION

Takes a given taxid (eg Nematoda is 6231) and the nodes.dmp file that stores taxonomy child-parent relationship, and prints a list of all taxids that are descendants of the given taxid, including itself.
Eg, if you want all taxids in the phylum Nematoda, run:
taxid_children.pl 6231 /home/blastdb/nodes.dmp >Nematoda.taxid
To get all taxids other than Nematoda (eg to check for contamination), run:
cut -f 1 /home/blastdb/nodes.dmp | fgrep -v -f Nematoda.taxid >ALLnNematoda.taxid 
(cut -f 1 gets all the taxids in the first col in nodes.dmp, grep -v -f returns all the taxids that are NOT in Nematoda.taxid)

=head1 ARGUMENTS

Arg [1] : taxnonomy_id (Integer). If it doesn't exist in nodes.dmp then nothing is printed (no error is given)

Arg [2] : nodes.dmp (filepath). SQL dump of NCBI's taxonomy nodes table. Extracted from ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

=head1 SEE ALSO

See the team wiki page http://scratchy.internal.sanger.ac.uk/wiki/index.php/Team_133

=head1 AUTHOR

Sujai Kumar, <sk13@sanger.ac.uk

=head1 TIME AND DATE

200908241030

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2009 Genome Research Limited. All Rights Reserved.              
                                                                               
This program is free software; you can redistribute it and/or                 
modify it under the terms of the GNU General Public License                   
as published by the Free Software Foundation; either version 2                
of the License, or (at your option) any later version.                        
                                                                               
This program is distributed in the hope that it will be useful,               
but WITHOUT ANY WARRANTY; without even the implied warranty of                
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                 
GNU General Public License for more details.                                  
                                                                               
You should have received a copy of the GNU General Public License             
along with this program; if not, write to the Free Software                   
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA. 

=cut

#code goes here

############################################################

use strict;
use warnings;

my $taxid = shift @ARGV;
my $nodefile = shift @ARGV;

# die "Usage: taxid_children.pl <taxonomy_id> <nodes.dmp>" unless @ARGV == 2;

open NODES, "<$nodefile" or die "Could not open nodes.dmp at the location specified";

## create hashes
my %childof;

while (my $line = <NODES>)
{
	# line in nodes.dmp should match the regexp below.
	# Change the regexp if NCBI changes their file format
	next if $line !~ /^(\d+)\s*\|\s*(\d+)/;
	next if $1 == $2; # if any node has itself as a parent (eg taxid 1) then the recursion doesn't stop and the script falls over
	push @{$childof{$2}}, $1;
}

&print_children($taxid);

############################################################

sub print_children 
{
	my $current_id = shift @_;
	print "$current_id\n";
	if (exists $childof{$current_id})
	{
		for my $child ( @{$childof{$current_id}} )
		{
			&print_children($child)
		}
	}
}
