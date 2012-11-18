#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case);

my $delim = "\t";
my $twoway;
GetOptions (
  "delim:s" => \$delim,
  "twoway" => \$twoway,
);

my %onewaylinks;
my %twowaylinks;
my $count = 0;
while (<>)
{
	chomp;
	my @a = split $delim;
	next unless scalar @a == 2;
	if ($twoway)
	{
		$onewaylinks{$a[0]}{$a[1]} = 1;
	} else
	{
		$twowaylinks{$a[0]}{$a[1]} = $twowaylinks{$a[1]}{$a[0]} = 1;
	}
	print STDERR ($count % 1000000 ? ".":$count) unless $count % 100000;$count++;
}

if ($twoway)
{
	print STDERR "Done making one way links hash\n";
	for my $node1 (keys %onewaylinks)
	{
		for my $node2 (keys %{$onewaylinks{$node1}})
		{
			if ($onewaylinks{$node1}{$node2} and $onewaylinks{$node2}{$node1})
			{
				$twowaylinks{$node1}{$node2} = $twowaylinks{$node2}{$node1} = 1;
			}
		}
	}
}
print STDERR "Done making two way links hash\n";

################################################################
# make clusters based on %twowaylinks

my (%cluster);
my $clusterid = 1;
for my $node (keys %twowaylinks)
{
	next if exists $cluster{$node};
	&visitnode($node); # updates global variable %cluster
	$clusterid++;
}

# each node maps to a clusterid, now map each clusterid to set of nodes
my %cluster2node;
for my $node (keys %cluster) { $cluster2node{$cluster{$node}}{$node} = 1 }

# print contents of each cluster
foreach (sort {scalar keys %{$cluster2node{$b}} <=> scalar keys %{$cluster2node{$a}}} keys %cluster2node) { print join(" ", keys %{$cluster2node{$_}}) . "\n" }

################################################################
################################################################

sub visitnode
{
	my $node = shift @_;
	return if exists $cluster{$node};
	$cluster{$node} = $clusterid;
	foreach (keys %{$twowaylinks{$node}})
	{
		next if exists $cluster{$_};
		&visitnode($_)
	}
}
