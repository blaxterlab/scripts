#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case);

my ($remove_intervals_file, $blast_file) = ("","");
GetOptions (
  "remove:s" => \$remove_intervals_file,
  "blast:s" => \$blast_file,
);

die "Usage: remove_overlapping_blast.pl -r remove_intervals_file -b blast_file" unless $remove_intervals_file and $blast_file;

#--------------------------------------------------------

my %intervals;

my $remove_intervals_fh = &read_fh ($remove_intervals_file);
while (<$remove_intervals_fh>)
{
    push @{$intervals{$1}}, [$2,$3] if /^(\S+)\t(\d+)\t(\d+)/;
}
for my $chrontig (keys %intervals)
{
    @{$intervals{$chrontig}} = sort { $a->[0] <=> $b->[0] } @{$intervals{$chrontig}};
}

#--------------------------------------------------------

my $blast_fh = &read_fh ($blast_file);
while (<$blast_fh>)
{
    chomp;
    next if /^#/;
    next unless /^(.+?)\t(.+?)\t(.+?)\t(.+?)\t(.+?)\t(.+?)\t(.+?)\t(.+?)\t(.+?)\t(.+?)\t(.+?)\t(.+)$/;
    my ($qid, $tid, $perc_identity, $length, $gaps, $mismatches, $qst, $qen, $tst, $ten, $evalue, $score) = ($1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12);
    ($qst, $qen) = ($qen, $qst) if $qst > $qen;
    ($tst, $ten) = ($ten, $tst) if $tst > $ten;
    if ($qid =~ /^(.+?)_(\d+)_(\d+)$/) { $qst = $qst + $2 - 1; $qen = $qen + $2 - 1; $qid = $1 }
    if ($tid =~ /^(.+?)_(\d+)_(\d+)$/) { $tst = $tst + $2 - 1; $ten = $ten + $2 - 1; $tid = $1 }

    next if exists $intervals{$qid} and scalar(&binary_search_intervals ( $intervals{$qid}, $qst, $qen ));
    next if exists $intervals{$tid} and scalar(&binary_search_intervals ( $intervals{$tid}, $tst, $ten ));
    
    print $_ . "\n";
}

#############################################################################

sub read_fh {
    my $filename = shift @_;
    my $filehandle;
    if ($filename =~ /gz$/) {
        open $filehandle, "gunzip -dc $filename |" or die $!;
    }
    else {
        open $filehandle, "<$filename" or die $!;
    }
    return $filehandle;
}
#############################################################################

sub binary_search_intervals
{
    my ($array_ref, $range_st, $range_en) = @_;
    my @sorted = @$array_ref;
    
    my $low_index = 0; my $high_index = $#sorted;
    my $found = 0;
    while ( $low_index <= $high_index and not $found )
    {
        my $mid_index = int ( ( $low_index + $high_index ) / 2 );
        if ( $sorted[$mid_index]->[0] <= $range_en and $sorted[$mid_index]->[1] >= $range_st )
        {
            $found = 1;
            return $mid_index + 1;
        }
        if ( $sorted[$mid_index]->[0] < $range_en )
        {
            $low_index = $mid_index + 1;
        } else
        {
            $high_index = $mid_index - 1;
        }
    }
    return $found;
}