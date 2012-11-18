#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case);

my $mapfile;
GetOptions (
  "file:s" => \$mapfile,
);

die "\nUsage: freplace.pl -f mapfile STDIN\n\nDescription: mapfile has a mapping from old to new charstrings. freplace.pl looks for old char string and replaces it with new one if they are unambiguous, else it prints ambiguous lines to stderr\n\n" unless $mapfile;

my %map;
my $fh = &read_fh($mapfile);

#load mapfile into memory
while (<$fh>) {
    chomp;
    @_ = split /\t/;
    next unless scalar @_ == 2;
    $map{$_[0]} = $_[1];
}

#process main stream
while (<>) {
    chomp;
    while (/([\=>\s;:,]*)([^\=>\s;:,]+)([\=>\s;:,]*)/g) {
        print $1 . ((exists $map{$2}) ? $map{$2} : $2) . $3;
    }
    print "\n";
}

############################################################
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

