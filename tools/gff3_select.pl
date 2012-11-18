#!/bin/env perl

# -g gff3 with embedded fasta as input
# -i include list of contigs wanted in output
# -s split up as separate files with .gff3 extension (one for each contig)

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case no_auto_abbrev);

my ($gff_file, $include_file, $split) = ("-","","");
GetOptions (
  "g|gfffile:s" => \$gff_file,
  "s|split" => \$split,
  "i|inc|include|includefile:s" => \$include_file,
);

#-----------------------

my @includes;

if ($include_file) {
    my $include_fh = &read_fh ( $include_file );
    while (<$include_fh>) {
       /^>?(\S+)/ and push @includes, $1;
    }
}

#-----------------------

my $gff_fh = &read_fh ( $gff_file );

my %contig;

while (<$gff_fh>) {
LN: next if /^#/;
    if (/^>(\S+)/) {
        my $contig_id = $1;
        $contig {$contig_id} {fasta} = $_;
        while ( $_ = <$gff_fh> ) {
            goto LN unless /^[a-z]+$/i;
            $contig {$contig_id} {fasta} .= $_;
        }
    }
    else {
        my @fields = split /\t/;
        next unless @fields == 9;
        $contig { $fields[0] } {gff} .= $_;
    }
}

@includes = sort keys %contig unless $include_file;

my $fh;
open $fh, ">&1";
foreach ( @includes ) {
    open $fh, ">$_.gff3" if $split;
    print $fh "##gff-version 3\n";
    if (!exists $contig {$_}) {
      print STDERR "$_\n";
      next;
    }
    print $fh $contig{$_}{gff};
    print $fh "###\n";
    print $fh "##FASTA\n";
    print $fh $contig{$_}{fasta};
    print $fh "###\n";
}














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
