#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case);

my $gi_taxid_file = "/exports/work/blast/gi_taxid.dmp.gz";

GetOptions (
    "gi_taxid_file=s" => \$gi_taxid_file,
);

die "Usage: taxid2gid.pl -g /exports/work/blast/gi_taxid.dmp.gz <taxid_list>\n" unless -r $gi_taxid_file;

my %selected_taxids;
while (<>) {
    chomp; $selected_taxids{$_} = 1
}

my $fh = &read_fh($gi_taxid_file);
while (my $line = <$fh>) {
    next unless $line =~ /^(\S+)\s+(\S+)$/;
    print "$1\n" if exists $selected_taxids{$2};
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
