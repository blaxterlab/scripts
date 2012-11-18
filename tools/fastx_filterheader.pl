#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case);

my ($fastxfile,$regexp,$includefile,$excludefile);
my $fastxlines = 2; # default for next gen fasta files
GetOptions (
  "fastxfile:s" => \$fastxfile,
  "regexp:s" => \$regexp,
  "includefile:s" => \$includefile,
  "excludefile:s" => \$excludefile,
  "lines:i" => \$fastxlines,
);

die "Usage: fastx_filterheader.pl -f <fastxfile> [-i <includeseqs>] [-e <excludeseqs>] [-r <regexp>] [-l <2 or 4 depending on fasta or fastq, default 2>]\n" unless $fastxfile;

my (%include_headers, %exclude_headers);
if ($includefile) {
    my $fh = &read_fh($includefile);
    while (<$fh>) {
        chomp;
        $include_headers{$1}=1 if /^>?(\S+)/;
    }
    close $fh;
}
if ($excludefile) {
    my $fh = &read_fh($excludefile);
    while (<$fh>) {
        chomp;
        $exclude_headers{$1}=1 if /^>?(\S+)/;
    }
    close $fh;
}

my $fh = &read_fh($fastxfile);
while (<$fh>) {
    next unless /^[>@](\S+)/;
    my $header = $1;
    my $toprint = $_;
    next if $excludefile and exists $exclude_headers{$header};
    next if $includefile and not exists $include_headers{$header};
    for (2..$fastxlines) {
        $toprint .= <$fh>;
    }
    next if $regexp and $toprint !~ /$regexp/;
    print $toprint;
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
