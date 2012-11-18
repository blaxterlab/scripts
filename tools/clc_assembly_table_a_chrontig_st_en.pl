#!/usr/bin/perl

use strict;
use warnings;

# takes clc's assembly_table -a file.cas output and gives reference CONTIGNAME START END

while (my $line = <>) {
    next unless $line =~ /^\s*(\S+)\s.+? has (\d+) match/;
    if ($2 == 0) {
        print STDERR "$1\n";
        next;
    }
    $line = <>;
    $line = <>;
    die "Unexpected input\n" unless $line =~ /(\d+)\s+(\S+)\s+(\d+)\s+(\S+)/;
    print "$4\t$1\t$3\n";
}
