#!/usr/bin/env perl

=head1 NAME

gff_remove_orphans.pl

=head1 SYNOPSIS

gff_remove_orphans.pl gfffile

=head1 DESCRIPTION

- Takes a gff file as input
- Assumes parents/children will be together with parent first (eg: mRNA feature will have exons as children)
- Maintains a hash of all IDs as it goes along.
- If any row has a parent that is not present in the ID table, that row is not output (is sent to STDERR if -e is used)
- Takes a subtract file (contig<tab>start<tab>end) with start and end coords of region to be removed

=head1 AUTHORS

sujai.kumar@ed.ac.uk 2011.08.29

=cut

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case);
use IO::File;

# by default, don't print removed lines to stderr
my $error = '';

GetOptions (
    "error" => \$error,
);

#-------------------------------
# read gff one by one

my %id; #stores all IDs seen in file
while (<>) {
    # if not a gff coord line, print to stdout
    if ( ! /^(\S+)\t(\S+)\t(\S+)\t(\d+)\t(\d+)\t(\S+)\t(.)\t(.)\t(.+)/ ) {
        print;
        next;
    }
    my $attributes = $9;
    $id {$1} ++ if $attributes =~ /ID=(.+?);/;
    print if $attributes !~ /Parent=/;
    print if $attributes =~ /Parent=(.+?);/ and exists $id {$1};
    print STDERR if $attributes =~ /Parent=(.+?);/ and not exists $id {$1} and $error;
}
