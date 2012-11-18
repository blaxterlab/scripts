#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my @files;
my $delimiter = "\t";
my $column = 1;
my $unpair;
GetOptions (
    "delimiter=s" => \$delimiter,
    "unpair=s" => \$unpair,
    "column=i" => \$column,
    "files=s{,}" => \@files,
);

#---------------------------------

die <<USAGE

  Usage: join.pl [-d \"\\t\"] [-u string] [-c 1] -f file1 file2 [...]
  Joins delimited files on first column (use -c to change)
  Tab delimited by default. Use -d to change
  Only prints rows where all files have an entry unless -u is given in which case the string given is inserted instead of every missing data. 
  Only takes the first row if a key is present in multiple rows
  
USAGE
unless @files;

#---------------------------------
# load all files into memory

my %filehashes;

for my $file (@files)
{
    my %filehash;
    open FH, "<$file" or die $!;
    while (<FH>)
    {
        chomp;
        my @row = split /$delimiter/;
        my $row_key = $row[$column -1];
        @row = @row[1..$#row] if $column==1;
        $filehash{$row_key} = \@row unless exists $filehash{$row_key}; # ensures only first entry per key
    }
    $filehashes{$file} = \%filehash;
}

#---------------------------------
# get all keys for first file and use that to search all the others

my $testfile = shift @files; # the first file is now removed so 

for my $key (keys %{$filehashes{$testfile}})
{
    my @toprint = ($key);
    push @toprint, @{$filehashes{$testfile}->{$key}};
    my $found_in_all = 1;
    for my $file (@files)
    {
        if (not defined $filehashes{$file}->{$key})
        {
            if (not defined $unpair)
            {
                $found_in_all = 0;
                last;
            }else
            {
                push @toprint, $unpair;
            }
        }
        else
        {
            push @toprint, @{$filehashes{$file}->{$key}};
        }
    }
    if ($found_in_all)
    {
        print join ($delimiter, @toprint) . "\n" if $found_in_all;
    }
}
