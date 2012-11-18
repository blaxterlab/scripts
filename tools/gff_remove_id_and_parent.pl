#!/usr/bin/env perl

=head1 NAME

gff_remove_id_and_parent.pl

=head1 SYNOPSIS

gff_remove_id_and_parent.pl -g gfffile -i idfile

=head1 DESCRIPTION

Removes lines from GFF file that have IDs in the ID file and removes their parents as well

- Takes a gff file as input (or gzipped)
- Takes an id file (tab delimited or first column)
- Stdout is gff file with lines subtracted (and edited, in case there were multiple parent lines)
- Stderr is gff file with the lines that were removed

=head1 AUTHORS

sujai.kumar@ed.ac.uk 2011.10.10

=cut

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case);
use IO::File;

my $gfffile;
my $idfile;
my $tmpdir = "/tmp/";

GetOptions (
    "gfffile=s" => \$gfffile,
    "idfile=s"  => \$idfile,
    "tmpdir=s"  => \$tmpdir,
);

my $tmpfile = $tmpdir . "gff_remove_id_and_parent.pl." . int(rand(10000));

die "Usage: gff_remove_id_and_parent.pl -g gfffile -i idfile [-t tmpdir]\n"
    unless -r $gfffile;

#-------------------------------
# load idfile into memory

my %id_hash;
my $id_fh = &read_fh($idfile);
while (<$id_fh>) {
    chomp;
    die "File with IDs to be removed should have the ID in the first column\n"
        unless /^([^\t]+)/;
    $id_hash {$1} = 1;
}

my $gff_fh;

#-------------------------------
# make multiple passes through gff file until the $id_hash size stops changing 
# (eg should take 3 passes if there are 3 levels of hierarchy)

my %current_id_hash = %id_hash;
my $i = 1;
do {
    %id_hash = %current_id_hash;
    
    print STDERR "Pass through gff: " . $i++ . "\n";

    $gff_fh = &read_fh($gfffile);
    open TMP, ">$tmpfile" or die "Please specify -t tmpdir if you cannot write to $tmpdir\n";

    while (my $current_gff = <$gff_fh>) {
        #process only if this is a gff 9 col entry, else write to tmpfile
        if ($current_gff !~ /^\S+\t\S+\t\S+\t\d+\t\d+\t\S+\t.\t.\t(.+)/) {
            print TMP $current_gff;
            next;
        }

        my $attributes = $1;
        
        $current_id_hash {$1} = 1 if $attributes =~ /ID=(.+?);/ and exists $id_hash {$1};

        # check each parent (exons can have multiple parents)
        if ($attributes =~ /ID=(.+?);.*?Parent=(.+?);/) {

            my $all_parents_in_id_hash = 1;
            my @curr_parents;    
            for my $parent (split /,/,$2) {

                # if ID of that line is in id_hash, then add the PARENT(s) (if any) to id_hash
                $current_id_hash {$parent} = 1 if exists $id_hash {$1}; 

                # if a parent is not in id_hash, then add to curr_parents list
                if (not exists $id_hash {$parent}) {
                    $all_parents_in_id_hash = 0;
                    push @curr_parents, $parent;
                }

            }

            # if all of the PARENTs of that line are in id_hash, then add the ID to id_hash
            $current_id_hash {$1} = 1 if $all_parents_in_id_hash;
            
            # update parents for this line:
            my $curr_parent_string = join (",", @curr_parents);
            $current_gff =~ s/Parent=(.+?);/Parent=$curr_parent_string;/;
        }
        # print curr_gff to tmpfile;
        print TMP $current_gff;
        
    }
    close $gff_fh;
    close TMP;
} until (scalar keys %id_hash == scalar keys %current_id_hash);

#-------------------------------
# FINAL pass though gff file - 
# if ID of that line is in id_hash, print to STDERR, else print to STDOUT

$gff_fh = &read_fh($tmpfile);
while (<$gff_fh>) {
    if ( /^\S+\t\S+\t\S+\t\d+\t\d+\t\S+\t.\t.\t.*?ID=(.+?);/ and exists $id_hash {$1} ) {
        print STDERR;
        next;
    }
    print;
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

