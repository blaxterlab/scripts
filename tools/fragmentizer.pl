#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case);
# This creates $fragment_len long fragments of contigs
use Data::Dumper;

my $fragment_len ;
my $len_threshold ;
my $assembly_file = '';
my $verbose;

GetOptions (
    "file=s" => \$assembly_file,
    "threshold=i" => \$len_threshold,
    "len=i" => \$fragment_len,
    "verbose" => \$verbose
);
$len_threshold = (0.5) unless $len_threshold;
$fragment_len = (1000) unless $fragment_len;

die <<USAGE
Usage: fragmentizer.pl -f assembly_file.fa [-l] [-t]
-f Assembly file to be split
-l length at which to split the sequences [default: 1000]
-t percentage of length at which the end of the sequence gets split into a new sequence [default: 0.5].
   e.g at l = 1000, t = 0.5 : a 1500 nt long sequence does not get split and a 1501 nt long sequence gets split into two sequences (1000nt and 501nt) 
-v verbose
# If the program encounters more than 10 consecutive N's it replaces them with 10 N's 
USAGE
unless ($assembly_file ne '');

open IN, "<$assembly_file" || die "Can't read\n";

my $header = '';
my $seq    = '';
my %hash;
my @array;

while ( my $line = <IN> ) {
    if ( $line =~ /^>(\S+)\n/ ) {
        if ( $header ne '' ) {
            $seq =~ s/N{11,}/NNNNNNNNNN/g; # replace ocurrence of more than 10 N's with 10 N's
            push @array,
                { 'header' => $header, 'seq' => $seq, 'len' => length($seq) };
        }
        $header = $1;
        $seq    = '';
    }
    else {
        chomp $line;
        $seq .= $line;
    }
}
$seq =~ s/N{11,}/NNNNNNNNNN/g;
push @array, { 'header' => $header, 'seq' => $seq, 'len' => length($seq) };
close IN;
my $out_file = $assembly_file . "_" . $fragment_len . ".fa";

# command line args and usage

open OUT, ">$out_file" || die "Can't write\n";
for my $i ( 0 .. $#array ) {
    my $header                = $array[$i]{'header'};
    my $seq                   = $array[$i]{'seq'};
    my $len                   = $array[$i]{'len'};
    my $modulo                = $len % $fragment_len;
    my $number_of_full_chunks = ( $len - $modulo ) / $fragment_len;
    my @subsequences          = ( $seq =~ /(.{1,$fragment_len})/g );
    my $print_string = $header." ".$len."\n";
    if ( $number_of_full_chunks == 0 ) {
        print OUT ">" . $header . "_1\n" . $subsequences[0] . "\n";
        $print_string .= "_1\t".length($subsequences[0])."\n";
        next;
    }
    if ( $modulo > 0 && ( $modulo / $fragment_len ) <= $len_threshold ) {
        $subsequences[-2] .= pop @subsequences;
    }
    for my $j ( 0 .. $#subsequences ) {
        print OUT ">"
            . $header . "_"
            . ( $j + 1 ) . "\n"
            . $subsequences[$j] . "\n";
        $print_string .= ($j+1)."\t".length($subsequences[$j])."\n";
    }
    if ($verbose){
        print $print_string;
    }
}
close OUT;
