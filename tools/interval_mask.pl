#!/usr/bin/perl

# interval_mask.pl takes
# -f fasta file
# and converts to -ci (case_input) if specified
# takes -i interval file which has chrontig st en in three columns (st is 1 based)
# and masks bases in the input with the intervals specified to
# lower/upper case (specified with -co or -case_output)
# and a table to stderr with counts of masked, unmasked bases

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case);

my ($fasta_file, $interval_file);
my ($case_output, $case_input) = ("L", "N");
GetOptions (
  "fastafile:s" => \$fasta_file,
  "intervalfile:s" => \$interval_file,
  "case_output|co:s" => \$case_output,
  "case_input|ci:s" => \$case_input,
);
$case_input  = uc(substr($case_input,0,1));
$case_output = uc(substr($case_output,0,1));
$case_output = "U" if $case_input eq "L";

die "Usage: interval_mask.pl -f <fasta_file> -i <interval_file> [-co {U|L}] [-ci {U|L}]" unless $fasta_file and $interval_file;

my $sequences = &fastafile2hash ( $fasta_file, $case_input ) ;

my $interval_fh = &read_fh ( $interval_file ) ;
while (<$interval_fh>)
{
    # check for errors
    /^>?(\S+)[\t_ ](\d+)[\t_ ](\d+)$/ or do { print STDERR; next };
    exists $$sequences{$1}            or do { print STDERR; next };
    length $$sequences{$1}{seq} >= $3 or do { print STDERR; next };

    if ( $case_output eq "U" )
    {
        substr($$sequences{$1}{seq}, $2-1, $3-$2+1, uc(substr($$sequences{$1}{seq}, $2-1, $3-$2+1)) );
    } elsif ( $case_output eq "L" )
    {
        substr($$sequences{$1}{seq}, $2-1, $3-$2+1, lc(substr($$sequences{$1}{seq}, $2-1, $3-$2+1)) );
    } elsif ( $case_output eq "N" )
    {
        substr($$sequences{$1}{seq}, $2-1, $3-$2+1, "n" x ($3-$2+1) );
    }
}

foreach (sort { $$sequences{$a}{order} <=> $$sequences{$b}{order} } keys %{$sequences} )
{
    print ">$_ $$sequences{$_}{desc}\n$$sequences{$_}{seq}\n";
}

#############################################################################################

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

#############################################################################################

sub fastafile2hash
{
	my $fastafile  = shift @_;
	my $changecase = "N";
	my $order      = "S"; # S = same as input, or R = random
	$changecase    = substr(uc(shift @_),0,1) if @_;
	$order         = substr(uc(shift @_),0,1) if @_;
	my %sequences;
	my $fh = &read_fh($fastafile);
	my $seqid;
	my $seq_counter;
	while (<$fh>)
	{
		if (/^>(\S+)(.*)/) {
		    $seqid = $1;
		    $sequences{$seqid}{desc} = $2;
		    $sequences{$seqid}{order} = $order eq "S" ? $seq_counter++ : rand;
		}
		else {
		    chomp($sequences{$seqid}{seq} .= lc($_)) if $changecase eq "L";
		    chomp($sequences{$seqid}{seq} .= uc($_)) if $changecase eq "U";
		    chomp($sequences{$seqid}{seq} .= $_    ) if $changecase eq "N";
		}
	}
	return \%sequences;
}
