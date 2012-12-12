#!/usr/bin/env perl

# sge_align
#
# Generate SNP calls from raw FASTQ files and reference to VCF files
# using an SGE cluster, allowing selection of steps, multiple lanes per
# sample and multiple samples per species
#
# Takes config file specifying FASTQ files, reference files, sample, lane and
# species IDs, and stages to run. Individual references, samples, lanes and
# species can be chosen at command line. Stages can also be specified at
# command line.
#
# Requires a configuration file with a line for each lane of sequencing
# in this format:
# Species Reference SampleX LaneX FileX1 (FileX2)
# where Species is an ID for the species in question,
#       Reference is the name of the reference sequence,
#       SampleX is the name of the sample for the read group,
#       LaneX is the name of the lane for the read group,
#       And FileX1 & optional FileX2 are the names of the FASTQ read files.
# Standard options for Stampy, GATK baked in to script

# Author: John Davey john.davey@ed.ac.uk
# Begun 1/7/11

#############################################################################
###                                                                       ###
###                                 CODE                                  ###
###                                                                       ###
#############################################################################

use strict;
use warnings;
use Carp;
use English;
use Getopt::Long;
use Data::Dumper;
use File::Basename;
use Cwd qw(getcwd abs_path);

eval {
    require Parallel::ForkManager;
    Parallel::ForkManager->import();
};
croak "Can't find Parallel::ForkManager: $OS_ERROR\n"
  if $@;

# Autoflush output so reporting on progress works
$| = 1;

my $in_filename = "";
my $tmp_dir     = "\$TMPDIR";    # Will use Eddie's temp dir
my $output_dir  = getcwd();

my %run;
$run{stampy}              = "stampy.py";
$run{stampy_options}      = "--gatkcigarworkaround --baq  --alignquals";
$run{stampy_user_options} = "";
$run{picard_sort}         = "java -Xmx1500m -jar \$PICARDPATH/SortSam.jar";
$run{picard_sort_options} =
  "TMP_DIR=$tmp_dir MAX_RECORDS_IN_RAM=50000 SORT_ORDER=coordinate";
$run{picard_merge} = "java -Xmx7g -jar \$PICARDPATH/MergeSamFiles.jar";
$run{picard_merge_options} =
  "TMP_DIR=$tmp_dir CREATE_INDEX=true MAX_RECORDS_IN_RAM=50000";
$run{picard_rmdup} = "java -Xmx2500m -jar \$PICARDPATH/MarkDuplicates.jar";
$run{picard_rmdup_options} =
"TMP_DIR=$tmp_dir REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 MAX_RECORDS_IN_RAM=50000";
$run{picard_dict} =
  "java -Xmx1500m -jar \$PICARDPATH/CreateSequenceDictionary.jar";
$run{gatk_ug} =
"java -Xmx14g -jar \$GATKPATH/GenomeAnalysisTK.jar -T UnifiedGenotyper";
$run{gatk_ug_options} =
  "-out_mode EMIT_ALL_CONFIDENT_SITES -baq CALCULATE_AS_NECESSARY";
$run{gatk_rtc} =
"java -Xmx7g -jar \$GATKPATH/GenomeAnalysisTK.jar -T RealignerTargetCreator";
$run{gatk_rtc_options} = "";
$run{gatk_ir} =
"java -Xmx14g -jar \$GATKPATH/GenomeAnalysisTK.jar -T IndelRealigner";
$run{gatk_ir_options} = "";
$run{gatk_combine} =
"java -Xmx14g -jar \$GATKPATH/GenomeAnalysisTK.jar -T CombineVariants";
$run{gatk_combine_options} = "-assumeIdenticalSamples";
$run{generate_consensus} =
"perl /exports/work/biology_ieb_bioinfx/scripts/generate_consensus_from_fasta_and_vcf.pl";
$run{samtools} = "samtools";

$run{qsub_options} = "#!/bin/bash\n";
$run{qsub_options} .= "#\$ -cwd -l h_rt=5:59:00\n";
$run{qsub_scripts} = ". /etc/profile\n";
$run{qsub_scripts} .=
  ". /exports/work/biology_ieb_mblaxter/software/.softwarerc\n";

$run{job_size}   = 100_000;
$run{unit_size}  = 20;
$run{sleep_secs} = 600;
$run{dry_run}    = 0;
$run{force}      = 1;
$run{max_lanes}  = 3;
$run{stage}      = "ALL";
$run{continue}   = 0;
$run{cleanup}    = 1;

my $max_species    = 20;
my $user_species   = "";
my $user_reference = "";
my $user_sample    = "";
my $user_lane      = "";

# Stages:
# 1s = Pass 1 Stampy
# 1m = Pass 1 Merge
# 1d = Pass 1 Remove Duplicates
# 1i = Pass 1 Indel Realign
# 1v = Pass 1 VCF
# sp = Species reference
# 2s = Pass 2 Stampy
# 2m = Pass 2 Merge
# 2d = Pass 2 Remove Duplicates
# 2i = Pass 2 Indel Realign
# 2v = Pass 2 VCF

my %stages = (
    "ALL" => 0,
    "1s"  => 1,
    "1d"  => 2,
    "1i"  => 3,
    "1v"  => 4,
#    "sp"  => 5,
#    "2s"  => 6,
#    "2d"  => 7,
#    "2i"  => 8,
#    "2v"  => 9,
);

my $options_okay = GetOptions(
    'input=s'          => \$in_filename,
    'species=s'        => \$user_species,
    'reference=s'      => \$user_reference,
    'sample=s'         => \$user_sample,
    'lane=s'           => \$user_lane,
    'output_dir=s'     => \$output_dir,
    'max_species=i'    => \$max_species,
    'max_lanes=i'      => \$run{max_lanes},
    'stampy_options=s' => \$run{stampy_user_options},
    'stage=s'          => \$run{stage},
    'jobsize=i'        => \$run{job_size},
    'unitsize=f'       => \$run{unit_size},
    'sleep=i'          => \$run{sleep_secs},
    'dry_run'          => \$run{dry_run},
    'continue'         => \$run{continue},
    'force!'           => \$run{force},
    'cleanup!'         => \$run{cleanup},
);

# Always continue if running everything
if ( $run{stage} eq "ALL" ) { $run{continue} = 1; }

sub usage {
    croak
"\nUsage: sge_align -i config_file -j job_size -u unit_size -max_species max_species -max_lanes max_lanes -stage stage_id -species species_name -r reference -u sample -sleep seconds -l lane -d dry_run\n";
}

usage() if !$options_okay;

croak "\nPlease specify a config file with -i\n" if ( $in_filename eq "" );

print STDERR localtime() . " | MAIN : sge_align\n";
print STDERR localtime() . " | MAIN : Config file      = $in_filename\n";
print STDERR localtime()
  . " | MAIN : Job size         = $run{job_size} reads\n";
print STDERR localtime() . " | MAIN : Unit size        = $run{unit_size} Mbp\n";
print STDERR localtime() . " | MAIN : Stage            = $run{stage}";
if ( $run{continue} ) {
    print STDERR " to end";
}
print STDERR "\n";
if ( $user_reference ne "" ) {
    print STDERR localtime() . " | MAIN : Reference        = $user_reference\n";
}
if ( $user_species ne "" ) {
    print STDERR localtime() . " | MAIN : Species          = $user_species\n";
}
if ( $user_sample ne "" ) {
    print STDERR localtime() . " | MAIN : Sample           = $user_sample\n";
}
if ( $user_lane ne "" ) {
    print STDERR localtime() . " | MAIN : Lane             = $user_lane\n";
}
print STDERR localtime() . " | MAIN : Output directory = $output_dir\n";
print STDERR localtime()
  . " | MAIN : Sleep time       = $run{sleep_secs} seconds\n";
if ( $run{force} ) {
    print STDERR localtime() . " | MAIN : Force recreation of output files\n";
}
if ( $run{dry_run} ) {
    print STDERR localtime() . " | MAIN : Dry run\n";
}

print STDERR localtime() . " | MAIN : Initialized\n";

print STDERR localtime() . " | MAIN : Loading config file...";

# Process config file
open my $in_file, '<', $in_filename
  or croak "Can't open $in_filename: $OS_ERROR!\n";

my %species;
my %reference;
while ( my $lane = <$in_file> ) {
    chomp $lane;
    my @fields = split /\t/, $lane;
    next if ( ( $user_species   ne "" ) && ( $fields[0] ne $user_species ) );
    next if ( ( $user_reference ne "" ) && ( $fields[1] ne $user_reference ) );
    next if ( ( $user_sample    ne "" ) && ( $fields[2] ne $user_sample ) );
    next if ( ( $user_lane      ne "" ) && ( $fields[3] ne $user_lane ) );

    $reference{ $fields[1] }++;

    $species{ $fields[0] }{alignments}{ $fields[1] }++;

    # Check input files exist
    croak "Input file $fields[4] not found!\n" unless ( -e $fields[4] );

    $species{ $fields[0] }{samples}{ $fields[2] }{ $fields[3] }{1} = $fields[4];

    if ( defined $fields[5] ) {
        croak "Input file $fields[5] not found!\n" unless ( -e $fields[5] );
        $species{ $fields[0] }{samples}{ $fields[2] }{ $fields[3] }{2} =
          $fields[5];
    }

    if ( defined $fields[6] ) {
        $species{ $fields[0] }{samples}{ $fields[2] }{ $fields[3] }{parts} =
          int( $fields[6] / $run{job_size} ) + 1;
    }
}

close $in_file;

# Check references
print STDERR "OK\n"
  . localtime()
  . " | MAIN : Check references and generate intervals...\n";

foreach my $ref_filename ( keys %reference ) {

    # Check reference exists
    my ( $refbasename, $refdir ) =
      fileparse( $ref_filename, ".fasta", ".fas", ".fa", ".fna" );
    croak "Reference, index and hash for $ref_filename not found!\n"
      unless ( ( -e $ref_filename )
        && ( -e "$refdir/$refbasename\.stidx" )
        && ( -e "$refdir/$refbasename\.sthash" ) );

    # Generate intervals
    try_system(
"partition_reference_into_intervals.pl -r $ref_filename -u $run{unit_size}",
        $run{dry_run}, undef, "| MAIN"
    );
    my @interval_files = glob "$refbasename\.intervals/*.intervals";
    $run{refjobs}{$refbasename} = scalar @interval_files;
}

print STDERR "OK\n"
  . localtime()
  . " | MAIN : Setting up directories and threads...";

# Set up Parallel::ForkManager and directories
my $align_procs = 0;

$output_dir = abs_path($output_dir);
if ( !$run{dry_run} ) { mkdir("$output_dir") unless -d "$output_dir"; }

foreach my $species ( sort keys %species ) {
    if ( !$run{dry_run} ) {
        mkdir("$output_dir/$species") unless -d "$output_dir/$species";
    }
    foreach my $reference ( sort keys %{ $species{$species}{alignments} } ) {
        $align_procs++;

        next if $run{dry_run};
        my ( $refbasename, $refdir ) =
          fileparse( $reference, ".fasta", ".fas", ".fa", ".fna" );
        mkdir("$output_dir/$species/$refbasename")
          unless -d "$output_dir/$species/$refbasename";

        foreach my $sample ( sort keys %{ $species{$species}{samples} } ) {
            mkdir("$output_dir/$species/$refbasename/$sample")
              unless -d "$output_dir/$species/$refbasename/$sample";
            foreach
              my $lane ( sort keys %{ $species{$species}{samples}{$sample} } )
            {
                mkdir("$output_dir/$species/$refbasename/$sample/$lane")
                  unless -d "$output_dir/$species/$refbasename/$sample/$lane";
            }
        }
    }
}

my $species_pm = new Parallel::ForkManager($align_procs);
$species_pm->set_max_procs($max_species);
my $success = 0;
$species_pm->run_on_finish(
    sub {
        my ( $pid, $exit_code ) = @_;
        $success += $exit_code;
    }
);

print STDERR "OK\n" . localtime() . " | MAIN : Generating VCFs...\n";

# Generate VCFs for all species
foreach my $species ( sort keys %species ) {
    foreach my $reference ( sort keys %{ $species{$species}{alignments} } ) {

        print STDERR localtime()
          . " | MAIN : Start thread for species $species, reference $reference\n";

        # Initialise
        $species_pm->start and next;

        my $error_id = " | $species | $reference";
        print STDERR localtime() . "$error_id : Thread started\n";
        my ( $refbasename, $refdir ) =
          fileparse( $reference, ".fasta", ".fas", ".fa", ".fna" );
        my $samples_ref = $species{$species}{samples};

        # Make Pass 1 VCF against original reference
        if (   ( $run{stage} =~ /^1/ )
            or ( ( $run{continue} ) and ( $stages{ $run{stage} } <= 5 ) ) )
        {
            print STDERR localtime() . "$error_id | 1  : Make VCF\n";

            # Alignment and SNP call to original reference
            make_vcf(
                {
                    samples_ref => $species{$species}{samples},
                    species     => $species,
                    reference   => $reference,
                    refdir      => $refdir,
                    output_dir  => $output_dir,
                    species_pm  => $species_pm,
                    pass        => 1,
                    run_ref     => \%run,
                    stages_ref  => \%stages,
                }
            );
        }

        foreach my $sample ( sort keys %{ $species{$species}{samples} } ) {
            foreach
              my $lane ( sort keys %{ $species{$species}{samples}{$sample} } )
            {
                $species_pm->finish(0)
                  if (
                    no_input(
"$output_dir/$species/$refbasename/$sample/$lane/$species.$sample.$lane.$refbasename.pass1.rmdup.realign.bam",
                        $run{dry_run},
                        $error_id
                    )
                  );
            }
        }

        print STDERR localtime() . "$error_id : Complete\n";
        $species_pm->finish(0);
        # Make species reference
        my $species_reference = "$species.$refbasename.speciesref.fasta";
        if (   ( $run{stage} eq "sp" )
            or ( ( $run{continue} ) and ( $stages{ $run{stage} } <= 5 ) ) )
        {
            print STDERR localtime() . "$error_id : Make species reference\n";
            $species_reference = make_species_reference(
                {
                    species     => $species,
                    refbasename => $refbasename,
                    reference   => $reference,
                    output_dir  => $output_dir,
                    run_ref     => \%run,
                    species_pm  => $species_pm,
                }
            );
        }

        $species_pm->finish(0)
          if (
            (
                no_input(
"$output_dir/$species/$refbasename/$species\.$refbasename\.speciesref.fasta",
                    $run{dry_run},
                    $error_id
                )
            )
            or (
                no_input(
"$output_dir/$species/$refbasename/$species\.$refbasename\.speciesref.stidx",
                    $run{dry_run},
                    $error_id
                )
            )
            or (
                no_input(
"$output_dir/$species/$refbasename/$species\.$refbasename\.speciesref.sthash",
                    $run{dry_run},
                    $error_id
                )
            )
            or (
                no_input(
"$output_dir/$species/$refbasename/$species\.$refbasename\.speciesref.dict",
                    $run{dry_run},
                    $error_id
                )
            )
            or (
                no_input(
"$output_dir/$species/$refbasename/$species\.$refbasename\.speciesref.fasta.fai",
                    $run{dry_run},
                    $error_id
                )
            )
          );

        # Make Pass 2 VCF against species reference
        if (   ( $run{stage} =~ /^2/ )
            or ( ( $run{continue} ) and ( $stages{ $run{stage} } <= 9 ) ) )
        {
            print STDERR localtime() . "$error_id | 2  : Make VCF\n";
            make_vcf(
                {
                    samples_ref => $species{$species}{samples},
                    species     => $species,
                    reference =>
                      "$output_dir/$species/$refbasename/$species_reference",
                    output_dir => $output_dir,
                    species_pm => $species_pm,
                    pass       => 2,
                    run_ref    => \%run,
                    stages_ref => \%stages,
                }
            );

        }
        print STDERR localtime() . "$error_id : Complete\n";
        $species_pm->finish(0);
    }
}
$species_pm->wait_all_children;

if ( $success == 0 ) {
    print STDERR localtime() . " | MAIN : Complete\n";
}
else {
    print STDERR localtime()
      . " | MAIN : Failed to create all VCF files; exiting\n";
}
exit;

# SUBROUTINES

sub make_species_reference {
    my ($arg_ref)   = @_;
    my $species     = $arg_ref->{species};
    my $refbasename = $arg_ref->{refbasename};
    my $reference   = $arg_ref->{reference};
    my $output_dir  = $arg_ref->{output_dir};
    my $run_ref     = $arg_ref->{run_ref};
    my $species_pm  = $arg_ref->{species_pm};

    # Create species reference, generate index, hash and copy to Eddie

    my $species_reference = "$species.$refbasename.speciesref.fasta";
    my $species_refbase   = "$species.$refbasename.speciesref";
    my $error_id          = "| $species | $refbasename | r";

    return $species_reference
      if check_recreation(
        "$output_dir/$species/$refbasename/$species_reference",
        $run_ref->{force}, $error_id );

    print STDERR localtime()
      . " $error_id : Generate consensus from reference and Pass 1 VCF $species.$refbasename.vcf\n";

    my $spref_commands =
"$run_ref->{generate_consensus} -r $reference -v $output_dir/$species/$refbasename/$species.$refbasename.pass1.vcf > $output_dir/$species/$refbasename/$species_reference\n";
    $spref_commands .=
"if [ -e \"$output_dir/$species/$refbasename/$species_refbase\.stidx\" ]\nthen\n";
    $spref_commands .=
      "\trm $output_dir/$species/$refbasename/$species_refbase\.stidx\nfi\n";
    $spref_commands .=
"$run_ref->{stampy} -G $output_dir/$species/$refbasename/$species_refbase $output_dir/$species/$refbasename/$species_reference\n";
    $spref_commands .=
"if [ -e \"$output_dir/$species/$refbasename/$species_refbase\.sthash\" ]\nthen\n";
    $spref_commands .=
      "\trm $output_dir/$species/$refbasename/$species_refbase\.sthash\nfi\n";
    $spref_commands .=
"$run_ref->{stampy} -g $output_dir/$species/$refbasename/$species_refbase -H $output_dir/$species/$refbasename/$species_refbase\n";
    $spref_commands .=
"$run_ref->{picard_dict} R=$output_dir/$species/$refbasename/$species_refbase\.fasta O=$output_dir/$species/$refbasename/$species_refbase\.dict\n";
    $spref_commands .=
"$run_ref->{samtools} faidx $output_dir/$species/$refbasename/$species_refbase\.fasta\n";

    # Open error log
    my $log_filename =
      "$output_dir/$species/$refbasename/$species.$refbasename\.r.log";
    print STDERR localtime() . " $error_id : Opening error log $log_filename\n";

    my $qsub_log;
    if ( !$run_ref->{dry_run} ) {
        open $qsub_log, ">", $log_filename
          or croak "Can't open $log_filename: $OS_ERROR\n";
    }

    do_qsub(
        {
            commands  => $spref_commands,
            run_ref   => $run_ref,
            directory => "$output_dir/$species/$refbasename",
            name      => "$species.$refbasename",
            pm        => $species_pm,
            qsub_log  => $qsub_log,
            error_id  => $error_id,
            slots     => 1,
            pass      => 1,
            stage     => "r",
        }
    );
    if ( !$run_ref->{dry_run} ) { close $qsub_log; }

    print STDERR localtime()
      . " $error_id : Consensus $output_dir/$species/$refbasename/$species_reference created\n";

    return $species_reference;
}

sub make_vcf {
    my ($arg_ref)   = @_;
    my $samples_ref = $arg_ref->{samples_ref};
    my $species     = $arg_ref->{species};
    my $reference   = $arg_ref->{reference};
    my $output_dir  = $arg_ref->{output_dir};
    my $species_pm  = $arg_ref->{species_pm};
    my $pass        = $arg_ref->{pass};
    my $run_ref     = $arg_ref->{run_ref};
    my $stages_ref  = $arg_ref->{stages_ref};

    my $speciesref = "";
    if ( $pass == 2 ) { $speciesref = ".speciesref"; }

    my ( $refbasename, $refdir ) = fileparse(
        $reference,        "$speciesref.fasta",
        "$speciesref.fas", "$speciesref.fa",
        "$speciesref.fna"
    );

    if ( ( $pass == 2 ) && ( $refbasename =~ /$species\.(.+)/ ) ) {
        $refbasename = $1;
    }

    my $error_id = "| $species | $refbasename";

    return
      if check_recreation(
"$output_dir/$species/$refbasename/$species\.$refbasename\.pass$pass\.vcf",
        $run_ref->{force}, $error_id
      );

    # Set up Parallel::ForkManager
    my $procs = 0;
    foreach my $sample ( sort keys %{$samples_ref} ) {
        foreach my $lane ( sort keys %{ $samples_ref->{$sample} } ) {
            $procs++;
        }
    }

    my $success = 0;
    my $lane_pm = new Parallel::ForkManager($procs);
    $lane_pm->set_max_procs( $run_ref->{max_lanes} );
    $lane_pm->run_on_finish(
        sub {
            my ( $pid, $exit_code ) = @_;
            $success += $exit_code;
        }
    );

    # Generate BAM files for all lanes
    my $stage_incr = ( $pass == 1 ) ? 0 : ( $pass == 2 ) ? 5 : 0;

    foreach my $sample ( sort keys %{$samples_ref} ) {
        foreach my $lane ( sort keys %{ $samples_ref->{$sample} } ) {

            $lane_pm->start and next;

            $error_id .= " | $sample | $lane";
            print STDERR localtime() . " $error_id | $pass  : Thread started\n";

            make_bam(
                {
                    lane_ref    => $samples_ref->{$sample}{$lane},
                    output_dir  => $output_dir,
                    species     => $species,
                    reference   => $reference,
                    refbasename => $refbasename,
                    refdir      => $refdir,
                    sample      => $sample,
                    lane        => $lane,
                    run_ref     => $run_ref,
                    stages_ref  => $stages_ref,
                    stage_incr  => $stage_incr,
                    lane_pm     => $lane_pm,
                    error_id    => "$error_id | $pass",
                    pass        => $pass,

                }
            );
            print STDERR localtime() . " $error_id | $pass  : Complete\n";
            $lane_pm->finish(0);
        }
    }

    # Run GATK on all input lanes combined
    $lane_pm->wait_all_children;

    if ( $success > 0 ) {
        $species_pm->finish(1);
        return;
    }

    if (
        ( $run_ref->{stage} =~ /v$/ )
        or (   ( $run_ref->{continue} )
            && ( $stages{ $run_ref->{stage} } <= ( 4 + $stage_incr ) ) )
      )
    {
        run_ug(
            {
                run_ref     => $run_ref,
                samples_ref => $samples_ref,
                species     => $species,
                reference   => $reference,
                refbasename => $refbasename,
                output_dir  => $output_dir,
                error_id    => "$error_id | $pass",
                pm          => $species_pm,
                pass        => $pass,
            }
        );
    }
    return;
}

sub make_bam {
    my ($arg_ref) = @_;
    my %lane;

    my $run_ref = $arg_ref->{run_ref};
    $lane{ref}         = $arg_ref->{lane_ref};
    $lane{output_dir}  = $arg_ref->{output_dir};
    $lane{species}     = $arg_ref->{species};
    $lane{reference}   = $arg_ref->{reference};
    $lane{refbasename} = $arg_ref->{refbasename};
    $lane{refdir}      = $arg_ref->{refdir};
    $lane{sample}      = $arg_ref->{sample};
    $lane{lane}        = $arg_ref->{lane};
    $lane{stages_ref}  = $arg_ref->{stages_ref};
    $lane{stage_incr}  = $arg_ref->{stage_incr};
    $lane{pm}          = $arg_ref->{lane_pm};
    $lane{error_id}    = $arg_ref->{error_id};
    $lane{pass}        = $arg_ref->{pass};

    print STDERR localtime() . " $lane{error_id}  : Start alignment\n";
    $lane{write_dir} =
"$lane{output_dir}/$lane{species}/$lane{refbasename}/$lane{sample}/$lane{lane}";
    $lane{bam_stub} =
"$lane{species}.$lane{sample}.$lane{lane}.$lane{refbasename}.pass$lane{pass}";

    return
      if check_recreation(
        "$lane{write_dir}/$lane{bam_stub}\.rmdup.realign.bam",
        $run_ref->{force}, $lane{error_id} );

    # Calculate number of Stampy jobs
    if ( !defined $lane{ref}->{parts} ) {
        $lane{ref}->{parts} =
          count_jobs( $lane{ref}->{1}, $run_ref->{job_size}, $lane{error_id} );
    }
    else {
        print STDERR localtime()
          . " $lane{error_id} : Read number provided, $lane{ref}->{parts} jobs to run\n";
    }

    # RUN STAMPY
    if ( stage_ok( "s", "1", $run_ref, \%lane ) ) {
        make_stampy(
            {
                run_ref  => $run_ref,
                lane_ref => \%lane,
                error_id => $lane{error_id} . "s",
            }
        );
    }

    # REMOVE DUPLICATES
    if ( stage_ok( "d", "2", $run_ref, \%lane ) ) {
        remove_duplicates(
            {
                run_ref  => $run_ref,
                lane_ref => \%lane,
                error_id => $lane{error_id} . "d",
            }
        );
    }

    # REALIGN INDELS
    if ( stage_ok( "i", "3", $run_ref, \%lane ) ) {
        realign_indels(
            {
                run_ref  => $run_ref,
                lane_ref => \%lane,
                error_id => $lane{error_id} . "i",
            }
        );

    }

    # CLEAN UP
    if (   ( $run_ref->{cleanup} )
        && ( -e "$lane{write_dir}/$lane{bam_stub}\.rmdup.realign.bam" ) )
    {
        clean_up(
            {
                run_ref  => $run_ref,
                lane_ref => \%lane,
                error_id => $lane{error_id},
            }
        );
    }

    return;

}

sub run_ug {
    my ($arg_ref)   = @_;
    my $run_ref     = $arg_ref->{run_ref};
    my $samples_ref = $arg_ref->{samples_ref};
    my $species     = $arg_ref->{species};
    my $reference   = $arg_ref->{reference};
    my $refbasename = $arg_ref->{refbasename};
    my $output_dir  = $arg_ref->{output_dir};
    my $error_id    = $arg_ref->{error_id};
    my $pm          = $arg_ref->{pm};
    my $pass        = $arg_ref->{pass};

    my $write_dir = "$output_dir/$species/$refbasename";
    my $vcf_stub  = "$species\.$refbasename\.pass$pass";

    my $gatk_input = "";
    foreach my $sample ( sort keys %{$samples_ref} ) {
        foreach my $lane ( sort keys %{ $samples_ref->{$sample} } ) {
            return
              if no_input(
"$write_dir/$sample/$lane/$species\.$sample\.$lane\.$refbasename\.pass$pass\.rmdup.realign.bam",
                $run_ref->{dry_run}, $error_id
              );
            $gatk_input .=
"-I $write_dir/$sample/$lane/$species\.$sample\.$lane\.$refbasename\.pass$pass\.rmdup.realign.bam ";
        }
    }
    print STDERR localtime() . " $error_id | $pass  : Call genotypes\n";

    # Open error log
    my $log_filename = "$write_dir/$vcf_stub\.v.log";
    print STDERR localtime() . " $error_id : Opening error log $log_filename\n";

    my $qsub_log;
    if ( !$run_ref->{dry_run} ) {
        open $qsub_log, ">", $log_filename
          or croak "Can't open $log_filename: $OS_ERROR\n";
    }
    my $ug_command =
"$run_ref->{gatk_ug} $run_ref->{gatk_ug_options} -L $refbasename\.intervals/\$SGE_TASK_ID.intervals -R $reference -o $write_dir/$vcf_stub\.\$SGE_TASK_ID\.vcf $gatk_input\n";

    do_qsub_array(
        {
            commands  => $ug_command,
            run_ref   => $run_ref,
            jobs      => $run_ref->{refjobs}{$refbasename},
            directory => $write_dir,
            name      => $vcf_stub,
            pm        => $pm,
            qsub_log  => $qsub_log,
            error_id  => $error_id,
            slots     => 8,
            pass      => $pass,
            stage     => "v",
        }
    );

    if ( !$run_ref->{dry_run} ) { close $qsub_log; }

    merge_files(
        {
            jobs      => $run_ref->{refjobs}{$refbasename},
            directory => $write_dir,
            file_stub => "$vcf_stub",
            in_ext    => "vcf",
            out_ext   => "vcf",
            reference => $reference,
            pass      => $pass,
            stage     => "vm",
            pm        => $pm,
            run_ref   => $run_ref,
            error_id  => $error_id . "vm",
        }
    );

    print STDERR localtime() . " $error_id  : Genotyping complete\n";

    return;
}

sub stage_ok {
    my ( $type, $offset, $run_ref, $lane_ref ) = @_;

    return (
        ( $run_ref->{stage} =~ /$type$/ )
          or (
            ( $run_ref->{continue} )
            && ( $lane_ref->{stages_ref}->{ $run_ref->{stage} } <=
                ( $offset + $lane_ref->{stage_incr} ) )
          )
    );
}

sub count_jobs {
    my ( $reads_file, $job_size, $error_id ) = @_;
    print STDERR localtime() . " $error_id  : Calculating number of jobs...\n";
    my ( $basename, $dir ) = fileparse($reads_file);
    my $reads = 0;
    if ( $basename =~ /gz/ ) {
        $reads = `gunzip -c $reads_file | wc -l` / 4;
    }
    else {
        $reads = `wc -l < $reads_file` / 4;
    }
    my $jobs = int( $reads / $job_size ) + 1;
    print STDERR localtime() . " $error_id  : $jobs jobs to run\n";
    return $jobs;
}

sub make_stampy {
    my ($arg_ref) = @_;
    my $run_ref   = $arg_ref->{run_ref};
    my $lane_ref  = $arg_ref->{lane_ref};
    my $error_id  = $arg_ref->{error_id};

    print STDERR localtime() . " $error_id : Run alignment\n";

    # Check if output file already exists
    return
      if check_recreation( "$lane_ref->{write_dir}/$lane_ref->{bam_stub}\.bam",
        $run_ref->{force}, $lane_ref->{error_id} );

    # If final BAM doesn't exist, but part BAMs do, skip to merge
    if (
        !(
            check_recreation_glob(
                $lane_ref->{write_dir}, $lane_ref->{bam_stub},
                "bam",                  $lane_ref->{ref}->{parts},
                $run_ref->{force},      $error_id
            )
        )
      )
    {

        # Open error log
        my $log_filename =
          "$lane_ref->{write_dir}/$lane_ref->{bam_stub}\.s.log";
        print STDERR localtime()
          . " $error_id : Opening error log $log_filename\n";

        my $qsub_log;
        if ( !$run_ref->{dry_run} ) {
            open $qsub_log, ">", $log_filename
              or croak "Can't open $log_filename: $OS_ERROR\n";
        }

        # Unzip to tmp dir if gzipped

        my $in1stub = "";
        my $in2stub = "";
        my ( $in1basename, $in1dir ) = fileparse( $lane_ref->{ref}->{1} );
        my $in2basename = "";
        my $in2dir      = "";
        if ( defined $lane_ref->{ref}->{2} ) {
            ( $in2basename, $in2dir ) = fileparse( $lane_ref->{ref}->{2} );
        }

        if ( ( $in1basename =~ /(.+).gz/ ) or ( $in2basename =~ /(.+).gz/ ) ) {
            my $gzip_commands = "";
            if ( $in1basename =~ /(.+).gz$/ ) {
                $in1stub = $1;
                $gzip_commands .=
"gunzip -c $lane_ref->{ref}->{1} > /exports/work/scratch/$in1stub\n";
            }
            if ( $in2basename =~ /(.+).gz$/ ) {
                $in2stub = $1;
                $gzip_commands .=
"gunzip -c $lane_ref->{ref}->{2} > /exports/work/scratch/$in2stub\n";
            }
            do_qsub(
                {
                    commands  => $gzip_commands,
                    run_ref   => $run_ref,
                    directory => $lane_ref->{write_dir},
                    name      => $lane_ref->{bam_stub},
                    pm        => $lane_ref->{pm},
                    qsub_log  => $qsub_log,
                    error_id  => $error_id,
                    slots     => 1,
                    pass      => $lane_ref->{pass},
                    stage     => "sg",
                    sleep     => 60,
                }
            );
        }

        my $qsub_commands =
"if [ -e \"$lane_ref->{write_dir}/$lane_ref->{bam_stub}.\$SGE_TASK_ID.bam\" ]\nthen\n";

        $qsub_commands .=
"\trm $lane_ref->{write_dir}/$lane_ref->{bam_stub}.\$SGE_TASK_ID.bam\nfi\n";

        $qsub_commands =
"if [ -e \"$lane_ref->{write_dir}/$lane_ref->{bam_stub}.\$SGE_TASK_ID.sam\" ]\nthen\n";

        $qsub_commands .=
"\trm $lane_ref->{write_dir}/$lane_ref->{bam_stub}.\$SGE_TASK_ID.sam\nfi\n";

        $qsub_commands .=
"$run_ref->{stampy} $run_ref->{stampy_options} $run_ref->{stampy_user_options} ";
        $qsub_commands .=
          "--processpart \$SGE_TASK_ID/$lane_ref->{ref}->{parts} ";
        $qsub_commands .=
"--readgroup=ID:$lane_ref->{species}.$lane_ref->{sample}.$lane_ref->{lane},SM:$lane_ref->{species}.$lane_ref->{sample},PL:Illumina ";

        my $stampy_stub = "$lane_ref->{refdir}/$lane_ref->{refbasename}";
        if ( $lane_ref->{pass} == 2 ) {
            $stampy_stub =
"$lane_ref->{refdir}/$lane_ref->{species}\.$lane_ref->{refbasename}.speciesref";
        }
        $qsub_commands .= "-g $stampy_stub ";
        $qsub_commands .= "-h $stampy_stub ";
        $qsub_commands .=
"-o $lane_ref->{write_dir}/$lane_ref->{bam_stub}.\$SGE_TASK_ID.sam -M";

        if ( $in1basename =~ /.gz/ ) {
            $qsub_commands .= " /exports/work/scratch/$in1stub";
        }
        else {
            $qsub_commands .= " $in1basename";
        }

        if ( defined $lane_ref->{ref}->{2} ) {
            if ( $in2basename =~ /.gz/ ) {
                $qsub_commands .= " /exports/work/scratch/$in2stub";
            }
            else {
                $qsub_commands .= " $in2basename";
            }
        }
        $qsub_commands .= "\n";

        $qsub_commands .=
"$run_ref->{picard_sort} $run_ref->{picard_sort_options} INPUT=$lane_ref->{write_dir}/$lane_ref->{bam_stub}.\$SGE_TASK_ID.sam OUTPUT=$lane_ref->{write_dir}/$lane_ref->{bam_stub}.\$SGE_TASK_ID.bam\n";

        $qsub_commands .=
          "rm $lane_ref->{write_dir}/$lane_ref->{bam_stub}.\$SGE_TASK_ID.sam\n";

        do_qsub_array(
            {
                commands  => $qsub_commands,
                run_ref   => $run_ref,
                jobs      => $lane_ref->{ref}->{parts},
                directory => $lane_ref->{write_dir},
                name      => $lane_ref->{bam_stub},
                pm        => $lane_ref->{pm},
                qsub_log  => $qsub_log,
                error_id  => $error_id,
                slots     => 1,
                pass      => $lane_ref->{pass},
                stage     => "s",
            }
        );
        
        
        if ( !$run_ref->{dry_run} ) {
            close $qsub_log;

            if ( $in1basename =~ /.gz/ ) {
                try_rm( "/exports/work/scratch/$in1stub",
                    $run_ref->{dry_run}, $lane_ref->{pm}, $error_id );
            }

            if ( defined $lane_ref->{ref}->{2} ) {
                if ( $in2basename =~ /.gz/ ) {
                    try_rm( "/exports/work/scratch/$in2stub",
                        $run_ref->{dry_run}, $lane_ref->{pm}, $error_id );
                }
            }

        }

    }

    merge_files(
        {
            jobs      => $lane_ref->{ref}->{parts},
            directory => $lane_ref->{write_dir},
            file_stub => $lane_ref->{bam_stub},
            in_ext    => "bam",
            out_ext   => "bam",
            reference => $lane_ref->{reference},
            pass      => $lane_ref->{pass},
            stage     => "sm",
            pm        => $lane_ref->{pm},
            run_ref   => $run_ref,
            error_id  => $lane_ref->{error_id} . "sm",
        }
    );

    print STDERR localtime() . " $error_id : Complete\n";

    return;
}

sub remove_duplicates {
    my ($arg_ref) = @_;
    my $run_ref   = $arg_ref->{run_ref};
    my $lane_ref  = $arg_ref->{lane_ref};
    my $error_id  = $arg_ref->{error_id};

    print STDERR localtime()
      . " $error_id : Removing duplicates from $lane_ref->{write_dir}/$lane_ref->{bam_stub}\.bam\n";

    return
      if no_input( "$lane_ref->{write_dir}/$lane_ref->{bam_stub}\.bam",
        $run_ref->{dry_run}, $error_id );

    return
      if check_recreation(
        "$lane_ref->{write_dir}/$lane_ref->{bam_stub}\.rmdup.bam",
        $run_ref->{force}, $error_id );

    my $rmdup_commands =
"$run_ref->{picard_rmdup} $run_ref->{picard_rmdup_options} INPUT=$lane_ref->{write_dir}/$lane_ref->{bam_stub}\.bam OUTPUT=$lane_ref->{write_dir}/$lane_ref->{bam_stub}\.rmdup.bam METRICS_FILE=$lane_ref->{write_dir}/$lane_ref->{bam_stub}\.rmdup.metrics\n";

    $rmdup_commands .=
"$run_ref->{samtools} index $lane_ref->{write_dir}/$lane_ref->{bam_stub}\.rmdup.bam\n";

    # Open error log
    my $log_filename = "$lane_ref->{write_dir}/$lane_ref->{bam_stub}.d.log";
    print STDERR localtime() . " $error_id : Opening error log $log_filename\n";

    my $qsub_log;
    if ( !$run_ref->{dry_run} ) {
        open $qsub_log, ">", $log_filename
          or croak "Can't open $log_filename: $OS_ERROR\n";
    }

    do_qsub(
        {
            commands  => $rmdup_commands,
            run_ref   => $run_ref,
            directory => $lane_ref->{write_dir},
            name      => $lane_ref->{bam_stub},
            pm        => $lane_ref->{pm},
            qsub_log  => $qsub_log,
            error_id  => $error_id,
            slots     => 2,
            pass      => $lane_ref->{pass},
            stage     => "d",
        }
    );
    if ( !$run_ref->{dry_run} ) { close $qsub_log; }

    print STDERR localtime()
      . " $error_id : Removing duplicates complete; indexed $lane_ref->{write_dir}/$lane_ref->{bam_stub}\.rmdup.bam\n";

    return;
}

sub realign_indels {
    my ($arg_ref) = @_;
    my $run_ref   = $arg_ref->{run_ref};
    my $lane_ref  = $arg_ref->{lane_ref};
    my $error_id  = $arg_ref->{error_id};

    print STDERR localtime()
      . " $error_id : Starting indel realignment for $lane_ref->{write_dir}/$lane_ref->{bam_stub}\.rmdup.bam\n";

    return
      if no_input( "$lane_ref->{write_dir}/$lane_ref->{bam_stub}\.rmdup.bam",
        $run_ref->{dry_run}, $error_id );

    return
      if check_recreation(
        "$lane_ref->{write_dir}/$lane_ref->{bam_stub}\.rmdup.realign.bam",
        $run_ref->{force}, $error_id );

    # Open error log
    my $log_filename = "$lane_ref->{write_dir}/$lane_ref->{bam_stub}.i.log";
    print STDERR localtime() . " $error_id : Opening error log $log_filename\n";

    my $qsub_log;
    if ( !$run_ref->{dry_run} ) {
        open $qsub_log, ">", $log_filename
          or croak "Can't open $log_filename: $OS_ERROR\n";
    }

    if (!(check_recreation("$lane_ref->{write_dir}/$lane_ref->{bam_stub}\.intervals", $run_ref->{force}, $error_id))) {
        my $rtc_command =
"$run_ref->{gatk_rtc} $run_ref->{gatk_rtc_options} -R $lane_ref->{reference} -I $lane_ref->{write_dir}/$lane_ref->{bam_stub}\.rmdup.bam -o $lane_ref->{write_dir}/$lane_ref->{bam_stub}\.intervals\n";

        do_qsub(
            {
                commands  => $rtc_command,
                run_ref   => $run_ref,
                directory => $lane_ref->{write_dir},
                name      => $lane_ref->{bam_stub},
                pm        => $lane_ref->{pm},
                qsub_log  => $qsub_log,
                error_id  => $error_id,
                slots     => 4,
                pass      => $lane_ref->{pass},
                stage     => "i",
            }
        );
    }

    return
      if no_input( "$lane_ref->{write_dir}/$lane_ref->{bam_stub}\.intervals",
        $run_ref->{dry_run}, $error_id );

    my $ir_command =
"$run_ref->{gatk_ir} $run_ref->{gatk_ir_options} -R $lane_ref->{reference}  -targetIntervals $lane_ref->{write_dir}/$lane_ref->{bam_stub}\.intervals -I $lane_ref->{write_dir}/$lane_ref->{bam_stub}\.rmdup.bam --out $lane_ref->{write_dir}/$lane_ref->{bam_stub}\.rmdup.realign.bam\n";

    do_qsub(
        {
            commands  => $ir_command,
            run_ref   => $run_ref,
            directory => $lane_ref->{write_dir},
            name      => $lane_ref->{bam_stub},
            pm        => $lane_ref->{pm},
            qsub_log  => $qsub_log,
            error_id  => $error_id,
            slots     => 8,
            pass      => $lane_ref->{pass},
            stage     => "i",
        }
    );

# Parallel version, for when GATK fixes -L bug
#    my $ir_command =
#"$run_ref->{gatk_ir} $run_ref->{gatk_ir_options} -R $lane_ref->{refdir}/$lane_ref->{reference} -L $lane_ref->{refbasename}\.intervals/\$SGE_TASK_ID.intervals -targetIntervals $lane_ref->{write_dir}/$lane_ref->{bam_stub}\.intervals -I $lane_ref->{write_dir}/$lane_ref->{bam_stub}\.rmdup.bam --out $lane_ref->{write_dir}/$lane_ref->{bam_stub}\.rmdup.realign.\$SGE_TASK_ID\.bam\n";
#    do_qsub_array(
#        {
#            commands  => $ir_command,
#            run_ref   => $run_ref,
#            jobs      => $run_ref->{refjobs}{ $lane_ref->{refbasename} },
#            directory => $lane_ref->{write_dir},
#            name      => $lane_ref->{bam_stub},
#            pm        => $lane_ref->{pm},
#            qsub_log  => $qsub_log,
#            error_id  => $error_id,
#            slots     => 1,
#            pass      => $lane_ref->{pass},
#            stage     => "i",
#        }
#    );
#    merge_files(
#        {
#            jobs      => $run_ref->{refjobs}{ $lane_ref->{refbasename} },
#            directory => $lane_ref->{write_dir},
#            file_stub => "$lane_ref->{bam_stub}.rmdup.realign",
#            in_ext    => "bam",
#            out_ext   => "bam",
#            reference => "$lane_ref->{refdir}/$lane_ref->{reference}",
#            pass      => $lane_ref->{pass},
#            stage     => "im",
#            pm        => $lane_ref->{pm},
#            run_ref   => $run_ref,
#            error_id  => $lane_ref->{error_id} . "im",
#        }
#    );

    if ( !$run_ref->{dry_run} ) { close $qsub_log; }

    print STDERR localtime()
      . " $error_id : Indel realignment complete for $lane_ref->{write_dir}/$lane_ref->{bam_stub}\.rmdup.bam\n";

    return;
}

sub clean_up {
    my ($arg_ref) = @_;
    my $run_ref   = $arg_ref->{run_ref};
    my $lane_ref  = $arg_ref->{lane_ref};
    my $error_id  = $arg_ref->{error_id};

    if ( -e "$lane_ref->{write_dir}/$lane_ref->{bam_stub}\.rmdup.bam" ) {
        try_rm( "$lane_ref->{write_dir}/$lane_ref->{bam_stub}\.bam",
            $run_ref->{dry_run}, $lane_ref->{pm}, $error_id );
        try_rm( "$lane_ref->{write_dir}/$lane_ref->{bam_stub}\.bai",
            $run_ref->{dry_run}, $lane_ref->{pm}, $error_id );
    }

    if ( -e "$lane_ref->{write_dir}/$lane_ref->{bam_stub}\.rmdup.realign.bam" )
    {
        try_rm( "$lane_ref->{write_dir}/$lane_ref->{bam_stub}\.rmdup.bam",
            $run_ref->{dry_run}, $lane_ref->{pm}, $error_id );
        try_rm( "$lane_ref->{write_dir}/$lane_ref->{bam_stub}\.rmdup.bam.bai",
            $run_ref->{dry_run}, $lane_ref->{pm}, $error_id );

    }

    return;
}

sub merge_files {
    my ($arg_ref) = @_;
    my $jobs      = $arg_ref->{jobs};
    my $directory = $arg_ref->{directory};
    my $file_stub = $arg_ref->{file_stub};
    my $in_ext    = $arg_ref->{in_ext};
    my $out_ext   = $arg_ref->{out_ext};
    my $reference = $arg_ref->{reference};
    my $pass      = $arg_ref->{pass};
    my $stage     = $arg_ref->{stage};
    my $pm        = $arg_ref->{pm};
    my $run_ref   = $arg_ref->{run_ref};
    my $error_id  = $arg_ref->{error_id};

    print STDERR localtime()
      . " $error_id : Merging output ($directory/$file_stub\.*\.$in_ext)\n";

    return
      if no_input( "$directory/$file_stub\.*\.$in_ext",
        $run_ref->{dry_run}, $error_id );

    return
      if check_recreation( "$directory/$file_stub\.$out_ext",
        $run_ref->{force}, $error_id );

    my $merge_command = "";
    if ( ( $out_ext eq "sam" ) or ( $out_ext eq "bam" ) ) {

        # Generate Picard MergeSamFiles command for this lane
        $merge_command =
"$run_ref->{picard_merge} $run_ref->{picard_merge_options} OUTPUT=$directory/$file_stub\.$out_ext";
        for ( my $i = 1 ; $i <= $jobs ; $i++ ) {
            $merge_command .= " INPUT=$directory/$file_stub\.$i\.$in_ext";
        }
        $merge_command .= "\n";
    }
    elsif ( $out_ext eq "vcf" ) {
        $merge_command =
"$run_ref->{gatk_combine} $run_ref->{gatk_combine_options} -R $reference -o $directory/$file_stub\.$out_ext";
        for ( my $i = 1 ; $i <= $jobs ; $i++ ) {
            $merge_command .= " -B:$i,VCF $directory/$file_stub\.$i\.$in_ext";
        }
        $merge_command .= "\n";
    }
    else {
        print STDERR localtime()
          . "$error_id : File extension $out_ext not recognised\n";
        return;
    }

    # Open error log
    my $log_filename = "$directory/$file_stub\.$stage\.log";
    print STDERR localtime() . " $error_id : Opening error log $log_filename\n";

    my $qsub_log;
    if ( !$run_ref->{dry_run} ) {
        open $qsub_log, ">", $log_filename
          or croak "Can't open $log_filename: $OS_ERROR\n";
    }

    do_qsub(
        {
            commands  => $merge_command,
            run_ref   => $run_ref,
            directory => $directory,
            name      => $file_stub,
            pm        => $pm,
            qsub_log  => $qsub_log,
            error_id  => $error_id,
            slots     => 8,
            pass      => $pass,
            stage     => $stage,
        }
    );
    if ( !$run_ref->{dry_run} ) { close $qsub_log; }

    # Clean up
    if (   ( $run_ref->{cleanup} )
        && ( -e "$directory/$file_stub\.$out_ext" ) )
    {
        for ( my $i = 1 ; $i <= $jobs ; $i++ ) {
            try_rm( "$directory/$file_stub\.$i\.$in_ext",
                $run_ref->{dry_run}, $pm, $error_id, 1 );
            if ( $out_ext eq "vcf" ) {
                try_rm( "$directory/$file_stub\.$i\.$in_ext\.idx",
                    $run_ref->{dry_run}, $pm, $error_id, 1 );
            }
        }
        print STDERR localtime()
          . " $error_id : Removed $directory/$file_stub\.*\.$in_ext\n";
    }
    print STDERR localtime() . " $error_id : Merging done\n";

    return;
}

sub check_recreation {
    my ( $output_file, $force, $error_id ) = @_;
    if ( -e $output_file ) {
        if ($force) {
            print STDERR localtime()
              . " $error_id : $output_file exists, but forcing recreation...\n";
            return 0;
        }
        else {
            print STDERR localtime()
              . " $error_id : $output_file exists, skipping creation\n";
            return 1;
        }
    }
    else {
        print STDERR localtime()
          . " $error_id : $output_file doesn't exist, creating...\n";
        return 0;
    }
    return 0;
}

sub check_recreation_glob {
    my ( $write_dir, $file_stub, $file_ext, $jobs, $force, $error_id ) = @_;

    my $existing = 0;
    for my $i ( 1 .. $jobs ) {
        if ( -e "$write_dir/$file_stub\.$i\.$file_ext" ) { $existing++; }
    }

    if ($existing) {
        if ($force) {
            print STDERR localtime()
              . " $error_id : $existing $write_dir/$file_stub\.\*\.$file_ext files exist, but forcing recreation...\n";
            return 0;
        }
        else {
            print STDERR localtime()
              . " $error_id : $existing $write_dir/$file_stub\.\*\.$file_ext  files exist, skipping creation\n";
            return 1;
        }
    }
    else {
        print STDERR localtime()
          . " $error_id : $existing $write_dir/$file_stub\.\*\.$file_ext files exist, creating...\n";
    }
    return 0;
}

sub do_qsub {
    my ($arg_ref) = @_;
    my $commands  = $arg_ref->{commands};
    my $run_ref   = $arg_ref->{run_ref};
    my $directory = $arg_ref->{directory};
    my $name      = $arg_ref->{name};
    my $qsub_log  = $arg_ref->{qsub_log};
    my $pm        = $arg_ref->{pm};
    my $error_id  = $arg_ref->{error_id};
    my $slots     = $arg_ref->{slots};
    my $pass      = $arg_ref->{pass};
    my $stage     = $arg_ref->{stage};
    my $sleep     = $arg_ref->{sleep};

    print STDERR localtime() . " $error_id : $commands\n";

    return if ( $run_ref->{dry_run} );

    if ( !defined($slots) ) {
        $slots = 1;
    }

    open my $qsub_sh, ">", "$directory/$name\.$stage\.qsub.sh"
      or croak
"Can't create qsub shell script $directory/$name\.$stage\.qsub.sh:$OS_ERROR\n";

    print $qsub_sh $run_ref->{qsub_options};
    print $qsub_sh "\#\$ -o $directory\n";
    print $qsub_sh "\#\$ -e $directory\n";
    print $qsub_sh "\#\$ -pe memory-2G $slots -R y\n";
    print $qsub_sh $run_ref->{qsub_scripts};
    print $qsub_sh $commands;
    print $qsub_sh "touch $directory/$name\.$stage\.done\n";
    close $qsub_sh;

    try_rm( "$directory/$name\.$stage\.done",
        $run_ref->{dry_run}, $pm, $error_id );
    my $qsub_output = `qsub $directory/$name\.$stage\.qsub.sh`;
    my $job_name    = "";
    if ( $qsub_output =~ /Your job (\d+) (.+) has been submitted/ ) {
        $job_name = $1;
        print STDERR localtime() . " $error_id : Launched qsub job $job_name\n";
    }
    else {
        croak "Can't launch qsub job $directory/$name\.$stage\.qsub.sh!\n";
    }

    while ( !( -e "$directory/$name\.$stage\.done" ) ) {
        my $is_job_running = `qstat -j $job_name 2>&1`;
        if ( $is_job_running =~ /do not exist/ ) {
            print STDERR localtime()
              . " $error_id : Qsub job $job_name failed before completion, exiting\n";
            if ( defined $pm ) { $pm->finish(1); }
            return;
        }
        if ( defined $sleep ) {
            sleep $sleep;
        }
        else {
            sleep $run_ref->{sleep_secs};
        }
    }

    print STDERR localtime()
      . " $error_id : Waiting for qsub job $job_name to complete...\n";

    my $qstat_cleared = 0;
    while ( $qstat_cleared == 0 ) {
        my $is_job_running = `qstat -j $job_name 2>&1`;
        if ( $is_job_running =~ /do not exist/ ) {
            $qstat_cleared = 1;
        }
        else {
            sleep 15;
        }
    }

    print STDERR localtime() . " $error_id : Qsub job $job_name complete\n";

    cat_log( $qsub_log, "$directory/$name\.$stage\.qsub.sh.e$job_name",
        $job_name );
    cat_log( $qsub_log, "$directory/$name\.$stage\.qsub.sh.o$job_name",
        $job_name );
    try_rm( "$directory/$name\.$stage\.qsub.sh\.e$job_name",
        $run_ref->{dry_run}, $pm, $error_id, 1 );
    try_rm( "$directory/$name\.$stage\.qsub.sh\.o$job_name",
        $run_ref->{dry_run}, $pm, $error_id, 1 );
    try_rm( "$directory/$name\.$stage\.qsub.sh\.po$job_name",
        $run_ref->{dry_run}, $pm, $error_id, 1 );
    try_rm( "$directory/$name\.$stage\.qsub.sh\.pe$job_name",
        $run_ref->{dry_run}, $pm, $error_id, 1 );
    try_rm( "$directory/$name\.$stage\.done",
        $run_ref->{dry_run}, $pm, $error_id );

    print STDERR localtime() . " $error_id : Qsub job $job_name cleaned up\n";

    # Wait for a minute so accounting works
    sleep 60;

    try_system( "get_qsub_cost.pl -j $job_name",
        $run_ref->{dry_run}, $pm, $error_id );

    return;
}

sub do_qsub_array {
    my ($arg_ref) = @_;
    my $commands  = $arg_ref->{commands};
    my $run_ref   = $arg_ref->{run_ref};
    my $jobs      = $arg_ref->{jobs};
    my $directory = $arg_ref->{directory};
    my $name      = $arg_ref->{name};
    my $qsub_log  = $arg_ref->{qsub_log};
    my $pm        = $arg_ref->{pm};
    my $error_id  = $arg_ref->{error_id};
    my $slots     = $arg_ref->{slots};
    my $pass      = $arg_ref->{pass};
    my $stage     = $arg_ref->{stage};

    print STDERR localtime() . " $error_id : $commands\n";

    return if ( $run_ref->{dry_run} );

    if ( !defined($slots) ) {
        $slots = 1;
    }

    open my $qsub_sh, ">", "$directory/$name\.$stage\.qsub.sh"
      or croak
"Can't create qsub shell script $directory/$name\.$stage\.qsub.sh:$OS_ERROR\n";

    print $qsub_sh $run_ref->{qsub_options};
    print $qsub_sh "\#\$ -o $directory\n";
    print $qsub_sh "\#\$ -e $directory\n";
    print $qsub_sh "\#\$ -pe memory-2G $slots -R y\n";
    print $qsub_sh $run_ref->{qsub_scripts};
    print $qsub_sh $commands;
    print $qsub_sh "touch $directory/$name\.$stage\.\$SGE_TASK_ID.done\n";
    close $qsub_sh;

    for my $i ( 1 .. $jobs ) {
        try_rm( "$directory/$name.\$stage\.$i.done",
            $run_ref->{dry_run}, $pm, $error_id );
    }
    my $qsub_output = `qsub -t 1-$jobs $directory/$name\.$stage\.qsub.sh`;
    my $array_name  = "";
    if ( $qsub_output =~ /Your job-array (\d+)\.(.+)has been submitted/ ) {
        $array_name = $1;
        print STDERR localtime()
          . " $error_id : Launched qsub job array $array_name\n";
    }
    else {
        croak "Can't launch qsub job $directory/$name.$stage\.qsub.sh!\n";
    }

    my $jobs_done = 0;
    while ( $jobs_done < $jobs ) {
        my $is_job_running = `qstat -j $array_name 2>&1`;
        if ( $is_job_running =~ /do not exist/ ) {
            print STDERR localtime()
              . " $error_id : Qsub job array $array_name failed before completion, exiting\n";
            if ( defined $pm ) { $pm->finish(1); }
            return;
        }
        sleep $run_ref->{sleep_secs};
        if ( -e "$directory/$name\.$stage\.1.done" ) {
            $jobs_done = `ls -l $directory/$name.$stage\.*.done | wc -l`;
        }
    }

    print STDERR localtime()
      . " $error_id : Waiting for qsub job array $array_name to complete...\n";

    my $qstat_cleared = 0;
    while ( $qstat_cleared == 0 ) {
        my $is_job_running = `qstat -j $array_name 2>&1`;
        if ( $is_job_running =~ /do not exist/ ) {
            $qstat_cleared = 1;
        }
        else {
            sleep 15;
        }
    }

    print STDERR localtime()
      . " $error_id : Qsub job array $array_name complete\n";

    cat_log_array( $qsub_log, $directory, "e", $array_name, $jobs );
    cat_log_array( $qsub_log, $directory, "o", $array_name, $jobs );
    try_rm( "$directory/*\.e$array_name\.*",
        $run_ref->{dry_run}, $pm, $error_id, 1 );
    try_rm( "$directory/*\.o$array_name\.*",
        $run_ref->{dry_run}, $pm, $error_id, 1 );
    try_rm( "$directory/*\.po$array_name\.*",
        $run_ref->{dry_run}, $pm, $error_id, 1 );
    try_rm( "$directory/*\.pe$array_name\.*",
        $run_ref->{dry_run}, $pm, $error_id, 1 );
    try_rm( "$directory/$name\.$stage\.*.done",
        $run_ref->{dry_run}, $pm, $error_id );

    print STDERR localtime()
      . " $error_id : Qsub job array $array_name cleaned up\n";

    # Wait for a minute so accounting works
    sleep 60;

    try_system( "get_qsub_cost.pl -j $array_name",
        $run_ref->{dry_run}, $pm, $error_id );

    return;
}

sub cat_log_array {
    my ( $log_file, $directory, $type, $array_name, $jobs ) = @_;

    for my $i ( 1 .. $jobs ) {
        my @job_filename = glob("$directory/*\.$type$array_name\.$i");
        cat_log( $log_file, $job_filename[0], "$array_name\.$i" );
    }

    return;
}

sub cat_log {
    my ( $log_file, $job_filename, $job_name ) = @_;

    if ( defined $log_file ) {
        print {$log_file} localtime() . "\t$job_name\t$job_filename\n";
        open my $job_file, "<", $job_filename
          or croak "Can't open $job_filename: $OS_ERROR\n";
        while (<$job_file>) { print {$log_file} $_; }
        close $job_file;
    }
    return;
}

# Returns true if no input files exist
sub no_input {
    my ( $filename, $dry_run, $error_id ) = @_;

    return if ($dry_run);

    # Handles multiple files ...
    if ( $filename =~ /\*/ ) {
        my @files = glob($filename);
        if ( @files > 0 ) {
            print STDERR localtime()
              . " $error_id : Input files $filename exist, continuing...\n";
            return 0;
        }
        else {

            # No input file, fall through to last return
        }
    }

    # ... or single files
    elsif ( -e $filename ) {
        print STDERR localtime()
          . " $error_id : Input file $filename exists, continuing...\n";
        return 0;
    }
    else {

        # No input file, fall through to last return
    }
    print STDERR localtime()
      . " $error_id : Input file $filename does not exist, exiting...\n";

    return 1;
}

sub try_rm {
    my ( $filename, $dry_run, $pm, $error_id, $silent ) = @_;

    if ( !defined($silent) ) { $silent = 0; }

    # Handles multiple files ...
    if ( $filename =~ /\*/ ) {
        my @files = glob($filename);
        foreach my $globbedfile (@files) {

            # Silent removal, to avoid hundreds of messages
            try_system( "rm $globbedfile", $dry_run, $pm, $error_id, 1 );
        }
        print STDERR localtime() . " $error_id : Removed $filename\n";
    }

    # ... or single files
    elsif ( -e $filename ) {
        try_system( "rm $filename", $dry_run, $pm, $error_id, $silent );
    }
    return;
}

sub try_system {
    my ( $command, $dry_run, $pm, $error_id, $silent ) = @_;

    if ( !defined($error_id) ) { $error_id = ""; }
    if ( !defined($silent) )   { $silent   = 0; }

    if ( !$silent ) {
        print STDERR localtime() . " $error_id : $command\n";
    }
    return if $dry_run;
    system($command);
    if ( $? == -1 ) {
        print STDERR localtime()
          . " $error_id : Failed to execute \"$command\": $!\n";
        if ( defined $pm ) { $pm->finish(1); }
    }
    elsif ( $? & 127 ) {
        my $error = "\"$command\" died with signal " . ( $? & 127 ) . "\n";
        print STDERR localtime() . " $error_id : $error";
        if ( defined $pm ) { $pm->finish(1); }
    }
    elsif ( ( $? >> 8 ) > 0 ) {
        print STDERR localtime()
          . " $error_id : \"$command\" failed, exiting\n";
        if ( defined $pm ) { $pm->finish(1); }
    }
    return;
}
