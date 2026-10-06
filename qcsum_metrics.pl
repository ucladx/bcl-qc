#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

my %cfg;
my @keys = qw(prefix qcfolder pipeline_version platform fail_min_align_pct covered
              fail_min_roi_pct fail_min_avgcov fail_min_reads capture);
GetOptions(\%cfg, map { "$_=s" } @keys) or die "Invalid qcsum arguments\n";
for my $key (@keys) {
    die "Missing --$key\n" unless defined $cfg{$key} && length $cfg{$key};
}
for my $key (qw(fail_min_align_pct covered fail_min_roi_pct fail_min_avgcov fail_min_reads)) {
    die "Invalid --$key\n" unless looks_like_number($cfg{$key}) && $cfg{$key} !~ /nan|inf/i;
}
die "Unsupported coverage threshold: $cfg{covered}\n" unless $cfg{covered} =~ /^(20|100|250|500)$/;

# Read the HsMetrics table by column name; comments and column order may vary by Picard version.
my $input = "$cfg{qcfolder}/$cfg{prefix}.hsm.txt";
open my $in, '<', $input or die "Cannot read $input: $!\n";
my %metrics;
while (my $line = <$in>) {
    chomp $line;
    my @header = split /\t/, $line, -1;
    next unless grep { $_ eq 'TOTAL_READS' } @header;
    my $values = <$in> // die "Missing metrics row in $input\n";
    chomp $values;
    my @values = split /\t/, $values, -1;
    die "Incomplete metrics row in $input\n" unless @header == @values;
    @metrics{@header} = @values;
    last;
}
close $in;
sub metric {
    my ($name) = @_;
    my $value = $metrics{$name};
    die "Missing or invalid $name in $input\n"
        unless defined $value && looks_like_number($value) && $value !~ /nan|inf/i;
    return $value;
}
my $reads = metric('TOTAL_READS');
my $aligned = metric('PCT_PF_UQ_READS_ALIGNED') * 100;
my $avg_cov = metric('MEAN_TARGET_COVERAGE');
my $covered = metric("PCT_TARGET_BASES_$cfg{covered}X") * 100;
my $align_qc = $aligned < $cfg{fail_min_align_pct} || $reads < $cfg{fail_min_reads} ? 'FAIL' : 'PASS';
my $cov_qc = $covered < $cfg{fail_min_roi_pct} || $avg_cov < $cfg{fail_min_avgcov} ? 'FAIL' : 'PASS';
my @header = qw(Sample Sequencing_Platform Pipeline_version Alignment_QC Coverage_QC
                Total_Reads %Reads_Aligned Capture Avg_Capture_Coverage %On/Near_Bait_Bases
                %On_Bait_Bases FOLD_80_BASE_PENALTY Avg_ROI_Coverage MEDIAN_ROI_COVERAGE
                MAX_ROI_COVERAGE %ROI_1x %ROI_20x %ROI_100x %ROI_250x %ROI_500x);
my @values = ($cfg{prefix}, $cfg{platform}, $cfg{pipeline_version}, $align_qc, $cov_qc,
              $reads, $aligned, $cfg{capture}, metric('MEAN_BAIT_COVERAGE'),
              metric('PCT_SELECTED_BASES') * 100, 100 - metric('PCT_OFF_BAIT') * 100,
              metric('FOLD_80_BASE_PENALTY'), $avg_cov, metric('MEDIAN_TARGET_COVERAGE'),
              metric('MAX_TARGET_COVERAGE'), map { metric("PCT_TARGET_BASES_${_}X") * 100 } (1, 20, 100, 250, 500));
my $output = "$cfg{qcfolder}/$cfg{prefix}.qcsum.txt";
open my $out, '>', $output or die "Cannot write $output: $!\n";
print $out join(',', @header), "\n", join(',', @values), "\n";
close $out or die "Cannot close $output: $!\n";
