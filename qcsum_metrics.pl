#!/usr/bin/perl
# last update: 5-25-2018
use Getopt::Long;

my ($prefix,
    $qcfolder,
    $qcsum,
    $pipeline_version,
    $platform,
    $pass_min_align_pct,
    $fail_min_align_pct,
    $covered,
    $pass_min_roi_pct,
    $fail_min_roi_pct,
    $pass_min_avgcov,
    $fail_min_avgcov,
    $pass_min_reads,
    $fail_min_reads,
    # $pass_min_contem,
    $capture,
    $capture_version,
  ) = ("") x 16;

GetOptions(
	"prefix=s" => \$prefix,
	"qcfolder=s"  => \$qcfolder,
    "pipeline_version=s" => \$pipeline_version,
    "platform=s" => \$platform,
    "pass_min_align_pct=s" =>\$pass_min_align_pct,
    "fail_min_align_pct=s" =>\$fail_min_align_pct,
    "covered=s" =>\$covered,
    "pass_min_roi_pct=s" =>\$pass_min_roi_pct,
    "fail_min_roi_pct=s" =>\$fail_min_roi_pct,
    "pass_min_avgcov=s" =>\$pass_min_avgcov,
    "fail_min_avgcov=s" =>\$fail_min_avgcov,
    "pass_min_reads=s" =>\$pass_min_reads,
    "fail_min_reads=s" =>\$fail_min_reads,
    "capture=s" => \$capture,
    "capture_version=s" =>\$capture_version,
          );

# DEFINE ACCEPTABLE CRITERIA FOR DIFFERENT METRICS
# only define pass and fail thresholds. Anything in between will issue a warning

my $hs = $qcfolder."/".$prefix.".hsm.txt";

# define picard hs metrics variables
my ($BAIT_SET,
$TOTAL_READS, 
$PCT_PF_UQ_READS_ALIGNED, 
$ON_BAIT_BASES, 
$NEAR_BAIT_BASES, 
$PCT_ON_BAIT, 
$ON_TARGET_BASES, 
$PCT_SELECTED_BASES, 
$MEAN_BAIT_COVERAGE, 
$MEAN_TARGET_COVERAGE, 
$MEDIAN_TARGET_COVERAGE, 
$MAX_TARGET_COVERAGE, 
$MIN_TARGET_COVERAGE, 
$FOLD_80_BASE_PENALTY, 
$PCT_TARGET_BASES_1X, 
$PCT_TARGET_BASES_20X, 
$PCT_TARGET_BASES_100X, 
$PCT_TARGET_BASES_500X) = ("") x 18;

# Open picard hs metrics file
open (DATA, "$hs");
readline(DATA);
$x = 0;
while (<DATA>){
    $rm = $_;
    chomp $rm;
    $x++;
    my @line = split(/\t/, $rm);

     # Check if we are at the line containing stats yet
     ## with Picard v2, the column indices have been changed.
     ## changes have been made accordingly below

    if (($x==7)){
    	$BAIT_SET = $line[0];
    	$TOTAL_READS = $line[22];
    	$PCT_PF_UQ_READS_ALIGNED = $line[32]*100;
    	$ON_BAIT_BASES = $line[13];
    	$NEAR_BAIT_BASES = $line[3];
        $PCT_ON_BAIT = 100-($line[7]*100);
        $ON_TARGET_BASES = $line[29];
        $PCT_SELECTED_BASES = $line[6]*100;
    	$MEAN_BAIT_COVERAGE = $line[9];
        $MEAN_TARGET_COVERAGE = $line[33];
        $MEDIAN_TARGET_COVERAGE = $line[34];
        $MAX_TARGET_COVERAGE =  $line[35];
        $MIN_TARGET_COVERAGE = $line[36];
        $FOLD_80_BASE_PENALTY = $line[44];
        $PCT_TARGET_BASES_1X = $line[45]*100;
        $PCT_TARGET_BASES_20X = $line[48]*100;
        $PCT_TARGET_BASES_100X = $line[52]*100;
        $PCT_TARGET_BASES_500X = $line[54]*100;
	}
}

# DETERMINE ALIGNMENT QC (PASS, WARN, FAIL)
my $alignqc = "PASS";

if($PCT_PF_UQ_READS_ALIGNED < $pass_min_align_pct)
{
    $alignqc = "WARN";
}
if($PCT_PF_UQ_READS_ALIGNED < $fail_min_align_pct)
{
    $alignqc = "FAIL";
}
if($TOTAL_READS < $pass_min_reads)
{
    $alignqc = "WARN";
}
if($TOTAL_READS < $fail_min_reads)
{
    $alignqc = "FAIL";
}


# DETERMINE COVERAGE QC (PASS, WARN, FAIL)
my $cov_th = "";
if ($covered == 20)
{
    $cov_th = $PCT_TARGET_BASES_20X;
}
elsif ($covered == 100)
{
    $cov_th = $PCT_TARGET_BASES_100X;
}
elsif ($covered == 500)
{
    $cov_th = $PCT_TARGET_BASES_500X;
}
my $covqc = "PASS";
if($cov_th < $pass_min_roi_pct)
{
    $covqc = "WARN";
}
if($cov_th < $fail_min_roi_pct)
{
    $covqc = "FAIL";
}
if($MEAN_TARGET_COVERAGE < $pass_min_avgcov)
{
    $covqc = "WARN";
}
if($MEAN_TARGET_COVERAGE < $fail_min_avgcov)
{
    $covqc = "FAIL";
}
# Print prefix (project,sample, etc. whatever is passed to the script as -p)
print MYFILE "$prefix\t";
# Now we have all data ready to be written in the qc summary file
my $qcsumfile = $qcfolder."/".$prefix.".qcsum.txt";

#open the qcsum file to write
open (MYFILE, ">". $qcsumfile);

# Print headers
print MYFILE "#Sample";
print MYFILE "\tSequencing_Platform\tPipeline_version";
print MYFILE "\tAlignment_QC\tCoverage_QC";
print MYFILE "\tTotal_Reads\t%Reads_Aligned\tCapture\tAvg_Capture_Coverage";
print MYFILE "\t%On/Near_Bait_Bases\t%On_Bait_Bases\tFOLD_80_BASE_PENALTY\tAvg_ROI_Coverage";
print MYFILE "\tMEDIAN_ROI_COVERAGE\tMAX_ROI_COVERAGE";
print MYFILE "\t%ROI_1x\t%ROI_20x\t%ROI_100x\t%ROI_500x";
print MYFILE "\n";

#print sample information
print MYFILE "$prefix";
print MYFILE "\t$platform\t$pipeline_version";
print MYFILE "\t$alignqc\t$covqc";
print MYFILE "\t$TOTAL_READS\t$PCT_PF_UQ_READS_ALIGNED\t$capture\t$MEAN_BAIT_COVERAGE";
print MYFILE "\t$PCT_SELECTED_BASES\t$PCT_ON_BAIT\t$FOLD_80_BASE_PENALTY\t$MEAN_TARGET_COVERAGE";
print MYFILE "\t$MEDIAN_TARGET_COVERAGE\t$MAX_TARGET_COVERAGE";
print MYFILE "\t$PCT_TARGET_BASES_1X\t$PCT_TARGET_BASES_20X\t$PCT_TARGET_BASES_100X\t$PCT_TARGET_BASES_500X";
print MYFILE "\n";

close MYFILE or warn $! ? "Error closing the qcsum file $!"
                   : "Exit status $? from the file";
