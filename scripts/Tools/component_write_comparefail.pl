#! /usr/bin/env perl
use strict;

if ($#ARGV == -1) {
    die " ERROR component_write_comparefail: must specify a caseroot and testlog input arguments";
}
my ($rundir, $testlog) = @ARGV;

open(fhout, '>>', $testlog) or die "Could not open output testlog $testlog' $!";

print fhout "---summarizing more details of test failures if any: --- \n\n" ;

my @cprnc_files = glob("$rundir/*cprnc.out");
foreach my $file (@cprnc_files) {

    my $ndiffs     = `grep RMS $file | wc -l`;
    my $nfilldiffs = `grep FILLDIFF $file | wc -l`;

    if (($ndiffs > 0) || ($nfilldiffs > 0)) {
	print fhout "$file had the following fields that are NOT b4b  \n";
	print fhout "\n";
    }
    
    if ($ndiffs > 0) {
	open(fhin, "<$file") or die "Could not open file $file to read";
	while (my $line = <fhin>) {
	    if ($line =~ /RMS\s+(\S+)\s+(\S+)/) {
		print fhout"  $line "; 
	    }
	}
	close(fhin);
    }
    
    if ($nfilldiffs > 0) {
	open(fhin, "<$file") or die "Could not open file $file to read";
	while (my $line = <fhin>) {
	    if ($line =~ /FILLDIFF\s+(\S+)/) {
		print fhout"  $line "; 
	    }
	}
	close(fhin);
    }
}

my $sdate = `date +"%Y-%m-%d %H:%M:%S"`;

print fhout "\n---finished summarizing more details of test failures: ---- \n\n ";
print fhout "test completed $sdate"; 

close (fhout);

exit(0);
