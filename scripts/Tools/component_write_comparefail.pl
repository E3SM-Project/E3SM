#! /usr/bin/env perl
use strict;

if ($#ARGV == -1) {
    die " ERROR component_write_comparefail: must specify a caseroot";
}
my $rundir = $ARGV[0];

print "---summarizing more details of test failures if any: --- \n\n" ;

my $rv = 0;
my @cprnc_files = glob("$rundir/*cprnc.out");
foreach my $file (@cprnc_files) {

    my $ndiffs     = `grep RMS $file | wc -l`;
    my $nfilldiffs = `grep FILLDIFF $file | wc -l`;

    if (($ndiffs > 0) || ($nfilldiffs > 0)) {
	print "$file had the following fields that are NOT b4b  \n\n";
        $rv = 1;
    }
    
    if ($ndiffs > 0) {
	open(fhin, "<$file") or die "Could not open file $file to read";
	while (my $line = <fhin>) {
	    if ($line =~ /RMS\s+(\S+)\s+(\S+)/) {
		print "  $line ";
	    }
	}
	close(fhin);
    }
    
    if ($nfilldiffs > 0) {
	open(fhin, "<$file") or die "Could not open file $file to read";
	while (my $line = <fhin>) {
	    if ($line =~ /FILLDIFF\s+(\S+)/) {
		print "  $line ";
	    }
	}
	close(fhin);
    }
}

my $sdate = `date +"%Y-%m-%d %H:%M:%S"`;

print "\n---finished summarizing more details of test failures: ---- \n\n ";
print "test completed $sdate\n"; 

exit($rv);
