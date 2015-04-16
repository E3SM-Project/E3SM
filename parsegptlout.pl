#!/usr/bin/perl

use strict;

our ($verbose) = 0;   # output verbosity
our ($maxval);        # max value across all processes/threads
our ($minval);        # min value across all processes/threads
our ($sum);           # total
our ($nval);          # number of entries found across all processes/threads
our ($totcalls);      # number of calls found across all processes/threads
our ($numthreads);    # number of threads

my ($fn);             # file name
my ($fnroot) = "timing";
my ($target);         # region to search for
my ($arg);            # element of @ARGV
my ($started);        # flag indicates initial "Stats for thread" found
my ($thread);         # thread number
my ($threadmax);      # thread number for max value
my ($threadmin);      # thread number for min value
my ($task);           # MPI task index (= ntask at loop completion)
my ($taskmax);        # task index for max value
my ($taskmin);        # task index for min value
my ($line);           # input line read in from file
my ($idx);            # index
my ($hidx);           # heading index
my ($mean);           # mean value
my ($found);          # flag indicates region name found
my ($totposs);        # number of threads * number of tasks
my ($heading) = "Wallclock"; # heading (column) to search for 

my (@vals);           # values for region
my (@headinglist);    # list of headings in input files

# Parse arg list

while ($arg = shift (@ARGV)) {
    if ($arg eq "-f") {
	$fnroot = shift (@ARGV);    # change root of file name
    } elsif ($arg eq "-h") {
	$heading = shift (@ARGV);   # change heading to search for
    } elsif ($arg eq "-v") {
	$verbose = 1;
    } else {
	if ( ! defined $target ) {
	    $target = "$arg";         # region name
	    chomp ($target);
	} else {
	    die_usemsg ("Unknown argument $arg\n");
	}
    }
}

die_usemsg ("Target region name not defined\n") if ( ! defined $target );
&initstats();     # Initialize stats
$found = 0;       # false

# Loop through output files

for ($task = 0; -e "${fnroot}.$task"; $task++) {
    $fn = "${fnroot}.$task";
    open (FILE, "<$fn") || die ("Can't open $fn for reading\n");
    $started = 0;

# Read all the lines in the file, looking for "Stats for thread", followed by
# thre region name

    while ($line = <FILE>) {
	chomp ($line);
	if ($line =~ /^Stats for thread (\d*):/) {
	    $started = 1;
	    $thread = $1;
	    $numthreads = $thread if ($thread > $numthreads);

# Next line contains the headings. Parse for later printing
# Chop off leading whitespace--in can foul up later parsing

	    $line = <FILE>;
	    chomp ($line);
	    if ($line =~ /^\s+(.*)$/) {
		$line = $1;
	    }
	    @headinglist = split (/\s+/, $line);
	    for ($hidx = 0; $hidx <= $#headinglist; $hidx++) {
		last if ($headinglist[$hidx] eq $heading);
	    }
	    if ($hidx > $#headinglist) {
		die ("Heading $heading not found in $fn. Giving up\n");
	    }
	} elsif ($started && ($line =~ /^[* ]\s*${target}\s+(.*)$/)) {
	    $found = 1;
	    @vals = split (/\s+/, $1);
	    ($#vals >= $hidx) || die ("$heading not found in input:\n$1\n");
	    $totcalls += $vals[0];
	    $sum += $vals[$hidx];
	    $nval++;

	    if ($vals[$hidx] > $maxval) {
		$maxval = $vals[$hidx];
		$taskmax = $task;
		$threadmax = $thread;
	    }
	    if ($vals[$hidx] < $minval) {
		$minval = $vals[$hidx];
		$taskmin = $task;
		$threadmin = $thread;
	    }
	    $started = 0;
	    next;   # Look for next "Stats for thread"
	}
    }
}

die ("Found no occurrences of $target in any of $task files\n") if ( ! $found );

print (STDOUT "Searched for region $target\n");
$numthreads++;  # convert from 0-based to 1-based
print (STDOUT "Found $totcalls calls across $task tasks and $numthreads threads per task\n");
$totposs = $numthreads * $task;
print (STDOUT "$nval of a possible $totposs tasks and threads had entries for $target\n");
print (STDOUT "Heading is $heading\n");
printf (STDOUT "Max   = %.3g on thread %d task %d\n", $maxval, $threadmax, $taskmax);
printf (STDOUT "Min   = %.3g on thread %d task %d\n", $minval, $threadmin, $taskmin);
$mean = $sum / $nval;
printf (STDOUT "Mean  = %.3g\n", $mean);
printf (STDOUT "Total = %.3g\n", $sum);

exit 0;

sub initstats {
    our ($verbose);
    our ($maxval);
    our ($minval);
    our ($sum);
    our ($nval);
    our ($totcalls);
    our ($numthreads);

    $totcalls = 0;
    $numthreads = 0;
    $minval = 9.99e19;
    $maxval = -9.99e19;
    $nval = 0;
    $sum = 0.;
    $taskmax = -1;
    $taskmin = -1;
    $threadmax = -1;
    $threadmin = -1;
}

sub die_usemsg {
    defined $_[0] && print (STDOUT "$_[0]");
    print (STDOUT "Usage: $0 [-v] [-f file-root] [-h heading] region\n",
	   " -f file-root => look for files named <file-root>.<taskid> (default file-root is 'timing')\n",
	   " -h heading   => use <heading> as the default search target (default is Wallclock)\n",
	   " -v           => verbose\n",
	   " region       => region name to search for\n");
    exit 1;
}
