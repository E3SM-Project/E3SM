#!/usr/bin/env perl
#BSUB -P P93300075            # project code
#BSUB -W 0:40                # wall-clock time (hrs:mins)
#BSUB -n 64            # number of tasks in job         
#BSUB -R "span[ptile=16]"    # run 16 MPI tasks per node
#BSUB -J testpio                # job name
#BSUB -o testpio.%J.out         # output file name in which %J is replaced by the job ID
#BSUB -e testpio.%J.err         # error file name in which %J is replaced by the job ID
#BSUB -q small             # queue
#BSUB -a poe
##BSUB -XF
##BSUB -a tv
#BSUB -N
#BSUB -x
use strict;
use File::Copy ;
require Utils;


my $host="yellowstone";
$host = Utils->host() unless(defined $host);
print "host = $host\n";
Utils->loadmodules("$host");



my $piosrc="$ENV{HOME}/pio2_0";
my $testdir="/glade/scratch/$ENV{USER}/pio2test/pio";

opendir(TNL,"$piosrc/testpio/namelists");
my @namelists = grep(/testpio_in.*\d$/,readdir(TNL));
closedir(TNL);
$ENV{LD_LIBRARY_PATH}="$ENV{LD_LIBRARY_PATH}:/glade/apps/opt/hdf5-mpi/1.8.11/intel/13.1.2/lib";
my $passcnt=0;
my $failcnt=0;

open(T,">$testdir/testpio/TestStatus");
chdir "$testdir/unittests";
copy("$piosrc/unittests/input.nl","input.nl");
print T "Running unittests ... ";
system("mpirun.lsf ./piotest 1> unittest.out 2>&1");
open(F,"unittest.out");
my $cnt = grep /PASSED unit testing/ , <F>;
close(F);
if($cnt>0){
    $passcnt++;
    print "PASS \n";
    print T "PASS \n";
}else{
    $failcnt++;
    print "FAIL \n";
    print T "FAIL \n";
}



foreach my $nl (sort @namelists){
    chdir "$testdir/testpio";
    $nl =~ /testpio_in\.(.*)/;
    my $test = "test.$1";

    next if($test =~ /\.b/);
    next if($test =~ /\.a/);
    next if($test =~ /\.n4/);

    print T "Running test $1 ... ";
    print "Running test $1 ... ";
    mkdir $test;
    chdir $test;
    open(F,"$piosrc/testpio/namelists/$nl");
    my @namelist = <F>;
    close(F);

    open(G,">testpio_in");
    foreach(@namelist){
	if(/nx_global/){
	    print G "  nx_global = 256\n";
	    next;
	}
	if(/ny_global/){
	    print G "  ny_global = 128\n";
	    next;
	}
	if(/maxiter/){
	    print G " maxiter=1\n";
	    next;
	}
	print G $_;
    }
    close(G);
#    copy("$piosrc/testpio/namelists/$nl","testpio_in");
    mkdir "none";
    system("mpirun.lsf ../testpio > $test.out");
    open(F,"$test.out");
    my $cnt = grep /testpio completed successfully/ , <F>;
    close(F);

    if($cnt>0){
	$passcnt++;
	print "PASS \n";
	print T "PASS \n";
    }else{
	$failcnt++;
	print "FAIL \n";
	print T "FAIL \n";
    }
    
}
close(T);
