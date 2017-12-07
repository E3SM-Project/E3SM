#!/usr/bin/env perl
use strict;
use File::Copy ;
#BSUB -P P93300606            # project code
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

my $piosrc="$ENV{HOME}/pio_trunk";
my $testdir="/glade/scratch/$ENV{USER}/piotest/pio.all/pio";

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
system("mpirun.lsf ./piotest > unittest.out");
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
    print T "Running test $1 ... ";
    mkdir $test;
    chdir $test;
    copy("$piosrc/testpio/namelists/$nl","testpio_in");
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
