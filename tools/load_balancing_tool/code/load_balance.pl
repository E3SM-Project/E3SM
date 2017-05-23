#!/usr/bin/perl

$cpus_array = $ARGV[0];
$rdir = $ARGV[1];
$current = $ARGV[2];

# Copy over the first case's env_mach_pes file to modify
# with new layouts

#my $case = 't'.$ARGV[3].'1_'.$ARGV[4];
#my $casebase =  $ARGV[5].'_'.$ARGV[6].'_'.$ARGV[7];
#my $caseDir_1 =  $ARGV[8].'/'.$case.'_'.$casebase;

my $case = $ARGV[3];
my $casebase = $ARGV[4];
my $caseDir_1 = $ARGV[5];
my $fv_constraints = $ARGV[6];
my $full_res = $ARGV[7];

system("cp $caseDir_1/env_mach_pes.xml $rdir \n");
system("cp $caseDir_1/env_mach_pes.xml . \n");
system("cp $caseDir_1/env_run.xml . \n");
system("cp $caseDir_1/env_build.xml . \n");
system("cp $caseDir_1/env_case.xml . \n");
system("cp $caseDir_1/env_archive.xml . \n");
system("cp $caseDir_1/xmlchange $rdir \n");
system("cp -r $caseDir_1/Tools $rdir \n");

# Copy over the timing files from the scaling curve runs
my $test_file = "$rdir/test_list.out";
open my $list_handle, '<', $test_file;
chomp(my @test_list = <$list_handle>);
foreach my $f (@test_list){
  my $file_string = $f."/timing/";
  if (-e $file_string){
    system("cp $f/timing/*_timing.* $rdir \n");
  }
}

system("perl $current/code/get_cesm_times.pl $rdir \n");
@cpus = split(/,/, $cpus_array);
foreach my $c (@cpus){
  system("perl $current/code/create_dataFile.pl $rdir $c $current \n");
  system("mv model.data $current/code/ \n");
  system("/usr/bin/python $current/code/merge.py $current/code/model.mod $current/code/model.data $current/code/model.run >& $rdir/job.xml \n");
  system("/usr/bin/python $current/code/neos.py $rdir/job.xml >& $rdir/job.out \n");
  system("tail -n 6 $rdir/job.out >& $rdir/minmax_times.txt \n");

  # If trying to find a FV layout, run twice to make sure the the ATM pe task count can be decomposed
  if ($fv_constraints == 1){
   my $minmax = "$rdir/minmax_times.txt";
   open (F1, $minmax);
   my $cpus;
   while(<F1>) {
    if (/atm \s+(\d+)/){
      $cpus = $1;
    }
   }
   my $res = substr($full_res,0,3);
   my $codeDir = "$current/code/";
   system("perl $current/code/fv_second_pass.pl $res $cpus $codeDir $current $rdir \n");
  }

  my $minmax = "$rdir/minmax_times.txt";

  open (F, $minmax);

  while(<F>) {

    if (/atm \s+(\d+)/){
      my $cpus = $1;
      system("$rdir/xmlchange -file env_mach_pes.xml -id NTASKS_ATM -val $cpus \n");
      system("$rdir/xmlchange -file env_mach_pes.xml -id ROOTPE_ATM -val 0 \n");

      system("$rdir/xmlchange -file env_mach_pes.xml -id NTASKS_CPL -val $cpus \n");
      system("$rdir/xmlchange -file env_mach_pes.xml -id ROOTPE_CPL -val 0 \n");

      system("$rdir/xmlchange -file env_mach_pes.xml -id ROOTPE_OCN -val $cpus \n");

      system("$rdir/xmlchange -file env_mach_pes.xml -id NTASKS_GLC -val 1 \n");
      system("$rdir/xmlchange -file env_mach_pes.xml -id ROOTPE_GLC -val 0 \n");
    }

    if (/lnd \s+(\d+)/){
      my $cpus = $1;
      system("$rdir/xmlchange -file env_mach_pes.xml -id NTASKS_LND -val $cpus \n");
      system("$rdir/xmlchange -file env_mach_pes.xml -id ROOTPE_LND -val 0 \n");

      system("$rdir/xmlchange -file env_mach_pes.xml -id NTASKS_ROF -val $cpus \n");
      system("$rdir/xmlchange -file env_mach_pes.xml -id ROOTPE_ROF -val 0 \n");

      system("$rdir/xmlchange -file env_mach_pes.xml -id ROOTPE_ICE -val $cpus \n");
      system("$rdir/xmlchange -file env_mach_pes.xml -id ROOTPE_WAV -val $cpus \n");
    }

    if (/ice \s+(\d+)/){
      my $cpus = $1;
      system("$rdir/xmlchange -file env_mach_pes.xml -id NTASKS_ICE -val $cpus \n");
      system("$rdir/xmlchange -file env_mach_pes.xml -id NTASKS_WAV -val $cpus \n");
    }

    if (/ocn \s+(\d+)/){
      my $cpus = $1;
      system("$rdir/xmlchange -file env_mach_pes.xml -id NTASKS_OCN -val $cpus \n");
    }

  }
  close(F);
  system("mv $rdir/minmax_times.txt $rdir/minmax_times_$c.txt");
  system("cp env_mach_pes.xml $rdir/env_mach_pes_$c.xml");
}

system("rm -f $current/env_* \n");
system("gnuplot $current/code/cesm_scaling.gplot \n");
system("mv *.gif $rdir \n");
system("mv *.dat $rdir \n");


