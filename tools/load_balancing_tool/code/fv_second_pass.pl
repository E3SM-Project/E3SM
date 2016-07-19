#!/usr/bin/env perl

my $res = $ARGV[0];
my $target = $ARGV[1];
my $codeDir = $ARGV[2];
my $current = $ARGV[3];
my $rdir = $ARGV[4];

my $count_file = $codeDir."/".$res."_peList.txt";
open my $list_handle, '<', $count_file;
chomp(my @test_list = <$list_handle>);

$index = 0;
I: foreach my $count (@test_list){
     if ($count >= $target){
       last I;
     } else {
       $index = $index + 1;
     }
   }

if ($test_list[$index] != $target){
  if ($index < 10){
    $lower = 0;
  } else {
    $lower = $index - 10;
  }
  $upper = $index + 9;

  my $dataFile = $codeDir."/model.data";
  open(MD, ">>$dataFile");

  print MD "\n";
  print MD "param AtmPart := \n";
  $counter = 0;
  for ($i=$lower; $i<=$upper; ++$i){
    $counter = $counter + 1 if exists $test_list[$i];
    print MD "         $counter $test_list[$i] \n" if exists $test_list[$i];
  }
  print MD ";\n";
  print MD "\n";
  print MD "param AtmMaxPart := ".$counter.";\n";
  print MD "\n";
  close(MD);

  system("/usr/bin/python $current/code/merge.py $current/code/fv_model.mod $current/code/model.data $current/code/fv_model.run >& $rdir/job.xml \n");
  system("/usr/bin/python $current/code/neos.py $rdir/job.xml >& $rdir/job.out \n");
  system("tail -n 6 $rdir/job.out >& $rdir/minmax_times.txt \n");
}

