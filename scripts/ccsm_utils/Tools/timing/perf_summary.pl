#!/usr/bin/env perl

###This is a program to extract timing information based on ccsm4_0_alpha versions.

if($#ARGV == 0 ){
   $fhead=@ARGV[0];
}
else {
   print "timing_summary.pl input ERROR\n";
   print "Usage: timing_summary.pl input_file_name \n";
   print "ie     timing_suumary.pl /tmp/timing \n";
   print "          all files that match /tmp/timing* will be used as input\n";
   exit;
}

@files=<$fhead.*>;
$ntask = -1;
$nthrd = -1;

foreach $fin (@files) {
@lines=`cat $fin`;
for($i=0;$i<=$#lines;$i++){
  $cline=$lines[$i];
  $fword=$cline;
###  $fword =~ s/(^\s*\w+)(\s*[-y]\s+.*$)/\1/;
  $fword =~ s/^\*?(\s*\w+)(\s*[-y]\s+.*$)/\1/;
  $fword =~ s/\n//;
  if ($cline =~ m/\*+\s*PROCESS\s*(\d*)\s*\(\s*(\d*)\s*\)\s*\*+/) {
      $ntask=$2;
#      print "tcxx $ntask $nthrd\n";
  }
  if ($cline =~ m/\s*Stats for thread\s*(\d*)\s*:/) {
      $nthrd=$1;
#      print "tcxy $ntask $nthrd\n";
  }
  if ($cline =~ m/^\s*\*?\s*(\w+)\s*[-y]\s+(\d*\.*\d*)\s+-+\s+(\d*\.*\d*)\s+(\d*\.*\d*)\s+(\d*\.*\d*)\s+(\d*\.*\d*)/) {
#      print "tcx $1 $2 $3 $4 $5 $6\n";
       $ccnt=$2;
       $ctim=$3;
       $found=0;
       for ($j=0;$j<=$#vname;$j++){
          if ($found == 0) {
          if ($vname[$j] =~ /^\s*$fword$/) {
             if ($ctim < $vminx[$j]) {
                @vminx[$j]=$ctim;
                @vminp[$j]=$ntask;
                @vmint[$j]=$nthrd;
             }
             if ($ctim > $vmaxx[$j]) {
                @vmaxx[$j]=$ctim;
                @vmaxp[$j]=$ntask;
                @vmaxt[$j]=$nthrd;
             }
             if ($ccnt > $vcntx[$j]) {
		 @vcntt[$j]=$ccnt;
             }
             $found=1;
	  }
          }
       }
       if ($found == 0) {
          $j=$#vname+1;
          @vname[$j]=$fword;
          @vminx[$j]=$ctim;
          @vminp[$j]=$ntask;
          @vmint[$j]=$nthrd;
          @vmaxx[$j]=$ctim;
          @vmaxp[$j]=$ntask;
          @vmaxt[$j]=$nthrd;
          $vcntx[$j]=$ccnt;
#          print "tcx 1st $j $fword\n";
       }         
    }	
}
}

print " \n";
print "     Timer_Name                                  Count       Max_Time (Task,Thread)     Min_Time (Task,Thread)\n";
print "  ----------------                             ---------     -------- -------------     -------- -------------\n";
for ($j=0;$j<=$#vname;$j++){
    printf ("%-43s %12u %12.4f (%7u,%3u) %12.4f (%7u,%3u)\n",@vname[$j],@vcntx[$j],@vmaxx[$j],$vmaxp[$j],$vmaxt[$j],@vminx[$j],$vminp[$j],$vmint[$j]);
}
print " \n";
