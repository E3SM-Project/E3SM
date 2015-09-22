#!/usr/bin/perl
use Getopt::Long;

my $match = '';
my $mem = '';
my $rounded_new = '';
my $rounded_old = '';
my $num_dates = 0;
my @dates = ''; 
my $date = '';
my $memleak = 0;
my $memincrfrombaseline =0; 
my @mem_highwater = ''; 
my $maxmemincrperday = 1; 
my $maxmempctincrbaseline = 20; 
my $maxmemleak = 0;

my $result = GetOptions("file1=s" => \$cpllog1,             # cpl.log file (required)  
                        "file2:s" => \$cpllog2,             # cpl.log file (optional) 
                        "m:f" => \$maxmemincrperday,        # max allowable memory increase per day (optional) 
                        "mbase:f" => \$maxmempctincrbaseline); # maximum allowable memory increase since baseline (optional)  


if (! defined $cpllog1  ) { 
   print "Usage: memory_check -file1 cpl.log -file2 baseline.cpl.log (optional) -m maxmemincr (optional) -mbase maxmemincrbase
(optional) \n";
   exit 1; 
}


#if baseline cpl.log given on command line then compare "pe max memory" for current run with baseline   
if ( $cpllog2 ) { 
   $memincrfrombaseline = &compare_maxmemory_baseline($cpllog1,$cpllog2,$maxmempctincrbaseline);  
}
else {
   $memleak = &check_memory_leak_current_run($cpllog1,$maxmemincrperday);  
}


if  ($memleak > 0) { 
    print "\nMemory leak. memory increased by $maxmemincrperday MB/day or more\n"; 
   }
if  ( $memincrfrombaseline > $maxmempctincrbaseline ) { 
    print "\npesmaxmem_incr = $memincrfrombaseline\n"; 
   } 

if ( ( $memincrfrombaseline > $maxmempctincrbaseline ) || ( $memleak > 0 ) )  { 
    print "\nFAIL\n"; 
}
else {
    print "not a memory leak";
    print "\nPASS\n"; 
}

#-----------------------------------------------------------------------------------
# check_memory_leak_current_run: Check for memory leak/creep in this test run.  
#                                For consecutive models days compare memory highwater    
#                                mark and determine whether there is an increase greater
#                                than $maxmemincreaseperday specified by the user on the 
#                                command line.  

sub check_memory_leak_current_run{
    my ($log1,$maxmemincreaseperday);                  # new private variables   
       ($log1,$$maxmemincreaseperday) = @_;            # assign input parameters names 

    open my $logfile_1, $log1 or die "Could not open $cpllog1: $!";

#find memory high water marks per day
    while( my $line = <$logfile_1>)  {   
       chomp($line);
       #  collect memory highwater marks per/day 
       if ( my ($date, $mem) =  $line =~ m/model date = (.{8}).*memory =(.*)MB.*highwater/ ) { 
       $rounded_new = sprintf("%.3f", $mem); 
       push (@dates, $date);
       push (@mem_highwater, $rounded_new);
       $num_dates++;
       }   
    }

# Check for memory leak in this run.  For consecutive model days determine  
# whether memory usage has increased by more than $maxmemincrperday 
    for my $x (1 .. ($num_dates-1)) 
    { 
      if ( ($dates[$x+1] - $dates[$x]) == 1 ) { 
        if ( $mem_highwater[$x+1] - $mem_highwater[$x] >= $maxmemincrperday ) {
         $memleak = 1;
         $memleakval = $mem_highwater[$x+1] - $mem_highwater[$x]; 
           if ($memleakval > $maxmemleak) { 
              $maxmemleak = $memleakval; 
           } 
         }
         else { 
           $memleak = 0;
         }  
        } 
      if ($memleak == 0 ) { 
        last; 
      }
    }
    
    if ($memleak == 1) { 
      print "memleak = $maxmemleak\n"; 
     }

    close $logfile_1;
    return $memleak; 
}
#-----------------------------------------------------------------------------------
# compare_maxmemory_baseline:  compare "pes max memory" value from current run
#                              for "pes max memory" from baseline.  If 'pes max memory' 
#                              usage has increased more than $mbase, specified
#                              by user on the command line, then test fails.  

sub compare_maxmemory_baseline { 
my ($log1,$log2,$memthreshold_for_base);           #  new private variables  
   ($log1,$log2,$memthreshold_for_base) =@_;       #  assign input parameters names 

my $memincr = 0; 

open my $logfile_1, $log1 or die "Could not open $cpllog1: $!";
open my $logfile_2, $log2 or die "Could not open $cpllog2: $!";

    
while( my $line = <$logfile_1>)  {   
   chomp($line);
   # find max highwater mark for current run
   if ( my ($match) =  $line =~ m/pes max memory highwater.*.MB. +(\d+.\d+) =*.*$/ ) {
    $rounded_new = sprintf("%.3f", $match);
     }
}

#compare 'pes max memory value' to baseline 'pes max memory' value 
  while( my $line = <$logfile_2>)  {   
     chomp($line);
     # find max highwater mark for baseline run
     if ( my ($match) =  $line =~ m/pes max memory highwater.*.MB. +(\d+.\d+) =*.*$/ ) { 
     $rounded_old = sprintf("%.3f", $match); 
       }   
  }


  if ($rounded_new > $rounded_old)  
   { 
    if  ( ( ( ( $rounded_new - $rounded_old ) / $rounded_old ) * 100) > $memthreshold_for_base ){
     my $mem_percent_incr = ( ( ( $rounded_new - $rounded_old ) / $rounded_old ) * 100);
     $mem_pnt_incr = sprintf("%.3f", $mem_percent_incr); 
     print "\mem_pnt_incr = $mem_pnt_incr\n";
   }

   } 

    close $logfile_2;
    return $mem_pnt_incr;
}
