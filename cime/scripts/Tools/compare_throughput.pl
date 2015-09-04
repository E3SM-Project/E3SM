#!/usr/bin/perl
# Compare throughtput of current run with throughput of previous baseline
use Getopt::Long;

my $match = '';
my $rounded_1 = '';
my $rounded_2 = '';
my $maxtputdecr_ypd =.2; 
my $maxtputdecr_pct = 2; 
my $performance_decrease = 0;

my $result = GetOptions("file1=s" => \$cpllog1,             # cpl.log file (required)  
                        "file2=s" => \$cpllog2,             # cpl.log file (required) 
                        "tput_ypd:f"  => \$maxtputdecr_ypd, # max allowable throughput decrease (optional) 
                        "tput_pct:f"  => \$maxtputdecr_pct # max allowable throughput decrease (optional) 
                        );   

if (! result ) { 
   die "Usage: compare_throughput -file1 cpl.log -file2 baseline.cpl.log -m_ypd maxtputdecr (optional) -m_pct maxtputdecrpct(optional)" 
}

if (! defined $cpllog1 || ! defined $cpllog2  ) { 
   print "Usage: compare_throughput -file1 cpl.log -file2 baseline.cpl.log -m_ypd maxtputdecr (optional) -m_pct maxtputdecrpct(optional) \n"; 
   exit 1; 
}

open my $logfile_1, $cpllog1 or die "Could not open $cpllog1: $!";
open my $logfile_2, $cpllog2 or die "Could not open $cpllog2: $!";
 
print "\n IN COMPARE_THROUGHPUT\n";

while( my $line = <$logfile_1>)  {   
   chomp($line);
   # find throughput for current run
   if ( my ($match) =  $line =~ m/# simulated years \/ cmp-day = +(\d+.\d+) =*.*$/ ) {
    print "\nline = $line\n";  
    $rounded_1 = sprintf("%.3f", $match);
    print "rounded_1 = $rounded_1\n";
     }
}

while( my $line = <$logfile_2>)  {   
   chomp($line);
   # find max highwater mark for baseline run
   if ( my ($match) =  $line =~ m/# simulated years \/ cmp-day = +(\d+.\d+) =*.*$/ ) { 
     $rounded_2 = sprintf("%.3f", $match); 
     print "rounded_2 = $rounded_2\n";
   }   
}

if ($rounded_2 > $rounded_1 ) { 
    if ($rounded_2 - $rounded_1 > $maxtputdecr_ypd)  
     {
     my $tput_decrease = $rounded_2 - $rounded_1; 
     print "Throughput decrease greater than $maxtputdecr_ypd ypd\n";
     print "\ntput_decr = $tput_decrease\n"; 
     $performance_decrease = 1;
    }
    if  ( ( ( ( $rounded_2 - $rounded_1 ) / $rounded_2 ) * 100) > $maxtputdecr_pct ){  
     print "Throughput decrease greater than $maxtputdecr_pct percent\n";
     my $percent_decrease = ( ( ( $rounded_2 - $rounded_1 ) / $rounded_2 ) * 100); 
     print "\ntput_percent_decr = $percent_decrease\n";
     $performance_decrease = 1;
   }
}
if  ($performance_decrease == 1) { 
    print "\nFAIL\n";
   }
else {
    print "throughput passes";
    print "\nPASS\n"; 
}
