#!/usr/bin/perl
use strict;
my @files = @ARGV;
my $ofile = $files[0];
$ofile =~ s/(\d+)\.nc/\.nc/;
my $stride = $ENV{SAMPLE_STRIDE};
$stride=1 unless($stride>1);

if(-e "combined/$ofile" ){
    my $a1 = -A "$files[0]"; 
    my $a2 = -A "combined/$ofile";
    if($a2 < $a1){ 
	print "File combined/$ofile is up to date\n";
        exit(0);
    }
}
system("mkdir combined") unless(-d "combined");


foreach(@files){

  system("ncpdq -a latlonp,time $_ tmp$_");
  $_ ="tmp".$_;
}
system("ncrcat -d latlonp,,,$stride @files out.nc");
unlink @files;
system("ncpdq -a time,latlonp out.nc combined/$ofile"); 
unlink "out.nc";
