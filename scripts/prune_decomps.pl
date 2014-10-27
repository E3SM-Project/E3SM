#!/usr/bin/perl
use strict;

my $rundir = shift;

opendir(F,$rundir);
my @decompfiles = grep(/^piodecomp/,readdir(F));
closedir(F);

for(my $i=0; $i< $#decompfiles; $i++){
    my $file  = $decompfiles[$i];
    for(my $j=$i+1;$j<$#decompfiles;$j++){
	my $nfile = $decompfiles[$j];
	my $cmp = `cmp $file $nfile`;
#	print ">$cmp<\n";
	unless( $cmp =~ /differ/){
	    print "Removing $file\n";
	    unlink($file);
	    last;
	}
    }
}
    
