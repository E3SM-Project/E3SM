#!/bin/env perl
use strict;

my $timingDir = $ARGV[0];

opendir(D,$timingDir);
my @files = grep /_timing/, readdir(D);
closedir(D);
my $component;
open(T,">tp.dat");
open(C,">cost.dat");
foreach my $file (@files){
    my $full_fn = $timingDir . $file;
    open(F,"$full_fn") or die "could not open $full_fn";
    my $tasks;
    my $threads;
    my $tp;
    my $mc;
    my $pes;
    foreach(<F>){
	if(/(\w+) = (\w+)\s+\d+\s+\d+\s+(\d+)\s+x\s+(\d+)/){
	    my $comp = $1;
	    $comp =~ tr/a-z/A-Z/;
	    $tasks->{$comp}=$3;
	    $threads->{$comp}=1;
	}
	if(/(\w+) Run Time:\s+(\d+\.\d+) seconds \s+(\d+\.\d+) seconds/){
	    my $comp = $1;
	    next if ($comp eq 'TOT' or $comp eq 'GLC');
            $component->{$comp}{$tasks->{$comp}}{$threads->{$comp}} = $3;
	}
        if(/ Model Throughput:\s+(\d+\.\d+)/){
	    $tp .= $1;
	}
        if(/Model Cost:\s+(\d+\.\d+)/){
            $mc .= $1;
        }
        if(/total pes active           :\s+(\d+)/){
	    $pes .= $1;
        }
    }
    print T "$pes $tp \n" ;
    print C "$pes $mc \n" ;
    close(F);
}

close(T);
close(C);


foreach my $comp (keys %$component){
    open(F,">$comp.dat");
    print F "Tasks seconds/model-day \n";
    foreach my $tasks (sort numerically keys %{$component->{$comp}}){
        my $nodes = $tasks;
	print F "$nodes $component->{$comp}{$tasks}{1} $component->{$comp}{$tasks}{2} \n";
    }
    close(F);
}

sub numerically{ $a <=> $b; }
