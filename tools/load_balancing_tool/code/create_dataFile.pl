#!/bin/env perl
use strict;

my $timingDir = $ARGV[0];
my $CPUS = $ARGV[1];
my $current = $ARGV[2];
my $codeDir = "$current/code/";

opendir(D,$timingDir);
my @files = grep /cesm_timing/, readdir(D);
closedir(D);
my $component;
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
    }
    close(F);
}


open(ATM_FILE, "<tp.dat") or die "Could not open file: $!";
my $Tcount = 0;
while (<ATM_FILE>) {
  $Tcount++;
}

my $x = 1;
my $TcountS;
while ($x <= $Tcount){
  $TcountS = $TcountS." ".$x;
  $x++;
}

open(S,">all.dat");
print S "Tasks: \n";

my $dataFile = "$codeDir/model.data";

open(I, ">model.data");

print I "data;\n\n";
print I "param D := $Tcount;\n";
print I "param CPUS := $CPUS;\n";
print I "param CPN := 16;\n";
print I "param Tsync := 3.0;\n";
print I "param Etarget := 0.5;\n";
print I "param MinNodes := 64;\n";
print I "param MaxNodes := 48160;\n";
print I " \n";
print I "param rawx:  $TcountS :=\n";

foreach my $comp (keys %$component){
    print S "$comp  ";
    if ($comp eq 'LND' or $comp eq 'ICE' or $comp eq 'ATM' or $comp eq 'OCN'){
      my $compLC = lc $comp;
      print I "          \'$compLC\' ";
    }
    foreach my $tasks (sort numerically keys %{$component->{$comp}}){
        my $nodes = $tasks;
        print S "$nodes  ";
        if ($comp eq 'LND' or $comp eq 'ICE' or $comp eq 'ATM' or $comp eq 'OCN'){
          print I " $nodes";
        }
    }
    print S "\n";
    if ($comp eq 'LND' or $comp eq 'ICE' or $comp eq 'ATM' or $comp eq 'OCN'){
      print I "\n";
    }
}

print I ";\n";
print I "\n";
print I "param rawy:  $TcountS :=\n";

print S "Timings: \n";
foreach my $comp (keys %$component){
    print S "$comp  ";
    if ($comp eq 'LND' or $comp eq 'ICE' or $comp eq 'ATM' or $comp eq 'OCN'){
      my $compLC = lc $comp;
      print I "          \'$compLC\' ";
    }
    foreach my $tasks (sort numerically keys %{$component->{$comp}}){
        my $nodes = $tasks;
        print S "$component->{$comp}{$tasks}{1}  ";
        if ($comp eq 'LND' or $comp eq 'ICE' or $comp eq 'ATM' or $comp eq 'OCN'){
          print I " $component->{$comp}{$tasks}{1}";
        }
    }
    print S "\n";
    if ($comp eq 'LND' or $comp eq 'ICE' or $comp eq 'ATM' or $comp eq 'OCN'){
      print I "\n";
    }
}
close(S);

print I ";\n";
print I "\n";

close(I);

sub numerically{ $a <=> $b; }
