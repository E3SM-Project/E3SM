#!/usr/bin/env perl
use strict;

my $nml		= shift;
my $basenml	= shift;
my $baseid	= shift;

open(F,$nml) or die "Could not open $nml";
my @nlfile = <F>;
close(F);

open(G,$basenml) or die "Could not open $basenml";
my @basenml = <G>;
close(G);

my @diffs1;
my @diffs2;
my @added;
my @removed;

my $line;
my $cnt1=$#nlfile;
my $cnt2=$#basenml;
my $shiftbl=1;
foreach $line (@nlfile){
    my $bline = shift(@basenml) if($shiftbl);
    $shiftbl=1;
    next if($line eq $bline);
    $bline=shift(@basenml) if(($bline =~ /^\s*[#!\[\/]/) or ($bline =~ /^\s*$/));

    next if(defined $baseid && $line =~ /$baseid/); 
    next if($line =~ /^\s*[#!\[\/]/);
    next if($line =~ /^\s*$/);
    next if($line =~ /runid/);    
    next if($line =~ /model_version/);
    next if($line =~ /logfile/);
    next if($line =~ /username/);
    my($name1,$name2);
    if($line =~ /(.*)=(.*)/){
	$name1 = $1;
    }
    if($bline =~ /(.*)=(.*)/){
        $name2 = $1;
    }
    chomp $line;
    chomp $bline;
    if($name1 eq $name2){
	push(@diffs1,$line);
	push(@diffs2,$bline);
    }elsif($cnt1>$cnt2){
	push(@added,$line);
	$shiftbl=0;
	$cnt1--;
    }elsif($cnt2>$cnt1){
	push(@removed,$bline);
	shift(@basenml);
	$cnt2--;
    }


}

my $fname = $nml;
if($nml=~ /.*\/(.*)$/){
    $fname = $1;
}

if($#diffs1>=0){
    print "\nFAIL namelist compare:  $fname differs\n";
    foreach my $line (@diffs1){
	my $bline = shift @diffs2;
	print "   NEW:        $line\n";
	print "   BASELINE: $bline\n";
    }
    exit -1;
}
   
print "\nPASS namelist compare: $nml and $basenml are the same\n";
if($#added>=0){
    print "   The following variables were added to $nml\n";
    foreach(@added){
	print "    $_\n";
    }
}

if($#removed>=0){
    print "   The following variables were removed from $nml\n";
    foreach(@removed){
	print "    $_\n";
    }
}

exit 0;
