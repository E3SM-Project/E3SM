#!/usr/bin/perl
use strict;

my $rundir = shift;

opendir(F,$rundir);
my @decompfiles = grep(/^piodecomp/,readdir(F));
closedir(F);
my $rmfile=0;
for(my $i=0; $i< $#decompfiles; $i++){
    my $file  = $decompfiles[$i];
    my $fsize = -s $file;
    for(my $j=$i+1;$j<$#decompfiles;$j++){
	my $nfile = $decompfiles[$j];
	my $f2size = -s $nfile;
	if($fsize == $f2size){
	    open(F1,$file);
	    my @file1 = <F1>;
	    open(F2,$nfile);
	    my @file2 = <F2>;
	    foreach my $line (@file1){
		my $nline = shift (@file2);
		if($line =~ /Obtained/){
		    print "Files $file and $nfile are the same\n";
		    $rmfile=1;
		}
		next if($line == $nline);
		last;
	    }
	    close(F1);
	    close(F2);
	    unlink($nfile) if ($rmfile==1);
	}
    }
}
opendir(F,$rundir);
my @decompfiles = grep(/^piodecomp/,readdir(F));
closedir(F);
for(my $i=0; $i<= $#decompfiles; $i++){
    my $file = $decompfiles[$i];
    open(F1,$file);
    my @file1 = <F1>;
    close(F1);
    open(F1,">$file");
    foreach(@file1){
	if(/\[(.*)\]/){
	    my $decode = `addr2line -e ../bld/cesm.exe $1`;
	    print F1 "$decode\n";
	    print  "$decode\n";
	}else{
	    print F1 $_;
	}
	    
    }
    close(F1);

}
    
