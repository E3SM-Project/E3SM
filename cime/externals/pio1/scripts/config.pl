#!/usr/bin/perl
use strict;
use Getopt::Long;
require Utils;

my $host;
my $compiler;
my $build;
my $result = GetOptions("host=s"=>\$host,
			                 "compiler=s"=>\$compiler,
                                         "build=s"=>\$build);


$host = Utils->host() unless(defined $host);
print "host = $host\n";
Utils->loadmodules("$host");

$build="all" unless(defined($build));

my ($scratch,$netcdf,$pnetcdf,$mpi,$cc,$fc,$filesystem);

($scratch,$netcdf,$pnetcdf,$cc,$fc,$filesystem) = Utils->hostmods($host,$mpi);

print "$scratch\n";
print "$cc\n";
print "$filesystem\n";
print "$fc\n";

# override pnetcdf with my own
#$pnetcdf =  "/glade/u/home/jedwards/pnetcdf/svn1544/intel/";


my $piosrc = `pwd`;
chomp $piosrc;
$piosrc.="/../";

$ENV{CC}=$cc;
$ENV{FC}=$fc;

my $cmake_opts;

if($build eq "netcdf" or $build eq "all"){
    $cmake_opts .= "  -DNETCDF_DIR=$netcdf ";
}
if($build eq "pnetcdf" or $build eq "all"){
    $cmake_opts .= "  -DPNETCDF_DIR=$pnetcdf ";
}
if(defined($filesystem)){
    $cmake_opts .= " -DPIO_FILESYSTEM_HINTS=$filesystem ";
}

$cmake_opts .= " -DPIO_BUILD_TIMING=ON --debug-trycompile ";


mkdir "$scratch";

chdir "$scratch" or die "Could not make directory $scratch";

system("cmake  $cmake_opts -DCMAKE_VERBOSE_MAKEFILE=1 $piosrc");
system("gmake -j 4");


