#!/usr/bin/env perl

# Test methods of the NMLTest::CompFiles object.
# NMLTest::CompFiles is a module to help you build a unit-tester
# for a build-namelist. It allows you to create namelist and output
# files from build-namelist with different options and then compare
# them afterwards.

#########################

#use Test::More qw(no_plan); # use "no_plan" until number of tests is determined
use Test::More tests => 15;

#########################

use strict;

use lib "..";
use Cwd;
use English;
use diagnostics;
BEGIN {use_ok('NMLTest::CompFiles')};

my @files = ( "lnd_in", "run.log" );
my $cwd   = `pwd`;
chomp( $cwd );
my $cfiles  = NMLTest::CompFiles->new( $cwd, @files );
my $options = "default";
my $mode    = "standard";
# Copy *.same files to the standard name
foreach my $file ( @files ) {
   system( "cp $file.same $file" );
}
$cfiles->checkfilesexist( "$options", $mode );
$cfiles->dodiffonfile( $files[0], "$options", $mode );
$cfiles->dodiffonfile( $files[1], "$options", $mode );
# Copy file to standard name and compare (to itself)
$cfiles->copyfiles(    "$options", "$mode" );
$cfiles->comparefiles( "$options", "$mode" );
$cfiles->shownmldiff(  "$options", "$mode" );

# OK, now do comparisons where one of the files is different
system( "cp lnd_in.different $files[0]" );
$cfiles->checkfilesexist( "$options", $mode );
$cfiles->doNOTdodiffonfile( $files[0], "$options", $mode );
$cfiles->dodiffonfile(      $files[1], "$options", $mode );
$cfiles->comparefiles( "$options", "$mode" );
$cfiles->shownmldiff(  "$options", "$mode" );
# Ensure that lnd_in.same and lnd_in.different are indeed different
my $stat = `cmp lnd_in.same lnd_in.different`;
isnt( $stat, 0, "lnd_in.same and lnd_in.different MUST be different");

# Check that things that should FAIL do...
my $baddir = "zztop";
eval{ $cfiles  = NMLTest::CompFiles->new( $baddir, @files ); };
like( $@, qr/ERROR/, "new CompFiles with a bad directory" );

my $mfiles  = NMLTest::CompFiles->new( $cwd, @files );
eval{ $mfiles->comparefiles( $options, $mode ); };
like( $@, qr/ERROR/, "comparefiles too early before setup" );
$mfiles  = NMLTest::CompFiles->new( $cwd, @files );
# get the diffs setup and then compare to an invalid directory
$mfiles->checkfilesexist( "$options", $mode );
eval{ $mfiles->comparefiles( $options, $mode, $baddir ); };
like( $@, qr/ERROR/, "bad directory name to comparefiles" );

# Cleanup
$mfiles->rmfiles( $options, $mode );

print "\nSuccessfully ran all tests\n";


