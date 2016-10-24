#!/usr/bin/env perl 

use strict;
use warnings;
use Data::Dumper;
use Test::More qw(no_plan);
use XML::LibXML;
my $scriptsdir = "../../../scripts";

my $banner = "===============================================================================";
print "$banner\nRUNNING UNIT TESTS\n$banner\n";
my $dir = '../';
unshift(@INC, $dir);

require_ok('Testing::CESMTest');
require_ok('Testing::TestLists');
require Testing::CESMTest;
require Testing::TestLists;
my $listobj = Testing::TestLists->new(scriptsdir => $scriptsdir);
my %case;
$case{'compset'} = 'B55TRWCN';
$case{'grid'} = 'f19_g16';
my $msg = $listobj->findTestsForCase(\%case);
like($msg, qr/bluewaters, eos, hopper, intrepid, titan, and yellowstone/);
like($msg, qr/ibm, intel, and pgi/);
like($msg, qr/aux_drv, aux_waccm, prebeta, and prerelease/);

my %case2;
$case2{'compset'} = 'BCN';
$case2{'grid'} = 'ne120_g16';
my $msg2 = $listobj->findTestsForCase(\%case2);

like($msg2, qr/hopper, intrepid, and titan/);
like($msg2, qr/ibm and pgi/);
like($msg2, qr/prerelease/);

my %case3;
$case3{'compset'} = 'B1850BPRP';
$case3{'grid'} = 'f09_g16';
my $msg3 = $listobj->findTestsForCase(\%case3);

like($msg3, qr/bluewaters, eos, goldbach, hopper, intrepid, janus, titan, and yellowstone/);
like($msg3, qr/gnu, ibm, intel, and pgi/);
like($msg3, qr/aux_science, prebeta, and prerelease/);

my %untestedcase;
$untestedcase{'compset'} = 'EGCN';
$untestedcase{'grid'} = 'T31_g37';
my $untestedmsg = $listobj->findTestsForCase(\%untestedcase);
like($untestedmsg, qr/Unfortunately/);
