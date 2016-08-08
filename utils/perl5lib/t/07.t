#!/usr/bin/env perl

# Test methods of the Testing::TestLists object.
#

#########################

use Test::More qw(no_plan); # use "no_plan" until number of tests is determined
#use Test::More tests => 4;

#########################

use strict;

use lib "..";
use Cwd;
use English;
use diagnostics;
use XML::LibXML;
BEGIN {use_ok('Testing::CESMTest')};
BEGIN {use_ok('Testing::TestLists')};

# Return for untested compset/grid
my %case;
$case{'compset'} = 'BC5';
$case{'grid'} = 'ne30_g16';
my $testlistobj = Testing::TestLists->new(scriptsdir => ".");
my $msg = $testlistobj->findTestsForCase(\%case);
like( $msg, qr/WARNING:: The following compset\/grid combination/, 'Untested combination');
print "Untested Msg: $msg\n";

# Return for well tested compset/grid
$case{'compset'} = 'I';
$case{'grid'} = 'f45_f45';
$testlistobj = Testing::TestLists->new(scriptsdir => ".");
$msg = $testlistobj->findTestsForCase(\%case);
like( $msg, qr/are tested on the following/, 'Tested combination');
print "Tested Msg: $msg\n";

# Failure testing
# directory not given
eval{ $testlistobj = Testing::TestLists->new(); };
like( $@, qr/the scripts dir must be provided to use this module/, 'Directory NOT given');
# bad directory
eval{ $testlistobj = Testing::TestLists->new(".."); };
like( $@, qr/the scripts dir must be provided to use this module/, 'Bad directory');
