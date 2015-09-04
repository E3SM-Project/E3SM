#!/usr/bin/env perl

#-# =========================================================================================

=head1 wrapClmTest.pl

=head1 Overview

This is a wrapper script that is called from test_driver.sh for either interactive or batch
tests.  It calls the CTOR for the xFail::expectedFail.pm module and also parses the td*.status
file to create a new file with xFails listed.

It takes the following arguments:

   numberOfTests ->  number of tests from test_driver.sh
   statusFile    ->  name of the td.<pid>.status file
   callingScript ->  name of script calling this.  For test_driver.sh it may be one of:
                        1) test_driver.sh-i for interactive tests
                        2) test_driver.sh   for batch tests

=head1 Notes

This script may be run standalone which is useful for testing purposes.
 
=cut

#-# =========================================================================================

use strict;
use Getopt::Long;
use English;
use Cwd;
use Scalar::Util qw(looks_like_number);

my $DEBUG=0;

sub usage {
    die <<EOF;
SYNOPSIS
     wrapClmTests [options]

     Usually called from test_driver.sh.  Scans the td.*.status file and checks for expected test failures.
OPTIONS
     -help [or -h]                        Print usage to STDOUT.                               
     -numberOfTests "numberOfTests"
     -statusFile    "statusFile"          Name of status file
     -callingScript "callingScript"       Name of calling script

EOF
}

#
# Process command-line options.
#
my %opts = ( help        => 0,
             numberOfTests  => undef,
             statusFile     => undef,
             callingScript  => undef,
            );

GetOptions(
    "h|help"           => \$opts{'help'},
    "numberOfTests=s"  => \$opts{'numberOfTests'},
    "statusFile=s"     => \$opts{'statusFile'},
    "callingScript=s"  => \$opts{'callingScript'},
)  or usage();

# Give usage message.
usage() if $opts{'help'};

my $statFoo = undef;
my $nTests = undef;
my $script= undef;

if (defined($opts{'statusFile'})) {
    $statFoo = $opts{'statusFile'};
}
if (defined($opts{'numberOfTests'})) {
    $nTests = $opts{'numberOfTests'};
}
if (defined($opts{'callingScript'})) {
    $script = $opts{'callingScript'};
}

my ( $self ) = @_;
 
#Figure out where configure directory is and where can use the XML/Lite module from
my $ProgName;
($ProgName = $PROGRAM_NAME) =~ s!(.*)/!!; # name of program
my $ProgDir = $1;                         # name of directory where program lives
 
my $cwd = getcwd();  # current working directory
my $cfgdir;
 
if ($ProgDir) { $cfgdir = $ProgDir; }
else { $cfgdir = $cwd; }
 
#-----------------------------------------------------------------------------------------------
# Add $cfgdir to the list of paths that Perl searches for modules
#-----------------------------------------------------------------------------------------------
my @dirs = ( $cfgdir,
             "$cfgdir/../",
             "$cfgdir/../../../../../../scripts/ccsm_utils/Tools/perl5lib");
unshift @INC, @dirs;
my $result = eval "require expectedFail";
if ( ! defined($result) ) {
   die <<"EOF";
** Cannot find perl module \"xFail/expectedFail.pm\" from directories: @dirs **
EOF
}

#_# ====================================
#_# setup work complete.  Now parse file
#_# ====================================

if ($DEBUG) {
   print (" wrapClmTests.pl:: calling script $script \n");
   print (" wrapClmTests.pl:: number of tests $nTests \n");
   print (" wrapClmTests.pl:: processing $statFoo \n");
}

#_# compGen not used for CLM batch or interactive tests, but we use "compare" as the default in this case
my $compGen="compare";  
my $xFail = xFail::expectedFail->new($script,$compGen,$nTests);

$xFail->parseOutputCLM($statFoo);

exit(0);
