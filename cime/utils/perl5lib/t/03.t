#!/usr/bin/env perl

# Test methods of the Decomp::Config object.

#########################

use Test::More tests => 36;
#use Test::More tests => "no_plan";

#########################

use strict;

use lib "..";
use Cwd;
use English;
use diagnostics;
BEGIN {use_ok('Decomp::Config')};

# check that XML::Lite is being found
BEGIN {use_ok('XML::Lite')};

# create Decomp::Config object
my %inputopts;
my $cwd = getcwd;
$inputopts{'printing'}  = 1;
$inputopts{'ProgName'}  = "$PROGRAM_NAME";
$inputopts{'ProgDir'}   = ".";
$inputopts{'cmdline'}   = "";
$inputopts{'res'}       = "gx1v5";
$inputopts{'nproc'}     = 256;
$inputopts{'platform'}  = "XT";
$inputopts{'model'}     = "cice";
my $dcmp= Decomp::Config->new( \%inputopts );
isa_ok($dcmp, "Decomp::Config", "created Decomp config object");

# Check creating new decomp and reading in a XML file on it

my $file = "configInfo.xml";
my %expect= ( 
              "256_gx1v5"   =>{maxblocks=>1, bsize_x=>20, bsize_y=>24, 
                              decomptype=>"cartesian", nlats=>384, nlons=>320},
              "480_tx0.1v2" =>{maxblocks=>1, bsize_x=>180, bsize_y=>100, 
                              decomptype=>"cartesian", nlats=>2400, nlons=>3600},
              "1024_tx0.1v2"=>{maxblocks=>3, bsize_x=>72, bsize_y=>48, 
                              decomptype=>"spacecurve", nlats=>2400, nlons=>3600},
              "2000_tx0.1v2"=>{maxblocks=>1, bsize_x=>72, bsize_y=>60, 
                              decomptype=>"cartesian", nlats=>2400, nlons=>3600},
              "2123_tx0.1v2"=>{maxblocks=>0, bsize_x=>0, bsize_y=>0, 
                              decomptype=>"", nlats=>2400, nlons=>3600},
              "1234_nomatch"=>{maxblocks=>0, bsize_x=>0, bsize_y=>0, 
                              decomptype=>"", nlats=>0, nlons=>0},
            );
foreach my $i( keys(%expect) ) {
   my %decomp = ( maxblocks=>0, bsize_x=>0, bsize_y=>0, decomptype=>"",
               nlats=>0, nlons=>0 );
   my $nproc;
   my $res;
   if ( $i =~ /^([0-9]+)_(gx[13]v[2345]|tx0.1v2|nomatch)$/ ) {
       $nproc = $1;
       $res   = $2;
   }
   my $model = "cice";
   $inputopts{'res'}   = $res;
   $inputopts{'nproc'} = $nproc;
   $inputopts{'model'} = $model;
   print "Test nproc=$nproc, model=$model, res=$res\n";
   $dcmp= Decomp::Config->new( \%inputopts );
   my $matches =  $dcmp->ReadXML( $file, \%decomp );
   if ( $nproc == 2123 ) { 
      is($matches,1,"Check that one matches when expected"); 
   } elsif ( $nproc == 1234 ) { 
      is($matches,0,"Check that no matches when expected"); 
   } elsif ( $nproc == 2000 ) { 
      is($matches,3,"Check that 3 matches when expected"); 
   } else {
      is($matches,2,"Check that two matches when expected"); 
   }
   is_deeply( \%decomp, $expect{$i}, "Make sure data files are as expected" );
   if ( $res eq "gx1v5" ) { $res = "gx1v3"; }
   $model = "pop";
   $inputopts{'res'}   = $res;
   $inputopts{'model'} = $model;
   print "Test nproc=$nproc, model=$model, res=$res\n";
   my $dcmp= Decomp::Config->new( \%inputopts );
   $matches = $dcmp->ReadXML( $file, \%decomp );
   if ( $nproc == 2123 ) { 
      is($matches,1,"Check that one matches when expected"); 
   } elsif ( $nproc == 1234 ) { 
      is($matches,0,"Check that no matches when expected"); 
   } elsif ( $nproc == 2000 ) { 
      is($matches,3,"Check that 3 matches when expected"); 
   } else {
      is($matches,2,"Check that two matches when expected"); 
   }
   if ( $nproc != 2000 ) {
      is_deeply( \%decomp, $expect{$i}, "Make sure data files are as expected" );
   } else {
      $expect{$i}{'decomptype'} = "spacecurve";
      is_deeply( \%decomp, $expect{$i}, "Make sure data files are as expected" );
   }
}

# Check for CAM stuff
my %camexpect= ( 
              "200_fv0.47x0.63" =>{npr_y=>50, npr_z=>4, modcomm_gatscat=>0,
                                   nlats=>384, nlons=>573, nlevs=>26},
              "100_fv0.47x0.63" =>{npr_y=>50, npr_z=>2, modcomm_gatscat=>0,
                                   nlats=>384, nlons=>573, nlevs=>26},
            );
foreach my $i( keys(%camexpect) ) {
   my %decomp = ( npr_y=>0, npr_z=>0, modcomm_gatscat=>0,
                  nlats=>0, nlons=>0, nlevs=>0 );
   my $nproc;
   my $res;
   if ( $i =~ /^([0-9]+)_([^_]+)$/ ) {
       $nproc = $1;
       $res   = $2;
   }
   my $model = "cam";
   $inputopts{'res'}   = $res;
   $inputopts{'nproc'} = $nproc;
   $inputopts{'model'} = $model;
   print "Test nproc=$nproc, model=$model, res=$res\n";
   $dcmp= Decomp::Config->new( \%inputopts );
   my $matches =  $dcmp->ReadXML( $file, \%decomp );
   is($matches,2,"Check that two matches when expected"); 
   is_deeply( \%decomp, $camexpect{$i}, "Make sure data files are as expected" );
}

# Check that doesn't work if input decomp hash is not correct
$inputopts{'nproc'} = 256;
$inputopts{'model'} = "cice";
$inputopts{'res'}   = "tx0.1v2";
my %decomp = ( maxblocks=>0, bsize_x=>0, bsize_y=>0, decomptype=>"",
               nlats=>0, nlons=>0 );
$file = "configInfo.xml";
my %newdecomp = ( stuff=>0, stuff2=>0, thing=>0, thingy=>0, thingy2=>0 );
eval{ $dcmp->ReadXML( $file, \%newdecomp ) };
like( $@, qr/is NOT a valid element for the decomp output hash/, 'Read with bad input hash');
# Check that doesn't work if input decomp is not a hash
$file = "configInfo.xml";
my @newdecomp = keys( %newdecomp );
eval{ $dcmp->ReadXML( $file, \@newdecomp ) };
like( $@, qr/input decomp is not a hash/, 'Read with input hash -- not a hash');
# Check opening file that does not exist
$file = "This_file_does_not_exist_ZZTop.txt";
eval{ $dcmp->ReadXML( $file, \%decomp ) };
like( $@, qr/Trouble opening or reading $file/, 'Read of invalid file');
# create Decomp::Config object without right input elements
my %newinputopts;
$newinputopts{'printing'}  = 1;
$newinputopts{'ProgName'}  = $PROGRAM_NAME;
eval{ Decomp::Config->new( \%newinputopts ); };
like( $@, qr/Required input variable .+ was not found/, 'Create new config without right hash elements');

# create Decomp::Config object without giving it an input hash
my @newinputopts = keys( %inputopts );
eval{ Decomp::Config->new( \@newinputopts ); };
like( $@, qr/ERROR:: input opts is not a hash!/, 'Create new config without input hash');


print "\nSuccessfully ran all tests\n";
