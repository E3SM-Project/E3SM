#!/usr/bin/env perl

# Test methods of the Stream::Template object.

#########################

#use Test::More qw(no_plan); # use "no_plan" until number of tests is determined
use Test::More tests => 27;

#########################

use strict;

use lib "..";
use Cwd;
use English;
use diagnostics;
BEGIN {use_ok('Streams::Template')};

# check that XML::Lite is being found
BEGIN {use_ok('XML::Lite')};

# create Streams::Template object
my $outfile  = "datm.out";
my $cmpfile  = "datm.streams.txt";
my $template = "datm.template.streams.xml";
my %inputopts;
my $cwd = getcwd;
$inputopts{'printing'} = 1;
$inputopts{'ProgName'}  = "$PROGRAM_NAME";
$inputopts{'ProgDir'}   = $cwd;
$inputopts{'yearfirst'} = 1;
$inputopts{'yearlast'}  = 2;
$inputopts{'filepath'}  = ".";
$inputopts{'domain'}    = "";
$inputopts{'domainpath'}= "";
$inputopts{'filenames'} = "";
$inputopts{'cmdline'}   = "";
$inputopts{'test'}      = 0;
$inputopts{'type'}      = "";
$inputopts{'res'}       = "";
$inputopts{'datasource'}= "CLMNCEP";
$inputopts{'case'}      = "casename";
$inputopts{'csmdata'}   = "csmdata";
my $strms= Streams::Template->new( \%inputopts );
isa_ok($strms, "Streams::Template", "created streams template object");

# Check reading in streams template file

$strms->Read( $template );

## check Write method
$strms->Write( $outfile, 1 );

# Check if files are the same.
# The evaluation is messy because it maybe that the pathname to the program will
# be different from one run of this program to another on a different machine/ directory.
#
my $rval = `diff -wbs $outfile $cmpfile`;
my $eval = "Files $outfile and $cmpfile are identical\n";
like($rval, 
qr/$eval|[0-9]+c[0-9]+[ <\n]+$cwd\/02.t -help[ \n->]+\S+\/02.t -help[ \n]*[0-9]+c[0-9]+[ <\n]+$cwd\/02.t[ \n->]+\S+\/02.t/, 'check writing output file');

# Check ability to get domainpath, domainfile, datapath and datafiles 
# when set 
$inputopts{'filepath'}   = "%d/filepath";
$inputopts{'domain'}     = "domain";
$inputopts{'domainpath'} = "%p";
$inputopts{'filenames'}  = "%c.%ym";
$strms= Streams::Template->new( \%inputopts );

$inputopts{'filepath'}   = $inputopts{'csmdata'} . "/filepath";
$inputopts{'domainpath'} = $inputopts{'filepath'};

$strms->Read( $template );
foreach my $lastmonth ( undef, 1 ) {
   my @datafiles;
   if ( $lastmonth ) {
      my $year = $inputopts{'yearfirst'}-1;
      my $month = 12;
      push( @datafiles, $inputopts{'filepath'} . "/" . $inputopts{'case'} . 
                              ".000$year-$month" );
   }
   for( my $year = $inputopts{'yearfirst'}; $year <= $inputopts{'yearlast'}; $year++ ) {
      my @months = ( "01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"
   );
      foreach my $month ( @months ) {
         push( @datafiles, $inputopts{'filepath'} . "/" . $inputopts{'case'} . 
                           ".000$year-$month" );
      }
   }
   my $expandEnv  = 0;
   my $datapath   = $strms->GetDataFilepath(  "data",   $expandEnv );
   my $domainpath = $strms->GetDataFilepath(  "domain", $expandEnv );
   my @data       = $strms->GetDataFilenames( "data",   $expandEnv, $lastmonth );
   my @domain     = $strms->GetDataFilenames( "domain", $expandEnv );

   is($datapath,   $inputopts{'filepath'},   'check getting data   filepath');
   is($domainpath, $inputopts{'domainpath'}, 'check getting domain filepath');
   is_deeply( \@datafiles, \@data, "Make sure data files are as expected" );
   my @domfile = ( $inputopts{'domainpath'} . "/" . $inputopts{'domain'} );
   is_deeply( \@domain, \@domfile, "Make sure domain file is as expected" );
}

# Check that testing of file existance works
$inputopts{'filepath'}   = ".";
$inputopts{'domain'}     = "$outfile";
$inputopts{'domainpath'} = ".";
$inputopts{'filenames'}  = "$cmpfile";
$strms= Streams::Template->new( \%inputopts );
$strms->Read( $template );
$strms->TestFilesExist( "data" );
$strms->TestFilesExist( "domain" );

#
# Check that tInterpAlgo/offset works correctly
#
$inputopts{'datasource'}= "CLMNCEP.Solar";
$strms= Streams::Template->new( \%inputopts );
$strms->Read( $template );
$strms->Write( $outfile );
my $fh = IO::File->new($outfile, '<' ) or die "ERROR: can't open $outfile\n";
my $found    = 0;
my $offfound = 0;
while( $_ = <$fh> ) {
  if ( $_ =~ /tInterpAlgo/ ) {
     $_ = <$fh>;
     $found = 1;
     like($_, qr/^\s*coszen\s*$/,  'Check that tInterpAlgo is correct');
     $_ = <$fh>;
  }
  if ( $_ =~ /offset/ ) {
     $_ = <$fh>;
     $offfound = 1;
     like($_, qr/^\s*-21600\s*$/,  'Check that offset is correct');
     $_ = <$fh>;
  }
}
is($found, 1,   "Check that found tInterpAlgo in $outfile");
is($offfound, 1,"Check that found offset in $outfile");

# Check that fails appropriately if file does not exist
$inputopts{'filepath'}  = ".";
$inputopts{'filenames'} = "xxx";
$strms= Streams::Template->new( \%inputopts );
$strms->Read( $template );
eval{ $strms->TestFilesExist( "data" ); };
like( $@, qr/local file .+ does NOT exist -- aborting/, 'Test for non-existant file');

# check that read of an invalid file works ok
eval { $strms->Read( "xxx" ); };
like( $@, qr/Trouble opening or reading xxx/, 'Read of invalid file');

# check that if try to do something before the Read -- that it will die
$strms= Streams::Template->new( \%inputopts );
eval { $strms->TestFilesExist( "data" ); };
like( $@, qr/a template has NOT been read in yet -- abort/, 'TestFilesExist before a read' );
eval { $strms->Write( "data" ); };
like( $@, qr/a template has NOT been read in yet -- abort/, 'Write before a read' );
eval { $strms->GetDataFilenames( "data" ); };
like( $@, qr/a template has NOT been read in yet -- abort/, 'GetDataFilenames before a read' );
eval { $strms->GetDataFilepath( "data" ); };
like( $@, qr/a template has NOT been read in yet -- abort/, 'GetDataFilepath before a read' );

# Check that if create a streams file and then use it as an input template -- 
# results are nearly the same (other than template name and command line)
$inputopts{'filepath'}   = "%dp";
$inputopts{'domain'}     = "domainfile";
$inputopts{'domainpath'} = ".";
$inputopts{'filenames'}  = "";
$inputopts{'datasource'} = "CAMHIST";
$strms= Streams::Template->new( \%inputopts );
$strms->Read(  $template );
$strms->Write( $outfile  );
my $outfile2 = "$outfile.tmp";
$strms= Streams::Template->new( \%inputopts );
$strms->Read(  $outfile );
$strms->Write( $outfile2  );
$rval = `diff -wbs $outfile $outfile2`;
$eval = "Files $outfile and $outfile2 are identical\n";
like($rval, qr/$eval/, 'Check using output stream files as input template');

#Check that if do not have the right stuff in the inputopts hash that dies appropriately
$inputopts{'domainpath'}= undef;
eval { $strms = Streams::Template->new( \%inputopts ); };
like( $@, qr/Required input variable domainpath was not found/, 'Test abort if inputopts hash does not include domainpath' );

# Check that if need something set and it is not will abort appriopriately
$inputopts{'domainpath'} = "";
$inputopts{'case'}       = "";
$inputopts{'filenames'}  = "%c.%ym";
$inputopts{'datasource'} = "CAMHIST";
$strms= Streams::Template->new( \%inputopts );
$strms->Read(  $template );
eval { $strms->GetDataFilenames( "data" ); };
like( $@, qr/case was NOT defined on command line and needs to be set/, 'Check will abort if need case set' );

# Check that if need something set and it is not will abort appriopriately
$inputopts{'domainpath'} = "";
$inputopts{'case'}       = "case";
$inputopts{'yearfirst'}  = -1;
$inputopts{'yearlast'}   = -1;
$inputopts{'filenames'}  = "%c.%ym";
$inputopts{'datasource'} = "CAMHIST";
$strms= Streams::Template->new( \%inputopts );
$strms->Read(  $template );
eval { $strms->GetDataFilenames( "data" ); };
like( $@, qr/yearfirst and yearlast  was NOT defined on command line and needs to be set/,
 'Check will abort if need year_first/last set' );

# Check for ymd generic CPLHIST case read/write
$cmpfile  = "datm.ymd.streams.txt";
$inputopts{'datasource'} = "CPLHIST";
$inputopts{'domain'}     = "domainfile.nc";
$inputopts{'filenames'}  = "";
$inputopts{'yearfirst'}  = 1;
$inputopts{'yearlast'}   = 1;
$strms= Streams::Template->new( \%inputopts );
$strms->Read(  $template );
$strms->Write( $outfile );

# Check if files are the same.
# The evaluation is messy because it maybe that the pathname to the program will
# be different from one run of this program to another on a different machine/ directory.
#
$rval = `diff -wbs $outfile $cmpfile`;
$eval = "Files $outfile and $cmpfile are identical\n";
like($rval, 
qr/$eval|[0-9]+c[0-9]+[ <\n]+$cwd\/02.t -help[ \n->]+\S+\/02.t -help[ \n]*[0-9]+c[0-9]+[
<\n]+$cwd\/02.t[ \n->]+\S+\/02.t/, 'check writing CPLHIST ymd output file');

#
# Cleanup
#
system( "/bin/rm -f $outfile" );
system( "/bin/rm -f $outfile2" );

print "\nSuccessfully ran all tests\n";
