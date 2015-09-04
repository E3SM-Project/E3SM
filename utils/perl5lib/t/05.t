#!/usr/bin/env perl

# Test methods of the Stream::TemplateGeneric object.

#########################

#use Test::More qw(no_plan); # use "no_plan" until number of tests is determined
use Test::More tests => 28;

#########################

use strict;

use lib "..";
use Cwd;
use English;
#use diagnostics;

# TemplateGeneric object that we are testing...

BEGIN {use_ok('Streams::TemplateGeneric')};

# Also get the Config, and namelistdefaults objects to use for some test settings
BEGIN {use_ok('Build::Config')};

BEGIN {use_ok('Build::NamelistDefaults')};

# check that XML::Lite is being found
BEGIN {use_ok('XML::Lite')};

# Some of the test settings require that we read namelists_defaults file
my $nl_definition_file  = "namelist_definition_cam.xml";
my $nl_defaults_file    = "namelist_defaults_datm.xml";
my $config_file         = "config_cache.xml";
my $inputdata_rootdir   = "/fs/cgd/csm/inputdata";

# Create a namelist defaults object.
my $cfg = Build::Config->new($config_file);
my $defaults = Build::NamelistDefaults->new($nl_defaults_file, $cfg);
my %default_namelist_opts;
# create Streams::TemplateGeneric object
my $outfile  = "datm.out";
my $cmpfile  = "datm.streams.txt";
my $template = "datm.template.streams.xml";
my %inputopts;
my $cwd = getcwd;
my %xmlvar;
$inputopts{'csmdata'}                  = "/fs/cgd/csm/inputdata";
$xmlvar{'DIN_LOC_ROOT'}                = $inputopts{'csmdata'};
$default_namelist_opts{'stream'}       = "CLM_QIAN.Precip";
$default_namelist_opts{'DIN_LOC_ROOT'} = $inputopts{'csmdata'};
$inputopts{'printing'}    = 1;
$inputopts{'ProgName'}    = "$PROGRAM_NAME";
$inputopts{'ProgDir'}     = $cwd;
$inputopts{'cmdline'}     = "";
$inputopts{'yearfirst'}   = &Streams::TemplateGeneric::expandXMLVar( $defaults->get_value( "strm_year_start", \%default_namelist_opts ), \%xmlvar );
$inputopts{'yearlast'}    = &Streams::TemplateGeneric::expandXMLVar( $defaults->get_value( "strm_year_end"  , \%default_namelist_opts ), \%xmlvar );
$inputopts{'filepath'}    = &Streams::TemplateGeneric::expandXMLVar( $defaults->get_value( "strm_datdir"    , \%default_namelist_opts ), \%xmlvar );
$inputopts{'domain'}      = &Streams::TemplateGeneric::expandXMLVar( $defaults->get_value( "strm_domfil"    , \%default_namelist_opts ), \%xmlvar );
$inputopts{'domainpath'}  = &Streams::TemplateGeneric::expandXMLVar( $defaults->get_value( "strm_domdir"    , \%default_namelist_opts ), \%xmlvar );
$inputopts{'datvarnames'} = &Streams::TemplateGeneric::expandXMLVar( $defaults->get_value( "strm_datvar"    , \%default_namelist_opts ), \%xmlvar );
$inputopts{'domvarnames'} = &Streams::TemplateGeneric::expandXMLVar( $defaults->get_value( "strm_domvar"    , \%default_namelist_opts ), \%xmlvar );
$inputopts{'filenames'}   = &Streams::TemplateGeneric::expandXMLVar( $defaults->get_value( "strm_datfil"    , \%default_namelist_opts ), \%xmlvar );
$inputopts{'offset'}      = &Streams::TemplateGeneric::expandXMLVar( $defaults->get_value( "strm_offset"    , \%default_namelist_opts ), \%xmlvar );



my $strms= Streams::TemplateGeneric->new( \%inputopts );
isa_ok($strms, "Streams::TemplateGeneric", "created streams template object");

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
$inputopts{'filepath'}   = "$xmlvar{'DIN_LOC_ROOT'}/filepath";
my $case                 = "b40.1850.track1.1deg.006a";
$xmlvar{'DATM_CPL_CASE'} = $case;
$inputopts{'filenames'}  = "${case}.%ym";

$inputopts{'filepath'}   = $inputopts{'csmdata'} . "/filepath";
$inputopts{'domainpath'} = $inputopts{'filepath'};
$strms= Streams::TemplateGeneric->new( \%inputopts );

$strms->Read( $template );
foreach my $lastmonth ( undef, 1 ) {
   my @datafiles;
   if ( $lastmonth ) {
      my $year = $inputopts{'yearfirst'}-1;
      my $month = 12;
      push( @datafiles, $inputopts{'filepath'} . "/${case}.$year-$month" );
   }
   for( my $year = $inputopts{'yearfirst'}; $year <= $inputopts{'yearlast'}; $year++ ) {
      my @months = ( "01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12" );
      foreach my $month ( @months ) {
         push( @datafiles, $inputopts{'filepath'} . "/${case}.$year-$month" );
      }
   }
   my $datapath   = $strms->GetDataFilepath(  "data",   );
   my $domainpath = $strms->GetDataFilepath(  "domain"  );
   my @data       = $strms->GetDataFilenames( "data",   $lastmonth );
   my @domain     = $strms->GetDataFilenames( "domain"  );

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
$strms= Streams::TemplateGeneric->new( \%inputopts );
$strms->Read( $template );

# Try two CLM1PT type scenarios
$xmlvar{'DATM_CLMNCEP_YR_ALIGN'} = 1;
$xmlvar{'DATM_CLMNCEP_YR_START'} = 1990;
$xmlvar{'DATM_CLMNCEP_YR_END'}   = 1999;
$xmlvar{'CLM_USRDAT_NAME'}       = "13x12pt_f19_alaskaUS";
foreach my $stream ( "1x1_mexicocityMEX", "CLM_USRDAT" ) {
   $default_namelist_opts{'stream'} = "CLM1PT.$stream";
   $inputopts{'yearfirst'}   = &Streams::TemplateGeneric::expandXMLVar( $defaults->get_value( "strm_year_start", \%default_namelist_opts ), \%xmlvar );
   $inputopts{'yearlast'}    = &Streams::TemplateGeneric::expandXMLVar( $defaults->get_value( "strm_year_end"  , \%default_namelist_opts ), \%xmlvar );
   $inputopts{'filepath'}    = &Streams::TemplateGeneric::expandXMLVar( $defaults->get_value( "strm_datdir"    , \%default_namelist_opts ), \%xmlvar );
   $inputopts{'domain'}      = &Streams::TemplateGeneric::expandXMLVar( $defaults->get_value( "strm_domfil"    , \%default_namelist_opts ), \%xmlvar );
   $inputopts{'domainpath'}  = &Streams::TemplateGeneric::expandXMLVar( $defaults->get_value( "strm_domdir"    , \%default_namelist_opts ), \%xmlvar );
   $inputopts{'datvarnames'} = &Streams::TemplateGeneric::expandXMLVar( $defaults->get_value( "strm_datvar"    , \%default_namelist_opts ), \%xmlvar );
   $inputopts{'domvarnames'} = &Streams::TemplateGeneric::expandXMLVar( $defaults->get_value( "strm_domvar"    , \%default_namelist_opts ), \%xmlvar );
   $inputopts{'filenames'}   = &Streams::TemplateGeneric::expandXMLVar( $defaults->get_value( "strm_datfil"    , \%default_namelist_opts ), \%xmlvar );
   $inputopts{'offset'}      = &Streams::TemplateGeneric::expandXMLVar( $defaults->get_value( "strm_offset"    , \%default_namelist_opts ), \%xmlvar );
   if ( $stream eq "CLM_USRDAT" ) {
      is( $inputopts{'yearfirst'},  $xmlvar{'DATM_CLMNCEP_YR_START'},      "yearfirst matches"   );
      is( $inputopts{'yearlast'},   $xmlvar{'DATM_CLMNCEP_YR_END'},        "yearlast matches"    );
      is( $inputopts{'domainpath'}, $xmlvar{'DIN_LOC_ROOT'}."/domainpath", "domainpath matches"  );
      is( $inputopts{'filepath'},   $xmlvar{'DIN_LOC_ROOT'}."/CLM1PT_data/".$xmlvar{'CLM_USRDAT_NAME'}, "domainpath matches"  );
      is( $inputopts{'domain'},     "domain.lnd.1x1pt-$xmlvar{'CLM_USRDAT_NAME'}.nc",                   "domainfile matches"  );
   }

   my $strms= Streams::TemplateGeneric->new( \%inputopts );
   isa_ok($strms, "Streams::TemplateGeneric", "created streams template object for: $stream");
   $strms->Read( $template );
   $strms->Write( $outfile, 1 );
}

# check that if try to do something before the Read -- that it will die
$strms= Streams::TemplateGeneric->new( \%inputopts );
eval { $strms->Write( "data" ); };
like( $@, qr/a template has NOT been read in yet -- abort/, 'Write before a read' );
eval { $strms->GetDataFilenames( "data" ); };
like( $@, qr/a template has NOT been read in yet -- abort/, 'GetDataFilenames before a read' );
eval { $strms->GetDataFilepath( "data" ); };
like( $@, qr/a template has NOT been read in yet -- abort/, 'GetDataFilepath before a read' );

#Check that if do not have the right stuff in the inputopts hash that dies appropriately
$inputopts{'domainpath'}= undef;
eval { $strms = Streams::TemplateGeneric->new( \%inputopts ); };
like( $@, qr/Required input variable domainpath was not found/, 'Test abort if inputopts hash does not include domainpath' );

# Check that if need something set and it is not will abort appriopriately
$default_namelist_opts{'stream'} = "CPLHIST3HrWx.Solar";
$inputopts{'domainpath'} = "";
$inputopts{'yearfirst'}  = -1;
$inputopts{'yearlast'}   = -1;
$inputopts{'filenames'}  = "case.%ym";
$strms= Streams::TemplateGeneric->new( \%inputopts );
$strms->Read(  $template );
eval { $strms->GetDataFilenames( "data" ); };
like( $@, qr/yearfirst and yearlast  was NOT defined on command line and needs to be set/,
 'Check will abort if need year_first/last set' );

# Check for ymd generic CPLHIST case read/write
$cmpfile  = "datm.ymd.streams.txt";
$default_namelist_opts{'stream'} = "CPLHIST3HrWx.Solar";
$inputopts{'yearfirst'}   = &Streams::TemplateGeneric::expandXMLVar( $defaults->get_value( "strm_year_start", \%default_namelist_opts ), \%xmlvar );
$inputopts{'yearlast'}    = &Streams::TemplateGeneric::expandXMLVar( $defaults->get_value( "strm_year_end"  , \%default_namelist_opts ), \%xmlvar );
$inputopts{'filepath'}    = &Streams::TemplateGeneric::expandXMLVar( $defaults->get_value( "strm_datdir"    , \%default_namelist_opts ), \%xmlvar );
$inputopts{'domain'}      = &Streams::TemplateGeneric::expandXMLVar( $defaults->get_value( "strm_domfil"    , \%default_namelist_opts ), \%xmlvar );
$inputopts{'domainpath'}  = &Streams::TemplateGeneric::expandXMLVar( $defaults->get_value( "strm_domdir"    , \%default_namelist_opts ), \%xmlvar );
$inputopts{'datvarnames'} = &Streams::TemplateGeneric::expandXMLVar( $defaults->get_value( "strm_datvar"    , \%default_namelist_opts ), \%xmlvar );
$inputopts{'domvarnames'} = &Streams::TemplateGeneric::expandXMLVar( $defaults->get_value( "strm_domvar"    , \%default_namelist_opts ), \%xmlvar );
$inputopts{'filenames'}   = &Streams::TemplateGeneric::expandXMLVar( $defaults->get_value( "strm_datfil"    , \%default_namelist_opts ), \%xmlvar );
$inputopts{'offset'}      = &Streams::TemplateGeneric::expandXMLVar( $defaults->get_value( "strm_offset"    , \%default_namelist_opts ), \%xmlvar );
$strms= Streams::TemplateGeneric->new( \%inputopts );
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


# Check %glc substitution
# Note that this uses dlnd-like stream information, but for simplicity
# it is still included in files with name datm
$cmpfile = "dlnd.sno.streams.txt";
$default_namelist_opts{'stream'} = "sno.cplhist";
$inputopts{'yearfirst'}   = &Streams::TemplateGeneric::expandXMLVar( $defaults->get_value( "strm_year_start", \%default_namelist_opts ), \%xmlvar );
$inputopts{'yearlast'}    = &Streams::TemplateGeneric::expandXMLVar( $defaults->get_value( "strm_year_end"  , \%default_namelist_opts ), \%xmlvar );
$inputopts{'filepath'}    = &Streams::TemplateGeneric::expandXMLVar( $defaults->get_value( "strm_datdir"    , \%default_namelist_opts ), \%xmlvar );
$inputopts{'domain'}      = &Streams::TemplateGeneric::expandXMLVar( $defaults->get_value( "strm_domfil"    , \%default_namelist_opts ), \%xmlvar );
$inputopts{'domainpath'}  = &Streams::TemplateGeneric::expandXMLVar( $defaults->get_value( "strm_domdir"    , \%default_namelist_opts ), \%xmlvar );
$inputopts{'datvarnames'} = &Streams::TemplateGeneric::expandXMLVar( $defaults->get_value( "strm_datvar"    , \%default_namelist_opts ), \%xmlvar );
$inputopts{'domvarnames'} = &Streams::TemplateGeneric::expandXMLVar( $defaults->get_value( "strm_domvar"    , \%default_namelist_opts ), \%xmlvar );
$inputopts{'filenames'}   = &Streams::TemplateGeneric::expandXMLVar( $defaults->get_value( "strm_datfil"    , \%default_namelist_opts ), \%xmlvar );
$inputopts{'offset'}      = &Streams::TemplateGeneric::expandXMLVar( $defaults->get_value( "strm_offset"    , \%default_namelist_opts ), \%xmlvar );
$inputopts{'glc_nec'} = 10;
$strms= Streams::TemplateGeneric->new( \%inputopts );
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
<\n]+$cwd\/02.t[ \n->]+\S+\/02.t/, 'check %glc substitution');


#
# Cleanup
#
system( "/bin/rm -f $outfile" );

print "\nSuccessfully ran all tests\n";
