#!/usr/bin/env perl

# Test methods of the datm build-namelist script.
# Try to test that all the different options at least work.
# Test that inconsistentcies are appropriately caught.

#########################

use Test::More;
use IO::File;

#########################

use strict;
use Getopt::Long;

sub usage {
    die <<EOF;
SYNOPSIS
     build-namelist_test.pl [options]

     Test the the datm build-namelist 
OPTIONS
     -help [or -h]                 Print usage to STDOUT.                               
     -compare <directory>          Compare namelists for this version to namelists
                                   created by another version.
     -generate                     Leave the namelists in place to do a later compare.

EOF
}
# The env variables that will need to be set
my %xmlenv = ( "DATM_CO2_TSERIES", "RUN_TYPE", "DIN_LOC_ROOT", "ATM_DOMAIN_FILE", "ATM_DOMAIN_PATH", "DATM_MODE", "DATM_PRESAERO", "ATM_GRID", "GRID" );
my $envxmlfile = "env_utrun.xml";   # Filename to write env settings to
#
# Process command-line options.
#
my %opts = ( help     => 0,
             generate => 0,
             compare  => undef,
            );

GetOptions(
    "h|help"     => \$opts{'help'},
    "compare=s"  => \$opts{'compare'},
    "generate"   => \$opts{'generate'},
)  or usage();

# Give usage message.
usage() if $opts{'help'};

# Check for unparsed arguments
if (@ARGV) {
    print "ERROR: unrecognized arguments: @ARGV\n";
    usage();
}

#
# Figure out number of tests that will run
#
my $ntests = 490;
if ( defined($opts{'compare'}) ) {
   $ntests += 50;
}
plan( tests=>$ntests );

# Set testcase to a subdirectory
my $cwd = `pwd`;
chomp( $cwd );
my $CASEROOT = "$cwd/testcase";
print "CASEROOT = $CASEROOT\n";
my $confdir = "$CASEROOT/Buildconf/datmconf";

# Check that SCRIPTSROOT set so can run
if ( $ENV{'SCRIPTSROOT'} eq "" ) {
   my $dir = "../../../../../scripts"; 
   if ( ! -d "$dir" ) {
      die "SCRIPTSROOT NOT set, and scripts directory $dir does not exist\n";
   }
   $ENV{'SCRIPTSROOT'} = $dir;
   print "SCRIPTSROOT = $ENV{'SCRIPTSROOT'}\n";
} else {
   if ( ! -d "$ENV{'SCRIPTSROOT'}" ) {
      die "SCRIPTSROOT dir does not exist\n";
   }
}
#
# 
my $bldnml = "../build-namelist -debug -caseroot $CASEROOT -scriptsroot $ENV{'SCRIPTSROOT'}";

my $tempfile = "temp_file.txt";
if ( -f $tempfile ) {
  system( "/bin/rm $tempfile" );
}
# Set some basic ENV vars
$xmlenv{'RUN_TYPE'}              = "startup";
if ( $ENV{'CSMDATA'} ne "" ) {
   $xmlenv{'DIN_LOC_ROOT'}          = $ENV{'CSMDATA'};
} else {
   $xmlenv{'DIN_LOC_ROOT'}          = "/fs/cgd/csm/inputdata";
}
$xmlenv{'DIN_LOC_ROOT_CLMFORC'}     = "/glade/proj2/cgd/tss";
$xmlenv{'ATM_DOMAIN_FILE'}          = "domain.lnd.T31_gx3v7.090928.nc";
$xmlenv{'ATM_DOMAIN_PATH'}          = "\$DIN_LOC_ROOT/share/domains";
$xmlenv{'ATM_GRID'}                 = "48x96";
$xmlenv{'GRID'}                     = "48x96_g37";
$xmlenv{'DATM_MODE'}                = "CLMQIAN";
$xmlenv{'DATM_PRESAERO'}            = "clim_2000";
$xmlenv{'DATM_CO2_TSERIES'}         = "none";
# Set some ENV vars needed for CLM_QIAN
$xmlenv{'DATM_CLMNCEP_YR_START'} = 1948;
$xmlenv{'DATM_CLMNCEP_YR_END'}   = 2004;
$xmlenv{'DATM_CLMNCEP_YR_ALIGN'} = 2000;

# Set some ENV vars needed for CPLHIST
$xmlenv{'DATM_CPLHIST_CASE'}     = "b40.1850.track1.1deg.006a";
$xmlenv{'DATM_CPLHIST_YR_START'} = 1948;
$xmlenv{'DATM_CPLHIST_YR_END'}   = 2004;
$xmlenv{'DATM_CPLHIST_YR_ALIGN'} = 2000;
foreach my $env ( "DATM_CLMNCEP_YR_START", "DATM_CLMNCEP_YR_END", "DATM_CLMNCEP_YR_ALIGN",
                  "DATM_CPLHIST_YR_START", "DATM_CPLHIST_YR_END", "DATM_CPLHIST_YR_ALIGN",
                  "DATM_CPLHIST_CASE",
                ) {
   print "$env = ".$xmlenv{$env}."\n";
}


# Simple test -- just run build-namelist with -help option
eval{ system( "$bldnml -help > $tempfile 2>&1 " ); };
is( $@, '', "help" );

# Simple test -- just run build-namelist
$xmlenv{'DATM_MODE'} = "CLM_QIAN";
print "DATM_MODE      = $xmlenv{'DATM_MODE'}\n";
mkdir( "$CASEROOT" );
&writeEnv( %xmlenv );
eval{ system( "$bldnml" ); };
ok( ! $?, "plain build-namelist" );
&checkfilesexist( );
my %diff;
&copydefaultfiles( "default", \%diff );

# Compare to baseline
if ( defined($opts{'compare'}) ) {
   &comparefiles( \%diff, "default", "$xmlenv{'DATM_MODE'}", $opts{'compare'} );
}
# Cleanup
system( "/bin/rm -rf $confdir"         );

my %wc;
# Run all the different options 
foreach my $option ( "-print 0", "-print 1", "-print 2", "-test -print 2", 
                     "-namelist \"&datm_exp taxmode='extend'/\"",
                     "-user_xml_dir .",
                     "-user_xml_dir $cwd/user_xml_dir -print 2"  ) {
   &writeEnv( %xmlenv );
   eval{ system( "$bldnml $option > $tempfile 2>&1 " ); };
   ok( ! $?, "$option" );
   $wc{$option} = `wc -l $tempfile | awk '{print \$1}'`;
   chomp( $wc{$option} );
   &checkfilesexist( );
   if ( $option =~ /extend/ ) {
      &doNOTdodiffonfile( \%diff, "datm_atm_in", "default" );
   } else {
      &dodiffonfile( \%diff,    "datm_atm_in", "default" );
   }
   if ( $option =~ /$cwd\/user_xml_dir/ ) {
      &doNOTdodiffonfile( \%diff, "datm_atm_in",                      "default" );
      &doNOTdodiffonfile( \%diff, "datm.streams.txt.CLM_QIAN.Precip", "default" );
      &doNOTdodiffonfile( \%diff, "datm.streams.txt.CLM_QIAN.Solar",  "default" );
      &doNOTdodiffonfile( \%diff, "datm.streams.txt.CLM_QIAN.TPQW",   "default" );
   }
   &comparefiles( \%diff, "default", "$xmlenv{'DATM_MODE'}" );
   system( "/bin/rm -rf $confdir"         );
}
ok( $wc{'-print 0'} < $wc{'-print 1'},     , "Silent printing smaller than normal" );
ok( $wc{'-print 1'} < $wc{'-print 2'},     , "Normal printing smaller than verbose" );
ok( $wc{'-print 2'} < $wc{'-test -print 2'}, "Verbose printing smaller than test/verbose" );
&dodiffonfile( \%diff,    "datm_atm_in", "default" );
# Run with user datm stream file
mkdir( "$CASEROOT" );
&writeEnv( %xmlenv );
my $streamfile = "datm.streams.txt.CLM_QIAN.Solar";
system( "cp user_streams/user_$streamfile $CASEROOT" );
eval{ system( "$bldnml" ); };
ok( ! $?, "user_stream build-namelist" );
&checkfilesexist( );
my %diff;

&comparefiles( \%diff, "default", "$xmlenv{'DATM_MODE'}", $opts{'compare'} );
my $diffstat = `diff -q user_streams/user_$streamfile $CASEROOT/Buildconf/datmconf/$streamfile`;
ok( ! $?, "stream file same as user version" );

# Cleanup
system( "/bin/rm -rf $confdir"         );

# Run different RUN_TYPE (all should be the same as default namelist)
foreach my $run_type ( "startup", "branch", "hybrid" ) {
   $xmlenv{'RUN_TYPE'} = $run_type;
   print "RUN_TYPE = $run_type\n";
   &writeEnv( %xmlenv );
   eval{ system( "$bldnml" ); };
   ok( ! $?, "type = $run_type" );
   &checkfilesexist( );
   &comparefiles( \%diff, "default", "$xmlenv{'DATM_MODE'}" );
   system( "/bin/rm -rf $confdir"         );
}
$xmlenv{'RUN_TYPE'} = "startup";
print "RUN_TYPE = $xmlenv{'RUN_TYPE'}\n";
# Run different presaero
foreach my $presaero ( 
    "1x1_brazil", 
    "1x1_camdenNJ",
    "1x1_tropicAtl", 
    "1x1_asphaltjungleNJ", 
    "1x1_vancouverCAN",
    "1x1_mexicocityMEX", 
    "1x1_urbanc_alpha",
    "1x1_numaIA",   
    "1x1_smallvilleIA",
    "5x5_amazon",   
    "clim_1850.1x1_tropicAtl",
    "clim_2000.1x1_tropicAtl",,
    "trans_1850-2000.1x1_tropicAtl",
    "clim_1850",
    "clim_2000",
    "trans_1850-2000",
    "rcp2.6",
    "rcp4.5",
    "rcp6.0",
    "rcp8.5",
 ) {
   if ( $presaero =~ /^(clim|trans)*([_0-9\.-]*)([0-9]+x[0-9]+_[a-zA-Z0-9_]+)$/ ) {
      $xmlenv{'GRID'}     = $3;
      $xmlenv{'ATM_GRID'} = $3;
      if ( $3 eq "1x1_mexicocityMEX" ) {
         $xmlenv{'DATM_MODE'} = "CLM1PT";
      } else {
         $xmlenv{'DATM_MODE'} = "CLM_QIAN";
      }
      $xmlenv{'DATM_PRESAERO'} = "pt1_pt1";
      &doNOTdodiffonfile( \%diff, "datm_atm_in", "default" );
      &doNOTdodiffonfile( \%diff, "datm.streams.txt.$xmlenv{'GRID'}", "default" );
   } else {
      $xmlenv{'GRID'}          = "1.9x2.5_g16";
      $xmlenv{'ATM_GRID'}      = "1.9x2.5";
      $xmlenv{'DATM_PRESAERO'} = "$presaero";
      if ( $presaero eq "clim_2000" ) {
         &dodiffonfile(      \%diff, "datm_atm_in",                "default" );
         &dodiffonfile(      \%diff, "datm.streams.txt.$presaero", "default" );
      } else {
         &doNOTdodiffonfile( \%diff, "datm_atm_in",                "default" );
         &doNOTdodiffonfile( \%diff, "datm.streams.txt.$presaero", "default" );
      }
   }
   print "ATM_GRID       = $xmlenv{'ATM_GRID'}";
   print "DATM_MODE      = $xmlenv{'DATM_MODE'}";
   print "DATM_PRESAERO  = $xmlenv{'DATM_PRESAERO'}";
   &writeEnv( %xmlenv );
   eval{ system( "$bldnml" ); };
   ok( ! $?, "presaero=$presaero" );
   &shownmldiff( "default", "$xmlenv{'DATM_MODE'}" );
   &checkfilesexist( );
   if ( $xmlenv{'DATM_MODE'} ne "CLM1PT" ) {
      &comparefiles( \%diff, "default", "$xmlenv{'DATM_MODE'}" );
   }
   system( "/bin/rm -rf $confdir"         );
}
$xmlenv{'DATM_MODE'}     = "CLM_QIAN";
$xmlenv{'DATM_PRESAERO'} = "clim_2000";
print "DATM_MODE    = $xmlenv{'DATM_MODE'}\n";
print "DATM_PRESAERO= $xmlenv{'DATM_PRESAERO'}\n";
# Check CLM_USRDAT_NAME
$xmlenv{'GRID'}               = "CLM_USRDAT";
$xmlenv{'ATM_GRID'}           = $xmlenv{'GRID'};
print "GRID           = $xmlenv{'GRID'}\n";
print "ATM_GRID       = $xmlenv{'ATM_GRID'}\n";
$xmlenv{'CLM_USRDAT_NAME'} = "13x12pt_f19_alaskaUSA";
$xmlenv{'ATM_DOMAIN_FILE'} ="domain.lnd.\${CLM_USRDAT_NAME}_navy.nc";
$xmlenv{'ATM_DOMAIN_PATH'} ="\$DIN_LOC_ROOT/share/domains/domain.clm";
foreach my $env ( "CLM_USRDAT_NAME", "ATM_DOMAIN_FILE", "ATM_DOMAIN_PATH" ) {
   print "$env\t\t= ".$xmlenv{$env}."\n";
}
&writeEnv( %xmlenv );
eval{ system( "$bldnml" ); };
ok( ! $?, "CLM_USRDAT_NAME" );
&doNOTdodiffonfile( \%diff, "datm_atm_in", "default" );
&checkfilesexist( );
&comparefiles( \%diff, "default", "$xmlenv{'DATM_MODE'}" );
&shownmldiff( "default", "$xmlenv{'DATM_MODE'}" );
# Turn transient CO2 on
$xmlenv{'DATM_CO2_TSERIES'} = "20tr";
&writeEnv( %xmlenv );
eval{ system( "$bldnml" ); };
ok( ! $?, "CLM_USRDAT_NAME" );
&doNOTdodiffonfile( \%diff, "datm_atm_in", "CO2_tseries" );
&checkfilesexist( );
&comparefiles( \%diff, "default", "$xmlenv{'DATM_MODE'}" );
&shownmldiff( "default", "$xmlenv{'DATM_MODE'}" );
$xmlenv{'DATM_CO2_TSERIES'} = "none";
# Check CLM_USRDAT_NAME for CLM1PT mode
$xmlenv{'GRID'}               = "CLM_USRDAT";
$xmlenv{'ATM_GRID'}           = $xmlenv{'GRID'};
$xmlenv{'DATM_MODE'}     = "CLM1PT";
$xmlenv{'DATM_PRESAERO'} = "pt1_pt1";
print "DATM_MODE     = $xmlenv{'DATM_MODE'}\n";
print "DATM_PRESAERO = $xmlenv{'DATM_PRESAERO'}\n";
$xmlenv{'CLM_USRDAT_NAME'} = "US-UMB";
$xmlenv{'ATM_DOMAIN_FILE'} ="domain.lnd.\${CLM_USRDAT_NAME}_navy.nc";
$xmlenv{'ATM_DOMAIN_PATH'} ="\$DIN_LOC_ROOT/share/domains/domain.clm";
foreach my $env ( "CLM_USRDAT_NAME", "ATM_DOMAIN_FILE", "ATM_DOMAIN_PATH" ) {
   print "$env\t\t= ".$xmlenv{$env}."\n";
}
&writeEnv( %xmlenv );
eval{ system( "$bldnml" ); };
ok( ! $?, "CLM_USRDAT_NAME with CLM1PT" );
&doNOTdodiffonfile( \%diff, "datm_atm_in", "default" );
&checkfilesexist( );
&shownmldiff( "default", "$xmlenv{'DATM_MODE'}" );
# Make sure there isn't any unresolved env variables in the output namelist ($)
my $stat = `fgrep \$ $confdir/datm_atm_in`;
isnt( $?, 0, "check for unresolved env variables" );
system( "/bin/rm -rf $confdir"         );
# Run different DATM_MODE
foreach my $mode ( "CORE2_NYF", "CORE2_IAF", "CPLHIST3HrWx", "CLMCRUNCEP","CLMCRUNCEP_V5" ) {
   $xmlenv{'DATM_MODE'} = $mode;
   print "DATM_MODE       = $xmlenv{'DATM_MODE'}\n";
   &writeEnv( %xmlenv );
   eval{ system( "$bldnml" ); };
   ok( ! $?, "mode=$mode" );
   &shownmldiff( "default", "CLMQIAN" );
   &checkfilesexist( );
   &copydefaultfiles( "default", \%diff );
   if ( defined($opts{'compare'}) ) {
      &comparefiles( \%diff, "default", "$xmlenv{'DATM_MODE'}", $opts{'compare'} );
   }
   system( "/bin/rm -rf $confdir"         );
}
# Run with presaero="none";
$xmlenv{'DATM_MODE'}     = "CORE2_NYF";
$xmlenv{'DATM_PRESAERO'} = "none";
print "DATM_MODE       = $xmlenv{'DATM_MODE'}\n";
print "DATM_PRESAERO   = $xmlenv{'DATM_PRESAERO'}\n";
&writeEnv( %xmlenv );
eval{ system( "$bldnml" ); };
ok( ! $?, "presaero=none" );
&shownmldiff( "default", "$xmlenv{'DATM_MODE'}" );
&checkfilesexist( );
system( "/bin/rm -rf $confdir"         );
$xmlenv{'DATM_MODE'}     = "CLM_QIAN";
$xmlenv{'DATM_PRESAERO'} = "clim_2000";
print "DATM_MODE       = $xmlenv{'DATM_MODE'}\n";
print "DATM_PRESAERO   = $xmlenv{'DATM_PRESAERO'}\n";

#
# Error testing -- verify that conditions that should FAIL -- do
#

# Bad arguments or bad namelist value to build-namelist
foreach my $option ( "-zztop", "-namelist '&datm_exp zztop=24/'", "-user_xml_dir zztop" ) {
  &writeEnv( %xmlenv );
  eval{ system( "$bldnml $option  > $tempfile 2>&1 " ); };
  isnt( $?, 0, "Bad argument to build-namelist: $option" );
  system( "cat $tempfile" );
  system( "/bin/rm -rf $confdir"         );
}

# Bad ENV settings
my %bad_env = (
   CLMQIAN_N_PAERONONE =>{ DATM_MODE    =>"CLM_QIAN",   DATM_PRESAERO=>"none"     },
   BAD_DATM_MODE       =>{ DATM_MODE    =>"zztop"       },
   BAD_PRESAERO        =>{ DATM_PRESAERO=>"zztop"       },
   BAD_CQYR_RANGE      =>{ DATM_CLMNCEP_YR_START=>2004, DATM_CLMNCEP_YR_END=>1948 },
   MISSING_CQYR_START  =>{ DATM_CLMNCEP_YR_START=>undef },
   MISSING_CQYR_END    =>{ DATM_CLMNCEP_YR_END  =>undef },
   MISSING_CQYR_ALIGN  =>{ DATM_CLMNCEP_YR_ALIGN=>undef },
   MISSING_CHYR_START  =>{ DATM_MODE    =>"CPLHIST3HrWx",    DATM_CPLHIST_YR_START=>undef },
   MISSING_CHYR_END    =>{ DATM_MODE    =>"CPLHIST3HrWx",    DATM_CPLHIST_YR_END  =>undef },
   MISSING_CHYR_ALIGN  =>{ DATM_MODE    =>"CPLHIST3HrWx",    DATM_CPLHIST_YR_ALIGN=>undef },
   BAD_CHYR_RANGE      =>{ DATM_MODE    =>"CPLHIST3HrWx", 
                           DATM_CPLHIST_YR_START=>1000, DATM_CPLHIST_YR_END=>500 },
   BAD_DATM_CO2TSERIES =>{ DATM_CO2_TSERIES=>"thing" },
   READONLY_USERSTRM   =>{ DATM_MODE    =>"CLM_QIAN" },
              );
foreach my $test ( keys(%bad_env) ) {
   my $envvarref = $bad_env{$test};
   # Set values
   my %def;
   print "Test: $test\n";
   foreach my $env ( keys(%$envvarref) ) {
      $def{$env} = $xmlenv{$env};
      $xmlenv{$env} = $bad_env{$test}{$env};
      print "$env = ".$xmlenv{$env}."\n";
   }
   &writeEnv( %xmlenv );
   # Readonly user stream file
   if ( $test =~ /READONLY_USERSTRM/ ) {
      my $streamfile = "datm.streams.txt.CLM_QIAN.Solar";
      system( "cp -p user_streams/user_$streamfile.readonly $CASEROOT/user_$streamfile" );
   }
   eval{ system( "$bldnml  > $tempfile 2>&1 " ); };
   isnt( $?, 0, "$test" );
   system( "cat $tempfile" );
   # Set values back
   foreach my $env ( keys(%$envvarref) ) {
      $xmlenv{$env} = $def{$env};
      print "$env back to ".$xmlenv{$env}."\n";
   }
   system( "/bin/rm -rf $confdir"         );
}
print "\nCompleted all error tests...\n\n";

# End and cleanup
print "\nCompleted all tests cleanup and quit...\n\n";
system( "/bin/rm -rf $CASEROOT" );
system( "/bin/rm $tempfile" );
if ( ! $opts{'generate'} ) {
   system( "/bin/rm *.default $CASEROOT/$envxmlfile" );
}

sub writeEnv {
#
# Write the input env hash to an env_*.xml file
#
  my %xmlenv = @_;

  my $fh = new IO::File;
  $fh->open(">$CASEROOT/$envxmlfile") or die "can't open file: $CASEROOT/$envxmlfile\n";
  print $fh <<EOF;
<?xml version="1.0"?>

<config_definition>

EOF
  foreach my $id ( keys(%xmlenv) ) {
     print $fh "<entry id=\"$id\"   value=\"$xmlenv{$id}\"  />\n\n";
  }
  print $fh "\n</config_definition>\n";
  $fh->close();
}

sub checkfilesexist {
#
# Check that files exist
#
   my @files = ( 
                 "atm_modelio.nml",  "../datm.input_data_list",
                 "config_cache.xml", "datm_atm_in",
                 "datm_in" 
               );
   foreach my $file ( glob("$confdir/datm.streams.txt.*") ) {
      $file =~ s|$confdir/||;
      push( @files, $file );
   }
   foreach my $file ( @files ) {
      ok( (-f "$confdir/$file" ), "$file file exists" );
   }
}

sub comparefiles {
#
# Compare the the resultant files to the default versions
#
   my $diffref   = shift;
   my $type      = shift;
   my $comp_mode = shift;
   my $compdir   = shift;

   if ( ! defined($type)    ) {
      $type = "default";
   }
   my $compare = "compare to previous tag";
   if ( ! defined($compdir) ) {
      $compdir = ".";
      $compare = undef;
   }
   if ( ! -d "$compdir" ) {
      die "Compare directory $compdir does NOT exist!\n";
   }
   print "Compare files for $type type DATM_MODE=$comp_mode $compare\n";
   my $diffstat;
   my %diffhas = %$diffref;
   my $same = "file the same as expected";
   my $diff = "file different as expected";
   foreach my $file ( "datm_in", "datm_atm_in", "atm_modelio.nml" ) {
      if ( ! -f "$compdir/${file}.$comp_mode.${type}" ) {
         die "File $compdir/${file}.$comp_mode.${type} does NOT exist!\n";
      }
      if ( ! exists($diffhas{$comp_mode}{$type}{$file}) ) {
         next;
      }
      $diffstat = `diff -q $confdir/${file}     $compdir/${file}.$comp_mode.${type}`;
      if ( $diffhas{$comp_mode}{$type}{$file} ) {
         ok( ! $?, "$file $same for $comp_mode" );
      } else {
         like( $diffstat, "/Files $confdir/$file and $compdir/${file}.$comp_mode.${type} differ/", "$file $diff" );
      }
   }
   foreach my $file ( glob("$confdir/datm.streams.txt.*") ) {
      $file =~ s|$confdir/||;
      if ( ! exists($diffhas{$comp_mode}{$type}{$file}) ) {
         next;
      }
      # Get rid of comment about command line options
      system( "sed 's/build-namelist -debug.*/build-namelist -debug   /' $confdir/$file > $confdir/${file}.tmp" );
      system( "/bin/mv $confdir/${file}.tmp  $confdir/${file}" );
      if ( ! -f "$compdir/${file}.$comp_mode.${type}" ) {
         die "File $compdir/${file}.$comp_mode.${type} does NOT exist!\n";
      }
      $diffstat = `diff -q $confdir/$file $compdir/${file}.$comp_mode.${type}`;
      if ( $diffhas{$comp_mode}{$type}{$file} ) {
         ok( ! $?, "$file $same" );
      } else {
         like( $diffstat, "/Files $confdir/$file and $compdir/${file}.$comp_mode.${type} differ/", "$file $diff" );
      }
   }
}

sub copydefaultfiles {
#
# Copy the namelist files to default names for comparisions
#
   my $type    = shift;
   my $diffref = shift;

   my $mode = $xmlenv{'DATM_MODE'};
   foreach my $file ( "datm_in", "datm_atm_in", "atm_modelio.nml" ) {
      system( "/bin/cp  $confdir/$file     ${file}.${mode}.${type}"     );
      $$diffref{${mode}}{${type}}{$file} = 1;
   }
   foreach my $file ( glob("$confdir/datm.streams.txt.*") ) {
      $file =~ s|$confdir/||;
      system( "/bin/cp  $confdir/$file ${file}.${mode}.${type}" );
      $$diffref{${mode}}{${type}}{$file} = 1;
   }
   print "$type namelists for $mode\n";
   system( "/bin/cat datm_in.${mode}.${type}"     );
   system( "/bin/cat datm_atm_in.${mode}.${type}" );
}

sub shownmldiff {
#
# Show the differences in the namelists
#
   my $type      = shift;
   my $comp_mode = shift;

   foreach my $file ( "datm_in", "datm_atm_in" ) {
      my $file1 = "$confdir/$file";
      if ( ! -f "$file1" ) {
         print "$file1 does NOT exist\n";
         return;
      }
      my $file2 = "${file}.${comp_mode}.${type}";
      if ( ! -f "$file2" ) {
         print "$file2 does NOT exist\n";
         return;
      }
      print "Diff in in $file to $type $comp_mode version";
      system( "diff $file1 $file2" ); 
   }

} 

sub dodiffonfile {
#
# Set it so that it does do a difference on the given input file
#
  my $diffref = shift;
  my $file    = shift;
  my $type    = shift;

  if ( ! defined($type) ) {
     $type = "default";
  }
  my $mode = $xmlenv{'DATM_MODE'};
  if ( exists($$diffref{$mode}{$type}{$file}) ) {
     $$diffref{$mode}{$type}{$file} = 1;
  } 
}

sub doNOTdodiffonfile {
#
# Set it so that it does NOT do a difference on the given input file
#
  my $diffref = shift;
  my $file    = shift;
  my $type    = shift;
  
  if ( ! defined($type) ) {
     $type = "default";
  }
  my $mode = $xmlenv{'DATM_MODE'};
  if ( exists($$diffref{$mode}{$type}{$file}) ) {
     $$diffref{$mode}{$type}{$file} = 0;
  } 
}
