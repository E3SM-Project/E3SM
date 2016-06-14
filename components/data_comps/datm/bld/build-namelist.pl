#!/usr/bin/env perl
#-----------------------------------------------------------------------------------------------
#
# build-namelist
#
# This is the build-namelist script for the CIME datm (Data Atmosphere Model).
#--------------------------------------------------------------------------------------------

use strict;
use Cwd qw(getcwd abs_path);
use English;
use Getopt::Long;
use IO::File;
#-----------------------------------------------------------------------------------------------

####################################
# Process command-line options.
####################################

sub usage {
    die <<EOF;
SYNOPSIS
     build-namelist [options]

     Build the dice namelist
OPTIONS
     -help [or -h]            	   Print usage to STDOUT.
     -infile "filepath"       	   Specify a file or list of files (comma delimited)
                              	   containing namelists to read values from.
     -print "level"           	   Print level for debugging:
                              	      0 = silent
                              	      1 = regular
                              	      2 = verbose
     -caseroot                     caseroot directory variable
     -cimeroot                     cimeroot directory variable
     -inst_string                  INST_STRING variable
     -user_xml_dir "directory"     Directory of where to look for user versions of
                                      namelist XML files (usually your SourceMods/src.dice directory)
                                      (such as namelist_definition_dice.xml, or namelist_defaults_dice.xml)
     -test                    	   Enable checking that input datasets exist on local filesystem.
     -debug                        Run script using default values for env variable settings
                                   for debugging.


Note: The precedence for setting the values of namelist variables is (highest to lowest):
      1. values read from the file specified by -infile,
      2. values from the namelist defaults file.

EOF
}

my %opts = ( help        => 0,
	     silent      => 0,
             debug       => 0,
	     test        => 0,
	     caseroot    => undef,
	     cimeroot    => undef,
	     inst_string => undef,
             user_xml_dir=> undef,
	    );

GetOptions(
    "h|help"         => \$opts{'help'},
    "infile=s"       => \$opts{'infile'},
    "print=i"        => \$opts{'print'},
    "debug"          => \$opts{'debug'},
    "test"           => \$opts{'test'},
    "caseroot=s"     => \$opts{'caseroot'},
    "cimeroot=s"     => \$opts{'cimeroot'},
    "inst_string=s"  => \$opts{'inst_string'},
    "user_xml_dir=s" => \$opts{'user_xml_dir'},
)  or usage();

# Give usage message.
usage() if $opts{'help'};

# Check for unparsed arguments
if (@ARGV) {
    print "ERROR: unrecognized arguments: @ARGV\n";
    usage();
}

# Define print levels:
#   0 - only issue fatal error messages
#   1 - only informs what files are created (default)
#   2 - verbose
my $print = $opts{'print'};

# Set caseroot for Debug mode
if ( $opts{'debug'} && ! defined($opts{'caseroot'}) ) {
    my $cwd = getcwd();
    $opts{'caseroot'} = "$cwd/unit_testers";
    if ( $print > 1 ) { print "caseroot = $opts{'caseroot'}\n"; }
}

# user_xml_dir
my $opt = 'user_xml_dir';
if (defined $opts{$opt}) {
   my $dir = $opts{$opt};
   if ( ! -d "$dir" ) {
       die "** $opt: $dir does NOT exist";
   }
   my @files = glob("$dir/*.xml");
   if ( $#files == -1 && $print > 0 ) {
       print "** Warning NO XML files exist in $opt directory $dir";
   }
}

####################################
# Create required objects
####################################

my $caseroot    = $opts{'caseroot'};
my $cimeroot    = $opts{'cimeroot'};
my $inst_string = $opts{'inst_string'};

my $confdir = "$caseroot/Buildconf/datmconf";

if ( $opts{'debug'} ) {
   my $cmd = "mkdir -p $confdir";
   print "Execute: $cmd\n";
   system( "$cmd" );
   chdir( "$confdir" );
}

# Add $cfgdir/perl5lib to the list of paths that Perl searches for modules
my @dirs = ( "$cimeroot/utils/perl5lib");
unshift @INC, @dirs;

require Build::NamelistDefinition;
require Build::NamelistDefaults;
require Build::Namelist;
require Streams::BuildNamelistUtils;

my ($definition, $defaults, $nl) = BuildNamelistUtils::create_namelist_objects($cimeroot, $caseroot, $confdir, $opts{'user_xml_dir'}, $opts{'infile'});

####################################
# Required xml variables           #
####################################

my %xmlvars = ();
SetupTools::getxmlvars(${caseroot},\%xmlvars);
foreach my $attr (keys %xmlvars) {
    $xmlvars{$attr} = SetupTools::expand_xml_var($xmlvars{$attr}, \%xmlvars);
}
foreach my $var ( "DIN_LOC_ROOT", "ATM_DOMAIN_FILE", "ATM_DOMAIN_PATH",
		  "DATM_MODE", "DATM_PRESAERO", "DATM_TOPO", "DATM_CO2_TSERIES",
                  "ATM_GRID", "GRID" ) {
    if ( $print > 1 ) { print "$var = $xmlvars{$var}\n"; }
    if ( ! defined($xmlvars{$var})  || $xmlvars{$var} =~ /^(UNSET|)$/ ) {
	die "** $var is NOT set  ** \n"
    }
}

my $DIN_LOC_ROOT    = $xmlvars{'DIN_LOC_ROOT'};
my $ATM_DOMAIN_FILE = $xmlvars{'ATM_DOMAIN_FILE'};
my $ATM_DOMAIN_PATH = $xmlvars{'ATM_DOMAIN_PATH'};
my $DATM_MODE       = $xmlvars{'DATM_MODE'};
my $DATM_PRESAERO   = $xmlvars{'DATM_PRESAERO'};
my $DATM_TOPO       = $xmlvars{'DATM_TOPO'};
my $DATM_CO2_TSERIES= $xmlvars{'DATM_CO2_TSERIES'};
my $ATM_GRID        = $xmlvars{'ATM_GRID'};
my $GRID            = $xmlvars{'GRID'};
my $CLM_USRDAT_NAME = $xmlvars{'CLM_USRDAT_NAME'};

# Validate some of the env values.
if ( $DATM_MODE =~ /CLM/ && $DATM_PRESAERO eq "none" ) {
   die "A DATM_MODE for CLM is incompatible with DATM_PRESAERO=none\n";
}
if ( $DATM_MODE =~ /CLM/ && $DATM_TOPO eq "none" ) {
   die "A DATM_MODE for CLM is incompatible with DATM_TOPO=none\n";
}

if ( $GRID eq "CLM_USRDAT" && $CLM_USRDAT_NAME =~ /^(UNSET|)$/ ) {
    die "** GRID=CLM_USRDAT and CLM_USRDAT_NAME is NOT set  **";
}

(-d $DIN_LOC_ROOT)  or mkdir $DIN_LOC_ROOT;
if ($print>=2) { print "Inputdata root directory: $DIN_LOC_ROOT\n"; }

if ($opts{'test'}) {
    (-d $DIN_LOC_ROOT)  or  die <<"EOF";
** Inputdata root is not a directory: \"$DIN_LOC_ROOT\" **
EOF
}
my $test_files = $opts{'test'};

my $var   = "DATM_MODE";
my $group = $definition->get_group_name($var);
$nl->set_variable_value( $group, $var, "\'$xmlvars{$var}\'" );
my $var   = "DATM_PRESAERO";
my $group = $definition->get_group_name($var);
$nl->set_variable_value( $group, $var, "\'$xmlvars{$var}\'" );

####################################
# Streams file(s)                  #
####################################

# Get defaults for data manipulation options associated with
# each stream (mapping, filling, time-interp etc.)

if ($print>=1) { print "  datm mode is $DATM_MODE \n"; }
if ($print>=1) { print "  datm presaero mode is $DATM_PRESAERO \n"; }
if ($print>=1) { print "  datm topo mode is $DATM_TOPO \n"; }
if ($print>=1) { print "  datm model grid is $GRID \n"; }
if ($print>=1) { print "  datm atm   grid is $ATM_GRID \n" };

# Hash for parsing default_namelist_datm.xml
my %default_namelist_opts;
$default_namelist_opts{'grid'}             = $GRID;
$default_namelist_opts{'atm_grid'}         = $ATM_GRID;
$default_namelist_opts{'datm_mode'}        = $DATM_MODE;
$default_namelist_opts{'presaero_mode'}    = $DATM_PRESAERO;
$default_namelist_opts{'datm_co2_tseries'} = $DATM_CO2_TSERIES;
if ($DATM_PRESAERO ne 'none') {
    $default_namelist_opts{'presaero_flag'} = "active";
} else {
    $default_namelist_opts{'presaero_flag'} = "none";
}

# Create streams template file(s) - loop over streams
my $streams = $defaults->get_value( "streamslist", \%default_namelist_opts );
$streams = SetupTools::expand_xml_var( $streams, \%xmlvars );
my @streams = split ",", $streams, -1;
if ($DATM_PRESAERO ne "none") {
    if ($DATM_PRESAERO eq "pt1_pt1") {
	push (@streams, "presaero.$DATM_PRESAERO.$ATM_GRID");
    } else {
	push (@streams, "presaero.$DATM_PRESAERO");
    }
 }
if ($DATM_TOPO ne "none") {
   push (@streams, "topo.$DATM_TOPO");
}
if ($DATM_CO2_TSERIES ne "none") {
    push (@streams, "co2tseries.$DATM_CO2_TSERIES");
}

# anomaly forcing: Check for bias correction streams
my $bias_correct = $nl->get_value( 'bias_correct' );
$bias_correct =~ s/[\'\"]//g;
if ( $bias_correct ) {
    push (@streams, $bias_correct); #from namelist_defaults_datm.xml
}

#Check for anomaly forcing streams
my $anomaly_forcing = $nl->get_value( 'anomaly_forcing' );
$anomaly_forcing =~ s/[\'\"]//g;
my @anomaly_forcing = split ",", $anomaly_forcing, -1;
if ( @anomaly_forcing ) {
    push (@streams, @anomaly_forcing);
}

my %streams_namelists = (
    ostreams  => undef,
    omapalgo  => undef,
    omapmask  => undef,
    otintalgo => undef,
    otaxmode  => undef,
    ofillalgo => undef,
    ofillmask => undef,
    odtlimit  => undef, 
    );

# Create input data list 
my $fh_out = new IO::File;
$fh_out->open(">$caseroot/Buildconf/datm.input_data_list") or
    die "** can't open filepath file: datm.input_data_list\n";

foreach my $stream ( @streams ) {

    # Exit if there is no prescribed aerosol stream
    if ($stream eq "presaero" && $DATM_PRESAERO eq "none" ) {
	next;
    }
    if ($print>=1) { print "  datm stream is $stream$inst_string \n";}

    # Create hash needed to parse namelist_defaults.xml file
    # Hash default_namelist_opts contains the attribute values needed for $defaults->get_value

    $default_namelist_opts{'stream'} = $stream;
    if ( $stream =~ /presaero/ ) {
	$default_namelist_opts{'ispresaerostream'} = "TRUE";
    } else {
	$default_namelist_opts{'ispresaerostream'} = "FALSE";
    }

    # Determine stream txt file
    my $outstream = "datm.streams.txt" . ".$stream" . "$inst_string";
    if (-e "$caseroot/user_$outstream") {
        if ( ! -w "$caseroot/user_$outstream" ) {
           print "Your user streams file is read-only: $caseroot/user_$outstream\n";
           die "Make it writable to continue\n";
        }
	my $command = "cp -p $caseroot/user_$outstream $confdir/$outstream";
	system($command) == 0  or die "system $command failed: $? \n";

    } else {
	BuildNamelistUtils::create_stream_file($fh_out, 'datm', $defaults, \%default_namelist_opts, 
					       \%xmlvars, $stream, $outstream, $opts{'test'});
    }
    BuildNamelistUtils::create_streams_namelists($defaults, \%default_namelist_opts, 
						 \%xmlvars, $stream, $outstream, \%streams_namelists);
}

$fh_out->close;

##########################################################
# namelist group: shr_strdata_nml  (in file datm_atm_in) #
##########################################################

my $datamode   = $defaults->get_value( "datamode", \%default_namelist_opts );
BuildNamelistUtils::add_default($definition, $defaults, $DIN_LOC_ROOT, 
				$nl, 'datamode', "$datamode");  

my $vectors    = $defaults->get_value( "vectors",  \%default_namelist_opts );
BuildNamelistUtils::add_default($definition, $defaults, $DIN_LOC_ROOT, 
				$nl, 'vectors', "$vectors");

BuildNamelistUtils::add_default($definition, $defaults, $DIN_LOC_ROOT, 
				$nl, 'domainfile', "${ATM_DOMAIN_PATH}/${ATM_DOMAIN_FILE}" );

BuildNamelistUtils::add_default($definition, $defaults, $DIN_LOC_ROOT, 
				$nl, 'streams', $streams_namelists{"ostreams"});

BuildNamelistUtils::add_default($definition, $defaults, $DIN_LOC_ROOT, 
				$nl, 'mapalgo', $streams_namelists{"omapalgo"});

BuildNamelistUtils::add_default($definition, $defaults, $DIN_LOC_ROOT, 
				$nl, 'mapmask', $streams_namelists{"omapmask"});

BuildNamelistUtils::add_default($definition, $defaults, $DIN_LOC_ROOT, 
				$nl, 'tintalgo', $streams_namelists{"otintalgo"});

BuildNamelistUtils::add_default($definition, $defaults, $DIN_LOC_ROOT, 
				$nl, 'taxmode', $streams_namelists{"otaxmode"});

BuildNamelistUtils::add_default($definition, $defaults, $DIN_LOC_ROOT, 
				$nl, 'fillalgo', $streams_namelists{"ofillalgo"});

BuildNamelistUtils::add_default($definition, $defaults, $DIN_LOC_ROOT, 
				$nl, 'fillmask', $streams_namelists{"ofillmask"});

BuildNamelistUtils::add_default($definition, $defaults, $DIN_LOC_ROOT, 
				$nl, 'dtlimit', $streams_namelists{"odtlimit"});

##########################################################
# namelist group: datm_nml  (in file datm_in)            #
##########################################################

BuildNamelistUtils::add_default($definition, $defaults, $DIN_LOC_ROOT, 
				$nl, 'atm_in' , "datm_atm_in${inst_string}");

BuildNamelistUtils::add_default($definition, $defaults, $DIN_LOC_ROOT, 
				$nl, 'iradsw');

if ($DATM_MODE =~ /^CORE/) {
    my $factorfn = "$DIN_LOC_ROOT/atm/datm7/CORE2/COREv2.correction_factors.T62.121007.nc";
    BuildNamelistUtils::add_default($definition, $defaults, $DIN_LOC_ROOT, 
				    $nl, 'factorfn', "$factorfn" );
}
if ($DATM_PRESAERO eq "none" ) {
    BuildNamelistUtils::add_default($definition, $defaults, $DIN_LOC_ROOT, 
				    $nl, 'presaero', '.false.');
} else {
    BuildNamelistUtils::add_default($definition, $defaults, $DIN_LOC_ROOT, 
				    $nl, 'presaero', '.true.');
}
BuildNamelistUtils::add_default($definition, $defaults, $DIN_LOC_ROOT, 
				$nl, 'decomp',    '1d');

BuildNamelistUtils::add_default($definition, $defaults, $DIN_LOC_ROOT, 
				$nl, 'force_prognostic_true', '.false.');

BuildNamelistUtils::add_default($definition, $defaults, $DIN_LOC_ROOT, 
				$nl, 'restfilm',  'undefined');

BuildNamelistUtils::add_default($definition, $defaults, $DIN_LOC_ROOT, 
				$nl, 'restfils',  'undefined');


##########################################################
# Write output files
##########################################################

# Validate that the entire resultant namelist is valid
$definition->validate($nl);

my $note = "";

# datm_atm_in
my @groups = qw(shr_strdata_nml);
my $outfile = "./datm_atm_in";
$nl->write($outfile, 'groups'=>\@groups, 'note'=>"$note" );
if ($print>=2) { print "Writing datm_dshr namelist to $outfile \n"; }

# datm_in
@groups = qw(datm_nml);
$outfile = "./datm_in";
$nl->write($outfile, 'groups'=>\@groups, 'note'=>"$note" );
if ($print>=2) { print "Writing datm_in namelist to $outfile \n"; }

# atm_modelio
@groups = qw(modelio);
$outfile = "./atm_modelio.nml";
$nl->set_variable_value( "modelio", "logfile", "'atm.log'" );
$nl->write($outfile, 'groups'=>\@groups, 'note'=>"$note" );
if ($print>=2) { print "Writing atm_modelio.nml namelist to $outfile \n"; }

# Test that input files exist locally.
BuildNamelistUtils::check_input_files($definition, 
				      $nl, $DIN_LOC_ROOT, "$caseroot/Buildconf/datm.input_data_list");



