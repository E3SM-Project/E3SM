#-----------------------------------------------------------------------------------------------
#
# build-namelist
#
# This script builds the namelists for CLM
#
# The simplest use of build-namelist is to execute it from the build directory where configure
# was run.  By default it will use the config_cache.xml file that was written by configure to
# determine the build time properties of the executable, and will write the files that contain
# the output namelists in that same directory.  But if multiple runs are to made using the
# same executable, successive invocations of build-namelist will overwrite previously generated
# namelist files.  So generally the best strategy is to invoke build-namelist from the run
# directory and use the -config option to provide the filepath of the config_cache.xml file.
#
#
# Date        Contributor      Modification
# -------------------------------------------------------------------------------------------
# 2009-01-20  Vertenstein      Original version
# 2010-04-27  Kluzek           Add ndep streams capability
# 2011-07-25  Kluzek           Add multiple ensemble's of namelists
# 2012-03-23  Kluzek           Add megan namelist and do checking on it
# 2012-07-01  Kluzek           Add some common CESM namelist options
# 2013-12     Andre            Refactor everything into subroutines
#--------------------------------------------------------------------------------------------

package CLMBuildNamelist;

require 5;

use strict;
#use warnings;
#use diagnostics;

use Cwd qw(getcwd abs_path);
use File::Basename qw(dirname);
use English;
use Getopt::Long;
use IO::File;
use File::Glob ':glob';

#-------------------------------------------------------------------------------
#
# Define a small number of global variables
#
#-------------------------------------------------------------------------------

(my $ProgName = $0) =~ s!(.*)/!!;      # name of this script
$ProgName = "CLM " . "$ProgName";

my $cwd = abs_path(getcwd());          # absolute path of the current working directory

my $verbosity = 1;   # Define print level
my $print_verbose = 2;

# Some regular expressions...
###my $TRUE  = qr/\.true\./i;
###my $FALSE = qr/\.false\./i;
# **N.B.** the use of qr// for precompiling regexps isn't supported until perl 5.005.
my $TRUE  = '\.true\.';
my $FALSE = '\.false\.';

#-------------------------------------------------------------------------------

sub usage {
    die <<EOF;
SYNOPSIS
     build-namelist [options]

     Create the namelist for CLM
OPTIONS
     -bgc "value"             Build CLM with BGC package [ sp | cn | bgc ]
                              (default is sp).
                                CLM Biogeochemistry mode
                                sp    = Satellite Phenology (SP)
                                cn    = Carbon Nitrogen model (CN)
                                        (or CLM45CN if phys=clm4_5, use_cn=true)
                                bgc   = Carbon Nitrogen with methane, nitrification, vertical soil C,
                                        CENTURY decomposition
                                        (or CLM45BGC if phys=clm4_5, use_cn=true, use_vertsoilc=true,
                                         use_century_decomp=true, use_nitrif_denitrif=true, and use_lch4=true)
                                    This toggles on the namelist variables:
                                          use_cn, use_lch4, use_nitrif_denitrif, use_vertsoilc, use_century_decomp

     -bgc_spinup "on|off"     CLM 4.5 Only. For CLM 4.0, spinup is controlled from configure.
                              Turn on given spinup mode for BGC setting of CN
                                  on : Turn on Accelerated Decomposition   (spinup_state = 1)
                                  off : run in normal mode                 (spinup_state = 0)

                              Default is off.

                              Spinup is now a two step procedure. First, run the model
                              with spinup = "on". Then run the model for a while with
                              spinup = "off". The exit spinup step happens automatically
                              on the first timestep when using a restart file from spinup
                              mode.

                              The spinup state is saved to the restart file.
                              If the values match between the model and the restart 
                              file it proceeds as directed. 

                              If the restart file is in spinup mode and the model is in
                              normal mode, then it performs the exit spinup step 
                              and proceeds in normal mode after that. 

                              If the restart file has normal mode and the model is in 
                              spinup, then it enters spinup. This is useful if you change
                              a parameter and want to rapidly re-equilibrate without doing
                              a cold start.

     -[no-]chk_res            Also check [do NOT check] to make sure the resolution and
                              land-mask is valid.
     -clm_demand "list"       List of variables to require on clm namelist besides the usuals.
                              "-clm_demand list" to list valid options.
                              (can include a list member "null" which does nothing)
     -clm_start_type "type"   Start type of simulation
                              (default, cold, arb_ic, startup, continue, or branch)
                              (default=do the default type for this configuration)
                              (cold=always start with arbitrary initial conditions)
                              (arb_ic=start with arbitrary initial conditions if
                               initial conditions don't exist)
                              (startup=ensure that initial conditions are being used)
     -clm_usr_name     "name" Dataset resolution/descriptor for personal datasets.
                              Default: not used
                              Example: 1x1pt_boulderCO_c090722 to describe location,
                                       number of pts, and date files created
     -co2_type "value"        Set CO2 the type of CO2 variation to use.
     -co2_ppmv "value"        Set CO2 concentration to use when co2_type is constant (ppmv).
     -config "filepath"       Read the given CLM configuration cache file.
                              Default: "config_cache.xml".
     -crop                    Toggle for prognostic crop model. (default is off)
                              (can ONLY be turned on when BGC type is CN or BGC)
                              This turns on the namelist variable: use_crop
     -csmdata "dir"           Root directory of CESM input data.
                              Can also be set by using the CSMDATA environment variable.
     -d "directory"           Directory where output namelist file will be written
                              Default: current working directory.
     -drydep                  Produce a drydep_inparm namelist that will go into the
                              "drv_flds_in" file for the driver to pass dry-deposition to the atm.
                              Default: -no-drydep
                              (Note: buildnml.csh copies the file for use by the driver)
     -dynamic_vegetation      Toggle for dynamic vegetation model. (default is off)
                              (can ONLY be turned on when BGC type is 'cn' or 'bgc')
                              This turns on the namelist variable: use_cndv
     -envxml_dir "directory"  Directory name of env_*.xml case files to read in.
                              (if read they allow user_nl_clm and CLM_BLDNML_OPTS to expand
                               variables [for example to use \$DIN_LOC_ROOT])
     -glc_grid "grid"         Glacier model grid and resolution when glacier model,
                              Only used if glc_nec > 0 for determining fglcmask
                              Default:  gland5UM
                              (i.e. gland20, gland10 etcetera)
     -glc_nec <name>          Glacier number of elevation classes [0 | 3 | 5 | 10 | 36]
                              (default is 0) (standard option with land-ice model is 10)
     -glc_smb <value>         Only used if glc_nec > 0
                              If .true., pass surface mass balance info to GLC
                              If .false., pass positive-degree-day info to GLC
                              Default: true
     -help [or -h]            Print usage to STDOUT.
     -ignore_ic_date          Ignore the date on the initial condition files
                              when determining what input initial condition file to use.
     -ignore_ic_year          Ignore just the year part of the date on the initial condition files
                              when determining what input initial condition file to use.
     -infile "filepath"       Specify a file (or list of files) containing namelists to
                              read values from.

                              If used with a CLM build with multiple ensembles (ninst_lnd>1)
                              and the filename entered is a directory to files of the
                              form filepath/filepath and filepath/filepath_\$n where \$n
                              is the ensemble member number. the "filepath/filepath"
                              input namelist file is the master input namelist file
                              that is applied to ALL ensemble members.

                              (by default for CESM this is setup for files of the
                               form \$CASEDIR/user_nl_clm/user_nl_clm_????)
     -inputdata "filepath"    Writes out a list containing pathnames for required input datasets in
                              file specified.
     -irrig "value"           If .true. turn irrigation on with namelist logical irrigate (for CLM4.5 physics)
                              (requires use_crop to be true in the clm configuration)
                              Seek surface datasets with irrigation turned on.  (for CLM4.0 physics)
                              Default: .false.
     -l_ncpl "LND_NCPL"       Number of CLM coupling time-steps in a day.
     -lnd_frac "domainfile"   Land fraction file (the input domain file)
     -mask "landmask"         Type of land-mask (default, navy, gx3v5, gx1v5 etc.)
                              "-mask list" to list valid land masks.
     -namelist "namelist"     Specify namelist settings directly on the commandline by supplying
                              a string containing FORTRAN namelist syntax, e.g.,
                                 -namelist "&clm_inparm dt=1800 /"
     -no-megan                DO NOT PRODUCE a megan_emis_nl namelist that will go into the
                              "drv_flds_in" file for the driver to pass VOCs to the atm.
                              MEGAN (Model of Emissions of Gases and Aerosols from Nature)
                              (Note: buildnml.csh copies the file for use by the driver)
     -[no-]note               Add note to output namelist  [do NOT add note] about the
                              arguments to build-namelist.
     -rcp "value"             Representative concentration pathway (rcp) to use for
                              future scenarios.
                              "-rcp list" to list valid rcp settings.
     -res "resolution"        Specify horizontal grid.  Use nlatxnlon for spectral grids;
                              dlatxdlon for fv grids (dlat and dlon are the grid cell size
                              in degrees for latitude and longitude respectively)
                              "-res list" to list valid resolutions.
     -s                       Turns on silent mode - only fatal messages issued.
     -sim_year "year"         Year to simulate for input datasets
                              (i.e. 1850, 2000, 1850-2000, 1850-2100)
                              "-sim_year list" to list valid simulation years
     -test                    Enable checking that input datasets exist on local filesystem.
     -use_case "case"         Specify a use case which will provide default values.
                              "-use_case list" to list valid use-cases.
     -verbose [or -v]         Turn on verbose echoing of informational messages.
     -version                 Echo the SVN tag name used to check out this CLM distribution.
     -vichydro                Toggle to turn on VIC hydrologic parameterizations (default is off)
                              This turns on the namelist variable: use_vichydro


Note: The precedence for setting the values of namelist variables is (highest to lowest):
      0. namelist values set by specific command-line options, like, -d, -sim_year
             (i.e.  compset choice and CLM_BLDNML_OPTS env_run variable)
     (NOTE: If you try to contradict these settings by methods below, an error will be triggered)
      1. values set on the command-line using the -namelist option,
             (i.e. CLM_NAMELIST_OPTS env_run variable)
      2. values read from the file(s) specified by -infile,
             (i.e.  user_nl_clm files)
      3. datasets from the -clm_usr_name option,
             (i.e.  CLM_USRDAT_NAME env_run variable)
      4. values set from a use-case scenario, e.g., -use_case
             (i.e.  CLM_NML_USE_CASE env_run variable)
      5. values from the namelist defaults file.
EOF
}

#-------------------------------------------------------------------------------

sub process_commandline {
  # Process command-line options and return the hash
  my ($nl_flags) = @_;

  # Save the command line arguments to the script. NOTE: this must be
  # before GetOptions() is called because items are removed from from
  # the array!
  $nl_flags->{'cmdline'} = "@ARGV";

  my %opts = ( config                => "config_cache.xml",
               csmdata               => undef,
               clm_usr_name          => undef,
               co2_type              => undef,
               co2_ppmv              => undef,
               clm_demand            => "null",
               help                  => 0,
               glc_nec               => "default",
               glc_grid              => "default",
               glc_smb               => "default",
               l_ncpl                => undef,
               lnd_frac              => undef,
               dir                   => "$cwd",
               rcp                   => "default",
               sim_year              => "default",
               bgc_spinup            => "default",
               chk_res               => undef,
               note                  => undef,
               drydep                => 0,
               megan                 => 1,
               irrig                 => "default",
               res                   => "default",
               silent                => 0,
               mask                  => "default",
               test                  => 0,
               bgc                   => "default",
               crop                  => 0,
               dynamic_vegetation    => 0,
               envxml_dir            => undef,
               vichydro              => 0,
               maxpft                => "default",
             );

  GetOptions(
             "clm_demand=s"              => \$opts{'clm_demand'},
             "co2_ppmv=f"                => \$opts{'co2_ppmv'},
             "co2_type=s"                => \$opts{'co2_type'},
             "config=s"                  => \$opts{'config'},
             "csmdata=s"                 => \$opts{'csmdata'},
             "clm_usr_name=s"            => \$opts{'clm_usr_name'},
             "envxml_dir=s"              => \$opts{'envxml_dir'},
             "drydep!"                   => \$opts{'drydep'},
             "chk_res!"                  => \$opts{'chk_res'},
             "note!"                     => \$opts{'note'},
             "megan!"                    => \$opts{'megan'},
             "glc_nec=i"                 => \$opts{'glc_nec'},
             "glc_grid=s"                => \$opts{'glc_grid'},
             "glc_smb=s"                 => \$opts{'glc_smb'},
             "irrig=s"                   => \$opts{'irrig'},
             "d:s"                       => \$opts{'dir'},
             "h|help"                    => \$opts{'help'},
             "ignore_ic_date"            => \$opts{'ignore_ic_date'},
             "ignore_ic_year"            => \$opts{'ignore_ic_year'},
             "infile=s"                  => \$opts{'infile'},
             "lnd_frac=s"                => \$opts{'lnd_frac'},
             "l_ncpl=i"                  => \$opts{'l_ncpl'},
             "inputdata=s"               => \$opts{'inputdata'},
             "mask=s"                    => \$opts{'mask'},
             "namelist=s"                => \$opts{'namelist'},
             "res=s"                     => \$opts{'res'},
             "rcp=s"                     => \$opts{'rcp'},
             "s|silent"                  => \$opts{'silent'},
             "sim_year=s"                => \$opts{'sim_year'},
             "bgc_spinup=s"              => \$opts{'bgc_spinup'},
             "clm_start_type=s"          => \$opts{'clm_start_type'},
             "test"                      => \$opts{'test'},
             "use_case=s"                => \$opts{'use_case'},
             "bgc=s"                     => \$opts{'bgc'},
             "crop"                      => \$opts{'crop'},
             "dynamic_vegetation"        => \$opts{'dynamic_vegetation'},
             "vichydro"                  => \$opts{'vichydro'},
             "maxpft=i"                  => \$opts{'maxpft'},
             "v|verbose"                 => \$opts{'verbose'},
             "version"                   => \$opts{'version'},
            )  or usage();

  # Give usage message.
  usage() if $opts{'help'};

  # Check for unparsed arguments
  if (@ARGV) {
    print "ERROR: unrecognized arguments: @ARGV\n";
    usage();
  }
  return %opts;
}

#-------------------------------------------------------------------------------

sub set_print_level {
  # Define print levels:
  # 0 - only issue fatal error messages
  # 1 - only informs what files are created (default)
  # 2 - verbose
  my %opts = %{shift()};
  if ($opts{'silent'})  { $verbosity = 0; }
  if ($opts{'verbose'}) { $verbosity = 2; }
}

#-------------------------------------------------------------------------------

sub check_for_perl_utils {
  my $cfgdir = shift;
  # Make sure we can find required perl modules, definition, and defaults files.
  # Look for them under the directory that contains the configure script.

  # The root directory for the input data files must be specified.

  #The root directory to cesm utils Tools
  my $cesm_tools = abs_path( "$cfgdir/../../../../scripts/ccsm_utils/Tools" );

  # The XML::Lite module is required to parse the XML files.
  (-f "$cesm_tools/perl5lib/XML/Lite.pm")  or
    fatal_error("Cannot find perl module \"XML/Lite.pm\" in directory\n" .
                "\"$cesm_tools/perl5lib\"");

  # The Build::Config module provides utilities to access the configuration information
  # in the config_cache.xml file
  (-f "$cesm_tools/perl5lib/Build/Config.pm")  or
    fatal_error("Cannot find perl module \"Build/Config.pm\" in directory\n" .
                "\"$cesm_tools/perl5lib\"");

  # The Build::NamelistDefinition module provides utilities to validate that the output
  # namelists are consistent with the namelist definition file
  (-f "$cesm_tools/perl5lib/Build/NamelistDefinition.pm")  or
    fatal_error("Cannot find perl module \"Build/NamelistDefinition.pm\" in directory\n" .
                "\"$cesm_tools/perl5lib\"");

  # The Build::NamelistDefaults module provides a utility to obtain default values of namelist
  # variables based on finding a best fit with the attributes specified in the defaults file.
  (-f "$cesm_tools/perl5lib/Build/NamelistDefaults.pm")  or
    fatal_error("Cannot find perl module \"Build/NamelistDefaults.pm\" in directory\n" .
                "\"$cesm_tools//perl5lib\"");

  # The Build::Namelist module provides utilities to parse input namelists, to query and modify
  # namelists, and to write output namelists.
  (-f "$cesm_tools/perl5lib/Build/Namelist.pm")  or
    fatal_error("Cannot find perl module \"Build/Namelist.pm\" in directory\n" .
                "\"$cesm_tools/perl5lib\"");

  #-----------------------------------------------------------------------------
  # Add $cesm_tools/perl5lib to the list of paths that Perl searches for modules
  my @dirs = ( $cfgdir, "$cesm_tools/perl5lib");
  unshift @INC, @dirs;

  # required cesm perl modules
  require XML::Lite;
  require Build::Config;
  require Build::NamelistDefinition;
  require Build::NamelistDefaults;
  require Build::Namelist;
  require Streams::TemplateGeneric;

}

#-------------------------------------------------------------------------------

sub read_configure_definition {
  # Read the configure definition and specific config_cache file for this case
  # configure are the build-time settings for CLM
  my ($cfgdir, $opts) = @_;

  verbose_message("Setting CLM configuration script directory to $cfgdir");

  # Create a configuration object from the default config_definition file
  my $configfile;
  if ( -f $opts->{'config'} ) {
    $configfile = $opts->{'config'};
  } else {
    $configfile = "$cfgdir/config_files/config_definition.xml";
  }

  # Check that configuration cache file exists.
  verbose_message("Using CLM configuration cache file $opts->{'config'}");
  if ( $configfile ne $opts->{'config'} ) {
    fatal_error("Cannot find configuration cache file: \"$opts->{'config'}\"\n");
  }

  my $cfg = Build::Config->new("$configfile");

  return $cfg;
}

#-----------------------------------------------------------------------------------------------

sub read_namelist_definition {
  my ($drvblddir, $opts, $nl_flags) = @_;

  # The namelist definition file contains entries for all namelist
  # variables that can be output by build-namelist.
  my @nl_definition_files = ( "$drvblddir/namelist_files/namelist_definition_drv.xml",
                              "$drvblddir/namelist_files/namelist_definition_modio.xml",
                              "$nl_flags->{'cfgdir'}/namelist_files/namelist_definition_$nl_flags->{'phys'}.xml" );
  foreach my $nl_defin_file  ( @nl_definition_files ) {
    (-f "$nl_defin_file")  or  fatal_error("Cannot find namelist definition file \"$nl_defin_file\"\n");

    verbose_message("Using namelist definition file $nl_defin_file");
  }

  # Create a namelist definition object.  This object provides a
  # method for verifying that the output namelist variables are in the
  # definition file, and are output in the correct namelist groups.
  my $definition = Build::NamelistDefinition->new( shift(@nl_definition_files) );
  foreach my $nl_defin_file ( @nl_definition_files ) {
    $definition->add( "$nl_defin_file" );
  }

  return $definition;
}

#-----------------------------------------------------------------------------------------------

sub read_envxml_case_files {
  # read the contents of the env*.xml files in the case directory
  my ($opts) = @_;
  
  my %envxml = ();
  if ( defined($opts->{'envxml_dir'}) ) {
      (-d $opts->{'envxml_dir'})  or  fatal_error( "envxml_dir is not a directory" );
      my @files = glob( $opts->{'envxml_dir'}."/env_*xml" );
      ($#files >= 0)              or  fatal_error( "there are no env_*xml files in the envxml_dir" );
      foreach my $file (@files) {
          verbose_message( "Open env.xml file: $file" );
          my $xml = XML::Lite->new( "$file" );
          my @e   = $xml->elements_by_name('entry');
          while ( my $e = shift @e ) {
              my %a = $e->get_attributes();
              $envxml{$a{'id'}} = $a{'value'};
          }
      }
      foreach my $attr (keys %envxml) {
          if ( $envxml{$attr} =~ m/\$/ ) {
             $envxml{$attr} = &Streams::TemplateGeneric::expandXMLVar( $envxml{$attr}, \%envxml );
          }
      }
  }
  return( %envxml );
}

#-----------------------------------------------------------------------------------------------

sub read_namelist_defaults {
  my ($drvblddir, $opts, $nl_flags, $cfg) = @_;

  # The namelist defaults file contains default values for all required namelist variables.
  my @nl_defaults_files = ( "$nl_flags->{'cfgdir'}/namelist_files/namelist_defaults_overall.xml",
                            "$nl_flags->{'cfgdir'}/namelist_files/namelist_defaults_$nl_flags->{'phys'}.xml",
                            "$drvblddir/namelist_files/namelist_defaults_drv.xml",
                            "$nl_flags->{'cfgdir'}/namelist_files/namelist_defaults_drydep.xml" );

  # Add the location of the use case defaults files to the options hash
  $opts->{'use_case_dir'} = "$nl_flags->{'cfgdir'}/namelist_files/use_cases";

  if (defined $opts->{'use_case'}) {
    if ( $opts->{'use_case'} ne "list" ) {
      unshift( @nl_defaults_files, "$opts->{'use_case_dir'}/$opts->{'use_case'}.xml" );
    }
  }

  foreach my $nl_defaults_file ( @nl_defaults_files ) {
    (-f "$nl_defaults_file")  or  fatal_error("Cannot find namelist defaults file \"$nl_defaults_file\"\n");

    verbose_message("Using namelist defaults file $nl_defaults_file");
  }

  # Create a namelist defaults object.  This object provides default
  # values for variables contained in the input defaults file.  The
  # configuration object provides attribute values that are relevent
  # for the CLM executable for which the namelist is being produced.
  my $defaults = Build::NamelistDefaults->new( shift( @nl_defaults_files ), $cfg);
  foreach my $nl_defaults_file ( @nl_defaults_files ) {
    $defaults->add( "$nl_defaults_file" );
  }
  return $defaults;
}

#-------------------------------------------------------------------------------

sub check_cesm_inputdata {
  # Check that the CESM inputdata root directory has been specified.  This must be
  # a local or nfs mounted directory.

  my ($opts, $nl_flags) = @_;

  $nl_flags->{'inputdata_rootdir'} = undef;
  if (defined($opts->{'csmdata'})) {
    $nl_flags->{'inputdata_rootdir'} = $opts->{'csmdata'};
  }
  elsif (defined $ENV{'CSMDATA'}) {
    $nl_flags->{'inputdata_rootdir'} = $ENV{'CSMDATA'};
  }
  else {
    fatal_error("CESM inputdata root directory must be specified by either -csmdata\n" .
                "argument or by the CSMDATA environment variable.\n");
  }
  if ( ! defined($ENV{'DIN_LOC_ROOT'}) ) {
    $ENV{'DIN_LOC_ROOT'} = $nl_flags->{'inputdata_rootdir'};
  }

  if ($opts->{'test'}) {
    (-d $nl_flags->{'inputdata_rootdir'})  or  fatal_error("CESM inputdata root is not a directory: \"$nl_flags->{'inputdata_rootdir'}\"\n");
  }

  verbose_message("CESM inputdata root directory: $nl_flags->{'inputdata_rootdir'}");
}

#-------------------------------------------------------------------------------

sub process_namelist_user_input {
  # Process the user input in general by order of precedence.  At each point
  # we'll only add new values to the namelist and not overwrite
  # previously specified specified values which have higher
  # precedence. The one exception to this rule are the specifc command-line
  # options which are done last as if the user contradicts these settings
  # CLM build-namelist will abort with an error.
  #
  # 1. values set on the command-line using the -namelist option,
  #         (i.e. CLM_NAMELIST_OPTS env_run variable)
  # 2. values read from the file(s) specified by -infile,
  #         (i.e.  user_nl_clm files)
  # After the above are done the command line options are processed and they
  # are made sure the user hasn't contradicted any of their settings with
  # anything above. Because of this they are condsidered to have the highest
  # precedence.
  # 0. namelist values set by specific command-line options, like, -d, -sim_year
  #         (i.e.  CLM_BLDNML_OPTS env_run variable)
  # The results of these are needed for the final two user input
  # 3. datasets from the -clm_usr_name option,
  #         (i.e.  CLM_USRDAT_NAME env_run variable)
  # 4. values set from a use-case scenario, e.g., -use_case
  #         (i.e.  CLM_NML_USE_CASE env_run variable)
  #
  # Finally after all the above is done, the defaults are found from the
  # namelist defaults file (outside of this routine).
  #


  my ($opts, $nl_flags, $definition, $defaults, $nl, $cfg, $envxml_ref) = @_;

  # Get the inputs that will be coming from the user...
  process_namelist_commandline_namelist($opts, $definition, $nl, $envxml_ref);
  process_namelist_commandline_infile($opts, $definition, $nl, $envxml_ref);

  # Apply the commandline options and make sure the user didn't change it above
  process_namelist_commandline_options($opts, $nl_flags, $definition, $defaults, $nl, $cfg);

  # The last two process command line arguments for usr_name and use_case
  # They require that process_namelist_commandline_options was called before this
  process_namelist_commandline_clm_usr_name($opts, $nl_flags, $definition, $defaults, $nl, $cfg, $envxml_ref);
  process_namelist_commandline_use_case($opts, $nl_flags, $definition, $defaults, $nl, $cfg, $envxml_ref);

  # Set the start_type by the command line setting for clm_start_type
  process_namelist_commandline_clm_start_type($opts->{'test'}, $nl_flags, $definition, $defaults, $nl);

}

#-------------------------------------------------------------------------------

sub process_namelist_commandline_options {
  # First process the commandline args that provide specific namelist values.
  #
  # First get the command-line specified overall values or their defaults
  # Obtain default values for the following build-namelist input arguments
  # : res, mask, rcp, sim_year, sim_year_range, and bgc_spinup.
  #
  # NOTE: cfg only needs to be passed to functions that work with
  # clm4_0 compile time functionality!

  my ($opts, $nl_flags, $definition, $defaults, $nl, $cfg) = @_;

  setup_cmdl_chk_res($opts, $defaults);
  setup_cmdl_resolution($opts, $nl_flags, $definition, $defaults);
  setup_cmdl_mask($opts, $nl_flags, $definition, $defaults, $nl);
  setup_cmdl_bgc($opts, $nl_flags, $definition, $defaults, $nl, $cfg);
  setup_cmdl_crop($opts, $nl_flags, $definition, $defaults, $nl, $cfg);
  setup_cmdl_maxpft($opts, $nl_flags, $definition, $defaults, $nl, $cfg);
  setup_cmdl_glc_nec($opts, $nl_flags, $definition, $defaults, $nl);
  setup_cmdl_irrigation($opts, $nl_flags, $definition, $defaults, $nl);
  setup_cmdl_rcp($opts, $nl_flags, $definition, $defaults, $nl);
  setup_cmdl_bgc_spinup($opts, $nl_flags, $definition, $defaults, $nl, $cfg);
  setup_cmdl_simulation_year($opts, $nl_flags, $definition, $defaults, $nl);
  setup_cmdl_run_type($opts, $nl_flags, $definition, $defaults, $nl);
  setup_cmdl_dynamic_vegetation($opts, $nl_flags, $definition, $defaults, $nl);
  setup_cmdl_vichydro($opts, $nl_flags, $definition, $defaults, $nl);
}

#-------------------------------------------------------------------------------

sub setup_cmdl_chk_res {
  my ($opts, $defaults) = @_;

  my $var = "chk_res";
  if ( ! defined($opts->{$var}) ) {
    $opts->{$var} = $defaults->get_value($var);
  }
}

sub setup_cmdl_resolution {
  my ($opts, $nl_flags, $definition, $defaults) = @_;

  my $var = "res";
  my $val;

  if ( $opts->{$var} ne "default" ) {
    $val = $opts->{$var};
  } else {
    $val= $defaults->get_value($var);
  }

  $nl_flags->{'res'} = $val;
  verbose_message("CLM atm resolution is $nl_flags->{'res'}");
  $opts->{$var} = $val;
  if ( $opts->{'chk_res'} ) {
    $val = &quote_string( $nl_flags->{'res'} );
    if (  ! $definition->is_valid_value( $var, $val ) ) {
      my @valid_values   = $definition->get_valid_values( $var );
      if ( ! defined($opts->{'clm_usr_name'}) || $nl_flags->{'res'} ne $opts->{'clm_usr_name'} ) {
        fatal_error("$var has a value ($val) that is NOT valid. Valid values are: @valid_values\n");
      }
    }
  }
}

#-------------------------------------------------------------------------------

sub setup_cmdl_mask {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  my $var = "mask";
  my $val;

  if ( $opts->{$var} ne "default" ) {
    $val = $opts->{$var};
  } else {
    my %tmp = ( 'hgrid'=>$nl_flags->{'res'} );
    $val = $defaults->get_value($var, \%tmp );
  }

  $nl_flags->{'mask'} = $val;
  $opts->{'mask'} = $nl_flags->{'mask'};
  if ( $opts->{'chk_res'} ) {
    $val = &quote_string( $val );
    my $group = $definition->get_group_name($var);
    $nl->set_variable_value($group, $var, $val);
    if (  ! $definition->is_valid_value( $var, $val ) ) {
      my @valid_values   = $definition->get_valid_values( $var );
      fatal_error("$var has a value ($val) that is NOT valid. Valid values are: @valid_values\n");
    }
  }
  verbose_message("CLM land mask is $nl_flags->{'mask'}");
}

#-------------------------------------------------------------------------------

sub setup_cmdl_bgc {
  # BGC - alias for group of biogeochemistry related use_XXX namelists

  my ($opts, $nl_flags, $definition, $defaults, $nl, $cfg) = @_;

  my $val;
  my $var = "bgc";

  $val = $opts->{$var};
  $nl_flags->{'bgc_mode'} = $val;
  if ($nl_flags->{'phys'} eq "clm4_0") {
    if ( $nl_flags->{'bgc_mode'} ne "default" ) {
      fatal_error("-bgc option used with clm4_0 physics. -bgc can ONLY be used with clm4_5 physics");
    }
    $nl_flags->{'bgc_mode'} = $cfg->get($var);
  } else {
    my $var = "bgc_mode";
    if ( $nl_flags->{$var} eq "default" ) {
       $nl_flags->{$var} = $defaults->get_value($var);
    }
    my $group = $definition->get_group_name($var);
    $nl->set_variable_value($group, $var, quote_string( $nl_flags->{$var} ) );
    if (  ! $definition->is_valid_value( $var, quote_string( $nl_flags->{$var}) ) ) {
      my @valid_values   = $definition->get_valid_values( $var );
      fatal_error("$var has a value (".$nl_flags->{$var}.") that is NOT valid. Valid values are: @valid_values\n");
    }
    verbose_message("Using $nl_flags->{$var} for bgc.");

    # now set the actual name list variables based on the bgc alias
    my $setting = ".false.";
    if ($nl_flags->{$var} eq "cn") {
      $nl_flags->{'use_cn'} = ".true.";
    } elsif ($nl_flags->{$var} eq "bgc") {
      $nl_flags->{'use_cn'} = ".true.";
      $setting = ".true.";
    } else {
      $nl_flags->{'use_cn'} = ".false.";
    }
    if ( defined($nl->get_value("use_cn")) && ($nl_flags->{'use_cn'} ne $nl->get_value("use_cn")) ) {
      fatal_error("The namelist variable use_cn is inconsistent with the -bgc option");
    }
    # If the variable has already been set use it, if not set to the value defined by the bgc_mode
    my @list  = (  "use_lch4", "use_nitrif_denitrif", "use_vertsoilc", "use_century_decomp" );
    my $ndiff = 0;
    foreach my $var ( @list ) {
       if ( ! defined($nl->get_value($var))  ) {
          $nl_flags->{$var} = $setting;
       } else {
          if ( $nl->get_value($var) ne $setting ) {
             $ndiff += 1;
          }
          $nl_flags->{$var} = $nl->get_value($var);
       }
       $val = $nl_flags->{$var};
       my $group = $definition->get_group_name($var);
       $nl->set_variable_value($group, $var, $val);
       if (  ! $definition->is_valid_value( $var, $val ) ) {
         my @valid_values   = $definition->get_valid_values( $var );
         fatal_error("$var has a value ($val) that is NOT valid. Valid values are: @valid_values\n");
       }
    }
    # If all the variables are different report it as an error
    if ( $ndiff == ($#list + 1) ) {
       fatal_error("You are contradicting the -bgc setting with the namelist variables: @list" );
    }

    # Now set use_cn
    $var = "use_cn";
    $val = $nl_flags->{'use_cn'};
    my $group = $definition->get_group_name($var);
    $nl->set_variable_value($group, $var, $val);
    if (  ! $definition->is_valid_value( $var, $val ) ) {
      my @valid_values   = $definition->get_valid_values( $var );
      fatal_error("$var has a value ($val) that is NOT valid. Valid values are: @valid_values\n");
    }
  }
} # end bgc

#-------------------------------------------------------------------------------

sub setup_cmdl_crop {
  my ($opts, $nl_flags, $definition, $defaults, $nl, $cfg) = @_;

  $nl_flags->{'use_crop'} = ".false.";
  my $val;
  my $var = "crop";
  if ($nl_flags->{'phys'} eq "clm4_0") {
    $nl_flags->{'crop'} = $cfg->get($var);
    if ( $nl_flags->{'crop'} eq "on" ) {
      $nl_flags->{'use_crop'} = ".true.";
    }
  } else {
    $val = $opts->{$var};
    $nl_flags->{'crop'} = $val;
    if ( $nl_flags->{'crop'} eq 1 ) {
      $nl_flags->{'use_crop'} = ".true.";
    }
    if ( defined($nl->get_value("use_crop")) && ($nl_flags->{'use_crop'} ne $nl->get_value("use_crop")) ) {
      fatal_error("Namelist item use_crop contradicts the command-line option -crop, use the command line option");
    }
    if ( ($nl_flags->{'crop'} eq 1 ) && ($nl_flags->{'bgc_mode'} eq "sp") ) {
      fatal_error("** Cannot turn crop mode on mode bgc=sp\n" .
                  "**\n" .
                  "** Set the bgc mode to 'cn' or 'bgc' by the following means from highest to lowest precedence:\n" .
                  "** * by the command-line options -bgc cn\n" .
                  "** * by a default configuration file, specified by -defaults\n");
    }

    $var = "use_crop";
    $val = ".false.";
    if ($nl_flags->{'crop'} eq 1) {
      $val = ".true.";
    }
    my $group = $definition->get_group_name($var);
    $nl->set_variable_value($group, $var, $val);
    if (  ! $definition->is_valid_value( $var, $val ) ) {
      my @valid_values   = $definition->get_valid_values( $var );
      fatal_error("$var has a value ($val) that is NOT valid. Valid values are: @valid_values\n");
    }
  }
}

#-------------------------------------------------------------------------------

sub setup_cmdl_maxpft {
  my ($opts, $nl_flags, $definition, $defaults, $nl, $cfg) = @_;

  my $val;
  my $var = "maxpft";
  if ($nl_flags->{'phys'} eq "clm4_0") {
    $nl_flags->{'maxpft'} = $cfg->get($var);
    # NOTE: maxpatchpft sizes already checked for clm4_0 by configure.
  } else {
    my %maxpatchpft;
    $maxpatchpft{'.true.'}   = 25;
    $maxpatchpft{'.false.'} = 17;
    if ( $opts->{$var} ne "default") {
      $val = $opts->{$var};
    } else {
      $val = $maxpatchpft{$nl_flags->{'use_crop'}};
    }
    $nl_flags->{'maxpft'} = $val;

    if ( ($nl_flags->{'bgc_mode'} ne "sp") && ($nl_flags->{'maxpft'} != $maxpatchpft{$nl_flags->{'use_crop'}}) ) {
      fatal_error("** For CN or BGC mode you MUST set max patch PFT's to $maxpatchpft{$nl_flags->{'use_crop'}}\n" .
                  "**\n" .
                  "** When the crop model is on then it must be set to $maxpatchpft{'crop'} otherwise to $maxpatchpft{'nocrop'}\n" .
                  "** Set the bgc mode, crop and maxpft by the following means from highest to lowest precedence:\n" .
                  "** * by the command-line options -bgc, -crop and -maxpft\n" .
                  "** * by a default configuration file, specified by -defaults\n" .
                  "**\n");
    }
    if ( $nl_flags->{'maxpft'} > $maxpatchpft{$nl_flags->{'use_crop'}} ) {
      fatal_error("** Max patch PFT's can NOT exceed $maxpatchpft{$nl_flags->{'use_crop'}}\n" .
                  "**\n" .
                  "** Set maxpft by the following means from highest to lowest precedence:\n" .
                  "** * by the command-line options -maxpft\n" .
                  "** * by a default configuration file, specified by -defaults\n" .
                  "**\n");
    }
    if ( $nl_flags->{'maxpft'} != $maxpatchpft{$nl_flags->{'use_crop'}} ) {
      warning("running with maxpft NOT equal to $maxpatchpft{$nl_flags->{'use_crop'}} is " .
              "NOT validated / scientifically supported.\n");
    }
    verbose_message("Using $nl_flags->{'maxpft'} for maxpft.");

    $var = "maxpatch_pft";
    $val = $nl_flags->{'maxpft'};
    my $group = $definition->get_group_name($var);
    $nl->set_variable_value($group, $var, $val);
    if (  ! $definition->is_valid_value( $var, $val ) ) {
      my @valid_values   = $definition->get_valid_values( $var );
      fatal_error("$var has a value ($val) that is NOT valid. Valid values are: @valid_values\n");
    }
  }
}

#-------------------------------------------------------------------------------

sub setup_cmdl_glc_nec {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  my $val;
  my $var = "glc_nec";

  if ( $opts->{$var} ne "default" ) {
    $val = $opts->{$var};
  } else {
    $val = $defaults->get_value($var);
  }

  $nl_flags->{'glc_nec'} = $val;
  $opts->{'glc_nec'} = $val;
  my $group = $definition->get_group_name($var);
  $nl->set_variable_value($group, $var, $val);
  if (  ! $definition->is_valid_value( $var, $val ) ) {
    my @valid_values   = $definition->get_valid_values( $var );
    fatal_error("$var has a value ($val) that is NOT valid. Valid values are: @valid_values\n");
  }
  verbose_message("Glacier number of elevation classes is $val");
}

#-------------------------------------------------------------------------------

sub setup_cmdl_irrigation {
  # Must be after setup_cmdl_crop
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  my $var   = "irrig";

  if ( $opts->{$var} eq "default" ) {
    $nl_flags->{$var} = $defaults->get_value($var);
  } else {
    $nl_flags->{$var} = $opts->{$var};
  }
  my $val   = $nl_flags->{$var};
  my $group = $definition->get_group_name($var);
  $nl->set_variable_value($group, $var, $val);
  if (  ! $definition->is_valid_value( $var, $val ) ) {
    my @valid_values   = $definition->get_valid_values( $var );
    fatal_error("$var has a value ($val) that is NOT valid. Valid values are: @valid_values\n");
  }
  verbose_message("Irrigation $val");
  if ($nl_flags->{'phys'} eq "clm4_0") {
    if ( $nl_flags->{'irrig'} =~ /$TRUE/i && $nl_flags->{'use_crop'} eq ".true." ) {
      fatal_error("You've turned on both irrigation and crop.\n" .
                  "Irrigation is only applied to generic crop currently,\n" .
                  "which negates it's practical usage.\n." .
                  "We also have a known problem when both are on " .
                  "(see bug 1326 in the models/lnd/clm/doc/KnownBugs file)\n" .
                  "both irrigation and crop can NOT be on.\n");
    }
  } else {
    if ( $nl_flags->{'irrig'} =~ /$TRUE/i && $nl_flags->{'use_crop'} =~ /$FALSE/i ) {
      fatal_error("The -irrig=.true. option requires -crop");
    }
    if ( defined($nl->get_value("irrigate")) && $nl->get_value("irrigate") ne $nl_flags->{'irrig'} ) {
      fatal_error("The namelist value irrigate contradicts the command line option -irrig");
    }
  }
}

#-------------------------------------------------------------------------------

sub setup_cmdl_rcp {
  # representative concentration pathway
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  my $val;
  my $var = "rcp";
  if ( $opts->{$var} ne "default" ) {
    $val = $opts->{$var};
  } else {
    $val = $defaults->get_value($var);
  }
  $nl_flags->{'rcp'} = $val;
  $opts->{'rcp'} = $nl_flags->{'rcp'};
  my $group = $definition->get_group_name($var);
  $nl->set_variable_value($group, $var, $val);
  if (  ! $definition->is_valid_value( $var, $val ) ) {
    my @valid_values   = $definition->get_valid_values( $var );
    fatal_error("$var has a value ($val) that is NOT valid. Valid values are: @valid_values\n");
  }
  verbose_message("CLM future scenario representative concentration is $nl_flags->{'rcp'}");
}

#-------------------------------------------------------------------------------

sub setup_cmdl_bgc_spinup {
  # CLM 4.0 --> spinup mode controlled from "spinup" in configure
  # CLM 4.5 --> spinup mode controlled from "bgc_spinup" in build-namelist
  my ($opts, $nl_flags, $definition, $defaults, $nl, $cfg) = @_;

  my $val;
  my $var;
  $nl_flags->{'spinup'} = undef;
  $nl_flags->{'bgc_spinup'} = undef;
  if ($nl_flags->{'phys'} eq "clm4_0") {
    $nl_flags->{'spinup'} = $cfg->get('spinup');
    if ($opts->{"bgc_spinup"} ne "default") {
      fatal_error("bgc_spinup can not be controlled from the namelist in CLM 4.0. Try configure using CLM_CONFIG_OPTS (-spinup).");
    }
  } elsif ($nl_flags->{'phys'} eq "clm4_5") {
    $var = "bgc_spinup";
    if ( $opts->{$var} ne "default" ) {
      $val = $opts->{$var};
    } else {
      $val = $defaults->get_value($var);
    }
    $nl_flags->{$var} = $val;
    my $group = $definition->get_group_name($var);
    $nl->set_variable_value($group, $var, quote_string($val) );
    if (  ! $definition->is_valid_value( $var, $val , 'noquotes' => 1) ) {
      my @valid_values   = $definition->get_valid_values( $var );
      fatal_error("$var has an invalid value ($val). Valid values are: @valid_values\n");
    }
    if ( $nl_flags->{'bgc_spinup'} eq "on" && $nl_flags->{'use_cn'} ne ".true.") {
      fatal_error("$var can not be '$nl_flags->{'bgc_spinup'}' if CN is turned off (use_cn=$nl_flags->{'use_cn'}).");
    }
    if ( $nl->get_value("spinup_state") eq 0 && $nl_flags->{'bgc_spinup'} eq "on" ) {
      fatal_error("Namelist spinup_state contradicts the command line option bgc_spinup" );
    }
    if ( $nl->get_value("spinup_state") eq 1 && $nl_flags->{'bgc_spinup'} eq "off" ) {
      fatal_error("Namelist spinup_state contradicts the command line option bgc_spinup" );
    }
  }

  if ($nl_flags->{'phys'} eq "clm4_0") {
    $val = $nl_flags->{'spinup'};
  } else {
    $val = $nl_flags->{'bgc_spinup'};
  }
  verbose_message("CLM CN bgc_spinup mode is $val");
}

#-------------------------------------------------------------------------------

sub setup_cmdl_simulation_year {
  my ($opts, $nl_flags, $definition, $defaults, $nl, $cfg) = @_;

  my $val;
  my $var = "sim_year";
  if ( $opts->{$var} ne "default" ) {
    $val = $opts->{$var};
  } else {
    $val = $defaults->get_value($var);
  }

  $nl_flags->{'sim_year_range'} = $defaults->get_value("sim_year_range");
  $nl_flags->{'sim_year'}       = $val;
  if ( $val =~ /([0-9]+)-([0-9]+)/ ) {
    $nl_flags->{'sim_year'}       = $1;
    $nl_flags->{'sim_year_range'} = $val;
  }
  $val = $nl_flags->{'sim_year'};
  my $group = $definition->get_group_name($var);
  $nl->set_variable_value($group, $var, $val );
  if (  ! $definition->is_valid_value( $var, $val, 'noquotes'=>1 ) ) {
    my @valid_values   = $definition->get_valid_values( $var );
    fatal_error("$var of $val is NOT valid. Valid values are: @valid_values\n");
  }
  $nl->set_variable_value($group, $var, $val );
  verbose_message("CLM sim_year is $nl_flags->{'sim_year'}");

  $var = "sim_year_range";
  $val = $nl_flags->{'sim_year_range'};
  if ( $val ne "constant" ) {
    $opts->{$var}   = $val;
    $group = $definition->get_group_name($var);
    $nl->set_variable_value($group, $var, $val );
    if (  ! $definition->is_valid_value( $var, $val, 'noquotes'=>1 ) ) {
      my @valid_values   = $definition->get_valid_values( $var );
      fatal_error("$var of $val is NOT valid. Valid values are: @valid_values\n");
    }
    $val = "'".$defaults->get_value($var)."'";
    $nl->set_variable_value($group, $var, $val );
    verbose_message("CLM sim_year_range is $nl_flags->{'sim_year_range'}");
  }
}

#-------------------------------------------------------------------------------

sub setup_cmdl_run_type {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  my $val;
  my $var = "clm_start_type";
  if (defined $opts->{$var}) {
    if ($opts->{$var} eq "default" ) {
      add_default($opts->{'test'}, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var );
    } else {
      my $group = $definition->get_group_name($var);
      $nl->set_variable_value($group, $var, quote_string( $opts->{$var} ) );
    }
  } else {
    add_default($opts->{'test'}, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var );
  }
  $nl_flags->{'clm_start_type'} = $nl->get_value($var);
}

#-------------------------------------------------------------------------------

sub setup_cmdl_dynamic_vegetation {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  my $val;
  my $var = "dynamic_vegetation";
  $val = $opts->{$var};
  $nl_flags->{'dynamic_vegetation'} = $val;
  if ($nl_flags->{'phys'} eq "clm4_0") {
    # not applicable
    if ( $nl_flags->{'dynamic_vegetation'}eq 1) {
       fatal_error("** Turn dynamic_vegetation mode on with CLM_CONFIG_OPTS (-bgc cndv) for clm4_0 physics.\n" );
    }
  } else {
    if ( ($nl_flags->{'dynamic_vegetation'} eq 1 ) && ($nl_flags->{'bgc_mode'} eq "sp") ) {
      fatal_error("** Cannot turn dynamic_vegetation mode on with bgc=sp.\n" .
                  "**\n" .
                  "** Set the bgc mode to 'cn' or 'bgc' by the following means from highest to lowest precedence:\n" .
                  "** * by the command-line options -bgc cn\n");
    }

    $var = "use_cndv";
    $nl_flags->{$var} = ".false.";
    if ($nl_flags->{'dynamic_vegetation'} eq 1) {
      $val = ".true.";
      $nl_flags->{$var} = $val;
    }
    if ( defined($nl->get_value($var)) && $nl->get_value($var) ne $val ) {
      fatal_error("$var is inconsistent with the commandline setting of -dynamic_vegetation");
    }
    if ( $nl_flags->{$var} eq ".true." ) {
      my $group = $definition->get_group_name($var);
      $nl->set_variable_value($group, $var, $val);
      if (  ! $definition->is_valid_value( $var, $val ) ) {
        my @valid_values   = $definition->get_valid_values( $var );
        fatal_error("$var has a value ($val) that is NOT valid. Valid values are: @valid_values\n");
      }
    }
  }
}

#-------------------------------------------------------------------------------

sub setup_cmdl_vichydro {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  my $val;
  my $var = "vichydro";
  $val = $opts->{$var};
  $nl_flags->{'vichydro'} = $val;
  if ($nl_flags->{'phys'} eq "clm4_0") {
    # not relevant in clm4_0
    if ( $nl_flags->{'vichydro'}eq 1) {
       fatal_error("** Cannot turn vichydro on with clm4_0 physics.\n" );
    }
  } else {
    if ($nl_flags->{'vichydro'} eq 1) {
      message("Using VIC hydrology for runoff calculations.");
    }

    $var = "use_vichydro";
    $val = $nl->get_value($var);
    if ($nl_flags->{'vichydro'} eq 1) {
      my $group = $definition->get_group_name($var);
      my $set = ".true.";
      if ( defined($val) && $set ne $val ) {
        fatal_error("$var contradicts the command-line -vichydro option" );
      }
      $nl->set_variable_value($group, $var, $set);
      if ( ! $definition->is_valid_value($var, $val) ) {
        my @valid_values   = $definition->get_valid_values( $var );
        fatal_error("$var has a value ($val) that is NOT valid. Valid values are: @valid_values\n");
      }
    }
  }
}

#-------------------------------------------------------------------------------

sub process_namelist_commandline_namelist {
  # Process the commandline '-namelist' arg.
  my ($opts, $definition, $nl, $envxml_ref) = @_;

  if (defined $opts->{'namelist'}) {
    # Parse commandline namelist
    my $nl_arg = Build::Namelist->new($opts->{'namelist'});

    # Validate input namelist -- trap exceptions
    my $nl_arg_valid;
    eval { $nl_arg_valid = $definition->validate($nl_arg); };
    if ($@) {
      fatal_error("Invalid namelist variable in commandline arg '-namelist'.\n $@");
    }
    # Go through all variables and expand any XML env settings in them
    expand_env_variables_in_namelist( $nl_arg_valid, $envxml_ref );

    # Merge input values into namelist.  Previously specified values have higher precedence
    # and are not overwritten.
    $nl->merge_nl($nl_arg_valid);
  }
}

#-------------------------------------------------------------------------------

sub process_namelist_commandline_infile {
  # Process the commandline '-infile' arg.
  my ($opts, $definition, $nl, $envxml_ref) = @_;

  if (defined $opts->{'infile'}) {
    my @infiles = split( /,/, $opts->{'infile'} );
    foreach my $infile ( @infiles ) {
      # Make sure a valid file was found
      if (    -f "$infile" ) {
        # Otherwise abort as a valid file doesn't exist
      } else {
        fatal_error("input namelist file does NOT exist $infile.\n $@");
      }
      # Parse namelist input from the next file
      my $nl_infile = Build::Namelist->new($infile);

      # Validate input namelist -- trap exceptions
      my $nl_infile_valid;
      eval { $nl_infile_valid = $definition->validate($nl_infile); };
      if ($@) {
        fatal_error("Invalid namelist variable in '-infile' $infile.\n $@");
      }
      # Go through all variables and expand any XML env settings in them
      expand_env_variables_in_namelist( $nl_infile_valid, $envxml_ref );

      # Merge input values into namelist.  Previously specified values have higher precedence
      # and are not overwritten.
      $nl->merge_nl($nl_infile_valid);
    }
  }
}

#-------------------------------------------------------------------------------

sub process_namelist_commandline_clm_usr_name {
  # Process the -clm_usr_name argument
  my ($opts, $nl_flags, $definition, $defaults, $nl, $cfg, $envxml_ref) = @_;

  if (defined $opts->{'clm_usr_name'}) {
    # The user files definition is contained in an xml file with the same format as the defaults file.

    # The one difference is that variables are expanded.
    # Create a new NamelistDefaults object.
    my $nl_defaults_file = "$nl_flags->{'cfgdir'}/namelist_files/namelist_defaults_usr_files.xml";
    my $uf_defaults = Build::NamelistDefaults->new("$nl_defaults_file", $cfg );
    # Loop over the variables specified in the user files
    # Add each one to the namelist.
    my @vars = $uf_defaults->get_variable_names();
    my %settings;
    $settings{'mask'}           = $nl_flags->{'mask'};
    $settings{'sim_year'}       = $nl_flags->{'sim_year'};
    $settings{'rcp'}            = $nl_flags->{'rcp'};
    $settings{'sim_year_range'} = $nl_flags->{'sim_year_range'};
    $settings{'bgc_spinup'}     = $nl_flags->{'bgc_spinup'};
    $settings{'clm_usr_name'}   = $opts->{'clm_usr_name'};

    if ( $nl_flags->{'inputdata_rootdir'} eq "\$DIN_LOC_ROOT" ) {
      $settings{'csmdata'}     = $ENV{'DIN_LOC_ROOT'};
    } else {
      $settings{'csmdata'}     = $nl_flags->{'inputdata_rootdir'};
    }

    my $nvars = 0;
    my $nl_usrfile = Build::Namelist->new();
    foreach my $var (@vars) {
      my $val = $uf_defaults->get_usr_file($var, $definition, \%settings);

      if ($val) {
        message("adding clm user file defaults for var $var with val $val");
        add_default($opts->{'test'}, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl_usrfile, $var, 'val'=>$val);
        $nvars++;
      }
    }
    if ( $nvars == 0 ) {
      warning("setting clm_usr_name -- but did NOT find any user datasets: $opts->{'clm_usr_name'}\n");
    }
    # Go through all variables and expand any XML env settings in them
    expand_env_variables_in_namelist( $nl_usrfile, $envxml_ref );
    # Merge input values into namelist.  Previously specified values have higher precedence
    # and are not overwritten.
    $nl->merge_nl($nl_usrfile);
  }
}

#-------------------------------------------------------------------------------

sub process_namelist_commandline_use_case {
  # Now process the -use_case arg.
  my ($opts, $nl_flags, $definition, $defaults, $nl, $cfg, $envxml_ref) = @_;

  if (defined $opts->{'use_case'}) {

    # The use case definition is contained in an xml file with the same format as the defaults file.
    # Create a new NamelistDefaults object.
    my $uc_defaults = Build::NamelistDefaults->new("$opts->{'use_case_dir'}/$opts->{'use_case'}.xml", $cfg);

    my %settings;
    $settings{'res'}            = $nl_flags->{'res'};
    $settings{'rcp'}            = $nl_flags->{'rcp'};
    $settings{'mask'}           = $nl_flags->{'mask'};
    $settings{'sim_year'}       = $nl_flags->{'sim_year'};
    $settings{'sim_year_range'} = $nl_flags->{'sim_year_range'};
    $settings{'phys'}           = $nl_flags->{'phys'};
    if ($nl_flags->{'phys'} eq "clm4_5") {
      $settings{'use_cn'}         = $nl_flags->{'use_cn'};
    } else {
      $settings{'bgc'}         = $nl_flags->{'bgc_mode'};
    }
    # Loop over the variables specified in the use case.
    # Add each one to the namelist.
    my @vars = $uc_defaults->get_variable_names();
    my $nl_usecase = Build::Namelist->new();
    foreach my $var (@vars) {
      my $val = $uc_defaults->get_value($var, \%settings );

      if ( defined($val) ) {
        message("CLM adding use_case $opts->{'use_case'} defaults for var '$var' with val '$val'");

        add_default($opts->{'test'}, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl_usecase, $var, 'val'=>$val);
      }
    }
    # Go through all variables and expand any XML env settings in them
    expand_env_variables_in_namelist( $nl_usecase, $envxml_ref );

    # Merge input values into namelist.  Previously specified values have higher precedence
    # and are not overwritten.
    $nl->merge_nl($nl_usecase);
  }
}

#-------------------------------------------------------------------------------

sub process_namelist_commandline_clm_start_type {
  # Set the start_type according to the command line clm_start_type option

  my ($test_files, $nl_flags, $definition, $defaults, $nl) = @_;

  # Run type for driver namelist - note that arb_ic implies that the run is startup
  my $var = "start_type";
  if ($nl_flags->{'clm_start_type'} eq "'cold'" || $nl_flags->{'clm_start_type'} eq "'arb_ic'") {
    # Add default is used here, but the value is explicitly set
    add_default($test_files, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var, 'val'=>'startup'   );
  } else {
    # Add default is used here, but the value is explicitly set
    add_default($test_files, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var, 'val'=>$nl_flags->{'clm_start_type'} );
  }
}

#-------------------------------------------------------------------------------

sub process_namelist_inline_logic {
  # Use the namelist default object to add default values for required
  # namelist variables that have not been previously set.
  my ($opts, $nl_flags, $definition, $defaults, $nl, $cfg, $envxml_ref) = @_;

  ##############################
  # namelist group: clm_inparm #
  ##############################
  setup_logic_site_specific($nl_flags, $definition, $nl);
  setup_logic_lnd_frac($opts, $nl_flags, $definition, $defaults, $nl, $envxml_ref);
  setup_logic_co2_type($opts, $nl_flags, $definition, $defaults, $nl);
  setup_logic_irrigate($opts, $nl_flags, $definition, $defaults, $nl);
  setup_logic_start_type($nl_flags, $nl);
  setup_logic_delta_time($opts, $nl_flags, $definition, $defaults, $nl);
  setup_logic_decomp_performance($opts->{'test'}, $nl_flags, $definition, $defaults, $nl);
  setup_logic_snow($opts->{'test_files'}, $nl_flags, $definition, $defaults, $nl);
  setup_logic_glacier($opts, $nl_flags, $definition, $defaults, $nl, $envxml_ref);
  setup_logic_params_file($opts->{'test'}, $nl_flags, $definition, $defaults, $nl);
  setup_logic_create_crop_landunit($opts->{'test'}, $nl_flags, $definition, $defaults, $nl);
  setup_logic_urban($opts->{'test'}, $nl_flags, $definition, $defaults, $nl);
  setup_logic_more_vertlayers($opts->{'test'}, $nl_flags, $definition, $defaults, $nl);
  setup_logic_demand($opts, $nl_flags, $definition, $defaults, $nl);
  setup_logic_surface_dataset($opts->{'test'}, $nl_flags, $definition, $defaults, $nl);
  setup_logic_initial_conditions($opts, $nl_flags, $definition, $defaults, $nl);
  setup_logic_bgc_spinup($opts->{'test'}, $nl_flags, $definition, $defaults, $nl);
  setup_logic_supplemental_nitrogen($opts->{'test'}, $nl_flags, $definition, $defaults, $nl);
  setup_logic_hydrology_switches($nl);

  ###############################
  # namelist group: ch4par_in   #
  ###############################
  setup_logic_methane($opts->{'test'}, $nl_flags, $definition, $defaults, $nl);
  setup_logic_c_isotope($nl_flags, $definition, $defaults, $nl);

  ###############################
  # namelist group: ndepdyn_nml #
  ###############################
  setup_logic_nitrogen_deposition($opts->{'test'}, $nl_flags, $definition, $defaults, $nl);

  #################################
  # namelist group: popd_streams  #
  #################################
  setup_logic_popd_streams($opts->{'test'}, $nl_flags, $definition, $defaults, $nl);

  ##################################
  # namelist group: light_streams  #
  ##################################
  setup_logic_lightning_streams($opts->{'test'}, $nl_flags, $definition, $defaults, $nl);

  #################################
  # namelist group: drydep_inparm #
  #################################
  setup_logic_dry_deposition($opts, $nl_flags, $definition, $defaults, $nl);

  #################################
  # namelist group: megan_emis_nl #
  #################################
  setup_logic_megan($opts, $nl_flags, $definition, $defaults, $nl);

}

#-------------------------------------------------------------------------------

sub setup_logic_site_specific {
  # site specific requirements
  my ($nl_flags, $definition, $nl) = @_;

  if ($nl_flags->{'phys'} eq "clm4_5") {
    # res check prevents polluting the namelist with an unnecessary
    # false variable for every run
    if ($nl_flags->{'res'} eq "1x1_vancouverCAN") {
      my $var = "use_vancouver";
      my $val = ".true.";
      my $group = $definition->get_group_name($var);
      $nl->set_variable_value($group, $var, $val);
    }
  }

  if ($nl_flags->{'phys'} eq "clm4_5") {
    # res check prevents polluting the namelist with an unnecessary
    # false variable for every run
    if ($nl_flags->{'res'} eq "1x1_mexicocityMEX") {
      my $var = "use_mexicocity";
      my $val = ".true.";
      my $group = $definition->get_group_name($var);
      $nl->set_variable_value($group, $var, $val);
    }
  }

  if ($nl_flags->{'phys'} eq "clm4_5" && $nl_flags->{'res'} eq "1x1_smallvilleIA") {
    if ($nl_flags->{'use_cn'} ne ".true." || $nl_flags->{'use_crop'} ne ".true.") {
      fatal_error("1x1_smallvilleIA grids must use a compset with CN and CROP turned on.\n");
    }
  }

  if ($nl_flags->{'phys'} eq "clm4_5" && $nl_flags->{'res'} eq "1x1_numaIA") {
    if ($nl_flags->{'use_cn'} ne ".true." || $nl_flags->{'use_crop'} ne ".true.") {
      fatal_error("1x1_numaIA grids must use a compset with CN and CROP turned on.\n");
    }
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_lnd_frac {

  my ($opts, $nl_flags, $definition, $defaults, $nl, $envxml_ref) = @_;

  my $var = "lnd_frac";
  if ( defined($opts->{$var}) ) {
    if ( defined($nl->get_value('fatmlndfrc')) ) {
      fatal_error("Can NOT set both -lnd_frac option (set via LND_DOMAIN_PATH/LND_DOMAIN_FILE " .
                  "env variables) AND fatmlndfrac on namelist\n");
    }
    my $lnd_frac = &Streams::TemplateGeneric::expandXMLVar( $opts->{$var}, $envxml_ref);
    add_default($opts->{'test'}, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'fatmlndfrc','val'=>$lnd_frac );
  }

  # Get the fraction file
  if (defined $nl->get_value('fatmlndfrc')) {
    # do nothing - use value provided by config_grid.xml and clm.cpl7.template
  } else {
    fatal_error("fatmlndfrc was NOT sent into CLM build-namelist.\n");
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_co2_type {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  my $var = "co2_type";
  if ( defined($opts->{$var}) ) {
    if ( ! defined($nl->get_value($var)) ) {
      add_default($opts->{'test'}, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'co2_type','val'=>"$opts->{'co2_type'}");
    } else {
      fatal_error("co2_type set on namelist as well as -co2_type option.\n");
    }
  }
  add_default($opts->{'test'}, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'co2_type');
  if ( $nl->get_value('co2_type') =~ /constant/ ) {
    my $var = 'co2_ppmv';
    if ( defined($opts->{$var}) ) {
      if ( $opts->{$var} <= 0.0 ) {
        fatal_error("co2_ppmv can NOT be less than or equal to zero.");
      }
      my $group = $definition->get_group_name($var);
      $nl->set_variable_value($group, $var, $opts->{$var});
    } else {
      add_default($opts->{'test'}, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var, 'sim_year'=>$nl_flags->{'sim_year'} );
    }
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_irrigate {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  if ($nl_flags->{'phys'} eq "clm4_5") {
    if ( $nl_flags->{'use_crop'} eq ".true." ) {
      add_default($opts->{'test'}, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'irrigate', 'val'=>$nl_flags->{'irrig'});
    }
    elsif ( defined($nl->get_value('irrigate')) ) {
      if ($nl->get_value('irrigate') =~ /$TRUE/i ) {
        fatal_error("irrigate TRUE needs crop TRUE but it is not\n");
      }
    }
    $nl_flags->{'irrigate'} = lc($nl->get_value('irrigate'));
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_start_type {
  my ($nl_flags, $nl) = @_;

  my $var = "start_type";
  my $drv_start_type = $nl->get_value($var);
  my $my_start_type  = $nl_flags->{'clm_start_type'};
  my $nsrest         = $nl->get_value('override_nsrest');

  if ( defined($nsrest) ) {
    if ( $nsrest == 0 ) { $my_start_type = "startup";  }
    if ( $nsrest == 1 ) { $my_start_type = "continue"; }
    if ( $nsrest == 3 ) { $my_start_type = "branch";   }
    if ( "$my_start_type" eq "$drv_start_type" ) {
      fatal_error("no need to set override_nsrest to same as start_type.\n");
    }
    if ( "$drv_start_type" !~ /startup/ ) {
      fatal_error("can NOT set override_nsrest if driver is NOT a startup type.\n");
    }
  }

  if ( $my_start_type =~ /branch/ ) {
    if (not defined $nl->get_value('nrevsn')) {
      fatal_error("nrevsn is required for a branch type.\n");
    }
  } else {
    if (defined $nl->get_value('nrevsn')) {
      fatal_error("nrevsn should ONLY be set for a branch type.\n");
    }
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_delta_time {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  if ( defined($opts->{'l_ncpl'}) ) {
    my $l_ncpl = $opts->{'l_ncpl'};
    if ( $l_ncpl <= 0 ) {
      fatal_error("bad value for -l_ncpl option.\n");
    }
    my $val = ( 3600 * 24 ) / $l_ncpl;
    my $dtime = $nl->get_value('dtime');
    if ( ! defined($dtime)  ) {
      add_default($opts->{'test'}, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'dtime', 'val'=>$val);
    } elsif ( $dtime ne $val ) {
      fatal_error("can NOT set both -l_ncpl option (via LND_NCPL env variable) AND dtime namelist variable.\n");
    }
  } else {
    add_default($opts->{'test'}, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'dtime', 'hgrid'=>$nl_flags->{'res'});
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_decomp_performance {
  my ($test_files, $nl_flags, $definition, $defaults, $nl) = @_;

  # Set the number of segments per clump
  add_default($test_files, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'nsegspc', 'hgrid'=>$nl_flags->{'res'});
}

#-------------------------------------------------------------------------------

sub setup_logic_snow {
  my ($test_files, $nl_flags, $definition, $defaults, $nl) = @_;

  add_default($test_files, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'fsnowoptics' );
  add_default($test_files, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'fsnowaging' );
}

#-------------------------------------------------------------------------------

sub setup_logic_glacier {
  #
  # Glacier multiple elevation class options
  #
  my ($opts, $nl_flags, $definition, $defaults, $nl, $envxml_ref) = @_;

  if ($nl_flags->{'phys'} eq "clm4_5") {
     # glc_do_dynglacier is set via CLM_UPDATE_GLC_AREAS; it cannot be set via
     # user_nl_clm (this is because we might eventually want the coupler and glc
     # to also respond to CLM_UPDATE_GLC_AREAS, by not bothering to send / map
     # these fields - so we want to ensure that CLM is truly listening to this
     # shared xml variable and not overriding it)
     my $var = "glc_do_dynglacier";
     my $val = logical_to_fortran($envxml_ref->{'CLM_UPDATE_GLC_AREAS'});
     add_default($opts->{'test'}, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var, 'val'=>$val);
     if (lc($nl->get_value($var)) ne lc($val)) {
        fatal_error("glc_do_dynglacier can only be set via the env variable CLM_UPDATE_GLC_AREAS: it can NOT be set in user_nl_clm\n");
     }

     # glc_dyn_runoff_routing is also set via CLM_UPDATE_GLC_AREAS by default,
     # but unlike glc_do_dynglacier, it can be overridden via user_nl_clm
     $var = "glc_dyn_runoff_routing";
     add_default($opts->{'test'}, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var, 'val'=>$val);     
  }

  my $var = "maxpatch_glcmec";
  add_default($opts->{'test'}, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var, 'val'=>$nl_flags->{'glc_nec'} );

  my $val = $nl->get_value($var);
  if ( $val != $nl_flags->{'glc_nec'} ) {
    fatal_error("$var set to $val does NOT agree with -glc_nec argument of $nl_flags->{'glc_nec'} (set with GLC_NEC env variable)\n");
  }
  if ( $nl_flags->{'glc_nec'} > 0 ) {
    foreach my $var ( "glc_grid", "glc_smb" ) {
      if ( $opts->{$var} ne "default" ) {
        add_default($opts->{'test'}, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var, 'val'=>$opts->{$var} );
        $val = $nl->get_value($var);
        $val =~ s/['"]//g;
        my $ucvar = $var;
        $ucvar =~ tr/a-z/A-Z/;
        if ( $val ne $opts->{$var} ) {
          fatal_error("$var set to $val does NOT agree with -$var argument of $opts->{$var} (set with $ucvar env variable)\n");
        }
      } else {
        add_default($opts->{'test'}, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var, 'glc_nec'=>$nl_flags->{'glc_nec'} );
      }
      $val = $nl->get_value($var);
      verbose_message("Glacier model $var is $val");
      if ( ! defined($val) ) {
        fatal_error("$var is NOT set, but glc_nec is positive");
      }
    }
    my $glc_grid = $nl->get_value('glc_grid');
    $glc_grid =~ s/['"]//g;
    add_default($opts->{'test'}, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'flndtopo'  , 'hgrid'=>$nl_flags->{'res'}, 'mask'=>$nl_flags->{'mask'} );
    add_default($opts->{'test'}, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'fglcmask'  , 'hgrid'=>$nl_flags->{'res'}, 'glc_grid'=>$glc_grid );

    if ($nl_flags->{'phys'} eq "clm4_5") {
      add_default($opts->{'test'}, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'glcmec_downscale_rain_snow_convert');
      add_default($opts->{'test'}, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'glcmec_downscale_longwave');
      add_default($opts->{'test'}, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'glc_snow_persistence_max_days');
    }

  } else {
    if (defined($nl->get_value('glc_grid'))) {
      my $glc_grid = $nl->get_value('glc_grid');
      $glc_grid =~ s/['"]//g;
      if ( defined($glc_grid) ) {
        fatal_error("glc_grid is set, but glc_nec is zero");
      }
    }
    # Error checking for glacier multiple elevation class options when glc_mec off
    # Make sure various glc_mec-specific logicals are not true, and fglcmask is not set
    my $create_glcmec = $nl->get_value('create_glacier_mec_landunit');
    if ( defined($create_glcmec) ) {
      if ( $create_glcmec =~ /$TRUE/i ) {
        fatal_error("create_glacer_mec_landunit is true, but glc_nec is equal to zero");
      }
    }
    my $glc_smb = $nl->get_value('glc_smb');
    if ( defined($glc_smb) ) {
      if ( $glc_smb =~ /$TRUE/i ) {
        fatal_error("glc_smb is true, but glc_nec is equal to zero");
      }
    }
    my $glc_dyntopo= $nl->get_value('glc_dyntopo');
    if ( defined($glc_dyntopo) ) {
      if ( $glc_dyntopo =~ /$TRUE/i ) {
        fatal_error("glc_dyntopo is true, but glc_nec is equal to zero");
      }
    }
    my $glc_do_dynglacier= $nl->get_value('glc_do_dynglacier');
    if ( defined($glc_do_dynglacier) ) {
      if ( $glc_do_dynglacier =~ /$TRUE/i ) {
        fatal_error("glc_do_dynglacier (set from CLM_UPDATE_GLC_AREAS env variable) is true, but glc_nec is equal to zero");
      }
    }
    my $glc_dyn_runoff_routing= $nl->get_value('glc_dyn_runoff_routing');
    if ( defined($glc_dyn_runoff_routing) ) {
      if ( $glc_dyn_runoff_routing =~ /$TRUE/i ) {
        fatal_error("glc_dyn_runoff_routing is true, but glc_nec is equal to zero");
      }
    }
    my $fglcmask = $nl->get_value('fglcmask');
    if ( defined($fglcmask) ) {
      fatal_error("fglcmask is set, but glc_nec is equal to zero");
    }
  }

  add_default($opts->{'test'}, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'albice', 'glc_nec'=>$nl_flags->{'glc_nec'});
}

#-------------------------------------------------------------------------------

sub setup_logic_params_file {
  # get param data. For 4_0, pft-physiology, for 4_5 old
  # pft-physiology was used but now now includes CN and BGC century
  # parameters.
  my ($test_files, $nl_flags, $definition, $defaults, $nl) = @_;

  if ($nl_flags->{'phys'} eq "clm4_5") {
    add_default($test_files, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'paramfile');
  } else {
    add_default($test_files, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'fpftcon');
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_create_crop_landunit {
  # Create crop land unit
  my ($test_files, $nl_flags, $definition, $defaults, $nl) = @_;

  if ($nl_flags->{'phys'} eq "clm4_0") {
    if ( $nl_flags->{'crop'} eq "on" ) {
      add_default($test_files, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'create_crop_landunit' );
    }
    add_default($test_files, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'create_crop_landunit', 'irrig'=>$nl_flags->{'irrig'} );
  } else {
    add_default($test_files, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'create_crop_landunit', 'use_crop'=>$nl_flags->{'use_crop'});
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_urban {
  my ($test_files, $nl_flags, $definition, $defaults, $nl) = @_;

  add_default($test_files, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'urban_hac');
  add_default($test_files, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'urban_traffic');
}

#-------------------------------------------------------------------------------

sub setup_logic_more_vertlayers {
  my ($test_files, $nl_flags, $definition, $defaults, $nl) = @_;

  if ($nl_flags->{'phys'} eq "clm4_5") {
    add_default($test_files, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'more_vertlayers', 'hgrid'=>$nl_flags->{'res'} );
    $nl_flags->{'more_vert'} = $nl->get_value('more_vertlayers');
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_demand {
  #
  # Deal with options that the user has said are required...
  #
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  my %settings;
  $settings{'hgrid'}          = $nl_flags->{'res'};
  $settings{'sim_year'}       = $nl_flags->{'sim_year'};
  $settings{'sim_year_range'} = $nl_flags->{'sim_year_range'};
  $settings{'mask'}           = $nl_flags->{'mask'};
  $settings{'crop'}           = $nl_flags->{'crop'};
  $settings{'irrig'}          = $nl_flags->{'irrig'};
  $settings{'rcp'}            = $nl_flags->{'rcp'};
  $settings{'glc_nec'}        = $nl_flags->{'glc_nec'};
  if ($nl_flags->{'phys'} eq "clm4_5") {
    # necessary for demand to be set correctly (fpftdyn requires
    # use_crop, maybe other options require other flags?)!
    $settings{'use_cn'}              = $nl_flags->{'use_cn'};
    $settings{'use_cndv'}            = $nl_flags->{'use_cndv'};
    $settings{'use_lch4'}            = $nl_flags->{'use_lch4'};
    $settings{'use_nitrif_denitrif'} = $nl_flags->{'use_nitrif_denitrif'};
    $settings{'use_vertsoilc'}       = $nl_flags->{'use_vertsoilc'};
    $settings{'use_century_decomp'}  = $nl_flags->{'use_century_decomp'};
    $settings{'use_crop'}            = $nl_flags->{'use_crop'};
  }

  my $demand = $nl->get_value('clm_demand');
  if (defined($demand)) {
    $demand =~ s/\'//g;   # Remove quotes
    if ( $demand =~ /.+/ ) {
      $opts->{'clm_demand'} .= ",$demand";
    }
  }

  $demand = $defaults->get_value('clm_demand', \%settings);
  if (defined($demand)) {
    $demand =~ s/\'//g;   # Remove quotes
    if ( $demand =~ /.+/ ) {
      $opts->{'clm_demand'} .= ",$demand";
    }
  }

  my @demandlist = split( ",", $opts->{'clm_demand'} );
  foreach my $item ( @demandlist ) {
    if ( $item eq "null" ) {
      next;
    }
    add_default($opts->{'test'}, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $item, %settings );
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_surface_dataset {
  #
  # Get surface dataset after fpftdyn so that we can get surface data
  # consistent with it
  # MUST BE AFTER: setup_logic_demand which is where fpftdyn is set
  #
  my ($test_files, $nl_flags, $definition, $defaults, $nl) = @_;

  $nl_flags->{'fpftdyn'} = "null";
  my $fpftdyn = $nl->get_value('fpftdyn');
  if (defined($fpftdyn)) {
    $fpftdyn =~ s!(.*)/!!;
    $fpftdyn =~ s/\'//;
    $fpftdyn =~ s/\"//;
    if ( $fpftdyn ne "" ) {
      $nl_flags->{'fpftdyn'} = $fpftdyn;
    }
  }
  $fpftdyn = $nl_flags->{'fpftdyn'};

  if ($nl_flags->{'phys'} eq "clm4_0") {
    if ($fpftdyn ne "null" && $nl_flags->{'bgc_mode'} eq "cndv" ) {
        fatal_error( "dynamic PFT's (setting fpftdyn) are incompatible with dynamic vegetation ('-bgc cndv' in CLM_CONFIG_OPTS)." );
    }

    add_default($test_files, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'fsurdat',
                'hgrid'=>$nl_flags->{'res'},
                'sim_year'=>$nl_flags->{'sim_year'}, 'irrig'=>$nl_flags->{'irrig'},
                'crop'=>$nl_flags->{'crop'}, 'glc_nec'=>$nl_flags->{'glc_nec'});
  } else{
    if ($fpftdyn ne "null" && $nl_flags->{'use_cndv'} =~ /$TRUE/i ) {
        fatal_error( "dynamic PFT's (setting fpftdyn) are incompatible with dynamic vegetation (use_cndv=.true)." );
    }
    add_default($test_files, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'fsurdat',
                'hgrid'=>$nl_flags->{'res'},
                'sim_year'=>$nl_flags->{'sim_year'}, 'irrig'=>$nl_flags->{'irrig'},
                'use_crop'=>$nl_flags->{'use_crop'}, 'glc_nec'=>$nl_flags->{'glc_nec'});
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_initial_conditions {
  # Initial conditions
  # The initial date is an attribute in the defaults file which should be matched unless
  # the user explicitly requests to ignore the initial date via the -ignore_ic_date option,
  # or just ignore the year of the initial date via the -ignore_ic_year option.
  #
  # MUST BE AFTER: setup_logic_demand   which is where fpftdyn is set
  #         AFTER: setup_logic_irrigate which is where irrigate is set
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  if ( $nl_flags->{'clm_start_type'} =~ /cold/ ) {
    if (defined $nl->get_value('finidat')) {
      fatal_error("setting finidat is incomptable with using start_type=cold.");
    }
    add_default($opts->{'test'}, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl,
                'finidat', 'val'=>"' '", 'no_abspath'=>1);
  }

  if (not defined $nl->get_value('finidat')) {
    my $ic_date = $nl->get_value('start_ymd');
    my $nofail = 1;
    my $var = "finidat";
    if ( $nl_flags->{'clm_start_type'} =~ /startup/  ) { $nofail = 0; }
    if ($opts->{'ignore_ic_date'}) {
      if ( $nl_flags->{'use_crop'} eq ".true." ) {
        fatal_error("using ignore_ic_date is incompatable with crop!");
      }
      if ($nl_flags->{'phys'} eq "clm4_0") {
        add_default($opts->{'test'}, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var,
                    'hgrid'=>$nl_flags->{'res'}, 'mask'=>$nl_flags->{'mask'},
                    'nofail'=>$nofail, 'fpftdyn'=>$nl_flags->{'fpftdyn'},
                    'sim_year'=>$nl_flags->{'sim_year'}, 'maxpft'=>$nl_flags->{'maxpft'},
                    'irrig'=>$nl_flags->{'irrig'}, 'glc_nec'=>$nl_flags->{'glc_nec'},
                    'crop'=>$nl_flags->{'crop'}, 'bgc'=>$nl_flags->{'bgc_mode'} );
      } else {
        add_default($opts->{'test'}, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var,
                    'hgrid'=>$nl_flags->{'res'}, 'mask'=>$nl_flags->{'mask'},
                    'nofail'=>$nofail, 'fpftdyn'=>$nl_flags->{'fpftdyn'},
                    'use_cn'=>$nl_flags->{'use_cn'}, 'use_cndv'=>$nl_flags->{'use_cndv'},
                    'use_nitrif_denitrif'=>$nl_flags->{'use_nitrif_denitrif'},
                    'use_vertsoilc'=>$nl_flags->{'use_vertsoilc'},
                    'use_century_decomp'=>$nl_flags->{'use_century_decomp'},
                    'sim_year'=>$nl_flags->{'sim_year'}, 'maxpft'=>$nl_flags->{'maxpft'},
                    'more_vertlayers'=>$nl_flags->{'more_vert'},
                    'glc_nec'=>$nl_flags->{'glc_nec'}, 'use_crop'=>$nl_flags->{'use_crop'},
                    'irrigate'=>$nl_flags->{'irrigate'} );
      }
    } elsif ($opts->{'ignore_ic_year'}) {
      if ($nl_flags->{'phys'} eq "clm4_0") {
        add_default($opts->{'test'}, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var,
                    'hgrid'=>$nl_flags->{'res'}, 'mask'=>$nl_flags->{'mask'},
                    'ic_md'=>$ic_date, 'nofail'=>$nofail, 'fpftdyn'=>$nl_flags->{'fpftdyn'},
                    'sim_year'=>$nl_flags->{'sim_year'}, 'maxpft'=>$nl_flags->{'maxpft'},
                    'irrig'=>$nl_flags->{'irrig'}, 'glc_nec'=>$nl_flags->{'glc_nec'},
                    'crop'=>$nl_flags->{'crop'}, 'bgc'=>$nl_flags->{'bgc_mode'} );
      } else {
        add_default($opts->{'test'}, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var,
                    'hgrid'=>$nl_flags->{'res'}, 'mask'=>$nl_flags->{'mask'},
                    'ic_md'=>$ic_date, 'nofail'=>$nofail, 'fpftdyn'=>$nl_flags->{'fpftdyn'},
                    'use_cn'=>$nl_flags->{'use_cn'}, 'use_cndv'=>$nl_flags->{'use_cndv'},
                    'use_nitrif_denitrif'=>$nl_flags->{'use_nitrif_denitrif'},
                    'use_vertsoilc'=>$nl_flags->{'use_vertsoilc'},
                    'use_century_decomp'=>$nl_flags->{'use_century_decomp'},
                    'sim_year'=>$nl_flags->{'sim_year'}, 'maxpft'=>$nl_flags->{'maxpft'},
                    'more_vertlayers'=>$nl_flags->{'more_vert'},
                    'glc_nec'=>$nl_flags->{'glc_nec'}, 'use_crop'=>$nl_flags->{'use_crop'},
                    'irrigate'=>$nl_flags->{'irrigate'} );
      }
    } else {
      if ($nl_flags->{'phys'} eq "clm4_0") {
        add_default($opts->{'test'}, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var,
                    'hgrid'=>$nl_flags->{'res'}, 'mask'=>$nl_flags->{'mask'},
                    'ic_ymd'=>$ic_date, 'nofail'=>$nofail, 'fpftdyn'=>$nl_flags->{'fpftdyn'},
                    'sim_year'=>$nl_flags->{'sim_year'}, 'maxpft'=>$nl_flags->{'maxpft'},
                    'irrig'=>$nl_flags->{'irrig'}, 'glc_nec'=>$nl_flags->{'glc_nec'},
                    'crop'=>$nl_flags->{'crop'}, 'bgc'=>$nl_flags->{'bgc_mode'} );
      } else {
        add_default($opts->{'test'}, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var,
                    'hgrid'=>$nl_flags->{'res'}, 'mask'=>$nl_flags->{'mask'},
                    'ic_ymd'=>$ic_date, 'nofail'=>$nofail, 'fpftdyn'=>$nl_flags->{'fpftdyn'},
                    'use_cn'=>$nl_flags->{'use_cn'}, 'use_cndv'=>$nl_flags->{'use_cndv'},
                    'use_nitrif_denitrif'=>$nl_flags->{'use_nitrif_denitrif'},
                    'use_vertsoilc'=>$nl_flags->{'use_vertsoilc'},
                    'use_century_decomp'=>$nl_flags->{'use_century_decomp'},
                    'sim_year'=>$nl_flags->{'sim_year'}, 'maxpft'=>$nl_flags->{'maxpft'},
                    'more_vertlayers'=>$nl_flags->{'more_vert'},
                    'glc_nec'=>$nl_flags->{'glc_nec'}, 'use_crop'=>$nl_flags->{'use_crop'},
                    'irrigate'=>$nl_flags->{'irrigate'} );
      }
    }
    my $finidat = $nl->get_value($var);
    if ( (not defined $finidat ) || $finidat =~ /null/ ) {
      my $group = $definition->get_group_name($var);
      $nl->set_variable_value($group, $var, "' '" );
    }
  }
} # end initial conditions

#-------------------------------------------------------------------------------

sub setup_logic_bgc_spinup {
  my ($test_files, $nl_flags, $definition, $defaults, $nl) = @_;

  if ($nl_flags->{'phys'} eq "clm4_5") {
    if ( $nl_flags->{'bgc_mode'} ne "sp" ) {
      # only set bgc_spinup state if CN is on.
      add_default($test_files, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'spinup_state', 'bgc_spinup'=>$nl_flags->{'bgc_spinup'} );
    }

    if ( $nl_flags->{'bgc_mode'} eq "sp" && defined($nl->get_value('override_bgc_restart_mismatch_dump'))) {
      fatal_error("CN must be on if override_bgc_restart_mismatch_dump is set.\n");
    }
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_supplemental_nitrogen {
  #
  # Supplemental Nitrogen for prognostic crop cases
  #
  my ($test_files, $nl_flags, $definition, $defaults, $nl) = @_;

  if ( $nl_flags->{'bgc_mode'} ne "sp" && $nl_flags->{'use_crop'} eq ".true." ) {
    add_default($test_files, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl,
                'suplnitro', 'use_cn'=>$nl_flags->{'use_cn'}, 'use_crop'=>$nl_flags->{'use_crop'});
  }

  #
  # Error checking for suplnitro
  #
  my $suplnitro = $nl->get_value('suplnitro');
  if ( defined($suplnitro) ) {
    if ( $nl_flags->{'bgc_mode'} eq "sp" ) {
      fatal_error("supplemental Nitrogen (suplnitro) is set, but neither CN nor CNDV is active!\n");
    }
    if ( $nl_flags->{'use_crop'} ne ".true." && $suplnitro =~ /PROG_CROP_ONLY/i ) {
      fatal_error("supplemental Nitrogen is set to run over prognostic crops, but prognostic crop is NOT active!\n");
    }
    
    if ( $suplnitro =~ /ALL/i ) {
      if ( $nl_flags->{'phys'} eq "clm4_0" && $nl_flags->{'spinup'} ne "normal" ) {
        fatal_error("There is no need to use a spinup mode when supplemental Nitrogen is on for all PFT's, as these modes spinup Nitrogen\n" .
                    "when spinup != normal you can NOT set supplemental Nitrogen (suplnitro) to ALL\n");
      }
      if ( $nl_flags->{'phys'} eq "clm4_5" && $nl_flags->{'bgc_spinup'} ne "off" ) {
        warning("There is no need to use a bgc_spinup mode when supplemental Nitrogen is on for all PFT's, as these modes spinup Nitrogen\n");
      }
    }
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_hydrology_switches {
  #
  # Check on Switches for hydrology
  #
  my ($nl) = @_;

  my $subgrid    = $nl->get_value('subgridflag' ) || 0;
  my $origflag   = $nl->get_value('origflag'    ) || 0;
  my $h2osfcflag = $nl->get_value('h2osfcflag'  ) || 0;
  if ( $origflag == 1 && $subgrid == 1 ) {
    fatal_error("if origflag is ON, subgridflag can NOT also be on!");
  }
  if ( $h2osfcflag == 1 && $subgrid != 1 ) {
    fatal_error("if h2osfcflag is ON, subgridflag can NOT be off!");
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_methane {
  #
  # CH4 model if bgc=CN or CNDV
  #
  my ($test_files, $nl_flags, $definition, $defaults, $nl) = @_;

  if ( $nl_flags->{'use_lch4'}  eq '.true.' ) {
    add_default($test_files, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'fin_use_fsat' );
    add_default($test_files, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'use_aereoxid_prog' );
    #
    # Check if use_aereoxid_prog is set.  If no, then read value of aereoxid from
    # parameters file
    #
    my $use_aereoxid_prog = $nl->get_value('use_aereoxid_prog');
    if ( defined($use_aereoxid_prog) && $use_aereoxid_prog =~ /$FALSE/i ) {
      warning("Using aereoxid value from parameters file.\n");
    }
  } else {
    my @vars = $nl->get_variable_names('ch4par_in');
    if ( $#vars >= 0 ) {
      fatal_error("ch4par_in namelist variables were set, but Methane model NOT defined in the configuration (use_lch4)");
    }
  }

  #
  # Ch4 namelist checking
  #
  if ( $nl_flags->{'use_lch4'}  eq ".true." ) {
    my $allowlakeprod = $nl->get_value('allowlakeprod');
    if ( ! defined($allowlakeprod) ||
         (defined($allowlakeprod) && $allowlakeprod =~ /$FALSE/i) ) {
      if ( defined($nl->get_value('lake_decomp_fact')) ) {
        fatal_error("lake_decomp_fact set without allowlakeprod=.true.\n");
      }
    }
    my $anoxia = $nl->get_value('anoxia');
    if ( ! defined($anoxia) ||
         (defined($anoxia) && $anoxia =~ /$FALSE/i) ) {
      if ( defined($nl->get_value('anoxia_wtsat')) ) {
        fatal_error("anoxia_wtsat set without anoxia=.true.\n");
      }
    }
    my $pftspec_rootprof = $nl->get_value('pftspecific_rootingprofile');
    if ( ! defined($pftspec_rootprof) ||
         (defined($pftspec_rootprof) && $pftspec_rootprof =~ /$TRUE/i) ) {
      if ( defined($nl->get_value('rootprof_exp')) ) {
        fatal_error("rootprof_exp set without pftspecific_rootingprofile=.false.\n");
      }
    }
  } else {
    my @vars = ( "allowlakeprod", "anoxia", "anoxia_wtsat", "pftspecific_rootingprofile" );
    foreach my $var ( @vars ) {
      if ( defined($nl->get_value($var)) ) {
        fatal_error("$var set without methane model configuration on (use_lch4)\n");
      }
    }
  }
} # end methane

#-------------------------------------------------------------------------------

sub setup_logic_c_isotope {
  #
  # Error checking for C-isotope options
  #
  my ($nl_flags, $definition, $defaults, $nl) = @_;

  my $use_c13 = $nl->get_value('use_c13');
  my $use_c14 = $nl->get_value('use_c14');
  if ( $nl_flags->{'bgc_mode'} ne "sp" ) {
    if ( $nl_flags->{'use_crop'} eq ".true." ) {
      if ( defined($use_c13) ||
           defined($use_c14) ||
           defined($nl->get_value('use_c14_bombspike')) ||
           defined($nl->get_value('atm_c14_filename')) ) {
        fatal_error("CROP is on and C isotope  namelist variables were set, both can't be used at the same time");
      }
    }
    if ( $nl_flags->{'bgc_mode'} ne "bgc" ) {
      if ( defined($use_c13) && $use_c13 =~ /$TRUE/i ) {
        warning("use_c13 is ONLY scientifically validated with the bgc=BGC configuration\n");
      }
      if ( defined($use_c14) && $use_c14 =~ /$TRUE/i ) {
        warning("use_c14 is ONLY scientifically validated with the bgc=BGC configuration\n");
      }
    }
    if ( defined($use_c14) ) {
      if ( $use_c14 =~ /$TRUE/i ) {
        my $use_c14_bombspike = $nl->get_value('use_c14_bombspike');
        if ( defined($use_c14_bombspike) && $use_c14_bombspike =~ /$TRUE/i &&
             ! defined($nl->get_value('atm_c14_filename')) ) {
          fatal_error("use_c14_bombspike TRUE but atm_c14_filename NOT set\n");
        }
      } else {
        if ( defined($nl->get_value('use_c14_bombspike')) ||
             defined($nl->get_value('atm_c14_filename')) ) {
          fatal_error("use_c14 is FALSE and use_c14_bombspike or atm_c14_filename set\n");
        }
      }
    } else {
      if ( defined($nl->get_value('use_c14_bombspike')) ||
           defined($nl->get_value('atm_c14_filename')) ) {
        fatal_error("use_c14 NOT set to .true., but use_c14_bompspike/atm_c14_filename defined.\n");
      }
    }
  } else {
    if ( defined($use_c13) ||
         defined($use_c14) ||
         defined($nl->get_value('use_c14_bombspike')) ||
         defined($nl->get_value('atm_c14_filename')) ) {
           fatal_error("bgc=sp and C isotope  namelist variables were set, both can't be used at the same time");
    }
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_nitrogen_deposition {
  my ($test_files, $nl_flags, $definition, $defaults, $nl) = @_;

  #
  # Nitrogen deposition for bgc=CN
  #

  if ($nl_flags->{'phys'} eq "clm4_0" && $nl_flags->{'bgc_mode'} ne "none") {
    add_default($test_files, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'ndepmapalgo', 'phys'=>$nl_flags->{'phys'}, 
                'bgc'=>$nl_flags->{'bgc_mode'}, 'hgrid'=>$nl_flags->{'res'} );
    add_default($test_files, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_year_first_ndep', 'phys'=>$nl_flags->{'phys'},
                'bgc'=>$nl_flags->{'bgc_mode'}, 'sim_year'=>$nl_flags->{'sim_year'},
                'sim_year_range'=>$nl_flags->{'sim_year_range'});
    add_default($test_files, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_year_last_ndep', 'phys'=>$nl_flags->{'phys'},
                'bgc'=>$nl_flags->{'bgc_mode'}, 'sim_year'=>$nl_flags->{'sim_year'},
                'sim_year_range'=>$nl_flags->{'sim_year_range'});

    # Set align year, if first and last years are different
    if ( $nl->get_value('stream_year_first_ndep') != $nl->get_value('stream_year_last_ndep') ) {
      add_default($test_files, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'model_year_align_ndep', 'sim_year'=>$nl_flags->{'sim_year'},
                  'sim_year_range'=>$nl_flags->{'sim_year_range'});
    }

    add_default($test_files, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_fldfilename_ndep', 'phys'=>$nl_flags->{'phys'},
                'bgc'=>$nl_flags->{'bgc_mode'}, 'rcp'=>$nl_flags->{'rcp'},
                'hgrid'=>"1.9x2.5" );

  } elsif ($nl_flags->{'phys'} eq "clm4_5" && $nl_flags->{'bgc_mode'} ne "sp") {
    add_default($test_files, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'ndepmapalgo', 'phys'=>$nl_flags->{'phys'},
                'use_cn'=>$nl_flags->{'use_cn'}, 'hgrid'=>$nl_flags->{'res'} );
    add_default($test_files, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_year_first_ndep', 'phys'=>$nl_flags->{'phys'},
                'use_cn'=>$nl_flags->{'use_cn'}, 'sim_year'=>$nl_flags->{'sim_year'},
                'sim_year_range'=>$nl_flags->{'sim_year_range'});
    add_default($test_files, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_year_last_ndep', 'phys'=>$nl_flags->{'phys'},
                'use_cn'=>$nl_flags->{'use_cn'}, 'sim_year'=>$nl_flags->{'sim_year'},
                'sim_year_range'=>$nl_flags->{'sim_year_range'});
    # Set align year, if first and last years are different
    if ( $nl->get_value('stream_year_first_ndep') != $nl->get_value('stream_year_last_ndep') ) {
      add_default($test_files, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'model_year_align_ndep', 'sim_year'=>$nl_flags->{'sim_year'},
                  'sim_year_range'=>$nl_flags->{'sim_year_range'});
    }
    add_default($test_files, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_fldfilename_ndep', 'phys'=>$nl_flags->{'phys'},
                'use_cn'=>$nl_flags->{'use_cn'}, 'rcp'=>$nl_flags->{'rcp'},
                'hgrid'=>"1.9x2.5" );
  } else {
    # If bgc is NOT CN/CNDV then make sure none of the ndep settings are set!
    if ( defined($nl->get_value('stream_year_first_ndep')) ||
         defined($nl->get_value('stream_year_last_ndep'))  ||
         defined($nl->get_value('model_year_align_ndep'))  ||
         defined($nl->get_value('stream_fldfilename_ndep'))
       ) {
      fatal_error("When bgc is NOT CN or CNDV none of: stream_year_first_ndep," .
                  "stream_year_last_ndep, model_year_align_ndep, nor stream_fldfilename_ndep" .
                  " can be set!\n");
    }
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_popd_streams {
  # population density streams require clm4_5 and CN/BGC
  my ($test_files, $nl_flags, $definition, $defaults, $nl) = @_;

  if ( $nl_flags->{'phys'} eq "clm4_5" ) {
    if ( $nl_flags->{'bgc_mode'} ne "sp" ) {
      add_default($test_files, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'popdensmapalgo', 'hgrid'=>$nl_flags->{'res'} );
      add_default($test_files, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_year_first_popdens', 'phys'=>$nl_flags->{'phys'},
                  'use_cn'=>$nl_flags->{'use_cn'}, 'sim_year'=>$nl_flags->{'sim_year'},
                  'sim_year_range'=>$nl_flags->{'sim_year_range'});
      add_default($test_files, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_year_last_popdens', 'phys'=>$nl_flags->{'phys'},
                  'use_cn'=>$nl_flags->{'use_cn'}, 'sim_year'=>$nl_flags->{'sim_year'},
                  'sim_year_range'=>$nl_flags->{'sim_year_range'});
      # Set align year, if first and last years are different
      if ( $nl->get_value('stream_year_first_popdens') != 
           $nl->get_value('stream_year_last_popdens') ) {
        add_default($test_files, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'model_year_align_popdens', 'sim_year'=>$nl_flags->{'sim_year'},
                    'sim_year_range'=>$nl_flags->{'sim_year_range'});
      }
      add_default($test_files, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_fldfilename_popdens', 'phys'=>$nl_flags->{'phys'},
                  'use_cn'=>$nl_flags->{'use_cn'}, 'hgrid'=>"0.5x0.5" );
    } else {
      # If bgc is NOT CN/CNDV then make sure none of the popdens settings are set
      if ( defined($nl->get_value('stream_year_first_popdens')) ||
           defined($nl->get_value('stream_year_last_popdens'))  ||
           defined($nl->get_value('model_year_align_popdens'))  ||
           defined($nl->get_value('stream_fldfilename_popdens'))   ) {
        fatal_error("When bgc is SP (NOT CN or BGC) none of: stream_year_first_popdens,\n" .
                    "stream_year_last_popdens, model_year_align_popdens, nor\n" .
                    "stream_fldfilename_popdens can be set!\n");
      }
    }
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_lightning_streams {
  # lightning streams require clm4_5 and CN/BGC
  my ($test_files, $nl_flags, $definition, $defaults, $nl) = @_;

  if ( $nl_flags->{'phys'} eq "clm4_5" ) {
    if ( $nl_flags->{'bgc_mode'} ne "sp" ) {
      add_default($test_files, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'lightngmapalgo', 'use_cn'=>$nl_flags->{'use_cn'},
                  'hgrid'=>$nl_flags->{'res'} );
      add_default($test_files, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_year_first_lightng', 'use_cn'=>$nl_flags->{'use_cn'},
                  'sim_year'=>$nl_flags->{'sim_year'},
                  'sim_year_range'=>$nl_flags->{'sim_year_range'});
      add_default($test_files, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_year_last_lightng', 'use_cn'=>$nl_flags->{'use_cn'},
                  'sim_year'=>$nl_flags->{'sim_year'},
                  'sim_year_range'=>$nl_flags->{'sim_year_range'});
      # Set align year, if first and last years are different
      if ( $nl->get_value('stream_year_first_lightng') != 
           $nl->get_value('stream_year_last_lightng') ) {
        add_default($test_files, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'model_year_align_lightng', 'sim_year'=>$nl_flags->{'sim_year'},
                    'sim_year_range'=>$nl_flags->{'sim_year_range'});
      }
      add_default($test_files, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_fldfilename_lightng', 'use_cn'=>$nl_flags->{'use_cn'},
                  'hgrid'=>"94x192" );
    } else {
      # If bgc is NOT CN/CNDV then make sure none of the Lightng settings are set
      if ( defined($nl->get_value('stream_year_first_lightng')) ||
           defined($nl->get_value('stream_year_last_lightng'))  ||
           defined($nl->get_value('model_year_align_lightng'))  ||
           defined($nl->get_value('stream_fldfilename_lightng'))   ) {
        fatal_error("When bgc is SP (NOT CN or BGC) none of: stream_year_first_lightng,\n" .
                    "stream_year_last_lightng, model_year_align_lightng, nor\n" .
                    "stream_fldfilename_lightng can be set!\n");
      }
    }
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_dry_deposition {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  if ($opts->{'drydep'} ) {
    add_default($opts->{'test'}, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'drydep_list');
    add_default($opts->{'test'}, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'drydep_method');
  } else {
    if ( defined($nl->get_value('drydep_list')) ||
         defined($nl->get_value('drydep_method')) ) {
      fatal_error("drydep_list or drydep_method defined, but drydep option NOT set\n");
    }
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_megan {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  if ($opts->{'megan'} ) {
    add_default($opts->{'test'}, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'megan_specifier');
    check_megan_spec( $nl, $definition );
    add_default($opts->{'test'}, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'megan_factors_file');
  } else {
    if ( defined($nl->get_value('megan_specifier')) ||
         defined($nl->get_value('megan_factors_file')) ) {
      fatal_error("megan_specifier or megan_factors_file defined, but megan option NOT set\n");
    }
  }
}

#-------------------------------------------------------------------------------

sub write_output_files {
  my ($opts, $nl_flags, $defaults, $nl) = @_;

  my $note = "";
  my $var = "note";
  if ( ! defined($opts->{$var}) ) {
    $opts->{$var} = $defaults->get_value($var);
  }
  if ( $opts->{$var} ) {
    $note = "Comment:\n" .
      "This namelist was created using the following command-line:\n" .
        "    $nl_flags->{'cfgdir'}/$ProgName $nl_flags->{'cmdline'}\n" .
          "For help on options use: $nl_flags->{'cfgdir'}/$ProgName -help";
  }

  # CLM component
  my @groups;
  if ($nl_flags->{'phys'} eq "clm4_0") {
    @groups = qw(clm_inparm ndepdyn_nml);
  } else {
    @groups = qw(clm_inparm ndepdyn_nml popd_streams light_streams clm_hydrology1_inparm clm_soilhydrology_inparm finidat_consistency_checks dynpft_consistency_checks);
    if ( $nl_flags->{'use_lch4'}  eq ".true." ) {
      push @groups, "ch4par_in";
    }
  }

  my $outfile;
  $outfile = "$opts->{'dir'}/lnd_in";
  $nl->write($outfile, 'groups'=>\@groups, 'note'=>"$note" );
  verbose_message("Writing clm namelist to $outfile");

  # Drydep or MEGAN namelist
  if ($opts->{'drydep'} || $opts->{'megan'} ) {
    @groups = qw(drydep_inparm megan_emis_nl);
    $outfile = "$opts->{'dir'}/drv_flds_in";
    $nl->write($outfile, 'groups'=>\@groups, 'note'=>"$note" );
    verbose_message("Writing @groups namelists to $outfile");
  }
}

#-------------------------------------------------------------------------------

sub add_default {

# Add a value for the specified variable to the specified namelist object.  The variables
# already in the object have the higher precedence, so if the specified variable is already
# defined in the object then don't overwrite it, just return.
#
# This method checks the definition file and adds the variable to the correct
# namelist group.
#
# The value can be provided by using the optional argument key 'val' in the
# calling list.  Otherwise a default value is obtained from the namelist
# defaults object.  If no default value is found this method throws an exception
# unless the 'nofail' option is set true.
#
# Example 1: Specify the default value $val for the namelist variable $var in namelist
#            object $nl:
#
#  add_default($inputdata_rootdir, $definition, $defaults, $nl, $var, 'val'=>$val)
#
# Example 2: Add a default for variable $var if an appropriate value is found.  Otherwise
#            don't add the variable
#
#  add_default($inputdata_rootdir, $definition, $defaults, $nl, $var, 'nofail'=>1)
#
#
# ***** N.B. ***** This routine assumes the following variables are in package main::
#  $definition        -- the namelist definition object
#  $defaults          -- the namelist defaults object
#  $inputdata_rootdir -- CESM inputdata root directory

  my $test_files = shift;
  my $inputdata_rootdir = shift;
  my $definition = shift;
  my $defaults = shift;
  my $nl = shift;
  my $var = shift;
  my %settings = @_;

  #my $nl = shift;     # namelist object
  #my $var = shift;    # name of namelist variable
  #my %settings = @_;      # options

  # If variable has quotes around it
  if ( $var =~ /'(.+)'/ ) {
    $var = $1;
  }
  # Query the definition to find which group the variable belongs to.  Exit if not found.
  my $group = $definition->get_group_name($var);
  unless ($group) {
    my $fname = $definition->get_file_name();
    fatal_error("variable \"$var\" not found in namelist definition file $fname.\n");
  }

  # check whether the variable has a value in the namelist object -- if so then skip to end
  my $val = $nl->get_variable_value($group, $var);
  if (! defined $val) {

    # Look for a specified value in the options hash

    if (defined $settings{'val'}) {
      $val = $settings{'val'};
    }
    # or else get a value from namelist defaults object.
    # Note that if the 'val' key isn't in the hash, then just pass anything else
    # in %settings to the get_value method to be used as attributes that are matched
    # when looking for default values.
    else {
      $val = $defaults->get_value($var, \%settings);

      # Truncate model_version appropriately

      if ( $var eq "model_version" ) {
        $val =~ /(URL: https:\/\/[a-zA-Z0-9._-]+\/)([a-zA-Z0-9\/._-]+)(\/bld\/.+)/;
        $val = $2;
      }
    }

    # if no value is found then exit w/ error (unless 'nofail' option set)
    unless ( defined($val) ) {
      unless ($settings{'nofail'}) {
        if ($var eq 'finidat') {
          warning("No default value found for $var.\n" .
                  "            Are defaults provided for this resolution and land mask?\n");
        } else {
          fatal_error("No default value found for $var.\n" .
                      "            Are defaults provided for this resolution and land mask?\n");
        }
      }
      else {
        return;
      }
    }

    # query the definition to find out if the variable is an input pathname
    my $is_input_pathname = $definition->is_input_pathname($var);

    # The default values for input pathnames are relative.  If the namelist
    # variable is defined to be an absolute pathname, then prepend
    # the CESM inputdata root directory.
    if (not defined $settings{'no_abspath'}) {
      if (defined $settings{'set_abspath'}) {
        $val = set_abs_filepath($val, $settings{'set_abspath'});
      } else {
        if ($is_input_pathname eq 'abs') {
          $val = set_abs_filepath($val, $inputdata_rootdir);
          if ( $test_files and ($val !~ /null/) and (! -f "$val") ) {
            fatal_error("file not found: $var = $val");
          }
        }
      }
    }

    # query the definition to find out if the variable takes a string value.
    # The returned string length will be >0 if $var is a string, and 0 if not.
    my $str_len = $definition->get_str_len($var);

    # If the variable is a string, then add quotes if they're missing
    if ($str_len > 0) {
      $val = quote_string($val);
    }

    # set the value in the namelist
    $nl->set_variable_value($group, $var, $val);
  }

}

#-------------------------------------------------------------------------------

sub expand_env_variables_in_namelist {
   # Go through all variables in the namelist and expand any XML env settings in them
   my ($nl, $envxml_ref) = @_;
   foreach my $group ( $nl->get_group_names() ) {
       foreach my $var ( $nl->get_variable_names($group) ) {
          my $val    = $nl->get_variable_value($group, $var);
          my $newval = &Streams::TemplateGeneric::expandXMLVar( $val, $envxml_ref );
          if ( $newval ne $val ) {
             $nl->set_variable_value($group, $var, $newval);
          }
       }
   }
}

#-------------------------------------------------------------------------------

sub check_input_files {

# For each variable in the namelist which is an input dataset, check to see if it
# exists locally.
#
# ***** N.B. ***** This routine assumes the following variables are in package main::
#  $definition        -- the namelist definition object
#  $nl                -- namelist object
#  $inputdata_rootdir -- if false prints test, else creates inputdata file

    my ($nl, $inputdata_rootdir, $outfile, $definition) = @_;

    open(OUTFILE, ">>$outfile") if defined $inputdata_rootdir;

    # Look through all namelist groups
    my @groups = $nl->get_group_names();
    foreach my $group (@groups) {

        # Look through all variables in each group
        my @vars = $nl->get_variable_names($group);
        foreach my $var (@vars) {

            # Is the variable an input dataset?
            my $input_pathname_type = $definition->is_input_pathname($var);

            # If it is, check whether it exists locally and print status
            if ($input_pathname_type) {

                # Get pathname of input dataset
                my $pathname = $nl->get_variable_value($group, $var);
                # Need to strip the quotes
                $pathname =~ s/['"]//g;

                if ($input_pathname_type eq 'abs') {
                    if ($inputdata_rootdir) {
                        #MV $pathname =~ s:$inputdata_rootdir::;
                        print OUTFILE "$var = $pathname\n";
                    }
                    else {
                        if (-e $pathname) {  # use -e rather than -f since the absolute pathname
                                             # might be a directory
                            print "OK -- found $var = $pathname\n";
                        }
                        else {
                            print "NOT FOUND:  $var = $pathname\n";
                        }
                    }
                }
                elsif ($input_pathname_type =~ m/rel:(.+)/o) {
                    # The match provides the namelist variable that contains the
                    # root directory for a relative filename
                    my $rootdir_var = $1;
                    my $rootdir = $nl->get_variable_value($group, $rootdir_var);
                    $rootdir =~ s/['"]//g;
                    if ($inputdata_rootdir) {
                        $pathname = "$rootdir/$pathname";
                        #MV $pathname =~ s:$inputdata_rootdir::;
                        print OUTFILE "$var = $pathname\n";
                    }
                    else {
                        if (-f "$rootdir/$pathname") {
                            print "OK -- found $var = $rootdir/$pathname\n";
                        }
                        else {
                            print "NOT FOUND:  $var = $rootdir/$pathname\n";
                        }
                    }
                }
            }
        }
    }
    close OUTFILE if defined $inputdata_rootdir;
    return 0 if defined $inputdata_rootdir;
}

#-------------------------------------------------------------------------------

sub set_abs_filepath {

# check whether the input filepath is an absolute path, and if it isn't then
# prepend a root directory

    my ($filepath, $rootdir) = @_;

    # strip any leading/trailing whitespace
    $filepath =~ s/^\s+//;
    $filepath =~ s/\s+$//;
    $rootdir  =~ s/^\s+//;
    $rootdir  =~ s/\s+$//;

    # strip any leading/trailing quotes
    $filepath =~ s/^['"]+//;
    $filepath =~ s/["']+$//;
    $rootdir =~ s/^['"]+//;
    $rootdir =~ s/["']+$//;

    my $out = $filepath;
    unless ( $filepath =~ /^\// ) {  # unless $filepath starts with a /
        $out = "$rootdir/$filepath"; # prepend the root directory
    }
    return $out;
}

#-------------------------------------------------------------------------------

sub valid_option {

    my ($val, @expect) = @_;

    my $expect;

    $val =~ s/^\s+//;
    $val =~ s/\s+$//;
    foreach $expect (@expect) {
        if ($val =~ /^$expect$/i) { return $expect; }
    }
    return undef;
}

#-------------------------------------------------------------------------------

sub check_use_case_name {
#
# Check the use-case name and ensure it follows the naming convention.
#
  my ($use_case) = @_;

  my $diestring = "bad use_case name $use_case, follow the conventions " .
                  "in namelist_files/use_cases/README\n";
  my $desc = "[a-zA-Z0-9]*";
  my $rcp  = "rcp[0-9\.]+";
  if (      $use_case =~ /^[0-9]+-[0-9]+([a-zA-Z0-9_\.]*)_transient$/ ) {
    my $string = $1;
    if (      $string =~ /^_($rcp)_*($desc)$/ ) {
       # valid name
    } elsif ( $string =~ /^_*($desc)$/ ) {
       # valid name
    } else {
      fatal_error($diestring);
    }
  } elsif ( $use_case =~ /^20thC([a-zA-Z0-9_\.]*)_transient$/ ) {
    my $string = $1;
    if (      $string =~ /^_($rcp)_*($desc)$/ ) {
       # valid name
    } elsif ( $string =~ /^_*($desc)$/ ) {
       # valid name
    } else {
      fatal_error($diestring);
    }
  } elsif ( $use_case =~ /^([0-9]+)_*($desc)_control$/   ) {
     # valid name
  } elsif ( $use_case =~ /^($desc)_pd$/   ) {
     # valid name
  } else {
      fatal_error($diestring);
  }
}

#-------------------------------------------------------------------------------

sub validate_options {

# $source -- text string declaring the source of the options being validated
# $cfg    -- configure object
# $opts   -- reference to hash that contains the options

    my ($source, $cfg, $opts) = @_;

    my ($opt, $old, @expect);

    # use_case
    $opt = 'use_case';
    if (defined $opts->{$opt}) {

        if ( $opts->{$opt} ne "list" ) {
           # create the @expect array by listing the files in $use_case_dir
           # and strip off the ".xml" part of the filename
           @expect = ();
           my @files = glob("$opts->{'use_case_dir'}/*.xml");
           foreach my $file (@files) {
               $file =~ m{.*/(.*)\.xml};
               &check_use_case_name( $1 );
               push @expect, $1;
           }

           $old = $opts->{$opt};
           $opts->{$opt} = valid_option($old, @expect)
               or fatal_error("invalid value of $opt ($old) specified in $source\n" .
                              "expected one of: @expect");
        } else {
           print "Use cases are:...\n\n";
           my @ucases;
           foreach my $file( sort( glob($opts->{'use_case_dir'}."/*.xml") ) ) {
              my $use_case;
              if ( $file =~ /\/([^\/]+)\.xml$/ ) {
                 &check_use_case_name( $1 );
                 $use_case = $1;
              } else {
                 fatal_error("Bad name for use case file = $file");
              }
              my $uc_defaults = Build::NamelistDefaults->new("$file", $cfg);
              printf "%15s = %s\n", $use_case, $uc_defaults->get_value("use_case_desc");
              push @ucases, $use_case;
           }
           exit_message("use cases : @ucases");
        }
    }
}

#-------------------------------------------------------------------------------

sub list_options {
#
# List the options for different command line values if asked for
#
    my ($opts_cmdl, $definition, $defaults) = @_;

    # options to list values that are in the defaults files
    my @opts_list = ( "res", "mask", "sim_year", "rcp" );
    my %opts_local;
    foreach my $var ( "res", "mask", "sim_year", "rcp" ) {
       my $val;
       if (      $opts_cmdl->{$var} eq "list" ) {
         $val = "default";
       } elsif ( $opts_cmdl->{$var} eq "default" ) {
         $val = $defaults->get_value($var, \%opts_local );
       } else {
         $val = $opts_cmdl->{$var};
       }
       my $vname = $var;
       if ( $vname eq "res" ) { $vname = "hgrid"; }
       $opts_local{$vname} = $val;
    }
    foreach my $var ( @opts_list ) {
       if (defined $opts_cmdl->{$var}) {

           if ( $opts_cmdl->{$var} eq "list" ) {
               my @valid_values   = $definition->get_valid_values( $var );
               if ( $var eq "sim_year" ) {
                   unshift( @valid_values,
                            $definition->get_valid_values( "sim_year_range" ) );
               }
               unshift( @valid_values, "default" );
               # Strip out quotes and the constant value
               for( my $i = 0; $i <= $#valid_values; $i++ ) {
                  $valid_values[$i] =~ s/('|')//g;
                  if ( $valid_values[$i] eq "constant" ) { $valid_values[$i] = undef; }
               }
               my $val= $defaults->get_value($var, \%opts_local);
               my $doc = $definition->get_var_doc( $var );
               $doc =~ s/\n//;
               chomp( $doc );
               exit_message("valid values for $var ($doc) :\n" .
                            "    Values: @valid_values\n" .
                            "    Default = $val\n" .
                            "    (NOTE: resolution and mask and other settings may influence what the default is)");
           }
       }
    }
    # clm_demand
    my $var = 'clm_demand';
    if (defined $opts_cmdl->{$var}) {

        if ( $opts_cmdl->{$var} eq "list" ) {
           my @vars = $definition->get_var_names( );
           my @demands = ( "null" );
           foreach my $var ( @vars ) {
              if ( $definition->get_group_name( $var ) ne "clm_inparm" ) { next; }
              if ( defined($defaults->get_value($var, $opts_cmdl ) ) ) {
                 push( @demands, $var );
              }
           }
           my $doc = $definition->get_var_doc( 'clm_demand' );
           $doc =~ s/\n//;
           chomp( $doc );
           exit_message("valid values for $var ($doc) :\n" .
                        "Namelist options to require: @demands\n" .
                        "any valid namelist item for clm_inparm can be set. However, not all are\n" .
                        "available in the clm defaults file. The defaults are also dependent on\n" .
                        "resolution and landmask, as well as other settings. Hence, the list above\n" .
                        "will vary depending on what you set for resolution and landmask.\n");
        }
    }
}

#-------------------------------------------------------------------------------

sub check_megan_spec {
#
# Check the megan specifier setting
#
    my ($nl, $definition) = @_;

    my $megan_spec      = $nl->get_value('megan_specifier');
    my @megan_spec_list = split( /\s*,\s*/, $megan_spec );
    foreach $megan_spec ( @megan_spec_list ) {
       if ( $megan_spec =~ /^['"]+[A-Za-z0-9]+\s*\=\s*([\sA-Za-z0-9+_-]+)["']+$/ ) {
          my $megan_list = $1;
          my @megan_cmpds = split( /\s*\+\s*/, $megan_list );
          my $var = "megan_cmpds";
          my $warn = 0;
          foreach my $megan_cmpd ( @megan_cmpds ) {
             if (  ! $definition->is_valid_value( $var, $megan_cmpd, 'noquotes'=>1 ) ) {
                warning("megan_compound $megan_cmpd NOT found in list");
                $warn++;
             }
          }
          if ( $warn > 0 ) {
             my @valid_values   = $definition->get_valid_values( $var, 'noquotes'=>1 );
             warning("list of megan compounds includes:\n" .
                     "@valid_values\n" .
                     "Does your megan_factors_file include more coumpounds?\n" .
                     "If NOT your simulation will fail.\n");
          }
       } else {
          fatal_error("Bad format for megan_specifier = $megan_spec");
       }
    }
}

#-------------------------------------------------------------------------------

sub quote_string {
    # Add quotes around a string, unless they are already there
    my ($str) = @_;
    $str =~ s/^\s+//;
    $str =~ s/\s+$//;
    unless ($str =~ /^['"]/) {        #"'
        $str = "\'$str\'";
    }
    return $str;
}

#-------------------------------------------------------------------------------

sub logical_to_fortran {
   # Given a logical variable ('true' / 'false'), convert it to a fortran-style logical ('.true.' / '.false.')
   # The result will be lowercase, regardless of the case of the input.
   my ($var) = @_;
   my $result;

   if (lc($var) eq 'true') {
      $result = ".true.";
   }
   elsif (lc($var) eq 'false') {
      $result = ".false.";
   }
   else {
      fatal_error("Unexpected value in logical_to_fortran: $var\n");
   }

   return $result;
}

#-------------------------------------------------------------------------------

sub version {
# The version is found in CLM ChangeLog file.
# $cfgdir is set by the configure script to the name of its directory.

    my ($cfgdir) = @_;

    my $logfile = "$cfgdir/../doc/ChangeLog";

    my $fh = IO::File->new($logfile, '<') or fatal_error("can't open ChangeLog file: $logfile");

    while (my $line = <$fh>) {

        if ($line =~ /^Tag name:\s*([a-zA-Z0-9_. -]*[clmcesm0-9_.-]+)$/ ) {
            exit_message("$1\n");
        }
    }
}

#-------------------------------------------------------------------------------
# Some simple subroutines to print messages out

sub message {
  my ($message) = @_;
  print "$message\n";
}

sub verbose_message {
  my ($message) = @_;
  if ($verbosity >= $print_verbose) {
    print "$message\n";
  }
}

#-------------------------------------------------------------------------------
# Some simple subroutines to do a clean exit, print warning, or a fatal error

sub exit_message {
  my ($message) = @_;
  print "${ProgName} : $message\n";
  exit;
}

#-------------------------------------------------------------------------------

sub warning {
  my ($message) = @_;
  my $func_name = (caller(1))[3];
  print "Warning : ${ProgName}::${func_name}() : $message\n";
}

#-------------------------------------------------------------------------------

sub fatal_error {
  my ($message) = @_;
  my $func_name = (caller(1))[3];
  die "ERROR : ${ProgName}::${func_name}() : $message\n";
}

#-------------------------------------------------------------------------------

sub main {
  my %nl_flags;
  $nl_flags{'cfgdir'} = dirname(abs_path($0));

  my %opts = process_commandline(\%nl_flags);

  version($nl_flags{'cfgdir'}) if $opts{'version'};
  set_print_level(\%opts);

  check_for_perl_utils($nl_flags{'cfgdir'});
  my $cfg = read_configure_definition($nl_flags{'cfgdir'}, \%opts);

  $nl_flags{'phys'} = $cfg->get('phys'); # This defines clm4_5 versus clm4_0 physics

  my $drvblddir  = abs_path( "$nl_flags{'cfgdir'}/../../../../models/drv/bld" );
  my $definition = read_namelist_definition($drvblddir, \%opts, \%nl_flags);
  my %env_xml    = read_envxml_case_files( \%opts );
  my $defaults   = read_namelist_defaults($drvblddir, \%opts, \%nl_flags, $cfg);

  # List valid values if asked for
  list_options(\%opts, $definition, $defaults);

  # Validate some of the commandline option values.
  validate_options("commandline", $cfg, \%opts);

  # Create an empty namelist object.
  my $nl = Build::Namelist->new();

  check_cesm_inputdata(\%opts, \%nl_flags);

  # Process the user inputs
  process_namelist_user_input(\%opts, \%nl_flags, $definition, $defaults, $nl, $cfg, \%env_xml );
  # Get any other defaults needed from the namelist defaults file
  process_namelist_inline_logic(\%opts, \%nl_flags, $definition, $defaults, $nl, $cfg, \%env_xml);

  # Validate that the entire resultant namelist is valid
  $definition->validate($nl);
  write_output_files(\%opts, \%nl_flags, $defaults, $nl);

  if ($opts{'inputdata'}) {
    check_input_files($nl, $nl_flags{'inputdata_rootdir'}, $opts{'inputdata'}, $definition);
  }
}

#-------------------------------------------------------------------------------

1;
