#!/usr/bin/env perl
#=======================================================================
#
#  This is a script to read the CLM namelist XML file
#
# Usage:
#
# queryDefaultNamelist.pl [options]
#
# To get help on options and usage:
#
# queryDefaultNamelist.pl -help
#
#=======================================================================

use Cwd;
use strict;
#use diagnostics;
use Getopt::Long;
use English;

#-----------------------------------------------------------------------------------------------

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
my @dirs = ( $cfgdir, "$cfgdir/perl5lib",
             "$cfgdir/../../../../scripts/ccsm_utils/Tools/perl5lib",
             "$cfgdir/../../../../models/utils/perl5lib",
           );
unshift @INC, @dirs;
my $result = eval "require XML::Lite";
if ( ! defined($result) ) {
   die <<"EOF";
** Cannot find perl module \"XML/Lite.pm\" from directories: @dirs **
EOF
}
require Build::Config;
require Build::NamelistDefinition;
require queryDefaultXML;

# Defaults
my $namelist = "clm_inparm";
my $config = "config_cache.xml";


sub usage {
    die <<EOF;
SYNOPSIS
     $ProgName [options]

     query default namelist values.
OPTIONS
     -config "file"                       CLM build configuration file created by configure.
     -cesm                                CESM mode set csmdata to \$DIN_LOC_ROOT.
     -usrname "name"                      Dataset resolution/descriptor for personal datasets.  
                                          Default : not used
                                          Example: 1x1pt_boulderCO to describe location,
                                          number of pts
     -csmdata "dir"                       Directory for head of csm inputdata.
     -demand                              Demand that something is returned.
     -filenameonly                        Only return the filename -- not the full path to it.
     -help  [or -h]                       Display this help.
     -justvalue                           Just display the values (NOT key = value).
     -var "varname"                       Variable name to match. Use "-var list" to
                                          list valid variable names from all namelists.
     -namelist "namelistname"             Namelist name to read in (by default $namelist).
     -onlyfiles                           Only output filenames.
     -options "item=value,item2=value2"   Set options to query for when matching.
                                          (comma delimited, with equality to set value).
     -phys "CLM-version" [or -p]          CLM version to use (clm4_0 or clm4_5)
     -res  "resolution"                   Resolution to use for files. Use "-res list" to
                                          list all valid resolutions. Use "-res any" to
                                          use any valid resolution.
     -silent [or -s]                      Don't do any extra printing.
     -test   [or -t]                      Test that files exists.
EXAMPLES

  To list all fsurdat files that match the resolution: 10x15:

  $ProgName -var "fsurdat" -res 10x15

  To only list files that match T42 resolution (or are for all resolutions)

  $ProgName  -onlyfiles  -res 64x128

  To test that all of the files exist on disk under a different default inputdata

  $ProgName  -onlyfiles -test -csmdata /spin/proj/ccsm/inputdata

  To query for namelist items that match particular configurations

  $ProgName  -namelist seq_infodata_inparm -options sim_year=2000,bgc=cn

  Only lists namelist items in the seqinfodata_inparm namelist with options for
  sim_year=2000 and BGC=cn.

EOF
}

#-----------------------------------------------------------------------------------------------

  my %opts = ( 
               namelist   => $namelist,
               model      => "clm4_5",
               var        => undef,
               hgrid      => undef,
               config     => undef,
               cesm       => undef, 
               csmdata    => undef,
               demand     => undef,
               test       => undef,
               onlyfiles  => undef,
               fileonly   => undef,
               silent     => undef,
               usrname    => undef,
               help       => undef,
               options    => undef,
             );

  my $cmdline = "@ARGV";
  GetOptions(
        "f|file=s"     => \$opts{'file'},
        "n|namelist=s" => \$opts{'namelist'},
        "v|var=s"      => \$opts{'var'},
        "p|phys=s"     => \$opts{'model'},
        "r|res=s"      => \$opts{'hgrid'},
        "config=s"     => \$opts{'config'},
        "cesm"         => \$opts{'cesm'},
        "csmdata=s"    => \$opts{'csmdata'},
        "demand"       => \$opts{'demand'},
        "options=s"    => \$opts{'options'},
        "t|test"       => \$opts{'test'},
        "onlyfiles"    => \$opts{'onlyfiles'},
        "filenameonly" => \$opts{'fileonly'},
        "justvalues"   => \$opts{'justvalues'},
        "usrname=s"    => \$opts{'usrname'},
        "s|silent"     => \$opts{'silent'},
        "h|elp"        => \$opts{'help'},
  ) or usage();

  # Check for unparsed arguments
  if (@ARGV) {
      print "ERROR: unrecognized arguments: @ARGV\n";
      usage();
  }
  if ( $opts{'help'} ) {
      usage();
  }
  # Set if should do extra printing or not (if silent mode is not set)
  my $printing = 1;
  if ( defined($opts{'silent'}) ) {
      $printing = 0;
  }
  # Get list of options from command-line into the settings hash
  my %settings;
  if ( defined($opts{'options'}) ) {
     $opts{'options'} =~ s/\s//g;    # Remove all white-space in options
     my @optionlist = split( ",", $opts{'options'} );
     foreach my $item ( @optionlist ) {
        my ($key,$value) = split( "=", $item );
        $settings{$key} = $value;
     }
  }
  my $csmdata = "";
  if ( defined($opts{'fileonly'}) ) {
     if ( ! defined($opts{'justvalues'}) ) { print "When -filenameonly option used, -justvalues is set as well\n" if $printing; }
     if ( ! defined($opts{'onlyfiles'}) )  { print "When -filenameonly option used, -onlyfiles is set as well\n"  if $printing; }
     $opts{'justvalues'} = 1;
     $opts{'onlyfiles'}  = 1;
  }
  # List of input options
  my %inputopts;
  $inputopts{empty_cfg_file} = "$cfgdir/config_files/config_definition.xml";
  my $datmblddir             = "$cfgdir/../../../../models/atm/datm/bld";
  my $drvblddir              = "$cfgdir/../../../../models/drv/bld";
  my $model                  = $opts{'model'};
  my @nl_definition_files    = ( "$datmblddir/namelist_files/namelist_definition_datm.xml",
                                 "$drvblddir/namelist_files/namelist_definition_drv.xml",
                                 "$cfgdir/namelist_files/namelist_definition_$model.xml" 
                               );
  $inputopts{nldef_files}    = \@nl_definition_files;
  $inputopts{namelist}       = $opts{namelist};
  $inputopts{printing}       = $printing;
  $inputopts{cfgdir}         = $cfgdir;
  $inputopts{ProgName}       = $ProgName;
  $inputopts{cmdline}        = $cmdline;

  my $exitearly = 0;
  my $definition = Build::NamelistDefinition->new( $nl_definition_files[0] );
  foreach my $nl_defin_file ( @nl_definition_files ) {
     if ( ! -f "$nl_defin_file" ) {
        die "($ProgName $cmdline) ERROR:: bad namelist definition filename: $nl_defin_file.\n";
     }
     $definition->add( "$nl_defin_file" );
  }

  if ( ! defined($opts{csmdata}) ) {
     $inputopts{csmdata} = "default";
  } else {
     $inputopts{csmdata} = $opts{csmdata};
  }
  if ( defined($opts{cesm}) ) {
     $inputopts{csmdata} = '$DIN_LOC_ROOT';
  }
  if ( ! defined($opts{config}) ) {
     $inputopts{config} = "noconfig";
  } else {
     $inputopts{config} = $opts{config};
  }
  if (      ! defined($opts{var}) ) {
    $settings{'var'}  = undef;
  } elsif ( $opts{var} eq "list" ) {
     print "Valid variables: " if $printing;
     my @vars = $definition->get_var_names( );
     print "@vars\n";
     $exitearly = 1;
  } else {
    $settings{'var'}  = $opts{'var'};
  }
  if (      ! defined($opts{hgrid}) ) {
     $inputopts{hgrid} = "any";
  } elsif ( $opts{hgrid} eq "list" ) {
     print "Valid resolutions: " if $printing;
     my @hgrids = $definition->get_valid_values( "res", 'noquotes'=>1 );
     print "@hgrids\n";
     $exitearly = 1;
  } else {
     if ( ! $definition->is_valid_value( "res", $opts{hgrid}, 'noquotes'=>1 ) ) {
        if ( $opts{'hgrid'} ne $opts{'usrname'} ) {
           die "($ProgName $cmdline) ERROR:: invalid resolution entered.\n";
        }
     }
     $inputopts{hgrid} = $opts{hgrid};
  }
  # The namelist defaults file contains default values for all required namelist variables.
  my @nl_defaults_files = ( "$cfgdir/namelist_files/namelist_defaults_overall.xml" );
  if ( defined($opts{'usrname'}) ) {
     my $nl_defaults_file =  "$cfgdir/namelist_files/namelist_defaults_usr_files.xml";
     push( @nl_defaults_files, $nl_defaults_file );
     $settings{'clm_usr_name'} = $opts{'usrname'};
     $settings{'notest'}       = ! $opts{'test'};
     $settings{'csmdata'}      = $inputopts{csmdata};
  } else {
     my @files = ( "$cfgdir/namelist_files/namelist_defaults_${model}.xml", 
                   "$cfgdir/namelist_files/namelist_defaults_${model}_tools.xml", 
                   "$drvblddir/namelist_files/namelist_defaults_drv.xml",
                   "$cfgdir/namelist_files/namelist_defaults_drydep.xml",
                   "$datmblddir/namelist_files/namelist_defaults_datm.xml",
                 );
     push( @nl_defaults_files, @files );
  }
  if ( ! $exitearly ) {
    $inputopts{files} = \@nl_defaults_files;

    my $defaults_ref = &queryDefaultXML::ReadDefaultXMLFile( \%inputopts, \%settings );
    my %defaults = %$defaults_ref;
    my @keys = keys(%defaults);
    if ( defined($opts{'demand'}) && ($#keys == -1) ) {
       die "($ProgName $cmdline) ERROR:: demand option is set and nothing was found.\n";
    }
    my $print;
    foreach my $var ( @keys ) {
       $print = 1;
       my $value   = $defaults{$var}{value};
       my $isadir  = $defaults{$var}{isdir};
       my $isafile = $defaults{$var}{isfile};
       my $isastr  = $defaults{$var}{isstr};
       # If onlyfiles option set do NOT print if is NOT a file
       if ( defined($opts{'onlyfiles'}) && (! $isafile) ) {
          $print = undef;
       }
       # If is a directory
       if ( $isadir ) {
          # Test that this directory exists
          if ( defined($opts{'test'})  && defined($print) ) {
             print "Test that directory $value exists\n" if $printing;
             if ( ! -d "$value" ) {
                die "($ProgName) ERROR:: directory $value does NOT exist!\n";
             }
          }
       }
       # If is a file
       if ( $isafile ) {
          # Test that this file exists
          if ( defined($opts{'test'})  && defined($print) ) {
             chomp( $value );
             print "Test that file $value exists\n" if $printing;
             if ( ! -f "$value" ) {
                die "($ProgName) ERROR:: file $value does NOT exist!\n";
             }
          }
       }
       # If a string
       if ( (! defined($opts{'justvalues'}) ) && ($isastr) ) {
         $value = "\'$value\'";
       }
       # if you just want the filename -- not the full path with the directory
       if ( defined($opts{'fileonly'}) ) {
         $value =~ s!(.*)/!!;
       }
       if ( defined($print) ) {
          if ( ! defined($opts{'justvalues'})  ) {
             print "$var = ";
          }
          print "$value\n";
       }
    }
  }
  if ( $printing && defined($opts{'test'}) ) {
     print "\n\nTesting was successful\n\n"
  }

