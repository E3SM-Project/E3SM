#!/usr/bin/env perl
#=======================================================================
#
#  This is a script to list the missing files in your CESM inputdata area
#  for a list of resolutions and model configurations. The list goes
#  out to the file: clm.input_data_list. The check_input_data script
#  can then be used to get this list of files from the SVN inputdata
#  repository.
#
# Usage:
#
# listDefaultNamelist.pl [options]
#
# To get help on options and usage:
#
# listDefaultNamelist.pl -help
#
# To then get the files from the CESM SVN repository:
#
# ../../../../scripts/ccsm_utils/Tools/check_input_data -datalistdir . -export
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

my $printTimes = 0;

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
require queryDefaultXML;

# Defaults
my $datmblddir  = "$cfgdir/../../../../models/atm/datm/bld";
my $drvblddir   = "$cfgdir/../../../../models/drv/bld";

# The namelist defaults file contains default values for all required namelist variables.
my @nl_defaults_files = ( "$cfgdir/namelist_files/namelist_defaults_overall.xml",
                          "$drvblddir/namelist_files/namelist_defaults_drv.xml",
                          "$datmblddir/namelist_files/namelist_defaults_datm.xml" );
my $list = "clm.input_data_list";
my %list_of_all_files;

sub usage {
    die <<EOF;
SYNOPSIS
     $ProgName [options]
OPTIONS
     -help  [or -h]                       Display this help.
     -csmdata [or -d]                     Path to CSMDATA.
     -phys "name"                         Model physics (either clm4_0 or clm4_5) (default clm4_0)
     -res  "resolution1,resolution2,..."  List of resolution to use for files.
                                          (At least one resolution is required)
                                          (If res is "all" will run over all resolutions)
     -usrdat "name"                       Allow resolution to be the given clm user-data name
     -silent [or -s]                      Do not do any extra printing.
EXAMPLES

  List all the files needed for resolution 10x15 to the file $list.

  $ProgName -res 10x15

  List all the files needed for resolutions 10x15,4x5, and 64x128 from jaguar.
  with silent mode on so no extra printing is done:

  $ProgName -res 10x15,4x5,64x128 -s

  to then read the resulting clm.input_data_list file and retreive the files

  ../../../../scripts/ccsm_utils/Tools/check_input_data -datalistdir . -export

EOF
}

#-----------------------------------------------------------------------------------------------

sub GetListofNeededFiles {
#
# Get list of files that are needed to be copied to disk from the XML file.
#
  my $inputopts_ref = shift;
  my $settings_ref  = shift;
  my $files_ref     = shift;
  
  my $defaults_ref = &queryDefaultXML::ReadDefaultXMLFile( $inputopts_ref, $settings_ref );
  my @keys     = keys(%$defaults_ref);
  my $csmdata  = $$inputopts_ref{'csmdata'};
  my $printing = $$inputopts_ref{'printing'};
  foreach my $var ( @keys ) {
     my $value   = $$defaults_ref{$var}{'value'};
     my $isafile = $$defaults_ref{$var}{'isfile'};
     # If is a file
     if ( $isafile ) {
        $value    =~ m#$csmdata/(.+?)/([^/]+)$#;
        my $dir   = $1;
        my $file  = $2;
        
        # If file is already in the list then do NOT do anything
        if ( defined($list_of_all_files{"$dir/$file"} ) ) {
        # Test that this file exists
        } elsif ( -f "$value" ) {
           print "File $value exists\n" if $printing;
           $list_of_all_files{"$dir/$file"} = 1;
        } else {
        # If doesn't exist add it to the list of files to copy
           my $cfile = $$inputopts_ref{'scpfrom'} . "$dir/$file";
           my @dfiles;
           if ( defined($$files_ref{$dir}) ) {
              my $dir_ref = $$files_ref{$dir};
              @dfiles     = @$dir_ref;
              my $match = 0;
              foreach my $i ( @dfiles ) {
                if ( $i eq $cfile ) { $match = 1; }
              }
              if ( $match == 0 ) { push( @dfiles, $cfile ); }
           } else {
              @dfiles = ( "$cfile" );
           }
           if ( ! defined($$files_ref{$dir}) ) {
              print "           ADD $cfile to list to copy\n";
           }
           $$files_ref{$dir} = \@dfiles;
           $list_of_all_files{"$dir/$file"} = 0;
        }
     }
  }
  $printTimes++;
}

#-----------------------------------------------------------------------------------------------

  my %opts = ( 
               res        => undef,
               silent     => undef,
               csmdata    => "default",
               list       => $list,
               usrdat     => undef,
               help       => undef,
               phys       => "clm4_0",
             );

  my $cmdline = "@ARGV";
  GetOptions(
        "d|csmdata=s"  => \$opts{'csmdata'},
        "r|res=s"      => \$opts{'res'},
        "s|silent"     => \$opts{'silent'},
        "u|usrdat=s"   => \$opts{'usrdat'},
        "h|elp"        => \$opts{'help'},
        "phys=s"       => \$opts{'phys'},
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
  #
  # Check for required arguments
  #
  foreach my $req ( "res", "list" ) {
     if ( ! defined($opts{$req}) ) {
         print "ERROR: $req NOT set and it is a required argument\n";
         usage();
     }
  }
  my %inputopts;
  my $datmblddir             = "$cfgdir/../../../../models/atm/datm/bld";
  my @nl_definition_files    = (
                                 "$datmblddir/namelist_files/namelist_definition_datm.xml",
                                 "$cfgdir/namelist_files/namelist_definition_$opts{'phys'}.xml"
                               );
  $inputopts{'nldef_files'}    = \@nl_definition_files;
  $inputopts{'empty_cfg_file'} = "$cfgdir/config_files/config_definition.xml";

  my $definition = Build::NamelistDefinition->new( $nl_definition_files[0] );
  foreach my $nl_defin_file ( @nl_definition_files ) {
     $definition->add( "$nl_defin_file" );
  }
  my $cfg = Build::Config->new( $inputopts{'empty_cfg_file'} );

  # Resolutions...
  my @resolutions;
  if ( $opts{'res'} eq "all" ) {
     @resolutions = $definition->get_valid_values( "res", 'noquotes'=>1 );
  } else {
     @resolutions = split( /,/, $opts{'res'} );
  }

  # Input options
  push @nl_defaults_files, "$cfgdir/namelist_files/namelist_defaults_$opts{'phys'}.xml";
  if ( defined($opts{'usrdat'}) ) {
      push @nl_defaults_files, "$cfgdir/namelist_files/namelist_defaults_usr_files.xml";
  }
  $inputopts{'files'}     = \@nl_defaults_files;
  $inputopts{'printing'}  = $printing;
  $inputopts{'ProgName'}  = $ProgName;
  $inputopts{'cmdline'}   = $cmdline;
  $inputopts{'cfgdir'}    = $cfgdir;
  if ( $opts{'csmdata'} eq "default" && $ENV{'CSMDATA'} ne "" ) {
     $opts{'csmdata'} = $ENV{'CSMDATA'};
  }
  $inputopts{'csmdata'}   = $opts{'csmdata'};
  $inputopts{'config'}    = "noconfig";
  my %files;
  #
  # Loop over all resolutions asked for: 1.9x2.5, 10x15, 64x128 etc.
  #
  foreach my $res ( @resolutions ) {
     if ( ! $definition->is_valid_value( "res", "'$res'" ) && $res ne $opts{'usrdat'} ) {
        die "ERROR: Input resolution: $res is NOT a valid resolution\n";
     }
     $inputopts{'hgrid'} = $res;
     print "Resolution = $res\n" if $printing;
     my %settings;
     if ( $res eq $opts{'usrdat'} ) {
        $settings{'clm_usr_name'} = $opts{'usrdat'};
        $settings{'csmdata'}      = $opts{'csmdata'};
        $settings{'notest'}       = 1;
     }
     #
     # Loop for all possible land masks: USGS, gx1v6, gx3v5 etc.
     #
     foreach my $mask ( $definition->get_valid_values( "mask", 'noquotes'=>1 ) ) {
        print "Mask = $mask \n" if $printing;
        $settings{'mask'} = $mask;
        #
        # Loop over all possible simulation year: 1890, 2000, 2100 etc.
        #
        $settings{'sim_year_range'} = "constant";
        my @rcps = $definition->get_valid_values( "rcp", 'noquotes'=>1 );
        $settings{'rcp'} = $rcps[0];   
YEAR:   foreach my $sim_year ( $definition->get_valid_values( "sim_year", 'noquotes'=>1 ) ) {
           print "sim_year = $sim_year\n" if $printing;
           $settings{'sim_year'} = $sim_year;   
           if ( $sim_year ne 1850 && $sim_year ne 2000 && $sim_year > 1800 ) { next YEAR; }

           my @bgcsettings   = $cfg->get_valid_values( "bgc" );
           #my @glc_meclasses = $cfg->get_valid_values( "glc_nec" );
           my @glc_meclasses = ( 0, 10 );
           my @glc_grids     = $definition->get_valid_values( "glc_grid", 'noquotes'=>1 );
           print "glc_nec = @glc_meclasses bgc=@bgcsettings glc_grids=@glc_grids\n" if $printing;
           #
           # Loop over all possible BGC settings
           #
	   my $phys = $opts{'phys'};
           foreach my $bgc ( @bgcsettings ) {
              $settings{'bgc'} = $bgc;
              my @crop_vals;
              if ( $bgc =~ /^cn/ ) {
                 @crop_vals = ( "on", "off" );
              } else {
                 @crop_vals = ( "off" );
              }
              #
              # Loop over all possible glc_nec settings
              #
              foreach my $glc_nec ( @glc_meclasses ) {
                 $settings{'glc_nec'} = $glc_nec;
                 #
                 # Loop over all possible glc_grid settings
                 #
                 my @glcgrds;
                 if ( $glc_nec == 0 ) { @glcgrds = ( "none" );
                 } else               { @glcgrds = @glc_grids; } 
                 foreach my $glc_grid ( @glcgrds ) {
                    $settings{'glc_grid'} = $glc_grid;
                    #
                    # Loop over all possible crop settings
                    #
                    foreach my $crop ( @crop_vals ) {
                       $settings{'crop'} = $crop;
                       if ( $crop eq "on" ) {
			   if ($phys eq "clm4_0") { 
			       $settings{'maxpft'} = 21;
			   } else {
			       $settings{'maxpft'} = 25;
			   }
		       } else {
                          $settings{'maxpft'} = 17;
                       }
                       my @irrigset;
                       if ( $glc_nec  == 0 && $sim_year == 2000 ) { 
                          @irrigset= ( ".true.", ".false." );
                       } else { 
                          @irrigset= ( ".false." );
                       }
                       #
                       # Loop over irrigation settings
                       #
                       foreach my $irrig ( @irrigset ) {
                          $settings{'irrig'}     = $irrig;
                          $inputopts{'namelist'} = "clm_inparm";
                          &GetListofNeededFiles( \%inputopts, \%settings, \%files );
                          if ( $printTimes >= 1 ) {
                             $inputopts{'printing'} = 0;
                          }
                       }
                    }
                 }
              }
           }
        }
        #
        # Now do sim-year ranges
        #
        $settings{'bgc'}       = "cn";
        $settings{'irrig'}     = ".false.";
        $settings{'glc_grid'}  = "none";
        $inputopts{'namelist'} = "clm_inparm";
        foreach my $sim_year_range ( $definition->get_valid_values( "sim_year_range", 'noquotes'=>1 ) ) {
           $settings{'sim_year_range'} = $sim_year_range;
           if ( $sim_year_range =~ /([0-9]+)-([0-9]+)/ ) {
              $settings{'sim_year'}  = $1;
           }
           #
           # Loop over all possible rcp's
           #
           print "sim_year_range=$sim_year_range rcp=@rcps\n" if $printing;
           foreach my $rcp ( @rcps ) {
              $settings{'rcp'}       = $rcp;
              &GetListofNeededFiles( \%inputopts, \%settings, \%files );
              if ( $printTimes >= 1 ) {
                 $inputopts{'printing'} = 0;
              }
           }
        }
     }
  }
  #
  # Loop over directories that need to have files copied into
  #
  my $hostname;
  my $csmdata = $inputopts{'csmdata'};
  open( OUT, ">$list" ) || die "ERROR: trouble opening output file: $list";
  foreach my $dir ( sort(keys(%files)) ) {
     if ( $dir eq "."     ) { next; }
     if ( $dir eq "/"     ) { next; }
     if ( $dir eq "\n"    ) { next; }
     if ( $dir eq ""      ) { next; }
     if ( ! defined($dir) ) { next; }
     my $files_ref = $files{$dir};
     my @files     = @$files_ref;
     foreach my $file ( @files ) {
        if ( $file !~ /\n$/ ) { $file = "$file\n"; }
        print OUT  "file = \$DIN_LOC_ROOT/$file";
     }
  }
  close( OUT );
  if ( $printing ) {
     print "\n\nSuccessful\n\n"
  }
