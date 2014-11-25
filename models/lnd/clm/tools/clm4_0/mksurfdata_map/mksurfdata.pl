#!/usr/bin/env perl
#
# Oct/30/2008                                         Erik Kluzek
#
# mksurfdata.pl Perl script to make surface datasets for all resolutions.
#
#
use Cwd;
use strict;
use English;
use IO::File;
use Getopt::Long;


#Figure out where configure directory is and where can use the XML/Lite module from
my $ProgName;
($ProgName = $PROGRAM_NAME) =~ s!(.*)/!!; # name of program
my $ProgDir = $1;                         # name of directory where program lives

my $cwd = getcwd();  # current working directory
my $scrdir;

if ($ProgDir) { $scrdir = $ProgDir; }
else { $scrdir = $cwd; }

#-----------------------------------------------------------------------------------------------
# Add $scrdir to the list of paths that Perl searches for modules
my @dirs = ( $scrdir, "$scrdir/../../../../../../scripts/ccsm_utils/Tools/perl5lib",
           );
unshift @INC, @dirs;
my $result = eval "require XML::Lite";
if ( ! defined($result) ) {
   die <<"EOF";
** Cannot find perl module \"XML/Lite.pm\" from directories: @dirs **
EOF
}
my $result = eval "require Build::NamelistDefinition";
if ( ! defined($result) ) {
   die <<"EOF";
** Cannot find perl module \"Build/NamelistDefinition.pm\" from directories: @dirs **
EOF
}
my $nldef_file     = "$scrdir/../../../bld/namelist_files/namelist_definition_clm4_0.xml";

my $definition = Build::NamelistDefinition->new( $nldef_file );

my $CSMDATA = "/glade/p/cesm/cseg/inputdata";

my %opts = ( 
               hgrid=>"all", 
               rcp=>"-999.9", 
               debug=>0,
               exedir=>undef,
               allownofile=>undef,
               crop=>undef,
               hires=>undef,
               irrig=>undef, 
               years=>"1850,2000",
               glc_nec=>0,
               help=>0,
               mv=>0,
               old_woodharv=>undef,
               irrig=>undef,
               pft_override=>undef,
               pft_frc=>undef,
               pft_idx=>undef,
               soil_override=>undef,
               soil_cly=>undef,
               soil_snd=>undef,
               soil_col=>undef,
               soil_fmx=>undef,
	       outnc_double=>undef,
	       outnc_dims=>"2",     
               usrname=>"",
               usr_mapdir=>"../../shared/mkmapdata",
               dynpft=>undef,
               csmdata=>$CSMDATA,
           );

my $numpft = 16;

#-----------------------------------------------------------------------------------------------
sub usage {
    die <<EOF;
SYNOPSIS 

     For supported resolutions:
     $ProgName -res <res>  [OPTIONS]
        -res [or -r] is the supported resolution(s) to use for files (by default $opts{'hgrid'} ).

      
     For unsupported, user-specified resolutions:	
     $ProgName -res usrspec -usr_gname <user_gname> -usr_gdate <user_gdate>  [OPTIONS]
        -usr_gname "user_gname"    User resolution name to find grid file with 
                                   (only used if -res is set to 'usrspec')
        -usr_gdate "user_gdate"    User map date to find mapping files with
                                   (only used if -res is set to 'usrspec')
                                   NOTE: all mapping files are assumed to be in mkmapdata
                                    - and the user needs to have invoked mkmapdata in 
                                      that directory first
        -usr_mapdir "mapdirectory" Directory where the user-supplied mapping files are
                                   Default: $opts{'usr_mapdir'}

OPTIONS
     -allownofile                  Allow the script to run even if one of the input files
                                   does NOT exist.
     -crop                         Add in crop datasets
     -dinlc [or -l]                Enter the directory location for inputdata 
                                   (default $opts{'csmdata'})
     -debug [or -d]                Do not actually run -- just print out what 
                                   would happen if ran.
     -dynpft "filename"            Dynamic PFT/harvesting file to use 
                                   (rather than create it on the fly) 
                                   (must be consistent with first year)
     -glc_nec "number"             Number of glacier elevation classes to use (by default $opts{'glc_nec'})
     -hires                        If you want to use high-resolution input datasets rather 
                                   than the default lower resolution datasets 
                                   (low resolution is typically at half-degree)
     -irrig                        If you want to include irrigated crop in the output file.
     -exedir "directory"           Directory where mksurfdata_map program is
                                   (by default assume it is in the current directory)
     -mv                           If you want to move the files after creation to the 
                                   correct location in inputdata
                                   (by default -nomv is assumed so files are NOT moved)
     -old_woodharv                 Use the old "bad" wood harvesting (for rcp6 and rcp8.5)
                                   (this is the version that was used for IPCC simulations)
                                   (by default use the new good wood harvest datasets)
     -years [or -y]                Simulation year(s) to run over (by default $opts{'years'}) 
                                   (can also be a simulation year range: i.e. 1850-2000)
     -help  [or -h]                Display this help.
     
     -rcp   [or -c] "rep-con-path" Representative concentration pathway(s) to use for 
                                   future scenarios 
                                   (by default $opts{'rcp'}, where -999.9 means historical ).
     -usrname "clm_usrdat_name"    CLM user data name to find grid file with.

      NOTE: years, res, and rcp can be comma delimited lists.


OPTIONS to override the mapping of the input gridded data with hardcoded input

     -pft_frc "list of fractions"  Comma delimited list of percentages for veg types
     -pft_idx "list of veg index"  Comma delimited veg index for each fraction
     -soil_cly "% of clay"         % of soil that is clay
     -soil_col "soil color"        Soil color (1 [light] to 20 [dark])
     -soil_fmx "soil fmax"         Soil maximum saturated fraction (0-1)
     -soil_snd "% of sand"         % of soil that is sand

EOF
}

sub check_soil {
#
# check that the soil options are set correctly
#
  foreach my $type ( "soil_cly", "soil_snd" ) {
     if ( ! defined($opts{$type} )  ) {
        die "ERROR: Soil variables were set, but $type was NOT set\n";
     }
  }
  #if ( $opts{'soil_col'} < 0 || $opts{'soil_col'} > 20 ) {
  #   die "ERROR: Soil color is out of range = ".$opts{'soil_col'}."\n";
  #}
  my $texsum = $opts{'soil_cly'} + $opts{'soil_snd'};
  my $loam   = 100.0 - $texsum;
  if ( $texsum < 0.0 || $texsum > 100.0 ) {
     die "ERROR: Soil textures are out of range: clay = ".$opts{'soil_cly'}.
         " sand = ".$opts{'soil_snd'}." loam = $loam\n";
  }
}

sub check_soil_col_fmx {
#
# check that the soil color or soil fmax option is set correctly
#
  if ( defined($opts{'soil_col'}) ) {
     if ( $opts{'soil_col'} < 0 || $opts{'soil_col'} > 20 ) {
        die "ERROR: Soil color is out of range = ".$opts{'soil_col'}."\n";
     }
  }
  if ( defined($opts{'soil_fmx'}) ) {
     if ( $opts{'soil_fmx'} < 0.0 || $opts{'soil_fmx'} > 1.0 ) {
        die "ERROR: Soil fmax is out of range = ".$opts{'soil_fmx'}."\n";
     }
  }
}

sub check_pft {
#
# check that the pft options are set correctly
#
  # Eliminate starting and ending square brackets
  $opts{'pft_idx'} =~ s/^\[//;
  $opts{'pft_idx'} =~ s/\]$//;
  $opts{'pft_frc'} =~ s/^\[//;
  $opts{'pft_frc'} =~ s/\]$//;
  foreach my $type ( "pft_idx", "pft_frc" ) {
     if ( ! defined($opts{$type} ) ) {
        die "ERROR: PFT variables were set, but $type was NOT set\n";
     }
  }
  my @pft_idx     = split( /,/, $opts{'pft_idx'} );
  my @pft_frc     = split( /,/, $opts{'pft_frc'} );
  if ( $#pft_idx != $#pft_frc ) {
     die "ERROR: PFT arrays are different sizes: pft_idx and pft_frc\n";
  }
  my $sumfrc = 0.0;
  for( my $i = 0; $i <= $#pft_idx; $i++ ) {
     # check index in range
     if ( $pft_idx[$i] < 0 || $pft_idx[$i] > $numpft ) {
         die "ERROR: pft_idx out of range = ".$opts{'pft_idx'}."\n";
     }
     # make sure there are no duplicates
     for( my $j = 0; $j < $i; $j++ ) {
        if ( $pft_idx[$i] == $pft_idx[$j] ) {
            die "ERROR: pft_idx has duplicates = ".$opts{'pft_idx'}."\n";
        }
     }
     # check fraction in range
     if ( $pft_frc[$i] <= 0.0 || $pft_frc[$i] > 100.0 ) {
         die "ERROR: pft_frc out of range (>0.0 and <=100.0) = ".$opts{'pft_frc'}."\n";
     }
     $sumfrc = $sumfrc + $pft_frc[$i];
  }
  # check that fraction sums up to 100%
  if ( abs( $sumfrc - 100.0) > 1.e-6 ) {
      die "ERROR: pft_frc does NOT add up to 100% = ".$opts{'pft_frc'}."\n";
  }
  
}
 
# Perl trim function to remove whitespace from the start and end of the string
sub trim($)
{
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}

#-----------------------------------------------------------------------------------------------

   my $cmdline = "@ARGV";
   GetOptions(
        "allownofile"  => \$opts{'allownofile'},
        "r|res=s"      => \$opts{'hgrid'},
        "usr_gname=s"  => \$opts{'usr_gname'},
        "usr_gdate=s"  => \$opts{'usr_gdate'},
        "usr_mapdir=s" => \$opts{'usr_mapdir'},
        "crop"         => \$opts{'crop'},
        "irrig"        => \$opts{'irrig'},
        "hires"        => \$opts{'hires'},
        "c|rcp=s"      => \$opts{'rcp'},
        "l|dinlc=s"    => \$opts{'csmdata'},
        "d|debug"      => \$opts{'debug'},
        "dynpft=s"     => \$opts{'dynpft'},
        "y|years=s"    => \$opts{'years'},
        "old_woodharv" => \$opts{'old_woodharv'},
        "exedir=s"     => \$opts{'exedir'},
        "h|help"       => \$opts{'help'},
        "usrname=s"    => \$opts{'usrname'},
        "glc_nec=i"    => \$opts{'glc_nec'},
        "irrig"        => \$opts{'irrig'},
        "mv"           => \$opts{'mv'},
        "pft_frc=s"    => \$opts{'pft_frc'},
        "pft_idx=s"    => \$opts{'pft_idx'},
        "soil_col=i"   => \$opts{'soil_col'},
        "soil_fmx=f"   => \$opts{'soil_fmx'},
        "soil_cly=f"   => \$opts{'soil_cly'},
        "soil_snd=f"   => \$opts{'soil_snd'},
   ) or usage();

   # Check for unparsed arguments
   if (@ARGV) {
       print "ERROR: unrecognized arguments: @ARGV\n";
       usage();
   }
   if ( $opts{'help'} ) {
       usage();
   }
   # If csmdata was changed from the default
   if ( $CSMDATA ne $opts{'csmdata'} ) {
      $CSMDATA = $opts{'csmdata'};
   }
   my $glc_nec = $opts{'glc_nec'};
   #
   # Set disk location to send files to, and list resolutions to operate over, 
   # set filenames, and short-date-name
   #
   my @hresols;
   my $mapdate; 
   if ( $opts{'hgrid'} eq "all" ) {
      my @all_hresols = $definition->get_valid_values( "res" );
      @hresols = @all_hresols;
   } elsif ( $opts{'hgrid'} eq "usrspec" ) {
      @hresols = $opts{'usr_gname'}; 
      $mapdate = $opts{'usr_gdate'}; 
   } else {
      @hresols = split( ",", $opts{'hgrid'} );
      # Check that resolutions are valid
      foreach my $res ( @hresols ) {
         if ( ! $definition->is_valid_value( "res", "'$res'" ) ) {
            if ( $opts{'usrname'} eq ""  || $res ne $opts{'usrname'} ) {
               print "** Invalid resolution: $res\n";
               usage();
            }
         }
      }
   }
   #
   # Set years to run over
   #
   my @years   = split( ",", $opts{'years'} );
   # Check that resolutions are valid
   foreach my $sim_year ( @years ) {
      if ( ! $definition->is_valid_value( "sim_year", $sim_year ) ) {
         if ( ! $definition->is_valid_value( "sim_year_range", "'$sim_year'" ) ) {
            print "** Invalid simulation year or simulation year range: $sim_year\n";
            usage();
         }
      }
   }
   #
   # Set rcp to use
   #
   my @rcpaths = split( ",", $opts{'rcp'} );
   # Check that rcp is valid
   foreach my $rcp ( @rcpaths  ) {
      if ( ! $definition->is_valid_value( "rcp", $rcp ) ) {
         if ( ! $definition->is_valid_value( "rcp", "$rcp" ) ) {
            print "** Invalid rcp: $rcp\n";
            usage();
         }
      }
   }
   # Check if soil set
   if ( defined($opts{'soil_cly'}) || 
        defined($opts{'soil_snd'}) ) {
       &check_soil( );
       $opts{'soil_override'} = 1;
   }
   # Check if pft set
   if ( defined($opts{'crop'}) ) { $numpft = 20; }   # First set numpft if crop is on
   if ( defined($opts{'pft_frc'}) || defined($opts{'pft_idx'}) ) {
       &check_pft( );
       $opts{'pft_override'} = 1;
   }
   # Check if dynpft set and is valid filename
   if ( defined($opts{'dynpft'}) ) {
       if ( ! -f $opts{'dynpft'} ) {
          print "** Dynamic PFT file does NOT exist: $opts{'dynpft'}\n";
          usage();
       }
   }

   my $nl = "namelist";
   my $sdate = "c" . `date +%y%m%d`;
   chomp( $sdate );

   my @ncfiles;
   my @lfiles;
   my @tsfiles;
   my $cfile = "clm.input_data_files";
   if ( -f "$cfile" ) {
      `/bin/mv -f $cfile ${cfile}.previous`;
   }
   my $cfh = IO::File->new;
   $cfh->open( ">$cfile" ) or die "** can't open file: $cfile\n";
   system( "\rm -f $cfile" );
   system( "touch $cfile" );
   print $cfh <<"EOF";
#! /bin/csh -f
set CSMDATA = $CSMDATA
EOF
   system( "chmod +x $cfile" );
   my $svnrepo = "https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata";
   my $svnmesg = "Update fsurdat files with mksurfdata_map";
   my $surfdir = "lnd/clm2/surfdata";

   #
   # Loop over all resolutions listed
   #
   foreach my $res ( @hresols ) {
      #
      # Query the XML default file database to get the appropriate files
      #
      my $queryopts, my $queryfilopts; 
      if ( $opts{'hgrid'} eq "usrspec" ) {
	  $queryopts = "-csmdata $CSMDATA -silent -justvalue -phys clm4_0";
      } else {
	  $queryopts = "-res $res -csmdata $CSMDATA -silent -justvalue -phys clm4_0";
      }
      $queryfilopts = "$queryopts -onlyfiles ";
      my $mkcrop = "";
      my $setnumpft = "";
      if ( defined($opts{'crop'}) ) {
         $mkcrop    = ",crop='on'";
         $setnumpft = "numpft = $numpft"
      }
      my $usrnam    = "";
      if ( $opts{'usrname'} ne "" && $res eq $opts{'usrname'} ) {
         $usrnam    = "-usrname ".$opts{'usrname'};
      }
      #
      # Mapping files
      #
      my %map; my %hgrd; my %lmsk; my %datfil;
      my $hires = "off";
      if ( defined($opts{'hires'}) ) {
         $hires = "on";
      }
      my $mopts  = "$queryopts -namelist default_settings $usrnam";
      my $mkopts = "-csmdata $CSMDATA -silent -justvalue -namelist clmexp $usrnam -phys clm4_0";
      foreach my $typ ( "lak", "veg", "voc", "top", "irr", "tex", "col", 
                        "fmx", "lai", "urb", "org", "glc", "utp", "wet" ) {
         my $lmask = `$scrdir/../../../bld/queryDefaultNamelist.pl $mopts -options type=$typ,glc_nec=$glc_nec,hires=$hires -var lmask`;
         $lmask = trim($lmask);
         my $hgrid = `$scrdir/../../../bld/queryDefaultNamelist.pl $mopts -options type=$typ,glc_nec=$glc_nec,hires=$hires -var hgrid`;
         $hgrid = trim($hgrid);
         my $filnm = `$scrdir/../../../bld/queryDefaultNamelist.pl $mopts -options type=$typ -var mksrf_filename`;
         $filnm = trim($filnm);
         $hgrd{$typ} = $hgrid;
         $lmsk{$typ} = $lmask;
	 if ( $opts{'hgrid'} eq "usrspec" ) {
	     $map{$typ} = $opts{'usr_mapdir'}."/map_${hgrid}_${lmask}_to_${res}_nomask_aave_da_c${mapdate}\.nc";
	 } else {
	     $map{$typ} = `$scrdir/../../../bld/queryDefaultNamelist.pl $queryfilopts -namelist clmexp -options frm_hgrid=$hgrid,frm_lmask=$lmask,to_hgrid=$res,to_lmask=nomask -var map`;
	 }	     
         $map{$typ} = trim($map{$typ});
         if ( $map{$typ} !~ /[^ ]+/ ) {
            die "ERROR: could NOT find a mapping file for this resolution: $res and type: $typ at $hgrid and $lmask.\n";
         }
         if ( ! defined($opts{'allownofile'}) && ! -f $map{$typ} ) {
            die "ERROR: mapping file for this resolution does NOT exist ($map{$typ}).\n";
         }
         $datfil{$typ} = `$scrdir/../../../bld/queryDefaultNamelist.pl $mkopts -options hgrid=$hgrid,lmask=$lmask,glc_nec=$glc_nec$mkcrop -var $filnm`;
         $datfil{$typ} = trim($datfil{$typ});
         if ( $datfil{$typ} !~ /[^ ]+/ ) {
            die "ERROR: could NOT find a $filnm data file for this resolution: $res and type: $typ at $hgrid and $lmask.\n";
         }
         if ( ! defined($opts{'allownofile'}) && ! -f $datfil{$typ} ) {
            die "ERROR: data file for this resolution does NOT exist ($datfil{$typ}).\n";
         }
         if ( $typ eq "irr" && ! defined($opts{'irrig'}) ) {
            $map{$typ}    = " ";
            $datfil{$typ} = " ";
         }
      }
      #
      # Grid file from the pft map file or grid if not found
      #
      my $griddata    = trim($map{'veg'});
      if ( $griddata eq "" ) {
         $griddata = `$scrdir/../../../bld/queryDefaultNamelist.pl $queryfilopts $usrnam -var fatmgrid`;
         if ( $griddata eq "" ) {
            die "ERROR: could NOT find a grid data file for this resolution: $res.\n";
         }
      }
      my $desc;
      my $desc_yr0;
      #
      # Check if all urban single point dataset
      #
      my @all_urb = ( "1x1_camdenNJ","1x1_vancouverCAN", "1x1_mexicocityMEX", 
                      "1x1_asphaltjungleNJ", "1x1_urbanc_alpha" );
      my $all_urb = ".false.";
      my $urb_pt  = 0;
      foreach my $urb_res ( @all_urb ) {
         if ( $res eq $urb_res ) {
            $all_urb = ".true.";
            $urb_pt  = 1;
         }
      }
      #
      # Always run at double precision for output
      #
      my $double = ".true.";
      #
      # Loop over each sim_year
      #
      RCP: foreach my $rcp ( @rcpaths ) {
         #
         # Loop over each sim_year
         #
         SIM_YEAR: foreach my $sim_year ( @years ) {
            #
            # Skip if urban unless sim_year=2000
            #
            if ( $urb_pt && $sim_year != 2000 ) {
               print "For urban -- skip this simulation year = $sim_year\n";
               next SIM_YEAR;
            }
            #
            # If year is 1850-2000 actually run 1850-2005
            #
            if ( $sim_year eq "1850-2000" ) {
               my $actual = "1850-2005";
               print "For $sim_year actually run $actual\n";
               $sim_year = $actual;
            }
            #
            # Irrigation dataset
            #
            my $irrdes = "";
            if ( defined($opts{'irrig'}) ) {
              $irrdes = "irrcr_";
            }
            #
            # Create namelist file
            #
            my $fh = IO::File->new;
            $fh->open( ">$nl" ) or die "** can't open file: $nl\n";
	    print "CSMDATA is $CSMDATA \n";
            print $fh <<"EOF";
&clmexp
 nglcec         = $glc_nec
 mksrf_fgrid    = '$griddata'
 map_fpft       = '$map{'veg'}'
 map_fglacier   = '$map{'glc'}'
 map_fsoicol    = '$map{'col'}'
 map_furban     = '$map{'urb'}'
 map_fmax       = '$map{'fmx'}'
 map_forganic   = '$map{'org'}'
 map_flai       = '$map{'lai'}'
 map_fharvest   = '$map{'lai'}'
 map_flakwat    = '$map{'lak'}'
 map_fwetlnd    = '$map{'wet'}'
 map_fvocef     = '$map{'voc'}'
 map_fsoitex    = '$map{'tex'}'
 map_firrig     = '$map{'irr'}'
 map_furbtopo   = '$map{'utp'}'
 map_flndtopo   = '$map{'top'}'
 mksrf_fsoitex  = '$datfil{'tex'}'
 mksrf_forganic = '$datfil{'org'}'
 mksrf_flakwat  = '$datfil{'lak'}'
 mksrf_fwetlnd  = '$datfil{'wet'}'
 mksrf_fmax     = '$datfil{'fmx'}'
 mksrf_fglacier = '$datfil{'glc'}'
 mksrf_fvocef   = '$datfil{'voc'}'
 mksrf_furbtopo = '$datfil{'utp'}'
 mksrf_flndtopo = '$datfil{'top'}'
 mksrf_firrig   = '$datfil{'irr'}'
 outnc_double   = $double
 all_urban      = $all_urb
EOF
            my $urbdesc = "urb3den";
            if ( ! $urb_pt ) {
               print $fh <<"EOF";
 mksrf_furban   = '$datfil{'urb'}'
EOF
            } else {
               my $urbdata = `$scrdir/../../../bld/queryDefaultNamelist.pl $queryfilopts -var fsurdat -filenameonly`;
               if ( $? != 0 ) {
                  die "ERROR:: furbinp file NOT found\n";
               }
               chomp( $urbdata );
               print $fh <<"EOF";
 mksrf_furban   = '$CSMDATA/lnd/clm2/surfdata/$urbdata'
EOF
            }

            my $resol = "";
            if ( $res ne "1x1_tropicAtl" ) {
               $resol = "-res $hgrd{'veg'}";
            }
            my $sim_yr0 = $sim_year;
            my $sim_yrn = $sim_year;
            if ( $sim_year =~ /([0-9]+)-([0-9]+)/ ) {
               $sim_yr0 = $1;
               $sim_yrn = $2;
            }
            my $cmd    = "$scrdir/../../../bld/queryDefaultNamelist.pl $queryfilopts $resol -options sim_year=${sim_yr0}$mkcrop -var mksrf_fvegtyp -namelist clmexp";
            my $vegtyp = `$cmd`;
            chomp( $vegtyp );
            if ( $vegtyp eq "" ) {
               die "** trouble getting vegtyp file with: $cmd\n";
            }
            if ( $rcp == -999.9 ) {
               $desc     = sprintf( "hist_simyr%4.4d-%4.4d", $sim_yr0, $sim_yrn );
               $desc_yr0 = sprintf( "simyr%4.4d",            $sim_yr0  );
            } else {
               $desc     = sprintf( "%s%2.1f_simyr%4.4d-%4.4d", "rcp", $rcp, $sim_yr0, $sim_yrn );
               $desc_yr0 = sprintf( "%s%2.1f_simyr%4.4d",       "rcp", $rcp, $sim_yr0  );
            }
            my $strlen = 195;
            my $dynpft_format = "%-${strlen}.${strlen}s %4.4d\n";
            my $options = "";
            my $crpdes  = "";
            if ( $mkcrop ne "" ) { 
               $options = "-options $mkcrop";
               $crpdes  = "mp20_";
            }
            my $landuse_timeseries_text_file;
	    if ( $sim_year ne $sim_yr0 ) {
		if ( ! defined($opts{'dynpft'}) && ! $opts{'pft_override'} ) {
		    $landuse_timeseries_text_file = "landuse_timeseries_$desc.txt";
		    my $fh_landuse_timeseries = IO::File->new;
		    $fh_landuse_timeseries->open( ">$landuse_timeseries_text_file" ) or die "** can't open file: $landuse_timeseries_text_file\n";
		    print "Writing out landuse_timeseries text file: $landuse_timeseries_text_file\n";

                    # use new good wood harvest by default:
                    my $ngwh = ",ngwh=on";
                    if ( defined($opts{'old_woodharv'}) ) {
                       $ngwh = "";
                    }

		    for( my $yr = $sim_yr0; $yr <= $sim_yrn; $yr++ ) {
                        my $vegtypyr = `$scrdir/../../../bld/queryDefaultNamelist.pl $queryfilopts $resol -options sim_year=$yr,rcp=${rcp}${mkcrop}${ngwh} -var mksrf_fvegtyp -namelist clmexp`;
			chomp( $vegtypyr );
			printf $fh_landuse_timeseries $dynpft_format, $vegtypyr, $yr;
			if ( $yr % 100 == 0 ) {
			    print "year: $yr\n";
			}
		    }
		    $fh_landuse_timeseries->close;
		    print "Done writing file\n";
		} elsif ( $opts{'pft_override'} && defined($opts{'dynpft'}) ) {
		    $landuse_timeseries_text_file = $opts{'dynpft'};
		} else {
		    $landuse_timeseries_text_file = "landuse_timeseries_override_$desc.txt";
		    my $fh_landuse_timeseries = IO::File->new;
		    $fh_landuse_timeseries->open( ">$landuse_timeseries_text_file" ) or die "** can't open file: $landuse_timeseries_text_file\n";
		    my $frstpft = "<pft_f>$opts{'pft_frc'}</pft_f>" . 
			"<pft_i>$opts{'pft_idx'}</pft_i>" .
			"<harv>0,0,0,0,0</harv><graz>0</graz>";
		    print "Writing out landuse_timeseries text file: $landuse_timeseries_text_file\n";
		    if ( (my $len = length($frstpft)) > $strlen ) {
			die "ERROR PFT line is too long ($len): $frstpft\n";
		    }
		    printf $fh_landuse_timeseries $dynpft_format, $frstpft, $sim_yr0;
		    $fh_landuse_timeseries->close;
		    print "Done writing file\n";
		}
	    }

            if ( defined($opts{'soil_override'}) ) {
               print $fh <<"EOF";
 soil_clay     = $opts{'soil_cly'}
 soil_sand     = $opts{'soil_snd'}
EOF
            }
            if ( defined($opts{'pft_override'}) ) {
               print $fh <<"EOF";
 pft_frc      = $opts{'pft_frc'}
 pft_idx      = $opts{'pft_idx'}
EOF
            }

            print $fh <<"EOF";
 mksrf_fvegtyp  = '$vegtyp'
 mksrf_fsoicol  = '$datfil{'col'}'
 mksrf_flai     = '$datfil{'lai'}'
EOF
	    my $ofile = "surfdata_${res}_${desc_yr0}_${sdate}";
	    print $fh <<"EOF";
 fsurdat        = '$ofile.nc'
 fsurlog        = '$ofile.log'     
EOF

	    my $ofile_ts = "landuse.timeseries_${res}_${desc}_${sdate}";
	    if ( $sim_year ne $sim_yr0 ) {
	        print $fh <<"EOF";
 mksrf_fdynuse  = '$landuse_timeseries_text_file'
 fdyndat        = '$ofile_ts.nc'
EOF
            } else {
		print $fh <<"EOF";
 mksrf_fdynuse  = ' '
 fdyndat        = ' '
EOF
            }
            print $fh <<"EOF";
 $setnumpft
/
EOF
            $fh->close;
	    print "resolution: $res rcp=$rcp sim_year = $sim_year\n";
 	    print "namelist: $nl\n";
	    # 
	    # Print namelist file 
            $fh->open( "<$nl" ) or die "** can't open file: $nl\n";
            while( $_ = <$fh> ) {
              print $_;
            }
            $fh->close;
            #
            # Delete previous versions of files that will be created
            #
            system( "/bin/rm -f $ofile.nc $ofile.log" );
            #
            # Run mksurfdata_map with the namelist file
            #
            my $exedir = $scrdir;
            if ( defined($opts{'exedir'}) ) {
               $exedir = $opts{'exedir'};
            }
            print "$exedir/mksurfdata_map < $nl\n";
            my $filehead;
            my $tsfilehead;
            if ( ! $opts{'debug'} ) {
               system( "$exedir/mksurfdata_map < $nl" );
               if ( $? ) { die "ERROR in mksurfdata_map: $?\n"; }
            } else {
               $filehead  = "surfdata_$res";
               $tsfilehead = "landuse.timeseries_testfile";
               system( "touch $filehead.nc" );
               system( "touch $tsfilehead.nc" );
               system( "touch $filehead.log" );
            }
	    print "\n===========================================\n\n";
            #
            # Check that files were created
            #
            @ncfiles  = glob( "$ofile.nc" );
            if ( $#ncfiles != 0 ) {
              die "ERROR surfdata netcdf file was NOT created!\n";
            }
            chomp( $ncfiles[0] );
            @lfiles = glob( "$ofile.log" );
            chomp( $lfiles[0] );
            @tsfiles = glob( "$ofile_ts.nc" );
            chomp( $tsfiles[0] );
            if ( $#tsfiles != 0 ) {
              die "ERROR surfdata landuse_timeseries netcdf file was NOT created!\n";
            }
            #
            # If urban point, append grid and frac file on top of surface dataset
            #
            if ( $urb_pt ) {
               my $cmd = "ncks -A $griddata $ncfiles[0]";
               print "$cmd\n";
               if ( ! $opts{'debug'} ) { system( $cmd ); }
               my $fracdata = `$scrdir/../../../bld/queryDefaultNamelist.pl $queryopts -var fatmlndfrc`;
               if ( $? != 0 ) {
                  die "ERROR:: fatmlndfrc file NOT found\n";
               }
               chomp( $fracdata );
               $cmd = "ncks -A $fracdata $ncfiles[0]";
               print "$cmd\n";
               if ( ! $opts{'debug'} ) { system( $cmd ); }
            }
            #
            # Rename files to CSMDATA
            #
            my $lsvnmesg = "'$svnmesg $urbdesc $desc'";
            if ( -f "$ncfiles[0]" && -f "$lfiles[0]" ) {
               my $outdir = "$CSMDATA/$surfdir";
               my $ofile = "surfdata_${res}_${crpdes}${desc_yr0}_${irrdes}${sdate}";
               my $mvcmd = "/bin/mv -f $ncfiles[0]  $outdir/$ofile.nc";
               if ( ! $opts{'debug'} && $opts{'mv'} ) {
                  print "$mvcmd\n";
                  system( "$mvcmd" );
                  chmod( 0444, "$outdir/$ofile.nc" );
               }
               my $mvcmd = "/bin/mv -f $lfiles[0] $outdir/$ofile.log";
               if ( ! $opts{'debug'} && $opts{'mv'} ) {
                  print "$mvcmd\n";
                  system( "$mvcmd" );
                  chmod( 0444, "$outdir/$ofile.log" );
               }
               if ( $opts{'mv'} ) {
                  print $cfh "# FILE = \$DIN_LOC_ROOT/$surfdir/$ofile.nc\n";
                  print $cfh "svn import -m $lsvnmesg \$CSMDATA/$surfdir/$ofile.nc " . 
                             "$svnrepo/$surfdir/$ofile.nc\n";
                  print $cfh "# FILE = \$DIN_LOC_ROOT/$surfdir/$ofile.log\n";
                  print $cfh "svn import -m $lsvnmesg \$CSMDATA/$surfdir/$ofile.log " .
                             "$svnrepo/$surfdir/$ofile.log\n";
               }
               # If running a transient case
               if ( $sim_year ne $sim_yr0 ) {
                  $ofile = "landuse.timeseries_${res}_${desc}_${sdate}";
                  $mvcmd = "/bin/mv -f $tsfiles[0] $outdir/$ofile.nc";
                  if ( ! $opts{'debug'} && $opts{'mv'} ) {
                     print "$mvcmd\n";
                     system( "$mvcmd" );
                     chmod( 0444, "$outdir/$ofile.nc" );
                  }
                  if ( $opts{'mv'} ) {
                     print $cfh "# FILE = \$DIN_LOC_ROOT/$surfdir/$ofile.nc\n";
                     print $cfh "svn import -m $lsvnmesg \$CSMDATA/$surfdir/$ofile.nc " .
                                "$svnrepo/$surfdir/$ofile.nc\n";
                  }
               }
   
            } elsif( ! $opts{'debug'} ) {
              die "ERROR files were NOT created: nc=$ncfiles[0] log=$lfiles[0]\n";
            }
            if ( (! $opts{'debug'}) && $opts{'mv'} && (-f "$ncfiles[0]" || -f "$lfiles[0]") ) {
              die "ERROR files were NOT moved: nc=$ncfiles[0] log=$lfiles[0]\n";
            }
            if ( ! $opts{'debug'} ) {
               system( "/bin/rm -f $filehead.nc $filehead.log $tsfilehead.nc" );
            }
         } # End of sim_year loop
      }    # End of rcp loop
   }
   close( $cfh );
   print "Successfully created fsurdat files\n";
