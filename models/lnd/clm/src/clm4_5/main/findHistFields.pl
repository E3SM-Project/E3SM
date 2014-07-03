#!/usr/bin/env perl
#
# This perl script reads in the histFldsMod.F90 file to find the total list of history 
# fields that can be added for this model version, regardless of namelist options, or
# CPP processing.
# 
use strict;
#use warnings;
#use diagnostics;

use Cwd;
use English;
use Getopt::Long;
use IO::File;
use File::Glob ':glob';

# Set the directory that contains the CLM configuration scripts.  If the command was
# issued using a relative or absolute path, that path is in $ProgDir.  Otherwise assume
# the
# command was issued from the current working directory.

(my $ProgName = $0) =~ s!(.*)/!!;      # name of this script
my $ProgDir = $1;                      # name of directory containing this script -- may be a
                                       # relative or absolute path, or null if the script
                                       # is in
                                       # the user's PATH
my $cmdline = "@ARGV";                 # Command line arguments to script
my $cwd = getcwd();                    # current working directory
my $cfgdir;                            # absolute pathname of directory that contains this script
my $nm = "${ProgName}::";              # name to use if script dies
if ($ProgDir) { 
    $cfgdir = $ProgDir;
} else {
    $cfgdir = $cwd;
}
# The namelist definition file contains entries for all namelist variables that
# can be output by build-namelist.
my $nl_definition_file = "$cfgdir/../../../bld/namelist_files/namelist_definition_clm4_5.xml";
(-f "$nl_definition_file")  or  die <<"EOF";
** $ProgName - Cannot find namelist definition file \"$nl_definition_file\" **
EOF
print "Using namelist definition file $nl_definition_file\n";

# The Build::NamelistDefinition module provides utilities to get the list of
# megan compounds

#The root directory to cesm utils Tools
my $cesm_tools = "$cfgdir/../../../../../../scripts/ccsm_utils/Tools";

(-f "$cesm_tools/perl5lib/Build/NamelistDefinition.pm")  or  die <<"EOF";
** $ProgName - Cannot find perl module \"Build/NamelistDefinition.pm\" in directory 
    \"$cesm_tools/perl5lib\" **
EOF
# Add $cfgdir/perl5lib to the list of paths that Perl searches for modules
my @dirs = ( $cfgdir, "$cesm_tools/perl5lib");
unshift @INC, @dirs;
require Build::NamelistDefinition;
# Create a namelist definition object.  This object provides a method for verifying that
# the
# output namelist variables are in the definition file, and are output in the correct
# namelist groups.
my $definition = Build::NamelistDefinition->new($nl_definition_file);


my $mxname  = 0;
my $mxlongn = 0;
my %fields;
my $fldnamevar = "fieldname_var";

sub matchKeyword {
#
# Match a keyword
#
  my $keyword = shift;
  my $line    = shift;
  my $fh      = shift;

  my $match = undef;
  if ( $line =~ /$keyword/ ) {
     if ( $line =~ /$keyword\s*=\s*['"]([^'"]+)['"]/ ) {
        $match = $1;
     } elsif ( $line =~ /$keyword\s*=\s*&\s*$/ ) {
        $line  = <$fh>;
        if ( $line =~ /^\s*['"]([^'"]+)['"]/ ) {
           $match = $1;
        } else {
           die "ERROR: Trouble getting keyword string\n Line: $line";
        }
     } else {
        if (      $line =~ /fname\s*=\s*fieldname/ ) {
           print STDERR "Found variable used for fieldname = $line\n";
           $match = $fldnamevar;
        } elsif ( $line =~ /fname\s*=\s*trim\(fname\)/ ) {
           $match = undef;
        } elsif ( $line =~ /units\s*=\s*units/ ) {
           $match = undef;
        } elsif ( $line =~ /long_name\s*=\s*long_name/ ) {
           $match = undef;
        } elsif ( $line =~ /long_name\s*=\s*longname/ ) {
           print STDERR "Found variable used for longname = $line\n";
           $match = "longname_var";
        } else {
          die "ERROR: Still have a match on $keyword\n Line: $line";
        }
     }
  }
  return( $match );
}

sub getFieldInfo {
#
# Get field Information
#
  my $fh   = shift;
  my $line = shift;

  my $fname = undef;
  my $units = undef;
  my $longn = undef;
  my $endin = undef;
  do {
    if ( $line =~ /MEG_/ ) {
       $line =~ s|'//'_'|_'|g;
       $line =~ s|'//trim\(meg_cmp\%name\)|megancmpd'|gi;
       if ( $line =~ /meg_cmp\%name/ ) {
          die "ERROR: Still have meg_cmp in a line\n";
       }
    }
    if ( ! defined($fname) ) {
       $fname = &matchKeyword( "fname",     $line, $fh );
    }
    if ( ! defined($units) ) {
       $units = &matchKeyword( "units",     $line, $fh );
    }
    if ( ! defined($longn) ) {
       $longn = &matchKeyword( "long_name", $line, $fh );
    }
    if ( $line =~ /\)\s*$/ ) {
       $endin = 1;
    }
    if ( ! defined($endin) ) { $line = <$fh>; }

  } until( (defined($fname) && defined($units) && defined($longn)) ||
           ! defined($line) || defined($endin) );
  if ( ! defined($fname) ) {
     die "ERROR: name undefined for field ending with: $line\n";
  }
  return( $fname, $longn, $units );
}

sub setField {
#
# Set the field
#
  my $name  = shift;
  my $longn = shift;
  my $units = shift;

  if ( defined($name) && $name ne $fldnamevar ) {
    if ( length($name)  > $mxname  ) { $mxname  = length($name);  }
    if ( length($longn) > $mxlongn ) { $mxlongn = length($longn); }
    my $len;
    if ( length($longn) > 90 ) {
       $len = 110;
    } elsif ( length($longn) > 60 ) {
       $len = 90;
    } else {
       $len = 60;
    }
    $fields{$name}{'field'} = sprintf( "%-${len}s\t(%s)", $longn, $units );
    $fields{$name}{'longn'} = $longn;
    $fields{$name}{'units'} = $units;
  }
}

sub XML_Header {
#
# Write out header to history fields file
#
  my $outfh       = shift;
  my $outfilename = shift;
  my $filename    = shift;

  print STDERR " Write out header to history fields file to: $outfilename\n";
  my $svnurl = '$URL: https://svn-ccsm-models.cgd.ucar.edu/clm2/trunk_tags/clm4_0_40/models/lnd/clm/src/main/findHistFields.pl $';
  my $svnid  = '$Id: findHistFields.pl 34757 2012-02-15 18:38:05Z erik $';
  print $outfh <<"EOF";
<?xml version="1.0"?>

\<\?xml-stylesheet type="text\/xsl" href="history_fields.xsl"\?\>

\<\!--
  List of history file field names, long-names and units for all the fields output
  by CLM. This was created by reading in the file: $filename
  SVN version information:
  $svnurl
  $svnid
--\>

\<history_fields\>
EOF
}

sub XML_Footer {
#
# Write out footer to history fields file
#
  my $outfh = shift;

  print STDERR " Write out footer to history fields file\n";
  print $outfh "\n</history_fields>\n";
}

my $pwd = `pwd`;
chomp( $pwd );
my @megcmpds  =  $definition->get_valid_values( "megan_cmpds", 'noquotes'=>1 );
my @filenames = ( "$pwd/histFldsMod.F90", "$pwd/../biogeochem/CNFireMod.F90" );

#
# Loop over all files that have hist_addfld calls in them
#
foreach my $filename ( @filenames ) {

   my $fh = IO::File->new($filename, '<') or die "** $ProgName - can't open history Fields file: $filename\n";
   #
   # Read in the list of fields from the source file
   #
   while (my $line = <$fh>) {
   
      # Comments
      if ($line =~ /(.*)\!/) {
        $line = $1;
      }
      if ($line =~ /end subroutine/) {
        last;
      }
      my $format = "\n<field name='%s' units='%s'\n long_name='%s'\n/>\n";
      if ($line =~ /call\s*hist_addfld/i ) {
         (my $name, my $longn, my $units) = &getFieldInfo( $fh, $line );
         if ( $name ne "MEG_megancmpd" ) {
            &setField( $name, $longn, $units );
            printf( <STDERR>, $format, $name, $units, $longn );
         } else {
            foreach my $megcmpd ( @megcmpds ) {
               my $name = "MEG_${megcmpd}";
               &setField( $name, $longn, $units );
               printf( <STDERR>, $format, $name, $units, $longn );
            }
         }
      }
   }
   close( $fh );
}
print STDERR " mxname  = $mxname\n";
print STDERR " mxlongn = $mxlongn\n";
my %pool_name = ( 
                     L1=> { hist=>'LITR1',      long=>'litter 1'            },
                     L2=> { hist=>'LITR2',      long=>'litter 2'            },
                     L3=> { hist=>'LITR3',      long=>'litter 3'            },
                     CWD=>{ hist=>'CWD',        long=>'coarse woody debris' },
                     S1=> { hist=>'SOIL1',      long=>'soil 1'              },
                     S2=> { hist=>'SOIL2',      long=>'soil 2'              },
                     S3=> { hist=>'SOIL3',      long=>'soil 3'              },
                     S4=> { hist=>'SOIL4',      long=>'soil 4'              },
                     atm=>{ hist=>'atmosphere', long=>'atmosphere'          },
                );

my %vrt_suffix = ( C=>" C", "C_vr"=>" C (vertically resolved)", C_1m=>" C to 1 meter", 
                   C_30cm=>" C to 30 cm", C_activelayer=>" C in active layer", 
                   N=>" C", "N_vr"=>" N (vertically resolved)", N_1m=>" N to 1 meter",
                   N_30cm=>" N to 30 cm", N_activelayer=>" N in active layer",
                );
my %firelist = (
                   C_TO_FIRE=>" C fire loss", C_TO_FIRE_vr=>" C fire loss", 
                   N_TO_FIRE=>" N fire loss", N_TO_FIRE_vr=>" N fire loss", 
               );
my %leechlist = (
                   C_TO_LEACHING=>" C leaching loss", C_TNDNCY_VERT_TRANSPORT=>" C tendency due to vertical transport",
                   N_TO_LEACHING=>" N leaching loss", N_TNDNCY_VERT_TRANSPORT=>" N tendency due to vertical transport",
               );
#
# Add fields that are looped over
#
my $name, my $longn, my $units;
foreach my $pool  ( keys(%pool_name) ) {
   my $fname = $pool_name{$pool}{'hist'};
   foreach my $fld ( keys(%vrt_suffix) ) {
      $name  = $fname . $fld;
      $longn = $pool_name{$pool}{'hist'} . $vrt_suffix{$fld};
      $units;
      if (      $fld eq "C_vr" ) {
         $units = "gC/m^3";
      } elsif ( $fld eq "N_vr" ) {
         $units = "gN/m^3";
      } elsif ( $fld =~ /^N/) {
         $units = "gN/m^2";
      } else {
         $units = "gC/m^2";
      }
      &setField( $name, $longn, $units );
      if ( $fld eq "C" || $fld eq "C_vr" ) {
         foreach my $ciso ( "C13", "C14" ) {
            $name  = $ciso."_".$fname . $fld;
            $longn = $ciso." ".$pool_name{$pool}{'long'} . $vrt_suffix{$fld};
            if ( $fld eq "C_vr" ) {
               $units = "g${ciso}m^3";
            } else {
               $units = "g${ciso}/m^2";
            }
            &setField( $name, $longn, $units );
         }
      }
      if ( $fld =~ "C_1m" || $fld eq "C_30m" || $fld eq "C_activelayer"  ) {
         foreach my $ciso ( "C14" ) {
            $name  = $ciso."_".$fname . $fld;
            $longn = $ciso." ".$pool_name{$pool}{'long'} . $vrt_suffix{$fld};
            $units = "g${ciso}/m^2";
            &setField( $name, $longn, $units );
         }
      }
   }
   # Fire list
   if ( $fname =~ /^CWD/ || $fname =~ /^LIT/ ) {
      foreach my $fld ( keys(%firelist) ) {
         $name  = "M_".$fname . $fld;
         $longn = $firelist{$fname};
         $units;
         if (      $fld =~ /_vr$/ ) {
            $units = "gC/m^3";
         } else {
            $units = "gC/m^2";
         }
         &setField( $name, $longn, $units );
         # Carbon isotopes (C13/C14)
         if ( $fld =~ /^C/ ) {
            foreach my $ciso ( "C13", "C14" ) {
               $name  = "${ciso}_M_".$fname . $fld;
               $longn = $ciso.$firelist{$fname};
               if (      $fld =~ /_vr$/ ) {
                  $units = "g${ciso}/m^3";
               } else {
                  $units = "g${ciso}/m^2";
               }
               &setField( $name, $longn, $units );
            }
         }
      }
   }
   # Potential loss coefficient
   $name  = "K_".$fname;
   $longn = $pool_name{$pool}{'long'} . " potential loss coefficient";
   $units = "1/s";
   &setField( $name, $longn, $units );
   #
   # Not CWD
   #
   if ( $fname !~ /^CWD/ ) {
      foreach my $fld ( keys(%leechlist) ) {
         $name  = "M_".$fname . $fld;
         $longn = $leechlist{$fname};
         my $elm;
         if ( $fld =~ /^N/ ) {
           $elm = "N";
         } else {
           $elm = "C";
         }
         if (      $fld =~ /VERT$/ ) {
            $units = "g${elm}/m^3";
         } else {
            $units = "g${elm}/m^2";
         }
         &setField( $name, $longn, $units );
      }
   }
}
my %translist = (
                   # CN transitions
                   L1S1 =>{d=>"L1",  r=>"S1"},  L2S2 =>{d=>"L2",  r=>"S2"}, 
                   L3S3 =>{d=>"L3",  r=>"S3"},  S1S2 =>{d=>"S1",  r=>"S2"}, 
                   S2S3 =>{d=>"S2",  r=>"S3"},  S3S4 =>{d=>"S3",  r=>"S4"}, 
                   S4   =>{d=>"S4",  r=>"atm"},
                   CWDL2=>{d=>"CWD", r=>"L2"},  CWDL3=>{d=>"CWD", r=>"L3"},
                   # CENTURY transitions NOT already given above
                   L2S1 =>{d=>"L2",  r=>"S1"},  L3S2 =>{d=>"L3",  r=>"S2"}, 
                   S1S3 =>{d=>"S1",  r=>"S3"},  S2S1 =>{d=>"S2",  r=>"S1"}, 
                   S3S1 =>{d=>"S3",  r=>"S1"},  
                );
#
# Transition list (NOT complete)
#
my $unitsvr;
foreach my $trans ( keys(%translist) ) {
   my $donor = $translist{$trans}{'d'};
   my $rcvr  = $translist{$trans}{'r'};
   if ( $trans ne "${donor}${rcvr}" && ($rcvr ne "atm" || $trans ne $donor) ) {
      die "ERROR: Either bad transition name: $trans or bad donor: $donor or receiver:
$rcvr\n";
   }
   # Carbon isotopes
   foreach my $ciso ( "", "C13", "C14" ) {
      if ( $ciso eq "" ) {
         $units   = "gC/m^2/s";
         $unitsvr = "gC/m^3/s";
      } else {
         $units   = "g${ciso}/m^2/s";
         $unitsvr = "g${ciso}/m^3/s";
      }
      if ( $donor ne "CWD" ) {
         my $ii = 0;
         foreach my $trans2 ( keys(%translist) ) {
            if ($donor eq $translist{$trans}{'d'} ) { $ii = $ii + 1; }
         }
         # HR
         if ( $ii == 1 ) {
            $name  = $pool_name{$donor}{'hist'}."_HR";
         } else {
            $name  = $pool_name{$donor}{'hist'}."_HR_$rcvr";
         }
         if ( $ciso ne "" ) {
            $name = "${ciso}$name";
         }
         $longn = 'Het. Resp. from '.$pool_name{$donor}{'long'};
         # vertically integrated fluxes
         &setField( $name,        $longn, $units   );
         # vertically resolved version
         &setField( "${name}_vr", $longn, $unitsvr );
      }
      if ( $rcvr ne "atm" ) {
         # transfer
         $name  = $pool_name{$donor}{'hist'}. "C_TO_" .
                  $pool_name{$rcvr}{'hist'}. "C";
         $longn = "decomp of " . $pool_name{$donor}{'long'}. " C to " .
                  $pool_name{$rcvr}{'long'}. " C";
         if ( $ciso ne "" ) {
            $name = "${ciso}$name";
         }
         # vertically integrated fluxes
         &setField( $name,        $longn, $units   );
         # vertically resolved version
         &setField( "${name}_vr", $longn, $unitsvr );
      }
   }

   #-- mineralization/immobilization fluxes (none from CWD)
   if ( $donor ne "CWD" ) {
      $units   = "gN/m^2/s";
      $unitsvr = "gN/m^3/s";
      if ( $rcvr ne "atm" ) {
         $name  = "SMINN_TO_".$pool_name{$rcvr}{'hist'}. "N_$donor";
      } else {
         $name  = $pool_name{$donor}{'hist'}. "N_TO_SMINN";
      }
      $longn = "mineral N flux for decomp. of " . $pool_name{$donor}{'hist'};
      # vertically integrated fluxes
      &setField( $name,        $longn, $units   );
      # vertically resolved fluxes
      &setField( "${name}_vr", $longn, $unitsvr );
      # transfer fluxes
      if ( $rcvr ne "atm" ) {
         $name  = $pool_name{$donor}{'hist'}. "N_TO_" .
                  $pool_name{$rcvr}{'hist'}. "N";
         $longn = "decomp of " . $pool_name{$donor}{'long'}. " N to " .
                  $pool_name{$rcvr}{'long'}. " N";
         # vertically integrated fluxes
         &setField( $name,        $longn, $units   );
         # vertically resolved fluxes
         &setField( "${name}_vr", $longn, $unitsvr );
      }
      # NITRIF_DENITRIF
      $name  = "SMINN_TO_DENIT_$trans";
      $longn = "denitrification for decomp. of " . $pool_name{$donor}{'long'} .
               "to ". $pool_name{$rcvr}{'hist'};
      &setField( $name,        $longn, "gN/m^2" );
      # vertically resolved fluxes
      &setField( "${name}_vr", $longn, "gN/m^3" );
   }
}

#
# List the fields in a neatly ordered list
# And Output to an XML file
#
my $outfilename = "$pwd/../../../bld/namelist_files/history_fields_clm4_5.xml";

my $outfh = IO::File->new($outfilename, '>') or die "** $ProgName - can't open output history Fields XML file: $outfilename\n";
foreach my $filename ( @filenames ) {
&XML_Header( $outfh, $outfilename, $filename );
foreach my $name ( sort(keys(%fields)) ) {
   my $len;
   if ( length($name) > 20 ) {
      $len = 40;
   } else {
      $len = 20;
   }
   printf( "%-${len}s = %s\n", $name, $fields{$name}{'field'} );
   printf( $outfh "\n<field name='%s' units='%s'\n long_name='%s'\n/>\n", 
           $name, $fields{$name}{'units'}, $fields{$name}{'longn'} );
}
}

&XML_Footer( $outfh );
close( $outfh );
