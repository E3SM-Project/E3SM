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
my $nl_definition_file = "$cfgdir/../../../bld/namelist_files/namelist_definition_clm4_0.xml";
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
    $fields{$name} = sprintf( "%-${len}s\t(%s)", $longn, $units );
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
  my $svnurl = '$URL: https://svn-ccsm-models.cgd.ucar.edu/clm2/trunk_tags/clm4_5_71/models/lnd/clm/src/clm4_0/main/findHistFields.pl $';
  my $svnid  = '$Id: findHistFields.pl 46656 2013-05-02 06:19:47Z erik $';
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
my $filename = "$pwd/histFldsMod.F90";

my $fh = IO::File->new($filename, '<') or die "** $ProgName - can't open history Fields file: $filename\n";
my @megcmpds =  $definition->get_valid_values( "megan_cmpds", 'noquotes'=>1 );
#
# Read in the list of fields from the source file
# And output to an XML file
#
my $outfilename = "$pwd/../../../bld/namelist_files/history_fields_clm4_0.xml";

my $outfh = IO::File->new($outfilename, '>') or die "** $ProgName - can't open output history Fields XML file: $outfilename\n";
&XML_Header( $outfh, $outfilename, $filename );
while (my $line = <$fh>) {

   # Comments
   if ($line =~ /(.*)\!/) {
     $line = $1;
   }
   my $format = "\n<field name='%s' units='%s'\n long_name='%s'\n/>\n";
   if ($line =~ /call\s*hist_addfld/i ) {
      (my $name, my $longn, my $units) = &getFieldInfo( $fh, $line );
      if ( $name ne "MEG_megancmpd" ) {
         &setField( $name, $longn, $units );
         printf( $outfh $format, $name, $units, $longn );
      } else {
         foreach my $megcmpd ( @megcmpds ) {
            my $name = "MEG_${megcmpd}";
            &setField( $name, $longn, $units );
            printf( $outfh $format, $name, $units, $longn );
         }
      }
   }
}
close( $fh );
&XML_Footer( $outfh );
close( $outfh );
print STDERR " mxname  = $mxname\n";
print STDERR " mxlongn = $mxlongn\n";

#
# List the fields in a neatly ordered list
#
foreach my $name ( sort(keys(%fields)) ) {
   my $len;
   if ( length($name) > 20 ) {
      $len = 40;
   } else {
      $len = 20;
   }
   printf( "%-${len}s = %s\n", $name, $fields{$name} );
}

