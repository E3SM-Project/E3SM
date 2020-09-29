#!/usr/bin/env perl
#
# July 18 2012                                         Muszala
#
# createMapEntry.pl - A simple script to dump a list of mappings for a specified resolution to then
# cut and paste into namelist_defaults_clm.xml.  A better way is to write the output of this script
# to a file and then directly insert that file into namelist_defaults_clm.xml (using :r foo in vim for
# example).
#
# Example usage:>> ./createMapEntry.pl 1x1_brazil   
#    will create XML entries for maps in ../lnd/clm2/mappingdata/maps/1x1_brazil  such as:
#
#    <map frm_hgrid="0.5x0.5"    frm_lmask="AVHRR"  to_hgrid="1x1_brazil"   to_lmask="nomask" 
#    >lnd/clm2/mappingdata/maps/1x1_brazil/map_0.5x0.5_AVHRR_to_1x1_brazil_nomask_aave_da_c120717.nc</map>
#
use Cwd;
use strict;
use English;
use IO::File;
use Getopt::Long;

   my $date = scalar localtime() ;
   my $scriptName;
   ($scriptName = $0) =~ s!(.*)/!!; # get name of script
   my $cwd = getcwd();
   my $CSMDATA = "/glade/p/cesm/cseg/inputdata";

   if ($#ARGV != 0 ) {
      usage();
	   exit;
   }
   my $grid=$ARGV[0];

   sub usage {
      die <<EOF;
         SYNOPSIS 
   
         $scriptName <res>  
            <res> is the resolution to use to dump text to paste into namelist_defaults_clm.xml
EOF
   }

   #~# set up directory paths
   my $pathStub="lnd/clm2/mappingdata/maps";
   my $partialPath="$pathStub/$grid";
   my $fullPath = "$CSMDATA/$partialPath";

   #~# open and read directory
   opendir DIR, $fullPath or die "Cannot read dir! $fullPath";
   my @list = readdir DIR;

   #~# print a unique start string in the XML comments 
   print "\n<!-- mapping files for $grid START added on $date-->";
   print "\n<!-- Created by lnd/clm/bld/namelist_files/$scriptName--> \n\n";

   foreach my $foo ( @list ) {
      next if ($foo =~ m/^\./);  #~# skip anything in the directory with a leading or stand alone 'dot'
      my @tokens = split(/_/, $foo); #~# split foo name by the underscore
      #~# write out lines for namelist_defaults_clm.xml
      print "<map frm_hgrid=\"$tokens[1]\"    frm_lmask=\"$tokens[2]\"  to_hgrid=\"$tokens[4]\"   to_lmask=\"$tokens[5]\" \n";
      print ">$partialPath/$foo</map>\n";
   }

   #~# print a unique end string in the XML comments 
   print "\n<!-- mapping files for $grid END --> \n";
   closedir(DIR);
   exit 0;
