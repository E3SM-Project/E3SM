#!/usr/bin/env perl

# Creates a file giving XML entries for all the mapping files in the
# current directory (mapping_entries.txt). Also creates another file
# giving commands to move these files to the inputdata space
# (mv_cmds.sh).
#
# Should be run with no arguments.
#
# See also bld/namelist_files/createMapEntry.pl, and mvNimport.sh in
# the current directory for scripts that share some of the
# functionality of this script.

# Bill Sacks
# March, 2013

use strict;

# ----------------------------------------------------------------------
# FUNCTIONS
# ----------------------------------------------------------------------

# Given a map filename, returns a hash giving the resolutions and
# masks implicit in that filename.
# Inputs:
#   - filename
# Output:
#   - hash containing:
#     - filename
#     - from_res
#     - from_mask
#     - to_res
#     - to_mask
#   Or does a bare return if the filename doesn't match the expected pattern
sub get_resolutions_and_masks {
   my $filename = shift;

   # The following match assumes that the destination mask is
   # "nomask".  This match will tolerate underscores in the
   # destination grid (e.g., 5x5_amazon), but be careful about
   # underscores in the source grid or source mask!
   if ($filename =~ m/^map_(.*)_(.*)_to_(.*)_nomask/) {
      my $from_res=$1;
      my $from_mask=$2;
      my $to_res=$3;
      my $to_mask="nomask";

      my %info = (filename  => $filename,
                  from_res  => $from_res,
                  from_mask => $from_mask,
                  to_res    => $to_res,
                  to_mask   => $to_mask);
         
      return %info;
   }
   else {
      return;
   }
}
                  

# ----------------------------------------------------------------------
# PARAMETERS DEFINED HERE
# ----------------------------------------------------------------------

my $CSMDATA = "/glade/p/cesm/cseg/inputdata";
my $maps_dir = "lnd/clm2/mappingdata/maps";   # directory where mapping files are stored within the inputdata directory

# ----------------------------------------------------------------------
# BEGIN MAIN PROGRAM
# ----------------------------------------------------------------------

my @files = glob "map*.nc";

# Make a hash containing all of the files at each destination resolution.
# The keys of the hash are destination resolutions; the values are
# references to arrays of hash references, where these low-level
# hashes are the return values of get_resolutions_and_masks.
my %dest_resols;
foreach my $file (@files) {
   my %info = get_resolutions_and_masks($file);
   if (%info) {
      my $to_res = $info{'to_res'};
      push @{$dest_resols{$to_res}}, \%info;
   }
   else {
      warn "WARNING: $file doesn't match expected mapping filename pattern; skipping\n";
   }
}      

open MAP_ENTRIES, ">", "mapping_entries.txt";
open MV_CMDS, ">", "mv_cmds.sh";

# Output xml entries (and mv commands) grouped by destination resolution
foreach my $to_res (sort keys %dest_resols) {
   my $full_maps_dir = "$maps_dir/$to_res";

   foreach my $info_ref (@{$dest_resols{$to_res}}) {
      my $filename  = ${$info_ref}{'filename'};
      my $from_res  = ${$info_ref}{'from_res'};
      my $from_mask = ${$info_ref}{'from_mask'};
      my $to_res    = ${$info_ref}{'to_res'};
      my $to_mask   = ${$info_ref}{'to_mask'};
      
      print MV_CMDS "mv $filename $CSMDATA/$full_maps_dir/$filename\n";
      print MAP_ENTRIES "<map frm_hgrid=\"$from_res\"    frm_lmask=\"$from_mask\"  to_hgrid=\"$to_res\"   to_lmask=\"$to_mask\" \n";
      print MAP_ENTRIES ">$full_maps_dir/$filename</map>\n";
   }

   # Print blank line between destination grids
   print MAP_ENTRIES "\n";
}

system "chmod", "755", "mv_cmds.sh";
close MAP_ENTRIES;
close MV_CMDS;
