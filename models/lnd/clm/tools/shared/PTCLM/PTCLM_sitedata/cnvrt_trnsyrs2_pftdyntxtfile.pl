#!/usr/bin/env perl
#
# cnvrt_trnsyrs2_pftdyntxtfile.pl                  Erik Kluzek
#                                                  Aug/5/2010
#
# Convert the transition years files to pftdyn text files.
#
use Cwd;
use strict;
use English;
use IO::File;
use Getopt::Long;

#
# Some global constants
#
my $maxlen    = 125;
my $numharv   = 5;
my $numgraz   = 1;
my $nbreak    = 5*2;
my $numarray  = 1 + $nbreak + $numharv + $numgraz + 2 - 1;
my @hd_pftarr;
my @hd_hrvarr;
my @hd_grzarr;

sub parse_header {
#
# Parse the header and make sure it's correct
#
  my $header = shift;

  my @harray    = split( /,/, $header );
  if ( $#harray != $numarray ) {
     die "** Number of elements in line is incorrect: $#harray should be: $numarray\n";
  }
  if ( (my $hyear = shift( @harray )) ne "trans_year" ) {
     die "** First header element is NOT trans_year as expected: $hyear\n";
  }
  foreach my $var ( "hold_graze\n", "hold_harv" ) {
     if ( (my $val = pop( @harray )) ne $var ) {
        die "** Last header elements are NOT $var as expected: $val\n";
     }
  }
  for( my $i = 0; $i < $nbreak; $i++ ) {
     push( @hd_pftarr, shift( @harray ) );
  }
  for( my $i = 0; $i < $numharv; $i++ ) {
     push( @hd_hrvarr, shift( @harray ) );
  }
  for( my $i = 0; $i < $numgraz; $i++ ) {
     push( @hd_grzarr, shift( @harray ) );
  }
}

sub parse_pft {
#
# Parse the PFT array
#
   my $frcname  = shift;
   my $idxname  = shift;
   my $head     = shift;
   my @pftarray = @_;

   my @header = split( /,/, $head );
   my $sum  = 0.0;
   my $frcline = "<$frcname>";
   my $idxline = "<$idxname>";
   my $n    = 1;
   my $endit= undef;
   for( my $i = 0; $i <= $#pftarray; $i=$i+2 ) {
       my $j    = $i + 1;
       my $expect = "pft_f$n";
       if ( $header[$i] ne $expect ) {
          die "**  PFT fraction header is wrong: $header[$i] expect $expect\n";
       }
       if ( ! defined($endit) ) { $frcline   .= $pftarray[$i]; }
       $sum     = $sum + $pftarray[$i];
       if ( $pftarray[$j] < 0 || $pftarray[$j] > 16 ) {
          die "** PFT index is out of range: $pftarray[$j]\n";
       }
       $expect = "pft_c$n";
       if ( $header[$j] ne $expect ) {
          die "**  PFT code header is wrong: $header[$j] expect $expect\n";
       }
       if ( ! defined($endit) ) { $idxline .= $pftarray[$j]; }
       if (      $sum >  100.0 ) { 
          die "** Sum of PFT fractions exceeds 100: $sum\n";
       } elsif ( $sum == 100.0 ) {
          $endit = 1;
       }
       if ( ! defined($endit) && ($j < $#pftarray) ) { $frcline .= ","; }
       if ( ! defined($endit) && ($j < $#pftarray) ) { $idxline .= ","; }
       $n++;
   }
   if ( $sum != 100.0 ) { 
      die "** Sum of PFT fractions does NOT go to 100: $sum\n";
   }
   $frcline .= "</$frcname>";
   $idxline .= "</$idxname>";

   return( "$frcline$idxline" );
}

sub parse_hrv {
#
# Parse the harvesting array
#
   my $name  = shift;
   my $head  = shift;
   my $exp   = shift;
   my @array = @_;

   my @exp    = split( /,/, $exp  );
   my @header = split( /,/, $head );
   my $line = "<$name>";
   for( my $i = 0; $i <= $#array; $i++ ) {
       if ( $header[$i] ne $exp[$i] ) {
          die "**  harvesting header is wrong: $header[$i], expecting $exp[$i]\n";
       }
       $line .= $array[$i];
       if ( $array[$i] < 0 || $array[$i] > 1 ) {
          die "** Bad value for harvest value: $array[$i]\n";
       }
       if ( $i < $#array ) { $line .= ","; }
   }
   $line .= "</$name>";

   return( $line );
}

sub printoutline {
#
# Print the line out with the fill as well
#
   my $pftline = shift;
   my $harline = shift;
   my $grzline = shift;
   my $year    = shift;

   my $outline = "$pftline$harline$grzline";
   my $length = length( $outline );
   if ( $length > $maxlen ) {
      die "** line length is too long = $length\n";
   }
   my $fill = "";
   for( my $i = 1; $i <= $maxlen - $length; $i++ ) { $fill .= " "; }

   print "$outline$fill $year\n";

}

if ( $#ARGV != 1 ) {
   die "** Wrong number of arguments: should just be filename and sim_year_range\n ".
       "$0 <filename> <sim_year_range>";
}
my $filename       = $ARGV[0];
my $sim_year_range = $ARGV[1];
chomp( $sim_year_range );
my $start_year;
my $end_year;
if ( $sim_year_range =~ /^([0-9]+)\-([0-9]+)$/ ) {
   $start_year = $1;
   $end_year   = $2;
} else {
  die "** bad format for sim_year_range (should be yyyy-yyyy): $sim_year_range\n";
}

my $fh = IO::File->new;
$fh->open( "<$filename" ) or die "** can't open file: $filename\n";

my $header    = <$fh>;
my $last_year = undef;
my $frst_year = undef;

&parse_header( $header );

my $pftline   = "";
my $harline   = "";
my $grzline   = "";
my $hold_h    = undef;
my $hold_g    = undef;
my $year      = undef;
while( my $line = <$fh> ) {
   my @array   = split( /,/, $line );
   if ( $#array != $numarray ) {
      die "** Number of elements in line is incorrect: $#array\n";
   }
   $year    = shift( @array );
   #
   # Write out the years from last year until current year
   #
   if ( defined($last_year) ) {
      for( my $yr = $last_year+1; $yr < $year; $yr++ ) {
         if ( $yr >= $start_year ) {
            &printoutline( $pftline, $harline, $grzline, $yr );
         }
      }
   }
   #
   # Last two parts of the array are harvesting and grazing hold values
   #
   $hold_g  = pop( @array );
   chomp( $hold_g );
   $hold_h  = pop( @array );
   #
   # Separate out the array into the different sections
   #
   my @pftarray;
   for( my $i = 0; $i < $nbreak; $i++ ) {
      push( @pftarray,  shift( @array  ) );
   }
   my @hrvarray;
   for( my $i = 0; $i < $numharv; $i++ ) {
      push( @hrvarray,  shift( @array ) );
   }
   my @grzarray;
   for( my $i = 0; $i < $numgraz; $i++ ) {
      push( @grzarray,  shift( @array ) );
   }
   #
   # Parse the different sections
   #
   $pftline = &parse_pft( "pft_f", "pft_i", join( ",", @hd_pftarr ), @pftarray );
   $harline = &parse_hrv( "harv",  join( ",", @hd_hrvarr ), "har_vh1,har_vh2,har_sh1,har_sh2,har_sh3", @hrvarray );
   $grzline = &parse_hrv( "graz",  join( ",", @hd_grzarr ), "graze", @grzarray );

   #
   # Write out the years from start year until first year
   #
   if ( ! defined($frst_year) && ($year > $start_year) ) {
      for( my $yr = $start_year; $yr < $year; $yr++ ) {
         &printoutline( $pftline, $harline, $grzline, $yr );
      }
   }
   #
   # Figure out the line length and the amount of fill to have and print it
   #
   if ( $year >= $start_year ) {
      &printoutline( $pftline, $harline, $grzline, $year );
   }
   # If NOT holding harvesting, set it to zero, for transition years
   if ( ! $hold_h ) {
      for( my $i = 0; $i <= $#hrvarray; $i++ ) { $hrvarray[$i] = 0; }
      $harline = &parse_hrv( "harv", join( ",", @hd_hrvarr ), "har_vh1,har_vh2,har_sh1,har_sh2,har_sh3", @hrvarray );
   }
   # If NOT holding grazing, set it to zero, for transition years
   if ( ! $hold_g ) {
      for( my $i = 0; $i <= $#grzarray; $i++ ) { $grzarray[$i] = 0; }
      $grzline = &parse_hrv( "graz", join( ",", @hd_grzarr ), "graze", @grzarray );
   }
   # Save last years value so can create transition years
   $last_year = $year;
   if ( ! defined($frst_year) ) { $frst_year = $year }
}
#
# Write out years from end to file to the end_year
#
for( my $yr = $year+1; $yr <= $end_year; $yr++ ) {
   &printoutline( $pftline, $harline, $grzline, $yr );
}
$fh->close();
