#!/usr/bin/env perl
#
# Limit the line length for output designed to go into the document.
#
use strict;
use Cwd;
use English;
use IO::File;
use Getopt::Long;
use IO::Handle;
#-----------------------------------------------------------------------------------------------

# Get the directory name and filename of this script.  If the command was
# issued using a relative or absolute path, that path is in $ProgDir.  Otherwise assume the
# command was issued from the current working directory.

(my $ProgName = $0) =~ s!(.*)/!!; # name of this script
my $ProgDir = $1;                 # name of directory containing this script -- may be a
                                  # relative or absolute path, or null if the script
                                  # is in
                                  # the user's PATH
my $nm = "$ProgName::";           # name to use if script dies
my $scrdir;
if ($ProgDir) { 
    $scrdir = $ProgDir;
} else {
    $scrdir = getcwd()
}
my $limitLen  = 99;

sub usage {
    my $msg = shift;

    print "ERROR:: $msg\n";
    die <<EOF;
SYNOPSIS
     $ProgName <input_file>
OPTIONS
    -l      = Limit line length to this value (default $limitLen)
EOF
}

sub LengthofwhiteSpaceNearLength {
  my $line = shift;
  my $leng = shift;

   my $l = $leng;
   while( substr( $line, $l, 1 ) !~ /\s|:|,|\// ) {
      # First search for white-space before desired length -- and then after
      if ( $l <= $leng ) {
         $l--;
      } else {
         $l++;
      }
      # Once reach beginning of line, go to the desired length+1 and increment
      if ( $l < 0 ) { $l = $leng+1; }
      # Once reach the very end of the line die as couldn't break it
      if ( $l >= length($line) ) {
         die "ERROR : went through entire line and did NOT find a place to break it\n";
      }
   }
   return( $l );
}

my %opts = ( limitLen => $limitLen );

GetOptions(
    "l=s"                => \$opts{'limitLen'},
) or usage();

if ( $#ARGV != 0 ) {
    &usage( "Wrong number of command line arguments" );
}

$limitLen = $opts{'limitLen'};
   
my $inputFile = $ARGV[0];

if ( ! -f $inputFile ) {
    &usage( "Input file does NOT exist : $inputFile" );
}

my $fh = IO::File->new($inputFile, '<') or die "** $nm - can't open input file: $inputFile\n";

while (my $line = <$fh>) {

   while( length($line) > $limitLen ) {
      print STDERR "Line length over $limitLen\n";
      my $lenlim = &LengthofwhiteSpaceNearLength( $line, $limitLen );
      if ( ($lenlim == length($line)) || $lenlim < 0 ) {
         print "Can NOT truncate long line: $line\n";
         die "ERROR : Having trouble breaking a long line\n";
      }
      my $substring = substr( $line, 0, $lenlim+1 );
      print "$substring \\ \n";
      my $newline = "   " . substr( $line, $lenlim+1, length($line) );
      $line = $newline;
   }
   print $line;

}
$fh->close;


