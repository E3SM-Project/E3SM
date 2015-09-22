#!/usr/bin/env perl
#
use strict;
use Cwd;
use English;
use IO::File;
use Getopt::Long;
use IO::Handle;
#-----------------------------------------------------------------------------------------------

# Get the directory name and filename of this script.  If the command was
# issued using a relative or absolute path, that path is in $ProgDir.  Otherwise assume
# the
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

sub usage {
    my $msg = shift;

    print "ERROR:: $msg\n";
    die <<EOF;
SYNOPSIS
     $ProgName <input_file>
OPTIONS
     NONE
EOF
}

my %opts = ( );

GetOptions(
) or usage();

if ( $#ARGV != 0 ) {
    &usage( "Wrong number of command line arguments" );
}

my $inputFile = $ARGV[0];

if ( ! -f $inputFile ) {
    &usage( "Input file does NOT exist : $inputFile" );
}

my $fh = IO::File->new($inputFile, '<') or die "** $nm - can't open input file:
$inputFile\n";

#
# Add in XML XHTML headers
#
print <<"EOF";
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE html>
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">
EOF
while (my $line = <$fh>) {
   if      ( $line =~ /^<html/ ) { 
       # Ignore this line as already have html line above
   } elsif ( $line =~ /^<hr>$/ ) { 
       print "<hr/>\n";
   } elsif ( $line =~ /^<meta/ ) { 
       print '<meta http-equiv="Content-Type" content="text/html; charset=UTF-8"/>'."\n";
   } else {
       print $line;
   }
}
$fh->close();
