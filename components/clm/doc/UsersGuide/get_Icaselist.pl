#!/usr/bin/env perl
#-----------------------------------------------------------------------------------------------
#
# get_Icaselist.pl
#
# This utility gets a list of the I cases from the CCSM compset database.
#
#-----------------------------------------------------------------------------------------------

use strict;
use Cwd;
use English;
use Getopt::Long;
use IO::File;
use IO::Handle;
#-----------------------------------------------------------------------------------------------

sub usage {
    die <<EOF;
SYNOPSIS
     get_Icaselist.pl  [options]
OPTIONS
EOF
}

#-----------------------------------------------------------------------------------------------
# Setting autoflush (an IO::Handle method) on STDOUT helps in debugging.  It forces the test
# descriptions to be printed to STDOUT before the error messages start.

*STDOUT->autoflush();                  

#-----------------------------------------------------------------------------------------------
my $cwd = getcwd();      # current working directory
my $cfgdir;              # absolute pathname of directory that contains this script
$cfgdir = $cwd;

#-----------------------------------------------------------------------------------------------
# Parse command-line options.
my %opts = (
	    );
GetOptions(
    "h|help"                    => \$opts{'help'},
)  or usage();

# Give usage message.
usage() if $opts{'help'};

# Check for unparsed argumentss
if (@ARGV) {
    print "ERROR: unrecognized arguments: @ARGV\n";
    usage();
}

# Check for manditory case input if not just listing valid values

my %cfg = ();           # build configuration

#-----------------------------------------------------------------------------------------------

# Check for the configuration definition file.
my $config_def_file = "config_definition.xml";
my $case_def_dir    = "$cfgdir/../../../../../scripts/ccsm_utils/Case.template";
(-f "$case_def_dir/$config_def_file")  or  die <<"EOF";
** Cannot find configuration definition file \"$config_def_file\" in directory 
    \"$case_def_dir\" **
EOF

# Compset definition file.
my $compset_file = 'config_compsets.xml';
(-f "$case_def_dir/$compset_file")  or  die <<"EOF";
** Cannot find compset parameters file \"$compset_file\" in directory 
    \"$case_def_dir\" **
EOF

my $xml_dir = "$cfgdir/../../../../../scripts/ccsm_utils/Tools/perl5lib";
# The XML::Lite module is required to parse the XML configuration files.
(-f "$xml_dir/XML/Lite.pm")  or  die <<"EOF";
** Cannot find perl module \"XML/Lite.pm\" in directory 
    \"$xml_dir\" **
EOF


#-----------------------------------------------------------------------------------------------
my @dirs = (  $cfgdir, $xml_dir, $case_def_dir );
unshift @INC, @dirs;
require XML::Lite;
require ConfigCase;

#-----------------------------------------------------------------------------------------------
my $cfg_ref = ConfigCase->new("$case_def_dir/$config_def_file"); 
print_compsets( "$case_def_dir/$compset_file" );

#-----------------------------------------------------------------------------------------------
# FINNISHED ####################################################################################
#-----------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------

sub print_compsets
{
    # Print all currently supported valid compsets

    my ($compset_file) = @_;
    my $xml = XML::Lite->new( $compset_file );
    my $root = $xml->root_element();

    # Check for valid root node
    my $name = $root->get_name();
    $name eq "config_compset" or die
	"file $compset_file is not a compset parameters file\n";

    # Read the compset parameters from $compset_file.
    my @e = $xml->elements_by_name( "compset" );
    my %a = ();
    my %data;
    while ( my $e = shift @e ) {
	%a        = $e->get_attributes();
        my $sname = $a{'SHORTNAME'};
	if ($a{GRID_MATCH} && exists($data{$sname}) && defined($data{$sname}{'DESC'} && defined($a{'DESC'}) )  ) {
            if ( $data{$sname}{'DESC'} =~ /^INVALID:/ ) {
                $data{$sname}{'DESC'} = $a{'DESC'};
            }
	} elsif ( $a{'SHORTNAME'} =~ /^I/ ) {
            $data{$sname}{'NAME'} = $a{'NAME'};
            $data{$sname}{'DESC'} = $a{'DESC'};
	}
    }
    print "<orderedlist>\n";
    foreach my $sname ( sort(keys(%data)) ) {
        print "<listitem><para><varname>$data{$sname}{'NAME'}</varname>" .
              "(<varname>$sname</varname>)\n";
	print "$data{$sname}{'DESC'}</para></listitem>\n";
    }
    print "</orderedlist>\n";
}

