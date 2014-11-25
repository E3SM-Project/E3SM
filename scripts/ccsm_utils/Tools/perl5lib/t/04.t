#!/usr/bin/env perl

# Test methods of the Namelist objects.

#########

use Test::More tests=>73;  # use "no_plan" until number of tests is determined

########


use strict;

use lib "..";

BEGIN {use_ok('Build::Config')};

BEGIN {use_ok('Build::Namelist')};

BEGIN {use_ok('Build::NamelistDefaults')};

BEGIN {use_ok('Build::NamelistDefinition')};

# check that XML::List is being found

BEGIN {use_ok('XML::Lite')};

#
my $nl_definition_file  = "namelist_definition_cam.xml";
my $nl_definition_file2 = "namelist_definition_clm.xml";
my $nl_defaults_file    = "namelist_defaults_cam.xml";
my $nl_defaults_file2   = "namelist_defaults_clm.xml";
my $config_file         = "config_cache.xml";
my $inputdata_rootdir   = "/fs/cgd/csm/inputdata";

# Create a namelist definition object.  This object provides a method for verifying that
# the
# output namelist variables are in the definition file, and are output in the correct
# namelist groups.
my $definition = Build::NamelistDefinition->new($nl_definition_file);
$definition->add( $nl_definition_file2 );

# Create a namelist defaults object.  This object provides default values for variables
# contained in the input defaults file.  The configuration object provides attribute
# values that are relevent for the CAM executable for which the namelist is being
# produced.
my $cfg = Build::Config->new($config_file);
my $defaults = Build::NamelistDefaults->new($nl_defaults_file, $cfg);
$defaults->add($nl_defaults_file2);

my $nl_defaults_usr_file = "namelist_defaults_usr_files.xml";
my $uf_defaults = Build::NamelistDefaults->new($nl_defaults_usr_file, $cfg);
my %settings;
$settings{'mask'}           = "USGS";
$settings{'sim_year'}       = "1850";
$settings{'sim_year_range'} = "1850-2000";
$settings{'clm_usr_name'}   = "1x1_boulderCO_c090804";
$settings{'csmdata'}        = "\$PWD";
my @nonexist_files     = ( "faerdep", "finidat", "fndepdat", "fndepdyn" );
my @nonexist_filenames = ( 
                           "\$PWD/usrfiles/aerosoldep_monthly_1850-2000_1x1_boulderCO_c090804.nc",
                           "\$PWD/usrfiles/clmi.1x1_boulderCO_c090804_USGS_simyr1850.nc", 
                           undef,
                           "\$PWD/usrfiles/fndep_clm_1850-2000_1x1_boulderCO_c090804.nc"
                         );
my %usr_files;
foreach my $var ( $uf_defaults->get_variable_names() ) {
   my $val = $uf_defaults->get_usr_file($var,$definition,\%settings);
   if ( defined($val) ) { 
       $usr_files{$var} = $val;
    } else {
       is( $var, shift(@nonexist_files), "Make sure that expected non-exist user files do NOT return true: $var" );
       $settings{'notest'} = 1;
       $val = $uf_defaults->get_usr_file($var,$definition,\%settings);
       is( $val, shift(@nonexist_filenames), "Expected non-exist user files do return if notest" );
       $settings{'notest'} = 0;
    }
}

my %expected_usr_files = ( 
                           'fatmgrid'   => "\$PWD/usrfiles/griddata_1x1_boulderCO_c090804.nc", 
                           'fatmlndfrc' => "\$PWD/usrfiles/fracdata_1x1_boulderCO_c090804_USGS.nc",
                           'fatmtopo'   => "\$PWD/usrfiles/topodata_1x1_boulderCO_c090804.nc",
                           'flndtopo'   => "\$PWD/usrfiles/topodata_1x1_boulderCO_c090804.nc",
                           'flanduse_timeseries'    => "\$PWD/usrfiles/landuse.timeseries_1x1_boulderCO_c090804_simyr1850-2000.nc",
                           'fsurdat'    => "\$PWD/usrfiles/surfdata_1x1_boulderCO_c090804_simyr1850.nc",
                          );
is_deeply( \%usr_files, \%expected_usr_files, "Make sure return expected user default files" );


# Create an empty namelist object.  Add values to it in order of precedence.
my $nl = Build::Namelist->new();

#-------------------------------------------------------------------------------

sub quote_string {
    my $str = shift;
    $str =~ s/^\s+//;
    $str =~ s/\s+$//;
    unless ($str =~ /^['"]/) {        #"'
        $str = "\'$str\'";
    }
    return $str;
}

#-----------------------------------------------------------------------------------------------

sub set_abs_filepath {

# check whether the input filepath is an absolute path, and if it isn't then
# prepend a root directory

    my ($filepath, $rootdir) = @_;

    # strip any leading/trailing whitespace
    $filepath =~ s/^\s+//;
    $filepath =~ s/\s+$//;
    $rootdir  =~ s/^\s+//;
    $rootdir  =~ s/\s+$//;

    # strip any leading/trailing quotes
    $filepath =~ s/^['"]+//;
    $filepath =~ s/["']+$//;
    $rootdir =~ s/^['"]+//;
    $rootdir =~ s/["']+$//;

    my $out = $filepath;
    unless ( $filepath =~ /^\// ) {  # unless $filepath starts with a /
  $out = "$rootdir/$filepath"; # prepend the root directory
    }
    return $out;
}

#-----------------------------------------------------------------------------------------------


sub absolute_path {
#
# Convert a pathname into an absolute pathname, expanding any . or .. characters.
# Assumes pathnames refer to a local filesystem.
# Assumes the directory separator is "/".
#
  my $path = shift;
  my $cwd = getcwd();  # current working directory
  my $abspath;         # resulting absolute pathname

# Strip off any leading or trailing whitespace.  (This pattern won't match if
# there's embedded whitespace.
  $path =~ s!^\s*(\S*)\s*$!$1!;

# Convert relative to absolute path.

  if ($path =~ m!^\.$!) {          # path is "."
      return $cwd;
  } elsif ($path =~ m!^\./!) {     # path starts with "./"
      $path =~ s!^\.!$cwd!;
  } elsif ($path =~ m!^\.\.$!) {   # path is ".."
      $path = "$cwd/..";
  } elsif ($path =~ m!^\.\./!) {   # path starts with "../"
      $path = "$cwd/$path";
  } elsif ($path =~ m!^[^/]!) {    # path starts with non-slash character
      $path = "$cwd/$path";
  }

  my ($dir, @dirs2);
  my @dirs = split "/", $path, -1;   # The -1 prevents split from stripping trailing nulls
                                     # This enables correct processing of the input "/".

  # Remove any "" that are not leading.
  for (my $i=0; $i<=$#dirs; ++$i) {
      if ($i == 0 or $dirs[$i] ne "") {
    push @dirs2, $dirs[$i];
      }
  }
  @dirs = ();

  # Remove any "."
  foreach $dir (@dirs2) {
      unless ($dir eq ".") {
    push @dirs, $dir;
      }
  }
  @dirs2 = ();

  # Remove the "subdir/.." parts.
  foreach $dir (@dirs) {
    if ( $dir !~ /\.\./ ) {
        push @dirs2, $dir;
    } else {
        pop @dirs2;   # remove previous dir when current dir is ..
    }
  }
  if ($#dirs2 == 0 and $dirs2[0] eq "") { return "/"; }
  $abspath = join '/', @dirs2;
  return( $abspath );
}

sub add_default {

# Add a value for the specified variable to the specified namelist object.  The variables
# already in the object have the higher precedence, so if the specified variable is
# already
# defined in the object then don't overwrite it, just return.
#
# This method checks the definition file and adds the variable to the correct
# namelist group.
#
# The value can be provided by using the optional argument key 'val' in the
# calling list.  Otherwise a default value is obtained from the namelist
# defaults object.  If no default value is found this method throws an exception
# unless the 'nofail' option is set true.
#
# Example 1: Specify the default value $val for the namelist variable $var in namelist
#            object $nl:
#
#  add_default($nl, $var, 'val'=>$val)
#
# Example 2: Add a default for variable $var if an appropriate value is found.  Otherwise
#            don't add the variable
#
#  add_default($nl, $var, 'nofail'=>1)
#
#
# ***** N.B. ***** This routine assumes the following variables are in package main::
#  $definition        -- the namelist definition object
#  $defaults          -- the namelist defaults object
#  $inputdata_rootdir -- CCSM inputdata root directory

    my $nl = shift;     # namelist object
    my $var = shift;    # name of namelist variable
    my %opts = @_;      # options

    # Query the definition to find which group the variable belongs to.  Exit if not
    # found.
    my $group = $definition->get_group_name($var);
    unless ($group) {
  my $fname = $definition->get_file_name();
  die "ERROR: variable \"$var\" not found in namelist definition file $fname.\n";
    }

    # check whether the variable has a value in the namelist object -- if so then skip to the end
    my $val = $nl->get_variable_value($group, $var);
    if ( ! defined( $val ) ) {

       # Look for a specified value in the options hash
       if (defined $opts{'val'}) {
           $val = $opts{'val'};
       }
       # or else get a value from namelist defaults object.
       # Note that if the 'val' key isn't in the hash, then just pass anything else
       # in %opts to the get_value method to be used as attributes that are matched
       # when looking for default values.
       else {
           $val = $defaults->get_value($var, \%opts);
       }

       # if no value is found then exit w/ error (unless 'nofail' option set)
       unless ($val) {
           unless ($opts{'nofail'}) {
               die "No default value found for $var";
           }
           else {
               return;
           }
       }

       # query the definition to find out if the variable is an input pathname
       my $is_input_pathname = $definition->is_input_pathname($var);

       # The default values for input pathnames are relative.  If the namelist
       # variable is defined to be an absolute pathname, then prepend
       # the CCSM inputdata root directory.
       if ($is_input_pathname eq 'abs') {
           $val = set_abs_filepath($val, $inputdata_rootdir);
       }

       # query the definition to find out if the variable takes a string value.
       # The returned string length will be >0 if $var is a string, and 0 if not.
       my $str_len = $definition->get_str_len($var);
   
       # If the variable is a string, then add quotes if they're missing
       if ($str_len > 0) {
           $val = quote_string($val);
       }

       # set the value in the namelist
       $nl->set_variable_value($group, $var, $val);

       # Get the type description hash for the variable
       # and validate the value with the datatype in the namelist object
       # This method throws an exception when an error is encountered.
       my %type_ref = $definition->_get_datatype($var);
       Build::Namelist::validate_variable_value($var,$val, \%type_ref);
    }
}


#-------------------------------------------------------------------------------

my $var = "scenario_prognostic_sulfur";
my @valid_values = $definition->get_valid_values( $var );
is( $#valid_values, 0, "Make sure can get valid values for a single value" );
is( $valid_values[0], "'RAMPED'", "Make sure can get valid values for a single valid option" );
my $val = "'RAMPED'";
my $group = $definition->get_group_name($var);
$nl->set_variable_value($group, $var, $val);
my %type_def = $definition->_get_datatype($var);
my $rc = Build::Namelist::validate_variable_value( $var, $val, \%type_def );
is( $rc, 1, "Test that validate_variable_value method works" );
my $var = "soil_erod";
my @valid_values = $definition->get_valid_values( $var );
is( $#valid_values, -1, "Make sure can get valid values when empty" );

my $var = "phys_alltoall";
my @valid_values = $definition->get_valid_values( $var );
my @exp_values = ( 0, 1, 2, 11, 12, 13 );
is( @valid_values, @exp_values, "Make sure can get valid values for integer array" );
my $val = 13;
my $group = $definition->get_group_name($var);
$nl->set_variable_value($group, $var, $val);
my %type_def = $definition->_get_datatype($var);
my $rc = Build::Namelist::validate_variable_value( $var, $val, \%type_def );
is( $rc, 1, "Test that validate_variable_value method works" );

my $var = "diag_cnst_conv_tend";
my @valid_values = $definition->get_valid_values( $var );
my $val = "'q_only'";
@exp_values=("'none'","'q_only'","'all'");
is( @valid_values, @exp_values, "Make sure can get valid values for string array" );
my $group = $definition->get_group_name($var);
$nl->set_variable_value($group, $var, $val);
my %type_def = $definition->_get_datatype($var);
my $rc = Build::Namelist::validate_variable_value( $var, $val, \%type_def );
is( $rc, 1, "Test that validate_variable_value method works" );

my $var = "stop_option";
my $val = "'nyears'";
$nl->set_variable_value( "generic_namelist", $var, $val);
my %type_def = $definition->_get_datatype($var);
my $rc = Build::Namelist::validate_variable_value( $var, $val, \%type_def );
is( $rc, 1, "Test that setting a namelist item in a generic_namelist works" );

my $var = "diag_cnst_conv_tend";
my @valid_values = $definition->get_valid_values( $var, 'noquotes'=>0 );
my $val = "q_only";
@exp_values=("none","q_only","all");
is( @valid_values, @exp_values, "Make sure can get non-quoted valid values for string array" );

my $var = "case_name";
my $group = $definition->get_group_name($var);
my $input = "'case_name_in_quotes'";
add_default($nl, $var, 'val'=>"$input" );
my $result = $nl->get_variable_value($group,$var);

is( $result, $input, "Checking that if set $var you get back what you set" );

# Add some defaults from the file
# Include some changes to the case of the variable name to make sure it will work.
add_default($nl, 'start_type');
add_default($nl, 'fPftcon');
add_default($nl, 'fsurDat');
add_default($nl, 'finIdAt',    'nofail'=>1);
add_default($nl, 'fatmgrid',   'nofail'=>1);
add_default($nl, 'fatMlndfrc', 'nofail'=>1);
add_default($nl, 'restART_option');
add_default($nl, 'START_YMD');
add_default($nl, 'stop_option');
add_default($nl, 'stop_N');
add_default($nl, 'orb_IYEAR_ad');
add_default($nl, 'dtime');
my $val = "'thingfirst'";
for( my $i=0; $i<99; $i++ ) {
  $val .= ", 'thing$i'";
}
add_default($nl, 'hist_fincl1', 'val'=>$val);
add_default($nl, 'hist_fincl2', 'val'=>"''");
add_default($nl, 'aTm_CPL_dt', 'val'=>3600 );
add_default($nl, 'lnd_cpl_DT', 'val'=>$nl->get_value('atm_cpl_dt'));
add_default($nl, 'ice_cpl_DT', 'val'=>$nl->get_value('atm_cpl_dt'));
add_default($nl, 'ocn_cpl_dt', 'val'=>$nl->get_value('atm_cpl_dt'));
add_default($nl, 'complex1', 'val'=>"(-5.2345d+23,+5.98765E-100), (5.2345,5.98765E-33)");
add_default($nl, 'ncdata');
add_default($nl, 'bnd_topo');
add_default($nl, 'absems_data');
add_default($nl, 'orb_eccen', 'val'=>"1.e-12" );

# Test what happens if input variable is unknown
eval { add_default($nl, "xxx" ); };
like( $@, qr/ERROR: variable "xxx" not found in namelist definition file/, 'Test trying to set an unknown name' );

# Test what happens if a default is not set
eval { add_default($nl, "irad" ); };
like( $@, qr/No default value found for irad/, 'Test trying to set an unset default' );

# Make sure can write the namelist
my $outnl = "output_test_namelist.nl";
ok( $nl->write( "$outnl.tmp", note=>"Add note to end" ), "test that can write namelist out" );

my $rval = `diff $outnl $outnl.tmp`;
is($rval,"","Check that output namelist is as expected");

# Make sure can read and validate and merge (a no-op) the output namelist

my $nl2 = Build::Namelist->new("$outnl");
ok( $definition->validate($nl2), "Test validate of resultant namelist" );
add_default($nl2, "case_name", 'val'=>"'newtest'" );
my $val_nl = $definition->validate($nl2);
$nl->merge_nl($val_nl);
ok( $definition->validate($nl2), "Test validate of resultant namelist" );

# Test the merge_nl 
my %cvals = ( "bnd_topo"=>  { "prev"=>undef, set=>"'bnd_topo'"   }, 
              "ncdata"=>    { "prev"=>undef, set=>"'ncdata'"     },
              "co2_type"=>  { "prev"=>undef, set=>"'prognostic'" },
              "start_type"=>{ "prev"=>undef, set=>"'branch'"     },
            );
foreach my $var ( keys(%cvals) ) {
   my $group = $definition->get_group_name($var);
   $cvals{$var}{'prev'} = $nl->get_variable_value($group, $var );
   $nl2->set_variable_value($group, $var, $cvals{$var}{'set'} );
}
# Do a standard merge make sure only the ones NOT set before were changed
$val_nl->merge_nl($nl2);
ok( $definition->validate($val_nl), "Test validate of merged namelist" );
foreach my $var ( keys(%cvals) ) {
   my $group = $definition->get_group_name($var);
   $val = $val_nl->get_variable_value($group, $var );
   if ( ! defined($cvals{$var}{'prev'})) {
      is($val,$cvals{$var}{'set'},"Check that merged namelist item is changed");
   } else {
      is($val,$cvals{$var}{'prev'},"Check that merged namelist item is unchanged");
   }
}
# Do a merge where you overwrite values on a namelist, make sure they all were changed
$val_nl->merge_nl($nl2, overwrite=>1 );
ok( $definition->validate($val_nl), "Test validate of merged overwritten namelist" );
foreach my $var ( keys(%cvals) ) {
   my $group = $definition->get_group_name($var);
   $val = $val_nl->get_variable_value($group, $var );
   is($val,$cvals{$var}{'set'},"Check that all merged namelist item is changed");
}
# Make sure reading in a simple namelist works
my $namelist = 
"&clm_inparm finidat=' ', hist_nhtfrq=1, create_crop_landunit=.false., hist_avgflag_pertape='A', 'I', 'X', hist_fincl1='A1','B2','C3' /";
my $nl2 = Build::Namelist->new( $namelist );
ok( $definition->validate($nl2), "Test validate of namelist read in on command line" );

# Check that validate properly finds problems when things are not properly set
my %vars = ( "aero_carbon"=>{'type'=>"logical" }, "bnd_topo"=>{'type'=>"string"},
             "orb_eccen"=>{'type'=>"real"}, "irad"=>{'type'=>"integer"},
             "complex1"=>{'type'=>"complex"}  );
foreach my $var ( keys(%vars) ) {
   my $group = $definition->get_group_name($var);
   my $prev_val = $nl->get_variable_value($group, $var);
   my %vals = ( 'integer'=>"1", 'real'=>"1.25", 'logical'=>".true.", 'string'=>"'1'",
                'complex'=>"(-5.2345d+23,+5.98765E-100)" );
   foreach my $type ( keys(%vals) ) {
     my $val = $vals{$type};
     if ( $vars{$var}{'type'} eq $type ) { next; }
     if ( $vars{$var}{'type'} eq "real" && $type eq "integer" ) { next; }
     #if ( $vars{$var}{'type'} eq "real" && $type eq "complex" ) { next; }

     $nl->set_variable_value($group, $var, $val);
     eval { my $nl3 = $definition->validate($nl); };
     my $regexpr;
     if ( $vars{$var}{'type'} eq "string" || $type eq "string" ) {
        $regexpr = qr/ERROR\: in _validate_pair \(package Build\:\:NamelistDefinition\)\: Variable name .+ has an input value of type .+/;
     } else {
        $regexpr = qr/ERROR\: in validate_variable_value \(package Build::Namelist\)\: Variable name .+ is not a valid FORTRAN .+/;
     }
     like( $@, $regexpr, "Test validating incorrect value: $var with $type"  );
   }
   $nl->set_variable_value($group, $var, $prev_val);
   my $nl3 = $definition->validate($nl);
}
# Test that can have End-Of-Line between integer arrays

my $var = "npr_yz";
my $group = $definition->get_group_name($var);
my $val = "1\n, 2,\n3,  \n4";
$nl->set_variable_value($group, $var, $val);
ok( $definition->validate($nl), "Test that can have end of line between integer arrays" );

# Test string too long

my $var = "scenario_prognostic_sulfur";
my $group = $definition->get_group_name($var);
my $prev_val = $nl->get_variable_value($group, $var);
my $val = "'string_too_long_for_a_short_variable_only_16_chars_long'";
$nl->set_variable_value($group, $var, $val);
eval { my $nl3 = $definition->validate($nl); };
my $regexpr;
$regexpr = qr/ERROR: in validate_variable_value \(package Build::Namelist\): Variable name .+ has a string element that is too long/;
like( $@, $regexpr, "Test validating string too long: $var"  );
$nl->set_variable_value($group, $var, $prev_val);
my $nl3 = $definition->validate($nl);

# Test that namelist with a mistake with quotes fails
eval { my $nl2 = Build::Namelist->new( "&clm_inparm hist_fincl1='a0', ''A1','B2','C3' /" ); };
my $regexpr;
$regexpr = qr/ERROR\(Build::Namelist::_parse_next\): expect a equal '=' sign instead got: ','B2','C3'/;
like( $@, $regexpr, "Test reading namelist with double quotes mistake fails"  );

# integer array too long
my $var = "npr_yz";
my $group = $definition->get_group_name($var);
my $val = "1,2,3,4,5,6,7";
$nl->set_variable_value($group, $var, $val);
eval { my $nl3 = $definition->validate($nl); };
my $regexpr;
$regexpr = qr/ERROR: in validate_variable_value \(package Build::Namelist\): Variable name .+ has exceeded the dimension size of the array/;
like( $@, $regexpr, "Test validating integer array too long: $var"  );
my $val = "1,2,3,4";
$nl->set_variable_value($group, $var, $val);
my $nl3 = $definition->validate($nl);

# string array too long
my $var = "avgflag_pertape";
my $group = $definition->get_group_name($var);
my $val = "'A','I','X'	,'M',  'M' ,'M', 'I','X'";
$nl->set_variable_value($group, $var, $val);
eval { my $nl3 = $definition->validate($nl); };
my $regexpr;
$regexpr = qr/ERROR: in validate_variable_value \(package Build::Namelist\): Variable name .+ has exceeded the dimension size of the array/;
like( $@, $regexpr, "Test validating string array too long: $var"  );
my $val = "'A','I','X','M','M','M'";
$nl->set_variable_value($group, $var, $val);
my $nl3 = $definition->validate($nl);

# Test get_defined_vars_in_group
# Check a namelist group that has some items:
my ($num_found, $defined_vars) = $nl->get_defined_vars_in_group('clm_inparm');
cmp_ok ($num_found, 'gt', 0, "Test get_defined_vars_in_group for non-empty group: num_found");
cmp_ok (length($defined_vars), 'gt', 0, "Test get_defined_vars_in_group for non-empty group: defined_vars");
# Check a namelist group that appears in the namelist definition, but does not have any items:
my ($num_found, $defined_vars) = $nl->get_defined_vars_in_group('phys_debug_nl');
is ($num_found, 0, "Test get_defined_vars_in_group for empty group: num_found");
is ($defined_vars, '', "Test get_defined_vars_in_group for empty group: defined_vars");
# Check a namelist group that doesn't exist at all:
my ($num_found, $defined_vars) = $nl->get_defined_vars_in_group('XXX_does_not_exist_XXX');
is ($num_found, 0, "Test get_defined_vars_in_group for non-existent group: num_found");
is ($defined_vars, '', "Test get_defined_vars_in_group for non-existent group: defined_vars");


# multi-dim array --- NOTE: THIS TEST MUST BE LAST -- ALWAYS DIES AFTER IT
my $var = "test_multi_dim_array";
my $group = $definition->get_group_name($var);
my $val = "'thing','thing'";
$nl->set_variable_value($group, $var, $val);
eval { my $nl3 = $definition->validate($nl); };
my $regexpr;
$regexpr = qr/Variable name .+ is defined as a multidimensional array -- which is invalid for a namelist/;
like( $@, $regexpr, "Test that multi-dim array is invalid: $var"  );

# Cleanup

system( "/bin/rm $outnl.tmp" );


print "\nSuccessfully ran all tests\n";
