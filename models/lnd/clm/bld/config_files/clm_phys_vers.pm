package config_files::clm_phys_vers;
my $pkg_nm = 'config_files::clm_phys_vers';
#-----------------------------------------------------------------------------------------------
#
# SYNOPSIS
#
# require config_files::clm_phys_vers;
#
# my $phys = config_files::clm_phys_vers->new("clm4_0");
# print $phys->as_float();
# print $phys->as_long();
# print $phys->as_string();
# print $phys->as_filename();
#
# DESCRIPTION
#
# Enter the physics version as a string, with a list of valid versions, and have the ability to convert it to 
# different formats.
#
# COLLABORATORS: None
# 
#-----------------------------------------------------------------------------------------------
#
# Date        Author                  Modification
# 03/06/2014  Erik Kluzek             creation
#
#--------------------------------------------------------------------------------------------

use strict;
use bigint;
#use warnings;
#use diagnostics;

my $major_mask      = 1000000;
my $minor_mask      =    1000;
my @version_strings = (       "clm4_0",                     "clm4_5",      "clm5_0" );
my @version_long    = (  4*$major_mask,  4*$major_mask+5*$minor_mask, 5*$major_mask );

#-------------------------------------------------------------------------------

sub new {
    # Constructor, enter version string as argument
    my $class       = shift;
    my $vers_string = shift;

    my $nm   = "$class\:\:new";
    my $self = {};
    bless($self, $class);
    $self->__validate_vers__( $vers_string );
    $self->{'vers_string'} = $vers_string;
    return( $self );
}

#-------------------------------------------------------------------------------

sub __validate_vers__ {
   # Make sure the version string is a valid one
   my $class       = shift;
   my $vers_string = shift;

   my $found = undef;
   foreach  my $i (0..$#version_strings) {
      if ( $vers_string eq $version_strings[$i] ) {
         $found = 1;
         last;
      }
   }
   if ( ! defined($found) ) {
      die "NOT a valid CLM version: $vers_string\n";
   }
}

#-------------------------------------------------------------------------------

sub as_long {
# Return the physics version as a long
  my $self = shift;
  my $vers = shift;

  if ( ! defined($vers) ) {
     $vers = $self->{'vers_string'};
  } else {
     $self->__validate_vers__( $vers );
  }
  my $phys = undef;
  for( my $i = 0; $i <= $#version_strings; $i++ ) {
     if ( $vers eq $version_strings[$i] ) {
        $phys = $version_long[$i];
        last;
     }
  }
  return( $phys );
}

#-------------------------------------------------------------------------------

sub as_float {
# Return the physics version as a float
  my $self = shift;

  my $long  = $self->as_long();
  my $major = int($long / $major_mask);
  my $minor = int(($long - $major*$major_mask)/ $minor_mask);
  my $rev   =  $long - $major*$major_mask - $minor*$minor_mask;
  {
     no  bigint;
     use bignum;

     my $phys  = $major*1.0 + $minor/10.0   + $rev / 10000.0;
     return( $phys );
  }
}

#-------------------------------------------------------------------------------

sub as_string {
# Return the physics version as a string
  my $self = shift;

  my $phys = $self->{'vers_string'};
  return( $phys );
}

#-------------------------------------------------------------------------------

sub as_filename {
# Return the physics version string with clm4_5 and clm5_0 pointing to the same name
  my $self = shift;

  my $phys = undef;
  if ( $self->as_long() < 5*$major_mask ) {
     $phys = $self->as_string();
  } else {
     $phys = "clm4_5";
  }
  return( $phys );
}

#-----------------------------------------------------------------------------------------------
# Unit testing of above
#-----------------------------------------------------------------------------------------------
if ( ! defined(caller) && $#ARGV == -1 ) {
   package phys_vers_unit_tester;

   require Test::More;
   Test::More->import( );

   plan( tests=>13 );

   sub testit {
      print "unit tester\n";
      my %lastv;
      my @vers_list = ( "clm4_0", "clm4_5", "clm5_0" );
      foreach my $vers ( @vers_list ) {
         my $phys = config_files::clm_phys_vers->new($vers);
         isa_ok($phys, "config_files::clm_phys_vers", "created clm_phys_vers object");
         print "$vers: long: ".$phys->as_long()." float: ".$phys->as_float()." string: ".$phys->as_string()." file: ".$phys->as_filename()."\n";
         if ( exists($lastv{"long"}) ) {
            is( $phys->as_long() > $lastv{'long'}, 1, "Definition of long is not increasing\n" );
         }
         if ( exists($lastv{"float"}) ) {
            is( $phys->as_float() > $lastv{'float'}, 1, "Definition of float is not increasing\n" );
         }
         # Check that also can get results of any valid value for long
         foreach my $chvers ( @vers_list ) {
            my $lvalue = $phys->as_long($chvers);
            print "Long value of $chvers = $lvalue\n";
         }
         # Check that a bad value gives an error
         eval { $phys->as_long('xxx'); };
         like( $@, qr/NOT a valid CLM version:/, "check that a bad version fails" );
         # Save last values to make sure increasing
         $lastv{'long'}   = $phys->as_long();
         $lastv{'float'}  = $phys->as_float();
      }
      my $phys = config_files::clm_phys_vers->new("clm4_0");
      is( 4.0, $phys->as_float(), "Make sure clm4_0 correct float value" );
      $phys = config_files::clm_phys_vers->new("clm4_5");
      no  bigint;
      use bignum;
      is( 4.5, $phys->as_float(), "Make sure clm4_5 correct float value" );
      no bignum;
      use bigint;
      $phys = config_files::clm_phys_vers->new("clm5_0");
      is( 5.0, $phys->as_float(), "Make sure clm5_0 correct float value" );
      print "\nSuccessfully ran all tests\n";
   }
}

#-----------------------------------------------------------------------------------------------
# Determine if you should run the unit test or if this is being called from a require statement
#-----------------------------------------------------------------------------------------------

if ( defined(caller) ) {
   1   # to make use or require happy
} elsif ( $#ARGV == -1 ) {
   &phys_vers_unit_tester::testit();
}
