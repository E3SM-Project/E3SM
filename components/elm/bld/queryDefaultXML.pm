#=======================================================================
#
#  This is a perl module to read in a list of namelist_default files.
#
#=======================================================================
use strict;
use Build::Config;
use Build::NamelistDefinition;
use Build::NamelistDefaults;
use Build::Namelist;

package queryDefaultXML;

#-------------------------------------------------------------------------------

sub read_cfg_file
#
# Read in the configuration cache XML file on the build-time configuration
#
{
    my ($file, $empty_cfg_file, $printing, $settings_ref) = @_;

    my $cfg;
    my %config;
    if ( $file eq "noconfig" ) {
       print "No configuration cache file to read in.\n" if $printing;
       $cfg = Build::Config->new( $empty_cfg_file );
    } elsif ( -f "$file" ) {
       $cfg = Build::Config->new($file);
    } else {
       die "Bad filename entered: $file does NOT exist or can not open it.\n";
    }
    #
    # Make sure variables are set to valid values
    #
    foreach my $key ( keys( %config ) ) {
       if ( $cfg->is_valid_name( $key ) ) {
          $cfg->set( $key, $config{$key} );
       }
    }
    foreach my $key ( $cfg->get_names( ) ) {
       if ( defined($$settings_ref{$key}) ) {
          if ( $cfg->is_valid_name( $key ) ) {
             $cfg->set( $key, $$settings_ref{$key} );
          }
       }
    }
    return( $cfg );
}

#-------------------------------------------------------------------------------

sub ReadDefaultXMLFile {
#
# Read in the default XML file for the default namelist settings
#
  my $opts_ref     = shift;
  my $settings_ref = shift;

  # Error check that input and opts hash has the expected variables
  my $ProgName     = $$opts_ref{'ProgName'};
  my $nm = "${ProgName}::ReadDefaultXMLFile";
  my @required_list = ( "files", "nldef_files", "empty_cfg_file", "config", "namelist", 
                        "csmdata", "hgrid", "printing", "ProgName", "cmdline",      
                        "cfgdir"  );
  foreach my $var ( @required_list ) {
     if ( ! defined($$opts_ref{$var}) ) {
        die "ERROR($nm): Required input variable $var was not found\n";
     }
  }
  my $printing = $$opts_ref{'printing'};
  my $cmdline  = $$opts_ref{'cmdline'};
  # Initialize some local variables
  my $files_ref          = $$opts_ref{'files'};
  my @files              = @$files_ref;
  my $nldef_ref          = $$opts_ref{'nldef_files'};
  my @nl_definition_files= @$nldef_ref;
  my $empty_config_file  = $$opts_ref{'empty_cfg_file'};
  my $namelist           = $$opts_ref{'namelist'};

  my $cfg = read_cfg_file( $$opts_ref{'config'}, $$opts_ref{'empty_cfg_file'},
                           $printing, $settings_ref );

  #
  # Set up options to send to namelist defaults object
  #
  my %nlopts;
  foreach my $var ( keys( %$settings_ref) ) {
     if ( $var ne "csmdata" ) {
        $nlopts{$var} = $$settings_ref{$var};
     }
  }
  if ( $$opts_ref{'hgrid'} ne "any" ) {
     $nlopts{'hgrid'} = $$opts_ref{'hgrid'};
  }
  #
  # Loop through all variables in files
  #
  print "($nm) Read: $files[0]\n" if $printing;
  my %defaults;
  my $nldefaults = Build::NamelistDefaults->new($files[0], $cfg);
  for ( my $i = 1; $i <= $#files; $i++ ) {
      print "($nm) Read: $files[$i]\n" if $printing;
      $nldefaults->add( $files[$i] );
  }
  my $definition = Build::NamelistDefinition->new( $nl_definition_files[0] );
  for ( my $i = 1; $i <= $#nl_definition_files; $i++ ) {
     print "($nm) Read: $nl_definition_files[$i]\n" if $printing;
     $definition->add(  $nl_definition_files[$i] );
  }
  if ( $$opts_ref{'csmdata'} eq "default" ) {
    $$opts_ref{'csmdata'} = $nldefaults->get_value( "csmdata", \%nlopts );
  } 
  $nlopts{'csmdata'} = $$opts_ref{'csmdata'};
  foreach my $name ( $nldefaults->get_variable_names() ) {
    my $value   = $nldefaults->get_value( $name, \%nlopts );
    if ( $value eq "null" ) { next; }
    if ( defined($$settings_ref{'var'}) ) {
       if ( $name ne $$settings_ref{'var'} ) { next; }
    }
    $value =~ s/\n//g;
    my $isafile = 0;
    if ( $definition->is_input_pathname($name) ) {

       if ( defined($$settings_ref{'clm_usr_name'}) ) {
          $value   = $nldefaults->get_usr_file( $name, $definition, \%nlopts );
       } 
       if ( $value && ($value !~ /^\/.+$/) ) {
          $value   = $$opts_ref{'csmdata'} . "/" . $value;
       }
       $isafile = 1;
    }
    my $isadir  = 0;
    my $isastr  = 0;
    if (  $definition->get_str_len($name) > 0 ) {
       $isastr = 1;
    }
    #
    # If is a directory (is a file and csmdata or a var with dir in name)
    #
    if ( $isafile  && (($name eq "csmdata") || ($name =~ /dir/)) ) {
      if ( $name eq "csmdata" ) {
         $value = $$opts_ref{'csmdata'};
         $isadir = 1;
      } else {
         $isadir = 1;
      }
    }
    # Return hash with the results
    my $group = $definition->get_group_name( $name );
    if ( $group eq $namelist && $value && (! exists($defaults{$name}{'value'}))  ) { 
       $defaults{$name}{'value'}  = $value;
       $defaults{$name}{'isfile'} = $isafile;
       $defaults{$name}{'isdir'}  = $isadir;
       $defaults{$name}{'isstr'}  = $isastr;
    }
  }
  return( \%defaults );
}

1 # To make use or require happy
