#=======================================================================
#
# NAME
#
#   Streams::Template - a perl module for reading/writing CCSM stream text files.
#
# SYNOPSIS
#
#   use Streams::Template;
#
#   # create a new stream template object
#   my %opts;
#   # Set the required hash elements given in the documentation on the new method
#   # all have to be set -- but most can be set to blank to indicate not used.
#   $opts{'ProgName'} = $0;
#   $opts{'case'}     = "";
#   $opts{'printing'} = 0;
#   .
#   .
#   .
#   $opts{'csmdata'}  = "/fs/cgd/csm/inputdata";
#
#   my $strms = Streams::Template->new( \%opts );
#
#   # Read in template
#   # NOTE: If the file streams.txt.readme exists in the same directory as the input
#   #       template -- a comment will be added to the top of the output stream text
#   #       file with the information from stream.txt.readme in it.
#   my $inputFile = "InputStreamsTemplateFile.txt";
#   $strms->Read( $inputFile );
#
#   # Write the stream text file out
#   my $outputFile = "OutputStreamsFile.txt";
#   $strms->Write( $outputFile );
#
#   # Test that files exist
#   $strms->TestFilesExist( "data" );
#   $strms->TestFilesExist( "domain" );
#
#   # Get filenames or filepath
#   my $expandEnv = 1;      # If should expand environment variables in get
#   my @filenames = $strms->GetDataFilenames( "data", $expandEnv );
#   my $filepath  = $strms->GetDataFilepath(  "data", $expandEnv );
#
# DESCRIPTION
#
#  This is a perl module to build a stream file
#
# Methods:
#
#     new ---------------- Constructor.
#     Read --------------- Read in input template. Add information from
#                          streams.txt.readme if it exists in the directory
#                          with the template file being used.
#     Write -------------- Write output stream text file based on input Read.
#     TestFilesExist ----- Test that the files exist.
#                          (can either test domain file or data files)
#     GetDataFilenames --- Get the full filepath to all of the filenames.
#                          (can either get domain file or data files)
#     GetDataFilepath ---- Get the directory where data or domain files exist.
#                          (can either get domain filepath or data filepath)
#
# COLLABORATORS
# 
# IO::File
# XML::Lite
#
# HISTORY
#
# Date        Author                  Modification
#----------------------------------------------------------------------------------------
#  2007-Sep   Erik Kluzek             Original version
#  2008-Jul   Erik Kluzek             Add tInterpAlgo
#  2008-Sep   Erik Kluzek             Add last-month option to write
#  2009-Jun   Tony Craig              Add offset
#  2009-Jun   Erik Kluzek             Add year-month-day indicators %ymd
#----------------------------------------------------------------------------------------
#
#=======================================================================
use strict;
use IO::File;

package Streams::Template;

$Streams::Template::name       = "streamstemplate";
$Streams::Template::Tmplt_name = "stream";
#-------------------------------------------------------------------------------

sub new {
#
# Create a Streams Template object
#
#  Required input to opts hash:
#
#   printing ----- Flag if should print out status information (files read etc.).
#   ProgName ----- Name of calling program.
#   ProgDir ------ Directory of calling program.
#   cmdline ------ Command-line give to calling program.
#   filepath ----- Input full filepath to data.
#   yearfirst ---- First year to loop data over.
#   yearlast ----- Last year to loop data over.
#   case --------- Input case name.
#   domain ------- Domain filename that should use.
#   csmdata ------ Directory to where standard CCSM inputdata is.
#   filenames ---- Input filename(s) or filename indicator.
#   datasource --- Data source type of where to get data.
#   res ---------- Input data resolution to use.
#   domainpath --- Directory to domain file.
#
  my $class        = shift;
  my $opts_ref     = shift;

  my %opts = %$opts_ref;

  my $ProgName     = $opts{'ProgName'};
  my $nm = "${ProgName}::new";
  # Error check that input and opts hash has the expected variables
  my @required_list = ( "printing",   "ProgName", "filepath", 
                        "yearfirst",  "yearlast",   "cmdline", 
                        "case",       "domain",   "csmdata",
                        "filenames",  "datasource", "ProgDir",
                        "res",        "domainpath"  );
  foreach my $var ( @required_list ) {
     if ( ! defined($opts{$var}) ) {
        die "ERROR($nm): Required input variable $var was not found\n";
     }
  }
  my $self = {};

  $self->{'ProgName'}   = $opts{'ProgName'};
  $self->{'ProgDir'}    = $opts{'ProgDir'};
  $self->{'cmdline'}    = $opts{'cmdline'};
  $self->{'RESOLUTION'} = $opts{'res'};
  $self->{'printing'}   = $opts{'printing'};
  $self->{'filepath'}   = $opts{'filepath'};
  $self->{'case'}       = $opts{'case'};
  $self->{'yearfirst'}  = $opts{'yearfirst'};
  $self->{'yearlast'}   = $opts{'yearlast'};
  $self->{'csmdata'}    = $opts{'csmdata'};
  $self->{'domain'}     = $opts{'domain'};
  $self->{'domainpath'} = $opts{'domainpath'};
  $self->{'filenames'}  = $opts{'filenames'};
  $self->{'datasource'} = $opts{'datasource'};
  # Set template to undefined so can abort if try to do something before read is done
  $self->{'template'}   = undef;

  bless( $self, $class );
  return( $self );
}

#-------------------------------------------------------------------------------

sub Read {
#
# Read in the stream template file and store in Streams::Template object.
#
# Only read an element if input datasource and RESOLUTION match value in file
# for that element. Will keep reading stream elements and overwrite what was
# stored as long as the attributes match. So the last valid match is what will
# be stored on output.
#
# Output is a reference to a hash contained in the Streams::Template object.
# The domainInfo and fieldInfo elements of that hash are also hashes themselves.
# The variableNames element of domainInfo or fieldInfo is an array reference.
#
# At this point if any input variables are not-blank they will replace the values read
# in.
#
  my $self      = shift;
  my $file      = shift;

  # Initialize some local variables
  my $ProgName  = $self->{'ProgName'};
  my $printing  = $self->{'printing'};
  my $filenames = $self->{'filenames'};
  my $domain    = $self->{'domain'};
  my $path      = $self->{'filepath'};
  my $dompath   = $self->{'domainpath'};
  my $nm        = "${ProgName}::Read";

  $self->{'template'} = $file;

  my %defaults;
  my $match = undef;

  # Open file
  my $xml = XML::Lite->new( $file );
  if ( ! defined($xml) ) {
    die "ERROR($nm): Trouble opening or reading $file\n";
  }
  #
  # Find the root streamstemplate for this file
  #
  my $rootname = "$Streams::Template::name";
  my $elm = $xml->root_element( );
  my @list = $xml->elements_by_name( $rootname );
  if ( $#list < 0 ) {
    die "ERROR($nm): could not find the main $rootname element in $file\n";
  }
  if ( $#list != 0 ) {
    die "ERROR($nm): $rootname root element in $file is duplicated, there should only be one\n";
  }
  $elm = $list[0];
  my @streams = $elm->get_children();
  if ( $#streams < 0 ) {
    die "ERROR($nm): There are no sub-elements to the $rootname element in $file\n";
  }
  #
  # Open and read readme file if it exists
  #
  $file   =~ s!(.*)/!!;
  my $Dir = $1;
  if ( ! defined($Dir) ) { $Dir = "."; }
  my $readmefile = "$Dir/streams.txt.readme";
  if ( -f $readmefile ) {
     print STDERR "($nm) Read template readme file: $readmefile\n" if $printing;
     my $fh = new IO::File;
     $fh->open("<$readmefile") or die "** can't open input stream readme file: " . 
                                      "$readmefile, $!\n";
     my $value = "";
     while( $_ = <$fh> ) {
        $value .= $_;
     }
     close( $fh );
     $defaults{'template_comment'} = $value;
  }
  #
  # Go through each stream template
  #
  foreach my $streams ( @streams ) {
    my $name  = $streams->get_name();
    my %atts  = $streams->get_attributes();
    #
    # Either read in a comment or a streams tag
    #
    if ( $name =~ /^[a-zA-Z0-9]+\_comment$/ ) {
       $defaults{$name} = $streams->get_text();
       next;
    } elsif ( $name ne $Streams::Template::Tmplt_name ) {
       die "ERROR($nm): The only valid sub-elements to the $rootname element can be" .
           " called $Streams::Template::Tmplt_name or *_comment not: $name\n";
    }
    #
    # Skip to next one if doesn't match attributes
    #
    if ( ! $self->__MatchesAtts__( \%atts ) ) { 
       next; 
    }
    #
    # Go through the sub-elements to the stream element
    #
    $match = 1;
    my @children = $streams->get_children();
    foreach my $child ( @children ) {
      # Name of element, and it's associated value
      my $child_name  = $child->get_name();
      my $child_value = $child->get_text();
      my %child_atts  = $child->get_attributes();
      my @keys  = keys( %child_atts );
      if ( $self->__MatchesAtts__( \%child_atts ) ) {
         #
         # If an element with sub-elements go through each of it's children
         #
         my @Grandchildren = $child->get_children();
         if ( $#Grandchildren > 0 ) {
            my $grandchild_defaults_ref;
            my %grandchild_defaults;
            # If already defined just overwrite existing hash
            if ( defined($defaults{$child_name}) ) {
               $grandchild_defaults_ref = $defaults{$child_name};
            } else {
               $grandchild_defaults_ref = \%grandchild_defaults;
            }
            foreach my $Grandchild ( @Grandchildren ) {
               my $name  = $Grandchild->get_name();
               my $value = $Grandchild->get_text();
               my %grandchild_atts  = $Grandchild->get_attributes();
               if ( $self->__MatchesAtts__( \%grandchild_atts ) ) {
                  $value =~  s/^[ \n]+//;          # remove leading spaces
                  $value =~  s/[ \n]+$//;          # remove ending spaces
                  if ( $child_name eq "domainInfo" && $name eq "fileNames" && ($domain ne "") ) { 
                     $value = "$domain";
                  #
                  # For variable-names -- split them up by returns and make sure they have two values
                  #
                  } elsif ( $name eq "variableNames" ) {
                     my @fields = split( /\n/, $value );
                     $value = "";
                     if ( $fields[0] =~ /^[\n ]*$/ ) {
                        shift( @fields );
                     }
                     if ( $fields[$#fields] =~ /^[\n ]*$/ ) {
                        pop( @fields );
                     }
                     for( my $i = 0; $i <= $#fields; $i++ ) {
                        # Check for 2 field values (non-white-space) separated by white-space
                        if ( $fields[$i] =~ /^[ ]*([^ ]+[ ]+[^ ]+)[ ]*$/ ) {
                           # Remove whitespace
                           $fields[$i] = $1;
                        } else {
                           die "ERROR($nm): Variablenames in template file $file NOT two-valued: $fields[$i]\n";
                        }
                     }
                     if ( $#fields < 0 ) {
                       die "ERROR($nm): No variable names -- template file must be bad\n";
                     }
                     $value = \@fields;
                  #
                  # For filepath replace with the input filepath if input path NOT blank
                  #
                  } elsif ( $name eq "filePath") {
                    if (      $child_name eq "fieldInfo"  && ($path    ne "") ) {
                      $value = "$path";
                    } elsif ( $child_name eq "domainInfo" && ($dompath ne ""))  {
                      $value = "$dompath";
                    }
                  #
                  # For file inidicator replace with the input filenames if input filenames NOT blank
                  #
                  } elsif ( $child_name eq "fieldInfo" && $name eq "fileNames" && ($filenames ne "") ) {
                    $value = "$filenames";
                  }
  
                  $$grandchild_defaults_ref{$name} = $value;
               }
            }
            $defaults{$child_name} = $grandchild_defaults_ref;
         } else {
            $child_value =~  s/^[ \n]+//;          # remove leading spaces
            $child_value =~  s/[ \n]+$//;          # remove ending spaces
  
            $defaults{$child_name} = $child_value;
         }
      }
    }
  }
  if ( ! $match ) {
     die "ERROR($nm): did NOT match any of the stream templates in input file template: $file\n";
  }
  $self->{defaults} = \%defaults;
  print STDERR "($nm) Successfully read template: $file\n" if $printing;
}

#-------------------------------------------------------------------------------

sub Write {
#
# Write out a stream text file
#
# This is where % indicators input from read or input from command line are
# substituted before written out to the file.
#
  my $self      = shift;
  my $outfile   = shift;
  my $lastmonth = shift;

  # Initialize some local variables
  my $ProgName     = $self->{'ProgName'};
  my $ProgDir      = $self->{'ProgDir'};
  my $printing     = $self->{'printing'};
  my $template     = $self->{'template'};
  my $nm           = "${ProgName}::Write";
  my $cmdline      = $self->{'cmdline'};

  if ( ! defined($template) ) {
     die "${nm}:: a template has NOT been read in yet -- abort.\n";
  }

  my $defaults_ref = $self->{'defaults'};
  my %defaults     = %$defaults_ref;
  #
  # Open output file
  #
  my $fh = new IO::File;
  if ( $outfile ne "" ) {
     $fh->open(">$outfile") or die "** can't open output stream file: $outfile, $!\n";
     print STDERR "($nm) Open output stream: $outfile\n" if $printing;
  } else {
    $fh->fdopen(fileno(STDOUT),"w") or die "** can't open STDOUT, $!\n";
  }
  #
  # Write any comments out to the top of the file
  #
  print $fh "<$Streams::Template::name>\n";
  my @keys = keys( %defaults );
  foreach my $com ( @keys ) {
     if ( $com =~ /^[a-zA-Z0-9]+_comment$/ ) {
        my $value = $defaults{$com};
        $value =~  s/^[ \n]+//;          # remove leading spaces
        $value =~  s/[ \n]+$//;          # remove ending spaces
        $value =~  s/\n[ ]+/\n/g;        # remove spaces at beginning of lines
        $self->__WriteValue__( $fh, $com, $value, 1 );
     }
  }
  # Begin stream
  print $fh "<$Streams::Template::Tmplt_name>\n";
  #
  # Write out comment and datasource
  #
  foreach my $i ( ("comment", "dataSource" ) ) {
     $self->__WriteValue__( $fh, $i, $defaults{$i}, 1 );
  }
  #
  # Write out the domainInfo and fieldInfo elements
  #
  my %comment;
  my %keys = ( domainInfo=>["variableNames", "filePath", "fileNames"],
               fieldInfo =>["variableNames", "filePath", "tInterpAlgo", "offset", "fileNames"] );
  foreach my $info ( ( "domainInfo", "fieldInfo" ) ) {
     my $level = 2;
     if ( defined( $defaults{$info} ) ) {
        my $spacing = $self->__Spacing__( 1 );
        print $fh     "$spacing<$info>\n";
        my $Info_ref = $defaults{$info};
        if ( ref($Info_ref) ne "HASH" ) {
           die "${nm}:: $info is NOT a hash -- something must be wrong from the Read\n";
        }
        my %Info     = %$Info_ref;
        my $keys_ref = $keys{$info};
        foreach my $name ( @$keys_ref ) {
           if ( ref($Info{$name}) eq "ARRAY" || ($name =~ /fileNames/) ) {
              $self->__WriteArray__( $fh, $name, $Info_ref, $level, $lastmonth );
           } else {
              # Make sure defined -- unless is tInterpAlgo or offset
              if ( defined($Info{$name}) || ($name ne "tInterpAlgo" && $name ne "offset")) {
                my $value = $self->__Sub__( $Info_ref, $name );
                $self->__WriteValue__( $fh, $name, $value, $level );
              }
           }
        }
        print $fh "$spacing</$info>\n";
     } else {
        die "${nm}:: $info was NOT defined from the Read of template: $template\n";
     }
  }
  #
  # Write out information on the program that created this file
  #
  my $spacing = $self->__Spacing__( 1 );
  print $fh "$spacing<!-- Information on the program that created this file -->\n";
  my $description="This CCSM stream text file was created by $ProgName using the command line:\n" . 
                  "      $ProgDir/$ProgName $cmdline   \n" .
                  "For more information on $ProgName:  \n".
                  "    $ProgDir/$ProgName -help          ";
  $self->__WriteValue__( $fh, "${ProgName}_documentation",  "$description", 1 );
  # Close stream and file
  print $fh "</$Streams::Template::Tmplt_name>\n";
  print $fh "</$Streams::Template::name>\n";
  
  close( $fh );
  print STDERR "($nm) Successfully created stream file: $outfile\n" if $printing;
}

#-------------------------------------------------------------------------------

sub TestFilesExist {
#
# Test that the files indicated in the stream exist
#
  my $self = shift;
  my $type = shift;  # Either data or domain

  # Initialize some local variables
  my $ProgName     = $self->{'ProgName'};
  my $printing     = $self->{'printing'};
  my $nm           = "${ProgName}::TestFilesExist";

  if ( ! defined($self->{'template'}) ) {
     die "${nm}:: a template has NOT been read in yet -- abort.\n";
  }
  #
  # Get hash to work on
  #
  my $defaults_ref = $self->{'defaults'};
  my %defaults     = %$defaults_ref;
  my $info;
  my $key = "fileNames";
  if (      $type eq "data" ) {
     $info = "fieldInfo";
  } elsif ( $type eq "domain" ) {
     $info = "domainInfo";
  } else {
     die "${nm}:: bad input type to method: $type should be data or domain\n";
  }
  my $Info_ref = $defaults{$info};
  if ( ref($Info_ref) ne "HASH" ) {
     die "${nm}:: $key is NOT a hash -- something must have went wrong in the Read\n";
  }
  #
  # Get filepath
  #
  my $expandEnv = 1;
  my $filepath  = $self->GetDataFilepath( $type, $expandEnv );

  #
  # Get local path -- if cp: or null: option set
  #
  if ( $filepath =~ /^[ ]*(cp:|null:|)([^:]+)[ ]*$/ ) {
     $filepath = $2;
  } else {
    die "${nm}:: can not test that files exist when filepath is not a local path: path=$filepath\n";
  }

  my $fileNames_ref = $self->__Sub__( $Info_ref, $key );
  my @filenames = @$fileNames_ref;
  foreach my $file ( @filenames ) {
     print STDERR "Test that the file $filepath/$file exists\n" if $printing;
     if ( ! -f "$filepath/$file" ) {
        die "${nm}:: local file $filepath/$file does NOT exist -- aborting\n";
     }
  }
}

#-------------------------------------------------------------------------------

sub GetDataFilenames {
#
# Just get the data filenames
#
  my $self      = shift;
  my $type      = shift;
  my $expandEnv = shift;
  my $lastmonth = shift;

  my $class        = ref( $self );
  my $nm           = "${class}::GetDataFilenames";

  if ( ! defined($self->{'template'}) ) {
     die "${nm}:: a template has NOT been read in yet -- abort.\n";
  }
  my $defaults_ref = $self->{'defaults'};
  my %defaults     = %$defaults_ref;
  #
  # Get path
  #
  my $filepath = $self->GetDataFilepath( $type, $expandEnv );
  my $key = "fileNames";
  my $info;
  if (      $type eq "data" ) {
     $info = "fieldInfo";
  } else {
     $info = "domainInfo";
  }
  #
  # Get filenames as an array and add path to them
  #
  my $Info_ref = $defaults{$info};
  if ( ref($Info_ref) ne "HASH" ) {
     die "${nm}:: $info is NOT a hash -- something must have went wrong in the Read\n";
  }
  my $fileNames_ref = $self->__Sub__( $Info_ref, $key, $lastmonth );
  my @files;
  foreach my $file ( @$fileNames_ref ) {
     push @files, "$filepath/$file";
  }
  return( @files );
}

#-------------------------------------------------------------------------------

sub GetDataFilepath {
#
# Just get the data filepath
#
  my $self      = shift;
  my $type      = shift;
  my $expandEnv = shift;

  my $class        = ref( $self );
  my $nm           = "${class}::GetDataFilepath";

  if ( ! defined($self->{'template'}) ) {
     die "${nm}:: a template has NOT been read in yet -- abort.\n";
  }
  my $defaults_ref = $self->{'defaults'};
  my %defaults     = %$defaults_ref;
  my $key;
  if (      $type eq "data" ) {
     $key = "fieldInfo";
  } elsif ( $type eq "domain" ) {
     $key = "domainInfo";
  } else {
     die "${nm}:: bad input type to method: $type should be data or domain\n";
  }
  my $Info_ref = $defaults{$key};
  if ( ref($Info_ref) ne "HASH" ) {
     die "${nm}:: $key is NOT a hash -- something must have went wrong in the Read\n";
  }
  my %Info = %$Info_ref;
  my $filepath = $self->__Sub__( $Info_ref, 'filePath' );
  #
  # Replace any env variables
  #
  if ( $expandEnv ) {
     while( $filepath =~ m/\$/ ) {
        if ( $filepath =~ m/(^[^\$]*)\$([^ \/]+)(.*$)/ ) {
           my $startpath = $1;
           my $envvarnm  = $2;
           my $endpath   = $3;
           if ( ! defined($ENV{$envvarnm}) ) {
              die "${nm}:: ENV variable $envvarnm is in pathname ($filepath) -- but NOT defined\n";
           }
           my $envvar = $ENV{$envvarnm};
           $filepath     = "${startpath}${envvar}${endpath}";
        } else {
           die "${nm}:: malformed ENV variable in pathname -- $filepath\n";
        }
     }
     if ( $filepath =~ m/(^.*)(DIN_LOC_ROOT)(.*$)/ ) {
        my $startpath = $1;
        my $envvarnm  = $2;
        my $endpath   = $3;
        if ( ! defined($ENV{$envvarnm}) ) {
           die "${nm}:: ENV variable $envvarnm is in pathname ($filepath) -- but NOT defined\n";
        }
        my $envvar = $ENV{$envvarnm};
        $filepath     = "${startpath}${envvar}${endpath}";
     }
  }

  return( $filepath );
}

#-------------------------------------------------------------------------------
# Private methods
#-------------------------------------------------------------------------------

sub __MatchesAtts__ {
#
# If this element matches the input attributes
#
  my $self = shift;
  my $atts_ref = shift;

  my $ProgName  = $self->{'ProgName'};
  my $nm        = "${ProgName}::__MatchesAtts__";
  my %atts = %$atts_ref;
  
  my @attributes = ( "datasource", "RESOLUTION" );
  foreach my $key ( @attributes ) {
     my $value= $self->{$key};
     if ( defined($atts{$key}) ) {
        if ( $key =~ /datasource|RESOLUTION/ && ($value eq "") ) {
            die "ERROR($nm): $key NOT set and is required -- set it on the command line\n";
        }
        if ($atts{$key} ne $value ) {
          return( 0 );
        }
     }
  }
  return( 1 );
}

#-------------------------------------------------------------------------------

sub __Sub__ {
#
# Substitute indicators with given values
#
# Replace any instance of the following substring indicators with the appropriate values:
#
#     %y    = year from the range yearfirst to yearlast
#     %ym   = year-month from the range yearfirst to yearlast with all 12 months
#     %ymd  = year-month-day from the range yearfirst to yearlast with all 12 months
#     %p    = input path
#     %d    = location of csmdata (inputdata for all CCSM models)
#     %dp   = domain file path
#     %do   = name of domain filename
#     %c    = input casename
#
# fileNames or anything with %y, %ym or %ymd will be returned as a reference to an array.
# Everything else is returned as a scalar.
#
   my $self      = shift;
   my $Info_ref  = shift;
   my $name      = shift;
   my $lastmonth = shift;

   my $ProgName = $self->{'ProgName'};
   my $nm       = "${ProgName}::__Sub__";

   my $value = $$Info_ref{$name};
   $value =~  s/^[ \n]+//;          # remove leading spaces
   $value =~  s/[ \n]+$//;          # remove ending spaces
 
   #
   # Replace % indicators appropriately
   #
   my $path       = $self->{'filepath'};
   my $case       = $self->{'case'};
   my $csmdata    = $self->{'csmdata'};
   my $domain     = $self->{'domain'};
   my $domainpath = $self->{'domainpath'};
   my $yearfirst  = $self->{'yearfirst'};
   my $yearlast   = $self->{'yearlast'};
   my %indic = ( "%p" =>{ replace=>"$path",        name=>"filepath"   }, 
                 "%c" =>{ replace=>"$case",        name=>"case"       },
                 "%do"=>{ replace=>"$domain",      name=>"domain"     },
                 "%dp"=>{ replace=>"$domainpath",  name=>"domainpath" },
                 "%d" =>{ replace=>"$csmdata",     name=>"csmdata"    } );
   foreach my $key ( keys(%indic) ) {
      my $val_ref = $indic{$key};
      my %vals    = %$val_ref;
      # If find an indicator in a field value -- replace it with the appropriate thing
      if ( $value =~ /^([^%]*)$key([^ ]*)[ ]*$/ ) {
         my $replace = $vals{'replace'};
         if ( $replace ne "" ) {
           $value = "${1}${replace}${2}";
         } else {
           my $name = $vals{'name'};
           die "${nm}:: $name was NOT defined on command line and needs to be set\n";
         }
      }
   }
   #
   # If year or year/month indicators exist
   #
   if ( $value =~ /%([1-9]*)y([m]?)([d]?)/ ) {
      my $digits = 4;
      if ( $1 ne "" ) { $digits = $2; }
      my $months = 1;
      if ( $2 eq "" ) { $months = 0; }
      my $days   = 0;
      if ( $3 ne "" ) { 
          if ( ! $months ) {
             die "${nm}:: Months NOT defined but days are? (\%yd is NOT valid indicator)\n";
          }
          $months = 0;
          $days   = 1; 
      }
      if ( ($yearfirst < 0)  || ($yearlast < 0) ) {
         die "${nm}:: yearfirst and yearlast  was NOT defined on command line and needs to be set\n";
      }
      #
      # Loop over year range
      #
      my @filenames;
      my $startfilename;
      my $endfilename;
      #
      # Include previous December if %ym form and lastmonth is true
      #
      if ( $lastmonth && $months ) {
          my $year = $yearfirst-1;
          my $filename = $value;
          my $month = 12;
          if ( $filename =~ /^([^%]*)%[1-9]?ym([^ ]*)$/ ) {
             $startfilename = $1;
             $endfilename   = $2;
             $filename = sprintf "%s%${digits}.${digits}d-%2.2d%s", $startfilename, 
                                 $year, $month, $endfilename;
             push @filenames, $filename;
          }
      }
      #
      # Include previous December/31 if %ymd form and lastmonth is true
      #
      if ( $lastmonth && $months ) {
          my $year = $yearfirst-1;
          my $filename = $value;
          my $month = 12;
          my $day   = 31;
          if ( $filename =~ /^([^%]*)%[1-9]?ymd([^ ]*)$/ ) {
             $startfilename = $1;
             $endfilename   = $2;
             $filename = sprintf "%s%${digits}.${digits}d-%2.2d-%2.2d%s", $startfilename, 
                                 $year, $month, $day, $endfilename;
             push @filenames, $filename;
          }
      }
      for ( my $year = $yearfirst; $year <= $yearlast; $year++ ) {
         #
         # If include year and months AND days
         #
         if ( $days ) {
            for ( my $month = 1; $month <= 12; $month++ ) {
               my $dpm = $self->__DaysPerMonth__( $month, $year );
               for ( my $day = 1; $day <= $dpm; $day++ ) {
                  my $filename = $value;
                  if ( $filename =~ /^([^%]*)%[1-9]?ymd([^ ]*)$/ ) {
                     $startfilename = $1;
                     $endfilename   = $2;
                     $filename = sprintf "%s%${digits}.${digits}d-%2.2d-%2.2d%s", 
                                    $startfilename, $year, $month, $day, $endfilename;
                     push @filenames, $filename;
                  }
               }
            }
         #
         # If include year and months
         #
         } elsif ( $months ) {
            for ( my $month = 1; $month <= 12; $month++ ) {
               my $filename = $value;
               if ( $filename =~ /^([^%]*)%[1-9]?ym([^ ]*)$/ ) {
                  $startfilename = $1;
                  $endfilename   = $2;
                  $filename = sprintf "%s%${digits}.${digits}d-%2.2d%s", $startfilename, 
                                        $year, $month, $endfilename;
                  push @filenames, $filename;
               }
            }
         #
         # If just years
         #
         } else {
            my $filename = $value;
            if ( $filename =~ /^([^%]*)%[1-9]?y([^ ]*)$/ ) {
               $startfilename = $1;
               $endfilename   = $2;
               $filename = sprintf "%s%${digits}.${digits}d%s", $startfilename, $year, 
                             $endfilename;
               push @filenames, $filename;
            }
         }
      }
      if ( $#filenames < 0 ) { 
         die "${nm}:: No output filenames -- must be something wrong in template or input filename indicator\n";
      }
      return( \@filenames );
   #
   # If fileNames then return as an array
   #
   } elsif( $name =~ /fileNames/ ) {
      my @filenames = split( /\n/, $value );
      $value = "";
      if ( $filenames[0] =~ /^[\n ]*$/ ) {
         shift( @filenames );
      }
      if ( $filenames[$#filenames] =~ /^[\n ]*$/ ) {
         pop( @filenames);
      }
      for( my $i = 0; $i <= $#filenames; $i++ ) {
         # Remove any whitespace before or after
         $filenames[$i] =~  s/^[ \n]+//;          # remove leading spaces
         $filenames[$i] =~  s/[ \n]+$//;          # remove ending spaces
      }
      if ( $#filenames < 0 ) { 
         die "${nm}:: No output filenames -- must be something wrong in template or input filename indicator\n";
      }
      return( \@filenames );
   #
   # Otherwise return a scalar value
   #
   } else {
      return( "$value" );
   }
}


#-------------------------------------------------------------------------------

sub __DaysPerMonth__ {
#
# Return the number of days per month for a given month
# (and in general year -- but right now just do a noleap calendar of 365 days/year)
#
  my $self  = shift;
  my $month = shift;
  my $year  = shift;

  my $ProgName = $self->{'ProgName'};
  my $nm = "${ProgName}::__DaysPerMonth___";
  my @dpm = ( 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 );
  if ( $month < 1 || $month > 12 ) {
     die "${nm}:: Input month is NOT valid = $month\n";
  }
  my $days = $dpm[$month-1];
  return( $days );
}

#-------------------------------------------------------------------------------

sub __WriteValue__ {
#
# Write a single value out
#
  my $self  = shift;
  my $fh    = shift;
  my $name  = shift;
  my $value = shift;
  my $level = shift;

  my $ProgName = $self->{'ProgName'};
  my $template = $self->{'file'};
  my $nm = "${ProgName}::__WriteValue__";

  my $spacing     = $self->__Spacing__( $level );
  my $val_spacing = $self->__Spacing__( $level+1 );
  if ( defined( $value ) && ($value =~ /\S/) ) {
     $value =~ s/\n/\n${val_spacing}/g;
     print $fh "$spacing<$name>\n${val_spacing}${value}\n${spacing}</$name>\n";
  } else {
     my $datasource = $self->{'datasource'};
     die "${nm}:: $name was NOT defined in the Read of template $template was datasource: " .
         " $datasource correct?\n";
  }
}

#-------------------------------------------------------------------------------

sub __Spacing__ {
#
# Figure out depth of spacing for writing output depending on the level of nesting.
#
  my $self  = shift;
  my $level = shift;

  if ( ! defined($level) ) { $level = 1; }
  my $spacing = "   ";
  for( my $i=$level; $i>0; $i-- ) {
      $spacing .= "   ";
  }
  return( $spacing );
}

#-------------------------------------------------------------------------------

sub __WriteArray__ {
#
# Write out an array
#
  my $self      = shift;
  my $fh        = shift;
  my $name      = shift;
  my $Info_ref  = shift;
  my $level     = shift;
  my $lastmonth = shift;

  my $spacing = $self->__Spacing__( $level );
  # Initialize some local variables
  my $ProgName     = $self->{'ProgName'};
  my $nm           = "${ProgName}::__WriteArray__";
  my $template     = $self->{'file'};
  print $fh "${spacing}<$name>\n";
  my $array_ref;
  if ( $name =~ /fileNames/ ) {
     $array_ref = $self->__Sub__( $Info_ref, $name, $lastmonth );
  } else {
     $array_ref = $$Info_ref{$name};
  }
  my $val_spacing = $self->__Spacing__( $level+1 );
  foreach my $value ( @$array_ref ) {
     print $fh "${val_spacing}$value\n";
  }
  print $fh "$spacing</$name>\n";
}

#-------------------------------------------------------------------------------

1 # to make use or require happy
