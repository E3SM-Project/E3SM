package Build::NamelistDefaults;
my $pkg = 'Build::NamelistDefaults';
#-----------------------------------------------------------------------------------------------
#
# SYNOPSIS
#
#   use Build::Config;
#   use Build::NamelistDefaults;
#
#   # Create an object containing the configuration from a previous run of configure
#   # (the config_cache.xml file written by configure contains all the parameters needed
#   # to build an executable).
#   my $cfg = Build::Config->new('config_cache.xml');
#
#   # Create a namelist defaults object (read the namelist defaults XML file).
#   my $defaults = Build::NamelistDefaults->new("namelist_defaults.xml", $cfg);
#
#   # Get the default value for the specified namelist variable
#   my $nl_var_default = $defaults->get_value('nl_variable');
#
#   # To add additional defaults from a separate file
#   $defaults->add('namelist_defaults.xml',$cfg);
#
# DESCRIPTION
#
# Build::NamelistDefaults objects are used to represent the default values of namelist
# variables that are stored in an XML file.  Default values may depend on a number of
# attributes that are listed in the XML file.  For example, default values may depend
# on specific parameters that were specified when building the executable code, such
# as the horizontal grid resolution.

#
# METHODS
#
# new() Reads xml file that contains the namelist defaults, and stores
#       A pointer to the object that contains the configuration information.
#
# add() Adds additional defaults from a separate file.
#
# get_value()
#
#
#
# COLLABORATORS
#
# IO::File
# XML::Lite
# Build::Config
#-----------------------------------------------------------------------------------------------
#
# Date        Author                  Modification
#-----------------------------------------------------------------------------------------------
# 2007-Sep    Brian Eaton             Original version
#-----------------------------------------------------------------------------------------------

use strict;
#use warnings;
#use diagnostics;

use XML::Lite;
use Build::Namelist;

#-----------------------------------------------------------------------------------------------
# Public methods
#-----------------------------------------------------------------------------------------------

sub new
{
    my $class = shift;
    my $defaults_filepath = shift;    # the xml file containing the defaults
    my $cfg               = shift;    # a reference to a config object

    # bless the object here so the initialization has access to object methods
    my $nl_defaults = {};
    bless( $nl_defaults, $class );

    # Add the filepath of the defaults file to the object attributes (this is for error reporting).
    $nl_defaults->{'defaults_filepath'} = $defaults_filepath;

    # Add the configuration object to the object attributes
    $nl_defaults->{'cfg_ref'} = $cfg;

    # Initialize the object with the namelist defaults.
    $nl_defaults->_initialize($defaults_filepath);

    return $nl_defaults;
}

#-----------------------------------------------------------------------------------------------

sub add
{
    my $self        = shift;
    my $defaults_filepath = shift;    # the xml file containing the defaults

    # Append the filepath of the defaults file.
    $self->{'defaults_filepath'} .= ", " . $defaults_filepath;

    # Add the additional namelist defaults.
    $self->_initialize($defaults_filepath);
}

#-----------------------------------------------------------------------------------------------

sub get_value
{
# Return a default value for the requested namelist variable.
# Return undef if no default found.

    my $self        = shift;
    my $var_name    = lc(shift);   # name of namelist variable (case insensitive interface)
    my $usr_att_ref = shift;       # reference to hash containing user supplied attributes

    # convenience variables
    my $cfg = $self->{'cfg_ref'};  # configuration object
    my %usr_att = ();              # hash of user supplied attributes

    if (defined $usr_att_ref) {
       ref($usr_att_ref) eq "HASH" or die
	   "ERROR: $pkg\:\:get_value -- user attributes arg must be a hash reference\n";
       %usr_att = %$usr_att_ref;
    }

    # Return now if this variable name is not in the list of names
    my $names_ref = $self->{'names'};  # list of names
    if ( ! exists($$names_ref{$var_name}) ) {
       return( undef );
    }

    # get all the elements that contain defaults for the requested namelist variable
    my $all_elements_ref = $self->{'elements'};  # list of elements
    my $elements_ref     = $$all_elements_ref{$var_name};
    my @elements         = @$elements_ref;

    # examine the attributes of each element to determine the "best fit"
    # keep track of the number of attributes that match the configuration
    my @fit = ();
    ELEMENT: for (my $i = 0; $i <= $#elements; $i++) {

	my $e = $elements[$i];
	my $matches = 0;

	# extract the element attributes
	my %attributes = $e->get_attributes();

	# Check each attribute; first against the configuration, then against values supplied
	# via the optional argument.
	# If an attribute doesn't match, then eliminate the element from further consideration
	foreach my $att_name (keys %attributes) {

	    # Is the attribute part of the configuration?
	    if ($cfg->is_valid_name($att_name)) {

		# Get the value for the attribute from the configuration
		my $cfg_val = $cfg->get($att_name);

		# Check for a match.  If the attributes don't match then skip this element
		# and move to the next after recording the no-match status
		if ($attributes{$att_name} eq $cfg_val) {
		    $matches++;
		}
		else {
		    $fit[$i] = -1;
		    next ELEMENT;
		}

	    }
	    else {

		# If the attribute isn't part of the configuration then do addition
		# checks here.
		#
		# Start with attributes that require special handling...

		if ($att_name eq "ic_ymd") {

		    # Has this attribute been supplied by the user?
		    if (defined $usr_att{'ic_ymd'}) {

			# Check for match (numeric)
			if ($attributes{$att_name} == $usr_att{'ic_ymd'}) {
			    $matches++;
			}
			else {
			    $fit[$i] = -1;
			    next ELEMENT;
			}
		    }
		    # Did user specify that only the month/day needs to match?
		    elsif (defined $usr_att{'ic_md'}) {

			# Check for match (numeric) against month/day part of ic_ymd
			my $ic_md = $attributes{$att_name} % 10000;
			if ($ic_md == $usr_att{'ic_md'} % 10000) {
			    $matches++;
			}
			else {
			    $fit[$i] = -1;
			    next ELEMENT;
			}

		    }

		}
		# Continue from here with generic checking of user supplied attributes
		else {
		    # Has this attribute been supplied by the user?
		    if (defined $usr_att{$att_name}) {

			# Check for a match.  If the attributes don't match then skip this element
			# and move to the next after recording the no-match status
			if ($attributes{$att_name} eq $usr_att{$att_name}) {
			    $matches++;
			}
			else {
			    $fit[$i] = -1;
			    next ELEMENT;
			}
		    }

		} # Finished checking user specified attributes

	    } # Finished attribute checks

	} # Finished loop over attributes

	# At this point the attribute checking has been successful.  Record the matches.
	$fit[$i] = $matches;

    } # Finished loop over elements in defaults file.

    # All elements have been examined.  Return the value from the best fit.  That's the
    # index of the max value of @fit.  In case of a tie it's the first one found.
    my $max_val = $fit[0];
    my $max_idx = 0;
    for (my $i = 1; $i <= $#elements; $i++) {
	if ($fit[$i] > $max_val) {
	    $max_val = $fit[$i];
	    $max_idx = $i;
	}
    }

    # If "best fit" is $max_val = -1, then no match was found.
    if ($max_val >= 0) {
	return $elements[$max_idx]->get_content();
    }
    else {
	return undef;
    }

}

#-----------------------------------------------------------------------------------------------

sub get_usr_file
{
# Return a default value for the requested namelist user file
# Return undef if no default found.
# Set the option 'notest'=>1 if you do NOT want to test the files for existance.

    my $self        = shift;
    my $var_name    = lc(shift);   # name of namelist variable (case insensitive interface)
    my $definition  = shift;
    my $usr_att_ref = shift;       # reference to hash containing user supplied attributes
    my $nm          = "get_usr_file";

    # convenience variables
    my $cfg = $self->{'cfg_ref'};  # configuration object
    my %usr_att = ();              # hash of user supplied attributes

    if (defined $usr_att_ref) {
       ref($usr_att_ref) eq "HASH" or die
	   "ERROR: $pkg\:\:$nm -- user attributes arg must be a hash reference\n";
       %usr_att = %$usr_att_ref;
    }

    # Return now if this variable name is not in the list of names
    my $names_ref = $self->{'names'};  # list of names
    if ( ! exists($$names_ref{$var_name}) ) {
       return( undef );
    }
    # Die if the definition variable entered isn't a definition object
    my $class = ref($definition);
    if ( $class ne "Build::NamelistDefinition" ) {
       die "ERROR: $pkg\:\:$nm -- definition object entered is NOT a Build::NamelistDefinition object\n";
    }
    # Return now if this variable name is not a file-type
    if ( ! defined($definition->is_input_pathname( $var_name )) ) {
       return( undef );
    }
    # get all the elements that contain defaults for the requested namelist variable
    my $all_elements_ref = $self->{'elements'};  # list of elements
    my $elements_ref     = $$all_elements_ref{$var_name};
    my @elements         = @$elements_ref;

    # examine the attributes of each element to determine the "best fit"
    # keep track of the number of attributes that match the configuration
    my @fit    = ();
    my @values = ();
    ELEMENT: for (my $i = 0; $i <= $#elements; $i++) {

	my $e = $elements[$i];
	my $matches = 0;

        #
        # Get the value and then expand any variables contained in it
        #
	my $val = $e->get_content();
        while( $val =~ /^([^\$]*)\$\{([^\}]+)\}(.*)$/ ) {
           my $first = $1;
           my $var   = $2;
           my $last  = $3;
           if ( exists($usr_att{$var}) ) {
              $val = "$first" . $usr_att{$var} . "$last";
           } else {
              die "ERROR: $pkg\:\:$nm -- variable $var NOT sent in for variable $val\n";
           }
        }
        #
        # Add the csmdata directory to the filepath
        #
        $val = $usr_att{'csmdata'} . "/$val";
        #
        # If want to test for file-existance -- then...
        # Make sure the file exists -- if not set fit to -1 and next ELEMENT
        #
        my $test = 1;
        if ( exists($usr_att{'notest'}) ) {
           if ( $usr_att{'notest'} ) { $test = 0; }
        }
        if ( $test ) {
           my $filepath = $val;
           if ( $filepath =~ /^\$([^\/]+)(.+)$/ ) { $filepath = "$ENV{$1}$2"; }
           if ( ! -f "$filepath" ) {
              $fit[$i] = -1;
              next ELEMENT;
           }
        }

	# extract the element attributes
	my %attributes = $e->get_attributes();

	# Check each attribute; first against the configuration, then against values supplied
	# via the optional argument.
	# If an attribute doesn't match, then eliminate the element from further consideration
	foreach my $att_name (keys %attributes) {

	    # Is the attribute part of the configuration?
	    if ($cfg->is_valid_name($att_name)) {

		# Get the value for the attribute from the configuration
		my $cfg_val = $cfg->get($att_name);

		# Check for a match.  If the attributes don't match then skip this element
		# and move to the next after recording the no-match status
		if ($attributes{$att_name} eq $cfg_val) {
		    $matches++;
		}
		else {
		    $fit[$i] = -1;
		    next ELEMENT;
		}

	    } else {

		# If the attribute isn't part of the configuration then do addition
		# checks here.
		#

		# Has this attribute been supplied by the user?
		if (defined $usr_att{$att_name}) {

                   # Check for a match.  If the attributes don't match then skip this element
                   # and move to the next after recording the no-match status
                   if ($attributes{$att_name} eq $usr_att{$att_name}) {
                      $matches++;
                   } else {
                      $fit[$i] = -1;
                      next ELEMENT;
                   }
                }

	    } # Finished attribute checks

	} # Finished loop over attributes

	# At this point the attribute checking has been successful.  Record the matches.
	$values[$i] = $val;
	$fit[$i]    = $matches;

    } # Finished loop over elements in defaults file.

    # All elements have been examined.  Return the value from the best fit.  That's the
    # index of the max value of @fit.  In case of a tie it's the first one found.
    my $max_val = $fit[0];
    my $max_idx = 0;
    for (my $i = 1; $i <= $#elements; $i++) {
	if ($fit[$i] > $max_val) {
	    $max_val = $fit[$i];
	    $max_idx = $i;
	}
    }

    # If "best fit" is $max_val = -1, then no match was found.
    if ($max_val >= 0) {
	return $values[$max_idx];
    }
    else {
	return undef;
    }

}

#-----------------------------------------------------------------------------------------------

sub get_variable_names
{
# Return a list of the varible names found in the defaults file.

    my $self = shift;

    my $names_ref = $self->{'names'};

    return( sort( keys( %$names_ref ) ) );
}

#-----------------------------------------------------------------------------------------------
# Private methods
#-----------------------------------------------------------------------------------------------

sub _initialize
{
# Extract the element objects from the XML namelist defaults file and store internally.  Each element
# from the file contains a default value for a namelist variable.  There may be multiple elements for
# a single variable since the default value may depend on certain attributes.  The xml elements contain
# the default values along with all the attributes.
#
# This method adds (or merges) the following attributes to a NamelistDefaults object:
# 'names'    -- Reference to a hash whose keys are the namelist variable names.  The names have
#               been converted to lower case to implement a case insensitive lookup.
# 'elements' -- Reference to a hash that contains all the element objects extracted from the xml
#               file(s).  The keys are the namelist variable names (lowercase) and the values are
#               references to arrays that contain all the element objects for that variable.

    my ($self, $defaults_file) = @_;

    # Process the definition file
    my $xml = XML::Lite->new( $defaults_file );
    my $root = $xml->root_element();

    # Check for valid root node
    my $name = $root->get_name();
    $name eq "namelist_defaults" or die
	"ERROR: $defaults_file is not a namelist defaults file\n";

    # The children of the XML root object are xml elements, each of which contains a default
    # value for a namelist variable along with attributes that determine the configuration for
    # which the default is appropriate.
    my @children = $root->get_children();

    # If this object has already been initialized then set some convenience variables.
    my %names;
    my %elements;
    if ( exists( $self->{'names'}) ) {
       my $names_ref    = $self->{'names'};
       my $elements_ref = $self->{'elements'};
       %names    = %$names_ref;
       %elements = %$elements_ref;
    }

    # Store the elements in the definition file in a hash using the namelist variable
    # name (converted to lowercase) as the key.

    foreach my $e (@children) {
       my $name = lc $e->get_name();
       $names{$name} = 1;

       # If there have been previous elements for this variable, then append it to
       # the existing array.  Otherwise create a new array.
       if (defined $elements{$name}) {
	   push @{$elements{$name}}, $e;
       }
       else {
	   $elements{$name} = [$e];
       }
    }

    $self->{'names'}    = \%names;
    $self->{'elements'} = \%elements;

}

#-----------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------

1; # to make use or require happy
