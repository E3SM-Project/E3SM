package Build::NamelistDefinition;
my $pkg = 'Build::NamelistDefinition';
#-----------------------------------------------------------------------------------------------
#
# SYNOPSIS
# 
#   use Build::Namelist;
#   use Build::NamelistDefinition;
# 
#   # Create a namelist definition object (read the namelist definition file).
#   my $nldef = Build::NamelistDefinition->new("namelist_definition.xml");
# 
#   # Create a namelist object from an input file that contains one or more namelist groups.
#   my $nl = Build::Namelist->new('namelist.in');
#
#   # Validate the namelist object.
#   my $nl_valid = $nldef->validate($nl);
#
#   # Query the definition object for the filepath of the definition file
#   my $definition_file = $nldef->get_file_name();
#
#   # Query the definition object to find which group a variable belongs to.
#   my $group = $nldef->get_group_name('variable_name');
#
#   # Query the definition object for the string length of a requested variable.
#   # get_str_len returns 0 if the requested variable isn't a character type
#   my $strlen = $nldef->get_str_len('variable_name');
#
#   # Query the definition object whether the requested variable is the pathname of
#   # an input dataset.  If it is return the input_pathname attribute.  If not
#   # return "".
#   my $pathname_type = $nldef->is_input_pathname('variable_name');
#
#   # Write validated namelist to output file.
#   $nl_valid->write('namelist.out');
# 
#   # If some of your namelist definition are in different files -- use Add
#   $nldef->Add("namelist_definition2.xml");
# 
# DESCRIPTION
# 
# Build::NamelistDefinition objects encapsulate a namelist definition.
# They provide a method used to validate a namelist object (created by the
# Build::Namelist module) against a namelist definition file.  The
# validation currently consists of making sure that each variable in the
# namelist object is defined in the definition file.
#
# All variables are checked and problems reported to stdout.  If any
# problems were encountered the object will throw an exception.
#
# If all variables are successfully validated, then the returned namelist
# object will have all variables in the namelist group(s) that are
# specified in the definition file.  This means that the namelist groups in
# the input namelist object are ignored.
#
# The module also provides accessor methods to extract data from the
# namelist definition file.  This is useful for scripts that produce
# namelist documentation from the definition file.
#
# METHODS
#
# new()
#       Reads xml file that contains the namelist definition.
#       The "namelist_definition.xml" file contains all the allowable
#       variables with a description of each one.  Where appropriate a
#       list of valid values of a variable is given.
#
# add()
#       Adds definitions from an additional file.
# 
# validate()
#       Validate a namelist object (created by the Build::Namelist module)
#       against a namelist definition file.  Each variable is checked to
#       verify that it is defined in the definition file.  Also the value
#       of each variable is checked to verify that it's a string or a numeric
#       type.
#
#       The returned namelist object will have all variables contained in
#       the correct namelist group as specified by the definition file.
#
# get_file_name()
#       Return the filepath of the file that contains the namelist definition.
#
# get_group_name()
#       If the requested variable name is contained in the namelist
#       definition then return its namelist group.  Otherwise return an
#       empty string.
#
# get_str_len()
#       Return 'str_len' attribute for requested variable.  This is the
#       length from the type declaration in the definition.  If the
#       variable is not a character type then the return value is 0.
#
# get_var_names()
#       Return an alphabetized list of all variable names that are
#       contained in the definition file.  If the catagory is specified
#       (via the optional argument pair 'catagory'=>'cat_name') then only
#       the variables whose 'catagory' attribute have the value 'cat_name'
#       are returned.  This option is designed for producing documentation.
#
# get_var_type()
#       Return the type declaration for the requested variable.
#
# get_var_doc()
#       Return the documentation for the requested variable.
#
# get_var_doc_html()
#       Return the documentation for the requested variable.  Insert html tags
#       for presentation in a table.
#
# get_valid_values()
#       Get list of valid values.
#
# is_input_pathname()
#       Return 'input_pathname' attribute for requested variable.  The
#       value is 'abs' if the variable contains the absolute pathname for
#       an input dataset.  A value of 'rel:var_name' means the variable
#       contains a relative pathname and var_name is the name of another
#       namelist variable that contains the root directory for the relative
#       pathname.  A value of '' is returned if the requested variable is
#       not the pathname of an input dataset.
#
# is_valid_value()
#       Check if a single (non-array, non-split) input value for a variable
#       is valid.
#
# COLLABORATORS
# 
# IO::File
# XML::Lite
# Build::Namelist
#-----------------------------------------------------------------------------------------------
#
# Date        Author                  Modification
#-----------------------------------------------------------------------------------------------
# 2007-Sep    Brian Eaton             Original version
# 2008-May    Erik Kluzek             Add methods to validate data types
# 2008-Sep    Erik Kluzek             Add add method
#-----------------------------------------------------------------------------------------------

use strict;
#use warnings;
#use diagnostics;

use IO::File;
use XML::Lite;
use Build::Namelist;

#-----------------------------------------------------------------------------------------------
# Public methods
#-----------------------------------------------------------------------------------------------

sub new
{
    my $class = shift;
    my ($definition_filepath) = @_;

    # bless the object here so the initialization has access to object methods
    my $nl_definition = {};
    bless( $nl_definition, $class );

    # Add the filepath of the definition file to the object attributes.
    $nl_definition->{'definition_filepath'} = $definition_filepath;

    # Initialize the object with the namelist definition.
    $nl_definition->_initialize($definition_filepath);

    return $nl_definition;
}

#-----------------------------------------------------------------------------------------------

sub add
{
    my $self  = shift;
    my ($definition_filepath) = @_;

    # Append the filepath of the definition file to the object attributes.
    $self->{'definition_filepath'} .= ", " . $definition_filepath;

    # Add additional namelist definitions.
    $self->_initialize($definition_filepath);
}

#-----------------------------------------------------------------------------------------------

sub validate
{

# Validate a namelist object (created by the Build::Namelist module)
# against a namelist definition file.  All variables are checked and
# problems reported to stdout.  If any problems were encountered the object
# will throw an exception.
#
# If all variables are successfully validated, then the returned namelist
# object will have all variables in the namelist group(s) that are
# specified in the definition file.  This means that the namelist groups
# in the input namelist object are ignored.

    my $self  = shift;
    my $nl    = shift;    # namelist to be validated
    
    # Create an empty namelist which will be populated with variables from the
    # input namelist as they are validated.
    my $nl_valid = Build::Namelist->new();

    # Loop over the groups in the input namelist
    my @groups = $nl->get_group_names();
    for my $group (@groups) {

	# Loop over the variables in the namelist group
	my @vars = $nl->get_variable_names($group);
	for my $var (@vars) {

	    # Get the variable's value
	    my $value = $nl->get_variable_value($group, $var);

	    # Validate the variable/value pair.  This method throws an exception
	    # when an error is encountered.  If the validation is successful, then
	    # the valid group name for the variable is returned.
	    my $valid_group = $self->_validate_pair($var, $value);

	    # Add the validated variable to the output namelist
	    $nl_valid->set_variable_value($valid_group, $var, $value);
	}
    }

    return $nl_valid;
}

#-----------------------------------------------------------------------------------------------

sub get_file_name
{
# Return the name of the file that contains the namelist definition.

    my $self = shift;

    return $self->{'definition_filepath'};
}

#-----------------------------------------------------------------------------------------------

sub get_group_name
{
# If the requested name is contained in the namelist definition then return its
# namelist group.  Otherwise return an empty string.

    my ($self, $name) = @_;
    my $lc_name = lc $name;

    return defined($self->{$lc_name}) ? $self->{$lc_name}->{'group'} : "";
}

#-----------------------------------------------------------------------------------------------

sub get_str_len
{
# Return 'str_len' attribute for requested variable

    my ($self, $name) = @_;
    my $lc_name = lc $name;

    return $self->{$lc_name}->{'str_len'};
}


#-----------------------------------------------------------------------------------------------

sub is_input_pathname
{
# Return 'input_pathname' attribute for requested variable

    my ($self, $name) = @_;
    my $lc_name = lc $name;

    return $self->{$lc_name}->{'input_pathname'};
}

#-----------------------------------------------------------------------------------------------

sub get_var_names
{
# Return alphabetized list of all variable names that are contained in the definition file.
# If the optional argument pair 'catagory'=>'cat_name' are supplied then only variables
# whose 'catagory' attribute has the value 'cat_name' will be returned.

    my $self = shift;
    my %opt  = @_;      # options

    # Put all keys from the definition object except for 'definition_filepath' into a new
    # hash.  Then return the sorted keys.
    my %var = ();
    foreach my $k (keys %$self) {
	unless ($k eq 'definition_filepath') {

	    # If a specific catagory has been requested then only add variables in 
	    # that catagory.
	    if ( defined $opt{'catagory'} ) {
		if ( $opt{'catagory'} eq $self->{$k}->{'catagory'} ) { $var{$k} = ''; }
	    }
	    else {
		$var{$k} = '';
	    }
	}
    }
    return sort keys %var;
}

#-----------------------------------------------------------------------------------------------

sub get_var_type
{
# Return 'type' attribute for requested variable

    my ($self, $name) = @_;
    my $lc_name = lc $name;

    return $self->{$lc_name}->{'type'};
}

#-----------------------------------------------------------------------------------------------

sub get_var_doc
{
# Return documentation for requested variable

    my ($self, $name) = @_;
    my $lc_name = lc $name;

    return $self->{$lc_name}->{'doc'};
}

#-----------------------------------------------------------------------------------------------

sub get_var_doc_html
{
# Return documentation for requested variable with html tags included.

    my ($self, $name) = @_;
    my $lc_name = lc $name;

    my $doc = $self->{$lc_name}->{'doc'};

    # Insert a line break in front of the 'Default:' token.
    $doc =~ s/(defaults?:)/<br\/>$1/i;

    return $doc;
}

#-----------------------------------------------------------------------------------------------

sub get_valid_values
{
# Return list of valid_values as an array for requested variable
# To return without quotes use the 'noquotes'=>1 option.
    my ($self, $name, %opts) = @_;
    my $lc_name = lc $name;

    my $valid_values = $self->{$lc_name}->{'valid_values'};
    my @values = split( /,/, $valid_values );
    my $str_len = $self->get_str_len( $lc_name );
    # if string type and NOT noquote option and have a list -- add quotes around values
    if ( ! defined($opts{'noquotes'}) || ! $opts{'noquotes'} ) {
       if ( $#values > -1 && ($str_len > 0) ) {
          for( my $i=0; $i <= $#values; $i++ ) {
             $values[$i] = "'$values[$i]'";
          }
       }
    }
    return( @values );
}

#-----------------------------------------------------------------------------------------------

sub is_valid_value {

# Check that a given single variable value matches the input list 
# NOTE: This only works for a single value entered (NOT a list of values)

  my $self           = shift;
  my $variable_name  = lc shift;
  my $variable_value = shift;
  my %opts           = @_;

  my @valid_values   = $self->get_valid_values( $variable_name, %opts );
  if ( $#valid_values > -1 ) {
     foreach my $val ( @valid_values ) {
        if ( $variable_value eq $val ) {
           return 1;
        }
     }
     return 0;
  } else {
     return 1;
  }
}

#-----------------------------------------------------------------------------------------------
# Private methods
#-----------------------------------------------------------------------------------------------

sub _initialize
{
# Read the namelist definition file.  Add a hash entry to the object for each
# namelist entry in the definition file.  The hash entries look like this:
#
#   id =>                               # id is the name of the namelist variable
#         {type           => "...",     # variable's type
#          str_len        => "...",     # length of a string type, 0 if not a string
#          arr_len        => "...",     # size of an array type, 0 if not an array
#          input_pathname => "...",     # if value is an input pathname then this attribute is
#                                       # set to either 'abs' for an absolute pathname, or 'rel:var_name'
#                                       # for a relative pathname.  If the pathname is relative, 
#                                       # var_name specifies the namelist variable that contains the
#                                       # root directory that the pathname is relative to.  If variable
#                                       # is not an input pathname this entry is set to ''.
#          catagory       => "...",     # catagory used in documentation
#          group          => "...",     # namelist group
#          valid_values   => "...",     # valid values (if easy to list)
#          doc            => "..."}     # documentation for variable
#

    my ($self, $definition_file) = @_;

    # Process the definition file
    my $xml = XML::Lite->new( $definition_file );
    my $root = $xml->root_element();

    # Check for valid root node
    my $name = $root->get_name();
    $name eq "namelist_definition" or die
	"ERROR: $definition_file is not a namelist definition file\n";

    # Each namelist variable is contained in an "entry" element.  Get all these elements.
    my @elements = $xml->elements_by_name('entry');

    # Loop over the elements...
    foreach my $e (@elements) {

        # and extract the attributes and element content.
	my %attributes = $e->get_attributes();
	my $content    = $e->get_content();

	# Look for the specific attributes that are contained in the namelist definition.
	my $id             = lc $attributes{'id'};    # make interfaces case insensitive
	my $type           = $attributes{'type'};
        my $input_pathname = defined $attributes{'input_pathname'} ? $attributes{'input_pathname'} : "";
        my $catagory       = $attributes{'catagory'};
        my $group          = $attributes{'group'};
	my $valid_values   = defined $attributes{'valid_values'} ? $attributes{'valid_values'} : "";

	# Parse the type specification for the following info:

        # Is the type string or numeric?  A string type will be indicated by $str_len > 0.
	# $str_len = 0 indicates a numeric type.
	my $str_len = 0;
	if ( $type =~ m/char\*(\d+)/ ) {
	    $str_len = $1;
	}

	# Is the type an array or a scalar?  An array will be indicated by $arr_len > 0 
	# where $arr_len is the size of the array.  $arr_len = 0 indicates a scalar.
	my $arr_len = 0;
	if ( $type =~ m{\(         # opening paren
                        (.*)       # capture everything between the parens in $1
                        \)         # closing paren
                       }x ) {
	    # split the dimensions between the parenthesis on "," and multiply the 
	    # dimensions together to get the array size
	    my @dims = split /,/, $1;  #/
            $arr_len = 1;
	    foreach my $dim (@dims) {
		$arr_len *= $dim;
	    }
	}

	# Now add the attributes and content to the object's internal data structure.
	$self->{$id} = {'type'           => $type,
			'str_len'        => $str_len,
			'arr_len'        => $arr_len,
			'input_pathname' => $input_pathname,
			'catagory'       => $catagory,
			'group'          => $group,
			'valid_values'   => $valid_values,
			'doc'            => $content,
	};
    }

}

#-----------------------------------------------------------------------------------------------

sub _validate_pair
{

# The validation consists of the following checks:
#
# . The variable must be defined in the namelist definition file.
# . The variable's value is verified as string or numeric.
#

    my $self   = shift;
    my $var_in = lc shift;   # namelist variable
    my $value  = shift;      # namelist variable's value

    my $def    = $self->{'definition_filepath'};

    # Is variable in namelist definition file?
    # If the variable is an array, and an array element has been specified in the namelist,
    # then we need to strip off the index specification before checking whether or not the
    # name is valid.
    my $var;
    if ($var_in =~ m/(\w[\w%]*)    # capture the array name (allow it to be an element in a derived type
                     \(            # opening paren for index spec
                     [\d,]+        # 1 or more digits and commas
                     \)            # closing paren of index spec
                    /ix) {
	$var = $1;
    }
    else {
	$var = $var_in;
    }

    $self->_is_valid_name($var) or die 
	"ERROR: in _validate_pair (package $pkg): Variable name $var not found in $def \n";

    # Parse the value string

    # Is the input value a string or numeric type?
    # A string type must start with a quote with optional leading whitespace.
    my $value_is_string = 0;
    if ( $value =~ m{^\s*['"]} ) { $value_is_string = 1; }

    # Get string length from definition file (returns 0 for a numeric type)
    my $str_len_def = $self->get_str_len($var);
    my $type_def    = $self->_get_type($var);

    # Check for mismatch
    if ( $value_is_string and ($str_len_def == 0) ) {
	die "ERROR: in _validate_pair (package $pkg): Variable name $var has an input value of type string, \n".
	    "$value \n".
	    "but is defined as type $type_def in $def \n";
    }
    elsif ( ! $value_is_string and ($str_len_def > 0) ) {
	die "ERROR: in _validate_pair (package $pkg): Variable name $var has an input value of type numeric, \n".
	    "$value \n".
	    "but is defined as type $type_def in $def \n";
    }	

    # 22 September 2007, bee
    # The intention is to include more rigorous checking for valid input values.
    # But postpone this for now.  It requires re-parsing the value that's currently 
    # only available as a string from the namelist parser.  The functionality of 
    # breaking input values into arrays of the input type belongs in the namelist 
    # parser and shouldn't be done here.
    # The following validations depend on breaking the input value string into
    # an array of elements.
    # . The variable's value is checked to verify that list values aren't specified for scalar
    #   variables.
    # . For string input check that the length of input strings doesn't exceed the declared
    #   string length
    # . The variable's value is checked against any valid values specified by the definition.

    # 6 Oct 2008, bee
    # Checking has been added for valid values that puts the validation in
    # the namelist object rather than here.  The reason is to make use of
    # the regexps that are defined in the namelist object to do the
    # parsing.  Rather than using those regexps to re-parse the values, the
    # namelist parser should be refactored to store that information during
    # the initial parsing.  Then the namelist object can then be queried
    # for the type of the input values, and the validation done here, i.e.,
    # in the module that knows the definition of the variable's type.


    # Get the type description hash for the variable
    # and validate the value with the datatype in the namelist object
    # This method throws an exception when an error is encountered.
    my %type_ref = $self->_get_datatype($var);
    Build::Namelist::validate_variable_value($var, $value, \%type_ref);

    # Checks all passed.  Return valid group name.
    return $self->get_group_name($var);
}

#-----------------------------------------------------------------------------------------------

sub _is_valid_name
{
# Return true if the requested name is contained in the namelist definition.

    my ($self, $name) = @_;
    my $lc_name = lc $name;

    return defined($self->{$lc_name}) ? 1 : 0;
}

#-----------------------------------------------------------------------------------------------

sub _get_type
{
# Return 'type' attribute for requested variable

    my ($self, $name) = @_;
    my $lc_name = lc $name;

    return $self->{$lc_name}->{'type'};
}

#-----------------------------------------------------------------------------------------------

sub _get_datatype 
#
# Return hash of description of data type read in from the file:
# Hash keys are:
#      type           type description (char, logical, integer, or real)       (string)
#      strlen         Length of string (if type char)                          (integer)
#      arrayNDims     Number of dimensions (0=scalar,1=array,2=2D array, etc.) (integer)
#      arrayDims      Reference to array of size of each dimension             (integer)
#      validValues    Reference to array of valid values                       (string)
#
{
    my ($self, $name) = @_;
    my $lc_name = lc $name;
    my $nm = "_get_datatype";

    my %datatype;
    my $type_def    = $self->_get_type($lc_name);
    $datatype{'strlen'} = $self->get_str_len($lc_name);
    if (      $type_def =~ /^(char|logical|integer|complex|real)/ ) {
      $datatype{'type'} = $1;
    } else {
	die "ERROR: in $nm (package $pkg): datatype $type_def is NOT valid\n";
    }
    # Arrays
    if ( $type_def =~ /\(/ ) {
      if ( $type_def =~ /\(([0-9, ]+)\)$/ ) {
        my @dimSizes = split( /,[ ]*/, $1 );
        $datatype{'arrayNDims'} = $#dimSizes + 1;
        $datatype{'arrayDims'}  = \@dimSizes;
      } else {
	die "ERROR: in $nm (package $pkg): datatype $type_def is NOT valid\n";
      }
    # Scalars
    } else {
      my @dimSizes;
      push( @dimSizes, "1" );
      $datatype{'arrayNDims'} = 0;
      $datatype{'arrayDims'}  = \@dimSizes;
    }
    my @valid_values = $self->get_valid_values( $lc_name );
    $datatype{'validValues'}  = \@valid_values;
    return( %datatype );
}


1; # to make use or require happy
