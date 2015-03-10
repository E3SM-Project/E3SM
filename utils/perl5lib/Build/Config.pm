package Build::Config;
#-----------------------------------------------------------------------------------------------
#
# SYNOPSIS
# 
#   use Build::Config;
# 
#   # read the configuration definition and defaults xml files.
#   my $cfg = Build::Config->new("config_definition.xml", "config_setup.xml");
# 
#   # set configuration parameters
#   $cfg->set('dyn', 'fv');
#   $cfg->set('hgrid', '2x2.5');
# 
#   # query configuration parameters
#   my $exe = $cfg->get('cam_exe');
# 
#   # write a configuration cache file
#   $cfg->write_file("config_cache.xml", "configure commandline");
# 
# DESCRIPTION
# 
# Build::Config objects are used to represent features of a model
# configuration that must be specified at build time.
#
# new() Reads xml files that contain the configuration definition and
#       default values of the configuration parameters.
#
#       The "config_definition.xml" file contains all the allowable
#       parameters with a description of each one.  Where appropriate a
#       list of valid values of a parameter is given.  Generic default
#       values that are useful across many model configurations are
#       provided for some parameters.
#
#       The generic default values that are provided in the definition file
#       may be overridden by values in a setup file ("config_setup.xml")
#       that is assumed to provide appropriate values for a specific model
#       configuration.  This setup file is optional.
#
#  ========================
# Query methods:
#  ========================
#
# is_valid_name() Returns true if the specified parameter name is contained in
#       the configuration definition file.
#
# is_valid_value() Returns true if the specified parameter name is contained in
#       the configuration definition file, and either 1) the specified value is
#       listed as a valid_value in the definition file, or 2) the definition file
#       doesn't specify the valid values.
#
# get() Return the value of the specified configuration parameter.  Triggers
#       an exception if the parameter name is not valid.
#       ***NOTE*** If you don't want to trap exceptions, then use the query
#                  functions before calling this routine.
#
# get_names() Return the list of valid names in the configuration definition file.
#
# get_valid_values( ) Returns list of valid values for the specified parameter name.
#
#  ========================
# Put and I/O methods:
#  ========================
#
# set() Sets values of the configuration parameters.  It takes the
#       parameter name and its value as arguments.  An invalid parameter
#       name (i.e., a name not present in the definition file) triggers an
#       exception.  If the definition file contains valid values of a
#       parameter, then the set method checks for a valid input value.  If
#       and invalid value is found then an exception is thrown.
#       ***NOTE*** If you don't want to trap exceptions, then use the query
#                  functions before calling this routine.
#
# write_file() Write a configuration xml file.  The first argument is the
#       filename.  The second argument, if present, is the commandline of the
#       configure command that was invoked to produce the output configuration
#       file.  It is written to the output file to help document the procedure
#       used to configure the executable.
#
# print() Print the configuration to STDOUT.
# 
# COLLABORATORS
# 
# IO::File
# XML::Lite
#-----------------------------------------------------------------------------------------------
#
# Date        Author                  Modification
#-----------------------------------------------------------------------------------------------
#  2007-Aug   Erik Kluzek             Add get_names method
#  2006-Nov   Brian Eaton             Original version
#-----------------------------------------------------------------------------------------------

use strict;
#use warnings;
#use diagnostics;

use IO::File;
use XML::Lite;

sub new
{
    my $class = shift;
    my ($definition_file, $default_file) = @_;

    # bless the object here so the initialization has access to object methods
    my $cfg = {};
    bless( $cfg, $class );

    # Initialize the object with the configuration definition and its initial setup.
    if ( defined($definition_file) ) {
       $cfg->_initialize($definition_file, $default_file);
    }
    else {
	die "ERROR: $class new method requires a definition file\n";
    }

    return $cfg;
}

#-----------------------------------------------------------------------------------------------

sub is_valid_name
{
# Return true if the requested name is contained in the configuration definition.

    my ($self, $name) = @_;

    return defined($self->{$name}) ? 1 : 0;
}

#-----------------------------------------------------------------------------------------------

sub get
{
# Return requested value.

    my ($self, $name) = @_;

    defined($self->{$name}) or die "ERROR: unknown parameter name: $name\n";

    return $self->{$name}->{'value'};
}

#-----------------------------------------------------------------------------------------------

sub get_names
{
# Return list of valid names.

    my ($self) = @_;


    my @names = sort( keys( %$self ) );

    return @names;
}


#-----------------------------------------------------------------------------------------------

sub is_valid_value
{
# Return true if the specified parameter name is contained in
# the configuration definition file, and either 1) the specified value is
# listed as a valid_value in the definition file, or 2) the definition file
# doesn't specify the valid values.

    my ($self, $id, $value) = @_;

    # Check that the parameter name is in the configuration definition
    unless ($self->is_valid_name($id)) { return 0; }

    # Check that a list value is not supplied when parameter takes a scalar value.
    my $is_list_value = $self->{$id}->{'list'};
    unless ($is_list_value) {  # this conditional is satisfied when the list attribute is false, i.e., for scalars
	if ($value =~ /.*,.*/) { return 0; }   # the pattern matches when $value contains a comma, i.e., is a list
    }

    # Check that the value is valid
    my @valid_values = $self->get_valid_values( $id );
    if ( @valid_values ) {  # if no valid values are specified, then $value is automatically valid
	if ($is_list_value) {
	    unless (_list_value_ok($value, @valid_values)) { return 0; }
	}
	else {
	    unless (_value_ok($value, @valid_values)) { return 0; }
	}

    }

    return 1;
}

#-----------------------------------------------------------------------------------------------

sub get_valid_values
{
# Return list of valid values for a given id

    my ($self, $id) = @_;

    my @return_array;
    # Check that the parameter name is in the configuration definition
    if ($self->is_valid_name($id)) { 

       my $valid_values = $self->{$id}->{'valid_values'};
       if ( $valid_values ) {
          @return_array = split /,/, $valid_values;             #/
       }
    }

    return( @return_array );

}

#-----------------------------------------------------------------------------------------------

sub set
{
# Set requested value.
#
# This routine handles errors by throwing exceptions.  It will report exactly what problem was
# found in either the parameter name or requested value.
#
# To avoid dealing with exceptions use the is_valid_name(), is_valid_value() methods to get a
# true/false return before calling the set method.

    my ($self, $id, $value) = @_;

    # Check that the parameter name is in the configuration definition
    $self->is_valid_name($id) or die
	"ERROR: parameter name $id is not in the configuration definition\n";

    # Check that the value is valid
    my $valid_values = $self->{$id}->{'valid_values'};
    $self->is_valid_value($id, $value) or die
	    "ERROR: $value is not a valid value for parameter $id: valid values are $valid_values\n";

    # Add the new value to the object's internal data structure.
    $self->{$id}->{'value'} = $value;

    return 1;
}

#-----------------------------------------------------------------------------------------------

sub write_file
{
# Write a configuration definition file.

    my ($self, $filename, $commandline) = @_;

    my $fh = IO::File->new($filename, '>') or die "** can't open file: $filename\n";

    # head of the xml file
    print $fh <<"EOD";
<?xml version="1.0"?>

<config_definition>

EOD

    # add commandline if present
    if (defined $commandline) {
	print $fh <<"EOD";
<commandline>
$commandline
</commandline>
EOD
    }

    # add the entry elements
    my @ids = keys %$self;
    foreach my $id (sort @ids) {
	print $fh <<"EOD";
<entry id="$id" value="$self->{$id}->{'value'}" list="$self->{$id}->{'list'}" valid_values="$self->{$id}->{'valid_values'}">
$self->{$id}->{'definition'}
</entry>
EOD
    }

    # tail of the xml file
    print $fh <<"EOD";

</config_definition>
EOD


}

#-----------------------------------------------------------------------------------------------

sub print
{
# Print the configuration to STDOUT.

    my ($self) = @_;

    my @ids = keys %$self;
    foreach my $id (sort @ids) {
	printf "%12s = %s\n", $id, $self->{$id}->{'value'};
    }
}

#-----------------------------------------------------------------------------------------------
# Private methods
#-----------------------------------------------------------------------------------------------

sub _initialize
{
# Read the configuration definition file.  Create an anonymous hash with the following
# structure:
# { id1 => {value => "xxx", list => "XXX", valid_values => "yyy", definition => "zzz"},
#   id2 => {value => "xxx", list => "XXX", valid_values => "yyy", definition => "zzz"},
#   ...
#   idn => {value => "xxx", list => "XXX", valid_values => "yyy", definition => "zzz"},
# }


    my ($self, $definition_file, $setup_file) = @_;

    # Process the definition file
    my $xml = XML::Lite->new( $definition_file );
    my $root = $xml->root_element();

    # Check for valid root node
    my $name = $root->get_name();
    $name eq "config_definition" or die
	"ERROR: $definition_file is not a configuration definition file\n";

    # Each parameter is contained in an "entry" element.  Get all these elements.
    my @elements = $xml->elements_by_name('entry');

    # Loop over the elements...
    foreach my $e (@elements) {

        # and extract the attributes and element content.
	my %attributes = $e->get_attributes();
	my $content    = $e->get_content();
	# if present strip initial and final newlines from content
	$content =~ s/^\n{1}//;
	$content =~ s/\n{1}$//;

	# Look for the specific attributes that are contained in the configuration definition.
	my $id    = $attributes{'id'};
	my $value = $attributes{'value'};
        my $list  = $attributes{'list'};
	my $valid_values = defined $attributes{'valid_values'} ? $attributes{'valid_values'} : "";

	# Now add the attributes and content to the object's internal data structure.
	$self->{$id} = {'value' => $value, 'list' => $list, 'valid_values' => $valid_values, 'definition' => $content};
    }

    # Process the setup file
    if (defined $setup_file) {
	my $xml = XML::Lite->new( $setup_file );
	my $root = $xml->root_element();

	# Check for valid root node
	my $name = $root->get_name();
	$name eq "config_definition" or die
	    "ERROR: $definition_file is not a configuration definition file\n";

	# Each parameter is contained in an "entry" element.  Get all these elements.
	@elements = ();
	@elements = $xml->elements_by_name('entry');

	# Loop over the elements...
	foreach my $e (@elements) {

	    # and extract the attributes
	    my %attributes = $e->get_attributes();

	    # just get the parameter name and value
	    my $id    = $attributes{'id'};
	    my $value = $attributes{'value'};

	    # set new value
	    $self->set($id, $value);
	}
    } # end processing setup file
}

#-----------------------------------------------------------------------------------------------

sub _list_value_ok
{
# Check that all input values ($values_in may be a comma separated list)
# are contained in the list of valid values (@valid_values).
# Return 1 (true) if all input values are valid, Otherwise return 0 (false).

    my ($values_in, @valid_values) = @_;

    my @values = split /,/, $values_in;            #/

    my $num_vals = scalar(@values);
    my $values_ok = 0;

    foreach my $value (@values) {

	if (_value_ok($value, @valid_values)) { ++$values_ok; }

    }

    ($num_vals == $values_ok) ? return 1 : return 0;
}

#-----------------------------------------------------------------------------------------------

sub _value_ok
{

# Check that the input value is contained in list of
# valid values (@valid_values).  Return 1 (true) if input value is valid,
# Otherwise return 0 (false).

    my ($value, @valid_values) = @_;

    # If the valid value list is null, all values are valid.
    unless (@valid_values) { return 1; }

    $value =~ s/^\s+//;
    $value =~ s/\s+$//;
    foreach my $expect (@valid_values) {
	if ($value =~ /^$expect$/i) { return 1; }
    }

    return 0;
}

#-----------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------

1; # to make use or require happy
