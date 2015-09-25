package ConfigCase;
my $pkg_nm = __PACKAGE__;

#-----------------------------------------------------------------------------------------------
# SYNOPSIS
# 
#   use ConfigCase;
# 
#   # create a new empty config case
#   my $cfg = ConfigCase->new("");
#
#   # set some parameters
#   $cfg->set($id, $value);
# 
#   # get some parameters
#   my $value = $cfg->get($id);
#
#   # Write an xml file out
#   $cfg->write_file("$caseroot/env_run.xml", "xml","$caseroot");
#
# DESCRIPTION
# 
# ConfigCase objects are used to represent features of a CIME model
# configuration that must be specified when a new case is created.
#
# new() Reads xml files that contain the configuration definition and
#       default values of the configuration parameters.
#
#       The "definition.xml" file contains all the allowable
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
# is_valid_name() 
#       Returns true if the specified parameter name is contained in
#       the configuration definition file.
#
# get() Return the value of the specified configuration parameter.  Triggers
#       an exception if the parameter name is not valid.
#       ***NOTE*** If you don't want to trap exceptions, then use the query
#                  functions before calling this routine.
#
# set() 
#       Sets values of the configuration parameters.  It takes the
#       parameter name and its value as arguments.  An invalid parameter
#       name (i.e., a name not present in the definition file) triggers an
#       exception.  If the definition file contains valid values of a
#       parameter, then the set method checks for a valid input value.  If
#       and invalid value is found then an exception is thrown.
#       ***NOTE*** If you don't want to trap exceptions, then use the query
#                  functions before calling this routine.
#
# write_file() 
#       Write an xml file.  The first argument is the
#       filename.  The second argument, if present, is the commandline of the
#       setup command that was invoked to produce the output configuration
#       file.  It is written to the output file to help document the procedure
#       used to petsetup the executable.
#
#-----------------------------------------------------------------------------------------------

use strict;
use English;
#use warnings;
#use diagnostics;
use IO::File;
use XML::LibXML;
use File::Basename;
use Log::Log4perl qw(get_logger);
use SetupTools;


my $logger;

BEGIN{
    $logger = get_logger();
}

#-----------------------------------------------------------------------------------------------
sub new
{
    my $class = shift;
    my ($definition_file, $default_file) = @_;

    # bless the object here so the initialization has access to object methods
    my $cfg = {};
    bless( $cfg, $class );
    return $cfg;
}

#-----------------------------------------------------------------------------------------------
sub add_config_variables
{
    # define variables in $cfg_ref from xml file

    my ($self, $file, $srcroot, $cimeroot, $model) = @_;

    (-f $file) or $logger->logdie ("ERROR ConfigCase.pm::add_config_variables file \'$file\' does not exist \n");
    my $xml = XML::LibXML->new( no_blanks => 1)->parse_file($file);
    my @nodes = $xml->findnodes(".//entry");
    if (! @nodes) {
	$logger->logdie( "ERROR add_config_variables: no variable elements in file $file \n"); 
    }
    foreach my $node (@nodes) 
    {
	my $id = $node->getAttribute('id');
	foreach my $define_node ($node->childNodes()) 
	{
	    my $node_name  = $define_node->nodeName();
	    my $node_value = $define_node->textContent();
	    if (defined $node_value) {
		$node_value =~ s/\$MODEL/$model/;
		$node_value =~ s/\$CIMEROOT/$cimeroot/;
		if (-d $srcroot) {
		    $node_value =~ s/\$SRCROOT/$srcroot/;
		}

		# now set the initial value to the default value - this can get overwritten
		if ($node_name eq 'default_value') {
		    $self->{$id}->{'value'} = $node_value;
		} else {
		    $self->{$id}->{$node_name} = $node_value;
		}
		$logger->debug("id= $id name = $node_name value = $node_value\n");
	    }
	}
	if (! defined $self->{$id}->{'value'} ) {
	    $logger->logdie( "ERROR add_config_variables: default_value must be set for $id in $file\n");
	}
    }
}

#-----------------------------------------------------------------------------------------------
sub set
{
    # Set requested value.
    # This routine handles errors by throwing exceptions.  It will report exactly what problem was
    # found in either the parameter name or requested value.
    # To avoid dealing with exceptions use the is_valid_name(), is_valid_value() methods to get a
    # true/false return before calling the set method.

    my ($self, $id, $value) = @_;

    # Check that the parameter name is in the configuration definition
    unless ($self->is_valid_name($id)) { 
	$logger->logdie ("ERROR ConfigCase::set: $id is not a valid name \n");
    }

    # Get the type description hash for the variable and check that the type is valid
    # This method throws an exception when an error is encountered.
    my %type_ref = $self->_get_typedesc($id);
    SetupTools::validate_variable_value($id, $value, \%type_ref);

    # Check that the value is valid
    my $valid_values = $self->{$id}->{'valid_values'};
    if ( defined $valid_values && $valid_values ne "" ) {
	my $value = _clean($value);
	my $is_list_value = $self->{$id}->{'list'};
	SetupTools::is_valid_value($id, $value, $valid_values, $is_list_value) 
	    or $logger->logdie(
		"ERROR: value of $value is not a valid value for parameter $id: valid values are $valid_values\n");
    }
    # Add the new value to the object's internal data structure.
    $self->{$id}->{'value'} = $value;

    return 1;
}

#-----------------------------------------------------------------------------------------------
sub get
{
    # Return requested value.
    my ($self, $name) = @_;

    defined($self->{$name}) or $logger->logde( "ERROR ConfigCase.pm::get: unknown parameter name: $name\n");
    $logger->debug("GET: $name $self->{$name}->{value}\n");
    return $self->{$name}->{'value'};
}

#-----------------------------------------------------------------------------------------------
sub get_valid_values
{
    # Return list of valid_values as an array for requested variable
    # To return without quotes use the 'noquotes'=>1 option.
    my ($self, $name, %opts) = @_;

    my $valid_values = $self->{$name}->{'valid_values'};
    my $type = $self->{$name}->{'type'};
    my @values;
    if(defined $valid_values){
	@values = split( /,/, $valid_values );

	# if string type and NOT noquote option and have a list -- add quotes around values
	if ( ! defined($opts{'noquotes'}) || ! $opts{'noquotes'} ) {
	    if ( $#values > -1 && ($type eq "char") ) {
		for( my $i=0; $i <= $#values; $i++ ) {
		    $values[$i] = "'$values[$i]'";
		}
	    }
	}
    }
    return( @values );
}

#-----------------------------------------------------------------------------------------------
sub is_valid_name
{
    # Return true if the requested name is contained in the configuration definition.

    my ($self, $name) = @_;
    return defined($self->{$name}) ? 1 : 0;
}

#-----------------------------------------------------------------------------------------------
sub write_file
{
    # Write an env_xxx.xml file (specified by $filename) in $caseroot
    my ($self, $output_xml_file, $caseheaders, $caseroot, $cimeroot, $input_xml_file) = @_;

    if ( -f $output_xml_file ) { unlink( $output_xml_file ); }
    my $fh = IO::File->new($output_xml_file, '>' ) or die "can't open file: $output_xml_file\n";

    print $fh "<?xml version=\"1.0\"?> \n";
    print $fh "\n";
    print $fh "<config_definition> \n";
    print $fh "\n";

    if ($caseheaders) {
	_print_file_header($fh, $caseheaders, $output_xml_file); 
    }

    if ($output_xml_file =~ /env_archive.xml/) {

	if (! $input_xml_file ) {
	    $logger->logdie ("ERROR write_file: must specify input_xml_file as argument for writing out $output_xml_file \n");
	} else {
	    open CONFIG_ARCHIVE, $input_xml_file or die $!;
	    while (<CONFIG_ARCHIVE>) {
		chomp;
		print $fh "$_\n";
	    }
	    close (CONFIG_ARCHIVE);
	}

    } else {

	my @groups;
	my $file_xml = basename($output_xml_file);
	foreach my $id (keys %$self) {
	    my $file_id  = $self->{$id}->{'file'};
	    next unless defined $file_id;

	    if ($file_id eq $file_xml) {
		my $group = $self->{$id}->{'group'};
		push (@groups, $group);
	    }
	}
	@groups = sort (@groups);
	@groups = _unique(@groups);


	# Write all the groups out to the target xml file
	print $fh "\n\n";
	print $fh "<groups>\n";   	    
	foreach my $group (@groups) {
	    print $fh "   <group>$group</group> \n"; 
	}
	print $fh "</groups>\n";   	    
	print $fh "\n";

	# Write out all necessary groups to the output xml file
	foreach my $group (@groups) {
	    foreach my $id (sort keys %$self) {
		if ($self->{$id}->{'group'} eq $group ) {
		    my $value         = $self->{$id}->{'value'};
		    my $type          = $self->{$id}->{'type'};
		    my $valid_values  = $self->{$id}->{'valid_values'};
		    my $desc          = $self->{$id}->{'desc'};
		    my $is_list_value = $self->{$id}->{'list'};
		    my $file          = $self->{$id}->{'file'};
		    if (defined $file) {
			if ($file eq $file_xml) {
			    write_xml_entry($fh, $id, $value, $type, $valid_values, $desc, $group, $is_list_value);
			}
		    } else {
			$logger->logdie("file attribute for variable $id is not defined \n");
		    }
		}
	    }
	}
    }
    print $fh "\n";
    print $fh "</config_definition> \n";
}

#-----------------------------------------------------------------------------------------------
sub write_xml_entry
{
    # Output xml file entry
    my ($fh, $id, $value, $type, $valid_values, $desc, $group, $is_list_value) = @_;

    $value =~ s/'/&apos;/g;
    $value =~ s/\</&lt;/g;
    $value =~ s/\</&gt;/g;
    if(defined $desc){
	$desc =~ s/^\n//;
	$desc =~ s/\n$//;
	$desc =~ s/^ *//;
	$desc =~ s/ *$//g;
	chomp $desc;
    }else{
	$desc = "no description available";
    }
    print $fh "\n";
    print $fh "<entry id=\"$id\"  value=\"$value\">\n";   	    
    print $fh "  <type>$type</type> \n"; 
    if (defined $valid_values && $valid_values  ne '') {print $fh "  <valid_values>$valid_values</valid_values> \n";}
    if (defined $is_list_value && $is_list_value ne '') {print $fh "  <list>$is_list_value</list> \n";}
    print $fh "  <group>$group</group> \n"; 
    print $fh "  <desc>$desc</desc> \n";
    print $fh "</entry> \n";
}

#-----------------------------------------------------------------------------------------------
#                               Private routines
#-----------------------------------------------------------------------------------------------
sub _get_type
{
# Return 'type' attribute for requested variable

    my ($self, $name) = @_;
    
    return $self->{$name}->{'type'};
}

#-----------------------------------------------------------------------------------------------
sub _get_typedesc
#
# Return hash of description of data type read in from the file:
# Hash keys are:
#      type           type description (char, logical, integer, or real) (string)
#      strlen         Length of string (if type char)                          (integer)
#      validValues    Reference to array of valid values                 (string)
#
{
    my ($self, $name) = @_;
    my $nm = "_get_typedesc";

    my %datatype;
    my $type_def = $self->_get_type($name);
    my $lc_name = lc $name;
    if ($type_def =~ /^(char|logical|integer|real)/ ) {
	$datatype{'type'} = $1;
    } else {
	$logger->logdie ("ERROR: in $nm (package $pkg_nm): datatype $type_def is NOT valid for $name \n");
    }
    if ( $datatype{'type'} eq "char" ) {
       if ($type_def =~ /^char\*([0-9]+)/ ) {
           $datatype{'strlen'} = $1;
       } else {
           $datatype{'strlen'} = 9999999;
       }
    } else {
       $datatype{'strlen'} = undef;
    }
    my @valid_values = $self->get_valid_values( $name );
    $datatype{'validValues'}  = \@valid_values;
    return( %datatype );
}

#-------------------------------------------------------------------------------
sub _get_group_names
{
    my ($headerfile, $filename) = @_;

    my $file = basename($filename);

    my $xml_nodes = XML::LibXML->new( no_blanks => 1)->parse_file($headerfile);
    my @group_nodes = $xml_nodes->findnodes(".//file[\@name=\"$file\"]/groups/group");
    my @groups;
    foreach my $group_node (@group_nodes) {
	my $group = $group_node->textContent();
	push (@groups, $group);
    }
    return (@groups);
}

#-------------------------------------------------------------------------------
sub _print_file_header
{
    my ($fh, $headerfile, $filename) = @_;

    my $outputfile = basename($filename);
    my $xml = XML::LibXML->new( no_blanks => 1)->parse_file($headerfile);
    my @nodes = $xml->findnodes(".//file[\@name=\"$outputfile\"]/header");
    if (! @nodes) {
	$logger->logdie (" ERROR: no header nodes found for file $outputfile \n");
    }
    my $text = $nodes[0]->textContent();
    chomp($text);
    print $fh "<header>\n";
    print $fh "$text";
    print $fh "</header>\n";
}


#-------------------------------------------------------------------------------
sub _clean
{
    my ($name) = @_;
    $name =~ s/^\s+//; # strip any leading whitespace 
    $name =~ s/\s+$//; # strip any trailing whitespace
    return ($name);
}

#-----------------------------------------------------------------------------------------------
#                               Public routines to be deprecated
#-----------------------------------------------------------------------------------------------
sub write_docbook_master
{
    # Write the documentation on the configuration to an output README file.

    my $self = shift;
    my $filename = shift;   # filepath for output namelist

    my $fh;
    if ( -f $filename ) { unlink( $filename ); }
    $fh = IO::File->new($filename, '>' ) or die "can't open file: $filename\n";

    my $gid;
    $logger->info("Writing $filename\n");
    if ($filename =~ "case") { 
        $gid = "case";
	print $fh "<table><title>env_case.xml variables</title>\n";
    } elsif($filename =~ "build") {
        $gid = "build";
	print $fh "<table><title>env_build.xml variables</title>\n";
    } elsif($filename =~ "mach_pes") {
        $gid = "mach_pes";
	print $fh "<table><title>env_mach_pes.xml variables</title>\n";
    } elsif($filename =~ "run") {
        $gid = "run";
	print $fh "<table><title>env_run.xml variables</title>\n";
    }
	print $fh "<tgroup cols=\"4\">\n";
	print $fh "<thead>\n";
	print $fh "<row> \n";
	print $fh "<entry>Name</entry>\n";
	print $fh "<entry>Type</entry>\n";
	print $fh "<entry>Default</entry>\n";
	print $fh "<entry>Description [Valid Values]</entry>\n";
	print $fh "</row> \n";
	print $fh "</thead>\n";
	print $fh "<tbody>\n";
    
    my @ids = keys %$self;
    foreach my $id (sort @ids) {

        my $desc = "";
        my $valid = "";
        my $value = "";
        my $type  = "";

        $desc = $self->{$id}->{'desc'};
  	$valid = $self->{$id}->{'valid_values'};
	$value = $self->{$id}->{'value'};
	$type  = $self->{$id}->{'type'};

        if ( ! defined($desc) ) { $desc = ""; }
        if ( ! defined($valid) ) { $valid = ""; }
        if ( ! defined($value) ) { $value = ""; }
        if ( ! defined($type) ) { $type = ""; }

	    if ( $self->{$id}->{'group'} =~ $gid ) {
		print $fh "<row> \n";
		print $fh "<entry>$id</entry>\n";
		print $fh "<entry>$type</entry>\n";
		print $fh "<entry>$value</entry>\n";
		if ($self->{$id}->{'valid_values'}) {
		    print $fh "<entry>$desc [$valid] </entry>\n";
		} else {
		    print $fh "<entry>$desc</entry>\n";
                }
		print $fh "</row>\n";
	    }	    
    }
    print $fh "</tbody>\n";
    print $fh "</tgroup>\n";
    print $fh "</table>\n";
}

#-------------------------------------------------------------------------------
sub _unique {
    my %seen;
    grep !$seen{$_}++, @_;
}


#-------------------------------------------------------------------------------
1; # to make use or require happy
