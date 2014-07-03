package ConfigCase;
my $pkg_nm = 'ConfigCase';
#-----------------------------------------------------------------------------------------------
#
# SYNOPSIS
# 
#   use ConfigCase;
# 
#   # read the configuration definition xml file
#   my $cfg = ConfigCase->new("$cfgdir/ccsm_utils/Case.template/config_definition.xml");
#
#   # set some parameters
#   $cfg->set($id, $value);
# 
#   # get some parameters
#   my $value = $cfg->get($id);
#
#   # Write an xml file out
#   $cfg->write_file("$caseroot/env_run.xml", "xml");
#
#   # Write out documentation in a readme file
#   $cfg->write_doc("$caseroot/README/readme_env");
#
#   # Reset the config definition file with all of the values from the xml files
#   # in the $caseroot directory
#   $cfg->reset_setup("$caseroot/env_build.xml");
# 
# DESCRIPTION
# 
# ConfigCase objects are used to represent features of a CCSM model
# configuration that must be specified when a new case is created.
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
# is_valid_name() Returns true if the specified parameter name is contained in
#       the configuration definition file.
#
# is_ignore_name() Returns true if the specified parameter name is a name to ignore.
#
# is_valid_value() Returns true if the specified parameter name is contained in
#       the configuration definition file, and either 1) the specified value is
#       listed as a valid_value in the definition file, or 2) the definition file
#       doesn't specify the valid values.
#
# is_char() Returns true is the specified parameter name is of character type.
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
# get() Return the value of the specified configuration parameter.  Triggers
#       an exception if the parameter name is not valid.
#       ***NOTE*** If you don't want to trap exceptions, then use the query
#                  functions before calling this routine.
#
# write_file() Write an xml file.  The first argument is the
#       filename.  The second argument, if present, is the commandline of the
#       setup command that was invoked to produce the output configuration
#       file.  It is written to the output file to help document the procedure
#       used to petsetup the executable.
#
# write_doc() Write documentation on the configuration to an output README file.
#
# reset_setup() Reset with all of the values from the xml files in the $caseroot directory
# 
# COLLABORATORS
# 
# IO::File
# XML::Lite
#-----------------------------------------------------------------------------------------------
#
# Date        Author                  Modification
#-----------------------------------------------------------------------------------------------
#  2008-Aug   Mariana Vertensten      Original version
#-----------------------------------------------------------------------------------------------

use strict;
use English;
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
    $cfg->_initialize($definition_file, $default_file);

    return $cfg;
}

#-----------------------------------------------------------------------------------------------

sub is_ignore_name
{
# Return true if the requested name is a name to ignore
# These are descriptive names put in the config files -- but NOT in config_definition.xml

    my ($self, $name) = @_;

    if (    $name eq "NAME" 
	    or $name eq "SHORTNAME" 
	    or $name eq "DESC" 
	    or $name eq "VALID_GRID_MATCH" 
	    or $name eq "GRID_MATCH" 
	    or $name eq "GEN_COMPSET_MATCH"  
	    or $name eq "BEG_COMPSET_MATCH"
	    or $name eq "RES_COMPSET_MATCH"
	    or $name eq "VALID_COMPSET_MATCH"
	    or $name eq "SSTICE_COMPSET_MATCH") {
        return( 1 );
    } else {
        return( 0 );
    }
}

#-----------------------------------------------------------------------------------------------

sub is_valid_name
{
# Return true if the requested name is contained in the configuration definition.

    my ($self, $name) = @_;

    return defined($self->{$name}) ? 1 : 0;
}

#-----------------------------------------------------------------------------------------------

sub is_char
{
# Return true if the requested name is of character type

    my ($self, $name) = @_;

    if ( $self->_get_type($name) eq "char" ) {
        return( 1 );
    } else {
        return( 0 );
    }
}

#-----------------------------------------------------------------------------------------------

sub get
{
# Return requested value.

    my ($self, $name) = @_;

    defined($self->{$name}) or die "ERROR: unknown parameter name: $name\n";

    return $self->{$name}->{'value'};
}

#
# returns the value set in name with all embeded parameters resolved.
#
sub getresolved
{
    my($self,$name) = @_;

    my $val = $self->get($name);
    
    my @vars = grep(/\$([\w_]+)/,$val);
    my $v1 = $val;

    while($v1 =~ /\$([\w_]+)(.*)$/){
	my $newvar=$1;
	$v1 = $2;
	if($self->is_valid_name($newvar)){
	    my $v2=$self->getresolved($newvar);
	    $val =~ s/\$$newvar/$v2/;
	}
    }
    return $val;
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
    my $valid_values = $self->{$id}->{'valid_values'};
    if ( $valid_values ne "" ) {  # if no valid values are specified, then $value is automatically valid
	if ($is_list_value) {
	    unless (_list_value_ok($value, $valid_values)) { return 0; }
	}
	else {
	    unless (_value_ok($value, $valid_values)) { return 0; }
	}

    }

    return 1;
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
    #$self->is_valid_value($id, $value) or die
	    #"ERROR: $value is not a valid value for parameter $id: valid values are $valid_values\n";

    # Get the type description hash for the variable and check that the type is valid
    # This method throws an exception when an error is encountered.
    my %type_ref = $self->_get_typedesc($id);
    $self->validate_variable_value($id, $value, \%type_ref);

    # Check that the value is valid
    if ( $valid_values ne "" ) {
       $self->is_valid_value($id, $value) or die
	       "ERROR: $value is not a valid value for parameter $id: valid values are $valid_values\n";
    }

    # Add the new value to the object's internal data structure.
    $self->{$id}->{'value'} = $value;

    return 1;
}

#-----------------------------------------------------------------------------------------------

sub write_file
{
# Write a configuration definition file.

    my $self = shift;
    my $filename = shift;   # filepath for output namelist
    my $format = shift;

    # determine what file to write
    my @groups;
    my $xmode;

    $xmode = $self->get("XMLMODE");
    my @comps_atm  = $self->get_valid_values("COMP_ATM"); 
    my @comps_lnd  = $self->get_valid_values("COMP_LND"); 
    my @comps_ocn  = $self->get_valid_values("COMP_OCN"); 
    my @comps_ice  = $self->get_valid_values("COMP_ICE"); 
    my @comps_glc  = $self->get_valid_values("COMP_GLC"); 
    my @comps_rof  = $self->get_valid_values("COMP_ROF"); 
    my @comps_wav  = $self->get_valid_values("COMP_WAV"); 
    my @comps_cpl  = $self->get_valid_values("COMP_CPL"); 
    my @comps      = (@comps_atm, @comps_lnd, @comps_ocn, @comps_ice, @comps_glc, @comps_rof, @comps_wav, @comps_cpl); 

    if ($filename =~ "env_case") {
	@groups = qw(case_def case_desc case_comp case_mach case_der case_last);
    } elsif ($filename =~ "env_mach_pes") {
	@groups = qw(mach_pes_def mach_pes_desc mach_pes_atm mach_pes_lnd 
		     mach_pes_ice mach_pes_ocn mach_pes_cpl mach_pes_glc mach_pes_rof
		     mach_pes_wav mach_pes_stride mach_pes_last);
    } elsif ($filename =~ "env_build") {
	@groups = qw(build_macros build_def build_component build_status 
		     build_grid build_derived);
    } elsif ($filename =~ "env_run") {
       @groups = qw(run_desc run_start run_stop run_rest 
                    run_coupling run_mpi run_pio run_flags
                    run_cplhist run_mach run_din run_dout 
                    run_component run_sstice run_cesm run_domain 
                    run_dirderv run_defpts); 
    } elsif ($filename =~ "env_test") {
       @groups = qw(case_test);
    }

    my $fh;
    if ( -f $filename ) { unlink( $filename ); }
    $fh = IO::File->new($filename, '>' ) or die "can't open file: $filename\n";

    if (($format eq "xml") || ($filename =~ m/xml/)) {
    print $fh <<"EOD";
<?xml version="1.0"?>

<config_definition>

EOD
    }

    if ($filename =~ "env_case.xml" && $xmode =~ "normal") {
    print $fh <<"EOD";
<!-- ========================================================================== -->
<!--                                                                            -->
<!--       These variables CANNOT BE CHANGED once a case has been created.      -->
<!--       Invoke create_newcase again if a different grid or component         -->
<!--       combination is required.                                             -->
<!--                                                                            -->
<!-- ========================================================================== -->

EOD
}

    if ($filename =~ "env_run.xml" && $xmode =~ "normal") {
    print $fh <<"EOD";
<!-- ========================================================================== -->
<!--                                                                            -->
<!--       These variables MAY BE CHANGED ANYTIME during a run.                 -->
<!--       Additional machine speific variables that can be changed             -->
<!--       during a run are contained in the env_mach_specific file             -->
<!--                                                                            -->
<!--       Note1: users SHOULD NOT modify BUILD_COMPETE in env_build.xml        -->
<!--              this is done automatically by the scripts                     -->
<!--                                                                            -->
<!-- ========================================================================== -->

EOD
}
    if ($filename =~ "env_build.xml" && $xmode =~ "normal") {
    print $fh <<"EOD";
<!-- ========================================================================== -->
<!--                                                                            -->
<!--       These variables SHOULD NOT be changed once the model has been built. -->
<!--       Currently, these variables are not cached.                           -->
<!--                                                                            -->
<!--       Note1: users SHOULD NOT modify BUILD_COMPETE below                   -->
<!--              this is done automatically by the scripts                     -->
<!--                                                                            -->
<!-- ========================================================================== -->
EOD
}
    if($filename =~ "env_test.xml" && $xmode =~ "normal") {
    print $fh <<"EOD";
<!-- ========================================================================== -->
<!--                                                                            -->
<!--       These are the variables specific to a test case.                     -->
<!--                                                                            -->
<!-- ========================================================================== -->
EOD
}
    if ($filename =~ "env_mach_pes.xml" && $xmode =~ "normal") {
    print $fh <<"EOD";
<!-- ========================================================================== -->
<!--                                                                            -->
<!--      These variables CANNOT be modified once cesm_setup has been           -->
<!--      invoked without first invoking cesm_setup -clean.                     -->
<!--                                                                            -->
<!-- component task/thread settings                                             -->
<!-- if the user wants to change the values below after ./cesm_setup, run       -->
<!--    ./cesm_setup -clean                                                     -->
<!--    ./cesm_setup                                                            -->
<!--  to reset the pes for the run                                              -->
<!--                                                                            -->
<!--  NTASKS are the total number of MPI tasks                                  -->
<!--  NTHRDS are the number of OpenMP threads per MPI task                      -->  
<!--  ROOTPE is the global mpi task associated with the root task               -->
<!--         of that component                                                  -->     
<!--  PSTRID is the stride of MPI tasks across the global                       -->
<!--         set of pes (for now this is set to 1)                              -->
<!--  NINST is the number of instances of the component (will be spread         -->
<!--        evenly across NTASKS)                                               -->
<!--                                                                            -->
<!--  for example, for a setting with                                           -->  
<!--    NTASKS = 8                                                              -->
<!--    NTHRDS = 2                                                              -->
<!--    ROOTPE = 32                                                             -->
<!--    NINST  = 2                                                              -->
<!--  the MPI tasks would be placed starting on global pe 32                    -->
<!--  and each pe would be threaded 2-ways for this component.                  -->  
<!--  These tasks will be divided amongst both instances (4 tasks each).        -->
<!--                                                                            -->
<!--  Note: PEs that support threading never have an MPI task associated        -->
<!--    with them for performance reasons.  As a result, NTASKS and ROOTPE      -->
<!--    are relatively independent of NTHRDS and they determine                 -->
<!--    the layout of mpi processors between components.  NTHRDS is used        -->
<!--    to determine how those mpi tasks should be placed across the machine.   -->
<!--                                                                            -->
<!--  The following values should not be set by the user since they'll be       --> 
<!--  overwritten by scripts.                                                   -->
<!--    TOTALPES                                                                -->
<!--    CCSM_PCOST                                                              -->
<!--    CCSM_ESTCOST                                                            -->
<!--    PES_LEVEL                                                               -->
<!--    MAX_TASKS_PER_NODE                                                      -->
<!--    PES_PER_NODE                                                            -->
<!--    CCSM_TCOST                                                              -->
<!--    CCSM_ESTCOST                                                            -->
<!--                                                                            -->
<!--  The user can copy env_mach_pes.xml from another run, but they'll need to  -->
<!--  do the following                                                          -->
<!--    ./cesm_setup -clean                                                     -->
<!--    ./cesm_setup                                                            -->
<!--    ./CASE.build                                                            -->
<!--                                                                            -->
<!-- ========================================================================== -->

EOD
}
    foreach my $group (@groups) {
	if ($group =~ /component/) {
	    foreach my $comp ( @comps ) {
		$comp =~ s/\'//g; # get rid of quotes in $comp
		foreach my $model qw(COMP_ATM COMP_LND COMP_ICE COMP_OCN COMP_GLC COMP_ROF COMP_WAV) {
		    if ($self->get($model) eq $comp) {
			my $groupname = $group . "_$comp";
			if (($format eq "xml") || ($filename =~ m/xml/)) {
			    if ($xmode =~ "expert") {
				$self->_write_xml2($fh, $groupname);
			    } else {
				$self->_write_xml($fh, $groupname);
				if ($group =~ 'build_component') {
				    # do nothing
				} else {
				    print $fh "\n<!-- ====================================== -->";
				}
			    }
			} else {
			    $self->_write_env($fh, $groupname);
			}
		    }
		}
	    }
	} else {
	    if (($format eq "xml") || ($filename =~ m/xml/)) {
		if ($group =~ m/mach_pes_/ || $xmode =~ "expert") {
		    $self->_write_xml2($fh, $group);
		} else {
		    $self->_write_xml($fh, $group);
		    if ($group =~ 'build_component') {
			# do nothing
		    } else {
			print $fh "\n<!-- ====================================== -->";
		    }
		}
	    } else {
		$self->_write_env($fh, $group);
	    }
	}
    }
    if (($format eq "xml") || ($filename =~ m/xml/)) {
    print $fh <<"EOD";

</config_definition>
EOD
    }
}

#-----------------------------------------------------------------------------------------------

sub write_doc
{
# Write the documentation on the configuration to an output README file.

    my $self = shift;
    my $filename = shift;   # filepath for output namelist

    my $fh;
    if ( -f $filename ) { unlink( $filename ); }
    $fh = IO::File->new($filename, '>' ) or die "can't open file: $filename\n";

    my @ids = keys %$self;
    foreach my $id (sort @ids) {

	print $fh "name: $id \n";
        my $ldesc = $self->{$id}->{'ldesc'};
        if ( ! defined($ldesc) ) { $ldesc = ""; }
	print $fh "description: " . $ldesc . "\n"; 
	if ($self->{$id}->{'valid_values'}) {
	    print $fh "valid values: $self->{$id}->{'valid_values'} \n"; 
	}
	print $fh "default value: $self->{$id}->{'value'} \n"; 
	if ( $self->{$id}->{'group'} =~ "case" ) {
	    print $fh "file: env_case.xml \n";
	    print $fh "locked_stage: create_newcase \n"; 
	}
	if ( $self->{$id}->{'group'} =~ "confrun" ) {
	    print $fh "file: env_run.xml \n";
	    print $fh "locked_stage: mach_pes none \n"; 
	}
	if ( $self->{$id}->{'group'} =~ "run" ) {
	    print $fh "file: env_run.xml \n";
	    print $fh "locked_stage: none \n"; 
	}
	if ( $self->{$id}->{'group'} =~ "build" ) {
	    print $fh "file: env_build.xml \n";
	    print $fh "locked_stage: build \n"; 
	}
	print $fh " \n";
    }
}

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

        my $ldesc = "";
        my $valid = "";
        my $value = "";
        my $type  = "";

        $ldesc = $self->{$id}->{'ldesc'};
  	$valid = $self->{$id}->{'valid_values'};
	$value = $self->{$id}->{'value'};
	$type  = $self->{$id}->{'type'};

        if ( ! defined($ldesc) ) { $ldesc = ""; }
        if ( ! defined($valid) ) { $valid = ""; }
        if ( ! defined($value) ) { $value = ""; }
        if ( ! defined($type) ) { $type = ""; }

	    if ( $self->{$id}->{'group'} =~ $gid ) {
		print $fh "<row> \n";
		print $fh "<entry>$id</entry>\n";
		print $fh "<entry>$type</entry>\n";
		print $fh "<entry>$value</entry>\n";
		if ($self->{$id}->{'valid_values'}) {
		    print $fh "<entry>$ldesc [$valid] </entry>\n";
		} else {
		    print $fh "<entry>$ldesc</entry>\n";
                }
		print $fh "</row>\n";
	    }	    
    }
    print $fh "</tbody>\n";
    print $fh "</tgroup>\n";
    print $fh "</table>\n";
}

#-----------------------------------------------------------------------------------------------

sub reset_setup
{

# Reset the config object from the setup file

    my ($self, $setup_file, $setup_id) = @_;

    # Process the setup file
    my $xml = XML::Lite->new( $setup_file );
    my $root = $xml->root_element();

    # Each parameter is contained in an "entry" element.  Get all these elements.
    my @elements = ();
    @elements = $xml->elements_by_name('entry');
    
    # Loop over the elements...
    foreach my $e (@elements) {
	
	# and extract the attributes
	my %attributes = $e->get_attributes();
	
	# just get the parameter name and value
	my $id    = $attributes{'id'};
	my $value = $attributes{'value'};
	
	# set new value
	if (defined($setup_id)) {
	    if ($id ne $setup_id) {$self->set($id, $value);}
	} else {
	    $self->set($id, $value);
	}
    } # end processing setup file
}

#-----------------------------------------------------------------------------------------------
# Private methods
#-----------------------------------------------------------------------------------------------

sub _initialize
{
# Read the configuration definition file.  Create an anonymous hash with the following
# structure:
# { id1 => {type=>'ttt', value=>"xxx", list =>"XXX", valid_values=>"yyy", 
#           group=>'vvv', sdesc=>'www', ldesc=>'WWW',definition=>"zzz"},
#   id2 => {type=>'ttt', value=>"xxx", list =>"XXX", valid_values=>"yyy", 
#           group=>'vvv', sdesc=>'www', ldesc=>'WWW',definition=>"zzz"},
#   ...												     										
#   idn => {type=>'ttt', value=>"xxx", list =>"XXX", valid_values=>"yyy", 
#           group=>'vvv', sdesc=>'www', ldesc=>'WWW',definition=>"zzz"},
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
    my $index = 0;
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
        my $sdesc = $attributes{'sdesc'};
        my $ldesc = $attributes{'ldesc'};
        my $list  = $attributes{'list'};
	my $group = $attributes{'group'};
	my $type  = $attributes{'type'};
	my $valid_values = defined $attributes{'valid_values'} ? $attributes{'valid_values'} : "";

	my $i = $index++; 
	# Now add the attributes and content to the object's internal data structure.
	$self->{$id} = {'type'         => $type, 
			'value'        => $value, 
			'list'         => $list, 
			'valid_values' => $valid_values, 
			'sdesc'        => $sdesc, 
			'ldesc'        => $ldesc, 
                        'definition'   => $content, 
			'group'        => $group,
                        'index'        => $i};
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

sub _write_env
{
# Write a single line for the input env variable.

    my $self = shift;
    my $fh = shift;   # filepath for output namelist
    my $group = shift;

    print $fh <<"EOD";

EOD

    # first determine all of the indices you will write out
    my %id_indices;  
    my @ids = keys %$self;
    foreach my $id (sort @ids) {
	if ($group eq $self->{$id}->{'group'}) {
	    my $i = $self->{$id}->{'index'};
	    $id_indices{$i} = $id;
	}
    }

    # add the entry elements
    my @notsorted = keys %id_indices;
    my @sorted = sort { $a <=> $b } @notsorted;
    foreach my $i (@sorted) {
	my $id = $id_indices{$i};
	my $var = $self->{$id}->{'value'};
	if ($var =~ /\s/ ) { 
	    $var = "\"$var\"";
	} elsif ($var =~ /apos/ ) { 
	    $var = "\"$var\"";
	}
        $var =~ s/\&apos\;/\'/g;
        $var =~ s/\&lt\;/\</g;
        $var =~ s/\*/\\*/g;
	if ($group eq $self->{$id}->{'group'}) {
	    my $comment;
	    if ($self->{$id}->{'valid_values'}) {
		$comment = $self->{$id}->{'valid_values'};
	    } else {
		$comment = $self->{$id}->{'sdesc'};
	    }
	    print $fh <<"EOD";
setenv $id  $var     # $comment
EOD
        }
    }
}

#-----------------------------------------------------------------------------------------------

sub _write_xml
{
# Write a single XML element out to a file

    my $self = shift;
    my $fh = shift;   # filepath for output namelist
    my $group = shift;

    # separate the groups with spaces
    print $fh <<"EOD";

EOD

    # first determine all of the indices you will write out
    my %id_indices;  
    my @ids = keys %$self;
    foreach my $id (sort @ids) {
	if ($group eq $self->{$id}->{'group'}) {
	    my $i = $self->{$id}->{'index'};
            $id_indices{$i} = $id;
	}
    }

    # add the entry elements
    my @notsorted = keys %id_indices;
    my @sorted = sort { $a <=> $b } @notsorted;
    foreach my $i (@sorted) {
	my $id = $id_indices{$i};
	if ($self->{$id}->{'valid_values'} ne '') {
	    print $fh <<"EOD";

<!--"$self->{$id}->{'sdesc'}, valid values: $self->{$id}->{'valid_values'} ($self->{$id}->{'type'}) " -->
<entry id="$id"   value="$self->{$id}->{'value'}"  />    
EOD
} else {
		print $fh <<"EOD";

<!--"$self->{$id}->{'sdesc'} ($self->{$id}->{'type'}) " -->
<entry id="$id"   value="$self->{$id}->{'value'}"  />    
EOD
}

    }
}


#-----------------------------------------------------------------------------------------------

sub _write_xml2
{
# Write a single XML element out to a file

    my $self = shift;
    my $fh = shift;   # filepath for output namelist
    my $group = shift;

    # separate the groups with spaces
    print $fh <<"EOD";

EOD

    # first determine all of the indices you will write out
    my %id_indices;  
    my @ids = keys %$self;
    foreach my $id (sort @ids) {
	if ($group eq $self->{$id}->{'group'}) {
	    my $i = $self->{$id}->{'index'};
            $id_indices{$i} = $id;
	}
    }

    # add the entry elements
    my @notsorted = keys %id_indices;
    my @sorted = sort { $a <=> $b } @notsorted;
    foreach my $i (@sorted) {
	my $id = $id_indices{$i};
	print $fh <<"EOD";
<entry id="$id"   value="$self->{$id}->{'value'}"  />    
EOD
    }
}


#-----------------------------------------------------------------------------------------------

sub _list_value_ok
{
# Check that all input values ($values_in may be a comma separated list)
# are contained in the comma separated list of valid values ($valid_values).
# Return 1 (true) if all input values are valid, Otherwise return 0 (false).

    my ($values_in, $valid_values) = @_;

    my @values = split /,/, $values_in;

    my $num_vals = scalar(@values);
    my $values_ok = 0;

    foreach my $value (@values) {

	if (_value_ok($value, $valid_values)) { ++$values_ok; }

    }

    ($num_vals == $values_ok) ? return 1 : return 0;
}

#-----------------------------------------------------------------------------------------------

sub _value_ok
{

# Check that the input value is contained in the comma separated list of
# valid values ($valid_values).  Return 1 (true) if input value is valid,
# Otherwise return 0 (false).

    my ($value, $valid_values) = @_;

    # If the valid value list is null, all values are valid.
    unless ($valid_values) { return 1; }

    my @expect = split /,/, $valid_values;

    $value =~ s/^\s+//;
    $value =~ s/\s+$//;
    foreach my $expect (@expect) {
	if ($value =~ /^$expect$/) { return 1; }
    }

    return 0;
}

#-----------------------------------------------------------------------------------------------

sub get_var_type
{
# Return 'type' attribute for requested variable

    my ($self, $name) = @_;

    return $self->{$name}->{'type'};
}

#-----------------------------------------------------------------------------------------------

sub get_valid_values
{
# Return list of valid_values as an array for requested variable
# To return without quotes use the 'noquotes'=>1 option.
    my ($self, $name, %opts) = @_;

    my $valid_values = $self->{$name}->{'valid_values'};
    my $type = $self->{$name}->{'type'};
    my @values = split( /,/, $valid_values );

    # if string type and NOT noquote option and have a list -- add quotes around values
    if ( ! defined($opts{'noquotes'}) || ! $opts{'noquotes'} ) {
       if ( $#values > -1 && ($type eq "char") ) {
          for( my $i=0; $i <= $#values; $i++ ) {
             $values[$i] = "'$values[$i]'";
          }
       }
    }
    return( @values );
}

#-----------------------------------------------------------------------------------------------

sub _get_type
{
# Return 'type' attribute for requested variable

    my ($self, $name) = @_;

    return $self->{$name}->{'type'};
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
	die "ERROR: in $nm (package $pkg_nm): datatype $type_def is NOT valid\n";
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

#-----------------------------------------------------------------------------------------------
# Perl regular expressions to match Fortran namelist tokens.

# Variable names.
# % for derived types, () for arrays
my $varmatch = "[A-Za-z_]{1,31}[A-Za-z0-9_]{0,30}[(0-9)%a-zA-Z_]*";

# Integer data.
my $valint = "[+-]?[0-9]+";
my $valint_repeat = "${valint}\\*$valint";

# Logical data.
my $vallogical1 = "[Tt][Rr][Uu][Ee]";
my $vallogical2 = "[Ff][Aa][Ll][Ss][Ee]";
my $vallogical = "$vallogical1|$vallogical2";
my $vallogical_repeat = "${valint}\\*$vallogical1|${valint}\\*$vallogical2";

# Real data.
# "_" are for f90 precision specification
my $valreal1 = "[+-]?[0-9]*\\.[0-9]+[EedDqQ]?[0-9+-]*";
my $valreal2 = "[+-]?[0-9]+\\.?[EedDqQ]?[0-9+-]*";
my $valreal = "$valreal1|$valreal2";
my $valreal_repeat = "${valint}\\*$valreal1|${valint}\\*$valreal2";

# Match for all valid data-types: integer, real or logical
# note: valreal MUST come before valint in this string to prevent integer portion of real 
#       being separated from decimal portion
my $valall = "$vallogical|$valreal|$valint";

# Match for all valid data-types with repeater: integer, real, logical, or string data
# note: valrepeat MUST come before valall in this string to prevent integer repeat specifier 
#       being accepted alone as a value
my $valrepeat = "$vallogical_repeat|$valreal_repeat|$valint_repeat";

# Match for all valid data-types with or without numberic repeater at the lead
my $valmatch = "$valrepeat|$valall";

# Same as above when a match isn't required
my $nrvalmatch = $valmatch. "||";

#-----------------------------------------------------------------------------------------------

sub validate_variable_value
#
# Validate that a given value matches the expected input type definition
# Expected description of keys for the input type hash is:
#      type           type description (char, logical, integer, or real)       (string)
#      strlen         Length of string (if type char)                          (integer)
#      validValues    Reference to array of valid values                       (string)
#
{
  my ($self, $var, $value, $type_ref) = @_;
  my $nm = "validate_variable_value";

  # Ensure type hash has required variables
  if ( ref($type_ref) !~ /HASH/ ) {
      die "ERROR: in $nm : Input type is not a HASH reference.\n";
  }
  foreach my $item ( "type", "validValues", "strlen" ) {
     if ( ! exists($$type_ref{$item}) ) {
         die "ERROR: in $nm: Variable name $item not defined in input type hash.\n";
     }
  }

  # If string check that less than defined string length
  my $str_len = 0;
  if ( $$type_ref{'type'} eq "char" ) {
      $str_len = $$type_ref{'strlen'};
      if ( length($value) > $str_len ) {
         die "ERROR: in $nm Variable name $var " .
             "has a string element that is too long: $value\n";
      }
  }
  # If not string -- check that array size is smaller than definition
  my @values;
  if ( $str_len == 0 ) {
      @values = split( /,/, $value );
      # Now check that values are correct for the given type
      foreach my $i ( @values ) {
	  my $compare;
	  if (      $$type_ref{'type'} eq "logical" ) {
	      $compare = $vallogical;
	  } elsif ( $$type_ref{'type'} eq "integer" ) {
	      $compare = $valint;
	  } elsif ( $$type_ref{'type'} eq "real" ) {
	      $compare = $valreal;
	  } else {
	      die "ERROR: in $nm (package $pkg_nm): Type of variable name $var is " . 
		  "not a valid FORTRAN type (logical, integer, real, or char).\n";
	  }
	  if ( $i !~ /^\s*(${compare})$/ ) {
	      die "ERROR: in $nm (package $pkg_nm): Variable name $var " .
		  "has a value ($i) that is not a valid type " . $$type_ref{'type'} . "\n";
	  }
      }
  }
}

sub set_machine
{
    # Set the parameters for the specified machine.  The
    # parameters are read from an input file, and if no machine matches are
    # found then issue error message.
    # This routine uses the configuration defined at the package level ($cfg_ref).

    my ($self, $machine_file, $machine, $print) = @_;
    my $xml = XML::Lite->new( $machine_file );
    my $root = $xml->root_element();

    # Check for valid root node
    my $name = $root->get_name();
    $name eq "config_machines" or die
	"file $machine_file is not a machine parameters file\n";

    # Read the machine parameters from $machine_file.
    my @e = $xml->elements_by_name( "machine" );
    my %a = ();

    my $mach=$machine;
    if($machine =~ /(.*)_/){
       $mach=$1;
    }
  # Search for matching compset.
    my $found = 0;
    my @mach_settings;
  MACHINE:
    while ( my $e = shift @e ) {
	%a = $e->get_attributes();
	if ( ($machine eq $a{'MACH'}) or ($mach eq $a{MACH})) {
	    $found = 1;
	    @mach_settings = $e->get_children();
	    last MACHINE;
	}
    }

    # Die unless search was successful.
    unless ($found) { 
	print "set_machine: no match for machine $machine - possible machine values are \n";
	print_machines( $machine_file );
	die "set_machine: exiting\n"; 
    }

    # Loop through all entry_ids of the $self object and if the corresponding 
    # attributed is defined in the compset hash, then reset the self object to
    # that value


    my $setting;
    foreach $setting (@mach_settings){
	my $sname = $setting->get_name();
	next if($self->is_ignore_name($sname)) ;
	if ( ! $self->is_valid_name($sname) ) { 
	    die "set_machine: invalid id $sname in machine $machine in file $machine_file exiting\n"; 
	}
        # allow for environment variables in the config_machines.xml file 
        # using $ENV{variablename} syntax
	my $text = $setting->get_text();
	if($text =~/^(.*)\$ENV{(.*)}(.*)$/){
	    $text = $1.$ENV{$2}.$3;
	}
	$self->set($sname,$text);
	

	print "cfg_ref: $sname set to ".$self->get($sname)."  $text\n" if($print==2);
    }
    $self->set('MACH',$machine);



}

#-------------------------------------------------------------------------------
sub print_machines
{
    # Print all currently supported machines
    my ($machine_file) = @_;
    my $xml = XML::Lite->new( $machine_file );
    my $root = $xml->root_element();

    # Check for valid root node
    my $name = $root->get_name();
    $name eq "config_machines" or die
	"file $machine_file is not a machine parameters file $name\n";

    print ("  \n");
    print ("  MACHINES:  name (description)\n");
    
    # Read the grid parameters from $grid_file.
    my @e = $xml->elements_by_name( "machine" );
    my %a = ();
    while ( my $e = shift @e ) {
	%a = $e->get_attributes();
	my @children = $e->get_children();
	my $child;
	foreach $child (@children){
	    if($child->get_name == "DESC"){
		my $desc = $child->get_text();
		print "    $a{'MACH'} ($desc) \n";		
		last;
	    }
	}
    }
}

#-----------------------------------------------------------------------------------------------
sub set_pes
{
    # Set the parameters for the pe layout.    
    my($self,$pes_file,$decomp_ref,$pecount_opts, $print) = @_;

#    print "tcx1: pecount_opts: $pecount_opts\n";

    # Initialize some local variables
    my $nm = "set_pes"; 
    my @matches = keys(%$self);
    if ( ref($decomp_ref) ne "HASH" ) { die "ERROR::($nm) input decomp is not a hash!\n"; }

    # Open and read the xml file
    my $xml = XML::Lite->new( $pes_file );
    if ( ! defined($xml) ) {
	die "ERROR($nm): Trouble opening or reading $pes_file\n";
    }
    my $elm  = $xml->root_element( );
    my $root = "pesinfo";
#    my $root = "config_definition";

    my @children = $xml->elements_by_name( $root );
    if ( $#children < 0 ) {
	die "ERROR($nm): could not find the main $root root element in $pes_file\n";
    }
    if ( $#children != 0 ) {
	die "ERROR($nm): $root root element in $pes_file is duplicated, there should only be one\n";
    }
    
    # examine the attributes of each element to determine the "best fit"
    my $possible_match ;
    my $pecount_vals = " ";
    for (my $i = 0; $i <= $#children; $i++) {
	#
	# Name of element, and it's associated attributes
	my $child = $children[$i];
	my $name = $child->get_name();
	my @children_level = $child->get_children();
	my $num_children = $#children_level+1;
	
	if ( $#children_level > -1 ) {
	    foreach my $child_level ( @children_level ) {
		
		# Check all the attributes for this element to determine if we have a complete match
		my %atts = $child_level->get_attributes;
		my @keys = keys(%atts);
		my $num_matches = 0;
		foreach my $key ( @keys ) {
		    foreach my $var ( @matches ) {
			my $match = $atts{$key};
                        # Match either the begining or the end of the line
			if ( ($key eq $var)){ 
			   if($var eq "GRID"  && (($self->get($var) =~ /$match/ ))){
			       $num_matches++; 
			   }elsif ($var eq "CCSM_LCOMPSET" &&
				   $self->get($var) =~ /$match/ ) {
			       # For CCSM_LCOMPSET, we consider it a match if it matches any part of
			       # the compset name (similar to GEN_COMPSET_MATCH in config_compsets.xml)
			       $num_matches++;
			   }elsif ($self->get($var) =~ /^$match/ ) {
			       $num_matches++; 
			   }
			   
			}
		    }

		    if ($key eq "PECOUNT") {
			if($pecount_opts =~ /^$atts{$key}$/){
			    $num_matches++;
			}
		    }

		}
		




		# Need all the attributes to match in order to read the element pes
		my $num_keys = $#keys + 1;

		if ($num_matches eq $num_keys) {

		    my %atts = $child_level->get_attributes;
		    my @keys = keys(%atts);


		    foreach my $key ( @keys ) {
			$possible_match->{$key} = $atts{$key};
		    }

		    my @decomp_children = $child_level->get_children();
		    if ( $#decomp_children > -1 ) {
			foreach my $decomp_child ( @decomp_children ) {
			    my $name  = $decomp_child->get_name();
			    my $value = $decomp_child->get_text();
                            # Use exists over defined, as MPILIB exists but set to undef by default
			    if ( ! exists($$decomp_ref{$name}) ) {
				die "ERROR($nm): $name is NOT a valid element for the decomp output hash\n";
			    }
			    $$decomp_ref{$name} = $value;
			}
		    } else {
			die "ERROR($nm): No sub-elements for $name \n";
		    }
		}  
	    }
	}
    }


    print "The PE layout for this case match these options:\n";
    foreach my $key (keys %$possible_match){
      print "$key =  $possible_match->{$key}\n";
  }
}

1; # to make use or require happy
