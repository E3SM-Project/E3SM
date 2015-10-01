package ConfigMachine;
my $pkg_nm = 'ConfigMachine';

use strict;
use English;
use Cwd qw( getcwd abs_path chdir);
use IO::File;
use XML::LibXML;
use File::Basename;
use Log::Log4perl qw(get_logger);
my $logger;

BEGIN{
    $logger = get_logger();
}

#-----------------------------------------------------------------------------------------------
sub setMachineFile
{
    my ($machine, $model, $cimeroot, $input_file) = @_;

    # Determine the machines specificaiton file and see if the
    # target machine is supported

    my $xml = XML::LibXML->new( no_blanks => 1)->parse_file($input_file);
    my @nodes = $xml->findnodes(".//entry[\@id=\"MACHINES_SPEC_FILE\"]/default_value");
    my $machines_file = $nodes[0]->textContent();
    $machines_file =~ s/\$CIMEROOT/$cimeroot/;
    $machines_file =~ s/\$MODEL/$model/;
    unless(-f "$machines_file"){
	$logger->logdie("*** Cannot find supported machines file $machines_file ***\n");
	exit 1;
    }
    if ($machine =~ /(.*)_(.*)/){
	$machine = $1;
    }

    my $machxml = XML::LibXML->new( no_blanks => 1)->parse_file($machines_file);
    my @machnodes = $machxml->findnodes(".//machine[\@MACH=\"$machine\"]");
    if (@machnodes) {
	$logger->info("Found machine \"$machine\" in $machines_file \n");
    } else {
	$logger->error( "ERROR ConfigMachine::setMachineFile: no match for machine $machine 
	                 - possible machine values are: \n");
	listMachines( "$machines_file" );
	return 0;
    }	    
    return $machines_file;
}

#-----------------------------------------------------------------------------------------------
sub setMachineValues
{
    # Set the parameters for the specified machine.  
    my ( $file_config, $primary_component, $machine, $config) = @_;

    my $model = $config->get('MODEL');
    my $machines_file = $config->get('MACHINES_SPEC_FILE');
    $machines_file =~ s/\$MODEL/$model/;
    if(! -f "$machines_file") {
	$logger->error("*** Cannot find supported machines file $machines_file ***\n");
	exit 1;
    }
    my $machines_dir  = dirname($machines_file);

    # First Check if the target machine is supported
    if ($machine =~ /(.*)_(.*)/){
	$machine  = $1;
	$config->set('COMPILER', "$2");
    }

    my $xml = XML::LibXML->new( no_blanks => 1)->parse_file($machines_file);
    my @nodes = $xml->findnodes(".//machine[\@MACH=\"$machine\"]");
    if (@nodes) {
	$logger->info( "Found machine \"$machine\" in $machines_file ");
    } else {
	$logger->fatal("ERROR ConfigMachine::setMachineValues: no match for machine $machine :");
	$logger->fatal("  - possible machine values are ");
	listMachines( "$machines_file" );
	$logger->logdie( "Exiting ");
    }	    
    my $compiler = $config->get('COMPILER');
    my $mpilib = $config->get('MPILIB');

    if (!defined $compiler) {
	my @nodes = $xml->findnodes(".//machine[\@MACH=\"$machine\"]/COMPILERS");
	my $compilers = $nodes[0]->textContent();
	my @compilers = split(/,/,$compilers);
	$compiler = $compilers[0];
	$config->set('COMPILER', "$compiler");
    }	
    if (!defined $mpilib) {
	my @nodes = $xml->findnodes(".//machine[\@MACH=\"$machine\"]/MPILIBS");
	my $mpilibs = $nodes[0]->textContent();
	my @mpilibs = split(/,/,$mpilibs);
	$mpilib = $mpilibs[0];
	$config->set('MPILIB', "$mpilib");
    }	
    $config->set('MACH'         ,  $machine);
    $config->set('MACHINES_FILE', "$machines_file");
    $config->set('MACHDIR'      , "$machines_dir");

    # Set the machine values obtained from the $machines_file
    _set_machine_values($config);

    # Check that compiler request for target machine matches a supported value
    # Or set default compiler - if not provided compiler request
    _check_machine_compilers($config);

    $logger->info( "Machine specifier: $machine.\n");

    # Determine pio settings for target machine
    # Note that any pio settings that are grid or compset dependent will be overwritten
    # by the config_pio.xml settings for the primary component
    _setPIOsettings($file_config, $primary_component, $config);

    $logger->info("Set pio settings for $machine.\n");
}

#-------------------------------------------------------------------------------
sub listMachines
{
    # Print the list of supported machines
    my ($machine_file) = @_;

    my $parser = XML::LibXML->new( no_blanks => 1);
    my $xml = $parser->parse_file($machine_file);

    $logger->warn ("  MACHINES:  name (description)\n");

    foreach my $node ($xml->findnodes(".//machine")) {
	my $name = $node->getAttribute('MACH');
	foreach my $child ($node->findnodes("./*")) {
	    if ($child->nodeName() eq 'DESC') {
		my $desc = $child->textContent();
		$logger->warn( "    $name ($desc) \n");		
	    }
	}
    }
}

#-----------------------------------------------------------------------------------------------
#                               Private routines
#-----------------------------------------------------------------------------------------------
sub _set_machine_values
{
    # open the specified xml file
    my ($config) = @_;

    my $machine       = $config->get('MACH'); 
    my $machines_file = $config->get('MACHINES_FILE');

    my $xml = XML::LibXML->new( no_blanks => 1)->parse_file($machines_file);
    my @machine_nodes = $xml->findnodes(".//machine[\@MACH=\"$machine\"]/*");
    if (@machine_nodes) 
    {
	foreach my $node (@machine_nodes) {
	    my $name  = $node->nodeName();
	    next if($name eq "mpirun");
	    next if($name eq "module_system");
	    my $value = $node->textContent();
	    if ( ! $config->is_valid_name($name) ) { 
		$logger->logdie("set_machine: invalid id $name in machine $machine file $machines_file exiting\n");
		return;
	    }
	    # allow for environment variables in the config_machines.xml file using $ENV{variablename} syntax
	    if ($value =~/^(.*)\$ENV{([^}]*)}(.*)$/){
		if(defined $ENV{$2}){
		    $value = $ENV{$2};
		    $value = $1.$value if (defined $1);
		    $value .= $3 if (defined $3);
		}else{
		    $logger->warn( "No environment setting found for $2 in $name=$value");
		}
	    }
	    if($machine eq "userdefined" && $value =~ /USERDEFINED/){
		$logger->warn("Value for $name not set");
	    }else{
		$config->set($name, $value);
	    }
	    $logger->debug( "config: $name set to ".$config->get($name)."  $value\n" );
	}
    } 
    else 
    {
	$logger->logdie( "ERROR: ConfigMachine::_set_machine_values: no specifications contained for machine $machine :\n");
	return; 
    }
}

#-------------------------------------------------------------------------------
sub _check_machine_compilers
{
    # Check that compiler request for target machine matches a supported value
    # Or set default compiler - if not provided compiler request

    my ( $config) = @_;

    my $machine   = $config->get('MACH'); 
    my $caseroot  = $config->get('CASEROOT');
    my $compiler  = $config->get('COMPILER');

    my $compilers;
    if ($machine =~ /userdefined/){
	$config->set('COMPILER', "USERDEFINED_required_build");
    } else { 
	$compilers = $config->get('COMPILERS');
	my @compilers = split ",", $compilers, -1;
	if ($compiler) {
	    if (! ($machine =~ "generic")){
		my $found = 0;
		foreach my $comp (@compilers) {
		    if ($compiler eq $comp) {
			$found = 1;
		    }
		}
		if (!$found) {
		    my $sysmod = "rm -rf $caseroot";
		    unless(system($sysmod) == 0) {
			$logger->fatalxs("ERROR: $sysmod failed: $?\n");
		    }
		    $logger->logdie( "ERROR: compiler setting of $compiler does not match supported values of $compilers \n");
		    return;
		}
	    }
	    $config->set('COMPILER', "$compiler");
	    $logger->info( "Machine compiler specifier: $compiler\n");
	} else {
	    $compiler = $compilers[0];   
	    $config->set('COMPILER', "$compiler");
	    $logger->info("Machine compiler specifier: $compiler\n");
	}
    }
}

#-------------------------------------------------------------------------------
sub _setPIOsettings
{
    # Set pio settings from config_machines.xml and config_pio.xml file

    my ($file_config, $primary_component, $config) = @_; 

    my $model    = $config->get('MODEL');
    my $cimeroot = $config->get('CIMEROOT');
 
    # read the PIO_SPEC_FILE file for grid and/or compset specific pio settings
    my $file = $config->get('PIO_SPEC_FILE');
    $file =~ s/\$MODEL/$model/;
    $file =~ s/\$CIMEROOT/$cimeroot/;
    my $xml = XML::LibXML->new( no_blanks => 1)->parse_file($file);

    foreach my $entry ($xml->findnodes(".//entry")) {
	my $id = $entry->getAttribute('id');
	my $default;

	# Loop over each value node for the given entry node
	my $index_match   = -1;
	my $value_counter = 0;
	my $max_matches   = 0;
	my @values = $entry->findnodes("./values/value");
	foreach my $value (@values) {
	    # Loop over each attribute in the value node and make sure they all match
	    # If they do increment the match counter
	    $value_counter++;
	    my $attr_value;
	    my $matches = 0;
	    my @attributes = $value->attributes();
            #
	    if($#attributes<0){
		$default = $value->textContent();
	    }
	    MATCH: foreach my $attr (@attributes) {
		$attr_value = $attr->value(); 
		my $attr_name  =uc $attr->name();
		my $target = $config->get( $attr_name);
		if (($attr_value =~ /^\!/) && ($target !~ m/$attr_value/)) {
		    $matches++;
		} elsif ($target =~ m/$attr_value/) {
		    $matches++;
		} else {
		    # Not all attributes match - reset match counter to 0 
		    # and go to next value node
		    $matches = 0;
		    last MATCH;
		}
	    }
	    # If the number of matches is greater than $max_matches - then
	    # reset $max_matches and set new value for $index_match
	    if ($matches > $max_matches) {
		$max_matches = $matches;
		$index_match = $value_counter-1;
	    }
	}

	# If $index_match is greater than -1, reset the value of $id
	if ($index_match > -1) {
	    my $newval = $values[$index_match]->textContent();
	    $config->set($id, $newval);
	}elsif(defined $default){
	    $config->set($id, $default);
	}

    }
}

1; # to make use or require happy
