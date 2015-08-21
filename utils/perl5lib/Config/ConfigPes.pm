package ConfigPes;
my $pkg_nm = 'ConfigPes';

use strict;
use English;
use IO::File;
use XML::LibXML;
use Data::Dumper;
use ConfigCase;

# Check for the existence of XML::LibXML in whatever perl distribution happens to be in use.
# If not found, print a warning message then exit.
eval {
    require XML::LibXML;
    XML::LibXML->import();
};
if($@)
{
    my $warning = <<END;
WARNING:
  The perl module XML::LibXML is needed for XML parsing in the CESM script system.
  Please contact your local systems administrators or IT staff and have them install it for
  you, or install the module locally.  

END
    print "$warning\n";
    exit(1);
}

#-----------------------------------------------------------------------------------------------
sub setPes {

    # Set the parameters for the pe layout.    

    my($pesize_opts, $primary_component, $opts_ref, $config) = @_;

    if (defined $$opts_ref{'pes_file'}) {

	# Reset the pes if a pes file is specified
	my $pes_file = $$opts_ref{'pes_file'};
	(-f "$pes_file")  or  die "** Cannot find pes_file \"$pes_file\" ***\n";
	$config->reset_setup("$pes_file");    

    } else {

	if ($pesize_opts =~ m!^([0-9]+)D?$!) {

	    my $ntasks = $1;
	    my $nthrds = 1;
	    my $rootpe =  0;
	    _setPESmatch1($pesize_opts, $ntasks, $nthrds, $rootpe, $config);

	} elsif ($pesize_opts =~ m!^([0-9]+)x([0-9]+)D?$!) {

	    my $ntasks = $1;
	    my $nthrds = $2;
	    my $rootpe =  0;
	    _setPESmatch2($pesize_opts, $ntasks, $nthrds, $rootpe, $config);

	} else {

	    # Determine pe layout settings
	    _setPESsettings($pesize_opts, $config);
	}
    }
}

#-----------------------------------------------------------------------------------------------
#                               Private routines
#-----------------------------------------------------------------------------------------------
sub _setPESmatch1
{
    my ($pesize_opts, $config) = @_; 

    my $ntasks = $1;
    my $nthrds = 1;
    my $root = 0;

    foreach my $comp ("ATM", "LND", "ICE", "OCN", "GLC", "WAV", "ROF", "CPL" ) {
	$config->set("NTASKS_$comp") = $ntasks; 
	$config->set("NTHRDS_$comp") = $nthrds; 
    }

    my $root;
    if ($pesize_opts =~ m!^([0-9]+)D$!) {
	$root = 0          ; $config->set('ROOTPE_ATM') = $root;
	$root = 1 * $ntasks; $config->set('ROOTPE_LND') = $root;
	$root = 2 * $ntasks; $config->set('ROOTPE_OCN') = $root;
	$root = 3 * $ntasks; $config->set('ROOTPE_ICE') = $root;
	$root = 4 * $ntasks; $config->set('ROOTPE_GLC') = $root;
	$root = 5 * $ntasks; $config->set('ROOTPE_WAV') = $root;
	$root = 6 * $ntasks; $config->set('ROOTPE_ROF') = $root;
	$root = 7 * $ntasks; $config->set('ROOTPE_CPL') = $root;
    }

}

#-------------------------------------------------------------------------------
sub _setPESmatch2
{
    my ($pesize_opts, $ntasks, $nthrds, $rootpe, $config) = @_; 

    foreach my $comp ("ATM", "LND", "ICE", "OCN", "GLC", "WAV", "ROF", "CPL") {
	$config->set("NTASKS_$comp") = $ntasks; 
	$config->set("NTHRDS_$comp") = $nthrds; 
    }
    
    my $root;
    if ($pesize_opts =~ m!^([0-9]+)x([0-9]+)D$!) {
	$root = 0          ; $config->set('ROOTPE_ATM') = $root;
	$root = 1 * $ntasks; $config->set('ROOTPE_LND') = $root;
	$root = 2 * $ntasks; $config->set('ROOTPE_OCN') = $root;
	$root = 3 * $ntasks; $config->set('ROOTPE_ICE') = $root;
	$root = 4 * $ntasks; $config->set('ROOTPE_GLC') = $root;
	$root = 5 * $ntasks; $config->set('ROOTPE_WAV') = $root;
	$root = 6 * $ntasks; $config->set('ROOTPE_ROF') = $root;
	$root = 7 * $ntasks; $config->set('ROOTPE_CPL') = $root;
    }
}    

#-------------------------------------------------------------------------------
sub _setPESsettings
{
    # Read xml file and obtain NTASKS, NTHRDS, ROOTPE and NINST for each component
    my ($pesize_opts, $config) = @_; 

    my $mach		= $config->get('MACH');
    my $grid_longname   = $config->get('GRID');
    my $compset_name	= $config->get('COMPSET');
    my $pes_file	= $config->get('PES_SPEC_FILE');

    # temporary hash
    my %decomp;
    $decomp{'NTASKS_ATM'} = 16;
    $decomp{'NTASKS_LND'} = 16;
    $decomp{'NTASKS_ICE'} = 16;
    $decomp{'NTASKS_OCN'} = 16;
    $decomp{'NTASKS_ROF'} = 16;
    $decomp{'NTASKS_GLC'} = 16;
    $decomp{'NTASKS_WAV'} = 16;
    $decomp{'NTASKS_CPL'} = 16;

    $decomp{'NTHRDS_ATM'} = 1;
    $decomp{'NTHRDS_LND'} = 1;
    $decomp{'NTHRDS_ICE'} = 1;
    $decomp{'NTHRDS_OCN'} = 1;
    $decomp{'NTHRDS_ROF'} = 1;
    $decomp{'NTHRDS_GLC'} = 1;
    $decomp{'NTHRDS_WAV'} = 1;
    $decomp{'NTHRDS_CPL'} = 1;

    $decomp{'ROOTPE_ATM'} = 0;
    $decomp{'ROOTPE_LND'} = 0;
    $decomp{'ROOTPE_ICE'} = 0;
    $decomp{'ROOTPE_OCN'} = 0;
    $decomp{'ROOTPE_ROF'} = 0;
    $decomp{'ROOTPE_GLC'} = 0;
    $decomp{'ROOTPE_WAV'} = 0;
    $decomp{'ROOTPE_CPL'} = 0;

    # --------------------------------------
    # look for model grid / model machine match
    # --------------------------------------

    # Parse the pes xml file
    my $xml = XML::LibXML->new( no_blanks => 1)->parse_file($pes_file);

    # First determine that the xml file has consistent attributes
    foreach my $grid ($xml->findnodes(".//grid")) {
	if (! defined $grid->getAttribute('name')) {
	    my $name = $grid->nodeName();
	    die "\n ERROR ConfigPes.pm: node <$name> does not have a required attribute of \"name\" in $pes_file \n\n";
	}
    }
    foreach my $mach ($xml->findnodes(".//mach")) {
	if (! defined $mach->getAttribute('name')) {
	    my $name = $mach->nodeName();
	    die "\n ERROR ConfigPes.pm: node <$name> does not have a required attribute of \"name\" in $pes_file \n\n";
	}
    }
    foreach my $pes ($xml->findnodes(".//pes")) {
	if (! defined $pes->getAttribute('pesize')) {
	    my $name = $pes->nodeName();
	    die "\n ERROR ConfigPes.pm: node <$name> does not have a required attribute of \"pesize\" in $pes_file \n\n";
	}
	if (! defined $pes->getAttribute('compset')) {
	    my $name = $pes->nodeName();
	    die "\n ERROR ConfigPes.pm: node <$name> does not have a required attribute of \"compset\" in $pes_file \n\n";
	}
    }

    # Set the match variables $mach_set and $grid_set
    # Determine possible matches in following order:
    # (1) model grid and model machine theh (2) model grid and 'any' machine
    
    my $mach_set;
    my $grid_set;

    # Determine if grid and machine node names that match
    my @pes = $xml->findnodes(".//grid/mach[not(\@name =\"any\")]");
    if (@pes) {
	foreach my $node (@pes) {
	    my $grid_match = $node->parentNode()->getAttribute('name');
	    my $mach_match = $node->getAttribute('name');
	    if ($grid_longname =~ m/$grid_match/ && $mach_match =~ m/$mach/) {
		$grid_set = $grid_match;
		$mach_set = $mach_match;
		last;
	    }
	}
    }
    if (! defined $grid_set && ! defined $mach_set) {
	$mach_set = 'any';
	@pes = $xml->findnodes(".//grid/mach[\@name=\"any\"]");
	if (@pes) {
	    foreach my $node (@pes) {
		my $grid_match = $node->parentNode()->getAttribute('name');
		if ($grid_longname =~ m/$grid_match/) {
		    $grid_set = $grid_match;
		    $mach_set = 'any';
		    last;
		}
	    }
	}
    }
    if (! defined $grid_set) {
	$grid_set = 'any';
    }

    print "ConfigPes: grid longname is $grid_longname \n";
    print "ConfigPES: grid match    is $grid_set \n";
    print "ConfigPES: machine match is $mach_set\n";

    if ($mach_set eq 'any' && $grid_set eq 'any') {
	# set the pes for the MAX_TASKS_PER_NODE in config_machines.xml
	my $max_tasks_per_node = $config->get('MAX_TASKS_PER_NODE');
	foreach my $name ('NTASKS_ATM','NTASKS_LND','NTASKS_ROF','NTASKS_ICE','NTASKS_OCN','NTASKS_GLC','NTASKS_GLC','NTASKS_WAV','NTASKS_CPL') {
	    $decomp{$name} = $max_tasks_per_node;
	}
	foreach my $name ('NTHRDS_ATM','NTHRDS_LND','NTHRDS_ROF','NTHRDS_ICE','NTHRDS_OCN','NTHRDS_GLC','NTHRDS_GLC','NTHRDS_WAV','NTHRDS_CPL') {
	    $decomp{$name} = 1;
	}
	foreach my $name ('ROOTPE_ATM','ROOTPE_LND','ROOTPE_ROF','ROOTPE_ICE','ROOTPE_OCN','ROOTPE_GLC','ROOTPE_GLC','ROOTPE_WAV','ROOTPE_CPL') {
	    $decomp{$name} = 0;
	}
    } else {
	# Now search the pes
	my $pe_select;
	my $compset_match;
	my $pesize_match;
	my @nodes = $xml->findnodes(".//grid[\@name=\"$grid_set\"]/mach[\@name=\"$mach_set\"]/*");
	if ((! defined $compset_match) && (! defined $pesize_match)) {
	    foreach my $node (@nodes) {
		my $compset_attr = $node->getAttribute('compset');
		my $pesize_attr  = $node->getAttribute('pesize');
		if (($compset_name =~ m/$compset_attr/) && ($pesize_opts =~ m/$pesize_attr/)) {
		    $pe_select = $node;
		    $compset_match = $compset_attr;
		    $pesize_match  = $pesize_attr;
		    last;
		} 
	    }
	}
	if ((! defined $compset_match) && (! defined $pesize_match)) {
	    foreach my $node (@nodes) {
		my $compset_attr = $node->getAttribute('compset');
		my $pesize_attr  = $node->getAttribute('pesize');
		if (($compset_name =~ m/$compset_attr/) && ($pesize_attr eq 'any')) {
		    $pe_select = $node;
		    $compset_match = $compset_attr;
		    $pesize_match  = $pesize_attr;
		    last;
		}
	    }
	} 
	if ((! defined $compset_match) && (! defined $pesize_match)) {
	    foreach my $node (@nodes) {
		my $compset_attr = $node->getAttribute('compset');
		my $pesize_attr  = $node->getAttribute('pesize');
		if (($compset_name eq 'any' ) && ($pesize_opts =~ m/$pesize_attr/)) {
		    $pe_select = $node;
		    $compset_match = $compset_attr;
		    $pesize_match  = $pesize_attr;
		    last;
		}
	    }
	}
	if ((! defined $compset_match) && (! defined $pesize_match)) {
	    foreach my $node (@nodes) {
		my $compset_attr = $node->getAttribute('compset');
		my $pesize_attr  = $node->getAttribute('pesize');
		if (($compset_attr eq 'any' ) && ($pesize_attr eq 'any')) {
		    $pe_select = $node;
		    $compset_match = $compset_attr;
		    $pesize_match  = $pesize_attr;
		    last;
		}
	    }
	}
	if (! defined $pe_select) {
	    die "ERROR: no pes match found in $pes_file \n";
	} else {
	    print "ConfigPES: compset_match is $compset_match\n"; 
	    print "ConfigPES: pesize match  is $pesize_match \n"; 
	}

	my @pes_ntasks = $pe_select->findnodes("./ntasks"); 
	my @pes_nthrds = $pe_select->findnodes("./nthrds"); 
	my @pes_rootpe = $pe_select->findnodes("./rootpe"); 

	foreach my $pes (@pes_ntasks, @pes_nthrds, @pes_rootpe) {
	    my @children = $pes ->childNodes();
	    foreach my $child (@children) {
		my $name  = uc $child->nodeName(); 
		my $value =    $child->textContent();
		$decomp{$name}  = $value;
	    }
	}
    }

    foreach my $comp ("ATM", "LND", "ICE", "OCN", "GLC", "WAV", "ROF", "CPL") {
	my $ntasks = _clean($decomp{"NTASKS_$comp"});
	my $nthrds = _clean($decomp{"NTHRDS_$comp"});
	my $rootpe = _clean($decomp{"ROOTPE_$comp"});

	$config->set("NTASKS_$comp" , $ntasks);
	$config->set("NTHRDS_$comp" , $nthrds);
	$config->set("ROOTPE_$comp" , $rootpe);
    }
}

#-------------------------------------------------------------------------------
sub _clean
{
    my ($name) = @_;
    $name =~ s/^\s+//; # strip any leading whitespace 
    $name =~ s/\s+$//; # strip any trailing whitespace
    return ($name);
}

1;



