package ConfigPes;
my $pkg_nm = __PACKAGE__;

use strict;
use English;
use IO::File;
use XML::LibXML;
use Data::Dumper;
use ConfigCase;
use Log::Log4perl qw(get_logger);
my $logger;

BEGIN{
    $logger = get_logger();
}


#-----------------------------------------------------------------------------------------------
sub setPes {

    # Set the parameters for the pe layout.    

    my($pesize_opts, $primary_component, $opts_ref, $config) = @_;

    if (defined $$opts_ref{'pes_file'}) {

	# Reset the pes if a pes file is specified
	my $pes_file = $$opts_ref{'pes_file'};
	(-f "$pes_file")  or  $logger->logdie("** Cannot find pes_file \"$pes_file\" ***");
	$config->reset_setup("$pes_file");    

    } else {

	if ($pesize_opts =~ m!^([0-9]+)D?$!) {
	    my $ntasks = $1;
	    my $nthrds =  1;
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
    my ($pesize_opts, $ntasks, $nthrds, $rootpe, $config) = @_; 

    foreach my $comp ("ATM", "LND", "ICE", "OCN", "GLC", "WAV", "ROF", "CPL" ) {
	$config->set("NTASKS_" . "${comp}", $ntasks); 
	$config->set("NTHRDS_" . "${comp}", $nthrds); 
	$config->set("ROOTPE_" . "${comp}", $rootpe); 
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
	$config->set("NTASKS_" . "$comp", $ntasks); 
	$config->set("NTHRDS_" . "$comp", $nthrds); 
	$config->set("ROOTPE_" . "$comp", $rootpe); 
    }
    
    my $root;
    if ($pesize_opts =~ m!^([0-9]+)x([0-9]+)D$!) {
	$root = 0          ; $config->set('ROOTPE_ATM',$root);
	$root = 1 * $ntasks; $config->set('ROOTPE_LND',$root);
	$root = 2 * $ntasks; $config->set('ROOTPE_OCN',$root);
	$root = 3 * $ntasks; $config->set('ROOTPE_ICE',$root);
	$root = 4 * $ntasks; $config->set('ROOTPE_GLC',$root);
	$root = 5 * $ntasks; $config->set('ROOTPE_WAV',$root);
	$root = 6 * $ntasks; $config->set('ROOTPE_ROF',$root);
	$root = 7 * $ntasks; $config->set('ROOTPE_CPL',$root);
    }
}    

#-------------------------------------------------------------------------------
sub _setPESsettings
{
    # Read xml file and obtain NTASKS, NTHRDS, ROOTPE and NINST for each component
    my ($pesize_opts, $config) = @_; 

    my $mach		 = $config->get('MACH');
    my $grid_longname    = $config->get('GRID');
    my $compset_longname = $config->get('COMPSET');
    my $pes_file	 = $config->get('PES_SPEC_FILE');

    my %decomp;

    # --------------------------------------
    # look for model grid / model machine match
    # --------------------------------------

    # Parse the pes xml file
    my $xml = XML::LibXML->new( no_blanks => 1)->parse_file($pes_file);

    # First determine that the xml file has consistent attributes
    foreach my $grid ($xml->findnodes(".//grid")) {
	next if($grid->nodeType == XML_COMMENT_NODE);

	if (! defined $grid->getAttribute('name')) {
	    my $name = $grid->nodeName();
	    $logger->logdie("\n ERROR ConfigPes.pm: node <$name> does not have a required attribute of \"name\" in $pes_file \n\n");
	}
    }
    foreach my $mach ($xml->findnodes(".//mach")) {
	next if($mach->nodeType == XML_COMMENT_NODE);
	if (! defined $mach->getAttribute('name')) {
	    my $name = $mach->nodeName();
	    $logger->logdie("\n ERROR ConfigPes.pm: node <$name> does not have a required attribute of \"name\" in $pes_file \n\n");
	}
    }
    foreach my $pes ($xml->findnodes(".//pes")) {
	next if($pes->nodeType == XML_COMMENT_NODE);
	if (! defined $pes->getAttribute('pesize')) {
	    my $name = $pes->nodeName();
	    $logger->logdie("\n ERROR ConfigPes.pm: node <$name> does not have a required attribute of \"pesize\" in $pes_file \n\n");
	}
	if (! defined $pes->getAttribute('compset')) {
	    my $name = $pes->nodeName();
	    $logger->logdie("\n ERROR ConfigPes.pm: node <$name> does not have a required attribute of \"compset\" in $pes_file \n\n");
	}
    }

    # Set the match variables $mach_set and $grid_set
    # Determine possible matches in following order:
    # (1) model grid and model machine theh (2) model grid and 'any' machine
    
    my $mach_set;
    my $grid_set;

    # Determine setting for grid="any" AND mach="any"
    # If value is < 0 then it refers to pes per node

    my $pe_select;
    my $grid_match = 'any';
    my $mach_match = 'any';
    my $compset_match;
    my $pesize_match;

    # Find default which will be overwritten by any match below 
    foreach my $node ($xml->findnodes(".//grid[(\@name =\"any\")]/mach[\@name =\"any\"]/pes[\@pesize =\"any\" and \@compset =\"any\"]")) {
	next if($node->nodeType == XML_COMMENT_NODE);
	$pe_select = $node;
    }


    
        # Find nodes with grid= not 'any' and mach='any' - this override the default for $grid_match and $mach_match
    foreach my $node ($xml->findnodes(".//grid[not(\@name =\"any\")]/mach[\@name =\"any\"]")) {
	next if($node->nodeType == XML_COMMENT_NODE);
	my $grid_attr = $node->parentNode()->getAttribute('name');
	if ($grid_longname =~ m/$grid_attr/) {
	    $grid_match = $grid_attr;
	    $mach_match = 'any';
	    last;
	}
    }

    # Find nodes with grid = 'any' and mach = not 'any'
    # The first match found will overwrite the above grid_match and mach_match default
    foreach my $node ($xml->findnodes(".//grid[\@name =\"any\"]/mach[not(\@name =\"any\")]")) {
	next if($node->nodeType == XML_COMMENT_NODE);
	my $mach_attr = $node->getAttribute('name');
	if ($mach =~ m/$mach_attr/) {
	    $grid_match = 'any';
	    $mach_match = $mach_attr;
	    last;
	}
    }

    # Find nodes with grid = not 'any' and mach = not 'any'
    # The first match found will overwrite the above grid_match and mach_match default
    foreach my $node ($xml->findnodes(".//grid[not(\@name =\"any\")]/mach[not(\@name =\"any\")]")) {
	next if($node->nodeType == XML_COMMENT_NODE);
	my $grid_attr = $node->parentNode()->getAttribute('name');
	my $mach_attr = $node->getAttribute('name');
	if ($grid_longname =~ m/$grid_attr/ && $mach =~ m/$mach_attr/) {
	    $grid_match = $grid_attr;
	    $mach_match = $mach_attr;
	    last;
	}
    }

    # Now $grid_match and $mach_match are set
    # Determine the values of $compset_match and $pe_size_match
    my @nodes = $xml->findnodes(".//grid[\@name=\"$grid_match\"]/mach[\@name=\"$mach_match\"]/*");
    if ((! defined $compset_match) && (! defined $pesize_match)) {
	foreach my $node (@nodes) {
	    next if($node->nodeType == XML_COMMENT_NODE);
	    my $compset_attr = $node->getAttribute('compset');
	    my $pesize_attr  = $node->getAttribute('pesize');
	    if (($compset_longname =~ m/$compset_attr/) && ($pesize_opts =~ m/$pesize_attr/)) {
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
	    if (($compset_longname =~ m/$compset_attr/) && ($pesize_attr eq 'any')) {
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
	    if (($compset_attr eq 'any' ) && ($pesize_opts =~ m/$pesize_attr/)) {
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
	$logger->logdie("ERROR: no pes match found in $pes_file ");
    } 
    $logger->info("ConfigPES: grid          is $grid_longname ");
    $logger->info("ConfigPES: compset       is $compset_longname ");
    $logger->info("ConfigPES: grid match    is $grid_match ");
    $logger->info("ConfigPES: machine match is $mach_match");
    $logger->info("ConfigPES: compset_match is $compset_match"); 
    $logger->info("ConfigPES: pesize match  is $pesize_match "); 
    
    my $pes_per_node = $config->get('PES_PER_NODE');
    my @pes_ntasks = $pe_select->findnodes("./ntasks"); 
    my @pes_nthrds = $pe_select->findnodes("./nthrds"); 
    my @pes_rootpe = $pe_select->findnodes("./rootpe"); 
    my %decomp;
    
    foreach my $pes (@pes_ntasks, @pes_nthrds, @pes_rootpe) {
	next if($pes->nodeType == XML_COMMENT_NODE);
	my @children = $pes ->childNodes();
	foreach my $child (@children) {
	    next if($child->nodeType == XML_COMMENT_NODE);
	    my $name  = uc $child->nodeName(); 
	    my $value =    $child->textContent();
            if ($value =~ /-?\d/){
            if ($value < 0) {
		$value = -($value * $pes_per_node);
	    }
	    $decomp{$name}  = $value;
           }
	}
    }
    if(!defined $decomp{NTASKS_ATM} || ! defined $decomp{NTHRDS_ATM} || !defined $decomp{ROOTPE_ATM}){
	die("NO pelayout found");
    }


    
    foreach my $comp ("ATM", "LND", "ICE", "OCN", "GLC", "WAV", "ROF", "CPL") {
	if($comp ne "ATM"){
	    if(!defined $decomp{"NTASKS_$comp"}){
		$decomp{"NTASKS_$comp"} = $decomp{NTASKS_ATM};
		warn "No Match found for NTASKS_$comp, using $decomp{NTASKS_ATM}";
	    }
	    if(!defined $decomp{"NTHRDS_$comp"}){
		$decomp{"NTHRDS_$comp"} = $decomp{NTHRDS_ATM};
		warn "No Match found for NTHRDS_$comp, using $decomp{NTHRDS_ATM}";
	    }
	    if(!defined $decomp{"ROOTPE_$comp"}){
		$decomp{"ROOTPE_$comp"} = $decomp{ROOTPE_ATM};
		warn "No Match found for ROOTPE_$comp, using $decomp{ROOTPE_ATM}";
	    }
	}
	    

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



