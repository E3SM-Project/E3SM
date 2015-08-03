package ConfigPes;
my $pkg_nm = 'ConfigPes';

use strict;
use English;
use IO::File;
use XML::LibXML;
use Data::Dumper;
use ConfigCase;
use ConfigCESM;

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

    my($file_config, $pesize_opts, $primary_component, $opts_ref, $config) = @_;

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

	    my $cimeroot = $config->get('CIMEROOT');

	    # Determine pes setup files
	    my $pes_file = _setPesSetupFile($file_config, $cimeroot, $primary_component);

	    # Determine pes override file
	    my $override_file = _setPesOverrideFile($file_config, $cimeroot, $primary_component);

	    # Determine target grid
	    my $target_grid = ConfigCESM::setTargetGridMatch($primary_component, $config );

	    # Determine pe layout settings
	    _setPESsettings($pes_file, $override_file, $target_grid, $pesize_opts, $config);
	}
    }
}

#-----------------------------------------------------------------------------------------------
#                               Private routines
#-----------------------------------------------------------------------------------------------
sub _setPesSetupFile
{
    my ($file_config, $cimeroot, $primary_component) = @_;

    my $xml = XML::LibXML->new( no_blanks => 1)->parse_file($file_config);
    my @files = $xml->findnodes(".//entry[\@id=\"PES_SPEC_FILE\"]/values/value[\@component=\"$primary_component\"]");
    if (! @files) {
	die " ERROR ConfigPes::_setPesSetupFile: no pes specification file found for $primary_component \n";
    }
    my $file = $files[0]->textContent();
    $file =~ s/\$CIMEROOT/$cimeroot/;
    return ($file);
}

#-------------------------------------------------------------------------------
sub _setPesOverrideFile
{
    my ($file_config, $cimeroot, $primary_component) = @_;

    my $xml = XML::LibXML->new( no_blanks => 1)->parse_file($file_config);
    my @files = $xml->findnodes(".//entry[\@id=\"OVERRIDE_SPEC_FILE\"]/values/value[\@component=\"$primary_component\"]");
    if (! @files) {
	die " ERROR: no pes specification file found for $primary_component \n";
    }
    my $file = $files[0]->textContent();
    $file =~ s/\$CIMEROOT/$cimeroot/;
    return ($file);
}

#-------------------------------------------------------------------------------
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
    my ($pes_file, $override_file, $target_grid, $pesize_opts, $config) = @_; 

    my $mach = $config->get('MACH');
    my $compset_name = $config->get('COMPSET');

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
    # look for target_grid / target_machine
    # --------------------------------------

    # Before looking for grid match, determine if there is any grid translation that has to be done
    # to match on the grid name
    my $xml = XML::LibXML->new( no_blanks => 1)->parse_file($override_file);
    my @nodes = $xml->findnodes("//pes_grid_translation/grid");
    if (defined @nodes) {
	foreach my $node (@nodes) {
	    my $attr  = $node->getAttribute('name');
	    my $value = $node->textContent();
	    if ($target_grid =~ /$attr/) {
		$target_grid = $value;
	    }
	}
    }

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
    # Determine if there are any settings for target grid AND target_machine
    # If not, look for match for target_grid and 'any' machine
    # If not, look for match for 'any' grid and target_machine
    # If not, look for match for 'any' grid and 'any' machine
    
    my $mach_set = $mach; 
    my $grid_set = $target_grid;
    my @pes = $xml->findnodes(".//grid[contains(\@name,\"$grid_set\")]/mach[contains(\@name,\"$mach_set\")]/pes");
    if (! @pes) 
    {
	$mach_set = 'any';
	$grid_set = $target_grid;
	@pes = $xml->findnodes(".//grid[contains(\@name,\"$grid_set\")]/mach[\@name='any']/pes");
	if (! @pes) 
	{
	    $mach_set = $mach;
	    $grid_set = 'any';
	    @pes = $xml->findnodes(".//grid[\@name='any']/mach[contains(\@name,\"$mach_set\")]/pes");
	    if (! @pes) 
	    {
		$mach_set = 'any';
		$grid_set = 'any';
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
	    }
	}
    }
    print "ConfigPES: target grid match is $target_grid \n";
    print "ConfigPES: grid match is $grid_set \n";
    print "ConfigPES: machine match is $mach_set\n";

    if (($mach_set ne 'any') || ($grid_set ne 'any')) {
	# Now search the pes
	my @pes = $xml->findnodes(".//grid[contains(\@name,\"$grid_set\")]/mach[contains(\@name,\"$mach_set\")]");
	my $pe_select;
	my $compset_match;
	my $pesize_match;
      SEARCH: foreach my $pe (@pes) {
	  my @nodes = $pe->findnodes("./*[contains(\@pesize,\"$pesize_opts\") and not(contains(\@compset,'any'))]");
	  if (defined @nodes) {
	      foreach my $node (@nodes) {
		  my $compset_attr = $node->getAttribute('compset');
		  print "DEBUG: i am HERE2 with compset_attr of $compset_attr \n";
		  if ($compset_name =~ m/$compset_attr/) {
		      $pe_select = $node;
		      $pesize_match = $pesize_opts;
		      $compset_match = 'any';
		      last SEARCH;
		  }
	      }
	  }
	  if ((! defined @nodes) || (! defined $pe_select)) {
              @nodes = $pe->findnodes("./*[\@pesize='any' and not(contains(\@compset,'any'))]");
	      if (defined @nodes) {
		  foreach my $node (@nodes) {
		      my $compset_attr = $node->getAttribute('compset');
		      if ($compset_name =~ m/$compset_attr/) {
			  $pe_select = $node;
			  $pesize_match = 'any';
			  $compset_match = $compset_attr;
			  last SEARCH;
		      }
		  }
	      }
	      if ((! defined @nodes) || (! defined $pe_select)) {
		  @nodes = $pe->findnodes("./*[\@pesize=\"$pesize_opts\" and \@compset='any']");
		  if (defined @nodes) {
		      $pe_select = $nodes[0];
		  }
		  if ((! defined @nodes) || (! defined $pe_select)) {
		      @nodes = $pe->findnodes("./*[\@pesize='any' and \@compset='any']");
		      if (defined @nodes) {
			  $pesize_match = 'any';
			  $compset_match = 'any';
			  $pe_select = $nodes[0];
		      }
		  }
	      }
	  }
        }
	if (! defined $pe_select) {
	    die "ERROR: no pes match found in $pes_file \n";
	} else {
	    print "ConfigPES: compset_match is $compset_match\n"; 
	    print "ConfigPES: pesize match is $pesize_match \n"; 
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

    # --------------------------------------
    # override pes settings
    # --------------------------------------
    
    my $xml = XML::LibXML->new( no_blanks => 1)->parse_file($override_file);
    foreach my $node_grid ($xml->findnodes(".//pes_override/*")) 
    {
	my $gridname = $node_grid->getAttribute('name');
	if (($gridname eq 'any') || ($gridname =~ /$grid_set/ )) 
	{
	    foreach my $node_mach ($node_grid->findnodes("./mach")) 
	    {
		my $machname = $node_mach->getAttribute('name');
		foreach my $node_pes ($node_mach->findnodes("./pes")) 
		{
		    my $pesize  = $node_pes->getAttribute('pesize');
		    my $compset = $node_pes->getAttribute('compset');
		    if (($pesize eq 'any' && $compset eq 'any') ||
			($pesize eq 'any' && $compset_name =~ /$compset/)) 
		    {
			foreach my $pes ($node_pes->findnodes("./ntasks"), 
					 $node_pes->findnodes("./nthrds"), 
					 $node_pes->findnodes("./rootpe")) {
			    my @children = $pes ->childNodes();
			    foreach my $child (@children) {
				my $name  = uc $child->nodeName(); 
				my $value =    $child->textContent();
				$decomp{$name}  = $value;
			    }
			}
		    }
		}
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



