package ConfigCompsetGrid;
my $pkg_nm = 'ConfigCompsetGrid';

use strict;
use English;
use Cwd qw( getcwd abs_path chdir);
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
  The perl module XML::LibXML is needed for XML parsing in the CIME script system.
  Please contact your local systems administrators or IT staff and have them install it for
  you, or install the module locally.  

END
    print "$warning\n";
    exit(1);
}

#-----------------------------------------------------------------------------------------------
sub setCompsetGrid {

    # Set the parameters for the specified compset and grid.  
    # The parameters are read from an input file, and if no compset/grid matches are
    # found then issue error message.
    # This routine uses the configuration defined at the package level ($cfg_ref).
    # input arguments

    my ($print_flag, $cfg_ref) = @_;
    
    my $cimeroot         = $cfg_ref->get('CASEROOT');
    my $caseroot         = $cfg_ref->get('CIMEROOT');		
    my $compset_longname = $cfg_ref->get('COMPSET');		
    my $grid_longname    = $cfg_ref->get('GRID');			
    my $grids_file       = $cfg_ref->get('GRIDS_SPEC_FILE');		
    
    (-f "$grids_file") or  die "Cannot find supported model grids file $grids_file ";

    # Note that DRV must be first in the list below
    my @setup_comp_files;
    foreach my $name ('CONFIG_DRV_FILE', 'CONFIG_ATM_FILE', 'CONFIG_ICE_FILE', 'CONFIG_GLC_FILE', 
		      'CONFIG_LND_FILE', 'CONFIG_ROF_FILE', 'CONFIG_OCN_FILE', 'CONFIG_WAV_FILE') {
	my $file = $cfg_ref->get($name);	
	push (@setup_comp_files, $file);
    }

    # ====================================================================
    # Determine general compset variables and
    # check that compset is supported for target grid
    # ====================================================================

    my %newxml;    # output in following call
    _setCompsetGeneralVars(\%newxml, $cfg_ref); 

    # ====================================================================
    # determine compgrid hash
    # ====================================================================

    # Assumes the following order for specifying a grid name
    #  a%aname_l%lname_oi%oiname_r%rname_m%mname_g%gname_w%wname

    my %compgrid;
    $grid_longname =~ /(a%)(.+)(_l%)/ ; $compgrid{'atm'}  = $2;
    $grid_longname =~ /(l%)(.+)(_oi%)/; $compgrid{'lnd'}  = $2; $compgrid{'glc'} = $2;
    $grid_longname =~ /(oi%)(.+)(_r%)/; $compgrid{'ocn'}  = $2; $compgrid{'ice'} = $2;
    $grid_longname =~ /(r%)(.+)(_m%)/ ; $compgrid{'rof'}  = $2; 
    $grid_longname =~ /(m%)(.+)(_g%)/ ; $compgrid{'mask'} = $2; 
    $grid_longname =~ /(g%)(.+)(_w%)/ ; $compgrid{'cism'} = $2; 
    $grid_longname =~ /(w%)(.+)$/     ; $compgrid{'wav'}  = $2; 

    # ========================================================
    # set grid component domains and related variables
    # ========================================================

    _setGridDomain($grids_file, \%compgrid, $cfg_ref);

    # ====================================================================
    # set grid mapping variables
    # ====================================================================

    _setGridMaps($grids_file, \%compgrid, $cfg_ref);

    # ====================================================================
    # Determine compset component configurations 
    # Loop over the component files (i.e. definitions_component.xml files)
    # determined by the compsets of the primary component
    # ====================================================================


    my $desc_comp = "";
    foreach my $setup_comp_file (@setup_comp_files) {
	_setComponent($setup_comp_file, \%compgrid, \%newxml, \$desc_comp, $cfg_ref);
    }

    # ====================================================================
    # Special case - if land and river grids are different AND there 
    # is no mapping file between  land and river then set the river mode 
    # to null if the river component is rtm
    # ====================================================================

    my $rof_comp = $cfg_ref->get('COMP_ROF');
    if ($rof_comp eq 'rtm') {
	my $lnd_grid    = $cfg_ref->get('LND_GRID');
	my $rof_grid    = $cfg_ref->get('ROF_GRID');
	my $map_lnd2rof = $cfg_ref->get('LND2ROF_FMAPNAME');

	if (($lnd_grid ne $rof_grid) && ($map_lnd2rof eq 'idmap')) {
	    print "No lnd2rof_fmapname exists - RTM mode set to null \n";
	    $cfg_ref->set('RTM_MODE', 'NULL');     
	    $newxml{"RTM_MODE"} = 'NULL'; 
	}
    }

    # ========================================================
    # Print compset/grid info
    # ========================================================

    _printGridCompsetInfo($grids_file, $desc_comp, \%newxml, $print_flag, $cfg_ref );

    if ($print_flag >= 2) { 
	my $compset = $cfg_ref->get('COMPSET');
	print "Compset specifier: $compset.\n"; 
	print "Grid is valid for this compset. \n"; 
    }
}

#-----------------------------------------------------------------------------------------------
sub getCompsetLongname
{
    # Determine compset longname, alias and support level
    my ($compsets_file, $compset_input,  $user_compset) = @_;

    my $compset_longname;
    my $compset_aliasname;
    my $support_level;

    my $xml_compsets = XML::LibXML->new( no_blanks => 1)->parse_file($compsets_file);

    if ($user_compset) {
	$compset_aliasname = ' ';
    } else {
	my $found = 0;
	my @nodes;
	if (! $found) {
	    @nodes = $xml_compsets->findnodes(".//COMPSET[lname=\"$compset_input\"]");
	    if (@nodes) {
		if ($#nodes > 1) {
		    die "ERROR: more than one node was found with compset name $compset_input \n";
		}
		$found = 1;
	    }
	}
	if (! $found) {
	    @nodes = $xml_compsets->findnodes(".//COMPSET[alias=\"$compset_input\"]");
	    if (@nodes) {
		if ($#nodes > 1) {
		    die "ERROR: more than one node was found with compset name $compset_input \n";
		}
		$found = 1;
	    }
	}
	unless ($found) { 
	    print "ERROR getCompsetLongname: no match for compset $compset_input \n";
	    print "  to see supported compsets issue \n";
	    print "  ./manage_case -list compsets -compsets_setby <target component name>\n";
	    die "setCompset: exiting\n"; 
	}
	my @lname_nodes = $nodes[0]->findnodes("./lname");
	my @alias_nodes = $nodes[0]->findnodes("./alias");
	my @support_nodes = $nodes[0]->findnodes("./support_level");

	my $lname = $lname_nodes[0]->textContent();
	my $alias = $alias_nodes[0]->textContent();
	my $support;
	if (@support_nodes) {
	    $support = $support_nodes[0]->textContent();
	}	
	if ( $support ) {$support_level .= "Compset ($alias): $support\n"; }

	return ($lname, $alias, $support);
    }
}

#-------------------------------------------------------------------------------
sub getGridLongname
{
    my ($cfg_ref, $grid_input) = @_;

    my ($grid_longname, $grid_shortname, $grid_aliasname);
    my $compset_match;

    my $grids_file = $cfg_ref->get('GRIDS_SPEC_FILE');

    my $xml = XML::LibXML->new( no_blanks => 1)->parse_file($grids_file);
    my @nodes_alias = $xml->findnodes(".//GRID[alias=\"$grid_input\"]");
    my @nodes_lname = $xml->findnodes(".//GRID[lname=\"$grid_input\"]");
    my @nodes_sname = $xml->findnodes(".//GRID[sname=\"$grid_input\"]");

    my $grid_node;
    if (@nodes_alias) {
	$grid_node = $nodes_alias[0];
    } elsif (@nodes_lname) {
	$grid_node = $nodes_lname[0];
    } elsif (@nodes_sname) {
	$grid_node = $nodes_sname[0];
    } else { 
	die " ERROR: no supported grid match for target grid $grid_input \n";
    }
    
    # set the compset grid alias and longname
    foreach my $node ($grid_node->childNodes()) {
	my $name = $node->nodeName();
	my $value = $node->textContent();
	if ($name eq 'lname') {$grid_longname   = $node->textContent();}
	if ($name eq 'sname') {$grid_shortname  = $node->textContent();}
	if ($name eq 'alias') {$grid_aliasname  = $node->textContent();}
    }	

    return ($grid_longname, $grid_shortname, $grid_aliasname);
}

#-----------------------------------------------------------------------------------------------
#                               Private routines
#-----------------------------------------------------------------------------------------------
sub _setGridDomain
{
    my ($grids_file, $compgrid_ref, $cfg_ref) = @_;
    
    my $parser = XML::LibXML->new( no_blanks => 1);
    my $xml = $parser->parse_file($grids_file);
  
    # set the component grid names
    foreach my $key (keys %$compgrid_ref) {
	my $comp = uc $key;
	my $grid = $$compgrid_ref{$key};
	my $mask = $$compgrid_ref{'mask'};
	
	if ($cfg_ref->is_valid_name("${comp}_GRID")) {
	    my @nodes = $xml->findnodes(".//domain[\@name=\"$grid\"]");
	    if (@nodes) {$cfg_ref->set("${comp}_GRID", $nodes[0]->getAttribute('name'));}
	}	
	
	if ($cfg_ref->is_valid_name("${comp}_NX")) {
	    my @nodes = $xml->findnodes(".//domain[\@name=\"$grid\"]/nx");
	    if (@nodes) {$cfg_ref->set("${comp}_NX", $nodes[0]->textContent());}
	}
	
	if ($cfg_ref->is_valid_name("${comp}_NY")) {
	    my @nodes = $xml->findnodes(".//domain[\@name=\"$grid\"]/ny");
	    if (@nodes) { $cfg_ref->set("${comp}_NY", $nodes[0]->textContent());}
	}
	
	# Note domain files are only specified for ATM, LND, ICE and OCN components 
	if ($comp eq 'ATM' || $comp eq 'LND' || $comp eq 'ICE' || $comp eq 'OCN') {
	    my @nodes;
	    my $var = "$comp" . "_DOMAIN_FILE";
	    if ($cfg_ref->is_valid_name($var)) {
		if ($comp eq 'ATM' || $comp eq 'LND') {	    
		    @nodes = $xml->findnodes(".//domain[\@name=\"$grid\"]/file[\@lnd_mask=\"$mask\"]");
		} else {
		    @nodes = $xml->findnodes(".//domain[\@name=\"$grid\"]/file[\@ocn_mask=\"$mask\"]");
		}
		if (@nodes) {$cfg_ref->set($var,$nodes[0]->textContent());}
	    }
	    my $var = "$comp" . "_DOMAIN_PATH";
	    if ($cfg_ref->is_valid_name($var)) {
		if ($comp eq 'ATM' || $comp eq 'LND') {	    
		    @nodes = $xml->findnodes(".//domain[\@name=\"$grid\"]/path[\@lnd_mask=\"$mask\"]");
		} else {
		    @nodes = $xml->findnodes(".//domain[\@name=\"$grid\"]/path[\@ocn_mask=\"$mask\"]");
		}
		if (@nodes) {$cfg_ref->set($var,$nodes[0]->textContent());}
	    }
	}
	
    }
    
    if ($$compgrid_ref{'cism'} ne 'null') {
	$cfg_ref->set('CISM_GRID',$$compgrid_ref{'cism'});
    }
}

#-------------------------------------------------------------------------------
sub _setGridMaps
{
    # set grid mapping variables
    my ($grids_file, $compgrid_ref, $cfg_ref) = @_;
    
    my $parser = XML::LibXML->new( no_blanks => 1);
    my $xml = $parser->parse_file($grids_file);

    # first set defaults
    my @nodes = $xml->findnodes(".//grid_mappings/defaults/*");
    foreach my $node (@nodes) {
	my $name  = $node->nodeName();
	my $value = $node->textContent();
	$cfg_ref->set($name, $value);
    }

    # overwrite the default values if appropriate
    my $atm_grid = $$compgrid_ref{'atm'};
    my $lnd_grid = $$compgrid_ref{'lnd'};
    my $rof_grid = $$compgrid_ref{'rof'};
    my $ocn_grid = $$compgrid_ref{'ocn'};
    my $ice_grid = $$compgrid_ref{'ice'};
    my $glc_grid = $$compgrid_ref{'glc'};
    my $wav_grid = $$compgrid_ref{'wav'};

    @nodes = $xml->findnodes(".//gridmap[\@atm_grid=\"$atm_grid\" and \@ocn_grid=\"$ocn_grid\"]/ATM2OCN_FMAPNAME");
    if (@nodes) {$cfg_ref->set('ATM2OCN_FMAPNAME',$nodes[0]->textContent())}

    @nodes = $xml->findnodes(".//gridmap[\@atm_grid=\"$atm_grid\" and \@ocn_grid=\"$ocn_grid\"]/ATM2OCN_SMAPNAME");
    if (@nodes) {$cfg_ref->set('ATM2OCN_SMAPNAME',$nodes[0]->textContent())}

    @nodes = $xml->findnodes(".//gridmap[\@atm_grid=\"$atm_grid\" and \@ocn_grid=\"$ocn_grid\"]/ATM2OCN_VMAPNAME");
    if (@nodes) {$cfg_ref->set('ATM2OCN_VMAPNAME',$nodes[0]->textContent())}

    @nodes = $xml->findnodes(".//gridmap[\@atm_grid=\"$atm_grid\" and \@ocn_grid=\"$ocn_grid\"]/OCN2ATM_FMAPNAME");
    if (@nodes) {$cfg_ref->set('OCN2ATM_FMAPNAME',$nodes[0]->textContent())}

    @nodes = $xml->findnodes(".//gridmap[\@atm_grid=\"$atm_grid\" and \@ocn_grid=\"$ocn_grid\"]/OCN2ATM_SMAPNAME");
    if (@nodes) {$cfg_ref->set('OCN2ATM_SMAPNAME',$nodes[0]->textContent())}

    @nodes = $xml->findnodes(".//gridmap[\@atm_grid=\"$atm_grid\" and \@lnd_grid=\"$lnd_grid\"]/ATM2LND_SMAPNAME");
    if (@nodes) {$cfg_ref->set('ATM2LND_SMAPNAME',$nodes[0]->textContent())}

    @nodes = $xml->findnodes(".//gridmap[\@atm_grid=\"$atm_grid\" and \@lnd_grid=\"$lnd_grid\"]/ATM2LND_FMAPNAME");
    if (@nodes) {$cfg_ref->set('ATM2LND_FMAPNAME',$nodes[0]->textContent())}

    @nodes = $xml->findnodes(".//gridmap[\@atm_grid=\"$atm_grid\" and \@lnd_grid=\"$lnd_grid\"]/LND2ATM_SMAPNAME");
    if (@nodes) {$cfg_ref->set('ATM2LND_SMAPNAME',$nodes[0]->textContent())}

    @nodes = $xml->findnodes(".//gridmap[\@atm_grid=\"$atm_grid\" and \@lnd_grid=\"$lnd_grid\"]/LND2ATM_FMAPNAME");
    if (@nodes) {$cfg_ref->set('LND2ATM_FMAPNAME',$nodes[0]->textContent())}

    @nodes = $xml->findnodes(".//gridmap[\@lnd_grid=\"$lnd_grid\" and \@rof_grid=\"$rof_grid\"]/LND2ROF_FMAPNAME");
    if (@nodes) {$cfg_ref->set('LND2ROF_FMAPNAME',$nodes[0]->textContent())}

    @nodes = $xml->findnodes(".//gridmap[\@lnd_grid=\"$lnd_grid\" and \@rof_grid=\"$rof_grid\"]/ROF2LND_FMAPNAME");
    if (@nodes) {$cfg_ref->set('ROF2LND_FMAPNAME',$nodes[0]->textContent())}

    @nodes = $xml->findnodes(".//gridmap[\@ocn_grid=\"$ocn_grid\" and \@rof_grid=\"$rof_grid\"]/ROF2OCN_FMAPNAME");
    if (@nodes) {$cfg_ref->set('ROF2OCN_FMAPNAME',$nodes[0]->textContent())}

    @nodes = $xml->findnodes(".//gridmap[\@ocn_grid=\"$ocn_grid\" and \@rof_grid=\"$rof_grid\"]/ROF2OCN_RMAPNAME");
    if (@nodes) {$cfg_ref->set('ROF2OCN_RMAPNAME',$nodes[0]->textContent())}

    @nodes = $xml->findnodes(".//gridmap[\@ocn_grid=\"$ocn_grid\" and \@rof_grid=\"$rof_grid\" and \@lnd_grid=\"$lnd_grid\"]/XROF_FLOOD_MODE");
    if (@nodes) {$cfg_ref->set('XROF_FLOOD_MODE',$nodes[0]->textContent())}
}

#-------------------------------------------------------------------------------
sub _setCompsetGeneralVars
{
    # Determine general compset variables
    my ($newxml_ref, $cfg_ref) = @_;

    my $compset_longname = $cfg_ref->get('COMPSET');
    my $grid_longname    = $cfg_ref->get('GRID');
    my $compsets_file    = $cfg_ref->get('COMPSETS_SPEC_FILE');	

    my $parser = XML::LibXML->new( no_blanks => 1);
    my $xml_compset = $parser->parse_file($compsets_file);

    # first check that compset is supported for target grid
    my @nodes = $xml_compset->findnodes(".//COMPSET[lname=\"$compset_longname\"]");
    if (! @nodes) {die " ERROR: $compset_longname not supported \n";}

    my $compset_grid = $nodes[0]->getAttribute('grid');
    if ($compset_grid) {
	$compset_grid = _clean($compset_grid);
	if ($grid_longname !~ /$compset_grid/) {
	    # TODO - this is CESM specific and should be removed
	    if ($grid_longname !~ /CLM_USRDAT/) {
		die "ERROR: $compset_longname \n is not supported for \n $grid_longname \n";
	    }
	}
    }

    foreach my $node ($xml_compset->findnodes(".//compset_variables/*")) {
	my $compset_match = $node->getAttribute('compset');
	my $grid_match    = $node->getAttribute('grid');

	my $set_cfg_ref;
	my $debug = $node->nodeName();
	if ($compset_match && $grid_match) {
	    if ($compset_longname =~ /$compset_match/ && $grid_longname =~ /$grid_match/) {	    
		$set_cfg_ref = 'yes';
	    }
	} elsif ($compset_match) {
	    if ($compset_longname =~ /$compset_match/) {
		$set_cfg_ref = 'yes';
	    }
	} elsif ($grid_match) {
	    if ($grid_longname =~ /$grid_match/) {
		$set_cfg_ref = 'yes';
	    }
	} else {
	    $set_cfg_ref = 'yes';
	}
	if ($set_cfg_ref) {
	    my $id      =  $node->nodeName();
	    my $new_val =  $node->textContent();
	    if ($cfg_ref->is_valid_name($id)) {
		$$newxml_ref{$id} = $new_val; 
		$cfg_ref->set($id, $new_val);
		
	    } else {
		die "ERROR: $id is not a valid name in $compsets_file \n";
	    }
	}
    }
}

#-------------------------------------------------------------------------------
sub _setComponent {

    # set variable value for compset component
    # first determine if there is a additive or merge values for the variable value 
    # for additive attributes on variables, new values are added to the end
    # of the string UNLESS there is already a match.  As an example:						
    # If initially CAM_CONFIG_OPTS is "-phys cam4 -chem none" and a new value
    # of CAM_CONFIG_OPTS is matched that has CAM_CONFIG_OPTS set to "-phys cam5"			
    # then the final CAM_CONFIG_OPTS  is "-phys cam5 -chem none"
    # All other xml variables take on the value of the last match. As an example
    # Assume there are two matches for xml variable RUN_STARTDATE, date1 and date2,
    # with date1 appearing first. The final value of RUN_STARTDATE will be date2.			

    my ($component_file, $compgrid_ref, $newxml_ref, $desc_comp_ref, $cfg_ref) = @_;

    my $grid_longname	   = $cfg_ref->get('GRID');
    my $compset_longname   = $cfg_ref->get('COMPSET');

    my $xml = XML::LibXML->new( no_blanks => 1)->parse_file($component_file);

    # Parse the component definition file and add all the variable definitions to $cfg_ref
    $cfg_ref->add_config_variables($component_file);

    # Now set the actual values of the variable
    my $new_val;
    my @variable_nodes = $xml->findnodes('.//entry');
    foreach my $variable_node (@variable_nodes) 
    {
	my $name     = $variable_node->getAttribute('id');
	my @additive = $variable_node->findnodes("./values[\@additive=\"yes\"]");
	my @merge    = $variable_node->findnodes("./values[\@merge=\"yes\"]");
	foreach my $value_node ($variable_node->findnodes('./values/value')) {
	    if ($value_node->nodeType() == XML_COMMENT_NODE) {
		# do nothing
	    } else {
		my $compset_match = $value_node->getAttribute('compset');
		my $grid_match    = $value_node->getAttribute('grid');
		my $match = 'no'; 
		if ($compset_match && $grid_match) {
		    if  (($compset_longname =~ /$compset_match/) && ($grid_longname =~ /$grid_match/)) {
			$match = 'yes';
		    }
		} elsif ($compset_match) {
		    if  ($compset_longname =~ /$compset_match/) {
			$match = 'yes';
		    }
		} elsif ($grid_match) {
		    if  ($grid_longname =~ /$grid_match/) {
			$match = 'yes';
		    }
		}
		
		if ($match eq 'yes') {
		    if ($cfg_ref->is_valid_name($name)) {
			# do nothing
		    } else {
			die "ERROR: $name is not a valid name in $component_file \n";
		    }

		    my $input_val = $value_node->textContent();
		    my $current_val = $cfg_ref->get($name);
		    if (@additive) {
			$new_val = _setComponentAdditiveNewVal($current_val, $input_val);
		    } elsif (@merge) {
			$new_val = _setComponentMergeNewVal($current_val, $input_val);
		    } else {
			$new_val = $input_val;
		    }
		    $$newxml_ref{$name} = $new_val; 
		    $cfg_ref->set($name, $new_val);
		}
	    }
	}
    }

    # set variable description
    my @desc_nodes = $xml->findnodes('.//description/desc');
    foreach my $desc_node (@desc_nodes) 
    {
	my $compset_match = $desc_node->getAttribute('compset');
	if  ($compset_longname =~ /$compset_match/) {
	    my $value = $desc_node->textContent();
	    if ($$desc_comp_ref =~ m/$value/) {
		# do nothing
	    } else {
		if($value =~ /\n/){
		    $$desc_comp_ref = $$desc_comp_ref . "\n$value\n";
		} else {
		    $$desc_comp_ref = $$desc_comp_ref . " $value";
		}
	    }
	}
    }

   
    # Special case for CLMUSRDAT in compset name
    if ($compset_longname =~ /(.+CLMUSRDAT%)(.*)/){
     	$cfg_ref->set('CLM_USRDAT_NAME',$2);
     	$$newxml_ref{"CLM_USR_DATNAME"} = $2;
    } 
    
    # Determine run_refcase and run_refdate
    my $run_refcase = _clean($cfg_ref->get('RUN_REFCASE'));
    my $run_refdate = _clean($cfg_ref->get('RUN_REFDATE'));
    if ($run_refcase ne 'case.std') {
     	$cfg_ref->set('RUN_TYPE','hybrid');
     	$cfg_ref->set('GET_REFCASE','TRUE');
     	$$newxml_ref{"RUN_TYPE"}    = 'hybrid'; 
     	$$newxml_ref{"GET_REFCASE"} = 'TRUE';
     }
}

#-------------------------------------------------------------------------------
sub _setComponentAdditiveNewVal 
{
    # input arguments
    my ($current_val, $input_val) = @_;

    my $new_val;
    if ($current_val !~ m/^\s*$/)  { 
	$new_val = $current_val;
	
	# Split into separate options, separated by '-'.
	# The regex is used to ensure '-' is only noticed if it is
	# either the first character or follows a space.
	# Note that the '-' will be stripped off.
	my @nameopts = split(/(?:^|\s)-/, $input_val);
	my @curopts  = split(/(?:^|\s)-/, $current_val);
	
	# First item in each array will be space or empty string, so
	# remove it with shift.
	shift @nameopts;
	shift @curopts;
	
	# Iterate through new options
	foreach my $nameopt (@nameopts) {
	    
	    # Grab option name.
	    my ($optname) = $nameopt =~ m/^(\w+)\s/;
	    my $name_found = 0;
	    
	    # Check current options for values to replace.
	    foreach my $curopt (@curopts) {
		if ($curopt =~ m/^$optname\s/) {
		    $name_found = 1;
		    # Substitute, adding one space just in case.
		    $new_val =~ s/$curopt/$nameopt /;
		}
	    }
	    # If the new option was not found in existing options, append it.
	    if ( ! $name_found) {
		$new_val = "$new_val -$nameopt";
	    }
	}
	# Get rid of extra spaces.
	$new_val =~ s/\s+/ /g; # spaces in middle
	$new_val =~ s/\s*$//; # spaces at end
    } else {
	$new_val = $input_val;
    }
    return $new_val
}	

#-------------------------------------------------------------------------------
sub _setComponentMergeNewVal 
{
    # input arguments
    my ($current_val, $input_val) = @_;
    
    my $new_val;
    if ($input_val =~ /$current_val/) {
	$new_val = $input_val;
    } else {
	$new_val = "$current_val" . " $input_val"; 
    }
    return $new_val;
}	

#-------------------------------------------------------------------------------
sub _printGridCompsetInfo 
{
    my ($grids_file, $desc_comp, $newxml, $print_flag, $cfg_ref) = @_;

    my $caseroot	   = $cfg_ref->get('CASEROOT');
    my $grid_longname	   = $cfg_ref->get('GRID');
    my $compset_longname   = $cfg_ref->get('COMPSET');

    my @grids;
    my $desc_grid = "";

    my $atm_nx   = $cfg_ref->get('ATM_NX'); 
    my $atm_ny   = $cfg_ref->get('ATM_NY'); 
    my $atm_grid = $cfg_ref->get('ATM_GRID');
    push (@grids, $atm_grid);

    my $lnd_nx   = $cfg_ref->get('LND_NX'); 
    my $lnd_ny   = $cfg_ref->get('LND_NY'); 
    my $lnd_grid = $cfg_ref->get('LND_GRID');
    push (@grids, $lnd_grid);

    my $ice_nx   = $cfg_ref->get('ICE_NX'); 
    my $ice_ny   = $cfg_ref->get('ICE_NY'); 
    my $ice_grid = $cfg_ref->get('ICE_GRID');
    push (@grids, $ice_grid);

    my $ocn_nx   = $cfg_ref->get('OCN_NX'); 
    my $ocn_ny   = $cfg_ref->get('OCN_NY'); 
    my $ocn_grid = $cfg_ref->get('OCN_GRID');
    push (@grids, $ocn_grid);

    my $rof_nx   = $cfg_ref->get('ROF_NX'); 
    my $rof_ny   = $cfg_ref->get('ROF_NY'); 
    my $rof_grid = $cfg_ref->get('ROF_GRID');
    push (@grids, $rof_grid);

    my $glc_nx   = $cfg_ref->get('GLC_NX'); 
    my $glc_ny   = $cfg_ref->get('GLC_NY'); 
    my $glc_grid = $cfg_ref->get('GLC_GRID');
    push (@grids, $glc_grid);

    my $wav_nx    = $cfg_ref->get('WAV_NX'); 
    my $wav_ny    = $cfg_ref->get('WAV_NY'); 
    my $wav_grid  = $cfg_ref->get('WAV_GRID');
    push (@grids, $wav_grid);

    my $cism_grid = $cfg_ref->get('CISM_GRID');
    push (@grids, $cism_grid);

    my $mask_grid = $cfg_ref->get('MASK_GRID');
    push (@grids, $mask_grid);

    my $parser = XML::LibXML->new( no_blanks => 1);
    my $xml_domain = $parser->parse_file($grids_file);

    my $support_level; 
    foreach my $grid (@grids)
    {
	my @nodes = $xml_domain->findnodes(".//domain[\@name=\"grid\"]/support_level");
	if (@nodes) {
	    my $value = @nodes[0]->textContent();
	    my $name  = @nodes[0]->nodeName();
	    $support_level .= "Grid ($name): $value\n";
	}
	my @nodes = $xml_domain->findnodes(".//domain[\@name=\"$grid\"]/desc");
	if (@nodes) {
	    my $value = @nodes[0]->textContent();
	    if ($desc_grid !~ m/$value/) {
		$desc_grid = $desc_grid . " $value" ;
	    }
	}
    }

    # ========================================================
    # Generate compset description - TODO - fill this in
    # ========================================================

    #my @compset_parts = split /[-_%]/, $compset_longname;
    # 	    if ($compset_longname =~ /$attr{'compset'}/)  {
    # 		if ($desc_comp =~ m/$val/) {
    # 		    # do nothing
    # 		} else {
    # 		    if($val =~ /\n/){
    # 			$desc_comp = $desc_comp . "\n$val\n";
    # 		    }else{
    # 			$desc_comp = $desc_comp . " $val";
    # 		    }
    # 		}
    # 	    }
    # 	    # Make sure that there is a description match for each part of the compset_name
    # 	    # this assures that there are no typos in the name.
    # 	    my $cnt=0; my $atrcnt=0;
    # 	    my @compset_tmp = @compset_parts;
    # 	    foreach my $part (@compset_tmp){
    # 		if(grep /$part/, $attr{'compset'}){
    # 		    splice(@compset_parts, $cnt, 1);
    # 		}else{
    # 		    $cnt++;
    # 		}
    # 	    }  
    # 	}
    # }
    # die "Could not find definition for part @compset_parts of compset name" if($#compset_parts >= 0) ;


    if ($print_flag) {
	my $compset_longname = $cfg_ref->get('COMPSET');

	my $fh_case = new IO::File;
	$fh_case->open(">>$caseroot/README.case") or die "can't open file: README.case\n";

	my $fh_stdout = *STDOUT;
	my @file_handles = ($fh_stdout, $fh_case);

	foreach my $fhandle (@file_handles) {
	    print $fhandle "Component set: longname \n";
	    print $fhandle "  $compset_longname \n";
	    print $fhandle "Component set Description: \n";
	    print $fhandle " $desc_comp \n";
	    print $fhandle "Grid: \n";
	    print $fhandle "  $grid_longname \n";
	    print $fhandle "  ATM_GRID = $atm_grid  NX_ATM=$atm_nx NY_ATM=$atm_ny \n";
	    print $fhandle "  LND_GRID = $lnd_grid  NX_LND=$lnd_nx NX_LND=$lnd_ny \n";
	    print $fhandle "  ICE_GRID = $ice_grid  NX_ICE=$ice_nx NX_ICE=$ice_ny \n";
	    print $fhandle "  OCN_GRID = $ocn_grid  NX_OCN=$ocn_nx NX_OCN=$ocn_ny \n";
	    print $fhandle "  ROF_GRID = $rof_grid  NX_ROF=$rof_nx NX_ROF=$rof_ny \n";
	    print $fhandle "  GLC_GRID = $glc_grid  NX_GLC=$glc_nx NX_GLC=$glc_ny \n";
	    print $fhandle "  WAV_GRID = $wav_grid  NX_WAV=$wav_nx NX_WAV=$wav_ny \n";
	    # TODO - the following is CESM specific
	    if ($compset_longname =~ /CISM/) {
		my $cism_grid = $cfg_ref->get('CISM_GRID');
		print $fhandle "  CISM_GRID = $cism_grid \n";
	    }
	    print $fhandle "Grid Description: \n";
	    print $fhandle " $desc_grid \n";
	    print $fhandle "Non-Default Options: \n";
	    my @ids = keys %$newxml;
	    foreach my $id (sort @ids) {
		my $value = $$newxml{$id};
		print $fhandle     "  $id: $value \n";
	    } 
	    print $fhandle "\n";
	}
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
