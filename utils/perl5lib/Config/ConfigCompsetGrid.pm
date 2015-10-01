package ConfigCompsetGrid;
my $pkg_nm = 'ConfigCompsetGrid';

use strict;
use English;
use Cwd qw( getcwd abs_path chdir);
use IO::File;
use XML::LibXML;
use Data::Dumper;
use ConfigCase;
use Log::Log4perl qw(get_logger);
my $logger;

BEGIN{
    $logger = get_logger();
}

# Global variable
my %compgrid;
my %newxml;    
my $desc_comp = "";

#-----------------------------------------------------------------------------------------------
sub getCompsetLongname
{
    # Determine compset longname, alias and support level
    my ($input_file, $input_compset, $config) = @_;

    $input_compset =~ s/^\s+//; # strip any leading whitespace 
    $input_compset =~ s/\s+$//; # strip any trailing whitespace

    # Note - the value of the config variable 'COMPSETS_SPEC_FILE' gives the full pathname of the
    # file containing the possible out of the box compsets that can be used by create_newcase

    # TODO: add logic for determining support level
    # TODO: Add support level to $config rather than an argument

    my $support_level;
    my $pes_setby;

    my $cimeroot = $config->get('CIMEROOT');
    my $srcroot  = $config->get('SRCROOT');
    my $model    = $config->get('MODEL');
    my $compset_longname;
    my $compset_aliasname;

    # First determine primary component (for now this is only CESM specific)
    # Each primary component is responsible for defining the compsets that turn of the
    # appropriate feedbacks for development of that component

    my $compsets_file;
    my $compset_longname;

    my $pes_setby;
    my $xml1 = XML::LibXML->new( no_blanks => 1)->parse_file("$input_file");

    # Loop through all of the files listed in COMPSETS_SPEC_FILE and find the file
    # that has a match for either the alias or the longname in that order
    my @nodes = $xml1->findnodes(".//entry[\@id=\"COMPSETS_SPEC_FILE\"]/values/value");
    foreach my $node_file (@nodes) {
	my $file = $node_file->textContent();
	$file =~ s/\$CIMEROOT/$cimeroot/;
	$file =~ s/\$SRCROOT/$srcroot/;
	$file =~ s/\$MODEL/$model/;
	if (! -f $file ) {next;}
	my $xml2 = XML::LibXML->new( no_blanks => 1)->parse_file("$file");

	# First determine if there is a match for the alias - if so stop
	my @alias_nodes = $xml2->findnodes(".//compset[alias=\"$input_compset\"]");
	if (@alias_nodes) {
	    if ($#alias_nodes > 0) {
		die "ERROR create_newcase: more than one match for alias element in file $file \n";
	    } else {
		my @name_nodes = $alias_nodes[0]->childNodes();
		foreach my $name_node (@name_nodes) {
		    my $debug = $name_node->nodeName();
		    if ($name_node->nodeName() eq 'lname') {
			$compset_longname = $name_node->textContent();
		    }		    
		}
	    }
	    $pes_setby = $node_file->getAttribute('component');
	    $compsets_file = $file;
	    $config->set('COMPSETS_SPEC_FILE', $compsets_file);
	    $config->set('COMPSET', "$compset_longname");
	    last;
	} 

	# If no alias match - then determine if there is a match for the longname
	my @lname_nodes = $xml2->findnodes(".//compset[lname=\"$input_compset\"]");
	if (@lname_nodes) {
	    if ($#lname_nodes > 0) {
		die "ERROR create_newcase: more than one match for lname element in file $file \n";
	    } else {
		my @name_nodes = $lname_nodes[0]->childNodes();
		foreach my $name_node (@name_nodes) {
		    my $debug = $name_node->nodeName();
		    if ($name_node->nodeName() eq 'lname') {
			$compset_longname = $name_node->textContent();
		    }		    
		}
	    }
	    $pes_setby = $node_file->getAttribute('component');
	    $compsets_file = $file;
	    $config->set('COMPSETS_SPEC_FILE', "$compsets_file");
	    $config->set('COMPSET', "$compset_longname");
	    last;
	} 

    }
    if (! defined $pes_setby) {
	my $outstr = "ERROR create_newcase: no compset match was found in any of the following files:";
	foreach my $node_file (@nodes) {
	    my $file = $node_file->textContent();
	    $file =~ s/\$CIMEROOT/$cimeroot/;
	    $file =~ s/\$SRCROOT/$srcroot/;
	    $file =~ s/\$MODEL/$model/;
	    $outstr .= "$file\n";
	}
	$logger->logdie ($outstr);
    } else {
	$logger->info( "File specifying possible compsets: $compsets_file ");
	$logger->info( "Primary component (specifies possible compsets, pelayouts and pio settings): $pes_setby ");
	$logger->info("Compset: $compset_longname ");
    }   

    return ($pes_setby, $support_level);
}


#-------------------------------------------------------------------------------
sub getGridLongname
{
    my ($grid_input, $config) = @_;

    my ($grid_longname, $grid_shortname, $grid_aliasname);
    my $compset_match;

    my $grids_file = $config->get('GRIDS_SPEC_FILE');

    my $xml = XML::LibXML->new( no_blanks => 1)->parse_file($grids_file);
    my @nodes_alias = $xml->findnodes(".//grid[alias=\"$grid_input\"]");
    my @nodes_sname = $xml->findnodes(".//grid[sname=\"$grid_input\"]");
    my @nodes_lname = $xml->findnodes(".//grid[lname=\"$grid_input\"]");

    my $grid_node;
    if (@nodes_alias) {
	$grid_node = $nodes_alias[0];
    } elsif (@nodes_lname) {
	$grid_node = $nodes_lname[0];
    } elsif (@nodes_sname) {
	$grid_node = $nodes_sname[0];
    } else { 
	die " ERROR: no supported grid match for target grid $grid_input ";
    }
    
    # set the compset grid alias and longname
    foreach my $node ($grid_node->childNodes()) {
	my $name = $node->nodeName();
	my $value = $node->textContent();
	if ($name eq 'lname') {$grid_longname   = $node->textContent();}
	if ($name eq 'sname') {$grid_shortname  = $node->textContent();}
	if ($name eq 'alias') {$grid_aliasname  = $node->textContent();}
    }	

    # determine compgrid hash (global variable)
    # Assume the following order for specifying a grid name
    #  a%aname_l%lname_oi%oiname_r%rname_m%mname_g%gname_w%wname

    $grid_longname =~ /(a%)(.+)(_l%)/ ; $compgrid{'atm'}  = $2;
    $grid_longname =~ /(l%)(.+)(_oi%)/; $compgrid{'lnd'}  = $2;
    $grid_longname =~ /(oi%)(.+)(_r%)/; $compgrid{'ocn'}  = $2; $compgrid{'ice'} = $2;
    $grid_longname =~ /(r%)(.+)(_m%)/ ; $compgrid{'rof'}  = $2; 
    $grid_longname =~ /(g%)(.+)(_w%)/ ; $compgrid{'glc'}  = $2; 
    $grid_longname =~ /(w%)(.+)$/     ; $compgrid{'wav'}  = $2; 
    $grid_longname =~ /(m%)(.+)(_g%)/ ; $compgrid{'mask'} = $2; 

    my @nodes = $xml->findnodes(".//grid[lname=\"$grid_longname\"]");
    if ($#nodes != 0) {
	die "ERROR ConfigCompsetGrid::checkGrid : no match found for $grid_longname \n";
    } 
    my $attr = $nodes[0]->getAttribute('compset');
    if (defined $attr) {
	my $compset = $config->get('COMPSET');
	if ($compset !~ m/$attr/) {
	    die "ERROR ConfigCompsetGrid::getGridLongame $grid_longname is not supported for $compset \n";
	}
    }

    return ($grid_longname);
}

#-----------------------------------------------------------------------------------------------
sub setGridDomain
{
    my ($grids_file, $config) = @_;
    
    my $parser = XML::LibXML->new( no_blanks => 1);
    my $xml = $parser->parse_file($grids_file);
  
    # set the component grid names
    foreach my $key (keys %compgrid) {
	my $comp = uc $key;
	my $grid = $compgrid{$key};
	my $mask = $compgrid{'mask'};
	
	if ($config->is_valid_name("${comp}_GRID")) {
	    my @nodes = $xml->findnodes(".//domain[\@name=\"$grid\"]");
	    if (@nodes) {$config->set("${comp}_GRID", $nodes[0]->getAttribute('name'));}
	}	
	
	if ($config->is_valid_name("${comp}_NX")) {
	    my @nodes = $xml->findnodes(".//domain[\@name=\"$grid\"]/nx");
	    if (@nodes) {$config->set("${comp}_NX", $nodes[0]->textContent());}
	}
	
	if ($config->is_valid_name("${comp}_NY")) {
	    my @nodes = $xml->findnodes(".//domain[\@name=\"$grid\"]/ny");
	    if (@nodes) { $config->set("${comp}_NY", $nodes[0]->textContent());}
	}
	
	# Note domain files are only specified for ATM, LND, ICE and OCN components 
	if ($comp eq 'ATM' || $comp eq 'LND' || $comp eq 'ICE' || $comp eq 'OCN') {
	    my @nodes;
	    my $var = "$comp" . "_DOMAIN_FILE";
	    if ($config->is_valid_name($var)) {
		if ($comp eq 'ATM' || $comp eq 'LND') {	    
		    @nodes = $xml->findnodes(".//domain[\@name=\"$grid\"]/file[\@lnd_mask=\"$mask\"]");
		} else {
		    @nodes = $xml->findnodes(".//domain[\@name=\"$grid\"]/file[\@ocn_mask=\"$mask\"]");
		}
		if (@nodes) {$config->set($var,$nodes[0]->textContent());}
	    }
	    my $var = "$comp" . "_DOMAIN_PATH";
	    if ($config->is_valid_name($var)) {
		if ($comp eq 'ATM' || $comp eq 'LND') {	    
		    @nodes = $xml->findnodes(".//domain[\@name=\"$grid\"]/path[\@lnd_mask=\"$mask\"]");
		} else {
		    @nodes = $xml->findnodes(".//domain[\@name=\"$grid\"]/path[\@ocn_mask=\"$mask\"]");
		}
		if (@nodes) {$config->set($var,$nodes[0]->textContent());}
	    }
	}
	
    }

# TODO - does this still need to be here ???    
#    if ($compgrid{'cism'} ne 'null') {
#	$config->set('CISM_GRID',$compgrid{'cism'});
#    }
}

#-------------------------------------------------------------------------------
sub setGridMaps
{
    # set grid mapping variables
    my ($grids_file, $config) = @_;
    
    my $parser = XML::LibXML->new( no_blanks => 1);
    my $xml = $parser->parse_file($grids_file);

    # first set defaults
    my @nodes = $xml->findnodes(".//grid_mappings/defaults/*");
    foreach my $node (@nodes) {
	my $name  = $node->nodeName();
	my $value = $node->textContent();
	$config->set($name, $value);
    }

    # overwrite the default values if appropriate
    my $atm_grid = $compgrid{'atm'};
    my $lnd_grid = $compgrid{'lnd'};
    my $rof_grid = $compgrid{'rof'};
    my $ocn_grid = $compgrid{'ocn'};
    my $ice_grid = $compgrid{'ice'};
    my $glc_grid = $compgrid{'glc'};
    my $wav_grid = $compgrid{'wav'};

    @nodes = $xml->findnodes(".//gridmap[\@atm_grid=\"$atm_grid\" and \@ocn_grid=\"$ocn_grid\"]/ATM2OCN_FMAPNAME");
    if (@nodes) {$config->set('ATM2OCN_FMAPNAME',$nodes[0]->textContent())}

    @nodes = $xml->findnodes(".//gridmap[\@atm_grid=\"$atm_grid\" and \@ocn_grid=\"$ocn_grid\"]/ATM2OCN_SMAPNAME");
    if (@nodes) {$config->set('ATM2OCN_SMAPNAME',$nodes[0]->textContent())}

    @nodes = $xml->findnodes(".//gridmap[\@atm_grid=\"$atm_grid\" and \@ocn_grid=\"$ocn_grid\"]/ATM2OCN_VMAPNAME");
    if (@nodes) {$config->set('ATM2OCN_VMAPNAME',$nodes[0]->textContent())}

    @nodes = $xml->findnodes(".//gridmap[\@atm_grid=\"$atm_grid\" and \@ocn_grid=\"$ocn_grid\"]/OCN2ATM_FMAPNAME");
    if (@nodes) {$config->set('OCN2ATM_FMAPNAME',$nodes[0]->textContent())}

    @nodes = $xml->findnodes(".//gridmap[\@atm_grid=\"$atm_grid\" and \@ocn_grid=\"$ocn_grid\"]/OCN2ATM_SMAPNAME");
    if (@nodes) {$config->set('OCN2ATM_SMAPNAME',$nodes[0]->textContent())}

    @nodes = $xml->findnodes(".//gridmap[\@atm_grid=\"$atm_grid\" and \@lnd_grid=\"$lnd_grid\"]/ATM2LND_SMAPNAME");
    if (@nodes) {$config->set('ATM2LND_SMAPNAME',$nodes[0]->textContent())}

    @nodes = $xml->findnodes(".//gridmap[\@atm_grid=\"$atm_grid\" and \@lnd_grid=\"$lnd_grid\"]/ATM2LND_FMAPNAME");
    if (@nodes) {$config->set('ATM2LND_FMAPNAME',$nodes[0]->textContent())}

    @nodes = $xml->findnodes(".//gridmap[\@atm_grid=\"$atm_grid\" and \@lnd_grid=\"$lnd_grid\"]/LND2ATM_SMAPNAME");
    if (@nodes) {$config->set('ATM2LND_SMAPNAME',$nodes[0]->textContent())}

    @nodes = $xml->findnodes(".//gridmap[\@atm_grid=\"$atm_grid\" and \@lnd_grid=\"$lnd_grid\"]/LND2ATM_FMAPNAME");
    if (@nodes) {$config->set('LND2ATM_FMAPNAME',$nodes[0]->textContent())}

    @nodes = $xml->findnodes(".//gridmap[\@lnd_grid=\"$lnd_grid\" and \@rof_grid=\"$rof_grid\"]/LND2ROF_FMAPNAME");
    if (@nodes) {$config->set('LND2ROF_FMAPNAME',$nodes[0]->textContent())}

    @nodes = $xml->findnodes(".//gridmap[\@lnd_grid=\"$lnd_grid\" and \@rof_grid=\"$rof_grid\"]/ROF2LND_FMAPNAME");
    if (@nodes) {$config->set('ROF2LND_FMAPNAME',$nodes[0]->textContent())}

    @nodes = $xml->findnodes(".//gridmap[\@ocn_grid=\"$ocn_grid\" and \@rof_grid=\"$rof_grid\"]/ROF2OCN_FMAPNAME");
    if (@nodes) {$config->set('ROF2OCN_FMAPNAME',$nodes[0]->textContent())}

    @nodes = $xml->findnodes(".//gridmap[\@ocn_grid=\"$ocn_grid\" and \@rof_grid=\"$rof_grid\"]/ROF2OCN_RMAPNAME");
    if (@nodes) {$config->set('ROF2OCN_RMAPNAME',$nodes[0]->textContent())}

    @nodes = $xml->findnodes(".//gridmap[\@ocn_grid=\"$ocn_grid\" and \@rof_grid=\"$rof_grid\" and \@lnd_grid=\"$lnd_grid\"]/XROF_FLOOD_MODE");
    if (@nodes) {$config->set('XROF_FLOOD_MODE',$nodes[0]->textContent())}

    @nodes = $xml->findnodes(".//gridmap[\@lnd_grid=\"$lnd_grid\" and \@glc_grid=\"$glc_grid\"]/LND2GLC_SMAPNAME");
    if (@nodes) {$config->set('LND2GLC_SMAPNAME',$nodes[0]->textContent())}

    @nodes = $xml->findnodes(".//gridmap[\@lnd_grid=\"$lnd_grid\" and \@glc_grid=\"$glc_grid\"]/LND2GLC_FMAPNAME");
    if (@nodes) {$config->set('LND2GLC_FMAPNAME',$nodes[0]->textContent())}

    @nodes = $xml->findnodes(".//gridmap[\@glc_grid=\"$glc_grid\" and \@lnd_grid=\"$lnd_grid\"]/GLC2LND_SMAPNAME");
    if (@nodes) {$config->set('GLC2LND_SMAPNAME',$nodes[0]->textContent())}

    @nodes = $xml->findnodes(".//gridmap[\@glc_grid=\"$glc_grid\" and \@lnd_grid=\"$lnd_grid\"]/GLC2LND_FMAPNAME");
    if (@nodes) {$config->set('GLC2LND_FMAPNAME',$nodes[0]->textContent())}

    @nodes = $xml->findnodes(".//gridmap[\@glc_grid=\"$glc_grid\" and \@ice_grid=\"$ice_grid\"]/GLC2ICE_RMAPNAME");
    if (@nodes) {$config->set('GLC2ICE_RMAPNAME',$nodes[0]->textContent())}

    @nodes = $xml->findnodes(".//gridmap[\@glc_grid=\"$glc_grid\" and \@ocn_grid=\"$ocn_grid\"]/GLC2OCN_RMAPNAME");
    if (@nodes) {$config->set('GLC2OCN_RMAPNAME',$nodes[0]->textContent())}
}

#-------------------------------------------------------------------------------
sub setCompsetGeneralVars
{
    # Determine general compset variables
    my ($config) = @_;

    my $compset_longname = $config->get('COMPSET');
    my $grid_longname    = $config->get('GRID');
    my $compsets_file    = $config->get('COMPSETS_SPEC_FILE');	

    my $parser = XML::LibXML->new( no_blanks => 1);
    my $xml_compset = $parser->parse_file($compsets_file);

    # first check that compset is supported for target grid
    my @nodes = $xml_compset->findnodes(".//compset[lname=\"$compset_longname\"]");
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

    foreach my $node ($xml_compset->findnodes(".//entry")) {
	my $id = $node->getAttribute('id');
	if ($config->is_valid_name($id)) {
	  VALUE: foreach my $child ($node->findnodes("./values/value")) {
	      my $compset_match = $child->getAttribute('compset');
	      my $grid_match    = $child->getAttribute('grid');
	      
	      if (! defined $compset_match) {
		  $compset_match = $compset_longname;
	      }
	      if (! defined $grid_match) {
		  $grid_match    = $grid_longname;
	      }
	      
	      if ($compset_longname =~ m/$compset_match/ && $grid_longname =~ m/$grid_match/) {	    
		  my $new_val =  $child->textContent();
		  $newxml{$id} = $new_val; 
		  $config->set($id, $new_val);
		  last VALUE;
	      }
	  }
	} else {
	    my $cimeroot = $config->get('CIMEROOT');
	    $logger->logdie( "ERROR : $id is not a valid name in $compsets_file 
	    *** See possible values in file $cimeroot/driver_cpl/cimeconfig/config_cime.xml ***");
	}
    }
}

#-------------------------------------------------------------------------------
sub setComponent {

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

    my ($setup_comp_file, $config) = @_;

    my $grid_longname	   = $config->get('GRID');
    my $compset_longname   = $config->get('COMPSET');

    my $xml = XML::LibXML->new( no_blanks => 1)->parse_file($setup_comp_file);

    # First make sure that there is a match for one entry in the description element that
    # matches the compset longname
    my $found_match;
    foreach my $node ($xml->findnodes(".//description/desc")) {
	my $compset_match = $node->getAttribute('compset');
	if ($compset_longname =~ m/$compset_match/) {
	    $found_match = 1;
	    last;
	}
    }
    if (! defined $found_match) {
	$logger->fatal( "ERROR ConfigCompsetGrid::setComponent: no match found in desc elements in file 
	         $setup_comp_file 
	      for $compset_longname  
	        Possible matches are: ");
	foreach my $node ($xml->findnodes(".//description/desc")) {
	    my $match = $node->getAttribute('compset');
	    my $desc  = $node->textContent();
	    $logger->fatal("     $match: $desc ");
	}
	$logger->logdie ("Exiting"); 
    }

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
		    unless ($config->is_valid_name($name)) {
			die "ERROR ConfigCompsetGrid::setComponent: $name is not a valid name \n";
		    }

		    my $input_val = $value_node->textContent();
		    my $current_val = $config->get($name);


		    if (@additive) {
			$new_val = _setComponentAdditiveNewVal($current_val, $input_val);
		    } elsif (@merge) {
			$new_val = _setComponentMergeNewVal($current_val, $input_val);
		    } else {
			$new_val = $input_val;
		    }
		    $newxml{$name} = $new_val; 
		    $config->set($name, $new_val);

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
	    if ($desc_comp =~ m/$value/) {
		# do nothing
	    } else {
		if($value =~ /\n/){
		    $desc_comp = $desc_comp . "\n$value\n";
		} else {
		    $desc_comp = $desc_comp . " $value";
		}
	    }
	}
    }

   
    # Special case for CLMUSRDAT in compset name
    if ($compset_longname =~ /(.+CLMUSRDAT%)(.*)/){
     	$config->set('CLM_USRDAT_NAME',$2);
     	$newxml{"CLM_USR_DATNAME"} = $2;
    } 
}

#-------------------------------------------------------------------------------
sub printGridCompsetInfo 
{
    my ($grids_file,  $config) = @_;

    my $caseroot	   = $config->get('CASEROOT');
    my $grid_longname	   = $config->get('GRID');
    my $compset_longname   = $config->get('COMPSET');

    my @grids;
    my $desc_grid = "";

    my $atm_nx   = $config->get('ATM_NX'); 
    my $atm_ny   = $config->get('ATM_NY'); 
    my $atm_grid = $config->get('ATM_GRID');
    push (@grids, $atm_grid);

    my $lnd_nx   = $config->get('LND_NX'); 
    my $lnd_ny   = $config->get('LND_NY'); 
    my $lnd_grid = $config->get('LND_GRID');
    push (@grids, $lnd_grid);

    my $ice_nx   = $config->get('ICE_NX'); 
    my $ice_ny   = $config->get('ICE_NY'); 
    my $ice_grid = $config->get('ICE_GRID');
    push (@grids, $ice_grid);

    my $ocn_nx   = $config->get('OCN_NX'); 
    my $ocn_ny   = $config->get('OCN_NY'); 
    my $ocn_grid = $config->get('OCN_GRID');
    push (@grids, $ocn_grid);

    my $rof_nx   = $config->get('ROF_NX'); 
    my $rof_ny   = $config->get('ROF_NY'); 
    my $rof_grid = $config->get('ROF_GRID');
    push (@grids, $rof_grid);

    my $glc_nx   = $config->get('GLC_NX'); 
    my $glc_ny   = $config->get('GLC_NY'); 
    my $glc_grid = $config->get('GLC_GRID');
    push (@grids, $glc_grid);

    my $wav_nx    = $config->get('WAV_NX'); 
    my $wav_ny    = $config->get('WAV_NY'); 
    my $wav_grid  = $config->get('WAV_GRID');
    push (@grids, $wav_grid);

    my $mask_grid = $config->get('MASK_GRID');
    push (@grids, $mask_grid);

    my $parser = XML::LibXML->new( no_blanks => 1);
    my $xml_domain = $parser->parse_file($grids_file);

    my $support; 
    foreach my $grid (@grids)
    {
	my @nodes = $xml_domain->findnodes(".//domain[\@name=\"grid\"]/support");
	if (@nodes) {
	    my $value = @nodes[0]->textContent();
	    my $name  = @nodes[0]->nodeName();
	    $ .= "Grid ($name): $value\n";
	}
	my @nodes = $xml_domain->findnodes(".//domain[\@name=\"$grid\"]/desc");
	if (@nodes) {
	    my $value = @nodes[0]->textContent();
	    if ($desc_grid !~ m/$value/) {
		$desc_grid = $desc_grid . " $value" ;
	    }
	}
    }

    my $compset_longname = $config->get('COMPSET');

    my $fh_case = new IO::File;
    $fh_case->open(">>$caseroot/README.case") or die "can't open file: README.case\n";

    my $fh_stdout = *STDOUT;
    my $outstr;
    $outstr .= "Component set: longname \n";
    $outstr .= "  $compset_longname \n";
    $outstr .= "Component set Description: \n";
    $outstr .= " $desc_comp \n";
    $outstr .= "Grid: \n";
    $outstr .= "  $grid_longname \n";
    $outstr .= "  ATM_GRID = $atm_grid  NX_ATM=$atm_nx NY_ATM=$atm_ny \n";
    $outstr .= "  LND_GRID = $lnd_grid  NX_LND=$lnd_nx NX_LND=$lnd_ny \n";
    $outstr .= "  ICE_GRID = $ice_grid  NX_ICE=$ice_nx NX_ICE=$ice_ny \n";
    $outstr .= "  OCN_GRID = $ocn_grid  NX_OCN=$ocn_nx NX_OCN=$ocn_ny \n";
    $outstr .= "  ROF_GRID = $rof_grid  NX_ROF=$rof_nx NX_ROF=$rof_ny \n";
    $outstr .= "  GLC_GRID = $glc_grid  NX_GLC=$glc_nx NX_GLC=$glc_ny \n";
    $outstr .= "  WAV_GRID = $wav_grid  NX_WAV=$wav_nx NX_WAV=$wav_ny \n";
    $outstr .= "Grid Description: \n";
    $outstr .= "  $desc_grid \n";
    $outstr .= "Non-Default Options: \n";
    my @ids = keys %newxml;
    foreach my $id (sort @ids) {
	my $value = $newxml{$id};
	$outstr .=     "  $id: $value \n";
    } 
    print $fh_case $outstr;
    $logger->info($outstr);
    $fh_case->close();
}	


#-----------------------------------------------------------------------------------------------
#                               Private routines
#-----------------------------------------------------------------------------------------------
sub _setComponentAdditiveNewVal 
{
    # input arguments
    my ($current_val, $input_val) = @_;

    my $new_val;
    # if $current_val is not empty
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
	    if (defined $optname) {
	    
		# Check current options for values to replace.
		foreach my $curopt (@curopts) {
		    if ($curopt =~ m/^$optname\s/) {
			$name_found = 1;
			# Substitute, adding one space just in case.
			$new_val =~ s/$curopt/$nameopt /;
		    }
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
sub _clean
{
    my ($name) = @_;
    
    $name =~ s/^\s+//; # strip any leading whitespace 
    $name =~ s/\s+$//; # strip any trailing whitespace
    return ($name);
}

1;
