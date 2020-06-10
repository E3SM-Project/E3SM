package SetupTools;
my $pkg_nm = 'SetupTools';

use strict;
use XML::LibXML;
use Data::Dumper;

use Log::Log4perl qw(get_logger);
my $logger;

BEGIN{
    $logger = get_logger();
}
#-----------------------------------------------------------------------------------------------
# SYNOPSIS
#
# The following routines are contained here
#
# is_valid_value()
#       Returns true if the specified parameter name is contained in
#       the configuration definition file, and either 1) the specified value is
#       listed as a valid_value in the definition file, or 2) the definition file
#       doesn't specify the valid values.
#

#-------------------------------------------------------------------------------
sub create_namelist_infile
{
    my ($caseroot, $user_nl_file, $namelist_infile, $infile_text) = @_;

    open( file_usernl,"<${user_nl_file}"   ) or $logger->logdie( "*** can't open file: $user_nl_file");
    open( file_infile,">${namelist_infile}") or $logger->logdie( "*** can't open file: $namelist_infile");

    print file_infile "\&comp_inparm \n";

    if ($infile_text) {
	print file_infile $infile_text;
    }

    while (my $line = <file_usernl>) {
	if ($line =~ /^[\&\/\!]/ ) {
	    next; # do nothing
	} elsif ($line =~ /\$([\w\_]+)/) {
	    my %xmlvars = ();
	    getxmlvars(${caseroot},\%xmlvars);
	    $line = expand_xml_var($line, \%xmlvars);
	}
	print file_infile "$line";
    }
    print file_infile "\/ \n";

    close(file_infile);
    close(file_usernl);
}

#-------------------------------------------------------------------------------
sub expand_env
{
    my ($value, $xmlvars_ref) = @_;

    if ($value =~ /\$\{*([\w_]+)}*(.*)$/) {
       my $subst = $xmlvars_ref->{$1};
       $subst = $ENV{$1} unless defined $subst;
       $value =~ s/\$\{*${1}\}*/$subst/g;
    }

    if ($value =~ /\$\{*[\w_]+\}*.*$/) {
	$value = expand_env($value, $xmlvars_ref)
    }
    return $value;
}

#-------------------------------------------------------------------------------
sub expand_xml_var
{
    my ($value, $xmlvars_ref) = @_;

    if($value =~ /\$ENV\{(.*)\}/){
	my $subst = $ENV{$1};
	if(defined $subst){
	    $value =~ s/\$ENV\{*${1}\}/$subst/g;
	}else{
	    $logger->warn ("No environment variable found for $1");
	}
    }
    if ($value =~ /\$\{*([\w_]+)}*(.*)$/) {
	my $found_xml;
	while ( my ($key, $subst) = each(%$xmlvars_ref) ) {
	    if ($key eq $1) {
		$value =~ s/\$\{*${1}\}*/$subst/g;
		$found_xml = 'yes';
		last;
	    }
	}
	if (! defined($found_xml) ) {
	    $value = expand_env($value, $xmlvars_ref) ;
	}
    }

    if ($value =~ /\$\{*[\w_]+\}*.*$/) {
	$value = expand_xml_var($value, $xmlvars_ref) ;
    }
    return $value;
}

#-------------------------------------------------------------------------------
sub getAllResolved
{
    # hash for all the parsers, and a hash for  all the config variables.
    my %parsers;
    my %masterconfig;

    # Get all the env*.xml files into an array...
    my @xmlfiles = qw( env_build.xml env_case.xml env_mach_pes.xml env_run.xml env_batch.xml);
    push(@xmlfiles, "env_test.xml") if(-e "./env_test.xml");
    push(@xmlfiles, "env_archive.xml") if(-e "./env_archive.xml");

    # Set up a new XML::LibXML parser for each xml file.
    foreach my $basefile(@xmlfiles)
    {
	my $xml = XML::LibXML->new()->parse_file($basefile);
	$parsers{$basefile} = $xml;
    }

    # find all the variable nodes.
    foreach my $basefile(@xmlfiles)
    {
	my $parser = $parsers{$basefile};
	my @nodes = $parser->findnodes("//entry");
	foreach my $node(@nodes)
	{
	    my $id = $node->getAttribute('id');
	    my $value = $node->getAttribute('value');
	    # if the variable value has an unresolved variable,
	    # we need to find it in whatever file it might be in.
	    $value = _resolveValues($value, \%parsers);
	    $masterconfig{$id} = $value;
	}
    }
    return %masterconfig;
}

#-------------------------------------------------------------------------------
# Resolve a single unresolved variable
#-------------------------------------------------------------------------------
sub getSingleResolved
{
    my $variable = shift;
    my $id_to_find = $variable;
    $id_to_find =~ s/\$//g;
    my %parsers;
    my $value;
    my @xmlfiles = qw (env_build.xml env_case.xml env_mach_pes.xml env_run.xml );
    push(@xmlfiles, "env_test.xml") if ( -e "./env_test.xml");
    push(@xmlfiles, "env_archive.xml") if ( -e "./env_archive.xml");

    foreach my $basefile(@xmlfiles)
    {
        my $xml = XML::LibXML->new()->parse_file($basefile);
        $parsers{$basefile} = $xml;
    }

    my @nodes;
    foreach my $basefile(@xmlfiles)
    {
        my $parser = $parsers{$basefile};
        my @nodes = $parser->findnodes("//entry[\@id=\'$id_to_find\']");
        if(@nodes)
        {
            my $node = $nodes[0];
            $value = $node->getAttribute('value');
            if($value =~/\$/)
            {
                $value = _resolValues($value, \%parsers);
            }

        }
        else
        {
            next;
        }
    }
    if(defined $value)
    {
        return $value;
    }
    else
    {
        return undef;
    }
}

#-------------------------------------------------------------------------------
sub getxmlvars
{
    # Read $caseroot xml files - put restuls in %xmlvars hash
    my ($caseroot, $xmlvars) = @_;

    my $parser = XML::LibXML->new( no_blanks => 1);

    my @files = <${caseroot}/*xml>;
    foreach my $file (@files) {
	my $xml_variables = $parser->parse_file($file);
	my @nodes = $xml_variables->findnodes(".//entry");
	foreach my $node (@nodes) {
	    my $id    = $node->getAttribute('id');
	    my $value = $node->getAttribute('value');
	    $xmlvars->{$id} = $value;
	}
    }
# These represent a workaround for a problem in resolving variables in perl

    if (defined ($xmlvars->{CIME_OUTPUT_ROOT})){
	$xmlvars->{CIME_OUTPUT_ROOT}=expand_xml_var($xmlvars->{CIME_OUTPUT_ROOT}, $xmlvars);
    }
    if (defined ($xmlvars->{DIN_LOC_ROOT})){
	$xmlvars->{DIN_LOC_ROOT}=expand_xml_var($xmlvars->{DIN_LOC_ROOT}, $xmlvars);
    }
}

#-----------------------------------------------------------------------------------------------
sub is_valid_value
{
    # Check if the input value is a valid value
    my ($id, $value, $valid_values, $is_list_value) = @_;

    # Check that a list value is not supplied when parameter takes a scalar value.
    unless ($is_list_value) {  # this conditional is satisfied when the list attribute is false, i.e., for scalars
	if ($value =~ /.*,.*/) {
	    # the pattern matches when $value contains a comma, i.e., is a list
 	    $logger->logdie( "Error is_valid_value; variable $id is a scalar but has a list value $value ");
	}
   }

    # Check that the value is valid
    # if no valid values are specified, then $value is automatically valid
    if ( $valid_values ne "" ) {
	if ($is_list_value) {
	    unless (_list_value_ok($value, $valid_values)) {
		$logger->logdie( "ERROR is_valid_value: $id has value $value which is not a valid value ");
	    }
	} else {
	    unless (_value_ok($value, $valid_values)) {
		$logger->logdie("ERROR is_valid_value: $id has value $value which is not a valid value ");
	    }
	}

    }
    return 1;
}

#-----------------------------------------------------------------------------------------------
sub validate_variable_value
{
    # Validate that a given value matches the expected input type definition
    # Expected description of keys for the input type hash is:
    #      type           type description (char, logical, integer, or real)       (string)
    #      strlen         Length of string (if type char)                          (integer)
    #      validValues    Reference to array of valid values                       (string)
    #
    my ($var, $value, $type_ref) = @_;
    my $nm = "validate_variable_value";

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

    # Ensure type hash has required variables
    if ( ref($type_ref) !~ /HASH/ ) {
	$logger->logdie("ERROR: in $nm : Input type is not a HASH reference.");
    }
    foreach my $item ( "type", "validValues", "strlen" ) {
	if ( ! exists($$type_ref{$item}) ) {
	    $logger->logdie( "ERROR: in $nm: Variable name $item not defined in input type hash.");
	}
    }
    # If string check that less than defined string length
    my $str_len = 0;
    if ( $$type_ref{'type'} eq "char" ) {
	$str_len = $$type_ref{'strlen'};
	if ( length($value) > $str_len ) {
	    $logger->logdie( "ERROR: in $nm Variable name $var " .
			     "has a string element that is too long: $value");
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
		$logger->logdie( "ERROR: in $nm (package $pkg_nm): Type of variable name $var is " .
		    "not a valid FORTRAN type (logical, integer, real, or char).");
	    }
	    if ( $i =~ /USERDEFINED/) {
		$logger->warn ("WARNING: in $nm (package $pkg_nm): Variable name $var " .
		    "has not been defined");
	    }elsif ( $i !~ /^\s*(${compare})$/  and ! $value =~ /^\$.*/ ) {   # give it a pass if set to another variable
		$logger->logdie( "ERROR: in $nm (package $pkg_nm): Variable name $var " .
		    "has a value ($i) that is not a valid type " . $$type_ref{'type'} );
	    }
	}
    }
}

#-----------------------------------------------------------------------------------------------
#                               Private routines
#-----------------------------------------------------------------------------------------------
sub _parse_hash
{
    my($href,$depth, $output_format, $cmakedebug) = @_;
    my @keys = keys %{$href};
    my $k1;

    my $width=2*$depth;
    foreach $k1 (@keys){
	if(ref($href->{$k1})){
	    my $k2;
	    if($output_format eq "make" or $k1 =~ /DEBUG/){
		foreach $k2 (keys %{$href->{$k1}}){
		    if($output_format eq "make"){
			printf(MACROS "%${width}s"," ") if($width>0);
			printf(MACROS "ifeq (\$(%s), %s) \n",$k1,$k2);
		    }

		    _parse_hash($href->{$k1}{$k2},$depth+1, $output_format, $k2);
		}
	    }
	}else{
	    if($output_format eq "make"){
		if($k1=~/ADD_(.*)/){
		    printf(MACROS "%${width}s %s +=%s\n"," ",$1,$href->{$k1});
		}else{
		    printf(MACROS "%${width}s %s :=%s\n"," ",$k1,$href->{$k1});
		}
	    }else{
		my $value = $href->{$k1};
		$value =~ s/\(/\{/g;
		$value =~ s/\)/\}/g;
		if($cmakedebug =~ /TRUE/){
		    if($k1 =~ /CFLAGS/){
			printf(MACROS "add_flags(CMAKE_C_FLAGS_DEBUG $value)\n\n");
		    }elsif($k1 =~ /FFLAGS/){
			printf(MACROS "add_flags(CMAKE_Fortran_FLAGS_DEBUG $value)\n\n");
		    }elsif($k1 =~ /CPPDEF/){
			printf(MACROS "add_config_definitions(DEBUG $value)\n\n");
		    }elsif($_ =~ /SLIBS/ or $_ =~ /LDFLAGS/){
			print MACROS "add_flags(CMAKE_EXE_LINKER_FLAGS_DEBUG $value)\n\n";
		    }
		}else{
		    if($k1 =~ /CFLAGS/){
			printf( MACROS "add_flags(CMAKE_C_FLAGS_RELEASE $value)\n\n");
		    }elsif($k1 =~ /FFLAGS/){
			printf( MACROS "add_flags(CMAKE_Fortran_FLAGS_RELEASE $value)\n\n");
		    }elsif($k1 =~ /CPPDEF/){
			printf(MACROS "add_config_definitions(RELEASE $value)\n\n");
		    }elsif($_ =~ /SLIBS/ or $_ =~ /LDFLAGS/){
			print MACROS "add_flags(CMAKE_EXE_LINKER_FLAGS_RELEASE $value)\n\n";
		    }
		}
	    }
	}
    }
    $width-=2;
    if($output_format eq "make"){
	printf(MACROS "%${width}s"," ") if($width>0);
	printf(MACROS "endif\n\n") if($depth>0) ;
    }
}

#-----------------------------------------------------------------------------------------------
sub _resolveValues
{
    # Recursively resolve the unresolved vars in an variable value.
    # Check the value passed in, and if it still has an unresolved var, keep calling the function
    # until all pieces of the variable are resolved.

    my $value = shift;
    my $parsers = shift;

    #print "in _resolveValues: value: $value\n";
    # We want to resolve $values from either tthe other xml files, or
    # the value can come from the
    if($value =~ /(\$[\w_]+)/)
    {
# too noisy
#	$logger->debug( "in _resolveValues: value: $value\n");
	my $unresolved = $1;

	#print "need to resolve: $unresolved\n";
	my $needed = $unresolved;
	$needed =~ s/\$//g;

	my $found = 0;
	foreach my $parser(values %$parsers)
	{
	    my @resolveplease = $parser->findnodes("//entry[\@id=\'$needed\']");
	    if(@resolveplease)
	    {
		$found = 1;
		foreach my $r(@resolveplease)
		{
		    my $rid = $r->getAttribute('id');
		    my $rvalue = $r->getAttribute('value');
		    $value =~ s/\$$needed/$rvalue/g;
# too noisy
#		    $logger->debug( "value after substitution: $value\n");
		}
	    }
	}
	# Check the environment if not found in the xml files.
	if(!$found)
	{
	    if(exists $ENV{$needed})
	    {
		$found = 1;
		my $rvalue = $ENV{$needed};
		$value =~ s/\$$needed/$rvalue/g;
	    }
	}
	#if the value is not found in any of the xml files or in the environment, then
	# return undefined.
	if(!$found)
	{
	    return undef;
	}
	_resolveValues($value, $parsers);
    }
    else
    {
# too verbose
#	$logger->debug( "returning $value\n");
	return $value;
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
#-------------------------------------------------------------------------------

1;
