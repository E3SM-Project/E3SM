package SetupTools;
my $pkg_nm = 'SetupTools';

use strict;
use XML::Lite;


sub expand_env_var
{
    my ($value,$vars) = @_;

    if ($value =~ /\$\{*([\w_]+)}*(.*)$/) {
       my $subst = $vars->{$1};
       $subst = $ENV{$1} unless defined $subst;
       $value =~ s/\$\{*${1}\}*/$subst/g;
    }
    
    $value = expand_env_var($value,$vars) if ($value =~ /\$\{*[\w_]+\}*.*$/) ;
    return $value; 
}       

sub getxmlvars
{
#-----------------------------------------------------------------------------------------------
# Read $caseroot xml files - put restuls in %xmlvars hash
#-----------------------------------------------------------------------------------------------
    my ($caseroot,$xmlvars) = @_;
    my @files = <${caseroot}/*xml>;

    foreach my $file (@files) {
	my $xml = XML::Lite->new( "$file" ) or die "Could not open file $file";
	my @e = $xml->elements_by_name('entry');
	while ( my $e = shift @e ) {
	    my %a = $e->get_attributes();
	    $xmlvars->{$a{'id'}} = $a{'value'};
	}
    }

}


#-------------------------------------------------------------------------------
sub set_compiler
{
    my ($os,$compiler_file, $compiler, $machine, $mpilib, $print, $macrosfile) = @_;

    print "$compiler_file $compiler $machine\n" ;#if($print>1);

#Parse the config_compiler.xml file into a Macros file for the given machine and compiler
    my $xml = XML::Lite->new( $compiler_file );


    my $root = $xml->root_element();
    # Check for valid root node
    die "no root defined" unless defined($root);

    my $name = $root->get_name();
    $name eq "config_compilers" or die
	"file $compiler_file is not a compiler parameters file\n";

#
# Read all settings for the given compiler and optionally the Machine
# more general settings are overwritten by more specific ones
# 
    my @elem = $xml->elements_by_name("compiler");
    my %a=();
    my $e;
    my @compiler_settings;
    foreach $e (@elem){
	%a = $e->get_attributes();
# Only pick up settings for which the defined attributes match
	next if(defined $a{COMPILER} && $a{COMPILER} ne $compiler);
	next if(defined $a{MACH} && $a{MACH} ne $machine);
	next if(defined $a{OS} && $a{OS} ne $os);
	next if(defined $a{MPILIB} && $a{MPILIB} ne $mpilib);

	print "compiler $compiler $a{COMPILER} $a{OS} $a{MACH}\n" if($print>1);
	push(@compiler_settings ,$e->get_children());
    }

    die "set_compiler: unrecognized compiler" unless($#compiler_settings);

    my $flag;
    my @allfkeys;
    my $macros;
    $macros->{_COND_}={};

#
# Parse the xml settings into the $macros hash structure
# put conditional settings in the _COND_ portion of the hash
# and handle them seperately
#

    foreach $flag (@compiler_settings){
	my $name = $flag->get_name();
	my $val = $flag->get_text();
	%a = $flag->get_attributes();
	my @keys =  keys %a;

	my $hash = $macros->{_COND_};
	if($#keys<0){
	    if($name =~ /^ADD_(.*)$/){
		my $basename = $1;
		if(defined $macros->{$basename}){
		    $macros->{$basename} .= $val ;
		}elsif(defined $macros->{$name}){
		    $macros->{$name}.=$val;
		}else{
		    $macros->{$name}=$val;
		}
		print "$basename+=$val\n" if($print>1);
	    }else{
		$macros->{$name}=$val;
		print "$name:=$val\n" if($print>1);
	    }
	}else{
	    my $key;
	    foreach $key (@keys){
		unless(defined $hash->{$key}{$a{$key}}){
		    $hash->{$key}{$a{$key}}={} ;
		}
		$hash=$hash->{$key}{$a{$key}};
	    }
	    $hash->{$name}=$val;

	}
    }
    my $compcpp = $compiler;
    $compcpp =~ tr/a-z/A-Z/;
    $macros->{ADD_CPPDEFS} .= " -D$os -DCPR$compcpp ";

    if(! defined($macros->{MPI_PATH}) && defined($ENV{MPI_PATH})){
      print "Setting MPI_PATH from Environment\n";
      $macros->{MPI_PATH}=$ENV{MPI_PATH};
    }
    if(! defined($macros->{NETCDF_PATH}) && defined($ENV{NETCDF_PATH})){
      print "Setting NETCDF_PATH from Environment\n";
      $macros->{NETCDF_PATH}=$ENV{NETCDF_PATH};
    }

#    print Dumper($macros);


#open Macros file to write
    open MACROS,">$macrosfile" or die "Could not open file $macrosfile to write";

    print MACROS "#\n# Makefile Macros generated from $compiler_file using\n";
    print MACROS "# COMPILER=$compiler\n";
    print MACROS "# OS=$os\n";
    print MACROS "# MACH=$machine\n#\n";

#    print Dumper($macros);
    my @keys =  sort keys %{$macros};

# print the settings out to the Macros file

    foreach (@keys){
	next if($_ eq '_COND_');
	if($_ =~ /^ADD_(.*)/){
	    print MACROS "$1+=".$macros->{$_}."\n\n";
	}else{
	    print MACROS "$_:=".$macros->{$_}."\n\n";
	}
    }
# Recursively print the conditionals, combining tests to avoid repetition    

    parse_hash($macros->{_COND_}, 0);
    close MACROS;
}

sub parse_hash
{
    my($href,$depth) = @_;
    my @keys = keys %{$href};
    my $k1;

    my $width=2*$depth;
    foreach $k1 (@keys){
	if(ref($href->{$k1})){
	    my $k2;
	    foreach $k2 (keys %{$href->{$k1}}){
		printf(MACROS "%${width}s"," ") if($width>0);
		printf(MACROS "ifeq (\$(%s), %s) \n",$k1,$k2);
		parse_hash($href->{$k1}{$k2},$depth+1);
	    }
	}else{
	    if($k1=~/ADD_(.*)/){
		printf(MACROS "%${width}s %s +=%s\n"," ",$1,$href->{$k1});
	    }else{
		printf(MACROS "%${width}s %s :=%s\n"," ",$k1,$href->{$k1});
	    }
	}
    }
    $width-=2;
    printf(MACROS "%${width}s"," ") if($width>0);
    printf(MACROS "endif\n\n") if($depth>0) ;    
   
}


#-------------------------------------------------------------------------------


1;
