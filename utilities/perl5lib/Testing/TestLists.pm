package Testing::TestLists;
use Exporter;
use Data::Dumper;
use Testing::CESMTest;
#my @ISA = qw(Exporter);
#my @EXPORTOK = qw(findTestsForCase);
BEGIN
{
	use vars qw( $VERSION @ISA );
	$VERSION = '0.01';
	@ISA     = qw();
}
	
#-----------------------------------------------------------------------------------------------
#
#
#
#-----------------------------------------------------------------------------------------------
use strict;

use XML::LibXML;

sub new
{
	my ($class, %params) = @_;

	my $self = {
		scriptsdir => $params{'scriptsdir'} || undef,
	};
	if(! defined $self->{scriptsdir})
	{
		die("the scripts dir must be provided to use this module!\n");
	}
	$self->{'testlistxml'} = $self->{scriptsdir} . "/ccsm_utils/Testlistxml/testlist.xml";
	stat($self->{'testlistxml'}) or die "cannot find $self->{'testlistxml'}";
	
	bless ($self, $class);
	return ($self);
}

sub findTestsForCase
{
	my $self = shift;
	my ($args) = @_;
	my $compset = $$args{'compset'} if defined $$args{'compset'};
	my $grid = $$args{'grid'} if defined $$args{'compset'};
	my $machine = $$args{'machine'} if defined $$args{'machine'};
	my $compiler = $$args{'compiler'} if defined $$args{'compiler'};
	my $cesmtest = new CESMTest(compset => $compset, grid => $grid, 
							    machine => $machine, compiler => $compiler);

	my $testxml = $self->_readTestListXML();

	my $root = $testxml->getDocumentElement();
	my @machinesforcase;
	my @compilersforcase;
	my @testsforcase;

	# find the matching compset..
	foreach my $compsetnode($root->findnodes('/testlist/compset'))
	{
		my $compsetname = $compsetnode->getAttribute('name');
		# skip unless the compset names match
		if(defined $cesmtest->{compset})
		{
			next unless $cesmtest->{compset} eq $compsetname;
		}
	
		my @gridnodes = $compsetnode->findnodes('./grid');
		foreach my $gridnode($compsetnode->findnodes('./grid'))
		{
			my $gridname = $gridnode->getAttribute('name');

			# skip unless the grid names match
			if(defined $cesmtest->{grid})
			{
				next unless $cesmtest->{grid} eq $gridname;
			}
			
			foreach my $testnode($gridnode->findnodes('./test'))
			{
				my $testname = $testnode->getAttribute('name');
				my @machnodes = $testnode->findnodes('./machine');
				foreach my $machnode(@machnodes)
				{
					push(@machinesforcase, $machnode->textContent());
					push(@compilersforcase, $machnode->getAttribute('compiler'));
					push(@testsforcase, $machnode->getAttribute('testtype'));
				}
			}
		}
	}

	# Now, make sure the machine, compiler names are unique..
	my @compilers;
	my @machines;
	my @testtypes;
	my %uniqcompilers;
	my %uniqmachines;
	my %uniqtesttypes;

	map { $uniqcompilers{$_} = 1 } @compilersforcase;
	@compilers = sort keys %uniqcompilers;
	map { $uniqmachines{$_} = 1 } @machinesforcase;
	@machines = sort keys %uniqmachines;
	map { $uniqtesttypes{$_} = 1 } @testsforcase;
	@testtypes = sort keys %uniqtesttypes;
	
	if (! @compilers && ! @machines && ! @testtypes)
	{
		my $msg = <<END;
WARNING:: The following compset/grid combination $cesmtest->{compset}/$cesmtest->{grid} is NOT 
tested during the standard CESM development process. Thus you may likely find that this configuration
will NOT work, and are on your own to figure out how to get it working.
END
		return $msg;
	}
		
my $msg = <<END;
The compset $cesmtest->{compset} and grid $cesmtest->{grid} are tested on the following
machines, compilers, and/or test categories: 
END
	my $line = '';
	if(@machines)
	{
		$msg .= "Machines: ";
		$line = commify(@machines);

    	$msg .= "$line\n";
	}
	
	if(@compilers)
	{
		$msg .= "Compilers: ";
		$line = commify(@compilers);
    	$msg .= "$line\n";
	}
	
	if(@testtypes)
	{
		$msg .= "Test Types: ";
		$line = commify(@testtypes);
    	$msg .= "$line\n";
	}
$msg .= <<END;
The closer the tests are to the machine and compiler you are using and the more tests
that are done, the more likely your case will work without trouble.
END
	
	return $msg;
	
}

sub _readTestListXML
{
	my $self = shift;
	my $parser = XML::LibXML->new( no_blanks => 1);
	my $testxml = $parser->parse_file($self->{'testlistxml'});
	return $testxml;

}

sub commify
{
	(@_ == 0) ? ''                :
	(@_ == 1) ? $_[0]             : 
	(@_ == 2) ?  join(" and ", @_) : 
				join(", ", @_[0 .. ($#_-1)], "and $_[-1]");
}

sub main
{
	my %case;
	$case{'compset'} = 'BC5';
	$case{'grid'} = 'ne30_g16';
	my $msg = Testing::TestLists->findTestsForCase(\%case);
}
main(@ARGV) unless caller();
1;
