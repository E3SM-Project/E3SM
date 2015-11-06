package Config::envBatch;
my  $pkg_nm = __PACKAGE__;
use strict;
use warnings;
use XML::LibXML;
use Log::Log4perl qw(get_logger);

my $logger;

BEGIN{
    $logger = get_logger();
}

sub new {
     my $class = shift();
     my $this = {};
     bless($this, $class);
     $this->_init(@_);
     return $this;
}

sub _init {
  my ($this) = @_;
#  $this->SUPER::_init();

}

sub read {
    my ($self, $file) = @_;

    my $xml = XML::LibXML->new( no_blanks=>1)->parse_file($file);

    my @jobs = $xml->findnodes(".//job");
    
    foreach my $job (@jobs){
	my $name = $job->getAttribute('name');
	foreach my $entry ($job->childNodes()){
	    if($entry->nodeName() eq "entry"){
		my $id = $entry->getAttribute('id');
		my $value = $entry->getAttribute('value');
		$logger->debug ("job= $name id=$id value=$value");
		$self->{$name}{$id}=$value;
	    }
	}
    }
    
}

sub get
{
    my ($self, $job) = @_;
    
    return($self->{$job});
}


sub set
{
    my($self, $id, $job, $value)  = @_;
  
    if(defined $value){
  
    $self->{$job}{$id} = $value;
   
    qx(./xmlchange -noecho -file env_batch.xml -id $id -val $value -subgroup $job);
  }
}



1;
