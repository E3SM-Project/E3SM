package Utils;

use strict;
BEGIN {
        use vars       qw( $VERSION @ISA );
        $VERSION = '0.10';
        @ISA         = qw();
} # end BEGIN
# non-exported package globals go here
use vars      qw();


sub host{
    my $host = `hostname`;

    if($host =~ /^fr\d+en/){
	$host = "frost";
    }elsif($host =~ /^be\d+en/){
	$host = "bluefire";
    }elsif($host =~ /^ja/ or $host =~ /^yo/){
	$host = "jaguar";
    }elsif($host =~ /^ath/ or $host =~ /^log/){
	$host = "athena";
    }elsif($host =~ /(\w+)\./){
	$host = $1;
    }
}

sub loadmodules{
    my ($mod,$host) = @_;
    print "host = $host\n";

    my $modpath = {bluefire => "/contrib/Modules/3.2.6/",
		   jaguar  => "/opt/modules/default/",
		   athena => "/opt/modules/default/"};

    return unless(defined $modpath->{$host});

    $ENV{MODULESHOME} = $modpath->{$host};

    if($modpath->{$host} =~ /([^\/]*)\/?$/){
	$ENV{MODULE_VERSION}=$1;
    }
    if (! defined $ENV{MODULEPATH} ) {
	open(F,"$modpath->{$host}/init/.modulespath") || die "could not open $modpath->{$host}/init/.modulespath";
	my @file = <F>;
	close(F);
	my $modulepath;
	foreach(@file){
	    if(/^([\/\w+]+)\s*/){
		if(defined $modulepath){
		    $modulepath = "$modulepath:$1";
		}else{
		    $modulepath = $1;
		}
	    }
	}
	$ENV{MODULEPATH} = $modulepath;
	}

    if (! defined $ENV{"LOADEDMODULES"} ) {
	$ENV{"LOADEDMODULES"} = "";
    }


    
    if($host eq "bluefire"){
#	module("load xlf12");
#        module("list");
    }elsif($host eq "jaguar"){
#	require "/opt/modules/default/init/perl";
	module(" purge");
	module(" load PrgEnv-pgi Base-opts");
	module(" load xtpe-barcelona");
	module(" load torque moab");
	module(" switch pgi pgi/7.1.6");
	module(" load netcdf/3.6.2");      
	module(" load p-netcdf/1.0.3");
	module(" swap xt-asyncpe xt-asyncpe/1.0c");
	module(" load xt-binutils-quadcore/2.0.1");
        module("list");
    }elsif($host eq "athena"){
#	require "/opt/modules/default/init/perl";
	module(" purge");
	module(" load PrgEnv-pgi Base-opts");
	module(" load xtpe-quadcore");
	module(" load torque moab");
        module(" load xt-mpt");
	module(" switch pgi pgi/7.1.6");
	module(" load netcdf/3.6.2");      
	module(" load p-netcdf/1.0.3");
	module(" swap xt-asyncpe xt-asyncpe/1.0c");
	module(" swap xt-binutils-quadcore xt-binutils-quadcore/2.0.1");
    }
}


sub module {
    my $exec_prefix = "$ENV{MODULESHOME}";

    if(-e "$exec_prefix/bin/modulecmd"){
	eval `$exec_prefix/bin/modulecmd perl @_`;
    }else{
	die "Could not find $exec_prefix/bin/modulecmd";
    }
}




1;
