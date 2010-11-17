package Utils;

use strict;
BEGIN {
        use vars       qw( $VERSION @ISA );
        $VERSION = '0.10';
        @ISA         = qw();
} # end BEGIN
# non-exported package globals go here
use vars      qw();
use POSIX qw(ceil);

sub host{
    my $host = `hostname -f`;
    $host = `hostname` if($?);
#HOST SPECIFIC START
    if($host =~ /intrepid/){
	$host = "intrepid";
    }elsif($host =~ /^fr\d+en/){
	$host = "frost";
    }elsif($host =~ /^be\d+en/){
	$host = "bluefire";
    }elsif($host =~ /^ja/ or $host =~ /^yo/){
	$host = "jaguar";
    }elsif($host =~ /^ath/ or $host =~ /^log/){
	$host = "athena";
    }elsif($host =~ /^kra/){
	$host = "kraken";
    }elsif($host =~ /^lynx/){
        $host = "lynx";
    }elsif($host =~ /(\w+)\./){
	$host = $1;
    }
#HOST SPECIFIC END
}

sub projectInfo{
   my ($mod,$host,$user) = @_;
   my $projectInfo;
   my $project;
#HOST SPECIFIC START
   if($host =~ "bluefire" or $host =~ "frost"){
      open(G,"/etc/project.ncar");
      foreach(<G>){
         if($_ =~ /^$user:(\d+),?/){
            $project = $1;
            last;
         }
      }
      close(G);
      if($host =~ "bluefire") {
        $projectInfo = "#BSUB -R \"span[ptile=64]\"\n#BSUB -P $project\n";
      }
   }elsif($host =~ "jaguar"){
     $project = `/sw/xt5/bin/showproj -s jaguar | tail -1`;
     $projectInfo ="#PBS -A $project\n";
   }elsif($host =~ "athena" or $host =~ "kraken"){
#    $project = `showproj -s athena | tail -1`;
     $projectInfo ="##PBS -A $project\n";
   }elsif($host =~ "columbia"){
     $project = "";
     $projectInfo ="##PBS -W group_list=$project\n";
   }
#HOST SPECIFIC END
}

sub preambleResource{
  my ($mod,$host,$pecount,$corespernode) = @_;
  my $nodes;
  my $preambleResource;
  if($host =~ "bluefire") {
     $preambleResource = "#BSUB -n $pecount\n";
  }elsif($host =~ "frost"){
     $preambleResource = "";
  }elsif($host =~ "edinburgh"){
     $nodes = ceil($pecount/$corespernode);
     $preambleResource = "#PBS -l nodes=$nodes:ppn=$corespernode\n"; 
  }elsif($host =~ "aum"){
     $nodes = ceil($pecount/$corespernode);
     $preambleResource = "#PBS -l nodes=$nodes:ppn=$corespernode\n"; 
  }elsif($host =~ "cyberstar" ){
     $nodes = ceil($pecount/$corespernode);
     $preambleResource = "#PBS -l nodes=$nodes:nehalem:ppn=$corespernode\n"; 
  }elsif($host =~ "lynx"){
     my $pecnt = $corespernode*ceil($pecount/$corespernode);
     $preambleResource = "#PBS -l mppwidth=$pecnt\n"; 
  }elsif($host =~ "athena"){
     my $pecnt = $corespernode*ceil($pecount/$corespernode);
     $preambleResource = "#PBS -l size=$pecnt\n"; 
  }elsif($host =~ "jaguar" or $host =~ "kraken"){
     my $pecnt = $corespernode*ceil($pecount/$corespernode);
     $preambleResource = "#PBS -l size=$pecnt\n"; 
  }elsif($host =~ "columbia"){
     $preambleResource = "#PBS -l ncpus=$pecount\n"; 
  }
}

sub runString{
  my ($mod,$host,$pecount,$run,$exename,$log)=@_;
  my $runString;
  if($host =~ "bluefire") {
    $runString = "$run $exename 1> $log 2>&1";
  }elsif($host eq "frost" ) {
    $runString = "$run $log -np $pecount $exename";
    #$runString = "$run -np $pecount $exename";
  }elsif($host eq "intrepid") {
    $runString = "$run $log -np $pecount $exename ";
  }elsif($host =~ "columbia"){
    $runString = "$run -np $pecount $exename 1> $log 2>&1";
#  } elsif($host =~ "kraken" or $host =~ "jaguar" or $host =~ "athena"){
# make this default
  }else{
   $runString = "$run -n $pecount $exename 1> $log 2>&1";
  }

}

sub submitString{
   my ($mod,$host,$pecount,$corespernode,$submit,$script)=@_;
   my $submitString;
   my $nodecnt;
   if($host =~ "frost"){
      $nodecnt=ceil($pecount/$corespernode);
      $submitString = "$submit -n $nodecnt $script";
   }else{
      $submitString = "$submit $script";
   }
}

sub loadmodules{
    my ($mod,$host) = @_;

#HOST SPECIFIC START
    my $modpath = {bluefire => "/contrib/Modules/3.2.6/",
		   jaguar  => "/opt/modules/default/",
		   athena => "/opt/modules/default/",
		   kraken => "/opt/modules/default/",
		   lynx => "/opt/modules/default/",
                   columbia => "/usr/share/modules/"};
#HOST SPECIFIC END

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


    
#HOST SPECIFIC START
    if($host =~ "bluefire"){
#	module("load xlf12");
#        module("list");
    }elsif($host =~ "jaguar"){
#	require "/opt/modules/default/init/perl";
#	module(" purge");
#        module(" load xt-mpt/4.0.0");
	module(" load PrgEnv-pgi Base-opts");
	module(" load xtpe-istanbul");
	module(" load torque moab");
	module(" switch pgi pgi/9.0.4");
	module(" load netcdf/3.6.2");      
	module(" load p-netcdf/1.1.1");
#	module(" swap xt-asyncpe xt-asyncpe/1.0c");
#	module(" load xt-binutils-quadcore/2.0.1");
        
        module("list");
    }elsif($host =~ "athena"){
#	require "/opt/modules/default/init/perl";
	module(" purge");
	module(" load PrgEnv-pgi Base-opts");
	module(" load xtpe-quadcore");
	module(" load torque moab");
        module(" load xt-mpt");
	module(" switch pgi pgi/7.1.6");
	module(" load netcdf/3.6.2");      
	module(" load p-netcdf/1.1.1");
	module(" swap xt-asyncpe xt-asyncpe/1.0c");
	module(" swap xt-binutils-quadcore xt-binutils-quadcore/2.0.1");
    }elsif($host =~ "kraken"){
	require "/opt/modules/default/init/perl";
	module(" load netcdf/3.6.2");      
	module(" load p-netcdf/1.1.1");
    }elsif($host =~ "columbia"){
        module(" load pd-netcdf.3.6.2");
!        module(" load pd-pnetcdf.1.1.1");
    }elsif($host =~ "lynx"){
	require "/opt/modules/default/init/perl";
	module(" load netcdf");
        module("list");
    }
	
#HOST SPECIFIC END
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
