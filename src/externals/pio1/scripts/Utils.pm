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
    }elsif($host =~ /^eos/){
	$host = "eos";
    }elsif($host =~ /^titan/){
	$host = "titan";
    }elsif($host =~ /^ath/ ){
	$host = "athena";
    }elsif($host =~ /^kra/){
	$host = "kraken";
    }elsif($host =~ /^lynx/){
        $host = "lynx";
    }elsif($host =~ /^hopp/){
	$host = "hopper";
    }elsif($host =~ /^cvrs/) {
        $host = "carver";
    }elsif($host =~/erlogin/) {
	$host="erebus";
    }elsif($host =~/yslogin/) {
	$host="yellowstone";
    }elsif( $host =~ /^login/){
	if(-d "/lustre/janus_scratch"){
	    $host="janus";
        }else{
            $host = "athena";
	}
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
   if($host =~ "erebus" or $host =~ "yellowstone"){
       if(defined $ENV{ACCOUNT}){
	   $project=$ENV{ACCOUNT};
       }else{
	   $project="P93300606";
       }
       $projectInfo = "#BSUB -R \"span[ptile=16]\"\n#BSUB -P $project\n";
   }elsif($host =~ "titan"){
     $project = `showproj -s $host | tail -1`;
     $projectInfo ="#PBS -A $project\n";
   }elsif($host =~ "athena" or $host =~ "kraken"){
#    $project = `showproj -s athena | tail -1`;
     $projectInfo ="##PBS -A $project\n";
   }elsif($host =~ "columbia" or $host =~ "pleiades" or $host=~ "carver"){
     $project = "";
     $projectInfo ="##PBS -W group_list=$project\n";
   }
#HOST SPECIFIC END
}

sub preambleResource{
  my ($mod,$host,$pecount,$corespernode) = @_;
  my $nodes;
  my $preambleResource;
  if($host =~ "bluefire" or $host =~ "erebus" or $host =~ "yellowstone") {
     $preambleResource = "#BSUB -n $pecount\n";
  }elsif($host =~ "frost"){
     $preambleResource = "";
  }elsif($host =~ "edinburgh" or $host =~ "carver"){
     $nodes = ceil($pecount/$corespernode);
     $preambleResource = "#PBS -l nodes=$nodes:ppn=$corespernode\n";
  }elsif($host =~ "aum"){
     $nodes = ceil($pecount/$corespernode);
     $preambleResource = "#PBS -l nodes=$nodes:ppn=$corespernode\n";
  }elsif($host =~ "cyberstar" ){
     $nodes = ceil($pecount/$corespernode);
     $preambleResource = "#PBS -l nodes=$nodes:nehalem:ppn=$corespernode\n";
  }elsif($host =~ "lynx" or $host =~ "hopper"){
     my $pecnt = $corespernode*ceil($pecount/$corespernode);
     $preambleResource = "#PBS -l mppwidth=$pecnt\n";
  }elsif($host =~ "athena"  or $host =~ /janus/){
     my $pecnt = $corespernode*ceil($pecount/$corespernode);
     $preambleResource = "#PBS -l size=$pecnt\n";
  }elsif($host =~ "kraken"){
     my $pecnt = $corespernode*ceil($pecount/$corespernode);
     $preambleResource = "#PBS -l size=$pecnt\n";
  }elsif($host =~ "titan"){
     my $nodecnt = ceil($pecount/$corespernode);
     $preambleResource = "#PBS -l nodes=$nodecnt\n";
  }elsif($host =~ "columbia" or $host =~ "pleiades"){
     $preambleResource = "#PBS -l ncpus=$pecount\n";
  }
}

sub runString{
  my ($mod,$host,$pecount,$run,$exename,$log)=@_;
  my $runString;
  if($host =~ "bluefire" || $host =~ "erebus" || $host =~ "yellowstone") {
    $runString = "$run $exename 1> $log 2>&1";
  }elsif($host eq "frost" ) {
    $runString = "$run $log -np $pecount $exename";
    #$runString = "$run -np $pecount $exename";
  }elsif($host eq "intrepid") {
    $runString = "$run $log -np $pecount $exename ";
  }elsif($host =~ "columbia" or $host =~ "pleiades"){
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
   $submit =~ s/&lt;/</;
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
    my $modpath = {
		   erebus => "/usr/share/Modules/",
		   yellowstone => "/glade/apps/opt/modulefiles",
		   yellowstone_pgi => "/glade/apps/opt/modulefiles",
		   yellowstone_gnu => "/glade/apps/opt/modulefiles",
		   titan  => "/opt/modules/default/",
		   eos  => "/opt/modules/default/",
		   athena => "/opt/modules/default/",
		   kraken => "/opt/modules/default/",
		   hopper => "/opt/modules/default/",
		   lynx => "/opt/modules/default/",
		   lynx_intel => "/opt/modules/default/",
		   pleiades => "/usr",
                   carver => "/usr/common/nsg/opt/Modules/default/",
                   columbia => "/usr/share/modules/"};
#HOST SPECIFIC END

    return unless(defined $modpath->{$host});


#HOST SPECIFIC START
    if($host =~ "titan"){
	require "/opt/modules/default/init/perl";
	module_check($modpath,$host);
        module("switch cray-mpich2    cray-mpich2/5.6.3");
        module(" switch xt-libsci xt-libsci/12.0.00");
        module(" swap xt-asyncpe xt-asyncpe/5.16");
        module("load szip/2.1");
	module(" switch pgi pgi/12.10.0");
	module(" load netcdf-hdf5parallel/4.2.0");
	module(" load parallel-netcdf/1.3.1");
        module("list");
    }elsif($host eq "eos"){
	require "/opt/modules/default/init/perl";
	module_check($modpath,$host);
	module("rm netcdf");
	module("rm cray-netcdf");
	module("rm cray-netcdf-hdf5parallel");
	module("rm pnetcdf");
	module("switch cray-mpich cray-mpich/6.0.2");
	module("switch intel      intel/13.1.3.192");
	module("load cray-netcdf-hdf5parallel/4.3.0");
	module("load  cray-parallel-netcdf/1.3.1.1");
	module("load cmake/2.8.11.2");
    }elsif($host =~ "athena"){
#	require "/opt/modules/default/init/perl";
	module_check($modpath,$host);
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
	module_check($modpath,$host);
	module(" load netcdf/3.6.3");
	module(" load p-netcdf/1.2.0");
    }elsif($host =~ "hopper"){
	require "/opt/modules/default/init/perl";
	module_check($modpath,$host);
	module(" load cray-netcdf-hdf5parallel/4.3.0");
	module(" load cray-parallel-netcdf/1.3.1.1");
        module("list");
    }elsif($host =~ "pleiades"){
        module(" load netcdf/4.0-i10.1");
    }elsif($host =~ "columbia"){
        module(" load pd-netcdf.3.6.2");
!        module(" load pd-pnetcdf.1.1.1");
    }elsif($host =~ "lynx_intel"){
	require "/opt/modules/default/init/perl";
	module_check($modpath,$host);
	module(" rm PrgEnv-pgi ");
	module(" load PrgEnv-intel");
	module(" switch intel intel/12.1.0");
        module(" load INTEL/netcdf4/4.1.3_seq");
        module(" load pnetcdf/1.2.0");
	module(" list");
    }elsif($host =~ "hopper_gnu"){
	require "/opt/modules/default/init/perl";
	module_check($modpath,$host);
	module(" rm PrgEnv-pgi ");
	module(" load PrgEnv-gnu");
	module(" switch gcc gcc/4.7.1");
        module(" load netcdf-hdf5parallel/4.2.0");
        module(" load parallel-netcdf/1.2.0");
	module(" list");
    }elsif($host eq "carver"){
	require "/usr/common/nsg/opt/Modules/default/init/perl";
	module_check($modpath,$host);
        module("load intel ");
	module("load openmpi-intel/1.6");
	module("load netcdf-intel/4.1.1");
#        module("load pnetcdf/1.3.0");
        module("list");
    }elsif($host eq "erebus"){
	require "/glade/apps/opt/lmod/lmod/init/perl";
	module_check($modpath,$host);
        module("load intel/13.1.2");
	module("load ncarcompilers/1.0");
        module("rm netcdf");
        module("load netcdf-mpi/4.2");
#        module("load netcdf/4.2");
        module("load pnetcdf/1.3.0");
        module("load ncarenv/1.0");
	module("load ncarbinlibs/1.1");
	module("list");
    }elsif($host eq "yellowstone_pgi"){
	print "Loading modules for $host\n";
	require "/glade/apps/opt/lmod/lmod/init/perl";
	module_check($modpath,"yellowstone");
        module("rm netcdf");
        module("rm intel");
        module("load pgi/13.9");
	module("load ncarcompilers/1.0");
        module("unload netcdf");
        module("load netcdf/4.2");
        module("load pnetcdf/1.3.0");
        module("load ncarenv/1.0");
	module("load ncarbinlibs/1.1");
	module("list");
    }elsif($host eq "yellowstone_gnu"){
	print "Loading modules for $host\n";
	require "/glade/apps/opt/lmod/lmod/init/perl";
	module_check($modpath,"yellowstone");
        module("rm netcdf");
        module("rm intel");
        module("load gnu/4.8.0");
	module("load ncarcompilers/1.0");
        module("unload netcdf");
        module("load netcdf/4.3.0-rc4");
        module("load pnetcdf/1.3.0");
        module("load ncarenv/1.0");
	module("load ncarbinlibs/0.0");
	module("list");
    }elsif($host eq "yellowstone"){
	require "/glade/apps/opt/lmod/lmod/init/perl";
	module_check($modpath,$host);
#        module("purge");
        module("load intel/13.1.2");
	module("load ncarcompilers/1.0");
        module("rm netcdf");
        module("load netcdf-mpi/4.3.0");
        module("load pnetcdf/1.3.0");
        module("load ncarenv/1.0");
	module("load ncarbinlibs/1.0");
	module("load cmake");
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

sub module_check{
    my($modpath,$host) = @_;

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
    print "module path = $ENV{MODULEPATH}\n";
}

sub hostmods{
    my($self,$host,$mpi) = @_;
    my ($scratch,$netcdf,$pnetcdf,$cc,$fc,$filesystem);
    print "host = $host\n";
    if($host =~ /yellowstone/){
	$scratch = "/glade/scratch/$ENV{USER}/piotest";
	$netcdf = $ENV{NETCDF};
	$pnetcdf = $ENV{PNETCDF};
	$filesystem = "gpfs";
	if($mpi ne "mpi-serial") {
	    $cc = "mpicc";
	    $fc  = "mpif90";
	}
    }
    if($host =~ /^eos/){
	$scratch = "$ENV{WORKDIR}/$ENV{PROJECT}/testpio";
	$netcdf = $ENV{NETCDF_DIR};
	$pnetcdf = $ENV{PARALLEL_NETCDF_DIR};
	$filesystem = "lustre";
	$cc = "cc";
	$fc = "ftn";
    }
    return ($scratch,$netcdf,$pnetcdf,$cc,$fc,$filesystem);
}

1;
