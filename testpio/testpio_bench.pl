#!/usr/bin/perl
use strict;
use Cwd;
use Getopt::Long;
use File::Copy;

my $preambleResource;
my $projectInfo;
my $project;
my $nodecount;
my $suites;
my $retry=0;
my $help=0;
my $host;
my $pecount;
my $bname;
my $iofmt;
my $rearr;
my $numIO;
my $stride;
my $maxiter;
my $dir;
my $debug=0;
my $numagg;
my $numvars;
my $iodecomp;
my $logfile = 'testpio.out';
my $result = GetOptions("suites=s@"=>\$suites,
                        "retry"=>\$retry,
                        "host=s"=>\$host,
			"pecount=i"=>\$pecount,
			"bench=s"=>\$bname,
                        "iofmt=s"=>\$iofmt,
                        "rearr=s"=>\$rearr,
			"numIO|numIOtasks=i"=>\$numIO,
                        "stride=i"=>\$stride,
 			"maxiter=i"=>\$maxiter,
                        "dir=s"=>\$dir,
                        "numagg=i"=>\$numagg,
			"numvars=i"=>\$numvars,
                        "decomp=s"=>\$iodecomp,
			"debug"=>\$debug,
		        "log=s"=>\$logfile,
		        "help"=>\$help);

usage() if($help);
sub usage{
    print "--suites                : Test only the listed suites (all, snet, pnet, mpiio, ant, bench)\n";
    print "--retry                 : Do not repeat tests that have already passed\n";
    print "--host                  : Force a hostname for testing\n";
    print "--pecount               : Select the processor count on which to run benchmark (defined in config_bench.xml) \n";
    print "--bench                 : Select the name of the benchmark to run (defined in config_bench.xml)\n";
    print "--iofmt                 : Selects the type of file to write (pnc,snc,bin)\n";
    print "--rearr                 : Selects the type of rearrangement (box,mct,none)\n";
    print "--numIOtasks (--numIO)  : Sets the number of IO tasks used by PIO\n";
    print "--stride                : Sets the stride between IO tasks, Note this is ignored on Blue Gene\n";
    print "--maxiter               : Sets the number of files to write\n";
    print "--dir                   : Sets the subdirectory for which to write files \n";
    print "--numagg                : Sets the number of MPI-IO aggregators to use \n";
    print "--numvars               : Sets the number of variables to write to each file \n";
    print "--decomp                : Sets the form of the IO-decomposition (x,y,z,xy,xye,xz,xze,yz,yze,xyz,xyze,setblk,cont1d,cont1dm)\n";
    print "--help                  : Print this message\n";
    print "--debug		   : Generate the runscript but do not submit it\n";
    print "--log                   : Manually sets the log output file for the benchmark\n";
    exit;
}

my $cfgdir = `pwd`;
chomp $cfgdir;
my $clean = 'yes';
my @valid_env = qw(NETCDF_PATH PNETCDF_PATH MPI_LIB MPI_INC F90 FC CC ALLCFLAGS FFLAGS
                   MPICC MPIF90 LDLIBS);

my @testsuites = qw(bench);

# The XML::Lite module is required to parse the XML configuration files.
(-f "$cfgdir/../testpio/perl5lib/XML/Lite.pm")  or  die <<"EOF";
** Cannot find perl module \"XML/Lite.pm\" in directory \"$cfgdir/../testpio/perl5lib\" **
EOF

unshift @INC, "$cfgdir/../testpio","$cfgdir/../testpio/perl5lib";
require XML::Lite;
require Utils;

$host = Utils->host() unless(defined $host);
Utils->loadmodules("$host");

my $xml = XML::Lite->new( "$cfgdir/../testpio/build_defaults.xml" );

my $root = $xml->root_element();
my $settings = $xml->elements_by_name($host);
my %attributes = $settings->get_attributes;


foreach(keys %attributes){
    if($attributes{$_} =~  /\$\{?(\w+)\}?/){
	my $envvar = $ENV{$1};
	$attributes{$_}=~ s/\$\{?$1\}?/$envvar/
    }
#    if(/ADDENV_(.*)/){
#	print F "\$ENV{$1}=\"$attributes{$_}:\$ENV{$1}\n\"";
#    }elsif(/ENV_(.*)/){
#        print "set $1 $attributes{$_}\n";
#	print F "\$ENV{$1}=\"$attributes{$_}\n\"";
#    }	
    
}

if(defined $suites){
    @testsuites = @$suites;
}elsif(defined $attributes{testsuites}){
    @testsuites = split(' ',$attributes{testsuites});
}



my $workdir = $attributes{workdir};

if(-d $workdir){
    print "Using existing directory $workdir\n";
}else{
    print "Creating directory $workdir\n";
    mkdir $workdir or die "Could not create directory"
}

my $config = XML::Lite->new("$cfgdir/../testpio/config_bench.xml");
my $elm = $config->root_element();
print "pecount is $pecount\n";

my $ldx=0;
my $ldy=0;
my $ldz=0;
my $nx_global=0;
my $ny_global=0;
my $nz_global=0;

my %configuration = ( ldx => 0,
                      ldy => 0,
                      ldz => 0,
                      iofmt => 'pnc', 
                      rearr => 'box',
                      numprocsIO => -10,
   		      stride => -1,
                      maxiter => 10,
                      dir => './none/',
		      iodecomp => 'yze',
                      numagg => -1,
                      numvars => 10);

#-------------------------------------------------
# Modify the configuration based on arguments
#-------------------------------------------------
if (defined $iofmt)   {$configuration{'iofmt'} = $iofmt;}
if (defined $rearr)   {$configuration{'rearr'} = $rearr;}
if (defined $numIO)   {$configuration{'numprocsIO'} = $numIO;}
if (defined $stride)  {$configuration{'stride'} = $stride;}
if (defined $maxiter) {$configuration{'maxiter'} = $maxiter;}
if (defined $numvars) {$configuration{'numvars'} = $numvars;}
if (defined $dir)     {$configuration{'dir'} = $dir;}
if (defined $numagg)  {$configuration{'numagg'} = $numagg;}

# See if you can find the general benchmark description first
my @blist = $config->elements_by_name("BenchConfig");
my $bchildren = $elm->get_children();
my $found=0;
foreach my $child (@blist) {
  my %atts = $child->get_attributes;
  my $bn = $atts{"bench_name"};
  if ($bn =~ $bname) {
     my $nx_global = $atts{"nx_global"}; $configuration{'nx_global'} = $nx_global;
     my $ny_global = $atts{"ny_global"}; $configuration{'ny_global'} = $ny_global;
     my $nz_global = $atts{"nz_global"}; $configuration{'nz_global'} = $nz_global;
     $found = 1;
  }
}
if(!$found) {
  printf "Could not find configuration for benchmark: %s\n" ,$bname;
  exit(-1);
} else {
  print "nx_global: $configuration{'nx_global'} ny_global: $configuration{'ny_global'} nz_global: $configuration{'nz_global'}\n";
}

my $root = "CompConfig";
my @list = $config->elements_by_name($root);
my $children = $elm->get_children();

my $found=0;
foreach my $child (@list ) {
  my %atts = $child->get_attributes;
  my $name = $child->get_name();
  my @keys = keys(%atts);
  my $np = $atts{"nprocs"};
  my $bn = $atts{"bench_name"};
#  printf "bench_name is $bn\n";
#  printf "bname is $bname\n";
  if(($np eq $pecount) & ($bn =~ $bname)) {
     my @gchildren = $child->get_children();
     foreach my $grandchild (@gchildren) {
	   my $name  = $grandchild->get_name();
           my $value = $grandchild->get_text();
	   $configuration{$name}=$value;
     }
     $found = 1;
  } 
}
my $suffix = $bname . "-" . $pecount;
my $testname = "bench." . $suffix;
printf "testname: %s\n",$testname;

if (defined $iodecomp) {$configuration{'iodecomp'} = $iodecomp;}

#print "ldx: $configuration{'ldx'} ldy: $ldy ldz: $ldz\n";
my $outfile;
if ($found) {
  print "ldx: $configuration{'ldx'} ldy: $configuration{'ldy'} ldz: $configuration{'ldz'}\n";
  $outfile = "$cfgdir/testpio_in." . $suffix;
  unlink("$outfile") if(-e "$outfile");

  open(F,"+> $outfile");
  gen_io_nml();       # Generate the io_nml namelist
  gen_compdof_nml();  # Generate the compdof_nml namelist
  if(defined $iodecomp) {
      gen_iodof_nml();    # Generate the iodof_nml namelist
  }
  gen_prof_inparm();  # Generate the prof_inparm namelist

  close(F);
} else {
  printf "Could not find configuration for benchmark: %s on %d MPI tasks \n" ,$bname, $pecount;
  exit(-1);
}

my $corespernode = $attributes{corespernode};

my $srcdir = "$workdir/src";
my $tstdir = "$srcdir/testpio";
copy("$outfile","$tstdir");
my $testpiodir = cwd();
my $piodir = "$testpiodir/..";
my $date = `date +%y%m%d-%H%M%S`;
my $user = $ENV{USER};
chomp $date;

my $outfile = "$testpiodir/testpio.out.$date";
my $script  = "$testpiodir/testpio.sub.$date";

open(F,">$script");
print F "#!/usr/bin/perl\n";
$preambleResource = Utils->preambleResource("$host","$pecount","$corespernode");
print F $preambleResource;
print F "$attributes{preamble}\n";

# Create a valid project string for this user
$projectInfo = Utils->projectInfo("$host","$user");
print F $projectInfo;

my @env;
foreach(keys %attributes){
#    if($attributes{$_} =~  /\$\{?(\w+)\}?/){
#	my $envvar = $ENV{$1};
#	$attributes{$_}=~ s/\$\{?$1\}?/$envvar/
#    }
    if(/ADDENV_(.*)/){
	print F "\$ENV{$1}=\"$attributes{$_}:\$ENV{$1}\"\;\n";
    }elsif(/ENV_(.*)/){
        print "set $1 $attributes{$_}\n";
	print F "\$ENV{$1}=\"$attributes{$_}\"\;\n";
    }	
    
}

my $log;
my $setlogfile;
my $run = $attributes{run};
my $exename = "./testpio";
if($logfile) {
  $log = "$logfile";
  $setlogfile = 1
} else {
  $log = "$cfgdir/testpio.log.";
  $setlogfile = 2
}
my $foo= Utils->runString($host,$pecount,$run,$exename,$log);

print "EXEC command: ($foo)\n";

print F << "EOF";
use strict;
use lib "$cfgdir";
use File::Copy;
use POSIX qw(ceil);
#unshift \@INC, "$cfgdir/../testpio";

#chmod 0755,"$cfgdir/../testpio/Utils.pm";
use Utils;

chdir ("$cfgdir");

mkdir "$srcdir" if(! -d "$srcdir");

my \$rc = 0xffff & system("rsync -rp $piodir $srcdir");
if(\$rc != 0) {
    system("cp -fr $piodir/pio $srcdir");
    system("cp -fr $piodir/mct $srcdir");
    system("cp -fr $piodir/timing $srcdir");
    system("cp -fr $piodir/testpio $srcdir");
}

#my \$confopts = {bench=>"--enable-pnetcdf --enable-mpiio --enable-netcdf --enable-timing"};
my \$confopts = {bench=>""};

my \$testlist = {bench=>["generated"]};

unlink("$workdir/wr01.dof.txt") if(-e "$workdir/wr01.dof.txt");
my \$suite;
my \$passcnt=0;
my \$failcnt=0;
my \$host   = "$host";
my \$pecount = $pecount;
my \$run     = "$attributes{run}";

foreach \$suite (qw(@testsuites)){
    my \$confopts = \$confopts->{\$suite};
#    my \@testlist = \@{\$testlist->{\$suite}};
    my \@testlist = \"$suffix";
#    unlink("../pio/Makefile.conf");
#    May want to uncomment this system call so that testpio gets build 
#    system("perl ./testpio_build.pl --conopts=\\"\$confopts\\" --host=$host");
    copy("testpio","$tstdir");     # copy executable into test directory
#    copy("testpio_in","$tstdir"); # copy the namelist file into test directory
    
    chdir ("$tstdir");
    my \$test;
    my \$run = "$attributes{run}";
    if(-e "../pio/Makefile.conf" && -e "testpio"){
	foreach \$test (\@testlist){
	    my \$casedir = "$workdir/\$suite.\$test";
	    mkdir \$casedir unless(-d \$casedir);
	    chdir(\$casedir) or die "Could not cd to \$casedir";
	    print "\$suite \$test    ";
	    if($retry && -e "TestStatus"){
		open(T,"TestStatus");
		my \$result = <T>;
		close(T);
		if(\$result =~ /PASS/){
		    \$passcnt++;
		    print "Test already PASSED\\n";
		    next;
		}
	    }

	    unlink("testpio") if(-e "testpio");
	    copy("$tstdir/testpio","testpio");
	    chmod 0755,"testpio";
#	    symlink("$tstdir/namelists/testpio_in.\$test","testpio_in");
#	    symlink("$tstdir/testpio_in.\$test","testpio_in");
	    symlink("$tstdir/testpio_in.\$test","testpio_in");
	    mkdir "none" unless(-d "none");
	    my \$exename = "./testpio";
            my \$log;
            if($setlogfile) {
                \$log = "$logfile";
		unlink("\$log") if(-e "\$log");
            } else {
	        \$log = "\$casedir/testpio.out.$date";
            }
            my \$sysstr =  Utils->runString(\$host,\$pecount,\$run,\$exename,\$log);
	    system(\$sysstr);
	    open(LOG,\$log);
	    my \@logout = <LOG>;
	    close(LOG);
	    
	    my \$cnt = grep /testpio completed successfully/ , \@logout;
            open(T,">TestStatus");
	    if(\$cnt>0){
		\$passcnt++;
		print "PASS \\n";
		print T "PASS \\n";
	    }else{
		\$failcnt++;
		print "FAIL \\n";
		print T "FAIL \\n";
	    }
	    close(T);
	}
    }else{
	print "suite \$suite FAILED to configure or build\\n";	
    }
}
print "test complete on $host \$passcnt tests PASS, \$failcnt tests FAIL\\n";
EOF
close(F);
chmod 0755, $script;
my $subsys = Utils->submitString($host,$pecount,$corespernode,$attributes{submit},$script);
print "submit: ($subsys)\n";
if($debug) {
    print "Created script ($script)\n";
} else {
#    exec("$subsys");
    my @foo2 = `$subsys`;
    my $jobid;
    foreach my $i (@foo2) {
     ($jobid) = ($i =~/([0-9]+)/);
#     print "jobid: ($jobid)\n";
    }
    exec("cqwait $jobid");
}


sub gen_compdof_nml{
  print F "&compdof_nml\n";
  print F "grddecomp = 'setblk'\n";
  print F "gdx = $configuration{'ldx'}\n";
  print F "gdy = $configuration{'ldy'}\n";
  print F "gdz = $configuration{'ldz'}\n";
  print F "/\n";
}
sub gen_iodof_nml {
  print F "&iodof_nml\n";
  print F "grddecomp = '$configuration{'iodecomp'}'\n";
  print F "/\n";
}
sub gen_prof_inparm {
  print F "&prof_inparm\n";
  print F "profile_disable = .false.\n";
  print F "profile_barrier = .true.\n";
  print F "profile_single_file = .false.\n";
  print F "profile_depth_limit = 10\n";
  print F "profile_detail_limit = 0\n";
  print F "/\n";
}
sub gen_io_nml {
  print F "&io_nml\n";
  print F "casename       = '$suffix'\n";
  print F "nx_global      = $configuration{'nx_global'}\n";
  print F "ny_global      = $configuration{'ny_global'}\n";
  print F "nz_global      = $configuration{'nz_global'}\n";
  print F "nvars        = $configuration{'numvars'}\n";
  print F "iofmt          = '$configuration{'iofmt'}'\n";
  print F "rearr          = '$configuration{'rearr'}'\n";
  print F "nprocsIO       = $configuration{'numprocsIO'}\n";
  print F "stride         = $configuration{'stride'}\n";
  print F "maxiter        = $configuration{'maxiter'}\n";
  print F "dir            = '$configuration{'dir'}'\n";
  print F "num_aggregator = $configuration{'numagg'}\n";
  print F "DebugLevel  = 0\n";
  print F "compdof_input = 'namelist'\n";
  if(defined $iodecomp) {
     print F "iodof_input = 'namelist'\n";
  }
  print F "compdof_output = 'none'\n";
  print F "/\n";
}
