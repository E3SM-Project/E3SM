#!/usr/bin/perl
use strict;
use Cwd;
use Getopt::Long;

my $suites;
my $retry=0;
my $help=0;
my $host;
my $result = GetOptions("suites=s@"=>\$suites,"retry"=>\$retry,"host=s"=>\$host,"help"=>\$help);

usage() if($help);
sub usage{
    print "--suites : Test only the listed suites (all, snet, pnet, mpiio, ant, bench)\n";
    print "--retry  : Do not repeat tests that have already passed\n";
    print "--host   : Force a hostname for testing\n";
    print "--help   : Print this message\n";
    exit;
}




my $cfgdir = `pwd`;
chomp $cfgdir;
my $clean = 'yes';
my @valid_env = qw(NETCDF_PATH PNETCDF_PATH MPI_LIB MPI_INC F90 FC CC FFLAGS
                   MPICC MPIF90 LDLIBS);


my @testsuites = qw(all snet pnet mpiio ant bench);



# The XML::Lite module is required to parse the XML configuration files.
(-f "$cfgdir/perl5lib/XML/Lite.pm")  or  die <<"EOF";
** Cannot find perl module \"XML/Lite.pm\" in directory \"$cfgdir/perl5lib\" **
EOF

unshift @INC, "$cfgdir/perl5lib";
require XML::Lite;
require Utils;

$host = Utils->host() unless(defined $host);
Utils->loadmodules("$host");

print "host = $host\n";

my $xml = XML::Lite->new( "build_defaults.xml" );

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

my $srcdir = "$workdir/src";
my $tstdir = "$srcdir/testpio";
my $testpiodir = cwd();
my $piodir = "$testpiodir/..";
my $date = `date +%y%m%d-%H%M%S`;
my $user = $ENV{USER};
chomp $date;

my $outfile = "$testpiodir/testpio.out.$date";
my $script  = "$testpiodir/testpio.sub.$date";

open(F,">$script");
print F "#!/usr/bin/perl\n";
print F "$attributes{preamble}\n";

# get a valid project number for this user
my $project;
if($host eq "bluefire" or $host eq "frost"){
    open(G,"/etc/project.ncar");
    foreach(<G>){
	if($_ =~ /^$user:(\d+),?/){
	    $project = $1;
	    last;
	}
    }
    close(G);
}elsif($host eq "jaguar"){
    $project = `/sw/xt5/bin/showproj -s jaguar | tail -1`;
    print F "#PBS -A $project\n";
}elsif($host eq "athena"){
#    $project = `showproj -s athena | tail -1`;
    print F "##PBS -A $project\n";
}
if($host eq "bluefire"){
    print F "#BSUB -R \"span[ptile=64]\"\n";
    print F "#BSUB -P $project\n";
}

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




    




print F << "EOF";
use strict;
use File::Copy;

chdir ("$cfgdir");

mkdir "$srcdir" if(! -d "$srcdir");

my \$rc = 0xffff & system("rsync -rp $piodir $srcdir");
if(\$rc != 0) {
    system("cp -fr $piodir/pio $srcdir");
    system("cp -fr $piodir/mct $srcdir");
    system("cp -fr $piodir/timing $srcdir");
    system("cp -fr $piodir/testpio $srcdir");
}

my \$confopts = {all=>"--enable-pnetcdf --enable-mpiio --enable-netcdf --enable-timing",
		snet=>"--disable-pnetcdf --disable-mpiio --enable-netcdf --enable-timing",
		pnet=>"--enable-pnetcdf --disable-mpiio --disable-netcdf --enable-timing",
		ant=>"--enable-pnetcdf --enable-mpiio --enable-netcdf --disable-timing",
		mpiio=>"--disable-pnetcdf --enable-mpiio --disable-netcdf --enable-timing",
		bench=>"--enable-pnetcdf --enable-mpiio --enable-netcdf --enable-timing"
	    };

my \$testlist = {all=>["sn01","sn02","sn03","sb01","sb02","sb03","sb04","sb05","sb06","sb07","sb08","pn01",
                      "pn02","pn03","pb01","pb02","pb03","pb04","pb05","pb06","pb07","pb08","bn01","bn02",
                      "bn03","bb01","bb02","bb03","bb04","bb05","bb06","bb07","bb08","wr01","rd01"],
		snet=>["sn01","sn02","sn03","sb01","sb02","sb03","sb04","sb05","sb06","sb07","sb08"],
		pnet=>["pn01","pn02","pn03","pb01","pb02","pb03","pb04","pb05","pb06","pb07","pb08"],
		ant=>["sn02","sb02","pn02","pb02","bn02","bb02"],
		mpiio=>["bn01","bn02","bn03","bb01","bb02","bb03","bb04","bb05","bb06","bb07","bb08"],
		bench=>["b14"]
	    };

unlink("$workdir/wr01.dof.txt") if(-e "$workdir/wr01.dof.txt");
my \$suite;
my \$passcnt=0;
my \$failcnt=0;


foreach \$suite (qw(@testsuites)){
    my \$confopts = \$confopts->{\$suite};
    my \@testlist = \@{\$testlist->{\$suite}};
    chdir ("$tstdir");
    unlink("../pio/Makefile.conf");
    system("perl ./testpio_build.pl --conopts=\\"\$confopts\\" --host=$host");
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
	    symlink("$tstdir/namelists/testpio_in.\$test","testpio_in");
	    mkdir "none" unless(-d "none");
	    my \$log = "testpio.out.$date";
	    
	    if("$host" eq "frost"){
		open(JC, "\$run -p $project -o \$log ./testpio < testpio_in |");
		my \$rv = <JC>;
		close(JC);
		chomp \$rv;
		my \$rv2 = system("cqwait \$rv");
		print "Job \$rv complete\\n";
	    }else{
		system("\$run ./testpio 1> \$log 2>&1");
	    }
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
my $submit = $attributes{submit};
exec("$submit $script");
