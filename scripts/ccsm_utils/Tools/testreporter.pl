#!/usr/bin/env perl
use Getopt::Long;
use Data::Dumper;
use LWP;
use HTTP::Request;
use XML::LibXML;
#-------------------------------------------------------------------------------
# testreporter.pl
# Perl script that watches the CESM tests as they progress, and sends the reports to 
# the testdb application at csegweb.cgd.ucar.edu

#-------------------------------------------------------------------------------
# Constants. Filenames we look in for test results, and the time we sleep between
# reporting test results.  
#-------------------------------------------------------------------------------
my $teststatusfilename = "TestStatus";
my $teststatusoutfilename = "TestStatus.out";
my $iopstatusfilename = "TestStatus.IOP";
my $casestatusfilename = "CaseStatus";
my $sleeptime = 120;
my $baselinetag;
# The URL we send test results to. 
my $posturl = "https://csegweb.cgd.ucar.edu/testdb/cgi-bin/processXMLtest.cgi";

# Options and global variables.
#-------------------------------------------------------------------------------
# root of the test suite currently running. 
my $testroot = undef;
# The tag name you are testing. 
my $tagname = undef;
# the testid parameter specified for ./create_test_suite.  This script uses this parameter to 
# find the tests you are running.  
my $testid = undef;
# the email address to send test reports to. 
my $email = undef;
# the test type: should be one of prealpha, prebeta, or prerelease. 
my $testtype = undef;
my $debug = 0;
my $dumpxml = 0;
my $printreport = 0;
my $username = undef;
my $password = undef;

#-------------------------------------------------------------------------------
# Main
# Get the options first. 
# Then, get the test directories, get the test suite info, get the test status for all the tests, 
# and send the results. 
#-------------------------------------------------------------------------------
opts();
authenticate();

my @testdirs;
my %suiteinfo;
while(1)
{
    @testdirs = &getTestDirs($testroot, $testid);
    %suiteinfo = &getTestSuiteInfo(\@testdirs);
    my $teststatus;
	my $nlfailreport;
    ($teststatus, $nlfailreport) = getTestStatus(\@testdirs, $tagname, $testid);
    &Debug( eval { Dumper $teststatus} );
    &Debug( eval { Dumper \%suiteinfo } );
    my $testxml = &makeResultsXml($teststatus, \%suiteinfo, $nlfailreport);
	&sendresults(\%teststatus, \%suiteinfo, $testxml);
	&printreport($teststatus, $nlfailreport) if $printreport;
    exit(0);
    sleep($sleeptime);
}
#-------------------------------------------------------------------------------
# End Main
#-------------------------------------------------------------------------------
  
#-------------------------------------------------------------------------------
# Get the options. 
#-------------------------------------------------------------------------------
sub opts
{
    my $opt_help;
    GetOptions("testroot=s" => \$testroot,
	          "tagname=s" => \$tagname,
	          "testid=s" => \$testid,
	          "debug|d" => \$debug,
	          "testtype=s" => \$testtype,
	          "help" => \$opt_help,
              "dumpxml|x" => \$dumpxml,
			  "printreport|p" => \$printreport,
	);
    # Show usage if the required options aren't specified. 
    &help if (defined $opt_help);
    &usage if ( (! defined $testroot) ||(! defined $tagname) || (! defined $testid) || (!defined $testtype) );
}

#-------------------------------------------------------------------------------
# Show the usage, and exit.  
#-------------------------------------------------------------------------------
sub usage
{
    print <<'END';
    Usage: 
./testreporter --testroot /glade/scratch/$user/tests/cesm1_1_alphaXX --tagname cesm1_1_alpha15c --testid testid --testtype prealpha|prebeta|prerelease
               [--dumpxml|--printreport]
END
    exit(1);
}

sub help
{
    
    print <<'END';
This is the CESM test reporter script, intended to be used to simplify the reporting of CESM test sets.  
Usage is as follows:
./testreporter --testroot /glade/scratch/$user/tests/cesm1_1_alphaXX 
--tagname cesm1_1_alpha15c --testid testid --testtype prealpha
It gathers all the test results found in the --testroot, gets the relevant test 
status fields, and sends the results to the test reporting system on csegweb. 

    Options:
--testroot This is the testroot you defined when running the test suite
--tagname  The name of the tag you are testing. 
           For example, if you are testing a sandbox that will eventilally 
           become cesm1_1_alpha16c, then put cesm1_1_alpha16c. If you are 
           running a suite that will become a beta tag, put the last alpha
       tag used.  
--testid   The testid you specified to create_test_suite.  
--testtype The type of test you are running: prealpha, prebeta, or prerelease. 
           
END
}

#-------------------------------------------------------------------------------
# Show debugging information if desired. 
#-------------------------------------------------------------------------------
sub Debug
{
    if($debug)
    {
	my ($msg) = @_;
	chomp $msg;
	print "Debug: $msg\n";
    }
}

#-------------------------------------------------------------------------------
# Using the testroot, find the test directories that end with the specified testid.
# If no matching directories are found, then exit. 
#-------------------------------------------------------------------------------
sub getTestDirs
{
    my ($testd, $tid) = @_;

    # Abort if the testroot does not exist.  
    if ( ! -d $testd)
    {
	    print STDERR "The testroot does not exist! Aborting...\n";
	    exit(1);
    }
    # open the testroot, find the test directories.  If no test directories ending in 
    # the testid exist, then abort. 
    opendir(my $DIR, $testd) or die "can't open $testd, error was $!";
    my @testdirs = grep { $_ =~ /($tid)$/ } readdir($DIR);
    closedir $DIR;
    &Debug("in gettestdirs: test directories: \n");
    &Debug( eval { Dumper \@testdirs} );
    if(@testdirs)
    {
	    return @testdirs;
    }
    else
    {
	    print STDERR "It appears that the test root exists, but there aren't any test directories\n";
	    print STDERR "in the testroot. Aborting...\n"; 
	    exit(1);
    }
    return @testdirs;
}

#-------------------------------------------------------------------------------
# Get the suite info from config_definition, and send it back.  
#-------------------------------------------------------------------------------
sub getTestSuiteInfo
{
    my $testlist = shift;
    my %caseinfo; 
    #my $firsttest = (@$testlist)[0];
	my $testpath;
	foreach my $testdir(@$testlist)
	{
		$testpath = $testroot . "/" . $testdir;
		if( -e "$testpath/env_case.xml" && -e "$testpath/env_run.xml" && "$testpath/env_build.xml")
		{
			last;
		}
	}
    #my $abspath = $testroot . "/" . $firsttest;
    &Debug("test path:  $testpath\n");
    my @dirs = ( $testpath, $testpath . "/Tools");
    unshift @INC, @dirs;
    require XML::Lite;
    require ConfigCase;
    my $caseenv = ConfigCase->new("$testpath/Tools/config_definition.xml", "$testpath/env_case.xml");
    my $runenv = ConfigCase->new("$testpath/Tools/config_definition.xml", "$testpath/env_run.xml");
    my $buildenv = ConfigCase->new("$testpath/Tools/config_definition.xml", "$testpath/env_build.xml");
    $caseinfo{'mach'} = $caseenv->get('MACH');
    $caseinfo{'compiler'} = $buildenv->get('COMPILER');
    $caseinfo{'mpilib'} = $buildenv->get('MPILIB');

    my $testspecfile = "$testroot/testspec.$testid.$caseinfo{'mach'}.xml";
	Debug( "testspecfile: $testspecfile\n");
	my $parser = XML::LibXML->new;
	my $spec = $parser->parse_file($testspecfile);
	my $root = $spec->getDocumentElement();
    my @bltagnodes = $root->findnodes('/testlist/baselinetag');
 	$caseinfo{'baselinetag'} = $bltagnodes[0]->textContent();
    Debug( "baselinetag: $caseinfo{'baselinetag'}\n") ;
    
    &Debug("caseinfo: " . eval { Dumper \%caseinfo} );
    return %caseinfo;
}

#-------------------------------------------------------------------------------
# Get the test status. open the $testroot , look for all the test directories, 
# then get the test status. 
#-------------------------------------------------------------------------------
sub getTestStatus
{
    my ($testdirs, $tag)  = @_;
    my %teststatushash;
	my %nlreporthash;
    my $time = localtime;
    print "$time\n";
    
    # Iterate through each of the test directories, and get the requisite test information. 
    foreach my $testcase(@$testdirs)
    {
	    # Get the test status 
	    my $statusfile = $testroot  . "/"  . $testcase .  "/" . $teststatusfilename; 
      if( ! -e $statusfile)
      {
	       warn("$statusfile does not exist, skipping to next test.");
         $teststatushash{$testcase}{'status'} = "TFAIL";
         $tetstatushash{$testcase}{'comment'} = "TestStatus file could not be found!";
         next;
      }
	    &Debug( "Status file: $statusfile\n");
	    open (my $teststatusfile, "<", $statusfile) or die "cannot open TestStatus file for $testcase, $!";
	    my $teststatus = <$teststatusfile>;
	    chomp $teststatus;
	    $teststatus = (split(/\s+/, $teststatus))[0];
	    &Debug("Testcase:   $testcase\n");
	    &Debug( "Teststatus: $teststatus\n"); 
	    $teststatushash{$testcase}{'status'} = $teststatus;

	    # Now go through the TestStats getting the memleak, compare, baseline tag, throughput, and comments if any. 
	    my @statuslines = <$teststatusfile>;

	    my @comparelines = grep { /compare_hist/} @statuslines;
	    my ($comparestatus,$comparetest)  = split(/\s+/, $comparelines[0]);
	    $teststatushash{$testcase}{'compare'} = $comparestatus;
	    my $comparetag = (split(/\./, $comparetest))[-1];
	    $baselinetag = $comparetag unless defined $baselinetag;

  
        # If this a normal test, ie NOT a set of SBN tests, then 
        # send the following fields. If this is a namelist test, then skip these fields, they
        # will never be filled in for SBN tests. 
        if($testtype ne 'namelist') 
		{
			my @memleaklines = grep { /memleak/ } @statuslines;
			my $memleakstatus = (split(/\s+/, $memleaklines[0]))[0];
  	  		$teststatushash{$testcase}{'memleak'} = $memleakstatus;

			my @memcomplines = grep { /memcomp/} @statuslines;
			my $memcompstatus = (split(/\s+/, $memcomplines[0]))[0];
			$teststatushash{$testcase}{'memcomp'} = $memcompstatus;

			my @tputcomplines = grep { /tputcomp/ } @statuslines;
			my $tputcompstatus = (split(/\s+/, $tputcomplines[0]))[0];
			$teststatushash{$testcase}{'tputcomp'} = $tputcompstatus;
		}

	    my @nlcomplines = grep { /nlcomp/i } @statuslines;
	    my $nlcompstatus = (split(/\s+/, $nlcomplines[0]))[0];
	    $teststatushash{$testcase}{'nlcomp'} = $nlcompstatus;

	    my @commentlines = grep { /COMMENT/ } @statuslines;
	    my $comment = (split(/\s+/, $commentlines[0], 2) )[1];
	    chomp $comment;
	    $teststatushash{$testcase}{'comment'} = $comment;
	
	    close $teststatusfile;
	
	    # Check the CaseStatus, and print out the last line...
	    my $casestatusfile = $testroot . "/"  . $testcase . "/" . $casestatusfilename;
      if( -e $casestatusfile)
      {
        open (my $casestatusfile, "<", $casestatusfile) or die "cannot open CaseStatusfile for $testcase, $!";

        my $lastline;
        while(<$casestatusfile>)
        {
          $lastline = $_ if eof;
        }
        close $casestatusfile;
        chomp $lastline;
        &Debug ("last line of CaseStatus: $lastline\n");
        $teststatushash{$testcase}{'casestatus'} = $lastline;
      }
      else
      {
        &Debug("Case status file $casestatus doesn't exist");
        $teststatushash{$testcase}{'casestatus'} = "CaseStatus file not found";
      }

	    # If the test is an IOP test, set a flag in the test status hash indicating it as such.
	    if($testcase =~ /IOP\./)
	    {
	      $teststatushash{$testcase}{'isioptest'} = "true";
	    }

  	  # Get the IOP test status if the file exists.   
	    # If so, create separate iop* entries in the teststatus hash for this test.  
	    my $iopstatusfile = $testroot . "/" . $testcase . "/" . $iopstatusfilename;
	    if( -e $iopstatusfile)
	    {
	      open (my $iopfh, "<", $iopstatusfile) or die " cannot open IOP status file for $testcase, $!";
	      my $iopstatus = <$iopfh>;
	      chomp $iopstatus;
	      $iopstatus = (split(/\s+/, $iopstatus))[0];
	      $teststatushash{$testcase}{'iopstatus'} = $iopstatus;
	      @statuslines = <$iopstatusfile>;
	    
	      @memleaklines = grep { /memleak/ } @statuslines;
	      $memleakstatus = (split(/\s+/, $memleaklines[0]))[0];
	      $teststatushash{$testcase}{'iopmemleak'} = $memleakstatus;

        @comparelines = grep { /compare_hist/} @statuslines;
	      ($comparestatus,$comparetest)  = split(/\s+/, $comparelines[0]);
	      $teststatushash{$testcase}{'compare'} = $comparestatus;
	      $comparetag = (split(/\./, $comparetest))[-1];
	      $teststatushash{$testcase}{'iopbaselinetag'} = $comparetag;

	      @memcomplines = grep { /memcomp/} @statuslines;
	      $memcompstatus = (split(/\s+/, $memcomplines[0]))[0];
	      $teststatushash{$testcase}{'iopmemcomp'} = $memcompstatus;

	      @tputcomplines = grep { /tputcomp/ } @statuslines;
	      $tputcompstatus = (split(/\s+/, $tputcomplines[0]))[0];
	      $teststatushash{$testcase}{'ioptputcomp'} = $tputcompstatus;

  	    @nlcomplines = grep { /nlcomp/i } @statuslines;
	      $nlcompstatus = (split(/\s+/, $nlcomplines[0]))[0];
	      $teststatushash{$testcase}{'iopnlcomp'} = $nlcompstatus;

	      @commentlines = grep { /COMMENT/ } @statuslines;
	      $comment = (split(/\s+/, $commentlines[0], 2) )[1];
	      chomp $comment;
	      $teststatushash{$testcase}{'iopcomment'} = $comment;
	    
	      $teststatushash{$testcase}{'iopcasestatus'} = $teststatushash{$testcase}{'casestatus'};
	    
	      close $iopfh;
	    }


        if($testtype eq 'namelist')
		{
			my $statusoutfile = $testroot  . "/"  . $testcase .  "/" . $teststatusoutfilename;
            if( -e $statusoutfile  && $teststatushash{$testcase}{'nlcomp'} eq 'FAIL')
			{
				open my $STATUSOUT, "<", $statusoutfile or warn "can't open $statusoutfile, $!";
				my @statusoutlines = <$STATUSOUT>;	
				my $commentflag = 0;
				foreach my $soline(@statusoutlines)
				{
					chomp $soline;
					next if ($soline =~ /^$/);
					$commentflag = 1 if ($soline =~ /^FAIL/);	
					$commentflag = 0 if ($soline =~ /^PASS/);
					if($commentflag)
					{
						$teststatushash{$testcase}{'comment'} .= "$soline\n";
						Debug( eval { Dumper \$teststatushash{$testcase}{'comment'} } );
					}
				}
				$nlreporthash{$teststatushash{$testcase}{'comment'}}{$testcase} = 1;
			}
			foreach my $blank(keys %nlreporthash)
			{
				delete $nlreporthash{$blank} if(length $nlreporthash{$blank} == 0);
			}
		}
    }
	
	foreach my $blank(keys %nlreporthash)
	{
		delete $nlreporthash{$blank} if(length $nlreporthash{$blank} == 0);
	}
	
    return \%teststatushash, \%nlreporthash;
}


# New sendresults, using XML::LibXML
# This is a work in progress, it is not yet functional.  
# TODO: finish implementing the programmatic text-based implementation. 
sub makeResultsXml
{
	my ($testresults, $suiteinfo, $nlfailreport) = @_;
	my $testxml = XML::LibXML::Document->new('1.0', 'UTF-8');
    my $root =  $testxml->createElement('testrecord');
	$root->appendTextChild('tag_name', $tagname);
	$root->appendTextChild('mach', $suiteinfo{'mach'});
    my $compilerelem = XML::LibXML::Element->new('compiler');
    $compilerelem->setAttribute('version', '');
    $compilerelem->appendText($suiteinfo{'compiler'});
    $root->appendChild($compilerelem);
    my $mpielem = XML::LibXML::Element->new('mpilib');
    $mpielem->setAttribute('version', '');
    $mpielem->appendText($suiteinfo{'mpilib'});
    $root->appendChild($mpielem);
	$root->appendTextChild('testroot', $testroot);
	$root->appendTextChild('testtype', $testtype);
	$root->appendTextChild('baselinetag', $suiteinfo{'baselinetag'});

	foreach my $test(reverse sort keys %$testresults)
	{
		my $testelem = $testxml->createElement('tests');
		$testelem->setAttribute('testname', $test);
		foreach my $detail(sort keys %{$$testresults{$test}})
		{
			my $catelem = $testxml->createElement('category');
			$catelem->setAttribute('name', $detail);
			if($detail eq 'comment')
			{
				$cdata = XML::LibXML::CDATASection->new($$testresults{$test}{$detail});
				$catelem->appendChild($cdata);
			}
			else
			{
				$catelem->appendText($$testresults{$test}{$detail});
			}
			$testelem->appendChild($catelem);
			#print Dumper $testelem;
			
		}
		$root->appendChild($testelem);
	}

	if($testtype eq 'namelist')
	{
		my $nlfailreportelem = XML::LibXML::Element->new('namelistfailuresreport');
		foreach my $nlfailkey(sort keys %$nlfailreport)
		{
			my $nlfailelem = $testxml->createElement('namelistfailure');
			my $nlfailtxt = $testxml->createElement('namelistfailuretext');
			my $nlfailcdata = XML::LibXML::CDATASection->new($nlfailkey);
			$nlfailtxt->appendChild($nlfailcdata);
			$nlfailelem->appendChild($nlfailtxt);

			foreach my $testname(sort keys %{$$nlfailreport{$nlfailkey}})
			{
				$nlfailelem->appendTextChild('test', $testname);
			}
			$nlfailreportelem->appendChild($nlfailelem);
		}
	 	$root->appendChild($nlfailreportelem);
	}

	$testxml->setDocumentElement($root);
	$xmlstring = $testxml->toString(1);
	Debug("testxml file: ");
    Debug( $xmlstring);
	if($debug || $dumpxml)
	{
		open my $xmldumpfile, ">", "./testreporter.dump.xml" or die $!;
		print $xmldumpfile $testxml->toString(1);
		close $xmldumpfile;
		print "wrote xml test report to testreporter.dump.xml\n";
    }	
	return $testxml;
}

sub sendresults
{
	my ($testresults, $suiteinfo, $testxml) = @_;
	my $resultsstr = $testxml->toString(1);

    my $useragent = LWP::UserAgent->new(ssl_opts => {verify_hostname => 0});
    my $req = HTTP::Request->new(POST => $posturl);
    $req->content_type('application/x-www-form-urlencoded');
    $req->content("username=$username&password=$password&testXML=$resultsstr");
    my $response = $useragent->request($req);

    if($response->is_success)
    {
        print "Test results successfully posted to csegweb\n";
    }
    elsif($response->code eq '401')
    {
      my $errmsg = "The server responded with '401 - Unauthorized'\n";
      $errmsg .=   "Your svn username & password is most likely incorrect!\n";
      $errmsg .=   "Please re-run the script, and provide the correct svn username & password\n";
      die $errmsg;

    }
    else
    {
        print "Posting the results to $posturl failed! \n";
        my $status = $response->status_line;
        print "$status\n";
        print "aborting!\n";
        exit(1);
    }
}

sub printreport
{
	my ($teststatus, $nlfailreport) = @_;
	foreach my $basefail(keys %$nlfailreport)
	{
		my @tests = sort keys %{$$nlfailreport{$basefail}};
    	print "====================================================================\n";
    	print "namelist difference: \n";
    	print "$basefail\n";
    	print "--------------------------------------------------------------------\n";
		print "Tests which had the above namelist diffs:\n";
    	map { print "$_\n"} sort @tests;
    	print "--------------------------------------------------------------------\n";
	}
}

#
#-------------------------------------------------------------------------------
# Send the test results to csegweb testdb
#-------------------------------------------------------------------------------
sub sendresultsOld
{
    my ($testresults, $suiteinfo) = @_;
    my $resultsstr = "";


#-------------------------------------------------------------------------------
# ASB - using the following DTD to send XML formatted results, comments 
# wrapped in CDATA to ensure proper escaping and encoding of anything that 
# could possibly 
#<testrecord>
#  <tag_name> </tag_name>
#  <machine> </machine>
#  <compiler version=' '> </compiler>
#  <mpilib version=' '> </mpilib>
#  <testroot> </testroot>
#  <testtype> </testtype>
#  <tests testname name=' '>
#      <category name=' '> </category>
#  </tests>
#</testrecord>
#-------------------------------------------------------------------------------
    $resultsstr .= "<?xml version='1.0'?>";
    $resultsstr .= "<testrecord>\n"; 
    $resultsstr .= "<tag_name>$tagname</tag_name>\n"; 
    $resultsstr .= "<mach>$suiteinfo{'mach'}</mach>\n";
    $resultsstr .= "<compiler version=\'$suiteinfo{'compiler_version'}\'>$suiteinfo{'compiler'}</compiler>\n";
    $resultsstr .= "<mpilib version=\'$suiteinfo{'mpilib_version'}\'>$suiteinfo{'mpilib'}</mpilib>\n";
    $resultsstr .= "<testroot>$testroot</testroot>\n";
    $resultsstr .= "<testtype>$testtype</testtype>\n";
    $resultsstr .= "<baselinetag>$baselinetag</baselinetag>\n";

    foreach my $test(sort keys %$testresults)
    {
	    $resultsstr .= "<tests testname='$test'>\n";
	    foreach my $detail(sort keys %{$$testresults{$test}})
	    {
		  if($detail eq 'comment')
		  {
			$resultsstr .= "<category name='$detail'><![CDATA[$$testresults{$test}{$detail}]]></category>\n";
		  }
		  else
		  {
	      	$resultsstr .= "<category name='$detail'>$$testresults{$test}{$detail}</category>\n";
		  }
	    }
	    $resultsstr .= "</tests>\n";
    }
    $resultsstr .= "</testrecord>\n"; 
    &Debug("resultsstr: \n $resultsstr\n");
    if($debug)
    {
        open my $xmldumpfile, ">", "./testreporter.dump.xml" or die $!;
        print $xmldumpfile $resultsstr;
        close $xmldumpfile;
    }

    my $useragent = LWP::UserAgent->new(ssl_opts => {verify_hostname => 0});
    my $req = HTTP::Request->new(POST => $posturl);
    $req->content_type('application/x-www-form-urlencoded');
    $req->content("username=$username&password=$password&testXML=$resultsstr");  
    my $response = $useragent->request($req);

    if($response->is_success)
    {
	    print "Test results successfully posted to csegweb\n";
    }
    elsif($response->code eq '401')
    {
      my $errmsg = "The server responded with '401 - Unauthorized'\n";
      $errmsg .=   "Your svn username & password is most likely incorrect!\n";
      $errmsg .=   "Please re-run the script, and provide the correct svn username & password\n";
      die $errmsg;
      
    }
    else 
    {
	    print "Posting the results to $posturl failed! \n";
	    my $status = $response->status_line;
	    print "$status\n";
	    print "aborting!\n";
	    exit(1);
    }
}

sub authenticate
{
  print "Enter your username: \n";
  $username = <STDIN>;
  print "Enter your password: \n";
  system('stty','-echo');
  chop($password = <STDIN>);
  system('stty','echo');
  chomp $username;
  chomp $password;
}
