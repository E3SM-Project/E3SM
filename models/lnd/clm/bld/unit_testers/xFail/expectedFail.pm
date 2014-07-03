=head1 expectedFail.pm

Documentation for expectedFail.pm

=head1 Overview

The module expectedFail.pm supplies the capability of checking if a failed test is expected to fail.
It is called directly from either test_driver.sh (for batch and interactive tests) or build-namelist_test.pl.  
Future plans involve integrating this module into cesm tests.

=head1 Use Case

This is a new feature being added to the existing CLM test infrastructure.  The use case would roughly be 
along the lines of:

   1) Run the test suite (CLM batch,interactive or namelist) 
   2) Search for test failures 
      a) Fix failed tests
      b) -or- Add new xFail entries to XML file if a test is supposed to fail (eg. due to some missing resolution).
   3) Check for new tests that now pass.  This is for modifying the ChangeLog.
   4) Update XML file by either adding new entries or removing old ones.
   5) update the ChangeLog to reflect important changes in test behavior (Tests that now pass that failed before, tests that 
         are now xFail, etc...

=head2 Public methods

There are two public methods needed.  The "new" ctor and one of the parseOutput* methods.
Everything else is private.

   xFail::expectedFail->new
   parseOutput
   parseOutputCLM

=head2 Private methods

   sub _searchExpectedFail
   sub _readXml
   sub _testNowPassing
   sub _printOutput
   sub _getTestType
   sub _getMachInfo
 
=cut

package xFail::expectedFail;

our $VERSION = '1.00';
 
use Cwd;
use strict;
use Getopt::Long;
use English;
use Scalar::Util qw(looks_like_number);

my @testList={};
my $DEBUG=0;

my $pass=" PASS";
my $fail=" FAIL";
my $xfail="xFAIL";

##############################################################################
#
##############################################################################

=head1 CTOR

Constructor for the class.  Reads in three arguments:
   _callingName     -> name of the script creating the new object
   _compareGenerate -> compare or generate option
   _totTests        -> total number of tests to run

Calls _readXml which reads the file expectedClmTestFails.xml and stores it memory
for later searches.

returns:  new object ($self)

=cut

##############################################################################
#
##############################################################################
sub new {
   my ($class_name) = @_;
   my $self = {
      _className             => shift,
      _callingName           => shift,
      _compareGenerate       => shift,
      _totTests              => shift,
      _foundList             => undef,
      _numericalTestId       => undef
   };

   if ($DEBUG) {
      print "$self->{_callingName}\n";
      print "$self->{_compareGenerate}\n";
   }

    bless ($self, $class_name);

    $self->{_numericalTestId}=0;
    $self->{_created} = 1;

    $self->_readXml();

    return $self;
}

##############################################################################
#
##############################################################################

=head1 parseOutput

parseOutput parsese the output from the build-namelist_test.pl script.  It is similar
to, but not interchangable with parseOutputCLM.

The only argument is that of the reference variable that contains the information dumped
by Test::More.

returns:  nothing

=cut

##############################################################################
#
##############################################################################
sub parseOutput
{


   my $report;
   my $testId;
   my @testName={};
   my $testReason;

   my ($self, $output) = @_ ;
    
   #_#===========================================
   #_# keep this in for logging
   #_#===========================================
   print ("captured output is :: \n $output \n");

   #_# split the output from Test::More output on newline
   my @refList = split('\n', $output);

   #_# process any buffered output which happens when a subroutine from build-namelist_test.pl 
   #_# itself calls some testing routines
   foreach my $refSplit (@refList) {

      #_# always look at the last element of refSplit since that will have the info. from the
      #_# last test run

      my @outArr=split(/ /,$refSplit);

      if ($DEBUG) {
         print ("\nxFail::expectedFail::parseOutput @outArr[0] \n");
         print ("xFail::expectedFail::parseOutput @outArr[1] \n");
         print ("xFail::expectedFail::parseOutput @outArr[2] \n");
         print ("xFail::expectedFail::parseOutput @outArr[3] \n");
         print ("xFail::expectedFail::parseOutput @outArr[4] \n");
      }

      my $size = @outArr-1;

      #_# first case, we have a passed (ok) test
      if (@outArr[0] eq "ok") {
         $self->{_numericalTestId}++;

         $report=$pass;
         $testId=@outArr[1];
         @testName=@outArr[3..$size];
         $testReason="";

         my ($retVal,$xFailText)=$self->_searchExpectedFail($testId);

         my $testReason=$self->_testNowPassing($testId,$retVal,$xFailText);

         if($DEBUG){
            print("$testReason \n");
         }

         $self->_printOutput($report,$testId,$testReason,@testName);


      #_# deal with the case of a failed (not ok) test
      } elsif  (@outArr[0] eq "not") {
         $self->{_numericalTestId}++;

         $testId=@outArr[2];
         my ($retVal,$xFailText)=$self->_searchExpectedFail($testId);
         
         if ($DEBUG) {
            print ("xFail::expectedFail::parseOutput Id $retVal,$xFailText \n");
         }
   
         @testName=@outArr[4..$size];
   
         if ($retVal eq "TRUE"){
            #_# found an expected FAIL (xFAIL)
            $report=$xfail;
            $testReason= "<Note: $xFailText>";
         } else {
            #_# print a regular FAIL
            $report=$fail;
            $testReason="";
         }

         $self->_printOutput($report,$testId,$testReason,@testName);

      } else {
         #_# skipping line.  Trying to parse error code from Test::More
      }

   }

   #_# this resets the reference that points to $output (\$captOut) on the caller side
   @_[1]="";
   
}
   
##############################################################################
#
##############################################################################

=head1 parseOutputCLM

parseOutputCLM parsese the output from the test_driver.sh script.  It is similar
to, but not interchangable with parseOutput.

parseOutputCLM takes one arguments:
   $statFoo->     the name of the td.<pid>.status file 

returns:  nothing

=cut

##############################################################################
#
##############################################################################
sub parseOutputCLM
{

   my $report;
   my $testId;
   my @testName={};
   my $testReason;

   my ($self, $statFoo) = @_ ;
    
   open(FOO, "< $statFoo"); # open for input
   open(FOO_OUT, "> $statFoo.xFail"); # open for input

   my(@reportLines);

   while (<FOO>) {

      my($line) = $_;

      my @outArr=split(/ /,$line);
      if (looks_like_number(@outArr[0])) {

         $self->{_numericalTestId}++;

         my $num=sprintf("%03d", $self->{_numericalTestId});
         my $totNum=sprintf("%03d", $self->{_totTests});

         #_# last element has the pass/fail info.
         chomp(@outArr[-1]);
         my $repPass=substr(@outArr[-1], -4, 4);

         if ($DEBUG) {
            print ("xFail::expectedFail::parseOutput @outArr[0] \n");
            print ("xFail::expectedFail::parseOutput @outArr[1] \n");
            print ("xFail::expectedFail::parseOutput @outArr[2] \n");
            print ("xFail::expectedFail::parseOutput @outArr[3] \n");
            print ("xFail::expectedFail::parseOutput @outArr[4] \n");
            print ("xFail::expectedFail::parseOutput @outArr[5] \n");
            print ("xFail::expectedFail::parseOutput @outArr[6] \n");
            print ("xFail::expectedFail::parseOutput @outArr[-1] \n");
            print ("xFail::expectedFail::parseOutput $repPass \n");
         }

         my $size = @outArr-1;
         if ($DEBUG) {
            print ("size of line $size \n");
         }
         my $endOfDesc=$size-1;
   
         if ($repPass eq "PASS") {
            $report=$pass;
            $testId=@outArr[1];
            @testName=@outArr[2..$endOfDesc];

            my ($retVal,$xFailText)=$self->_searchExpectedFail($testId);

            my $testReason=$self->_testNowPassing($testId,$retVal,$xFailText);

            #_# print out the test results
            print FOO_OUT ("$num/$totNum <$report> <Test Id: $testId> <Desc: @testName> $testReason \n"); 

         } else {
            $testId=@outArr[1];
            my ($retVal,$xFailText)=$self->_searchExpectedFail($testId);
            
            if ($DEBUG) {
               print ("xFail::expectedFail::parseOutput Id $retVal,$xFailText \n");
            }
         
            @testName=@outArr[2..$endOfDesc];
         
            if ($retVal eq "TRUE"){
               #_# found an expected FAIL (xFAIL)
               $report=$xfail;
               $testReason= "<Note: $xFailText>";
            } else {
               #_# print a regular FAIL
               $report=$fail;
               $testReason="";
            }

            #_# print out the test results
            print FOO_OUT ("$num/$totNum <$report> <Test Id: $testId> <Desc: @testName> $testReason \n"); 

         }

      } else {
         print FOO_OUT $line;
      }
   }
   close(FOO);
   close(FOO_OUT);
}

##############################################################################
#
##############################################################################

=head1 _searchExpectedFail

searches the list of expected fails for a match with testId.

_searchExpectedFail takes one arguments:
   $testId->    the test id (numerical or string) that we want to search for

returns: $retVal (TRUE or FALSE) if id was found
         $text   text from XML file

=cut

##############################################################################
#
##############################################################################
sub _searchExpectedFail
{
   my ( $self,$testId) = @_;

   #search through list for test ID
   my $retVal="FALSE";

   if ($DEBUG) {
      print ("here 2 Id $self->{_foundList} \n");
   }
   if ($self->{_foundList} eq "FALSE"){
      if ($DEBUG) {
         print ("returning early Id \n");
      }
      return $retVal;
   }

   my $failType;
   my $text;
   foreach my $tL (@testList) {
      my %tAtts = $tL->get_attributes();
      my $tid=$tAtts{'testId'};
      if ($DEBUG) {
         print ("_seachExpectedFail Id $tid $testId \n");
      }
      if ($tid eq $testId) {
         if ($DEBUG) {
            print ("here Id \n");
         }  
         #~# found the test we're looking for
         $text=$tL->get_text();
         $failType=$tAtts{'failType'};
         if ($failType eq "xFail"){
            $retVal="TRUE";
         }
      }
   }
   return ($retVal,$text);
}

##############################################################################
#
##############################################################################

=head1 _readXml

reads the xml file for a particular machine, compiler, test type and (compare
| generate) setup and saves it in memory for searching by _searchExpectedFail.

_readXml takes no arguments

returns: nothing

=cut

##############################################################################
#
##############################################################################
sub _readXml
{
   my ( $self ) = @_;

   #Figure out where configure directory is and where can use the XML/Lite module from
   my $ProgName;
   ($ProgName = $PROGRAM_NAME) =~ s!(.*)/!!; # name of program
   my $ProgDir = $1;                         # name of directory where program lives
   
   my $cwd = getcwd();  # current working directory
   my $cfgdir;
   
   if ($ProgDir) { $cfgdir = $ProgDir; }
   else { $cfgdir = $cwd; }
   
   #-----------------------------------------------------------------------------------------------
   # Add $cfgdir to the list of paths that Perl searches for modules
   my @dirs = ( $cfgdir, "$cfgdir/perl5lib",
               "$cfgdir/../../../../../scripts/ccsm_utils/Tools/perl5lib",
               "$cfgdir/../../../../../models/utils/perl5lib",
            );
   unshift @INC, @dirs;
   my $result = eval "require XML::Lite";
   if ( ! defined($result) ) {
      die <<"EOF";
   ** Cannot find perl module \"XML/Lite.pm\" from directories: @dirs **
EOF
   }

   #-----------------------------------------------------------------------------------------------
   
   my ($machine,$compiler)=_getMachInfo();

   my $testType=$self->_getTestType($self->{_callingName});


   my $xmlFile=undef;
   if ($testType eq "clmInteractive" || $testType eq "clmBatch") { 
      $xmlFile         = "$cfgdir/expectedClmTestFails.xml";
   } elsif ($testType eq "namelistTest") {
      $xmlFile         = "xFail/expectedClmTestFails.xml";
   } else {
      $xmlFile         = "xFail/expectedClmTestFails.xml";
   }
   my $xml = XML::Lite->new($xmlFile);

   my $root = $xml->root_element();
   
   if ($DEBUG) {
      print "_readXml $self->{_callingName}\n";
      print "_readXml $self->{_compareGenerate}\n";
      print "_readXml $xmlFile \n";
      print ("_readXml Debug testType $testType \n");
      print ("_readXml Debug machine $machine \n");
      print ("_readXml Debug compiler $compiler \n");
   }

   # Check for valid root node
   my $name = $root->get_name();
   $name eq "expectedFails" or die
      "readExpectedFail.pm::_readXml :: $xmlFile is not a file that contains expected test failures\n";
   
   my @e = $xml->elements_by_name($testType);

   $self->{_foundList}="FALSE";
   
   ### populate list of tests for a specfic test type, machine and compiler
   ### there's got to be a better way to write this
   while ( my $e = shift @e ) {
      my @mChildren = $e->get_children();
      foreach my $mChild (@mChildren) {
         my $mName=$mChild->get_name();
         if ($mName eq $machine){
            my @cChildren = $mChild->get_children();
            foreach my $cChild (@cChildren) {
               my $cName=$cChild->get_name();
               if ($cName eq $compiler) {
                  my @cgChildren=$cChild->get_children();
                  foreach my $cgChild (@cgChildren) {
                     my $cgName=$cgChild->get_name();
                     if($cgName eq $self->{_compareGenerate}){
                        @testList=$cgChild->get_children();
                        $self->{_foundList}="TRUE";
                        last;
                     }
                  }
               }
            }
         }
      }
   }
   if ($DEBUG) {
      print ("here 1 $self->{_foundList} \n");
   }
}

##############################################################################
#
##############################################################################

=head1 _testNowPassing

reads the xml file for a particular machine, compiler, test type and (compare
| generate) setup and saves it in memory for searching by _searchExpectedFail.

_testNowPassing takes three arguments:
   $id - test id to print out
   $retVal - TRUE or FALSE.  Was the id found in the expected fail list
   $xmlText - text from the XML notes section of the file.  (Currently not used,
               may be used in future).

   returns:  a text string 

=cut

##############################################################################
#
##############################################################################
sub _testNowPassing
{

   my ($self, $id, $retVal, $xmlText) = @_ ;
   my $text=undef;

   if ($retVal eq "TRUE") {
      #_# found a test that passes now, but is listed as an xFail 
      $text = "<NOTE: $id is a new PASS; was xFAIL>\n";
      
   } else {
      #_# this test passes and was not previously listed as an xFail
      #_# noOp
   }

   return $text;
}

##############################################################################
#
##############################################################################

=head1 _printOutput

method that prints output for status files.

_printOutput takes four arguments:
   $report     - PASS,FAIL,xFAIL
   $testId     - test id to print out
   $testReason - for xFAIL and new PASSES, additional reporting
   @testName   - test description from original test

   returns:  a text string 

=cut

##############################################################################
#
##############################################################################
sub _printOutput
{

   my ($self, $report, $testId, $testReason, @testName) = @_ ;

   #_# print out the test results
   my $num=sprintf("%03d", $self->{_numericalTestId});
   my $totNum=sprintf("%03d", $self->{_totTests});
   print ("$num/$totNum <$report> <Test Id: $testId> <Desc: @testName> $testReason \n"); 

}

##############################################################################
#
##############################################################################

=head1 _getTestType

method that takes the name of the calling script and returns the type of 
test.  Used for searching the expected fail list.

_getTestType takes four arguments:
   $name       - name of calling script

   returns:  $type, the type of test

=cut

##############################################################################
#
##############################################################################
sub _getTestType
{

   my ($self, $name) = @_ ;

   if($DEBUG){
      print ("_getTestType $name");
   }

   my %testTypes = (
      "build-namelist_test.pl"  => "namelistTest",
      "test_driver.sh-i"   => "clmInteractive",
      "test_driver.sh"   => "clmBatch",
      "clm-cesm.sh"   => "cesm"
   );
   
   my $type = $testTypes {lc $name} || "unknown";
   return $type;

}

##############################################################################
#
##############################################################################

=head1 _getMachInfo

method that figures out on what platform this is running and returns a 2 digit
machine identifier and the compiler.  This will eventually contain multiple
compiler for various machines.

_getMachInfo takes no arguments

   returns:  $mach - the machine I'm running on
             $comp - the compiler being used

=cut

##############################################################################
#
##############################################################################
sub _getMachInfo
{

   my $name=`uname -n`;
   $name = substr($name, 0, 2);

   my %machNames = (
      "ys"  => "yellowstone",
      "ja"   => "janus",
      "mi"   => "mirage",
      "ly"   => "lynx",
      "ys"   => "yellowstone"
   );

   my %compNames = (
      "ys"  => "INTEL",
      "ja"   => "intel",
      "mi"   => "PGI",
      "ly"   => "intel",
      "ys"   => "intel"
   );
   
   my $mach = $machNames {lc $name} || "unknown";
   my $comp = $compNames {lc $name} || "unknown";

   print"SPM DEBUG :: $mach , $comp \n";

   return ($mach,$comp);

}

# A Perl module must end with a true value or else it is considered not to
# have loaded.  By convention this value is usually 1 though it can be
# any true value.  A module can end with false to indicate failure but
# this is rarely used and it would instead die() (exit with an error).
1;
