#!/usr/bin/env perl 
use strict;
use warnings;
use Data::Dumper;
use Test::More tests => 28;
use XML::LibXML;
#push(@INC, "..");

#my $testfilesdir = "./testfiles/";
my $testfilesdir = ".";
my $banner = "===============================================================================";
print "$banner\nRUNNING UNIT TESTS\n$banner\n";

ok (require('../manage_xml_entries'), "loaded manage_xml_entries ok..");
#ok (require('manage_xml_entries'), "loaded manage_xml_entries ok..");

# Can we create CESMTest objects??
my $cesmTest = new CESMTest(compset => "B", grid => 'f19_g16', testname => 'ERS', machine => 'frankfurt', 
							compiler => 'intel', 'testmods' => 'clm/Default', comment => 'always test your code', testtype => 'prebeta');
isa_ok($cesmTest, "CESMTest", "we should be able to create CESMTest objects ok..");
# Object creation tests..
ok($cesmTest->{compset} eq 'B', "CESMTest compset should match..");
ok($cesmTest->{grid} eq 'f19_g16', "CESMTest grid should match..");
ok($cesmTest->{testname} eq 'ERS', "CESMTest testname should match..");
ok($cesmTest->{machine} eq 'frankfurt', "CESMTest frankfurt should match..");
ok($cesmTest->{compiler} eq 'intel', "CESMTest compiler should match..");
ok($cesmTest->{testmods} eq 'clm/Default', "CESMTest testmods should match..");
ok($cesmTest->{comment} eq 'always test your code', "CESMTest comment should match..");
ok($cesmTest->{testtype} eq 'prebeta', "CESMTest testtype should match..");


# We should be able to parse the testfiles/frankfurt_nag_prealpha test list...
my @expectedffnag;
my $n1 = new CESMTest(compset => 'A', grid => 'f45_g37_rx1', testname => 'SMS_D_Mmpich', 
					  machine => 'frankfurt', compiler => 'nag', comment => 'This is a comment');

my $n2 = new CESMTest(compset => 'F', grid => 'f19_f19', testname => 'ERS_D_Mmpich', 
					  machine => 'frankfurt', compiler => 'nag' );

my $n3 = new CESMTest(compset => 'FSCM5A97', grid => 'T42_T42', testname => 'SMS_D_Mmpi-serial', 
					  machine => 'frankfurt', compiler => 'nag' );

my $n4 = new CESMTest(compset => 'I20TRCLM45BGC', grid => 'f10_f10', testname => 'ERS_D_Mmpich', 
					  machine => 'frankfurt', compiler => 'nag' );

my $n5 = new CESMTest(compset => 'ICLM45BGC', grid => 'ne30_g16', testname => 'ERS_D_Mmpich', 
					  machine => 'frankfurt', compiler => 'nag' );
push(@expectedffnag, $n1);
push(@expectedffnag, $n2);
push(@expectedffnag, $n3);
push(@expectedffnag, $n4);
push(@expectedffnag, $n5);
my $ffnagprealpha = "$testfilesdir/frankfurt_nag_prealpha";
my $actualffnag = parseTextList("$testfilesdir/frankfurt_nag_prealpha");
#print Dumper \@expectedffnag;
#print Dumper $actualffnag;

is_deeply(\@expectedffnag, $actualffnag, "frankfurt nag prealpha tests should parse correctly..");

# We should be able to parse the aux_clm_short testlist...
my @expectedauxclmshort;
my $cs1 = new CESMTest(compset => 'I1850CLM45BGC', grid => 'f10_f10', testname => 'PET_P15x2_D', 
                       machine => 'yellowstone', compiler => 'pgi', testmods => 'clm/ciso', 
                       comment => 'This test fails miserably');
push(@expectedauxclmshort, $cs1);
my $cs2 = new CESMTest(compset => 'I1850CLM45BGC', grid => 'f10_f10', testname => 'PET_P15x2_D', machine => 'yellowstone', 
                       compiler => 'intel', testmods => 'clm/ciso');
push(@expectedauxclmshort, $cs2);
my $cs3 = new CESMTest(compset => 'I1850CLM45BGC', grid => 'f10_f10', testname => 'PET_P16x2_D', machine => 'frankfurt', 
                       compiler => 'pgi', testmods => 'clm/ciso');
push(@expectedauxclmshort, $cs3);
my $cs4 = new CESMTest(compset => 'I1850CLM45BGC', grid => 'f10_f10', testname => 'PET_P16x2_D', machine => 'frankfurt', 
                       compiler => 'intel', testmods => 'clm/ciso');
push(@expectedauxclmshort, $cs4);
my $cs5 = new CESMTest(compset => 'I1850CLM45BGC', grid => 'f10_f10', testname => 'PET_P16x2_D_Mmpich', machine => 'frankfurt', 
                       compiler => 'nag', testmods => 'clm/ciso');
push(@expectedauxclmshort, $cs5);


my $cs6 = new CESMTest(compset => 'I20TRCLM45BGC', grid => 'f10_f10', testname => 'ERS', machine => 'frankfurt', 
                       compiler => 'pgi');
push(@expectedauxclmshort, $cs6);
my $cs7 = new CESMTest(compset => 'I20TRCLM45BGC', grid => 'f10_f10', testname => 'ERS', machine => 'frankfurt', 
                       compiler => 'intel');
push(@expectedauxclmshort, $cs7);
my $cs8 = new CESMTest(compset => 'I20TRCLM45BGC', grid => 'f10_f10', testname => 'ERS', machine => 'yellowstone', 
                       compiler => 'pgi');
push(@expectedauxclmshort, $cs8);
my $cs9 = new CESMTest(compset => 'I20TRCLM45BGC', grid => 'f10_f10', testname => 'ERS', machine => 'yellowstone', 
                       compiler => 'intel');
push(@expectedauxclmshort, $cs9);
my $cs10 = new CESMTest(compset => 'I20TRCLM45BGC', grid => 'f10_f10', testname => 'ERS_D', machine => 'frankfurt', 
                       compiler => 'pgi');
push(@expectedauxclmshort, $cs10);
my $cs11 = new CESMTest(compset => 'I20TRCLM45BGC', grid => 'f10_f10', testname => 'ERS_D', machine => 'frankfurt', 
                       compiler => 'intel');
push(@expectedauxclmshort, $cs11);
my $cs12 = new CESMTest(compset => 'I20TRCLM45BGC', grid => 'f10_f10', testname => 'ERS_D', machine => 'yellowstone', 
                       compiler => 'pgi');
push(@expectedauxclmshort, $cs12);
my $cs13 = new CESMTest(compset => 'I20TRCLM45BGC', grid => 'f10_f10', testname => 'ERS_D', machine => 'yellowstone', 
                       compiler => 'intel');
push(@expectedauxclmshort, $cs13);
my $cs14 = new CESMTest(compset => 'I20TRCLM45BGC', grid => 'f10_f10', testname => 'ERS_D_Mmpich', machine => 'frankfurt', 
                       compiler => 'nag');
push(@expectedauxclmshort, $cs14);
my $cs15 = new CESMTest(compset => 'I20TRCLM45BGC', grid => 'f10_f10', testname => 'ERS_Mmpich', machine => 'frankfurt', 
                       compiler => 'nag');
push(@expectedauxclmshort, $cs15);
my $cs16 = new CESMTest(compset => 'ICLM45BGCCROP', grid => 'f19_g16', testname => 'SMS_Ly1', machine => 'frankfurt', 
                       compiler => 'pgi', testmods => 'clm/reduceOutput');
push(@expectedauxclmshort, $cs16);
my $cs17 = new CESMTest(compset => 'ICLM45BGCCROP', grid => 'f19_g16', testname => 'SMS_Ly1', machine => 'frankfurt', 
                       compiler => 'intel', testmods => 'clm/reduceOutput');
push(@expectedauxclmshort, $cs17);
my $cs18 = new CESMTest(compset => 'ICLM45BGCCROP', grid => 'f19_g16', testname => 'SMS_Ly1', machine => 'yellowstone', 
                       compiler => 'pgi');
push(@expectedauxclmshort, $cs18);
my $cs19 = new CESMTest(compset => 'ICLM45BGCCROP', grid => 'f19_g16', testname => 'SMS_Ly1', machine => 'yellowstone', 
                       compiler => 'intel');
push(@expectedauxclmshort, $cs19);

my $cs20 = new CESMTest(compset => 'ICLM45BGCCROP', grid => 'f19_g16', testname => 'SMS_Ly1_Mmpich', machine => 'frankfurt',
                        compiler => 'nag', testmods => 'clm/reduceOutput');
push(@expectedauxclmshort, $cs20);

my $cs21 = new CESMTest(compset => 'ICN', grid => 'f10_f10', testname => 'ERS_D', machine => 'frankfurt', 
                       compiler => 'pgi');
push(@expectedauxclmshort, $cs21);
my $cs22 = new CESMTest(compset => 'ICN', grid => 'f10_f10', testname => 'ERS_D', machine => 'frankfurt', 
                       compiler => 'intel');
push(@expectedauxclmshort, $cs22);
my $cs23 = new CESMTest(compset => 'ICN', grid => 'f10_f10', testname => 'ERS_D', machine => 'yellowstone', 
                       compiler => 'pgi');
push(@expectedauxclmshort, $cs23);
my $cs24 = new CESMTest(compset => 'ICN', grid => 'f10_f10', testname => 'ERS_D', machine => 'yellowstone', 
                       compiler => 'intel');
push(@expectedauxclmshort, $cs24);
my $cs25 = new CESMTest(compset => 'ICN', grid => 'f10_f10', testname => 'ERS_D_Mmpich', machine => 'frankfurt', 
                       compiler => 'nag');
push(@expectedauxclmshort, $cs25);

my $actualauxclmshort = parseTextList("$testfilesdir/aux_clm_short.07Oct2013");

my $sexpectedauxclmshort = sortCESMTests(\@expectedauxclmshort);
my $sactualauxclmshort = sortCESMTests($actualauxclmshort);
#print "sorted expected clm short: \n";
#print Dumper $sexpectedauxclmshort;
#print "sorted actual clm short: \n";
#print Dumper $sactualauxclmshort;
#print scalar @$sexpectedauxclmshort . "\n";;
is_deeply($sactualauxclmshort, $sexpectedauxclmshort, "aux_clm_short should parse correctly..");

my @expectedtestmodsdupes;
my $d1 = new CESMTest(compset => 'B', grid => 'f10_f10', testname => 'ERS', machine => 'frankfurt',
                      compiler => 'pgi', testmods => 'clm/default', comment => 'Testmods comment');
push(@expectedtestmodsdupes, $d1);
my $d2 = new CESMTest(compset => 'B', grid => 'f10_f10', testname => 'ERS', machine => 'frankfurt',
                      compiler => 'pgi',  comment => 'No Testmods comment');
push(@expectedtestmodsdupes, $d2);
my $d3 = new CESMTest(compset => 'B', grid => 'f10_f10', testname => 'ERS', machine => 'frankfurt',
                      compiler => 'pgi',  testmods => 'clm/decStart', comment => 'clm decStart');
push(@expectedtestmodsdupes, $d3);
my $d4 = new CESMTest(compset => 'B', grid => 'f10_f10', testname => 'ERS', machine => 'frankfurt',
                      compiler => 'pgi',  testmods => 'clm/Foo1', comment => 'Foo1 comment');
push(@expectedtestmodsdupes, $d4);
my $d5 = new CESMTest(compset => 'B', grid => 'f10_f10', testname => 'ERS', machine => 'frankfurt',
                      compiler => 'pgi',  testmods => 'clm/Foo2', comment => 'Foo2 comment');
push(@expectedtestmodsdupes, $d5);

my $actualtestmodsdupes = parseTextList("$testfilesdir/testmods_dupes_test");

is_deeply($actualtestmodsdupes, \@expectedtestmodsdupes, "We should be able to add duplicate testmods to the same test combination..");

# Now we need to test the removeXML functionality.  

# First, can we remove the aux_clm_short category??

my $testparser = XML::LibXML->new(no_blanks => 1);
my $origxml = $testparser->parse_file("$testfilesdir/testlist.original.15Oct2013.xml");
my $expected_auxclmshort_removed = $testparser->parse_file("$testfilesdir/testlist.auxclmshortremoved.15Oct2013.xml");
my $actual_auxclmshort_removed  = removeXMLTests($origxml, undef, undef, undef, undef, undef, \"aux_clm_short", undef);
my $expectedstring = $expected_auxclmshort_removed->toString(1);
my $actualstring   = $actual_auxclmshort_removed->toString(1);
#open my $ACT, ">", "./auxclm.xml" or die $!;
#print $ACT $actualstring;
#close $ACT;
is_deeply($expectedstring , $actualstring, "We should be able to remove the aux_clm_short category successfully...") or diag explain $actualstring;

# can we remove the X compset?
$origxml = $testparser->parse_file("$testfilesdir/testlist.original.15Oct2013.xml");
my $expected_xcompset_removed  = $testparser->parse_file("$testfilesdir/testlist.xcompsetremoved.15Oct2013.xml");
my $actual_xcompset_removed =    removeXMLTests($origxml, \"X", undef, undef, undef, undef, undef, undef);
$expectedstring = $expected_xcompset_removed->toString(1);
$actualstring = $actual_xcompset_removed->toString(1);
is_deeply($expectedstring , $actualstring, "We should be able to remove X compsets successfully...") or diag explain $actualstring;

# can we remove the 5x5amazon grid?
$origxml = $testparser->parse_file("$testfilesdir/testlist.original.15Oct2013.xml");
my $expected_5by5amazon_removed  = $testparser->parse_file("$testfilesdir/testlist.5by5amazonremoved.15Oct2013.xml");
my $actual_5by5amazon_removed =    removeXMLTests($origxml, undef, \"5x5_amazon", undef, undef, undef, undef, undef);
$expectedstring = $expected_5by5amazon_removed->toString(1);
$actualstring = $actual_5by5amazon_removed->toString(1);
is_deeply($expectedstring , $actualstring, "We should be able to remove the 5x5_amazon grid successfully...") or diag explain $actualstring;

# can we remove the SEQ_PFC test?
$origxml = $testparser->parse_file("$testfilesdir/testlist.original.15Oct2013.xml");
my $expected_seqpfc_removed  = $testparser->parse_file("$testfilesdir/testlist.SEQ_PFCremoved.15Oct2013.xml");
my $actual_seqpfc_removed =    removeXMLTests($origxml, undef, undef, \"SEQ_PFC", undef, undef, undef, undef);
$expectedstring = $expected_seqpfc_removed->toString(1);
$actualstring = $actual_seqpfc_removed->toString(1);
is_deeply($expectedstring , $actualstring, "We should be able to remove the SEQ_PFC test successfully...") or diag explain $actualstring;

# can we remove the olympus machine?
$origxml = $testparser->parse_file("$testfilesdir/testlist.original.15Oct2013.xml");
my $expected_olympus_removed  = $testparser->parse_file("$testfilesdir/testlist.olympusremoved.xml");
my $actual_olympus_removed =    removeXMLTests($origxml, undef, undef, undef, \"olympus", undef, undef, undef);
$expectedstring = $expected_olympus_removed->toString(1);
$actualstring = $actual_olympus_removed->toString(1);
is_deeply($expectedstring , $actualstring, "We should be able to remove the olympus machine successfully...") or diag explain $actualstring;

# can we remove the nag compiler?
$origxml = $testparser->parse_file("$testfilesdir/testlist.original.15Oct2013.xml");
my $expected_nag_removed  = $testparser->parse_file("$testfilesdir/testlist.nagremoved.xml");
my $actual_nag_removed =    removeXMLTests($origxml, undef, undef, undef, undef, \"nag", undef, undef);
$expectedstring = $expected_nag_removed->toString(1);
$actualstring = $actual_nag_removed->toString(1);
is_deeply($expectedstring , $actualstring, "We should be able to remove the nag compiler successfully...") or diag explain $actualstring;

# can we remove the clm decStart testmods?
$origxml = $testparser->parse_file("$testfilesdir/testlist.original.15Oct2013.xml");
my $expected_clmdecStart_removed  = $testparser->parse_file("$testfilesdir/testlist.testmods.clmdecStartremoved.xml");
my $clmdecStart = "clm/decStart";
my $actual_clmdecStart_removed =    removeXMLTests($origxml, undef, undef, undef, undef, undef, undef, \"clm/decStart");
$expectedstring = $expected_clmdecStart_removed->toString(1);
$actualstring = $actual_clmdecStart_removed->toString(1);
is_deeply($expectedstring , $actualstring, "We should be able to remove the clm/decStart testmods successfully...") or diag explain $actualstring;

$origxml = $testparser->parse_file("$testfilesdir/testlist.original.15Oct2013.xml");
my @actualadds;
my $xa1 = new CESMTest(compset => 'XA', grid => 'f10_f10', testname => 'SMS', 
                       machine => 'frankfurt', compiler => 'intel');
my $xa2 = new CESMTest(compset => 'XA', grid => 'f10_f10', testname => 'SMS', 
                       machine => 'frankfurt', compiler => 'pgi');
push(@actualadds, $xa1);
push(@actualadds, $xa2);
my $expected_XAfake_added = $testparser->parse_file("$testfilesdir/testlist.fakecompsetXA.added.xml");
my $actual_XAfake_added = addXMLTests(\@actualadds, $origxml, \"faketype");
$expectedstring = $expected_XAfake_added->toString(1);
$actualstring = $actual_XAfake_added->toString(1);
is_deeply($expectedstring, $actualstring, "adding fake compsets should work..");

my @testmods2add;
my $tm1 = new CESMTest(compset => 'A', grid => 'T31_g37_rx1', testname => 'NCK', machine => 'eos', 
                       compiler => 'intel', testmods => 'clm/default', comment => 'A compset clm default');
my $tm2 = new CESMTest(compset => 'A', grid => 'T31_g37_rx1', testname => 'NCK', machine => 'eos', 
                       compiler => 'intel', testmods => 'clm/ciso', comment => 'A compset clm ciso');
my $tm3 = new CESMTest(compset => 'A', grid => 'T31_g37_rx1', testname => 'NCK', machine => 'eos', 
                       compiler => 'intel', testmods => 'clm/foo1', comment => 'A compset clm foo1');
my $tm4 = new CESMTest(compset => 'A', grid => 'T31_g37_rx1', testname => 'NCK', machine => 'eos', 
                       compiler => 'intel', testmods => 'clm/foo2', comment => 'A compset clm foo2');
push(@testmods2add, $tm1);
push(@testmods2add, $tm2);
push(@testmods2add, $tm3);
push(@testmods2add, $tm4);
$origxml = $testparser->parse_file("$testfilesdir/testlist.original.15Oct2013.xml");
my $expected_testmods_added = $testparser->parse_file("$testfilesdir/testlist.testmodstesting.xml");
my $actual_testmods_added = addXMLTests(\@testmods2add, $origxml, \"prebeta");
$expectedstring = $expected_testmods_added->toString(1);
$actualstring =   $actual_testmods_added->toString(1);
is_deeply($expectedstring, $actualstring, "adding duplicate testmods to the same compset/grid/test/machine/compiler combination should work..");

# Final adds test, let's add 10 eos tests to the list, cloned from titan pgi
my @expeos;
my $e1 = new CESMTest(compset => 'A', grid => 'T31_g37_rx1', testname => 'ERS_D', machine => 'eos', 
					  compiler => 'intel');
my $e2 = new CESMTest(compset => 'A', grid => 'f19_g16_rx1', testname => 'ERS_N2_D', machine => 'eos', 
					  compiler => 'intel');
my $e3 = new CESMTest(compset => 'A', grid => 'f45_g37_rx1', testname => 'NCK', machine => 'eos', 
					  compiler => 'intel');
my $e4 = new CESMTest(compset => 'A', grid => 'ne30_f19_g16_rx1', testname => 'ERS', machine => 'eos', 
					  compiler => 'intel');
my $e5 = new CESMTest(compset => 'B1850BPRP', grid => 'f09_g16', testname => 'CME', machine => 'eos', 
					  compiler => 'intel');
my $e6 = new CESMTest(compset => 'B1850C5', grid => 'f19_g16', testname => 'PET_PT', machine => 'eos', 
					  compiler => 'intel');
my $e7 = new CESMTest(compset => 'B1850C5CN', grid => 'ne30_g16', testname => 'PFS', machine => 'eos', 
					  compiler => 'intel');
my $e8 = new CESMTest(compset => 'B1850CN', grid => 'f09_g16', testname => 'ERI_PT', machine => 'eos', 
					  compiler => 'intel');
my $e9 = new CESMTest(compset => 'B1850CN', grid => 'ne30_g16', testname => 'ERS', machine => 'eos', 
					  compiler => 'intel');
my $e10 = new CESMTest(compset => 'B1850RMCN', grid => 'f19_g16', testname => 'ERS', machine => 'eos', 
					  compiler => 'intel');
push(@expeos, $e1);
push(@expeos, $e2);
push(@expeos, $e3);
push(@expeos, $e4);
push(@expeos, $e5);
push(@expeos, $e6);
push(@expeos, $e7);
push(@expeos, $e8);
push(@expeos, $e9);
push(@expeos, $e10);

$origxml = $testparser->parse_file("$testfilesdir/testlist.original.15Oct2013.xml");
my $expected_eos = $testparser->parse_file("$testfilesdir/testlist.eosadded.xml");
my $actual_eos = addXMLTests(\@expeos, $origxml, \"prebeta");
$expectedstring = $expected_eos->toString(1);
$actualstring = $actual_eos->toString(1);
#open my $ACT, ">", "./eos.xml" or die $!;
#print $ACT $actualstring;
#close $ACT;
is_deeply($expectedstring, $actualstring, "should be able to add eos tests successfully");

# ------------------------------------------------------------------------
# Test convertXMLToCESMTests
# ------------------------------------------------------------------------

my $aux_glc_subset_xml = $testparser->parse_file("$testfilesdir/testlist.aux_glc_subset.xml");
my @expected_aux_glc_subset;
my $tst;
$tst = new CESMTest(compset => 'TGIS2', grid => 'f09_g16_gl10', testname => 'PEA_P1_M_Ly2', machine => 'yellowstone', compiler => 'intel', testtype => 'aux_glc', testmods => 'cism/no_trilinos', comment => 'needs to be no_trilinos because trilinos cannot build with mpi-serial');
push(@expected_aux_glc_subset, $tst);
$tst = new CESMTest(compset => 'TGIS2', grid => 'f09_g16_gl10', testname => 'SMS_D_Ly1', machine => 'yellowstone', compiler => 'intel', testtype => 'aux_glc');
push(@expected_aux_glc_subset, $tst);
$tst = new CESMTest(compset => 'TGIS2', grid => 'f09_g16_gl10', testname => 'SMS_D_Ly1', machine => 'yellowstone', compiler => 'intel', testmods => 'cism/trilinos', testtype => 'aux_glc');
push(@expected_aux_glc_subset, $tst);
$tst = new CESMTest(compset => 'TGIS2', grid => 'f09_g16_gl4', testname => 'SMS_Ly1', machine => 'yellowstone', compiler => 'intel', testmods => 'cism/gradient_margin_2', comment => 'include one short test of the typical production resolution for CISM2', testtype => 'aux_glc');
push(@expected_aux_glc_subset, $tst);

my $actual_aux_glc_subset = convertXMLToCESMTests($aux_glc_subset_xml);

is_deeply(\@expected_aux_glc_subset, $actual_aux_glc_subset, "convertXMLToCESMTests should work correctly for a simple input file");

# ------------------------------------------------------------------------
# Test cleanTestlistXML
# ------------------------------------------------------------------------

# If we read in a mangled xml file (with things out of order and duplicate
# nodes), and run it through the cleanTestlistXML routine, the resulting list of
# CESMTests should be the same as before. i.e., nothing should be added or
# removed. 

my $aux_glc_subset_mangled_xml = $testparser->parse_file("$testfilesdir/testlist.aux_glc_subset_mangled.xml");

my $aux_glc_subset_cleaned_xml = cleanTestlistXML($aux_glc_subset_mangled_xml);

my $cesm_tests_orig = convertXMLToCESMTests($aux_glc_subset_mangled_xml);
$cesm_tests_orig = sortCESMTests($cesm_tests_orig);
my $cesm_tests_new = convertXMLToCESMTests($aux_glc_subset_cleaned_xml);
$cesm_tests_new = sortCESMTests($cesm_tests_new);
is_deeply($cesm_tests_orig, $cesm_tests_new, "cleanTestlistXML should result in the same test list");

# ------------------------------------------------------------------------
# Now, let's do a bit of end-to-end testing:
# Add a test category or two, or compset, or grid, or testmods, then remove them, 
# do we get the original test list?? 
# ------------------------------------------------------------------------

print "$banner\nRUNNING INTEGRATION TESTS\n$banner\n";

# can we remove and re-add the frankfurt nag prealpha list, and will the resulting xml match the original?
my $ffnagprealphaorig = parseTextList("$testfilesdir/frankfurtnagprealpha.orig.txt");
$origxml = $testparser->parse_file("$testfilesdir/testlist.original.15Oct2013.xml");
my $modxml = removeXMLTests($origxml, undef, undef, undef, \"frankfurt", \"nag", \"prealpha", undef);
$modxml = addXMLTests($ffnagprealphaorig, $modxml, \"prealpha");
$origxml = $testparser->parse_file("$testfilesdir/testlist.original.15Oct2013.xml");
$expectedstring = $origxml->toString(1);
$actualstring = $modxml->toString(1);
#open my $ACT, ">", "./prealpha.xml" or die $!;
#print $ACT $actualstring;
#close $ACT;
is_deeply($expectedstring, $actualstring, 
		  "removing/adding frankfurt nag should result in the modified xml matching the original test xml");

# can we remove and re-add the aux_clm_short list, and will the resulting xml match the original?
my $auxclmshortorig = parseTextList("$testfilesdir/aux_clm_short.original.txt");
$origxml = $testparser->parse_file("$testfilesdir/testlist.original.15Oct2013.xml");
$modxml = removeXMLTests($origxml, undef, undef, undef, undef, \"aux_clm_short", undef);
$modxml = addXMLTests($auxclmshortorig, $modxml, \"aux_clm_short");
$origxml = $testparser->parse_file("$testfilesdir/testlist.original.15Oct2013.xml");
$expectedstring = $origxml->toString(1);
$actualstring = $modxml->toString(1);
is_deeply($expectedstring, $actualstring, "removing/adding aux_clm_short should result in the modified xml matching the original..");

# ------------------------------------------------------------------------
# End-to-end testing involving convertXMLToCESMTests
#
# If we start with a complex file, then run it through convertXMLToCESMTests,
# then create a new xml file based on this, we should be back to the original.
# Note that this round-trip workflow is implemented by the cleanTestlistXMl
# function, so that is what we call here.
# ------------------------------------------------------------------------

# In contrast to testlist.original.15Oct2013.xml, this 'cleaned' version has
# gone through this round-trip already, resulting in some cleanup, such as
# consolidating duplicated entries in the xml. (If we tried to simply run this
# test on the non-cleaned version, we would end up with inconsequential
# differences between origxml and modxml.)
$origxml = $testparser->parse_file("$testfilesdir/testlist.original.15Oct2013.cleaned.xml");
$modxml = cleanTestlistXML($origxml);
$expectedstring = $origxml->toString(1);
$actualstring = $modxml->toString(1);
is_deeply($expectedstring, $actualstring, "converting xml to a list of CESMTests then back to xml should match the original");

#done_testing();
