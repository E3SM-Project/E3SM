#----------------------------------------------------------------------------------
# Class to read in a test list XML file for PTCLM
#----------------------------------------------------------------------------------
from   xml.sax.handler import ContentHandler
class PTCLMtestlistXML( ContentHandler ):

   def startDocument(self):
     self.testlist     = [];
     self.failtestlist = [];

   def startElement(self, name, attrs):

     attributes = [ 'id', 'type', 'site', 'opts', 'resultfile', 'compdir' ]
     testmap    = {}
     for key in attributes:
       testmap[key] = str( attrs.get(key,"") )

     if name == 'test':
       for test in self.testlist:
          if ( testmap['compdir'] != "" and testmap['compdir'] == test['compdir'] ):
             print "compdir is repeated: "+test['compdir']
             sys.exit(100)
          if ( testmap['id'] == test['id'] and testmap['site'] == test['site']):
             print "id and site is duplicated: "+test['id']+" "+test['site']
             sys.exit(100)

       self.testlist.append( testmap )

     if name == 'failtest':
       self.failtestlist.append( testmap )

#----------------------------------------------------------------------------------
# Class to read in a test list for PTCLM and do operations on it
#----------------------------------------------------------------------------------
from   xml.sax         import make_parser
import os
import sys
class PTCLMtestlist:
   # Class data
   testlist     = "PTCLMtestlist.xml"
   testing_dir  = "testing_dir"
   inputdatadir = "/glade/p/cesmdata/cseg/inputdata"
    
   # Construct the class
   def Setup( self, cesmdir ):
     self.testXML = make_parser()
     self.xml     = PTCLMtestlistXML()
     self.testXML.setContentHandler(self.xml)
     # Get the CLM root directory
     self.cesmdir = cesmdir
     if ( not os.path.exists(self.cesmdir) or not os.path.isdir(self.cesmdir) or not os.path.exists(self.cesmdir+"/ChangeLog") ):
         self.error("CESM_ROOT does NOT exist or NOT a directory (set with CESM_ROOT env variable):"+self.cesmdir)

   # Set the input-dir
   def SetDataDir( self, inputdatadir ):
     self.inputdatadir = inputdatadir
     if ( not os.path.exists(self.inputdatadir) or \
          not os.path.isdir(self.inputdatadir)    ):
         self.error("Inputdatadir does NOT exist or NOT a directory:"+self.inputdatadir)

   # Read in the testlist file
   def Read( self ):
     if ( not os.path.exists( self.testlist) ): self.error("File does NOT exist:"+self.testlist)
     print "Open file: "+self.testlist
     self.testXML.parse( self.testlist )

   # --  Error function ---------------------------------
   def error( self, desc ):
       "error function"
       print "ERROR(PTCLMtestlist):: "+desc
       sys.exit(100)

   # Get the list of tests to do
   def get_testlist( self ):
     return( self.xml.testlist )

   # Get the list of fail tests to do
   def get_failtestlist( self ):
     return( self.xml.failtestlist )

   # Get an Identifier name for this test
   def get_testID( self, test ):
     tid = ""
     for att in ["id","opts","site"]:
        val  = test[att].replace(" ","+")
        if ( val != "" ):
           tid  = tid + val + "."

     return(  tid )

   # Fail test or not...
   def IsNOTFailTest( self, test ):
     if ( test['type'] != "Fail" ):
        return( True  )
     else:
        return( False )
     
   # Get the PTCLM command line options to use
   def get_PTCLMoptions( self, test ):
     opts = test['opts']
     opts = opts + " --cesm_root "+self.cesmdir
     if ( self.IsNOTFailTest( test ) ):
        opts = opts + " -s " + test['site']
        opts = opts + " -d "+self.inputdatadir
        opts = opts + " --debug"
        if ( test['type'] == "RUN" ):
           opts = opts + " --mydatadir "+self.testing_dir+"/"+test['id']


     return( opts )

   # run the test
   def run_PTCLMtest( self, test, redo_comp_files ):
     opts = self.get_PTCLMoptions( test )
     tid  = self.get_testID( test )
     tlog = "run.log"
     if ( self.IsNOTFailTest( test ) ):
        errcode =  os.system( "../PTCLMmkdata "+opts+" > "+tlog  )
     else:
        errcode =  os.system( "../PTCLMmkdata "+opts )
     if ( errcode == 0 ):
        stat = True
     else:
        stat = False

     teststatus = ""
     if ( stat == self.IsNOTFailTest( test ) ):
        teststatus = "PASS"
     else:
        teststatus = "Fail"

     print teststatus+" "+tid

     overallcompstatus = "no-comparisons-done"
     if ( test['resultfile'] != "" ):
        dstat = os.system( "diff "+tlog+" "+test['resultfile']  )
        if ( dstat == 0 ):
           compstatus = "PASS"
           desc       = " compare to result file"
        else:
           compstatus = "compare.Fail"
           desc       = " different from result file: "+test['resultfile']

        if ( redo_comp_files ):
           os.system( "cp "+tlog+" "+test['resultfile']  )
        overallcompstatus = compstatus

        print compstatus+" "+tid+" "+desc

     if ( test['compdir'] != "" ):
        testdir  = os.path.abspath(self.testing_dir)+"/"+test['id']
        stdout   = os.popen("cd "+testdir+"/*"+test['site']+"; pwd")
        testdir  = os.path.abspath( stdout.read().rstrip( ) )
        os.system( "mv "+tlog+" "+testdir )
        filelist = ["README.PTCLM", tlog, "user_nl_clm", "xmlchange_cmnds" ]
        for file in filelist:
           srcfile = testdir+"/"+file
           cmpfile = "compdirs/"+test['compdir']+"/"+file
           if ( not os.path.exists( srcfile ) ):
                 compstatus = "compare.Fail"
                 desc       = "source compare file does NOT exist: "+srcfile
           elif ( not os.path.exists( cmpfile ) ):
                 compstatus = "compare.Fail"
                 desc       = "compare file does NOT exist: "+cmpfile
           else:
              dstat = os.system( "diff "+srcfile+" "+cmpfile )

              if ( dstat == 0 ):
                 compstatus = "PASS"
                 desc       = "same as comp directory file: "+cmpfile
              else:
                 compstatus = "compare.Fail"
                 desc       = "different from comp directory file: "+cmpfile

           if ( redo_comp_files ):
              os.system( "cp "+srcfile+" "+cmpfile )
                 
           if ( overallcompstatus != "compare.Fail" ): overallcompstatus = compstatus

           print compstatus+" "+tid+" "+desc+" "+file

     if ( teststatus == "PASS" and overallcompstatus == "PASS" and os.path.exists( tlog ) ): 
        os.system( "/bin/rm "+tlog )

     return( [teststatus, overallcompstatus] )
#
# Unit testing for above classes
#
import unittest

class test_PTCLMtestlist(unittest.TestCase):

   def setUp( self ):
     self.test = PTCLMtestlist()
     cesmdir = os.getenv("CESM_ROOT", "../../../../../../.." )
     self.test.Setup(cesmdir)

   def test_read( self ):
     inpd_orig = self.test.inputdatadir
     self.assertRaises(SystemExit, self.test.SetDataDir, "nodir" )
     self.assertRaises(SystemExit, self.test.SetDataDir, "PTCLMtestlist.py" )
     self.test.Read()
     self.test.SetDataDir( inpd_orig )
     self.test.Read()
     print "\nValid Tests: "
     for test in self.test.get_testlist():
        print str(test)
        self.assertTrue( self.test.IsNOTFailTest( test ) )
     print "\n\n";
     print "\nValid Fail Tests: "
     for test in self.test.get_failtestlist():
        print str(test)
        self.assertTrue( not self.test.IsNOTFailTest( test ) )
     print "\n\n";

   def test_getopts( self ):
     self.test.Read()
     for test in self.test.get_testlist():
       print self.test.get_testID( test )+" = "+ self.test.get_PTCLMoptions( test )

     for test in self.test.get_failtestlist():
       print self.test.get_testID( test )+" = "+ self.test.get_PTCLMoptions( test )

   def test_run( self ):
     self.test.Read()
     testlist     = self.test.get_testlist()
     failtestlist = self.test.get_failtestlist()

     self.test.run_PTCLMtest( testlist[0],     False )
     self.test.run_PTCLMtest( testlist[5],     False )
     self.test.run_PTCLMtest( failtestlist[0], False )

if __name__ == '__main__':
     unittest.main()
