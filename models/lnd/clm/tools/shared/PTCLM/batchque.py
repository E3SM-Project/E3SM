#########################################################################################
#
# batchque.py
#
# Python class to handle batch submission of single-processor command-line jobs.
#
#########################################################################################
import os, sys

class batchque:
#----------------------------------------------------------------------------------------
# Class to handle batch queue submission
#----------------------------------------------------------------------------------------
   # Class data
   setup     = False
   mach      = ""
   account   = ""
   submit    = False
   jobscript = ""
   #
   # hash's keyed off the list of machines known
   #
   # Basic options giving queue name, number of processors (1) and walltime
   # yellowstone(LSF):   -n 1=1 task, -R=Number of tasks on node, -q=queue name, -N=, -a=process type, -W=wallclock time
   # edison/hopper(PBS): -l=tasks, processors per node, and wallclock time, -q=queue name, -V=use ALL env variables, 
   #                     -m=mail options (ae send mail on submit and exit)
   # -j oe on edison/hopper and -oo on yellowstone (without -e/-eo means combine stderr and stdout
   opts    = { 'yellowstone':"-n 1 -R 'span[ptile=15]' -q caldera -N -a poe ", \
               'edison'     :"-l nodes=1:ppn=1 -q regular -V -m ae -j oe ", \
               'hopper'     :"-l nodes=1:ppn=1 -q regular -V -m ae -j oe ", \
               'janus'      :"-l nodes=1:ppn=1 -q janus-short -V -m ae -j oe " }
   # batch submission command
   bsub      = { 'yellowstone':"bsub",   'edison':"qsub", 'hopper':"qsub", 'janus':"qsub" }
   # Option to give file for standard output
   bs_stdout = { 'yellowstone':" -oo ",  'edison':" -o ", 'hopper':" -o ", 'janus':" -o "   }
   # Option to give job name to use
   bs_jobnam = { 'yellowstone':" -J ",   'edison':" -N ", 'hopper':" -N ", 'janus':" -N "   }
   # Option to give current directory to use
   bs_curdir = { 'yellowstone':" -cwd ", 'edison':" -d ", 'hopper':" -d ", 'janus':" -d "   }
   # Option to give account name to use
   bs_accnt  = { 'yellowstone':" -P ",   'edison':"",     'hopper':"",     'janus':""     }
   # If jobcommand needs to be script file
   bs_script = { 'yellowstone':False,    'edison':True,   'hopper':True,   'janus':True   }
   # Option to give wallclock time to use
   bs_wtime  = { 'yellowstone':" -W ",   'edison':" -l walltime=", \
                 'hopper':" -l walltime=", 'janus':" -l walltime="     }

   def Initialize( self, prog, mach="yellowstone", account="" ):
      "Initialize the batchque"
      if ( not self.bsub.has_key(mach) ):
         print "List of valid machines: "+str(self.bsub.keys())
         prog.error( "Machine NOT in list of valid machines for batch queue: "+mach )

      self.mach    = mach
      if ( self.bs_accnt[mach] == "" and account != "" ):
         prog.error( "Account entered but this machine does NOT have an account option: "+mach )

      self.account = account

      self.setup   = True
      self.submit  = False
 
   def Get_OutFilename( self, prog ):
       "Get the output log filename"
       if ( not self.setup ):
          prog.error( "Trying to get the output filename and Initialize was NOT run first!" )
       if ( not self.submit ):
          prog.error( "Trying to get the output filename and Submit was NOT run first!" )

       return( self.stdout )

   def Submit( self, prog, jobcommand, curdir=os.getcwd(), jobname="batchjob", wall="4:00", submit=True ):
       "Get the command to submit the job to the batch queue"
       if ( not self.setup ):
          prog.error( "Initialize was NOT run first!" )

       if ( not os.path.exists(curdir) ):
          prog.error( "Input current directory does NOT exist: "+curdir )

       cmd    =  self.bsub[self.mach]
       opts   =  ""
       stdout = str(jobname)+".stdout.out"
       self.stdout = stdout
       opts += self.bs_stdout[self.mach]+stdout+" "
       opts += self.bs_jobnam[self.mach]+str(jobname)+" "
       opts += self.bs_curdir[self.mach]+curdir+" "
       opts += self.bs_wtime[self.mach]+wall+" "
       if ( self.account != ""  and self.bs_accnt[self.mach] != "" ):
          opts += self.bs_accnt[self.mach]+self.account+" "
       opts +=  self.opts[self.mach]+" "
       if ( self.bs_script[self.mach] ):
          self.jobscript = jobname+".job"
          if ( os.path.exists( self.jobscript ) ):
             os.system( "/bin/rm -rf "+self.jobscript )
          js = open(self.jobscript,"w")
          js.write( jobcommand+"\n" )
          js.close()
          os.chmod(self.jobscript,0555)
          cmd  += " "+opts+self.jobscript
       else:
          cmd  += " "+opts+jobcommand

       if ( os.path.exists( self.stdout ) ):
         os.system( "/bin/rm "+self.stdout )
       if ( submit ):
          status = os.system( cmd )
          if ( status != 0 ):
             prog.error( "Batch submit returns an error" )

       self.submit = True

       return( cmd )

   def SubmitCleanup( self, prog, rmout=False ):
       "Cleanup any files made in submit and reset output filename -- only DO AFTER BATCH HAS RUN!"
       if ( not self.setup  ):
          prog.error( "Initialize was NOT run first!" )
       if ( not self.submit ):
          prog.error( "Submit was NOT run first!" )
       outfile = self.Get_OutFilename( prog )
       if ( not os.path.exists(outfile) ):
          prog.error( "SubmitCleanup called before batch output was returned" )

       if ( self.bs_script[self.mach] ):
         os.system( "/bin/rm -rf "+self.jobscript )
       if ( rmout ):
         os.system( "/bin/rm "+outfile )

       self.submit = False

#
# Unit testing for above classes
#
import unittest

class error_prog:
     def error( self, desc ):
         print desc
         sys.exit( 100 )

class test_batchque(unittest.TestCase):

   def setUp( self ):
       "Setup tests"
       self.prog = error_prog()
       self.que  = batchque()

   def test_badinit( self ):
       "test bad initialization"
       # Bad machine name
       self.assertRaises(SystemExit, self.que.Initialize, self.prog, mach="zztop" )
       # account given on machine without account
       self.assertRaises(SystemExit, self.que.Initialize, self.prog, mach="edison", account="thing" )
       # test using submit and Get_OutFilename before Initialization
       self.assertRaises(SystemExit, self.que.Submit, self.prog, "ls" )
       self.assertRaises(SystemExit, self.que.Get_OutFilename, self.prog )
       # Test that all hashes have the same list of keys
       keylist = str(self.que.opts.keys())
       self.assertTrue(keylist == str(self.que.bsub.keys())      )
       self.assertTrue(keylist == str(self.que.bs_stdout.keys()) )
       self.assertTrue(keylist == str(self.que.bs_jobnam.keys()) )
       self.assertTrue(keylist == str(self.que.bs_accnt.keys())  )
       self.assertTrue(keylist == str(self.que.bs_curdir.keys())  )

   def test_init( self ):
       "test initialization and submit"

       machlist = self.que.opts.keys()
       for mach in machlist:
          print "Test initialization for: "+mach
          self.que.Initialize( self.prog, mach=mach )
          cmd = self.que.Submit( self.prog, "ls", jobname=mach, submit=False )
          print cmd+"\n"
          outfile = self.que.Get_OutFilename( self.prog )
          os.system( "touch "+outfile )
          self.que.SubmitCleanup( self.prog, rmout=True )

       mach = "yellowstone"
       self.que.Initialize( self.prog, mach=mach, account="account" )
       cmd = self.que.Submit( self.prog, "ls", jobname=mach, submit=False )
       print cmd+"\n"
       outfile = self.que.Get_OutFilename( self.prog )
       os.system( "touch "+outfile )
       print "outfile: "+outfile
       self.que.SubmitCleanup( self.prog, rmout=True )

   def test_bad_submit( self ):
       "test bad submit"
       mach = "yellowstone"
       self.que.Initialize( self.prog, mach=mach, account="account" )
       self.assertRaises(SystemExit, self.que.Submit, self.prog, "ls", curdir="zztop", jobname=mach, submit=False )
       self.assertRaises(SystemExit, self.que.Get_OutFilename, self.prog )

   def test_submit( self ):
       "test submitting to local machine if on list"
       stdout    = os.popen("hostname")
       host      = stdout.read().rstrip( )
       startname = { 'ys':'yellowstone', 'edison':'edison', 'hopper':'hopper', 'login':'janus' }
       mach      = ""
       for sname in startname:
          if ( host.startswith(sname) ):
             mach = startname[sname]
       if ( mach != "" ):
          if ( mach == "yellowstone" ):
             account = "P93300606"
          else:
             account = ""
          self.que.Initialize( self.prog, mach=mach, account=account )
       else:
          print "Machine not known, so NOT trying a test submit"
          return

       print "Submit ls to batch queue"
       # Submit and get the output filename
       scmd = "ls PTCLMmkdata"
       wall = "0:01"
       cmd = self.que.Submit( self.prog, scmd, jobname=mach, wall=wall, submit=False )
       print "submit: "+cmd
       outfile = self.que.Get_OutFilename( self.prog )
       # make sure SubmitCleanup will fail since it wasn't submitted yet and output not returned
       self.assertRaises(SystemExit, self.que.SubmitCleanup, self.prog )
       status = os.system( cmd )
       self.assertTrue( status == 0 )
       iter = 0
       exists = os.path.exists(outfile)
       while( status == 0 and iter < 100 and not exists ):
          iter += 1
          exists = os.path.exists(outfile)
          if ( not exists ):
             print "Sleep for a bit to check if outfile was created yet"
             os.system( "sleep 20" )
          else:
             print "Out file created cat it... should be ls of PTCLMmkdata"
             self.assertTrue( os.path.exists(outfile) )
             os.system( "cat "+outfile )

       # Cleanup after the submital
       self.que.SubmitCleanup( self.prog, rmout=True )

       # make sure files were deleted and submit status changed
       self.assertTrue( not self.que.submit )
       self.assertRaises(SystemExit, self.que.Get_OutFilename, self.prog )
       self.assertTrue( not os.path.exists(outfile) )
       if ( self.que.jobscript != "" ):
          self.assertTrue( not os.path.exists(self.que.jobscript) )
       # Now test that a bad submit returns an error (give bad walltime)
       self.assertRaises(SystemExit, self.que.Submit, self.prog, scmd, jobname=mach, wall="--zztop" )

if __name__ == '__main__':
     unittest.main()
