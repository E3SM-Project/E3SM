#! /usr/bin/env python 
import sys, getopt, os 
import numpy as np 
import Nio 
import time
import utility

#
# Command options 
#
def disp_usage(callType):
    print '\n Sets up CESM cases for either an ensemble of runs or a small (size 3)\n'
    print 'test set (default). Then Use pyCECT utilities to create an ensemble \n'
    print 'summary file or to evaluate the small test set of runs against the ensemble.\n'
    print '  '
    print '  ----------------------------'
    print '   Args for ensemble.py :'
    print '  ----------------------------'
    print '   ensemble.py'
    print '   -h                      : prints out this usage message'


    print 'Required flags:'
    print '  --case <name>    Case name passed on to create_newcase'
    print '  --mach <name>    Machine name passed on to create_newcase'
    print ' ' 
    print 'Optional flags (+ all "--" options to create_newcase): '
    print '  --project <num>  Project number to charge in job scripts'
    if callType == 'single_run.py':
       print '  --pertlim <num>     Run CAM with non-zero pertlim'
   print '  --walltime <hr:mn> Amount of walltime requested (default = 4:30, or 0:10 with --uf enabled)'
   print '  --compiler <name>  Compiler to use (default = same as Machine default) '
   print '  --compset <name>   Compset to use (default = F2000)'
   print '  --res <name>       Resolution to run (default = f19_f19)'
   print '  --uf               Enable ninth time step runs (ultra-fast mode) - otherwise the default is 12-month runs'
   if callType == 'ensemble.py': 
       print '  --nb               Disables auto building the root case of the ensemble.'
       print '  --ns               Disables auto submitting any members of the ensemble.'
       print '  --ensemble <size>  Build the ensemble (instead of building 3 cases with random pertlim values for verification),'
       print '                     and specify the number of ensemble members to generate (e.g.: 151 for annual averages or 350 for ultra-fast mode)'
   else:
       print '  --nb               Disables building (and submitting) the single case.'
       print '  --ns               Disables submitting the single case.'


def process_args_dict(caller, caller_argv):


   # Pull in and analyze the command line arguements
    s='case= mach= project= compiler= compset= res= uf nb ns ensemble= verbose silent test multi-driver pecount= nist= mpilib= pesfile= gridfile= srcroot= output-root= script-root= queue= user-modes-dir= input-dir= pertlim= walltime= h'

    optkeys=s.split()

    try: 
        opts, args = getopt.getopt(caller_argv,"hf:",optkeys)

    except getopt.GetoptError:
        disp_usage(caller)
        sys.exit(2)

    #check for help   
    for opt,arg in opts:
        if opt == '-h':
            disp_usage(caller)
            sys.exit()

    #opts_dict and defaults    
    opts_dict={}
    opts_dict['res']='f19_f19'
    opts_dict['compset']='F2000'
    opts_dict['walltime']='04:30'
    opts_dict['pertlim']= 0
    opts_dict['nb'] = False
    opts_dict['ns'] = False
    opts_dict['uf'] = False
    opts_dict['ensemble'] = 0
    #for create newcase
    opts_dict['verbose'] = False
    opts_dict['silent'] = False
    opts_dict['test'] = False
    opts_dict['multi-driver'] = False


    #opts_dict = utility.getopt_parseconfig(opts, optkeys, caller, opts_dict)
    for opt, arg in opts:
        if opt == '--case':
            opts_dict['case'] = arg
        elif opt == '--mach':
            opts_dict['mach'] = arg
        elif opt == '--project':
            opts_dict['project'] = arg
        elif opt == '--compset':
            opts_dict['compset'] = arg
        elif opt == '--res':
            opts_dict['res'] = arg
        elif opt == '--ensemble':
            opts_dict['ensemble'] = arg
        elif opt == '--mach':
            opts_dict['mach'] = arg
        elif opt == '--pertlim':
            opts_dict['pertlim'] = arg
            if caller == 'ensemble.py':
                print "WARNING: pertlim ignored for ensemble.py."
                opts_dict['pertlim'] = 0
        elif opt == '--project':
            opts_dict['project'] = arg
        elif opt == '--uf':
            opts_dict['uf'] = True
        elif opt == '--nb':
            opts_dict['nb'] = True
        elif opt == '--ns':
            opts_dict['ns'] = True
        elif opt == '--verbose':
            opts_dict['verbose'] = True
        elif opt == '--silent':
            opts_dict['silent'] = True
        elif opt == '--test':
            opts_dict['test'] = True
        elif opt == '--multi-driver':
            opts_dict['multi-driver'] = True        
        elif opt == '--nist':
            opts_dict['nist'] = arg
        elif opt == '--pecount':
            opts_dict['pecount'] = arg     
        elif opt == '--mpilib':
            opts_dict['mpilib'] = arg
        elif opt == '--pesfile':
            opts_dict['pesfile'] = arg     
        elif opt == '--srcroot':
            opts_dict['srcroot'] = arg
        elif opt == '--output-root':
            opts_dict['output-root'] = arg     
        elif opt == '--script-root':
            opts_dict['script-root'] = arg     
        elif opt == '--queue':
            opts_dict['queue'] = arg     
        elif opt == '--input-dir':
            opts_dict['input-dir'] = arg     
        elif opt == '--user-modes-dir':
            opts_dict['user-modes-dir'] = arg     
        elif opt == '--walltime':
            opts_dict['walltime'] = arg     
        
    return opts_dict

def single_case(opts_dict):
    
    #make sure we're in directory with single_run.py and ensemble.py
    root_dir = os.path.basename(__file__)

    #change to rootdir then to scripts

    os.system("cd ../../scripts")
    os.fchdir(fd)




def main(argv):

    caller = 'single_run.py'

    opts_dict = process_args_dict(caller, argv)


