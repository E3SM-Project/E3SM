"""The run CISM commands.

   A set of functions to run the tests on personal computers (PCs) or high 
   performance computers (HPCs). 
"""

import os
import sys
import time
import subprocess

from util import paths
from util import dicts

def personal(args, cism_driver, data_dir, test_dict):
    """
    Run commands for personal computers (PCs)
    """
    
    test_run = {}
    for case in test_dict:
        case_dir = str.split(case," ")[0]
        case_data_dir = os.path.normpath(data_dir+os.sep+case_dir)
        run_script, mod_dict = test_dict[case]
        
        mod_list = paths.file_modifier_list(args)
        mod_arg = ""
        if mod_list:
            mod_arg = " -m "+str.join("-", mod_list)

        if args.tmod and mod_dict:
            # run the default case that also have performance runs for timing
            case_data_dir += os.sep+'timing'
            paths.mkdir_p(case_data_dir)


        # run default test
        test_commands = ["cd "+os.path.normpath(args.cism_dir+os.sep+'tests'+os.sep+case_dir),
                         "./"+run_script+" -q -e "+cism_driver+" -o "+case_data_dir+mod_arg+' -n 1' ]
        
        if args.tmod:
            test_commands.extend(["cd "+os.path.normpath(case_data_dir),
                'find ./ -not -iname "*.results" -not -iname "*.cism_timing*" -type f -exec rm -f {} \\;', 
                "exit"])
        else:
            test_commands.append("exit")
        
        #print(str.join(" ; ",test_commands))
        if not args.tmod or (args.tmod and mod_dict):
            print("   Spawning "+case+" test default...")
            test_run[case] = subprocess.Popen(str.join(" ; ",test_commands),executable='/bin/bash',shell=True)


        
        # run performance test if specified
        if args.performance and mod_dict:
            for mod in mod_dict:
                print("   Spawning "+case+" test "+mod+"...")
                
                test_commands = ["cd "+os.path.normpath(args.cism_dir+os.sep+'tests'+os.sep+case_dir)+" ",
                                 "./"+run_script+" -q -e "+cism_driver+" -o "+case_data_dir+mod_arg+" "+mod_dict[mod] ]

                if args.tmod:
                    test_commands.extend(["cd "+os.path.normpath(case_data_dir),
                        'find ./ -not -iname "*.results" -not -iname "*.cism_timing*" -type f -exec rm -f {} \\;',
                        "exit"])
                else:
                    test_commands.append("exit")
                
                #print(str.join(" ; ",test_commands))
                test_run[case+' '+mod] = subprocess.Popen(str.join(" ; ",test_commands),executable='/bin/bash',shell=True)

    print("\n   All tests spawned.\n")

    print(  "   Waiting for processes to finish.\n")
    waiting = True
    waited = 0
    while waiting:
        time.sleep(args.sleep)
        waited += args.sleep


        print("\n   Total wait time: "+str(waited)+" seconds.")
        print(  "   Checking processes:")
        running = 0
        for pros in test_run:
            if test_run[pros].poll() is None:
                print("      Still waiting on "+pros)
                running += 1

        if (running == 0):
            waiting = False
            print("      Processes finished.\n")
        else:
            print("      Will check again in "+str(args.sleep)+" seconds...")
    
def create_job(args, job_name, p_replace, run_commands):
    """
    Create the job script for the HPC queues.
    """

    base_job = open('util/job.template','r') 
    job_file = open(job_name, 'w')
    for line in base_job:
        for src, target in p_replace.iteritems():
            line = line.replace(src, target)
        job_file.write(line)

    job_file.write("\n")
    job_file.write("# THE RUN COMMANDS:\n")
    
    for command in run_commands:
        job_file.write(command)

    job_file.write("\nwait \n# FINISH\n")

    base_job.close()
    job_file.close()



def hpc(args, cism_driver, data_dir, test_dict):
    """
    Run commands for high performance computers (HPCs).
    """
   
    # ----------------------
    # Setup all run commands
    # ----------------------
    platform_dict = dicts.hpc_dict[args.platform]
    perf_large_dict = dicts.perf_dict

    jobs_dir = data_dir+os.sep+'jobs'
    paths.mkdir_p(jobs_dir)

    # get file name modifier
    mod_list = paths.file_modifier_list(args)
    mod_arg = ""
    if mod_list:
        mod_arg = " -m "+str.join("-", mod_list)

    timing_exports = set()
    timing_dirs = set()
    timing_commands = []

    for case in test_dict:
        case_dir = str.split(case," ")[0]
        case_data_dir = os.path.normpath(data_dir+os.sep+case_dir)
        cism_test_dir = os.path.normpath(args.cism_dir+os.sep+'tests'+os.sep+case_dir)
        run_script, mod_dict = test_dict[case]
        

        print("   Setting up "+case+" tests")
        test_commands = ["cd "+cism_test_dir,
                "export PYTHONPATH=$PYTHONPATH:"+cism_test_dir,
               "./"+run_script+" -q -e "+cism_driver+" -o "+case_data_dir+mod_arg+" -s -n 1 --hpc",
               "exit"] 
        
        subprocess.check_call(str.join(" ; ",test_commands),executable='/bin/bash',shell=True)
        
        # run performance tests (always do this for hpc systems)
        if mod_dict:
            for mod in mod_dict:
                
                test_commands = ["cd "+cism_test_dir,
                        "export PYTHONPATH=$PYTHONPATH:"+cism_test_dir,
                        "./"+run_script+" -q -e "+cism_driver+" -o "+case_data_dir+mod_arg+" "+mod_dict[mod]+" -s --hpc",
                        "exit"]
                
                subprocess.check_call(str.join(" ; ",test_commands),executable='/bin/bash',shell=True)
        

        # get info to setup timing runs.
        if args.timing and mod_dict:
            timing_dir = case_data_dir+os.sep+'timing'
            for rnd in range(10):
                print("   Setting up "+case+" small timing test "+str(rnd))
                if mod_arg:
                    timing_mod = mod_arg+'-t'+str(rnd)
                else:
                    timing_mod = " -m t"+str(rnd)
                
                timing_dirs.add(timing_dir)
                timing_exports.add("export PYTHONPATH=$PYTHONPATH:"+cism_test_dir)
                timing_commands.extend(["cd "+cism_test_dir, 
                        "./"+run_script+" -q -e "+cism_driver+" -o "+timing_dir+timing_mod+" -s -n 1 --hpc"])

                for mod in mod_dict:
                    timing_commands.extend(["cd "+cism_test_dir, 
                            "./"+run_script+" -q -e "+cism_driver+" -o "+timing_dir+timing_mod+" "+mod_dict[mod]+" -s --hpc"])
       

    # -------------------------
    # Setup the small batch job
    # -------------------------
    # Get the default and small perf run files
    small_run_files = paths.recursive_glob(data_dir,"*.run")
    
    # get the default and small perf run commands
    small_run_commands = []
    for rf in small_run_files:
        rfo = open(rf,'r')
        rfo.next() # skip shebang
        for command in rfo:
            small_run_commands.append(command)

    ## set all aprun commands to background
    #small_run_commands = [command.replace('\n',' & \n') if 'aprun' in command else command for command in small_run_commands ]

    # create the default and small perf job script.
    platform_dict['PBS_N'] = 'small'

    small_job_name = jobs_dir+os.sep+args.platform+'_job.small'
    create_job(args, small_job_name, platform_dict, small_run_commands)
    
    # ----------------------------------------
    # setup the large performance run commands
    # ----------------------------------------
    large_timing_commands = []
    
    for case in perf_large_dict:
        case_dir = str.split(case," ")[0]
        case_data_dir = os.path.normpath(data_dir+os.sep+case_dir)
        cism_test_dir = os.path.normpath(args.cism_dir+os.sep+'tests'+os.sep+case_dir)
        run_script, mod_dict = perf_large_dict[case]

        if mod_dict:
            for mod in mod_dict:
                
                test_commands = ["cd "+cism_test_dir,
                        "export PYTHONPATH=$PYTHONPATH:"+cism_test_dir,
                        "./"+run_script+" -q -e "+cism_driver+" -o "+case_data_dir+mod_arg+" "+mod_dict[mod]+" -s --hpc",
                        "exit"]
               
                subprocess.check_call(str.join(" ; ",test_commands),executable='/bin/bash',shell=True)
        
        # get info to setup timing runs.
        if args.timing and mod_dict:
            timing_dir = case_data_dir+os.sep+'timing'
            for rnd in range(10):
                print("   Setting up "+case+" large timing test "+str(rnd))
                if mod_arg:
                    timing_mod = mod_arg+'-t'+str(rnd)
                else:
                    timing_mod = " -m t"+str(rnd)
                
                timing_dirs.add(timing_dir)
                timing_exports.add("export PYTHONPATH=$PYTHONPATH:"+cism_test_dir)
                large_timing_commands.extend(["cd "+cism_test_dir, 
                        "./"+run_script+" -q -e "+cism_driver+" -o "+timing_dir+timing_mod+" -s -n 1 --hpc"])

                for mod in mod_dict:
                    large_timing_commands.extend(["cd "+cism_test_dir, 
                            "./"+run_script+" -q -e "+cism_driver+" -o "+timing_dir+timing_mod+" "+mod_dict[mod]+" -s --hpc"])
 
    # -------------------------
    # Setup the large batch job
    # -------------------------
    
    # Get large perf run files
    all_run_files = paths.recursive_glob(data_dir,"*.run")
    large_run_files = list( set(small_run_files) ^ set(all_run_files) )   # get the new run files

    # get the  large perf run commands
    large_run_commands = []
    for rf in large_run_files:
        rfo = open(rf,'r')
        rfo.next() # skip shebang
        for command in rfo:
            large_run_commands.append(command)

    # set all aprun commands to background
    #large_run_commands = [command.replace('\n',' & \n') if 'aprun' in command else command for command in large_run_commands ]

    # create the default job script.
    platform_dict['PBS_N'] = 'large'
    platform_dict['PBS_walltime'] = '01:00:00'
    if args.platform.lower() == 'hopper':
        platform_dict['RES_NUM'] = str(11*24)
    else:
        platform_dict['RES_NUM'] = '16'
    
    large_job_name = jobs_dir+os.sep+args.platform+'_job.large'
    create_job(args, large_job_name, platform_dict, large_run_commands)

    if args.timing:
        # ----------------------
        # setup small timing job
        # ----------------------
        timing_exports_all = ['# Setup the environment variables \n'] 
        timing_exports_all.extend(timing_exports)
        timing_exports_all.append("\n")
        
        timing_commands.append("exit")
        subprocess.check_call(str.join(" ; ",timing_commands),executable='/bin/bash',shell=True)

        all_run_files = paths.recursive_glob(data_dir,"*.run")
        small_timing_run_files = list( set(small_run_files + large_run_files) ^ set(all_run_files) )   # get the new run files
      
        small_timing_jobs = set()
        for rnd in range(10):
            subset_run_files = [f for f in small_timing_run_files if '-t'+str(rnd) in f]
        
            # get the small timing run commands
            small_timing_run_commands = []
            for rf in subset_run_files:
                rfo = open(rf,'r')
                rfo.next() # skip shebang
                for command in rfo:
                    small_timing_run_commands.append(command)

            # set all aprun commands to background
            #small_timing_run_commands = [command.replace('\n',' & \n') if 'aprun' in command else command for command in small_timing_run_commands ]

            # create the default job script.
            platform_dict['PBS_N'] = 'small_timing_'+str(rnd)
            platform_dict['PBS_walltime'] = '1:00:00'
            if args.platform.lower() == 'hopper':
                platform_dict['RES_NUM'] = str(1*24)
            else:
                platform_dict['RES_NUM'] = '1'
            
            small_timing_job_name = jobs_dir+os.sep+args.platform+'_job.small_timing_'+str(rnd)
            create_job(args, small_timing_job_name, platform_dict, small_timing_run_commands)
            small_timing_jobs.add(small_timing_job_name)

        # ----------------------
        # setup large timing job
        # ----------------------
        large_timing_commands.append("exit")
        subprocess.check_call(str.join(" ; ",large_timing_commands),executable='/bin/bash',shell=True)

        all_run_files = paths.recursive_glob(data_dir,"*.run")
        large_timing_run_files = list( set(small_run_files + large_run_files + small_timing_run_files) ^ set(all_run_files) )   # get the new run files
       
        large_timing_jobs = set()
        for rnd in range(10):
            subset_run_files = [f for f in large_timing_run_files if '-t'+str(rnd) in f]

            # get the large timing run commands
            large_timing_run_commands = []
            for rf in subset_run_files:
                rfo = open(rf,'r')
                rfo.next() # skip shebang
                for command in rfo:
                    large_timing_run_commands.append(command)

            # set all aprun commands to background
            #large_timing_run_commands = [command.replace('\n',' & \n') if 'aprun' in command else command for command in large_timing_run_commands ]

            # create the default job script.
            platform_dict['PBS_N'] = 'large_timing_'+str(rnd)
            if args.platform.lower() == 'hopper':
                platform_dict['PBS_walltime'] = '20:00'
                platform_dict['RES_NUM'] = str(11*24)
            else:
                platform_dict['PBS_walltime'] = '00:20:00'
                platform_dict['RES_NUM'] = '16'
            
            large_timing_job_name = jobs_dir+os.sep+args.platform+'_job.large_timing_'+str(rnd)
            create_job(args, large_timing_job_name, platform_dict, large_timing_run_commands)
            large_timing_jobs.add(large_timing_job_name)


        # create a script to submit all batch jobs
        sub_script_script = data_dir+os.sep+"submit_all_jobs.bash" 
        sub_script_file = open(sub_script_script,'w') 
        
        sub_script_file.write('#!/bin/bash \n \n')
        sub_script_file.write('qsub '+small_job_name+'\n \n')
        sub_script_file.write('qsub '+large_job_name+'\n \n')
        if args.timing:
            for sm_jb in small_timing_jobs:
                sub_script_file.write('qsub '+sm_jb+'\n \n')
            for lg_jb in large_timing_jobs:
                sub_script_file.write('qsub '+lg_jb+'\n \n')

        sub_script_file.close()
        os.chmod(sub_script_script, 0o755)   # uses an octal number!
        
        
        # create a script to clean out the timing directory.
        clean_script = data_dir+os.sep+"clean_timing.bash" 
        clean_file = open(clean_script,'w') 
        
        clean_file.write('#!/bin/bash \n')
        for dr in timing_dirs:
            clean_file.write("cd "+dr+" \n")
            clean_file.write('find ./ -not -iname "*.results" -not -iname "*.cism_timing*" -type f -exec rm -f {} \\; \n')
            clean_file.write(" \n")

        clean_file.close()
        os.chmod(clean_script, 0o755)   # uses an octal number!


    # -----
    # DONE!
    # -----
    print("\n   Created batch job scripts:")
    print(  "      "+small_job_name)
    if args.timing:
        for sm_jb in small_timing_jobs:
            print(  "      "+sm_jb)
    
    print("\n      "+large_job_name)
    if args.timing:
        for lg_jb in large_timing_jobs:
            print(  "      "+lg_jb)
        print("\n   Submit all jobs with this script:")
        print(  "      "+sub_script_script)
        print("\n   Created script to clean out timing directory:")
        print(  "      "+clean_script)
        print("\n      Run this script after ALL jobs finish to remove every unneeded file in the timing directories.")

