#!/usr/bin/env python

# A quick hack to submit jobs in a multi-node SMP runtime environment.
# Written in March 2006 by Helgi Adalsteinsson and modified for our
# purpose in April 2006 by Bert Debusschere.
# Modified further by Khachik Sargsyan 2008-09

import os
import thread
import time
import shutil
import sys
import string
import getopt




# Get an array containing names for all of the CPUs we plan on using. We assume
# that there are a given number of cpus available for us in the SMP machine. We
# give them fictitional names cpu0, cpu1, ..., cpun. The number of cpus that is
# specified determines the number of parallel threads we can have going.
def avail_cpus(ncpus = 3):
    """The optional argument specified how many cpus we plan on using."""
    cpus = []
    for icpu in range(ncpus):
       cpus.append('cpu' + str(icpu))
    return cpus

# Use fork/exec to run the given command with arguments,
# returns control when the shell command completes
def run_command(cmd, *args):
    """Given a command and its argument list, use fork/execvp to run the
    command.  Note that this is a direct wrapper for the C library function,
    so the first argument in args should (by convention) be the command name"""
    if len(args) == 0:
        # execvp does not permit an empty tuple as an argument list
        args = ("",)
    # Fork off
    pid = os.fork()
    if pid < 0:
        raise RuntimeError("run_command:  Failed to fork")
    elif pid == 0:
        # I am the child
        os.execvp(cmd, args)
    else:
        # I am the parent
        (retid, status) = os.wait()
    return status

# Use fork/exec to run the given command with arguments in specified directory
# returns control when the shell command completes
def run_command_in_dir(dir, cmd, *args):
    """Given a command and its argument list, use fork/execvp to run the
    command.  Note that this is a direct wrapper for the C library function,
    so the first argument in args should (by convention) be the command name"""
    if len(args) == 0:
        # execvp does not permit an empty tuple as an argument list
        args = ("",)
    # Fork off
    pid = os.fork()
    if pid < 0:
        raise RuntimeError("run_command:  Failed to fork")
    elif pid == 0:
        # I am the child
	os.chdir(dir)
        os.execvp(cmd, args)
    else:
        # I am the parent
        (retid, status) = os.wait()
    return status

# This is what each thread does.
def thread_command(running, cpu, lock, tasks):
    """A task entity for each of the running threads.  All arguments are
    passed by reference (in Python, objects are passed by reference while
    literals are passed by value)."""
    while True:
        lock.acquire_lock()
        if len(tasks) > 0:
            my_task = tasks.pop()
            lock.release_lock()
            #print "Running",my_task,"on cpu",cpu
        else:
            lock.release_lock()
            #print "Task queue on cpu",cpu,"is done"
            break
	# This does not work...
	# Need to first change to the directory where the task script is
	# and then run the script without the path in it...
        #run_command(my_task,my_task)
	dir = my_task[0]
	script = my_task[1]
        #print cpu, dir, script
	starttime = time.time()
	run_command_in_dir(dir,script,script)
	stoptime = time.time()
	print "CPU",cpu,"finished task in",dir,"in",(stoptime-starttime),"seconds"
    #    print "=============================================================================="
    lock.acquire_lock()
    running[0] -= 1
    lock.release_lock()




#######################################################################
def get_tasks(list_args):
    
    param_file='args.in'
    script='./srun.x'
    tasks = []

    #os.system('rm -rf task_*')

    #f = open(param_file,'r')

    for it in range(len(list_args)):

        dir = 'task_' + str(it+1) #str(list_args[it])
        if not os.path.exists(dir):
           print "Creating directory ", dir 
           os.mkdir(dir,0755)


        
        fo = open(dir+os.sep+param_file,'w')
        #print >>fo," ".join(f.readline().split(' '))
        print >>fo,list_args[it]
        fo.close()


        shutil.copy(os.path.dirname(os.path.realpath(__file__))+os.sep+script,dir)


        os.popen("chmod +x " + dir + os.sep + script)
        tasks.append( [dir,script] )

#f.close()

    tasks.reverse()  # to convert 'pop()' into 'shift()'

    return tasks;



# Main routine.




def main(args):
   
   
    
    tasks_file=args[0]
    ncpus=int(args[1])
    #print tasks_file
    print "Running tasks in file", tasks_file,"on", ncpus,"CPUs"



    list_args=open(tasks_file).read().splitlines()
    
    #print list_args
    tasks = get_tasks(list_args)
    
    if (ncpus==1):
        for it in range(len(list_args)):
            dir = 'task_' + str(it+1)
            os.chdir(dir)
            os.system('./srun.x')
            os.chdir('../')
    else:
        running = [0]  # only arrays and lists are passed by reference
        cpus = avail_cpus(ncpus)	# can give optional argument with number of cpus available in the system
        lock = thread.allocate_lock()
        # Make the same number of threads as there are cpus.
        for cpu in cpus:
           lock.acquire_lock()
           running[0] += 1
           lock.release_lock()
           thread.start_new_thread(thread_command, (running, cpu, lock, tasks))
    # Wait for threads to finish (I don't think there is a wait_threads)


        while running[0] > 0:
           time.sleep(2) # sleep for 10 seconds before checking if the threads have finished.
        		   # to avoid spending too much cpu time waiting

    # All done.
        print "All threads have exited"

    
    

if __name__ == "__main__":
    main(sys.argv[1:])

