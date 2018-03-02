#!/usr/bin/env python
#
# Magnus Hagdorn
#
# Create test configuration and run test

import time, sys, optparse, os.path, os,fcntl

# test setups
EXCLUDE = ['name']
TEST_CONFIGURATION = {
    'B' : { 'name'   : 'verifBC',
            'lambda' : 0.,
            'H0'     : 3600.,
            'R0'     : 750,
            'tDel'   : 25000.
            },
    'C' : { 'name'   : 'verifBC',
            'lambda' : 5.,
            'H0'     : 3600.,
            'R0'     : 750,
            'tDel'   : -1.
            },
    'D' : { 'name'   : 'verifD',
            'Cp'     : 200.,
            'Tp'     : 5000.,
            'H0'     : 3600.,
            'R0'     : 750,
            'tf'     : 25000.
            }
    }
SOLVERS = {'lin' : 0, 'non-lin' : 2, 'ADI' : 1 }
TOKENS = ['#GRIDSIZE#', '#SOLVER#', '#DELTATIME#', '#TITLE#', '#OUTNAME#']
MODEL_BINARY = 'verif_glide'

def find_model(path):
    """Find the model.

    path: path to model specified on command line."""

    # first, check if there is a path given on command line
    if path != None:
        model = os.path.join(path,MODEL_BINARY)
        if os.path.exists(model):
            return model
    # second, check if the model is in the PATH environment
    pathe = os.environ['PATH']
    for p in pathe.split(os.path.pathsep):
        model = os.path.join(p,MODEL_BINARY)
        if os.path.exists(model):
            return model
    # third, check if we can find model from commandline used
    comm = sys.argv[0]
    if comm != '':
        comm = os.path.dirname(comm.split()[0])
        model = os.path.join(comm,'..','src',MODEL_BINARY)
        if os.path.exists(model):
            return model
    raise RuntimeError, 'Could not find model binary'

def create_config(experiment, solver, gridsize, timestep, tname):
    """Create configuration file.

    experiment: name test configuration
    solver:     one of 'lin', 'non-lin', 'ADI'
    gridsize:   gridsize in km
    timestep:   timestep in a
    tname:      name of template file"""


    # construct base name
    base_name = 'test_%s_%s_%dkm_%.2fa'%(experiment, solver, gridsize,timestep)
    config_name = '%s.config'%base_name

    tfile = open(tname,'r')
    f = open(config_name,'w')

    # print warning etc
    f.write('# automatically generated on %s\n'%(time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())))

    # write experiment config section
    f.write('[%s]\n'%TEST_CONFIGURATION[experiment]['name'])
    for c in TEST_CONFIGURATION[experiment]:
        if c not in EXCLUDE:
            f.write('%s = %s\n'%(c,TEST_CONFIGURATION[experiment][c]))
    f.write('\n')
                    
    # loop over template
    for l in tfile.readlines():
        for t in TOKENS:
            found = False
            if l.find(t) != -1:
                found = True
                break
        if found:
            if t == '#GRIDSIZE#':
                f.write('%s'% ( l.replace(t,'%f'%(gridsize*1000.)) ))
            elif t=='#SOLVER#':
                f.write('%s'% ( l.replace(t,'%d'%SOLVERS[solver]) ))
            elif t=='#DELTATIME#':
                f.write('%s'% ( l.replace(t,'%f'%timestep) ))
            elif t=='#TITLE#':
                f.write('%s'% ( l.replace(t,'Test %s, %s, %dkm, %.2fa'%(experiment, solver, gridsize,timestep)) ))
            elif t=='#OUTNAME#':
                    f.write('%s'% ( l.replace(t,'%s.nc'%base_name) ))
        else:
            f.write('%s'%l)
    tfile.close()
    f.close()

    return config_name

if __name__ == '__main__':

    # setup options
    parser = optparse.OptionParser(usage = "usage: %prog [options] template_file")
    parser.add_option('--experiment',default='B',metavar='TEST',type='choice',choices=TEST_CONFIGURATION.keys(),help='select test scenario (default: B)')
    parser.add_option('--solver',default='lin',metavar='SOLVER',type='choice',choices=SOLVERS.keys(),help='select solver (default: lin)')
    parser.add_option('--gridsize',default=20,metavar='DELTAX',type='int',help='grid spacing in km (default: 20)')
    parser.add_option('--timestep',default=10.,metavar='DELTAT',type='float',help='time step in a (default: 10)')
    parser.add_option('-f','--file',metavar='SPEC',help='an alternative way of setting up the model. The SPEC string takes the following format: test_EXP_SOLVER_Xkm_Ta where EXP is the experiment (see -e), SOLVER the solver (see -s), X the gridspacing in km (see -g) and T the time step in years (see -t)')
    group = optparse.OptionGroup(parser,"Options used for running the model.")
    group.add_option('--only-configure',action="store_true",default=False,help="only produce model configuration file.")
    group.add_option("-m", "--model", help="name of model binary to be launched", metavar="BINARY")
    group.add_option("-r", "--results", help="name of file where timing info is stored (default: results)", metavar="RESULTS", default="results")
    group.add_option("-s", "--submit-sge",action="store_true",default=False,help="submit job to Sun Grid Engine")
    group.add_option("-o", "--submit-options",default="",help="set additional options for cluster submission")
    parser.add_option_group(group)
    
    (options, args) = parser.parse_args()

    if len(args)!=1:
        parser.error("no template file name given")

    # create configuration file
    if options.file == None:
        configname = create_config(options.experiment, options.solver, options.gridsize, options.timestep, args[0])
    else:
        config = options.file.split('_')
        try:
            exp = config[1]
            solver = config[2]
            dx = int(config[3][:-2])
            dt = float(config[4][:-1])
        except:
            parser.error("Cannot parse SPEC string <%s>"%options.file)
        if exp not in TEST_CONFIGURATION.keys():
            parser.error("No such experiment %s, select one of %s"%(exp, str(TEST_CONFIGURATION.keys())))
        if solver not in SOLVERS.keys():
            parser.error("No such solver %s, select one of %s"%(solver, str(SOLVERS.keys())))
        configname = create_config(exp,solver,dx,dt,args[0])

    if options.only_configure:
        print 'Create configuration file %s.config'%base_name
        sys.exit(0)

    prefix = os.path.abspath(os.path.dirname(sys.argv[0]))
    p = prefix.split(os.sep)[-1]
    sge_script = os.path.abspath(os.path.join(prefix,'..','..','..','scripts'))
    if p == 'bin':
        sge_script = os.path.abspath(os.path.join(prefix,'..','share','glimmer'))
    sge_script = os.path.join(sge_script,'qsub_glide.sh')
    if options.model == None:
        model = os.path.join(prefix,MODEL_BINARY)
    else:
        model = options.model

    if not os.path.isfile(model):
        sys.stderr.write("Cannot find model executable %s"%model)
        sys.exit(0)
    prog = "%s -r %s %s"%(model,options.results,configname)

    if options.submit_sge:
        if not os.path.isfile(sge_script):
            sys.stderr.write("Cannot find model submission script %s\n"%sge_script)
            sys.exit(0)
        prog = "qsub %s %s %s"%(options.submit_options,sge_script,prog)

#    try:
#        retcode = subprocess.call(prog,shell=True)
#        if retcode < 0:
#            sys.stderr.write("glide model %s was terminated by signal %d\n"%(model,-retcode))
#    except OSError, e:
#        sys.stderr.write("Execution failed: %s\n"%e)
    print os.popen(prog,'r').read()


