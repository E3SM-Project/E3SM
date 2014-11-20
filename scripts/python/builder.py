from __future__ import generators
import abc
import os
import sys
# imports of NCAR scripts
lib_path = os.path.join('scripts/python/contrib/unit_testing')
sys.path.append(lib_path)
import environment as lmod
import subprocess
import shutil

class platformBuilder(object):
    __metaclass__ = abc.ABCMeta
    """ class to extend for building on various platforms. implements
        interfaces and a factory pattern.  creates a relevant platform
        class (darwin_gnu, yellowstone_intel, goldbach_nag, etc...
        that configures cmake, builds pio and tests, and runs the unit
        tests.
    """

    def __init__(self, compiler,test,mpi,debug):
        """ user defined ctor so we can put stuff in a class instead of as
            class attributes.  Override this in platform specific classes
        """
        self.test = test
        self.CMAKE_EXE = ''
        if mpi is True:
            self.FC = 'mpif90'
            self.CC = 'mpicc'
            self.CXX = 'mpiCC'
        else:
            self.CC = 'cc'
            self.FC = 'f90'
            self.CXX=''

        self.LDFLAGS=''

        if debug is True:
            bldtype = "PIO_DEBUG"
        else:
            bldtype = "PIO"


        self.OFLAGS = ('-D CMAKE_BUILD_TYPE:STRING={0} '
                                  '-D PIO_BUILD_TESTS:LOGICAL=TRUE '
                                  '-D PIO_BUILD_TIMING:LOGICAL=TRUE '.format(bldtype))
        self.MPIEXEC = ''
        self.EXECCA = ''
        self.TEST_CMD = 'ctest --verbose'
        self.MAKE_CMD = 'make all'
        self.envMod = {}

    @classmethod
    def _raise_not_implemented(cls, method_name):
        raise NotImplementedError(cls.__name__ +
                                  " does not implement method " +
                                  method_name+".")

    @abc.abstractmethod
    def runModuleCmd(self, modname):
        """Method not implemented."""
        self._raise_not_implemented("runModuleCmd")

    def metaBuild(self):
        """ routine where everything gets kicked off from
        """
        shutil.rmtree(self.BUILD_DIR, True)
        self.runModuleCmd()
        self.cmakeCmd()
        self.buildCmd()
        if self.test:
          self.testCmd()

    def buildCmd(self):
        """ run build
        """
        p = subprocess.Popen(self.MAKE_CMD,
                             shell=True, env=self.envMod)
        p.wait()

    def testCmd(self):
        """ run tests
        """
        p = subprocess.Popen(self.TEST_CMD,
                             shell=True, env=self.envMod)
        p.wait()

    def cmakeCmd(self):
        """ cmake command to run
        """
        # ~# make build directory and move to it.
        if not os.path.exists(self.BUILD_DIR):
            os.makedirs(self.BUILD_DIR)

        os.chdir(self.BUILD_DIR)

        # ~# change environment, first get existing env
        self.envMod = dict(os.environ)
        # ~# add to env-        
        self.envMod['FC'] = self.FC
        self.envMod['CC'] = self.CC
        if not self.CXX == '':
            self.envMod['CXX'] = self.CXX
        if not self.LDFLAGS == '':
            self.envMod['LDFLAGS'] = self.LDFLAGS

        cmakeString = (self.CMAKE_EXE +' '+ self.OFLAGS + ' '+ self.EXECCA + ' '+self.MPIEXEC + ' ..')
        print(cmakeString)
        p = subprocess.Popen(cmakeString,
                             shell=True, env=self.envMod)
        p.wait()

    @staticmethod
    def factory(platform,compiler,test,mpi,debug):
        """ factory method for instantiating the appropriate class
        """

        if platform == "darwin":
            return darwin(compiler,test,mpi,debug)
        if platform == "goldbach":
            return goldbach(compiler,test,mpi,debug)
        if platform == "yellowstone":
            return yellowstone(compiler,test,mpi,debug)
        if platform == "caldera":
            return caldera(compiler,test,mpi,debug)
        if platform == "mira":
            return cetus(compiler,test,mpi,debug)
        if platform == "cetus":
            return cetus(compiler,test,mpi,debug)
#        return platformBuilder(compiler)


""" these subclasses should probably be in their own files.
    each platform then needs one class to extend testCmd and
    runModuleCmd.  That in turn should be extended only for __init__
    for each compiler...all to reduce duplicated code.
"""
    
class darwin(platformBuilder):

    def __init__(self, compiler, test,mpi,debug):
        """ user defined ctor so we can put stuff in a class instead of as
            class attributes
        """
        super(darwin,self).__init__(compiler, test,mpi,debug)
        
        self.CMAKE_EXE = '/opt/local/bin/cmake'
        self.OFLAGS += ('-D PLATFORM:STRING=darwin ')
        if mpi is True:
            self.MPIEXEC = ('-D  MPIEXEC:FILEPATH='
                            '/opt/local/bin/mpiexec-mpich-gcc48 ')
        self.EXECCA = ''
    
    def runModuleCmd(self):
        """ implement ABC...give pass in this case...run module cmds
        """
        # ~# not implemented for a system without lmod (or
        # ~# somthing similar)
        pass

class elm(darwin):

    def __init__(self, compiler, test,mpi,debug):
        """ user defined ctor so we can put stuff in a class instead of as
            class attributes
        """
        super(elm,self).__init__(compiler, test,mpi,debug)

class goldbach(platformBuilder):
    
    def __init__(self, compiler, test,mpi,debug):
        """ user defined ctor so we can put stuff in a class instead of as
            class attributes
        """
        
        super(goldbach, self).__init__(compiler, test,mpi, debug)
        if compiler == 'nag':
            self.moduleList = ['compiler/nag/5.3.1-907']
        if compiler == 'intel':
            self.moduleList = ['compiler/intel/14.0.2']

        self.BUILD_DIR = "build_goldbach_" + compiler
        self.runModuleCmd()
        
        self.CMAKE_EXE = '/usr/bin/cmake --debug-trycompile'
   
        self.OFLAGS += ( '-D PLATFORM:STRING=goldbach ')
        if mpi is True:
            self.MPIEXEC = ('mpirun ')
        self.EXECCA = ''
        self.LDFLAGS = '-lcurl'
    def runModuleCmd(self):
        """ implement ABC...run module cmds
        """
        # ~# not implemented for a system without lmod (or
        # ~# somthing similar)
        self.lmod = lmod.ModuleInterface()
        self.lmod.python_init("/usr/share/Modules/init/python.py")

        self.lmod.purge()

        for cmd in self.moduleList:
            self.lmod.load(cmd)
        self.lmod.list()
  
class yellowstone(platformBuilder):

    def __init__(self, compiler, test, mpi, debug):
        """ user defined ctor so we can put stuff in a class instead of as
            class attributes
        """
        super(yellowstone,self).__init__( compiler, test, mpi,debug)

        self.moduleList = ['ncarenv/1.0 ',
                                           'cmake ',
                                           'python ',
                           'ncarbinlibs/1.1 ']
        if compiler == 'intel':
            self.moduleList += ['intel/15.0.0',
                           'ncarcompilers/1.0']
            if mpi is True:
                self.moduleList += ['netcdf-mpi/4.3.2']
            else:
                self.moduleList += ['netcdf/4.3.2']
                self.CC = 'icc'
                self.FC = 'ifort'

        if compiler == 'pgi':
            self.moduleList += ['pgi/14.7',
                           'netcdf/4.3.0',
                            'ncarcompilers/1.0']

        if compiler == 'gnu':
            self.moduleList += ['gnu/4.8.0',
                           'ncarcompilers/1.0',
                           'netcdf/4.3.0']

        if mpi is True:
            self.moduleList += ['pnetcdf/1.4.1']
            self.FC = 'mpif90'
            self.CXX = 'mpiCC'

        self.BUILD_DIR = "build_yellowstone_" + compiler
        self.runModuleCmd()

        self.CMAKE_EXE = 'cmake'

        self.NUMPE = '4'

        self.OFLAGS += ('-D PLATFORM:STRING=yellowstone ')

        self.MPIEXEC = ('-D MPIEXEC:FILEPATH="mpirun.lsf " ')
        self.TEST_CMD = ('execca ctest --verbose')

    def testCmd(self):
        """ override testCmd s.t. on yellowstone we open a caldera interactive
            node, run the tests (same as the base class)
            and then exit the queue.
        """
        self.envMod['DAV_CORES'] = self.NUMPE

        p = subprocess.Popen(self.TEST_CMD,
                             shell=True, env=self.envMod)
        p.wait()

    def runModuleCmd(self):
        """ implement ABC...add the lmod commands for yellowstone
        """
        self.lmod = lmod.ModuleInterface()
        self.lmod.python_init("/glade/apps/opt/lmod/lmod/init/"
                              "env_modules_python.py")
        self.lmod.purge()

        for cmd in self.moduleList:
            self.lmod.load(cmd)


class caldera(yellowstone):
    def __init__(self, compiler, test,mpi,debug):
        """ user defined ctor so we can put stuff in a class instead of as
        class attributes
        """
        self.test = test
        super(caldera,self).__init__(compiler, test,mpi,debug)
        self.EXECCA = ''
        self.TEST_CMD = ('ctest --verbose')


class cetus(platformBuilder):

    def __init__(self, compiler, test, mpi, debug):
        """ user defined ctor so we can put stuff in a class instead of as
            class attributes
        """
        super(cetus,self).__init__( compiler, test, mpi, debug)

        self.moduleList = ['+mpiwrapper-xl ',
                           '@ibm-compilers-2014-02 ',
                           '+cmake ']

        self.BUILD_DIR = "build_cetus_" + compiler
        self.runModuleCmd()

        self.CMAKE_EXE = 'cmake '

        self.FC = 'mpixlf2003_r'
        self.CC = 'mpixlc_r'
        self.CXX = 'mpixlcxx_r'

        self.LDFLAGS = '-L/soft/libraries/netcdf/4.3.0-f4.2/cnk-xl/V1R2M0-20131211/lib -L/soft/libraries/hdf5/1.8.10/cnk-xl/V1R2M0-20130405/lib -L/soft/libraries/alcf/current/xl/ZLIB/lib -lnetcdf -lhdf5_hl -lhdf5 -lm -lz'

        self.NUMPE = '4'

        self.OFLAGS += (' -D PLATFORM:STRING=cetus ')

        """ qsub on Cetus does not allow specifying scripts or executables
            when submitting interactive jobs. So we only have two options,
            1. Open an interactive shell and ask user to type "ctest" from
              the interactive shell
            2. Submit the test as a script job
            Option 1 is not feasible since we don't know how long the user 
            might have to wait for the interactive command prompt. So we
            go with Option 2
        """
        self.TEST_CMD = ('qsub -t 15 -n 1 --mode script ../scripts/cetus_test.sh')
        self.MAKE_CMD = ("/bin/csh -c \"" + "source ../scripts/cetus_env.sh && " +
                          "make all " + "\"")

    def buildCmd(self):
        """ run build
        """
        p = subprocess.Popen(self.MAKE_CMD,
                             shell=True, env=self.envMod)
        p.wait()

    def testCmd(self):

        p = subprocess.Popen(self.TEST_CMD,
                             shell=True, env=self.envMod)
        p.wait()

    def runModuleCmd(self):
        """ Using the soft environment requires sourcing the soft
            environment script file - this is not possible right now
            with the current framework. So the next best option is
            to source an environment file everytime we run a command
            (cmake, make) - cetus_env.sh
            This module sources the cetus environment script and grabs
            the PNETCDF and NETCDF installation directories set by
            the script. Note that there are no soft environments for
            these libraries on Cetus
        """
        cmd = ("/bin/csh -c \"" + "source ./scripts/cetus_env.sh && " +
                          "env " + "\"")
        p = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE)
        p.wait()
        for line in p.stdout:
          line = line.decode("UTF-8")
          (key, _, val) = line.partition("=")
          if key == "PNETCDF" or key == "NETCDF":
            os.environ[key] = val.rstrip()
          else:
            pass
				
    def cmakeCmd(self):
        """ cmake command to run
            For cetus the cetus environment script, cetus_env.sh,
            is sourced before running cmake. Overriding this function
            makes this workflow easier.
        """
        # ~# make build directory and move to it.
        if not os.path.exists(self.BUILD_DIR):
            os.makedirs(self.BUILD_DIR)

        os.chdir(self.BUILD_DIR)

        # ~# change environemnt, first get existing env
        self.envMod = dict(os.environ)
        # ~# add to env
        self.envMod['FC'] = self.FC
        self.envMod['CC'] = self.CC
        self.envMod['CXX'] = self.CXX

        cmakeString = (self.CMAKE_EXE + self.OFLAGS + " ..")
        cmakeString = ("/bin/csh -c \"" + "source ../scripts/cetus_env.sh && " +
                        cmakeString.replace('"', r'\"') +
                       "\"")
#        print(cmakeString)

        p = subprocess.Popen(cmakeString,
                             shell=True, env=self.envMod)
        p.wait()

