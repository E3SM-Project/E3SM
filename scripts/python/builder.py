from __future__ import generators
from __future__ import print_function
import abc
import os
import sys
import fileinput

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

    def __init__(self, compiler,test,mpilib,debug):
        """ user defined ctor so we can put stuff in a class instead of as
            class attributes.  Override this in platform specific classes
        """
        self.test = test
        self.CMAKE_EXE = ''
        if mpilib == 'mpi-serial':
            self.CC = 'cc'
            self.FC = 'f90'
            self.CXX=''
        else:
            self.FC = 'mpif90'
            self.CC = 'mpicc'
            self.CXX = 'mpiCC'

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
        self.TEST_CMD = 'ctest '
        self.MAKE_CMD = 'make all'
        self.envMod = dict(os.environ)
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
        # ~# change environment, first get existing env

        # ~# add to env-        
        self.envMod['FC'] = self.FC
        self.envMod['CC'] = self.CC
        if not self.CXX == '':
            self.envMod['CXX'] = self.CXX
        if not self.LDFLAGS == '':
            self.envMod['LDFLAGS'] = self.LDFLAGS

        self.cmakeCmd()
        self.buildCmd()
        if self.test:
          self.testCmd()

    def buildCmd(self):
        """ run build
        """
        p = subprocess.Popen(self.MAKE_CMD,shell=True,env=self.envMod)
        p.wait()

    def testCmd(self):
        """ run tests
        """
        p = subprocess.Popen(self.TEST_CMD, shell=True, env=self.envMod)
        p.wait()

    def cmakeCmd(self):
        """ cmake command to run
        """
        # ~# make build directory and move to it.
        if not os.path.exists(self.BUILD_DIR):
            os.makedirs(self.BUILD_DIR)

        os.chdir(self.BUILD_DIR)

        cmakeString = (self.CMAKE_EXE +' '+ self.OFLAGS + ' '+ self.EXECCA + ' '+self.MPIEXEC + ' ..')

        p = subprocess.Popen(cmakeString,
                             shell=True, env=self.envMod)
        p.wait()

    @staticmethod
    def factory(platform,compiler,test,mpilib,debug):
        """ factory method for instantiating the appropriate class
        """

        if platform == "darwin":
            return darwin(compiler,test,mpilib,debug)
        if platform == "goldbach":
            return goldbach(compiler,test,mpilib,debug)
        if platform == "yellowstone":
            return yellowstone(compiler,test,mpilib,debug)
        if platform == "caldera":
            return caldera(compiler,test,mpilib,debug)
        if platform == "mira":
            return cetus(compiler,test,mpilib,debug)
        if platform == "cetus":
            return cetus(compiler,test,mpilib,debug)
#        return platformBuilder(compiler)


""" these subclasses should probably be in their own files.
    each platform then needs one class to extend testCmd and
    runModuleCmd.  That in turn should be extended only for __init__
    for each compiler...all to reduce duplicated code.
"""
    
class darwin(platformBuilder):

    def __init__(self, compiler, test,mpilib,debug):
        """ user defined ctor so we can put stuff in a class instead of as
            class attributes
        """
        super(darwin,self).__init__(compiler, test,mpilib,debug)
        
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

    def __init__(self, compiler, test,mpilib,debug):
        """ user defined ctor so we can put stuff in a class instead of as
            class attributes
        """
        super(elm,self).__init__(compiler, test,mpilib,debug)

class cray(platformBuilder):
    def __init__(self, compiler, test,mpilib,debug):
        """ user defined ctor so we can put stuff in a class instead of as
            class attributes
        """
        
        super(cray, self).__init__(compiler, test,mpilib, debug)
        
        self.BUILD_DIR = "build_cray_" + compiler
        self.OFLAGS += ( '-D CMAKE_SYSTEM_NAME:STRING=Catamount ')


class bluewaters(cray):

    def __init__(self, compiler, test,mpilib,debug):
        """ user defined ctor so we can put stuff in a class instead of as
            class attributes
        """
        super(bluewaters,self).__init__(compiler, test,mpilib,debug)
        if compiler == 'cray':            
            self.moduleList = ['PrgEnv-cray/5.2.40','cce/8.3.8']
        if compiler == 'pgi':
            self.moduleList = ['PrgEnv-pgi/5.2.40','pgi/14.2.0']
        self.moduleList += ['cray-netcdf-hdf5parallel/4.3.2',
                            'cmake']


class goldbach(platformBuilder):
    
    def __init__(self, compiler, test,mpilib,debug):
        """ user defined ctor so we can put stuff in a class instead of as
            class attributes
        """
        
        super(goldbach, self).__init__(compiler, test,mpilib, debug)
        if compiler == 'nag':
            self.moduleList = ['compiler/nag/5.3.1-907']
        if compiler == 'intel':
            self.moduleList = ['compiler/intel/14.0.2']

        self.BUILD_DIR = "build_goldbach_" + compiler
        self.runModuleCmd()
        
        self.CMAKE_EXE = '/usr/bin/cmake  '
   
        self.OFLAGS += ( '-D PLATFORM:STRING=goldbach ')
        if mpilib is not "mpi-serial":
            self.MPIEXEC = ('mpirun ')
        self.EXECCA = ''
        self.LDFLAGS = '-lcurl'
        os.environ["MODULESHOME"] = "/usr/share/Modules"
        os.environ["MODULEPATH"]="/usr/share/Modules/modulefiles:/etc/modulefiles"
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

    def __init__(self, compiler, test, mpilib, debug):
        """ user defined ctor so we can put stuff in a class instead of as
            class attributes
        """
        super(yellowstone,self).__init__( compiler, test, mpilib,debug)
        os.environ["LMOD_DEFAULT_MODULEPATH"]="/glade/apps/opt/modulefiles/ca/compilers:/glade/apps/opt/modulefiles/ca/idep" 

        self.moduleList = ['ncarenv/1.0 ',
                                           'cmake ',
                                           'python ',
                           'ncarbinlibs/1.1 ']
        if compiler == 'intel':
            self.moduleList += ['intel/15.0.1',
                           'ncarcompilers/1.0']
            if mpilib is not "mpi-serial":
                self.moduleList += ['netcdf-mpi/4.3.3-rc3']
                os.environ["PNETCDF"]="/glade/u/home/jedwards/pnetcdf/svn2013/"
            else:
                self.moduleList += ['netcdf/4.3.2']
                self.CC = 'icc'
                self.FC = 'ifort'

        if compiler == 'pgi':
            self.moduleList += ['pgi/14.10',
                           'netcdf/4.3.0',
                            'ncarcompilers/1.0']
            if mpilib is not "mpi-serial":
                os.environ["PNETCDF"]="/glade/u/home/jedwards/pnetcdf/svn1920/pgi"


        if compiler == 'gnu':
            self.moduleList += ['gnu/4.9.2',
                           'ncarcompilers/1.0',
                           'netcdf/4.3.0']
            if mpilib is not "mpi-serial":
                os.environ["PNETCDF"]="/glade/u/home/jedwards/pnetcdf/svn1920/gnu"

        if mpilib is not 'mpi-serial':
#            self.moduleList += ['pnetcdf/1.4.1']
            self.FC = 'mpif90'
            self.CXX = 'mpiCC'

        self.BUILD_DIR = "build_yellowstone_" + compiler
#        os.environ["LMOD_CMD"]="/glade/apps/opt/lmod/lmod/libexec/lmod"
#        os.environ["LMOD_DIR"]="/glade/apps/opt/lmod/lmod/libexec/"

#        for key in os.environ.keys():
#            print("%30s %s\n" % (key,os.environ[key]))



        self.runModuleCmd()

        self.CMAKE_EXE = 'cmake'

        self.NUMPE = '4'

        self.OFLAGS += ('-D PLATFORM:STRING=yellowstone ')

        self.MPIEXEC = (' -D MPIEXEC:FILEPATH="mpirun.lsf " ')
        #self.TEST_CMD = ('execca ctest --verbose -D Experimental')
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
        self.lmod.python_init("/glade/apps/opt/lmod/lmod/init/env_modules_python.py")
        self.lmod.purge()

        for cmd in self.moduleList:
            print("Loading module "+cmd)
            self.lmod.load(cmd)
        self.lmod.list()


class caldera(yellowstone):
    def __init__(self, compiler, test,mpilib,debug):
        """ user defined ctor so we can put stuff in a class instead of as
        class attributes
        """
        self.test = test
        super(caldera,self).__init__(compiler, test,mpilib,debug)
        self.EXECCA = ''
        #self.TEST_CMD = ('ctest --verbose -D Experimental')
        self.TEST_CMD = ('execca ctest --verbose ')
#        os.environ["LMOD_DEFAULT_MODULEPATH"]="/glade/apps/opt/modulefiles/ca/compilers:/glade/apps/opt/modulefiles/ca/idep" 


class cetus(platformBuilder):

    def __init__(self, compiler, test, mpilib, debug):
        """ user defined ctor so we can put stuff in a class instead of as
            class attributes
        """
        super(cetus,self).__init__( compiler, test, mpilib, debug)

        self.moduleList = ['+mpiwrapper-xl ',
                           '@ibm-compilers-2015-02 ',
                           '+cmake ']

        self.srcroot = os.getcwd()
        self.BUILD_DIR = self.srcroot+"/build_cetus_" + compiler

        self.CMAKE_EXE = 'cmake '

        self.FC = '/home/pkcoff/mpich-sandboxes/onesidedromio/install-gpfsbgq-xl/bin/mpixlf2003_r'
        self.CC = '/home/pkcoff/mpich-sandboxes/onesidedromio/install-gpfsbgq-xl/bin/mpixlc_r'
        self.CXX = '/home/pkcoff/mpich-sandboxes/onesidedromio/install-gpfsbgq-xl/bin/mpixlcxx'
        self.LDFLAGS = '-Wl,--relax -Wl,--allow-multiple-definition -Wl,--whole-archive -L/soft/libraries/hdf5/1.8.14/cnk-xl/V1R2M2-20150213/lib -lhdf5_hl -lhdf5 -L /soft/libraries/alcf/current/xl/ZLIB/lib -lz  -Wl,--no-whole-archive '        
        self.MPIEXEC = (' -D  MPIEXEC:FILEPATH=/usr/bin/runjob')
        MPIEXEC_PREFLAGS =('GPFSMPIO_NAGG_PSET=16:ROMIO_HINTS=/home/pkcoff/public/romio_hints:GPFSMPIO_BALANCECONTIG=1:GPFSMPIO_AGGMETHOD=2:PAMID_TYPED_ONESIDED=1:PAMID_RMA_PENDING=1M:GPFSMPIO_BRIDGERINGAGG=1 ')
#        self.envMod['MPIEXEC_PREFLAGS'] = MPIEXEC_PREFLAGS
        # We use a sh wrapper script on cetus so we need to escape any quotes in the cmake command line    
#        self.OFLAGS += ('-D MPIEXEC_PREFLAGS=\$ENV{MPIEXEC_PREFLAGS}')
        self.NUMPE = ''

        self.OFLAGS += (' -D PLATFORM:STRING=cetus -DCMAKE_C_COMPILER='+self.CC)
        self.OFLAGS += (' -DCMAKE_Fortran_COMPILER='+self.FC)
        self.OFLAGS += (' -DCMAKE_CXX_COMPILER='+self.CXX)
        self.TEST_CMD = ('qsub -o pio2build.out  -t 30 -n 1 --mode script '+self.srcroot+'/scripts/cetus_test.sh ')
        self.MAKE_CMD = ("/bin/sh"+" ./cetus_env.sh"+" make all ")
        self.runModuleCmd()

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
        # ~# make build directory and move to it.
        if not os.path.exists(self.BUILD_DIR):
            os.makedirs(self.BUILD_DIR)

        os.chdir(self.BUILD_DIR)
        f = open("cetus_env.sh", 'w')
        f.write("#!/bin/sh -x\n")
        f.write(". /etc/profile.d/00softenv.sh\n")
        for line in self.moduleList: 
            f.write("soft add "+line+"\n")
        f.write("export LDFLAGS=\""+self.LDFLAGS+"\"\n")
        f.write("echo $@\n")
        f.write("$@\n")
        f.close()
      				
    def cmakeCmd(self):
        """ cmake command to run
            For cetus the cetus environment script, cetus_env.sh,
            is sourced before running cmake. Overriding this function
            makes this workflow easier.
        """

        # ~# change environemnt, first get existing env
        self.envMod = dict(os.environ)
        # ~# add to env
        self.envMod['FC'] = self.FC
        self.envMod['CC'] = self.CC
        self.envMod['CXX'] = self.CXX
#        self.envMod['LDFLAGS'] = self.LDFLAGS

        cmakeString = (self.CMAKE_EXE + self.OFLAGS + self.MPIEXEC+" "+self.srcroot)
        cmakeString = ("/bin/sh"+" ./cetus_env.sh " +
                       cmakeString )

        print(cmakeString)

        p = subprocess.Popen(cmakeString,
                             shell=True, env=self.envMod)
        p.wait()
# The cmake generated CTestTestfile.cmake has incorrect formating 
# I havent found a way to fix it in cmake (where it should be fixed)
# replace all occurrences of '\$' with '$' 

        for i, line in enumerate(fileinput.input(self.BUILD_DIR+'/unittests/CTestTestfile.cmake', inplace=1)):
            sys.stdout.write(line.replace('\$', '$'))  

        for i, line in enumerate(fileinput.input(self.BUILD_DIR+'/test/CTestTestfile.cmake', inplace=1)):
            sys.stdout.write(line.replace('\$', '$'))  




