from __future__ import generators
import os
import subprocess
import environment as lmod
import abc

class platformBuilder(object):
    __metaclass__  = abc.ABCMeta
    """ abstract base class that implements interfaces and
        a factory pattern.  creates a relevant
        platform class (darwin_gnu, yellowstone_intel, goldbach_nag, etc... that
        configures cmake, builds pio and tests, and runs the unit tests
    """
    
    def __init__(self):
        """ user defined ctor so we can put stuff in a class instead of as
            class attributes.  Override this in platform specific classes
        """
        self.setInvariantClassAttr()
    
        self.CMAKE_EXE = ''
    
        self.FC = ''
        self.CC = ''
        self.LDFLAGS = ''
    
        self.FFLAGS = ''
        self.CFLAGS = ''
        self.OFLAGS = ''
        self.QUEUE = ' '
        self.MPIEXEC = ''
    
        self.envMod = {}
    
    @classmethod
    def _raise_not_implemented(cls, method_name):
        raise NotImplementedError(
                                  cls.__name__ +" does not implement method " +
                                  method_name+".")
    
    @abc.abstractmethod
    def runModuleCmd(self, modname):
        """Method not implemented."""
        self._raise_not_implemented("runModuleCmd")
    
    def metaBuild(self):
        """ routine where everything gets kicked off from
        """
        self.runModuleCmd()
        self.cmakeCmd()
        self.buildCmd()
        self.testCmd()
    
    def setInvariantClassAttr(self):
        """ figure out some things that shouldn't change in subclasses
        """
        self.className = self.__class__.__name__
        self.BUILD_DIR = "build_" + self.className
        self.TEST_CMD = 'ctest'
        self.MAKE_CMD = 'make all'
    
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
    
        # ~# change environemnt, first get existing env
        self.envMod = dict(os.environ)
        # ~# add to env
        self.envMod['FC'] = self.FC
        self.envMod['CC'] = self.CC
        self.envMod['LDFLAGS'] = self.LDFLAGS
    
        cmakeString = (self.CMAKE_EXE + self.FFLAGS + self.CFLAGS +
                       self.OFLAGS + self.QUEUE + self.MPIEXEC + ' ..')
        print cmakeString
        p = subprocess.Popen(cmakeString,
                             shell=True, env=self.envMod)
        p.wait()

    def factory(type):
        """ factory method for instantiating the appropriate class
        """
        if type == "darwin_gnu":
            return darwin_gnu()
        if type == "yellowstone_intel":
            return yellowstone_intel()
        assert 0, "build platform not supported: " + type
    factory = staticmethod(factory)


class darwin_gnu(platformBuilder):

    def __init__(self):
        """ user defined ctor so we can put stuff in a class instead of as
            class attributes
        """
        self.setInvariantClassAttr()

        self.CMAKE_EXE = '/opt/local/bin/cmake'
    
        self.FC = '/opt/local/bin/mpifort-mpich-gcc48'
        self.CC = '/opt/local/bin/mpicc-mpich-mp'
        self.LDFLAGS = '-lcurl'
    
        self.FFLAGS = (' -D CMAKE_Fortran_FLAGS:STRING="-O -fconvert=big-endian '
                       '-ffree-line-length-none -ffixed-line-length-none '
                       '-fno-range-check '
                       '-g -Wall  -DDarwin  -DMCT_INTERFACE -DNO_MPI2 -DNO_MPIMOD '
                       '-DFORTRANUNDERSCORE -DNO_R16 -DSYSDARWIN  -DDarwin '
                       '-DCPRGNU -I. " ')
        self.CFLAGS = ('-D CMAKE_C_FLAGS:STRING=" -DDarwin  -DMCT_INTERFACE -DNO_MPI2 '
                       '-DNO_MPIMOD -DFORTRANUNDERSCORE -DNO_R16 -DSYSDARWIN  -DDarwin '
                       '-DCPRGNU -I. " ')
        self.OFLAGS = ('-D CMAKE_VERBOSE_MAKEFILE:BOOL=ON '
                       '-D NETCDF_DIR:STRING=/opt/local '
                       '-D PNETCDF_DIR:STRING=/opt/local '
                       '-D PIO_BUILD_TESTS:LOGICAL=TRUE ')
                        
                        ### '-D WITH_PNETCDF:LOGICAL=FALSE
                                  
        self.QUEUE = (' ')
        self.MPIEXEC = ('-D  MPIEXEC:FILEPATH=/opt/local/bin/mpiexec-mpich-gcc48 ')

    def runModuleCmd(self):
        """ run module cmds
        """
        # ~# not implemented for a system without lmod (or
        # ~# somthing similar)
        pass


class yellowstone_intel(platformBuilder):

    def __init__(self):
        """ user defined ctor so we can put stuff in a class instead of as
            class attributes
        """
        self.setInvariantClassAttr()

        self.moduleList = ['intel/14.0.2',
                           'ncarcompilers/1.0',
                           'netcdf-mpi/4.3.0',
                           'pnetcdf/1.4.1',
                           'ncarenv/1.0',
                           'cmake',
                           'ncarbinlibs/1.1']

        self.CMAKE_EXE = 'cmake'
    
        self.FC = 'mpif90'
        self.CC = 'mpicc'
        self.LDFLAGS = ''

        self.FFLAGS = (' -D CMAKE_Fortran_FLAGS:STRING="-fp-model source '
                       '-convert '
                       'big_endian -assume byterecl -ftz -traceback -assume '
                       'realloc_lhs '
                       '-xHost  -O2   -DLINUX  -DNDEBUG -DMCT_INTERFACE '
                       '-DHAVE_MPI '
                       '-DFORTRANUNDERSCORE -DNO_R16 -DHAVE_NANOTIME  -DLINUX '
                       '-DCPRINTEL '
                       '-DHAVE_SLASHPROC -I. " ')
        self.CFLAGS = (' -D CMAKE_C_FLAGS:STRING="-O2 -fp-model precise -xHost '
                       '-DLINUX  -DNDEBUG -DMCT_INTERFACE -DHAVE_MPI '
                       '-DFORTRANUNDERSCORE -DNO_R16 -DHAVE_NANOTIME  -DLINUX '
                       '-DCPRINTEL  -DHAVE_SLASHPROC -I. " ')
        self.OFLAGS = ('-D CMAKE_VERBOSE_MAKEFILE:BOOL=ON '
                       '-D NETCDF_DIR:STRING='
                       '/glade/apps/opt/netcdf-mpi/4.3.2/intel/default '
                       '-D PIO_FILESYSTEM_HINTS:STRING=gpfs '
                       '-D PIO_BUILD_TESTS:LOGICAL=TRUE ')

        self.QUEUE = ('-D QUEUE:FILEPATH=execca ')
        self.MPIEXEC = ('-D MPIEXEC:FILEPATH=mpirun.lsf ')

    def runModuleCmd(self):
        """ run module cmds
        """
        self.lmod = lmod.ModuleInterface()
        self.lmod.python_init("/glade/apps/opt/lmod/lmod/init/"
                              "env_modules_python.py")
        self.lmod.purge()

        for cmd in self.moduleList:
            self.lmod.load(cmd)
