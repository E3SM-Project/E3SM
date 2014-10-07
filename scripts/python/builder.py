from __future__ import generators
import os
import subprocess

module = 'eval `/glade/apps/opt/lmod/lmod/libexec/lmod tcsh !*`'

class platformBuilder(object):
    """ class that implements a factory pattern.  creates a relevant
        platform class (darwin, yellowstone, goldbach, edison, etc... that
        configures cmake, builds pio and tests, and runs the unit tests
    """
    def factory(type):
        """ factory method for instantiating the appropriate class
        """
        if type == "darwin":
            return darwin()
        if type == "yellowstone":
            return yellowstone()
        assert 0, "build platform not supported: " + type
    factory = staticmethod(factory)


class darwin(platformBuilder):

    CMAKE_EXE = '/opt/local/bin/cmake'
    BUILD_DIR = 'build'
    MAKE_CMD = 'make all'
    TEST_CMD = 'ctest'

    FC = '/opt/local/bin/mpifort-mpich-gcc48'
    CC = '/opt/local/bin/mpicc-mpich-mp'
    LDFLAGS = '-lcurl'

    FFLAGS = (' -D CMAKE_Fortran_FLAGS:STRING="-O -fconvert=big-endian '
              '-ffree-line-length-none -ffixed-line-length-none '
              '-fno-range-check '
              '-g -Wall  -DDarwin  -DMCT_INTERFACE -DNO_MPI2 -DNO_MPIMOD '
              '-DFORTRANUNDERSCORE -DNO_R16 -DSYSDARWIN  -DDarwin '
              '-DCPRGNU -I. " ')
    CFLAGS = ('-D CMAKE_C_FLAGS:STRING=" -DDarwin  -DMCT_INTERFACE -DNO_MPI2 '
              '-DNO_MPIMOD -DFORTRANUNDERSCORE -DNO_R16 -DSYSDARWIN  -DDarwin '
              '-DCPRGNU -I. " ')
    OFLAGS = ('-D CMAKE_VERBOSE_MAKEFILE:BOOL=ON -D '
              'NETCDF_DIR:STRING=/opt/local '
              '-D WITH_PNETCDF:LOGICAL=FALSE -D '
              'PIO_BUILD_TESTS:LOGICAL=TRUE ')
    MPIEXEC = ' -D  MPIEXEC:FILEPATH=/opt/local/bin/mpiexec-mpich-gcc48 '

    envMod = {}

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
                       self.OFLAGS + self.MPIEXEC + ' ..')
        p = subprocess.Popen(cmakeString,
                             shell=True, env=self.envMod)
        p.wait()

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


class yellowstone(platformBuilder):

    moduleList = [ 'module load intel/14.0.2',
                   'module load netcdf-mpi/4.3.0',
                   'module load pnetcdf/1.4.1',
                   'module load ncarenv/1.0',
                   'module load ncarbinlibs/1.1',
                   'module load cmake',
                   'module list']

    CMAKE_EXE = 'make'
    BUILD_DIR = 'build'
    MAKE_CMD = 'make all'
    TEST_CMD = 'ctest'
    
    MODULE_PRE  = '[\'/bin/bash\', \'-i\', \'-c\', \''
    MODULE_POST = '\']'

    FC = 'mpif90'
    CC = 'mpicc'
    LDFLAGS = ''

    FFLAGS = (' -D CMAKE_Fortran_FLAGS:STRING="-O -fconvert=big-endian '
              '-ffree-line-length-none -ffixed-line-length-none '
              '-fno-range-check '
              '-g -Wall  -DDarwin  -DMCT_INTERFACE -DNO_MPI2 -DNO_MPIMOD '
              '-DFORTRANUNDERSCORE -DNO_R16 -DSYSDARWIN  -DDarwin '
              '-DCPRGNU -I. " ')
    CFLAGS = ('-D CMAKE_C_FLAGS:STRING=" -DDarwin  -DMCT_INTERFACE -DNO_MPI2 '
              '-DNO_MPIMOD -DFORTRANUNDERSCORE -DNO_R16 -DSYSDARWIN  -DDarwin '
              '-DCPRGNU -I. " ')
    OFLAGS = ('-D CMAKE_VERBOSE_MAKEFILE:BOOL=ON -D '
              'NETCDF_DIR:STRING=/opt/local '
              '-D WITH_PNETCDF:LOGICAL=FALSE -D '
              'PIO_BUILD_TESTS:LOGICAL=TRUE ')
    MPIEXEC = ' -D  MPIEXEC:FILEPATH=/opt/local/bin/mpiexec-mpich-gcc48 '

    envMod = {}

    def runModuleCmd(self):
        """ run module cmds
        """
        for cmd in self.moduleList:
            print cmd
            p = subprocess.Popen(['/bin/bash', '-i', '-c', cmd])
            p.wait()


    def cmakeCmd(self):
        """ cmake command to run
        """
        print("yellostone cmake")
        # ~# make build directory and move to it.
        if not os.path.exists(self.BUILD_DIR):
            os.makedirs(self.BUILD_DIR)
        
        os.chdir(self.BUILD_DIR)
        # ~#
        self.runModuleCmd()

    def buildCmd(self):
        """ run build
        """
        print("yellowstone build")

    def testCmd(self):
        """ run tests
        """
        print("yellowstone ctest")


