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

    def __init__(self, compiler,test):
        """ user defined ctor so we can put stuff in a class instead of as
            class attributes.  Override this in platform specific classes
        """
        self.test = test
        self.CMAKE_EXE = ''

        self.FC = ''
        self.CC = ''
        self.CXX = ''
        self.LDFLAGS = ''


        self.CXXFLAGS = ''
        self.CFLAGS = '-I. '
        self.OFLAGS = ('-D PIO_BUILD_TESTS:LOGICAL=TRUE '
                                  '-D PIO_BUILD_TIMING:LOGICAL=TRUE ')
        self.MPIEXEC = ''
        self.EXECCA = ''
        self.TEST_CMD = 'ctest --verbose'
        self.MAKE_CMD = 'make all'
        self.envMod = {}
        if compiler == 'pgi':
            self.LDFLAGS = ('-time -Wl,--allow-multiple-definition  -nomp ')
            self.FFLAGS = (' -i4 -gopt -Mlist '
                           '-time -Mextend -byteswapio -Mflushz -Kieee -O '
                           '-nomp   -DLINUX  -DNDEBUG  -DHAVE_MPI '
                           '-DFORTRANUNDERSCORE -DNO_SHR_VMATH -DNO_R16   -DLINUX '
                           '-DCPRPGI  -DHAVE_SLASHPROC -I.  ')
            self.CFLAGS = (' -gopt -Mlist -time  -O  '
                           '-nomp   -DLINUX  -DNDEBUG  -DHAVE_MPI '
                           '-DFORTRANUNDERSCORE -DNO_SHR_VMATH -DNO_R16 '
                           '-DLINUX -DCPRPGI  -DHAVE_SLASHPROC -I.  ')
            self.CXXFLAGS = ''
        if compiler == 'intel':
            self.FFLAGS = (' -fp-model source '
                            '-convert big_endian -assume byterecl -ftz -traceback -assume '
                            'realloc_lhs '
                            '-xHost  -O2   -DLINUX  -DNDEBUG  '
                            '-DHAVE_MPI '
                            '-DFORTRANUNDERSCORE -DNO_R16 -DHAVE_NANOTIME  -DLINUX '
                            '-DCPRINTEL '
                            '-DHAVE_SLASHPROC -I.  ')
            self.CFLAGS = (' -O2 -fp-model precise '
                            '-xHost '
                            '-DLINUX  -DNDEBUG -DHAVE_MPI '
                            '-DFORTRANUNDERSCORE -DNO_R16 -DHAVE_NANOTIME  '
                            '-DLINUX -DBIT64 -DHAVE_VPRINTF  '
                            '-DHAVE_TIMES -DHAVE_GETTIMEOFDAY  '
                            '-DCPRINTEL  -DHAVE_SLASHPROC -I.  ')
            self.CXXFLAGS = (' -O2 -fp-model precise '
                              '-xHost '
                              '  -DNDEBUG -DHAVE_MPI '
                              '-DFORTRANUNDERSCORE -DNO_R16 -DHAVE_NANOTIME  '
                              '-DLINUX '
                              '-DCPRINTEL  -DHAVE_SLASHPROC -I.  ')
        if compiler == 'gnu':
            self.FFLAGS = (' -O '
                           '-fconvert=big-endian -ffree-line-length-none '
                           '-ffixed-line-length-none -DLINUX -DNDEBUG '
                           '-D_NETCDF -D_PNETCDF '
                           ' -DHAVE_MPI -DFORTRANUNDERSCORE '
                           ' -DLINUX -DCPRGNU -DHAVE_SLASHPROC -I.  ')
            self.CFLAGS = ('-DLINUX  -DNDEBUG '
                           ' -DHAVE_MPI -DFORTRANUNDERSCORE '
                           '-D_NETCDF -D_PNETCDF '
                           ' -DLINUX -DCPRGNU -DHAVE_SLASHPROC -I.  ')
            self.CXXFLAGS = ('-DLINUX  -DNDEBUG '
                             '-D_NETCDF -D_PNETCDF '
                             ' -DHAVE_MPI -DFORTRANUNDERSCORE '
                             ' -DLINUX -DCPRGNU -DHAVE_SLASHPROC -I.  ')
        if compiler == 'nag':
            self.FFLAGS = ('-Wp,-macro=no_com '
                           '-kind=byte -wmismatch=mpi_send,mpi_recv,mpi_bcast,'
                           'mpi_allreduce,mpi_reduce,mpi_isend,mpi_irecv,'
                           'mpi_irsend,mpi_rsend,mpi_gatherv,mpi_gather,'
                           'mpi_scatterv,'
                           'mpi_allgather,mpi_alltoallv,mpi_file_read_all,'
                           'mpi_file_write_all,mpibcast,mpiscatterv '
                           '-convert=BIG_ENDIAN '
                           '-gline  -C=all -g -time -f2003 -ieee=stop -DLINUX '
                           ' -DHAVE_MPI -DFORTRANUNDERSCORE '
                           '-DNO_CRAY_POINTERS -DNO_SHR_VMATH -DNO_C_SIZEOF '
                           '-DLINUX '
                           '-DCPRNAG -DHAVE_SLASHPROC -I. '
                           '-I/usr/local/netcdf-gcc-nag/include '
                           '-I/usr/local/openmpi-gcc-nag/include  ')
            self.CFLAGS = (' -g -Wl,--as-needed,'
                           '--allow-shlib-undefined -DLINUX  '
                           '-DHAVE_MPI -DFORTRANUNDERSCORE -DNO_CRAY_POINTERS '
                           '-DNO_SHR_VMATH -DNO_C_SIZEOF -DLINUX -DCPRNAG  '
                           '-DHAVE_SLASHPROC -I. '
                           '-I/usr/local/openmpi-gcc-nag/include  ')
            self.CXXFLAGS = (' -g -Wl,--as-needed,'
                             '--allow-shlib-undefined -DLINUX  '
                             '-DHAVE_MPI -DFORTRANUNDERSCORE -DNO_CRAY_POINTERS '
                             '-DNO_SHR_VMATH -DNO_C_SIZEOF -DLINUX -DCPRNAG  '
                             '-DHAVE_SLASHPROC -I. '
                             '-I/usr/local/openmpi-gcc-nag/include  ')
        if compiler == 'ibm':
            self.FFLAGS = ('-g')
            self.CFLAGS = ('-g')
            self.CXXFLAGS = ('-g')

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

        # ~# change environemnt, first get existing env
        self.envMod = dict(os.environ)
        # ~# add to env
        self.envMod['FC'] = self.FC
        self.envMod['CC'] = self.CC
        self.envMod['CXX'] = self.CXX
        self.envMod['LDFLAGS'] = self.LDFLAGS

        fflags = (' -D CMAKE_Fortran_FLAGS:STRING="{0}" '.format(self.FFLAGS))
        cflags = (' -D CMAKE_C_FLAGS:STRING="{0}"  '.format(self.CFLAGS))
        cxxflags =(' -D CMAKE_CXX_FLAGS:STRING="{0}" '.format(self.CXXFLAGS))

        cmakeString = (self.CMAKE_EXE + fflags + cflags +
                       cxxflags + self.OFLAGS + self.EXECCA + self.MPIEXEC + ' ..')

        p = subprocess.Popen(cmakeString,
                             shell=True, env=self.envMod)
        p.wait()

    @staticmethod
    def factory(platform,compiler,test):
        """ factory method for instantiating the appropriate class
        """

#        type = platform + '_' + compiler
        if platform == "darwin":
            return darwin(compiler,test)
        if platform == "goldbach":
            return goldbach(compiler,test)
        if platform == "yellowstone":
            return yellowstone(compiler,test)
        if platform == "caldera":
            return caldera(compiler,test)
        if platform == "mira":
            return cetus(compiler,test)
        if platform == "cetus":
            return cetus(compiler,test)
#        return platformBuilder(compiler)


""" these subclasses should probably be in their own files.
    each platform then needs one class to extend testCmd and
    runModuleCmd.  That in turn should be extended only for __init__
    for each compiler...all to reduce duplicated code.
"""
    
class darwin(platformBuilder):

    def __init__(self, compiler, test):
        """ user defined ctor so we can put stuff in a class instead of as
            class attributes
        """
        super(darwin,self).__init__(compiler, test)
        
        self.CMAKE_EXE = '/opt/local/bin/cmake'
        self.FC = '/opt/local/bin/mpifort-mpich-gcc48'
        self.CC = '/opt/local/bin/mpicc-mpich-mp'
        self.CXX = '/opt/local/bin/mpicxx-mpich-mp'
        self.LDFLAGS = '-lcurl'

        self.FFLAGS = (' -O '
                       '-fconvert=big-endian '
                       '-ffree-line-length-none -ffixed-line-length-none '
                       '-fno-range-check '
                       '-g -Wall  -DDarwin   -DNO_MPI2 '
                       '-DNO_MPIMOD '
                       '-DFORTRANUNDERSCORE -DNO_R16 -DSYSDARWIN  -DDarwin '
                       '-DCPRGNU -I.  ')
        self.CFLAGS += (' -D CMAKE_C_FLAGS:STRING=" -DDarwin  '
                       '-DNO_MPI2 '
                       '-DNO_MPIMOD -DFORTRANUNDERSCORE '
                       '-DNO_R16 -DSYSDARWIN -DDarwin '
                       '-DCPRGNU  -I." ')
        self.CXXFLAGS = (' -D CMAKE_CXX_FLAGS:STRING=" -DDarwin  '
                         '-DNO_MPI2 '
                         '-DNO_MPIMOD -DFORTRANUNDERSCORE '
                         '-DNO_R16 -DSYSDARWIN -DDarwin '
                         '-DCPRGNU  -I." ')
        self.OFLAGS += (' -D CMAKE_VERBOSE_MAKEFILE:BOOL=ON '
                       '-D NETCDF_DIR:STRING=/opt/local '
                       '-D PNETCDF_DIR:STRING=/opt/local '
                       '-D PLATFORM:STRING=darwin ')
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

    def __init__(self, compiler, test):
        """ user defined ctor so we can put stuff in a class instead of as
            class attributes
        """
        super(elm,self).__init__(compiler, test)

class goldbach(platformBuilder):
    
    def __init__(self, compiler, test):
        """ user defined ctor so we can put stuff in a class instead of as
            class attributes
        """
        
        super(goldbach, self).__init__(compiler, test)
        if compiler == 'nag':
            self.moduleList = ['compiler/nag/5.3.1-907']
        if compiler == 'intel':
            self.moduleList = ['compiler/intel/14.0.2']
            self.FFLAGS += (' -DNO_MPIMOD ')
        self.BUILD_DIR = "build_goldbach_" + compiler
        self.runModuleCmd()

        self.CMAKE_EXE = '/usr/bin/cmake '
        self.FC = 'mpif90'
        self.CC = 'mpicc'
        self.CXX = 'mpicxx'
   
        self.OFLAGS += ('-D NETCDF_DIR:STRING={netcdf} '
                        '-D WITH_PNETCDF:LOGICAL=FALSE '
                        '-D PLATFORM:STRING=goldbach '.format(
                        netcdf=os.environ['NETCDF_PATH']))

        self.MPIEXEC = ('mpirun ')
        self.EXECCA = ''
        self.LDFLAGS = '-lcurl'

    def runModuleCmd(self):
        """ implement ABC...run module cmds
        """
        # ~# not implemented for a system without lmod (or
        # ~# somthing similar)
        self.lmod = lmod.ModuleInterface()
        self.lmod.python_init("scripts/python/contrib/standAlone/"
                              "env_modules_python_goldbach.py")
        self.lmod.purge()

        for cmd in self.moduleList:
            self.lmod.load(cmd)

  
class yellowstone(platformBuilder):

    def __init__(self, compiler, test):
        """ user defined ctor so we can put stuff in a class instead of as
            class attributes
        """
        super(yellowstone,self).__init__( compiler, test)

        self.moduleList = ['ncarenv/1.0 ',
                                           'cmake ',
                                           'python ',
                           'ncarbinlibs/1.1 ']
        if compiler == 'intel':
            self.moduleList += ['intel/15.0.0',
                           'ncarcompilers/1.0',
                           'netcdf-mpi/4.3.2',
                           'pnetcdf/1.4.1']
        if compiler == 'pgi':
            self.moduleList += ['pgi/14.7',
                           'netcdf/4.3.0',
                           'pnetcdf/1.4.1',
                            'ncarcompilers/1.0']
        if compiler == 'gnu':
            self.moduleList += ['gnu/4.8.0',
                           'ncarcompilers/1.0',
                           'netcdf/4.3.0',
                           'pnetcdf/1.4.1']

        self.BUILD_DIR = "build_yellowstone_" + compiler
        self.runModuleCmd()

        self.OFLAGS += ('-D PNETCDF_DIR:STRING={pnetcdf} '
                       '-D NETCDF_DIR:STRING={netcdf} '.format(
                        netcdf=os.environ['NETCDF'],
                        pnetcdf=os.environ['PNETCDF']))

        self.CMAKE_EXE = 'cmake'

        self.FC = 'mpif90'
        self.CC = 'mpicc'
        self.CXX = 'mpicxx'
        self.LDFLAGS = ''
        self.NUMPE = '4'

        self.OFLAGS += (' -D CMAKE_VERBOSE_MAKEFILE:BOOL=ON '
                       '-D PIO_FILESYSTEM_HINTS:STRING=gpfs '
                       '-D PLATFORM:STRING=yellowstone ')
        self.CXXFLAGS += (' -DLINUX -DHAVE_MPI -DFORTRANUNDERSCORE ')

        self.MPIEXEC = ('-D MPIEXEC:FILEPATH="mpirun.lsf " ')
#        self.EXECCA = ('-D EXECCA:FILEPATH="execca " ')
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
    def __init__(self, compiler, test):
        """ user defined ctor so we can put stuff in a class instead of as
        class attributes
        """
        self.test = test
        super(caldera,self).__init__(compiler, test)
        self.EXECCA = ''
        self.TEST_CMD = ('ctest --verbose')


class cetus(platformBuilder):

    def __init__(self, compiler, test):
        """ user defined ctor so we can put stuff in a class instead of as
            class attributes
        """
        super(cetus,self).__init__( compiler, test)

        self.moduleList = ['+mpiwrapper-xl ',
                           '@ibm-compilers-2014-02 ',
                           '+cmake ']

        self.BUILD_DIR = "build_cetus_" + compiler
        self.runModuleCmd()

        self.OFLAGS += ('-D PNETCDF_DIR:STRING={pnetcdf} '
                       '-D NETCDF_DIR:STRING={netcdf} '.format(
                        netcdf=os.environ['NETCDF'],
                        pnetcdf=os.environ['PNETCDF']))

        self.CMAKE_EXE = 'cmake'

        self.FC = 'mpixlf2003_r'
        self.CC = 'mpixlc_r'
        self.CXX = 'mpixlcxx_r'
        self.LDFLAGS = ''
        self.NUMPE = '4'

        self.OFLAGS += (' -D CMAKE_VERBOSE_MAKEFILE:BOOL=ON '
                       '-D PIO_FILESYSTEM_HINTS:STRING=gpfs '
                       '-D PLATFORM:STRING=cetus '
                        '-D ADDITIONAL_LIBS=\"-Wl,--relax -Wl,--allow-multiple-definition -L/soft/libraries/hdf5/1.8.10/cnk-xl/current/lib -lhdf5_hl -lhdf5 -L/soft/libraries/alcf/current/xl/ZLIB/lib -lz\"')
#        self.CXXFLAGS += (' -DLINUX -DHAVE_MPI -DFORTRANUNDERSCORE ')

#        self.MPIEXEC = ('-D MPIEXEC:FILEPATH="runjob " ')
#        self.EXECCA = ('-D EXECCA:FILEPATH="execca " ')
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
				
#        self.lmod = lmod.SoftEnvInterface()

#        for module in self.moduleList:
#            self.lmod.load_str(module)

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
        self.envMod['LDFLAGS'] = self.LDFLAGS

        fflags = (' -D CMAKE_Fortran_FLAGS:STRING="{0}" '.format(self.FFLAGS))
        cflags = (' -D CMAKE_C_FLAGS:STRING="{0}"  '.format(self.CFLAGS))
        cxxflags =(' -D CMAKE_CXX_FLAGS:STRING="{0}" '.format(self.CXXFLAGS))

        cmakeString = (self.CMAKE_EXE + fflags + cflags +
                       cxxflags + self.OFLAGS + " ..")
        cmakeString = ("/bin/csh -c \"" + "source ../scripts/cetus_env.sh && " +
                        cmakeString.replace('"', r'\"') +
                       "\"")
#        print(cmakeString)

        p = subprocess.Popen(cmakeString,
                             shell=True, env=self.envMod)
        p.wait()

