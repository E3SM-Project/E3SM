from __future__ import generators
import abc
import os
import sys
# imports of NCAR scripts
#~#lib_path = os.path.join('scripts/python/contrib/unit_testing')
#~#sys.path.append(lib_path)
sys.path.insert(1, 'scripts/python/contrib/unit_testing')
import environment as lmod
import subprocess
from pprint import pprint

class platformBuilder(object):
    __metaclass__ = abc.ABCMeta
    """ class to extend for building on various platforms. implements
        interfaces and a factory pattern.  creates a relevant platform
        class (darwin_gnu, yellowstone_intel, goldbach_nag, etc...
        that configures cmake, builds pio and tests, and runs the unit
        tests.
    """

    def __init__(self):
        """ user defined ctor so we can put stuff in a class instead of as
            class attributes.  Override this in platform specific classes
        """
        self.setInvariantClassAttr()

        self.CMAKE_EXE = ''

        self.FC = ''
        self.CC = ''
        self.CXX = ''
        self.LDFLAGS = ''

        self.FFLAGS = ''
        self.CXXFLAGS = ''
        self.CFLAGS = '-I. '
        self.OFLAGS = ('-D PIO_BUILD_TESTS:LOGICAL=TRUE '
                                  '-D PIO_BUILD_TIMING:LOGICAL=TRUE ')
        self.MPIEXEC = ''
        self.EXECCA = ''

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
        self.runModuleCmd()
        self.cmakeCmd()
        self.buildCmd()
        self.testCmd()

    def setInvariantClassAttr(self):
        """ figure out some things that shouldn't change in subclasses
        """
        self.className = self.__class__.__name__
        self.BUILD_DIR = "build_" + self.className
        self.TEST_CMD = 'ctest --verbose'
        self.MAKE_CMD = 'make all'
        self.envMod = {}

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
    def factory(type):
        """ factory method for instantiating the appropriate class
        """
        if type == "darwin_gnu":
            return darwin_gnu()
        if type == "goldbach_nag":
            return goldbach_nag()
        if type == "goldbach_intel":
            return goldbach_intel()
        if type == "yellowstone_intel":
            return yellowstone_intel()
        if type == "yellowstone_pgi":
            return yellowstone_pgi()
        if type == "yellowstone_gnu":
            return yellowstone_gnu()
        if type == "caldera_intel":
            return caldera_intel()

        assert 0, "build platform not supported: " + type

""" these subclasses should probably be in their own files.
    each platform then needs one class to extend testCmd and
    runModuleCmd.  That in turn should be extended only for __init__
    for each compiler...all to reduce duplicated code.
"""
    
class darwin_gnu(platformBuilder):

    def __init__(self):
        """ user defined ctor so we can put stuff in a class instead of as
            class attributes
        """
        platformBuilder.__init__(self)
        self.setInvariantClassAttr()

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

class goldbach_nag(platformBuilder):
    
    def __init__(self):
        """ user defined ctor so we can put stuff in a class instead of as
            class attributes
        """
        
        platformBuilder.__init__(self)
        self.setInvariantClassAttr()

        self.moduleList = ['compiler/nag/5.3.1-907']

        self.CMAKE_EXE = '/usr/bin/cmake'
        self.COMPILE_PATH = ('/usr/local/'
                             'openmpi-1.6.5-gcc-g++-4.4.7-3-nag-5.3.1-907'
                             '/bin/')
        self.FC = self.COMPILE_PATH + 'mpif90'
        self.CC = self.COMPILE_PATH + 'mpicc'
        self.CXX = self.COMPILE_PATH + 'mpicxx'
        self.LDFLAGS = '-lcurl'

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
        self.OFLAGS += (' -D CMAKE_VERBOSE_MAKEFILE:BOOL=ON '
                       '-D NETCDF_DIR:STRING=/usr/local/netcdf-gcc-nag '
                       '-D WITH_PNETCDF:LOGICAL=FALSE '
                       '-D PLATFORM:STRING=goldbach ')

        self.MPIEXEC = ('-D MPIEXEC:FILEPATH='
                        '/usr/local/'
                        'openmpi-1.6.5-gcc-g++-4.4.7-3-nag-5.3.1-907'
                        '/bin/mpirun ')
        self.EXECCA = ''

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

        

class goldbach_intel(platformBuilder):

    def __init__(self):
        """ user defined ctor so we can put stuff in a class instead of as
            class attributes
        """
        platformBuilder.__init__(self)
        self.setInvariantClassAttr()

        self.moduleList = ['compiler/intel/14.0.2']

        self.CMAKE_EXE = '/usr/bin/cmake'
        self.COMPILE_PATH = ('/usr/mpi/intel/openmpi-1.4.3-qlc/bin/')
        self.FC = self.COMPILE_PATH + 'mpif90'
        self.CC = self.COMPILE_PATH + 'mpicc'
        self.CXX = self.COMPILE_PATH + 'mpicxx'
        self.LDFLAGS = '-lcurl'

        self.FFLAGS = (' -g '
                       '-DHAVE_MPI -DFORTRANUNDERSCORE '
                       '-DCPRINTEL -DHAVE_SLASHPROC -I. '
                       '-I/usr/local/netcdf-4.3.0-intel-cluster-2013.4.183/include '
                       '-I/usr/mpi/intel/openmpi-1.4.3-qlc/include  ')
        self.CFLAGS = (' -g '
                       ' -DLINUX  '
                       '-DHAVE_MPI -DFORTRANUNDERSCORE '
                       '-DLINUX -DCPRINTEL  '
                       '-DHAVE_SLASHPROC -I. '
                       '-I/usr/mpi/intel/openmpi-1.4.3-qlc/include  ')
        self.CXXFLAGS = (' -g '
                         '-DLINUX  '
                         '-DHAVE_MPI -DFORTRANUNDERSCORE  '
                         ' -DLINUX -DCPRINTEL '
                         '-DHAVE_SLASHPROC -I. '
                         '-I/usr/mpi/intel/openmpi-1.4.3-qlc/include  ')
        self.OFLAGS += (' -D CMAKE_VERBOSE_MAKEFILE:BOOL=ON '
                       '-D NETCDF_DIR:STRING=/usr/local/netcdf-4.3.0-intel-cluster-2013.4.183 '
                       '-D WITH_PNETCDF:LOGICAL=FALSE '
                       '-D PLATFORM:STRING=goldbach ')

        self.MPIEXEC = ('-D MPIEXEC:FILEPATH='
                        '/usr/mpi/intel/openmpi-1.4.3-qlc'
                        '/bin/mpirun ')
        self.EXECCA = ''

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

    def __init__(self):
        """ user defined ctor so we can put stuff in a class instead of as
            class attributes
        """
        platformBuilder.__init__(self)
        self.setInvariantClassAttr()

        self.moduleList = ['ncarenv/1.0 ',
                           'cmake ',
                           'python ',
                           'ncarbinlibs/1.1 ']

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
        self.EXECCA = ('-D EXECCA:FILEPATH="execca " ')

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


class yellowstone_intel(yellowstone):

    def __init__(self):
        """ user defined ctor so we can put stuff in a class instead of as
            class attributes
        """
        yellowstone.__init__(self)
        self.setInvariantClassAttr()

        self.moduleList += ['intel/15.0.0',
                           'ncarcompilers/1.0',
                           'netcdf-mpi/4.3.2',
                           'pnetcdf/1.4.1']

        self.runModuleCmd()
        
        

        self.FFLAGS += (' -fp-model source '
                       '-convert big_endian -assume byterecl -ftz -traceback -assume '
                       'realloc_lhs '
                       '-xHost  -O2   -DLINUX  -DNDEBUG  '
                       '-DHAVE_MPI '
                       '-DFORTRANUNDERSCORE -DNO_R16 -DHAVE_NANOTIME  -DLINUX '
                       '-DCPRINTEL '
                       '-DHAVE_SLASHPROC -I.  ')
        self.CFLAGS += (' -O2 -fp-model precise '
                       '-xHost '
                       '-DLINUX  -DNDEBUG -DHAVE_MPI '
                       '-DFORTRANUNDERSCORE -DNO_R16 -DHAVE_NANOTIME  '
                       '-DLINUX -DBIT64 -DHAVE_VPRINTF  '
                       '-DHAVE_TIMES -DHAVE_GETTIMEOFDAY  '
                       '-DCPRINTEL  -DHAVE_SLASHPROC -I.  ')
        self.CXXFLAGS += (' -O2 -fp-model precise '
                         '-xHost '
                         '  -DNDEBUG -DHAVE_MPI '
                         '-DFORTRANUNDERSCORE -DNO_R16 -DHAVE_NANOTIME  '
                         '-DLINUX '
                         '-DCPRINTEL  -DHAVE_SLASHPROC -I.  ')
        self.OFLAGS += ('-D PNETCDF_DIR:STRING={pnetcdf} '
                       '-D NETCDF_DIR:STRING={netcdf} '.format(
                        netcdf=os.environ['NETCDF'],
                        pnetcdf=os.environ['PNETCDF']))


class caldera_intel(yellowstone_intel):

    def __init__(self):
        """ user defined ctor so we can put stuff in a class instead of as
            class attributes
        """
        yellowstone_intel.__init__(self)
        self.setInvariantClassAttr()
        self.MPIEXEC = ('-D MPIEXEC:FILEPATH="mpirun.lsf " ')
        self.EXECCA = ''

class yellowstone_pgi(yellowstone):

    def __init__(self):
        """ user defined ctor so we can put stuff in a class instead of as
            class attributes
        """
        yellowstone.__init__(self)
        self.setInvariantClassAttr()

        self.moduleList += ['pgi/14.7',
                           'netcdf/4.3.0',
                           'pnetcdf/1.4.1',
                            'ncarcompilers/1.0']

# Need to do this here so that we update the environment
        self.runModuleCmd()

        self.CXX = ''
        self.LDFLAGS = ('-time -Wl,--allow-multiple-definition  -nomp ')


        self.FFLAGS += (' -i4 -gopt -Mlist '
                       '-time -Mextend -byteswapio -Mflushz -Kieee -O '
                       '-nomp   -DLINUX  -DNDEBUG  -DHAVE_MPI '
                       '-DFORTRANUNDERSCORE -DNO_SHR_VMATH -DNO_R16   -DLINUX '
                       '-DCPRPGI  -DHAVE_SLASHPROC -I.  ')
        self.CFLAGS += (' -gopt -Mlist -time  -O  '
                       '-nomp   -DLINUX  -DNDEBUG  -DHAVE_MPI '
                       '-DFORTRANUNDERSCORE -DNO_SHR_VMATH -DNO_R16 '
                       '-DLINUX -DCPRPGI  -DHAVE_SLASHPROC -I.  ')
        self.CXXFLAGS = ''
        self.OFLAGS += ('-D PNETCDF_DIR:STRING={pnetcdf} '
                       '-D NETCDF_DIR:STRING={netcdf} '.format(
                        netcdf=os.environ['NETCDF'],
                        pnetcdf=os.environ['PNETCDF']))



class yellowstone_gnu(yellowstone):

    def __init__(self):
        """ user defined ctor so we can put stuff in a class instead of as
            class attributes
        """
        yellowstone.__init__(self)
        self.setInvariantClassAttr()

        self.moduleList += ['gnu/4.8.0',
                           'ncarcompilers/1.0',
                           'netcdf/4.3.0',
                           'pnetcdf/1.4.1']

        self.FFLAGS += (' -O '
                       '-fconvert=big-endian -ffree-line-length-none '
                       '-ffixed-line-length-none -DLINUX -DNDEBUG '
                       '-D_NETCDF -D_PNETCDF '
                       ' -DHAVE_MPI -DFORTRANUNDERSCORE '
                       ' -DLINUX -DCPRGNU -DHAVE_SLASHPROC -I.  ')
        self.CFLAGS += ('-DLINUX  -DNDEBUG '
                       ' -DHAVE_MPI -DFORTRANUNDERSCORE '
                       '-D_NETCDF -D_PNETCDF '
                       ' -DLINUX -DCPRGNU -DHAVE_SLASHPROC -I.  ')
        self.CXXFLAGS += ('-DLINUX  -DNDEBUG '
                         '-D_NETCDF -D_PNETCDF '
                         ' -DHAVE_MPI -DFORTRANUNDERSCORE '
                         ' -DLINUX -DCPRGNU -DHAVE_SLASHPROC -I.  ')
        self.OFLAGS += ('-D PNETCDF_DIR:STRING={pnetcdf} '
                       '-D NETCDF_DIR:STRING={netcdf} '.format(
                        netcdf=os.environ['NETCDF'],
                        pnetcdf=os.environ['PNETCDF']))

