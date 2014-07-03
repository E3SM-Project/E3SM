#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_INIT(PIO, 1.7.2)
AC_CONFIG_MACRO_DIR([m4])
#
# take input from Makefile.conf.in
# output to Makefile.conf.out
#
AC_CONFIG_FILES(Makefile.conf)
AC_CONFIG_HEADER(config.h)

#
# configure options
#

# example:
#   AC_ARG_ENABLE(myvar,
#   AC_HELP_STRING([--enable-myvar],[turn on myvar]),
#   [USEMYVAR=$enableval],[USEMYVAR=default])
#   AC_SUBST(USEMYVAR)

AC_ARG_ENABLE(pnetcdf,
AC_HELP_STRING( [--enable-pnetcdf],
  [use pnetcdf (Uses PNETCDF_PATH) (default auto detect)]),
  [USE_PNETCDF=$enableval],[USE_PNETCDF=default])

AC_ARG_ENABLE(netcdf,
AC_HELP_STRING( [--enable-netcdf],
[use serial netcdf (Uses NETCDF_PATH) (default yes)]), 
[USE_NETCDF=$enableval],[USE_NETCDF=yes])

AC_ARG_ENABLE(netcdf4,
AC_HELP_STRING( [--enable-netcdf4],
[use netcdf4 with HDF5 capability (Uses NETCDF_PATH) (default no)]),
[USE_NETCDF4=$enableval],[USE_NETCDF4=no])

AC_ARG_ENABLE(mpiio,
AC_HELP_STRING([--enable-mpiio],
[Enable support for MPI-IO (default auto detect)]),
[USE_MPIIO=$enableval],[USE_MPIIO=default])

AC_ARG_ENABLE(debug,
AC_HELP_STRING([--enable-debug],
[Turn off optimization and turn on -g (same as OPT=-g)]),
[DEBUG=$enableval],[DEBUG=no])
AC_SUBST(DEBUG)

AC_ARG_ENABLE(compression,
AC_HELP_STRING( [--enable-compression],
  [use compression (Uses VDC2 compression) (default no)]),
  [USE_COMPRESSION=$enableval],[USE_COMPRESSION=no])

AC_ARG_ENABLE(timing,
AC_HELP_STRING([--enable-timing],
[Turn on the gptl timing library interface]),
[USE_TIMING=$enableval],[USE_TIMING=no])

AC_ARG_ENABLE(mpiserial,
AC_HELP_STRING([--enable-mpiserial],
[Build for use with mpi-serial lib (default no)]),
[USE_MPISERIAL=$enableval],[USE_MPISERIAL=no])

AC_ARG_ENABLE(mpi2,
AC_HELP_STRING([--enable-mpi2],
[Enable use of mpi2 calls (default yes)]),
[USE_MPI2=$enableval],[USE_MPI2=yes])

AC_ARG_ENABLE(c_sizeof,
AC_HELP_STRING([--enable-c_sizeof],
[Enable use of c_sizeof() in Fortran (default auto detect)]),
[FC_USE_C_SIZEOF=$enableval],[FC_USE_C_SIZEOF=default])

AC_ARG_ENABLE(memmon,
AC_HELP_STRING([--enable-memmon],
[use memmon library (Must set MEMMON_PATH)]),
[USEMEMMON=$enableval],[USEMEMMON=no])

AC_ARG_ENABLE(stackmon,
AC_HELP_STRING([--enable-stackmon],
[use stackmonitor library (Must set STACKMON_PATH)]),
[USESTACKMON=$enableval],[USESTACKMON=no])

AC_ARG_ENABLE(exithandler,
AC_HELP_STRING([--enable-exithandler],
[use exithandler library (Must set EXITHANDLER_PATH)]),
[USEEXITHANDLER=$enableval],[USEEXITHANDLER=no])


AC_ARG_ENABLE(filesystem-hints,
AC_HELP_STRING([--enable-filesystem-hints],
[set hints for gpfs or lustre parallel file systems]),
[FSHINTS=$enableval],[FSHINTS=no])

AC_ARG_WITH([piovdc],[--with-piovdc  Use the piovdc interface library (required for compression)], [with_piovdc="$withval"],[with_piovdc=no])  


##############################################################################

#
# declare "precious" environment variables
# (implicit AC_SUBST)
#


AC_ARG_VAR(MPI_INC,[MPI include flags])
AC_ARG_VAR(MPI_LIB,[MPI lib flags])

AC_ARG_VAR(FC,The Fortran compiler)
AC_ARG_VAR(MPIFC,The MPI Fortran compiler)

AC_ARG_VAR(FFLAGS,User-defined Fortran compiler flags)
AC_ARG_VAR(CFLAGS,Customized C source compilation flags)

AC_ARG_VAR(CC,The C compiler)
AC_ARG_VAR(MPICC,The MPI C compiler)

# fortran flags

AC_ARG_VAR(OPT,Fortran compiler flag for optimization level) 
AC_ARG_VAR(REAL8,[Fortran compiler flag for setting the default REAL size to REAL(KIND=8)])
AC_ARG_VAR(BIT64,Fortran compiler flag for generating 64-bit objects)
AC_ARG_VAR(ENDIAN,Fortran compiler flag for converting big-endian to little-endian)
AC_ARG_VAR(MODINC,Fortran compiler flag for specifying module search path)

# end fortran flags

AC_ARG_VAR(MCT_PATH,Path to the MCT root directory)

AC_ARG_VAR(INCLUDEPATH,Additional library and module paths with MODINC)

AC_ARG_VAR(AR,Archive command)
AC_ARG_VAR(ARFLAGS,Archive command flags)
AC_ARG_VAR(RM,rm command)
AC_ARG_VAR(AWK,awk command)

AC_ARG_VAR(PNETCDF_PATH,Path to pnetcdf)
AC_ARG_VAR(NETCDF_PATH,[Path to (serial) netcdf])

AC_ARG_VAR(MEMMON_PATH,Path to memmon library)
AC_ARG_VAR(STACKMON_PATH,Path to stackmonitor library)
AC_ARG_VAR(EXITHANDLER_PATH,Path to exithandler library)

AC_ARG_VAR(SIZEOF_CHAR,[Size of char])
AC_ARG_VAR(SIZEOF_INT,[Size of int])
AC_ARG_VAR(SIZEOF_FLOAT,[Size of float])
AC_ARG_VAR(SIZEOF_DOUBLE,[Size of double])
AC_ARG_VAR(SIZEOF_VOID_P,[Size of void*])


############################################################################


# SET TEMPORARY VARIABLES

# OS AND PLATFORM NAME
test "$osname"=NONE && osname=`uname -s`
test "$machinename"=NONE && machinename=`uname -m`
fullhostname=`hostname -f 2>/dev/null | sed 1q`

# CHECK FOR THE COMPILERS
AC_PROG_CC(xlc cc gcc)
#AX_PROG_CC_MPI()


# CHECK FOR BYTE ORDERING
AC_C_BIGENDIAN

AC_PROG_FC
#AX_PROG_FC_MPI()

AC_LANG_PUSH(Fortran)
AC_FC_SRCEXT(F90)

AC_FC_FREEFORM()
AC_FC_LINE_LENGTH()

# When using a cache file it appears that sometimes ac_ext does not get set
if test -z "$ac_ext"; then
   ac_ext=".F90"
fi

# How to preprocess fortran files
AC_FC_PP_DEFINE()

# CHECK HOW TO GET THE COMPILER VERSION.

#AX_FC_VERSION()

#fc_version="$ac_fc_version_output"



# CHECK MPI BY DEFAULT
if test "$USE_MPISERIAL" = "yes"
  then CPPDEFS="$CPPDEFS -D_MPISERIAL -DNO_MPIMOD "
       echo "building for use with mpi-serial"
  MPICC=$CC

  if test -z "$MPI_INC"; then
    if test -n "$MCT_PATH"
         then 
           echo 'using MCT_PATH from environment'
	   MPI_PATH=$MCT_PATH/mpi-serial
         else MPI_PATH=`pwd`/../mct/mpi-serial    # note: will not work for VPATH build
              echo 'using default mpi-serial from the pio link'
    fi

    MPI_LIB="-L$MPI_PATH -lmpi-serial"
    MPI_INC="-I$MPI_PATH" 
  fi
  USE_MPIIO="no"
  USE_PNETCDF="no"
  MPIFC=$FC
else
  CPPDEFS="$CPPDEFS -DSPMD -DHAVE_MPI"
  ACX_MPI()
  # Check for mpi.mod
  AC_LANG_PUSH(Fortran)
  AC_FC_SRCEXT(F90)
  AC_MSG_CHECKING([for mpi.mod])
  OLDFC=$FC
  FC=$MPIFC
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[       use mpi])],[AC_MSG_RESULT(yes)],
    [CPPDEFS="$CPPDEFS -DNO_MPIMOD"; AC_MSG_RESULT(no)])
  FC=$OLDFC

fi

AC_LANG_POP(Fortran)



###########################################################################
if test "$FC_USE_C_SIZEOF" = "default"
  then AC_LANG_PUSH(Fortran)
       AC_MSG_CHECKING([Fortran support for c_sizeof()])
       AC_COMPILE_IFELSE( [[
        use iso_c_binding, only : c_sizeof 
        implicit none
        integer :: a,b 
        b = c_sizeof(a)
        end]],
        FC_USE_C_SIZEOF="yes", FC_USE_C_SIZEOF="no" )
       AC_LANG_POP(Fortran)

       if test "$FC_USE_C_SIZEOF" = "yes"
         then AC_MSG_RESULT(yes)
         else 
            AC_MSG_RESULT(no)
            AC_MSG_CHECKING([Fortran support for sizeof()])
            AC_COMPILE_IFELSE( [[
              implicit none
              integer :: a,b 
              b = sizeof(a)
              end]],
              FC_USE_SIZEOF="yes", FC_USE_SIZEOF="no" )
         AC_LANG_POP(Fortran)

         if test "$FC_USE_SIZEOF" = "yes"
           then AC_MSG_RESULT(yes)
           else AC_MSG_RESULT(no)
         fi
    fi
fi

if test "$FC_USE_C_SIZEOF" = "no"
  then CPPDEFS=" $CPPDEFS -DNO_C_SIZEOF "
  if test "$FC_USE_SIZEOF" = "no"
    then CPPDEFS=" $CPPDEFS -DNO_SIZEOF "
       echo "disabling Fortran sizeof()"
  fi
fi
#####

# A HACK TO FIX ACX_MPI TO GET MPI_LIB TO BE AN EMPTY STRING
if test "$MPI_LIB" = " "; then
   MPI_LIB=""
fi


if test "$USE_MPIIO" = "default"
  then AC_MSG_CHECKING(MPI-IO support in MPI implementation)
       FFLAGS="$FFLAGS $MPI_INC"
       AC_LANG_PUSH(Fortran)
       AC_COMPILE_IFELSE( [program mpicheck
                           include "mpif.h"
                           integer,parameter:: a = MPI_MODE_RDWR
                           end program],
                USE_MPIIO="yes", USE_MPIIO="no" )
       AC_LANG_POP(Fortran)
       if test "$USE_MPIIO" = "yes"
         then AC_MSG_RESULT(yes)
         else AC_MSG_RESULT(no)
       fi
fi


if test "$USE_MPIIO" = "yes"
  then
       CPPDEFS="$CPPDEFS -DUSEMPIIO"
fi

#####


# DEFINE VARIABLES ACCORDING TO OS AND COMPILER

echo "Full hostname= $fullhostname"
echo "Hostname=$ac_hostname"
echo "Machine=$machinename"
echo "OS=$osname"

# CHECK OS NAME
if echo $osname | grep -i aix >/dev/null 2>&1; then
   SYSDEF="AIX"
fi
if echo $osname | grep -i darwin >/dev/null 2>&1; then
   SYSDEF="DARWIN"
fi
if echo $osname | grep -i unix_system_v >/dev/null 2>&1; then
   SYSDEF="UNIXSYSTEMV"
fi
if echo $osname | grep -i irix >/dev/null 2>&1; then
   SYSDEF="IRIX"
fi
if echo $osname | grep -i irix64 >/dev/null 2>&1; then
   SYSDEF="IRIX64"
fi
if echo $osname | grep -i linux >/dev/null 2>&1; then
   SYSDEF="LINUX"
fi
if echo $osname | grep -i osf1 >/dev/null 2>&1; then
   SYSDEF="OSF1"
fi
if echo $osname | grep -i super >/dev/null 2>&1; then
   SYSDEF="SUPERUX"
fi
if echo $osname | grep -i sun >/dev/null 2>&1; then
   SYSDEF="SUNOS"
fi
if echo $osname | grep -i t3e >/dev/null 2>&1; then
   SYSDEF="T3E"
fi
if echo $osname | grep -i unicos >/dev/null 2>&1; then
   SYSDEF="UNICOS"
fi
if test -z "$SYSDEF"; then
   AC_MSG_WARN([OPERATING SYSTEM UNKNOWN])
   SYSDEF="UNKNOWNOS"
fi

CPPDEFS=" $CPPDEFS -DSYS$SYSDEF "

AC_FC_MODULE_EXTENSION

AC_FC_MODULE_FLAG


# ARCHIVE COMMAND SIMILAR ACROSS ALL PLATFORMS 
if test -z "$AR"; then
   AR="ar"
   ARFLAGS="cq"
fi

if test -z "$RM"; then
   RM="rm -f"
fi

if test -z "$AWK"; then
   AWK="awk"
fi


##########################################################

#
# optional features for PIO
#
# these must be handled in the order they need to be linked
#
# need to handle the -I and the -D separately so can
# construct an xlf fortran version of the -D args if necessary
#

INCLUDEPATH=" $INCLUDEPATH$MPI_INC"
PLIBS=""

if test "$USE_NETCDF4" = "yes"
  then
   USE_NETCDF=yes
fi

if test "$USE_NETCDF" = "yes"
  then  
       if test -n "$NETCDF_PATH"
         then echo "using NETCDF_PATH from environment"
         else NETCDF_PATH=/usr/local/netcdf
              AC_MSG_WARN( [NETCDF_PATH not found in environment, \ 
                           defaulting to $NETCDF_PATH])
       fi             
       AC_LANG_PUSH(Fortran)
       AC_FC_SRCEXT(F90)
       ACX_NETCDF()

       AC_LANG_POP(Fortran)

       CPPDEFS=" $CPPDEFS -D_NETCDF $NETCDF_CPPFLAGS"
       FCFLAGS=" $FCFLAGS $NETCDF_FCFLAGS"
       LDLIBS=" $LDLIBS $NETCDF_LDFLAGS $NETCDF_LIBS "

       if test ! -f $NETCDF_PATH/include/netcdf.mod && \
          test ! -f $NETCDF_PATH/include/NETCDF.mod
         then AC_MSG_ERROR( [netcdf.mod not found in NETCDF_PATH/include \
                            check the environment variable NETCDF_PATH])
       fi

  else
       CPPDEFS="$CPPDEFS -D_NONETCDF"
       USE_NETCDF4="no"
fi

if test "$USE_NETCDF4" = "yes"
  then
   CPPDEFS="$CPPDEFS -D_NETCDF4"
   if test -e "$NETCDF_PATH/lib/pkgconfig/"
     then
     export PKG_CONFIG_PATH="$PKG_CONFIG_PATH:$NETCDF_PATH/lib/pkgconfig/"
     PKG_CHECK_MODULES([NETCDF4], netcdf)
     LDLIBS=" $LDLIBS $NETCDF4_LIBS"
   else
     # Did the user specify a location for the HDF5 library?
     AC_MSG_CHECKING([whether a location for the HDF5 library was specified])
     AC_ARG_WITH([hdf5],
                 [AS_HELP_STRING([--with-hdf5=<directory>],
                                 [Specify location of HDF5 library. HDF5 is required for netCDF-4. Configure will expect to find subdir lib.])],
               [HDF5DIR=$with_hdf5])
     AC_MSG_RESULT([$HDF5DIR])
     AC_SUBST(HDF5DIR, [$HDF5DIR])

     # Did the user specify a location for the ZLIB library?
     AC_MSG_CHECKING([whether a location for the ZLIB library was specified])
     AC_ARG_WITH([zlib],
                 [AS_HELP_STRING([--with-zlib=<directory>],
                                 [Specify location of ZLIB library. ZLIB is required for netCDF-4. Configure will expect to find subdir lib.])],
               [ZLIBDIR=$with_zlib])
     AC_MSG_RESULT([$ZLIBDIR])
     AC_SUBST(ZLIBDIR, [$ZLIBDIR])

     # Did the user specify a location for the SZLIB library?
     AC_MSG_CHECKING([whether a location for the SZLIB library was specified])
     AC_ARG_WITH([szlib],
                 [AS_HELP_STRING([--with-szlib=<directory>],
                                 [Specify location of SZLIB library. SZLIB is not required for netCDF-4, but you must specify its location \
   			      if HDF5 was built with SZLIB. Configure will expect to find subdir lib.])],
               [SZLIBDIR=$with_szlib])
     AC_MSG_RESULT([$SZLIBDIR])
     AC_SUBST(SZLIBDIR, [$SZLIBDIR])
	
     CPPDEFS="$CPPDEFS -D_NETCDF4"
     if test -n "$HDF5DIR"
       then LDLIBS=" $LDLIBS -L$HDF5DIR/lib"
     fi 
     if test -n "$ZLIBDIR"
       then LDLIBS=" $LDLIBS -L$ZLIBDIR/lib"
     fi 
     if test -n "$SZLIBDIR"
       then LDLIBS=" $LDLIBS -L$SZLIBDIR/lib"
     fi 
     LDLIBS=" $LDLIBS -lhdf5 -lz -lsz"   
   fi
fi




if test ! "$USE_PNETCDF" = "no"
  then
       if test -n "$PNETCDF_PATH"
         then echo "using PNETCDF_PATH from environment"
         else PNETCDF_PATH=/usr/local/pnetcdf
              AC_MSG_WARN( [PNETCDF_PATH not found in environment, \
                           defaulting to $PNETCDF_PATH])
       fi
       USE_PNETCDF="yes"
       AC_MSG_CHECKING([for pnetcdf.inc in $PNETCDF_PATH/include])
       if test -f $PNETCDF_PATH/include/pnetcdf.inc
       then AC_MSG_RESULT([yes])
       else
         AC_MSG_WARN( [pnetcdf.inc not found in $PNETCDF_PATH/include \
                             disabling pnetcdf support ])
         USE_PNETCDF=no
       fi

       AC_MSG_CHECKING([for libpnetcdf.a in $PNETCDF_PATH/lib])
       if test -f $PNETCDF_PATH/lib/libpnetcdf.a
         then
	 AC_MSG_RESULT(yes) 
         else
         AC_MSG_WARN( [libpnetcdf.a not found in $PNETCDF_PATH/lib \
                             disabling pnetcdf support ])
         USE_PNETCDF=no
       fi
fi
if test "$USE_PNETCDF" = "yes"
  then
     #  AC_MSG_CHECKING([for pnetcdf])
        CPPDEFS="$CPPDEFS -D_PNETCDF"
        INCLUDEPATH=" $INCLUDEPATH -I$PNETCDF_PATH/include"
        LDLIBS=" $LDLIBS -L$PNETCDF_PATH/lib -lpnetcdf"
  else
       CPPDEFS="$CPPDEFS -D_NOPNETCDF"
       USE_PNETCDF=no
fi

if test "$USE_COMPRESSION" = "yes"
   then
	AC_MSG_CHECKING([for libpiovdc])
	if test "$with_piovdc" != "no"
        then
           LDLIBS=" $LDLIBS -L$with_piovdc/lib -lpiovdc -lpiocommon "
        fi
#        AC_SEARCH_LIBS([createvdf], piovdc, AC_MSG_RESULT([yes]),
#           AC_MSG_ERROR([Could not find libpiovdc required for compression]),[ -lpiocommon ])
#        fi
        AC_FC_LIBRARY_LDFLAGS


        AC_PROG_CXX
	AC_LANG_PUSH(C)
	AM_WITH_EXPAT()
	if test "$with_expat" != "no"
        then
          LDLIBS=" $LDLIBS -L$with_expat/lib -lexpat"
        else
          LDLIBS=" $LDLIBS -lexpat"
        fi	
	AC_LANG_POP
	CPPDEFS="$CPPDEFS -D_COMPRESSION"
        AC_LANG_PUSH(Fortran)
        AC_MSG_CHECKING([how to link f90 with cxx])
        LDFLAGS_TMP=$LDFLAGS
	LDFLAGS="$LDFLAGS -cxxlib"
        AC_LINK_IFELSE([AC_LANG_PROGRAM([],[  integer :: a     ])],[LD=$MPIFC ;AC_MSG_RESULT([ using $MPIFC])],
          [LD=$MPICXX; LDFLAGS=$LDFLAGS_TMP; AC_MSG_RESULT([ using $MPICXX])])
else 
        LD=$MPIFC
fi


if test "$USE_TIMING" = "yes"
  then
    CPPDEFS=" $CPPDEFS -DTIMING"
    INCLUDEPATH=" $INCLUDEPATH -I../timing"
    LDLIBS=" $LDLIBS -L../timing -ltiming"
fi



#####


if test "$USE_MPI2" = "no"
  then CPPDEFS="$CPPDEFS -DNO_MPI2"
       echo "disabling MPI-2 calls"
fi

#####

if test "$USEMEMMON" = "yes"
  then
        CPPDEFS="$CPPDEFS -D_TESTMEM -D_MEMMON"
        if test -n "$MEMMON_PATH"
          then FCFLAGS="$FCFLAGS -qnullterm"
               LDLIBS=" $LDLIBS -L$MEMMON_PATH -lmemmon.rts"
          else AC_MSG_WARN([MEMMON_PATH is not set])
       fi
fi


if test "$USESTACKMON" = "yes"
  then
       CPPDEFS="$CPPDEFS -D_TESTMEM -D_STACKMON" 
       if test -n "$STACKMON_PATH"
         then LDLIBS="$LDIBS -L$STACKMON_PATH -lstackmonitor"
         else AC_MSG_WARN([STACKMON_PATH is not set])
       fi
fi


if test "$USEEXITHANDLER" = "yes"
  then
       if test -n "$EXITHANDLER_PATH"
         then LDLIBS=" $LDLIBS -L$EXITHANDLER_PATH -lexithandler"
         else AC_MSG_WARN([EXITHANDLER_PATH is not set])
       fi
fi

if test "$FSHINTS" != "no"; then
  fshints_set="no"
  if test "$FSHINTS" = "gpfs"; then
    CPPDEFS="$CPPDEFS -DPIO_GPFS_HINTS"
    fshints_set="yes"
  fi
  if test "$FSHINTS" = "lustre"; then
    CPPDEFS="$CPPDEFS -DPIO_LUSTRE_HINTS"
    fshints_set="yes"
  fi
  if test "$fshints_set" = "no"; then
    AC_MSG_WARN([$FSHINTS is not a valid option for --enable-filesystem-hints.  Valid Options are gpfs and lustre.])
  fi 
fi

if test "$DEBUG" = "yes"
  then OPT="-g"
fi


PIO_PATH=`pwd`


#############################################################################

echo
echo Output Variables: \
{CC=$CC} \
{FC=$FC} \
{CXX=$CXX} \
{MPICC=$MPICC} \
{MPIFC=$MPIFC}\
{MPICXX=$MPICXX} \
{FPPFLAGS=$FPPFLAGS} \
{LD=$LD} \
{OPT=$OPT} \
{BIT64=$BIT64} \
{REAL8=$REAL8} \
{ENDIAN=$ENDIAN} \
{FFLAGS=$FFLAGS} \
{CFLAGS=$CFLAGS} \
{FCFLAGS=$FCFLAGS} \
{FC_DEFINE=$FC_DEFINE} \
{MPI_INC=$MPI_INC} \
{MPI_LIB=$MPI_LIB} \
{MODINC=$MODINC}  \
{AR=$AR} \
{ARFLAGS=$ARFLAGS} \
{USE_PNETCDF=$USE_PNETCDF} \
{PNETCDF_PATH=$PNETCDF_PATH} \
{USE_NETCDF=$USE_NETCDF} \
{NETCDF_PATH=$NETCDF_PATH} \
{USE_COMPRESSION=$USE_COMPRESSION} \
{CPPDEFS=$CPPDEFS} \
{INCLUDEPATH=$INCLUDEPATH} \
{LDLIBS=$LDLIBS} \
{DEBUG=$DEBUG}

echo
AC_SUBST(PIO_PATH)
AC_SUBST(USE_COMPRESSION)
AC_SUBST(CPPDEFS)
AC_SUBST(LDLIBS)
AC_SUBST(LD)
AC_OUTPUT

echo Please check the Makefile.conf
echo Have a nice day!

# test -z is true for empty strings
# test -n is true for non-empty strings




