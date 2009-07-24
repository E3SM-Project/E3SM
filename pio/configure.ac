#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.


AC_INIT(PIO, 1.0)


#
# take input from Makefile.conf.in
# output to Makefile.conf.out
#
AC_CONFIG_FILES(Makefile.conf)


#
# configure options
#

# example:
#   AC_ARG_ENABLE(myvar,
#   AC_HELP_STRING([--enable-myvar],[turn on myvar]),
#   [USEMYVAR=$enableval],[USEMYVAR=default])
#   AC_SUBST(USEMYVAR)

AC_ARG_ENABLE(fpp,
AC_HELP_STRING([--enable-fpp],[explicit preprocessing of fortran source]),
[EXPLICIT_FPP=$enableval],[EXPLICIT_FPP=default])
AC_SUBST(EXPLICIT_FPP)

AC_ARG_ENABLE(pnetcdf,
AC_HELP_STRING( [--enable-pnetcdf],
  [use pnetcdf (Uses PNETCDF_PATH) (default auto detect)]),
  [USE_PNETCDF=$enableval],[USE_PNETCDF=default])

AC_ARG_ENABLE(netcdf,
AC_HELP_STRING( [--enable-netcdf],
[use serial netcdf (Uses NETCDF_PATH) (default yes)]),
[USE_NETCDF=$enableval],[USE_NETCDF=yes])

AC_ARG_ENABLE(mpiio,
AC_HELP_STRING([--enable-mpiio],
[Enable support for MPI-IO (default auto detect)]),
[USE_MPIIO=$enableval],[USE_MPIIO=default])

AC_ARG_ENABLE(mct,
AC_HELP_STRING([--enable-mct],[use mct (Uses MCT_PATH) (default yes)]),
[USE_MCT=$enableval],[USE_MCT=yes])

AC_ARG_ENABLE(box,
AC_HELP_STRING([--enable-box],[use box rearranger (default yes)]),
[USEBOX=$enableval],[USEBOX=yes])

AC_ARG_ENABLE(debug,
AC_HELP_STRING([--enable-debug],
[Turn off optimization and turn on -g (same as OPT=-g)]),
[DEBUG=$enableval],[DEBUG=no])
AC_SUBST(DEBUG)

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

AC_ARG_ENABLE(sizeof,
AC_HELP_STRING([--enable-sizeof],
[Enable use of sizeof() in Fortran (default auto detect)]),
[USE_SIZEOF=$enableval],[USE_SIZEOF=default])

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


##############################################################################

#
# declare "precious" environment variables
# (implicit AC_SUBST)
#


AC_ARG_VAR(MPI_INC,[MPI include flags])
AC_ARG_VAR(MPI_LIB,[MPI lib flags])

AC_ARG_VAR(FPP,C-preprocessor for Fortran source code)
AC_ARG_VAR(FPPFLAGS,C-preprocessing flags for Fortran source code)

AC_ARG_VAR(MACHDEFS,MCT compiler and OS flags)

AC_ARG_VAR(FC,The Fortran compiler)
AC_ARG_VAR(MPIF90,The MPI Fortran compiler)

AC_ARG_VAR(FFLAGS,User-defined Fortran compiler flags)
AC_ARG_VAR(ALLCFLAGS,Customized C source compilation flags)

AC_ARG_VAR(CC,The C compiler)
AC_ARG_VAR(MPICC,The MPI C compiler)

# fortran flags

AC_ARG_VAR(OPT,Fortran compiler flag for optimization level) 

AC_ARG_VAR(REAL8,[Fortran compiler flag for setting the default REAL size to REAL(KIND=8)])
AC_ARG_VAR(BIT64,Fortran compiler flag for generating 64-bit objects)
AC_ARG_VAR(ENDIAN,Fortran compiler flag for converting big-endian to little-endian)

# end fortran flags

AC_ARG_VAR(MCT_PATH,Path to the MCT root directory)

AC_ARG_VAR(MODINC,Fortran compiler flag for specifying module search path)
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

# HARDCODE SPECIFIC MACHINES FOR EXTRAORDINARY CIRCUMSTANCES

# CHECK IF WE ARE ON THE EARTH SIMULATOR
ES="NO"
if echo $osname | grep -i esos >/dev/null 2>&1; then
   ES="YES"
fi
if echo $osname | grep -i hp-ux >/dev/null 2>&1; then   
   if test "$ac_hostname" = "moon"; then 
      ES="YES"
      # TELLS CONFIGURE NOT TO RUN ANY TESTS THAT REQUIRE EXECUTION
      cross_compiling="yes"
   fi
fi
if test "$ES" = "YES"; then   
   echo "Using preset configuration values for the Earth Simulator"
   if test -z "$CC"; then
      CC="escc"
   fi
   if test -z "$F90"; then
      F90="esf90"
   fi
   if test -z "$MPIF90"; then
      MPIF90="esmpif90"
   fi
   if test -z "$AR"; then
      AR="esar cqs"
   fi
   if test -z "FPP"; then
      FPPFLAGS=" "
   fi
   if test -z "$F90FLAGS"; then
      F90FLAGS="-EP -Wf'-pvctl fullmsg -L fmtlist transform map'"
   fi
   if test -z "$OPT"; then
      OPT="-C vopt"
   fi      
   if test -z "$MACHDEFS"; then
      MACHDEFS="-DESVEC"
   fi
fi

# CHECK IF WE ARE ON THE NCSA IA-64 CLUSTER TITAN
if test "$ac_hostname" = "user02"; then
   if echo $osname | grep -i linux >/dev/null 2>&1; then
      if test -z "$MPI_LIB"; then
         if test -z "$MPI_INC"; then
            echo "Using preset MPI configuration for titan"
	    MPI_LIB="-L/usr/local/vmi/mpich/lib/intel -lmpich -lvmi -ldl -lpthread"
	    MPI_INC="-I/usr/local/vmi/mpich/include"
         fi
      fi
   fi
fi



# START TESTS

# CHECK FOR THE C COMPILER

AC_PROG_CC(xlc_r cc)

# CHECK MPI BY DEFAULT
if test "$USE_MPISERIAL" = "yes"
  then PDEFS="$PDEFS -D_MPISERIAL"
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
else
  ACX_MPI()
fi

#
# determine sizes
#
if test -n "$SIZEOF_CHAR"
then
  echo "Using SIZEOF_CHAR from environment = $SIZEOF_CHAR"
else
  if test -z "$host_alias"
    then
      AC_CHECK_SIZEOF([char])
      SIZEOF_CHAR=$ac_cv_sizeof_char
    else
      SIZEOF_CHAR=1
    fi
fi

if test -n "$SIZEOF_INT"
then
  echo "Using SIZEOF_INT from environment = $SIZEOF_INT"
else
  if test -z "$host_alias"
    then
      AC_CHECK_SIZEOF([int])
      SIZEOF_INT=$ac_cv_sizeof_int
    else
      SIZEOF_INT=4
    fi
fi

if test -n "$SIZEOF_FLOAT"
then
  echo "Using SIZEOF_FLOAT from environment = $SIZEOF_FLOAT"
else
  if test -z "$host_alias"
    then
      AC_CHECK_SIZEOF([float])
      SIZEOF_FLOAT=$ac_cv_sizeof_float
    else
      SIZEOF_FLOAT=4
    fi
fi

if test -n "$SIZEOF_DOUBLE"
then
  echo "Using SIZEOF_DOUBLE from environment = $SIZEOF_DOUBLE"
else
  if test -z "$host_alias"
    then
      AC_CHECK_SIZEOF([double])
      SIZEOF_DOUBLE=$ac_cv_sizeof_double
    else
      SIZEOF_DOUBLE=8
    fi
fi

if test -n "$SIZEOF_VOID_P"
then
  echo "Using SIZEOF_VOID_P from environment = $SIZEOF_VOID_P"
else
  if test -z "$host_alias"
    then
      AC_CHECK_SIZEOF([void *])
      SIZEOF_VOID_P=$ac_cv_sizeof_void_p
    else
      SIZEOF_VOID_P=8
    fi
fi

###########################################################################

AC_SUBST(PIO_PATH)

AC_SUBST(PDEFS)
AC_SUBST(PDEFSF)
AC_SUBST(PINCS)
AC_SUBST(LDLIBS)

# CHECK FOR THE FORTRAN COMPILER
AC_PROG_F90

AC_LANG_PUSH(Fortran 90)

# CHECK MPI BY DEFAULT
if test "$USE_MPISERIAL" = "yes"; then
  USE_PNETCDF="no"
  MPIF90=$F90
else
  ACX_MPI()
fi

# CHECK FOR MPI LIBRARIES

OLDF90FLAGS="$F90FLAGS"

if test -n "MPI_INC"; then
   F90FLAGS="$F90FLAGS $MPI_INC"
fi



F90FLAGS="$OLDF90FLAGS"

# A HACK TO FIX ACX_MPI TO GET MPI_LIB TO BE AN EMPTY STRING
if test "$MPI_LIB" = " "; then
   MPI_LIB=""
fi

# FOR SANITY, CHECK THAT FILENAME EXTENSION FOR FC IS CONSISTENT WITH F90
OLDF90="$F90"
F90="$MPIF90"

AC_COMPILE_IFELSE(
  [subroutine oof(); return; end], [],
  [AC_MSG_WARN([$FC FAILED TO COMPILE FILENAME EXTENSION $ac_ext])
  ])
#PDEFS=""
if test "$USE_MPIIO" = "default"
  then AC_MSG_CHECKING(MPI-IO support in MPI implementation)
       FFLAGS="$FFLAGS $MPI_INC"
       AC_LANG_PUSH(Fortran 90)
       AC_COMPILE_IFELSE( [program mpicheck
                           include "mpif.h"
                           integer,parameter:: a = MPI_MODE_RDWR
                           end program],
                USE_MPIIO="yes", USE_MPIIO="no" )
       AC_LANG_POP(Fortran 90)
       if test "$USE_MPIIO" = "yes"
         then AC_MSG_RESULT(yes)
         else AC_MSG_RESULT(no)
       fi
fi


if test "$USE_MPIIO" = "yes"
  then
       PDEFS="$PDEFS -DUSEMPIIO"
fi

#####

F90="$OLDF90"


# CHECK HOW TO GET THE COMPILER VERSION.

_AC_PROG_F90_VERSION

f90_version="$ac_f90_version_output"


AC_LANG_POP(Fortran 90)

# CHECK FOR BYTE ORDERING
AC_C_BIGENDIAN


# CHECK COMPILER NAME
if echo $F90 | grep xlf >/dev/null 2>&1; then
   echo "Fortran Compiler is XLF"
   CPRDEF="XLF"
#   if test -z "$REAL8"; then
#      REAL8="-qrealsize=8"
#   fi
   if test -z "$OPT"; then
      OPT="-O4 "
   fi
   if test -z "$F90FLAGS"; then
      if test "$ac_cv_f90_fpp_compiler" = "yes"; then
         F90FLAGS="-qarch=auto -qsuffix=f=F90:cpp=F90"
      fi
      if test "$ac_cv_f90_fpp_compiler" = "no"; then
         if test "$ac_cv_f90_fpp_external" = "yes"; then
	    F90FLAGS="-qarch=auto -qsuffix=f=f90"
         fi
      fi
   fi
fi

if echo $F90 | grep pgf90 >/dev/null 2>&1; then
   echo "Fortran Compiler is Portland Group"
   CPRDEF="PGI"
   if test -z "$REAL8"; then
      REAL8="-r8"
   fi
   if test -z "$BIT64"; then
      BIT64="-pc 64"
   fi
   if test -z "$ENDIAN"; then
      ENDIAN="-byteswapio"
   fi
   if test -z "$OPT"; then
      OPT="-O2"
   fi
fi

if echo $F90 | grep ftn >/dev/null 2>&1; then
 if echo $f90_version | grep -i Portland >/dev/null 2>&1; then
   echo "Fortran Compiler is Portland Group, XT3"
   CPRDEF="PGI"
   SYSDEF="CATAMOUNT"
   if test -z "$REAL8"; then
      REAL8="-r8"
   fi
   if test -z "$BIT64"; then
      BIT64="-pc 64"
   fi
   if test -z "$ENDIAN"; then
      ENDIAN="-byteswapio"
   fi
   if test -z "$OPT"; then
      OPT="-O2"
   fi
 else
   echo "Fortran Compiler is Cray"
   CPRDEF="CRAY" 
#   if test -z "$F90FLAGS"; then 
#      F90FLAGS="-c -em -dy -rm -Omsgs,negmsgs"
#   fi
   if test -z "$MODINC"; then
      MODINC="-p"
   fi
#   if test -z "$OPT"; then
#      OPT="-Oscalar3 -Ovector1 -Ostream1"
#   fi
 fi
fi

if echo $F90 | grep frt >/dev/null 2>&1; then
   echo "Fortran Compiler is UXP/V"
   echo "Suggested additional vectorization flags: -Wv,-s5,-t3,-noalias,-ilfunc,-md"
   CPRDEF="FUJITSU"
   if test -z "$F90FLAGS"; then
      F90FLAGS="-Am -X9"
   fi
   if test -z "$BIT64"; then
      BIT64="-KA64"
   fi
   if test -z "$REAL8"; then
      REAL8="-Ad"
   fi
fi

if echo $F90 | grep epcf90 >/dev/null 2>&1; then
   echo "Everest has never tested the epcf90 compiler :("
   CPRDEF="EPC"
fi

if echo $F90 | grep lf95 >/dev/null 2>&1; then
   echo "Fortran Compiler is Lahey"
   CPRDEF="LAHEY"
fi

if echo $F90 | grep ifc >/dev/null 2>&1; then
   echo "Fortran Compiler is Intel 7.x or earlier"
   echo "Intel ifc compiler must set the environment variable F_UFMTENDIAN=big to do endian conversion"
   CPRDEF="INTEL"
   if test -z "$REAL8"; then
      REAL8="-r8"
   fi
   if test -z "$F90FLAGS"; then
      F90FLAGS="-w"
   fi
   if test -z "$OPT"; then
      OPT="-O2"
   fi
fi

if echo $F90 | grep efc >/dev/null 2>&1; then
   echo "Fortran Compiler is Intel 7.x or earlier for IA-64"
   echo "Intel efc compiler must set the environment variable F_UFMTENDIAN=big to do endian conversion"
   CPRDEF="INTEL"
   if test -z "$REAL8"; then
      REAL8="-r8"
   fi
   if test -z "$F90FLAGS"; then
      F90FLAGS="-w -ftz"
   fi
   if test -z "$OPT"; then
      OPT="-O2"
   fi
fi

if echo $F90 | grep ifort >/dev/null 2>&1; then
   echo "Fortran Compiler is Intel ifort"
   CPRDEF="INTEL"
   if test -z "$REAL8"; then
      REAL8="-r8"
   fi
   if test -z "$F90FLAGS"; then
      F90FLAGS="-w -ftz -assume byterecl"
   fi
   if test -z "$ENDIAN"; then
      ENDIAN="-convert big_endian"
   fi
   if test -z "$OPT"; then
      OPT="-O2"
   fi
fi

if echo $F90 | grep pathf90 >/dev/null 2>&1; then
   echo "Fortran Compiler is PathScale"
   CPRDEF="PATHSC"
   if test -z "$REAL8"; then
      REAL8="-r8"
   fi
   if test -z "$BIT64"; then
      BIT64="-m64"
   fi
   if test -z "$OPT"; then
      OPT="-O2"
   fi
fi

if echo $F90 | grep g95 >/dev/null 2>&1; then
   echo "Fortran Compiler is GNU"
   CPRDEF="GNU"
fi

if echo $f90_version | grep -i nag >/dev/null 2>&1; then
   echo "Fortran Compiler is NAG"
   CPRDEF="NAG"
   if test -z "$F90FLAGS"; then
      F90FLAGS="-dusty -kind=byte"
   fi
   if test -z "$OPT"; then
      OPT="-O2"
   fi
fi

if echo $f90_version | grep -i absoft >/dev/null 2>&1; then
   echo "Fortran Compiler is Absoft"
   CPRDEF="ABSOFT"
   if test -z "$REAL8"; then
      REAL8="-N113"
   fi    
   if test -z "$MODINC"; then
      MODINC="-p"
   fi
   if test -z "$OPT"; then
      OPT="-O2"
   fi
fi

if echo $f90_version | grep -i workshop >/dev/null 2>&1; then
   echo "Fortran Compiler is Workshop"
   CPRDEF="WORKSHOP"
   if test -z "$MODINC"; then
      MODINC="-M"
   fi
fi

if echo $f90_version | grep -i mipspro >/dev/null 2>&1; then
   echo "Fortran Compiler is MIPSPro"
   CPRDEF="MIPSPRO"
   EXTRACFLAGS="-64"
   if test -z "$OPT"; then
      OPT="-O3"
   fi
   if test -z "$REAL8"; then
      REAL8="-r8"
   fi
   if test -z "$BIT64"; then
      BIT64="-64"
   fi
fi

if echo $f90_version | grep -i compaq >/dev/null 2>&1; then
   echo "Fortran Compiler is Compaq"
   CPRDEF="COMPAQ"
   MPI_LIB="$MPI_LIB -lelan"
   if test -z "$OPT"; then
      OPT="-fast"
   fi
   if test -z "$REAL8"; then
      REAL8="-real_size 64"
   fi
   if test -z "$ENDIAN"; then
      ENDIAN="-convert big_endian"
   fi
fi

# Compaq Fortran changed its name to HP Fortran.
# Lets support both versions for now.
if echo $f90_version | grep HP >/dev/null 2>&1; then
   echo "Fortran Compiler is HP"
   CPRDEF="COMPAQ"
   MPI_LIB="$MPI_LIB -lelan"
   if test -z "$OPT"; then
      OPT="-fast"
   fi
   if test -z "$REAL8"; then
      REAL8="-real_size 64"
   fi
   if test -z "$ENDIAN"; then
      ENDIAN="-convert big_endian"
   fi
fi

if echo $f90_version | grep -i sx >/dev/null 2>&1; then
   echo "Fortran Compiler is SX"
   CPRDEF="SX"
   if test -z "$F90FLAGS"; then
      F90FLAGS="-EP -Wf'-pvctl noassoc'"
   fi
   if test -z "$OPT"; then
      OPT="-Chopt"
   fi
fi

if test -z "$CPRDEF"; then
   AC_MSG_WARN([UNKNOWN FORTRAN 90 COMPILER])
   CPRDEF="UNKNOWNCPR"
fi


# CHECK IF FORTRAN COMPILER CREATES MODULE FILES IN UPPERCASE OR LOWERCASE
OLDF90="$F90"
OLDF90FLAGS="$F90FLAGS"

F90="$MPIF90" # use MPIF90, since differences may exist between MPIF90 and F90 
if test "$F90" = "frt"; then # -Am flag needed for frt to gerate mod files
   F90FLAGS="$F90FLAGS -Am"
fi
#if test "$F90" = "ftn"; then # -em flag needed for ftn to gerate mod files
#   F90FLAGS="$F90FLAGS -em"
#fi

AC_F90_MOD_UPPERCASE

# CHECK HOW TO NAME MANGLE C FUNCTIONS SO THAT IT CAN BE CALLED FROM FORTRAN

AC_F90_C_NAME_MANGLING

case $ac_cv_f90_mangling in
  "lower case, underscore") 
	 name_mangling="-DFORTRAN_UNDERSCORE_ -DFORTRANUNDERSCORE";;
  "lower case, no underscore")
	 name_mangling="-DFORTRAN_SAME";;
  "upper case, no underscore")  
	 name_mangling="-DFORTRAN_CAPS_ -DFORTRANCAPS";;
  "lower case, double underscore") 
	 name_mangling="-DFORTRAN_GNUF2C";;
  *) 
         name_mangling="-DFORTRAN_MANGLING_UNKNOWN" 
	 AC_MSG_WARN([UNKNOWN FORTRAN 90 NAME-MANGLING SCHEME]);;
esac

F90="$OLDF90"

# FPP AND FPPFLAGS ARE HARDCODED AS A LAST RESORT:
# DO NOT USE THE GNU C PREPROCESESSOR- 
# IT DOES NOT PROCESS FORTRAN SOURCES CORRECTLY.
# LOOK FOR MACHINE DEPENDENT CPP.
# THE HARCODED FPPFLAGS WORKS WITH MOST CPP COMMANDS.
if test -z "$FPP"; then
   AC_CHECK_PROGS(FPP,[cpp],$CPP)
fi
if test -z "$FPPFLAGS"; then
   FPPFLAGS="-P -C -N -traditional"
fi

# CHECK THAT THE FORTRAN COMPILER CAN CORRECTLY PROCESS THESE DIRECTIVES
# IF NOT, USE THE EXTERNAL C PREPROCESSOR
OLDF90="$F90"
F90="$MPIF90" # use MPIF90, since differences may exist between MPIF90 and F90

defineflag="-Daardvark"
if test "$CPRDEF" = "XLF"; then
   defineflag="-WF,-Daardvark"
fi   
if test "$OLDF90" = "frt"; then
   defineflag="-Wp,-Daardvark"
fi

AC_CHECK_FPP_COMPILER(
[
  implicit none
#if zebra || aardvark
  character(len=25) :: ch  ! '
#endif
  ch="Does this test &
!Comment-Line
          &work? "// &
          "YES!"
], [$defineflag])

if test -n "$FPP"; then
  defineflag="-Daardvark"
AC_CHECK_FPP_EXTERNAL(
[
  implicit none
#if zebra || aardvark
  character(len=25) :: ch  ! '
#endif
  ch="Does this test &
!Comment-Line
          &work? "// &
          "YES!"
], [$defineflag])

fi

F90="$OLDF90"


if test "$EXPLICIT_FPP" = "default"
then
  if test "$ac_cv_f90_fpp_compiler" = "yes" 
  then
     EXPLICIT_FPP="no"
  else
    if test "$ac_cv_f90_fpp_external" = "no"
    then
      EXPLICIT_FPP="no"
      AC_MSG_WARN([FORTRAN SOURCE PREPROCESSING HAS FAILED])
      if test -z "$FPP"; then
        AC_MSG_WARN([TRY DEFINING CPP PATH IN CONFIGURE INVOCATION])
      fi
    else
      EXPLICIT_FPP="yes"
    fi
  fi
fi


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


# SET THE PREPROCESSOR DEFINE FLAGS
if test -z "$MACHDEFS"; then
   MACHDEFS="-DSYS$SYSDEF -D$SYSDEF -DCPR$CPRDEF" # DEFAULT FLAGS
   if test "$CPRDEF" = "XLF"; then # SPECIAL XLF FLAGS
     if test "$ac_cv_f90_fpp_compiler" = "yes"; then
       MACHDEFS="-WF,-DSYS$SYSDEF,-D$SYSDEF,-DCPR$CPRDEF"
     fi
   fi
   if test "$F90" = "frt"; then # SPECIAL FRT FLAGS
      if test "$ac_cv_f90_fpp_compiler" = "yes"; then
         MACHDEFS="-Wp,-DSYS$SYSDEF,-D$SYSDEF,-DCPR$CPRDEF"
      fi
   fi
fi

if test -z "$ALLCFLAGS"; then # SET C COMPILE FLAGS
   ALLCFLAGS="$name_mangling -DSYS$SYSDEF -D$SYSDEF -DCPR$CPRDEF"
   if test -n "$EXTRACFLAGS"; then
      ALLCFLAGS="$name_mangling -O $EXTRACFLAGS"
   fi
fi

# SET HARDCODED VARIABLES AS A LAST RESORT


# INCLUDE FLAG IF NOT ALREADY SET IS MOST LIKELY -I
if test -z "$MODINC"; then
   MODINC="-I"
fi

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

PINCS="$MPI_INC"
PLIBS=""

if test "$USE_NETCDF" = "yes"
  then
       if test -n "$NETCDF_PATH"
         then echo "using NETCDF_PATH from environment"
         else NETCDF_PATH=/usr/local/netcdf
              AC_MSG_WARN( [NETCDF_PATH not found in environment, \ 
                           defaulting to $NETCDF_PATH])
       fi             

       PDEFS="$PDEFS -D_NETCDF"
       PINCS="$PINCS -I$NETCDF_PATH/include"
       PLIBS="$PLIBS -L$NETCDF_PATH/lib -lnetcdf"

       if test ! -f $NETCDF_PATH/include/netcdf.mod && \
          test ! -f $NETCDF_PATH/include/NETCDF.mod
         then AC_MSG_ERROR( [netcdf.mod not found in NETCDF_PATH/include \
                            check the environment variable NETCDF_PATH])
       fi

       if test ! -f $NETCDF_PATH/lib/libnetcdf.a
         then AC_MSG_ERROR( [libnetcdf.a not found in NETCDF_PATH/lib \
                            check the environment variable NETCDF_PATH ])
       fi
  else
       PDEFS="$PDEFS -D_NONETCDF"
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
       if test ! -f $PNETCDF_PATH/include/pnetcdf.inc
         then AC_MSG_WARN( [pnetcdf.inc not found in PNETCDF_PATH/include \
                             disabling pnetcdf support ])
         USE_PNETCDF=no
       fi

       if test ! -f $PNETCDF_PATH/lib/libpnetcdf.a
         then AC_MSG_WARN( [libpnetcdf.a not found in PNETCDF_PATH/lib \
                             disabling pnetcdf support ])
         USE_PNETCDF=no
       fi
fi
if test "$USE_PNETCDF" = "yes"
  then
        PDEFS="$PDEFS -D_PNETCDF"
        PINCS="$PINCS -I$PNETCDF_PATH/include"
        PLIBS="$PLIBS -L$PNETCDF_PATH/lib -lpnetcdf"
  else
       PDEFS="$PDEFS -D_NOPNETCDF"
       USE_PNETCDF=no
fi


if test "$USE_TIMING" = "yes"
  then
    PDEFS="$PDEFS -DTIMING"
    PINCS="$PINCS -I../timing"
    PLIBS="$PLIBS -L../timing -ltiming"
fi



#####

if test "$USE_MCT" = "yes"
  then
       PDEFS="$PDEFS -D_USEMCT"
       if test -n "$MCT_PATH"
         then echo 'using MCT_PATH from environment'
         else MCT_PATH=`pwd`/../mct    # note: will not work for VPATH build
              echo 'using default mct from the pio link'
       fi
       PDEFS="$PDEFS -D_USEMCT"
       PINCS="$PINCS -I$MCT_PATH/mct -I$MCT_PATH/mpeu"
       PLIBS="$PLIBS -L$MCT_PATH/mct -lmct -L$MCT_PATH/mpeu -lmpeu"

       if test ! -f $MCT_PATH/mct/libmct.a
         then AC_MSG_ERROR( [libmct.a not found in MCT_PATH/mct \
                             check the environment variable MCT_PATH ])
       fi


  else
       PDEFS="$PDEFS -D_NOUSEMCT"
fi


if test "$USEBOX" = "yes"
  then PDEFS="$PDEFS -D_USEBOX"
  else PDEFS="$PDEFS -D_NOUSEBOX"
fi


if test "$USE_MPI2" = "no"
  then PDEFS="$PDEFS -DNO_MPI2"
       echo "disabling MPI-2 calls"
fi

#####

if test "$USE_SIZEOF" = "default"
  then AC_LANG_PUSH(Fortran 90)
       AC_MSG_CHECKING([Fortran 90 support for sizeof()])
       AC_TRY_COMPILE( [], [integer :: a,b = sizeof(a)],
                USE_SIZEOF="yes", USE_SIZEOF="no" )
       AC_LANG_POP(Fortran 90)

       if test "$USE_SIZEOF" = "yes"
         then AC_MSG_RESULT(yes)
         else AC_MSG_RESULT(no)
       fi
fi

if test "$USE_SIZEOF" = "no"
  then PDEFS="$PDEFS -DNO_SIZEOF"
       echo "disabling Fortran sizeof()"
fi

#####

if test "$USEMEMMON" = "yes"
  then
        PDEFS="$PDEFS -D_MEMMON"
        if test -n "$MEMMON_PATH"
          then F90FLAGS="$F90FLAGS -qnullterm"
               PLIBS="$PLIBS -L$MEMMON_PATH -lmemmon.rts"
          else AC_MSG_WARN([MEMMON_PATH is not set])
       fi
fi


if test "$USESTACKMON" = "yes"
  then
       PDEFS="$PDEFS -D_STACKMON" 
       if test -n "$STACKMON_PATH"
         then PLIBS="$DLIBS -L$STACKMON_PATH -lstackmonitor"
         else AC_MSG_WARN([STACKMON_PATH is not set])
       fi
fi


if test "$USEEXITHANDLER" = "yes"
  then
       if test -n "$EXITHANDLER_PATH"
         then PLIBS="$PLIBS -L$EXITHANDLER_PATH -lexithandler"
         else AC_MSG_WARN([EXITHANDLER_PATH is not set])
       fi
fi

if test "$FSHINTS" != "no"; then
  fshints_set="no"
  if test "$FSHINTS" = "gpfs"; then
    PDEFS="$PDEFS -DPIO_GPFS_HINTS"
    fshints_set="yes"
  fi
  if test "$FSHINTS" = "lustre"; then
    PDEFS="$PDEFS -DPIO_LUSTRE_HINTS"
    fshints_set="yes"
  fi
  if test "$fshints_set" = "no"; then
    AC_MSG_WARN([$FSHINTS is not a valid option for --enable-filesystem-hints.  Valid Options are gpfs and lustre.])
  fi 
fi


if test "$EXPLICIT_FPP" = "yes"; then
   PDEFSF=""
else
  if test "$CPRDEF" = "XLF"
    then PDEFSF="-WF,"`echo $PDEFS | sed 's/[ ]/,/g'`
    else PDEFSF="$PDEFS"
  fi
fi

LDLIBS="$LDLIBS $PLIBS $MPI_LIB"



if test "$DEBUG" = "yes"
  then OPT="-g"
fi


PIO_PATH=`pwd`


#############################################################################

echo
echo Temp variables: {F90=$F90} {MPIF90=$MPIF90}

echo
echo Output Variables: {CC=$CC} \
{MPICC=$MPICC} \
{ALLCFLAGS=$ALLCFLAGS} \
{EXPLICIT_FPP=$EXPLICIT_FPP} {FPP=$FPP} {FPPFLAGS=$FPPFLAGS} \
{FC=$F90} \
{MPIF90=$MPIF90}\
{OPT=$OPT} \
{BIT64=$BIT64} \
{REAL8=$REAL8} \
{ENDIAN=$ENDIAN} \
{FFLAGS=$FFLAGS} \
{F90FLAGS=$F90FLAGS} \
{MACHDEFS=$MACHDEFS} \
{MPI_INC=$MPI_INC} {MPI_LIB=$MPI_LIB} \
{MODINC=$MODINC} {INCLUDEPATH=$INCLUDEPATH} \
{AR=$AR} \
{ARFLAGS=$ARFLAGS} \
{USE_PNETCDF=$USE_PNETCDF} \
{PNETCDF_PATH=$PNETCDF_PATH} \
{USE_NETCDF=$USE_NETCDF} \
{NETCDF_PATH=$NETCDF_PATH} \
{USE_MCT=$USE_MCT} \
{USEBOX=$USEBOX} \
{PDEFS=$PDEFS} \
{PDEFSF=$PDEFSF} \
{PINCS=$PINCS} \
{DEBUG=$DEBUG}

echo

AC_OUTPUT

echo Please check the Makefile.conf
echo Have a nice day!

# test -z is true for empty strings
# test -n is true for non-empty strings




