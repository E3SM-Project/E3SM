c_ext#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT(MCT, 2.0)

# PROCESS THE FOLLOWING MAKEFILES
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_FILES(Makefile.conf)
AC_CONFIG_HEADER(config.h)

# DECLARE PACKAGE OPTIONS

AC_ARG_ENABLE(mpiserial,
AC_HELP_STRING([--enable-mpiserial],
[Use the included MPI replacement library for single processor]),
[DONOTCHECKMPI="DONOTCHECKMPI"]
)

AC_ARG_ENABLE(selectedrealkind,
AC_HELP_STRING([--enable-selectedrealkind],
[define single precision and double precision numbers using the selected_real_kind function. Default uses the kind inquiry function.]),
[SRKDEF="SELECTEDREALKIND"]
)

AC_ARG_ENABLE(sequence,
AC_HELP_STRING([--enable-sequence],[Modify MCT types to make them contiguous in memory.]),
[SRKDEF="SEQUENCE"],)

AC_ARG_ENABLE(babel,
AC_HELP_STRING([--enable-babel],[Supply this option if you plan on building the Babel bindings to MCT]),
[SRKDEF="SEQUENCE"],)



# DECLARE THE FOLLOWING PRECIOUS VARIABLES

AC_ARG_VAR(MPILIBS,[MPI library command line invocation])
AC_ARG_VAR(MPIHEADER,[MPI header include path with INCLUDEFLAG])
AC_ARG_VAR(FPP,C-preprocessor for Fortran source code)
AC_ARG_VAR(FPPFLAGS,C-preprocessing flags for Fortran source code)
AC_ARG_VAR(MACHDEFS,MCT compiler and OS flags)
AC_ARG_VAR(FC,The Fortran compiler)
AC_ARG_VAR(FCFLAGS,User-defined Fortran compiler flags)
AC_ARG_VAR(CFLAGS,Customized C source compilation flags)
AC_ARG_VAR(OPT,Fortran compiler flag for optimization level) 
AC_ARG_VAR(REAL8,[Fortran compiler flag for setting the default REAL size to REAL(KIND=8)])
AC_ARG_VAR(BIT64,Fortran compiler flag for generating 64-bit objects)
AC_ARG_VAR(ENDIAN,Fortran compiler flag for converting big-endian to little-endian)
AC_ARG_VAR(INCLUDEFLAG,Fortran compiler flag for specifying module search path)
AC_ARG_VAR(INCLUDEPATH,Additional library and module paths with INCLUDEFLAG)
AC_ARG_VAR(AR,Archive command)
AC_ARG_VAR(BABELROOT,Root directory of your Babel installation. i.e.: $BABELROOT/bin/babel $BABELROOT/lib/libsidl.so)
AC_ARG_VAR(COMPILER_ROOT,Root directory of your FORTRAN compiler)
AC_ARG_VAR(FORT_SIZE, Number of bits in Fortran real and double kind)

# INCLUDE BABELROOT and COMPILER_ROOT in Makefile.conf(autoconf output)
AC_SUBST(BABELROOT)
AC_SUBST(COMPILER_ROOT)
AC_SUBST(PYTHON)
AC_SUBST(PYTHONOPTS)

# SET TEMPORARY VARIABLES

# OS AND PLATFORM NAME
test "$osname"=NONE && osname=`uname -s`
test "$machinename"=NONE && machinename=`uname -m`
fullhostname=`hostname -f`


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
      if test -z "$MPILIBS"; then
         if test -z "$MPIHEADER"; then
            echo "Using preset MPI configuration for titan"
	    MPILIBS="-L/usr/local/vmi/mpich/lib/intel -lmpich -lvmi -ldl -lpthread"
	    MPIHEADER="-I/usr/local/vmi/mpich/include"
         fi
      fi
   fi
fi


# Check if we are on the ANL BG/P

if echo $fullhostname | egrep -q '.\.(surveyor|intrepid)\.alcf\.anl\.gov' 
  then if test -z "$FC"; then
         FC=mpixlf90_r
       fi
       if test -z "$CC"; then
         CC=mpixlc_r
       fi
fi



# START TESTS

# CHECK FOR THE C COMPILER
AC_PROG_CC([cc])

# CHECK FOR BYTE ORDERING
AC_C_BIGENDIAN

# CHECK FOR THE FORTRAN COMPILER
AC_PROG_FC

# CHECK FOR MPI LIBRARIES
AC_LANG_PUSH(Fortran)

AC_FC_SRCEXT(F90)

OLDFCFLAGS="$FCFLAGS"

if test -n "$MPIHEADER"; then
   FCFLAGS="$FCFLAGS $MPIHEADER"
fi

# CHECK MPI BY DEFAULT
if test -z "$DONOTCHECKMPI"; then
   ACX_MPI
fi

# DONT CHECK MPI IF SERIALMPI OPTION IS ENABLED
if test -n "$DONOTCHECKMPI"; then
   echo "MPISERIAL ENABLED: BYPASSING MPI CHECK"
   if test -z "$MPIFC"; then
      MPIFC=$FC
   fi
   if test -z "$FORT_SIZE"; then
      FORT_SIZE="real4double8"
      echo "FORT_SIZE IS PRESET TO $FORT_SIZE"
   fi   
   abs_top_builddir=`pwd`
   MPISERPATH=$abs_top_builddir/mpi-serial
   AC_SUBST(MPISERPATH)
   MPIHEADER=-I$MPISERPATH
   MPILIBS="-L$MPISERPATH -lmpi-serial"
fi

FCFLAGS="$OLDFCFLAGS"

# A HACK TO FIX ACX_MPI TO GET MPILIBS TO BE AN EMPTY STRING
if test "$MPILIBS" = " "; then
   MPILIBS=""
fi

# SET FC TO MPIFC. IF MPILIBS IS PRESENT, SET FC TO FC.
if test -z "$FC"; then
   FC=$MPIFC
   if test "$FC" != "$MPIFC";  then
      if test -n "$MPILIBS"; then
	 FC=$FC
      fi
   fi
fi

# FOR SANITY, CHECK THAT FILENAME EXTENSION FOR FC IS CONSISTENT WITH FC
OLDFC="$FC"
FC="$FC"

AC_COMPILE_IFELSE(
  [       subroutine oof()
             return 
          end], [],
  [AC_MSG_WARN([$FC FAILED TO COMPILE FILENAME EXTENSION $ac_ext])
  ])



FC="$OLDFC"

# CHECK HOW TO GET THE COMPILER VERSION.
AX_FC_VERSION()

AC_LANG_POP(Fortran)

# Check how to use the cpp with fortran

AC_FC_PP_DEFINE()

# CHECK IF FORTRAN COMPILER CREATES MODULE FILES IN UPPERCASE OR LOWERCASE
OLDFC="$FC"
OLDFCFLAGS="$FCFLAGS"

FC="$FC" # use MPIFC, since differences may exist between MPIFC and FC 


FC="$OLDFC"
FCFLAGS="$OLDFCFLAGS"

# CHECK HOW TO NAME MANGLE C FUNCTIONS SO THAT IT CAN BE CALLED FROM FORTRAN
OLDFC="$FC"
FC="$FC" # use MPIFC, since differences may exist between MPIFC and FC

AC_FC_WRAPPERS()

FC="$OLDFC"

# CHECK THAT THE FORTRAN COMPILER CAN CORRECTLY PROCESS THESE DIRECTIVES
# IF NOT, USE THE EXTERNAL C PREPROCESSOR
OLDFC="$FC"
FC="$FC" # use MPIFC, since differences may exist between MPIFC and FC

defineflag="-Daardvark"
if test "$OLDFC" = "xlf90"; then
   defineflag="-WF,-Daardvark"
fi   
if test "$OLDFC" = "frt"; then
   defineflag="-Wp,-Daardvark"
fi

# DEFINE VARIABLES ACCORDING TO OS AND COMPILER

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
CPPDEFS+=" -DSYS$SYSDEF"

# SET HARDCODED VARIABLES AS A LAST RESORT

# ALWAYS ENABLE CRULE IN MAKEFILE
AC_SUBST(CRULE,[.c.o])

AC_SUBST(CPPDEFS)

# INCLUDE FLAG IF NOT ALREADY SET IS MOST LIKELY -I
if test -z "$INCLUDEFLAG"; then
   INCLUDEFLAG="-I"
fi

# ARCHIVE COMMAND SIMILAR ACROSS ALL PLATFORMS 
if test -z "$AR"; then
   AR="ar cq"
fi

echo
echo Output Variables: {CC=$CC} {CFLAGS=$CFLAGS} \
{FC=$FC} {FCFLAGS=$FCFLAGS} \
{MACHDEFS=$MACHDEFS} {OPT=$OPT} {REAL8=$REAL8} \
{BIT64=$BIT64} {ENDIAN=$ENDIAN} {MPIFC=$MPIFC} \
{MPILIBS=$MPILIBS} {MPIHEADER=$MPIHEADER}  \
{INCLUDEFLAG=$INCLUDEFLAG} {INCLUDEPATH=$INCLUDEPATH} \
{AR=$AR} {BABELROOT=$BABELROOT} {COMPILER_ROOT=$COMPILER_ROOT} \
{PYTHON=$PYTHON} {PYTHONOPTS=$PYTHONOPTS} {FORT_SIZE=$FORT_SIZE} {prefix=$prefix} \
{SRCDIR=$SRCDIR}
echo

AC_OUTPUT

echo Please check the Makefile.conf
echo Have a nice day!

# test -z is true for empty strings
# test -n is true for non-empty strings




