#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT(MCT, 2.0)

# PROCESS THE FOLLOWING MAKEFILES

AC_CONFIG_FILES(Makefile.conf)

# DECLARE PACKAGE OPTIONS

AC_ARG_ENABLE(selectedrealkind,
AC_HELP_STRING([--enable-selectedrealkind],
[define single precision and double precision numbers using the selected_real_kind function. Default uses the kind inquiry function.]),
[SRKDEF="SELECTEDREALKIND"]
)

# DECLARE THE FOLLOWING PRECIOUS VARIABLES

AC_ARG_VAR(MPILIBS,[MPI library command line invocation])
AC_ARG_VAR(MPIHEADER,[MPI header include path])
AC_ARG_VAR(FPP,C-preprocessor for Fortran source code)
AC_ARG_VAR(FPPFLAGS,C-preprocessing flags for Fortran source code)
AC_ARG_VAR(MACHDEFS,MCT compiler and OS flags)
AC_ARG_VAR(FC,The Fortran compiler)
AC_ARG_VAR(FCFLAGS,User-defined Fortran compiler flags)
AC_ARG_VAR(ALLCFLAGS,Customized C source compilation flags)
AC_ARG_VAR(OPT,Fortran compiler flag for optimization level) 
AC_ARG_VAR(REAL8,[Fortran compiler flag for setting the default REAL size to REAL(KIND=8)])
AC_ARG_VAR(BIT64,Fortran compiler flag for generating 64-bit objects)
AC_ARG_VAR(ENDIAN,Fortran compiler flag for converting big-endian to little-endian)
AC_ARG_VAR(INCLUDEFLAG,Fortran compiler flag for specifying module search path)
AC_ARG_VAR(INCLUDEPATH,Additional library and module paths)
AC_ARG_VAR(AR,Archive command)

# SET TEMPORARY VARIABLES

# OS AND PLATFORM NAME
test "$osname"=NONE && osname=`uname -s`
test "$machinename"=NONE && machinename=`uname -m`

# SAVE F90FLAGS IN CASE ITS SET IN THE COMMAND LINE
OLDF90FLAGS="$F90FLAGS"

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
	    MPIHEADER="/usr/local/vmi/mpich/include"
         fi
      fi
   fi
fi

# START TESTS

# CHECK FOR THE C COMPILER
AC_PROG_CC([cc])

# CHECK FOR THE FORTRAN COMPILER
AC_PROG_F90

# CHECK FOR MPI LIBRARIES
AC_LANG_PUSH(Fortran 90)

if test -n "$MPIHEADER"; then
   F90FLAGS="$F90FLAGS $MPIHEADER"
fi

ACX_MPI

F90FLAGS="$OLDF90FLAGS"

# A HACK TO FIX ACX_MPI TO GET MPILIBS TO BE AN EMPTY STRING
if test "$MPILIBS" = " "; then
   MPILIBS=""
fi

# SET FC TO MPIF90. IF MPILIBS IS PRESENT, SET FC TO F90.
if test -z "$FC"; then
   FC=$MPIF90
   if test "$F90" != "$MPIF90";  then
      if test -n "$MPILIBS"; then
	 FC=$F90
      fi
   fi
fi

# CHECK HOW TO GET THE COMPILER VERSION.
# THIS TEST DOES NOT WORK FOR XLF COMPILERS.
if test "$F90" != "xlf90"; then
   _AC_PROG_F90_VERSION
fi
f90_version="$ac_f90_version_output"

AC_LANG_POP(Fortran 90)

# CHECK FOR BYTE ORDERING
AC_C_BIGENDIAN

# CHECK IF FORTRAN COMPILER CREATES MODULE FILES IN UPPERCASE OR LOWERCASE
if test "$F90" = "frt"; then # -Am flag needed for frt to gerate mod files
   F90FLAGS="$F90FLAGS -Am"
fi
if test "$F90" = "ftn"; then # -em flag needed for ftn to gerate mod files
   F90FLAGS="$F90FLAGS -em"
fi

AC_F90_MOD_UPPERCASE

F90FLAGS="$OLDF90FLAGS"

# CHECK HOW TO NAME MANGLE C FUNCTIONS SO THAT IT CAN BE CALLED FROM FORTRAN
AC_F90_C_NAME_MANGLING
case $ac_cv_f90_mangling in
  "lower case, underscore") 
	 name_mangling="FORTRAN_UNDERSCORE_";;
  "lower case, no underscore")
	 name_mangling="FORTRAN_SAME";;
  "upper case, no underscore")  
	 name_mangling="FORTRAN_CAPS_";;
  *) 
         name_mangling="FORTRAN_MANGLING_UNKNOWN" 
	 AC_MSG_WARN([unknown Fortran 90 name-mangling scheme]);;
esac

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
defineflag="-Daardvark"
if test "$F90" = "xlf90"; then
   defineflag="-WF,-Daardvark"
fi   
if test "$F90" = "frt"; then
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

else

AC_MSG_WARN(cpp is not defined in your path)

fi

# ACTIVATE F90RULE IF FORTRAN COMPILER SUCCEEDS AT PREPROCESSING
# ACTIVATE F90RULECPP IF FOTRAN COMPILER FAILS AT PREPROCESSING
if test "$ac_cv_f90_fpp_compiler" = "yes"; then
   AC_SUBST(F90RULE,[.F90.o])
   AC_SUBST(F90RULECPP,[.F90RULECPP])
fi

if test "$ac_cv_f90_fpp_compiler" = "no"; then
if test "$ac_cv_f90_fpp_external" = "no"; then
   AC_MSG_WARN(Fortran Source Preprocessing has failed)
fi
if test "$ac_cv_f90_fpp_external" = "yes"; then
   AC_SUBST(F90RULE,[.F90RULE])
   AC_SUBST(F90RULECPP,[.F90.o])
fi
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

# CHECK COMPILER NAME
case $F90 in

xlf90)
   echo "Fortran Compiler is XLF"
   CPRDEF="XLF"
   if test -z "$REAL8"; then
      REAL8="-qrealsize=8"
   fi
   if test -z "$OPT"; then
      OPT="-O2"
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
   fi;;

pgf90)
   CPRDEF="PGI"
   echo "Fortran Compiler is Portland Group"
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
   fi;;

ftn)
   echo "Fortran Compiler is Cray"
   CPRDEF="CRAY" 
   if test -z "$F90FLAGS"; then 
      F90FLAGS="-em -dy -m0 -M1068,1132,399"
   fi
   if test -z "$INCLUDEFLAG"; then
      INCLUDEFLAG="-p"
   fi
   if test -z "$OPT"; then
      OPT="-Oscalar3 -Ovector3 -Ostream3"
   fi;;

frt)
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
   fi;;

epcf90)
   echo "Everest has never tested the epcf90 compiler :("
   CPRDEF="EPC";;

lf95)
   echo "Fortran Compiler is Lahey"
   CPRDEF="LAHEY";;

ifc)
   echo "Fortran Compiler is Intel IA-32"
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
   fi;;

efc)
   echo "Fortran Compiler is Intel Itanium"
   CPRDEF="INTEL"
   echo "Intel efc compiler must set the environment variable F_UFMTENDIAN=big to do endian conversion"
   if test -z "$REAL8"; then
      REAL8="-r8"
   fi
   if test -z "$F90FLAGS"; then
      F90FLAGS="-w -ftz"
   fi
   if test -z "$OPT"; then
      OPT="-O2"
   fi;;

ifort)
   echo "Fortran Compiler is Intel Itanium"
   CPRDEF="INTEL"
   if test -z "$REAL8"; then
      REAL8="-r8"
   fi
   if test -z "$F90FLAGS"; then
      F90FLAGS="-w -ftz -assume byterecl -convert big_endian"
   fi
   if test -z "$OPT"; then
      OPT="-O2"
   fi;;

g95)
   echo "Fortran Compiler is GNU"
   CPRDEF="GNU";;

esac

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
   if test -z "$INCLUDEFLAG"; then
      INCLUDEFLAG="-p"
   fi
   if test -z "$OPT"; then
      OPT="-O2"
   fi
fi

if echo $f90_version | grep -i workshop >/dev/null 2>&1; then
   echo "Fortran Compiler is Workshop"
   CPRDEF="WORKSHOP"
   if test -z "$INCLUDEFLAG"; then
      INCLUDEFLAG="-M"
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
   MPILIBS="$MPILIBS -lelan"
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
   MPILIBS="$MPILIBS -lelan"
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

# SET THE PREPROCESSOR DEFINE FLAGS
if test -z "$MACHDEFS"; then

   if test -z "$SYSDEF"; then
      echo "WARNING: OPERATING SYSTEM UNKNOWN"
      SYSDEF="UNKNOWNOS"
   fi

   MACHDEFS="-DSYS$SYSDEF -DCPR$CPRDEF"
   if test -n "$SRKDEF"; then
      MACHDEFS="$MACHDEFS -D$SRKDEF"
   fi

   if test "$F90" = "xlf90"; then
      if test "$ac_cv_f90_fpp_compiler" = "yes"; then
         MACHDEFS="-WF,-DSYS$SYSDEF,-DCPR$CPRDEF"
	 if test -n "$SRKDEF"; then
            MACHDEFS="$MACHDEFS,-D$SRKDEF"
         fi
      fi
   fi

   if test "$F90" = "frt"; then
      if test "$ac_cv_f90_fpp_compiler" = "yes"; then
         MACHDEFS="-Wp,-DSYS$SYSDEF,-DCPR$CPRDEF"
         if test -n "$SRKDEF"; then
            MACHDEFS="$MACHDEFS,-D$SRKDEF"
         fi
      fi
   fi

fi

if test -z "$ALLCFLAGS"; then
   ALLCFLAGS="-D$name_mangling -O"
   if test -n "$EXTRACFLAGS"; then
      ALLCFLAGS="-D$name_mangling -O $EXTRACFLAGS"
   fi
fi

# SET HARDCODED VARIABLES AS A LAST RESORT

# ALWAYS ENABLE CRULE IN MAKEFILE
AC_SUBST(CRULE,[.c.o])

# INCLUDE FLAG IF NOT ALREADY SET IS MOST LIKELY -I
if test -z "$INCLUDEFLAG"; then
   INCLUDEFLAG="-I"
fi

# ARCHIVE COMMAND SIMILAR ACROSS ALL PLATFORMS 
if test -z "$AR"; then
   AR="ar cq"
fi

echo
echo Output Variables: {CC=$CC} {ALLCFLAGS=$ALLCFLAGS} \
{FPP=$FPP} {FPPFLAGS=$FPPFLAGS} {FC=$FC} \
{F90=$F90} {FCFLAGS=$FCFLAGS} {F90FLAGS=$F90FLAGS} \
{MACHDEFS=$MACHDEFS} {OPT=$OPT} {REAL8=$REAL8} \
{BIT64=$BIT64} {ENDIAN=$ENDIAN} {MPIF90=$MPIF90} \
{MPILIBS=$MPILIBS} {MPIHEADER=$MPIHEADER}  \
{INCLUDEFLAG=$INCLUDEFLAG} {INCLUDEPATH=$INCLUDEPATH} \
{AR=$AR} {prefix=$prefix}
echo

AC_OUTPUT

echo Please check the Makefile.conf
echo Have a nice day!

# test -z is true for empty strings
# test -n is true for non-empty strings




