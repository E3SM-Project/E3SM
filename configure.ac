#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT(MCT, 1.8)

# Process the following makefiles
AC_CONFIG_FILES(Makefile.conf)

# Hardcode values for specific hostnames

test "$osname"=NONE && osname=`uname -s`
test "$machinename"=NONE && machinename=`uname -m`

# Check if we are on the Earth Simulator
ES="NO"
if echo $osname | grep -i esos >/dev/null 2>&1; then
   ES="YES"
fi
if echo $osname | grep -i hp-ux >/dev/null 2>&1; then   
   if test "$ac_hostname" = "moon"; then 
      ES="YES"
      # Tells configure not to run any tests that require execution
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
   if test -z "$FCFLAGS"; then
      FCFLAGS="-EP -Wf'-pvctl fullmsg -L fmtlist transform map'"
   fi
   if test -z "$OPT"; then
      OPT="-C vopt"
   fi      
   MACHDEFS="-DESVEC $MACHDEFS"
fi

# Check if we are on NCSA IA-64 cluster titan
if test "$ac_hostname" = "user02"; then
   if test -z "$MPILIBS"; then
   if test -z "$MPIHEADER"; then
      echo "Using preset MPI configuration for titan"
      MPILIBS="-L/usr/local/vmi/mpich/lib/intel -lmpich -lvmi -ldl -lpthread"
      MPIHEADER="/usr/local/vmi/mpich/include"
   fi
   fi
fi

# Check for the C compiler
AC_PROG_CC([cc])

# Preset FCFLAGS and F90FLAGS
if test -n "$F90FLAGS"; then
   F90TEMPFLAGS="$F90FLAGS"
   F90FLAGS=""
if test -n "$FCFLAGS"; then
   F90FLAGS="$FCFLAGS"
fi
else
if test -n "$FCFLAGS"; then
   F90TEMPFLAGS=""
   F90FLAGS="$FCFLAGS"
fi
fi
OLDF90FLAGS="$F90FLAGS"

# Check for the Fortran compiler. 
AC_PROG_F90

# Check for MPI libraries. 
AC_LANG_PUSH(Fortran 90)

AC_ARG_VAR(MPILIBS,[MPI library command line invocation])
AC_ARG_VAR(MPIHEADER,[MPI header include path])
if test -n "$MPIHEADER"; then
   F90FLAGS="$MPIHEADER"
fi
ACX_MPI
F90FLAGS="$OLDF90FLAGS"

# A hack to fix ACX_MPI to get MPILIBS to be an empty string
if test "$MPILIBS" = " "; then
   MPILIBS=""
fi

# Check how to get the compiler version.
# This test does not work for xlf compilers.
if test "$F90" != "xlf90"; then
   _AC_PROG_F90_VERSION
fi
f90_version="$ac_f90_version_output"

AC_LANG_POP(Fortran 90)

# Checks for compiler characteristics.

# Check for endianness
AC_C_BIGENDIAN

# Check if fortran compiler creates modules files
# in uppercase or lowercase
if test "$F90" = "frt"; then # -Am flag needed for frt to gerate mod files
   F90FLAGS="-Am"
fi
if test "$F90" = "ftn"; then # -em flag needed for ftn to gerate mod files
   F90FLAGS="-em"
fi
AC_F90_MOD_UPPERCASE
F90FLAGS="$OLDF90FLAGS"

# Check the fortran name mangling scheme
AC_F90_C_NAME_MANGLING
case $ac_cv_f90_mangling in
  "lower case, underscore") 
	 name_mangling="FORTRAN_UNDERSCORE_" ;;
  "lower case, no underscore")
	 name_mangling="FORTRAN_SAME" ;;
  "upper case, no underscore")  
	 name_mangling="FORTRAN_CAPS_" ;;
  *) name_mangling="UNKNOWN" ;;
esac

# Set FPP and FPPFLAGS
AC_ARG_VAR(FPP,C-preprocessor for Fortran source code)
AC_ARG_VAR(FPPFLAGS,C-preprocessing flags for Fortran source code)

if test -z "$FPP"; then
   AC_CHECK_PROGS(FPP,[cpp],$CPP)
fi
if test -z "$FPPFLAGS"; then
   FPPFLAGS="-P -C -N -traditional"
fi

# Check that the Fortran compiler can correctly process these directives
# If not, use the C preprocessor.

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

if test "$ac_cv_f90_fpp_compiler" = "yes"; then
   AC_SUBST(FRULE,[.F.o])
   AC_SUBST(F90RULE,[.F90.o])
   AC_SUBST(FRULECPP,[.FRULECPP])
   AC_SUBST(F90RULECPP,[.F90RULECPP])
fi

if test "$ac_cv_f90_fpp_compiler" = "no"; then
if test "$ac_cv_f90_fpp_external" = "no"; then
   AC_MSG_WARN(Fortran Source Preprocessing has failed)
fi
if test "$ac_cv_f90_fpp_external" = "yes"; then
   AC_SUBST(FRULE,[.FRULE])
   AC_SUBST(F90RULE,[.F90RULE])
   AC_SUBST(FRULECPP,[.F.o])
   AC_SUBST(F90RULECPP,[.F90.o])
fi
fi

# Always enable the CRULE

AC_SUBST(CRULE,[.c.o])

# Checks OS, compiler and machine type

echo "Hostname=$ac_hostname"
echo "Machine=$machinename"
echo "OS=$osname"

# Set F90FLAGS and MACHDEFS

F90FLAGS="$F90TEMPFLAGS"
AC_ARG_VAR(MACHDEFS,MCT compiler and OS flags)

# Check OS name

if echo $osname | grep -i aix >/dev/null 2>&1; then
   if test "$ac_cv_f90_fpp_compiler" = "yes"; then
      MACHDEFS="-WF,-DSYSAIX"
   else
      MACHDEFS="-DSYSAIX $MACHDEFS"
   fi
fi
if echo $osname | grep -i darwin >/dev/null 2>&1; then
   MACHDEFS="-DSYSDARWIN $MACHDEFS"
fi
if echo $osname | grep -i unix_system_v >/dev/null 2>&1; then
   if test "$ac_cv_f90_fpp_compiler" = "yes"; then
      MACHDEFS="-Wp,-DSYSUNIXSYSTEMV"
   else
      MACHDEFS="-DSYSUNIXSYSTEMV $MACHDEFS"
   fi
fi
if echo $osname | grep -i irix64 >/dev/null 2>&1; then
   MACHDEFS="-DSYSIRIX64 $MACHDEFS"
fi
if echo $osname | grep -i irix >/dev/null 2>&1; then
if echo $MACHDEFS | grep -i irix64 >/dev/null 2>&1; then
MACHDEFS="$MACHDEFS" # do nothing
else
MACHDEFS="-DSYSIRIX $MACHDEFS"
fi
fi
if echo $osname | grep -i linux >/dev/null 2>&1; then
   MACHDEFS="-DSYSLINUX $MACHDEFS"
fi
if echo $osname | grep -i osf1 >/dev/null 2>&1; then
   MACHDEFS="-DSYSOSF1 $MACHDEFS"
fi
if echo $osname | grep -i super >/dev/null 2>&1; then
   MACHDEFS="-DSYSSUPERUX $MACHDEFS"
fi
if echo $osname | grep -i sun >/dev/null 2>&1; then
   MACHDEFS="-DSYSSUNOS $MACHDEFS"
fi
if echo $osname | grep -i t3e >/dev/null 2>&1; then
   MACHDEFS="-DSYST3E $MACHDEFS"
fi
if echo $osname | grep -i unicos >/dev/null 2>&1; then
   MACHDEFS="-DSYSUNICOS $MACHDEFS"
fi

# Set FC to MPIF90. If MPILIBS is present, set FC to F90.
AC_ARG_VAR(FC,The Fortran compiler)
if test -z "$FC"; then
   FC=$MPIF90
   if test "$F90" != "$MPIF90";  then
      if test -n "$MPILIBS"; then
	 FC=$F90
      fi
   fi
fi

# Set the following variables by checking the compiler name
AC_ARG_VAR(FCFLAGS,Fortran compiler flags for all source files)
AC_ARG_VAR(ALLCFLAGS,Customized C source compilation flags)
AC_ARG_VAR(OPT,Fortran compiler flag for optimization level) 
AC_ARG_VAR(REAL8,[Fortran compiler flag for setting the default REAL size to REAL(KIND=8)])
AC_ARG_VAR(BIT64,Fortran compiler flag for generating 64-bit objects)
AC_ARG_VAR(ENDIAN,Fortran compiler flag for converting big-endian to little-endian)
AC_ARG_VAR(INCLUDEFLAG,Fortran compiler flag for specifying module search path)
AC_ARG_VAR(INCLUDEPATH,Additional library and module paths)
AC_ARG_VAR(AR,Archive command)

# Set AR 
if test -z "$AR"; then
   AR="ar cq"
fi

# Check compiler name
case $F90 in
xlf90) 
echo "Fortran Compiler is XLF"
if test -z "$REAL8"; then
   REAL8="-qrealsize=8"
fi
if test "$MACHDEFS" = "-WF,-DSYSAIX"; then
   MACHDEFS="$MACHDEFS,-DCPRXLF"
else
   MACHDEFS="-DCPRXLF $MACHDEFS"
fi
if test -z "$FCFLAGS"; then
   FCFLAGS="-qarch=auto"
fi
if test -z "$OPT"; then
   OPT="-O2"
fi
if test -z "$F90FLAGS"; then
  if test "$ac_cv_f90_fpp_compiler" = "yes"; then
     F90FLAGS="-qsuffix=f=F90:cpp=F90"
  fi
  if test "$ac_cv_f90_fpp_compiler" = "no"; then
     if test "$ac_cv_f90_fpp_external" = "yes"; then
	F90FLAGS="-qsuffix=f=f90"
     fi
  fi
fi;;	
pgf90) 
echo "Fortran Compiler is Portland Group"
MACHDEFS="-DCPRPGI $MACHDEFS"
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
MACHDEFS="-DCPRCRAY $MACHDEFS"
if test -z "$FCFLAGS"; then 
   FCFLAGS="-em -dy -m0 -M1068,1132,399"
fi
if test -z "$INCLUDEFLAG"; then
   INCLUDEFLAG="-p"
fi
if test -z "$OPT"; then
   OPT="-Oscalar3 -Ovector3 -Ostream3"
fi;;
frt)
echo "Fortran Compiler is UXP/V"
if test "$MACHDEFS" = "-Wp,-DSYSUNIXSYSTEMV"; then
   MACHDEFS="$MACHDEFS,-DSYSFUJITSU"
else
   MACHDEFS="-DSYSFUJITSU $MACHDEFS"
fi
if test -z "$FCFLAGS"; then
   FCFLAGS="-Am -X9"
fi
if test -z "$BIT64"; then
   BIT64="-KA64"
fi
if test -z "$REAL8"; then
   REAL8="-Ad"
fi
echo "Suggested additional vectorization flags: -Wv,-s5,-t3,-noalias,-ilfunc,-md";;
epcf90) 
"Everest has never tested this fortran compiler :("
MACHDEFS="-DCPREPC $MACHDEFS";;
lf95) 
echo "Fortran Compiler is Lahey"
MACHDEFS="-DCPRLAHEY $MACHDEFS";;
ifc | efc)
echo "Fortran Compiler is Intel"
MACHDEFS="-DCPRINTEL $MACHDEFS"
if test -z "$REAL8"; then
   REAL8="-r8"
fi
if test -z "$FCFLAGS"; then
   FCFLAGS="-w"
fi
if test -z "$OPT"; then
   OPT="-O2"
fi
echo "Intel Compilers must set the environment variable F_UFMTENDIAN=big to do endian conversion";;
g95)
echo "Fortran Compiler is GNU"
MACHDEFS="-DCPRGNU $MACHDEFS";;
esac

if echo $f90_version | grep -i nag >/dev/null 2>&1; then
   echo "Fortran Compiler is NAG"
   MACHDEFS="-DCPRNAG $MACHDEFS"
   if test -z "$FCFLAGS"; then
      FCFLAGS="-dusty -kind=byte"
   fi
   if test -z "$OPT"; then
      OPT="-O2"
   fi
fi
if echo $f90_version | grep -i absoft >/dev/null 2>&1; then
   echo "Fortran Compiler is Absoft"
   MACHDEFS="-DCPRABSOFT $MACHDEFS"
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
   MACHDEFS="-DCPRWORKSHOP $MACHDEFS"
   if test -z "$INCLUDEFLAG"; then
      INCLUDEFLAG="-M"
   fi
fi
if echo $f90_version | grep -i mipspro >/dev/null 2>&1; then
   echo "Fortran Compiler is MIPSPro"
   MACHDEFS="-DCPRMIPSPRO $MACHDEFS"
   if test -z "$OPT"; then
      OPT="-O3"
   fi
   if test -z "$ALLCFLAGS"; then
      ALLCFLAGS="-64"
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
   MACHDEFS="-DCPRCOMPAQ $MACHDEFS"
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
   MACHDEFS="-DCPRCOMPAQ $MACHDEFS"
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
   MACHDEFS="-DCPRSX $MACHDEFS"
   if test -z "$FCFLAGS"; then
      FCFLAGS="-EP -Wf'-pvctl noassoc'"
   fi
   if test -z "$OPT"; then
      OPT="-Chopt"
   fi
fi

case $name_mangling in
  FORTRAN_CAPS_) ALLCFLAGS="-DFORTRAN_CAPS_ -O $ALLCFLAGS";;
  FORTRAN_SAME)  ALLCFLAGS="-DFORTRAN_SAME -O $ALLCFLAGS";;
  FORTRAN_UNDERSCORE_) ALLCFLAGS="-O $ALLCFLAGS";;
  *) AC_MSG_WARN([unusual Fortran 90 name-mangling scheme])
     ALLCFLAGS="-O $ALLCFLAGS";;
esac

if test -z "$INCLUDEFLAG"; then
   INCLUDEFLAG="-I"
fi

AC_OUTPUT

echo
echo Output Variables: {CC=$CC} {ALLCFLAGS=$ALLCFLAGS} \
{FPP=$FPP} {FPPFLAGS=$FPPFLAGS} {FC=$FC} {F90=$F90} \
{FCFLAGS=$FCFLAGS} {F90FLAGS=$F90FLAGS} {OPT=$OPT} \
{REAL8=$REAL8} {BIT64=$BIT64} {MPIF90=$MPIF90} \
{MPILIBS=$MPILIBS} {MPIHEADER=$MPIHEADER} {MACHDEFS=$MACHDEFS} \
{ENDIAN=$ENDIAN} {INCLUDEFLAG=$INCLUDEFLAG} \
{INCLUDEPATH=$INCLUDEPATH} {AR=$AR} \
{prefix=$prefix}

echo
echo Please check the Makefile.conf
echo Have a nice day!

# test -z is true for empty strings
# test -n is true for non-empty strings
