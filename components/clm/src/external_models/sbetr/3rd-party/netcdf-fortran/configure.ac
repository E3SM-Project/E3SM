#                                               -*- Autoconf -*-
## Process this file with autoconf to produce a configure script.

# This is part of Unidata's netCDF package. Copyright 2011, see the
# COPYRIGHT file for more information.

# Running autoconf on this file will trigger a warning if
# autoconf is not at least the specified version.
AC_PREREQ([2.59])

# Initialize with name, version, and support email address.
AC_INIT([netCDF-Fortran], [4.4.3], [support-netcdf@unidata.ucar.edu])

# Create the VERSION file, which contains the package version from
# AC_INIT.
echo -n AC_PACKAGE_VERSION>VERSION
AC_SUBST(PACKAGE_VERSION)

AC_MSG_NOTICE([netCDF-Fortran AC_PACKAGE_VERSION])

# Keep libtool macros in an m4 directory.
AC_CONFIG_MACRO_DIR([m4])

# Find out about the host we're building on.
AC_CANONICAL_HOST

# Find out about the target we're building for.
AC_CANONICAL_TARGET

# This call is required by automake.
AM_INIT_AUTOMAKE([foreign dist-zip subdir-objects])
AM_MAINTAINER_MODE()

# Check for the existence of this file before proceeding.
AC_CONFIG_SRCDIR([libsrc/fort-attio.c])

AC_MSG_NOTICE([checking user options])

# Does the user want to run extra tests with valgrind?
AC_MSG_CHECKING([whether extra valgrind tests should be run])
AC_ARG_ENABLE([valgrind-tests],
              [AS_HELP_STRING([--enable-valgrind-tests],
                              [build with valgrind-tests (valgrind is required, static builds only)])])
test "x$enable_valgrind_tests" = xyes || enable_valgrind_tests=no
AC_MSG_RESULT($enable_valgrind_tests)

# Does the user want to run extra parallel tests when parallel netCDF-4 is built?
AC_MSG_CHECKING([whether parallel IO tests should be run])
AC_ARG_ENABLE([parallel-tests],
              [AS_HELP_STRING([--enable-parallel-tests],
                              [Run extra parallel IO tests. Ignored if \
                               netCDF-4 is not enabled, or built on a system \
                               without parallel I/O support.])])
test "x$enable_parallel_tests" = xyes || enable_parallel_tests=no
AC_MSG_RESULT($enable_parallel_tests)

# Does the user want to do some extra tests?
AC_MSG_CHECKING([whether netCDF extra tests should be run (developers only)])
AC_ARG_ENABLE([extra-tests],
              [AS_HELP_STRING([--enable-extra-tests],
                              [run some extra tests that may not pass because of known issues])])
test "x$enable_extra_tests" = xyes || enable_extra_tests=no
AC_MSG_RESULT($enable_extra_tests)
if test "x$enable_extra_tests" = xyes; then
   AC_DEFINE([EXTRA_TESTS], [1], [if true, run extra tests which may not work yet])
fi
AM_CONDITIONAL(EXTRA_TESTS, [test x$enable_extra_tests = xyes])

#####
# Determine if we want to enable doxygen-generated documentation.
#####
AC_ARG_ENABLE([doxygen],
	[AS_HELP_STRING([--enable-doxygen],
	[Enable generation of documentation with doxygen.])])
    test "x$enable_doxygen" = xyes || enable_doxygen=no
AM_CONDITIONAL([BUILD_DOCS], [test "x$enable_doxygen" = xyes])

# Does the user want to generate dot-based documentation?
AC_ARG_ENABLE([dot],
	[AS_HELP_STRING([--enable-dot],
	[Use dot (provided by graphviz) to generate charts and graphs in the doxygen-based documentation.])])
	test "x$enable_dot" = xyes || enable_dot=no

# Does the user want to generate internal documentation?
AC_ARG_ENABLE([internal-docs],
  [AS_HELP_STRING([--enable-internal-docs],
    [Include documentation of library internals. This is of interest only to those developing the netCDF library.])])
test "x$enable_internal_docs" = xyes || enable_internal_docs=no
AC_SUBST([BUILD_INTERNAL_DOCS], [$enable_internal_docs])

# Is doxygen installed? If so, have configure construct the Doxyfile.
AC_CHECK_PROGS([DOXYGEN], [doxygen])
if test -z "$DOXYGEN"; then
   AC_MSG_WARN([Doxygen not found - documentation will not be built])
fi

# Is graphviz/dot installed? If so, we'll use dot to create
# graphs in the documentation.
AC_CHECK_PROGS([DOT], [dot])
if test -z "$DOT"; then
   AC_MSG_WARN([dot not found - will use simple charts in documentation])
   HAVE_DOT=NO
elif test "x$enable_dot" = xno; then
   HAVE_DOT=NO
else
   HAVE_DOT=YES
fi
# If we have doxygen, and it's enabled, then process the file.
if test "x$enable_doxygen" != xno; then
   if test -n "$DOXYGEN"; then
        AC_SUBST(HAVE_DOT)
        AC_CONFIG_FILES([docs/Doxyfile])
   fi
# Note: the list of files to input to doxygen
# has been moved to docs/Doxyfile.in so
# that make distcheck works correctly.
# Any new inputs should be inserted into
# docs/Doxyfile.in and possibley docs/Makefile.am
fi

#####
# End doxygen-related documentation block.
#####

# Find the C compiler.
AC_MSG_NOTICE([finding C compiler])

AC_PROG_CC
AM_PROG_CC_C_O
AC_C_CONST

# If F03 is set, FC isn't, then set FC to whatever F03 was set
if test "${F03+set}" = set -a "${FC+set}" != set; then
   FC=$F03
fi

# If F90 is set, FC isn't, then set FC to whatever F90 was set to.
if test "${F90+set}" = set -a "${FC+set}" != set; then
   FC=$F90
fi

# Find fortran compiler. If FC isn't set after AC_PROG_FC, we couldn't
# find a fortran compiler. If the user requested the fortran API, and
# has told us not to recover from a missing compiler, then bail out
# right here.
AC_MSG_NOTICE([finding Fortran compiler])
AC_PROG_FC
AC_PROG_F77($FC)

# If we couldn't find an F90 compiler, deal with it.
if test "${FC+set}" = set && test "x$FC" = x; then
   AC_MSG_FAILURE([Can't find F90 compiler.])
fi

# Set the default fortran builds; default is to build f03
nc_build_f90=no
nc_build_f03=yes

# Does the user want to test for Fortran 2003 compiler.
# Added to get around NAG evaluation compiler error turning off F03
AC_MSG_CHECKING([whether Fortran compiler is checked for ISO_C_BINDING support])
AC_ARG_ENABLE([f03-compiler-check],
              [AS_HELP_STRING([--disable-f03-compiler-check],
                              [disable check of ISO_C_BINDING support in Fortran compiler])])
test "x$enable_f03_compiler_check" = xno || enable_f03_compiler_check=yes
AC_MSG_RESULT($enable_f03_compiler_check)

# next check to see if the Fortran compiler will support
# ISO_C_BINDING
nc_iso_c_binding=yes

# Some f90 compilers change the case of the mod file names. Others
# require special flags to be set to deal with .f90 files. Find out
# about this compiler.
AC_LANG_PUSH(Fortran)
AC_FC_SRCEXT(f90)
AC_LANG_POP(Fortran)
AX_F90_MODULE_FLAG
AC_SUBST(MOD_FLAG, ${ax_cv_f90_modflag})
if test "x$enable_f03_compiler_check" = xyes; then
  AC_MSG_CHECKING([if Fortran compiler supports Fortran 2003 ISO_C_BINDING])
  cat <<EOF >conftest.f90
module conftest_module
USE ISO_C_BINDING
end module conftest_module
EOF
  doit='$FC -c ${FCFLAGS} ${FCFLAGS_f90} conftest.f90'
  if AC_TRY_EVAL(doit); then
     nc_iso_c_binding=yes
  else
     nc_iso_c_binding=no
  fi
  AC_MSG_RESULT([$nc_iso_c_binding])

  AC_MSG_CHECKING([if Fortran compiler supports Fortran 2008 ISO_FORTRAN_ENV additions])
  cat <<EOF >conftest.f90
Program test
USE ISO_FORTRAN_ENV, ONLY: REAL32, REAL64, INT8, INT16, INT32, INT64
End Program
EOF
  doit='$FC -c ${FCFLAGS} ${FCFLAGS_f90} conftest.f90'
  if AC_TRY_EVAL(doit); then
     nc_f2008_is0_env_additions=yes
     AC_DEFINE([HAVE_F2008], 1, [if true, use Fortran 2008 ISO_FORTRAN_ENV additions])
  else
     nc_f2008_is0_env_additions=no
  fi
  AC_MSG_RESULT([$nc_f2008_is0_env_additions])

  AC_MSG_CHECKING([if Fortran compiler supports TS29113 standard extension])
  cat <<EOF >conftest.f90
Program test
USE ISO_C_BINDING, ONLY: C_PTRDIFF_T
End Program
EOF
  doit='$FC -c ${FCFLAGS} ${FCFLAGS_f90} conftest.f90'
  if AC_TRY_EVAL(doit); then
     nc_ts29113_support=yes
     AC_DEFINE([HAVE_TS29113_SUPPORT], 1, [if true, Fortran compiler supports TS29113 extension])
  else
     nc_ts29113_support=no
  fi
  AC_MSG_RESULT([$nc_ts29113_support])

  if test "x$nc_iso_c_binding" = xno; then
    if test "x$enable_compiler_recover" = xyes; then
      nc_build_f90=yes
      AC_MSG_WARN([ISO_C_BINDING not supported - defaulting to F90 interfaces])
    else
      AC_MSG_FAILURE([Can't use build F03 interfaces and cannot do compiler recovery])
    fi
  fi
fi


## Do we want to use the fortran 90 wrappers instead of 2003 native code?
AC_MSG_CHECKING([whether F03 native code is desired])
AC_ARG_ENABLE([f03],
              [AS_HELP_STRING([--disable-f03],
                              [suppress netCDF Fortran 2003 native code])])
test "x$enable_f03" = xno || enable_f03=yes
AC_MSG_RESULT([$enable_f03])

if test "x$enable_f03" = xyes ; then nc_build_f03=yes; fi

# Decide on f90 versus 2003
if test "x$nc_iso_c_binding" = xno; then
  #Cannot build f03 native code
  AC_MSG_NOTICE([Fortran 2003 interface not supported; Using fortran 90 interface])
  nc_build_f90=yes
fi

# Guarantee build state; f90 currently overrides f03
if test "x$nc_build_f90" = xyes ; then
  nc_build_f03=no
else
  nc_build_f03=yes
fi

# Does the user want to check into fortran type information?
AC_MSG_CHECKING([whether fortran type sizes should be checked])
AC_ARG_ENABLE([fortran-type-check],
              [AS_HELP_STRING([--disable-fortran-type-check],
                              [cause the Fortran type sizes checks to be skipped])])
test "x$enable_fortran_type_check" = xno || enable_fortran_type_check=yes
AC_MSG_RESULT($enable_fortran_type_check)

# Does the user want to run tests for large files (> 2GiB)?
AC_MSG_CHECKING([whether large file (> 2GB) tests should be run])
AC_ARG_ENABLE([large-file-tests],
              [AS_HELP_STRING([--enable-large-file-tests],
                              [Run tests which create very large data files (~13 GB disk space
                              required, but it will be recovered when tests are complete). See
                              option --with-temp-large to specify temporary directory])])
test "x$enable_large_file_tests" = xyes || enable_large_file_tests=no
AC_MSG_RESULT($enable_large_file_tests)
AM_CONDITIONAL(LARGE_FILE_TESTS, [test x$enable_large_file_tests = xyes])
if test "x$enable_large_file_tests" = xyes; then
   AC_DEFINE([LARGE_FILE_TESTS], [1], [do large file tests])
fi

# Does the user want to run benchmarks?
AC_MSG_CHECKING([whether benchmarks should be run (experimental)])
AC_ARG_ENABLE([benchmarks],
              [AS_HELP_STRING([--enable-benchmarks],
                              [Run benchmarks. This is an experimental feature.])])
test "x$enable_benchmarks" = xyes || enable_benchmarks=no
AC_MSG_RESULT($enable_benchmarks)
AM_CONDITIONAL(BUILD_BENCHMARKS, [test x$enable_benchmarks = xyes])

# If the env. variable TEMP_LARGE is set, or if
# --with-temp-large=<directory>, use it as a place for the large
# (i.e. > 2 GiB) files created during the large file testing.
AC_MSG_CHECKING([where to put large temp files if large file tests are run])
AC_ARG_WITH([temp-large],
            [AS_HELP_STRING([--with-temp-large=<directory>],
                            [specify directory where large files (i.e. >2 GB) \
                            will be written, if large files tests are run with
                            --enable-large-file-tests])],
            [TEMP_LARGE=$with_temp_large])
TEMP_LARGE=${TEMP_LARGE-.}
AC_MSG_RESULT($TEMP_LARGE)
#AC_SUBST(TEMP_LARGE)
AC_DEFINE_UNQUOTED([TEMP_LARGE], ["$TEMP_LARGE"], [Place to put very large netCDF test files.])


# Turn off building old F77 interfaces if nc_build_f03 is still yes

test "x$nc_build_f03" = xyes && nc_build_f77=no

if test "x$nc_build_f90" = xno; then
   F90=
   # If F77 is set, and FC isn't, then set FC to whatever F77 was set to.
   AS_IF([test x${F77+set} = xset -a x${FC+set} = x],[FC=$F77])

fi

# See if the fortran 90 build is desired. If so, set some stuff
if test "x$nc_build_f90" = xyes; then
   AC_MSG_NOTICE([setting up Fortran 90])
   if test "${F90+set}" != set; then
      F90=$FC
   fi
   AC_MSG_CHECKING([what FC is set to])
   AC_MSG_RESULT([$FC])

   if test "x$nc_build_f03" = xno; then
     # Some f90 compilers change the case of the mod file names. Others
     # require special flags to be set to deal with .f90 files. Find out
     # about this compiler.
     AC_LANG_PUSH(Fortran)
     AC_FC_SRCEXT(f90)
     AC_LANG_POP(Fortran)
     AX_F90_MODULE_FLAG
     AC_SUBST(MOD_FLAG, ${ax_cv_f90_modflag})
   fi

   # Set the FCLIBS flag to help with shared libraries.
#   AC_FC_LIBRARY_LDFLAGS
#   AC_F77_LIBRARY_LDFLAGS
fi

AM_CONDITIONAL(BUILD_F03, [test "x$nc_build_f03" = xyes])

# Set up libtool.
AC_MSG_NOTICE([setting up libtool])
LT_PREREQ([2.2])
LT_INIT

# Valgrind tests don't work with shared builds because of some libtool
# weirdness.
if test "x$enable_shared" = xyes; then
   if test $enable_valgrind_tests = yes; then
      AC_MSG_WARN([No valgrind tests with shared libraries])
      enable_valgrind_tests=no
   fi
fi

AC_MSG_NOTICE([finding other utilities])

# Find the install program.
AC_PROG_INSTALL

AC_MSG_NOTICE([trying to set fortran flags for this platform])
AC_MSG_CHECKING([what FC is set to])
AC_MSG_RESULT([$FC])
# Set fortran flag if the user has not already set it in CPPFLAGS.
case "$CPPFLAGS" in
 *pgiFortran*|*NAGf90Fortran*|*f2cFortran*|*hpuxFortran*|*apolloFortran*|*sunFortran*|*IBMR2Fortran*|*CRAYFortran*|*PATHSCALE_COMPILER*|*gFortran*|*mipsFortran*|*DECFortran*|*vmsFortran*|*CONVEXFortran*|*PowerStationFortran*|*AbsoftUNIXFortran*|*AbsoftProFortran*|*SXFortran*)
    ;;
   *)
      case "$FC" in
         pgf95|pgf90|pgf77|ifort*|lf90|lf95)
            AC_DEFINE(pgiFortran, [1], [Turned on by netCDF configure.])
            ;;
          *gfortran)
             if (gfortran --version | grep '(GCC 4.1.' || gfortran --version | grep '(GCC 4.0.'); then
                AC_DEFINE(gFortran, [1], [Turned on by netCDF configure.])
             else
                AC_DEFINE(pgiFortran, [1], [Turned on by netCDF configure.])
             fi
             ;;
          *g77)
             AC_DEFINE(f2cFortran, [1], [Turned on by netCDF configure.])
             ;;
          *)
             # If it's a gnu compiler, guess f2c.
             if test "x$ac_cv_fc_compiler_gnu" = xyes; then
                AC_DEFINE(f2cFortran, [1], [Turned on by netCDF configure.])
             fi
             ;;
      esac
esac

# Check to see if any macros must be set to enable large (>2GB) files.
AC_SYS_LARGEFILE

#AC_MSG_NOTICE([displaying some results])

## This next macro just prints some results for debugging
## support issues.
#UD_DISPLAY_RESULTS

# For nightly build testing, output CC, FC, etc.
#echo "CPPFLAGS=$CPPFLAGS CC=$CC CFLAGS=$CFLAGS CXX=$CXX CXXFLAGS=$CXXFLAGS FC=$FC FCFLAGS=$FCFLAGS F77=$F77 FFLAGS=$FFLAGS LDFLAGS=$LDFLAGS LIBS=$LIBS FLIBS=$FLIBS F90LIBS=$F90LIBS" >> comps.txt

ac_cv_prog_f90_uppercase_mod=no

# This will test the f77 compiler, as well as check some fortran types.
if test "$enable_fortran_type_check" = yes; then
   UD_FORTRAN_TYPES
   if test "x$ud_fc_failed" = xyes; then
      AC_MSG_ERROR([F77 compiler doesn't work, and compiler recovery disabled.])
   fi
else
   AC_DEFINE(NCBYTE_T, byte)
   AC_DEFINE(NCSHORT_T, integer*2)
   AC_DEFINE(NF_INT1_T, byte)
   AC_DEFINE(NF_INT2_T, integer*2)
   AC_DEFINE(NF_INT1_IS_C_SIGNED_CHAR, 1, [default])
   AC_DEFINE(NF_INT2_IS_C_SHORT, 1, [default])
   AC_DEFINE(NF_INT_IS_C_INT, 1, [default])
   AC_DEFINE(NF_REAL_IS_C_FLOAT, 1, [default])
   AC_DEFINE(NF_DOUBLEPRECISION_IS_C_DOUBLE, 1, [default])
fi

# Test whether F90 module names get capitalized.
if test "x$nc_build_f90" = xyes; then
   AC_PROG_FC_UPPERCASE_MOD
fi
if test "x$ac_cv_prog_f90_uppercase_mod" = xyes ; then
AC_MSG_FAILURE([Upper Case Modules names not supported])
fi
AM_CONDITIONAL(UPPER_CASE_MOD, [test "x$ac_cv_prog_f90_uppercase_mod" = xyes])

# Does the user want to run extra example tests
AC_MSG_CHECKING([whether extra example tests should be run])
AC_ARG_ENABLE([extra-example-tests],
              [AS_HELP_STRING([--enable-extra-example-tests],
                              [Run extra example tests; requires GNU sed. Ignored if \
                               netCDF-4 is not enabled.])])
test "x$enable_extra_example_tests" = xyes || enable_extra_example_tests=no
AC_MSG_RESULT($enable_extra_example_tests)
AM_CONDITIONAL(EXTRA_EXAMPLE_TESTS, [test "x$enable_extra_example_tests" = xyes])

# If --enable-dll is specified the DLL will be built. This only works
# on mingw.
AC_MSG_CHECKING([whether a win32 DLL is desired])
AC_ARG_ENABLE([dll],
              [AS_HELP_STRING([--enable-dll],
                              [build a win32 DLL (only works on mingw)])])
test "x$enable_dll" = xyes || enable_dll=no
AC_MSG_RESULT([$enable_dll])
AM_CONDITIONAL(BUILD_DLL, [test x$enable_dll = xyes])
if test "x$enable_dll" = xyes; then
   AC_DEFINE(DLL_NETCDF, 1, [set this only when building a DLL under MinGW])
fi

AC_MSG_NOTICE([checking types, headers, and functions])

# Find the netCDF header and library.
AC_CHECK_HEADERS([netcdf.h], [], [AC_MSG_ERROR([netcdf.h could not be found. Please set CPPFLAGS.])])
AC_SEARCH_LIBS([nc_open], [netcdf])

# See if various libraries are available
AC_CHECK_FUNCS([nc_def_opaque nccreate nc_set_log_level oc_open nc_use_parallel_enabled])
#test "x$ac_cv_func_oc_open" = xyes && nc_has_dap=yes || nc_has_dap=no
# Whether we build/test some functionality depends on what we found in
# the C library.
#nc_build_v2=$ac_cv_func_nccreate
#nc_build_v4=$ac_cv_func_nc_def_opaque

nc_build_v2=no
nc_build_v4=no
nc_has_logging=no
nc_has_dap=no
nc_has_parallel=no

test "x$ac_cv_func_nc_def_opaque" = xyes && nc_build_v4=yes
test "x$ac_cv_func_nccreate" = xyes && nc_build_v2=yes
test "x$ac_cv_func_nc_set_log_level" = xyes && nc_has_logging=yes
test "x$ac_cv_func_oc_open" = xyes && nc_has_dap=yes
test "x$ac_cv_func_nc_use_parallel_enabled" = xyes && nc_has_parallel=yes

#AC_CHECK_FUNC([nc_def_opaque],[nc_build_v4=yes])
#AC_CHECK_FUNC([nccreate],[nc_build_v2=yes])
#AC_CHECK_FUNC([nc_set_log_level],[nc_has_logging=yes])
#AC_CHECK_FUNC([oc_open],[nc_has_dap=yes])
#AC_CHECK_FUNC([nc_use_parallel_enabled],[nc_has_parallel=yes])

if test "x$nc_build_v4" = xyes; then
   AC_DEFINE([USE_NETCDF4], [1], [if true, build netCDF-4])
fi

AC_MSG_CHECKING([netCDF v2 API present])
AC_MSG_RESULT([$nc_build_v2])
AC_MSG_CHECKING([netCDF-4 present])
AC_MSG_RESULT([$nc_build_v4])

AM_CONDITIONAL([USE_NETCDF4], [test "x$nc_build_v4" = xyes])
AM_CONDITIONAL([BUILD_V2], [test "x$nc_build_v2" = xyes])
AM_CONDITIONAL([USE_LOGGING], [test "x$nc_has_logging" = xyes])
AM_CONDITIONAL([BUILD_DAP], [test "x$nc_has_dap" = xyes])
AM_CONDITIONAL([BUILD_PARALLEL], [test "x$nc_has_parallel" = xyes])
AC_CHECK_HEADER(stdlib.h, ,AC_DEFINE([NO_STDLIB_H], [], [no stdlib.h]))
AC_CHECK_HEADER(sys/types.h, ,AC_DEFINE([NO_SYS_TYPES_H], [], [no sys_types.h]))
AC_CHECK_HEADERS([sys/param.h])
#AC_CHECK_HEADERS([locale.h])
AC_CHECK_HEADERS([locale.h stdio.h stdarg.h errno.h ctype.h fcntl.h malloc.h stdlib.h string.h \
			   strings.h unistd.h sys/stat.h getopt.h sys/time.h sys/resource.h])
AC_HEADER_DIRENT
AC_HEADER_STDC
AC_HEADER_SYS_WAIT
AC_HEADER_TIME

AC_FUNC_VPRINTF

# Check for <stdbool.h> that conforms to C99 requirements
AC_HEADER_STDBOOL

# Check for these functions...
#AC_CHECK_FUNCS([strlcat strerror snprintf strchr strrchr mktemp strcat strcpy strdup \
#			strcasecmp strtod strtoll strtoull getrlimit gettimeofday fsync MPI_Comm_f2c])
AC_FUNC_ALLOCA
#AC_CHECK_DECLS([isnan, isinf, isfinite, signbit],,,[#include <math.h>])
#AC_STRUCT_ST_BLKSIZE
#UD_CHECK_IEEE
AC_TYPE_SIZE_T
AC_TYPE_OFF_T
AC_CHECK_TYPES([ssize_t, ptrdiff_t, uchar, longlong])
AC_C_CHAR_UNSIGNED
#AC_C_BIGENDIAN
AC_CHECK_SIZEOF(short)
AC_CHECK_SIZEOF(int)
AC_CHECK_SIZEOF(long)
AC_CHECK_SIZEOF(float)
AC_CHECK_SIZEOF(double)
AC_CHECK_SIZEOF(off_t)
AC_CHECK_SIZEOF(size_t)
AC_CHECK_SIZEOF(ptrdiff_t)
AC_COMPUTE_INT([longlong_is_long],[sizeof(long long) == sizeof(long)],[],[longlong_is_long=0])
if test $longlong_is_long = 1; then
   CFLAGS="$CFLAGS -DLONGLONG_IS_LONG"
fi

# Automake conditionals need to be called, whether the answer is yes
# or no.
AM_CONDITIONAL(TEST_PARALLEL, [test "x$enable_parallel_tests" = xyes])
AM_CONDITIONAL(USE_VALGRIND_TESTS, [test "x$enable_valgrind_tests" = xyes])

# Flags for nf-config script; by design $prefix, $includir, $libdir,
# etc.  are left as shell variables in the script so as to facilitate
# relocation
NC_LIBS="-lnetcdf"
NC_LIBS="$LDFLAGS $NC_LIBS $LIBS"

case "x$target_os" in
xsolaris*)
  NEWNCLIBS=""
  for x in $NC_LIBS ; do
    case "$x" in
    -L*) r=`echo "$x" | sed -e 's|^-L|-R|'`
	 NEWNCLIBS="$NEWNCLIBS $x $r"
	 ;;
    *)	 NEWNCLIBS="$NEWNCLIBS $x" ;;
    esac
  done
  NC_LIBS="$NEWNCLIBS"
  ;;
*);;
esac

NC_FLIBS="-lnetcdff $NC_LIBS"

AC_SUBST(NC_LIBS,[$NC_LIBS])
AC_SUBST(NC_FLIBS,[$NC_FLIBS])
AC_SUBST(HAS_F90,[$nc_build_f90])
AC_SUBST(HAS_F03,[$nc_build_f03])
AC_SUBST(HAS_DAP,[$nc_has_dap])
AC_SUBST(HAS_NC2,[$nc_build_v2])
AC_SUBST(HAS_NC4,[$nc_build_v4])

AC_MSG_NOTICE([generating header files and makefiles])

# We will output a fortran inc file, nfconfig.inc.
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_HEADERS([libsrc/nfconfig1.inc:libsrc/nfconfig.in],
   [sed '1d;s%^/\* \(.*\) \*/%\1%' libsrc/nfconfig1.inc >conftest.tmp && mv conftest.tmp libsrc/nfconfig.inc])

AC_CONFIG_FILES([Makefile
                 nf-config
                 netcdf-fortran.pc
                 docs/Makefile
                 libsrc/Makefile
                 fortran/Makefile
                 nf_test/Makefile
                 nf03_test/Makefile
                 CMakeExtras/Makefile
                 examples/Makefile
                 examples/F90/Makefile
                 examples/F77/Makefile],
		 [test -f  nf-config && chmod 755 nf-config ])

AC_OUTPUT()
