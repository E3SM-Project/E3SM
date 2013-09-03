dnl Process this file with autoconf to produce a configure script.
AC_PREREQ(2.59)
AC_INIT(HOMME, 2.0.0, jgarcia@ucar.edu)
AC_CONFIG_AUX_DIR([./helpers]) 
AC_CONFIG_SRCDIR([src/main.F90])
AM_INIT_AUTOMAKE

dnl Check for *FLAGS environments
dnl AX_COMPILER_MAXOPT_FLAGS

dnl Checks for programs.
AC_PROG_F77
AC_PROG_FC
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_RANLIB

AC_ARG_ENABLE([blas],
            [AS_HELP_STRING([--enable-blas],
              [builds a version of BLAS that comes with this source distribution])],
            [enable_blas=yes],
            [enable_blas=no])


AC_ARG_ENABLE([lapack],
            [AS_HELP_STRING([--enable-lapack],
              [builds a version of LAPACK that comes with this source distribution])],
            [enable_lapack=yes],
            [enable_lapack=no])


# defines if we want blas in this build
AS_IF([test "x$enable_blas" = "xyes"], [
  dnl Do the stuff needed for enabling the feature
            echo 'Local BLAS enabled'
            MAYBE_BLAS=libs/blas
],[MAYBE_BLAS=])

# defines if we want lapack in this build
AS_IF([test "x$enable_lapack" = "xyes"], [
  dnl Do the stuff needed for enabling the feature
            echo 'Local LAPACK enabled'
            MAYBE_LAPACK=libs/lapack
],[MAYBE_LAPACK=])


AC_SUBST([MAYBE_BLAS])
AC_SUBST([MAYBE_LAPACK])

AC_CONFIG_SUBDIRS([utils/pio libs/blas libs/lapack build/preqx build/sweqx])


AC_CONFIG_FILES([Makefile])

AC_OUTPUT
