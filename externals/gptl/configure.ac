dnl Process this file with autoconf to produce a configure script.
AC_INIT(private.h)
AC_LANG(C)
echo "This script provides suggestions for settings to apply in macros.make"
echo "You can pass things like FC=gfortran or CC=pathcc to it to override defaults."
AC_PROG_MAKE_SET
AC_PROG_CC()
echo "Assuming C compiler is $CC"
AC_PROG_FC()
echo "Assuming Fortran compiler is $FC"
dnl This barfs when CC is pgcc and F77 pgf90
AC_FC_FUNC(z_zz,z_zz)
case $z_zz in
    z_zz__)     echo "Fortran name mangling: UNDERSCORING=-DFORTRANDOUBLEUNDERSCORE" ;;
    z_zz_)      echo "Fortran name mangling: UNDERSCORING=-DFORTRANUNDERSCORE" ;;
    Z_ZZ)       echo "Fortran name mangling: UNDERSCORING=-DFORTRANCAPS" ;;
    *)          echo "Fortran name mangling: UNDERSCORING=" ;;
esac

dnl For proper handling of const and inline
AC_C_INLINE
echo "Inlining: -Dinline=$ac_cv_c_inline"
dnl AC_FUNC_VPRINTF
AC_CHECK_LIB(papi,PAPI_library_init,[echo "PAPI library found: OK to set HAVE_PAPI=yes"],
                                    [echo "PAPI library not found: HAVE_PAPI=no"])

AC_CHECK_FUNC(backtrace_symbols,[echo "backtrace_symbols found: OK to set HAVE_BACKTRACE=yes"],
				[echo "backtrace_symbols NOT found: HAVE_BACKTRACE=no"])

unset usempich;
unset usempi;
AC_CHECK_LIB(mpich,MPI_Init,[echo "libmpich.a found: OK to set HAVE_MPI=yes";usempich=yes],
	     AC_CHECK_LIB(mpi,MPI_Init,[echo "libmpi.a found: OK to set HAVE_MPI=yes OK";usempi=yes],[echo "MPI library not found: HAVE_MPI=no"]))
AC_CHECK_FILE(/proc,[echo "/proc found: HAVE_SLASHPROC=yes"],
		    [echo "/proc not found: HAVE_SLASHPROC=no"])
AC_CHECK_LIB(pthread,pthread_mutex_init,[echo "pthreads library found: OK to set PTHREADS=yes"],
					[echo "pthreads library not found: PTHREADS=no"])

if test -n "${usempich}" ; then
AC_CHECK_LIB(mpich,iargc,[echo "iargc found in libmpich.a: OK to set HAVE_IARGCGETARG=yes"],	                   [echo "iargc not found in libmpich.a: HAVE_IARGCGETARG=no"])
fi

if test -n "${usempi}" ; then
AC_CHECK_LIB(mpi,iargc,[echo "iargc found in libmpi.a: OK to set HAVE_IARGCGETARG=yes"],
		       [echo "iargc not found in libmpi.a: HAVE_IARGCGETARG=no"])
fi

AC_CHECK_SIZEOF([void *])
if test "$ac_cv_sizeof_void_p" = 8; then
  echo "Pointer size = 8 so BIT64=yes";
elif test "$ac_cv_sizeof_void_p" = 4; then
  echo "Pointer size = 4 so BIT64=no";
fi
