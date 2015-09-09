# SYNOPSIS
#
#   ACX_HDF5([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
# This macro looks for the HDF5 library
#
#
#   ACTION-IF-FOUND is a list of shell commands to run if the HDF5 library
#   is found, and ACTION-IF-NOT-FOUND is a list of commands to run it if it
#   is not found. If ACTION-IF-FOUND is not specified, the default action
#   will define HAVE_HDF5.
#
#   This macro will set the following variables
#   HDF5_LIBS HDF5_CPPFLAGS HDF5_LDFLAGS
#
# LICENSE
#
#   Copyright (c) 2009 Magnus Hagdorn <Magnus.Hagdorn@ed.ac.uk>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

AC_DEFUN([ACX_HDF5], [

acx_hdf5_ok=no

CPPFLAGSsave=$CPPFLAGS
LIBSsave=$LIBS
LDFLAGSsave=$LDFLAGS
HDF5_LIBS=""
HDF5_LDFLAGS=""
HDF5_CPPFLAGS=""

AC_ARG_WITH(hdf5,
  [AS_HELP_STRING([--with-hdf5],[root directory path where HDF5 is installed. (defaults to /usr/local or /usr if not found in /usr/local)])],
  [],
  [with_hdf5=yes])


# check whether hdf5 is disabled
AC_MSG_CHECKING([for hdf5])
AC_MSG_RESULT($with_hdf5)
AS_IF([test "x$with_hdf5" != xno],
      [ 
	hdf5_inc_path=""
	hdf5_lib_path=""
        # check if with_hdf5 is a path and if so setup various search paths
        if test -d "$with_hdf5"; then
          hdf5_inc_path="$with_hdf5"/include
          hdf5_lib_path="$with_hdf5"/lib
        fi
        # check whether we should use a non standard include path
        AC_ARG_WITH(hdf5-include,
          [AS_HELP_STRING([--with-hdf5-include],[path to where hdf5 header files can be found])],
          [
           if test -d "$withval"; then
               hdf5_inc_path=$withval
           else
               AC_MSG_ERROR([Cannot find directory "$withval"])
           fi
          ])

        # check whether we should use a non standard library path
        AC_ARG_WITH(hdf5-lib,
          [AS_HELP_STRING([--with-hdf5-lib],[path to where hdf5 library files can be found])],
          [
           if test -d "$withval"; then
               hdf5_lib_path="$withval"
           else
               AC_MSG_ERROR([Cannot find directory "$withval"])
           fi
          ])
        
        if test x"$hdf5_inc_path"x != xx ; then
           HDF5_CPPFLAGS=-I"$hdf5_inc_path"
        fi
        if test x"$hdf5_lib_path"x != xx ; then
           HDF5_LDFLAGS=-L"$hdf5_lib_path"
        fi

	# check for libraries
	LDFLAGS="$LDFLAGS $HDF5_LDFLAGS"
        # we always need to link to the C libraries, so let's look for them
        AC_LANG_PUSH([C])
        CPPFLAGS="$CPPFLAGS $HDF5_CPPFLAGS"
	AC_CHECK_LIB(hdf5,H5open,[acx_hdf5_ok=yes; HDF5_LIBS="-lhdf5"],[acx_hdf5_ok=no])
	AC_CHECK_LIB(hdf5_hl,H5Dclose,[acx_hdf5_ok=yes; HDF5_LIBS="$HDF5_LIBS -lhdf5_hl"],[acx_hdf5_ok=no],[$HDF5_LIBS])
        AC_LANG_POP([C])

        # figure out how to use HDF5 from various languages
        AC_LANG_CASE(
        [C],[
	AC_CHECK_HEADER([hdf5.h],[acx_hdf5_ok=yes],[acx_hdf5_ok=no])
	AC_CHECK_HEADER([hdf5_hl.h],[acx_hdf5_ok=yes],[acx_hdf5_ok=no],[AC_INCLUDES_DEFAULT 
                                                                        #include <hdf5.h>])
        ],
        [C++],[AC_MSG_NOTICE([C++ not checked for yet])],
        [Fortran 77],[AC_MSG_NOTICE([F77 not checked for yet])],
        [Fortran],[AC_MSG_NOTICE([Fortran not checked for yet])])
      ])

AC_SUBST(HDF5_LIBS)
AC_SUBST(HDF5_LDFLAGS)
AC_SUBST(HDF5_CPPFLAGS)

# reset variables to original values
CPPFLAGS=$CPPFLAGSsave
LIBS=$LIBSsave
LDFLAGS=$LDFLAGSsave
  
# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_hdf5_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_HDF5,1,[Define if you have the HDF5 library.]),[$1])
        :
else
        acx_hdf5_ok=no
        $2
fi
])dnl ACX_HDF5
