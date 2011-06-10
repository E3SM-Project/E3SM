# SYNOPSIS
#
#   ACX_NETCDF([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
# This macro looks for the NETCDF library
#
#
#   ACTION-IF-FOUND is a list of shell commands to run if the NETCDF library
#   is found, and ACTION-IF-NOT-FOUND is a list of commands to run it if it
#   is not found. If ACTION-IF-FOUND is not specified, the default action
#   will define HAVE_NETCDF. If the netCDF 4 API is found the symbol
#   HAVE_NETCDF4 is also defined.
#
#   This macro will set the following variables
#   NETCDF_LIBS NETCDF_CPPFLAGS NETCDF_FCFLAGS NETCDF_LDFLAGS
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

AC_DEFUN([ACX_NETCDF], [
AC_REQUIRE([ACX_HDF5])

acx_netcdf_ok=no

CPPFLAGSsave=$CPPFLAGS
FCFLAGSsave=$FCFLAGS
LIBSsave=$LIBS
LDFLAGSsave=$LDFLAGS
# initialise outputs
NETCDF_LIBS=""
NETCDF_LDFLAGS=""
NETCDF_CPPFLAGS=""
NETCDF_FCFLAGS=""

AC_ARG_WITH(netcdf,
  [AS_HELP_STRING([--with-netcdf],[root directory path where NETCDF is installed. (defaults to /usr/local or /usr if not found in /usr/local)])],
  [],
  [with_netcdf=yes])

# Check for variable NETCDF_PATH in environment
if test -d "$NETCDF_PATH"; then
  with_netcdf=$NETCDF_PATH
fi


# check whether netcdf is disabled
AC_MSG_CHECKING([for netcdf])
AC_MSG_RESULT($with_netcdf)
AS_IF([test "x$with_netcdf" != xno],
      [ 
	netcdf_inc_path=""
	netcdf_lib_path=""
        # check if with_netcdf is a path and if so setup various search paths
        if test -d "$with_netcdf"; then
          netcdf_inc_path="$with_netcdf"/include
          netcdf_lib_path="$with_netcdf"/lib
        fi
        # check whether we should use a non standard include path
        AC_ARG_WITH(netcdf-include,
          [AS_HELP_STRING([--with-netcdf-include],[path to where netcdf header files can be found])],
          [
           if test -d "$withval"; then
               netcdf_inc_path=$withval
           else
               AC_MSG_ERROR([Cannot find directory "$withval"])
           fi
          ])

        # check whether we should use a non standard library path
        AC_ARG_WITH(netcdf-lib,
          [AS_HELP_STRING([--with-netcdf-lib],[path to where netcdf library files can be found])],
          [
           if test -d "$withval"; then
               netcdf_lib_path="$withval"
           else
               AC_MSG_ERROR([Cannot find directory "$withval"])
           fi
          ])
        
        if test x"$netcdf_inc_path"x != xx ; then
           NETCDF_CPPFLAGS=-I"$netcdf_inc_path"
        fi
        if test x"$netcdf_lib_path"x != xx ; then
           NETCDF_LDFLAGS=-L"$netcdf_lib_path"
        fi

	# check for libraries
	LDFLAGS="$LDFLAGS $NETCDF_LDFLAGS $HDF5_LDFLAGS"
        LIBS="$HDF5_LIBS"
        # we always need to link to the C libraries, so let's look for them
        AC_LANG_PUSH(C)
        CPPFLAGS="$CPPFLAGS $NETCDF_CPPFLAGS"
        AC_SEARCH_LIBS(nc_inq_libvers,netcdf,[acx_netcdf_ok=yes],[acx_netcdf_ok=no;AC_MSG_ERROR(cannot find netCDF C library)])
        AC_LANG_POP([C])

        # figure out how to use netcdf from various languages
        AC_LANG_CASE(
        [C],[
	# check for header files
	AC_CHECK_HEADER([netcdf.h],[acx_netcdf_ok=yes],[acx_netcdf_ok=no])
        ],
        [C++],[
	# check for C++ header files
	AC_CHECK_HEADER([netcdfcpp.h],[acx_netcdf_ok=yes],[acx_netcdf_ok=no])
	LIBS="-lnetcdf_c++ $LIBS"
	AC_MSG_CHECKING([for netCDF C++ library])
	AC_LINK_IFELSE([
                AC_LANG_PROGRAM(
                    [[
@%:@include <netcdfcpp.h>
                    ]],
                    [[
NcError err_handler;
                    ]]
                )],
                [
                acx_netcdf_ok=yes
                AC_MSG_RESULT([yes])
                ],
                [
                acx_netcdf_ok=no
                AC_MSG_RESULT([no])                ]
            )
	],
        [Fortran 77],[
        AC_SEARCH_LIBS(NF_INQ_LIBVERS,netcdff netcdf,[acx_netcdf_ok=yes],[acx_netcdf_ok=no;AC_MSG_ERROR(cannot find netCDF fortran library)])
        ],
        [Fortran],[
        AC_SEARCH_LIBS(NF_INQ_LIBVERS,netcdff netcdf,[acx_netcdf_ok=yes],[acx_netcdf_ok=no; AC_MSG_ERROR(cannot find netCDF fortran library)])
        AC_REQUIRE([AC_FC_MODULE_FLAG])
        if test x"$netcdf_inc_path"x != xx ; then
           NETCDF_FCFLAGS="$ac_cv_fc_module_flag$netcdf_inc_path"
        else
           NETCDF_FCFLAGS="$ac_cv_fc_module_flag/usr/include"
        fi
	FCFLAGS="$FCFLAGS $NETCDF_FCFLAGS"
        AC_MSG_CHECKING([for f90 netCDF interface])
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[use netcdf])],[acx_netcdf_ok=yes; AC_MSG_RESULT([yes])],
                                                             [acx_netcdf_ok=no; AC_MSG_RESULT([no])])
        ])

	# check if netCDF4 API is available
        acx_netcdf4_ok=no
	if test x"$acx_netcdf_ok" = xyes; then
           AC_MSG_CHECKING([for netcdf 4 API])
           AC_MSG_RESULT()
           AC_LANG_CASE(
           [C],[ AC_CHECK_FUNC(nc_inq_grps,[acx_netcdf4_ok=yes;AC_MSG_RESULT([yes])],[AC_MSG_RESULT([no])]) ], 
           [C++],[AC_MSG_NOTICE([C++ not checked for yet])], 
           [Fortran 77],[AC_CHECK_FUNC(nf_inq_grps,[acx_netcdf4_ok=yes;AC_MSG_RESULT([yes])],[AC_MSG_RESULT([no])]) ],
           [Fortran],[
            AC_MSG_CHECKING([for nf90_inq_grps])
            AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[use netcdf;i=nf90_inq_grps(id1,id2,id3)])],[acx_netcdf4_ok=yes; AC_MSG_RESULT([yes])],
                                                                 [acx_netcdf4_ok=no; AC_MSG_RESULT([no])])
           ])
        fi
        if test x"$acx_netcdf4_ok" = xyes; then
          AC_DEFINE(HAVE_NETCDF4,1,[Define if the netCDF library is compiled with netCDF 4 API.])
        fi
      ])
NETCDF_LIBS=$LIBS
AC_SUBST(NETCDF_LIBS)
AC_SUBST(NETCDF_LDFLAGS)
AC_SUBST(NETCDF_CPPFLAGS)
AC_SUBST(NETCDF_FCFLAGS)

# reset variables to original values
CPPFLAGS=$CPPFLAGSsave
LIBS=$LIBSsave
LDFLAGS=$LDFLAGSsave
FCFLAGS=$FCFLAGSsave

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_netcdf_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_NETCDF,1,[Define if you have the NETCDF library.]),[$1])
        :
else
        acx_netcdf_ok=no
        $2
fi
])dnl ACX_NETCDF
