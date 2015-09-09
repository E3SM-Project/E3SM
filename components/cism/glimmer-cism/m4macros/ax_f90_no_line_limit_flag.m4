# SYNOPSIS
#
# AX_F90_NO_LINE_LIMIT_FLAG
#
# DESCRIPTION
#
#   Find Fortran compiler flag which disables line length limitations.
#   Standard free-form fortran source code is limited to 132 characters per
#   line.
#
# LICENSE
#
#   Copyright (c) 2010 Magnus Hagdorn <Magnus.Hagdorn@ed.ac.uk>
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

AC_DEFUN([AX_F90_NO_LINE_LIMIT_FLAG],[
AC_CACHE_CHECK([fortran 90 flag for removing limit of number of columns in free form source files],
ax_cv_f90_ff_nolimit,
[AC_LANG_PUSH(Fortran)
ax_cv_f90_ff_nolimit="not found"
for ax_flag in "" "-ffree-line-length-none" "-Mextend"; do
  if test "$ax_cv_f90_ff_nolimit" = "not found" ; then
    ax_save_FCFLAGS="$FCFLAGS"
    FCFLAGS="$ax_save_FCFLAGS ${ax_flag}"
    AC_COMPILE_IFELSE([
program conftest_program
write(*,*) "0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789"
write(*,*) "0         1         2         3         4         5         6         7         8         9        10        11        12        13        14        15         "
end program conftest_program
],[ax_cv_f90_ff_nolimit="$ax_flag"],[])
    FCFLAGS="$ax_save_FCFLAGS"
  fi
done
if test "$ax_cv_f90_ff_nolimit" = "not found" ; then
  AC_MSG_ERROR([unable to find compiler flag for disabling column limit])
fi
AC_LANG_POP(Fortran)
])])

