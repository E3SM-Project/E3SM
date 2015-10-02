dnl Process this file with autoconf to produce a configure script.

dnl set version
m4_define([glimmer_major_version], [1])
m4_define([glimmer_minor_version], [7])
m4_define([glimmer_micro_version], [1])
m4_define([glimmer_version],
	[glimmer_major_version.glimmer_minor_version.glimmer_micro_version])

AC_INIT([GLIMMER-CISM],[glimmer_version],[glimmer-cism-devel@lists.berlios.de])
AC_CONFIG_HEADERS(config.inc)
AC_CONFIG_MACRO_DIR([m4macros])
AM_INIT_AUTOMAKE

GLIMMER_MAJOR_VERSION=glimmer_major_version
GLIMMER_MINOR_VERSION=glimmer_minor_version
GLIMMER_MICRO_VERSION=glimmer_micro_version
GLIMMER_VERSION=glimmer_version
AC_SUBST(GLIMMER_MAJOR_VERSION)
AC_SUBST(GLIMMER_MINOR_VERSION)
AC_SUBST(GLIMMER_MICRO_VERSION)
AC_SUBST(GLIMMER_VERSION)

AC_CONFIG_SRCDIR([glimmer-cism-config.in])

# save compiler flags, autoconf will change them later
glimmer_FCFLAGS=$FCFLAGS
glimmer_FFLAGS=$FFLAGS
AC_DEFINE_UNQUOTED(GLIMMER_FCFLAGS,"${glimmer_FCFLAGS}",[fortran compiler flags the model is configured with])

if test x"$F77"x == xx ; then
  export F77=$FC
fi
# we are building fortran libraries so, let's make them static
AC_PROG_LIBTOOL

AC_PROG_CC
AC_PROG_CXX
AC_PROG_FC
AC_PROG_F77
AC_PROG_EGREP


dnl Checking for python
build_python=1
AC_MSG_CHECKING(if we should build python bindings)
AC_ARG_ENABLE(python, AS_HELP_STRING([--disable-python], [Do not build Python modules]),
 [case "${enableval}" in
   yes) build_python=1  
        AC_MSG_RESULT([yes]);;
   no)  build_python=0
        AC_MSG_RESULT([no])  ;;
   *) AC_MSG_ERROR(bad value ${enableval} for --disable-python) ;;
 esac],[AC_MSG_RESULT([yes])])
if test "$build_python" == 1; then
  AC_PYTHON_DEVEL([>= '2.3'])
  proj_ok=no
  AC_LANG_PUSH([C])
  AC_CHECK_HEADER(projects.h,proj_ok=yes,proj_ok=no)
  if test x"$proj_ok"x != xyesx; then
    build_python=0
  fi	
  AC_SEARCH_LIBS(pj_fwd,proj,proj_ok=yes,proj_ok=no)
  if test x"$proj_ok"x != xyesx; then
    build_python=0
  fi
  if test "$build_python" == "0"; then
    AC_MSG_WARN(libproj unusable)
  fi
  AC_LANG_POP([C])
fi
if test "$build_python" == "1"; then
  AM_PATH_PYTHON([2.3],,[AC_MSG_RESULT(need at least python 2.3)])
  AX_PYTHON_MODULE(numpy,[build_python=0])
else
  AC_MSG_WARN(not building python extension)
fi
AM_CONDITIONAL([BUILD_PYTHON], [test "$build_python" == "1"])




FCFLAGS=$glimmer_FCFLAGS
FFLAGS=$glimmer_FFLAGS
AC_LANG(Fortran)
AC_FC_SRCEXT(f90)
# Define FCLIBS for fortran libraries needed for CXX linking
AC_FC_LIBRARY_LDFLAGS
AC_FC_PP_DEFINE

AC_FC_MODULE_EXTENSION
# Define FLIBS for fortran libraries needed for CXX linking 
AC_F77_LIBRARY_LDFLAGS
# get the correct flag to set module include path
#AX_F90_MODULE_FLAG
#F90_MODFLAG=$ax_cv_f90_modflag
#AC_SUBST(F90_MODFLAG)
AC_FC_MODULE_FLAG

#AX_F90_NO_LINE_LIMIT_FLAG
#FCFLAGS="$FCFLAGS $ax_cv_f90_ff_nolimit"

# checking for f2003 command line parser
# 
AC_FC_SRCEXT(f90)
AC_MSG_CHECKING([for f2003 command line parser])
AC_LINK_IFELSE([AC_LANG_PROGRAM([],
               [i=command_argument_count()
                call get_command_argument(i)  ])],
               [AC_DEFINE(HAVE_2003ARGS,,[defined if you have fortran 2003 intrinsics for parsing command line])
                AC_MSG_RESULT(yes)], [AC_MSG_RESULT(no)])

# which Fortran compiler are we running?

AC_LANG_PUSH(Fortran)
AX_FC_VERSION()

f95string="$ac_fc_version_output"
f95_preprocessor="-fpp"

fs=`echo $f95string |$EGREP "XL Fortran"`
if test "$fs" != "" ; then
  f95v="IBM xlf"
  CPPFLAGS="$CPPFLAGS -DCPRIBM"
fi
# testing for SUN f95 compiler
fs=`echo $f95string |$EGREP "Sun|Forte"`
if test "$fs" != "" ; then
  f95v="SUN f95"
  CPPFLAGS="$CPPFLAGS -DsunFortran"
fi
# testing for NAG f95 compiler
fs=`echo $f95string |$EGREP "NAGWare"`
if test "$fs" != "" ; then
  f95v="NAGWare f95"
  CPPFLAGS="$CPPFLAGS -DNAGf90Fortran"
fi
# testing for INTEL ifort compiler
fs=`echo $f95string |$EGREP "Intel"`
if test "$fs" != "" ; then
  f95v="Intel ifort"
  CPPFLAGS="$CPPFLAGS -DpgiFortran"
fi
# testing for gfortran
fs=`echo $f95string |$EGREP "gcc"`
if test "$fs" != "" ; then
  f95v="gfortran"
  f95_preprocessor="-x f95-cpp-input"
  CPPFLAGS="$CPPFLAGS -DpgiFortran"
fi
# testing for Portland group
fs=`echo $f95string |$EGREP "Portland"`
if test "$fs" != "" ; then
  f95v="Portland pgf90"
  f95_preprocessor="-Mpreprocess"
  CPPFLAGS="$CPPFLAGS -DpgiFortran"
fi


AC_MSG_RESULT([$f95v])
AC_FC_WRAPPERS
# checking for BLAS library
ac_gl_have_blas=0
ACX_BLAS([ac_gl_have_blas=1])
AM_CONDITIONAL(EXT_BLAS,test "$ac_gl_have_blas" == 1)

# checking for netCDF
AC_LANG_PUSH([C++])
ac_gl_have_netcdf_cxx=0
ACX_NETCDF([ac_gl_have_netcdf_cxx=1],[])
NETCDF_CXX_LIBS=$NETCDF_LIBS
AC_SUBST(NETCDF_CXX_LIBS)
if test "$ac_gl_have_netcdf_cxx" != 1; then
  AC_MSG_WARN([Cannot build netCDF compare utility because C++ bindings to netCDF are not available])
fi
AM_CONDITIONAL(BUILD_COMPARE, test "$ac_gl_have_netcdf_cxx" == 1)
AC_LANG_POP([C++])
ACX_NETCDF([],[AC_MSG_ERROR(cannot find netCDF fortran library)])

# checking for ENMABAL
AC_ARG_WITH(enmabal,AS_HELP_STRING([--with-enmabal],[location ENMABAL library]),[ac_cv_enmabal_prefix=$withval])
if test "$ac_cv_enmabal_prefix" != ""; then
  LIBS="$LIBS -lsmb"
  AC_MSG_RESULT(enmabal location: $ac_cv_enmabal_prefix)
  AC_DEFINE(USE_ENMABAL,1,[Define to 1 if you want to use ENMABAL])
  FCFLAGS="$FCFLAGS $ax_cv_f90_modflag$ac_cv_enmabal_prefix"
  LDFLAGS="$LDFLAGS -L$ac_cv_enmabal_prefix"
  AC_MSG_CHECKING([for f90 ENMABAL interface])
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[use smb_mecons])],AC_MSG_RESULT([yes]),AC_MSG_FAILURE([no]))
fi

# building the documentation
have_latex=1
AC_CHECK_PROGS(LATEX,[latex])
if test x"$LATEX"x == xx; then
  have_latex=0
fi
AC_CHECK_PROGS(BIBTEX,[bibtex])
if test x"$BIBTEX"x == xx; then
  have_latex=0
fi
AC_CHECK_PROGS(DVIPS,[dvips]) 
if test x"$DVIPS"x == xx; then
  have_latex=0
fi
AC_CHECK_PROGS(GNUPLOT,[gnuplot]) 
if test x"$GNUPLOT"x == xx; then
  have_latex=0
fi
AC_CHECK_PROGS(DIA,[dia]) 
if test x"$DIA"x == xx; then
  have_latex=0
fi
AC_CHECK_PROGS(FIG2DEV,[fig2dev]) 
if test x"$FIG2DEV"x == xx; then
  have_latex=0
fi
AC_CHECK_PROGS(PS2EPSI,[ps2epsi]) 
if test x"$PS2EPSI"x == xx; then
  have_latex=0
fi
AC_CHECK_PROGS(INKSCAPE, [inkscape])
if test x"$INKSCAPE"x == xx; then
  have_latex=0
fi
AC_CHECK_PROGS(SED,[sed])
if test x"$SED"x == xx; then
  have_latex=0
fi
AC_CHECK_PROGS(PDF2PS,[pdf2ps])
if test x"PDF2PS"x == xx; then
  have_latex=0
fi
AC_CHECK_PROGS(CONVERT,[convert])
if test x"convert"x == xx; then
  have_latex=0
fi
AC_CHECK_PROGS(DVIPDFT,[dvipdft dvipdfm])

# build postscript documentation
do_doc_ps=0
AC_MSG_CHECKING([whether postscript documentation should be built])
AC_ARG_ENABLE(doc-ps, AS_HELP_STRING([--enable-doc-ps],[build postscript documentation]),
[case "${enableval}" in
 yes)  if test $have_latex == 0; then
          AC_MSG_RESULT([no])
          AC_MSG_WARN([did not find all the required tools])
       else
          AC_MSG_RESULT([yes])
          do_doc_ps=1
       fi;;
 no)   AC_MSG_RESULT([no]);;
 *) AC_MSG_ERROR(bad value ${enableval} for --enable-doc-ps) ;;
esac],[AC_MSG_RESULT([no])])
AM_CONDITIONAL([DO_DOC_PS],  [test $do_doc_ps == 1])

# build pdf documentation
do_doc_pdf=0
AC_MSG_CHECKING([whether PDF documentation should be built])
AC_ARG_ENABLE(doc-pdf, AS_HELP_STRING([--enable-doc-pdf],[build PDF documentation]),
[case "${enableval}" in
 yes)  if test $have_latex == 0 || test x"$DVIPDFT"x == xx; then
          AC_MSG_RESULT([no])
          AC_MSG_WARN([did not find all the required tools])
       else
          AC_MSG_RESULT([yes])
          do_doc_pdf=1
       fi;;
 no)   AC_MSG_RESULT([no]);;
 *) AC_MSG_ERROR(bad value ${enableval} for --enable-pdf) ;;
esac],[AC_MSG_RESULT([no])])
AM_CONDITIONAL([DO_DOC_PDF],  [test $do_doc_pdf == 1])

# doxygen documentation
DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(OFF)
DX_PS_FEATURE(OFF)
DX_INIT_DOXYGEN([Glimmer-CISM],glimmer-cism.cfg)

# some options
AC_MSG_CHECKING([whether profilng is enabled])
AC_ARG_ENABLE(profile, AS_HELP_STRING([--enable-profile],[Turn on profiling]),
[case "${enableval}" in
 yes)  AC_DEFINE(PROFILING,1,[Define to 1 if you want to enable extra profiling])
       FCFLAGS="$FCFLAGS -pg"
       AC_MSG_RESULT([yes]);;	
 no)   AC_MSG_RESULT([no]);;
 *) AC_MSG_ERROR(bad value ${enableval} for --enable-profile) ;;
esac],[AC_MSG_RESULT([no])])

AC_MSG_CHECKING([whether command-line arguments are enabled])
AC_ARG_ENABLE(command-line, AS_HELP_STRING([--enable-command-line],[Turn on command-line arguments]),
[case "${enableval}" in
 yes)  AC_DEFINE(COMMAND_LINE,1,[Define to 1 if you want to enable command-line arguments])
       AC_MSG_RESULT([yes]);;	
 no)   AC_MSG_RESULT([no]);;
 *) AC_MSG_ERROR(bad value ${enableval} for --enable-command-line) ;;
esac],[AC_MSG_RESULT([no])])

#AC_MSG_CHECKING([whether extended restarts are enabled])
#AC_ARG_ENABLE(restarts, AS_HELP_STRING([--enable-restarts],[Build extended restart code]),
#[case "${enableval}" in
# yes)  AC_DEFINE(RESTARTS,1,[Define to 1 if you want to enable extended restarts])
#       AC_MSG_RESULT([yes]);;	
# no)   AC_MSG_RESULT([no]);;
# *) AC_MSG_ERROR(bad value ${enableval} for --enable-restarts) ;;
#esac],[AC_MSG_RESULT([no])])
#AM_CONDITIONAL([BUILD_RESTARTS],[test "$enable_restarts" == "yes"])

AC_MSG_CHECKING([whether the build is single precision])
AC_ARG_ENABLE(sp, AS_HELP_STRING([--enable-sp],[Build code as single precision]),
[case "${enableval}" in
 yes)  AC_DEFINE(GLIMMER_SP,1,[Define to 1 if you want to enable extended restarts])
       AC_MSG_RESULT([yes]);;	
 no)   AC_MSG_RESULT([no]);;
 *) AC_MSG_ERROR(bad value ${enableval} for --enable-sp) ;;
esac],[AC_MSG_RESULT([no])])

AC_MSG_CHECKING([whether the build is parallel MPI])
AC_ARG_ENABLE(mpi, AS_HELP_STRING([--enable-mpi],[Build code as parallel MPI]),
[case "${enableval}" in
 yes)  AC_MSG_RESULT([yes]);;    
 no)   AC_MSG_RESULT([no]);;
 *) AC_MSG_ERROR(bad value ${enableval} for --enable-mpi) ;;
esac],[AC_MSG_RESULT([no])])
AM_CONDITIONAL(DISTRIBUTED, test "$enable_mpi" == "yes")

# should we run the tests
#AC_MSG_CHECKING([whether we should run the tests])
#AC_ARG_ENABLE(tests, AS_HELP_STRING([--enable-tests],[Run short tests]),
#[case "${enableval}" in
# yes)  AC_MSG_RESULT([yes]);;	
# no)   AC_MSG_RESULT([no]);;
# *) AC_MSG_ERROR(bad value ${enableval} for --enable-tests) ;;
#esac],[AC_MSG_RESULT([no])])
#AM_CONDITIONAL([RUN_TESTS], [test "$enable_tests" == "yes"])

#specify the directory where the gold standard files are located
AC_MSG_CHECKING([whether location of gold standard files are specified])
AC_ARG_WITH(gold_std,AS_HELP_STRING([--with-gold-std],[location of gold_std]),
              [GOLD_STD_PREFIX=$withval])
AC_MSG_RESULT(gold_std location: $GOLD_STD_PREFIX)
AC_SUBST(GOLD_STD_PREFIX)

dnl checking if model runs should be submitted to SGE
m4_pattern_allow([^AM_SGE_OPTS$])
AM_SGE_OPTS=""
AC_MSG_CHECKING([whether model runs should be submitted to SGE])
AC_ARG_WITH(sge-jobs,AS_HELP_STRING([--with-sge-jobs],[use Sun Grid Engine to run models]),[
  AM_SGE_OPTS="-s"
  SGE_OPTS=""
  if test x"$with_sge_jobs"x != xyesx; then
    AM_SGE_OPTS="$AM_SGE_OPTS -o"
    SGE_OPTS="\"$with_sge_jobs\""
  fi 
  AC_MSG_RESULT([$withval])],[AC_MSG_RESULT([no])])
AC_SUBST(AM_SGE_OPTS)
AC_SUBST(SGE_OPTS)

dnl checking which extensions should be built
# erosion module
build_extension_erosion="false"
AC_MSG_CHECKING([whether extension: erosion should be built])
AC_ARG_ENABLE(extension-erosion,AS_HELP_STRING([--enable-extension-erosion],[Build Glimmer-CISM with extension erosion]),
[case "${enableval}" in
 yes)  build_extension_erosion="true"
       AC_MSG_RESULT([yes]);;	
 no)   AC_MSG_RESULT([no]);;
 *) AC_MSG_ERROR(bad value ${enableval} for --enable-extension-erosion) ;;
esac],[AC_MSG_RESULT([no])])
AM_CONDITIONAL([BUILD_EXT_EROSION],[test "$build_extension_erosion" == "true"])
AC_SUBST(build_extension_erosion)

# checking for PHAML
build_extension_libphaml="false"
AC_MSG_CHECKING([whether extension: libphaml should be built])
PHAML_PREFIX=""
AC_ARG_WITH(phaml,AS_HELP_STRING([--with-phaml],[location of PHAML library]),[PHAML_PREFIX=$withval])
PHAML_GRAPHICS_PREFIX=""
AC_ARG_WITH(phaml-graphics,AS_HELP_STRING([--with-phaml-graphics],[location F90GL library]),[PHAML_GRAPHICS_PREFIX=$withval])
if test "$PHAML_PREFIX" != ""; then
  AC_MSG_RESULT([yes])
  ACX_MPI()
  LIBS_save=$LIBS
  LDFLAGS_save=$LDFLAGS
  FCFLAGS_save=$FCFLAGS
  LIBS="$LIBS -lphaml"
  LDFLAGS="$LDFLAGS -L$PHAML_PREFIX/lib"
  FCFLAGS="$FCFLAGS -I$PHAML_PREFIX/modules -I$PHAML_PREFIX/src -L$PHAML_PREFIX/lib"
  AC_MSG_RESULT(phaml location: $PHAML_PREFIX)
  AC_MSG_CHECKING([for PHAML interface])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],[use phaml])],AC_MSG_RESULT([yes]),AC_MSG_FAILURE([Can not link to PHAML Library]))

  AC_MSG_CHECKING([whether libphaml extension should be built with graphics support])
  if test "$PHAML_GRAPHICS_PREFIX" != ""; then
    AC_MSG_RESULT([yes])
    LIBS="$LIBS -L$PHAML_GRAPHICS_PREFIX/lib -L/usr/X11R6/lib -lXaw -lXt -lXmu -lXi -lX11 -lXext -lm -lf90glut -lf90GLU -lf90GL -lglut -lGLU -lGL"
    FCFLAGS="$FCFLAGS -I$PHAML_GRAPHICS_PREFIX/include/GL -L$PHAML_GRAPHICS_PREFIX/lib -L/usr/X11R6/lib -I/usr/include/GL"
    AC_MSG_RESULT(F90GL location: $PHAML_GRAPHICS_PREFIX)
    AC_MSG_CHECKING([for F90GL interface])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],[use opengl_gl])],AC_MSG_RESULT([yes]),AC_MSG_FAILURE([Can not link to PHAML Library]))
    AC_SUBST([phaml_graphics], ['phaml_graphics${EXEEXT}'])
  else
    AC_MSG_RESULT([no])
  fi

  LIBS=$LIBS_save
  LDFLAGS=$LDFLAGS_save
  FCFLAGS=$FCFLAGS_save
  build_extension_libphaml="true"
else
 AC_MSG_RESULT([no])
fi
AM_CONDITIONAL(BUILD_EXT_PHAML, test "$PHAML_PREFIX" != "")
AM_CONDITIONAL(USE_PHAML_GRAPHICS, test "$PHAML_GRAPHICS_PREFIX" != "")
AC_SUBST(PHAML_PREFIX)
AC_SUBST(PHAML_GRAPHICS_PREFIX)
AC_SUBST(build_extension_libphaml)

# Extra libraries can be added by user at configure time 
# This is used on Jaguar to add   f90main.o  for successful linking with CXX
EXTRA_LDLIBS_SUBST="$EXTRA_LDLIBS"
AC_SUBST(EXTRA_LDLIBS_SUBST)

#Checking for Trilinos
  AC_ARG_WITH(trilinos,AS_HELP_STRING([--with-trilinos],[Optional location of trilinos solvers]),[ac_cv_trilinos_prefix=$withval])
if test "$ac_cv_trilinos_prefix" != ""; then
        AC_MSG_RESULT([Using Trilinos to solve sparse linear systems, (located at $ac_cv_trilinos_prefix)])

  # automate findingof trilinos libraries, used by Makefile.Trilinos.export
  TRILINOS_PREFIX="$ac_cv_trilinos_prefix"
  AC_SUBST(TRILINOS_PREFIX)

  CPPFLAGS="$CPPFLAGS -DTRILINOS -DGLIMMER_MPI "
  CXXFLAGS="$CXXFLAGS -I$ac_cv_trilinos_prefix/include"
  LDFLAGS="$LDFLAGS -L$ac_cv_trilinos_prefix/lib"
fi
  AM_CONDITIONAL(TRILINOS, test "$ac_cv_trilinos_prefix" != "")
# end trilinos section

# Checking if sources are from subversion
AC_CHECK_FILE($srcdir/libglimmer/.svn/entries,have_svn_entries="yes",have_svn_entries="no")
GLIMMER_SVN_VERS=""
if test "$have_svn_entries" == "yes"; then
   GLIMMER_SVN_VERS=" svn_vers_string"
fi
AC_SUBST(GLIMMER_SVN_VERS)
AM_CONDITIONAL([HAVE_SVN],[test "$have_svn_entries" == "yes"])



FCCPPFLAGS=$(echo $CPPFLAGS | sed 's/-D/ $(FC_DEFINE)/g')
FCDEFS=$(echo $DEFS | sed 's/-D/ $(FC_DEFINE)/g')
#echo "FCCPPFLAGS=$FCCPPFLAGS"
#echo "FCDEFS=$FCDEFS"
AC_SUBST(FCCPPFLAGS)
AC_SUBST(FCDEFS)





AC_CONFIG_FILES([Makefile m4macros/Makefile scripts/Makefile \
                libglimmer-trilinos/Makefile.Trilinos.export libglimmer-trilinos/Makefile \
                utils/Makefile libglimmer/Makefile \
                libglimmer-solve/Makefile libglimmer-solve/SLAP/Makefile \
                libglimmer-solve/blas/Makefile libglint/Makefile \
                libglide/Makefile tests/Makefile example-drivers/Makefile extensions/Makefile \
                extensions/erosion/Makefile extensions/erosion/src/Makefile \
                extensions/erosion/scripts/Makefile \
                extensions/libphaml/Makefile \
                tests/EISMINT/Makefile tests/EISMINT/benchmarks/Makefile \
                tests/EISMINT/EISMINT-1/Makefile  tests/EISMINT/EISMINT-2/Makefile \
                tests/EISMINT/isos/Makefile  tests/EISMINT/lin-non_lin/Makefile \
                tests/EISMINT/erosion/Makefile \
                tests/exact-isothermal/Makefile tests/exact-isothermal/src/Makefile \
                tests/exact-isothermal/matlab/Makefile tests/exact-isothermal/scripts/Makefile \
                tests/exact-isothermal/tests/Makefile \
                tests/higher-order/Makefile tests/higher-order/dome/Makefile tests/higher-order/shelf/Makefile \
                PyGC/Makefile PyGC/PyGC/Makefile PyGC/PyGC/IO/Makefile PyGC/src/Makefile PyGC/progs/Makefile \
                PyGC/data/Makefile \
                example-drivers/eis/Makefile example-drivers/eis/src/Makefile \
                example-drivers/eis/fenscan/Makefile example-drivers/EISMINT-3/EISMINT-3/Makefile \
                example-drivers/EISMINT-3/Makefile example-drivers/EISMINT-3/src/Makefile \
                example-drivers/glint/Makefile example-drivers/glint/src/Makefile \
                example-drivers/simple_glide/Makefile example-drivers/simple_glide/src/Makefile \
                doc/Makefile doc/common/Makefile doc/dg/Makefile doc/dg/figures/Makefile \
                doc/ug/Makefile doc/ug/figs/Makefile \
                doc/num/Makefile doc/num/figs/Makefile doc/num/gnu/Makefile \
                doc/doclets/Makefile doc/doclets/figs/Makefile \
                doc/tut/Makefile doc/ext/Makefile \
                doc/ext/erosion/Makefile doc/ext/erosion/figs/Makefile doc/ext/erosion/gnu/Makefile \
                doc/ext/libphaml/Makefile doc/ext/libphaml/figures/Makefile \
                doc/ext/libphaml/figures/png/Makefile doc/ext/libphaml/figures/pdf/Makefile])
AC_CONFIG_FILES([doc/common/version.tex])
AC_CONFIG_FILES([glimmer-cism-config libglimmer/glimmer_vers.F90.in])
AC_OUTPUT
