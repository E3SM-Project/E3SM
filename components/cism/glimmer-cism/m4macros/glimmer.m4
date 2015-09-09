# Configure paths for GLIMMER
# Magnus Hagdorn
# Large parts are shamelessly stolen from the equivalent gimp-2.0 macro
# written by Manish Singh and Sven Neumann

dnl AM_PATH_GLIMMER([MINIMUM-VERSION, [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl Test for GLIMMER, abd define 
dnl GLIMMER_PREFIX, GLIMMER_FFLAGS and GLIMMER_LIBS
dnl
AC_DEFUN([AM_PATH_GLIMMER],
  [dnl
  dnl Find out where GLIMMER is installed

  no_glimmer="no"

  AC_ARG_WITH(glimmer-prefix,[  --with-glimmer-prefix=PFX   Prefix where GLIMMER is installed (optional)],
              glimmer_config_prefix="$withval", glimmer_config_prefix="")

  if test x$glimmer_config_prefix != x ; then
       glimmer_config_args="$glimmer_config_args --prefix=$glimmer_config_prefix"
       if test x${GLIMMER_CONFIG+set} != xset ; then
          GLIMMER_CONFIG=$glimmer_config_prefix/bin/glimmer-cism-config
       fi
  fi

  AC_PATH_PROG(GLIMMER_CONFIG, glimmer-cism-config, no)
  if ! test -x "$GLIMMER_CONFIG"; then
    GLIMMER_CONFIG="no"
  fi
  if test "$GLIMMER_CONFIG" = "no" ; then
       no_glimmer="yes"
  else

       dnl Checking for version of glimmer
       min_glimmer_version=$1
       min_glimmer_major_version=`echo $min_glimmer_version | \
              sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\1/'`
       min_glimmer_minor_version=`echo $min_glimmer_version | \
              sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\2/'`
       min_glimmer_micro_version=`echo $min_glimmer_version | \
              sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\3/'`

       AC_MSG_CHECKING(for GLIMMER - version >= $min_glimmer_version)
       glimmer_config_major_version=`$GLIMMER_CONFIG $glimmer_config_args --version | \
              sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\1/'`
       glimmer_config_minor_version=`$GLIMMER_CONFIG $glimmer_config_args --version | \
              sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\2/'`
       glimmer_config_micro_version=`$GLIMMER_CONFIG $glimmer_config_args --version | \
              sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\3/'`
       if test  "$min_glimmer_micro_version" -gt "$glimmer_config_micro_version"; then
           no_glimmer="yes"
       fi
       if test  "$min_glimmer_minor_version" -gt "$glimmer_config_minor_version"; then
           no_glimmer="yes"
       elif test  "$min_glimmer_minor_version" -lt "$glimmer_config_minor_version"; then
           no_glimmer="no"
       fi
       if test  "$min_glimmer_major_version" -gt "$glimmer_config_major_version"; then
           no_glimmer="yes"
       elif test "$min_glimmer_major_version" -lt "$glimmer_config_major_version"; then
           no_glimmer="no"
       fi
  fi

  if test "$no_glimmer" = "no" ; then
      GLIMMER_PREFIX=`$GLIMMER_CONFIG $glimmer_config_args --prefix`
      GLIMMER_FFLAGS=`$GLIMMER_CONFIG $glimmer_config_args --fcflags`
      GLIMMER_LIBS=`$GLIMMER_CONFIG $glimmer_config_args --libs`

      AC_MSG_RESULT([yes $glimmer_config_major_version.$glimmer_config_minor_version.$glimmer_config_micro_version])

      ifelse([$2], , :, [$2])
  else
      AC_MSG_RESULT([no])
      if test "$GLIMMER_CONFIG" = "no" ; then
          echo "*** Could not find glimmer-cism-config. Either add the directory where glimmer-cism-config is installed"
          echo "*** in your PATH or run configure with the --glimmer_prefix option."
      else
          echo "*** Found glimmer version $glimmer_config_major_version.$glimmer_config_minor_version.$glimmer_config_micro_version which is too old"
          echo "*** Get the latest version of glimmer from https://developer.berlios.de/projects/glimmer-cism/"
      fi   
      GLIMMER_PREFIX=""
      GLIMMER_FFLAGS=""
      GLIMMER_LIBS=""
      ifelse([$3], , :, [$3])
  fi

  AC_SUBST(GLIMMER_PREFIX)
  AC_SUBST(GLIMMER_FFLAGS)
  AC_SUBST(GLIMMER_LIBS)

])
