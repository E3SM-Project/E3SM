#!/bin/csh
#===============================================================================
# SVN $Id: build.yellowstone.csh 46158 2013-04-19 18:41:34Z mlevy@ucar.edu $
# SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/tools/mapping/trunk_tags/mapping_140422a/gen_mapping_files/runoff_to_ocn/build.yellowstone.csh $
#===============================================================================
# 
# Notes:
# - will build the CCSM runoff correcting/smoothing code in ./obj
# - must specify location of src code, Makefile, Macros file, dependancy generator 
#===============================================================================

setenv SRCDIR `pwd`/src
setenv TOOLDIR `pwd`/tools

echo source dir: $SRCDIR

if !(-d obj) mkdir obj
cd obj

cc -o makdep $TOOLDIR/makdep.c

echo $SRCDIR  >! Filepath

gmake VPFILE=Filepath THREAD=TRUE -f $SRCDIR/Makefile MACFILE=$SRCDIR/Macros.yellowstone || exit -1

cd ..
rm              runoff_map
ln -s obj/a.out runoff_map

