#! /bin/csh -f 

cd $OBJROOT/lnd/obj

if (-f $CASEBUILD/clmconf/Filepath) then
   cp $CASEBUILD/clmconf/Filepath ./tmp_filepath 
else
   echo "clm.buildexe.csh ERROR - missing $CASEBUILD/clmconf/Filepath"
   exit -1
endif
if (-f Filepath) then
  cmp -s tmp_filepath Filepath || mv -f tmp_filepath Filepath 
else
  mv -f tmp_filepath Filepath 
endif

set clmdefs = "`cat $CASEBUILD/clmconf/CESM_cppdefs`"
$GMAKE complib -j $GMAKE_J MODEL=clm COMPLIB=$LIBROOT/liblnd.a USER_CPPDEFS="$clmdefs" -f $CASETOOLS/Makefile || exit 2

wait



