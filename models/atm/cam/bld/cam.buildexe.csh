#! /bin/csh -f

cd $OBJROOT/atm/obj

cp $CASEBUILD/camconf/Filepath ./tmp_filepath 
if (-f Filepath) then
  cmp -s tmp_filepath Filepath || mv -f tmp_filepath Filepath 
else
  mv -f tmp_filepath Filepath 
endif

set camdefs = "`cat $CASEBUILD/camconf/CCSM_cppdefs`"
gmake complib -j $GMAKE_J MODEL=cam COMPLIB=$LIBROOT/libatm.a USER_CPPDEFS="$camdefs" -f $CASETOOLS/Makefile || exit 2

wait

