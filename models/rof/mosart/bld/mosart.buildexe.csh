#! /bin/csh -f 

cd $OBJROOT/rof/obj

#$CODEROOT/rof/mosart/src/wrm

cat >! tmp_filepath << EOF1
$CASEROOT/SourceMods/src.mosart
$CODEROOT/rof/mosart/src/riverroute
$CODEROOT/rof/mosart/src/cpl
$CODEROOT/rof/mosart/src/cpl_share
EOF1

if (-f Filepath) then
  cmp -s tmp_filepath Filepath || mv -f tmp_filepath Filepath 
else
  mv -f tmp_filepath Filepath 
endif

#$GMAKE complib -j $GMAKE_J MODEL=mosart COMPLIB=$LIBROOT/librof.a USER_CPPDEFS="-DINCLUDE_WRM" -f $CASETOOLS/Makefile || exit 2
$GMAKE complib -j $GMAKE_J MODEL=mosart COMPLIB=$LIBROOT/librof.a -f $CASETOOLS/Makefile || exit 2

wait



