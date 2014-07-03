#! /bin/csh -f 

cd $OBJROOT/rof/obj

set comp = 'unknown'
if ($COMP_INTERFACE == 'MCT' ) set comp = mct
if ($COMP_INTERFACE == 'ESMF') set comp = esmf

cat >! tmp_filepath << EOF1
$CASEROOT/SourceMods/src.rtm
$CODEROOT/rof/rtm/src/riverroute
$CODEROOT/rof/rtm/src/cpl_$comp
$CODEROOT/rof/rtm/src/cpl_share
EOF1

if (-f Filepath) then
  cmp -s tmp_filepath Filepath || mv -f tmp_filepath Filepath 
else
  mv -f tmp_filepath Filepath 
endif

$GMAKE complib -j $GMAKE_J MODEL=rtm COMPLIB=$LIBROOT/librof.a -f $CASETOOLS/Makefile || exit 2

wait



