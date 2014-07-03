#! /bin/csh -f 

cd $OBJROOT/ice/obj

set comp = 'unknown'
if ($COMP_INTERFACE == 'MCT' ) set comp = mct
if ($COMP_INTERFACE == 'ESMF') set comp = esmf

cat >! Filepath << EOF1
$CASEROOT/SourceMods/src.dice
$CODEROOT/ice/dice
$CODEROOT/ice/dice/cpl_$comp
EOF1

gmake complib -j $GMAKE_J MODEL=dice COMPLIB=$LIBROOT/libice.a -f $CASETOOLS/Makefile MACFILE=$CASEROOT/Macros.$MACH || exit 2

