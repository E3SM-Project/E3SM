#! /bin/csh -f 

cd $OBJROOT/ice/obj

set comp = 'unknown'
if ($COMP_INTERFACE == 'MCT' ) set comp = mct
if ($COMP_INTERFACE == 'ESMF') set comp = esmf

cat >! Filepath << EOF
$CASEROOT/SourceMods/src.sice
$CODEROOT/ice/sice
$CODEROOT/ice/sice/cpl_$comp
EOF

gmake complib -j $GMAKE_J MODEL=sice COMPLIB=$LIBROOT/libice.a -f $CASETOOLS/Makefile MACFILE=$CASEROOT/Macros.$MACH || exit 2

