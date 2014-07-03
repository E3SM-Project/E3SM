#! /bin/csh -f 

cd $OBJROOT/ice/obj

set comp = 'unknown'
if ($COMP_INTERFACE == 'MCT' ) set comp = mct
if ($COMP_INTERFACE == 'ESMF') set comp = esmf

cat >! Filepath << EOF
$CASEROOT/SourceMods/src.xice
$CODEROOT/ice/xice
$CODEROOT/ice/xice/cpl_$comp
EOF

gmake complib -j $GMAKE_J MODEL=xice COMPLIB=$LIBROOT/libice.a -f $CASETOOLS/Makefile MACFILE=$CASEROOT/Macros.$MACH || exit 2
