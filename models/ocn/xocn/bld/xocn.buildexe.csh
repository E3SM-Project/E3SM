#! /bin/csh -f 

cd $OBJROOT/ocn/obj

set comp = 'unknown'
if ($COMP_INTERFACE == 'MCT' ) set comp = mct
if ($COMP_INTERFACE == 'ESMF') set comp = esmf

cat >! Filepath << EOF
$CASEROOT/SourceMods/src.xocn
$CODEROOT/ocn/xocn
$CODEROOT/ocn/xocn/cpl_$comp
EOF

gmake complib -j $GMAKE_J MODEL=xocn COMPLIB=$LIBROOT/libocn.a -f $CASETOOLS/Makefile MACFILE=$CASEROOT/Macros.$MACH || exit 2

