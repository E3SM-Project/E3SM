#! /bin/csh -f 

cd  $OBJROOT/glc/obj

set comp = 'unknown'
if ($COMP_INTERFACE == 'MCT' ) set comp = mct
if ($COMP_INTERFACE == 'ESMF') set comp = esmf

cat >! Filepath << EOF
$CASEROOT/SourceMods/src.xglc
$CODEROOT/glc/xglc
$CODEROOT/glc/xglc/cpl_$comp
EOF

gmake complib -j $GMAKE_J MODEL=xglc COMPLIB=$LIBROOT/libglc.a -f $CASETOOLS/Makefile MACFILE=$CASEROOT/Macros.$MACH || exit 2
