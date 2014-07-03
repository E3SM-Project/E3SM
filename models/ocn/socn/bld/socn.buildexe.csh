#! /bin/csh -f 

cd $OBJROOT/ocn/obj
set comp = 'unknown'
if ($COMP_INTERFACE == 'MCT' ) set comp = mct
if ($COMP_INTERFACE == 'ESMF') set comp = esmf

cat >! Filepath << EOF
$CASEROOT/SourceMods/src.socn
$CODEROOT/ocn/socn
$CODEROOT/ocn/socn/cpl_$comp
EOF

gmake complib -j $GMAKE_J MODEL=socn COMPLIB=$LIBROOT/libocn.a -f $CASETOOLS/Makefile MACFILE=$CASEROOT/Macros.$MACH || exit 2



