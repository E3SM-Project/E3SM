#! /bin/csh -f 

cd $OBJROOT/ocn/obj

set comp = 'unknown'
if ($COMP_INTERFACE == 'MCT' ) set comp = mct
if ($COMP_INTERFACE == 'ESMF') set comp = esmf

cat >! Filepath << EOF1
$CASEROOT/SourceMods/src.docn
$CODEROOT/ocn/docn
$CODEROOT/ocn/docn/cpl_$comp
EOF1

gmake complib -j $GMAKE_J MODEL=docn COMPLIB=$LIBROOT/libocn.a -f $CASETOOLS/Makefile MACFILE=$CASEROOT/Macros.$MACH || exit 2


