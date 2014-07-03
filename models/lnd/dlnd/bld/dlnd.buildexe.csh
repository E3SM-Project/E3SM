#! /bin/csh -f 

cd $OBJROOT/lnd/obj

set comp = 'unknown'
if ($COMP_INTERFACE == 'MCT' ) set comp = mct
if ($COMP_INTERFACE == 'ESMF') set comp = esmf

cat >! Filepath << EOF1
$CASEROOT/SourceMods/src.dlnd
$CODEROOT/lnd/dlnd
$CODEROOT/lnd/dlnd/cpl_$comp
EOF1

gmake complib -j $GMAKE_J MODEL=dlnd COMPLIB=$LIBROOT/liblnd.a -f $CASETOOLS/Makefile MACFILE=$CASEROOT/Macros.$MACH || exit 2

