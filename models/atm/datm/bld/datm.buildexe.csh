#! /bin/csh -f 

cd $OBJROOT/atm/obj

set comp = 'unknown'
if ($COMP_INTERFACE == 'MCT' ) set comp = mct
if ($COMP_INTERFACE == 'ESMF') set comp = esmf

cat >! Filepath << EOF1
$CASEROOT/SourceMods/src.datm
$CODEROOT/atm/datm
$CODEROOT/atm/datm/cpl_$comp
EOF1

gmake complib -j $GMAKE_J MODEL=datm COMPLIB=$LIBROOT/libatm.a -f $CASETOOLS/Makefile MACFILE=$CASEROOT/Macros.$MACH || exit 2


