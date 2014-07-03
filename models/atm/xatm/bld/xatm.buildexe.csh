#! /bin/csh -f 

cd $OBJROOT/atm/obj

set comp = 'unknown'
if ($COMP_INTERFACE == 'MCT' ) set comp = mct
if ($COMP_INTERFACE == 'ESMF') set comp = esmf

cat >! Filepath << EOF
$CASEROOT/SourceMods/src.xatm
$CODEROOT/atm/xatm
$CODEROOT/atm/xatm/cpl_$comp
EOF

gmake complib -j $GMAKE_J MODEL=xatm COMPLIB=$LIBROOT/libatm.a -f $CASETOOLS/Makefile MACFILE=$CASEROOT/Macros.$MACH || exit 2


