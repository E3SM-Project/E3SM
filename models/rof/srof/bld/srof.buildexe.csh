#! /bin/csh -f 

cd $OBJROOT/rof/obj

set comp = 'unknown'
if ($COMP_INTERFACE == 'MCT' ) set comp = mct
if ($COMP_INTERFACE == 'ESMF') set comp = esmf

cat >! Filepath << EOF
$CASEROOT/SourceMods/src.srof
$CODEROOT/rof/srof
$CODEROOT/rof/srof/cpl_$comp
EOF

gmake complib -j $GMAKE_J MODEL=srof COMPLIB=$LIBROOT/librof.a -f $CASETOOLS/Makefile MACFILE=$CASEROOT/Macros.$MACH || exit 2

