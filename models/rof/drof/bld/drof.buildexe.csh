#! /bin/csh -f 

cd $OBJROOT/rof/obj

set comp = 'unknown'
if ($COMP_INTERFACE == 'MCT' ) set comp = mct
if ($COMP_INTERFACE == 'ESMF') set comp = esmf

cat >! Filepath << EOF1
$CASEROOT/SourceMods/src.drof
$CODEROOT/rof/drof
$CODEROOT/rof/drof/cpl_$comp
EOF1

gmake complib -j $GMAKE_J MODEL=drof COMPLIB=$LIBROOT/librof.a -f $CASETOOLS/Makefile MACFILE=$CASEROOT/Macros.$MACH || exit 2

