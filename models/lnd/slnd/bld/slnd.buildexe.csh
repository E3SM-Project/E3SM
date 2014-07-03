#! /bin/csh -f 

cd $OBJROOT/lnd/obj

set comp = 'unknown'
if ($COMP_INTERFACE == 'MCT' ) set comp = mct
if ($COMP_INTERFACE == 'ESMF') set comp = esmf

cat >! Filepath << EOF
$CASEROOT/SourceMods/src.slnd
$CODEROOT/lnd/slnd
$CODEROOT/lnd/slnd/cpl_$comp
EOF

gmake complib -j $GMAKE_J MODEL=slnd COMPLIB=$LIBROOT/liblnd.a -f $CASETOOLS/Makefile MACFILE=$CASEROOT/Macros.$MACH || exit 2

