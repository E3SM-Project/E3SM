#! /bin/csh -f 

cd $OBJROOT/wav/obj

set comp = 'unknown'
if ($COMP_INTERFACE == 'MCT' ) set comp = mct
if ($COMP_INTERFACE == 'ESMF') set comp = esmf

cat >! Filepath << EOF
$CASEROOT/SourceMods/src.swav
$CODEROOT/wav/swav
$CODEROOT/wav/swav/cpl_$comp
EOF

gmake complib -j $GMAKE_J MODEL=swav COMPLIB=$LIBROOT/libwav.a -f $CASETOOLS/Makefile MACFILE=$CASEROOT/Macros.$MACH || exit 2
