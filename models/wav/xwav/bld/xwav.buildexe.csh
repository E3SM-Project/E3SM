#! /bin/csh -f 

cd $OBJROOT/wav/obj

set comp = 'unknown'
if ($COMP_INTERFACE == 'MCT' ) set comp = mct
if ($COMP_INTERFACE == 'ESMF') set comp = esmf

cat >! Filepath << EOF
$CASEROOT/SourceMods/src.xwav
$CODEROOT/wav/xwav
$CODEROOT/wav/xwav/cpl_$comp
EOF

gmake complib -j $GMAKE_J MODEL=xwav COMPLIB=$LIBROOT/libwav.a -f $CASETOOLS/Makefile MACFILE=$CASEROOT/Macros.$MACH || exit 2


