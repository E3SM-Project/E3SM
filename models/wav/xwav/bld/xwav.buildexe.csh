#! /bin/csh -f 

cd $OBJROOT/wav/obj

cat >! Filepath << EOF
$CASEROOT/SourceMods/src.xwav
$CODEROOT/wav/xwav
$CODEROOT/wav/xwav/cpl
EOF

gmake complib -j $GMAKE_J MODEL=xwav COMPLIB=$LIBROOT/libwav.a -f $CASETOOLS/Makefile MACFILE=$CASEROOT/Macros.$MACH || exit 2


