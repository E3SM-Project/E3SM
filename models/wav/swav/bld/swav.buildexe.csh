#! /bin/csh -f 

cd $OBJROOT/wav/obj

cat >! Filepath << EOF
$CASEROOT/SourceMods/src.swav
$CODEROOT/wav/swav
$CODEROOT/wav/swav/cpl
EOF

gmake complib -j $GMAKE_J MODEL=swav COMPLIB=$LIBROOT/libwav.a -f $CASETOOLS/Makefile MACFILE=$CASEROOT/Macros.$MACH || exit 2
