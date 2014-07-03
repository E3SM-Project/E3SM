#! /bin/csh -f 

cd $OBJROOT/ice/obj

cat >! Filepath << EOF
$CASEROOT/SourceMods/src.xice
$CODEROOT/ice/xice
$CODEROOT/ice/xice/cpl
EOF

gmake complib -j $GMAKE_J MODEL=xice COMPLIB=$LIBROOT/libice.a -f $CASETOOLS/Makefile MACFILE=$CASEROOT/Macros.$MACH || exit 2
