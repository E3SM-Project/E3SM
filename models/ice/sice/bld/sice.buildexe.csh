#! /bin/csh -f 

cd $OBJROOT/ice/obj

cat >! Filepath << EOF
$CASEROOT/SourceMods/src.sice
$CODEROOT/ice/sice
$CODEROOT/ice/sice/cpl
EOF

gmake complib -j $GMAKE_J MODEL=sice COMPLIB=$LIBROOT/libice.a -f $CASETOOLS/Makefile MACFILE=$CASEROOT/Macros.$MACH || exit 2

