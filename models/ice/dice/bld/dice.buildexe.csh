#! /bin/csh -f 

cd $OBJROOT/ice/obj

cat >! Filepath << EOF1
$CASEROOT/SourceMods/src.dice
$CODEROOT/ice/dice
EOF1

gmake complib -j $GMAKE_J MODEL=dice COMPLIB=$LIBROOT/libice.a -f $CASETOOLS/Makefile MACFILE=$CASEROOT/Macros.$MACH || exit 2

