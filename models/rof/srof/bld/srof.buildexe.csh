#! /bin/csh -f 

cd $OBJROOT/rof/obj

cat >! Filepath << EOF
$CASEROOT/SourceMods/src.srof
$CODEROOT/rof/srof
$CODEROOT/rof/srof/cpl
EOF

gmake complib -j $GMAKE_J MODEL=srof COMPLIB=$LIBROOT/librof.a -f $CASETOOLS/Makefile MACFILE=$CASEROOT/Macros.$MACH || exit 2

