#! /bin/csh -f 

cd $OBJROOT/rof/obj

cat >! Filepath << EOF
$CASEROOT/SourceMods/src.xrof
$CODEROOT/rof/xrof
$CODEROOT/rof/xrof/cpl
EOF

gmake complib -j $GMAKE_J MODEL=xrof COMPLIB=$LIBROOT/librof.a -f $CASETOOLS/Makefile MACFILE=$CASEROOT/Macros.$MACH || exit 2
