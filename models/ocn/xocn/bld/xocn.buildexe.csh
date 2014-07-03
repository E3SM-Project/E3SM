#! /bin/csh -f 

cd $OBJROOT/ocn/obj

cat >! Filepath << EOF
$CASEROOT/SourceMods/src.xocn
$CODEROOT/ocn/xocn
$CODEROOT/ocn/xocn/cpl
EOF

gmake complib -j $GMAKE_J MODEL=xocn COMPLIB=$LIBROOT/libocn.a -f $CASETOOLS/Makefile MACFILE=$CASEROOT/Macros.$MACH || exit 2

