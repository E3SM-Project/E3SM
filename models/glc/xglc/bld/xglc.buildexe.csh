#! /bin/csh -f 

cd  $OBJROOT/glc/obj

cat >! Filepath << EOF
$CASEROOT/SourceMods/src.xglc
$CODEROOT/glc/xglc
$CODEROOT/glc/xglc/cpl
EOF

gmake complib -j $GMAKE_J MODEL=xglc COMPLIB=$LIBROOT/libglc.a -f $CASETOOLS/Makefile MACFILE=$CASEROOT/Macros.$MACH || exit 2
