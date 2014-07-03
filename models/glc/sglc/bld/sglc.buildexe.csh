#! /bin/csh -f 

cd $OBJROOT/glc/obj

cat >! Filepath << EOF
$CASEROOT/SourceMods/src.sglc
$CODEROOT/glc/sglc
$CODEROOT/glc/sglc/cpl
EOF

gmake complib -j $GMAKE_J MODEL=sglc COMPLIB=$LIBROOT/libglc.a -f $CASETOOLS/Makefile MACFILE=$CASEROOT/Macros.$MACH || exit 2

