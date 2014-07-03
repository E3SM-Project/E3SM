#! /bin/csh -f 

cd $OBJROOT/ocn/obj

cat >! Filepath << EOF
$CASEROOT/SourceMods/src.socn
$CODEROOT/ocn/socn
$CODEROOT/ocn/socn/cpl
EOF

gmake complib -j $GMAKE_J MODEL=socn COMPLIB=$LIBROOT/libocn.a -f $CASETOOLS/Makefile MACFILE=$CASEROOT/Macros.$MACH || exit 2



