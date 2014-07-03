#! /bin/csh -f 

cd $OBJROOT/ocn/obj


cat >! Filepath << EOF1
$CASEROOT/SourceMods/src.docn
$CODEROOT/ocn/docn
EOF1

gmake complib -j $GMAKE_J MODEL=docn COMPLIB=$LIBROOT/libocn.a -f $CASETOOLS/Makefile MACFILE=$CASEROOT/Macros.$MACH || exit 2


