#! /bin/csh -f 

cd $OBJROOT/lnd/obj

cat >! Filepath << EOF1
$CASEROOT/SourceMods/src.dlnd
$CODEROOT/lnd/dlnd
EOF1

gmake complib -j $GMAKE_J MODEL=dlnd COMPLIB=$LIBROOT/liblnd.a -f $CASETOOLS/Makefile MACFILE=$CASEROOT/Macros.$MACH || exit 2

