#! /bin/csh -f 

cd $OBJROOT/rof/obj

cat >! Filepath << EOF1
$CASEROOT/SourceMods/src.drof
$CODEROOT/rof/drof
EOF1

gmake complib -j $GMAKE_J MODEL=drof COMPLIB=$LIBROOT/librof.a -f $CASETOOLS/Makefile MACFILE=$CASEROOT/Macros.$MACH || exit 2

