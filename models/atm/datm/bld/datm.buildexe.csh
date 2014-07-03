#! /bin/csh -f 

cd $OBJROOT/atm/obj

cat >! Filepath << EOF1
$CASEROOT/SourceMods/src.datm
$CODEROOT/atm/datm
EOF1

gmake complib -j $GMAKE_J MODEL=datm COMPLIB=$LIBROOT/libatm.a -f $CASETOOLS/Makefile MACFILE=$CASEROOT/Macros.$MACH || exit 2


