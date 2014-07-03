#! /bin/csh -f 

cd $OBJROOT/atm/obj

cat >! Filepath << EOF
$CASEROOT/SourceMods/src.satm
$CODEROOT/atm/satm
$CODEROOT/atm/satm/cpl
EOF

gmake complib -j $GMAKE_J MODEL=satm COMPLIB=$LIBROOT/libatm.a -f $CASETOOLS/Makefile MACFILE=$CASEROOT/Macros.$MACH || exit 2
