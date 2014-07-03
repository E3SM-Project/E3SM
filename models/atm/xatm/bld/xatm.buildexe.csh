#! /bin/csh -f 

cd $OBJROOT/atm/obj

cat >! Filepath << EOF
$CASEROOT/SourceMods/src.xatm
$CODEROOT/atm/xatm
$CODEROOT/atm/xatm/cpl
EOF

gmake complib -j $GMAKE_J MODEL=xatm COMPLIB=$LIBROOT/libatm.a -f $CASETOOLS/Makefile MACFILE=$CASEROOT/Macros.$MACH || exit 2


