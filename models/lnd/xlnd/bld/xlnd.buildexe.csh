#! /bin/csh -f 

cd $OBJROOT/lnd/obj

cat >! Filepath << EOF
$CASEROOT/SourceMods/src.xlnd
$CODEROOT/lnd/xlnd
$CODEROOT/lnd/xlnd/cpl
EOF

gmake complib -j $GMAKE_J MODEL=xlnd COMPLIB=$LIBROOT/liblnd.a -f $CASETOOLS/Makefile MACFILE=$CASEROOT/Macros.$MACH || exit 2

