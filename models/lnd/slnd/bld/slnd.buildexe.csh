#! /bin/csh -f 

cd $OBJROOT/lnd/obj

cat >! Filepath << EOF
$CASEROOT/SourceMods/src.slnd
$CODEROOT/lnd/slnd
$CODEROOT/lnd/slnd/cpl
EOF

gmake complib -j $GMAKE_J MODEL=slnd COMPLIB=$LIBROOT/liblnd.a -f $CASETOOLS/Makefile MACFILE=$CASEROOT/Macros.$MACH || exit 2

