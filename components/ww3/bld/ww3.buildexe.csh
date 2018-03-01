#! /bin/csh -f

cd $OBJROOT/wav/obj

set comp = 'unknown'
if ($COMP_INTERFACE == 'MCT' ) set comp = mct
if ($COMP_INTERFACE == 'ESMF') set comp = esmf

#------------------------------------------------------------------------------
# Filepath: List of source code directories (in order of importance).
#------------------------------------------------------------------------------

\cat >! Filepath << EOF1
$CASEROOT/SourceMods/src.ww3
$CODEROOT/wav/ww3/src/source
$CODEROOT/wav/ww3/src/cpl_share
$CODEROOT/wav/ww3/src/cpl_$comp
EOF1

#------------------------------------------------------------------------------
# run make
#------------------------------------------------------------------------------

gmake complib -j $GMAKE_J MODEL=ww3 COMPLIB=$LIBROOT/libwav.a -f $CASETOOLS/Makefile MACFILE=$CASEROOT/Macros.$MACH || exit 2


