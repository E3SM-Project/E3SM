#! /bin/csh -f

cd $EXEROOT/cesm/obj

echo -------------------------------------------------------------------------
echo " Building a single executable version of CESM "
echo -------------------------------------------------------------------------

cat >! Filepath << EOF
$CASEROOT/SourceMods/src.drv
$CCSMROOT/models/drv/driver
EOF

## Set up multiple component instances
## WJS (2-9-11): I think these changes actually are not required here
## any more, since seq_multiinst_mod is built as part of csm_share;
## but I am leaving them in to be safe

set multiinst_cppdefs = "-DNUM_COMP_INST_ATM=$NINST_ATM -DNUM_COMP_INST_LND=$NINST_LND -DNUM_COMP_INST_OCN=$NINST_OCN -DNUM_COMP_INST_ICE=$NINST_ICE -DNUM_COMP_INST_GLC=$NINST_GLC -DNUM_COMP_INST_WAV=$NINST_WAV"

gmake exec_se -j $GMAKE_J EXEC_SE=$EXEROOT/cesm.exe MODEL=driver \
                 USER_CPPDEFS="$multiinst_cppdefs" -f $CASETOOLS/Makefile || exit 2

exit 0
