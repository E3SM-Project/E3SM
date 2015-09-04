#!/bin/csh -f

#------------------------------------------------------------------------------------
# For now, set streams manually. You must only set as many streams as are declared
#  in the tavg_nml section. For example, if there are three streams:
#  @ s1 = $my_stream
#  @ s2 = $s1 + 1
#  @ s3 = $s2 + 1
#------------------------------------------------------------------------------------

@ my_stream = $1
if ($my_stream < 1) then
   echo invalid my_stream number  ($my_stream)
   exit 5
endif

@ s1 = 1   # use base-model stream 1

cat >! $CASEROOT/Buildconf/pop2conf/cfc_tavg_contents << EOF
$s1  CFC_IFRAC
$s1  CFC_XKW
$s1  CFC_ATM_PRESS
$s1  STF_CFC11
$s1  STF_CFC12
$s1  CFC11
$s1  CFC12
EOF

if ($OCN_TAVG_TRACER_BUDGET == TRUE) then
cat >> $CASEROOT/Buildconf/pop2conf/cfc_tavg_contents << EOF
$s1  KPP_SRC_CFC11
$s1  KPP_SRC_CFC12
$s1  DIA_IMPVF_CFC11
$s1  DIA_IMPVF_CFC12
$s1  HDIFE_CFC11
$s1  HDIFE_CFC12
$s1  HDIFN_CFC11
$s1  HDIFN_CFC12
$s1  HDIFB_CFC11
$s1  HDIFB_CFC12
$s1  UE_CFC11
$s1  UE_CFC12
$s1  VN_CFC11
$s1  VN_CFC12
$s1  WT_CFC11
$s1  WT_CFC12
EOF
endif

#===============================================================================
# The following are fields computed by the CFC modules that are not placed in
# the tavg file by default.
#
#1  pCFC11
#1  pCFC12
#1  CFC11_SCHMIDT
#1  CFC12_SCHMIDT
#1  CFC11_PV
#1  CFC11_surf_sat
#1  CFC12_PV
#1  CFC12_surf_sat
#===============================================================================
