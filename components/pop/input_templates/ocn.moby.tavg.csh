#!/bin/csh -f

cat >! $CASEROOT/Buildconf/pop2conf/moby.tavg.nml << EOF
tavg_freq_opt             = 'nday'          'nyear'
tavg_freq                 =  1              1
tavg_stream_filestrings   = 'darwin.nday1'  'darwin.nyear1'
tavg_file_freq_opt        = 'nmonth'        'nyear'
tavg_file_freq            =  1              1
tavg_start_opt            = 'nstep'         'nstep'
tavg_start                =  0              0
tavg_fmt_in               = 'nc'            'nc'
tavg_fmt_out              = 'nc'            'nc'
ltavg_has_offset_date     = .false.         .false.
tavg_offset_years         =  1              1
tavg_offset_months        =  1              1
tavg_offset_days          =  2              2
ltavg_one_time_header     = .false.         .false.
EOF

#------------------------------------------------------------------------------------
# For now, set streams manually. You must only set AS MANY STREAMS AS ARE DEFINED
#  in the tavg_nml section. For example, if there are three streams, set:
#  @ s1 = $my_stream
#  @ s2 = $s1 + 1
#  @ s3 = $s2 + 1
#------------------------------------------------------------------------------------

    @ my_stream = $1
    if ($my_stream < 1) then
       echo invalid my_stream number  ($my_stream)
       exit 5
    endif
    @ s1 = 1             # use the base-model stream 1
    @ s2 = $my_stream    # use a moby-defined stream
    @ s3 = $s2 + 1       # use a moby-defined stream

cat >! $CASEROOT/Buildconf/pop2conf/moby_tavg_contents << EOF
$s1  PO4
$s1  NO3
$s1  FeT
$s1  SiO2
$s1  DOP
$s1  DON
$s1  DOFe
$s1  ZOO1P
$s1  ZOO1N
$s1  ZOO1Fe
$s1  ZOO1Si
$s1  ZOO2P
$s1  ZOO2N
$s1  ZOO2Fe
$s1  ZOO2Si
$s1  POP
$s1  PON
$s1  POFe
$s1  POSi
$s1  NH4
$s1  NO2
$s1  Phy01
$s1  Phy02
$s1  Phy03
$s1  Phy04
$s1  Phy05
$s1  Phy06
$s1  Phy07
$s1  Phy08
$s1  Phy09
$s1  Chl01
$s1  Chl02
$s1  Chl03
$s1  Chl04
$s1  Chl05
$s1  Chl06
$s1  Chl07
$s1  Chl08
$s1  Chl09
$s1  DIC
$s1  DOC
$s1  POC
$s1  PIC
$s1  ALK
$s1  O2
$s1  ZOOC1
$s1  ZOOC2
#s1  MOBY_CO2_FLUX # same as FG_CO2
$s1  FG_CO2
$s1  STF_O2
$s1  FvPER_DIC
$s1  FvICE_DIC
$s1  FvPER_ALK
$s1  FvICE_ALK
$s1  J_DIC
$s1  Jint_100m_PO4
$s1  Jint_100m_NO3
$s1  Jint_100m_FeT
$s1  Jint_100m_SiO2
$s1  Jint_100m_NH4
$s1  Jint_100m_NO2
$s1  Jint_100m_DIC
$s1  Jint_100m_DOC
$s1  Jint_100m_ALK
$s1  Jint_100m_O2
$s1  Jint_100m_ZOOC1
$s1  Jint_100m_ZOOC2
$s1  tend_zint_100m_PO4
$s1  tend_zint_100m_NO3
$s1  tend_zint_100m_FeT
$s1  tend_zint_100m_SiO2
$s1  tend_zint_100m_NH4
$s1  tend_zint_100m_NO2
$s1  tend_zint_100m_DIC
$s1  tend_zint_100m_DOC
$s1  tend_zint_100m_POC
$s1  tend_zint_100m_PIC
$s1  tend_zint_100m_ALK
$s1  tend_zint_100m_O2
$s2  Phy01_zint_100m
$s2  Phy02_zint_100m
$s2  Phy03_zint_100m
$s2  Phy04_zint_100m
$s2  Phy05_zint_100m
$s2  Phy06_zint_100m
$s2  Phy07_zint_100m
$s2  Phy08_zint_100m
$s2  Phy09_zint_100m
$s2  ZOO1P_zint_100m
$s2  ZOO1N_zint_100m
$s2  ZOO1Fe_zint_100m
$s2  ZOO1Si_zint_100m
$s2  ZOO2P_zint_100m
$s2  ZOO2N_zint_100m
$s2  ZOO2Fe_zint_100m
$s2  ZOO2Si_zint_100m
$s2  Chl01_SURF
$s2  Chl02_SURF
$s2  Chl03_SURF
$s2  Chl04_SURF
$s2  Chl05_SURF
$s2  Chl06_SURF
$s2  Chl07_SURF
$s2  Chl08_SURF
$s2  Chl09_SURF
$s3  J_PO4
$s3  J_NO3
$s3  J_FeT
$s3  J_SiO2
$s3  J_NH4
$s3  J_NO2
$s3  J_ALK
$s3  UE_FeT
$s3  VN_FeT
$s3  WT_FeT
$s3  UE_DIC
$s3  VN_DIC
$s3  WT_DIC
$s3  KPP_SRC_DIC
$s3  DIA_IMPVF_DIC
$s3  HDIFE_DIC
$s3  HDIFN_DIC
$s3  HDIFB_DIC
$s3  UE_DOC
$s3  VN_DOC
$s3  WT_DOC
$s3  DIA_IMPVF_DOC
$s3  HDIFE_DOC
$s3  HDIFN_DOC
$s3  HDIFB_DOC
$s3  UE_POC
$s3  VN_POC
$s3  WT_POC
$s3  DIA_IMPVF_POC
$s3  HDIFE_POC
$s3  HDIFN_POC
$s3  HDIFB_POC
$s3  UE_PIC
$s3  VN_PIC
$s3  WT_PIC
$s3  DIA_IMPVF_PIC
$s3  HDIFE_PIC
$s3  HDIFN_PIC
$s3  HDIFB_PIC
$s3  KPP_SRC_FeT
$s3  DIA_IMPVF_FeT
$s3  HDIFE_FeT
$s3  HDIFN_FeT
$s3  HDIFB_FeT
$s3  UE_O2
$s3  VN_O2
$s3  WT_O2
$s3  KPP_SRC_O2
$s3  DIA_IMPVF_O2
$s3  HDIFE_O2
$s3  HDIFN_O2
$s3  HDIFB_O2
# NOTE: the following fields are commented out; to uncomment, delete the leading #
#$s1  Jint_PO4
#$s1  Jint_NO3
#$s1  Jint_Fe
#$s1  Jint_SiO2
#$s1  Jint_ZOO1P
#$s1  Jint_ZOO1N
#$s1  Jint_ZOO1Fe
#$s1  Jint_ZOO1Si
#$s1  Jint_ZOO2P
#$s1  Jint_ZOO2N
#$s1  Jint_ZOO2Fe
#$s1  Jint_ZOO2Si
#$s1  Jint_NH4
#$s1  Jint_NO2
#$s1  Jint_Phy01
#$s1  Jint_Phy02
#$s1  Jint_Phy03
#$s1  Jint_Phy04
#$s1  Jint_Phy05
#$s1  Jint_Phy06
#$s1  Jint_Phy07
#$s1  Jint_Phy08
#$s1  Jint_Phy09
#$s1  Jint_Chl01
#$s1  Jint_Chl02
#$s1  Jint_Chl03
#$s1  Jint_Chl04
#$s1  Jint_Chl05
#$s1  Jint_Chl06
#$s1  Jint_Chl07
#$s1  Jint_Chl08
#$s1  Jint_Chl09
#$s1  Jint_DIC
#$s1  Jint_DOC
#$s1  Jint_POC
#$s1  Jint_PIC
#$s1  Jint_ALK
#$s1  Jint_O2
#$s1  Jint_ZOOC1
#$s1  Jint_ZOOC2
EOF

#1  Jint_zooC
