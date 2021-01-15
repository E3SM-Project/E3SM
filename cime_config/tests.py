# Here are the tests belonging to e3sm suites. Format is
# <test>.<grid>.<compset>[.<testmod>]
#
# suite_name : {
#     "inherit" : (suite1, suite2, ...), # Optional. Suites to inherit tests from. Default is None. Tuple, list, or str.
#     "time"    : "HH:MM:SS",            # Optional. Recommended upper-limit on test time.
#     "share"   : True|False,            # Optional. If True, all tests in this suite share a build. Default is False.
#     "tests"   : (test1, test2, ...)    # Optional. The list of tests for this suite. See above for format. Tuple, list, or str.
# }

_TESTS = {

    "e3sm_land_developer" : {
        "time"  : "0:45:00",
        "tests" : (
            "ERS.f19_f19.IELM",
            "ERS.f19_f19.I1850ELMCN",
            "ERS.f09_g16.I1850ELMCN",
            "ERS.f19_f19.I20TRELMCN",
            "SMS_Ld1.hcru_hcru.I1850CRUELMCN",
            "ERS.f19_g16.I1850CNECACNTBC.elm-eca",
            "ERS.f19_g16.I1850CNECACTCBC.elm-eca",
            "SMS_Ly2_P1x1.1x1_smallvilleIA.IELMCNCROP.elm-force_netcdf_pio",
            "ERS_Ld20.f45_f45.IELMED.elm-fates",
            "ERS.f19_g16.I1850ELM.elm-betr",
            "ERS.f19_g16.I1850ELM.elm-vst",
            "ERS.f09_g16.I1850ELMCN.elm-bgcinterface",
            "ERS.ne11_oQU240.I20TRELM",
            "ERS.f19_g16.I1850CNRDCTCBC.elm-rd",
            "ERS.f19_g16.I1850GSWCNPECACNTBC.elm-eca_f19_g16_I1850GSWCNPECACNTBC",
            "ERS.f19_g16.I20TRGSWCNPECACNTBC.elm-eca_f19_g16_I20TRGSWCNPECACNTBC",
            "ERS.f19_g16.I1850GSWCNPRDCTCBC.elm-ctc_f19_g16_I1850GSWCNPRDCTCBC",
            "ERS.f19_g16.I20TRGSWCNPRDCTCBC.elm-ctc_f19_g16_I20TRGSWCNPRDCTCBC",
            "ERS.f09_g16.IELMBC",
            "SMS.r05_r05.I1850ELMCN.elm-qian_1948",
            "SMS_Ly2_P1x1.1x1_smallvilleIA.IELMCNCROP.elm-lulcc_sville",
            "ERS.r05_r05.RMOSGPCC.mosart-gpcc_1972",
            "ERS.MOS_USRDAT.RMOSGPCC.mosart-mos_usrdat",
            "SMS.MOS_USRDAT.RMOSGPCC.mosart-unstructure",
            "ERS.ELM_USRDAT.I1850ELM.elm-usrdat"
            )
        },

    "e3sm_atm_developer" : {
        "inherit" : ("eam_theta_pg2"),
        "tests"   : (
            "ERP_Ln9.ne4_ne4.FC5AV1C-L",
            "SMS_Ln9.ne4_ne4.FC5AV1C-L.eam-outfrq9s",
            "SMS.ne4_ne4.FC5AV1C-L.eam-cosplite",
            "SMS_R_Ld5.ne4_ne4.FSCM5A97.eam-scm",
            "SMS_D_Ln5.ne4_ne4.FC5AV1C-L",
            "SMS_Ln5.ne4pg2_ne4pg2.FC5AV1C-L"
            )
        },

    "e3sm_atm_integration" : {
        "inherit" : ("eam_preqx", "eam_theta"),
        "tests" : (
            "ERP_Ln9.ne4_ne4.F-EAMv1-AQP1",
            "SMS_Ld1.ne4_ne4.F-EAMv1-AQP1.eam-clubb_only",
            "ERP_Ln9.ne4_ne4.F-EAMv1-RCEMIP",
            "PET_Ln5.ne4_ne4.FC5AV1C-L.allactive-mach-pet",
            "PEM_Ln5.ne4_ne4.FC5AV1C-L",
            "SMS_D_Ln5.ne4_ne4.FC5AV1C-L.eam-cosplite_nhtfrq5",
            "SMS_Ln1.ne4_ne4.FC5AV1C-L.eam-chem_pp",
            "ERS_Ld5.ne4_ne4.FC5AV1C-L.eam-rrtmgp",
            "ERS_Ld5.ne4_ne4.FC5AV1C-L.eam-gust_param",
            "REP_Ln5.ne4_ne4.FC5AV1C-L",
            "SMS_Ld9.ne4pg2_ne4pg2.FC5AV1C-04P2.eam-thetahy_sl_pg2_mass",
            )
        },

    #atmopheric tests for extra coverage
    "e3sm_atm_extra_coverage" : {
        "tests" : (
            "SMS_Lm1.ne4_ne4.FC5AV1C-L",
            "ERS_Ld31.ne4_ne4.FC5AV1C-L",
            "ERP_Lm3.ne4_ne4.FC5AV1C-L",
            "SMS_D_Ln5.ne30_ne30.FC5AV1C-L",
            "ERP_Ln7.ne30_ne30.FC5AV1C-L",
            "SMS_Ly1.ne4_ne4.FC5AV1C-L",
	    "SMS_D_Ln5.ne45pg2_ne45pg2.F-EAMv1-AQP1",
            )
        },

    #atmopheric tests for hi-res
    "e3sm_atm_hi_res" : {
        "time" : "01:30:00",
        "tests" : "SMS.ne120_ne120.FC5AV1C-H01A"
        },

    #atmopheric tests to mimic low res production runs
    "e3sm_atm_prod" : {
        "tests" : (
            "SMS_Ln5.ne30pg2_r05_oECv3.F2010SC5-CMIP6.eam-wcprod",
            "SMS.ne30pg2_r05_oECv3.F20TRC5-CMIP6.eam-wcprod",
            )
        },

    #atmopheric nbfb tests
    "e3sm_atm_nbfb" : {
        "tests" : (
            "PGN_P1x1.ne4_ne4.FC5AV1C-L",
            "TSC.ne4_ne4.FC5AV1C-L",
            "MVK_PL.ne4_ne4.FC5AV1C-L",
            )
        },

    "e3sm_ocnice_extra_coverage" : {
        "tests" : (
            "ERS_P480_Ld5.T62_oEC60to30v3wLI.GMPAS-DIB-IAF-ISMF",
            "PEM_P480_Ld5.T62_oEC60to30v3wLI.GMPAS-DIB-IAF-ISMF",
            )
        },

    "e3sm_developer" : {
        "inherit" : ("e3sm_land_developer", "e3sm_atm_developer"),
        "time"    : "0:45:00",
        "tests"   : (
            "ERS.f19_g16_rx1.A",
            "ERS.ne30_g16_rx1.A",
            "SEQ.f19_g16.X",
            "ERIO.ne30_g16_rx1.A",
            "NCK.f19_g16_rx1.A",
            "SMS.ne30_f19_g16_rx1.A",
            "ERS_Ld5.T62_oQU120.CMPASO-NYF",
            "ERS.f09_g16_g.MALISIA",
            "SMS.T62_oQU120_ais20.MPAS_LISIO_TEST",
            "SMS.f09_g16_a.IGELM_MLI",
            "SMS_P12x2.ne4_oQU240.A_WCYCL1850.allactive-mach_mods",
            "SMS_B.ne4_ne4.F-EAMv1-AQP1.eam-hommexx",
            )
        },

    "homme_integration" : {
        "time"    : "0:120:00",
        "tests"   : (
            "HOMME_P24.f19_g16_rx1.A",
            "HOMMEBFB_P24.f19_g16_rx1.A",
            )
        },

    "e3sm_integration" : {
        "inherit" : ("e3sm_developer", "e3sm_atm_integration"),
        "time"    : "03:00:00",
        "tests"   : (
            "ERS.ne11_oQU240.A_WCYCL1850",
            "SMS_D_Ld1.ne30_oECv3_ICG.A_WCYCL1850S_CMIP6.allactive-v1cmip6",
            "ERS_Ln9.ne4_ne4.FC5AV1C-L",
            #"ERT_Ld31.ne16_g37.B1850C5",#add this line back in with the new correct compset
            "NCK.ne11_oQU240.A_WCYCL1850",
            "PET.f19_g16.X.allactive-mach-pet",
            "PET.f45_g37_rx1.A.allactive-mach-pet",
            "PET_Ln9_PS.ne30_oECv3_ICG.A_WCYCL1850S.allactive-mach-pet",
            "PEM_Ln9.ne30_oECv3_ICG.A_WCYCL1850S",
            "ERP_Ld3.ne30_oECv3_ICG.A_WCYCL1850S.allactive-pioroot1",
            "SMS_D_Ln5.conusx4v1_conusx4v1.FC5AV1C-L",
            "SMS.ne30_oECv3.BGCEXP_BCRC_CNPECACNT_1850.elm-bgcexp",
            "SMS.ne30_oECv3.BGCEXP_BCRC_CNPRDCTC_1850.elm-bgcexp",
            "SMS_D_Ld1.T62_oEC60to30v3.DTESTM",
            "SMS_D_Ld1.ne30_r05_oECv3.A_WCYCL1850",
            "ERS_Ln9_P96x1.ne4pg2_ne4pg2.F-MMF1.eam-crmout",
            "ERS_Ln9_P96x1.ne4pg2_ne4pg2.F-MMFXX.eam-genmmf",
            "ERS_Ln9_P96x1.ne4pg2_ne4pg2.F-MMF1-RCEMIP.eam-genmmf",
            )
        },

    #e3sm tests for extra coverage
    "e3sm_extra_coverage" : {
        "inherit" : ("e3sm_atm_extra_coverage", "e3sm_ocnice_extra_coverage"),
        "tests"   : (
            "SMS_D_Ln5.enax4v1_enax4v1.FC5AV1C-L",
            "SMS_D_Ln5.twpx4v1_twpx4v1.FC5AV1C-L",
            )
        },

    #e3sm tests for hi-res
    "e3sm_hi_res" : {
        "inherit" : "e3sm_atm_hi_res",
        "tests"   : (
            "SMS.ne120_oRRS18v3_ICG.A_WCYCL2000_H01AS.eam-cosplite",
            "SMS.T62_oRRS30to10v3wLI.GMPAS-IAF",
            )
        },

    #e3sm tests for RRM grids
    "e3sm_rrm" : {
        "tests" : (
            "SMS_D_Ln5.conusx4v1_conusx4v1.FC5AV1C-L",
            "SMS_D_Ln5.enax4v1_enax4v1.FC5AV1C-L",
            "SMS_D_Ln5.twpx4v1_twpx4v1.FC5AV1C-L",
            )
        },

    #e3sm MMF tests for development
    "e3sm_mmf" : {
        "time" : "02:00:00",
        "tests" : (
            # MMF tests
            "SMS_D_Ln3_P96x1.ne4pg2_ne4pg2.F-MMF1",
            "SMS_D_Ln3_P96x1.ne4pg2_ne4pg2.F-MMF2",
            "ERS_Ln9_P96x1.ne4pg2_ne4pg2.F-MMFXX",
            "ERS_Ln9_P96x1.ne4pg2_ne4pg2.F-MMF1.eam-crmout",
            "ERS_Ln9_P96x1.ne4pg2_ne4pg2.F-MMF2",
            "ERS_Ln9_P96x1.ne4pg2_ne4pg2.F-MMF2-ECPP",
            # non-MMF tests with RRTMGP
            "ERP_Ln9.ne4pg2_ne4pg2.FC5AV1C-L.eam-rrtmgp",
            )
        },

    #e3sm tests to mimic production runs
    "e3sm_prod" : {
        "inherit" : "e3sm_atm_prod",
        "tests"   : "SMS_Ld1.ne30pg2_r05_EC30to60E2r2-1900_ICG.A_WCYCL1850S_CMIP6.allactive-wcprod"
        },

    #e3sm tests to mimic BGC production runs
    "e3sm_bgcprod" : {
        "tests"   :  (
               "SMS_Ld2.ne30_oECv3.BGCEXP_BCRC_CNPRDCTC_1850.allactive-v1bgc_1850",
               "SMS_Ld2.ne30_oECv3.BGCEXP_BCRD_CNPRDCTC_20TR.allactive-v1bgc",
               "SMS_Ld2.ne30_oECv3_ICG.BGCEXP_BCRC_CNPECACNT_1850S.allactive-v1bgceca_1850",
               "SMS_Ld2.ne30_oECv3_ICG.BGCEXP_BDRD_CNPECACNT_20TRS.allactive-v1bgceca",
               )
        },

    "fates" : {
        "tests" : (
            "ERS_Ld9.1x1_brazil.IELMED",
            "ERS_D_Ld9.1x1_brazil.IELMED",
            "SMS_D_Lm6.1x1_brazil.IELMED",
            )
        },


    #atmopheric tests for ftypes with 2 builds only
    #ftype2 is a default and tested in other suites for preqx
    # preqx ftype0
    # preqx ftype1
    # preqx ftype4
    # theta-l hy ftype0
    # theta-l hy ftype1
    # theta-l hy ftype2
    # theta-l hy ftype4
    # theta-l nh ftype0
    # theta-l nh ftype1
    # theta-l nh ftype2
    # theta-l nh ftype4
    # theta-l hy SL
    "eam_preqx" : {
        "share"    : True,
        "time"     : "01:00:00",
        "tests"    : (
                 "SMS.ne4_ne4.FC5AV1C-L.eam-preqx_ftype0",
                 "SMS.ne4_ne4.FC5AV1C-L.eam-preqx_ftype1",
                 "SMS.ne4_ne4.FC5AV1C-L.eam-preqx_ftype4",
                 )
    },
    "eam_theta" : {
        "share"    : True,
        "time"     : "02:00:00",
        "tests"    : (
                 "SMS.ne4_ne4.FC5AV1C-L.eam-thetahy_ftype0",
                 "SMS.ne4_ne4.FC5AV1C-L.eam-thetahy_ftype1",
                 "SMS.ne4_ne4.FC5AV1C-L.eam-thetahy_ftype2",
                 "SMS.ne4_ne4.FC5AV1C-L.eam-thetahy_ftype4",
                 "SMS.ne4_ne4.FC5AV1C-L.eam-thetanh_ftype0",
                 "SMS.ne4_ne4.FC5AV1C-L.eam-thetanh_ftype1",
                 "SMS.ne4_ne4.FC5AV1C-L.eam-thetanh_ftype2",
                 "SMS.ne4_ne4.FC5AV1C-L.eam-thetanh_ftype4",
                 "SMS.ne4_ne4.FC5AV1C-L.eam-thetahy_sl",
                 "ERS.ne4_ne4.FC5AV1C-L.eam-thetahy_ftype2",
                 "ERS.ne4_ne4.FC5AV1C-L.eam-thetanh_ftype2",
                 )
    },
    "eam_theta_pg2" : {
        "share"    : True,
        "time"     : "02:00:00",
        "tests"    : (
                 "SMS_Ln5.ne4pg2_ne4pg2.FC5AV1C-L.eam-thetahy_pg2",
                 "SMS_Ln5.ne4pg2_ne4pg2.FC5AV1C-L.eam-thetahy_sl_pg2",
                 "ERS_Ln5.ne4pg2_ne4pg2.FC5AV1C-L.eam-thetahy_sl_pg2"
                 )
    },
    "e3sm_bench_hires_g" : {
        "share"    : True,
        "time"     : "01:00:00",
        "tests"    : (
                 "PFS_P2560.T62_oRRS18to6v3.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P2792.T62_oRRS18to6v3.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P3072.T62_oRRS18to6v3.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P3200.T62_oRRS18to6v3.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P4096.T62_oRRS18to6v3.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P4800.T62_oRRS18to6v3.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P5120.T62_oRRS18to6v3.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P5200.T62_oRRS18to6v3.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P5584.T62_oRRS18to6v3.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P6400.T62_oRRS18to6v3.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P7200.T62_oRRS18to6v3.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P8192.T62_oRRS18to6v3.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P9600.T62_oRRS18to6v3.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P11168.T62_oRRS18to6v3.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P12000.T62_oRRS18to6v3.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P12800.T62_oRRS18to6v3.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P16000.T62_oRRS18to6v3.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P16384.T62_oRRS18to6v3.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P19200.T62_oRRS18to6v3.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P21600.T62_oRRS18to6v3.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P22400.T62_oRRS18to6v3.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P24000.T62_oRRS18to6v3.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P25600.T62_oRRS18to6v3.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P26000.T62_oRRS18to6v3.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P28000.T62_oRRS18to6v3.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P28800.T62_oRRS18to6v3.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P30000.T62_oRRS18to6v3.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P32000.T62_oRRS18to6v3.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P36000.T62_oRRS18to6v3.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P48000.T62_oRRS18to6v3.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P64000.T62_oRRS18to6v3.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P96000.T62_oRRS18to6v3.GMPAS-IAF.bench-gmpas_noio",
                 )
    },
    "e3sm_bench_hires_f" : {
        "time"     : "01:00:00",
        "tests"    : (
                 "PFS_P7200.ne120_ne120.FC5AV1C-H01A.eam-bench-noio",
                 "PFS_P8640.ne120_ne120.FC5AV1C-H01A.eam-bench-noio",
                 "PFS_P10800.ne120_ne120.FC5AV1C-H01A.eam-bench-noio",
                 "PFS_P14400.ne120_ne120.FC5AV1C-H01A.eam-bench-noio",
                 "PFS_P21600.ne120_ne120.FC5AV1C-H01A.eam-bench-noio",
                 "PFS_P43200.ne120_ne120.FC5AV1C-H01A.eam-bench-noio",
                 "PFS_P86400.ne120_ne120.FC5AV1C-H01A.eam-bench-noio",
                 )
    },
    "e3sm_bench_hires" : {
        "inherit" : ("e3sm_bench_hires_g", "e3sm_bench_hires_f"),
        "time"    : "01:00:00",
        "tests"   : (
                 "PFS_PS.ne120_oRRS18v3_ICG.A_WCYCL1950S_CMIP6_HR.bench-wcycl-hires",
                 "PFS_PM.ne120_oRRS18v3_ICG.A_WCYCL1950S_CMIP6_HR.bench-wcycl-hires",
                 "PFS_PL.ne120_oRRS18v3_ICG.A_WCYCL1950S_CMIP6_HR.bench-wcycl-hires",
                 )
    },
    "e3sm_bench_lores_g" : {
        "share"    : True,
        "time"     : "01:00:00",
        "tests"    : (
                 "PFS_P320.T62_oEC60to30v3.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P480.T62_oEC60to30v3.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P640.T62_oEC60to30v3.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P960.T62_oEC60to30v3.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P1280.T62_oEC60to30v3.GMPAS-IAF.bench-gmpas_noio",
                 )
    },
    "e3sm_bench_lores_f" : {
        "time"     : "01:00:00",
        "tests"    : (
                 "PFS_P1350.ne30_ne30.FC5AV1C-L.eam-bench-noio",
                 "PFS_P2700.ne30_ne30.FC5AV1C-L.eam-bench-noio",
                 "PFS_P5400.ne30_ne30.FC5AV1C-L.eam-bench-noio",
                 )
    },
    "e3sm_bench_lores" : {
        "inherit" : ("e3sm_bench_lores_g", "e3sm_bench_lores_f"),
        "time"    : "01:00:00",
        "tests"   : (
                 "PFS_PS.ne30_oECv3_ICG.A_WCYCL1850S_CMIP6.bench-wcycl-lores",
                 "PFS_PM.ne30_oECv3_ICG.A_WCYCL1850S_CMIP6.bench-wcycl-lores",
                 "PFS_PL.ne30_oECv3_ICG.A_WCYCL1850S_CMIP6.bench-wcycl-lores",
                 )
    },
    "e3sm_bench_all" : {
        "inherit" : ("e3sm_bench_hires", "e3sm_bench_lores"),
        "time"    : "01:00:00",
    },

    "e3sm_scream" : {
        "time"  : "03:00:00",
        "tests" : (
            "SMS_D.ne4pg2_ne4pg2.F2010-SCREAM-HR",
            "SMS_D.ne4pg2_ne4pg2.F2010-SCREAM-LR",
            "ERS.ne4pg2_ne4pg2.F2010-SCREAM-HR",
            "ERS.ne4pg2_ne4pg2.F2010-SCREAM-LR",
            "ERP.ne4pg2_ne4pg2.F2010-SCREAM-HR.eam-double_memleak_tol",
            "ERP.ne4pg2_ne4pg2.F2010-SCREAM-LR.eam-double_memleak_tol",
            "REP.ne4pg2_ne4pg2.F2010-SCREAM-HR",
            "REP.ne4pg2_ne4pg2.F2010-SCREAM-LR",
            "PEM.ne4pg2_ne4pg2.F2010-SCREAM-HR",
            "PEM.ne4pg2_ne4pg2.F2010-SCREAM-LR",
            )
    },

    "e3sm_gpu" : {
        "tests"    : (
                 "SMS_P36x1_Ld1.T62_oEC60to30v3.CMPASO-NYF",
                 "SMS_P36x1_Ld1.T62_oEC60to30v3.DTESTM",
                 )
    },

}
