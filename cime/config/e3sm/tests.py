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
            "ERS.f19_f19.ICLM45",
            "ERS.f19_f19.I1850CLM45CN",
            "ERS.f09_g16.I1850CLM45CN",
            "ERS.f19_f19.I20TRCLM45CN",
            "SMS_Ld1.hcru_hcru.I1850CRUCLM45CN",
            "ERS.f19_g16.I1850CNECACNTBC.clm-eca",
            "ERS.f19_g16.I1850CNECACTCBC.clm-eca",
            "SMS_Ly2_P1x1.1x1_smallvilleIA.ICLM45CNCROP.clm-force_netcdf_pio",
            "ERS_Ld3.f45_f45.ICLM45ED.clm-fates",
            "ERS.f19_g16.I1850CLM45.clm-betr",
            "ERS.f19_g16.I1850CLM45.clm-vst",
            "ERS.f09_g16.I1850CLM45CN.clm-bgcinterface",
            "ERS.ne11_oQU240.I20TRCLM45",
            "ERS.f19_g16.I1850CNRDCTCBC.clm-rd",
            "ERS.f19_g16.I1850GSWCNPECACNTBC.clm-eca_f19_g16_I1850GSWCNPECACNTBC",
            "ERS.f19_g16.I20TRGSWCNPECACNTBC.clm-eca_f19_g16_I20TRGSWCNPECACNTBC",
            "ERS.f19_g16.I1850GSWCNPRDCTCBC.clm-ctc_f19_g16_I1850GSWCNPRDCTCBC",
            "ERS.f19_g16.I20TRGSWCNPRDCTCBC.clm-ctc_f19_g16_I20TRGSWCNPRDCTCBC",
            "ERS.f09_g16.ICLM45BC",
            "SMS.r05_r05.I1850CLM45CN",
            "SMS_Ly2_P1x1.1x1_smallvilleIA.ICLM45CNCROP.clm-lulcc_sville",
            )
        },

    "e3sm_atm_developer" : {
        "tests"   : (
            "ERP_Ln9.ne4_ne4.FC5AV1C-L",
            "SMS_Ln9.ne4_ne4.FC5AV1C-L.cam-outfrq9s",
            "SMS.ne4_ne4.FC5AV1C-L.cam-cosplite",
            "SMS_R_Ld5.ne4_ne4.FSCM5A97",
            "SMS_D_Ln5.ne4_ne4.FC5AV1C-L",
            "SMS_Ln5.ne4pg2_ne4pg2.FC5AV1C-L",
            "SMS_Ln5.ne4pg2_ne4pg2.FC5AV1C-L.cam-thetahy_pg2",
            "SMS_Ln5.ne4pg2_ne4pg2.FC5AV1C-L.cam-thetahy_sl_pg2",
            )
        },

    "e3sm_atm_integration" : {
        "inherit" : ("eam_preqx", "eam_theta"),
        "tests" : (
            "ERP_Ln9.ne4_ne4.F-EAMv1-AQP1",
            "SMS_Ld1.ne4_ne4.F-EAMv1-AQP1.cam-clubb_only",
            "ERP_Ln9.ne4_ne4.F-EAMv1-RCEMIP",
            "PET_Ln5.ne4_ne4.FC5AV1C-L.allactive-mach-pet",
            "PEM_Ln5.ne4_ne4.FC5AV1C-L",
            "SMS_D_Ln5.ne4_ne4.FC5AV1C-L.cam-cosplite_nhtfrq5",
            "ERS_Ld5.ne4_ne4.FC5AV1C-L.cam-rrtmgp",
            "ERS_Ld5.ne4_ne4.FC5AV1C-L.cam-gust_param",
            "REP_Ln5.ne4_ne4.FC5AV1C-L",
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
            "SMS_Ln5.ne30_ne30.FC5AV1C-L.cam-cosplite",
            "SMS.ne30_r05_ne30.F20TRC5-CMIP6",
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
            "SMS.f09_g16_a.IGCLM45_MLI",
            "SMS_P12x2.ne4_oQU240.A_WCYCL1850.allactive-mach_mods",
            "SMS_B.ne4_ne4.F-EAMv1-AQP1.cam-hommexx",
            )
        },

    "homme_integration" : {
        "time"    : "0:45:00",
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
            "ERP_Ld3.ne30_oECv3_ICG.A_WCYCL1850S",
            "SMS.f09_g16_a.MALI",
            "SMS_D_Ln5.conusx4v1_conusx4v1.FC5AV1C-L",
            "SMS.ne30_oECv3.BGCEXP_BCRC_CNPECACNT_1850.clm-bgcexp",
            "SMS.ne30_oECv3.BGCEXP_BCRC_CNPRDCTC_1850.clm-bgcexp",
            "SMS_D_Ld1.T62_oEC60to30v3.DTESTM",
            "SMS_D_Ld1.ne30_r05_oECv3.A_WCYCL1850",
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
            "SMS.ne120_oRRS18v3_ICG.A_WCYCL2000_H01AS.cam-cosplite",
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
            "ERP_Ln9_P96.ne4_ne4.F-MMF1-TEST.cam-crmout",
            "ERP_Ln9_P96.ne4pg2_ne4pg2.F-MMF2-TEST",
            "ERP_Ln9_P96.ne4_ne4.F-MMF2-ECPP-TEST",
            "SMS_D_Ln3_P96.ne4_ne4.F-MMF1-TEST",
            "SMS_D_Ln3_P96.ne4pg2_ne4pg2.F-MMF2-TEST",
            # non-MMF tests with RRTMGP
            "ERP_Ln9.ne4pg2_ne4pg2.FC5AV1C-L.cam-rrtmgp",
            )
        },

    #e3sm tests to mimic production runs
    "e3sm_prod" : {
        "inherit" : "e3sm_atm_prod",
        "tests"   : "SMS_Ld2.ne30_oECv3_ICG.A_WCYCL1850S_CMIP6.allactive-v1cmip6"
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
            "ERS_Ld9.1x1_brazil.ICLM45ED",
            "ERS_D_Ld9.1x1_brazil.ICLM45ED",
            "SMS_D_Lm6.1x1_brazil.ICLM45ED",
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
                 "SMS.ne4_ne4.FC5AV1C-L.cam-preqx_ftype0",
                 "SMS.ne4_ne4.FC5AV1C-L.cam-preqx_ftype1",
                 "SMS.ne4_ne4.FC5AV1C-L.cam-preqx_ftype4",
                 )
    },
    "eam_theta" : {
        "share"    : True,
        "time"     : "02:00:00",
        "tests"    : (
                 "SMS.ne4_ne4.FC5AV1C-L.cam-thetahy_ftype0",
                 "SMS.ne4_ne4.FC5AV1C-L.cam-thetahy_ftype1",
                 "SMS.ne4_ne4.FC5AV1C-L.cam-thetahy_ftype2",
                 "SMS.ne4_ne4.FC5AV1C-L.cam-thetahy_ftype4",
                 "SMS.ne4_ne4.FC5AV1C-L.cam-thetanh_ftype0",
                 "SMS.ne4_ne4.FC5AV1C-L.cam-thetanh_ftype1",
                 "SMS.ne4_ne4.FC5AV1C-L.cam-thetanh_ftype2",
                 "SMS.ne4_ne4.FC5AV1C-L.cam-thetanh_ftype4",
                 "SMS.ne4_ne4.FC5AV1C-L.cam-thetahy_sl",
                 "ERS.ne4_ne4.FC5AV1C-L.cam-thetahy_ftype2",
                 "ERS.ne4_ne4.FC5AV1C-L.cam-thetanh_ftype2",
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
                 "PFS_P7200.ne120_ne120.FC5AV1C-H01A.cam-bench-noio",
                 "PFS_P8640.ne120_ne120.FC5AV1C-H01A.cam-bench-noio",
                 "PFS_P10800.ne120_ne120.FC5AV1C-H01A.cam-bench-noio",
                 "PFS_P14400.ne120_ne120.FC5AV1C-H01A.cam-bench-noio",
                 "PFS_P21600.ne120_ne120.FC5AV1C-H01A.cam-bench-noio",
                 "PFS_P43200.ne120_ne120.FC5AV1C-H01A.cam-bench-noio",
                 "PFS_P86400.ne120_ne120.FC5AV1C-H01A.cam-bench-noio",
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
                 "PFS_P1350.ne30_ne30.FC5AV1C-L.cam-bench-noio",
                 "PFS_P2700.ne30_ne30.FC5AV1C-L.cam-bench-noio",
                 "PFS_P5400.ne30_ne30.FC5AV1C-L.cam-bench-noio",
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

}
