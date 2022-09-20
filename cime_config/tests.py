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

    "e3sm_mosart_developer" : {
        "share" : True,
        "time"  : "0:45:00",
        "tests" : (
            "ERS.r05_r05.RMOSGPCC.mosart-gpcc_1972",
            "ERS.MOS_USRDAT.RMOSGPCC.mosart-mos_usrdat",
            "SMS.MOS_USRDAT.RMOSGPCC.mosart-unstructure",
            )
        },

    "e3sm_mosart_exenoshare": {
        "time"  : "0:45:00",
        "tests" : (
            "ERS.ne30pg2_r05_EC30to60E2r2.GPMPAS-JRA.mosart-rof_ocn_2way",
            )
        },

    "e3sm_land_exeshare" : {
        "share" : True,
        "time"  : "0:45:00",
        "tests" : (
            "ERS.f09_g16.IELMBC",
            "ERS.f19_g16.I1850CNECACNTBC.elm-eca",
            "ERS.f19_g16.I1850CNECACTCBC.elm-eca",
            "ERS.f19_g16.I1850CNRDCTCBC.elm-rd",
            "ERS.f09_g16.I1850GSWCNPRDCTCBC.elm-vstrd",
            "ERS.f19_g16.I1850GSWCNPECACNTBC.elm-eca_f19_g16_I1850GSWCNPECACNTBC",
            "ERS.f19_g16.I20TRGSWCNPECACNTBC.elm-eca_f19_g16_I20TRGSWCNPECACNTBC",
            "ERS.f19_g16.I1850GSWCNPRDCTCBC.elm-ctc_f19_g16_I1850GSWCNPRDCTCBC",
            "ERS.f19_g16.I20TRGSWCNPRDCTCBC.elm-ctc_f19_g16_I20TRGSWCNPRDCTCBC",
            )
        },

    "e3sm_land_exenoshare" : {
        "time"  : "0:45:00",
        "tests" : (
            "ERS_Ld20.f45_f45.IELMFATES.elm-fates",
            "ERS.hcru_hcru.I20TRGSWCNPRDCTCBC.elm-erosion",
            )
        },

    "e3sm_land_developer" : {
        "share" : True,
        "time"  : "0:45:00",
        "inherit" : ("e3sm_mosart_developer", "e3sm_mosart_exenoshare", "e3sm_land_exeshare", "e3sm_land_exenoshare"),
        "tests" : (
            "ERS.f19_f19.IELM",
            "ERS.f19_f19.I1850ELMCN",
            "ERS.f09_g16.I1850ELMCN",
            "ERS.f19_f19.I20TRELMCN",
            "SMS_Ld1.hcru_hcru.I1850CRUELMCN",
            "SMS_Ly2_P1x1.1x1_smallvilleIA.IELMCNCROP.elm-force_netcdf_pio",
            "SMS_Ld20.f45_f45.IELMFATES.elm-fates_rd",
            "SMS_Ld20.f45_f45.IELMFATES.elm-fates_eca",
            "SMS_Ld30.f45_f45.IELMFATES.elm-fates_satphen",
            "ERS.f19_g16.I1850ELM.elm-betr",
            "ERS.f19_g16.I1850ELM.elm-vst",
            "ERS.f09_g16.I1850ELMCN.elm-bgcinterface",
            "ERS.ne11_oQU240.I20TRELM",
            "SMS.r05_r05.I1850ELMCN.elm-qian_1948",
            "SMS_Ly2_P1x1.1x1_smallvilleIA.IELMCNCROP.elm-lulcc_sville",
            "SMS_Ly2_P1x1.1x1_smallvilleIA.IELMCNCROP.elm-per_crop",
            "SMS.r05_r05.IELM.elm-topounit",
            "ERS.ELM_USRDAT.I1850ELM.elm-usrdat",
            "ERS.r05_r05.IELM.elm-V2_ELM_MOSART_features",
            "ERS.f09_f09.IELM.elm-solar_rad",
            "ERS.f09_f09.IELM.elm-lnd_rof_2way",
            "ERS.f09_f09.IELM.elm-koch_snowflake"
            )
        },

    "e3sm_atm_developer" : {
        "inherit" : ("eam_theta_pg2"),
        "tests"   : (
            "ERP_Ln18.ne4_oQU240.F2010",
            "SMS_Ln9.ne4_oQU240.F2010.eam-outfrq9s",
            "SMS.ne4_oQU240.F2010.eam-cosplite",
            "SMS_R_Ld5.ne4_ne4.FSCM-ARM97.eam-scm",
            "SMS_D_Ln5.ne4_oQU240.F2010",
            "SMS_Ln5.ne4pg2_oQU480.F2010",
            "ERS.ne4_oQU240.F2010.eam-hommexx"
            )
        },

    "e3sm_atm_integration" : {
        "inherit" : ("eam_preqx", "eam_theta"),
        "tests" : (
            "ERP_Ln9.ne4_ne4.FAQP",
            "SMS_Ld1.ne4_ne4.FAQP.eam-clubb_only",
            "ERP_Ln9.ne4_ne4.FRCE",
            "PET_Ln5.ne4_oQU240.F2010.allactive-mach-pet",
            "PEM_Ln5.ne4_oQU240.F2010",
            "SMS_D_Ln5.ne4_oQU240.F2010.eam-cosplite_nhtfrq5",
            "SMS_Ln1.ne4_oQU240.F2010.eam-chem_pp",
            "SMS_Ln5.ne30pg2_r05_EC30to60E2r2.BGCEXP_LNDATM_CNPRDCTC_20TR",
            "SMS_Ln5.ne30pg2_r05_EC30to60E2r2.BGCEXP_LNDATM_CNPRDCTC_1850",
            "SMS_D_Ln5.ne4_oQU240.F2010.eam-clubb_sp",
            "ERS_Ld5.ne4_oQU240.F2010.eam-rrtmgp",
            "ERS_Ld5.ne4_oQU240.F2010.eam-rrtmgpxx",
            "REP_Ln5.ne4_oQU240.F2010",
            "SMS_Ld9.ne4pg2_oQU480.F2010.eam-thetahy_sl_pg2_mass",
            "ERP_Ld9.ne4_ne4.FIDEAL",
            )
        },

    #atmopheric tests for extra coverage
    "e3sm_atm_extra_coverage" : {
        "tests" : (
            "SMS_Lm1.ne4_oQU240.F2010",
            "ERS_Ld31.ne4_oQU240.F2010",
            "ERP_Lm3.ne4_oQU240.F2010",
            "SMS_D_Ln5.ne30_oECv3.F2010",
            "ERP_Ld3.ne30_oECv3.F2010.allactive-pioroot1",
            "SMS_Ly1.ne4_oQU240.F2010",
	    "SMS_D_Ln5.ne45pg2_ne45pg2.FAQP",
            "SMS_D_Ln5.ne4_oQU240.F2010.eam-implicit_stress",
            "ERS_Ld5.ne30_oECv3.F2010.eam-implicit_stress",
            )
        },

    #atmopheric tests for hi-res
    "e3sm_atm_hi_res" : {
        "time" : "01:30:00",
        "tests" : "SMS.ne120pg2_r0125_oRRS18to6v3.F2010"
        },

    #atmopheric tests to mimic low res production runs
    "e3sm_atm_prod" : {
        "tests" : (
            "SMS_Ln5.ne30pg2_r05_oECv3.F2010.eam-wcprod",
            "SMS.ne30pg2_r05_oECv3.F20TR.eam-wcprod",
            )
        },

    #atmopheric nbfb tests
    "e3sm_atm_nbfb" : {
        "tests" : (
            "PGN_P1x1.ne4_oQU240.F2010",
            "TSC.ne4_oQU240.F2010-CICE",
            "MVK_PS.ne4_oQU240.F2010",
            )
        },

    "e3sm_ocnice_extra_coverage" : {
        "tests" : (
            "ERS_P480_Ld5.T62_oEC60to30v3wLI.GMPAS-DIB-IAF-ISMF",
            "PEM_P480_Ld5.T62_oEC60to30v3wLI.GMPAS-DIB-IAF-ISMF",
            "SMS.ne30_oECv3_gis.IGELM_MLI.elm-extrasnowlayers",
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
            "SMS_P12x2.ne4_oQU240.WCYCL1850NS.allactive-mach_mods",
            "ERS_Ln9.ne4pg2_ne4pg2.F2010-MMF1.eam-mmf_crmout",
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
        "inherit" : ("e3sm_developer", "e3sm_atm_integration", "e3sm_mmf_integration"),
        "time"    : "03:00:00",
        "tests"   : (
            "ERS.ne11_oQU240.WCYCL1850NS",
            "SMS_D_Ld1.ne30pg2_EC30to60E2r2.WCYCL1850.allactive-wcprod",
            "SMS_D_Ld1.ne30pg2_EC30to60E2r2.WCYCLSSP370.allactive-wcprodssp",
            "ERS_Ld3.ne4_oQU240.F2010",
            #"ERT_Ld31.ne16_g37.B1850C5",#add this line back in with the new correct compset
            "NCK.ne11_oQU240.WCYCL1850NS",
            "PET.f19_g16.X.allactive-mach-pet",
            "PET.f45_g37_rx1.A.allactive-mach-pet",
            "PET_Ln9_PS.ne30pg2_EC30to60E2r2.WCYCL1850.allactive-mach-pet",
            "PEM_Ln9.ne30pg2_EC30to60E2r2.WCYCL1850",
            "ERP_Ld3.ne30pg2_EC30to60E2r2.WCYCL1850.allactive-pioroot1",
            "SMS_D_Ln5.conusx4v1_r05_oECv3.F2010",
            "SMS_Ld2.ne30pg2_r05_EC30to60E2r2.BGCEXP_CNTL_CNPECACNT_1850.elm-bgcexp",
            "SMS_Ld2.ne30pg2_r05_EC30to60E2r2.BGCEXP_CNTL_CNPRDCTC_1850.elm-bgcexp",
            "SMS_D_Ld1.T62_oEC60to30v3.DTESTM",
            "SMS_D_Ld3.T62_oQU120.CMPASO-IAF",
            "SMS_D_Ld1.ne30pg2_r05_EC30to60E2r2.WCYCL1850",
            )
        },

    #e3sm tests for extra coverage
    "e3sm_extra_coverage" : {
        "inherit" : ("e3sm_atm_extra_coverage", "e3sm_ocnice_extra_coverage"),
        "tests"   : (
            "SMS_D_Ln5.enax4v1_enax4v1.F2010-CICE",
            "SMS_D_Ln5.twpx4v1_twpx4v1.F2010-CICE",
            )
        },

    #e3sm tests for hi-res
    "e3sm_hi_res" : {
        "inherit" : "e3sm_atm_hi_res",
        "tests"   : (
            "SMS.ne120pg2_r0125_oRRS18to6v3.WCYCL1950.eam-cosplite",
            "SMS.T62_oRRS30to10v3wLI.GMPAS-IAF",
            )
        },

    #e3sm tests for RRM grids
    "e3sm_rrm" : {
        "tests" : (
            "SMS_D_Ln5.conusx4v1_r05_oECv3.F2010",
            "SMS_D_Ln5.enax4v1_enax4v1.F2010-CICE",
            "SMS_D_Ln5.twpx4v1_twpx4v1.F2010-CICE",
            )
        },

    #e3sm MMF tests for development
    "e3sm_mmf_integration" : {
        "tests" : (
            "ERP_Ln9.ne4pg2_ne4pg2.F2010-MMF1.eam-mmf_fixed_subcycle",
            "ERS_Ln9.ne4pg2_ne4pg2.FRCE-MMF1.eam-cosp_nhtfrq9",
            "SMS_Ln5.ne4_ne4.FSCM-ARM97-MMF1",
            )
        },

    #e3sm tests to mimic production runs
    "e3sm_prod" : {
        "inherit" : "e3sm_atm_prod",
        "tests"   : (
            "SMS_Ld1.ne30pg2_r05_EC30to60E2r2.WCYCL1850.allactive-wcprod",
            "SMS_Ld1.ne30pg2_EC30to60E2r2.WCYCL1850.allactive-wcprod",
            "SMS_Ld1.ne30pg2_EC30to60E2r2.WCYCLSSP370.allactive-wcprodssp",
            "SMS_Ld1.ne30pg2_EC30to60E2r2.WCYCLSSP585.allactive-wcprodssp",
            "SMS_PS.northamericax4v1pg2_WC14to60E2r3.WCYCL1850.allactive-wcprodrrm",
            )
        },

    #e3sm tests to mimic BGC production runs
    "e3sm_bgcprod" : {
        "tests"   :  (
               "SMS_Ld2.ne30_oECv3.BGCEXP_CNTL_CNPRDCTC_1850.allactive-v1bgc_1850",
               "SMS_Ld2.ne30_oECv3.BGCEXP_BCRD_CNPRDCTC_20TR.allactive-v1bgc",
               "SMS_Ld2.ne30_oECv3.BGCEXP_CNTL_CNPECACNT_1850S.allactive-v1bgceca_1850",
               "SMS_Ld2.ne30_oECv3.BGCEXP_BDRD_CNPECACNT_20TRS.allactive-v1bgceca",
               )
        },

    #e3sm performance-benching of production-like runs
    "e3sm_prod_bench" : {
        "tests"   : (
            "PFS.ne30pg2_r05_oECv3.F2010.bench-noio",
            "PFS.ne30pg2_r05_oECv3.F20TR.bench-noio",
            "PFS.ne30pg2_r05_EC30to60E2r2.WCYCL1850.bench-noio",
            "PFS.ne30pg2_EC30to60E2r2.WCYCL1850.bench-noio",
            "PFS_PS.northamericax4v1pg2_WC14to60E2r3.WCYCL1850.bench-noio",
            )
        },

    "fates" : {
        "tests" : (
            "ERS_Ld9.1x1_brazil.IELMFATES",
            "ERS_D_Ld9.1x1_brazil.IELMFATES",
            "SMS_D_Lm6.1x1_brazil.IELMFATES",
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
                 "SMS.ne4_oQU240.F2010.eam-preqx_ftype0",
                 "SMS.ne4_oQU240.F2010.eam-preqx_ftype1",
                 "SMS.ne4_oQU240.F2010.eam-preqx_ftype4",
                 )
    },
    "eam_theta" : {
        "share"    : True,
        "time"     : "02:00:00",
        "tests"    : (
                 "SMS.ne4_oQU240.F2010.eam-thetahy_ftype0",
                 "SMS.ne4_oQU240.F2010.eam-thetahy_ftype1",
                 "SMS.ne4_oQU240.F2010.eam-thetahy_ftype2",
                 "SMS.ne4_oQU240.F2010.eam-thetahy_ftype4",
                 "SMS.ne4_oQU240.F2010.eam-thetanh_ftype0",
                 "SMS.ne4_oQU240.F2010.eam-thetanh_ftype1",
                 "SMS.ne4_oQU240.F2010.eam-thetanh_ftype2",
                 "SMS.ne4_oQU240.F2010.eam-thetanh_ftype4",
                 "SMS.ne4_oQU240.F2010.eam-thetahy_sl",
                 "ERS.ne4_oQU240.F2010.eam-thetahy_ftype2",
                 "ERS.ne4_oQU240.F2010.eam-thetanh_ftype2",
                 )
    },
    "eam_theta_pg2" : {
        "share"    : True,
        "time"     : "02:00:00",
        "tests"    : (
                 "SMS_Ln5.ne4pg2_oQU480.F2010.eam-thetahy_pg2",
                 "SMS_Ln5.ne4pg2_oQU480.F2010.eam-thetahy_sl_pg2",
                 "ERS_Ld3.ne4pg2_oQU480.F2010.eam-thetahy_sl_pg2",
                 "SMS_Ln5.ne4pg2_oQU480.F2010.eam-thetahy_sl_pg2_ftype0",
                 "ERS_Ld3.ne4pg2_oQU480.F2010.eam-thetahy_sl_pg2_ftype0",
                 )
    },
    "e3sm_bench_hires_g" : {
        "share"    : True,
        "time"     : "03:00:00",
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
        "share"    : True,
        "time"     : "03:00:00",
        "tests"    : (
                 "PFS_P7200.ne120pg2_r05_EC30to60E2r2.F2010.eam-bench-noio",
                 "PFS_P8640.ne120pg2_r05_EC30to60E2r2.F2010.eam-bench-noio",
                 "PFS_P10800.ne120pg2_r05_EC30to60E2r2.F2010.eam-bench-noio",
                 "PFS_P14400.ne120pg2_r05_EC30to60E2r2.F2010.eam-bench-noio",
                 "PFS_P21600.ne120pg2_r05_EC30to60E2r2.F2010.eam-bench-noio",
                 "PFS_P43200.ne120pg2_r05_EC30to60E2r2.F2010.eam-bench-noio",
                 "PFS_P86400.ne120pg2_r05_EC30to60E2r2.F2010.eam-bench-noio",
                 )
    },
    "e3sm_bench_hires" : {
        "share"    : True,
        "inherit" : ("e3sm_bench_hires_g", "e3sm_bench_hires_f"),
        "time"    : "03:00:00",
        "tests"   : (
                 "PFS_PS.ne120pg2_r0125_oRRS18to6v3.WCYCL1950.bench-wcycl-hires",
                 "PFS_PM.ne120pg2_r0125_oRRS18to6v3.WCYCL1950.bench-wcycl-hires",
                 "PFS_PL.ne120pg2_r0125_oRRS18to6v3.WCYCL1950.bench-wcycl-hires",
                 )
    },
    "e3sm_bench_lores_g" : {
        "share"    : True,
        "time"     : "03:00:00",
        "tests"    : (
                 "PFS_P320.T62_EC30to60E2r2.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P480.T62_EC30to60E2r2.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P640.T62_EC30to60E2r2.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P960.T62_EC30to60E2r2.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P1280.T62_EC30to60E2r2.GMPAS-IAF.bench-gmpas_noio",
                 )
    },
    "e3sm_bench_lores_f" : {
        "share"    : True,
        "time"     : "03:00:00",
        "tests"    : (
                 "PFS_P1350.ne30pg2_EC30to60E2r2.F2010.eam-bench-noio",
                 "PFS_P2700.ne30pg2_EC30to60E2r2.F2010.eam-bench-noio",
                 "PFS_P5400.ne30pg2_EC30to60E2r2.F2010.eam-bench-noio",
                 )
    },
    "e3sm_bench_lores" : {
        "share"    : True,
        "inherit" : ("e3sm_bench_lores_g", "e3sm_bench_lores_f"),
        "time"    : "03:00:00",
        "tests"   : (
                 "PFS_PS.ne30pg2_EC30to60E2r2.WCYCL1850.bench-wcycl-lores",
                 "PFS_PM.ne30pg2_EC30to60E2r2.WCYCL1850.bench-wcycl-lores",
                 "PFS_PL.ne30pg2_EC30to60E2r2.WCYCL1850.bench-wcycl-lores",
                 )
    },
    "e3sm_bench_all" : {
        "inherit" : ("e3sm_bench_hires", "e3sm_bench_lores"),
        "time"    : "03:00:00",
    },


    "e3sm_gpuacc" : {
        "tests"    : (
                 "SMS_Ld1.T62_oEC60to30v3.CMPASO-NYF",
                 "SMS_Ld1.T62_oEC60to30v3.DTESTM",
                 )
    },

    "e3sm_gpuomp" : {
        "tests"    : (
                 "SMS_Ld1.T62_oEC60to30v3.DTESTM",
                 )
    },

    "e3sm_gpucxx" : {
        "tests"    : (
                 "SMS_Ln9.ne4pg2_ne4pg2.F2010-MMF1",
                 )
    },

    "e3sm_gpuall" : {
        "inherit" : ("e3sm_gpuacc", "e3sm_gpuomp", "e3sm_gpucxx"),
    },

    "eam_nl" : {
        "tests"    : (
            "SBN.ne4_oQU240.F2010.eam-thetadycore",
            "SBN.ne11_ne11.F2010-CICE.eam-thetadycore",
            "SBN.ne16_ne16.F2010-CICE.eam-thetadycore",
            "SBN.ne30_ne30.F2010-CICE.eam-thetadycore",
            "SBN.ne45pg2_r05_oECv3.F2010-CICE.allactive-thetadycore",
            "SBN.ne120_ne120.F2010-CICE.eam-thetadycore",
            "SBN.ne240_ne240.F2010-CICE.eam-thetadycore",
            "SBN.ne512np4_360x720cru_ne512np4.F2010-CICE.allactive-thetadycore",
            "SBN.ne1024np4_360x720cru_ne1024np4.F2010-CICE.allactive-thetadycore",
            "SBN.conusx4v1_conusx4v1.F2010-CICE.eam-thetadycore",
            "SBN.enax4v1_enax4v1.F2010-CICE.eam-thetadycore",
            "SBN.northamericax4v1_r0125_northamericax4v1.F2010-CICE.eam-thetadycore",
            "SBN.antarcticax4v1_r0125_antarcticax4v1.F2010-CICE.allactive-thetadycore",
            "SBN.antarcticax4v1pg2_r0125_antarcticax4v1pg2.F2010-CICE.allactive-thetadycore",
            "SBN.ne4_oQU240.F2010.eam-preqxdycore",
            "SBN.ne11_ne11.F2010-CICE.eam-preqxdycore",
            "SBN.ne16_ne16.F2010-CICE.eam-preqxdycore",
            "SBN.ne30_ne30.F2010-CICE.eam-preqxdycore",
            "SBN.ne45pg2_r05_oECv3.F2010-CICE.allactive-preqxdycore",
            "SBN.ne120_ne120.F2010-CICE.eam-preqxdycore",
            "SBN.ne240_ne240.F2010-CICE.eam-preqxdycore",
            "SBN.ne512np4_360x720cru_ne512np4.F2010-CICE.allactive-preqxdycore",
            "SBN.ne1024np4_360x720cru_ne1024np4.F2010-CICE.allactive-preqxdycore",
            "SBN.conusx4v1_conusx4v1.F2010-CICE.eam-preqxdycore",
            "SBN.enax4v1_enax4v1.F2010-CICE.eam-preqxdycore",
            "SBN.northamericax4v1_r0125_northamericax4v1.F2010-CICE.eam-preqxdycore",
            "SBN.antarcticax4v1_r0125_antarcticax4v1.F2010-CICE.allactive-preqxdycore",
            "SBN.antarcticax4v1pg2_r0125_antarcticax4v1pg2.F2010-CICE.allactive-preqxdycore",
            )
    },

}
