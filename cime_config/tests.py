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
        "inherit" : ("e3sm_mosart_sediment"),
        "tests" : (
            "ERS.r05_r05.RMOSGPCC.mosart-gpcc_1972",
            "ERS.MOS_USRDAT.RMOSGPCC.mosart-mos_usrdat",
            "SMS.MOS_USRDAT.RMOSGPCC.mosart-unstructure",
            "ERS.r05_r05.RMOSGPCC.mosart-heat",
            )
        },

    "e3sm_mosart_exenoshare": {
        "time"  : "0:45:00",
        "tests" : (
            "ERS.ne30pg2_r05_IcoswISC30E3r5.GPMPAS-JRA.mosart-rof_ocn_2way",
            )
        },

    "e3sm_mosart_sediment" : {
        "time"  : "0:45:00",
        "tests" : (
            "ERS.MOS_USRDAT.RMOSNLDAS.mosart-sediment",
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
            "ERS.f19_g16.I20TRGSWCNPRDCTCBC.elm-ctc_f19_g16_I20TRGSWCNPRDCTCBC",
            "ERS.r05_r05.ICNPRDCTCBC.elm-cbudget",
            "ERS.ELM_USRDAT.I1850CNPRDCTCBC.elm-snowveg_arctic",
            "ERS.ELM_USRDAT.I1850CNPRDCTCBC.elm-usrpft_default_I1850CNPRDCTCBC",
            "ERS.ELM_USRDAT.I1850CNPRDCTCBC.elm-usrpft_codetest_I1850CNPRDCTCBC",
            )
        },

    "e3sm_land_exenoshare" : {
        "time"  : "0:45:00",
        "tests" : (
            "ERS.f19_g16.IERA5ELM",
            "ERS.f19_g16.IERA56HRELM",
            "ERS_Ld20.f45_f45.IELMFATES.elm-fates",
            "ERS.f09_g16.IELMBC.elm-simple_decomp",
            "ERS.hcru_hcru.IELM.elm-multi_inst",
            )
        },

    "e3sm_land_debug" : {
        "time"  : "0:45:00",
        "tests" : (
            "ERS_D.f19_f19.IELM.elm-ic_f19_f19_ielm",
            "ERS_D.f09_g16.I1850ELMCN",
            "ERS_D.ne4pg2_oQU480.I20TRELM.elm-disableDynpftCheck",
            "SMS_Ly2_P1x1_D.1x1_smallvilleIA.IELMCNCROP.elm-lulcc_sville",
            "ERS_D.f19_g16.I1850GSWCNPRDCTCBC.elm-ctc_f19_g16_I1850GSWCNPRDCTCBC",
            "ERS_D.f09_f09.IELM.elm-solar_rad",
            "ERS_D.f09_f09.IELM.elm-koch_snowflake",
            )
        },

    "e3sm_land_developer" : {
        "share" : True,
        "time"  : "0:45:00",
        "inherit" : ("e3sm_mosart_developer", "e3sm_land_exeshare", "e3sm_land_exenoshare", "e3sm_land_debug", "fates_elm_developer"),
        "tests" : (
            "ERS.f19_f19.I1850ELMCN",
            "ERS.f19_f19.I20TRELMCN",
            "SMS_Ld1.hcru_hcru.I1850CRUELMCN",
            "SMS_Ly2_P1x1.1x1_smallvilleIA.IELMCNCROP.elm-force_netcdf_pio",
            "ERS.f19_g16.I1850ELM.elm-betr",
            "ERS.f19_g16.I1850ELM.elm-vst",
            "ERS.f09_g16.I1850ELMCN.elm-bgcinterface",
            "SMS.r05_r05.I1850ELMCN.elm-qian_1948",
            "SMS_Ly2_P1x1.1x1_smallvilleIA.IELMCNCROP.elm-per_crop",
            "SMS_Ly2_P1x1.1x1_smallvilleIA.IELMCNCROP.elm-fan",
            "SMS.r05_r05.IELM.elm-topounit",
            "ERS.ELM_USRDAT.I1850ELM.elm-usrdat",
            "ERS.r05_r05.IELM.elm-lnd_rof_2way",
            "ERS.r05_r05.IELM.elm-V2_ELM_MOSART_features",
            "ERS.ELM_USRDAT.IELM.elm-surface_water_dynamics"
            )
        },

    "e3sm_atm_developer" : {
        "inherit" : ("eam_theta_pg2"),
        "tests"   : (
            "ERP_Ld3.ne4pg2_oQU480.F2010",
            "SMS_Ln9.ne4pg2_oQU480.F2010.eam-outfrq9s",
            "SMS.ne4pg2_oQU480.F2010.eam-cosplite",
            "SMS_R_Ld5.ne4_ne4.FSCM-ARM97.eam-scm",
            "SMS_D_Ln5.ne4pg2_oQU480.F2010",
            "SMS_Ln5.ne4pg2_oQU480.F2010",
            "ERS_D.ne4pg2_oQU480.F2010.eam-hommexx",
            "SMS_Ln9_P24x1.ne4_ne4.FDPSCREAM-ARM97",
            )
        },

    "e3sm_ice_developer" : {
        "tests"   : (
            "SMS_D_Ld1.TL319_IcoswISC30E3r5.DTESTM-JRA1p5.mpassi-jra_1958",
            "ERS_Ld5.T62_oQU240.DTESTM",
            "PEM_Ln5.T62_oQU240wLI.DTESTM",
            "PET_Ln5.T62_oQU240.DTESTM",
            )
        },

    "e3sm_cryo_developer" : {
        "tests"   : (
            "SMS_D_Ld1.TL319_IcoswISC30E3r5.GMPAS-JRA1p5-DIB-PISMF.mpaso-jra_1958",
            "ERS_Ld5.T62_oQU240wLI.GMPAS-DIB-IAF-PISMF",
            "PEM_Ln5.T62_oQU240wLI.GMPAS-DIB-IAF-PISMF",
            "PET_Ln5.T62_oQU240wLI.GMPAS-DIB-IAF-PISMF",
            "ERS_Ld5.T62_oQU240wLI.GMPAS-DIB-IAF-DISMF",
            "PEM_Ln5.T62_oQU240wLI.GMPAS-DIB-IAF-DISMF",
            "PET_Ln5.T62_oQU240wLI.GMPAS-DIB-IAF-DISMF",
            )
        },

    "e3sm_landice_developer" : {
        "tests"   : (
            "SMS.ne30pg2_r05_IcoswISC30E3r5_gis20.IGELM_MLI.mali-gis20km",
            "ERS.ne30pg2_r05_IcoswISC30E3r5_gis20.IGELM_MLI.mali-gis20km",
            "SMS.ne30pg2_r05_IcoswISC30E3r5_gis20.BGWCYCL1850.allactive-gis20km",
            "SMS.ne30_oECv3_gis.IGELM_MLI.elm-extrasnowlayers",
            )
        },

    "eam_condidiag" : {
        "tests"   : (
            "SMS_D_Ln5.ne4pg2_oQU480.F2010",
            "SMS_D_Ln5.ne4pg2_oQU480.F2010.eam-condidiag_dcape",
            "ERP_Ld3.ne4pg2_oQU480.F2010.eam-condidiag_dcape",
            "ERP_Ld3.ne4pg2_oQU480.F2010.eam-condidiag_rhi",
            )
        },

    "e3sm_zm_developer" : {
        "tests"   : (
            "ERP.ne4pg2_oQU480.F2010.eam-zm_enhancements",
            "REP_Ln5.ne4pg2_oQU480.F2010.eam-zm_enhancements",
            "PET.ne4pg2_oQU480.F2010.eam-zm_enhancements",
            "PEM_Ln18.ne4pg2_oQU480.F2010.eam-zm_enhancements",
            "SMS_Ln5.ne30pg2_EC30to60E2r2.F2010.eam-zm_enhancements",
            "SMS_D_Ln5.ne4_oQU240.F2010.eam-zm_enhancements",
            "SMS_Ln5.ne4pg2_oQU480.F2010.eam-zm_enhancements",
            "ERS.ne4pg2_oQU480.F2010.eam-zm_enhancements"
            )
        },

    "e3sm_atm_stealth" : {
        "tests"   : (
            "ERP_Ln18.ne4_oQU240.F2010.eam-cflx_cpl_2",
            "SMS_D_Ln5.ne4_oQU240.F2010.eam-cflx_cpl_2",
            "ERS.ne4pg2_oQU480.F2010.eam-p3",
            "SMS_D_Ln5.ne4pg2_oQU480.F2010.eam-p3",
            )
        },

        "e3sm_p3_developer" : {
        "tests"   : (
            "ERP.ne4pg2_oQU480.F2010.eam-p3",
            "REP_Ln5.ne4pg2_oQU480.F2010.eam-p3",
            "PET.ne4pg2_oQU480.F2010.eam-p3",
            "PEM_Ln18.ne4pg2_oQU480.F2010.eam-p3",
            "SMS_Ln5.ne30pg2_EC30to60E2r2.F2010.eam-p3",
            "SMS_D_Ln5.ne4pg2_oQU480.F2010.eam-p3",
            "SMS_Ln5.ne4pg2_oQU480.F2010.eam-p3",
            "ERS.ne4pg2_oQU480.F2010.eam-p3"
            )
        },

    "e3sm_atm_integration" : {
        "inherit" : ("eam_preqx", "eam_theta"),
        "tests" : (
            "ERP_Ln9.ne4pg2_ne4pg2.FAQP",
            "SMS_Ld1.ne4pg2_ne4pg2.FAQP.eam-clubb_only",
            "ERP_Ln9.ne4pg2_ne4pg2.FRCE",
            "PET_Ln5.ne4pg2_oQU480.F2010.allactive-mach-pet",
            "PEM_Ln5.ne4pg2_oQU480.F2010",
            "SMS_D_Ln5.ne4pg2_oQU480.F2010.eam-cosplite_nhtfrq5",
            "SMS_Ln1.ne4pg2_oQU480.F2010.eam-chem_pp",
            "SMS_Ln5.ne30pg2_r05_IcoswISC30E3r5.BGCEXP_LNDATM_CNPRDCTC_20TR",
            "SMS_Ln5.ne30pg2_r05_IcoswISC30E3r5.BGCEXP_LNDATM_CNPRDCTC_1850",
            "SMS_D_Ln5.ne4pg2_oQU480.F2010.eam-clubb_sp",
            "ERS_Ld5.ne4pg2_oQU480.F2010.eam-rrtmgp",
            "ERS_Ld5.ne4pg2_oQU480.F2010.eam-rrtmgpxx",
            "REP_Ln5.ne4pg2_oQU480.F2010",
            "SMS_Ld3.ne4pg2_oQU480.F2010.eam-thetahy_sl_pg2_mass",
            "ERP_Ld3.ne4pg2_ne4pg2.FIDEAL.allactive-pioroot1",
            "ERS_Ld5.ne4pg2_oQU480.F2010.eam-sathist_F2010",
            )
        },

    #atmopheric tests for extra coverage
    "e3sm_atm_extra_coverage" : {
        "tests" : (
            "SMS_Lm1.ne4pg2_oQU480.F2010",
            "ERS_Ld31.ne4pg2_oQU480.F2010",
            "ERP_Lm3.ne4pg2_oQU480.F2010",
            "SMS_D_Ln5.ne30pg2_r05_IcoswISC30E3r5.F2010",
            "ERP_Ld3.ne30pg2_r05_IcoswISC30E3r5.F2010.allactive-pioroot1",
            "SMS_Ly1.ne4pg2_oQU480.F2010",
	    "SMS_D_Ln5.ne45pg2_ne45pg2.FAQP",
            "SMS_D_Ln5.ne4pg2_oQU480.F2010.eam-implicit_stress",
            "ERS_Ld5.ne30pg2_r05_IcoswISC30E3r5.F2010.eam-implicit_stress",
            "ERP_Ld3.ne4pg2_oQU480.F2010.eam-condidiag_dcape",
            "ERP_Ld3.ne4pg2_oQU480.F2010.eam-condidiag_rhi",
            )
        },

    #atmopheric tests for hi-res
    "e3sm_atm_hi_res" : {
        "time" : "01:30:00",
        "tests" : "SMS.ne120pg2_r025_RRSwISC6to18E3r5.F2010"
        },

    #atmopheric tests to mimic low res production runs
    "e3sm_atm_prod" : {
        "tests" : (
            "SMS_Ln5.ne30pg2_r05_IcoswISC30E3r5.F2010.eam-wcprod_F2010",
            "SMS_Ld1.ne30pg2_r05_IcoswISC30E3r5.F20TR.eam-wcprod_F20TR",
            )
        },

    #atmopheric nbfb tests
    "e3sm_atm_nbfb" : {
        "tests" : (
            "PGN_P1x1.ne4pg2_oQU480.F2010",
            "TSC_PS.ne4pg2_oQU480.F2010",
            "MVK_PS.ne4pg2_oQU480.F2010",
            )
        },

    #ocean non bit-for-bit test
    "e3sm_ocn_nbfb": {
        "tests": (
            "MVKO_PS.T62_oQU240.GMPAS-NYF",
            )
        },

    "e3sm_nbfb": {
        "inherit": ("e3sm_atm_nbfb", "e3sm_ocn_nbfb")
    },

    "e3sm_ocnice_stealth_features" : {
        "tests" : (
            "SMS_D_Ld1.T62_oQU240wLI.GMPAS-IAF-PISMF.mpaso-impl_top_drag",
            "SMS_D_Ld1.T62_oQU240.GMPAS-IAF.mpaso-harmonic_mean_drag",
            "SMS_D_Ld1.T62_oQU240.GMPAS-IAF.mpaso-upwind_advection",
            "ERS_Ld5_D.T62_oQU240.GMPAS-IAF.mpaso-conservation_check",
            )
        },

    "e3sm_ocnice_extra_coverage" : {
        "inherit" : ("e3sm_ocnice_stealth_features"),
        "tests" : (
            "ERS_P480_Ld5.TL319_IcoswISC30E3r5.GMPAS-JRA1p5-DIB-PISMF.mpaso-jra_1958",
            "PEM_P480_Ld5.TL319_IcoswISC30E3r5.GMPAS-JRA1p5-DIB-PISMF.mpaso-jra_1958",
            "SMS_P480_Ld5.TL319_IcoswISC30E3r5.GMPAS-JRA1p5-DIB-PISMF-TMIX.mpaso-jra_1958",
            "PET_P480_Ld2.TL319_IcoswISC30E3r5.GMPAS-JRA1p5-DIB-PISMF-DSGR.mpaso-jra_1958",
            )
        },

    "e3sm_atm_dustemis" : {
        "time"  : "1:45:00",
        "tests"   : (
            "ERP.ne4pg2_oQU480.F2010.eam-v3atm_dustemis",
            "REP.ne4pg2_oQU480.F2010.eam-v3atm_dustemis",
            "SMS.ne30pg2_IcoswISC30E3r5.F2010.eam-v3atm_dustemis",
            "SMS_D_Ln5.ne4pg2_oQU480.F2010.eam-v3atm_dustemis",
            "PET_Ln5.ne4pg2_oQU480.F2010.eam-v3atm_dustemis",
            "PEM_Ln5.ne4pg2_oQU480.F2010.eam-v3atm_dustemis",
            "ERS_D.ne4pg2_oQU480.F2010.eam-v3atm_dustemis"
            )
        },

    "e3sm_developer" : {
        "inherit" : ("e3sm_land_developer", "e3sm_atm_developer", "e3sm_ice_developer", "e3sm_cryo_developer"),
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
            "ERS_Ld5.TL319_oQU240wLI_ais8to30.MPAS_LISIO_JRA1p5.mpaso-ocn_glcshelf",
            "SMS_P12x2.ne4pg2_oQU480.WCYCL1850NS.allactive-mach_mods",
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
        "inherit" : ("e3sm_developer", "e3sm_atm_integration", "e3sm_mmf_integration", "e3sm_rrm"),
        "time"    : "03:00:00",
        "tests"   : (
            "ERS.ne4pg2_oQU480.WCYCL1850NS",
            "ERS_Vmoab.ne4pg2_oQU480.WCYCL1850NS",
            "SMS_D_Ld1.ne30pg2_r05_IcoswISC30E3r5.WCYCL1850.allactive-wcprod",
            "SMS_D_Ld1.ne30pg2_r05_IcoswISC30E3r5.WCYCLSSP370.allactive-wcprodssp",
            "ERS_Ld3.ne4pg2_oQU480.F2010",
            "NCK.ne4pg2_oQU480.WCYCL1850NS",
            "PET.f19_g16.X.allactive-mach-pet",
            "PET.f45_g37_rx1.A.allactive-mach-pet",
            "PET_Ln9_PS.ne30pg2_r05_IcoswISC30E3r5.WCYCL1850.allactive-mach-pet",
            "PEM_Ln9.ne30pg2_r05_IcoswISC30E3r5.WCYCL1850",
            "ERP_Ld3.ne30pg2_r05_IcoswISC30E3r5.WCYCL1850.allactive-pioroot1",
            "SMS_Ld2.ne30pg2_r05_IcoswISC30E3r5.BGCEXP_CNTL_CNPECACNT_1850.elm-bgcexp",
            "SMS_Ld2.ne30pg2_r05_IcoswISC30E3r5.BGCEXP_CNTL_CNPRDCTC_1850.elm-bgcexp",
            "SMS_D_Ld3.T62_oQU120.CMPASO-IAF",
            "SMS_Ln5.ne30pg2_ne30pg2.F2010-SCREAM-LR-DYAMOND2",
            "ERS_Ld3.ne30pg2_r05_IcoswISC30E3r5.WCYCL1850.allactive-nlmaps",
            "SMS_D_Ld1.ne30pg2_r05_IcoswISC30E3r5.CRYO1850-DISMF",
            "ERS.hcru_hcru.I20TRGSWCNPRDCTCBC.elm-erosion",
            "ERS.ne30pg2_r05_IcoswISC30E3r5.GPMPAS-JRA.mosart-rof_ocn_2way",
            )
        },

    #e3sm tests for extra coverage
    "e3sm_extra_coverage" : {
        "inherit" : ("e3sm_atm_extra_coverage", "e3sm_ocnice_extra_coverage"),
        "tests"   : (
            "SMS_D_Ln3.TL319_EC30to60E2r2_wQU225EC30to60E2r2.GMPAS-JRA1p5-WW3.ww3-jra_1958",
            )
        },

    #e3sm tests for hi-res
    "e3sm_hi_res" : {
        "inherit" : "e3sm_atm_hi_res",
        "tests"   : (
            "SMS_Ld3.ne120pg2_r025_RRSwISC6to18E3r5.WCYCL1850NS.eam-cosplite",
            "SMS.T62_SOwISC12to60E2r4.GMPAS-IAF",
            )
        },

    #e3sm tests for RRM grids
    "e3sm_rrm" : {
        "tests" : (
            "SMS_D_Ln5.conusx4v1pg2_r05_IcoswISC30E3r5.F2010",
            )
        },

    #e3sm MMF tests for development
    "e3sm_mmf_integration" : {
        "tests" : (
            "ERP_Ln9.ne4pg2_oQU480.WCYCL20TRNS-MMF1.allactive-mmf_fixed_subcycle",
            "ERS_Ln9.ne4pg2_ne4pg2.FRCE-MMF1.eam-cosp_nhtfrq9",
            "SMS_Ln5.ne4_ne4.FSCM-ARM97-MMF1",
            "SMS_Ln3.ne4pg2_oQU480.F2010-MMF2",
            )
        },

    #e3sm tests to mimic production runs
    "e3sm_prod" : {
        "inherit" : "e3sm_atm_prod",
        "tests"   : (
            "SMS_Ld1.ne30pg2_r05_IcoswISC30E3r5.WCYCL1850-1pctCO2.allactive-wcprod_1850_1pctCO2",
            "SMS_Ld1.ne30pg2_r05_IcoswISC30E3r5.WCYCL1850-4xCO2.allactive-wcprod_1850_4xCO2",
            "SMS_Ld1.ne30pg2_r05_IcoswISC30E3r5.WCYCL1850.allactive-wcprod_1850",
            "SMS_Ld1.ne30pg2_r05_IcoswISC30E3r5.WCYCLSSP370.allactive-wcprodssp",
            "SMS_Ld1.ne30pg2_r05_IcoswISC30E3r5.WCYCLSSP585.allactive-wcprodssp",
            "SMS_Ld1_PS.northamericax4v1pg2_WC14to60E2r3.WCYCL1850.allactive-wcprodrrm_1850",
            "SMS_D_Ld1.ne30pg2_r05_IcoswISC30E3r5.CRYO1850",
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
            "PFS.ne30pg2_r05_IcoswISC30E3r5.F2010.bench-noio",
            "PFS.ne30pg2_r05_IcoswISC30E3r5.WCYCL1850.bench-noio",
            )
        },

    #fates debug tests included in e3sm land developer test runs
    "fates_elm_debug" : {
        "tests" : (
            "SMS_D_Ld20.f45_f45.IELMFATES.elm-fates_rd",
            "ERS_D_Ld15.f45_g37.IELMFATES.elm-fates_cold_treedamage",
            )
        },

    #fates non-debug tests included in e3sm land developer test runs
    "fates_elm_developer" : {
        "inherit" : ("fates_elm_debug"),
        "tests" : (
            "ERS_Ld30.f45_f45.IELMFATES.elm-fates_satphen",
            "ERS_Ld30.f45_g37.IELMFATES.elm-fates_cold_sizeagemort",
            "SMS_Ld20.f45_f45.IELMFATES.elm-fates_eca",
            "SMS_Ld5_PS.f19_g16.IELMFATES.elm-fates_cold",
            )
        },

    #fates long duration tests runs
    "fates_long_tests" : {
        "time"    : "00:40:00",
        "tests"   : (
            "SMS_D_Lm6.f45_g37.IELMFATES.elm-fates_cold",
            "ERS_D_Lm13.ne4pg2_ne4pg2.IELMFATES.elm-fates_long",
            "ERS_Lm25.ne4pg2_ne4pg2.IELMFATES.elm-fates_cold_nocomp",
            )
        },

    #fates testmod coverage
    "fates_landuse" : {
        "time"    : "00:40:00",
        "tests"   : (
            "ERS_Ld60.f45_g37.IELMFATES.elm-fates_cold_logging",
            "ERS_D_Ld30.f45_g37.IELMFATES.elm-fates_cold_landuse",
            "ERS_D_Ld30.f45_g37.IELMFATES.elm-fates_cold_luh2",
            "ERS_D_Ld30.f45_g37.IELMFATES.elm-fates_cold_luh2harvestarea",
            "ERS_D_Ld30.f45_g37.IELMFATES.elm-fates_cold_luh2harvestmass",
            )
        },
    #fates testmod coverage
    "fates" : {
        "inherit" : ("fates_long_tests", "fates_elm_developer", "fates_landuse"),
        "tests" : (
            "ERP_Ld15.ne4pg2_ne4pg2.IELMFATES.elm-fates_cold_allvars",
            "ERP_Ld3.f09_g16.IELMFATES.elm-fates_cold",
            "ERP_D_Ld3.f19_g16.IELMFATES.elm-fates_cold",
            "ERS_D_Ld3_PS.f09_g16.IELMFATES.elm-fates_cold",
            "ERS_D_Ld5.f45_g37.IELMFATES.elm-fates_cold",
            "ERS_Ld30.f45_g37.IELMFATES.elm-fates_satphen",
            "ERS_Ld30.f45_g37.IELMFATES.elm-fates_cold_fixedbiogeo",
            "ERS_Ld30.f45_g37.IELMFATES.elm-fates_cold_nocomp",
            "ERS_Ld30.f45_g37.IELMFATES.elm-fates_cold_nocomp_fixedbiogeo",
            "ERS_Ld60.f45_g37.IELMFATES.elm-fates",
            "ERS_Ld60.f45_g37.IELMFATES.elm-fates_cold_nofire",
            "ERS_Ld60.f45_g37.IELMFATES.elm-fates_cold_st3",
            "ERS_Ld60.f45_g37.IELMFATES.elm-fates_cold_pphys",
            "SMS_D_Ld15.f45_g37.IELMFATES.elm-fates_cold_twostream",
            )
        },

    #e3sm v3atm related tests for development
    "e3sm_v3atm_integration" : {
        "tests" : (
            "ERP_Ld3.ne4pg2_oQU480.F2010",
            "ERS_Ld3.ne4pg2_oQU480.F20TR",
            "SMS_Ld1.ne30pg2_EC30to60E2r2.WCYCL1850.allactive-wcprod",
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
                 "SMS.ne4pg2_oQU480.F2010.eam-preqx_ftype0",
                 "SMS.ne4pg2_oQU480.F2010.eam-preqx_ftype1",
                 "SMS.ne4pg2_oQU480.F2010.eam-preqx_ftype4",
                 )
    },
    "eam_theta" : {
        "share"    : True,
        "time"     : "02:00:00",
        "tests"    : (
                 "SMS.ne4pg2_oQU480.F2010.eam-thetahy_ftype0",
                 "SMS.ne4pg2_oQU480.F2010.eam-thetahy_ftype1",
                 "SMS.ne4pg2_oQU480.F2010.eam-thetahy_ftype2",
                 "SMS.ne4pg2_oQU480.F2010.eam-thetahy_ftype2_energy",
                 "SMS.ne4pg2_oQU480.F2010.eam-thetahy_ftype4",
                 "SMS.ne4pg2_oQU480.F2010.eam-thetanh_ftype0",
                 "SMS.ne4pg2_oQU480.F2010.eam-thetanh_ftype1",
                 "SMS.ne4pg2_oQU480.F2010.eam-thetanh_ftype2",
                 "SMS.ne4pg2_oQU480.F2010.eam-thetanh_ftype4",
                 "SMS.ne4pg2_oQU480.F2010.eam-thetahy_sl",
                 "ERS.ne4pg2_oQU480.F2010.eam-thetahy_ftype2",
                 "ERS.ne4pg2_oQU480.F2010.eam-thetanh_ftype2",
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
                 "PFS_P2560.T62_RRSwISC6to18E3r5.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P2792.T62_RRSwISC6to18E3r5.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P3072.T62_RRSwISC6to18E3r5.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P3200.T62_RRSwISC6to18E3r5.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P4096.T62_RRSwISC6to18E3r5.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P4800.T62_RRSwISC6to18E3r5.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P5120.T62_RRSwISC6to18E3r5.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P5200.T62_RRSwISC6to18E3r5.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P5584.T62_RRSwISC6to18E3r5.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P6400.T62_RRSwISC6to18E3r5.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P7200.T62_RRSwISC6to18E3r5.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P8192.T62_RRSwISC6to18E3r5.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P9600.T62_RRSwISC6to18E3r5.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P11168.T62_RRSwISC6to18E3r5.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P12000.T62_RRSwISC6to18E3r5.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P12800.T62_RRSwISC6to18E3r5.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P16000.T62_RRSwISC6to18E3r5.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P16384.T62_RRSwISC6to18E3r5.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P19200.T62_RRSwISC6to18E3r5.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P21600.T62_RRSwISC6to18E3r5.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P22400.T62_RRSwISC6to18E3r5.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P24000.T62_RRSwISC6to18E3r5.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P25600.T62_RRSwISC6to18E3r5.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P26000.T62_RRSwISC6to18E3r5.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P28000.T62_RRSwISC6to18E3r5.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P28800.T62_RRSwISC6to18E3r5.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P30000.T62_RRSwISC6to18E3r5.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P32000.T62_RRSwISC6to18E3r5.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P36000.T62_RRSwISC6to18E3r5.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P48000.T62_RRSwISC6to18E3r5.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P64000.T62_RRSwISC6to18E3r5.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P96000.T62_RRSwISC6to18E3r5.GMPAS-IAF.bench-gmpas_noio",
                 )
    },
    "e3sm_bench_hires_f" : {
        "share"    : True,
        "time"     : "03:00:00",
        "tests"    : (
                 "PFS_P7200.ne120pg2_r025_RRSwISC6to18E3r5.F2010.eam-bench-noio",
                 "PFS_P8640.ne120pg2_r025_RRSwISC6to18E3r5.F2010.eam-bench-noio",
                 "PFS_P10800.ne120pg2_r025_RRSwISC6to18E3r5.F2010.eam-bench-noio",
                 "PFS_P14400.ne120pg2_r025_RRSwISC6to18E3r5.F2010.eam-bench-noio",
                 "PFS_P21600.ne120pg2_r025_RRSwISC6to18E3r5.F2010.eam-bench-noio",
                 "PFS_P43200.ne120pg2_r025_RRSwISC6to18E3r5.F2010.eam-bench-noio",
                 "PFS_P86400.ne120pg2_r025_RRSwISC6to18E3r5.F2010.eam-bench-noio",
                 )
    },
    "e3sm_bench_hires" : {
        "share"    : True,
        "inherit" : ("e3sm_bench_hires_g", "e3sm_bench_hires_f"),
        "time"    : "03:00:00",
        "tests"   : (
                 "PFS_PS.ne120pg2_r025_RRSwISC6to18E3r5.WCYCL1850NS.bench-wcycl-hires",
                 "PFS_PM.ne120pg2_r025_RRSwISC6to18E3r5.WCYCL1850NS.bench-wcycl-hires",
                 "PFS_PL.ne120pg2_r025_RRSwISC6to18E3r5.WCYCL1850NS.bench-wcycl-hires",
                 )
    },
    "e3sm_bench_hires_tiny" : {
        "time"  : "03:00:00",
        "tests" : (
                "PFS_P16384.T62_RRSwISC6to18E3r5.GMPAS-IAF.bench-gmpas_noio",
                "PFS_P16384.ne120pg2_r025_RRSwISC6to18E3r5.F2010.eam-bench-noio",
                "PFS.ne120pg2_r025_RRSwISC6to18E3r5.WCYCL1850NS.bench-wcycl-hires",
                )
    },
    "e3sm_bench_lores_g" : {
        "share"    : True,
        "time"     : "03:00:00",
        "tests"    : (
                 "PFS_P320.T62_IcoswISC30E3r5.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P480.T62_IcoswISC30E3r5.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P640.T62_IcoswISC30E3r5.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P960.T62_IcoswISC30E3r5.GMPAS-IAF.bench-gmpas_noio",
                 "PFS_P1280.T62_IcoswISC30E3r5.GMPAS-IAF.bench-gmpas_noio",
                 )
    },
    "e3sm_bench_lores_f" : {
        "share"    : True,
        "time"     : "03:00:00",
        "tests"    : (
                 "PFS_P1350.ne30pg2_r05_IcoswISC30E3r5.F2010.eam-bench-noio",
                 "PFS_P2700.ne30pg2_r05_IcoswISC30E3r5.F2010.eam-bench-noio",
                 "PFS_P5400.ne30pg2_r05_IcoswISC30E3r5.F2010.eam-bench-noio",
                 )
    },
    "e3sm_bench_lores" : {
        "share"    : True,
        "inherit" : ("e3sm_bench_lores_g", "e3sm_bench_lores_f"),
        "time"    : "03:00:00",
        "tests"   : (
                 "PFS_PS.ne30pg2_r05_IcoswISC30E3r5.WCYCL1850.bench-wcycl-lores",
                 "PFS_PM.ne30pg2_r05_IcoswISC30E3r5.WCYCL1850.bench-wcycl-lores",
                 "PFS_PL.ne30pg2_r05_IcoswISC30E3r5.WCYCL1850.bench-wcycl-lores",
                 )
    },
    "e3sm_bench_all" : {
        "inherit" : ("e3sm_bench_hires", "e3sm_bench_lores"),
        "time"    : "03:00:00",
    },

    "e3sm_scream" : {
        "time"    : "03:00:00",
        "inherit" : ("e3sm_scream_v0"),
    },

    "e3sm_scream_v0" : {
        "time"    : "03:00:00",
        "inherit" : ("e3sm_scream_v0_lowres"),
    },

    "e3sm_scream_v0_lowres" : {
        "time"  : "03:00:00",
        "tests" : (
            "SMS_D.ne4pg2_ne4pg2.F2010-SCREAM-HR",
            "SMS_D.ne4pg2_ne4pg2.F2010-SCREAM-LR",
            "ERP.ne4pg2_ne4pg2.F2010-SCREAM-HR.eam-double_memleak_tol",
            "ERP.ne4pg2_ne4pg2.F2010-SCREAM-LR.eam-double_memleak_tol",
            "ERS_R_Ln10.ne4_ne4.FDPSCREAM-ARM97",
            )
    },

    "e3sm_scream_v1" : {
        "time"    : "03:00:00",
        "inherit" : ("e3sm_scream_v1_lowres", "e3sm_scream_v1_medres", "e3sm_scream_v1_mpassi"),
    },


    "e3sm_scream_v1_lowres" : {
        "time"  : "01:00:00",
        "inherit" : ("e3sm_scream_mam4xx_v1_lowres"),
        "tests" : (
            "ERP_D_Lh4.ne4_ne4.F2010-SCREAMv1.scream-output-preset-1",
            "ERS_Ln9.ne4_ne4.F2000-SCREAMv1-AQP1.scream-output-preset-2",
            "SMS_D_Ln9.ne4_ne4.F2010-SCREAMv1-noAero.scream-output-preset-3",
            "ERP_Ln22.ne4pg2_ne4pg2.F2010-SCREAMv1.scream-output-preset-4",
            "ERS_D_Ln22.ne4pg2_ne4pg2.F2010-SCREAMv1.scream-rad_frequency_2--scream-output-preset-5",
            "ERS_Ln22.ne4pg2_ne4pg2.F2010-SCREAMv1.scream-small_kernels--scream-output-preset-5",
            "ERS_Ln22.ne4pg2_ne4pg2.F2010-SCREAMv1.scream-small_kernels_p3--scream-output-preset-5",
            "ERS_Ln22.ne4pg2_ne4pg2.F2010-SCREAMv1.scream-small_kernels_shoc--scream-output-preset-5",
            "SMS_D_Ln5.ne4pg2_oQU480.F2010-SCREAMv1-MPASSI.scream-mam4xx-all_mam4xx_procs",
            )
    },

    "e3sm_scream_v1_dp-eamxx" : {
        "time"  : "01:00:00",
        # each test runs with 225 dynamics and 100 physics columns, roughly size of ne2
        "tests" : (
            "ERS_P16_Ln22.ne30pg2_ne30pg2.FIOP-SCREAMv1-DP.scream-dpxx-dycomsrf01",
            "ERS_P16_Ln22.ne30pg2_ne30pg2.FIOP-SCREAMv1-DP.scream-dpxx-arm97",
            "ERS_P16_Ln22.ne30pg2_ne30pg2.FIOP-SCREAMv1-DP.scream-dpxx-comble",
            "ERS_P16_Ln22.ne30pg2_ne30pg2.FRCE-SCREAMv1-DP",
            )
    },

    # Tests run on exclusively on mappy for scream AT testing. These tests
    # should be fast, so we limit it to low res and add some thread tests
    # specifically for mappy.
    "e3sm_scream_v1_at" : {
        "inherit" : ("e3sm_scream_v1_lowres", "e3sm_scream_v1_dp-eamxx"),
        "tests"   : ("PET_Ln9_P32x2.ne4pg2_ne4pg2.F2010-SCREAMv1.scream-output-preset-1")
    },

    "e3sm_scream_v1_medres" : {
        "time"  : "02:00:00",
        "tests" : (
            #  "SMS_D_Ln2.ne30_ne30.F2000-SCREAMv1-AQP1", # Uncomment once IC file for ne30 is ready
            "ERS_Ln22.ne30_ne30.F2010-SCREAMv1.scream-internal_diagnostics_level--scream-output-preset-3",
            "PEM_Ln90.ne30pg2_ne30pg2.F2010-SCREAMv1.scream-spa_remap--scream-output-preset-4",
            "ERS_Ln90.ne30pg2_ne30pg2.F2010-SCREAMv1.scream-small_kernels--scream-output-preset-5",
            "ERP_Ln22.conusx4v1pg2_r05_oECv3.F2010-SCREAMv1-noAero.scream-bfbhash--scream-output-preset-6",
            "ERS_Ln22.ne30pg2_ne30pg2.F2010-SCREAMv1.scream-L128--scream-output-preset-4",
            "REP_Ld5.ne30pg2_ne30pg2.F2010-SCREAMv1.scream-L128--scream-output-preset-6"
            )
    },

    # Used to track performance
    "e3sm_scream_v1_hires" : {
        "time"  : "01:00:00",
        "tests" : (
            "SMS_Ln300.ne30pg2_ne30pg2.F2010-SCREAMv1.scream-perf_test--scream-output-preset-1"
            )
    },

    "e3sm_scream_v1_mpassi" : {
        "time"  : "01:00:00",
        "tests" : (
         #  "ERP_D_Ln9.ne4_oQU240.F2010-SCREAMv1-MPASSI.atmlndactive-rtm_off",
         #  "SMS_D_Ln9.ne4_oQU240.F2010-SCREAMv1-MPASSI-noAero",
         # Disable the two 111422-commented tests b/c they fail on pm-gpu and
         # we're not using MPASSI right now.
         #111422 "ERP_Ln22.ne4pg2_oQU480.F2010-SCREAMv1-MPASSI.atmlndactive-rtm_off",
         "ERS_D_Ln22.ne4pg2_oQU480.F2010-SCREAMv1-MPASSI.atmlndactive-rtm_off--scream-output-preset-1",
         #  "ERS_Ln22.ne30_oECv3.F2010-SCREAMv1-MPASSI.atmlndactive-rtm_off",
         #111422 "PEM_Ln90.ne30pg2_EC30to60E2r2.F2010-SCREAMv1-MPASSI",
         #  "ERS_Ln22.ne30pg2_EC30to60E2r2.F2010-SCREAMv1-MPASSI.atmlndactive-rtm_off",
            )
    },

    "e3sm_scream_v1_long" : {
        "time"  : "01:00:00",
        "tests" : (
            "ERP_D_Lh182.ne4pg2_ne4pg2.F2010-SCREAMv1",
            "ERS_Ln362.ne30pg2_ne30pg2.F2010-SCREAMv1"
            )
    },

    "e3sm_scream_v1_long_crusher" : {
        # _D builds take a long longer on crusher than ascent or pm-gpu, so
        # don't run the long _D test.
        "time"  : "01:00:00",
        "tests" : (
            "ERS_Ln362.ne30pg2_ne30pg2.F2010-SCREAMv1"
            )
    },

    "e3sm_scream_mam4xx_v1_lowres" : {
        "time"  : "01:00:00",
        "tests" : (
            "SMS_D_Ln5.ne4pg2_oQU480.F2010-SCREAMv1-MPASSI.scream-mam4xx-optics",
            "SMS_D_Ln5.ne4pg2_oQU480.F2010-SCREAMv1-MPASSI.scream-mam4xx-aci",
            "SMS_D_Ln5.ne4pg2_oQU480.F2010-SCREAMv1-MPASSI.scream-mam4xx-wetscav",
            "SMS_D_Ln5.ne4pg2_oQU480.F2010-SCREAMv1-MPASSI.scream-mam4xx-drydep",
        )
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
                 "ERP_Ln9.ne4pg2_ne4pg2.F2010-SCREAMv1",
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

    "e3sm_wav_developer" : {
        "time"    : "1:00:00",
        "tests"   : (
            "SMS_D_Ln3.TL319_EC30to60E2r2_wQU225EC30to60E2r2.GMPAS-JRA1p5-WW3.ww3-jra_1958",
            "ERS.ne30pg2_IcoswISC30E3r5_wQU225Icos30E3r5.WCYCL1850-WW3",
            "PEM_P480.ne30pg2_IcoswISC30E3r5_wQU225Icos30E3r5.WCYCL1850-WW3",
            "PET.ne30pg2_IcoswISC30E3r5_wQU225Icos30E3r5.WCYCL1850-WW3",
            "SMS_D_Ln3.ne30pg2_IcoswISC30E3r5_wQU225Icos30E3r5.WCYCL1850-WW3",
            )
    },

    # super-BFB OCN
    "e3sm_superbfb_ocn_opt" : { # opt + pureMPI
        "share"   : True,
        "time"    : "00:30:00",
        "tests"   : (
            "ERS_Ld3.T62_IcoswISC30E3r5.CMPASO-NYF.pemod-omp1",
            "PEM_Lh3.T62_IcoswISC30E3r5.CMPASO-NYF.pemod-omp1",
        )
    },

    "e3sm_superbfb_ocn_dbg" : { # dbg + pureMPI
        "share"   : True,
        "time"    : "01:00:00",
        "tests"   : (
            "ERS_Ld3_D.T62_IcoswISC30E3r5.CMPASO-NYF.pemod-omp1",
            "PEM_Lh3_D.T62_IcoswISC30E3r5.CMPASO-NYF.pemod-omp1",
        )
    },

    "e3sm_superbfb_ocn_opt_thrd" : { # opt + threads
        "share"   : True,
        "time"    : "00:30:00",
        "tests"   : (
            "ERS_Ld3.T62_IcoswISC30E3r5.CMPASO-NYF.pemod-omp2",
            "PET_Lh3.T62_IcoswISC30E3r5.CMPASO-NYF.pemod-ompfull",
        )
    },

    "e3sm_superbfb_ocn_dbg_thrd" : { # dbg + threads
        "share"   : True,
        "time"    : "01:00:00",
        "tests"   : (
            "ERS_Ld3_D.T62_IcoswISC30E3r5.CMPASO-NYF.pemod-omp2",
            "PET_Lh3_D.T62_IcoswISC30E3r5.CMPASO-NYF.pemod-ompfull",
        )
    },

    "e3sm_superbfb_ocn" : {
        "inherit" : ("e3sm_superbfb_ocn_dbg", "e3sm_superbfb_ocn_opt",
                     "e3sm_superbfb_ocn_dbg_thrd", "e3sm_superbfb_ocn_opt_thrd"),
    },

    # super-BFB ICE
    "e3sm_superbfb_ice_opt" : { # opt + pureMPI
        "share"   : True,
        "time"    : "00:30:00",
        "tests"   : (
            "ERS_Lh3.T62_IcoswISC30E3r5.DTESTM.pemod-omp1",
            "PEM_Lh3.T62_IcoswISC30E3r5.DTESTM.pemod-omp1",
        )
    },

    "e3sm_superbfb_ice_dbg" : { # dbg + pureMPI
        "share"   : True,
        "time"    : "01:00:00",
        "tests"   : (
            "ERS_Lh3_D.T62_IcoswISC30E3r5.DTESTM.pemod-omp1",
            "PEM_Lh3_D.T62_IcoswISC30E3r5.DTESTM.pemod-omp1",
        )
    },

    "e3sm_superbfb_ice_opt_thrd" : { # opt + threads
        "share"   : True,
        "time"    : "00:30:00",
        "tests"   : (
            "ERS_Lh3.T62_IcoswISC30E3r5.DTESTM.pemod-ompfull",
            "PET_Lh3.T62_IcoswISC30E3r5.DTESTM.pemod-omp2",
        )
    },

    "e3sm_superbfb_ice_dbg_thrd" : { # dbg + threads
        "share"   : True,
        "time"    : "01:00:00",
        "tests"   : (
            "ERS_Lh3_D.T62_IcoswISC30E3r5.DTESTM.pemod-ompfull",
            "PET_Lh3_D.T62_IcoswISC30E3r5.DTESTM.pemod-omp2",
        )
    },

    "e3sm_superbfb_ice" : {
        "inherit" : ("e3sm_superbfb_ice_dbg", "e3sm_superbfb_ice_opt",
                     "e3sm_superbfb_ice_dbg_thrd", "e3sm_superbfb_ice_opt_thrd"),
    },

    # super-BFB LND
    "e3sm_superbfb_lnd_opt" : { # opt + pureMPI
        "share"   : True,
        "time"    : "00:30:00",
        "tests"   : (
            "ERS_Lh3.r05_r05.IELMTEST.pemod-omp1",
            "PEM_Lh3.r05_r05.IELMTEST.pemod-omp1",
        )
    },

    "e3sm_superbfb_lnd_dbg" : { # dbg + pureMPI
        "share"   : True,
        "time"    : "01:00:00",
        "tests"   : (
            "ERS_Lh3_D.r05_r05.IELMTEST.pemod-omp1",
            "PEM_Lh3_D.r05_r05.IELMTEST.pemod-omp1",
        )
    },

    "e3sm_superbfb_lnd_opt_thrd" : { # opt + threads
        "share"   : True,
        "time"    : "00:30:00",
        "tests"   : (
            "PET_Lh3.r05_r05.IELMTEST.pemod-omp2",
            "ERS_Lh3.r05_r05.IELMTEST.pemod-ompfull",
        )
    },

    "e3sm_superbfb_lnd_dbg_thrd" : { # dbg + threads
        "share"   : True,
        "time"    : "01:00:00",
        "tests"   : (
            "PET_Lh3_D.r05_r05.IELMTEST.pemod-omp2",
            "ERS_Lh3_D.r05_r05.IELMTEST.pemod-ompfull",
        )
    },

    "e3sm_superbfb_lnd" : {
        "inherit" : ("e3sm_superbfb_lnd_dbg", "e3sm_superbfb_lnd_opt",
                     "e3sm_superbfb_lnd_dbg_thrd", "e3sm_superbfb_lnd_opt_thrd"),
    },

    # super-BFB ROF
    "e3sm_superbfb_rof_opt" : { # opt + pureMPI
        "share"   : True,
        "time"    : "00:30:00",
        "tests"   : (
            "ERS_Ld3.r05_r05.RMOSGPCC.pemod-omp1",
            "PEM_Ld3.r05_r05.RMOSGPCC.pemod-omp1",
        )
    },

    "e3sm_superbfb_rof_dbg" : { # dbg + pureMPI
        "share"   : True,
        "time"    : "01:00:00",
        "tests"   : (
            "ERS_Ld3_D.r05_r05.RMOSGPCC.pemod-omp1",
            "PEM_Ld3_D.r05_r05.RMOSGPCC.pemod-omp1",
        )
    },

    "e3sm_superbfb_rof_opt_thrd" : { # opt + threads
        "share"   : True,
        "time"    : "00:30:00",
        "tests"   : (
            "PET_Ld3.r05_r05.RMOSGPCC.pemod-omp2",
            "ERS_Ld3.r05_r05.RMOSGPCC.pemod-ompfull",
        )
    },

    "e3sm_superbfb_rof_dbg_thrd" : { # dbg + threads
        "share"   : True,
        "time"    : "01:00:00",
        "tests"   : (
            "PET_Ld3_D.r05_r05.RMOSGPCC.pemod-omp2",
            "ERS_Ld3_D.r05_r05.RMOSGPCC.pemod-ompfull",
        )
    },

    "e3sm_superbfb_rof" : {
        "inherit" : ("e3sm_superbfb_rof_dbg", "e3sm_superbfb_rof_opt",
                     "e3sm_superbfb_rof_dbg_thrd", "e3sm_superbfb_rof_opt_thrd"),
    },

    # super-BFB ATM
    "e3sm_superbfb_atm_opt" : { # opt + pureMPI
        "share"   : True,
        "time"    : "00:30:00",
        "tests"   : (
            "ERS_Lh3.ne30pg2_ne30pg2.FAQP.pemod-omp1",
            "PEM_Lh3.ne30pg2_ne30pg2.FAQP.pemod-omp1",
        )
    },

    "e3sm_superbfb_atm_dbg" : { # dbg + pureMPI
        "share"   : True,
        "time"    : "01:00:00",
        "tests"   : (
            "ERS_Lh3_D.ne30pg2_ne30pg2.FAQP.pemod-omp1",
            "PEM_Lh3_D.ne30pg2_ne30pg2.FAQP.pemod-omp1",
        )
    },

    "e3sm_superbfb_atm_opt_thrd" : { # opt + threads
        "share"   : True,
        "time"    : "00:30:00",
        "tests"   : (
            "PET_Lh3.ne30pg2_ne30pg2.FAQP.pemod-omp2",
            "ERS_Lh3.ne30pg2_ne30pg2.FAQP.pemod-ompfull",
        )
    },

    "e3sm_superbfb_atm_dbg_thrd" : { # dbg + threads
        "share"   : True,
        "time"    : "01:00:00",
        "tests"   : (
            "PET_Lh3_D.ne30pg2_ne30pg2.FAQP.pemod-omp2",
            "ERS_Lh3_D.ne30pg2_ne30pg2.FAQP.pemod-ompfull",
        )
    },

    "e3sm_superbfb_atm" : {
        "inherit" : ("e3sm_superbfb_atm_dbg", "e3sm_superbfb_atm_opt",
                     "e3sm_superbfb_atm_dbg_thrd", "e3sm_superbfb_atm_opt_thrd"),
    },

    # super-BFB all-active
    "e3sm_superbfb_wcycl_opt" : { # opt + pureMPI
        "share"   : True,
        "time"    : "01:00:00",
        "tests"   : (
            "ERS_Ld3.ne30pg2_r05_IcoswISC30E3r5.WCYCL1850.pemod-omp1",
            "PEM_Ld3.ne30pg2_r05_IcoswISC30E3r5.WCYCL1850.pemod-omp1",
        )
    },

    "e3sm_superbfb_wcycl_dbg" : { # dbg + pureMPI
        "share"   : True,
        "time"    : "02:00:00",
        "tests"   : (
            "ERS_Ld3_D.ne30pg2_r05_IcoswISC30E3r5.WCYCL1850.pemod-omp1",
            "PEM_Ld3_D.ne30pg2_r05_IcoswISC30E3r5.WCYCL1850.pemod-omp1",
        )
    },

    "e3sm_superbfb_wcycl_opt_thrd" : { # opt + threads
        "share"   : True,
        "time"    : "01:00:00",
        "tests"   : (
            "PET_Ld3.ne30pg2_r05_IcoswISC30E3r5.WCYCL1850.pemod-omp2",
            "ERS_Ld3.ne30pg2_r05_IcoswISC30E3r5.WCYCL1850.pemod-ompfull",
        )
    },

    "e3sm_superbfb_wcycl_dbg_thrd" : { # dbg + threads
        "share"   : True,
        "time"    : "02:00:00",
        "tests"   : (
            "PET_Ld3_D.ne30pg2_r05_IcoswISC30E3r5.WCYCL1850.pemod-omp2",
            "ERS_Ld3_D.ne30pg2_r05_IcoswISC30E3r5.WCYCL1850.pemod-omp2",
        )
    },

    "e3sm_superbfb_wcycl" : {
        "inherit" : ("e3sm_superbfb_wcycl_dbg", "e3sm_superbfb_wcycl_opt",
                     "e3sm_superbfb_wcycl_dbg_thrd", "e3sm_superbfb_wcycl_opt_thrd"),
    },

    "e3sm_superbfb" : {
        "inherit" : ("e3sm_superbfb_ocn", "e3sm_superbfb_ice",
                     "e3sm_superbfb_lnd", "e3sm_superbfb_rof",
                     "e3sm_superbfb_atm", "e3sm_superbfb_wcycl"),
    },
}
