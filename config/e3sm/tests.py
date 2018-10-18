# Here are the tests belonging to e3sm suites. Format is
# <test>.<grid>.<compset>.
# suite_name -> (inherits_from, timelimit, [test [, mods[, machines]]])
#   To elaborate, if no mods are needed, a string representing the testname is all that is needed.
#   If testmods are needed, a 2-ple must be provided  (test, mods)
#   If you want to restrict the test mods to certain machines, than a 3-ple is needed (test, mods, [machines])
_TESTS = {

    "e3sm_land_developer" : (None, "0:45:00",
                             ("ERS.f19_f19.ICLM45",
                              "ERS.f19_f19.I1850CLM45CN",
                              "ERS.f09_g16.I1850CLM45CN",
                              "ERS.f19_f19.I20TRCLM45CN",
                              "SMS_Ld1.hcru_hcru.I1850CRUCLM45CN",
                             ("ERS.f19_g16.I1850CNECACNTBC" ,"clm-eca"),
                             ("ERS.f19_g16.I1850CNECACTCBC" ,"clm-eca"),
                             ("SMS_Ly2_P1x1.1x1_smallvilleIA.ICLM45CNCROP", "clm-force_netcdf_pio"),
                             ("ERS_Ld3.f45_f45.ICLM45ED","clm-fates"),
                             ("ERS.f19_g16.I1850CLM45","clm-betr"),
                             ("ERS.f19_g16.I1850CLM45","clm-vst"),
                             ("ERS.f09_g16.I1850CLM45CN","clm-bgcinterface"),
                              "ERS.ne11_oQU240.I20TRCLM45",
                             ("ERS.f19_g16.I1850CNRDCTCBC","clm-rd"),
                             ("ERS.f19_g16.I1850GSWCNPECACNTBC","clm-eca_f19_g16_I1850GSWCNPECACNTBC"),
                             ("ERS.f19_g16.I20TRGSWCNPECACNTBC","clm-eca_f19_g16_I20TRGSWCNPECACNTBC"),
                             ("ERS.f19_g16.I1850GSWCNPRDCTCBC","clm-ctc_f19_g16_I1850GSWCNPRDCTCBC"),
                             ("ERS.f19_g16.I20TRGSWCNPRDCTCBC","clm-ctc_f19_g16_I20TRGSWCNPRDCTCBC"),
                              "ERS.f09_g16.ICLM45BC")
                             ),

    "e3sm_atm_developer" : (None, None,
                            ("ERP_Ln9.ne4_ne4.FC5AV1C-L",
                             ("SMS_Ln9.ne4_ne4.FC5AV1C-L", "cam-outfrq9s"),
                             ("SMS.ne4_ne4.FC5AV1C-L", "cam-cosplite"),
                             "SMS_R_Ld5.T42_T42.FSCM5A97",
                             "SMS_D_Ln5.ne4_ne4.FC5AV1C-L")
                            ),

    "e3sm_atm_integration" : (None, None,
                              ("ERP_Ln9.ne4_ne4.FC5AV1C-L-AQUAP",
                              ("SMS_Ld1.ne4_ne4.FC5AV1C-L-AQUAP","cam-clubb_only"),
                               ("PET_Ln5.ne4_ne4.FC5AV1C-L","allactive-mach-pet"),
                               "PEM_Ln5.ne4_ne4.FC5AV1C-L",
                               ("SMS_D_Ln5.ne4_ne4.FC5AV1C-L", "cam-cosplite_nhtfrq5"),
                               ("ERS_Ld5.ne4_ne4.FC5AV1C-L", "cam-rrtmgp"),
                               "REP_Ln5.ne4_ne4.FC5AV1C-L")
                              ),
    #atmopheric tests for extra coverage
    "e3sm_atm_extra_coverage" : (None, None,
                         ("SMS_Lm1.ne4_ne4.FC5AV1C-L",
                          "ERS_Ld31.ne4_ne4.FC5AV1C-L",
                          "ERP_Lm3.ne4_ne4.FC5AV1C-L",
                          "SMS_D_Ln5.ne30_ne30.FC5AV1C-L",
                          ("ERP_Ln5.ne30_ne30.FC5AV1C-L"),
                          "SMS_Ly1.ne4_ne4.FC5AV1C-L")
                         ),
    #atmopheric tests for hi-res
    "e3sm_atm_hi_res" : (None, "01:30:00",
                         (
                          "SMS.ne120_ne120.FC5AV1C-H01A",
                         )),
    #atmopheric tests to mimic low res production runs
    "e3sm_atm_prod" : (None, None,
                       (("SMS_Ln5.ne30_ne30.FC5AV1C-L", "cam-cosplite"),
                        )
                       ),

    #atmopheric nbfb tests
    "e3sm_atm_nbfb" : (None, None,
                                 ("PGN_P1x1.ne4_ne4.FC5AV1C-L",
                                  "TSC.ne4_ne4.FC5AV1C-L")
                                 ),

    "e3sm_developer" : (("e3sm_land_developer","e3sm_atm_developer"), "0:45:00",
                        ("ERS.f19_g16_rx1.A",
                         "ERS.ne30_g16_rx1.A",
                         "SEQ.f19_g16.X",
                         "ERIO.ne30_g16_rx1.A",
                         "HOMME_P24.f19_g16_rx1.A",
                         "NCK.f19_g16_rx1.A",
                         "SMS.ne30_f19_g16_rx1.A",
                         "ERS_Ld5.T62_oQU120.CMPASO-NYF",
                         "ERS.f09_g16_g.MALISIA",
                         "SMS.T62_oQU120_ais20.MPAS_LISIO_TEST",
                         "SMS.f09_g16_a.IGCLM45_MLI",
                        ("SMS_P12x2.ne4_oQU240.A_WCYCL1850","allactive-mach_mods")
                        )),

    "e3sm_integration" : (("e3sm_developer", "e3sm_atm_integration"),"03:00:00",
                          ("ERS.ne11_oQU240.A_WCYCL1850",
		           ("SMS_D_Ld1.ne30_oECv3_ICG.A_WCYCL1850S_CMIP6","allactive-v1cmip6"),
                           "ERS_Ln9.ne4_ne4.FC5AV1C-L",
                          #"ERT_Ld31.ne16_g37.B1850C5",#add this line back in with the new correct compset
                           "NCK.ne11_oQU240.A_WCYCL1850",
                           ("PET.f19_g16.X","allactive-mach-pet"),
                           ("PET.f45_g37_rx1.A","allactive-mach-pet"),
                           ("PET_Ln9_PS.ne30_oECv3_ICG.A_WCYCL1850S","allactive-mach-pet"),
                           "PEM_Ln9.ne30_oECv3_ICG.A_WCYCL1850S",
                           "ERP_Ld3.ne30_oECv3_ICG.A_WCYCL1850S",
                           "SMS.f09_g16_a.MALI",
                           "SMS_D_Ln5.conusx4v1_conusx4v1.FC5AV1C-L",
                           ("SMS.ne30_oECv3.BGCEXP_BCRC_CNPECACNT_1850","clm-bgcexp"),
                           ("SMS.ne30_oECv3.BGCEXP_BCRC_CNPRDCTC_1850","clm-bgcexp"))
                          ),
    #e3sm tests for extra coverage
    "e3sm_extra_coverage" : (("e3sm_atm_extra_coverage",), None,
                             ("SMS_D_Ln5.enax4v1_enax4v1.FC5AV1C-L",
                              "SMS_D_Ln5.twpx4v1_twpx4v1.FC5AV1C-L")),

    #e3sm tests for hi-res
    "e3sm_hi_res" : (("e3sm_atm_hi_res",),None,
                     (
                      ("SMS.ne120_oRRS18v3_ICG.A_WCYCL2000_H01AS", "cam-cosplite"),
                       "SMS.T62_oRRS30to10v3wLI.GMPAS-IAF",
                     )),

    #e3sm tests for RRM grids
    "e3sm_rrm" : (None, None,
                  ("SMS_D_Ln5.conusx4v1_conusx4v1.FC5AV1C-L",
                   "SMS_D_Ln5.enax4v1_enax4v1.FC5AV1C-L",
                   "SMS_D_Ln5.twpx4v1_twpx4v1.FC5AV1C-L")
                 ),

    #e3sm tests to mimic production runs
    "e3sm_prod" : (("e3sm_atm_prod",),None,
                     (
		      ("SMS_Ld2.ne30_oECv3_ICG.A_WCYCL1850S_CMIP6","allactive-v1cmip6"),
		      )),

    "fates" : (None, None,
                         ("ERS_Ld9.1x1_brazil.ICLM45ED",
                          "ERS_D_Ld9.1x1_brazil.ICLM45ED",
                          "SMS_D_Lm6.1x1_brazil.ICLM45ED")
               ),

}
