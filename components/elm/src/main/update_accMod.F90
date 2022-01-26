module update_accMod

  public :: update_acc_variables
contains


subroutine update_acc_variables()
  use pftvarcon
  use elm_varctl
  use elm_varcon
  use elm_varpar
  use soilorder_varcon
  use LakeCon
  use SoilWaterMovementMod
  use SharedParamsMod
  use MaintenanceRespMod
  use NitrifDenitrifMod
  use CNStateType
  use AllocationMod
  !---------- MaintenanceRespMod ------------!
  !$acc update device(br_mr_Inst)
  !----------- NitrifDenitrifMod -----------!
  !$acc update device (no_frozen_nitrif_denitrif )
  !---------- CNStateType ---------------!
  !$acc update device(fert_type  &
  !$acc               ,fert_continue &
  !$acc               ,fert_dose     &
  !$acc               ,fert_start    &
  !$acc               ,fert_end      &
  !---------- soilorder_varcon -------------!
  !$acc               ,smax(:)       &
  !$acc               ,ks_sorption(:)&
  !$acc               ,r_weather(:)  &
  !$acc               ,r_adsorp(:)   &
  !$acc               ,r_desorp(:)   &
  !$acc               ,r_occlude(:)  &
  !$acc               ,k_s1_biochem(:) &
  !$acc               ,k_s2_biochem(:) &
  !$acc               ,k_s3_biochem(:) &
  !$acc               ,k_s4_biochem(:) &
  !$acc               ,r_mort_soilorder(:) )
  !---------- SharedParamsMod -------------- !
  !$acc update device(anoxia_wtsat, nlev_soildecomp_standard )
  !--------- elm_varpar --------------------!
  !$acc update device( &
  !$acc       nlevsoi         &
  !$acc       ,nlevsoifl       &
  !$acc       ,nlevurb         &
  !$acc       ,nlevlak         &
  !$acc       ,nlevdecomp      &
  !$acc       ,nlevdecomp_full &
  !$acc       ,nlevtrc_soil    &
  !$acc       ,nlevtrc_full    &
  !$acc       ,nlevgrnd        &
  !$acc       ,natpft_lb       &
  !$acc       ,natpft_ub       &
  !$acc       ,natpft_size     &
  !$acc       ,cft_lb          &
  !$acc       ,cft_ub          &
  !$acc       ,cft_size        &
  !$acc       ,i_met_lit       &
  !$acc       ,i_cel_lit       &
  !$acc       ,i_lig_lit       &
  !$acc       ,i_cwd           &
  !$acc       ,maxpatch_glcmec   &
  !$acc       ,max_patch_per_col &
  !$acc       ,mach_eps         &
  !$acc       ,maxpatch_pft )
  !--------- elm_varcon --------------------!
  !$acc update device(&
  !$acc       zlak(:)        &
  !$acc       ,dzlak(:)       &
  !$acc       ,zsoi(:)        &
  !$acc       ,dzsoi(:)       &
  !$acc       ,zisoi(:)       &
  !$acc       ,dzsoi_decomp(:) &
  !$acc       ,nlvic(:)       &
  !$acc       ,dzvic(:)       &
  !$acc       ,zsoifl(:)      &
  !$acc       ,zisoifl(:)     &
  !$acc       ,dzsoifl(:)     &
  !$acc       ,denh2o         &
  !-------- elm_varctl ---------------------!
  !$acc       ,glcmec_downscale_rain_snow_convert &
  !$acc       ,glcmec_downscale_longwave &
  !$acc       ,subgridflag &
  !$acc       ,use_nofire          &
  !$acc       ,use_lch4            &
  !$acc       ,use_vertsoilc       &
  !$acc       ,use_extralakelayers &
  !$acc       ,use_vichydro        &
  !$acc       ,use_century_decomp  &
  !$acc       ,use_cn              )
  ! -------- pftvarcon ------------------!
  !$acc update device( &
  !$acc        dleaf(:)       &
  !$acc       ,c3psn(:)       &
  !$acc       ,xl(:)          &
  !$acc       ,rhol(:,:)      &
  !$acc       ,rhos(:,:)      &
  !$acc       ,taul(:,:)      &
  !$acc       ,taus(:,:)      &
  !$acc       ,z0mr(:)        &
  !$acc       ,displar(:)     &
  !$acc       ,roota_par(:)   &
  !$acc       ,rootb_par(:)   &
  !$acc       ,crop(:)        &
  !$acc       ,irrigated(:)   &
  !$acc       ,smpso(:)       &
  !$acc       ,smpsc(:)       &
  !$acc       ,fnitr(:)       &
  !$acc       ,slatop(:)      &
  !$acc       ,dsladlai(:)    &
  !$acc       ,leafcn(:)      &
  !$acc       ,flnr(:)        &
  !$acc       ,woody(:)       &
  !$acc       ,lflitcn(:)     &
  !$acc       ,frootcn(:)     &
  !$acc       ,livewdcn(:)    &
  !$acc       ,deadwdcn(:)    &
  !$acc       ,grperc(:)      &
  !$acc       ,grpnow(:)      &
  !$acc       ,rootprof_beta(:) &
  !$acc       ,leafcp(:)   &
  !$acc       ,lflitcp(:)  &
  !$acc       ,frootcp(:)  &
  !$acc       ,livewdcp(:) &
  !$acc       ,deadwdcp(:) &
  !$acc    ,mergetoelmpft (:) &
  !$acc    ,is_pft_known_to_model (:) &
  !$acc        ,graincn(:)       &
  !$acc        ,graincp(:)       &
  !$acc        ,mxtmp(:)         &
  !$acc        ,baset(:)         &
  !$acc        ,declfact(:)      &
  !$acc        ,bfact(:)         &
  !$acc        ,aleaff(:)        &
  !$acc        ,arootf(:)        &
  !$acc        ,astemf(:)        &
  !$acc        ,arooti(:)        &
  !$acc        ,fleafi(:)        &
  !$acc        ,allconsl(:)      &
  !$acc        ,allconss(:)      &
  !$acc        ,ztopmx(:)        &
  !$acc        ,laimx(:)         &
  !$acc        ,gddmin(:)        &
  !$acc        ,hybgdd(:)        &
  !$acc        ,lfemerg(:)       &
  !$acc        ,grnfill(:)       &
  !$acc        ,mxmat(:)         &
  !$acc        ,mnNHplantdate(:) &
  !$acc        ,mxNHplantdate(:) &
  !$acc        ,mnSHplantdate(:) &
  !$acc        ,mxSHplantdate(:) &
  !$acc        ,planttemp(:)     &
  !$acc        ,minplanttemp(:)  &
  !$acc        ,froot_leaf(:)    &
  !$acc        ,stem_leaf(:)     &
  !$acc        ,croot_stem(:)    &
  !$acc        ,flivewd(:)       &
  !$acc        ,fcur(:)          &
  !$acc        ,lf_flab(:)       &
  !$acc        ,lf_fcel(:)       &
  !$acc        ,lf_flig(:)       &
  !$acc        ,fr_flab(:)       &
  !$acc        ,fr_fcel(:)       &
  !$acc        ,fr_flig(:)       &
  !$acc        ,leaf_long(:)     &
  !$acc        ,froot_long(:)    &
  !$acc        ,evergreen(:)     &
  !$acc        ,stress_decid(:)  &
  !$acc        ,season_decid(:)  &
  !$acc        ,pconv(:)         &
  !$acc        ,pprod10(:)       &
  !$acc        ,pprod100(:)      &
  !$acc        ,pprodharv10(:)   &
  !$acc       ,cc_leaf(:)  &
  !$acc       ,cc_lstem(:) &
  !$acc       ,cc_dstem(:) &
  !$acc       ,cc_other(:) &
  !$acc       ,fm_leaf(:)  &
  !$acc       ,fm_lstem(:) &
  !$acc       ,fm_dstem(:) &
  !$acc       ,fm_other(:) &
  !$acc       ,fm_root(:)  &
  !$acc       ,fm_lroot(:) &
  !$acc       ,fm_droot(:) &
  !$acc       ,fsr_pft(:)  &
  !$acc       ,fd_pft(:)   &
  !$acc       ,fertnitro(:) &
  !$acc       ,fleafcn(:)   &
  !$acc       ,ffrootcn(:)  &
  !$acc       ,fstemcn(:)   &
  !$acc       ,presharv(:)  &
  !$acc       ,convfact(:)  &
  !$acc       ,fyield(:)    &
  !$acc       ,root_dmx(:)  &
  !$acc     ,VMAX_PLANT_NH4(:)   &
  !$acc     ,VMAX_PLANT_NO3(:)   &
  !$acc     ,VMAX_PLANT_P(:)     &
  !$acc     ,VMAX_MINSURF_P_vr(:,:) &
  !$acc     ,KM_PLANT_NH4(:)      &
  !$acc     ,KM_PLANT_NO3(:)      &
  !$acc     ,KM_PLANT_P(:)        &
  !$acc     ,KM_MINSURF_P_vr(:,:) )
  !$acc update device( &
  !$acc     KM_DECOMP_NH4        &
  !$acc     ,KM_DECOMP_NO3        &
  !$acc     ,KM_DECOMP_P          &
  !$acc     ,KM_NIT               &
  !$acc     ,KM_DEN               &
  !$acc     ,decompmicc_patch_vr(:,:) &
  !$acc     ,alpha_nfix(:)            &
  !$acc     ,alpha_ptase(:)           &
  !$acc     ,ccost_nfix(:)            &
  !$acc     ,pcost_nfix(:)            &
  !$acc     ,ccost_ptase(:)           &
  !$acc     ,ncost_ptase(:)           &
  !$acc     ,VMAX_NFIX(:)       &
  !$acc     ,KM_NFIX(:)         &
  !$acc     ,VMAX_PTASE(:)      &
  !$acc     ,KM_PTASE           &
  !$acc     ,lamda_ptase        &
  !$acc     ,i_vc(:)            &
  !$acc     ,s_vc(:)            &
  !$acc     ,leafcn_obs(:)          &
  !$acc     ,frootcn_obs(:)         &
  !$acc     ,livewdcn_obs(:)        &
  !$acc     ,deadwdcn_obs(:)        &
  !$acc     ,leafcp_obs(:)          &
  !$acc     ,frootcp_obs(:)         &
  !$acc     ,livewdcp_obs(:)        &
  !$acc     ,deadwdcp_obs(:)        &
  !$acc     ,leafcn_obs_flex(:,:)   &
  !$acc     ,frootcn_obs_flex(:,:)  &
  !$acc     ,livewdcn_obs_flex(:,:) &
  !$acc     ,deadwdcn_obs_flex(:,:) &
  !$acc     ,leafcp_obs_flex(:,:)   &
  !$acc     ,frootcp_obs_flex(:,:)  &
  !$acc     ,livewdcp_obs_flex(:,:) &
  !$acc     ,deadwdcp_obs_flex(:,:) &
  !$acc     ,fnr(:)        &
  !$acc     ,act25(:)      &
  !$acc     ,kcha(:)       &
  !$acc     ,koha(:)       &
  !$acc     ,cpha(:)       &
  !$acc     ,vcmaxha(:)    &
  !$acc     ,jmaxha(:)     &
  !$acc     ,tpuha(:)      &
  !$acc     ,lmrha(:)      &
  !$acc     ,vcmaxhd(:)    &
  !$acc     ,jmaxhd(:)     &
  !$acc     ,tpuhd(:)      &
  !$acc     ,lmrhd(:)      &
  !$acc     ,lmrse(:)      &
  !$acc     ,qe(:)         &
  !$acc     ,theta_cj(:)   &
  !$acc     ,bbbopt(:)     &
  !$acc     ,mbbopt(:)     &
  !$acc     ,nstor(:)      &
  !$acc     ,br_xr(:)      &
  !$acc     ,tc_stress     &
  !$acc     ,vcmax_np1(:)  &
  !$acc     ,vcmax_np2(:)  &
  !$acc     ,vcmax_np3(:)  &
  !$acc     ,vcmax_np4(:)  &
  !$acc     ,jmax_np1      &
  !$acc     ,jmax_np2      &
  !$acc     ,jmax_np3      &
  !$acc     ,laimax        &
  !$acc    ,rsub_top_globalmax &
  !------------- LakeCon ------------------!
  !$acc    ,fcrit      &
  !$acc    ,minz0lake  &
  !$acc     ,pudz &
  !$acc     ,depthcrit &
  !$acc     ,mixfact &
  !$acc     ,betavis &
  !$acc     ,lakepuddling &
  !$acc     ,lake_no_ed )
  !------------ SoilWaterMovementMod ---------------- !
  !$acc update device(soilroot_water_method )
  !----------- AllocationMod ------------------- !
  !$acc update device( nu_com_leaf_physiology &
  !$acc    , nu_com_root_kinetics   &
  !$acc    , nu_com_phosphatase     &
  !$acc    , nu_com_nfix            &
  !$acc    , bdnr                   &
  !$acc    , dayscrecover           &
  !$acc    , arepr(:)               &
  !$acc    , aroot(:)               &
  !$acc    , decompmicc(:)          &
  !$acc    , crop_supln             &
  !$acc     )

end subroutine update_acc_variables


end module update_accMod
