module dynConsBiogeochemMod
   
   #include "shr_assert.h"
   
   !---------------------------------------------------------------------------
   !
   ! !DESCRIPTION:
   ! Handle conservation of biogeochemical quantities (C & N) with dynamic land cover.
   !
   ! !USES:
   use shr_kind_mod             , only : r8 => shr_kind_r8
   use shr_log_mod              , only : errMsg => shr_log_errMsg
   use decompMod                , only : bounds_type
   use abortutils               , only : endrun
   use elm_varctl               , only : iulog, use_c13, use_c14
   use VegetationPropertiesType , only : veg_vp
   use CanopyStateType          , only : canopystate_type
   use PhotosynthesisType       , only : photosyns_type
   use CNStateType              , only : cnstate_type
   use GridcellDataType         , only : grc_cf, c13_grc_cf, c14_grc_cf
   use GridcellDataType         , only : grc_nf, grc_pf
   use LandunitType             , only : lun_pp
   use ColumnType               , only : col_pp
   use ColumnDataType           , only : column_carbon_state, column_nitrogen_state
   use ColumnDataType           , only : column_phosphorus_state
   use ColumnDataType           , only : col_cf, c13_col_cf, c14_col_cf
   use ColumnDataType           , only : col_nf, col_pf
   use VegetationType           , only : veg_pp
   use VegetationDataType       , only : vegetation_carbon_state, vegetation_carbon_flux
   use VegetationDataType       , only : vegetation_nitrogen_state
   use VegetationDataType       , only : vegetation_phosphorus_state
   use VegetationDataType       , only : veg_cf, c13_veg_cf, c14_veg_cf
   use VegetationDataType       , only : veg_nf, veg_pf
   use elm_varcon               , only : c14ratio
   use dynPatchStateUpdaterMod  , only : patch_state_updater_type
   use dynSubgridAdjustmentsMod , only : dyn_veg_cs_Adjustments, dyn_col_cs_Adjustments
   use dynSubgridAdjustmentsMod , only : dyn_veg_ns_Adjustments, dyn_col_ns_Adjustments
   use dynSubgridAdjustmentsMod , only : dyn_veg_ps_Adjustments, dyn_col_ps_Adjustments
   
   
   !
   ! !PUBLIC MEMBER FUNCTIONS:
   implicit none
   private
   
   save
   
   public :: dyn_cnbal_patch
   public :: dyn_cnbal_column
   !-----------------------------------------------------------------------
   
   contains
   
   !-----------------------------------------------------------------------
   subroutine dyn_cnbal_patch(bounds, &
      num_soilp_with_inactive, filter_soilp_with_inactive, &
      num_soilc_with_inactive, filter_soilc_with_inactive, &
      prior_weights, &
      patch_state_updater, &
      canopystate_vars, photosyns_vars, cnstate_vars, &
      veg_cs, c13_veg_cs, c14_veg_cs, &
      veg_ns, veg_ps, dt)
      !
      ! !DESCRIPTION:
      ! Modify pft-level state and flux variables to maintain carbon and nitrogen balance with
      ! dynamic pft-weights.
      !
      ! !USES:
      use shr_const_mod      , only : SHR_CONST_PDB
      use landunit_varcon    , only : istsoil, istcrop
      use elm_varpar         , only : numveg, nlevdecomp, max_patch_per_col
      use pftvarcon          , only : pconv, pprod10, pprod100
      use elm_varcon         , only : c13ratio, c14ratio
      use dynPriorWeightsMod , only : prior_weights_type
      !
      ! !ARGUMENTS:
      type(bounds_type)              , intent(in)    :: bounds
      integer                        , intent(in)    :: num_soilp_with_inactive ! number of points in filter
      integer                        , intent(in)    :: filter_soilp_with_inactive(:) ! soil patch filter that includes inactive points
      integer                        , intent(in)    :: num_soilc_with_inactive ! number of points in filter
      integer                        , intent(in)    :: filter_soilc_with_inactive(:) ! soil column filter that includes inactive points
      type(prior_weights_type)       , intent(in)    :: prior_weights ! weights prior to the subgrid weight updates
      type(patch_state_updater_type) , intent(in)    :: patch_state_updater
      type(canopystate_type)         , intent(inout) :: canopystate_vars
      type(photosyns_type)           , intent(inout) :: photosyns_vars
      type(cnstate_type)             , intent(inout) :: cnstate_vars
      type(vegetation_carbon_state)  , intent(inout) :: veg_cs
      type(vegetation_carbon_state)  , intent(inout) :: c13_veg_cs
      type(vegetation_carbon_state)  , intent(inout) :: c14_veg_cs
      type(vegetation_nitrogen_state), intent(inout) :: veg_ns
      type(vegetation_phosphorus_state),intent(inout) :: veg_ps
      real(r8)                         ,intent(in)    :: dt                            ! land model time step (sec)
      
      !
      ! !LOCAL VARIABLES:
      integer   :: p,c,l,g,j,fp,fc               ! indices
      integer   :: ier                           ! error code
      real(r8)  :: dwt                           ! change in pft weight (relative to column)
      character(len=32)             :: subname='dyn_cbal'            ! subroutine name
      
      !! ACTUAL VARIABLES that will be re-used for each species 
      real(r8)  :: dwt_leaf_seed 
      real(r8)  :: dwt_deadstem_seed  ! pft-level mass gain due to seeding of new area
      real(r8)  :: dwt_pool_seed      ! pft-level mass gain due to seeding of new area
      real(r8)  :: dwt_froot_to_litter(1:num_soilp_with_inactive)! pft-level mass loss due to weight shift
      real(r8)  :: dwt_livecroot_to_litter(1:num_soilp_with_inactive)   ! pft-level mass loss due to weight shift
      real(r8)  :: dwt_deadcroot_to_litter(1:num_soilp_with_inactive)   ! pft-level mass loss due to weight shift
      real(r8)  :: conv_flux(1:num_soilp_with_inactive)                 ! pft-level mass loss due to weight shift
      real(r8)  :: prod10_flux(1:num_soilp_with_inactive)               ! pft-level mass loss due to weight shift
      real(r8)  :: prod100_flux(1:num_soilp_with_inactive)              ! pft-level mass loss due to weight shift
      real(r8)  :: crop_product_flux(1:num_soilp_with_inactive)         ! pft-level mass loss due to weight shift
      integer   :: patch_to_soil_filter(bounds%begp:bounds%endp)
      
      
      !! C13
      real(r8), allocatable :: dwt_leafc13_seed(:)           ! pft-level mass gain due to seeding of new area
      real(r8), allocatable :: dwt_deadstemc13_seed(:)       ! pft-level mass gain due to seeding of new area
      real(r8), allocatable :: dwt_frootc13_to_litter(:)     ! pft-level mass loss due to weight shift
      real(r8), allocatable :: dwt_livecrootc13_to_litter(:) ! pft-level mass loss due to weight shift
      real(r8), allocatable :: dwt_deadcrootc13_to_litter(:) ! pft-level mass loss due to weight shift
      real(r8), allocatable :: conv_c13flux(:)               ! pft-level mass loss due to weight shift
      real(r8), allocatable :: prod10_c13flux(:)             ! pft-level mass loss due to weight shift
      real(r8), allocatable :: prod100_c13flux(:)            ! pft-level mass loss due to weight shift
      real(r8), allocatable :: crop_product_c13flux(:)       ! pft-level mass loss due to weight shift
      !! C14
      real(r8), allocatable :: dwt_leafc14_seed(:)           ! pft-level mass gain due to seeding of new area
      real(r8), allocatable :: dwt_deadstemc14_seed(:)       ! pft-level mass gain due to seeding of new area
      real(r8), allocatable :: dwt_frootc14_to_litter(:)     ! pft-level mass loss due to weight shift
      real(r8), allocatable :: dwt_livecrootc14_to_litter(:) ! pft-level mass loss due to weight shift
      real(r8), allocatable :: dwt_deadcrootc14_to_litter(:) ! pft-level mass loss due to weight shift
      real(r8), allocatable :: conv_c14flux(:)               ! pft-level mass loss due to weight shift
      real(r8), allocatable :: prod10_c14flux(:)             ! pft-level mass loss due to weight shift
      real(r8), allocatable :: prod100_c14flux(:)            ! pft-level mass loss due to weight shift
      real(r8), allocatable :: crop_product_c14flux(:)       ! pft-level mass loss due to weight shift
      real(r8) :: froot, croot
      real(r8) :: fr_flab, fr_fcel, fr_flig
      real(r8) :: startt, stopt
      real(r8) :: sum1, sum2, sum3, sum4 ,sum5
      !-----------------------------------------------------------------------
      
      if ( use_c13 ) then
         allocate(dwt_leafc13_seed           (num_soilp_with_inactive), stat=ier)
         allocate(dwt_deadstemc13_seed       (num_soilp_with_inactive), stat=ier)
         allocate(dwt_frootc13_to_litter     (num_soilp_with_inactive), stat=ier)
         allocate(dwt_livecrootc13_to_litter (num_soilp_with_inactive), stat=ier)
         allocate(dwt_deadcrootc13_to_litter (num_soilp_with_inactive), stat=ier)
         allocate(conv_c13flux               (num_soilp_with_inactive), stat=ier)
         allocate(prod10_c13flux             (num_soilp_with_inactive), stat=ier)
         allocate(prod100_c13flux            (num_soilp_with_inactive), stat=ier)
         allocate(crop_product_c13flux       (num_soilp_with_inactive), stat=ier)
      endif
      if ( use_c14 ) then
         allocate(dwt_leafc14_seed           (num_soilp_with_inactive), stat=ier)
         allocate(dwt_deadstemc14_seed       (num_soilp_with_inactive), stat=ier)
         allocate(dwt_frootc14_to_litter     (num_soilp_with_inactive), stat=ier)
         allocate(dwt_livecrootc14_to_litter (num_soilp_with_inactive), stat=ier)
         allocate(dwt_deadcrootc14_to_litter (num_soilp_with_inactive), stat=ier)
         allocate(conv_c14flux               (num_soilp_with_inactive), stat=ier)
         allocate(prod10_c14flux             (num_soilp_with_inactive), stat=ier)
         allocate(prod100_c14flux            (num_soilp_with_inactive), stat=ier)
         allocate(crop_product_c14flux       (num_soilp_with_inactive), stat=ier)
      endif
      
      call cpu_time(startt) 

      !$acc enter data create(dwt_leaf_seed,&
      !$acc dwt_deadstem_seed    ,&
      !$acc dwt_pool_seed       ,&
      !$acc dwt_froot_to_litter(:) ,&
      !$acc dwt_livecroot_to_litter(:),&
      !$acc dwt_deadcroot_to_litter(:),&
      !$acc conv_flux(:)              ,&
      !$acc prod10_flux(:)            ,&
      !$acc prod100_flux(:)           ,&
      !$acc crop_product_flux(:)     ,&
      !$acc  patch_to_soil_filter(:)  )
       !$acc enter data create(sum1, sum2, sum3, sum4 ,sum5)  
      !$acc parallel loop independent gang vector default(present) 
      do fp = 1, num_soilp_with_inactive
         ! initialize all the pft-level local flux arrays
         dwt_pool_seed         = 0.0_r8
         dwt_froot_to_litter(fp)    = 0.0_r8
         dwt_livecroot_to_litter(fp)= 0.0_r8
         dwt_deadcroot_to_litter(fp)= 0.0_r8
         conv_flux(fp)              = 0.0_r8
         prod10_flux(fp)            = 0.0_r8
         prod100_flux(fp)           = 0.0_r8
         crop_product_flux(fp)      = 0.0_r8
      enddo
      
      
      if(use_c13) then
         do fp = 1, num_soilp_with_inactive
            dwt_leafc13_seed(fp)           = 0._r8
            dwt_deadstemc13_seed(fp)       = 0._r8
            dwt_frootc13_to_litter(fp)     = 0._r8
            dwt_livecrootc13_to_litter(fp) = 0._r8
            dwt_deadcrootc13_to_litter(fp) = 0._r8
            conv_c13flux(fp)               = 0._r8
            prod10_c13flux(fp)             = 0._r8
            prod100_c13flux(fp)            = 0._r8
            crop_product_c13flux(fp)       = 0._r8
            
         enddo
      end if
      
      if ( use_c14 ) then
         do fp = 1, num_soilp_with_inactive
            dwt_leafc14_seed(fp)           = 0._r8
            dwt_deadstemc14_seed(fp)       = 0._r8
            dwt_frootc14_to_litter(fp)     = 0._r8
            dwt_livecrootc14_to_litter(fp) = 0._r8
            dwt_deadcrootc14_to_litter(fp) = 0._r8
            conv_c14flux(fp)               = 0._r8
            prod10_c14flux(fp)             = 0._r8
            prod100_c14flux(fp)            = 0._r8
            crop_product_c14flux(fp)       = 0._r8
         enddo
      endif
      
      !$acc parallel loop independent gang vector default(present) present(veg_cs,veg_ns,veg_ps) private(p,c,l,dwt)
      do fp = 1, num_soilp_with_inactive
         p = filter_soilp_with_inactive(fp)
         c = veg_pp%column(p)
         l = veg_pp%landunit(p)
         patch_to_soil_filter(p) = fp
         
         ! calculate the change in weight for the timestep
         dwt = veg_pp%wtcol(p)-prior_weights%pwtcol(p)
         cnstate_vars%lfpftd_patch(p) = -dwt
         
         ! Patches for which weight increases on this timestep
         if (dwt > 0._r8) then
            
            ! first identify Patches that are initiating on this timestep
            ! and set all the necessary state and flux variables
            if (prior_weights%pwtcol(p) == 0._r8) then
               
               ! set initial conditions for PFT that is being initiated
               ! in this time step.  Based on the settings in cnIniTimeVar.
               
               ! pft-level carbon state variables
               
               veg_cs%leafc(p)              = 0._r8
               veg_cs%leafc_storage(p)      = 0._r8
               veg_cs%leafc_xfer(p)         = 0._r8
               veg_cs%frootc(p)             = 0._r8
               veg_cs%frootc_storage(p)     = 0._r8
               veg_cs%frootc_xfer(p)        = 0._r8
               veg_cs%livestemc(p)          = 0._r8
               veg_cs%livestemc_storage(p)  = 0._r8
               veg_cs%livestemc_xfer(p)     = 0._r8
               veg_cs%deadstemc(p)          = 0._r8
               veg_cs%deadstemc_storage(p)  = 0._r8
               veg_cs%deadstemc_xfer(p)     = 0._r8
               veg_cs%livecrootc(p)         = 0._r8
               veg_cs%livecrootc_storage(p) = 0._r8
               veg_cs%livecrootc_xfer(p)    = 0._r8
               veg_cs%deadcrootc(p)         = 0._r8
               veg_cs%deadcrootc_storage(p) = 0._r8
               veg_cs%deadcrootc_xfer(p)    = 0._r8
               veg_cs%gresp_storage(p)      = 0._r8
               veg_cs%gresp_xfer(p)         = 0._r8
               veg_cs%cpool(p)              = 0._r8
               veg_cs%xsmrpool(p)           = 0._r8
               veg_cs%ctrunc(p)             = 0._r8
               veg_cs%dispvegc(p)           = 0._r8
               veg_cs%storvegc(p)           = 0._r8
               veg_cs%totvegc(p)            = 0._r8
               veg_cs%totpftc(p)            = 0._r8
               
               veg_ns%leafn(p)              = 0._r8
               veg_ns%leafn_storage(p)      = 0._r8
               veg_ns%leafn_xfer(p)         = 0._r8
               veg_ns%frootn(p)             = 0._r8
               veg_ns%frootn_storage(p)     = 0._r8
               veg_ns%frootn_xfer(p)        = 0._r8
               veg_ns%livestemn(p)          = 0._r8
               veg_ns%livestemn_storage(p)  = 0._r8
               veg_ns%livestemn_xfer(p)     = 0._r8
               veg_ns%deadstemn(p)          = 0._r8
               veg_ns%deadstemn_storage(p)  = 0._r8
               veg_ns%deadstemn_xfer(p)     = 0._r8
               veg_ns%livecrootn(p)         = 0._r8
               veg_ns%livecrootn_storage(p) = 0._r8
               veg_ns%livecrootn_xfer(p)    = 0._r8
               veg_ns%deadcrootn(p)         = 0._r8
               veg_ns%deadcrootn_storage(p) = 0._r8
               veg_ns%deadcrootn_xfer(p)    = 0._r8
               veg_ns%retransn(p)           = 0._r8
               veg_ns%npool(p)              = 0._r8
               veg_ns%ntrunc(p)             = 0._r8
               veg_ns%dispvegn(p)           = 0._r8
               veg_ns%storvegn(p)           = 0._r8
               veg_ns%totvegn(p)            = 0._r8
               veg_ns%totpftn(p)            = 0._r8
               
               veg_ps%leafp(p)              = 0._r8
               veg_ps%leafp_storage(p)      = 0._r8
               veg_ps%leafp_xfer(p)         = 0._r8
               veg_ps%frootp(p)             = 0._r8
               veg_ps%frootp_storage(p)     = 0._r8
               veg_ps%frootp_xfer(p)        = 0._r8
               veg_ps%livestemp(p)          = 0._r8
               veg_ps%livestemp_storage(p)  = 0._r8
               veg_ps%livestemp_xfer(p)     = 0._r8
               veg_ps%deadstemp(p)          = 0._r8
               veg_ps%deadstemp_storage(p)  = 0._r8
               veg_ps%deadstemp_xfer(p)     = 0._r8
               veg_ps%livecrootp(p)         = 0._r8
               veg_ps%livecrootp_storage(p) = 0._r8
               veg_ps%livecrootp_xfer(p)    = 0._r8
               veg_ps%deadcrootp(p)         = 0._r8
               veg_ps%deadcrootp_storage(p) = 0._r8
               veg_ps%deadcrootp_xfer(p)    = 0._r8
               veg_ps%retransp(p)           = 0._r8
               veg_ps%ppool(p)              = 0._r8
               veg_ps%ptrunc(p)             = 0._r8
               veg_ps%dispvegp(p)           = 0._r8
               veg_ps%storvegp(p)           = 0._r8
               veg_ps%totvegp(p)            = 0._r8
               veg_ps%totpftp (p)           = 0._r8
               
               canopystate_vars%laisun_patch(p) = 0._r8
               canopystate_vars%laisha_patch(p) = 0._r8
               
               cnstate_vars%dormant_flag_patch(p)          = 1._r8
               cnstate_vars%days_active_patch(p)           = 0._r8
               cnstate_vars%onset_flag_patch(p)            = 0._r8
               cnstate_vars%onset_counter_patch(p)         = 0._r8
               cnstate_vars%onset_gddflag_patch(p)         = 0._r8
               cnstate_vars%onset_fdd_patch(p)             = 0._r8
               cnstate_vars%onset_gdd_patch(p)             = 0._r8
               cnstate_vars%onset_swi_patch(p)             = 0._r8
               cnstate_vars%offset_flag_patch(p)           = 0._r8
               cnstate_vars%offset_counter_patch(p)        = 0._r8
               cnstate_vars%offset_fdd_patch(p)            = 0._r8
               cnstate_vars%offset_swi_patch(p)            = 0._r8
               cnstate_vars%lgsf_patch(p)                  = 0._r8
               cnstate_vars%bglfr_patch(p)                 = 0._r8
               cnstate_vars%bglfr_leaf_patch(p)            = 0._r8
               cnstate_vars%bglfr_froot_patch(p)           = 0._r8
               cnstate_vars%bgtr_patch(p)                  = 0._r8
               cnstate_vars%annavg_t2m_patch(p)            = cnstate_vars%annavg_t2m_col(c)
               cnstate_vars%tempavg_t2m_patch(p)           = 0._r8
               cnstate_vars%alloc_pnow_patch(p)            = 1._r8
               cnstate_vars%c_allometry_patch(p)           = 0._r8
               cnstate_vars%n_allometry_patch(p)           = 0._r8
               cnstate_vars%p_allometry_patch(p)           = 0._r8
               cnstate_vars%tempsum_potential_gpp_patch(p) = 0._r8
               cnstate_vars%annsum_potential_gpp_patch(p)  = 0._r8
               cnstate_vars%tempmax_retransn_patch(p)      = 0._r8
               cnstate_vars%annmax_retransn_patch(p)       = 0._r8
               cnstate_vars%downreg_patch(p)               = 0._r8
               
               cnstate_vars%tempmax_retransp_patch(p)      = 0._r8
               cnstate_vars%annmax_retransp_patch(p)       = 0._r8
               
               if ( use_c14 ) then
                  cnstate_vars%rc14_atm_patch(p) = c14ratio
                  cnstate_vars%rc14_atm_patch(p) = 0._r8
               endif
               veg_cf%xsmrpool_recover(p)      = 0._r8
               veg_cf%plant_calloc(p)          = 0._r8
               veg_cf%excess_cflux(p)          = 0._r8
               veg_cf%prev_leafc_to_litter(p)  = 0._r8
               veg_cf%prev_frootc_to_litter(p) = 0._r8
               veg_cf%availc(p)                = 0._r8
               veg_cf%gpp_before_downreg(p)    = 0._r8
               veg_cf%tempsum_npp(p)           = 0._r8
               veg_cf%annsum_npp(p)            = 0._r8
      
               if ( use_c13 ) then
                  veg_cf%xsmrpool_c13ratio(p) = c13ratio
               end if
               veg_nf%plant_ndemand(p)         = 0._r8
               veg_nf%avail_retransn(p)        = 0._r8
               veg_nf%plant_nalloc(p)          = 0._r8
               
               veg_pf%plant_pdemand(p)         = 0._r8
               veg_pf%avail_retransp(p)        = 0._r8
               veg_pf%plant_palloc(p)          = 0._r8
               
               ! if ( use_c13 ) then
               !    call CarbonStateVarsInit(c13_veg_cs, p)
               ! endif
               ! if ( use_c14 ) then
               !    call CarbonStateVarsInit(c14_veg_cs, p)
               ! endif
               
               ! add phosphorus related variables
               photosyns_vars%psnsun_patch(p) = 0._r8
               photosyns_vars%psnsha_patch(p) = 0._r8
               !if ( use_c13 ) then
               !   photosyns_vars%alphapsnsun_patch(p) = 0._r8
               !   photosyns_vars%alphapsnsha_patch(p) = 0._r8
               !   photosyns_vars%rc13_canair_patch(p) = 0._r8
               !   photosyns_vars%rc13_psnsun_patch(p) = 0._r8
               !   photosyns_vars%rc13_psnsha_patch(p) = 0._r8
               !   photosyns_vars%c13_psnsun_patch(p) = 0._r8
               !   photosyns_vars%c13_psnsha_patch(p) = 0._r8
               !   
               !endif
               !if ( use_c14 ) then
               !   photosyns_vars%c14_psnsun_patch(p) = 0._r8
               !   photosyns_vars%c14_psnsha_patch(p) = 0._r8
               !end if
               
            end if  ! end initialization of new pft
         end if  ! weight decreasing
      end do     ! patch loop
      
      call cpu_time(stopt) 
      write(iulog,*) "TIMING dyn_cnbal_patch::create and init ", (stopt-startt)*1.E+3,"ms"
      
      call cpu_time(startt) 

      !$acc parallel loop independent gang vector default(present) private(p,c,l,g) &
      !$acc present(conv_flux(:),dwt_froot_to_litter(:), &
      !$acc dwt_livecroot_to_litter(:),dwt_deadcroot_to_litter(:),prod10_flux(:), prod100_flux(:), &
      !$acc crop_product_flux(:)      )  
      do fp = 1, num_soilp_with_inactive
         p = filter_soilp_with_inactive(fp)
         c = veg_pp%column(p)
         l = veg_pp%landunit(p)
         g = veg_pp%gridcell(p) 
         dwt_leaf_seed   = 0.0_r8
         dwt_deadstem_seed = 0.0_r8

         
         call dyn_veg_cs_Adjustments(    &
            l, c, p,        &
            prior_weights,                 &
            patch_state_updater,           &
            dwt_leaf_seed,             &
            dwt_deadstem_seed,         &
            conv_flux(fp),                 &
            dwt_froot_to_litter(fp),       &
            dwt_livecroot_to_litter(fp),   &
            dwt_deadcroot_to_litter(fp),   &
            prod10_flux(fp),               &
            prod100_flux(fp),              &
            crop_product_flux(fp),         &
            veg_cs                         &
            )
         
         veg_cf%dwt_seedc_to_leaf(p) = dwt_leaf_seed/dt
         veg_cf%dwt_seedc_to_deadstem(p) = dwt_deadstem_seed/dt
         veg_cf%dwt_conv_cflux(p) = -conv_flux(fp)/dt
         veg_cf%dwt_prod10c_gain(p) = -prod10_flux(fp)/dt
         veg_cf%dwt_prod100c_gain(p) = - prod100_flux(fp)/dt
         veg_cf%dwt_crop_productc_gain(p) = - crop_product_flux(fp)/dt
         
         ! Note that patch-level fluxes are stored per unit GRIDCELL area - thus, we don't
         ! need to multiply by the patch's gridcell weight when translating patch-level
         ! fluxes into gridcell-level fluxes.
         !$acc atomic update
         grc_cf%dwt_seedc_to_leaf(g) =  grc_cf%dwt_seedc_to_leaf(g) + veg_cf%dwt_seedc_to_leaf(p)
         !$acc end atomic 
         
         !$acc atomic update 
         grc_cf%dwt_seedc_to_deadstem(g) = grc_cf%dwt_seedc_to_deadstem(g) +  veg_cf%dwt_seedc_to_deadstem(p)
         !$acc end atomic 
         
         !$acc atomic update
         grc_cf%dwt_conv_cflux(g) = grc_cf%dwt_conv_cflux(g) +  veg_cf%dwt_conv_cflux(p)
         !$acc end atomic
         
         !$acc atomic update 
         grc_cf%dwt_prod10c_gain(g) = grc_cf%dwt_prod10c_gain(g) + veg_cf%dwt_prod10c_gain(p)
         !$acc end atomic
         
         !$acc atomic update 
         grc_cf%dwt_prod100c_gain(g) = grc_cf%dwt_prod100c_gain(g) + veg_cf%dwt_prod100c_gain(p)
         !$acc end atomic 
      end do
      call cpu_time(stopt) 
      write(iulog,*) "dyn_cnbal_patch::veg_cs_Adjustment ",(stopt-startt)*1.E+3, "ms"

      call cpu_time(startt) 
      ! calculate pft-to-column for fluxes into litter and CWD pools 
      !$acc parallel loop independent gang worker collapse(2) default(present) private(sum1,sum2,sum3,sum4,sum5,c) &
      !$acc present(dwt_froot_to_litter(:),dwt_livecroot_to_litter(:),dwt_deadcroot_to_litter(:)) 
      do j = 1, nlevdecomp
         do fc = 1, num_soilc_with_inactive
            c = filter_soilc_with_inactive(fc)
            l = col_pp%landunit(c) 
            
            sum1 = 0.0_r8; sum2 = 0.0_r8 
            sum3 = 0.0_r8; sum4 = 0.0_r8;
            sum5 = 0.0_r8;
            !
            !$acc loop vector reduction(+:sum1,sum2,sum3,sum4,sum5) private(froot,croot,fr_flab,fr_fcel,fr_flig,fp)
            do p = col_pp%pfti(c), col_pp%pftf(c) 
               if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then
                  fp = patch_to_soil_filter(p)
                  froot   = cnstate_vars%froot_prof_patch(p,j)
                  croot   = cnstate_vars%croot_prof_patch(p,j)
                  fr_flab = veg_vp%fr_flab(veg_pp%itype(p))
                  fr_fcel = veg_vp%fr_fcel(veg_pp%itype(p))
                  fr_flig = veg_vp%fr_flig(veg_pp%itype(p))
                  
                  ! fine root litter carbon fluxes
                  sum1 = sum1 + (dwt_froot_to_litter(fp)* fr_flab)/dt * froot
                  sum2 = sum2 + (dwt_froot_to_litter(fp)* fr_fcel)/dt * froot
                  sum3 = sum3 + (dwt_froot_to_litter(fp)* fr_flig)/dt * froot
                  !
                  sum4 = sum4 + (dwt_livecroot_to_litter(fp))/dt * croot
                  sum5 = sum5 + (dwt_deadcroot_to_litter(fp))/dt * croot
               end if 
            end do 
            ! 
            col_cf%dwt_frootc_to_litr_met_c(c,j) = col_cf%dwt_frootc_to_litr_met_c(c,j) + sum1 
            col_cf%dwt_frootc_to_litr_cel_c(c,j) = col_cf%dwt_frootc_to_litr_cel_c(c,j) + sum2 
            col_cf%dwt_frootc_to_litr_lig_c(c,j) = col_cf%dwt_frootc_to_litr_lig_c(c,j) + sum3 
            !
            col_cf%dwt_livecrootc_to_cwdc(c,j) = col_cf%dwt_livecrootc_to_cwdc(c,j) + sum4
            !
            col_cf%dwt_deadcrootc_to_cwdc(c,j) = col_cf%dwt_deadcrootc_to_cwdc(c,j) + sum5 
         end do 
      end do 
      
      !$acc parallel loop independent gang worker default(present) private(sum1,sum2,sum3,sum4,sum5) &
      !$acc present(dwt_froot_to_litter(:),dwt_livecroot_to_litter(:),dwt_deadcroot_to_litter(:),&
      !$acc  conv_flux(:), crop_product_flux(:), prod10_flux(:),prod100_flux(:)) 
      do fc = 1, num_soilc_with_inactive 
         c = filter_soilc_with_inactive(fc) 
         l = col_pp%landunit(c) 
         sum1 = 0.0_r8; sum2 = 0.0_r8; 
         sum3 = 0.0_r8; sum4 = 0.0_r8 
         sum5 = 0.0_r8;
         
         !$acc loop vector reduction(+:sum1,sum2,sum3,sum4,sum5) private(fp)
         do p = col_pp%pfti(c), col_pp%pftf(c) 
            if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then
               fp = patch_to_soil_filter(p)
               ! column-level fluxes are accumulated as positive fluxes.
               ! column-level C flux updates
               sum1 = sum1 - prod10_flux(fp)/dt
               sum2 = sum2 - prod100_flux(fp)/dt
               sum3 = sum3 - crop_product_flux(fp)/dt
               sum4 = sum4 - conv_flux(fp)/dt
               sum5 = sum5 + dwt_froot_to_litter(fp)/dt + dwt_livecroot_to_litter(fp)/dt + dwt_deadcroot_to_litter(fp)/dt
               
            endif 
         end do 
         col_cf%dwt_prod10c_gain(c) = col_cf%dwt_prod10c_gain(c) + sum1
         col_cf%dwt_prod100c_gain(c)= col_cf%dwt_prod100c_gain(c)+ sum2
         col_cf%dwt_crop_productc_gain(c) = col_cf%dwt_crop_productc_gain(c) + sum3 
         col_cf%dwt_conv_cflux(c) = col_cf%dwt_conv_cflux(c) + sum4 
         col_cf%dwt_slash_cflux(c) = col_cf%dwt_slash_cflux(c) + sum5 
      end do 
      
      call cpu_time(stopt) 
      write(iulog,*) "cnbal_patch::col_cf reductions ",(stopt-startt)*1.E+3,"ms"
      
      call cpu_time(startt) 
      !$acc parallel loop independent gang vector default(present)
      do fp = 1, num_soilp_with_inactive
         p = filter_soilp_with_inactive(fp)
         c = veg_pp%column(p)
         l = veg_pp%landunit(p)
         dwt_leaf_seed   = 0.0_r8
         dwt_deadstem_seed = 0.0_r8
         dwt_pool_seed = 0._r8
         call dyn_veg_ns_Adjustments(    &
            l,c,p,              &
            prior_weights,                 &
            patch_state_updater,           &
            dwt_leaf_seed,             &
            dwt_deadstem_seed,         &
            dwt_pool_seed,             &
            conv_flux(fp),                 &
            dwt_froot_to_litter(fp),       &
            dwt_livecroot_to_litter(fp),   &
            dwt_deadcroot_to_litter(fp),   &
            prod10_flux(fp),               &
            prod100_flux(fp),              &
            crop_product_flux(fp),         &
            veg_ns                         &
            )
         
         veg_nf%dwt_seedn_to_leaf(p)   = dwt_leaf_seed/dt
         veg_nf%dwt_seedn_to_npool(p) = dwt_pool_seed/dt
         veg_nf%dwt_seedn_to_deadstem(p) = dwt_deadstem_seed/dt
         veg_nf%dwt_conv_nflux(p) = -conv_flux(fp)/dt
         veg_nf%dwt_prod10n_gain(p) = -prod10_flux(fp)/dt
         veg_nf%dwt_prod100n_gain(p)= -prod100_flux(fp)/dt
         veg_nf%dwt_crop_productn_gain(p) = -crop_product_flux(fp)/dt
         ! N fluxes
         
         !$acc atomic update 
         grc_nf%dwt_seedn_to_leaf(g) = grc_nf%dwt_seedn_to_leaf(g) + veg_nf%dwt_seedn_to_leaf(p)
         !$acc end atomic 
         
         !$acc atomic update 
         grc_nf%dwt_seedn_to_deadstem(g) = grc_nf%dwt_seedn_to_deadstem(g) + veg_nf%dwt_seedn_to_deadstem(p)
         !$acc end atomic 
         
         !$acc atomic update 
         grc_nf%dwt_seedn_to_npool(g) = grc_nf%dwt_seedn_to_npool(g) + veg_nf%dwt_seedn_to_npool(p)
         !$acc end atomic 
         
         !$acc atomic update 
         grc_nf%dwt_conv_nflux(g) =  grc_nf%dwt_conv_nflux(g) + veg_nf%dwt_conv_nflux(p)
         !$acc end atomic 
         
         !$acc atomic update 
         grc_nf%dwt_prod10n_gain(g) = grc_nf%dwt_prod10n_gain(g) + veg_nf%dwt_prod10n_gain(p)
         !$acc end atomic
         
         !$acc atomic update 
         grc_nf%dwt_prod100n_gain(g) = grc_nf%dwt_prod100n_gain(g) + veg_nf%dwt_prod100n_gain(p)
         !$acc end atomic 
      enddo
      call cpu_time(stopt) 
      write(iulog,*) "dyn_cnbal_patch::veg_ns_adjustments ",(stopt-startt)*1.E+3,"ms"
      
      call cpu_time(startt) 
      ! calculate pft-to-column for fluxes into litter and CWD pools 
      !$acc parallel loop independent gang worker collapse(2) default(present) private(sum1,sum2,sum3,sum4,sum5,c,l) &
      !$acc present(dwt_froot_to_litter(:),dwt_livecroot_to_litter(:),dwt_deadcroot_to_litter(:)) 
      do j = 1, nlevdecomp
         do fc = 1, num_soilc_with_inactive
            c = filter_soilc_with_inactive(fc)
            l = col_pp%landunit(c) 
            
            sum1 = 0.0_r8
            sum2 = 0.0_r8 
            sum3 = 0.0_r8 
            sum4 = 0.0_r8 
            sum5 = 0.0_r8 
            
            !$acc loop vector reduction(+:sum1,sum2,sum3,sum4,sum5) private(froot,croot,fr_flab,fr_fcel,fr_flig,fp)
            do p = col_pp%pfti(c), col_pp%pftf(c) 
               if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then
                  fp = patch_to_soil_filter(p)

                  froot   = cnstate_vars%froot_prof_patch(p,j)
                  croot   = cnstate_vars%croot_prof_patch(p,j)
                  fr_flab = veg_vp%fr_flab(veg_pp%itype(p))
                  fr_fcel = veg_vp%fr_fcel(veg_pp%itype(p))
                  fr_flig = veg_vp%fr_flig(veg_pp%itype(p))
                  
                  ! fine root litter carbon fluxes
                  sum1 = sum1 + (dwt_froot_to_litter(fp)* fr_flab)/dt * froot
                  sum2 = sum2 + (dwt_froot_to_litter(fp)* fr_fcel)/dt * froot
                  sum3 = sum3 + (dwt_froot_to_litter(fp)* fr_flig)/dt * froot
                  
                  ! livecroot fluxes to cwd
                  sum4 = sum4 + (dwt_livecroot_to_litter(fp))/dt * croot
                  ! 
                  sum5 = sum5 + (dwt_deadcroot_to_litter(fp))/dt * croot
               end if 
            end do 
            ! 
            col_nf%dwt_frootn_to_litr_met_n(c,j) = col_nf%dwt_frootn_to_litr_met_n(c,j) + sum1 
            col_nf%dwt_frootn_to_litr_cel_n(c,j) = col_nf%dwt_frootn_to_litr_cel_n(c,j) + sum2 
            col_nf%dwt_frootn_to_litr_lig_n(c,j) = col_nf%dwt_frootn_to_litr_lig_n(c,j) + sum3 
            ! 
            col_nf%dwt_livecrootn_to_cwdn(c,j)   = col_nf%dwt_livecrootn_to_cwdn(c,j) + sum4 
            !
            col_nf%dwt_deadcrootn_to_cwdn(c,j) = col_nf%dwt_deadcrootn_to_cwdn(c,j) + sum5 
         end do 
      end do 
      
      !$acc parallel loop independent gang worker default(present) private(c,l,sum1,sum2,sum3,sum4,sum5) &
      !$acc present(dwt_froot_to_litter(:),dwt_livecroot_to_litter(:),dwt_deadcroot_to_litter(:),&
      !$acc  conv_flux(:), crop_product_flux(:), prod10_flux(:),prod100_flux(:))
      do fc = 1, num_soilc_with_inactive 
         c = filter_soilc_with_inactive(fc) 
         l = col_pp%landunit(c) 
         sum1 = 0.0_r8; sum2 = 0.0_r8; 
         sum3 = 0.0_r8; sum4 = 0.0_r8 
         sum5 = 0.0_r8;
         
         !$acc loop vector reduction(+:sum1,sum2,sum3,sum4,sum5) private(fp)
         do p = col_pp%pfti(c), col_pp%pftf(c) 
            if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then
               ! column-level fluxes are accumulated as positive fluxes.
               ! column-level C flux updates
               fp = patch_to_soil_filter(p)

               sum1 = sum1 - prod10_flux(fp)/dt
               sum2 = sum2 - prod100_flux(fp)/dt
               sum3 = sum3 - crop_product_flux(fp)/dt
               sum4 = sum4 - conv_flux(fp)/dt
               sum5 = sum5 + dwt_froot_to_litter(fp)/dt + &
               dwt_livecroot_to_litter(fp)/dt + &
               dwt_deadcroot_to_litter(fp)/dt
               
            endif 
         end do 
         col_nf%dwt_prod10n_gain(c) = col_nf%dwt_prod10n_gain(c) + sum1
         col_nf%dwt_prod100n_gain(c)= col_nf%dwt_prod100n_gain(c)+ sum2
         col_nf%dwt_crop_productn_gain(c) = col_nf%dwt_crop_productn_gain(c) + sum3 
         col_nf%dwt_conv_nflux(c) = col_nf%dwt_conv_nflux(c) + sum4 
         col_nf%dwt_slash_nflux(c) = col_nf%dwt_slash_nflux(c) + sum5 
      end do 
      
      call cpu_time(stopt) 
      write(iulog,*) "cnbal_patch::col_nf reductions ",(stopt-startt)*1.E+3,"ms"
      
      call cpu_time(startt) 
      !$acc parallel loop independent gang vector default(present) private(p,c,l)
      do fp = 1, num_soilp_with_inactive
         p = filter_soilp_with_inactive(fp)
         c = veg_pp%column(p)
         l = veg_pp%landunit(p)
         dwt_leaf_seed   = 0.0_r8
         dwt_deadstem_seed = 0.0_r8
         dwt_pool_seed = 0._r8
         
         call dyn_veg_ps_Adjustments( &
         l,c,p,                  &
         prior_weights,                &
         patch_state_updater,          &
         dwt_leaf_seed,            &
         dwt_deadstem_seed,        &
         dwt_pool_seed,            &
         conv_flux(fp),                &
         dwt_froot_to_litter(fp),      &
         dwt_livecroot_to_litter(fp),  &
         dwt_deadcroot_to_litter(fp),  &
         prod10_flux(fp),              &
         prod100_flux(fp),             &
         crop_product_flux(fp),        &
         veg_ps                        &
         )
         
         ! P fluxes
         veg_pf%dwt_seedp_to_leaf(p)   = dwt_leaf_seed/dt
         veg_pf%dwt_seedp_to_ppool(p) = dwt_pool_seed/dt
         veg_pf%dwt_seedp_to_deadstem(p) = dwt_deadstem_seed/dt
         veg_pf%dwt_conv_pflux(p) = -conv_flux(fp)/dt
         veg_pf%dwt_prod10p_gain(p)  = -prod10_flux(fp)/dt
         veg_pf%dwt_prod100p_gain(p) = -prod100_flux(fp)/dt
         veg_pf%dwt_crop_productp_gain(p) = -crop_product_flux(fp)/dt
         
         !$acc atomic update 
         grc_pf%dwt_seedp_to_leaf(g) = grc_pf%dwt_seedp_to_leaf(g) + veg_pf%dwt_seedp_to_leaf(p)
         !$acc end atomic 
         
         !$acc atomic update 
         grc_pf%dwt_seedp_to_deadstem(g) = grc_pf%dwt_seedp_to_deadstem(g) + veg_pf%dwt_seedp_to_deadstem(p)
         !$acc end atomic 
         
         !$acc atomic update 
         grc_pf%dwt_seedp_to_ppool(g) = grc_pf%dwt_seedp_to_ppool(g) + veg_pf%dwt_seedp_to_ppool(p)
         !$acc end atomic 
         
         !$acc atomic update 
         grc_pf%dwt_conv_pflux(g) =  grc_pf%dwt_conv_pflux(g) + veg_pf%dwt_conv_pflux(p)
         !$acc end atomic 
         
         !$acc atomic update 
         grc_pf%dwt_prod10p_gain(g) = grc_pf%dwt_prod10p_gain(g) + veg_pf%dwt_prod10p_gain(p)
         !$acc end atomic
         
         !$acc atomic update 
         grc_pf%dwt_prod100p_gain(g) = grc_pf%dwt_prod100p_gain(g) + veg_pf%dwt_prod100p_gain(p)
         !$acc end atomic 
      end do

      call cpu_time(stopt) 
      write(iulog,*) "dyn_cnbal_patch::veg_ps_adjustments ",(stopt-startt)*1.E+3,"ms"
      
      call cpu_time(startt) 
      ! calculate pft-to-column for fluxes into litter and CWD pools 
      !$acc parallel loop independent gang worker collapse(2) default(present) private(sum1,sum2,sum3,sum4,sum5,c,l) &
      !$acc present(dwt_froot_to_litter(:),dwt_livecroot_to_litter(:),dwt_deadcroot_to_litter(:) )
      do j = 1, nlevdecomp
         do fc = 1, num_soilc_with_inactive
            c = filter_soilc_with_inactive(fc)
            l = col_pp%landunit(c) 
            
            ! c = veg_pp%column(p)
            sum1 = 0.0_r8; sum2 = 0.0_r8 
            sum3 = 0.0_r8; sum4 = 0.0_r8 
            sum5 = 0.0_r8 
            
            !$acc loop vector reduction(+:sum1,sum2,sum3,sum4,sum5) private(froot,croot,fr_flab,fr_fcel,fr_flig,fp)
            do p = col_pp%pfti(c), col_pp%pftf(c) 
               if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then
                  fp = patch_to_soil_filter(p)
                  froot   = cnstate_vars%froot_prof_patch(p,j)
                  croot   = cnstate_vars%croot_prof_patch(p,j)
                  fr_flab = veg_vp%fr_flab(veg_pp%itype(p))
                  fr_fcel = veg_vp%fr_fcel(veg_pp%itype(p))
                  fr_flig = veg_vp%fr_flig(veg_pp%itype(p))
                  
                  ! fine root litter carbon fluxes
                  sum1 = sum1 + (dwt_froot_to_litter(fp)* fr_flab)/dt * froot
                  sum2 = sum2 + (dwt_froot_to_litter(fp)* fr_fcel)/dt * froot
                  sum3 = sum3 + (dwt_froot_to_litter(fp)* fr_flig)/dt * froot
                  
                  ! livecroot fluxes to cwd
                  sum4 = sum4 + (dwt_livecroot_to_litter(fp))/dt * croot
                  ! 
                  sum5 = sum5 + (dwt_deadcroot_to_litter(fp))/dt * croot
               end if 
            end do 
            ! 
            col_pf%dwt_frootp_to_litr_met_p(c,j) = col_pf%dwt_frootp_to_litr_met_p(c,j) + sum1 
            col_pf%dwt_frootp_to_litr_cel_p(c,j) = col_pf%dwt_frootp_to_litr_cel_p(c,j) + sum2 
            col_pf%dwt_frootp_to_litr_lig_p(c,j) = col_pf%dwt_frootp_to_litr_lig_p(c,j) + sum3 
            ! 
            col_pf%dwt_livecrootp_to_cwdp(c,j) = col_pf%dwt_livecrootp_to_cwdp(c,j) + sum4 
            !
            col_pf%dwt_deadcrootp_to_cwdp(c,j) = col_pf%dwt_deadcrootp_to_cwdp(c,j) + sum5 
         end do 
      end do 
      
      !$acc parallel loop independent gang worker default(present) private(sum1,sum2,sum3,sum4,sum5,c,l) &
      !$acc present(dwt_froot_to_litter(:),dwt_livecroot_to_litter(:),dwt_deadcroot_to_litter(:),&
      !$acc  conv_flux(:), crop_product_flux(:), prod10_flux(:),prod100_flux(:))
      do fc = 1, num_soilc_with_inactive 
         c = filter_soilc_with_inactive(fc) 
         l = col_pp%landunit(c) 
         sum1 = 0.0_r8; sum2 = 0.0_r8; 
         sum3 = 0.0_r8; sum4 = 0.0_r8 
         sum5 = 0.0_r8;
         
         !$acc loop vector reduction(+:sum1,sum2,sum3,sum4,sum5) private(fp)
         do p = col_pp%pfti(c), col_pp%pftf(c) 
            if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then
               ! column-level fluxes are accumulated as positive fluxes.
               ! column-level C flux updates
               fp = patch_to_soil_filter(p)
               
               sum1 = sum1 - prod10_flux(fp)/dt
               sum2 = sum2 - prod100_flux(fp)/dt
               sum3 = sum3 - crop_product_flux(fp)/dt
               sum4 = sum4 - conv_flux(fp)/dt
               sum5 = sum5 + dwt_froot_to_litter(fp)/dt + &
                           dwt_livecroot_to_litter(fp)/dt + &
                           dwt_deadcroot_to_litter(fp)/dt
               
            endif 
         end do 
         col_pf%dwt_prod10p_gain(c) = col_pf%dwt_prod10p_gain(c) + sum1
         col_pf%dwt_prod100p_gain(c)= col_pf%dwt_prod100p_gain(c)+ sum2
         col_pf%dwt_crop_productp_gain(c) = col_pf%dwt_crop_productp_gain(c) + sum3 
         col_pf%dwt_conv_pflux(c) = col_pf%dwt_conv_pflux(c) + sum4 
         col_pf%dwt_slash_pflux(c) = col_pf%dwt_slash_pflux(c) + sum5 
      end do 
      
      call cpu_time(stopt) 
      write(iulog,*) "cnbal_patch::col_pf reductions ",(stopt-startt)*1.E+3,"ms"
      
      
      !$acc exit data delete(dwt_leaf_seed,&
      !$acc dwt_deadstem_seed , &
      !$acc dwt_pool_seed     , &
      !$acc dwt_froot_to_litter(:) ,&
      !$acc dwt_livecroot_to_litter(:),&
      !$acc dwt_deadcroot_to_litter(:),&
      !$acc conv_flux(:)              ,&
      !$acc prod10_flux(:)            ,&
      !$acc prod100_flux(:)           ,&
      !$acc crop_product_flux(:)      ,&
      !$acc  patch_to_soil_filter(:)  )
      !$acc exit data delete(sum1,sum2,sum3,sum4,sum5)  
      ! Deallocate pft-level flux arrays
      if ( use_c13 ) then
         deallocate(dwt_leafc13_seed)
         deallocate(dwt_deadstemc13_seed)
         deallocate(dwt_frootc13_to_litter)
         deallocate(dwt_livecrootc13_to_litter)
         deallocate(dwt_deadcrootc13_to_litter)
         deallocate(conv_c13flux)
         deallocate(prod10_c13flux)
         deallocate(prod100_c13flux)
         deallocate(crop_product_c13flux)
      endif
      
      if ( use_c14 ) then
         deallocate(dwt_leafc14_seed)
         deallocate(dwt_deadstemc14_seed)
         deallocate(dwt_frootc14_to_litter)
         deallocate(dwt_livecrootc14_to_litter)
         deallocate(dwt_deadcrootc14_to_litter)
         deallocate(conv_c14flux)
         deallocate(prod10_c14flux)
         deallocate(prod100_c14flux)
         deallocate(crop_product_c14flux)
      endif
      
   end subroutine dyn_cnbal_patch
   
   !-----------------------------------------------------------------------
   subroutine CarbonStateVarsInit(cs, p)
      !
      ! !DESCRIPTION:
      ! Initializes p-th patch of carbonstate_type
      !
      !$acc routine seq 
      implicit none
      !
      ! !ARGUMENT
      type(vegetation_carbon_state), intent(inout) :: cs
      integer                    , intent(in)    :: p
      
      cs%leafc(p)              = 0._r8
      cs%leafc_storage(p)      = 0._r8
      cs%leafc_xfer(p)         = 0._r8
      cs%frootc(p)             = 0._r8
      cs%frootc_storage(p)     = 0._r8
      cs%frootc_xfer(p)        = 0._r8
      cs%livestemc(p)          = 0._r8
      cs%livestemc_storage(p)  = 0._r8
      cs%livestemc_xfer(p)     = 0._r8
      cs%deadstemc(p)          = 0._r8
      cs%deadstemc_storage(p)  = 0._r8
      cs%deadstemc_xfer(p)     = 0._r8
      cs%livecrootc(p)         = 0._r8
      cs%livecrootc_storage(p) = 0._r8
      cs%livecrootc_xfer(p)    = 0._r8
      cs%deadcrootc(p)         = 0._r8
      cs%deadcrootc_storage(p) = 0._r8
      cs%deadcrootc_xfer(p)    = 0._r8
      cs%gresp_storage(p)      = 0._r8
      cs%gresp_xfer(p)         = 0._r8
      cs%cpool(p)              = 0._r8
      cs%xsmrpool(p)           = 0._r8
      cs%ctrunc(p)             = 0._r8
      cs%dispvegc(p)           = 0._r8
      cs%storvegc(p)           = 0._r8
      cs%totvegc(p)            = 0._r8
      cs%totpftc(p)            = 0._r8
      
   end subroutine CarbonStateVarsInit
   
   !-----------------------------------------------------------------------
   subroutine NitrogenStateVarsInit(veg_ns, p)
      !
      ! !DESCRIPTION:
      ! Initializes p-th patch of nitrogenstate_type
      !
      !$acc routine seq 
      implicit none
      !
      ! !ARGUMENT
      type(vegetation_nitrogen_state), intent(inout) :: veg_ns
      integer                 , intent(in)    :: p
      
      veg_ns%leafn(p)              = 0._r8
      veg_ns%leafn_storage(p)      = 0._r8
      veg_ns%leafn_xfer(p)         = 0._r8
      veg_ns%frootn(p)             = 0._r8
      veg_ns%frootn_storage(p)     = 0._r8
      veg_ns%frootn_xfer(p)        = 0._r8
      veg_ns%livestemn(p)          = 0._r8
      veg_ns%livestemn_storage(p)  = 0._r8
      veg_ns%livestemn_xfer(p)     = 0._r8
      veg_ns%deadstemn(p)          = 0._r8
      veg_ns%deadstemn_storage(p)  = 0._r8
      veg_ns%deadstemn_xfer(p)     = 0._r8
      veg_ns%livecrootn(p)         = 0._r8
      veg_ns%livecrootn_storage(p) = 0._r8
      veg_ns%livecrootn_xfer(p)    = 0._r8
      veg_ns%deadcrootn(p)         = 0._r8
      veg_ns%deadcrootn_storage(p) = 0._r8
      veg_ns%deadcrootn_xfer(p)    = 0._r8
      veg_ns%retransn(p)           = 0._r8
      veg_ns%npool(p)              = 0._r8
      veg_ns%ntrunc(p)             = 0._r8
      veg_ns%dispvegn(p)           = 0._r8
      veg_ns%storvegn(p)           = 0._r8
      veg_ns%totvegn(p)            = 0._r8
      veg_ns%totpftn(p)            = 0._r8
      
   end subroutine NitrogenStateVarsInit
   
   !-----------------------------------------------------------------------
   subroutine PhosphorusStateVarsInit(veg_ps, p)
      !
      ! !DESCRIPTION:
      ! Initializes p-th patch of phosphorusstate_type
      !
      !$acc routine seq 
      implicit none
      !
      ! !ARGUMENT
      type(vegetation_phosphorus_state), intent(inout) :: veg_ps
      integer                   , intent(in)    :: p
      
      veg_ps%leafp(p)              = 0._r8
      veg_ps%leafp_storage(p)      = 0._r8
      veg_ps%leafp_xfer(p)         = 0._r8
      veg_ps%frootp(p)             = 0._r8
      veg_ps%frootp_storage(p)     = 0._r8
      veg_ps%frootp_xfer(p)        = 0._r8
      veg_ps%livestemp(p)          = 0._r8
      veg_ps%livestemp_storage(p)  = 0._r8
      veg_ps%livestemp_xfer(p)     = 0._r8
      veg_ps%deadstemp(p)          = 0._r8
      veg_ps%deadstemp_storage(p)  = 0._r8
      veg_ps%deadstemp_xfer(p)     = 0._r8
      veg_ps%livecrootp(p)         = 0._r8
      veg_ps%livecrootp_storage(p) = 0._r8
      veg_ps%livecrootp_xfer(p)    = 0._r8
      veg_ps%deadcrootp(p)         = 0._r8
      veg_ps%deadcrootp_storage(p) = 0._r8
      veg_ps%deadcrootp_xfer(p)    = 0._r8
      veg_ps%retransp(p)           = 0._r8
      veg_ps%ppool(p)              = 0._r8
      veg_ps%ptrunc(p)             = 0._r8
      veg_ps%dispvegp(p)           = 0._r8
      veg_ps%storvegp(p)           = 0._r8
      veg_ps%totvegp(p)            = 0._r8
      veg_ps%totpftp (p)           = 0._r8
      
   end subroutine PhosphorusStateVarsInit
   
   !-----------------------------------------------------------------------
   subroutine CanopyStateVarsInit(canopystate_vars, p)
      !
      ! !DESCRIPTION:
      ! Initializes p-th patch of canopystate_type
      !
      !$acc routine seq 
      implicit none
      !
      ! !ARGUMENT
      type(canopystate_type), intent(inout) :: canopystate_vars
      integer               , intent(in)    :: p
      
      canopystate_vars%laisun_patch(p) = 0._r8
      canopystate_vars%laisha_patch(p) = 0._r8
      
   end subroutine CanopyStateVarsInit
   
   !-----------------------------------------------------------------------
   subroutine CNStateVarsInit(cnstate_vars, p, c)
      !
      ! !DESCRIPTION:
      ! Initializes p-th patch of cnstate_type
      !
      !$acc routine seq 
      use elm_varcon, only : c14ratio
      implicit none
      !
      ! !ARGUMENT
      type(cnstate_type), intent(inout) :: cnstate_vars
      integer           , intent(in)    :: p
      integer           , intent(in)    :: c
      
      cnstate_vars%dormant_flag_patch(p)          = 1._r8
      cnstate_vars%days_active_patch(p)           = 0._r8
      cnstate_vars%onset_flag_patch(p)            = 0._r8
      cnstate_vars%onset_counter_patch(p)         = 0._r8
      cnstate_vars%onset_gddflag_patch(p)         = 0._r8
      cnstate_vars%onset_fdd_patch(p)             = 0._r8
      cnstate_vars%onset_gdd_patch(p)             = 0._r8
      cnstate_vars%onset_swi_patch(p)             = 0._r8
      cnstate_vars%offset_flag_patch(p)           = 0._r8
      cnstate_vars%offset_counter_patch(p)        = 0._r8
      cnstate_vars%offset_fdd_patch(p)            = 0._r8
      cnstate_vars%offset_swi_patch(p)            = 0._r8
      cnstate_vars%lgsf_patch(p)                  = 0._r8
      cnstate_vars%bglfr_patch(p)                 = 0._r8
      cnstate_vars%bglfr_leaf_patch(p)            = 0._r8
      cnstate_vars%bglfr_froot_patch(p)           = 0._r8
      cnstate_vars%bgtr_patch(p)                  = 0._r8
      cnstate_vars%annavg_t2m_patch(p)            = cnstate_vars%annavg_t2m_col(c)
      cnstate_vars%tempavg_t2m_patch(p)           = 0._r8
      cnstate_vars%alloc_pnow_patch(p)            = 1._r8
      cnstate_vars%c_allometry_patch(p)           = 0._r8
      cnstate_vars%n_allometry_patch(p)           = 0._r8
      cnstate_vars%p_allometry_patch(p)           = 0._r8
      cnstate_vars%tempsum_potential_gpp_patch(p) = 0._r8
      cnstate_vars%annsum_potential_gpp_patch(p)  = 0._r8
      cnstate_vars%tempmax_retransn_patch(p)      = 0._r8
      cnstate_vars%annmax_retransn_patch(p)       = 0._r8
      cnstate_vars%downreg_patch(p)               = 0._r8
      
      cnstate_vars%tempmax_retransp_patch(p)      = 0._r8
      cnstate_vars%annmax_retransp_patch(p)       = 0._r8
      
      if ( use_c14 ) then
         cnstate_vars%rc14_atm_patch(p) = c14ratio
         cnstate_vars%rc14_atm_patch(p) = 0._r8
      endif
      
   end subroutine CNStateVarsInit
   
   !-----------------------------------------------------------------------
   subroutine CarbonFluxVarsInit(cf, p)
      !
      ! !DESCRIPTION:
      ! Initializes p-th patch of carbonflux_type
      !
      !$acc routine seq 
      use elm_varcon, only : c13ratio
      !
      implicit none
      !
      ! !ARGUMENT
      type(vegetation_carbon_flux), intent(inout) :: cf
      integer              , intent(in)    :: p
      
      cf%xsmrpool_recover(p)      = 0._r8
      cf%plant_calloc(p)          = 0._r8
      cf%excess_cflux(p)          = 0._r8
      cf%prev_leafc_to_litter(p)  = 0._r8
      cf%prev_frootc_to_litter(p) = 0._r8
      cf%availc(p)                = 0._r8
      cf%gpp_before_downreg(p)    = 0._r8
      cf%tempsum_npp(p)           = 0._r8
      cf%annsum_npp(p)            = 0._r8
      
      if ( use_c13 ) then
         cf%xsmrpool_c13ratio(p) = c13ratio
      end if
      
   end subroutine CarbonFluxVarsInit
   
   !-----------------------------------------------------------------------
   subroutine NitrogenFluxVarsInit(p)
      !
      ! !DESCRIPTION:
      ! Initializes p-th patch of nitrogenflux_type
      !
      !$acc routine seq 
      implicit none
      !
      ! !ARGUMENT
      integer                , intent(in)    :: p
      
      veg_nf%plant_ndemand(p)         = 0._r8
      veg_nf%avail_retransn(p)        = 0._r8
      veg_nf%plant_nalloc(p)          = 0._r8
      
   end subroutine NitrogenFluxVarsInit
   
   !-----------------------------------------------------------------------
   subroutine PhosphorusFluxVarsInit( p)
      !
      ! !DESCRIPTION:
      ! Initializes p-th patch of phosphorusflux_type
      !
      !$acc routine seq 
      implicit none
      !
      ! !ARGUMENT
      integer                  , intent(in)    :: p
      
      veg_pf%plant_pdemand(p)         = 0._r8
      veg_pf%avail_retransp(p)        = 0._r8
      veg_pf%plant_palloc(p)          = 0._r8
      
   end subroutine PhosphorusFluxVarsInit
   
   !-----------------------------------------------------------------------
   subroutine dyn_cnbal_column( bounds, clump_index, column_state_updater, &
      col_cs, c13_col_cs, c14_col_cs, &
      col_ns, col_ps)
      !
      ! !DESCRIPTION:
      ! Modify column-level state variables to maintain carbon, nitrogen
      ! and phosphorus balance with dynamic column weights.
      !
      ! !USES:
      !$acc routine seq 
      use dynColumnStateUpdaterMod, only : column_state_updater_type
      use dynPriorWeightsMod      , only : prior_weights_type
      use elm_varctl              , only : use_lch4
      !
      ! !ARGUMENTS:
      type(bounds_type)               , intent(in)    :: bounds
      integer                         , intent(in)    :: clump_index
      type(column_state_updater_type) , intent(in)    :: column_state_updater
      type(column_carbon_state)       , intent(inout) :: col_cs
      type(column_carbon_state)       , intent(inout) :: c13_col_cs
      type(column_carbon_state)       , intent(inout) :: c14_col_cs
      type(column_nitrogen_state)     , intent(inout) :: col_ns
      type(column_phosphorus_state)   , intent(inout) :: col_ps
      !
      ! !LOCAL VARIABLES:
      
      !    character(len=*), parameter :: subname = 'dyn_cnbal_col'
      !-----------------------------------------------------------------------
      !call dyn_col_cs_Adjustments(bounds, clump_index, column_state_updater, col_cs)
      ! if (use_c13) then
      !    call dyn_col_cs_Adjustments(bounds, clump_index, column_state_updater, c13_col_cs)
      ! end if
      ! if (use_c14) then
      !    call dyn_col_cs_Adjustments(bounds, clump_index, column_state_updater, c14_col_cs)
      ! end if
      
      ! call dyn_col_ns_Adjustments(bounds, clump_index, column_state_updater, col_ns)
      ! call dyn_col_ps_Adjustments(bounds, clump_index, column_state_updater, col_ps)
      
      ! DynamicColumnAdjustments for CH4 needs to be implemented
      
   end subroutine dyn_cnbal_column
   
end module dynConsBiogeochemMod
