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
  use CNCarbonFluxType         , only : carbonflux_type
  use CNCarbonStateType        , only : carbonstate_type
  use CNNitrogenFluxType       , only : nitrogenflux_type
  use CNNitrogenStateType      , only : nitrogenstate_type
  use PhosphorusFluxType       , only : phosphorusflux_type
  use PhosphorusStateType      , only : phosphorusstate_type
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
  use elm_varcon               , only : c3_r2, c4_r2, c14ratio
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
       carbonflux_vars, c13_carbonflux_vars, c14_carbonflux_vars, &
       nitrogenstate_vars, nitrogenflux_vars, &
       veg_ns, &
       phosphorusstate_vars,phosphorusflux_vars, &
       veg_ps)
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
    use clm_time_manager   , only : get_step_size
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
    type(carbonflux_type)          , intent(inout) :: carbonflux_vars
    type(carbonflux_type)          , intent(inout) :: c13_carbonflux_vars
    type(carbonflux_type)          , intent(inout) :: c14_carbonflux_vars
    type(nitrogenstate_type)       , intent(inout) :: nitrogenstate_vars
    type(nitrogenflux_type)        , intent(inout) :: nitrogenflux_vars
    type(vegetation_nitrogen_state), intent(inout) :: veg_ns
    type(phosphorusstate_type)     , intent(inout) :: phosphorusstate_vars
    type(phosphorusflux_type)      , intent(inout) :: phosphorusflux_vars
    type(vegetation_phosphorus_state),intent(inout) :: veg_ps
    !
    ! !LOCAL VARIABLES:
    integer                       :: pi,p,c,l,g,j                  ! indices
    integer                       :: ier                           ! error code
    real(r8)                      :: dwt                           ! change in pft weight (relative to column)
    real(r8)                      :: dt                            ! land model time step (sec)
    real(r8)                      :: init_h2ocan                   ! initial canopy water mass
    real(r8)                      :: new_h2ocan                    ! canopy water mass after weight shift
    real(r8), allocatable         :: dwt_leafc_seed(:)             ! pft-level mass gain due to seeding of new area
    real(r8), allocatable         :: dwt_leafn_seed(:)             ! pft-level mass gain due to seeding of new area
    real(r8), allocatable         :: dwt_deadstemc_seed(:)         ! pft-level mass gain due to seeding of new area
    real(r8), allocatable         :: dwt_deadstemn_seed(:)         ! pft-level mass gain due to seeding of new area
    real(r8), allocatable         :: dwt_npool_seed(:)             ! pft-level mass gain due to seeding of new area
    real(r8), allocatable         :: dwt_frootc_to_litter(:)       ! pft-level mass loss due to weight shift
    real(r8), allocatable         :: dwt_livecrootc_to_litter(:)   ! pft-level mass loss due to weight shift
    real(r8), allocatable         :: dwt_deadcrootc_to_litter(:)   ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_frootn_to_litter(:)       ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_livecrootn_to_litter(:)   ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_deadcrootn_to_litter(:)   ! pft-level mass loss due to weight shift
    real(r8), allocatable         :: conv_cflux(:)                 ! pft-level mass loss due to weight shift
    real(r8), allocatable         :: prod10_cflux(:)               ! pft-level mass loss due to weight shift
    real(r8), allocatable         :: prod100_cflux(:)              ! pft-level mass loss due to weight shift
    real(r8), allocatable         :: crop_product_cflux(:)         ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: conv_nflux(:)                 ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod10_nflux(:)               ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod100_nflux(:)              ! pft-level mass loss due to weight shift
    real(r8), allocatable         :: crop_product_nflux(:)         ! pft-level mass loss due to weight shift
    real(r8)                      :: t1,t2,wt_new,wt_old
    real(r8)                      :: init_state, change_state, new_state
    real(r8)                      :: tot_leaf, pleaf, pstor, pxfer
    real(r8)                      :: leafc_seed, leafn_seed
    real(r8)                      :: deadstemc_seed, deadstemn_seed, npool_seed
    real(r8), pointer             :: dwt_ptr0, dwt_ptr1, dwt_ptr2, dwt_ptr3, ptr
    character(len=32)             :: subname='dyn_cbal'            ! subroutine name

    ! ! add phosphorus local variables
    real(r8), allocatable         :: dwt_leafp_seed(:)             ! pft-level mass gain due to seeding of new area    
    real(r8), allocatable         :: dwt_deadstemp_seed(:)         ! pft-level mass gain due to seeding of new area
    real(r8), allocatable         :: dwt_ppool_seed(:)             ! pft-level mass gain due to seeding of new area
    real(r8), allocatable, target :: dwt_frootp_to_litter(:)       ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_livecrootp_to_litter(:)   ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_deadcrootp_to_litter(:)   ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: conv_pflux(:)                 ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod10_pflux(:)               ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod100_pflux(:)              ! pft-level mass loss due to weight shift
    real(r8), allocatable         :: crop_product_pflux(:)         ! pft-level mass loss due to weight shift
    real(r8)                      :: leafp_seed
    real(r8)                      :: deadstemp_seed, ppool_seed

    !! C13
    real(r8), allocatable         :: dwt_leafc13_seed(:)           ! pft-level mass gain due to seeding of new area
    real(r8), allocatable         :: dwt_deadstemc13_seed(:)       ! pft-level mass gain due to seeding of new area
    real(r8), allocatable, target :: dwt_frootc13_to_litter(:)     ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_livecrootc13_to_litter(:) ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_deadcrootc13_to_litter(:) ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: conv_c13flux(:)               ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod10_c13flux(:)             ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod100_c13flux(:)            ! pft-level mass loss due to weight shift
    real(r8), allocatable         :: crop_product_c13flux(:)       ! pft-level mass loss due to weight shift
    real(r8)                      :: leafc13_seed, deadstemc13_seed
    !! C14
    real(r8), allocatable         :: dwt_leafc14_seed(:)           ! pft-level mass gain due to seeding of new area
    real(r8), allocatable         :: dwt_deadstemc14_seed(:)       ! pft-level mass gain due to seeding of new area
    real(r8), allocatable, target :: dwt_frootc14_to_litter(:)     ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_livecrootc14_to_litter(:) ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_deadcrootc14_to_litter(:) ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: conv_c14flux(:)               ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod10_c14flux(:)             ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod100_c14flux(:)            ! pft-level mass loss due to weight shift
    real(r8), allocatable         :: crop_product_c14flux(:)       ! pft-level mass loss due to weight shift
    real(r8)                      :: leafc14_seed, deadstemc14_seed
    real(r8) :: froot, croot
    real(r8) :: fr_flab, fr_fcel, fr_flig
    !-----------------------------------------------------------------------
    
    associate(&
         cs     => veg_cs    , &
         c13_cs => c13_veg_cs, &
         c14_cs => c14_veg_cs, &
         cf     => carbonflux_vars     , &
         c13_cf => c13_carbonflux_vars , &
         c14_cf => c14_carbonflux_vars , &
         ns     => veg_ns  , &
         nf     => nitrogenflux_vars   ,  &
         ps     => veg_ps  , &
         pf     => phosphorusflux_vars     &
         )


    ! Allocate pft-level mass loss arrays
    allocate(dwt_leafc_seed           (bounds%begp:bounds%endp), stat=ier)
    allocate(dwt_leafn_seed           (bounds%begp:bounds%endp), stat=ier)
    allocate(dwt_deadstemc_seed       (bounds%begp:bounds%endp), stat=ier)
    allocate(dwt_deadstemn_seed       (bounds%begp:bounds%endp), stat=ier)
    allocate(dwt_npool_seed           (bounds%begp:bounds%endp), stat=ier)
    allocate(dwt_frootc_to_litter     (bounds%begp:bounds%endp), stat=ier)
    allocate(dwt_livecrootc_to_litter (bounds%begp:bounds%endp), stat=ier)
    allocate(dwt_deadcrootc_to_litter (bounds%begp:bounds%endp), stat=ier)
    allocate(dwt_frootn_to_litter     (bounds%begp:bounds%endp), stat=ier)
    allocate(dwt_livecrootn_to_litter (bounds%begp:bounds%endp), stat=ier)
    allocate(dwt_deadcrootn_to_litter (bounds%begp:bounds%endp), stat=ier)
    allocate(conv_cflux               (bounds%begp:bounds%endp), stat=ier)
    allocate(prod10_cflux             (bounds%begp:bounds%endp), stat=ier)
    allocate(prod100_cflux            (bounds%begp:bounds%endp), stat=ier)
    allocate(crop_product_cflux       (bounds%begp:bounds%endp), stat=ier)
    allocate(conv_nflux               (bounds%begp:bounds%endp), stat=ier)
    allocate(prod10_nflux             (bounds%begp:bounds%endp), stat=ier)
    allocate(prod100_nflux            (bounds%begp:bounds%endp), stat=ier)
    allocate(crop_product_nflux       (bounds%begp:bounds%endp), stat=ier)

    ! Allocate P arrays
    allocate(dwt_leafp_seed           (bounds%begp:bounds%endp), stat=ier)
    allocate(dwt_deadstemp_seed       (bounds%begp:bounds%endp), stat=ier)
    allocate(dwt_ppool_seed           (bounds%begp:bounds%endp), stat=ier)
    allocate(dwt_frootp_to_litter     (bounds%begp:bounds%endp), stat=ier)
    allocate(dwt_livecrootp_to_litter (bounds%begp:bounds%endp), stat=ier)
    allocate(dwt_deadcrootp_to_litter (bounds%begp:bounds%endp), stat=ier)
    allocate(conv_pflux               (bounds%begp:bounds%endp), stat=ier)
    allocate(prod10_pflux             (bounds%begp:bounds%endp), stat=ier)
    allocate(prod100_pflux            (bounds%begp:bounds%endp), stat=ier)
    allocate(crop_product_pflux       (bounds%begp:bounds%endp), stat=ier)

    if ( use_c13 ) then
       allocate(dwt_leafc13_seed           (bounds%begp:bounds%endp), stat=ier)
       allocate(dwt_deadstemc13_seed       (bounds%begp:bounds%endp), stat=ier)
       allocate(dwt_frootc13_to_litter     (bounds%begp:bounds%endp), stat=ier)
       allocate(dwt_livecrootc13_to_litter (bounds%begp:bounds%endp), stat=ier)
       allocate(dwt_deadcrootc13_to_litter (bounds%begp:bounds%endp), stat=ier)
       allocate(conv_c13flux               (bounds%begp:bounds%endp), stat=ier)
       allocate(prod10_c13flux             (bounds%begp:bounds%endp), stat=ier)
       allocate(prod100_c13flux            (bounds%begp:bounds%endp), stat=ier)
       allocate(crop_product_c13flux       (bounds%begp:bounds%endp), stat=ier)
    endif
    if ( use_c14 ) then
       allocate(dwt_leafc14_seed           (bounds%begp:bounds%endp), stat=ier)
       allocate(dwt_deadstemc14_seed       (bounds%begp:bounds%endp), stat=ier)
       allocate(dwt_frootc14_to_litter     (bounds%begp:bounds%endp), stat=ier)
       allocate(dwt_livecrootc14_to_litter (bounds%begp:bounds%endp), stat=ier)
       allocate(dwt_deadcrootc14_to_litter (bounds%begp:bounds%endp), stat=ier)
       allocate(conv_c14flux               (bounds%begp:bounds%endp), stat=ier)
       allocate(prod10_c14flux             (bounds%begp:bounds%endp), stat=ier)
       allocate(prod100_c14flux            (bounds%begp:bounds%endp), stat=ier)
       allocate(crop_product_c14flux       (bounds%begp:bounds%endp), stat=ier)
    endif
    
    ! Get time step
    dt = real( get_step_size(), r8 )
    
    do p = bounds%begp,bounds%endp
       c = veg_pp%column(p)
       ! initialize all the pft-level local flux arrays
       dwt_leafc_seed(p)           = 0._r8
       dwt_deadstemc_seed(p)       = 0._r8
       dwt_frootc_to_litter(p)     = 0._r8
       dwt_livecrootc_to_litter(p) = 0._r8
       dwt_deadcrootc_to_litter(p) = 0._r8
       conv_cflux(p)               = 0._r8
       prod10_cflux(p)             = 0._r8
       prod100_cflux(p)            = 0._r8
       crop_product_cflux(p)       = 0._r8
       
       dwt_leafn_seed(p)           = 0._r8
       dwt_deadstemn_seed(p)       = 0._r8
       dwt_npool_seed(p)           = 0._r8
       dwt_frootn_to_litter(p)     = 0._r8
       dwt_livecrootn_to_litter(p) = 0._r8
       dwt_deadcrootn_to_litter(p) = 0._r8
       conv_nflux(p)               = 0._r8
       prod10_nflux(p)             = 0._r8
       prod100_nflux(p)            = 0._r8
       crop_product_nflux(p)       = 0._r8
       
       dwt_leafp_seed(p)           = 0._r8
       dwt_deadstemp_seed(p)       = 0._r8
       dwt_ppool_seed(p)           = 0._r8
       dwt_frootp_to_litter(p)     = 0._r8
       dwt_livecrootp_to_litter(p) = 0._r8
       dwt_deadcrootp_to_litter(p) = 0._r8
       conv_pflux(p)               = 0._r8
       prod10_pflux(p)             = 0._r8
       prod100_pflux(p)            = 0._r8
       crop_product_pflux(p)       = 0._r8

       if ( use_c13 ) then
          dwt_leafc13_seed(p)           = 0._r8
          dwt_deadstemc13_seed(p)       = 0._r8
          dwt_frootc13_to_litter(p)     = 0._r8
          dwt_livecrootc13_to_litter(p) = 0._r8
          dwt_deadcrootc13_to_litter(p) = 0._r8
          conv_c13flux(p)               = 0._r8
          prod10_c13flux(p)             = 0._r8
          prod100_c13flux(p)            = 0._r8
          crop_product_c13flux(p)       = 0._r8
       endif
       
       if ( use_c14 ) then
          dwt_leafc14_seed(p)           = 0._r8
          dwt_deadstemc14_seed(p)       = 0._r8
          dwt_frootc14_to_litter(p)     = 0._r8
          dwt_livecrootc14_to_litter(p) = 0._r8
          dwt_deadcrootc14_to_litter(p) = 0._r8
          conv_c14flux(p)               = 0._r8
          prod10_c14flux(p)             = 0._r8
          prod100_c14flux(p)            = 0._r8
          crop_product_c14flux(p)       = 0._r8
       endif
       
       l = veg_pp%landunit(p)
       if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then
          
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
                call CarbonStateVarsInit     (cs, p)
                call NitrogenStateVarsInit   (veg_ns, p)
                call PhosphorusStateVarsInit (veg_ps, p)
                call CanopyStateVarsInit     (canopystate_vars, p)
                call CNStateVarsInit         (cnstate_vars, p, c)
                call CarbonFluxVarsInit      (veg_cf, p)
                call NitrogenFluxVarsInit    (nf, p)
                call PhosphorusFluxVarsInit  (pf, p)
                
                if ( use_c13 ) then
                   call CarbonStateVarsInit(c13_cs, p)
                endif
                
                if ( use_c14 ) then
                   call CarbonStateVarsInit(c14_cs, p)
                endif
                                 
                ! add phosphorus related variables

                call photosyns_vars%NewPatchinit(p)
                
             end if  ! end initialization of new pft
                                                    
          end if       ! weight decreasing
       end if           ! is soil
    end do               ! patch loop
    
    call dyn_veg_cs_Adjustments(    &
         bounds,                        &
         num_soilp_with_inactive,       &
         filter_soilp_with_inactive,    &
         prior_weights,                 &
         patch_state_updater,           &
         dwt_leafc_seed,                &
         dwt_deadstemc_seed,            &
         conv_cflux,                    &
         dwt_frootc_to_litter,          &
         dwt_livecrootc_to_litter,      &
         dwt_deadcrootc_to_litter,      &
         prod10_cflux,                  &
         prod100_cflux,                 &
         crop_product_cflux,            &
         veg_cs                         &
         )

    if (use_c13) then
       call dyn_veg_cs_Adjustments( &
            bounds,                        &
            num_soilp_with_inactive,       &
            filter_soilp_with_inactive,    &
            prior_weights,                 &
            patch_state_updater,           &
            dwt_leafc13_seed,              &
            dwt_deadstemc13_seed,          &
            conv_c13flux,                  &
            dwt_frootc13_to_litter,        &
            dwt_livecrootc13_to_litter,    &
            dwt_deadcrootc13_to_litter,    &
            prod10_c13flux,                &
            prod100_c13flux,               &
            crop_product_c13flux,          &
            c13_veg_cs                     &
            )
    endif

    if (use_c14) then
       call dyn_veg_cs_Adjustments( &
            bounds,                        &
            num_soilp_with_inactive,       &
            filter_soilp_with_inactive,    &
            prior_weights,                 &
            patch_state_updater,           &
            dwt_leafc14_seed,              &
            dwt_deadstemc14_seed,          &
            conv_c14flux,                  &
            dwt_frootc14_to_litter,        &
            dwt_livecrootc14_to_litter,    &
            dwt_deadcrootc14_to_litter,    &
            prod10_c14flux,                &
            prod100_c14flux,               &
            crop_product_c14flux,          &
            c14_veg_cs                     &
            )
    endif

    call dyn_veg_ns_Adjustments(    &
         bounds,                        &
         num_soilp_with_inactive,       &
         filter_soilp_with_inactive,    &
         prior_weights,                 &
         patch_state_updater,           &
         dwt_leafn_seed,                &
         dwt_deadstemn_seed,            &
         dwt_npool_seed,                &
         conv_nflux,                    &
         dwt_frootn_to_litter,          &
         dwt_livecrootn_to_litter,      &
         dwt_deadcrootn_to_litter,      &
         prod10_nflux,                  &
         prod100_nflux,                 &
         crop_product_nflux,            &
         veg_ns                         &
         )

    call dyn_veg_ps_Adjustments(    &
         bounds,                        &
         num_soilp_with_inactive,       &
         filter_soilp_with_inactive,    &
         prior_weights,                 &
         patch_state_updater,           &
         dwt_leafp_seed,                &
         dwt_deadstemp_seed,            &
         dwt_ppool_seed,                &
         conv_pflux,                    &
         dwt_frootp_to_litter,          &
         dwt_livecrootp_to_litter,      &
         dwt_deadcrootp_to_litter,      &
         prod10_pflux,                  &
         prod100_pflux,                 &
         crop_product_pflux,            &
         veg_ps                         &
         )

    ! calculate column-level seeding fluxes
    do p = bounds%begp, bounds%endp
       g = veg_pp%gridcell(p)

       ! C fluxes
       veg_cf%dwt_seedc_to_leaf(p) = dwt_leafc_seed(p)/dt
       grc_cf%dwt_seedc_to_leaf(g)   = &
            grc_cf%dwt_seedc_to_leaf(g) + &
            veg_cf%dwt_seedc_to_leaf(p)

       veg_cf%dwt_seedc_to_deadstem(p) = dwt_deadstemc_seed(p)/dt
       grc_cf%dwt_seedc_to_deadstem(g)   = &
            grc_cf%dwt_seedc_to_deadstem(g) + &
            veg_cf%dwt_seedc_to_deadstem(p)

       if ( use_c13 ) then
          c13_veg_cf%dwt_seedc_to_leaf(p) = dwt_leafc_seed(p)/dt
          c13_grc_cf%dwt_seedc_to_leaf(g)   = &
               c13_grc_cf%dwt_seedc_to_leaf(g) + &
               c13_veg_cf%dwt_seedc_to_leaf(p)

          c13_veg_cf%dwt_seedc_to_deadstem(p) = dwt_deadstemc_seed(p)/dt
          c13_grc_cf%dwt_seedc_to_deadstem(g)   = &
               c13_grc_cf%dwt_seedc_to_deadstem(g) + &
               c13_veg_cf%dwt_seedc_to_deadstem(p)
       endif

       if ( use_c14 ) then
          c14_veg_cf%dwt_seedc_to_leaf(p) = dwt_leafc_seed(p)/dt
          c14_grc_cf%dwt_seedc_to_leaf(g)   = &
               c14_grc_cf%dwt_seedc_to_leaf(g) + &
               c14_veg_cf%dwt_seedc_to_leaf(p)

          c14_veg_cf%dwt_seedc_to_deadstem(p) = dwt_deadstemc_seed(p)/dt
          c14_grc_cf%dwt_seedc_to_deadstem(g)   = &
               c14_grc_cf%dwt_seedc_to_deadstem(g) + &
               c14_veg_cf%dwt_seedc_to_deadstem(p)
       endif

       ! N fluxes
       veg_nf%dwt_seedn_to_leaf(p)   = dwt_leafn_seed(p)/dt
       grc_nf%dwt_seedn_to_leaf(g)     = &
            grc_nf%dwt_seedn_to_leaf(g) + &
            veg_nf%dwt_seedn_to_leaf(p)

       veg_nf%dwt_seedn_to_deadstem(p) = dwt_deadstemn_seed(p)/dt
       grc_nf%dwt_seedn_to_deadstem(g)   = &
            grc_nf%dwt_seedn_to_deadstem(g) + &
            veg_nf%dwt_seedn_to_deadstem(p)


       veg_nf%dwt_seedn_to_npool(p) = dwt_npool_seed(p)/dt
       grc_nf%dwt_seedn_to_npool(g)   = &
            grc_nf%dwt_seedn_to_npool(g) + &
            veg_nf%dwt_seedn_to_npool(p)

       ! P fluxes
       veg_pf%dwt_seedp_to_leaf(p)   = dwt_leafp_seed(p)/dt
       grc_pf%dwt_seedp_to_leaf(g)     = &
            grc_pf%dwt_seedp_to_leaf(g) + &
            veg_pf%dwt_seedp_to_leaf(p)

       veg_pf%dwt_seedp_to_deadstem(p) = dwt_deadstemp_seed(p)/dt
       grc_pf%dwt_seedp_to_deadstem(g)   = &
            grc_pf%dwt_seedp_to_deadstem(g) + &
            veg_pf%dwt_seedp_to_deadstem(p)


       veg_pf%dwt_seedp_to_ppool(p) = dwt_npool_seed(p)/dt
       grc_pf%dwt_seedp_to_ppool(g)   = &
            grc_pf%dwt_seedp_to_ppool(g) + &
            veg_pf%dwt_seedp_to_ppool(p)

    end do
    
    ! calculate patch-to-column slash fluxes into litter and CWD pools
    do p = bounds%begp, bounds%endp
       c = veg_pp%column(p)

       ! fine and coarse root to litter and CWD slash carbon fluxes
       col_cf%dwt_slash_cflux(c) =            &
            col_cf%dwt_slash_cflux(c)       + &
            dwt_frootc_to_litter(p)     /dt + &
            dwt_livecrootc_to_litter(p) /dt + &
            dwt_deadcrootc_to_litter(p) /dt

       if ( use_c13 ) then
          c13_col_cf%dwt_slash_cflux(c) =          &
               c13_col_cf%dwt_slash_cflux(c)     + &
               dwt_frootc13_to_litter(p)     /dt + &
               dwt_livecrootc13_to_litter(p) /dt + &
               dwt_deadcrootc13_to_litter(p) /dt
       endif

       if ( use_c14 ) then
          c14_col_cf%dwt_slash_cflux(c) =          &
               c14_col_cf%dwt_slash_cflux(c)     + &
               dwt_frootc14_to_litter(p)     /dt + &
               dwt_livecrootc14_to_litter(p) /dt + &
               dwt_deadcrootc14_to_litter(p) /dt
       endif

       col_nf%dwt_slash_nflux(c) =            &
            col_nf%dwt_slash_nflux(c)       + &
            dwt_frootn_to_litter(p)     /dt + &
            dwt_livecrootn_to_litter(p) /dt + &
            dwt_deadcrootn_to_litter(p) /dt

       col_pf%dwt_slash_pflux(c) =            &
            col_pf%dwt_slash_pflux(c)       + &
            dwt_frootp_to_litter(p)     /dt + &
            dwt_livecrootp_to_litter(p) /dt + &
            dwt_deadcrootp_to_litter(p) /dt

    end do

    ! calculate pft-to-column for fluxes into litter and CWD pools
    do j = 1, nlevdecomp
       do p = bounds%begp, bounds%endp
          c = veg_pp%column(p)

          froot   = cnstate_vars%froot_prof_patch(p,j)
          croot   = cnstate_vars%croot_prof_patch(p,j)
          fr_flab = veg_vp%fr_flab(veg_pp%itype(p))
          fr_fcel = veg_vp%fr_fcel(veg_pp%itype(p))
          fr_flig = veg_vp%fr_flig(veg_pp%itype(p))

          ! fine root litter carbon fluxes
          col_cf%dwt_frootc_to_litr_met_c(c,j) = &
               col_cf%dwt_frootc_to_litr_met_c(c,j) + &
               (dwt_frootc_to_litter(p)* fr_flab)/dt * froot

          col_cf%dwt_frootc_to_litr_cel_c(c,j) = &
               col_cf%dwt_frootc_to_litr_cel_c(c,j) + &
               (dwt_frootc_to_litter(p)* fr_fcel)/dt * froot

          col_cf%dwt_frootc_to_litr_lig_c(c,j) = &
               col_cf%dwt_frootc_to_litr_lig_c(c,j) + &
               (dwt_frootc_to_litter(p)* fr_flig)/dt * froot


          ! fine root litter nitrogen fluxes
          col_nf%dwt_frootn_to_litr_met_n(c,j) = &
               col_nf%dwt_frootn_to_litr_met_n(c,j) + &
               (dwt_frootn_to_litter(p)* fr_flab)/dt * froot
          col_nf%dwt_frootn_to_litr_cel_n(c,j) = &

               col_nf%dwt_frootn_to_litr_cel_n(c,j) + &
               (dwt_frootn_to_litter(p)* fr_fcel)/dt * froot

          col_nf%dwt_frootn_to_litr_lig_n(c,j) = &
               col_nf%dwt_frootn_to_litr_lig_n(c,j) + &
               (dwt_frootn_to_litter(p)* fr_flig)/dt * froot
               
          ! fine root litter phosphorus fluxes
          col_pf%dwt_frootp_to_litr_met_p(c,j) = &
               col_pf%dwt_frootp_to_litr_met_p(c,j) + &
               (dwt_frootp_to_litter(p)* fr_flab)/dt * froot
          col_pf%dwt_frootp_to_litr_cel_p(c,j) = &

               col_pf%dwt_frootp_to_litr_cel_p(c,j) + &
               (dwt_frootp_to_litter(p)* fr_fcel)/dt * froot

          col_pf%dwt_frootp_to_litr_lig_p(c,j) = &
               col_pf%dwt_frootp_to_litr_lig_p(c,j) + &
               (dwt_frootp_to_litter(p)* fr_flig)/dt * froot

          ! livecroot fluxes to cwd
          col_cf%dwt_livecrootc_to_cwdc(c,j) = &
               col_cf%dwt_livecrootc_to_cwdc(c,j) + &
               (dwt_livecrootc_to_litter(p))/dt * croot

          col_nf%dwt_livecrootn_to_cwdn(c,j) = &
               col_nf%dwt_livecrootn_to_cwdn(c,j) + &
               (dwt_livecrootn_to_litter(p))/dt * croot

          col_pf%dwt_livecrootp_to_cwdp(c,j) = &
               col_pf%dwt_livecrootp_to_cwdp(c,j) + &
               (dwt_livecrootp_to_litter(p))/dt * croot

          ! deadcroot fluxes to cwd
          col_cf%dwt_deadcrootc_to_cwdc(c,j) = &
               col_cf%dwt_deadcrootc_to_cwdc(c,j) + &
               (dwt_deadcrootc_to_litter(p))/dt * croot

          col_nf%dwt_deadcrootn_to_cwdn(c,j) = &
               col_nf%dwt_deadcrootn_to_cwdn(c,j) + &
               (dwt_deadcrootn_to_litter(p))/dt * croot

          col_pf%dwt_deadcrootp_to_cwdp(c,j) = &
               col_pf%dwt_deadcrootp_to_cwdp(c,j) + &
               (dwt_deadcrootp_to_litter(p))/dt * croot
               
          if ( use_c13 ) then
             ! C13 fine root litter fluxes
             c13_col_cf%dwt_frootc_to_litr_met_c(c,j) = &
                  c13_col_cf%dwt_frootc_to_litr_met_c(c,j) + &
                  (dwt_frootc13_to_litter(p)* fr_flab)/dt * froot

             c13_col_cf%dwt_frootc_to_litr_cel_c(c,j) = &
                  c13_col_cf%dwt_frootc_to_litr_cel_c(c,j) + &
                  (dwt_frootc13_to_litter(p)* fr_fcel)/dt * froot

             c13_col_cf%dwt_frootc_to_litr_lig_c(c,j) = &
                  c13_col_cf%dwt_frootc_to_litr_lig_c(c,j) + &
                  (dwt_frootc13_to_litter(p)* fr_flig)/dt * froot

             ! livecroot fluxes to cwd
             c13_col_cf%dwt_livecrootc_to_cwdc(c,j) = &
                  c13_col_cf%dwt_livecrootc_to_cwdc(c,j) + &
                  (dwt_livecrootc13_to_litter(p))/dt * croot

             ! deadcroot fluxes to cwd
             c13_col_cf%dwt_deadcrootc_to_cwdc(c,j) = &
                  c13_col_cf%dwt_deadcrootc_to_cwdc(c,j) + &
                  (dwt_deadcrootc13_to_litter(p))/dt * croot

          endif
          
          if ( use_c14 ) then
             ! C14 fine root litter fluxes
             c14_col_cf%dwt_frootc_to_litr_met_c(c,j) = &
                  c14_col_cf%dwt_frootc_to_litr_met_c(c,j) + &
                  (dwt_frootc14_to_litter(p)* fr_flab)/dt * froot

             c14_col_cf%dwt_frootc_to_litr_cel_c(c,j) = &
                  c14_col_cf%dwt_frootc_to_litr_cel_c(c,j) + &
                  (dwt_frootc14_to_litter(p)* fr_fcel)/dt * froot

             c14_col_cf%dwt_frootc_to_litr_lig_c(c,j) = &
                  c14_col_cf%dwt_frootc_to_litr_lig_c(c,j) + &
                  (dwt_frootc14_to_litter(p)* fr_flig)/dt * froot

             ! livecroot fluxes to cwd
             c14_col_cf%dwt_livecrootc_to_cwdc(c,j) = &
                  c14_col_cf%dwt_livecrootc_to_cwdc(c,j) + &
                  (dwt_livecrootc14_to_litter(p))/dt * croot

             ! deadcroot fluxes to cwd
             c14_col_cf%dwt_deadcrootc_to_cwdc(c,j) = &
                  c14_col_cf%dwt_deadcrootc_to_cwdc(c,j) + &
                  (dwt_deadcrootc14_to_litter(p))/dt * croot
          endif
       end do
    end do
    ! calculate pft-to-column for fluxes into product pools and conversion flux
    do p = bounds%begp, bounds%endp 
       c = veg_pp%column(p)
       g = veg_pp%gridcell(p) 

       ! column-level fluxes are accumulated as positive fluxes.
       ! column-level C flux updates
       col_cf%dwt_conv_cflux(c) = col_cf%dwt_conv_cflux(c) - conv_cflux(p)/dt
       col_cf%dwt_prod10c_gain(c) = col_cf%dwt_prod10c_gain(c) - prod10_cflux(p)/dt
       col_cf%dwt_prod100c_gain(c) = col_cf%dwt_prod100c_gain(c) - prod100_cflux(p)/dt
       col_cf%dwt_crop_productc_gain(c) = col_cf%dwt_crop_productc_gain(c) - crop_product_cflux(p)/dt

       veg_cf%dwt_prod10c_gain(p) = - prod10_cflux(p)/dt
       grc_cf%dwt_prod10c_gain(g)   = grc_cf%dwt_prod10c_gain(g) + veg_cf%dwt_prod10c_gain(p)

       veg_cf%dwt_prod100c_gain(p) = - prod100_cflux(p)/dt
       grc_cf%dwt_prod100c_gain(g)   = grc_cf%dwt_prod100c_gain(g) + veg_cf%dwt_prod100c_gain(p)

       veg_cf%dwt_crop_productc_gain(p) = - crop_product_cflux(p)/dt

       if ( use_c13 ) then
          ! C13 column-level flux updates
          c13_col_cf%dwt_conv_cflux(c) = c13_col_cf%dwt_conv_cflux(c) - conv_c13flux(p)/dt
          c13_col_cf%dwt_prod10c_gain(c) = c13_col_cf%dwt_prod10c_gain(c) - prod10_c13flux(p)/dt
          c13_col_cf%dwt_prod100c_gain(c) = c13_col_cf%dwt_prod100c_gain(c) - prod100_c13flux(p)/dt
          c13_col_cf%dwt_crop_productc_gain(c) = c13_col_cf%dwt_crop_productc_gain(c) - crop_product_c13flux(p)/dt

          c13_veg_cf%dwt_prod10c_gain(p) = - prod10_c13flux(p)/dt
          c13_grc_cf%dwt_prod10c_gain(g)   = c13_grc_cf%dwt_prod10c_gain(g) + c13_veg_cf%dwt_prod10c_gain(p)

          c13_veg_cf%dwt_prod100c_gain(p) = - prod100_c13flux(p)/dt
          c13_grc_cf%dwt_prod100c_gain(g)   = c13_grc_cf%dwt_prod100c_gain(g) + c13_veg_cf%dwt_prod100c_gain(p)

          c13_veg_cf%dwt_crop_productc_gain(p) = - crop_product_c13flux(p)/dt

       endif
       
       if ( use_c14 ) then
          ! C14 column-level flux updates
          c14_col_cf%dwt_conv_cflux(c) = c14_col_cf%dwt_conv_cflux(c) - conv_c14flux(p)/dt
          c14_col_cf%dwt_prod10c_gain(c) = c14_col_cf%dwt_prod10c_gain(c) - prod10_c14flux(p)/dt
          c14_col_cf%dwt_prod100c_gain(c) = c14_col_cf%dwt_prod100c_gain(c) - prod100_c14flux(p)/dt
          c14_col_cf%dwt_crop_productc_gain(c) = c14_col_cf%dwt_crop_productc_gain(c) - crop_product_c14flux(p)/dt

          c14_veg_cf%dwt_prod10c_gain(p) = - prod10_c14flux(p)/dt
          c14_grc_cf%dwt_prod10c_gain(g)   = c14_grc_cf%dwt_prod10c_gain(g) + c14_veg_cf%dwt_prod10c_gain(p)

          c14_veg_cf%dwt_prod100c_gain(p) = - prod100_c14flux(p)/dt
          c14_grc_cf%dwt_prod100c_gain(g)   = c14_grc_cf%dwt_prod100c_gain(g) + c14_veg_cf%dwt_prod100c_gain(p)

          c14_veg_cf%dwt_crop_productc_gain(p) = - crop_product_c14flux(p)/dt

       endif

       ! column-level N flux updates
       col_nf%dwt_conv_nflux(c) = col_nf%dwt_conv_nflux(c) - conv_nflux(p)/dt
       col_nf%dwt_prod10n_gain(c) = col_nf%dwt_prod10n_gain(c) - prod10_nflux(p)/dt
       col_nf%dwt_prod100n_gain(c) = col_nf%dwt_prod100n_gain(c) - prod100_nflux(p)/dt
       col_nf%dwt_crop_productn_gain(c) = col_nf%dwt_crop_productn_gain(c) - crop_product_nflux(p)/dt

       veg_nf%dwt_prod10n_gain(p) = -prod10_nflux(p)/dt
       grc_nf%dwt_prod10n_gain(g)   = grc_nf%dwt_prod10n_gain(g) + veg_nf%dwt_prod10n_gain(p)

       veg_nf%dwt_prod100n_gain(p)= -prod100_nflux(p)/dt
       grc_nf%dwt_prod100n_gain(g)  = grc_nf%dwt_prod100n_gain(g) + veg_nf%dwt_prod100n_gain(p)

       veg_nf%dwt_crop_productn_gain(p) = -crop_product_nflux(p)/dt
       
       ! column-level P flux updates

       col_pf%dwt_conv_pflux(c) = col_pf%dwt_conv_pflux(c) - conv_pflux(p)/dt
       col_pf%dwt_prod10p_gain(c) = col_pf%dwt_prod10p_gain(c) - prod10_pflux(p)/dt
       col_pf%dwt_prod100p_gain(c) = col_pf%dwt_prod100p_gain(c) - prod100_pflux(p)/dt
       col_pf%dwt_crop_productp_gain(c) = col_pf%dwt_crop_productp_gain(c) - crop_product_pflux(p)/dt

       veg_pf%dwt_prod10p_gain(p) = -prod10_pflux(p)/dt
       grc_pf%dwt_prod10p_gain(g)   = grc_pf%dwt_prod10p_gain(g) + veg_pf%dwt_prod10p_gain(p)

       veg_pf%dwt_prod100p_gain(p)= -prod100_pflux(p)/dt
       grc_pf%dwt_prod100p_gain(g)  = grc_pf%dwt_prod100p_gain(g) + veg_pf%dwt_prod100p_gain(p)

       veg_pf%dwt_crop_productp_gain(p) = -crop_product_pflux(p)/dt
             
    end do
    
    do p = bounds%begp, bounds%endp
       g = veg_pp%gridcell(p)

       ! Note that patch-level fluxes are stored per unit GRIDCELL area - thus, we don't
       ! need to multiply by the patch's gridcell weight when translating patch-level
       ! fluxes into gridcell-level fluxes.

       veg_cf%dwt_conv_cflux(p) = -conv_cflux(p)/dt
       grc_cf%dwt_conv_cflux(g) = &
            grc_cf%dwt_conv_cflux(g) + &
            veg_cf%dwt_conv_cflux(p)

       if ( use_c13 ) then
          ! C13 column-level flux updates
          c13_veg_cf%dwt_conv_cflux(p) = -conv_c13flux(p)/dt
          c13_grc_cf%dwt_conv_cflux(g) = &
               c13_grc_cf%dwt_conv_cflux(g) + &
               c13_veg_cf%dwt_conv_cflux(p)
       endif

       if ( use_c14 ) then
          ! C14 column-level flux updates
          c14_veg_cf%dwt_conv_cflux(p) = -conv_c14flux(p)/dt
          c14_grc_cf%dwt_conv_cflux(g) = &
               c14_grc_cf%dwt_conv_cflux(g) + &
               c14_veg_cf%dwt_conv_cflux(p)
       endif

       veg_nf%dwt_conv_nflux(p) = -conv_nflux(p)/dt
       grc_nf%dwt_conv_nflux(g) = &
            grc_nf%dwt_conv_nflux(g) + &
            veg_nf%dwt_conv_nflux(p)

       veg_pf%dwt_conv_pflux(p) = -conv_pflux(p)/dt
       grc_pf%dwt_conv_pflux(g) = &
            grc_pf%dwt_conv_pflux(g) + &
            veg_pf%dwt_conv_pflux(p)

    end do

    ! Deallocate pft-level flux arrays
    deallocate(dwt_leafc_seed)
    deallocate(dwt_leafn_seed)
    deallocate(dwt_deadstemc_seed)
    deallocate(dwt_deadstemn_seed)
    deallocate(dwt_npool_seed)
    deallocate(dwt_frootc_to_litter)
    deallocate(dwt_livecrootc_to_litter)
    deallocate(dwt_deadcrootc_to_litter)
    deallocate(dwt_frootn_to_litter)
    deallocate(dwt_livecrootn_to_litter)
    deallocate(dwt_deadcrootn_to_litter)
    deallocate(conv_cflux)
    deallocate(prod10_cflux)
    deallocate(prod100_cflux)
    deallocate(crop_product_cflux)
    deallocate(conv_nflux)
    deallocate(prod10_nflux)
    deallocate(prod100_nflux)
    deallocate(crop_product_nflux)

    deallocate(dwt_leafp_seed)
    deallocate(dwt_deadstemp_seed)
    deallocate(dwt_ppool_seed)
    deallocate(dwt_frootp_to_litter)
    deallocate(dwt_livecrootp_to_litter)
    deallocate(dwt_deadcrootp_to_litter)
    deallocate(conv_pflux)
    deallocate(prod10_pflux)
    deallocate(prod100_pflux)
    deallocate(crop_product_pflux)

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
    
   end associate
 end subroutine dyn_cnbal_patch
 
 !-----------------------------------------------------------------------
 subroutine CarbonStateVarsInit(cs, p)
   !
   ! !DESCRIPTION:
   ! Initializes p-th patch of carbonstate_type
   !
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
 subroutine NitrogenFluxVarsInit(nf, p)
   !
   ! !DESCRIPTION:
   ! Initializes p-th patch of nitrogenflux_type
   !
   implicit none
   !
   ! !ARGUMENT
   type(nitrogenflux_type), intent(inout) :: nf
   integer                , intent(in)    :: p

   veg_nf%plant_ndemand(p)         = 0._r8
   veg_nf%avail_retransn(p)        = 0._r8
   veg_nf%plant_nalloc(p)          = 0._r8

 end subroutine NitrogenFluxVarsInit

 !-----------------------------------------------------------------------
 subroutine PhosphorusFluxVarsInit(pf, p)
   !
   ! !DESCRIPTION:
   ! Initializes p-th patch of phosphorusflux_type
   !
   implicit none
   !
   ! !ARGUMENT
   type(phosphorusflux_type), intent(inout) :: pf
   integer                  , intent(in)    :: p

   veg_pf%plant_pdemand(p)         = 0._r8
   veg_pf%avail_retransp(p)        = 0._r8
   veg_pf%plant_palloc(p)          = 0._r8

 end subroutine PhosphorusFluxVarsInit

 !-----------------------------------------------------------------------
 subroutine dyn_cnbal_column( bounds, clump_index, column_state_updater, &
       col_cs, c13_col_cs, c14_col_cs, &
       phosphorusstate_vars, col_ns, col_ps)
   !
   ! !DESCRIPTION:
   ! Modify column-level state variables to maintain carbon, nitrogen
   ! and phosphorus balance with dynamic column weights.
   !
   ! !USES:
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
   type(phosphorusstate_type)      , intent(inout) :: phosphorusstate_vars
   type(column_nitrogen_state)     , intent(inout) :: col_ns
   type(column_phosphorus_state)   , intent(inout) :: col_ps
   !
   ! !LOCAL VARIABLES:

   character(len=*), parameter :: subname = 'dyn_cnbal_col'
   !-----------------------------------------------------------------------

   associate(&
         cs     => col_cs,     &
         c13_cs => c13_col_cs, &
         c14_cs => c14_col_cs, &
         ns     => col_ns  , &
         ps     => phosphorusstate_vars  &
         )

    call dyn_col_cs_Adjustments(bounds, clump_index, column_state_updater, col_cs)

   if (use_c13) then
      call dyn_col_cs_Adjustments(bounds, clump_index, column_state_updater, c13_col_cs)
   end if

   if (use_c14) then
      call dyn_col_cs_Adjustments(bounds, clump_index, column_state_updater, c14_col_cs)
   end if

   call dyn_col_ns_Adjustments(bounds, clump_index, column_state_updater, col_ns)

   call dyn_col_ps_Adjustments(bounds, clump_index, column_state_updater, col_ps)

   ! DynamicColumnAdjustments for CH4 needs to be implemented

 end associate

end subroutine dyn_cnbal_column

 end module dynConsBiogeochemMod
