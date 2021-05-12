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
    type(vegetation_nitrogen_state), intent(inout) :: veg_ns
    type(vegetation_phosphorus_state),intent(inout) :: veg_ps
    real(r8)                         ,intent(in)    :: dt                            ! land model time step (sec)

    !
    ! !LOCAL VARIABLES:
    integer   :: p,c,l,g,j,fp               ! indices
    integer   :: ier                           ! error code
    real(r8)  :: dwt                           ! change in pft weight (relative to column)
    real(r8)  :: dwt_leafc_seed(1:num_soilp_with_inactive) ! pft-level mass gain due to seeding of new area
    real(r8)  :: dwt_leafn_seed(1:num_soilp_with_inactive) ! pft-level mass gain due to seeding of new area
    real(r8)  :: dwt_deadstemc_seed(1:num_soilp_with_inactive) ! pft-level mass gain due to seeding of new area
    real(r8)  :: dwt_deadstemn_seed(1:num_soilp_with_inactive) ! pft-level mass gain due to seeding of new area
    real(r8)  :: dwt_npool_seed(1:num_soilp_with_inactive)     ! pft-level mass gain due to seeding of new area
    real(r8)  :: dwt_frootc_to_litter(1:num_soilp_with_inactive)       ! pft-level mass loss due to weight shift
    real(r8)  :: dwt_livecrootc_to_litter(1:num_soilp_with_inactive)   ! pft-level mass loss due to weight shift
    real(r8)  :: dwt_deadcrootc_to_litter(1:num_soilp_with_inactive)   ! pft-level mass loss due to weight shift
    real(r8)  :: dwt_frootn_to_litter(1:num_soilp_with_inactive)       ! pft-level mass loss due to weight shift
    real(r8)  :: dwt_livecrootn_to_litter(1:num_soilp_with_inactive)   ! pft-level mass loss due to weight shift
    real(r8)  :: dwt_deadcrootn_to_litter(1:num_soilp_with_inactive)   ! pft-level mass loss due to weight shift
    real(r8)  :: conv_cflux(1:num_soilp_with_inactive)                 ! pft-level mass loss due to weight shift
    real(r8)  :: prod10_cflux(1:num_soilp_with_inactive)               ! pft-level mass loss due to weight shift
    real(r8)  :: prod100_cflux(1:num_soilp_with_inactive)              ! pft-level mass loss due to weight shift
    real(r8)  :: crop_product_cflux(1:num_soilp_with_inactive)         ! pft-level mass loss due to weight shift
    real(r8)  :: conv_nflux(1:num_soilp_with_inactive)                 ! pft-level mass loss due to weight shift
    real(r8)  :: prod10_nflux(1:num_soilp_with_inactive)               ! pft-level mass loss due to weight shift
    real(r8)  :: prod100_nflux(1:num_soilp_with_inactive)              ! pft-level mass loss due to weight shift
    real(r8)  :: crop_product_nflux(1:num_soilp_with_inactive)         ! pft-level mass loss due to weight shift
    character(len=32)             :: subname='dyn_cbal'            ! subroutine name

    ! ! add phosphorus local variables
    real(r8)  :: dwt_leafp_seed(1:num_soilp_with_inactive)             ! pft-level mass gain due to seeding of new area
    real(r8)  :: dwt_deadstemp_seed(1:num_soilp_with_inactive)         ! pft-level mass gain due to seeding of new area
    real(r8)  :: dwt_ppool_seed(1:num_soilp_with_inactive)             ! pft-level mass gain due to seeding of new area
    real(r8)  :: dwt_frootp_to_litter(1:num_soilp_with_inactive)       ! pft-level mass loss due to weight shift
    real(r8)  :: dwt_livecrootp_to_litter(1:num_soilp_with_inactive)   ! pft-level mass loss due to weight shift
    real(r8)  :: dwt_deadcrootp_to_litter(1:num_soilp_with_inactive)   ! pft-level mass loss due to weight shift
    real(r8)  :: conv_pflux(1:num_soilp_with_inactive)                 ! pft-level mass loss due to weight shift
    real(r8)  :: prod10_pflux(1:num_soilp_with_inactive)               ! pft-level mass loss due to weight shift
    real(r8)  :: prod100_pflux(1:num_soilp_with_inactive)              ! pft-level mass loss due to weight shift
    real(r8)  :: crop_product_pflux(1:num_soilp_with_inactive)         ! pft-level mass loss due to weight shift

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

    do fp = 1, num_soilp_with_inactive
      ! initialize all the pft-level local flux arrays
      dwt_leafc_seed(fp)           = 0._r8
      dwt_deadstemc_seed(fp)       = 0._r8
      dwt_frootc_to_litter(fp)     = 0._r8
      dwt_livecrootc_to_litter(fp) = 0._r8
      dwt_deadcrootc_to_litter(fp) = 0._r8
      conv_cflux(fp)               = 0._r8
      prod10_cflux(fp)             = 0._r8
      prod100_cflux(fp)            = 0._r8
      crop_product_cflux(fp)       = 0._r8

      dwt_leafn_seed(fp)           = 0._r8
      dwt_deadstemn_seed(fp)       = 0._r8
      dwt_npool_seed(fp)           = 0._r8
      dwt_frootn_to_litter(fp)     = 0._r8
      dwt_livecrootn_to_litter(fp) = 0._r8
      dwt_deadcrootn_to_litter(fp) = 0._r8
      conv_nflux(fp)               = 0._r8
      prod10_nflux(fp)             = 0._r8
      prod100_nflux(fp)            = 0._r8
      crop_product_nflux(fp)       = 0._r8

      dwt_leafp_seed(fp)           = 0._r8
      dwt_deadstemp_seed(fp)       = 0._r8
      dwt_ppool_seed(fp)           = 0._r8
      dwt_frootp_to_litter(fp)     = 0._r8
      dwt_livecrootp_to_litter(fp) = 0._r8
      dwt_deadcrootp_to_litter(fp) = 0._r8
      conv_pflux(fp)               = 0._r8
      prod10_pflux(fp)             = 0._r8
      prod100_pflux(fp)            = 0._r8
      crop_product_pflux(fp)       = 0._r8
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

    do fp = 1, num_soilp_with_inactive
       p = filter_soilp_with_inactive(fp)
       c = veg_pp%column(p)
       l = veg_pp%landunit(p)
       if (.not.(lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) ) Then
         print *, "istsoil,istcrop",istsoil, istcrop
         print *, lun_pp%itype(l)
       end if
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
            call CarbonStateVarsInit     (veg_cs, p)
            call NitrogenStateVarsInit   (veg_ns, p)
            call PhosphorusStateVarsInit (veg_ps, p)
            call CanopyStateVarsInit     (canopystate_vars, p)
            call CNStateVarsInit         (cnstate_vars, p, c)
            call CarbonFluxVarsInit      (veg_cf, p)
            call NitrogenFluxVarsInit    ( p)
            call PhosphorusFluxVarsInit  ( p)

            if ( use_c13 ) then
               call CarbonStateVarsInit(c13_veg_cs, p)
            endif

            if ( use_c14 ) then
               call CarbonStateVarsInit(c14_veg_cs, p)
            endif

            ! add phosphorus related variables

            if ( use_c13 ) then
               photosyns_vars%alphapsnsun_patch(p) = 0._r8
               photosyns_vars%alphapsnsha_patch(p) = 0._r8
               photosyns_vars%rc13_canair_patch(p) = 0._r8
               photosyns_vars%rc13_psnsun_patch(p) = 0._r8
               photosyns_vars%rc13_psnsha_patch(p) = 0._r8
               photosyns_vars%c13_psnsun_patch(p) = 0._r8
               photosyns_vars%c13_psnsha_patch(p) = 0._r8

            endif
            photosyns_vars%psnsun_patch(p) = 0._r8
            photosyns_vars%psnsha_patch(p) = 0._r8
            if ( use_c14 ) then
               photosyns_vars%c14_psnsun_patch(p) = 0._r8
               photosyns_vars%c14_psnsha_patch(p) = 0._r8
            end if

        end if  ! end initialization of new pft
      end if       ! weight decreasing
    end do     ! patch loop

    do fp = 1, num_soilp_with_inactive
       p = filter_soilp_with_inactive(fp)
       c = veg_pp%column(p)
       l = veg_pp%landunit(p)

       call dyn_veg_cs_Adjustments(    &
            l, c, p,        &
            prior_weights,                 &
            patch_state_updater,           &
            dwt_leafc_seed(fp),                &
            dwt_deadstemc_seed(fp),            &
            conv_cflux(fp),                    &
            dwt_frootc_to_litter(fp),          &
            dwt_livecrootc_to_litter(fp),      &
            dwt_deadcrootc_to_litter(fp),      &
            prod10_cflux(fp),                  &
            prod100_cflux(fp),                 &
            crop_product_cflux(fp),            &
            veg_cs                         &
            )

    enddo

    if (use_c13) then
      do fp = 1, num_soilp_with_inactive
         p = filter_soilp_with_inactive(fp)
         c = veg_pp%column(p)
         l = veg_pp%landunit(p)
         call dyn_veg_cs_Adjustments( &
              l, c, p,        &
              prior_weights,                 &
              patch_state_updater,           &
              dwt_leafc13_seed(fp),              &
              dwt_deadstemc13_seed(fp),          &
              conv_c13flux(fp),                  &
              dwt_frootc13_to_litter(fp),        &
              dwt_livecrootc13_to_litter(fp),    &
              dwt_deadcrootc13_to_litter(fp),    &
              prod10_c13flux(fp),                &
              prod100_c13flux(fp),               &
              crop_product_c13flux(fp),          &
              c13_veg_cs                     &
              )
      enddo
    endif

    if (use_c14) then
      do fp = 1, num_soilp_with_inactive
         p = filter_soilp_with_inactive(fp)
         c = veg_pp%column(p)
         l = veg_pp%landunit(p)
         call dyn_veg_cs_Adjustments( &
              l, c, p,        &
              prior_weights,                 &
              patch_state_updater,           &
              dwt_leafc14_seed(fp),              &
              dwt_deadstemc14_seed(fp),          &
              conv_c14flux(fp),                  &
              dwt_frootc14_to_litter(fp),        &
              dwt_livecrootc14_to_litter(fp),    &
              dwt_deadcrootc14_to_litter(fp),    &
              prod10_c14flux(fp),                &
              prod100_c14flux(fp),               &
              crop_product_c14flux(fp),          &
              c14_veg_cs                     &
            )
      enddo
    endif


    do fp = 1, num_soilp_with_inactive
       p = filter_soilp_with_inactive(fp)
       c = veg_pp%column(p)
       l = veg_pp%landunit(p)
       call dyn_veg_ns_Adjustments(    &
            l,c,p,              &
            prior_weights,                 &
            patch_state_updater,           &
            dwt_leafn_seed(fp),             &
            dwt_deadstemn_seed(fp),         &
            dwt_npool_seed(fp),             &
            conv_nflux(fp),                 &
            dwt_frootn_to_litter(fp),       &
            dwt_livecrootn_to_litter(fp),   &
            dwt_deadcrootn_to_litter(fp),   &
            prod10_nflux(fp),               &
            prod100_nflux(fp),              &
            crop_product_nflux(fp),         &
            veg_ns                         &
            )
    enddo

    do fp = 1, num_soilp_with_inactive
       p = filter_soilp_with_inactive(fp)
       c = veg_pp%column(p)
       l = veg_pp%landunit(p)
       call dyn_veg_ps_Adjustments(    &
           l,c,p,                  &
           prior_weights,                 &
           patch_state_updater,           &
           dwt_leafp_seed(fp),            &
           dwt_deadstemp_seed(fp),        &
           dwt_ppool_seed(fp),            &
           conv_pflux(fp),                &
           dwt_frootp_to_litter(fp),      &
           dwt_livecrootp_to_litter(fp),  &
           dwt_deadcrootp_to_litter(fp),  &
           prod10_pflux(fp),              &
           prod100_pflux(fp),             &
           crop_product_pflux(fp),        &
           veg_ps                         &
           )
    end do

    ! calculate column-level seeding fluxes
    do fp = 1, num_soilp_with_inactive
       p = filter_soilp_with_inactive(fp)
       g = veg_pp%gridcell(p)

       ! C fluxes
       veg_cf%dwt_seedc_to_leaf(p) = dwt_leafc_seed(fp)/dt
       grc_cf%dwt_seedc_to_leaf(g)   = &
            grc_cf%dwt_seedc_to_leaf(g) + &
            veg_cf%dwt_seedc_to_leaf(p)

       veg_cf%dwt_seedc_to_deadstem(p) = dwt_deadstemc_seed(fp)/dt
       grc_cf%dwt_seedc_to_deadstem(g)   = &
            grc_cf%dwt_seedc_to_deadstem(g) + &
            veg_cf%dwt_seedc_to_deadstem(p)

       if ( use_c13 ) then
          c13_veg_cf%dwt_seedc_to_leaf(p) = dwt_leafc_seed(fp)/dt
          c13_grc_cf%dwt_seedc_to_leaf(g)   = &
               c13_grc_cf%dwt_seedc_to_leaf(g) + &
               c13_veg_cf%dwt_seedc_to_leaf(p)

          c13_veg_cf%dwt_seedc_to_deadstem(p) = dwt_deadstemc_seed(fp)/dt
          c13_grc_cf%dwt_seedc_to_deadstem(g)   = &
               c13_grc_cf%dwt_seedc_to_deadstem(g) + &
               c13_veg_cf%dwt_seedc_to_deadstem(p)
       endif

       if ( use_c14 ) then
          c14_veg_cf%dwt_seedc_to_leaf(p) = dwt_leafc_seed(fp)/dt
          c14_grc_cf%dwt_seedc_to_leaf(g)   = &
               c14_grc_cf%dwt_seedc_to_leaf(g) + &
               c14_veg_cf%dwt_seedc_to_leaf(p)

          c14_veg_cf%dwt_seedc_to_deadstem(p) = dwt_deadstemc_seed(fp)/dt
          c14_grc_cf%dwt_seedc_to_deadstem(g)   = &
               c14_grc_cf%dwt_seedc_to_deadstem(g) + &
               c14_veg_cf%dwt_seedc_to_deadstem(p)
       endif

       ! N fluxes
       veg_nf%dwt_seedn_to_leaf(p)   = dwt_leafn_seed(fp)/dt
       grc_nf%dwt_seedn_to_leaf(g)     = &
            grc_nf%dwt_seedn_to_leaf(g) + &
            veg_nf%dwt_seedn_to_leaf(p)

       veg_nf%dwt_seedn_to_deadstem(p) = dwt_deadstemn_seed(fp)/dt
       grc_nf%dwt_seedn_to_deadstem(g)   = &
            grc_nf%dwt_seedn_to_deadstem(g) + &
            veg_nf%dwt_seedn_to_deadstem(p)


       veg_nf%dwt_seedn_to_npool(p) = dwt_npool_seed(fp)/dt
       grc_nf%dwt_seedn_to_npool(g)   = &
            grc_nf%dwt_seedn_to_npool(g) + &
            veg_nf%dwt_seedn_to_npool(p)

       ! P fluxes
       veg_pf%dwt_seedp_to_leaf(p)   = dwt_leafp_seed(fp)/dt
       grc_pf%dwt_seedp_to_leaf(g)     = &
            grc_pf%dwt_seedp_to_leaf(g) + &
            veg_pf%dwt_seedp_to_leaf(p)

       veg_pf%dwt_seedp_to_deadstem(p) = dwt_deadstemp_seed(fp)/dt
       grc_pf%dwt_seedp_to_deadstem(g)   = &
            grc_pf%dwt_seedp_to_deadstem(g) + &
            veg_pf%dwt_seedp_to_deadstem(p)


       veg_pf%dwt_seedp_to_ppool(p) = dwt_npool_seed(fp)/dt
       grc_pf%dwt_seedp_to_ppool(g)   = &
            grc_pf%dwt_seedp_to_ppool(g) + &
            veg_pf%dwt_seedp_to_ppool(p)

    end do

    ! calculate patch-to-column slash fluxes into litter and CWD pools
    do fp = 1, num_soilp_with_inactive
       p = filter_soilp_with_inactive(fp)
       c = veg_pp%column(p)

       ! fine and coarse root to litter and CWD slash carbon fluxes
       col_cf%dwt_slash_cflux(c) =            &
            col_cf%dwt_slash_cflux(c)       + &
            dwt_frootc_to_litter(fp)     /dt + &
            dwt_livecrootc_to_litter(fp) /dt + &
            dwt_deadcrootc_to_litter(fp) /dt

       if ( use_c13 ) then
          c13_col_cf%dwt_slash_cflux(c) =          &
               c13_col_cf%dwt_slash_cflux(c)     + &
               dwt_frootc13_to_litter(fp)     /dt + &
               dwt_livecrootc13_to_litter(fp) /dt + &
               dwt_deadcrootc13_to_litter(fp) /dt
       endif

       if ( use_c14 ) then
          c14_col_cf%dwt_slash_cflux(c) =          &
               c14_col_cf%dwt_slash_cflux(c)     + &
               dwt_frootc14_to_litter(fp)     /dt + &
               dwt_livecrootc14_to_litter(fp) /dt + &
               dwt_deadcrootc14_to_litter(fp) /dt
       endif

       col_nf%dwt_slash_nflux(c) =            &
            col_nf%dwt_slash_nflux(c)       + &
            dwt_frootn_to_litter(fp)     /dt + &
            dwt_livecrootn_to_litter(fp) /dt + &
            dwt_deadcrootn_to_litter(fp) /dt

       col_pf%dwt_slash_pflux(c) =            &
            col_pf%dwt_slash_pflux(c)       + &
            dwt_frootp_to_litter(fp)     /dt + &
            dwt_livecrootp_to_litter(fp) /dt + &
            dwt_deadcrootp_to_litter(fp) /dt

    end do

    ! calculate pft-to-column for fluxes into litter and CWD pools
    do j = 1, nlevdecomp
      do fp = 1, num_soilp_with_inactive
         p = filter_soilp_with_inactive(fp)
         c = veg_pp%column(p)

          froot   = cnstate_vars%froot_prof_patch(p,j)
          croot   = cnstate_vars%croot_prof_patch(p,j)
          fr_flab = veg_vp%fr_flab(veg_pp%itype(p))
          fr_fcel = veg_vp%fr_fcel(veg_pp%itype(p))
          fr_flig = veg_vp%fr_flig(veg_pp%itype(p))

          ! fine root litter carbon fluxes
          col_cf%dwt_frootc_to_litr_met_c(c,j) = &
               col_cf%dwt_frootc_to_litr_met_c(c,j) + &
               (dwt_frootc_to_litter(fp)* fr_flab)/dt * froot

          col_cf%dwt_frootc_to_litr_cel_c(c,j) = &
               col_cf%dwt_frootc_to_litr_cel_c(c,j) + &
               (dwt_frootc_to_litter(fp)* fr_fcel)/dt * froot

          col_cf%dwt_frootc_to_litr_lig_c(c,j) = &
               col_cf%dwt_frootc_to_litr_lig_c(c,j) + &
               (dwt_frootc_to_litter(fp)* fr_flig)/dt * froot


          ! fine root litter nitrogen fluxes
          col_nf%dwt_frootn_to_litr_met_n(c,j) = &
               col_nf%dwt_frootn_to_litr_met_n(c,j) + &
               (dwt_frootn_to_litter(fp)* fr_flab)/dt * froot
          col_nf%dwt_frootn_to_litr_cel_n(c,j) = &

               col_nf%dwt_frootn_to_litr_cel_n(c,j) + &
               (dwt_frootn_to_litter(fp)* fr_fcel)/dt * froot

          col_nf%dwt_frootn_to_litr_lig_n(c,j) = &
               col_nf%dwt_frootn_to_litr_lig_n(c,j) + &
               (dwt_frootn_to_litter(fp)* fr_flig)/dt * froot

          ! fine root litter phosphorus fluxes
          col_pf%dwt_frootp_to_litr_met_p(c,j) = &
               col_pf%dwt_frootp_to_litr_met_p(c,j) + &
               (dwt_frootp_to_litter(fp)* fr_flab)/dt * froot
          col_pf%dwt_frootp_to_litr_cel_p(c,j) = &

               col_pf%dwt_frootp_to_litr_cel_p(c,j) + &
               (dwt_frootp_to_litter(fp)* fr_fcel)/dt * froot

          col_pf%dwt_frootp_to_litr_lig_p(c,j) = &
               col_pf%dwt_frootp_to_litr_lig_p(c,j) + &
               (dwt_frootp_to_litter(fp)* fr_flig)/dt * froot

          ! livecroot fluxes to cwd
          col_cf%dwt_livecrootc_to_cwdc(c,j) = &
               col_cf%dwt_livecrootc_to_cwdc(c,j) + &
               (dwt_livecrootc_to_litter(fp))/dt * croot

          col_nf%dwt_livecrootn_to_cwdn(c,j) = &
               col_nf%dwt_livecrootn_to_cwdn(c,j) + &
               (dwt_livecrootn_to_litter(fp))/dt * croot

          col_pf%dwt_livecrootp_to_cwdp(c,j) = &
               col_pf%dwt_livecrootp_to_cwdp(c,j) + &
               (dwt_livecrootp_to_litter(fp))/dt * croot

          ! deadcroot fluxes to cwd
          col_cf%dwt_deadcrootc_to_cwdc(c,j) = &
               col_cf%dwt_deadcrootc_to_cwdc(c,j) + &
               (dwt_deadcrootc_to_litter(fp))/dt * croot

          col_nf%dwt_deadcrootn_to_cwdn(c,j) = &
               col_nf%dwt_deadcrootn_to_cwdn(c,j) + &
               (dwt_deadcrootn_to_litter(fp))/dt * croot

          col_pf%dwt_deadcrootp_to_cwdp(c,j) = &
               col_pf%dwt_deadcrootp_to_cwdp(c,j) + &
               (dwt_deadcrootp_to_litter(fp))/dt * croot

          if ( use_c13 ) then
             ! C13 fine root litter fluxes
             c13_col_cf%dwt_frootc_to_litr_met_c(c,j) = &
                  c13_col_cf%dwt_frootc_to_litr_met_c(c,j) + &
                  (dwt_frootc13_to_litter(fp)* fr_flab)/dt * froot

             c13_col_cf%dwt_frootc_to_litr_cel_c(c,j) = &
                  c13_col_cf%dwt_frootc_to_litr_cel_c(c,j) + &
                  (dwt_frootc13_to_litter(fp)* fr_fcel)/dt * froot

             c13_col_cf%dwt_frootc_to_litr_lig_c(c,j) = &
                  c13_col_cf%dwt_frootc_to_litr_lig_c(c,j) + &
                  (dwt_frootc13_to_litter(fp)* fr_flig)/dt * froot

             ! livecroot fluxes to cwd
             c13_col_cf%dwt_livecrootc_to_cwdc(c,j) = &
                  c13_col_cf%dwt_livecrootc_to_cwdc(c,j) + &
                  (dwt_livecrootc13_to_litter(fp))/dt * croot

             ! deadcroot fluxes to cwd
             c13_col_cf%dwt_deadcrootc_to_cwdc(c,j) = &
                  c13_col_cf%dwt_deadcrootc_to_cwdc(c,j) + &
                  (dwt_deadcrootc13_to_litter(fp))/dt * croot

          endif

          if ( use_c14 ) then
             ! C14 fine root litter fluxes
             c14_col_cf%dwt_frootc_to_litr_met_c(c,j) = &
                  c14_col_cf%dwt_frootc_to_litr_met_c(c,j) + &
                  (dwt_frootc14_to_litter(fp)* fr_flab)/dt * froot

             c14_col_cf%dwt_frootc_to_litr_cel_c(c,j) = &
                  c14_col_cf%dwt_frootc_to_litr_cel_c(c,j) + &
                  (dwt_frootc14_to_litter(fp)* fr_fcel)/dt * froot

             c14_col_cf%dwt_frootc_to_litr_lig_c(c,j) = &
                  c14_col_cf%dwt_frootc_to_litr_lig_c(c,j) + &
                  (dwt_frootc14_to_litter(fp)* fr_flig)/dt * froot

             ! livecroot fluxes to cwd
             c14_col_cf%dwt_livecrootc_to_cwdc(c,j) = &
                  c14_col_cf%dwt_livecrootc_to_cwdc(c,j) + &
                  (dwt_livecrootc14_to_litter(fp))/dt * croot

             ! deadcroot fluxes to cwd
             c14_col_cf%dwt_deadcrootc_to_cwdc(c,j) = &
                  c14_col_cf%dwt_deadcrootc_to_cwdc(c,j) + &
                  (dwt_deadcrootc14_to_litter(fp))/dt * croot
          endif
       end do
    end do
    ! calculate pft-to-column for fluxes into product pools and conversion flux
    do fp = 1, num_soilp_with_inactive
       p = filter_soilp_with_inactive(fp)
       c = veg_pp%column(p)
       g = veg_pp%gridcell(p)

       ! column-level fluxes are accumulated as positive fluxes.
       ! column-level C flux updates
       col_cf%dwt_conv_cflux(c) = col_cf%dwt_conv_cflux(c) - conv_cflux(fp)/dt
       col_cf%dwt_prod10c_gain(c) = col_cf%dwt_prod10c_gain(c) - prod10_cflux(fp)/dt
       col_cf%dwt_prod100c_gain(c) = col_cf%dwt_prod100c_gain(c) - prod100_cflux(fp)/dt
       col_cf%dwt_crop_productc_gain(c) = col_cf%dwt_crop_productc_gain(c) - crop_product_cflux(fp)/dt

       veg_cf%dwt_prod10c_gain(p) = - prod10_cflux(fp)/dt
       grc_cf%dwt_prod10c_gain(g)   = grc_cf%dwt_prod10c_gain(g) + veg_cf%dwt_prod10c_gain(p)

       veg_cf%dwt_prod100c_gain(p) = - prod100_cflux(fp)/dt
       grc_cf%dwt_prod100c_gain(g)   = grc_cf%dwt_prod100c_gain(g) + veg_cf%dwt_prod100c_gain(p)

       veg_cf%dwt_crop_productc_gain(p) = - crop_product_cflux(fp)/dt

       if ( use_c13 ) then
          ! C13 column-level flux updates
          c13_col_cf%dwt_conv_cflux(c) = c13_col_cf%dwt_conv_cflux(c) - conv_c13flux(fp)/dt
          c13_col_cf%dwt_prod10c_gain(c) = c13_col_cf%dwt_prod10c_gain(c) - prod10_c13flux(fp)/dt
          c13_col_cf%dwt_prod100c_gain(c) = c13_col_cf%dwt_prod100c_gain(c) - prod100_c13flux(fp)/dt
          c13_col_cf%dwt_crop_productc_gain(c) = c13_col_cf%dwt_crop_productc_gain(c) - crop_product_c13flux(fp)/dt

          c13_veg_cf%dwt_prod10c_gain(p) = - prod10_c13flux(fp)/dt
          c13_grc_cf%dwt_prod10c_gain(g)   = c13_grc_cf%dwt_prod10c_gain(g) + c13_veg_cf%dwt_prod10c_gain(p)

          c13_veg_cf%dwt_prod100c_gain(p) = - prod100_c13flux(fp)/dt
          c13_grc_cf%dwt_prod100c_gain(g)   = c13_grc_cf%dwt_prod100c_gain(g) + c13_veg_cf%dwt_prod100c_gain(p)

          c13_veg_cf%dwt_crop_productc_gain(p) = - crop_product_c13flux(fp)/dt

       endif

       if ( use_c14 ) then
          ! C14 column-level flux updates
          c14_col_cf%dwt_conv_cflux(c) = c14_col_cf%dwt_conv_cflux(c) - conv_c14flux(fp)/dt
          c14_col_cf%dwt_prod10c_gain(c) = c14_col_cf%dwt_prod10c_gain(c) - prod10_c14flux(fp)/dt
          c14_col_cf%dwt_prod100c_gain(c) = c14_col_cf%dwt_prod100c_gain(c) - prod100_c14flux(fp)/dt
          c14_col_cf%dwt_crop_productc_gain(c) = c14_col_cf%dwt_crop_productc_gain(c) - crop_product_c14flux(fp)/dt

          c14_veg_cf%dwt_prod10c_gain(p) = - prod10_c14flux(fp)/dt
          c14_grc_cf%dwt_prod10c_gain(g)   = c14_grc_cf%dwt_prod10c_gain(g) + c14_veg_cf%dwt_prod10c_gain(p)

          c14_veg_cf%dwt_prod100c_gain(p) = - prod100_c14flux(fp)/dt
          c14_grc_cf%dwt_prod100c_gain(g)   = c14_grc_cf%dwt_prod100c_gain(g) + c14_veg_cf%dwt_prod100c_gain(p)

          c14_veg_cf%dwt_crop_productc_gain(p) = - crop_product_c14flux(fp)/dt

       endif

       ! column-level N flux updates
       col_nf%dwt_conv_nflux(c) = col_nf%dwt_conv_nflux(c) - conv_nflux(fp)/dt
       col_nf%dwt_prod10n_gain(c) = col_nf%dwt_prod10n_gain(c) - prod10_nflux(fp)/dt
       col_nf%dwt_prod100n_gain(c) = col_nf%dwt_prod100n_gain(c) - prod100_nflux(fp)/dt
       col_nf%dwt_crop_productn_gain(c) = col_nf%dwt_crop_productn_gain(c) - crop_product_nflux(fp)/dt

       veg_nf%dwt_prod10n_gain(p) = -prod10_nflux(fp)/dt
       grc_nf%dwt_prod10n_gain(g)   = grc_nf%dwt_prod10n_gain(g) + veg_nf%dwt_prod10n_gain(p)

       veg_nf%dwt_prod100n_gain(p)= -prod100_nflux(fp)/dt
       grc_nf%dwt_prod100n_gain(g)  = grc_nf%dwt_prod100n_gain(g) + veg_nf%dwt_prod100n_gain(p)

       veg_nf%dwt_crop_productn_gain(p) = -crop_product_nflux(fp)/dt

       ! column-level P flux updates

       col_pf%dwt_conv_pflux(c) = col_pf%dwt_conv_pflux(c) - conv_pflux(fp)/dt
       col_pf%dwt_prod10p_gain(c) = col_pf%dwt_prod10p_gain(c) - prod10_pflux(fp)/dt
       col_pf%dwt_prod100p_gain(c) = col_pf%dwt_prod100p_gain(c) - prod100_pflux(fp)/dt
       col_pf%dwt_crop_productp_gain(c) = col_pf%dwt_crop_productp_gain(c) - crop_product_pflux(fp)/dt

       veg_pf%dwt_prod10p_gain(p) = -prod10_pflux(fp)/dt
       grc_pf%dwt_prod10p_gain(g)   = grc_pf%dwt_prod10p_gain(g) + veg_pf%dwt_prod10p_gain(p)

       veg_pf%dwt_prod100p_gain(p)= -prod100_pflux(fp)/dt
       grc_pf%dwt_prod100p_gain(g)  = grc_pf%dwt_prod100p_gain(g) + veg_pf%dwt_prod100p_gain(p)

       veg_pf%dwt_crop_productp_gain(p) = -crop_product_pflux(fp)/dt

    end do

    do fp = 1, num_soilp_with_inactive
       p = filter_soilp_with_inactive(fp)
       g = veg_pp%gridcell(p)

       ! Note that patch-level fluxes are stored per unit GRIDCELL area - thus, we don't
       ! need to multiply by the patch's gridcell weight when translating patch-level
       ! fluxes into gridcell-level fluxes.

       veg_cf%dwt_conv_cflux(p) = -conv_cflux(fp)/dt
       grc_cf%dwt_conv_cflux(g) = &
            grc_cf%dwt_conv_cflux(g) + &
            veg_cf%dwt_conv_cflux(p)

       if ( use_c13 ) then
          ! C13 column-level flux updates
          c13_veg_cf%dwt_conv_cflux(p) = -conv_c13flux(fp)/dt
          c13_grc_cf%dwt_conv_cflux(g) = &
               c13_grc_cf%dwt_conv_cflux(g) + &
               c13_veg_cf%dwt_conv_cflux(p)
       endif

       if ( use_c14 ) then
          ! C14 column-level flux updates
          c14_veg_cf%dwt_conv_cflux(p) = -conv_c14flux(fp)/dt
          c14_grc_cf%dwt_conv_cflux(g) = &
               c14_grc_cf%dwt_conv_cflux(g) + &
               c14_veg_cf%dwt_conv_cflux(p)
       endif

       veg_nf%dwt_conv_nflux(p) = -conv_nflux(fp)/dt
       grc_nf%dwt_conv_nflux(g) = &
            grc_nf%dwt_conv_nflux(g) + &
            veg_nf%dwt_conv_nflux(p)

       veg_pf%dwt_conv_pflux(p) = -conv_pflux(fp)/dt
       grc_pf%dwt_conv_pflux(g) = &
            grc_pf%dwt_conv_pflux(g) + &
            veg_pf%dwt_conv_pflux(p)

    end do

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
 subroutine NitrogenFluxVarsInit(p)
   !
   ! !DESCRIPTION:
   ! Initializes p-th patch of nitrogenflux_type
   !
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

   character(len=*), parameter :: subname = 'dyn_cnbal_col'
   !-----------------------------------------------------------------------

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

end subroutine dyn_cnbal_column

 end module dynConsBiogeochemMod
