module dynSubgridAdjustmentsMod

#include "shr_assert.h"

  !---------------------------------------------------------------------------
  !
  ! !DESCRIPTION:
  ! Holds the methods for adjusting state variables at each sub-grid level,
  ! based on weight changes and other information stored in the 
  ! state_updater structures.
  ! These subroutines were originally contained in the sub-grid level
  ! data type definitions, but moved here to keep all the dynamic subgrid
  ! functionality together.
  !
  ! !USES:
  use shr_kind_mod           , only : r8 => shr_kind_r8
  use shr_infnan_mod         , only : nan => shr_infnan_nan, assignment(=)
  use decompMod              , only : bounds_type
  use clm_varpar             , only : ndecomp_pools, nlevdecomp
  use clm_varctl             , only : use_crop
  use clm_varcon             , only : dzsoi_decomp
  use dynPatchStateUpdaterMod, only : patch_state_updater_type
  use dynColumnStateUpdaterMod, only : column_state_updater_type
  use ColumnDataType         , only : column_carbon_state
  use VegetationDataType     , only : vegetation_carbon_state
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  save
  
  public :: dynColumnAdjustments
  public :: dynVegetationAdjustments
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine dynVegetationAdjustments(      &
       bounds,                              &
       num_filterp_with_inactive,           &
       filterp_with_inactive,               &
       prior_weights,                       &
       patch_state_updater,                 &
       dwt_leafc_seed,                      &
       dwt_deadstemc_seed,                  &
       conv_cflux,                          &
       dwt_frootc_to_litter,                &
       dwt_livecrootc_to_litter,            &
       dwt_deadcrootc_to_litter,            &
       prod10_cflux,                        &
       prod100_cflux,                       &
       crop_product_cflux,                  &
       veg_cs                               &
       )
    !
    ! !DESCRIPTION:
    ! Adjust state variables and compute associated fluxes when patch areas change due to
    ! dynamic landuse
    !
    ! !USES:
    use pftvarcon          , only : pconv, pprod10, pprod100
    use dynPriorWeightsMod , only : prior_weights_type
    use landunit_varcon    , only : istsoil, istcrop
    use CNComputeSeedMod   , only : ComputeSeedAmounts
    !
    ! !ARGUMENTS:
    type(bounds_type)              , intent(in)    :: bounds
    integer                        , intent(in)    :: num_filterp_with_inactive ! number of points in filterp_with_inactive
    integer                        , intent(in)    :: filterp_with_inactive(:) ! patch filter that includes inactive points
    type(prior_weights_type)       , intent(in)    :: prior_weights
    type(patch_state_updater_type) , intent(in)    :: patch_state_updater
    real(r8)                       , intent(inout) :: dwt_leafc_seed           (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_deadstemc_seed       (bounds%begp:)
    real(r8)                       , intent(inout) :: conv_cflux               (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_frootc_to_litter     (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_livecrootc_to_litter (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_deadcrootc_to_litter (bounds%begp:)
    real(r8)                       , intent(inout) :: prod10_cflux             (bounds%begp:)
    real(r8)                       , intent(inout) :: prod100_cflux            (bounds%begp:)
    real(r8)                       , intent(inout) :: crop_product_cflux       (bounds%begp:)
    type(vegetation_carbon_state)  , intent(inout) :: veg_cs
    !
    ! !LOCAL VARIABLES:
    integer                     :: begp, endp
    integer                     :: l, c, p
    logical                     :: old_weight_was_zero(bounds%begp:bounds%endp)
    logical                     :: patch_grew(bounds%begp:bounds%endp)

    ! The following are only set for growing patches:
    real(r8)                    :: seed_leafc_patch(bounds%begp:bounds%endp)
    real(r8)                    :: seed_leafc_storage_patch(bounds%begp:bounds%endp)
    real(r8)                    :: seed_leafc_xfer_patch(bounds%begp:bounds%endp)
    real(r8)                    :: seed_deadstemc_patch(bounds%begp:bounds%endp)

    real(r8)                    :: wood_product_cflux(bounds%begp:bounds%endp)
    real(r8)                    :: deadstemc_patch_temp(bounds%begp:bounds%endp)

    character(len=*), parameter :: subname = 'CStateDynamicPatchAdjustments'
    !-----------------------------------------------------------------------

    begp = bounds%begp
    endp = bounds%endp

    SHR_ASSERT_ALL((ubound(dwt_leafc_seed           ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_deadstemc_seed       ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(conv_cflux               ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_frootc_to_litter     ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_livecrootc_to_litter ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_deadcrootc_to_litter ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(prod10_cflux             ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(prod100_cflux            ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(crop_product_cflux       ) == (/endp/)), errMsg(__FILE__, __LINE__))
   
    old_weight_was_zero = patch_state_updater%old_weight_was_zero(bounds)
    patch_grew          = patch_state_updater%patch_grew(bounds)

    call ComputeSeedAmounts(bounds                                        , &
         species                    = veg_cs%species                          , &
         leaf_patch                 = veg_cs%leafc(begp:endp)           , &
         leaf_storage_patch         = veg_cs%leafc_storage(begp:endp)   , &
         leaf_xfer_patch            = veg_cs%leafc_xfer(begp:endp)      , &

         ! Calculations only needed for patches that grew:
         compute_here_patch         = patch_grew(begp:endp)               , &

         ! For patches that previously had zero area, ignore the current state for the
         ! sake of computing leaf proportions:
         ignore_current_state_patch = old_weight_was_zero(begp:endp)      , &

         seed_leaf_patch            = seed_leafc_patch(begp:endp)         , &
         seed_leaf_storage_patch    = seed_leafc_storage_patch(begp:endp) , &
         seed_leaf_xfer_patch       = seed_leafc_xfer_patch(begp:endp)    , &
         seed_deadstem_patch        = seed_deadstemc_patch(begp:endp))

    ! 1) LEAFC_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%leafc   (begp:endp)           , &
         flux_out_grc_area = conv_cflux       (begp:endp)           , &
         seed              = seed_leafc_patch (begp:endp)           , &
         seed_addition     = dwt_leafc_seed   (begp:endp))


    ! 2) LEAFC_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%leafc_storage   (begp:endp)   , &
         flux_out_grc_area = conv_cflux               (begp:endp)   , &
         seed              = seed_leafc_storage_patch (begp:endp)   , &
         seed_addition     = dwt_leafc_seed           (begp:endp))

    ! 3) LEAF_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%leafc_xfer   (begp:endp)      , &
         flux_out_grc_area = conv_cflux            (begp:endp)      , &
         seed              = seed_leafc_xfer_patch (begp:endp)      , &
         seed_addition     = dwt_leafc_seed        (begp:endp))

    ! 4) FROOTC_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%frootc(begp:endp)             , &
         flux_out_col_area = dwt_frootc_to_litter(begp:endp))
    
    ! 5) FROOTC_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%frootc_storage(begp:endp)     , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 6) FROOTC_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%frootc_xfer(begp:endp)        , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 7) LIVESTEMC_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%livestemc(begp:endp)          , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 8) LIVESTEMC_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%livestemc_storage(begp:endp)  , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 9) LIVESTEMC_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%livestemc_xfer(begp:endp)     , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 10) PROD10_CFLUX
    wood_product_cflux(begp:endp)      = 0._r8
    deadstemc_patch_temp(begp:endp)    = veg_cs%deadstemc(begp:endp)
    call patch_state_updater%update_patch_state_partition_flux_by_type( &
         bounds                                                        , &
         num_filterp_with_inactive                                     , &
         filterp_with_inactive                                         , &
         flux1_fraction_by_pft_type = pprod10                          , &
         var                        = deadstemc_patch_temp    (begp:endp) , &
         flux1_out                  = prod10_cflux            (begp:endp) , &
         flux2_out                  = wood_product_cflux      (begp:endp) , &
         seed                       = seed_deadstemc_patch    (begp:endp) )

    ! 11) PROD100_CFLUX
    wood_product_cflux(begp:endp)      = 0._r8
    deadstemc_patch_temp(begp:endp)    = veg_cs%deadstemc(begp:endp)
    call patch_state_updater%update_patch_state_partition_flux_by_type( &
         bounds                                                        , &
         num_filterp_with_inactive                                     , &
         filterp_with_inactive                                         , &
         flux1_fraction_by_pft_type = pprod100                         , &
         var                        = deadstemc_patch_temp    (begp:endp) , &
         flux1_out                  = prod100_cflux            (begp:endp) , &
         flux2_out                  = wood_product_cflux      (begp:endp) , &
         seed                       = seed_deadstemc_patch    (begp:endp))

    ! 12) DEADSTEMC_PATCH
    wood_product_cflux(begp:endp)      = 0._r8
    call patch_state_updater%update_patch_state_partition_flux_by_type( &
         bounds                                                        , &
         num_filterp_with_inactive                                     , &
         filterp_with_inactive                                         , &
         flux1_fraction_by_pft_type = pconv                          , &
         var                        = veg_cs%deadstemc(begp:endp) , &
         flux1_out                  = conv_cflux              (begp:endp) , &
         flux2_out                  = wood_product_cflux       (begp:endp) , &
         seed                       = seed_deadstemc_patch    (begp:endp) , &
         seed_addition              = dwt_deadstemc_seed     (begp:endp))

    ! 13) DEADSTEMC_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%deadstemc_storage(begp:endp)  , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 14) DEADSTEMC_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%deadstemc_xfer(begp:endp)     , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 15) LIVECROOTC_PATCH 
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%livecrootc(begp:endp)         , &
         flux_out_col_area = dwt_livecrootc_to_litter(begp:endp))

    ! 16) LIVECROOTC_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%livecrootc_storage(begp:endp) , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 17) LIVECROOTC_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%livecrootc_xfer(begp:endp)    , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 18) DEADCROOTC_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%deadcrootc(begp:endp)         , &
         flux_out_col_area = dwt_deadcrootc_to_litter(begp:endp))

    ! 19) DEADCROOTC_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%deadcrootc_storage(begp:endp) , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 20) DEADCROOT_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%deadcrootc_xfer(begp:endp)    , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 21) GRESP_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%gresp_storage(begp:endp)      , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 22) GRESP_XFER_STORAGE
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%gresp_xfer(begp:endp)         , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 23) CPOOL_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%cpool(begp:endp)              , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 24) XSMRPOOL_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%xsmrpool(begp:endp)           , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 25) CTRUNC_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%ctrunc(begp:endp)             , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 26) DISPVEGC_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%dispvegc(begp:endp))

    ! 27) STORVEGC_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%storvegc(begp:endp))

    ! 28) TOTVEGC_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%totvegc(begp:endp))

    ! 29) TOTPFTC_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%totpftc(begp:endp))

    if (use_crop) then
       ! This is a negative pool. So any deficit that we haven't repaid gets sucked out
       ! of the atmosphere.
       call patch_state_updater%update_patch_state(         &
            bounds                                        , &
            num_filterp_with_inactive                     , &
            filterp_with_inactive                         , &
            var = veg_cs%cropseedc_deficit(begp:endp) , &
            flux_out_grc_area = conv_cflux(begp:endp))
    end if

    ! These fluxes are computed as negative quantities, but are expected to be positive,
    ! so flip the signs
    do p = begp, endp
       dwt_frootc_to_litter(p)     = -1._r8 * dwt_frootc_to_litter(p)
       dwt_livecrootc_to_litter(p) = -1._r8 * dwt_livecrootc_to_litter(p)
       dwt_deadcrootc_to_litter(p) = -1._r8 * dwt_deadcrootc_to_litter(p)
    end do

  end subroutine dynVegetationAdjustments
  
  !-----------------------------------------------------------------------
  subroutine dynColumnAdjustments(bounds, clump_index, column_state_updater, col_cs)
    !
    ! !DESCRIPTION:
    ! Adjust state variables and compute associated fluxes when patch areas change due to
    ! dynamic landuse
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(bounds_type)               , intent(in)    :: bounds
    integer                         , intent(in)    :: clump_index
    type(column_state_updater_type) , intent(in)    :: column_state_updater
    type(column_carbon_state)       , intent(in)    :: col_cs
    !
    ! !LOCAL VARIABLES:
    integer                     :: l, j
    integer                     :: begc, endc
    real(r8)                    :: adjustment_one_level(bounds%begc:bounds%endc)

    character(len=*), parameter :: subname = 'CStateDynamicColumnAdjustments'
    !-----------------------------------------------------------------------

    begc = bounds%begc
    endc = bounds%endc

    col_cs%dyn_cbal_adjustments(begc:endc) = 0._r8

    do l = 1, ndecomp_pools
       do j = 1, nlevdecomp
          call column_state_updater%update_column_state_no_special_handling( &
               bounds      = bounds,                                         &
               clump_index = clump_index,                                    &
               var         = col_cs%decomp_cpools_vr(begc:endc, j, l),     &
               adjustment  = adjustment_one_level(begc:endc))

          col_cs%dyn_cbal_adjustments(begc:endc) = &
               col_cs%dyn_cbal_adjustments(begc:endc) + &
               adjustment_one_level(begc:endc) * dzsoi_decomp(j)

       end do
    end do

    do j = 1, nlevdecomp
       call column_state_updater%update_column_state_no_special_handling( &
            bounds      = bounds,                                         &
            clump_index = clump_index,                                    &
            var         = col_cs%ctrunc_vr(begc:endc,j),     &
            adjustment  = adjustment_one_level(begc:endc))

       col_cs%dyn_cbal_adjustments(begc:endc) = &
            col_cs%dyn_cbal_adjustments(begc:endc) + &
            adjustment_one_level(begc:endc) * dzsoi_decomp(j)

    end do

  end subroutine dynColumnAdjustments
  
end module dynSubgridAdjustmentsMod
