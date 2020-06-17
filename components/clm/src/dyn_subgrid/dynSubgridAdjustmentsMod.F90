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
  use decompMod              , only : bounds_type
  use clm_varpar             , only : ndecomp_pools, nlevdecomp
  use clm_varctl             , only : use_crop
  use clm_varcon             , only : dzsoi_decomp
  use dynPatchStateUpdaterMod, only : patch_state_updater_type
  use dynColumnStateUpdaterMod, only : column_state_updater_type
  use ColumnDataType         , only : column_carbon_state, column_nitrogen_state
  use ColumnDataType         , only : column_phosphorus_state
  use VegetationDataType     , only : vegetation_carbon_state, vegetation_nitrogen_state
  use VegetationDataType     , only : vegetation_phosphorus_state
  use SpeciesMod             , only : CN_SPECIES_N, CN_SPECIES_P
  !---Currently contains module class methods for gpu
  use dynUpdateModAcc
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  save

  real(r8) , parameter :: npool_seed_param     = 0.1_r8
  real(r8) , parameter :: ppool_seed_param     = 0.01_r8
  !$acc declare copyin(npool_seed_param, ppool_seed_param)

  public :: dyn_veg_cs_Adjustments
  public :: dyn_col_cs_Adjustments
  public :: dyn_veg_ns_Adjustments
  public :: dyn_col_ns_Adjustments
  public :: dyn_veg_ps_Adjustments
  public :: dyn_col_ps_Adjustments
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine dyn_veg_cs_Adjustments(        &
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
    ! Adjust veg-level state variables and compute associated fluxes when patch areas change due to
    ! dynamic landuse
    !
    ! !USES:
      !$acc routine seq
    use pftvarcon          , only : pconv, pprod10, pprod100
    use dynPriorWeightsMod , only : prior_weights_type
    use landunit_varcon    , only : istsoil, istcrop
    use ComputeSeedMod   , only : ComputeSeedAmounts
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
    !-----------------------------------------------------------------------
    associate(&
      leafc         => veg_cs%leafc         , &
      leafc_storage => veg_cs%leafc_storage , &
      leafc_xfer    => veg_cs%leafc_xfer    , &
      frootc        => veg_cs%frootc                ,&
      frootc_storage => veg_cs%frootc_storage  ,&
      frootc_xfer    => veg_cs%frootc_xfer    , &
      livestemc      => veg_cs%livestemc                ,&
      livestemc_storage => veg_cs%livestemc_storage  ,&
      livestemc_xfer    => veg_cs%livestemc_xfer    , &
      deadstemc         => veg_cs%deadstemc                ,&
      deadstemc_storage => veg_cs%deadstemc_storage  ,&
      deadstemc_xfer    => veg_cs%deadstemc_xfer    , &
      livecrootc         => veg_cs%livecrootc          ,&
      livecrootc_storage => veg_cs%livecrootc_storage  ,&
      livecrootc_xfer    => veg_cs%livecrootc_xfer    , &
      deadcrootc         => veg_cs%deadcrootc                ,&
      deadcrootc_storage => veg_cs%deadcrootc_storage  ,&
      deadcrootc_xfer    => veg_cs%deadcrootc_xfer    , &
      gresp_storage      => veg_cs%gresp_storage  ,&
      gresp_xfer         => veg_cs%gresp_xfer    , &
      cpool              => veg_cs%cpool    , &
      xsmrpool           => veg_cs%xsmrpool    , &
      ctrunc             => veg_cs%ctrunc      , &
      dispvegc           => veg_cs%dispvegc  , &
      storvegc           => veg_cs%storvegc  , &
      totvegc            => veg_cs%totvegc   , &
      totpftc            => veg_cs%totpftc   , &
      cropseedc_deficit  => veg_cs%cropseedc_deficit &
      )
    begp = bounds%begp
    endp = bounds%endp


    old_weight_was_zero = old_weight_was_zeroAcc(patch_state_updater, bounds)
    patch_grew          = patch_grewAcc(patch_state_updater,bounds)

    call ComputeSeedAmounts(bounds                               , &
         species                    = veg_cs%species             , &
         leaf_patch                 = leafc        (begp:endp)   , &
         leaf_storage_patch         = leafc_storage(begp:endp)   , &
         leaf_xfer_patch            = leafc_xfer   (begp:endp)   , &
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
    call update_patch_stateAcc(patch_state_updater,          &
         bounds                                              , &
         num_filterp_with_inactive                           , &
         filterp_with_inactive                               , &
         var               = leafc(begp:endp)       , &
         flux_out_grc_area = conv_cflux(begp:endp)           , &
         seed              = seed_leafc_patch(begp:endp)     , &
         seed_addition     = dwt_leafc_seed(begp:endp)   )


    ! 2) LEAFC_STORAGE_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = leafc_storage(begp:endp)      , &
         flux_out_grc_area = conv_cflux(begp:endp)                  , &
         seed              = seed_leafc_storage_patch(begp:endp)    , &
         seed_addition     = dwt_leafc_seed(begp:endp)           )

    ! 3) LEAF_XFER_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = leafc_xfer(begp:endp)         , &
         flux_out_grc_area = conv_cflux(begp:endp)                  , &
         seed              = seed_leafc_xfer_patch(begp:endp)       , &
         seed_addition     = dwt_leafc_seed(begp:endp)        )

    ! 4) FROOTC_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = frootc(begp:endp)             , &
         flux_out_col_area = dwt_frootc_to_litter(begp:endp))

    ! 5) FROOTC_STORAGE_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = frootc_storage(begp:endp)     , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 6) FROOTC_XFER_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = frootc_xfer(begp:endp)        , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 7) LIVESTEMC_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = livestemc(begp:endp)          , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 8) LIVESTEMC_STORAGE_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = livestemc_storage(begp:endp)  , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 9) LIVESTEMC_XFER_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = livestemc_xfer(begp:endp)     , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 10) PROD10_CFLUX
    wood_product_cflux(begp:endp)      = 0._r8
    deadstemc_patch_temp(begp:endp)    = deadstemc(begp:endp)
    call update_patch_state_partition_flux_by_typeAcc(&
         patch_state_updater    ,&
         bounds                                                        , &
         num_filterp_with_inactive                                     , &
         filterp_with_inactive                                         , &
         flux1_fraction_by_pft_type = pprod10                           , &
         var                        = deadstemc_patch_temp(begp:endp)     , &
         flux1_out                  = prod10_cflux(begp:endp)             , &
         flux2_out                  = wood_product_cflux(begp:endp)       , &
         seed                       = seed_deadstemc_patch(begp:endp)     )

    ! 11) PROD100_CFLUX
    wood_product_cflux(begp:endp)      = 0._r8
    deadstemc_patch_temp(begp:endp)    = deadstemc(begp:endp)
    call update_patch_state_partition_flux_by_typeAcc(&
         patch_state_updater    ,&
         bounds                                                        , &
         num_filterp_with_inactive                                     , &
         filterp_with_inactive                                         , &
         flux1_fraction_by_pft_type = pprod100                            , &
         var                        = deadstemc_patch_temp (begp:endp)      , &
         flux1_out                  = prod100_cflux        (begp:endp)       , &
         flux2_out                  = wood_product_cflux   (begp:endp)      , &
         seed                       = seed_deadstemc_patch (begp:endp)     )

    ! 12) DEADSTEMC_PATCH
    wood_product_cflux(begp:endp)      = 0._r8
    call update_patch_state_partition_flux_by_typeAcc(&
         patch_state_updater    ,&
         bounds                                                        , &
         num_filterp_with_inactive                                     , &
         filterp_with_inactive                                         , &
         flux1_fraction_by_pft_type = pconv                          , &
         var                        = deadstemc(begp:endp) , &
         flux1_out                  = conv_cflux(begp:endp)               , &
         flux2_out                  = wood_product_cflux(begp:endp)        , &
         seed                       = seed_deadstemc_patch(begp:endp)     , &
         seed_addition              = dwt_deadstemc_seed(begp:endp)     )

    ! 13) DEADSTEMC_STORAGE_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = deadstemc_storage(begp:endp)  , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 14) DEADSTEMC_XFER_PATCH
     call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = deadstemc_xfer(begp:endp)     , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 15) LIVECROOTC_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = livecrootc(begp:endp)         , &
         flux_out_col_area = dwt_livecrootc_to_litter(begp:endp))

    ! 16) LIVECROOTC_STORAGE_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = livecrootc_storage(begp:endp) , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 17) LIVECROOTC_XFER_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = livecrootc_xfer(begp:endp)    , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 18) DEADCROOTC_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = deadcrootc (begp:endp)        , &
         flux_out_col_area = dwt_deadcrootc_to_litter(begp:endp))

    ! 19) DEADCROOTC_STORAGE_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = deadcrootc_storage(begp:endp) , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 20) DEADCROOT_XFER_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = deadcrootc_xfer(begp:endp)    , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 21) GRESP_STORAGE_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = gresp_storage(begp:endp)      , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 22) GRESP_XFER_STORAGE
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = gresp_xfer(begp:endp)         , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 23) CPOOL_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = cpool(begp:endp)              , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 24) XSMRPOOL_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = xsmrpool(begp:endp)           , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 25) CTRUNC_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ctrunc(begp:endp)             , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 26) DISPVEGC_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = dispvegc(begp:endp))

    ! 27) STORVEGC_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = storvegc(begp:endp))

    ! 28) TOTVEGC_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = totvegc(begp:endp))

    ! 29) TOTPFTC_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = totpftc(begp:endp))

    if (use_crop) then
       ! This is a negative pool. So any deficit that we haven't repaid gets sucked out
       ! of the atmosphere.
       call update_patch_stateAcc(patch_state_updater,                 &
            bounds                                        , &
            num_filterp_with_inactive                     , &
            filterp_with_inactive                         , &
            var = cropseedc_deficit(begp:endp) , &
            flux_out_grc_area = conv_cflux(begp:endp))
    end if

    ! These fluxes are computed as negative quantities, but are expected to be positive,
    ! so flip the signs
    do p = begp, endp
       dwt_frootc_to_litter(p)     = -1._r8 * dwt_frootc_to_litter(p)
       dwt_livecrootc_to_litter(p) = -1._r8 * dwt_livecrootc_to_litter(p)
       dwt_deadcrootc_to_litter(p) = -1._r8 * dwt_deadcrootc_to_litter(p)
    end do
   end associate
  end subroutine dyn_veg_cs_Adjustments

  !-----------------------------------------------------------------------
  subroutine dyn_col_cs_Adjustments(bounds, clump_index, column_state_updater, col_cs)
    !
    ! !DESCRIPTION:
    ! Adjust column-level carbon state variables and compute associated fluxes
    ! when patch areas change due to dynamic landuse
    !
    ! !USES:
    !
    ! !ARGUMENTS:
      !$acc routine seq
    type(bounds_type)               , intent(in)    :: bounds
    integer                         , intent(in)    :: clump_index
    type(column_state_updater_type) , intent(in)    :: column_state_updater
    type(column_carbon_state)       , intent(inout) :: col_cs
    !
    ! !LOCAL VARIABLES:
    integer                     :: l, j
    integer                     :: begc, endc
    real(r8)                    :: adjustment_one_level(bounds%begc:bounds%endc)
    !-----------------------------------------------------------------------
    associate(&
      decomp_cpools_vr    => col_cs%decomp_cpools_vr , &
      ctrunc_vr           => col_cs%ctrunc_vr       &
      )
    begc = bounds%begc
    endc = bounds%endc

    col_cs%dyn_cbal_adjustments(begc:endc) = 0._r8

    do l = 1, ndecomp_pools
       do j = 1, nlevdecomp
          call update_column_state_no_special_handling_acc( column_state_updater , &
               bounds      = bounds,                                         &
               clump_index = clump_index,                                    &
               var         = decomp_cpools_vr(begc:endc, j, l),     &
               adjustment  = adjustment_one_level(begc:endc))

          col_cs%dyn_cbal_adjustments(begc:endc) = &
               col_cs%dyn_cbal_adjustments(begc:endc) + &
               adjustment_one_level(begc:endc) * dzsoi_decomp(j)

       end do
    end do

    do j = 1, nlevdecomp
       call update_column_state_no_special_handling_acc( column_state_updater , &
            bounds      = bounds,                                         &
            clump_index = clump_index,                                    &
            var         = ctrunc_vr(begc:endc,j),     &
            adjustment  = adjustment_one_level(begc:endc))

       col_cs%dyn_cbal_adjustments(begc:endc) = &
            col_cs%dyn_cbal_adjustments(begc:endc) + &
            adjustment_one_level(begc:endc) * dzsoi_decomp(j)

    end do
    end associate
  end subroutine dyn_col_cs_Adjustments

  !-----------------------------------------------------------------------
  subroutine dyn_veg_ns_Adjustments(        &
       bounds,                              &
       num_filterp_with_inactive,           &
       filterp_with_inactive,               &
       prior_weights,                       &
       patch_state_updater,                 &
       dwt_leafn_seed,                      &
       dwt_deadstemn_seed,                  &
       dwt_npool_seed,                      &
       conv_nflux,                          &
       dwt_frootn_to_litter,                &
       dwt_livecrootn_to_litter,            &
       dwt_deadcrootn_to_litter,            &
       prod10_nflux,                        &
       prod100_nflux,                       &
       crop_product_nflux,                  &
       veg_ns                               &
       )
    !
    ! !DESCRIPTION:
    ! Adjust state variables and compute associated fluxes when patch areas change due to
    ! dynamic landuse
    !
    ! !USES:
      !$acc routine seq
    use pftvarcon          , only : pconv, pprod10, pprod100
    use dynPriorWeightsMod , only : prior_weights_type
    use landunit_varcon    , only : istsoil, istcrop
    use ComputeSeedMod   , only : ComputeSeedAmounts
    !
    ! !ARGUMENTS:
    type(bounds_type)              , intent(in)    :: bounds
    integer                        , intent(in)    :: num_filterp_with_inactive
    integer                        , intent(in)    :: filterp_with_inactive(:)
    type(prior_weights_type)       , intent(in)    :: prior_weights
    type(patch_state_updater_type) , intent(in)    :: patch_state_updater
    real(r8)                       , intent(inout) :: dwt_leafn_seed           (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_deadstemn_seed       (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_npool_seed           (bounds%begp:)
    real(r8)                       , intent(inout) :: conv_nflux               (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_frootn_to_litter     (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_livecrootn_to_litter (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_deadcrootn_to_litter (bounds%begp:)
    real(r8)                       , intent(inout) :: prod10_nflux             (bounds%begp:)
    real(r8)                       , intent(inout) :: prod100_nflux            (bounds%begp:)
    real(r8)                       , intent(inout) :: crop_product_nflux       (bounds%begp:)
    type(vegetation_nitrogen_state), intent(inout) :: veg_ns
    !
    ! !LOCAL VARIABLES:
    integer                     :: begp, endp
    integer                     :: l, c, p
    logical                     :: old_weight_was_zero      (bounds%begp:bounds%endp)
    logical                     :: patch_grew               (bounds%begp:bounds%endp)

    ! The following are only set for growing patches:
    real(r8)                    :: seed_leafn_patch         (bounds%begp:bounds%endp)
    real(r8)                    :: seed_leafn_storage_patch (bounds%begp:bounds%endp)
    real(r8)                    :: seed_leafn_xfer_patch    (bounds%begp:bounds%endp)
    real(r8)                    :: seed_deadstemn_patch     (bounds%begp:bounds%endp)
    real(r8)                    :: seed_npool_patch         (bounds%begp:bounds%endp)

    real(r8)                    :: wood_product_nflux       (bounds%begp:bounds%endp)
    real(r8)                    :: deadstemn_patch_temp     (bounds%begp:bounds%endp)

    !character(len=*), parameter :: subname = 'NStateDynamicPatchAdjustments'
    !-----------------------------------------------------------------------
    associate(&
      leafn          => veg_ns%leafn         , &
      leafn_storage  => veg_ns%leafn_storage , &
      leafn_xfer     => veg_ns%leafn_xfer    , &
      frootn         => veg_ns%frootn                ,&
      frootn_storage => veg_ns%frootn_storage  ,&
      frootn_xfer    => veg_ns%frootn_xfer    , &
      livestemn         => veg_ns%livestemn                ,&
      livestemn_storage => veg_ns%livestemn_storage  ,&
      livestemn_xfer    => veg_ns%livestemn_xfer    , &
      deadstemn         => veg_ns%deadstemn                ,&
      deadstemn_storage => veg_ns%deadstemn_storage  ,&
      deadstemn_xfer    => veg_ns%deadstemn_xfer    , &
      livecrootn         => veg_ns%livecrootn          ,&
      livecrootn_storage => veg_ns%livecrootn_storage  ,&
      livecrootn_xfer    => veg_ns%livecrootn_xfer    , &
      deadcrootn         => veg_ns%deadcrootn                ,&
      deadcrootn_storage => veg_ns%deadcrootn_storage  ,&
      deadcrootn_xfer    => veg_ns%deadcrootn_xfer    , &
      retransn           => veg_ns%retransn  ,&
      npool              => veg_ns%npool    , &
      ntrunc             => veg_ns%ntrunc      , &
      dispvegn           => veg_ns%dispvegn  , &
      storvegn           => veg_ns%storvegn  , &
      totvegn            => veg_ns%totvegn   , &
      totpftn            => veg_ns%totpftn   , &
      cropseedn_deficit  => veg_ns%cropseedn_deficit &
      )
    begp = bounds%begp
    endp = bounds%endp


    old_weight_was_zero = old_weight_was_zeroAcc(patch_state_updater,bounds)
    patch_grew          = patch_grewAcc(patch_state_updater,bounds)

    call ComputeSeedAmounts(bounds                                        , &
         species                    = CN_SPECIES_N                        , &
         leaf_patch                 = leafn         (begp:endp)  , &
         leaf_storage_patch         = leafn_storage (begp:endp)  , &
         leaf_xfer_patch            = leafn_xfer    (begp:endp)  , &

         ! Calculations only needed for patches that grew:
         compute_here_patch         = patch_grew(begp:endp)               , &

         ! For patches that previously had zero area, ignore the current state for the
         ! sake of computing leaf proportions:
         ignore_current_state_patch = old_weight_was_zero(begp:endp)      , &

         seed_leaf_patch            = seed_leafn_patch(begp:endp)         , &
         seed_leaf_storage_patch    = seed_leafn_storage_patch(begp:endp) , &
         seed_leaf_xfer_patch       = seed_leafn_xfer_patch(begp:endp)    , &
         seed_deadstem_patch        = seed_deadstemn_patch(begp:endp)     , &
         pool_seed_param            = npool_seed_param                    , &
         pool_seed_patch            = seed_npool_patch(begp:endp))

    ! 1) LEAFN_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                           , &
         num_filterp_with_inactive                        , &
         filterp_with_inactive                            , &
         var               = leafn            (begp:endp)   , &
         flux_out_grc_area = conv_nflux       (begp:endp) , &
         seed              = seed_leafn_patch (begp:endp) , &
         seed_addition     = dwt_leafn_seed   (begp:endp))

    ! 2) LEAFN_STORAGE_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                   , &
         num_filterp_with_inactive                                , &
         filterp_with_inactive                                    , &
         var               = leafn_storage(begp:endp)    , &
         flux_out_grc_area = conv_nflux(begp:endp)                , &
         seed              = seed_leafn_storage_patch(begp:endp)  , &
         seed_addition     = dwt_leafn_seed(begp:endp)           )

    ! 3) LEAFN_XFER_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                        , &
         num_filterp_with_inactive                                     , &
         filterp_with_inactive                                         , &
         var               = leafn_xfer(begp:endp)   , &
         flux_out_grc_area = conv_nflux(begp:endp)            , &
         seed              = seed_leafn_xfer_patch(begp:endp) , &
         seed_addition     = dwt_leafn_seed(begp:endp)        )

    ! 4) FROOTN_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = frootn(begp:endp)             , &
         flux_out_col_area = dwt_frootn_to_litter(begp:endp))

    ! 5) FROOTN_STORAGE_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = frootn_storage(begp:endp)     , &
         flux_out_grc_area = conv_nflux(begp:endp)  )

    ! 6) FROOTN_XFER_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                       , &
         num_filterp_with_inactive                    , &
         filterp_with_inactive                        , &
         var               = frootn_xfer(begp:endp)        , &
         flux_out_grc_area = conv_nflux (begp:endp)   )

    ! 7) LIVESTEMN_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = livestemn(begp:endp)          , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 8) LIVESTEMN_STORAGE_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = livestemn_storage(begp:endp)  , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 9) LIVESTEMN_XFER_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = livestemn_xfer(begp:endp)     , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 10) PROD10_NFLUX
    wood_product_nflux(begp:endp)      = 0._r8
    deadstemn_patch_temp(begp:endp)    = deadstemn(begp:endp)
    call update_patch_state_partition_flux_by_typeAcc(patch_state_updater , &
         bounds                                                           , &
         num_filterp_with_inactive                                        , &
         filterp_with_inactive                                            , &
         flux1_fraction_by_pft_type = pprod10                  , &
         var                        = deadstemn_patch_temp(begp:endp)     , &
         flux1_out                  = prod10_nflux(begp:endp)             , &
         flux2_out                  = wood_product_nflux(begp:endp)       , &
         seed                       = seed_deadstemn_patch(begp:endp)     )

    ! 11) PROD100_NFLUX
    wood_product_nflux(begp:endp)      = 0._r8
    deadstemn_patch_temp(begp:endp)    = deadstemn(begp:endp)
    call update_patch_state_partition_flux_by_typeAcc(patch_state_updater , &
         bounds                                                           , &
         num_filterp_with_inactive                                        , &
         filterp_with_inactive                                            , &
         flux1_fraction_by_pft_type = pprod100                 , &
         var                        = deadstemn_patch_temp(begp:endp)     , &
         flux1_out                  = prod100_nflux       (begp:endp)     , &
         flux2_out                  = wood_product_nflux  (begp:endp)     , &
         seed                       = seed_deadstemn_patch(begp:endp)    )

    ! 12) DEADSTEMN_PATCH
    wood_product_nflux(begp:endp)      = 0._r8
    call update_patch_state_partition_flux_by_typeAcc(patch_state_updater , &
         bounds                                                           , &
         num_filterp_with_inactive                                        , &
         filterp_with_inactive                                            , &
         flux1_fraction_by_pft_type = pconv                               , &
         var                        = deadstemn            (begp:endp) , &
         flux1_out                  = conv_nflux           (begp:endp)    , &
         flux2_out                  = wood_product_nflux   (begp:endp)    , &
         seed                       = seed_deadstemn_patch (begp:endp)    , &
         seed_addition              = dwt_deadstemn_seed   (begp:endp))

    ! 13) DEADSTEMN_STORAGE_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = deadstemn_storage(begp:endp)  , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 14) DEADSTEMN_XFER_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = deadstemn_xfer (begp:endp)    , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 15) LIVECROOTN_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = livecrootn (begp:endp)        , &
         flux_out_col_area = dwt_livecrootn_to_litter(begp:endp))

    ! 16) LIVECROOTN_STORAGE_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = livecrootn_storage(begp:endp) , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 17) LIVECROOTN_XFER_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = livecrootn_xfer(begp:endp)    , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 18) DEADCROOTN_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = deadcrootn (begp:endp)        , &
         flux_out_col_area = dwt_deadcrootn_to_litter(begp:endp))

    ! 19) DEADCROOTN_STORAGE_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = deadcrootn_storage(begp:endp) , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 20) DEADCROOT_XFER_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = deadcrootn_xfer(begp:endp)    , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 21) RETRANSN_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = retransn(begp:endp)           , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 22) NTRUNC_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ntrunc    (begp:endp)         , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 23) CPOOL_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = npool   (begp:endp)           , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 24) DISPVEGN_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = dispvegn(begp:endp))

    ! 25) STORVEGN_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = storvegn(begp:endp))

    ! 26) TOTVEGN_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = totvegn(begp:endp))

    ! 27) TOTPFTN_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = totpftn(begp:endp))

    if (use_crop) then
       ! This is a negative pool. So any deficit that we haven't repaid gets sucked out
       ! of the atmosphere.
       call update_patch_stateAcc(patch_state_updater,       &
            bounds                                        , &
            num_filterp_with_inactive                     , &
            filterp_with_inactive                         , &
            var = cropseedn_deficit(begp:endp) , &
            flux_out_grc_area = conv_nflux(begp:endp))
    end if

    ! These fluxes are computed as negative quantities, but are expected to be positive,
    ! so flip the signs
    do p = begp,endp
       dwt_frootn_to_litter(p)     = -1._r8 * dwt_frootn_to_litter(p)
       dwt_livecrootn_to_litter(p) = -1._r8 * dwt_livecrootn_to_litter(p)
       dwt_deadcrootn_to_litter(p) = -1._r8 * dwt_deadcrootn_to_litter(p)
    end do
   end associate
  end subroutine dyn_veg_ns_Adjustments

  !-----------------------------------------------------------------------
  subroutine dyn_col_ns_Adjustments(bounds, clump_index, column_state_updater, col_ns)
    !
    ! !DESCRIPTION:
    ! Adjust state variables and compute associated fluxes when patch areas change due to
    ! dynamic landuse
    !
    ! !USES:
      !$acc routine seq
    use dynPriorWeightsMod       , only : prior_weights_type
    use landunit_varcon          , only : istsoil, istcrop
    use dynColumnStateUpdaterMod , only : column_state_updater_type
        !
    ! !ARGUMENTS:
    type(bounds_type)               , intent(in)    :: bounds
    integer                         , intent(in)    :: clump_index
    type(column_state_updater_type) , intent(in)    :: column_state_updater
    type(column_nitrogen_state)     , intent(inout) :: col_ns
    !
    ! !LOCAL VARIABLES:
    integer                     :: l, j
    integer                     :: begc, endc
    real(r8)                    :: adjustment_one_level(bounds%begc:bounds%endc)
    !-----------------------------------------------------------------------
    associate(&
      decomp_npools_vr    => col_ns%decomp_npools_vr , &
      ntrunc_vr           => col_ns%ntrunc_vr      , &
      sminn_vr            => col_ns%sminn_vr &
      )

    begc = bounds%begc
    endc = bounds%endc
    col_ns%dyn_nbal_adjustments(begc:endc) = 0._r8
    do l = 1, ndecomp_pools
       do j = 1, nlevdecomp
          call update_column_state_no_special_handling_acc(column_state_updater, &
               bounds      = bounds,                                         &
               clump_index = clump_index,                                    &
               var         = decomp_npools_vr(begc:endc, j, l),     &
               adjustment  = adjustment_one_level(begc:endc))

          col_ns%dyn_nbal_adjustments(begc:endc) = &
               col_ns%dyn_nbal_adjustments(begc:endc) + &
               adjustment_one_level(begc:endc) * dzsoi_decomp(j)

       end do
    end do

    do j = 1, nlevdecomp
       call update_column_state_no_special_handling_acc(column_state_updater, &
            bounds      = bounds,                                         &
            clump_index = clump_index,                                    &
            var         = ntrunc_vr(begc:endc,j),     &
            adjustment  = adjustment_one_level(begc:endc))

       col_ns%dyn_nbal_adjustments(begc:endc) = &
            col_ns%dyn_nbal_adjustments(begc:endc) + &
            adjustment_one_level(begc:endc) * dzsoi_decomp(j)


       call update_column_state_no_special_handling_acc(column_state_updater, &
           bounds      = bounds                          , &
           clump_index = clump_index                     , &
           var         = sminn_vr(begc:endc, j), &
           adjustment  = adjustment_one_level(begc:endc))

       col_ns%dyn_nbal_adjustments(begc:endc) = &
           col_ns%dyn_nbal_adjustments(begc:endc) + &
           adjustment_one_level(begc:endc) * dzsoi_decomp(j)

    end do
    end associate
  end subroutine dyn_col_ns_Adjustments

  !-----------------------------------------------------------------------
  subroutine dyn_veg_ps_Adjustments(        &
       bounds,                              &
       num_filterp_with_inactive,           &
       filterp_with_inactive,               &
       prior_weights,                       &
       patch_state_updater,                 &
       dwt_leafp_seed,                      &
       dwt_deadstemp_seed,                  &
       dwt_ppool_seed,                      &
       conv_pflux,                          &
       dwt_frootp_to_litter,                &
       dwt_livecrootp_to_litter,            &
       dwt_deadcrootp_to_litter,            &
       prod10_pflux,                        &
       prod100_pflux,                       &
       crop_product_pflux,                  &
       veg_ps                               &
       )
    !
    ! !DESCRIPTION:
    ! Adjust state variables and compute associated fluxes when patch areas change due to
    ! dynamic landuse
    !
    ! !USES:
      !$acc routine seq
    use pftvarcon          , only : pconv, pprod10, pprod100
    use dynPriorWeightsMod , only : prior_weights_type
    use landunit_varcon    , only : istsoil, istcrop
    use ComputeSeedMod   , only : ComputeSeedAmounts
    !
    ! !ARGUMENTS:
    type(bounds_type)                , intent(in)    :: bounds
    integer                          , intent(in)    :: num_filterp_with_inactive
    integer                          , intent(in)    :: filterp_with_inactive(:)
    type(prior_weights_type)         , intent(in)    :: prior_weights
    type(patch_state_updater_type)   , intent(in)    :: patch_state_updater
    real(r8)                         , intent(inout) :: dwt_leafp_seed           (bounds%begp:)
    real(r8)                         , intent(inout) :: dwt_deadstemp_seed       (bounds%begp:)
    real(r8)                         , intent(inout) :: dwt_ppool_seed           (bounds%begp:)
    real(r8)                         , intent(inout) :: conv_pflux               (bounds%begp:)
    real(r8)                         , intent(inout) :: dwt_frootp_to_litter     (bounds%begp:)
    real(r8)                         , intent(inout) :: dwt_livecrootp_to_litter (bounds%begp:)
    real(r8)                         , intent(inout) :: dwt_deadcrootp_to_litter (bounds%begp:)
    real(r8)                         , intent(inout) :: prod10_pflux             (bounds%begp:)
    real(r8)                         , intent(inout) :: prod100_pflux            (bounds%begp:)
    real(r8)                         , intent(inout) :: crop_product_pflux       (bounds%begp:)
    type(vegetation_phosphorus_state), intent(inout) :: veg_ps
    !
    ! !LOCAL VARIABLES:
    integer                     :: begp, endp
    integer                     :: l, c, p
    logical                     :: old_weight_was_zero      (bounds%begp:bounds%endp)
    logical                     :: patch_grew               (bounds%begp:bounds%endp)

    ! The following are only set for growing patches:
    real(r8)                    :: seed_leafp_patch         (bounds%begp:bounds%endp)
    real(r8)                    :: seed_leafp_storage_patch (bounds%begp:bounds%endp)
    real(r8)                    :: seed_leafp_xfer_patch    (bounds%begp:bounds%endp)
    real(r8)                    :: seed_deadstemp_patch     (bounds%begp:bounds%endp)
    real(r8)                    :: seed_ppool_patch         (bounds%begp:bounds%endp)

    real(r8)                    :: wood_product_pflux       (bounds%begp:bounds%endp)
    real(r8)                    :: deadstemp_patch_temp     (bounds%begp:bounds%endp)

    !character(len=*), parameter :: subname = 'dyn_veg_ps_Adjustments'
    !-----------------------------------------------------------------------
    associate(&
      leafp          => veg_ps%leafp         , &
      leafp_storage  => veg_ps%leafp_storage , &
      leafp_xfer     => veg_ps%leafp_xfer    , &
      frootp         => veg_ps%frootp                ,&
      frootp_storage => veg_ps%frootp_storage  ,&
      frootp_xfer    => veg_ps%frootp_xfer    , &
      livestemp         => veg_ps%livestemp                ,&
      livestemp_storage => veg_ps%livestemp_storage  ,&
      livestemp_xfer    => veg_ps%livestemp_xfer    , &
      deadstemp         => veg_ps%deadstemp                ,&
      deadstemp_storage => veg_ps%deadstemp_storage  ,&
      deadstemp_xfer    => veg_ps%deadstemp_xfer    , &
      livecrootp         => veg_ps%livecrootp          ,&
      livecrootp_storage => veg_ps%livecrootp_storage  ,&
      livecrootp_xfer    => veg_ps%livecrootp_xfer    , &
      deadcrootp         => veg_ps%deadcrootp                ,&
      deadcrootp_storage => veg_ps%deadcrootp_storage  ,&
      deadcrootp_xfer    => veg_ps%deadcrootp_xfer    , &
      retransp           => veg_ps%retransp  ,&
      ppool              => veg_ps%ppool    , &
      ptrunc             => veg_ps%ptrunc      , &
      dispvegp           => veg_ps%dispvegp  , &
      storvegp           => veg_ps%storvegp  , &
      totvegp            => veg_ps%totvegp   , &
      totpftp            => veg_ps%totpftp   , &
      cropseedp_deficit  => veg_ps%cropseedp_deficit &
      )
    begp = bounds%begp
    endp = bounds%endp


    old_weight_was_zero = old_weight_was_zeroAcc(patch_state_updater,bounds)
    patch_grew          = patch_grewAcc(patch_state_updater,bounds)

    call ComputeSeedAmounts(bounds                                        , &
         species                    = CN_SPECIES_P                        , &
         leaf_patch                 = leafp        (begp:endp)     , &
         leaf_storage_patch         = leafp_storage(begp:endp)     , &
         leaf_xfer_patch            = leafp_xfer   (begp:endp)     , &
         ! Calculations only needed for patches that grew:
         compute_here_patch         = patch_grew(begp:endp)               , &
         ! For patches that previously had zero area, ignore the current state for the
         ! sake of computing leaf proportions:
         ignore_current_state_patch = old_weight_was_zero(begp:endp)      , &

         seed_leaf_patch            = seed_leafp_patch(begp:endp)         , &
         seed_leaf_storage_patch    = seed_leafp_storage_patch(begp:endp) , &
         seed_leaf_xfer_patch       = seed_leafp_xfer_patch(begp:endp)    , &
         seed_deadstem_patch        = seed_deadstemp_patch(begp:endp)     , &
         pool_seed_param            = ppool_seed_param                    , &
         pool_seed_patch            = seed_ppool_patch(begp:endp))

    ! 1) LEAFP_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                           , &
         num_filterp_with_inactive                        , &
         filterp_with_inactive                            , &
         var               = leafp           (begp:endp)  , &
         flux_out_grc_area = conv_pflux      (begp:endp)  , &
         seed              = seed_leafp_patch(begp:endp)  , &
         seed_addition     = dwt_leafp_seed  (begp:endp) )

    ! 2) LEAFP_STORAGE_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                   , &
         num_filterp_with_inactive                                , &
         filterp_with_inactive                                    , &
         var               = leafp_storage           (begp:endp)  , &
         flux_out_grc_area = conv_pflux              (begp:endp)  , &
         seed              = seed_leafp_storage_patch(begp:endp)  , &
         seed_addition     = dwt_leafp_seed          (begp:endp) )

    ! 3) LEAFP_XFER_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                        , &
         num_filterp_with_inactive                                     , &
         filterp_with_inactive                                         , &
         var               = leafp_xfer    (begp:endp) , &
         flux_out_grc_area = conv_pflux           (begp:endp) , &
         seed              = seed_leafp_xfer_patch(begp:endp) , &
         seed_addition     = dwt_leafp_seed       (begp:endp) )

    ! 4) FROOTP_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = frootp       (begp:endp)      , &
         flux_out_col_area = dwt_frootp_to_litter(begp:endp))

    ! 5) FROOTP_STORAGE_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = frootp_storage(begp:endp)     , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 6) FROOTP_XFER_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = frootp_xfer(begp:endp)        , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 7) PROD10_PFLUX
    wood_product_pflux(begp:endp)      = 0._r8
    deadstemp_patch_temp(begp:endp)    = deadstemp(begp:endp)
    call update_patch_state_partition_flux_by_typeAcc(patch_state_updater, &
         bounds                                                           , &
         num_filterp_with_inactive                                        , &
         filterp_with_inactive                                            , &
         flux1_fraction_by_pft_type = pprod10                             , &
         var                        = deadstemp_patch_temp(begp:endp)     , &
         flux1_out                  = prod10_pflux        (begp:endp)     , &
         flux2_out                  = wood_product_pflux  (begp:endp)     , &
         seed                       = seed_deadstemp_patch(begp:endp)     )

    ! 8) PROD100_PFLUX
    wood_product_pflux(begp:endp)      = 0._r8
    deadstemp_patch_temp(begp:endp)    = deadstemp(begp:endp)
    call update_patch_state_partition_flux_by_typeAcc(patch_state_updater, &
         bounds                                                           , &
         num_filterp_with_inactive                                        , &
         filterp_with_inactive                                            , &
         flux1_fraction_by_pft_type = pprod100                            , &
         var                        = deadstemp_patch_temp(begp:endp)     , &
         flux1_out                  = prod100_pflux       (begp:endp)     , &
         flux2_out                  = wood_product_pflux  (begp:endp)     , &
         seed                       = seed_deadstemp_patch(begp:endp)    )

    ! 9) DEADSTEMP_PATCH
    wood_product_pflux(begp:endp)      = 0._r8
    call update_patch_state_partition_flux_by_typeAcc(patch_state_updater, &
         bounds                                                           , &
         num_filterp_with_inactive                                        , &
         filterp_with_inactive                                            , &
         flux1_fraction_by_pft_type = pconv                               , &
         var                        = deadstemp           (begp:endp)  , &
         flux1_out                  = conv_pflux          (begp:endp)     , &
         flux2_out                  = wood_product_pflux  (begp:endp)     , &
         seed                       = seed_deadstemp_patch(begp:endp)     , &
         seed_addition              = dwt_deadstemp_seed  (begp:endp) )

    ! 10) DEADSTEMP_STORAGE_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = deadstemp_storage(begp:endp)  , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 11) DEADSTEMP_XFER_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = deadstemp_xfer(begp:endp)     , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 12) LIVESTEMP_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = livestemp(begp:endp)          , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 13) LIVESTEMP_STORAGE_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = livestemp_storage(begp:endp)  , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 14) LIVESTEMP_XFER_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = livestemp_xfer(begp:endp)     , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 15) LIVECROOTP_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = livecrootp(begp:endp)         , &
         flux_out_col_area = dwt_livecrootp_to_litter(begp:endp))

    ! 16) LIVECROOTP_STORAGE_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = livecrootp_storage(begp:endp) , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 17) LIVECROOTP_XFER_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = livecrootp_xfer(begp:endp)    , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 18) DEADCROOTP_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = deadcrootp(begp:endp)         , &
         flux_out_col_area = dwt_deadcrootp_to_litter(begp:endp))

    ! 19) DEADCROOTP_STORAGE_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = deadcrootp_storage(begp:endp) , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 20) DEADCROOT_XFER_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = deadcrootp_xfer (begp:endp)   , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 21) RETRANSP_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = retransp (begp:endp)          , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 22) PTRUNC_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ptrunc(begp:endp)             , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 23) PPOOL_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ppool(begp:endp)              , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 24) DISPVEGP_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = dispvegp(begp:endp))

    ! 25) STORVEGP_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = storvegp(begp:endp))

    ! 26) TOTVEGP_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = totvegp(begp:endp))

    ! 27) TOTPFTP_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = totpftp(begp:endp))

    if (use_crop) then
       ! This is a negative pool. So any deficit that we haven't repaid gets sucked out
       ! of the atmosphere.
       call update_patch_stateAcc(patch_state_updater,       &
            bounds                                        , &
            num_filterp_with_inactive                     , &
            filterp_with_inactive                         , &
            var = cropseedp_deficit(begp:endp) , &
            flux_out_grc_area = conv_pflux(begp:endp))
    end if

    ! These fluxes are computed as negative quantities, but are expected to be positive,
    ! so flip the signs
    do p = begp,endp
       dwt_frootp_to_litter(p)     = -1._r8 * dwt_frootp_to_litter(p)
       dwt_livecrootp_to_litter(p) = -1._r8 * dwt_livecrootp_to_litter(p)
       dwt_deadcrootp_to_litter(p) = -1._r8 * dwt_deadcrootp_to_litter(p)
    end do
    end associate
  end subroutine dyn_veg_ps_Adjustments

  !-----------------------------------------------------------------------
  subroutine dyn_col_ps_Adjustments(bounds, clump_index, column_state_updater, col_ps)
    !
    ! !DESCRIPTION:
    ! Adjust state variables and compute associated fluxes when patch areas change due to
    ! dynamic landuse
    !
    ! !USES:
      !$acc routine seq
    use dynPriorWeightsMod       , only : prior_weights_type
    use landunit_varcon          , only : istsoil, istcrop
    use dynColumnStateUpdaterMod , only : column_state_updater_type
    !
    ! !ARGUMENTS:
    type(bounds_type)               , intent(in)    :: bounds
    integer                         , intent(in)    :: clump_index
    type(column_state_updater_type) , intent(in)    :: column_state_updater
    type(column_phosphorus_state)   , intent(inout) :: col_ps
    !
    ! !LOCAL VARIABLES:
    integer           :: l, j
    integer           :: begc, endc
    real(r8)          :: adjustment_one_level(bounds%begc:bounds%endc)
    !-----------------------------------------------------------------------
    associate(&
      decomp_ppools_vr    => col_ps%decomp_ppools_vr , &
      ptrunc_vr           => col_ps%ptrunc_vr      , &
      solutionp_vr        => col_ps%solutionp_vr , &
      labilep_vr          => col_ps%labilep_vr , &
      secondp_vr          => col_ps%secondp_vr &
      )
    begc = bounds%begc
    endc = bounds%endc

    col_ps%dyn_pbal_adjustments(begc:endc) = 0._r8
    do l = 1, ndecomp_pools
       do j = 1, nlevdecomp

         !col_slice = col_ps%decomp_ppools_vr(begc:endc,j,l)

          call update_column_state_no_special_handling_acc( column_state_updater, &
               bounds      = bounds,                                         &
               clump_index = clump_index,                                    &
               var         = decomp_ppools_vr(begc:endc, j, l),     &
               adjustment  = adjustment_one_level(begc:endc) )

          col_ps%dyn_pbal_adjustments(begc:endc) =      &
               col_ps%dyn_pbal_adjustments(begc:endc) + &
               adjustment_one_level(begc:endc) * dzsoi_decomp(j)

       end do
    end do

    do j = 1, nlevdecomp
      !col_slice = col_ps%ptrunc_vr(begc:endc,j)
       call update_column_state_no_special_handling_acc( column_state_updater, &
            bounds      = bounds,                                         &
            clump_index = clump_index,                                    &
            var         = ptrunc_vr(begc:endc,j),                &
            adjustment  = adjustment_one_level(begc:endc))

       col_ps%dyn_pbal_adjustments(begc:endc) =      &
            col_ps%dyn_pbal_adjustments(begc:endc) + &
            adjustment_one_level(begc:endc) * dzsoi_decomp(j)
       !!
       !col_slice = col_ps%solutionp_vr(begc:endc,j)
       call update_column_state_no_special_handling_acc( column_state_updater, &
            bounds      = bounds,                                         &
            clump_index = clump_index,                                    &
            var         = solutionp_vr(begc:endc,j),             &
            adjustment  = adjustment_one_level(begc:endc))

       col_ps%dyn_pbal_adjustments(begc:endc) =      &
            col_ps%dyn_pbal_adjustments(begc:endc) + &
            adjustment_one_level(begc:endc) * dzsoi_decomp(j)
        !!
       call update_column_state_no_special_handling_acc( column_state_updater, &
            bounds      = bounds,                                         &
            clump_index = clump_index,                                    &
            var         = labilep_vr(begc:endc,j),               &
            adjustment  = adjustment_one_level(begc:endc))

       col_ps%dyn_pbal_adjustments(begc:endc) =      &
            col_ps%dyn_pbal_adjustments(begc:endc) + &
            adjustment_one_level(begc:endc) * dzsoi_decomp(j)
       !!
       call update_column_state_no_special_handling_acc( column_state_updater, &
            bounds      = bounds,                                         &
            clump_index = clump_index,                                    &
            var         = secondp_vr(begc:endc,j),               &
            adjustment  = adjustment_one_level(begc:endc))

       col_ps%dyn_pbal_adjustments(begc:endc) =      &
            col_ps%dyn_pbal_adjustments(begc:endc) + &
            adjustment_one_level(begc:endc) * dzsoi_decomp(j)

    end do
    end associate

  end subroutine dyn_col_ps_Adjustments

  !-----------------------------------------------------------------------
  subroutine dyn_col_ps_Adjustments_acc(bounds, clump_index, column_state_updater, col_ps)
    !
    ! !DESCRIPTION:
    ! Adjust state variables and compute associated fluxes when patch areas change due to
    ! dynamic landuse
    !$acc routine seq
    ! !USES:
    use dynPriorWeightsMod       , only : prior_weights_type
    use landunit_varcon          , only : istsoil, istcrop
    use dynColumnStateUpdaterMod , only : column_state_updater_type
    !
    ! !ARGUMENTS:
    type(bounds_type)               , intent(in)    :: bounds
    integer                         , intent(in)    :: clump_index
    type(column_state_updater_type) , intent(in)    :: column_state_updater
    type(column_phosphorus_state)   , intent(inout) :: col_ps
    !
    ! !LOCAL VARIABLES:
    integer                     :: l, j
    integer                     :: begc, endc
    real(r8)                    :: adjustment_one_level(bounds%begc:bounds%endc)
    !-----------------------------------------------------------------------
    associate(&
      decomp_ppools_vr    => col_ps%decomp_ppools_vr , &
      ptrunc_vr           => col_ps%ptrunc_vr      , &
      solutionp_vr        => col_ps%solutionp_vr , &
      labilep_vr          => col_ps%labilep_vr , &
      secondp_vr          => col_ps%secondp_vr &
      )
    begc = bounds%begc
    endc = bounds%endc
    !allocate(col_slice(begc:endc))

    col_ps%dyn_pbal_adjustments(begc:endc) = 0._r8
    do l = 1, ndecomp_pools
       do j = 1, nlevdecomp

         !col_slice = col_ps%decomp_ppools_vr(begc:endc,j,l)

          call update_column_state_no_special_handling_acc( column_state_updater, &
               bounds      = bounds,                                         &
               clump_index = clump_index,                                    &
               var         = decomp_ppools_vr(begc:endc, j, l),     &
               adjustment  = adjustment_one_level(begc:endc) )

          col_ps%dyn_pbal_adjustments(begc:endc) =      &
               col_ps%dyn_pbal_adjustments(begc:endc) + &
               adjustment_one_level(begc:endc) * dzsoi_decomp(j)

       end do
    end do

    do j = 1, nlevdecomp
      !col_slice = col_ps%ptrunc_vr(begc:endc,j)
       call update_column_state_no_special_handling_acc( column_state_updater, &
            bounds      = bounds,                                         &
            clump_index = clump_index,                                    &
            var         = ptrunc_vr(begc:endc,j),                &
            adjustment  = adjustment_one_level(begc:endc))

       col_ps%dyn_pbal_adjustments(begc:endc) =      &
            col_ps%dyn_pbal_adjustments(begc:endc) + &
            adjustment_one_level(begc:endc) * dzsoi_decomp(j)
       !!
       !col_slice = col_ps%solutionp_vr(begc:endc,j)
       call update_column_state_no_special_handling_acc( column_state_updater, &
            bounds      = bounds,                                         &
            clump_index = clump_index,                                    &
            var         = solutionp_vr(begc:endc,j),             &
            adjustment  = adjustment_one_level(begc:endc))

       col_ps%dyn_pbal_adjustments(begc:endc) =      &
            col_ps%dyn_pbal_adjustments(begc:endc) + &
            adjustment_one_level(begc:endc) * dzsoi_decomp(j)
        !!
       call update_column_state_no_special_handling_acc( column_state_updater, &
            bounds      = bounds,                                         &
            clump_index = clump_index,                                    &
            var         = labilep_vr(begc:endc,j),               &
            adjustment  = adjustment_one_level(begc:endc))

       col_ps%dyn_pbal_adjustments(begc:endc) =      &
            col_ps%dyn_pbal_adjustments(begc:endc) + &
            adjustment_one_level(begc:endc) * dzsoi_decomp(j)
       !!
       call update_column_state_no_special_handling_acc( column_state_updater, &
            bounds      = bounds,                                         &
            clump_index = clump_index,                                    &
            var         = secondp_vr(begc:endc,j),               &
            adjustment  = adjustment_one_level(begc:endc))

       col_ps%dyn_pbal_adjustments(begc:endc) =      &
            col_ps%dyn_pbal_adjustments(begc:endc) + &
            adjustment_one_level(begc:endc) * dzsoi_decomp(j)

    end do
    end associate

  end subroutine dyn_col_ps_Adjustments_acc


end module dynSubgridAdjustmentsMod
