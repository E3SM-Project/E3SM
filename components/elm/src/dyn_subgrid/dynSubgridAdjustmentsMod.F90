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
  use elm_varpar             , only : ndecomp_pools, nlevdecomp
  use elm_varctl             , only : use_crop, iulog
  use elm_varcon             , only : dzsoi_decomp
  use dynPatchStateUpdaterMod, only : patch_state_updater_type
  use dynPatchStateUpdaterMod, only : update_patch_state, update_patch_state_partition_flux_by_type
  use dynColumnStateUpdaterMod,only : column_state_updater_type, update_column_state_no_special_handling
  use ColumnDataType         , only : column_carbon_state, column_nitrogen_state
  use ColumnDataType         , only : column_phosphorus_state
  use VegetationDataType     , only : vegetation_carbon_state, vegetation_nitrogen_state
  use VegetationDataType     , only : vegetation_phosphorus_state
  use SpeciesMod             , only : CN_SPECIES_N, CN_SPECIES_P

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
       l,c,p,                       &
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
    integer           , value      , intent(in)    :: l, c, p
    type(prior_weights_type)       , intent(in)    :: prior_weights
    type(patch_state_updater_type) , intent(in)    :: patch_state_updater
    real(r8)                       , intent(inout) :: dwt_leafc_seed
    real(r8)                       , intent(inout) :: dwt_deadstemc_seed
    real(r8)                       , intent(inout) :: conv_cflux
    real(r8)                       , intent(inout) :: dwt_frootc_to_litter
    real(r8)                       , intent(inout) :: dwt_livecrootc_to_litter
    real(r8)                       , intent(inout) :: dwt_deadcrootc_to_litter
    real(r8)                       , intent(inout) :: prod10_cflux
    real(r8)                       , intent(inout) :: prod100_cflux
    real(r8)                       , intent(inout) :: crop_product_cflux
    type(vegetation_carbon_state)  , intent(inout) :: veg_cs
    !
    ! !LOCAL VARIABLES:
    integer  :: begp, endp
    logical  :: old_weight_was_zero
    logical  :: patch_grew

    ! The following are only set for growing patches:
    real(r8) :: seed_leafc_patch
    real(r8) :: seed_leafc_storage_patch
    real(r8) :: seed_leafc_xfer_patch
    real(r8) :: seed_deadstemc_patch
    real(r8) :: wood_product_cflux
    real(r8) :: deadstemc_patch_temp
    !-----------------------------------------------------------------------
    associate(&
      leafc          => veg_cs%leafc         (p), &
      leafc_storage  => veg_cs%leafc_storage (p), &
      leafc_xfer     => veg_cs%leafc_xfer    (p), &
      frootc         => veg_cs%frootc        (p)   ,&
      frootc_storage => veg_cs%frootc_storage(p)  ,&
      frootc_xfer    => veg_cs%frootc_xfer   (p)  , &
      livestemc      => veg_cs%livestemc     (p)        ,&
      livestemc_storage => veg_cs%livestemc_storage(p)  ,&
      livestemc_xfer    => veg_cs%livestemc_xfer   (p)  ,&
      deadstemc         => veg_cs%deadstemc        (p)        ,&
      deadstemc_storage => veg_cs%deadstemc_storage(p)  ,&
      deadstemc_xfer    => veg_cs%deadstemc_xfer   (p) , &
      livecrootc         => veg_cs%livecrootc      (p)    ,&
      livecrootc_storage => veg_cs%livecrootc_storage (p) ,&
      livecrootc_xfer    => veg_cs%livecrootc_xfer    (p) ,&
      deadcrootc         => veg_cs%deadcrootc         (p) ,&
      deadcrootc_storage => veg_cs%deadcrootc_storage (p) ,&
      deadcrootc_xfer    => veg_cs%deadcrootc_xfer    (p), &
      gresp_storage      => veg_cs%gresp_storage (p) ,&
      gresp_xfer         => veg_cs%gresp_xfer    (p), &
      cpool              => veg_cs%cpool         (p), &
      xsmrpool           => veg_cs%xsmrpool      (p), &
      ctrunc             => veg_cs%ctrunc        (p), &
      dispvegc           => veg_cs%dispvegc      (p), &
      storvegc           => veg_cs%storvegc      (p), &
      totvegc            => veg_cs%totvegc       (p), &
      totpftc            => veg_cs%totpftc       (p), &
      grainc             => veg_cs%grainc(p)        , &
      grainc_storage     => veg_cs%grainc_storage(p), &
      grainc_xfer        => veg_cs%grainc_xfer(p)   , &
      cropseedc_deficit  => veg_cs%cropseedc_deficit(p) &
      )

    patch_grew   = (patch_state_updater%dwt(p) > 0._r8)
    old_weight_was_zero = (patch_state_updater%pwtgcell_old(p) == 0._r8)

    call ComputeSeedAmounts(p                  , &
         species                    = veg_cs%species     , &
         leaf_patch                 = leafc           , &
         leaf_storage_patch         = leafc_storage   , &
         leaf_xfer_patch            = leafc_xfer      , &
         ! Calculations only needed for patches that grew:
         compute_here_patch         = patch_grew      , &
         ! For patches that previously had zero area, ignore the current state for the
         ! sake of computing leaf proportions:
         ignore_current_state_patch = old_weight_was_zero      , &
         seed_leaf_patch            = seed_leafc_patch         , &
         seed_leaf_storage_patch    = seed_leafc_storage_patch , &
         seed_leaf_xfer_patch       = seed_leafc_xfer_patch    , &
         seed_deadstem_patch        = seed_deadstemc_patch )

    ! 1) LEAFC_PATCH
    call update_patch_state(patch_state_updater,          &
         p,c                              , &
         var               = leafc        , &
         flux_out_grc_area = conv_cflux            , &
         seed              = seed_leafc_patch      , &
         seed_addition     = dwt_leafc_seed    )

    ! 2) LEAFC_STORAGE_PATCH
    call update_patch_state(patch_state_updater,                 &
         p, c                            , &
         var               = leafc_storage       , &
         flux_out_grc_area = conv_cflux                   , &
         seed              = seed_leafc_storage_patch     , &
         seed_addition     = dwt_leafc_seed            )

    ! 3) LEAF_XFER_PATCH
    call update_patch_state(patch_state_updater,                 &
         p , c                           , &
         var               = leafc_xfer          , &
         flux_out_grc_area = conv_cflux                   , &
         seed              = seed_leafc_xfer_patch        , &
         seed_addition     = dwt_leafc_seed         )

    ! 4) FROOTC_PATCH
    call update_patch_state(patch_state_updater,                 &
         p ,c               , &
         var               = frootc              , &
         flux_out_col_area = dwt_frootc_to_litter )

    ! 5) FROOTC_STORAGE_PATCH
    call update_patch_state(patch_state_updater,  &
         p, c                    , &
         var               = frootc_storage  , &
         flux_out_grc_area = conv_cflux )

    ! 6) FROOTC_XFER_PATCH
    call update_patch_state(patch_state_updater,                 &
         p,c                       , &
         var               = frootc_xfer   , &
         flux_out_grc_area = conv_cflux )

    ! 7) LIVESTEMC_PATCH
    call update_patch_state(patch_state_updater,                 &
         p,c                             , &
         var               = livestemc           , &
         flux_out_grc_area = conv_cflux )

    ! 8) LIVESTEMC_STORAGE_PATCH
    call update_patch_state(patch_state_updater,                 &
         p, c                   , &
         var               = livestemc_storage   , &
         flux_out_grc_area = conv_cflux )

    ! 9) LIVESTEMC_XFER_PATCH
    call update_patch_state(patch_state_updater,    &
         p,c                , &
         var               = livestemc_xfer      , &
         flux_out_grc_area = conv_cflux )

    ! 10) PROD10_CFLUX
    wood_product_cflux       = 0._r8
    deadstemc_patch_temp     = deadstemc
    call update_patch_state_partition_flux_by_type(&
         patch_state_updater    ,&
         p,c             ,&
         flux1_out_dest             = 'c'   , &
         flux1_fraction_by_pft_type = pprod10                   , &
         var                        = deadstemc_patch_temp      , &
         flux1_out                  = prod10_cflux              , &
         flux2_out                  = wood_product_cflux        , &
         seed                       = seed_deadstemc_patch      )

    ! 11) PROD100_CFLUX
    wood_product_cflux       = 0._r8
    deadstemc_patch_temp     = deadstemc
    call update_patch_state_partition_flux_by_type(&
         patch_state_updater    ,&
         p,c                   , &
         flux1_out_dest             = 'c'   , &
         flux1_fraction_by_pft_type = pprod100      , &
         var                        = deadstemc_patch_temp    , &
         flux1_out                  = prod100_cflux           , &
         flux2_out                  = wood_product_cflux      , &
         seed                       = seed_deadstemc_patch       )

    ! 12) DEADSTEMC_PATCH
    wood_product_cflux       = 0._r8
    call update_patch_state_partition_flux_by_type(&
         patch_state_updater    ,&
         p,c                   , &
         flux1_out_dest             = 'g'   , &
         flux1_fraction_by_pft_type = pconv               , &
         var                        = deadstemc  , &
         flux1_out                  = conv_cflux          , &
         flux2_out                  = wood_product_cflux  , &
         seed                       = seed_deadstemc_patch, &
         seed_addition              = dwt_deadstemc_seed      )

    ! 13) DEADSTEMC_STORAGE_PATCH
    call update_patch_state(patch_state_updater,                 &
         p,c                    , &
         var               = deadstemc_storage   , &
         flux_out_grc_area = conv_cflux )

    ! 14) DEADSTEMC_XFER_PATCH
     call update_patch_state(patch_state_updater,                 &
         p,c                              , &
         var               = deadstemc_xfer      , &
         flux_out_grc_area = conv_cflux )

    ! 15) LIVECROOTC_PATCH
    call update_patch_state(patch_state_updater,                 &
         p,c          , &
         var               = livecrootc          , &
         flux_out_col_area = dwt_livecrootc_to_litter )

    ! 16) LIVECROOTC_STORAGE_PATCH
    call update_patch_state(patch_state_updater,                 &
         p,c      , &
         var               = livecrootc_storage  , &
         flux_out_grc_area = conv_cflux )

    ! 17) LIVECROOTC_XFER_PATCH
    call update_patch_state(patch_state_updater,                 &
         p,c                                   , &
         var               = livecrootc_xfer     , &
         flux_out_grc_area = conv_cflux )

    ! 18) DEADCROOTC_PATCH
    call update_patch_state(patch_state_updater,                 &
         p,c                                     , &
         var               = deadcrootc          , &
         flux_out_col_area = dwt_deadcrootc_to_litter )

    ! 19) DEADCROOTC_STORAGE_PATCH
    call update_patch_state(patch_state_updater,                 &
         p,c                              , &
         var               = deadcrootc_storage  , &
         flux_out_grc_area = conv_cflux )

    ! 20) DEADCROOT_XFER_PATCH
    call update_patch_state(patch_state_updater,                 &
         p,c                              , &
         var               = deadcrootc_xfer     , &
         flux_out_grc_area = conv_cflux )

    ! 21) GRESP_STORAGE_PATCH
    call update_patch_state(patch_state_updater,                 &
         p,c                              , &
         var               = gresp_storage       , &
         flux_out_grc_area = conv_cflux )

    ! 22) GRESP_XFER_STORAGE
    call update_patch_state(patch_state_updater,                 &
         p,c                             , &
         var               = gresp_xfer          , &
         flux_out_grc_area = conv_cflux )

    ! 23) CPOOL_PATCH
    call update_patch_state(patch_state_updater,                 &
         p,c                            , &
         var               = cpool               , &
         flux_out_grc_area = conv_cflux )

    ! 24) XSMRPOOL_PATCH
    call update_patch_state(patch_state_updater,                 &
         p,c                              , &
         var               = xsmrpool            , &
         flux_out_grc_area = conv_cflux )

    ! 25) CTRUNC_PATCH
    call update_patch_state(patch_state_updater,                 &
         p,c                         , &
         var               = ctrunc              , &
         flux_out_grc_area = conv_cflux )

    ! 26) DISPVEGC_PATCH
    call update_patch_state(patch_state_updater,                 &
         p,c                   , &
         var               = dispvegc )

    ! 27) STORVEGC_PATCH
    call update_patch_state(patch_state_updater,                 &
         p,c                     , &
         var               = storvegc )

    ! 28) TOTVEGC_PATCH
    call update_patch_state(patch_state_updater,                 &
         p,c                               , &
         var               = totvegc )

    ! 29) TOTPFTC_PATCH
    call update_patch_state(patch_state_updater,                 &
         p,c                     , &
         var               = totpftc )

    if (use_crop) then

      call update_patch_state(patch_state_updater, &
           p,c               ,&
           var = grainc               , &
           flux_out_col_area = crop_product_cflux)

      call update_patch_state(patch_state_updater, &
           p,c               ,&
           var = grainc_storage        , &
           flux_out_grc_area = conv_cflux)

      call update_patch_state(patch_state_updater, &
           p,c               ,&
           var = grainc_xfer , &
           flux_out_grc_area = conv_cflux)
       !============================================================!

       ! This is a negative pool. So any deficit that we haven't repaid gets sucked out
       ! of the atmosphere.
       call update_patch_state(patch_state_updater,                 &
            p,c               ,&
            var = cropseedc_deficit  , &
            flux_out_grc_area = conv_cflux )

    end if

    ! These fluxes are computed as negative quantities, but are expected to be positive,
    ! so flip the signs
       dwt_frootc_to_litter     = -1._r8 * dwt_frootc_to_litter
       dwt_livecrootc_to_litter = -1._r8 * dwt_livecrootc_to_litter
       dwt_deadcrootc_to_litter = -1._r8 * dwt_deadcrootc_to_litter
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
    integer         :: l, j, c
    integer         :: begc, endc
    real(r8)        :: adjustment_one_level(bounds%begc:bounds%endc)
    !-----------------------------------------------------------------------
    associate(&
      decomp_cpools_vr    => col_cs%decomp_cpools_vr , &
      ctrunc_vr           => col_cs%ctrunc_vr      , &
      prod1c    => col_cs%prod1c , &
      prod10c   => col_cs%prod10c, &
      prod100c  => col_cs%prod100c &
      )
    begc = bounds%begc
    endc = bounds%endc

    col_cs%dyn_cbal_adjustments(begc:endc) = 0._r8

    do l = 1, ndecomp_pools
       do j = 1, nlevdecomp
          call update_column_state_no_special_handling( column_state_updater , &
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
       call update_column_state_no_special_handling( column_state_updater , &
            bounds      = bounds,                                         &
            clump_index = clump_index,                                    &
            var         = ctrunc_vr(begc:endc,j),     &
            adjustment  = adjustment_one_level(begc:endc))

       col_cs%dyn_cbal_adjustments(begc:endc) = &
            col_cs%dyn_cbal_adjustments(begc:endc) + &
            adjustment_one_level(begc:endc) * dzsoi_decomp(j)

    end do

    call update_column_state_no_special_handling( column_state_updater, &
         bounds      = bounds,                                         &
         clump_index = clump_index,                                    &
         var         = prod1c(begc:endc),     &
         adjustment  = adjustment_one_level(begc:endc))

    col_cs%dyn_cbal_adjustments(begc:endc) = &
         col_cs%dyn_cbal_adjustments(begc:endc) + &
         adjustment_one_level(begc:endc)

    call update_column_state_no_special_handling( column_state_updater, &
         bounds      = bounds,                                         &
         clump_index = clump_index,                                    &
         var         = prod10c(begc:endc),     &
         adjustment  = adjustment_one_level(begc:endc))

    col_cs%dyn_cbal_adjustments(begc:endc) = &
         col_cs%dyn_cbal_adjustments(begc:endc) + &
         adjustment_one_level(begc:endc)

    call update_column_state_no_special_handling( column_state_updater, &
         bounds      = bounds,                                         &
         clump_index = clump_index,                                    &
         var         = prod100c(begc:endc),     &
         adjustment  = adjustment_one_level(begc:endc))

    col_cs%dyn_cbal_adjustments(begc:endc) = &
         col_cs%dyn_cbal_adjustments(begc:endc) + &
         adjustment_one_level(begc:endc)

    end associate
  end subroutine dyn_col_cs_Adjustments

  !-----------------------------------------------------------------------
  subroutine dyn_veg_ns_Adjustments(  &
       l,c,p,                  &
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
    integer           , value      , intent(in)    :: l, c, p
    type(prior_weights_type)       , intent(in)    :: prior_weights
    type(patch_state_updater_type) , intent(in)    :: patch_state_updater
    real(r8)                       , intent(inout) :: dwt_leafn_seed
    real(r8)                       , intent(inout) :: dwt_deadstemn_seed
    real(r8)                       , intent(inout) :: dwt_npool_seed
    real(r8)                       , intent(inout) :: conv_nflux
    real(r8)                       , intent(inout) :: dwt_frootn_to_litter
    real(r8)                       , intent(inout) :: dwt_livecrootn_to_litter
    real(r8)                       , intent(inout) :: dwt_deadcrootn_to_litter
    real(r8)                       , intent(inout) :: prod10_nflux
    real(r8)                       , intent(inout) :: prod100_nflux
    real(r8)                       , intent(inout) :: crop_product_nflux
    type(vegetation_nitrogen_state), intent(inout) :: veg_ns
    !
    ! !LOCAL VARIABLES:
    logical    :: old_weight_was_zero
    logical    :: patch_grew
    ! The following are only set for growing patches:
    real(r8)   :: seed_leafn_patch
    real(r8)   :: seed_leafn_storage_patch
    real(r8)   :: seed_leafn_xfer_patch
    real(r8)   :: seed_deadstemn_patch
    real(r8)   :: seed_npool_patch
    real(r8)   :: wood_product_nflux
    real(r8)   :: deadstemn_patch_temp
    !-----------------------------------------------------------------------
    associate(&
      leafn          => veg_ns%leafn         (p), &
      leafn_storage  => veg_ns%leafn_storage (p), &
      leafn_xfer     => veg_ns%leafn_xfer    (p), &
      frootn         => veg_ns%frootn        (p), &
      frootn_storage => veg_ns%frootn_storage(p), &
      frootn_xfer    => veg_ns%frootn_xfer   (p), &
      livestemn         => veg_ns%livestemn  (p), &
      livestemn_storage => veg_ns%livestemn_storage (p) , &
      livestemn_xfer    => veg_ns%livestemn_xfer    (p) , &
      deadstemn         => veg_ns%deadstemn         (p) , &
      deadstemn_storage => veg_ns%deadstemn_storage (p) , &
      deadstemn_xfer    => veg_ns%deadstemn_xfer    (p) , &
      livecrootn         => veg_ns%livecrootn       (p)   , &
      livecrootn_storage => veg_ns%livecrootn_storage (p) , &
      livecrootn_xfer    => veg_ns%livecrootn_xfer    (p) , &
      deadcrootn         => veg_ns%deadcrootn         (p) , &
      deadcrootn_storage => veg_ns%deadcrootn_storage (p) , &
      deadcrootn_xfer    => veg_ns%deadcrootn_xfer    (p) , &
      retransn           => veg_ns%retransn (p) , &
      npool              => veg_ns%npool    (p) , &
      ntrunc             => veg_ns%ntrunc   (p)   , &
      dispvegn           => veg_ns%dispvegn (p) , &
      storvegn           => veg_ns%storvegn (p) , &
      totvegn            => veg_ns%totvegn  (p) , &
      totpftn            => veg_ns%totpftn  (p) , &
      grainn             => veg_ns%grainn(p)        , &
      grainn_storage     => veg_ns%grainn_storage(p), &
      grainn_xfer        => veg_ns%grainn_xfer(p)   , &
      cropseedn_deficit  => veg_ns%cropseedn_deficit(p) &
      )

    patch_grew   = (patch_state_updater%dwt(p) > 0._r8)
    old_weight_was_zero = (patch_state_updater%pwtgcell_old(p) == 0._r8)

    call ComputeSeedAmounts(p           , &
         species                    = CN_SPECIES_N   , &
         leaf_patch                 = leafn          , &
         leaf_storage_patch         = leafn_storage  , &
         leaf_xfer_patch            = leafn_xfer     , &
         ! Calculations only needed for patches that grew:
         compute_here_patch         = patch_grew     , &
         ! For patches that previously had zero area, ignore the current state for the
         ! sake of computing leaf proportions:
         ignore_current_state_patch = old_weight_was_zero      , &
         seed_leaf_patch            = seed_leafn_patch         , &
         seed_leaf_storage_patch    = seed_leafn_storage_patch , &
         seed_leaf_xfer_patch       = seed_leafn_xfer_patch    , &
         seed_deadstem_patch        = seed_deadstemn_patch     , &
         pool_seed_param            = npool_seed_param         , &
         pool_seed_patch            = seed_npool_patch )

    ! 1) LEAFN_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                                , &
         var               = leafn             , &
         flux_out_grc_area = conv_nflux        , &
         seed              = seed_leafn_patch  , &
         seed_addition     = dwt_leafn_seed   )

    ! 2) LEAFN_STORAGE_PATCH
    call update_patch_state(patch_state_updater,  &
         p,c                             , &
         var               = leafn_storage      , &
         flux_out_grc_area = conv_nflux         , &
         seed              = seed_leafn_storage_patch    , &
         seed_addition     = dwt_leafn_seed             )

    ! 3) LEAFN_XFER_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c     , &
         var               = leafn_xfer     , &
         flux_out_grc_area = conv_nflux              , &
         seed              = seed_leafn_xfer_patch   , &
         seed_addition     = dwt_leafn_seed          )

    ! 4) FROOTN_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                    , &
         var               = frootn          , &
         flux_out_col_area = dwt_frootn_to_litter  )

    ! 5) FROOTN_STORAGE_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c           , &
         var               = frootn_storage       , &
         flux_out_grc_area = conv_nflux    )

    ! 6) FROOTN_XFER_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                             , &
         var               = frootn_xfer          , &
         flux_out_grc_area = conv_nflux      )

    ! 7) LIVESTEMN_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                                       , &
         var               = livestemn            , &
         flux_out_grc_area = conv_nflux  )

    ! 8) LIVESTEMN_STORAGE_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                                   , &
         var               = livestemn_storage    , &
         flux_out_grc_area = conv_nflux  )

    ! 9) LIVESTEMN_XFER_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                               , &
         var               = livestemn_xfer       , &
         flux_out_grc_area = conv_nflux  )

    ! 10) PROD10_NFLUX
    wood_product_nflux        = 0._r8
    deadstemn_patch_temp      = deadstemn
    call update_patch_state_partition_flux_by_type(patch_state_updater , &
         p,c                  , &
         flux1_out_dest = 'c'                                , &
         flux1_fraction_by_pft_type = pprod10                  , &
         var                        = deadstemn_patch_temp       , &
         flux1_out                  = prod10_nflux               , &
         flux2_out                  = wood_product_nflux         , &
         seed                       = seed_deadstemn_patch       )

    ! 11) PROD100_NFLUX
    wood_product_nflux        = 0._r8
    deadstemn_patch_temp      = deadstemn
    call update_patch_state_partition_flux_by_type(patch_state_updater , &
         p,c              , &
         flux1_out_dest = 'c'                                , &
         flux1_fraction_by_pft_type = pprod100                 , &
         var                        = deadstemn_patch_temp       , &
         flux1_out                  = prod100_nflux              , &
         flux2_out                  = wood_product_nflux         , &
         seed                       = seed_deadstemn_patch      )

    ! 12) DEADSTEMN_PATCH
    wood_product_nflux        = 0._r8
    call update_patch_state_partition_flux_by_type(patch_state_updater , &
         p,c                         , &
         flux1_out_dest = 'g'                           , &
         flux1_fraction_by_pft_type = pconv               , &
         var                        = deadstemn           , &
         flux1_out                  = conv_nflux                 , &
         flux2_out                  = wood_product_nflux         , &
         seed                       = seed_deadstemn_patch       , &
         seed_addition              = dwt_deadstemn_seed     )

    ! 13) DEADSTEMN_STORAGE_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c               , &
         var               = deadstemn_storage    , &
         flux_out_grc_area = conv_nflux  )

    ! 14) DEADSTEMN_XFER_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c  , &
         var               = deadstemn_xfer       , &
         flux_out_grc_area = conv_nflux  )

    ! 15) LIVECROOTN_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c              , &
         var               = livecrootn           , &
         flux_out_col_area = dwt_livecrootn_to_litter  )

    ! 16) LIVECROOTN_STORAGE_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c  , &
         var               = livecrootn_storage   , &
         flux_out_grc_area = conv_nflux  )

    ! 17) LIVECROOTN_XFER_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c    , &
         var               = livecrootn_xfer      , &
         flux_out_grc_area = conv_nflux  )

    ! 18) DEADCROOTN_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c    , &
         var               = deadcrootn           , &
         flux_out_col_area = dwt_deadcrootn_to_litter  )

    ! 19) DEADCROOTN_STORAGE_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c          , &
         var               = deadcrootn_storage   , &
         flux_out_grc_area = conv_nflux  )

    ! 20) DEADCROOT_XFER_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c   , &
         var               = deadcrootn_xfer      , &
         flux_out_grc_area = conv_nflux  )

    ! 21) RETRANSN_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c     , &
         var               = retransn             , &
         flux_out_grc_area = conv_nflux  )

    ! 22) NTRUNC_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c    , &
         var               = ntrunc               , &
         flux_out_grc_area = conv_nflux  )

    ! 23) CPOOL_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c   , &
         var               = npool                , &
         flux_out_grc_area = conv_nflux  )

    ! 24) DISPVEGN_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c              , &
         var               = dispvegn  )

    ! 25) STORVEGN_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                         , &
         var               = storvegn  )

    ! 26) TOTVEGN_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c   , &
         var               = totvegn  )

    ! 27) TOTPFTN_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c, &
         var    = totpftn  )

    if (use_crop) then
      call update_patch_state(patch_state_updater, &
           p,c     , &
           var = grainn     , &
           flux_out_col_area = crop_product_nflux)

      call update_patch_state(patch_state_updater, &
           p,c           , &
           var = grainn_storage   , &
           flux_out_grc_area = conv_nflux)

      call update_patch_state(patch_state_updater, &
           p,c      , &
           var = grainn_xfer , &
           flux_out_grc_area = conv_nflux)

       ! This is a negative pool. So any deficit that we haven't repaid gets sucked out
       ! of the atmosphere.
       call update_patch_state(patch_state_updater,       &
            p,c, &
            var = cropseedn_deficit   , &
            flux_out_grc_area = conv_nflux  )
    end if

    ! These fluxes are computed as negative quantities, but are expected to be positive,
    ! so flip the signs
    dwt_frootn_to_litter      = -1._r8 * dwt_frootn_to_litter
    dwt_livecrootn_to_litter  = -1._r8 * dwt_livecrootn_to_litter
    dwt_deadcrootn_to_litter  = -1._r8 * dwt_deadcrootn_to_litter

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
    integer                     :: l, j, c
    integer                     :: begc, endc
    real(r8)                    :: adjustment_one_level(bounds%begc:bounds%endc)
    !-----------------------------------------------------------------------
    associate(&
      decomp_npools_vr    => col_ns%decomp_npools_vr , &
      ntrunc_vr           => col_ns%ntrunc_vr      , &
      sminn_vr            => col_ns%sminn_vr ,  &
      smin_nh4_vr         => col_ns%smin_nh4_vr  , &
      smin_no3_vr         => col_ns%smin_no3_vr  , &
      prod1n              => col_ns%prod1n       , &
      prod10n             => col_ns%prod10n      , &
      prod100n            => col_ns%prod100n     &
      )

    begc = bounds%begc
    endc = bounds%endc
    col_ns%dyn_nbal_adjustments(begc:endc) = 0._r8
    do l = 1, ndecomp_pools
       do j = 1, nlevdecomp
          call update_column_state_no_special_handling(column_state_updater, &
               bounds      = bounds,                                         &
               clump_index = clump_index,                                    &
               var         = decomp_npools_vr(begc:endc, j, l),     &
               adjustment  = adjustment_one_level(begc:endc))
          ! TRS Check the first element to see if it is NaN - this is
          ! causing problems for us while debugging.  I am not sure if
          ! col_ns%dyn_nbal_adjustments being NaN is an actual problem
          ! or not, but it bombs out while debugging so we have to trap
          ! for it
          if (.not. isnan(col_ns%dyn_nbal_adjustments(begc))) then
             col_ns%dyn_nbal_adjustments(begc:endc) = &
                  col_ns%dyn_nbal_adjustments(begc:endc) + &
                  adjustment_one_level(begc:endc) * dzsoi_decomp(j)
          endif

       end do
    end do

    do j = 1, nlevdecomp
       call update_column_state_no_special_handling(column_state_updater, &
            bounds      = bounds,                                         &
            clump_index = clump_index,                                    &
            var         = ntrunc_vr(begc:endc,j),     &
            adjustment  = adjustment_one_level(begc:endc))

       ! TRS NaN trap, for debugging
       if (.not. isnan(col_ns%dyn_nbal_adjustments(begc))) then
          col_ns%dyn_nbal_adjustments(begc:endc) = &
               col_ns%dyn_nbal_adjustments(begc:endc) + &
               adjustment_one_level(begc:endc) * dzsoi_decomp(j)
       endif


       call update_column_state_no_special_handling(column_state_updater, &
           bounds      = bounds                          , &
           clump_index = clump_index                     , &
           var         = sminn_vr(begc:endc, j), &
           adjustment  = adjustment_one_level(begc:endc))

       ! TRS NaN trap, for debugging
       if (.not. isnan(col_ns%dyn_nbal_adjustments(begc))) then 
          col_ns%dyn_nbal_adjustments(begc:endc) = &
               col_ns%dyn_nbal_adjustments(begc:endc) + &
               adjustment_one_level(begc:endc) * dzsoi_decomp(j)
       endif
       call update_column_state_no_special_handling(column_state_updater, &
            bounds      = bounds                          , &
            clump_index = clump_index                     , &
            var         = smin_nh4_vr(begc:endc, j) , &
            adjustment  = adjustment_one_level(begc:endc))

       call update_column_state_no_special_handling(column_state_updater, &
            bounds      = bounds                          , &
            clump_index = clump_index                     , &
            var         = smin_no3_vr(begc:endc, j) , &
            adjustment  = adjustment_one_level(begc:endc))
    end do
    call update_column_state_no_special_handling(column_state_updater, &
         bounds      = bounds,                                         &
         clump_index = clump_index,                                    &
         var         = prod1n(begc:endc),     &
         adjustment  = adjustment_one_level(begc:endc))

    if (.not. isnan(col_ns%dyn_nbal_adjustments(begc))) then
       col_ns%dyn_nbal_adjustments(begc:endc) = &
            col_ns%dyn_nbal_adjustments(begc:endc) + &
            adjustment_one_level(begc:endc)
    endif

    call update_column_state_no_special_handling(column_state_updater, &
         bounds      = bounds,                                         &
         clump_index = clump_index,                                    &
         var         = prod10n(begc:endc),     &
         adjustment  = adjustment_one_level(begc:endc))

    if (.not. isnan(col_ns%dyn_nbal_adjustments(begc))) then
       col_ns%dyn_nbal_adjustments(begc:endc) = &
            col_ns%dyn_nbal_adjustments(begc:endc) + &
            adjustment_one_level(begc:endc)
    endif

    call update_column_state_no_special_handling(column_state_updater, &
         bounds      = bounds,                                         &
         clump_index = clump_index,                                    &
         var         = prod100n(begc:endc),     &
         adjustment  = adjustment_one_level(begc:endc))

    ! TRS NaN trap, for debugging
    if (.not. isnan(col_ns%dyn_nbal_adjustments(begc))) then
       col_ns%dyn_nbal_adjustments(begc:endc) = &
            col_ns%dyn_nbal_adjustments(begc:endc) + &
            adjustment_one_level(begc:endc)
    endif
    !=======================================================!
    end associate

  end subroutine dyn_col_ns_Adjustments

  !-----------------------------------------------------------------------
  subroutine dyn_veg_ps_Adjustments(        &
       l,c,p,                              &
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
    use ComputeSeedMod   , only : ComputeSeedAmounts
    !
    ! !ARGUMENTS:
    integer,  value                  , intent(in)    :: l,c,p
    type(prior_weights_type)         , intent(in)    :: prior_weights
    type(patch_state_updater_type)   , intent(in)    :: patch_state_updater
    real(r8)                         , intent(inout) :: dwt_leafp_seed
    real(r8)                         , intent(inout) :: dwt_deadstemp_seed
    real(r8)                         , intent(inout) :: dwt_ppool_seed
    real(r8)                         , intent(inout) :: conv_pflux
    real(r8)                         , intent(inout) :: dwt_frootp_to_litter
    real(r8)                         , intent(inout) :: dwt_livecrootp_to_litter
    real(r8)                         , intent(inout) :: dwt_deadcrootp_to_litter
    real(r8)                         , intent(inout) :: prod10_pflux
    real(r8)                         , intent(inout) :: prod100_pflux
    real(r8)                         , intent(inout) :: crop_product_pflux
    type(vegetation_phosphorus_state), intent(inout) :: veg_ps
    !
    ! !LOCAL VARIABLES:
    logical   :: old_weight_was_zero
    logical   :: patch_grew

    ! The following are only set for growing patches:
    real(r8)  :: seed_leafp_patch
    real(r8)  :: seed_leafp_storage_patch
    real(r8)  :: seed_leafp_xfer_patch
    real(r8)  :: seed_deadstemp_patch
    real(r8)  :: seed_ppool_patch

    real(r8)  :: wood_product_pflux
    real(r8)  :: deadstemp_patch_temp
    !-----------------------------------------------------------------------
    associate(&
      leafp          => veg_ps%leafp         (p), &
      leafp_storage  => veg_ps%leafp_storage (p), &
      leafp_xfer     => veg_ps%leafp_xfer    (p), &
      frootp         => veg_ps%frootp        (p)        ,&
      frootp_storage => veg_ps%frootp_storage(p)  ,&
      frootp_xfer    => veg_ps%frootp_xfer   (p) , &
      livestemp         => veg_ps%livestemp  (p)              ,&
      livestemp_storage => veg_ps%livestemp_storage  (p) ,&
      livestemp_xfer    => veg_ps%livestemp_xfer     (p), &
      deadstemp         => veg_ps%deadstemp          (p)       ,&
      deadstemp_storage => veg_ps%deadstemp_storage  (p) ,&
      deadstemp_xfer    => veg_ps%deadstemp_xfer     (p), &
      livecrootp         => veg_ps%livecrootp        (p)   ,&
      livecrootp_storage => veg_ps%livecrootp_storage(p) ,&
      livecrootp_xfer    => veg_ps%livecrootp_xfer   (p), &
      deadcrootp         => veg_ps%deadcrootp        (p)       ,&
      deadcrootp_storage => veg_ps%deadcrootp_storage(p) ,&
      deadcrootp_xfer    => veg_ps%deadcrootp_xfer   (p), &
      retransp           => veg_ps%retransp (p) ,&
      ppool              => veg_ps%ppool    (p), &
      ptrunc             => veg_ps%ptrunc   (p)   , &
      dispvegp           => veg_ps%dispvegp (p) , &
      storvegp           => veg_ps%storvegp (p) , &
      totvegp            => veg_ps%totvegp  (p) , &
      totpftp            => veg_ps%totpftp  (p) , &
      grainp             => veg_ps%grainp(p)        , &
      grainp_storage     => veg_ps%grainp_storage(p), &
      grainp_xfer        => veg_ps%grainp_xfer(p)   , &
      cropseedp_deficit  => veg_ps%cropseedp_deficit(p) &
      )
      patch_grew   = (patch_state_updater%dwt(p) > 0._r8)
      old_weight_was_zero = (patch_state_updater%pwtgcell_old(p) == 0._r8)

    call ComputeSeedAmounts(p           , &
         species                    = CN_SPECIES_P , &
         leaf_patch                 = leafp               , &
         leaf_storage_patch         = leafp_storage       , &
         leaf_xfer_patch            = leafp_xfer          , &
         ! Calculations only needed for patches that grew:
         compute_here_patch         = patch_grew                 , &
         ! For patches that previously had zero area, ignore the current state for the
         ! sake of computing leaf proportions:
         ignore_current_state_patch = old_weight_was_zero        , &
         seed_leaf_patch            = seed_leafp_patch           , &
         seed_leaf_storage_patch    = seed_leafp_storage_patch   , &
         seed_leaf_xfer_patch       = seed_leafp_xfer_patch      , &
         seed_deadstem_patch        = seed_deadstemp_patch       , &
         pool_seed_param            = ppool_seed_param           , &
         pool_seed_patch            = seed_ppool_patch  )

    ! 1) LEAFP_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                           , &
         var               = leafp               , &
         flux_out_grc_area = conv_pflux          , &
         seed              = seed_leafp_patch    , &
         seed_addition     = dwt_leafp_seed     )

    ! 2) LEAFP_STORAGE_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                   , &
         var               = leafp_storage               , &
         flux_out_grc_area = conv_pflux                  , &
         seed              = seed_leafp_storage_patch    , &
         seed_addition     = dwt_leafp_seed             )

    ! 3) LEAFP_XFER_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                  , &
         var               = leafp_xfer       , &
         flux_out_grc_area = conv_pflux              , &
         seed              = seed_leafp_xfer_patch   , &
         seed_addition     = dwt_leafp_seed          )

    ! 4) FROOTP_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                           , &
         var               = frootp               , &
         flux_out_col_area = dwt_frootp_to_litter  )

    ! 5) FROOTP_STORAGE_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                             , &
         var               = frootp_storage       , &
         flux_out_grc_area = conv_pflux  )

    ! 6) FROOTP_XFER_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                         , &
         var               = frootp_xfer          , &
         flux_out_grc_area = conv_pflux  )

    ! 7) PROD10_PFLUX
    wood_product_pflux        = 0._r8
    deadstemp_patch_temp      = deadstemp
    call update_patch_state_partition_flux_by_type(patch_state_updater, &
         p,c                            , &
         flux1_out_dest = 'c'                                , &
         flux1_fraction_by_pft_type = pprod10                    , &
         var                        = deadstemp_patch_temp       , &
         flux1_out                  = prod10_pflux               , &
         flux2_out                  = wood_product_pflux         , &
         seed                       = seed_deadstemp_patch       )

    ! 8) PROD100_PFLUX
    wood_product_pflux        = 0._r8
    deadstemp_patch_temp      = deadstemp
    call update_patch_state_partition_flux_by_type(patch_state_updater, &
         p,c                          , &
         flux1_out_dest = 'c'                                , &
         flux1_fraction_by_pft_type = pprod100                   , &
         var                        = deadstemp_patch_temp       , &
         flux1_out                  = prod100_pflux              , &
         flux2_out                  = wood_product_pflux         , &
         seed                       = seed_deadstemp_patch      )

    ! 9) DEADSTEMP_PATCH
    wood_product_pflux        = 0._r8
    call update_patch_state_partition_flux_by_type(patch_state_updater, &
         p,c       , &
         flux1_out_dest = 'g'                                , &
         flux1_fraction_by_pft_type = pconv                   , &
         var                        = deadstemp               , &
         flux1_out                  = conv_pflux                 , &
         flux2_out                  = wood_product_pflux         , &
         seed                       = seed_deadstemp_patch       , &
         seed_addition              = dwt_deadstemp_seed     )

    ! 10) DEADSTEMP_STORAGE_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                     , &
         var               = deadstemp_storage    , &
         flux_out_grc_area = conv_pflux  )

    ! 11) DEADSTEMP_XFER_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c            , &
         var               = deadstemp_xfer       , &
         flux_out_grc_area = conv_pflux  )

    ! 12) LIVESTEMP_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c           , &
         var               = livestemp            , &
         flux_out_grc_area = conv_pflux  )

    ! 13) LIVESTEMP_STORAGE_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                 , &
         var               = livestemp_storage    , &
         flux_out_grc_area = conv_pflux  )

    ! 14) LIVESTEMP_XFER_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c       , &
         var               = livestemp_xfer       , &
         flux_out_grc_area = conv_pflux  )

    ! 15) LIVECROOTP_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                 , &
         var               = livecrootp           , &
         flux_out_col_area = dwt_livecrootp_to_litter  )

    ! 16) LIVECROOTP_STORAGE_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                          , &
         var               = livecrootp_storage   , &
         flux_out_grc_area = conv_pflux  )

    ! 17) LIVECROOTP_XFER_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                      , &
         var               = livecrootp_xfer      , &
         flux_out_grc_area = conv_pflux  )

    ! 18) DEADCROOTP_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                             , &
         var               = deadcrootp           , &
         flux_out_col_area = dwt_deadcrootp_to_litter  )

    ! 19) DEADCROOTP_STORAGE_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                       , &
         var               = deadcrootp_storage   , &
         flux_out_grc_area = conv_pflux  )

    ! 20) DEADCROOT_XFER_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                             , &
         var               = deadcrootp_xfer      , &
         flux_out_grc_area = conv_pflux  )

    ! 21) RETRANSP_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                     , &
         var               = retransp             , &
         flux_out_grc_area = conv_pflux  )

    ! 22) PTRUNC_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                           , &
         var               = ptrunc               , &
         flux_out_grc_area = conv_pflux  )

    ! 23) PPOOL_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                         , &
         var               = ppool                , &
         flux_out_grc_area = conv_pflux  )

    ! 24) DISPVEGP_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                      , &
         var               = dispvegp  )

    ! 25) STORVEGP_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                    , &
         var               = storvegp  )

    ! 26) TOTVEGP_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c             , &
         var               = totvegp  )

    ! 27) TOTPFTP_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                   , &
         var               = totpftp  )

    if (use_crop) then
      call update_patch_state(patch_state_updater, &
           p,c     , &
           var = grainp     , &
           flux_out_col_area = crop_product_pflux)

      call update_patch_state(patch_state_updater, &
           p,c           , &
           var = grainp_storage   , &
           flux_out_grc_area = conv_pflux)

      call update_patch_state(patch_state_updater, &
           p,c      , &
           var = grainp_xfer , &
           flux_out_grc_area = conv_pflux)

       ! This is a negative pool. So any deficit that we haven't repaid gets sucked out
       ! of the atmosphere.
       call update_patch_state(patch_state_updater,       &
            p,c     , &
            var = cropseedp_deficit   , &
            flux_out_grc_area = conv_pflux  )
    end if

    ! These fluxes are computed as negative quantities, but are expected to be positive,
    ! so flip the signs
    dwt_frootp_to_litter     = -1._r8 * dwt_frootp_to_litter
    dwt_livecrootp_to_litter = -1._r8 * dwt_livecrootp_to_litter
    dwt_deadcrootp_to_litter = -1._r8 * dwt_deadcrootp_to_litter
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
      secondp_vr          => col_ps%secondp_vr , &
      occlp_vr            => col_ps%occlp_vr , &
      primp_vr            => col_ps%primp_vr , &
      prod1p              => col_ps%prod1p  , &
      prod10p             => col_ps%prod10p , &
      prod100p            => col_ps%prod100p &
      )

    begc = bounds%begc
    endc = bounds%endc

    col_ps%dyn_pbal_adjustments(begc:endc) = 0._r8
    do l = 1, ndecomp_pools
       do j = 1, nlevdecomp

          call update_column_state_no_special_handling( column_state_updater, &
               bounds      = bounds,                                         &
               clump_index = clump_index,                                    &
               var         = decomp_ppools_vr(begc:endc, j, l),     &
               adjustment  = adjustment_one_level(begc:endc) )

          if (.not. isnan(col_ps%dyn_pbal_adjustments(begc))) then
             col_ps%dyn_pbal_adjustments(begc:endc) =      &
                  col_ps%dyn_pbal_adjustments(begc:endc) + &
                  adjustment_one_level(begc:endc) * dzsoi_decomp(j)
          endif

       end do
    end do

    do j = 1, nlevdecomp
       call update_column_state_no_special_handling( column_state_updater, &
            bounds      = bounds,                                         &
            clump_index = clump_index,                                    &
            var         = ptrunc_vr(begc:endc,j),                &
            adjustment  = adjustment_one_level(begc:endc))

       if (.not. isnan(col_ps%dyn_pbal_adjustments(begc))) then
          col_ps%dyn_pbal_adjustments(begc:endc) =      &
               col_ps%dyn_pbal_adjustments(begc:endc) + &
               adjustment_one_level(begc:endc) * dzsoi_decomp(j)
       endif
       call update_column_state_no_special_handling( column_state_updater, &
            bounds      = bounds,                                         &
            clump_index = clump_index,                                    &
            var         = solutionp_vr(begc:endc,j),             &
            adjustment  = adjustment_one_level(begc:endc))

       if (.not. isnan(col_ps%dyn_pbal_adjustments(begc))) then
          col_ps%dyn_pbal_adjustments(begc:endc) =      &
               col_ps%dyn_pbal_adjustments(begc:endc) + &
               adjustment_one_level(begc:endc) * dzsoi_decomp(j)
       endif

       call update_column_state_no_special_handling( column_state_updater, &
            bounds      = bounds,                                         &
            clump_index = clump_index,                                    &
            var         = labilep_vr(begc:endc,j),               &
            adjustment  = adjustment_one_level(begc:endc))

       if (.not. isnan(col_ps%dyn_pbal_adjustments(begc))) then
          col_ps%dyn_pbal_adjustments(begc:endc) =      &
               col_ps%dyn_pbal_adjustments(begc:endc) + &
               adjustment_one_level(begc:endc) * dzsoi_decomp(j)
       endif
       !!
       call update_column_state_no_special_handling( column_state_updater, &
            bounds      = bounds,                                         &
            clump_index = clump_index,                                    &
            var         = secondp_vr(begc:endc,j),               &
            adjustment  = adjustment_one_level(begc:endc))

       col_ps%dyn_pbal_adjustments(begc:endc) =      &
            col_ps%dyn_pbal_adjustments(begc:endc) + &
            adjustment_one_level(begc:endc) * dzsoi_decomp(j)

       call update_column_state_no_special_handling( column_state_updater, &
            bounds      = bounds,                                         &
            clump_index = clump_index,                                    &
            var         = occlp_vr(begc:endc,j),               &
            adjustment  = adjustment_one_level(begc:endc))

       call update_column_state_no_special_handling( column_state_updater, &
            bounds      = bounds,                                         &
            clump_index = clump_index,                                    &
            var         = primp_vr(begc:endc,j),               &
            adjustment  = adjustment_one_level(begc:endc))

    end do
    call update_column_state_no_special_handling( column_state_updater, &
         bounds      = bounds,                                         &
         clump_index = clump_index,                                    &
         var         = prod1p(begc:endc),     &
         adjustment  = adjustment_one_level(begc:endc))

    if (.not. isnan(col_ps%dyn_pbal_adjustments(begc))) then
       col_ps%dyn_pbal_adjustments(begc:endc) = &
            col_ps%dyn_pbal_adjustments(begc:endc) + &
            adjustment_one_level(begc:endc)
    endif

    call update_column_state_no_special_handling( column_state_updater, &
         bounds      = bounds,                                         &
         clump_index = clump_index,                                    &
         var         = prod10p(begc:endc),     &
         adjustment  = adjustment_one_level(begc:endc))

    if (.not. isnan(col_ps%dyn_pbal_adjustments(begc))) then
       col_ps%dyn_pbal_adjustments(begc:endc) = &
            col_ps%dyn_pbal_adjustments(begc:endc) + &
            adjustment_one_level(begc:endc)
    endif

    call update_column_state_no_special_handling( column_state_updater, &
         bounds      = bounds,                                         &
         clump_index = clump_index,                                    &
         var         = prod100p(begc:endc),     &
         adjustment  = adjustment_one_level(begc:endc))

    if (.not. isnan(col_ps%dyn_pbal_adjustments(begc))) then
       col_ps%dyn_pbal_adjustments(begc:endc) = &
            col_ps%dyn_pbal_adjustments(begc:endc) + &
            adjustment_one_level(begc:endc)
    endif

    end associate

  end subroutine dyn_col_ps_Adjustments


end module dynSubgridAdjustmentsMod
