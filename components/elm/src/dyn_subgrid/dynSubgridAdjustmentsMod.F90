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
  use elm_varctl             , only : use_crop
  use elm_varcon             , only : dzsoi_decomp
  use dynPatchStateUpdaterMod, only : patch_state_updater_type
  use dynPatchStateUpdaterMod, only : update_patch_state, update_patch_state_partition_flux_by_type
  use dynColumnStateUpdaterMod , only : column_state_updater_type, update_column_state_no_special_handling
  use ColumnDataType         , only : column_carbon_state, column_nitrogen_state
  use ColumnDataType         , only : column_phosphorus_state
  use VegetationDataType     , only : vegetation_carbon_state, vegetation_nitrogen_state
  use VegetationDataType     , only : vegetation_phosphorus_state
  use SpeciesMod             , only : CN_SPECIES_N, CN_SPECIES_P
  use decompMod, only :  get_clump_bounds_gpu
  use elm_varctl , only : iulog 

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
      leafc          => veg_cs%leafc        , &
      leafc_storage  => veg_cs%leafc_storage, &
      leafc_xfer     => veg_cs%leafc_xfer   , &
      frootc         => veg_cs%frootc          ,&
      frootc_storage => veg_cs%frootc_storage  ,&
      frootc_xfer    => veg_cs%frootc_xfer    , &
      livestemc      => veg_cs%livestemc            ,&
      livestemc_storage => veg_cs%livestemc_storage  ,&
      livestemc_xfer    => veg_cs%livestemc_xfer     ,&
      deadstemc         => veg_cs%deadstemc                ,&
      deadstemc_storage => veg_cs%deadstemc_storage  ,&
      deadstemc_xfer    => veg_cs%deadstemc_xfer    , &
      livecrootc         => veg_cs%livecrootc          ,&
      livecrootc_storage => veg_cs%livecrootc_storage  ,&
      livecrootc_xfer    => veg_cs%livecrootc_xfer     ,&
      deadcrootc         => veg_cs%deadcrootc          ,&
      deadcrootc_storage => veg_cs%deadcrootc_storage  ,&
      deadcrootc_xfer    => veg_cs%deadcrootc_xfer    , &
      gresp_storage      => veg_cs%gresp_storage ,&
      gresp_xfer         => veg_cs%gresp_xfer   , &
      cpool              => veg_cs%cpool        , &
      xsmrpool           => veg_cs%xsmrpool     , &
      ctrunc             => veg_cs%ctrunc       , &
      dispvegc           => veg_cs%dispvegc     , &
      storvegc           => veg_cs%storvegc     , &
      totvegc            => veg_cs%totvegc      , &
      totpftc            => veg_cs%totpftc      , &
      grainc             => veg_cs%grainc       , &
      grainc_storage     => veg_cs%grainc_storage, &
      grainc_xfer        => veg_cs%grainc_xfer   , &
      cropseedc_deficit  => veg_cs%cropseedc_deficit &
      )

    patch_grew   = (patch_state_updater%dwt(p) > 0._r8)
    old_weight_was_zero = (patch_state_updater%pwtgcell_old(p) == 0._r8)

    call ComputeSeedAmounts(p                  , &
         species                    = veg_cs%species     , &
         leaf_patch                 = leafc(p)           , &
         leaf_storage_patch         = leafc_storage(p)   , &
         leaf_xfer_patch            = leafc_xfer(p)      , &
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
         p,c                                 , &
         var               = leafc(p)        , &
         flux_out_grc_area = conv_cflux      , &
         seed              = seed_leafc_patch      , &
         seed_addition     = dwt_leafc_seed    )

    ! 2) LEAFC_STORAGE_PATCH
    call update_patch_state(patch_state_updater ,&
         p, c                            , & 
         var               = leafc_storage(p)  , &
         flux_out_grc_area = conv_cflux     , &
         seed              = seed_leafc_storage_patch , &
         seed_addition     = dwt_leafc_seed            )

    ! 3) LEAF_XFER_PATCH
    call update_patch_state(patch_state_updater, &
         p , c                           , &
         var               = leafc_xfer(p) , &
         flux_out_grc_area = conv_cflux, &
         seed              = seed_leafc_xfer_patch , &
         seed_addition     = dwt_leafc_seed         )

    ! 4) FROOTC_PATCH
    call update_patch_state(patch_state_updater,                 &
         p ,c               , &
         var               = frootc(p)              , &
         flux_out_col_area = dwt_frootc_to_litter )

    ! 5) FROOTC_STORAGE_PATCH
    call update_patch_state(patch_state_updater,  &
         p, c                    , &
         var               = frootc_storage(p) , &
         flux_out_grc_area = conv_cflux )

    ! 6) FROOTC_XFER_PATCH
    call update_patch_state(patch_state_updater,                 &
         p,c                       , &
         var               = frootc_xfer(p)   , &
         flux_out_grc_area = conv_cflux )

    ! 7) LIVESTEMC_PATCH
    call update_patch_state(patch_state_updater,                 &
         p,c                             , &
         var               = livestemc(p)           , &
         flux_out_grc_area = conv_cflux )

    ! 8) LIVESTEMC_STORAGE_PATCH
    call update_patch_state(patch_state_updater,  &
         p, c                   , &
         var               = livestemc_storage(p)   , &
         flux_out_grc_area = conv_cflux )

    ! 9) LIVESTEMC_XFER_PATCH
    call update_patch_state(patch_state_updater,    &
         p,c                , &
         var               = livestemc_xfer(p)      , &
         flux_out_grc_area = conv_cflux )

    ! 10) PROD10_CFLUX
    wood_product_cflux       = 0._r8
    deadstemc_patch_temp     = deadstemc(p)
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
    deadstemc_patch_temp     = deadstemc(p)
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
         var                        = deadstemc(p)  , &
         flux1_out                  = conv_cflux          , &
         flux2_out                  = wood_product_cflux  , &
         seed                       = seed_deadstemc_patch, &
         seed_addition              = dwt_deadstemc_seed      )

    ! 13) DEADSTEMC_STORAGE_PATCH
    call update_patch_state(patch_state_updater,                 &
         p,c                    , &
         var               = deadstemc_storage(p)   , &
         flux_out_grc_area = conv_cflux )

    ! 14) DEADSTEMC_XFER_PATCH
     call update_patch_state(patch_state_updater,                 &
         p,c                              , &
         var               = deadstemc_xfer(p)      , &
         flux_out_grc_area = conv_cflux )

    ! 15) LIVECROOTC_PATCH
    call update_patch_state(patch_state_updater,                 &
         p,c          , &
         var               = livecrootc(p)          , &
         flux_out_col_area = dwt_livecrootc_to_litter )

    ! 16) LIVECROOTC_STORAGE_PATCH
    call update_patch_state(patch_state_updater,                 &
         p,c      , &
         var               = livecrootc_storage(p)  , &
         flux_out_grc_area = conv_cflux )

    ! 17) LIVECROOTC_XFER_PATCH
    call update_patch_state(patch_state_updater,                 &
         p,c                                   , &
         var               = livecrootc_xfer(p)     , &
         flux_out_grc_area = conv_cflux )

    ! 18) DEADCROOTC_PATCH
    call update_patch_state(patch_state_updater,                 &
         p,c                                     , &
         var               = deadcrootc(p)          , &
         flux_out_col_area = dwt_deadcrootc_to_litter )

    ! 19) DEADCROOTC_STORAGE_PATCH
    call update_patch_state(patch_state_updater,                 &
         p,c                              , &
         var               = deadcrootc_storage(p)  , &
         flux_out_grc_area = conv_cflux )

    ! 20) DEADCROOT_XFER_PATCH
    call update_patch_state(patch_state_updater,                 &
         p,c                              , &
         var               = deadcrootc_xfer(p)     , &
         flux_out_grc_area = conv_cflux )

    ! 21) GRESP_STORAGE_PATCH
    call update_patch_state(patch_state_updater,                 &
         p,c                              , &
         var               = gresp_storage(p)       , &
         flux_out_grc_area = conv_cflux )

    ! 22) GRESP_XFER_STORAGE
    call update_patch_state(patch_state_updater,                 &
         p,c                             , &
         var               = gresp_xfer(p)          , &
         flux_out_grc_area = conv_cflux )

    ! 23) CPOOL_PATCH
    call update_patch_state(patch_state_updater ,&
         p,c                            , &
         var               = cpool(p)   , &
         flux_out_grc_area = conv_cflux )

    ! 24) XSMRPOOL_PATCH
    call update_patch_state(patch_state_updater,                 &
         p,c                             , &
         var               = xsmrpool(p) , &
         flux_out_grc_area = conv_cflux )

    ! 25) CTRUNC_PATCH
    call update_patch_state(patch_state_updater,                 &
         p,c                         , &
         var               = ctrunc(p) , &
         flux_out_grc_area = conv_cflux )

    ! 26) DISPVEGC_PATCH
    call update_patch_state(patch_state_updater,                 &
         p,c                   , &
         var               = dispvegc(p) )

    ! 27) STORVEGC_PATCH
    call update_patch_state(patch_state_updater,                 &
         p,c                     , &
         var               = storvegc(p) )

    ! 28) TOTVEGC_PATCH
    call update_patch_state(patch_state_updater,                 &
         p,c                               , &
         var               = totvegc(p) )

    ! 29) TOTPFTC_PATCH
    call update_patch_state(patch_state_updater,                 &
         p,c                     , &
         var               = totpftc(p) )

    if (use_crop) then

      call update_patch_state(patch_state_updater, &
           p,c               ,&
           var = grainc(p)               , &
           flux_out_col_area = crop_product_cflux)

      call update_patch_state(patch_state_updater, &
           p,c               ,&
           var = grainc_storage(p)  , &
           flux_out_grc_area = conv_cflux)

      call update_patch_state(patch_state_updater, &
           p,c               ,&
           var = grainc_xfer(p) , &
           flux_out_grc_area = conv_cflux)
       !============================================================!

       ! This is a negative pool. So any deficit that we haven't repaid gets sucked out
       ! of the atmosphere.
       call update_patch_state(patch_state_updater,                 &
            p,c               ,&
            var = cropseedc_deficit(p)  , &
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
  subroutine dyn_col_cs_Adjustments(proc_begc,proc_endc, nclumps, column_state_updater, col_cs)
    !
    ! !DESCRIPTION:
    ! Adjust column-level carbon state variables and compute associated fluxes
    ! when patch areas change due to dynamic landuse
    !
    ! !USES:
    use decompMod, only :  get_clump_bounds_gpu
    use elm_varctl, only : iulog 
    use dynUpdateModAcc , only : update_column_state_no_special_handling_acc
    !
    ! !ARGUMENTS:
    integer                         , intent(in)   :: proc_begc, proc_endc, nclumps
    type(column_state_updater_type) , intent(in)    :: column_state_updater
    type(column_carbon_state)       , intent(inout) :: col_cs
    !
    ! !LOCAL VARIABLES:
    type(bounds_type)        :: bounds
    integer         :: l, j,c, nc 
    integer         :: begc, endc
    real(r8)        :: adjustment_one_level(proc_begc:proc_endc,1:nlevdecomp, 1:ndecomp_pools)
    real(r8)        :: sum1 
    real :: startt, stopt 
    !-----------------------------------------------------------------------
    associate(&
      decomp_cpools_vr    => col_cs%decomp_cpools_vr , &
      ctrunc_vr           => col_cs%ctrunc_vr      , &
      prod1c    => col_cs%prod1c , &
      prod10c   => col_cs%prod10c, &
      prod100c  => col_cs%prod100c &
      )

    call cpu_time(startt) 
    !$acc enter data create(adjustment_one_level(:,:,:),sum1)
    !$acc parallel loop independent gang vector default(present) 
    do c = proc_begc, proc_endc
     col_cs%dyn_cbal_adjustments(c) = 0._r8
    end do 
    call cpu_time(stopt)
    !write(iulog,*) "col_cs_Adjustment::init ",(stopt-startt)*1.E+3,"ms"

    call cpu_time(startt)
    !$acc parallel loop gang independent default(present)
    do l = 1, ndecomp_pools
      !$acc loop worker independent 
       do j = 1, nlevdecomp
          !$acc loop vector independent private(bounds, begc,endc)
          do nc =1, nclumps
            if (column_state_updater%any_changes(nc)) then
               call get_clump_bounds_gpu(nc, bounds)
               begc = bounds%begc; endc = bounds%endc 
               call update_column_state_no_special_handling_acc( column_state_updater , &
                    bounds      = bounds,                                         &
                    clump_index = nc,                                    &
                    var         = decomp_cpools_vr(begc:endc, j, l),     &
                    adjustment  = adjustment_one_level(begc:endc,j,l))
            end if 
          end do 
       end do
    end do
    call cpu_time(stopt) 
    !write(iulog,*) "col_cs_Adjustment::3D loop ",(stopt-startt)*1.E+3,"ms"

    call cpu_time(startt) 
     !$acc parallel loop independent gang worker collapse(2) default(present) private(sum1)
     do l = 1, ndecomp_pools
          do c = proc_begc, proc_endc
               sum1 = 0._r8
               !$acc loop vector reduction(+:sum1)
               do j = 1, nlevdecomp 
                    sum1 = sum1 + adjustment_one_level(c,j,l) * dzsoi_decomp(j)
               end do
               col_cs%dyn_cbal_adjustments(c) = col_cs%dyn_cbal_adjustments(c) + sum1
          end do 
     end do
     call cpu_time(stopt)
    !write(iulog,*) "col_cs_Adjustment::Reduction ",(stopt-startt)*1.E+3,"ms"

     call cpu_time(startt) 
     !$acc parallel loop independent gang default(present) present(ctrunc_vr(:,:),adjustment_one_level(:,:,:)) 
     do j = 1, nlevdecomp
          !$acc loop vector independent private(begc,endc,bounds)
          do nc =1, nclumps
            if (column_state_updater%any_changes(nc)) then
               call get_clump_bounds_gpu(nc, bounds)
               begc = bounds%begc; endc = bounds%endc 

               call update_column_state_no_special_handling_acc( column_state_updater , &
                    bounds      = bounds,                                         &
                    clump_index = nc,                                    &
                    var         = ctrunc_vr(begc:endc,j),     &
                    adjustment  = adjustment_one_level(begc:endc,j, 1))
            end if 
          end do 
    end do 
    call cpu_time(stopt) 
    !write(iulog,*) "col_cs_Adjustment::2D loop ",(stopt-startt)*1.E+3,"ms"

    call cpu_time(startt) 
     !$acc parallel loop independent gang worker default(present) private(sum1)
    do c = proc_begc, proc_endc
          sum1 = 0._r8
          !$acc loop vector reduction(+:sum1)
          do j = 1, nlevdecomp 
               sum1 = sum1 + adjustment_one_level(c,j,1) * dzsoi_decomp(j)
          end do
          col_cs%dyn_cbal_adjustments(c) = col_cs%dyn_cbal_adjustments(c) + sum1
     end do 
     call cpu_time(stopt)
    !write(iulog,*) "col_cs_Adjustment::2D reduction ",(stopt-startt)*1.E+3,"ms"

     call cpu_time(startt) 
     !$acc parallel default(present)
     !$acc loop gang vector independent private(begc,endc,bounds)
     do nc =1, nclumps
         if (column_state_updater%any_changes(nc)) then
          call get_clump_bounds_gpu(nc, bounds)
          begc = bounds%begc; endc = bounds%endc 

          call update_column_state_no_special_handling_acc(column_state_updater, &
               bounds      = bounds,                                         &
               clump_index = nc,                                    &
               var         = prod1c(begc:endc),     &
               adjustment  = adjustment_one_level(begc:endc,1,1))
         
          !$acc loop seq 
          do c = begc, endc     
               col_cs%dyn_cbal_adjustments(c) = &
                    col_cs%dyn_cbal_adjustments(c) + &
                    adjustment_one_level(c,1,1)
          end do
         end if  
     end do 

    !$acc loop gang vector independent private(begc,endc,bounds)
     do nc =1, nclumps
      if (column_state_updater%any_changes(nc)) then 
         call get_clump_bounds_gpu(nc, bounds)
          begc = bounds%begc; endc = bounds%endc 

          call update_column_state_no_special_handling_acc(column_state_updater, &
               bounds      = bounds,                                         &
               clump_index = nc,                                    &
               var         = prod10c(begc:endc),     &
               adjustment  = adjustment_one_level(begc:endc,1,1))

          !$acc loop seq 
          do c = begc, endc     
               col_cs%dyn_cbal_adjustments(c) = &
                    col_cs%dyn_cbal_adjustments(c) + &
                    adjustment_one_level(c,1,1)
          end do
      end if 
     end do 

     !$acc loop gang vector independent private(begc,endc,bounds)
     do nc =1, nclumps
      if (column_state_updater%any_changes(nc)) then
          call get_clump_bounds_gpu(nc, bounds)
          begc = bounds%begc; endc = bounds%endc 

          call update_column_state_no_special_handling_acc(column_state_updater, &
               bounds      = bounds,                                         &
               clump_index = nc,                                    &
               var         = prod100c(begc:endc),     &
               adjustment  = adjustment_one_level(begc:endc,1,1))

          !$acc loop seq 
          do c = begc, endc     
               col_cs%dyn_cbal_adjustments(c) = &
                    col_cs%dyn_cbal_adjustments(c) + &
                    adjustment_one_level(c,1,1)
          end do
         end if 
     end do
     !$acc end parallel 
     call cpu_time(stopt) 
     
     !$acc exit data delete(adjustment_one_level(:,:,:),sum1)
     !write(iulog,*) "col_cs_Adjustment::Final parallel region ",(stopt-startt)*1.E+3,"ms"
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
      leafn          => veg_ns%leafn         , &
      leafn_storage  => veg_ns%leafn_storage , &
      leafn_xfer     => veg_ns%leafn_xfer    , &
      frootn         => veg_ns%frootn        , &
      frootn_storage => veg_ns%frootn_storage, &
      frootn_xfer    => veg_ns%frootn_xfer   , &
      livestemn         => veg_ns%livestemn  , &
      livestemn_storage => veg_ns%livestemn_storage  , &
      livestemn_xfer    => veg_ns%livestemn_xfer     , &
      deadstemn         => veg_ns%deadstemn          , &
      deadstemn_storage => veg_ns%deadstemn_storage  , &
      deadstemn_xfer    => veg_ns%deadstemn_xfer     , &
      livecrootn         => veg_ns%livecrootn          , &
      livecrootn_storage => veg_ns%livecrootn_storage  , &
      livecrootn_xfer    => veg_ns%livecrootn_xfer     , &
      deadcrootn         => veg_ns%deadcrootn          , &
      deadcrootn_storage => veg_ns%deadcrootn_storage  , &
      deadcrootn_xfer    => veg_ns%deadcrootn_xfer     , &
      retransn           => veg_ns%retransn  , &
      npool              => veg_ns%npool     , &
      ntrunc             => veg_ns%ntrunc    , &
      dispvegn           => veg_ns%dispvegn  , &
      storvegn           => veg_ns%storvegn  , &
      totvegn            => veg_ns%totvegn   , &
      totpftn            => veg_ns%totpftn   , &
      grainn             => veg_ns%grainn    , &
      grainn_storage     => veg_ns%grainn_storage, &
      grainn_xfer        => veg_ns%grainn_xfer   , &
      cropseedn_deficit  => veg_ns%cropseedn_deficit &
      )

    patch_grew   = (patch_state_updater%dwt(p) > 0._r8)
    old_weight_was_zero = (patch_state_updater%pwtgcell_old(p) == 0._r8)

    call ComputeSeedAmounts(p           , &
         species                    = CN_SPECIES_N   , &
         leaf_patch                 = leafn(p)          , &
         leaf_storage_patch         = leafn_storage(p)  , &
         leaf_xfer_patch            = leafn_xfer(p)     , &
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
         var               = leafn(p)       , &
         flux_out_grc_area = conv_nflux        , &
         seed              = seed_leafn_patch  , &
         seed_addition     = dwt_leafn_seed   )

    ! 2) LEAFN_STORAGE_PATCH
    call update_patch_state(patch_state_updater,  &
         p,c                             , &
         var               = leafn_storage(p)   , &
         flux_out_grc_area = conv_nflux         , &
         seed              = seed_leafn_storage_patch    , &
         seed_addition     = dwt_leafn_seed             )

    ! 3) LEAFN_XFER_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c     , &
         var               = leafn_xfer(p)  , &
         flux_out_grc_area = conv_nflux              , &
         seed              = seed_leafn_xfer_patch   , &
         seed_addition     = dwt_leafn_seed          )

    ! 4) FROOTN_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                    , &
         var               = frootn(p)       , &
         flux_out_col_area = dwt_frootn_to_litter  )

    ! 5) FROOTN_STORAGE_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c           , &
         var               = frootn_storage(p)    , &
         flux_out_grc_area = conv_nflux    )

    ! 6) FROOTN_XFER_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                             , &
         var               = frootn_xfer(p)       , &
         flux_out_grc_area = conv_nflux      )

    ! 7) LIVESTEMN_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                                       , &
         var               = livestemn(p)         , &
         flux_out_grc_area = conv_nflux  )

    ! 8) LIVESTEMN_STORAGE_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                                   , &
         var               = livestemn_storage(p) , &
         flux_out_grc_area = conv_nflux  )

    ! 9) LIVESTEMN_XFER_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                               , &
         var               = livestemn_xfer(p)    , &
         flux_out_grc_area = conv_nflux  )

    ! 10) PROD10_NFLUX
    wood_product_nflux        = 0._r8
    deadstemn_patch_temp      = deadstemn(p)
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
    deadstemn_patch_temp      = deadstemn(p)
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
         flux1_fraction_by_pft_type = pconv             , &
         var                        = deadstemn(p)      , &
         flux1_out                  = conv_nflux           , &
         flux2_out                  = wood_product_nflux   , &
         seed                       = seed_deadstemn_patch , &
         seed_addition              = dwt_deadstemn_seed     )

    ! 13) DEADSTEMN_STORAGE_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c               , &
         var               = deadstemn_storage(p)    , &
         flux_out_grc_area = conv_nflux  )

    ! 14) DEADSTEMN_XFER_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c  , &
         var               = deadstemn_xfer(p)       , &
         flux_out_grc_area = conv_nflux  )

    ! 15) LIVECROOTN_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c              , &
         var               = livecrootn(p)           , &
         flux_out_col_area = dwt_livecrootn_to_litter  )

    ! 16) LIVECROOTN_STORAGE_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c  , &
         var               = livecrootn_storage(p)   , &
         flux_out_grc_area = conv_nflux  )

    ! 17) LIVECROOTN_XFER_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c    , &
         var               = livecrootn_xfer(p)      , &
         flux_out_grc_area = conv_nflux  )

    ! 18) DEADCROOTN_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c    , &
         var               = deadcrootn(p)           , &
         flux_out_col_area = dwt_deadcrootn_to_litter  )

    ! 19) DEADCROOTN_STORAGE_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c          , &
         var               = deadcrootn_storage(p)   , &
         flux_out_grc_area = conv_nflux  )

    ! 20) DEADCROOT_XFER_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c   , &
         var               = deadcrootn_xfer(p)     , &
         flux_out_grc_area = conv_nflux  )

    ! 21) RETRANSN_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c     , &
         var               = retransn (p)            , &
         flux_out_grc_area = conv_nflux  )

    ! 22) NTRUNC_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c    , &
         var               = ntrunc(p)               , &
         flux_out_grc_area = conv_nflux  )

    ! 23) CPOOL_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c   , &
         var               = npool(p)                , &
         flux_out_grc_area = conv_nflux  )

    ! 24) DISPVEGN_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c              , &
         var               = dispvegn(p)  )

    ! 25) STORVEGN_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                         , &
         var               = storvegn(p)  )

    ! 26) TOTVEGN_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c   , &
         var               = totvegn(p)  )

    ! 27) TOTPFTN_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c, &
         var    = totpftn(p)  )

    if (use_crop) then
      call update_patch_state(patch_state_updater, &
           p,c     , &
           var = grainn(p)     , &
           flux_out_col_area = crop_product_nflux)

      call update_patch_state(patch_state_updater, &
           p,c           , &
           var = grainn_storage(p)   , &
           flux_out_grc_area = conv_nflux)

      call update_patch_state(patch_state_updater, &
           p,c      , &
           var = grainn_xfer(p) , &
           flux_out_grc_area = conv_nflux)

       ! This is a negative pool. So any deficit that we haven't repaid gets sucked out
       ! of the atmosphere.
       call update_patch_state(patch_state_updater, &
            p,c, &
            var = cropseedn_deficit(p)   , &
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
  subroutine dyn_col_ns_Adjustments(proc_begc,proc_endc, nclumps, column_state_updater, col_ns)
    !
    ! !DESCRIPTION:
    ! Adjust state variables and compute associated fluxes when patch areas change due to
    ! dynamic landuse
    !
    ! !USES:
    use landunit_varcon          , only : istsoil, istcrop
    use dynUpdateModAcc , only : update_column_state_no_special_handling_acc
    use dynColumnStateUpdaterMod , only : column_state_updater_type
    !
    ! !ARGUMENTS:
    integer                         , intent(in)  :: proc_begc, proc_endc, nclumps
    type(column_state_updater_type) , intent(in)    :: column_state_updater
    type(column_nitrogen_state)     , intent(inout) :: col_ns
    !
    ! !LOCAL VARIABLES:
    type(bounds_type) :: bounds
    integer           :: l, j,c,nc
    integer           :: begc, endc
    real(r8)          :: adjustment_one_level(proc_begc:proc_endc,1:nlevdecomp, 1:ndecomp_pools)
    real :: startt, stopt 
    real(r8)        :: sum1 
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

      call cpu_time(startt) 
    !$acc enter data create(adjustment_one_level(:,:,:),sum1)

    !$acc parallel loop independent gang vector default(present) 
    do c = proc_begc, proc_endc
     col_ns%dyn_nbal_adjustments(c) = 0._r8
    end do 
    call cpu_time(stopt)
    !write(iulog,*) "col_ps_Adjustment::init ",(stopt-startt)*1.E+3,"ms"

    call cpu_time(startt)
    !$acc parallel loop gang independent default(present) present(adjustment_one_level(:,:,:)) 
    do l = 1, ndecomp_pools
      !$acc loop worker independent 
       do j = 1, nlevdecomp
          !$acc loop vector independent private(bounds, begc,endc)
          do nc =1, nclumps
            if (column_state_updater%any_changes(nc)) then
               call get_clump_bounds_gpu(nc, bounds)
               begc = bounds%begc; endc = bounds%endc 
               
               call update_column_state_no_special_handling_acc( column_state_updater, &
                  bounds      = bounds,                                         &
                  clump_index = nc,                                    &
                  var         = decomp_npools_vr(begc:endc, j, l),     &
                  adjustment  = adjustment_one_level(begc:endc,j,l) )
            end if 
         end do 

       end do
    end do

    call cpu_time(stopt) 
    !write(iulog,*) "col_ps_Adjustment::3D loop ",(stopt-startt)*1.E+3,"ms"

    call cpu_time(startt) 
     !$acc parallel loop independent gang worker collapse(2) default(present) private(sum1)
     do l = 1, ndecomp_pools
          do c = proc_begc, proc_endc
               sum1 = 0._r8
               !$acc loop vector reduction(+:sum1)
               do j = 1, nlevdecomp 
                    sum1 = sum1 + adjustment_one_level(c,j,l) * dzsoi_decomp(j)
               end do
               col_ns%dyn_nbal_adjustments(c) = col_ns%dyn_nbal_adjustments(c) + sum1
          end do 
     end do
     call cpu_time(stopt)
    !write(iulog,*) "col_ns_Adjustment::Reduction ",(stopt-startt)*1.E+3,"ms"

    call cpu_time(startt) 

   !$acc parallel loop independent gang default(present) present(adjustment_one_level(:,:,:))
    do j = 1, nlevdecomp
      !$acc loop vector independent private(begc,endc,bounds)
      do nc =1, nclumps
        if (column_state_updater%any_changes(nc)) then
           call get_clump_bounds_gpu(nc, bounds)
           begc = bounds%begc; endc = bounds%endc       
           
           
           call update_column_state_no_special_handling_acc(column_state_updater, &
            bounds      = bounds,                                         &
            clump_index = nc,                                    &
            var         = ntrunc_vr(begc:endc,j),     &
            adjustment  = adjustment_one_level(begc:endc,j,1))


       call update_column_state_no_special_handling_acc(column_state_updater, &
           bounds      = bounds                          , &
           clump_index = nc                     , &
           var         = sminn_vr(begc:endc, j), &
           adjustment  = adjustment_one_level(begc:endc,j,2))

       call update_column_state_no_special_handling_acc(column_state_updater, &
            bounds      = bounds                          , &
            clump_index = nc                     , &
            var         = smin_nh4_vr(begc:endc, j) , &
            adjustment  = adjustment_one_level(begc:endc,j,3))

       call update_column_state_no_special_handling_acc(column_state_updater, &
            bounds      = bounds                          , &
            clump_index = nc                     , &
            var         = smin_no3_vr(begc:endc, j) , &
            adjustment  = adjustment_one_level(begc:endc,j,4))
        end if 
      end do 
    end do
    call cpu_time(stopt)
    !write(iulog,*) "col_ns_Adjustment::2D loop ",(stopt-startt)*1.E+3,"ms"

    call cpu_time(startt) 
    !$acc parallel loop independent gang worker collapse(2) default(present) private(sum1)
    do l = 1, 4
         do c = proc_begc, proc_endc
              sum1 = 0._r8
              !$acc loop vector reduction(+:sum1)
              do j = 1, nlevdecomp 
                   sum1 = sum1 + adjustment_one_level(c,j,l) * dzsoi_decomp(j)
              end do
              col_ns%dyn_nbal_adjustments(c) = col_ns%dyn_nbal_adjustments(c) + sum1
         end do 
    end do
    call cpu_time(stopt)
   !write(iulog,*) "col_ns_Adjustment::2D vr Reduction ",(stopt-startt)*1.E+3,"ms"

   call cpu_time(startt) 
     !$acc parallel default(present) present(adjustment_one_level(:,:,:))
     !$acc loop gang vector independent private(begc,endc,bounds)
     do nc =1, nclumps
         if (column_state_updater%any_changes(nc)) then
          call get_clump_bounds_gpu(nc, bounds)
          begc = bounds%begc; endc = bounds%endc 

          call update_column_state_no_special_handling_acc(column_state_updater, &
               bounds      = bounds,                                         &
               clump_index = nc,                                    &
               var         = prod1n(begc:endc),     &
               adjustment  = adjustment_one_level(begc:endc,1,1))
         
          !$acc loop seq 
          do c = begc, endc     
               col_ns%dyn_nbal_adjustments(c) = &
                    col_ns%dyn_nbal_adjustments(c) + &
                    adjustment_one_level(c,1,1)
          end do
         end if  
     end do 

    !$acc loop gang vector independent private(begc,endc,bounds)
     do nc =1, nclumps
      if (column_state_updater%any_changes(nc)) then 
         call get_clump_bounds_gpu(nc, bounds)
          begc = bounds%begc; endc = bounds%endc 

          call update_column_state_no_special_handling_acc(column_state_updater, &
               bounds      = bounds,                                         &
               clump_index = nc,                                    &
               var         = prod10n(begc:endc),     &
               adjustment  = adjustment_one_level(begc:endc,1,1))

          !$acc loop seq 
          do c = begc, endc     
               col_ns%dyn_nbal_adjustments(c) = &
                    col_ns%dyn_nbal_adjustments(c) + &
                    adjustment_one_level(c,1,1)
          end do
      end if 
     end do 

     !$acc loop gang vector independent private(begc,endc,bounds)
     do nc =1, nclumps
      if (column_state_updater%any_changes(nc)) then
          call get_clump_bounds_gpu(nc, bounds)
          begc = bounds%begc; endc = bounds%endc 

          call update_column_state_no_special_handling_acc(column_state_updater, &
               bounds      = bounds,                                         &
               clump_index = nc,                                    &
               var         = prod100n(begc:endc),     &
               adjustment  = adjustment_one_level(begc:endc,1,1))

          !$acc loop seq 
          do c = begc, endc     
               col_ns%dyn_nbal_adjustments(c) = &
                    col_ns%dyn_nbal_adjustments(c) + &
                    adjustment_one_level(c,1,1)
          end do
      end if 
     end do
     !$acc end parallel 
     call cpu_time(stopt) 
     
     !$acc exit data delete(adjustment_one_level(:,:,:), sum1)
     !write(iulog,*) "col_ns_Adjustment::Final parallel region ",(stopt-startt)*1.E+3,"ms"
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
      leafp          => veg_ps%leafp         , &
      leafp_storage  => veg_ps%leafp_storage , &
      leafp_xfer     => veg_ps%leafp_xfer    , &
      frootp         => veg_ps%frootp        , &
      frootp_storage => veg_ps%frootp_storage, &
      frootp_xfer    => veg_ps%frootp_xfer    , &
      livestemp         => veg_ps%livestemp   , &
      livestemp_storage => veg_ps%livestemp_storage  , &
      livestemp_xfer    => veg_ps%livestemp_xfer     , &
      deadstemp         => veg_ps%deadstemp          , &
      deadstemp_storage => veg_ps%deadstemp_storage  , &
      deadstemp_xfer    => veg_ps%deadstemp_xfer     , &
      livecrootp         => veg_ps%livecrootp        , &
      livecrootp_storage => veg_ps%livecrootp_storage, &
      livecrootp_xfer    => veg_ps%livecrootp_xfer   , &
      deadcrootp         => veg_ps%deadcrootp        , &
      deadcrootp_storage => veg_ps%deadcrootp_storage ,&
      deadcrootp_xfer    => veg_ps%deadcrootp_xfer   , &
      retransp           => veg_ps%retransp  ,&
      ppool              => veg_ps%ppool    , &
      ptrunc             => veg_ps%ptrunc   , &
      dispvegp           => veg_ps%dispvegp  , &
      storvegp           => veg_ps%storvegp  , &
      totvegp            => veg_ps%totvegp   , &
      totpftp            => veg_ps%totpftp   , &
      grainp             => veg_ps%grainp    , &
      grainp_storage     => veg_ps%grainp_storage , &
      grainp_xfer        => veg_ps%grainp_xfer    , &
      cropseedp_deficit  => veg_ps%cropseedp_deficit &
      )
      patch_grew   = (patch_state_updater%dwt(p) > 0._r8)
      old_weight_was_zero = (patch_state_updater%pwtgcell_old(p) == 0._r8)

    call ComputeSeedAmounts(p           , &
         species                    = CN_SPECIES_P , &
         leaf_patch                 = leafp(p)               , &
         leaf_storage_patch         = leafp_storage(p)       , &
         leaf_xfer_patch            = leafp_xfer(p)          , &
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
         var               = leafp(p)            , &
         flux_out_grc_area = conv_pflux          , &
         seed              = seed_leafp_patch    , &
         seed_addition     = dwt_leafp_seed     )

    ! 2) LEAFP_STORAGE_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                   , &
         var               = leafp_storage(p)            , &
         flux_out_grc_area = conv_pflux                  , &
         seed              = seed_leafp_storage_patch    , &
         seed_addition     = dwt_leafp_seed             )

    ! 3) LEAFP_XFER_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                  , &
         var               = leafp_xfer(p)       , &
         flux_out_grc_area = conv_pflux              , &
         seed              = seed_leafp_xfer_patch   , &
         seed_addition     = dwt_leafp_seed          )

    ! 4) FROOTP_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                           , &
         var               = frootp(p)               , &
         flux_out_col_area = dwt_frootp_to_litter  )

    ! 5) FROOTP_STORAGE_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                             , &
         var               = frootp_storage(p)       , &
         flux_out_grc_area = conv_pflux  )

    ! 6) FROOTP_XFER_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                         , &
         var               = frootp_xfer (p)         , &
         flux_out_grc_area = conv_pflux  )

    ! 7) PROD10_PFLUX
    wood_product_pflux        = 0._r8
    deadstemp_patch_temp      = deadstemp(p)
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
    deadstemp_patch_temp      = deadstemp(p)
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
         var                        = deadstemp(p)               , &
         flux1_out                  = conv_pflux                 , &
         flux2_out                  = wood_product_pflux         , &
         seed                       = seed_deadstemp_patch       , &
         seed_addition              = dwt_deadstemp_seed     )

    ! 10) DEADSTEMP_STORAGE_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                     , &
         var               = deadstemp_storage(p)    , &
         flux_out_grc_area = conv_pflux  )

    ! 11) DEADSTEMP_XFER_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c            , &
         var               = deadstemp_xfer (p)      , &
         flux_out_grc_area = conv_pflux  )

    ! 12) LIVESTEMP_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c           , &
         var               = livestemp(p)            , &
         flux_out_grc_area = conv_pflux  )

    ! 13) LIVESTEMP_STORAGE_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                 , &
         var               = livestemp_storage(p)    , &
         flux_out_grc_area = conv_pflux  )

    ! 14) LIVESTEMP_XFER_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c       , &
         var               = livestemp_xfer(p)      , &
         flux_out_grc_area = conv_pflux  )

    ! 15) LIVECROOTP_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                 , &
         var               = livecrootp(p)           , &
         flux_out_col_area = dwt_livecrootp_to_litter  )

    ! 16) LIVECROOTP_STORAGE_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                          , &
         var               = livecrootp_storage(p)   , &
         flux_out_grc_area = conv_pflux  )

    ! 17) LIVECROOTP_XFER_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                      , &
         var               = livecrootp_xfer(p)      , &
         flux_out_grc_area = conv_pflux  )

    ! 18) DEADCROOTP_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                             , &
         var               = deadcrootp(p)           , &
         flux_out_col_area = dwt_deadcrootp_to_litter  )

    ! 19) DEADCROOTP_STORAGE_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                       , &
         var               = deadcrootp_storage(p)  , &
         flux_out_grc_area = conv_pflux  )

    ! 20) DEADCROOT_XFER_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                             , &
         var               = deadcrootp_xfer(p)     , &
         flux_out_grc_area = conv_pflux  )

    ! 21) RETRANSP_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                     , &
         var               = retransp(p)             , &
         flux_out_grc_area = conv_pflux  )

    ! 22) PTRUNC_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                           , &
         var               = ptrunc (p)              , &
         flux_out_grc_area = conv_pflux  )

    ! 23) PPOOL_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                         , &
         var               = ppool(p)                , &
         flux_out_grc_area = conv_pflux  )

    ! 24) DISPVEGP_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                      , &
         var               = dispvegp(p)  )

    ! 25) STORVEGP_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                    , &
         var               = storvegp(p)  )

    ! 26) TOTVEGP_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c             , &
         var               = totvegp(p)  )

    ! 27) TOTPFTP_PATCH
    call update_patch_state(patch_state_updater,       &
         p,c                   , &
         var               = totpftp(p)  )

    if (use_crop) then
      call update_patch_state(patch_state_updater, &
           p,c     , &
           var = grainp (p)    , &
           flux_out_col_area = crop_product_pflux)

      call update_patch_state(patch_state_updater, &
           p,c           , &
           var = grainp_storage (p)  , &
           flux_out_grc_area = conv_pflux)

      call update_patch_state(patch_state_updater, &
           p,c      , &
           var = grainp_xfer(p) , &
           flux_out_grc_area = conv_pflux)

       ! This is a negative pool. So any deficit that we haven't repaid gets sucked out
       ! of the atmosphere.
       call update_patch_state(patch_state_updater,       &
            p,c     , &
            var = cropseedp_deficit(p)   , &
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
  subroutine dyn_col_ps_Adjustments(proc_begc,proc_endc, nclumps,column_state_updater, col_ps)
    !
    ! !DESCRIPTION:
    ! Adjust state variables and compute associated fluxes when patch areas change due to
    ! dynamic landuse
    !
    ! !USES:
    use landunit_varcon          , only : istsoil, istcrop
    use dynColumnStateUpdaterMod , only : column_state_updater_type
    use dynUpdateModAcc , only : update_column_state_no_special_handling_acc

    !
    ! !ARGUMENTS:
    integer                         , intent(in)  :: proc_begc, proc_endc, nclumps
    type(column_state_updater_type) , intent(in)    :: column_state_updater
    type(column_phosphorus_state)   , intent(inout) :: col_ps
    !
    ! !LOCAL VARIABLES:
    type(bounds_type)        :: bounds
    integer           :: l, j,c,nc
    integer           :: begc, endc
    real(r8)          :: adjustment_one_level(proc_begc:proc_endc,1:nlevdecomp, 1:ndecomp_pools)
    real :: startt, stopt 
    real(r8)        :: sum1 
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

     !$acc enter data create( adjustment_one_level(:,:,:), sum1)

    !$acc parallel loop independent gang vector default(present) 
    do c = proc_begc, proc_endc
     col_ps%dyn_pbal_adjustments(c) = 0._r8
    end do 
    
    !$acc parallel loop gang independent default(present)
    do l = 1, ndecomp_pools
      !$acc loop worker independent 
       do j = 1, nlevdecomp
          !$acc loop vector independent private(bounds, begc,endc)
          do nc =1, nclumps
            if (column_state_updater%any_changes(nc)) then
               call get_clump_bounds_gpu(nc, bounds)
               begc = bounds%begc; endc = bounds%endc 
               
               call update_column_state_no_special_handling_acc( column_state_updater, &
                  bounds      = bounds,                                         &
                  clump_index = nc,                                    &
                  var         = decomp_ppools_vr(begc:endc, j, l),     &
                  adjustment  = adjustment_one_level(begc:endc,j,l) )
            end if 
         end do 

       end do
    end do

     !$acc parallel loop independent gang worker collapse(2) default(present) private(sum1)
     do l = 1, ndecomp_pools
          do c = proc_begc, proc_endc
               sum1 = 0._r8
               !$acc loop vector reduction(+:sum1)
               do j = 1, nlevdecomp 
                    sum1 = sum1 + adjustment_one_level(c,j,l) * dzsoi_decomp(j)
               end do
               col_ps%dyn_pbal_adjustments(c) = col_ps%dyn_pbal_adjustments(c) + sum1
          end do 
     end do

    !$acc parallel loop independent gang default(present)
    do j = 1, nlevdecomp
      !$acc loop vector independent private(begc,endc,bounds)
      do nc =1, nclumps
        if (column_state_updater%any_changes(nc)) then
           call get_clump_bounds_gpu(nc, bounds)
           begc = bounds%begc; endc = bounds%endc

            call update_column_state_no_special_handling_acc( column_state_updater, &
               bounds      = bounds,                                         &
               clump_index = nc,                                    &
               var         = ptrunc_vr(begc:endc,j),                &
               adjustment  = adjustment_one_level(begc:endc,j,1))

            call update_column_state_no_special_handling_acc( column_state_updater, &
                 bounds      = bounds,                                         &
                 clump_index = nc,                                    &
                 var         = solutionp_vr(begc:endc,j),             &
                 adjustment  = adjustment_one_level(begc:endc,j,2))

            call update_column_state_no_special_handling_acc( column_state_updater, &
                 bounds      = bounds,                                         &
                 clump_index = nc,                                    &
                 var         = labilep_vr(begc:endc,j),               &
                 adjustment  = adjustment_one_level(begc:endc,j,3))

            !!
            call update_column_state_no_special_handling_acc( column_state_updater, &
                 bounds      = bounds,                                         &
                 clump_index = nc,                                    &
                 var         = secondp_vr(begc:endc,j),               &
                 adjustment  = adjustment_one_level(begc:endc,j,4))


            call update_column_state_no_special_handling_acc( column_state_updater, &
                 bounds      = bounds,                                         &
                 clump_index = nc,                                    &
                 var         = occlp_vr(begc:endc,j),               &
                 adjustment  = adjustment_one_level(begc:endc,j,5))

            call update_column_state_no_special_handling_acc( column_state_updater, &
                 bounds      = bounds,                                         &
                 clump_index = nc,                                    &
                 var         = primp_vr(begc:endc,j),               &
                 adjustment  = adjustment_one_level(begc:endc,j,5))
        end if 
      end do 
    end do

    !$acc parallel loop independent gang worker collapse(2) default(present) private(sum1)
    do l = 1, 4
         do c = proc_begc, proc_endc
              sum1 = 0._r8
              !$acc loop vector reduction(+:sum1)
              do j = 1, nlevdecomp 
                   sum1 = sum1 + adjustment_one_level(c,j,l) * dzsoi_decomp(j)
              end do
              col_ps%dyn_pbal_adjustments(c) = col_ps%dyn_pbal_adjustments(c) + sum1
         end do 
    end do
     
    !$acc parallel default(present)
     !$acc loop gang vector independent private(begc,endc,bounds)
     do nc =1, nclumps
         if (column_state_updater%any_changes(nc)) then
          call get_clump_bounds_gpu(nc, bounds)
          begc = bounds%begc; endc = bounds%endc 

          call update_column_state_no_special_handling_acc(column_state_updater, &
               bounds      = bounds,                                         &
               clump_index = nc,                                    &
               var         = prod1p(begc:endc),     &
               adjustment  = adjustment_one_level(begc:endc,1,1))
         
          !$acc loop seq 
          do c = begc, endc     
               col_ps%dyn_pbal_adjustments(c) = &
                    col_ps%dyn_pbal_adjustments(c) + &
                    adjustment_one_level(c,1,1)
          end do
         end if  
     end do 

    !$acc loop gang vector independent private(begc,endc,bounds)
     do nc =1, nclumps
      if (column_state_updater%any_changes(nc)) then 
         call get_clump_bounds_gpu(nc, bounds)
          begc = bounds%begc; endc = bounds%endc 

          call update_column_state_no_special_handling_acc(column_state_updater, &
               bounds      = bounds,                                         &
               clump_index = nc,                                    &
               var         = prod10p(begc:endc),     &
               adjustment  = adjustment_one_level(begc:endc,1,1))

          !$acc loop seq 
          do c = begc, endc     
               col_ps%dyn_pbal_adjustments(c) = &
                    col_ps%dyn_pbal_adjustments(c) + &
                    adjustment_one_level(c,1,1)
          end do
      end if 
     end do 

     !$acc loop gang vector independent private(begc,endc,bounds)
     do nc =1, nclumps
      if (column_state_updater%any_changes(nc)) then
          call get_clump_bounds_gpu(nc, bounds)
          begc = bounds%begc; endc = bounds%endc 

          call update_column_state_no_special_handling_acc(column_state_updater, &
               bounds      = bounds,                                         &
               clump_index = nc,                                    &
               var         = prod100p(begc:endc),     &
               adjustment  = adjustment_one_level(begc:endc,1,1))

          !$acc loop seq 
          do c = begc, endc     
               col_ps%dyn_pbal_adjustments(c) = &
                    col_ps%dyn_pbal_adjustments(c) + &
                    adjustment_one_level(c,1,1)
          end do
      end if 
     end do
     !$acc end parallel 
     
     !$acc exit data delete(adjustment_one_level(:,:,:),sum1)


    end associate

  end subroutine dyn_col_ps_Adjustments


end module dynSubgridAdjustmentsMod
