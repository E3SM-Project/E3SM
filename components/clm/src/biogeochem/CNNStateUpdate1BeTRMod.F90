module CNNStateUpdate1BeTRMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for nitrogen state variable updates, non-mortality fluxes.
  !
  ! !USES:
  use shr_kind_mod           , only: r8 => shr_kind_r8
  use clm_time_manager       , only : get_step_size
  use clm_varpar             , only : nlevdecomp, ndecomp_pools, ndecomp_cascade_transitions
  use clm_varpar             , only : crop_prog, i_met_lit, i_cel_lit, i_lig_lit, i_cwd
  use clm_varctl             , only : iulog, use_nitrif_denitrif
  use clm_varcon             , only : nitrif_n2o_loss_frac
  use pftvarcon              , only : npcropmin, nc3crop
  use VegetationPropertiesType      , only : veg_vp
  use CNDecompCascadeConType , only : decomp_cascade_con
  use CNStateType            , only : cnstate_type
  use CNNitrogenFluxType     , only : nitrogenflux_type
  use CNNitrogenStateType    , only : nitrogenstate_type
  use VegetationType                , only : veg_pp
  use tracer_varcon          , only : is_active_betr_bgc
  !! bgc interface & pflotran:
  use clm_varctl             , only : use_pflotran, pf_cmode
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: NStateUpdate1
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine NStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
       cnstate_vars, nitrogenflux_vars, nitrogenstate_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update all the prognostic nitrogen state
    ! variables (except for gap-phase mortality and fire fluxes)
    !
    ! !ARGUMENTS:
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                  , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnstate_type)       , intent(in)    :: cnstate_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    !
    ! !LOCAL VARIABLES:
    integer :: c,p,j,l,k ! indices
    integer :: fp,fc     ! lake filter indices
    real(r8):: dt        ! radiation time step (seconds)
    real(r8), parameter :: frootc_nfix_thc = 10._r8  !threshold fine root carbon for nitrogen fixation gC/m2
    !-----------------------------------------------------------------------

    associate(                                                                                           &
         ivt                   => veg_pp%itype                                , & ! Input:  [integer  (:)     ]  pft vegetation type

         woody                 => veg_vp%woody                         , & ! Input:  [real(r8) (:)     ]  binary flag for woody lifeform (1=woody, 0=not woody)

         cascade_donor_pool    => decomp_cascade_con%cascade_donor_pool    , & ! Input:  [integer  (:)     ]  which pool is C taken from for a given decomposition step
         cascade_receiver_pool => decomp_cascade_con%cascade_receiver_pool , & ! Input:  [integer  (:)     ]  which pool is C added to for a given decomposition step

         ndep_prof             => cnstate_vars%ndep_prof_col               , & ! Input:  [real(r8) (:,:)   ]  profile over which N deposition is distributed through column (1/m)
         nfixation_prof        => cnstate_vars%nfixation_prof_col          , & ! Input:  [real(r8) (:,:)   ]  profile over which N fixation is distributed through column (1/m)

         nf                    => nitrogenflux_vars                        , &
         ns                    => nitrogenstate_vars &
         )

      ! set time steps
      dt = real( get_step_size(), r8 )

      ! column-level fluxes

      ! seeding fluxes from dynamic landcover are now not accounted for in
      ! NStateUpdate1 (see CNNStateUpdate1Mod)

      ! patch loop

      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! phenology: transfer growth fluxes
         ns%leaf_patch(p)       = ns%leaf_patch(p)       + nf%leaf_xfer_to_leaf_patch(p)*dt
         ns%leaf_xfer_patch(p)  = ns%leaf_xfer_patch(p)  - nf%leaf_xfer_to_leaf_patch(p)*dt
         ns%froot_patch(p)      = ns%froot_patch(p)      + nf%froot_xfer_to_froot_patch(p)*dt
         ns%froot_xfer_patch(p) = ns%froot_xfer_patch(p) - nf%froot_xfer_to_froot_patch(p)*dt

         if (woody(ivt(p)) == 1.0_r8) then
            ns%livestem_patch(p)       = ns%livestem_patch(p)       + nf%livestem_xfer_to_livestem_patch(p)*dt
            ns%livestem_xfer_patch(p)  = ns%livestem_xfer_patch(p)  - nf%livestem_xfer_to_livestem_patch(p)*dt
            ns%deadstem_patch(p)       = ns%deadstem_patch(p)       + nf%deadstem_xfer_to_deadstem_patch(p)*dt
            ns%deadstem_xfer_patch(p)  = ns%deadstem_xfer_patch(p)  - nf%deadstem_xfer_to_deadstem_patch(p)*dt
            ns%livecroot_patch(p)      = ns%livecroot_patch(p)      + nf%livecroot_xfer_to_livecroot_patch(p)*dt
            ns%livecroot_xfer_patch(p) = ns%livecroot_xfer_patch(p) - nf%livecroot_xfer_to_livecroot_patch(p)*dt
            ns%deadcroot_patch(p)      = ns%deadcroot_patch(p)      + nf%deadcroot_xfer_to_deadcroot_patch(p)*dt
            ns%deadcroot_xfer_patch(p) = ns%deadcroot_xfer_patch(p) - nf%deadcroot_xfer_to_deadcroot_patch(p)*dt
         end if

         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            ! lines here for consistency; the transfer terms are zero
            ns%livestem_patch(p)       = ns%livestem_patch(p)      + nf%livestem_xfer_to_livestem_patch(p)*dt
            ns%livestem_xfer_patch(p)  = ns%livestem_xfer_patch(p) - nf%livestem_xfer_to_livestem_patch(p)*dt
            ns%grain_patch(p)          = ns%grain_patch(p)         + nf%grain_xfer_to_grain_patch(p)*dt
            ns%grain_xfer_patch(p)     = ns%grain_xfer_patch(p)    - nf%grain_xfer_to_grain_patch(p)*dt
         end if

         ! phenology: litterfall and retranslocation fluxes
         ns%leaf_patch(p)    = ns%leaf_patch(p)    - nf%leaf_to_litter_patch(p)*dt
         ns%froot_patch(p)   = ns%froot_patch(p)   - nf%froot_to_litter_patch(p)*dt
         ns%leaf_patch(p)    = ns%leaf_patch(p)    - nf%leafn_to_retransn_patch(p)*dt
         ns%retransn_patch(p) = ns%retransn_patch(p) + nf%leafn_to_retransn_patch(p)*dt

         ! live wood turnover and retranslocation fluxes
         if (woody(ivt(p)) == 1._r8) then
            ns%livestem_patch(p)  = ns%livestem_patch(p)  - nf%livestem_to_deadstem_patch(p)*dt
            ns%deadstem_patch(p)  = ns%deadstem_patch(p)  + nf%livestem_to_deadstem_patch(p)*dt
            ns%livestem_patch(p)  = ns%livestem_patch(p)  - nf%livestemn_to_retransn_patch(p)*dt
            ns%retransn_patch(p)   = ns%retransn_patch(p)   + nf%livestemn_to_retransn_patch(p)*dt
            ns%livecroot_patch(p) = ns%livecroot_patch(p) - nf%livecroot_to_deadcroot_patch(p)*dt
            ns%deadcroot_patch(p) = ns%deadcroot_patch(p) + nf%livecroot_to_deadcroot_patch(p)*dt
            ns%livecroot_patch(p) = ns%livecroot_patch(p) - nf%livecrootn_to_retransn_patch(p)*dt
            ns%retransn_patch(p)   = ns%retransn_patch(p)   + nf%livecrootn_to_retransn_patch(p)*dt
         end if
         if (ivt(p) >= npcropmin) then ! Beth adds retrans from froot
            ns%froot_patch(p)     = ns%froot_patch(p)     - nf%frootn_to_retransn_patch(p)*dt
            ns%retransn_patch(p)   = ns%retransn_patch(p)   + nf%frootn_to_retransn_patch(p)*dt
            ns%livestem_patch(p)  = ns%livestem_patch(p)  - nf%livestem_to_litter_patch(p)*dt
            ns%livestem_patch(p)  = ns%livestem_patch(p)  - nf%livestemn_to_retransn_patch(p)*dt
            ns%retransn_patch(p)   = ns%retransn_patch(p)   + nf%livestemn_to_retransn_patch(p)*dt
            ns%grain_patch(p)     = ns%grain_patch(p)     - nf%grain_to_food_patch(p)*dt
         end if

         ! uptake from soil mineral N pool
         ns%pool_patch(p) = &
              ns%pool_patch(p) + nf%sminn_to_npool_patch(p)*dt
         !write(*,*)'sminn uptake',p,nf%sminn_to_npool_patch(p)*dt
         ! deployment from retranslocation pool
         ns%pool_patch(p)    = ns%pool_patch(p)    + nf%retransn_to_npool_patch(p)*dt
         ns%retransn_patch(p) = ns%retransn_patch(p) - nf%retransn_to_npool_patch(p)*dt

         ! allocation fluxes
         ns%pool_patch(p)           = ns%pool_patch(p)          - nf%pool_to_leaf_patch(p)*dt
         ns%leaf_patch(p)           = ns%leaf_patch(p)          + nf%pool_to_leaf_patch(p)*dt
         ns%pool_patch(p)           = ns%pool_patch(p)          - nf%pool_to_leaf_storage_patch(p)*dt
         ns%leaf_storage_patch(p)   = ns%leaf_storage_patch(p)  + nf%pool_to_leaf_storage_patch(p)*dt
         ns%pool_patch(p)           = ns%pool_patch(p)          - nf%pool_to_froot_patch(p)*dt
         ns%froot_patch(p)          = ns%froot_patch(p)         + nf%pool_to_froot_patch(p)*dt
         ns%pool_patch(p)           = ns%pool_patch(p)          - nf%pool_to_froot_storage_patch(p)*dt
         ns%froot_storage_patch(p)  = ns%froot_storage_patch(p) + nf%pool_to_froot_storage_patch(p)*dt

         if (woody(ivt(p)) == 1._r8) then
            ns%pool_patch(p)              = ns%pool_patch(p)              - nf%pool_to_livestem_patch(p)*dt
            ns%livestem_patch(p)          = ns%livestem_patch(p)          + nf%pool_to_livestem_patch(p)*dt
            ns%pool_patch(p)              = ns%pool_patch(p)              - nf%pool_to_livestem_storage_patch(p)*dt
            ns%livestem_storage_patch(p)  = ns%livestem_storage_patch(p)  + nf%pool_to_livestem_storage_patch(p)*dt
            ns%pool_patch(p)              = ns%pool_patch(p)              - nf%pool_to_deadstem_patch(p)*dt
            ns%deadstem_patch(p)          = ns%deadstem_patch(p)          + nf%pool_to_deadstem_patch(p)*dt
            ns%pool_patch(p)              = ns%pool_patch(p)              - nf%pool_to_deadstem_storage_patch(p)*dt
            ns%deadstem_storage_patch(p)  = ns%deadstem_storage_patch(p)  + nf%pool_to_deadstem_storage_patch(p)*dt
            ns%pool_patch(p)              = ns%pool_patch(p)              - nf%pool_to_livecroot_patch(p)*dt
            ns%livecroot_patch(p)         = ns%livecroot_patch(p)         + nf%pool_to_livecroot_patch(p)*dt
            ns%pool_patch(p)              = ns%pool_patch(p)              - nf%pool_to_livecroot_storage_patch(p)*dt
            ns%livecroot_storage_patch(p) = ns%livecroot_storage_patch(p) + nf%pool_to_livecroot_storage_patch(p)*dt
            ns%pool_patch(p)              = ns%pool_patch(p)              - nf%pool_to_deadcroot_patch(p)*dt
            ns%deadcroot_patch(p)         = ns%deadcroot_patch(p)         + nf%pool_to_deadcroot_patch(p)*dt
            ns%pool_patch(p)              = ns%pool_patch(p)              - nf%pool_to_deadcroot_storage_patch(p)*dt
            ns%deadcroot_storage_patch(p) = ns%deadcroot_storage_patch(p) + nf%pool_to_deadcroot_storage_patch(p)*dt
         end if

         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            ns%pool_patch(p)              = ns%pool_patch(p)              - nf%pool_to_livestem_patch(p)*dt
            ns%livestem_patch(p)          = ns%livestem_patch(p)          + nf%pool_to_livestem_patch(p)*dt
            ns%pool_patch(p)              = ns%pool_patch(p)              - nf%pool_to_livestem_storage_patch(p)*dt
            ns%livestem_storage_patch(p)  = ns%livestem_storage_patch(p)  + nf%pool_to_livestem_storage_patch(p)*dt
            ns%pool_patch(p)              = ns%pool_patch(p)              - nf%pool_to_grain_patch(p)*dt
            ns%grain_patch(p)             = ns%grain_patch(p)             + nf%pool_to_grain_patch(p)*dt
            ns%pool_patch(p)              = ns%pool_patch(p)              - nf%pool_to_grain_storage_patch(p)*dt
            ns%grain_storage_patch(p)     = ns%grain_storage_patch(p)     + nf%pool_to_grain_storage_patch(p)*dt
         end if

         ! move storage pools into transfer pools
         ns%leaf_storage_patch(p)  = ns%leaf_storage_patch(p)  - nf%leaf_storage_to_xfer_patch(p)*dt
         ns%leaf_xfer_patch(p)     = ns%leaf_xfer_patch(p)     + nf%leaf_storage_to_xfer_patch(p)*dt
         ns%froot_storage_patch(p) = ns%froot_storage_patch(p) - nf%froot_storage_to_xfer_patch(p)*dt
         ns%froot_xfer_patch(p)    = ns%froot_xfer_patch(p)    + nf%froot_storage_to_xfer_patch(p)*dt

         if (woody(ivt(p)) == 1._r8) then
            ns%livestem_storage_patch(p)  = ns%livestem_storage_patch(p)  - nf%livestem_storage_to_xfer_patch(p)*dt
            ns%livestem_xfer_patch(p)     = ns%livestem_xfer_patch(p)     + nf%livestem_storage_to_xfer_patch(p)*dt
            ns%deadstem_storage_patch(p)  = ns%deadstem_storage_patch(p)  - nf%deadstem_storage_to_xfer_patch(p)*dt
            ns%deadstem_xfer_patch(p)     = ns%deadstem_xfer_patch(p)     + nf%deadstem_storage_to_xfer_patch(p)*dt
            ns%livecroot_storage_patch(p) = ns%livecroot_storage_patch(p) - nf%livecroot_storage_to_xfer_patch(p)*dt
            ns%livecroot_xfer_patch(p)    = ns%livecroot_xfer_patch(p)    + nf%livecroot_storage_to_xfer_patch(p)*dt
            ns%deadcroot_storage_patch(p) = ns%deadcroot_storage_patch(p) - nf%deadcroot_storage_to_xfer_patch(p)*dt
            ns%deadcroot_xfer_patch(p)    = ns%deadcroot_xfer_patch(p)    + nf%deadcroot_storage_to_xfer_patch(p)*dt
         end if

         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            ! lines here for consistency; the transfer terms are zero
            ns%livestem_storage_patch(p)  = ns%livestem_storage_patch(p) - nf%livestem_storage_to_xfer_patch(p)*dt
            ns%livestem_xfer_patch(p)     = ns%livestem_xfer_patch(p)    + nf%livestem_storage_to_xfer_patch(p)*dt
            ns%grain_storage_patch(p)     = ns%grain_storage_patch(p)    - nf%grain_storage_to_xfer_patch(p)*dt
            ns%grain_xfer_patch(p)        = ns%grain_xfer_patch(p)       + nf%grain_storage_to_xfer_patch(p)*dt
         end if

      end do

    end associate

  end subroutine NStateUpdate1

end module CNNStateUpdate1BeTRMod
