module CNNStateUpdate1BeTRMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for nitrogen state variable updates, non-mortality fluxes.
  !
  ! !USES:
  use shr_kind_mod           , only: r8 => shr_kind_r8
  use clm_time_manager       , only : get_step_size
  use elm_varpar             , only : nlevdecomp, ndecomp_pools, ndecomp_cascade_transitions
  use elm_varpar             , only : crop_prog, i_met_lit, i_cel_lit, i_lig_lit, i_cwd
  use clm_varctl             , only : iulog, use_nitrif_denitrif
  use elm_varcon             , only : nitrif_n2o_loss_frac
  use pftvarcon              , only : npcropmin, nc3crop
  use VegetationPropertiesType      , only : veg_vp
  use CNDecompCascadeConType , only : decomp_cascade_con
  use CNStateType            , only : cnstate_type
  use CNNitrogenFluxType     , only : nitrogenflux_type
  use CNNitrogenStateType    , only : nitrogenstate_type
  use VegetationType         , only : veg_pp
  use VegetationDataType     , only : veg_ns, veg_nf
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
         veg_ns%leafn(p)       = veg_ns%leafn(p)       + veg_nf%leafn_xfer_to_leafn(p)*dt
         veg_ns%leafn_xfer(p)  = veg_ns%leafn_xfer(p)  - veg_nf%leafn_xfer_to_leafn(p)*dt
         veg_ns%frootn(p)      = veg_ns%frootn(p)      + veg_nf%frootn_xfer_to_frootn(p)*dt
         veg_ns%frootn_xfer(p) = veg_ns%frootn_xfer(p) - veg_nf%frootn_xfer_to_frootn(p)*dt

         if (woody(ivt(p)) == 1.0_r8) then
            veg_ns%livestemn(p)       = veg_ns%livestemn(p)       + veg_nf%livestemn_xfer_to_livestemn(p)*dt
            veg_ns%livestemn_xfer(p)  = veg_ns%livestemn_xfer(p)  - veg_nf%livestemn_xfer_to_livestemn(p)*dt
            veg_ns%deadstemn(p)       = veg_ns%deadstemn(p)       + veg_nf%deadstemn_xfer_to_deadstemn(p)*dt
            veg_ns%deadstemn_xfer(p)  = veg_ns%deadstemn_xfer(p)  - veg_nf%deadstemn_xfer_to_deadstemn(p)*dt
            veg_ns%livecrootn(p)      = veg_ns%livecrootn(p)      + veg_nf%livecrootn_xfer_to_livecrootn(p)*dt
            veg_ns%livecrootn_xfer(p) = veg_ns%livecrootn_xfer(p) - veg_nf%livecrootn_xfer_to_livecrootn(p)*dt
            veg_ns%deadcrootn(p)      = veg_ns%deadcrootn(p)      + veg_nf%deadcrootn_xfer_to_deadcrootn(p)*dt
            veg_ns%deadcrootn_xfer(p) = veg_ns%deadcrootn_xfer(p) - veg_nf%deadcrootn_xfer_to_deadcrootn(p)*dt
         end if

         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            ! lines here for consistency; the transfer terms are zero
            veg_ns%livestemn(p)       = veg_ns%livestemn(p)      + veg_nf%livestemn_xfer_to_livestemn(p)*dt
            veg_ns%livestemn_xfer(p)  = veg_ns%livestemn_xfer(p) - veg_nf%livestemn_xfer_to_livestemn(p)*dt
            veg_ns%grainn(p)          = veg_ns%grainn(p)         + veg_nf%grainn_xfer_to_grainn(p)*dt
            veg_ns%grainn_xfer(p)     = veg_ns%grainn_xfer(p)    - veg_nf%grainn_xfer_to_grainn(p)*dt
         end if

         ! phenology: litterfall and retranslocation fluxes
         veg_ns%leafn(p)    = veg_ns%leafn(p)    - veg_nf%leafn_to_litter(p)*dt
         veg_ns%frootn(p)   = veg_ns%frootn(p)   - veg_nf%frootn_to_litter(p)*dt
         veg_ns%leafn(p)    = veg_ns%leafn(p)    - veg_nf%leafn_to_retransn(p)*dt
         veg_ns%retransn(p) = veg_ns%retransn(p) + veg_nf%leafn_to_retransn(p)*dt

         ! live wood turnover and retranslocation fluxes
         if (woody(ivt(p)) == 1._r8) then
            veg_ns%livestemn(p)  = veg_ns%livestemn(p)  - veg_nf%livestemn_to_deadstemn(p)*dt
            veg_ns%deadstemn(p)  = veg_ns%deadstemn(p)  + veg_nf%livestemn_to_deadstemn(p)*dt
            veg_ns%livestemn(p)  = veg_ns%livestemn(p)  - veg_nf%livestemn_to_retransn(p)*dt
            veg_ns%retransn(p)   = veg_ns%retransn(p)   + veg_nf%livestemn_to_retransn(p)*dt
            veg_ns%livecrootn(p) = veg_ns%livecrootn(p) - veg_nf%livecrootn_to_deadcrootn(p)*dt
            veg_ns%deadcrootn(p) = veg_ns%deadcrootn(p) + veg_nf%livecrootn_to_deadcrootn(p)*dt
            veg_ns%livecrootn(p) = veg_ns%livecrootn(p) - veg_nf%livecrootn_to_retransn(p)*dt
            veg_ns%retransn(p)   = veg_ns%retransn(p)   + veg_nf%livecrootn_to_retransn(p)*dt
         end if
         if (ivt(p) >= npcropmin) then ! Beth adds retrans from froot
            veg_ns%frootn(p)     = veg_ns%frootn(p)     - veg_nf%frootn_to_retransn(p)*dt
            veg_ns%retransn(p)   = veg_ns%retransn(p)   + veg_nf%frootn_to_retransn(p)*dt
            veg_ns%livestemn(p)  = veg_ns%livestemn(p)  - veg_nf%livestemn_to_litter(p)*dt
            veg_ns%livestemn(p)  = veg_ns%livestemn(p)  - veg_nf%livestemn_to_retransn(p)*dt
            veg_ns%retransn(p)   = veg_ns%retransn(p)   + veg_nf%livestemn_to_retransn(p)*dt
            veg_ns%grainn(p)     = veg_ns%grainn(p)     - veg_nf%grainn_to_food(p)*dt
         end if

         ! uptake from soil mineral N pool
         veg_ns%npool(p) = &
              veg_ns%npool(p) + veg_nf%sminn_to_npool(p)*dt
         !write(*,*)'sminn uptake',p,veg_nf%sminn_to_npool(p)*dt
         ! deployment from retranslocation pool
         veg_ns%npool(p)    = veg_ns%npool(p)    + veg_nf%retransn_to_npool(p)*dt
         veg_ns%retransn(p) = veg_ns%retransn(p) - veg_nf%retransn_to_npool(p)*dt

         ! allocation fluxes
         veg_ns%npool(p)           = veg_ns%npool(p)          - veg_nf%npool_to_leafn(p)*dt
         veg_ns%leafn(p)           = veg_ns%leafn(p)          + veg_nf%npool_to_leafn(p)*dt
         veg_ns%npool(p)           = veg_ns%npool(p)          - veg_nf%npool_to_leafn_storage(p)*dt
         veg_ns%leafn_storage(p)   = veg_ns%leafn_storage(p)  + veg_nf%npool_to_leafn_storage(p)*dt
         veg_ns%npool(p)           = veg_ns%npool(p)          - veg_nf%npool_to_frootn(p)*dt
         veg_ns%frootn(p)          = veg_ns%frootn(p)         + veg_nf%npool_to_frootn(p)*dt
         veg_ns%npool(p)           = veg_ns%npool(p)          - veg_nf%npool_to_frootn_storage(p)*dt
         veg_ns%frootn_storage(p)  = veg_ns%frootn_storage(p) + veg_nf%npool_to_frootn_storage(p)*dt

         if (woody(ivt(p)) == 1._r8) then
            veg_ns%npool(p)              = veg_ns%npool(p)              - veg_nf%npool_to_livestemn(p)*dt
            veg_ns%livestemn(p)          = veg_ns%livestemn(p)          + veg_nf%npool_to_livestemn(p)*dt
            veg_ns%npool(p)              = veg_ns%npool(p)              - veg_nf%npool_to_livestemn_storage(p)*dt
            veg_ns%livestemn_storage(p)  = veg_ns%livestemn_storage(p)  + veg_nf%npool_to_livestemn_storage(p)*dt
            veg_ns%npool(p)              = veg_ns%npool(p)              - veg_nf%npool_to_deadstemn(p)*dt
            veg_ns%deadstemn(p)          = veg_ns%deadstemn(p)          + veg_nf%npool_to_deadstemn(p)*dt
            veg_ns%npool(p)              = veg_ns%npool(p)              - veg_nf%npool_to_deadstemn_storage(p)*dt
            veg_ns%deadstemn_storage(p)  = veg_ns%deadstemn_storage(p)  + veg_nf%npool_to_deadstemn_storage(p)*dt
            veg_ns%npool(p)              = veg_ns%npool(p)              - veg_nf%npool_to_livecrootn(p)*dt
            veg_ns%livecrootn(p)         = veg_ns%livecrootn(p)         + veg_nf%npool_to_livecrootn(p)*dt
            veg_ns%npool(p)              = veg_ns%npool(p)              - veg_nf%npool_to_livecrootn_storage(p)*dt
            veg_ns%livecrootn_storage(p) = veg_ns%livecrootn_storage(p) + veg_nf%npool_to_livecrootn_storage(p)*dt
            veg_ns%npool(p)              = veg_ns%npool(p)              - veg_nf%npool_to_deadcrootn(p)*dt
            veg_ns%deadcrootn(p)         = veg_ns%deadcrootn(p)         + veg_nf%npool_to_deadcrootn(p)*dt
            veg_ns%npool(p)              = veg_ns%npool(p)              - veg_nf%npool_to_deadcrootn_storage(p)*dt
            veg_ns%deadcrootn_storage(p) = veg_ns%deadcrootn_storage(p) + veg_nf%npool_to_deadcrootn_storage(p)*dt
         end if

         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            veg_ns%npool(p)              = veg_ns%npool(p)              - veg_nf%npool_to_livestemn(p)*dt
            veg_ns%livestemn(p)          = veg_ns%livestemn(p)          + veg_nf%npool_to_livestemn(p)*dt
            veg_ns%npool(p)              = veg_ns%npool(p)              - veg_nf%npool_to_livestemn_storage(p)*dt
            veg_ns%livestemn_storage(p)  = veg_ns%livestemn_storage(p)  + veg_nf%npool_to_livestemn_storage(p)*dt
            veg_ns%npool(p)              = veg_ns%npool(p)              - veg_nf%npool_to_grainn(p)*dt
            veg_ns%grainn(p)             = veg_ns%grainn(p)             + veg_nf%npool_to_grainn(p)*dt
            veg_ns%npool(p)              = veg_ns%npool(p)              - veg_nf%npool_to_grainn_storage(p)*dt
            veg_ns%grainn_storage(p)     = veg_ns%grainn_storage(p)     + veg_nf%npool_to_grainn_storage(p)*dt
         end if

         ! move storage pools into transfer pools
         veg_ns%leafn_storage(p)  = veg_ns%leafn_storage(p)  - veg_nf%leafn_storage_to_xfer(p)*dt
         veg_ns%leafn_xfer(p)     = veg_ns%leafn_xfer(p)     + veg_nf%leafn_storage_to_xfer(p)*dt
         veg_ns%frootn_storage(p) = veg_ns%frootn_storage(p) - veg_nf%frootn_storage_to_xfer(p)*dt
         veg_ns%frootn_xfer(p)    = veg_ns%frootn_xfer(p)    + veg_nf%frootn_storage_to_xfer(p)*dt

         if (woody(ivt(p)) == 1._r8) then
            veg_ns%livestemn_storage(p)  = veg_ns%livestemn_storage(p)  - veg_nf%livestemn_storage_to_xfer(p)*dt
            veg_ns%livestemn_xfer(p)     = veg_ns%livestemn_xfer(p)     + veg_nf%livestemn_storage_to_xfer(p)*dt
            veg_ns%deadstemn_storage(p)  = veg_ns%deadstemn_storage(p)  - veg_nf%deadstemn_storage_to_xfer(p)*dt
            veg_ns%deadstemn_xfer(p)     = veg_ns%deadstemn_xfer(p)     + veg_nf%deadstemn_storage_to_xfer(p)*dt
            veg_ns%livecrootn_storage(p) = veg_ns%livecrootn_storage(p) - veg_nf%livecrootn_storage_to_xfer(p)*dt
            veg_ns%livecrootn_xfer(p)    = veg_ns%livecrootn_xfer(p)    + veg_nf%livecrootn_storage_to_xfer(p)*dt
            veg_ns%deadcrootn_storage(p) = veg_ns%deadcrootn_storage(p) - veg_nf%deadcrootn_storage_to_xfer(p)*dt
            veg_ns%deadcrootn_xfer(p)    = veg_ns%deadcrootn_xfer(p)    + veg_nf%deadcrootn_storage_to_xfer(p)*dt
         end if

         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            ! lines here for consistency; the transfer terms are zero
            veg_ns%livestemn_storage(p)  = veg_ns%livestemn_storage(p) - veg_nf%livestemn_storage_to_xfer(p)*dt
            veg_ns%livestemn_xfer(p)     = veg_ns%livestemn_xfer(p)    + veg_nf%livestemn_storage_to_xfer(p)*dt
            veg_ns%grainn_storage(p)     = veg_ns%grainn_storage(p)    - veg_nf%grainn_storage_to_xfer(p)*dt
            veg_ns%grainn_xfer(p)        = veg_ns%grainn_xfer(p)       + veg_nf%grainn_storage_to_xfer(p)*dt
         end if

      end do

    end associate

  end subroutine NStateUpdate1

end module CNNStateUpdate1BeTRMod
