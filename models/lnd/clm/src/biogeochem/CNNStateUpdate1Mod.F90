module CNNStateUpdate1Mod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for nitrogen state variable updates, non-mortality fluxes.
  !
  ! !USES:
  use shr_kind_mod                    , only: r8 => shr_kind_r8
  use clm_time_manager                , only : get_step_size
  use clm_varpar                      , only : nlevdecomp, ndecomp_pools, ndecomp_cascade_transitions
  use clm_varpar                      , only : crop_prog, i_met_lit, i_cel_lit, i_lig_lit, i_cwd
  use clm_varctl                      , only : iulog, use_nitrif_denitrif
  use clm_varcon                      , only : nitrif_n2o_loss_frac
  use pftconMod                       , only : npcropmin, pftcon
  use CNVegNitrogenStateType          , only : cnveg_nitrogenstate_type
  use CNVegNitrogenFluxType           , only : cnveg_nitrogenflux_type
  use SoilBiogeochemNitrogenFluxType  , only : soilbiogeochem_nitrogenflux_type
  use PatchType                       , only : patch                
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: NStateUpdate1
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine NStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
       cnveg_nitrogenflux_inst, cnveg_nitrogenstate_inst, soilbiogeochem_nitrogenflux_inst) 
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update all the prognostic nitrogen state
    ! variables (except for gap-phase mortality and fire fluxes)
    !
    ! !ARGUMENTS:
    integer                                 , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                                 , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                                 , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                                 , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnveg_nitrogenflux_type)           , intent(in)    :: cnveg_nitrogenflux_inst
    type(cnveg_nitrogenstate_type)          , intent(inout) :: cnveg_nitrogenstate_inst
    type(soilbiogeochem_nitrogenflux_type)  , intent(inout) :: soilbiogeochem_nitrogenflux_inst
    !
    ! !LOCAL VARIABLES:
    integer :: c,p,j,l,k ! indices
    integer :: fp,fc     ! lake filter indices
    real(r8):: dt        ! radiation time step (seconds)
    !-----------------------------------------------------------------------

    associate(                                                                   & 
         ivt                   => patch%itype                                    , & ! Input:  [integer  (:)     ]  patch vegetation type                                

         woody                 => pftcon%woody                                 , & ! Input:  binary flag for woody lifeform (1=woody, 0=not woody)

         nf_veg                => cnveg_nitrogenflux_inst                      , & ! Input:
         ns_veg                => cnveg_nitrogenstate_inst                     , & ! Output:
         nf_soil               => soilbiogeochem_nitrogenflux_inst               & ! Output:
         )

      ! set time steps
      dt = real( get_step_size(), r8 )


      ! soilbiogeochemistry fluxes TODO - this should be moved elsewhere
      ! plant to litter fluxes -  phenology and dynamic landcover fluxes
      do j = 1, nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)

            nf_soil%decomp_npools_sourcesink_col(c,j,i_met_lit) = &
                 ( nf_veg%phenology_n_to_litr_met_n_col(c,j) + nf_veg%dwt_frootn_to_litr_met_n_col(c,j) ) * dt

            nf_soil%decomp_npools_sourcesink_col(c,j,i_cel_lit) = &
                 ( nf_veg%phenology_n_to_litr_cel_n_col(c,j) + nf_veg%dwt_frootn_to_litr_cel_n_col(c,j) ) * dt

            nf_soil%decomp_npools_sourcesink_col(c,j,i_lig_lit) = &
                 ( nf_veg%phenology_n_to_litr_lig_n_col(c,j) + nf_veg%dwt_frootn_to_litr_lig_n_col(c,j) ) * dt

            nf_soil%decomp_npools_sourcesink_col(c,j,i_cwd)     = &
                 ( nf_veg%dwt_livecrootn_to_cwdn_col(c,j)    + nf_veg%dwt_deadcrootn_to_cwdn_col(c,j) )   * dt

         end do
      end do

      ! seeding fluxes, from dynamic landcover
      do fc = 1,num_soilc
         c = filter_soilc(fc)

         ns_veg%seedn_col(c) = ns_veg%seedn_col(c) - nf_veg%dwt_seedn_to_leaf_col(c) * dt
         ns_veg%seedn_col(c) = ns_veg%seedn_col(c) - nf_veg%dwt_seedn_to_deadstem_col(c) * dt
      end do

      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! phenology: transfer growth fluxes
         ns_veg%leafn_patch(p)       = ns_veg%leafn_patch(p)       + nf_veg%leafn_xfer_to_leafn_patch(p)*dt
         ns_veg%leafn_xfer_patch(p)  = ns_veg%leafn_xfer_patch(p)  - nf_veg%leafn_xfer_to_leafn_patch(p)*dt
         ns_veg%frootn_patch(p)      = ns_veg%frootn_patch(p)      + nf_veg%frootn_xfer_to_frootn_patch(p)*dt
         ns_veg%frootn_xfer_patch(p) = ns_veg%frootn_xfer_patch(p) - nf_veg%frootn_xfer_to_frootn_patch(p)*dt

         if (woody(ivt(p)) == 1.0_r8) then
            ns_veg%livestemn_patch(p)       = ns_veg%livestemn_patch(p)       + nf_veg%livestemn_xfer_to_livestemn_patch(p)*dt
            ns_veg%livestemn_xfer_patch(p)  = ns_veg%livestemn_xfer_patch(p)  - nf_veg%livestemn_xfer_to_livestemn_patch(p)*dt
            ns_veg%deadstemn_patch(p)       = ns_veg%deadstemn_patch(p)       + nf_veg%deadstemn_xfer_to_deadstemn_patch(p)*dt
            ns_veg%deadstemn_xfer_patch(p)  = ns_veg%deadstemn_xfer_patch(p)  - nf_veg%deadstemn_xfer_to_deadstemn_patch(p)*dt
            ns_veg%livecrootn_patch(p)      = ns_veg%livecrootn_patch(p)      + nf_veg%livecrootn_xfer_to_livecrootn_patch(p)*dt
            ns_veg%livecrootn_xfer_patch(p) = ns_veg%livecrootn_xfer_patch(p) - nf_veg%livecrootn_xfer_to_livecrootn_patch(p)*dt
            ns_veg%deadcrootn_patch(p)      = ns_veg%deadcrootn_patch(p)      + nf_veg%deadcrootn_xfer_to_deadcrootn_patch(p)*dt
            ns_veg%deadcrootn_xfer_patch(p) = ns_veg%deadcrootn_xfer_patch(p) - nf_veg%deadcrootn_xfer_to_deadcrootn_patch(p)*dt
         end if

         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            ! lines here for consistency; the transfer terms are zero
            ns_veg%livestemn_patch(p)       = ns_veg%livestemn_patch(p)      + nf_veg%livestemn_xfer_to_livestemn_patch(p)*dt
            ns_veg%livestemn_xfer_patch(p)  = ns_veg%livestemn_xfer_patch(p) - nf_veg%livestemn_xfer_to_livestemn_patch(p)*dt
            ns_veg%grainn_patch(p)          = ns_veg%grainn_patch(p)         + nf_veg%grainn_xfer_to_grainn_patch(p)*dt
            ns_veg%grainn_xfer_patch(p)     = ns_veg%grainn_xfer_patch(p)    - nf_veg%grainn_xfer_to_grainn_patch(p)*dt
         end if

         ! phenology: litterfall and retranslocation fluxes
         ns_veg%leafn_patch(p)    = ns_veg%leafn_patch(p)    - nf_veg%leafn_to_litter_patch(p)*dt
         ns_veg%frootn_patch(p)   = ns_veg%frootn_patch(p)   - nf_veg%frootn_to_litter_patch(p)*dt
         ns_veg%leafn_patch(p)    = ns_veg%leafn_patch(p)    - nf_veg%leafn_to_retransn_patch(p)*dt
         ns_veg%retransn_patch(p) = ns_veg%retransn_patch(p) + nf_veg%leafn_to_retransn_patch(p)*dt

         ! live wood turnover and retranslocation fluxes
         if (woody(ivt(p)) == 1._r8) then
            ns_veg%livestemn_patch(p)    = ns_veg%livestemn_patch(p)  - nf_veg%livestemn_to_deadstemn_patch(p)*dt
            ns_veg%deadstemn_patch(p)    = ns_veg%deadstemn_patch(p)  + nf_veg%livestemn_to_deadstemn_patch(p)*dt
            ns_veg%livestemn_patch(p)    = ns_veg%livestemn_patch(p)  - nf_veg%livestemn_to_retransn_patch(p)*dt
            ns_veg%retransn_patch(p)     = ns_veg%retransn_patch(p)   + nf_veg%livestemn_to_retransn_patch(p)*dt
            ns_veg%livecrootn_patch(p)   = ns_veg%livecrootn_patch(p) - nf_veg%livecrootn_to_deadcrootn_patch(p)*dt
            ns_veg%deadcrootn_patch(p)   = ns_veg%deadcrootn_patch(p) + nf_veg%livecrootn_to_deadcrootn_patch(p)*dt
            ns_veg%livecrootn_patch(p)   = ns_veg%livecrootn_patch(p) - nf_veg%livecrootn_to_retransn_patch(p)*dt
            ns_veg%retransn_patch(p)     = ns_veg%retransn_patch(p)   + nf_veg%livecrootn_to_retransn_patch(p)*dt
         end if
         if (ivt(p) >= npcropmin) then ! Beth adds retrans from froot
            ns_veg%frootn_patch(p)       = ns_veg%frootn_patch(p)     - nf_veg%frootn_to_retransn_patch(p)*dt
            ns_veg%retransn_patch(p)     = ns_veg%retransn_patch(p)   + nf_veg%frootn_to_retransn_patch(p)*dt
            ns_veg%livestemn_patch(p)    = ns_veg%livestemn_patch(p)  - nf_veg%livestemn_to_litter_patch(p)*dt
            ns_veg%livestemn_patch(p)    = ns_veg%livestemn_patch(p)  - nf_veg%livestemn_to_retransn_patch(p)*dt
            ns_veg%retransn_patch(p)     = ns_veg%retransn_patch(p)   + nf_veg%livestemn_to_retransn_patch(p)*dt
            ns_veg%grainn_patch(p)       = ns_veg%grainn_patch(p)     - nf_veg%grainn_to_food_patch(p)*dt
         end if

         ! uptake from soil mineral N pool
         ns_veg%npool_patch(p) = ns_veg%npool_patch(p) + nf_veg%sminn_to_npool_patch(p)*dt

         ! deployment from retranslocation pool
         ns_veg%npool_patch(p)    = ns_veg%npool_patch(p)    + nf_veg%retransn_to_npool_patch(p)*dt
         ns_veg%retransn_patch(p) = ns_veg%retransn_patch(p) - nf_veg%retransn_to_npool_patch(p)*dt

         ! allocation fluxes
         ns_veg%npool_patch(p)           = ns_veg%npool_patch(p)          - nf_veg%npool_to_leafn_patch(p)*dt
         ns_veg%leafn_patch(p)           = ns_veg%leafn_patch(p)          + nf_veg%npool_to_leafn_patch(p)*dt
         ns_veg%npool_patch(p)           = ns_veg%npool_patch(p)          - nf_veg%npool_to_leafn_storage_patch(p)*dt
         ns_veg%leafn_storage_patch(p)   = ns_veg%leafn_storage_patch(p)  + nf_veg%npool_to_leafn_storage_patch(p)*dt
         ns_veg%npool_patch(p)           = ns_veg%npool_patch(p)          - nf_veg%npool_to_frootn_patch(p)*dt
         ns_veg%frootn_patch(p)          = ns_veg%frootn_patch(p)         + nf_veg%npool_to_frootn_patch(p)*dt
         ns_veg%npool_patch(p)           = ns_veg%npool_patch(p)          - nf_veg%npool_to_frootn_storage_patch(p)*dt
         ns_veg%frootn_storage_patch(p)  = ns_veg%frootn_storage_patch(p) + nf_veg%npool_to_frootn_storage_patch(p)*dt

         if (woody(ivt(p)) == 1._r8) then
            ns_veg%npool_patch(p)              = ns_veg%npool_patch(p)              - nf_veg%npool_to_livestemn_patch(p)*dt
            ns_veg%livestemn_patch(p)          = ns_veg%livestemn_patch(p)          + nf_veg%npool_to_livestemn_patch(p)*dt
            ns_veg%npool_patch(p)              = ns_veg%npool_patch(p)              - nf_veg%npool_to_livestemn_storage_patch(p)*dt
            ns_veg%livestemn_storage_patch(p)  = ns_veg%livestemn_storage_patch(p)  + nf_veg%npool_to_livestemn_storage_patch(p)*dt
            ns_veg%npool_patch(p)              = ns_veg%npool_patch(p)              - nf_veg%npool_to_deadstemn_patch(p)*dt
            ns_veg%deadstemn_patch(p)          = ns_veg%deadstemn_patch(p)          + nf_veg%npool_to_deadstemn_patch(p)*dt
            ns_veg%npool_patch(p)              = ns_veg%npool_patch(p)              - nf_veg%npool_to_deadstemn_storage_patch(p)*dt
            ns_veg%deadstemn_storage_patch(p)  = ns_veg%deadstemn_storage_patch(p)  + nf_veg%npool_to_deadstemn_storage_patch(p)*dt
            ns_veg%npool_patch(p)              = ns_veg%npool_patch(p)              - nf_veg%npool_to_livecrootn_patch(p)*dt
            ns_veg%livecrootn_patch(p)         = ns_veg%livecrootn_patch(p)         + nf_veg%npool_to_livecrootn_patch(p)*dt
            ns_veg%npool_patch(p)              = ns_veg%npool_patch(p)              - nf_veg%npool_to_livecrootn_storage_patch(p)*dt
            ns_veg%livecrootn_storage_patch(p) = ns_veg%livecrootn_storage_patch(p) + nf_veg%npool_to_livecrootn_storage_patch(p)*dt
            ns_veg%npool_patch(p)              = ns_veg%npool_patch(p)              - nf_veg%npool_to_deadcrootn_patch(p)*dt
            ns_veg%deadcrootn_patch(p)         = ns_veg%deadcrootn_patch(p)         + nf_veg%npool_to_deadcrootn_patch(p)*dt
            ns_veg%npool_patch(p)              = ns_veg%npool_patch(p)              - nf_veg%npool_to_deadcrootn_storage_patch(p)*dt
            ns_veg%deadcrootn_storage_patch(p) = ns_veg%deadcrootn_storage_patch(p) + nf_veg%npool_to_deadcrootn_storage_patch(p)*dt
         end if

         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            ns_veg%npool_patch(p)              = ns_veg%npool_patch(p)              - nf_veg%npool_to_livestemn_patch(p)*dt
            ns_veg%livestemn_patch(p)          = ns_veg%livestemn_patch(p)          + nf_veg%npool_to_livestemn_patch(p)*dt
            ns_veg%npool_patch(p)              = ns_veg%npool_patch(p)              - nf_veg%npool_to_livestemn_storage_patch(p)*dt
            ns_veg%livestemn_storage_patch(p)  = ns_veg%livestemn_storage_patch(p)  + nf_veg%npool_to_livestemn_storage_patch(p)*dt
            ns_veg%npool_patch(p)              = ns_veg%npool_patch(p)              - nf_veg%npool_to_grainn_patch(p)*dt
            ns_veg%grainn_patch(p)             = ns_veg%grainn_patch(p)             + nf_veg%npool_to_grainn_patch(p)*dt
            ns_veg%npool_patch(p)              = ns_veg%npool_patch(p)              - nf_veg%npool_to_grainn_storage_patch(p)*dt
            ns_veg%grainn_storage_patch(p)     = ns_veg%grainn_storage_patch(p)     + nf_veg%npool_to_grainn_storage_patch(p)*dt
         end if

         ! move storage pools into transfer pools
         ns_veg%leafn_storage_patch(p)  = ns_veg%leafn_storage_patch(p)  - nf_veg%leafn_storage_to_xfer_patch(p)*dt
         ns_veg%leafn_xfer_patch(p)     = ns_veg%leafn_xfer_patch(p)     + nf_veg%leafn_storage_to_xfer_patch(p)*dt
         ns_veg%frootn_storage_patch(p) = ns_veg%frootn_storage_patch(p) - nf_veg%frootn_storage_to_xfer_patch(p)*dt
         ns_veg%frootn_xfer_patch(p)    = ns_veg%frootn_xfer_patch(p)    + nf_veg%frootn_storage_to_xfer_patch(p)*dt

         if (woody(ivt(p)) == 1._r8) then
            ns_veg%livestemn_storage_patch(p)  = ns_veg%livestemn_storage_patch(p)  - nf_veg%livestemn_storage_to_xfer_patch(p)*dt
            ns_veg%livestemn_xfer_patch(p)     = ns_veg%livestemn_xfer_patch(p)     + nf_veg%livestemn_storage_to_xfer_patch(p)*dt
            ns_veg%deadstemn_storage_patch(p)  = ns_veg%deadstemn_storage_patch(p)  - nf_veg%deadstemn_storage_to_xfer_patch(p)*dt
            ns_veg%deadstemn_xfer_patch(p)     = ns_veg%deadstemn_xfer_patch(p)     + nf_veg%deadstemn_storage_to_xfer_patch(p)*dt
            ns_veg%livecrootn_storage_patch(p) = ns_veg%livecrootn_storage_patch(p) - nf_veg%livecrootn_storage_to_xfer_patch(p)*dt
            ns_veg%livecrootn_xfer_patch(p)    = ns_veg%livecrootn_xfer_patch(p)    + nf_veg%livecrootn_storage_to_xfer_patch(p)*dt
            ns_veg%deadcrootn_storage_patch(p) = ns_veg%deadcrootn_storage_patch(p) - nf_veg%deadcrootn_storage_to_xfer_patch(p)*dt
            ns_veg%deadcrootn_xfer_patch(p)    = ns_veg%deadcrootn_xfer_patch(p)    + nf_veg%deadcrootn_storage_to_xfer_patch(p)*dt
         end if

         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            ! lines here for consistency; the transfer terms are zero
            ns_veg%livestemn_storage_patch(p)  = ns_veg%livestemn_storage_patch(p) - nf_veg%livestemn_storage_to_xfer_patch(p)*dt
            ns_veg%livestemn_xfer_patch(p)     = ns_veg%livestemn_xfer_patch(p)    + nf_veg%livestemn_storage_to_xfer_patch(p)*dt
            ns_veg%grainn_storage_patch(p)     = ns_veg%grainn_storage_patch(p)    - nf_veg%grainn_storage_to_xfer_patch(p)*dt
            ns_veg%grainn_xfer_patch(p)        = ns_veg%grainn_xfer_patch(p)       + nf_veg%grainn_storage_to_xfer_patch(p)*dt
         end if

      end do

    end associate

  end subroutine NStateUpdate1

end module CNNStateUpdate1Mod
