module CNNStateUpdate1Mod
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
  use EcophysconType         , only : ecophyscon
  use CNDecompCascadeConType , only : decomp_cascade_con
  use CNStateType            , only : cnstate_type
  use CNNitrogenFluxType     , only : nitrogenflux_type
  use CNNitrogenStateType    , only : nitrogenstate_type
  use PatchType              , only : pft
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
         ivt                   => pft%itype                                , & ! Input:  [integer  (:)     ]  pft vegetation type

         woody                 => ecophyscon%woody                         , & ! Input:  [real(r8) (:)     ]  binary flag for woody lifeform (1=woody, 0=not woody)

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

      ! seeding fluxes, from dynamic landcover
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         ns%seedn_col(c) = ns%seedn_col(c) - nf%dwt_seedn_to_leaf_col(c) * dt
         ns%seedn_col(c) = ns%seedn_col(c) - nf%dwt_seedn_to_deadstem_col(c) * dt
      end do

      if (.not. (is_active_betr_bgc) .and. .not.(use_pflotran .and. pf_cmode)) then

         do j = 1, nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)

               if (.not. use_nitrif_denitrif) then

                  ! N deposition and fixation
                  ns%sminn_vr_col(c,j) = ns%sminn_vr_col(c,j) + nf%ndep_to_sminn_col(c)*dt * ndep_prof(c,j)
                  ns%sminn_vr_col(c,j) = ns%sminn_vr_col(c,j) + nf%nfix_to_sminn_col(c)*dt * nfixation_prof(c,j)

               else

                  ! N deposition and fixation (put all into NH4 pool)
                  ns%smin_nh4_vr_col(c,j) = ns%smin_nh4_vr_col(c,j) + nf%ndep_to_sminn_col(c)*dt * ndep_prof(c,j)
                  ns%smin_nh4_vr_col(c,j) = ns%smin_nh4_vr_col(c,j) + nf%nfix_to_sminn_col(c)*dt * nfixation_prof(c,j)

               end if

               ! plant to litter fluxes
               ! phenology and dynamic landcover fluxes
               nf%decomp_npools_sourcesink_col(c,j,i_met_lit) = &
                    ( nf%phenology_n_to_litr_met_n_col(c,j) + nf%dwt_frootn_to_litr_met_n_col(c,j) ) * dt

               nf%decomp_npools_sourcesink_col(c,j,i_cel_lit) = &
                    ( nf%phenology_n_to_litr_cel_n_col(c,j) + nf%dwt_frootn_to_litr_cel_n_col(c,j) ) * dt

               nf%decomp_npools_sourcesink_col(c,j,i_lig_lit) = &
                    ( nf%phenology_n_to_litr_lig_n_col(c,j) + nf%dwt_frootn_to_litr_lig_n_col(c,j) ) * dt

               nf%decomp_npools_sourcesink_col(c,j,i_cwd)     = &
                    ( nf%dwt_livecrootn_to_cwdn_col(c,j)    + nf%dwt_deadcrootn_to_cwdn_col(c,j) )   * dt

            end do
         end do

         ! repeating N dep and fixation for crops
         if ( crop_prog )then
            do j = 1, nlevdecomp

               ! column loop
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  if (.not. use_nitrif_denitrif) then

                     ! N deposition and fixation
                     ns%sminn_vr_col(c,j) = ns%sminn_vr_col(c,j) + nf%fert_to_sminn_col(c)*dt * ndep_prof(c,j)
                     ns%sminn_vr_col(c,j) = ns%sminn_vr_col(c,j) + nf%soyfixn_to_sminn_col(c)*dt * nfixation_prof(c,j)
                  else

                     ! N deposition and fixation (put all into NH4 pool)
                     ns%smin_nh4_vr_col(c,j) = ns%smin_nh4_vr_col(c,j) + nf%fert_to_sminn_col(c)*dt * ndep_prof(c,j)
                     ns%smin_nh4_vr_col(c,j) = ns%smin_nh4_vr_col(c,j) + nf%soyfixn_to_sminn_col(c)*dt * nfixation_prof(c,j)

                  end if
               end do
            end do
         end if

         ! decomposition fluxes
         do k = 1, ndecomp_cascade_transitions
            do j = 1, nlevdecomp
               ! column loop
               do fc = 1,num_soilc
                  c = filter_soilc(fc)

                  nf%decomp_npools_sourcesink_col(c,j,cascade_donor_pool(k)) = &
                       nf%decomp_npools_sourcesink_col(c,j,cascade_donor_pool(k)) - &
                       nf%decomp_cascade_ntransfer_vr_col(c,j,k) * dt
               end do
            end do
         end do

         do k = 1, ndecomp_cascade_transitions
            if ( cascade_receiver_pool(k) /= 0 ) then  ! skip terminal transitions
               do j = 1, nlevdecomp
                  ! column loop
                  do fc = 1,num_soilc
                     c = filter_soilc(fc)

                     nf%decomp_npools_sourcesink_col(c,j,cascade_receiver_pool(k)) = &
                          nf%decomp_npools_sourcesink_col(c,j,cascade_receiver_pool(k)) + &
                          (nf%decomp_cascade_ntransfer_vr_col(c,j,k) + nf%decomp_cascade_sminn_flux_vr_col(c,j,k)) * dt
                  end do
               end do
            else  ! terminal transitions
               do j = 1, nlevdecomp
                  ! column loop
                  do fc = 1,num_soilc
                     c = filter_soilc(fc)
                     nf%decomp_npools_sourcesink_col(c,j,cascade_donor_pool(k)) = &
                          nf%decomp_npools_sourcesink_col(c,j,cascade_donor_pool(k)) - &
                          nf%decomp_cascade_sminn_flux_vr_col(c,j,k) * dt
                  end do
               end do
            end if
         end do

         if (.not. use_nitrif_denitrif) then

            !--------------------------------------------------------
            !-------------    NITRIF_DENITRIF OFF -------------------
            !--------------------------------------------------------

            ! immobilization/mineralization in litter-to-SOM and SOM-to-SOM fluxes and denitrification fluxes
            do k = 1, ndecomp_cascade_transitions
               if ( cascade_receiver_pool(k) /= 0 ) then  ! skip terminal transitions
                  do j = 1, nlevdecomp
                     ! column loop
                     do fc = 1,num_soilc
                        c = filter_soilc(fc)
                        ns%sminn_vr_col(c,j)  = ns%sminn_vr_col(c,j) - &
                             (nf%sminn_to_denit_decomp_cascade_vr_col(c,j,k) + nf%decomp_cascade_sminn_flux_vr_col(c,j,k))* dt
                     end do
                  end do
               else
                  do j = 1, nlevdecomp
                     ! column loop
                     do fc = 1,num_soilc
                        c = filter_soilc(fc)
                        ns%sminn_vr_col(c,j)  = ns%sminn_vr_col(c,j) - nf%sminn_to_denit_decomp_cascade_vr_col(c,j,k)* dt

                        ns%sminn_vr_col(c,j)  = ns%sminn_vr_col(c,j) + nf%decomp_cascade_sminn_flux_vr_col(c,j,k)* dt

                     end do
                  end do
               endif
            end do

            do j = 1, nlevdecomp
               ! column loop
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  ! "bulk denitrification"
                  ns%sminn_vr_col(c,j) = ns%sminn_vr_col(c,j) - nf%sminn_to_denit_excess_vr_col(c,j) * dt

                  ! total plant uptake from mineral N
                  ns%sminn_vr_col(c,j) = ns%sminn_vr_col(c,j) - nf%sminn_to_plant_vr_col(c,j)*dt

                  ! flux that prevents N limitation (when Carbon_only is set)
                  ns%sminn_vr_col(c,j) = ns%sminn_vr_col(c,j) + nf%supplement_to_sminn_vr_col(c,j)*dt
               end do
            end do

         else

            !--------------------------------------------------------
            !-------------    NITRIF_DENITRIF ON --------------------
            !--------------------------------------------------------

            do j = 1, nlevdecomp
               ! column loop
               do fc = 1,num_soilc
                  c = filter_soilc(fc)

                  ! mineralization fluxes (divert a fraction of this stream to nitrification flux, add the rest to NH4 pool)
                  ns%smin_nh4_vr_col(c,j) = ns%smin_nh4_vr_col(c,j) + nf%gross_nmin_vr_col(c,j)*dt

                  ! immobilization fluxes
                  ns%smin_nh4_vr_col(c,j) = ns%smin_nh4_vr_col(c,j) - nf%actual_immob_nh4_vr_col(c,j)*dt

                  ns%smin_no3_vr_col(c,j) = ns%smin_no3_vr_col(c,j) - nf%actual_immob_no3_vr_col(c,j)*dt

                  ! plant uptake fluxes
                  ns%smin_nh4_vr_col(c,j) = ns%smin_nh4_vr_col(c,j) - nf%smin_nh4_to_plant_vr_col(c,j)*dt

                  ns%smin_no3_vr_col(c,j) = ns%smin_no3_vr_col(c,j) - nf%smin_no3_to_plant_vr_col(c,j)*dt

                  ! Account for nitrification fluxes
                  ns%smin_nh4_vr_col(c,j) = ns%smin_nh4_vr_col(c,j) - nf%f_nit_vr_col(c,j) * dt

                  ns%smin_no3_vr_col(c,j) = ns%smin_no3_vr_col(c,j) + nf%f_nit_vr_col(c,j) * dt * (1._r8 - nitrif_n2o_loss_frac)

                  ! Account for denitrification fluxes
                  ns%smin_no3_vr_col(c,j) = ns%smin_no3_vr_col(c,j) - nf%f_denit_vr_col(c,j) * dt

                  ! flux that prevents N limitation (when Carbon_only is set; put all into NH4)
                  ns%smin_nh4_vr_col(c,j) = ns%smin_nh4_vr_col(c,j) + nf%supplement_to_sminn_vr_col(c,j)*dt

                  ! update diagnostic total
                  ns%sminn_vr_col(c,j) = ns%smin_nh4_vr_col(c,j) + ns%smin_no3_vr_col(c,j)

               end do ! end of column loop
            end do

         end if
      endif  !end if

      ! patch loop

      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! phenology: transfer growth fluxes
         ns%leafn_patch(p)       = ns%leafn_patch(p)       + nf%leafn_xfer_to_leafn_patch(p)*dt
         ns%leafn_xfer_patch(p)  = ns%leafn_xfer_patch(p)  - nf%leafn_xfer_to_leafn_patch(p)*dt
         ns%frootn_patch(p)      = ns%frootn_patch(p)      + nf%frootn_xfer_to_frootn_patch(p)*dt
         ns%frootn_xfer_patch(p) = ns%frootn_xfer_patch(p) - nf%frootn_xfer_to_frootn_patch(p)*dt

         if (woody(ivt(p)) == 1.0_r8) then
            ns%livestemn_patch(p)       = ns%livestemn_patch(p)       + nf%livestemn_xfer_to_livestemn_patch(p)*dt
            ns%livestemn_xfer_patch(p)  = ns%livestemn_xfer_patch(p)  - nf%livestemn_xfer_to_livestemn_patch(p)*dt
            ns%deadstemn_patch(p)       = ns%deadstemn_patch(p)       + nf%deadstemn_xfer_to_deadstemn_patch(p)*dt
            ns%deadstemn_xfer_patch(p)  = ns%deadstemn_xfer_patch(p)  - nf%deadstemn_xfer_to_deadstemn_patch(p)*dt
            ns%livecrootn_patch(p)      = ns%livecrootn_patch(p)      + nf%livecrootn_xfer_to_livecrootn_patch(p)*dt
            ns%livecrootn_xfer_patch(p) = ns%livecrootn_xfer_patch(p) - nf%livecrootn_xfer_to_livecrootn_patch(p)*dt
            ns%deadcrootn_patch(p)      = ns%deadcrootn_patch(p)      + nf%deadcrootn_xfer_to_deadcrootn_patch(p)*dt
            ns%deadcrootn_xfer_patch(p) = ns%deadcrootn_xfer_patch(p) - nf%deadcrootn_xfer_to_deadcrootn_patch(p)*dt
         end if

         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            ! lines here for consistency; the transfer terms are zero
            ns%livestemn_patch(p)       = ns%livestemn_patch(p)      + nf%livestemn_xfer_to_livestemn_patch(p)*dt
            ns%livestemn_xfer_patch(p)  = ns%livestemn_xfer_patch(p) - nf%livestemn_xfer_to_livestemn_patch(p)*dt
            ns%grainn_patch(p)          = ns%grainn_patch(p)         + nf%grainn_xfer_to_grainn_patch(p)*dt
            ns%grainn_xfer_patch(p)     = ns%grainn_xfer_patch(p)    - nf%grainn_xfer_to_grainn_patch(p)*dt
         end if

         ! phenology: litterfall and retranslocation fluxes
         ns%leafn_patch(p)    = ns%leafn_patch(p)    - nf%leafn_to_litter_patch(p)*dt
         ns%frootn_patch(p)   = ns%frootn_patch(p)   - nf%frootn_to_litter_patch(p)*dt
         ns%leafn_patch(p)    = ns%leafn_patch(p)    - nf%leafn_to_retransn_patch(p)*dt
         ns%retransn_patch(p) = ns%retransn_patch(p) + nf%leafn_to_retransn_patch(p)*dt

         ! live wood turnover and retranslocation fluxes
         if (woody(ivt(p)) == 1._r8) then
            ns%livestemn_patch(p)  = ns%livestemn_patch(p)  - nf%livestemn_to_deadstemn_patch(p)*dt
            ns%deadstemn_patch(p)  = ns%deadstemn_patch(p)  + nf%livestemn_to_deadstemn_patch(p)*dt
            ns%livestemn_patch(p)  = ns%livestemn_patch(p)  - nf%livestemn_to_retransn_patch(p)*dt
            ns%retransn_patch(p)   = ns%retransn_patch(p)   + nf%livestemn_to_retransn_patch(p)*dt
            ns%livecrootn_patch(p) = ns%livecrootn_patch(p) - nf%livecrootn_to_deadcrootn_patch(p)*dt
            ns%deadcrootn_patch(p) = ns%deadcrootn_patch(p) + nf%livecrootn_to_deadcrootn_patch(p)*dt
            ns%livecrootn_patch(p) = ns%livecrootn_patch(p) - nf%livecrootn_to_retransn_patch(p)*dt
            ns%retransn_patch(p)   = ns%retransn_patch(p)   + nf%livecrootn_to_retransn_patch(p)*dt
         end if
         if (ivt(p) >= npcropmin) then ! Beth adds retrans from froot
            ns%frootn_patch(p)     = ns%frootn_patch(p)     - nf%frootn_to_retransn_patch(p)*dt
            ns%retransn_patch(p)   = ns%retransn_patch(p)   + nf%frootn_to_retransn_patch(p)*dt
            ns%livestemn_patch(p)  = ns%livestemn_patch(p)  - nf%livestemn_to_litter_patch(p)*dt
            ns%livestemn_patch(p)  = ns%livestemn_patch(p)  - nf%livestemn_to_retransn_patch(p)*dt
            ns%retransn_patch(p)   = ns%retransn_patch(p)   + nf%livestemn_to_retransn_patch(p)*dt
            ns%grainn_patch(p)     = ns%grainn_patch(p)     - nf%grainn_to_food_patch(p)*dt
         end if

         ! uptake from soil mineral N pool
         ns%npool_patch(p) = &
              ns%npool_patch(p) + nf%sminn_to_npool_patch(p)*dt

         ! deployment from retranslocation pool
         ns%npool_patch(p)    = ns%npool_patch(p)    + nf%retransn_to_npool_patch(p)*dt
         ns%retransn_patch(p) = ns%retransn_patch(p) - nf%retransn_to_npool_patch(p)*dt

         ! allocation fluxes
         ns%npool_patch(p)           = ns%npool_patch(p)          - nf%npool_to_leafn_patch(p)*dt
         ns%leafn_patch(p)           = ns%leafn_patch(p)          + nf%npool_to_leafn_patch(p)*dt
         ns%npool_patch(p)           = ns%npool_patch(p)          - nf%npool_to_leafn_storage_patch(p)*dt
         ns%leafn_storage_patch(p)   = ns%leafn_storage_patch(p)  + nf%npool_to_leafn_storage_patch(p)*dt
         ns%npool_patch(p)           = ns%npool_patch(p)          - nf%npool_to_frootn_patch(p)*dt
         ns%frootn_patch(p)          = ns%frootn_patch(p)         + nf%npool_to_frootn_patch(p)*dt
         ns%npool_patch(p)           = ns%npool_patch(p)          - nf%npool_to_frootn_storage_patch(p)*dt
         ns%frootn_storage_patch(p)  = ns%frootn_storage_patch(p) + nf%npool_to_frootn_storage_patch(p)*dt

         if (woody(ivt(p)) == 1._r8) then
            ns%npool_patch(p)              = ns%npool_patch(p)              - nf%npool_to_livestemn_patch(p)*dt
            ns%livestemn_patch(p)          = ns%livestemn_patch(p)          + nf%npool_to_livestemn_patch(p)*dt
            ns%npool_patch(p)              = ns%npool_patch(p)              - nf%npool_to_livestemn_storage_patch(p)*dt
            ns%livestemn_storage_patch(p)  = ns%livestemn_storage_patch(p)  + nf%npool_to_livestemn_storage_patch(p)*dt
            ns%npool_patch(p)              = ns%npool_patch(p)              - nf%npool_to_deadstemn_patch(p)*dt
            ns%deadstemn_patch(p)          = ns%deadstemn_patch(p)          + nf%npool_to_deadstemn_patch(p)*dt
            ns%npool_patch(p)              = ns%npool_patch(p)              - nf%npool_to_deadstemn_storage_patch(p)*dt
            ns%deadstemn_storage_patch(p)  = ns%deadstemn_storage_patch(p)  + nf%npool_to_deadstemn_storage_patch(p)*dt
            ns%npool_patch(p)              = ns%npool_patch(p)              - nf%npool_to_livecrootn_patch(p)*dt
            ns%livecrootn_patch(p)         = ns%livecrootn_patch(p)         + nf%npool_to_livecrootn_patch(p)*dt
            ns%npool_patch(p)              = ns%npool_patch(p)              - nf%npool_to_livecrootn_storage_patch(p)*dt
            ns%livecrootn_storage_patch(p) = ns%livecrootn_storage_patch(p) + nf%npool_to_livecrootn_storage_patch(p)*dt
            ns%npool_patch(p)              = ns%npool_patch(p)              - nf%npool_to_deadcrootn_patch(p)*dt
            ns%deadcrootn_patch(p)         = ns%deadcrootn_patch(p)         + nf%npool_to_deadcrootn_patch(p)*dt
            ns%npool_patch(p)              = ns%npool_patch(p)              - nf%npool_to_deadcrootn_storage_patch(p)*dt
            ns%deadcrootn_storage_patch(p) = ns%deadcrootn_storage_patch(p) + nf%npool_to_deadcrootn_storage_patch(p)*dt
         end if

         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            ns%npool_patch(p)              = ns%npool_patch(p)              - nf%npool_to_livestemn_patch(p)*dt
            ns%livestemn_patch(p)          = ns%livestemn_patch(p)          + nf%npool_to_livestemn_patch(p)*dt
            ns%npool_patch(p)              = ns%npool_patch(p)              - nf%npool_to_livestemn_storage_patch(p)*dt
            ns%livestemn_storage_patch(p)  = ns%livestemn_storage_patch(p)  + nf%npool_to_livestemn_storage_patch(p)*dt
            ns%npool_patch(p)              = ns%npool_patch(p)              - nf%npool_to_grainn_patch(p)*dt
            ns%grainn_patch(p)             = ns%grainn_patch(p)             + nf%npool_to_grainn_patch(p)*dt
            ns%npool_patch(p)              = ns%npool_patch(p)              - nf%npool_to_grainn_storage_patch(p)*dt
            ns%grainn_storage_patch(p)     = ns%grainn_storage_patch(p)     + nf%npool_to_grainn_storage_patch(p)*dt
         end if

         ! move storage pools into transfer pools
         ns%leafn_storage_patch(p)  = ns%leafn_storage_patch(p)  - nf%leafn_storage_to_xfer_patch(p)*dt
         ns%leafn_xfer_patch(p)     = ns%leafn_xfer_patch(p)     + nf%leafn_storage_to_xfer_patch(p)*dt
         ns%frootn_storage_patch(p) = ns%frootn_storage_patch(p) - nf%frootn_storage_to_xfer_patch(p)*dt
         ns%frootn_xfer_patch(p)    = ns%frootn_xfer_patch(p)    + nf%frootn_storage_to_xfer_patch(p)*dt

         if (woody(ivt(p)) == 1._r8) then
            ns%livestemn_storage_patch(p)  = ns%livestemn_storage_patch(p)  - nf%livestemn_storage_to_xfer_patch(p)*dt
            ns%livestemn_xfer_patch(p)     = ns%livestemn_xfer_patch(p)     + nf%livestemn_storage_to_xfer_patch(p)*dt
            ns%deadstemn_storage_patch(p)  = ns%deadstemn_storage_patch(p)  - nf%deadstemn_storage_to_xfer_patch(p)*dt
            ns%deadstemn_xfer_patch(p)     = ns%deadstemn_xfer_patch(p)     + nf%deadstemn_storage_to_xfer_patch(p)*dt
            ns%livecrootn_storage_patch(p) = ns%livecrootn_storage_patch(p) - nf%livecrootn_storage_to_xfer_patch(p)*dt
            ns%livecrootn_xfer_patch(p)    = ns%livecrootn_xfer_patch(p)    + nf%livecrootn_storage_to_xfer_patch(p)*dt
            ns%deadcrootn_storage_patch(p) = ns%deadcrootn_storage_patch(p) - nf%deadcrootn_storage_to_xfer_patch(p)*dt
            ns%deadcrootn_xfer_patch(p)    = ns%deadcrootn_xfer_patch(p)    + nf%deadcrootn_storage_to_xfer_patch(p)*dt
         end if

         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            ! lines here for consistency; the transfer terms are zero
            ns%livestemn_storage_patch(p)  = ns%livestemn_storage_patch(p) - nf%livestemn_storage_to_xfer_patch(p)*dt
            ns%livestemn_xfer_patch(p)     = ns%livestemn_xfer_patch(p)    + nf%livestemn_storage_to_xfer_patch(p)*dt
            ns%grainn_storage_patch(p)     = ns%grainn_storage_patch(p)    - nf%grainn_storage_to_xfer_patch(p)*dt
            ns%grainn_xfer_patch(p)        = ns%grainn_xfer_patch(p)       + nf%grainn_storage_to_xfer_patch(p)*dt
         end if

      end do

    end associate

  end subroutine NStateUpdate1

end module CNNStateUpdate1Mod
