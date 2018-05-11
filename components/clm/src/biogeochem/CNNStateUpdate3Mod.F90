module CNNStateUpdate3Mod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for nitrogen state variable update, mortality fluxes.
  ! Also, sminn leaching flux.
  !
  ! !USES:
  use shr_kind_mod        , only: r8 => shr_kind_r8
  use clm_varpar          , only: nlevdecomp, ndecomp_pools
  use clm_time_manager    , only : get_step_size
  use clm_varctl          , only : iulog, use_nitrif_denitrif
  use clm_varpar          , only : i_cwd, i_met_lit, i_cel_lit, i_lig_lit
  use CNNitrogenStateType , only : nitrogenstate_type
  use CNNitrogenFLuxType  , only : nitrogenflux_type
  ! bgc interface & pflotran:
  use clm_varctl          , only : use_pflotran, pf_cmode
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: NStateUpdate3
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine NStateUpdate3(num_soilc, filter_soilc, num_soilp, filter_soilp, &
       nitrogenflux_vars, nitrogenstate_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update all the prognostic nitrogen state
    ! variables affected by gap-phase mortality fluxes. Also the Sminn leaching flux.
    ! NOTE - associate statements have been removed where there are
    ! no science equations. This increases readability and maintainability.
    !
    use tracer_varcon, only : is_active_betr_bgc      
    ! !ARGUMENTS:
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                  , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    !
    ! !LOCAL VARIABLES:
    integer :: c,p,j,l,k        ! indices
    integer :: fp,fc      ! lake filter indices
    real(r8):: dt         ! radiation time step (seconds)
    !-----------------------------------------------------------------------

    associate(                      & 
         nf => nitrogenflux_vars  , &
         ns => nitrogenstate_vars   &
         )

      ! set time steps
      dt = real( get_step_size(), r8 )

      if (.not. is_active_betr_bgc) then
         do j = 1, nlevdecomp
            ! column loop
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               
               if (.not. use_nitrif_denitrif) then
                  ! mineral N loss due to leaching
                  ns%sminn_vr_col(c,j) = ns%sminn_vr_col(c,j) - nf%sminn_leached_vr_col(c,j) * dt
               else
                  ! mineral N loss due to leaching and runoff
                  ns%smin_no3_vr_col(c,j) = max( ns%smin_no3_vr_col(c,j) - &
                       ( nf%smin_no3_leached_vr_col(c,j) + nf%smin_no3_runoff_vr_col(c,j) ) * dt, 0._r8)
                  
                  ns%sminn_vr_col(c,j) = ns%smin_no3_vr_col(c,j) + ns%smin_nh4_vr_col(c,j)
                  if (use_pflotran .and. pf_cmode) then 
                        ns%sminn_vr_col(c,j) = ns%sminn_vr_col(c,j) + ns%smin_nh4sorb_vr_col(c,j)
                  end if
               end if
               
               if (.not.(use_pflotran .and. pf_cmode)) then
                   ! column level nitrogen fluxes from fire
                   ! pft-level wood to column-level CWD (uncombusted wood)
                   ns%decomp_pools_vr_col(c,j,i_cwd) = ns%decomp_pools_vr_col(c,j,i_cwd) &
                        + nf%fire_mortality_n_to_cwdn_col(c,j) * dt

                   ! pft-level wood to column-level litter (uncombusted wood)
                   ns%decomp_pools_vr_col(c,j,i_met_lit) = ns%decomp_pools_vr_col(c,j,i_met_lit) &
                        + nf%m_n_to_litr_met_fire_col(c,j)* dt
                   ns%decomp_pools_vr_col(c,j,i_cel_lit) = ns%decomp_pools_vr_col(c,j,i_cel_lit) &
                        + nf%m_n_to_litr_cel_fire_col(c,j)* dt
                   ns%decomp_pools_vr_col(c,j,i_lig_lit) = ns%decomp_pools_vr_col(c,j,i_lig_lit) &
                        + nf%m_n_to_litr_lig_fire_col(c,j)* dt
               end if !(.not.(use_pflotran .and. pf_cmode))
            end do ! end of column loop
         end do
         
         ! litter and CWD losses to fire
         do l = 1, ndecomp_pools
            do j = 1, nlevdecomp
               ! column loop
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  ns%decomp_pools_vr_col(c,j,l) = ns%decomp_pools_vr_col(c,j,l) - nf%m_decomp_npools_to_fire_vr_col(c,j,l) * dt
               end do
            end do
         end do

      endif
      ! patch-level nitrogen fluxes 
      
      do fp = 1,num_soilp
         p = filter_soilp(fp)
         
         !from fire displayed pools
         ns%leaf_patch(p)              =  ns%leaf_patch(p)      - nf%m_leafn_to_fire_patch(p)      * dt
         ns%froot_patch(p)             =  ns%froot_patch(p)     - nf%m_frootn_to_fire_patch(p)     * dt
         ns%livestem_patch(p)          =  ns%livestem_patch(p)  - nf%m_livestemn_to_fire_patch(p)  * dt
         ns%deadstem_patch(p)          =  ns%deadstem_patch(p)  - nf%m_deadstemn_to_fire_patch(p)  * dt
         ns%livecroot_patch(p)         =  ns%livecroot_patch(p) - nf%m_livecrootn_to_fire_patch(p) * dt
         ns%deadcroot_patch(p)         =  ns%deadcroot_patch(p) - nf%m_deadcrootn_to_fire_patch(p) * dt

         ns%leaf_patch(p)              =  ns%leaf_patch(p)      - nf%m_leafn_to_litter_fire_patch(p)           * dt
         ns%froot_patch(p)             =  ns%froot_patch(p)     - nf%m_frootn_to_litter_fire_patch(p)          * dt
         ns%livestem_patch(p)          =  ns%livestem_patch(p)  - nf%m_livestemn_to_litter_fire_patch(p)       * dt
         ns%deadstem_patch(p)          =  ns%deadstem_patch(p)  - nf%m_deadstemn_to_litter_fire_patch(p)       * dt
         ns%livecroot_patch(p)         =  ns%livecroot_patch(p) - nf%m_livecrootn_to_litter_fire_patch(p)      * dt
         ns%deadcroot_patch(p)         =  ns%deadcroot_patch(p) - nf%m_deadcrootn_to_litter_fire_patch(p)      * dt

         ! storage pools
         ns%leaf_storage_patch(p)      =  ns%leaf_storage_patch(p)      - nf%m_leafn_storage_to_fire_patch(p)      * dt
         ns%froot_storage_patch(p)     =  ns%froot_storage_patch(p)     - nf%m_frootn_storage_to_fire_patch(p)     * dt
         ns%livestem_storage_patch(p)  =  ns%livestem_storage_patch(p)  - nf%m_livestemn_storage_to_fire_patch(p)  * dt
         ns%deadstem_storage_patch(p)  =  ns%deadstem_storage_patch(p)  - nf%m_deadstemn_storage_to_fire_patch(p)  * dt
         ns%livecroot_storage_patch(p) =  ns%livecroot_storage_patch(p) - nf%m_livecrootn_storage_to_fire_patch(p) * dt
         ns%deadcroot_storage_patch(p) =  ns%deadcroot_storage_patch(p) - nf%m_deadcrootn_storage_to_fire_patch(p) * dt

         ns%leaf_storage_patch(p)      =  ns%leaf_storage_patch(p)      - nf%m_leafn_storage_to_litter_fire_patch(p)      * dt
         ns%froot_storage_patch(p)     =  ns%froot_storage_patch(p)     - nf%m_frootn_storage_to_litter_fire_patch(p)     * dt
         ns%livestem_storage_patch(p)  =  ns%livestem_storage_patch(p)  - nf%m_livestemn_storage_to_litter_fire_patch(p)  * dt
         ns%deadstem_storage_patch(p)  =  ns%deadstem_storage_patch(p)  - nf%m_deadstemn_storage_to_litter_fire_patch(p)  * dt
         ns%livecroot_storage_patch(p) =  ns%livecroot_storage_patch(p) - nf%m_livecrootn_storage_to_litter_fire_patch(p) * dt
         ns%deadcroot_storage_patch(p) =  ns%deadcroot_storage_patch(p) - nf%m_deadcrootn_storage_to_litter_fire_patch(p) * dt


         ! transfer pools
         ns%leaf_xfer_patch(p)         =  ns%leaf_xfer_patch(p)      - nf%m_leafn_xfer_to_fire_patch(p)      * dt
         ns%froot_xfer_patch(p)        =  ns%froot_xfer_patch(p)     - nf%m_frootn_xfer_to_fire_patch(p)     * dt
         ns%livestem_xfer_patch(p)     =  ns%livestem_xfer_patch(p)  - nf%m_livestemn_xfer_to_fire_patch(p)  * dt
         ns%deadstem_xfer_patch(p)     =  ns%deadstem_xfer_patch(p)  - nf%m_deadstemn_xfer_to_fire_patch(p)  * dt
         ns%livecroot_xfer_patch(p)    =  ns%livecroot_xfer_patch(p) - nf%m_livecrootn_xfer_to_fire_patch(p) * dt
         ns%deadcroot_xfer_patch(p)    =  ns%deadcroot_xfer_patch(p) - nf%m_deadcrootn_xfer_to_fire_patch(p) * dt

         ns%leaf_xfer_patch(p)         =  ns%leaf_xfer_patch(p)      - nf%m_leafn_xfer_to_litter_fire_patch(p)      * dt
         ns%froot_xfer_patch(p)        =  ns%froot_xfer_patch(p)     - nf%m_frootn_xfer_to_litter_fire_patch(p)     * dt
         ns%livestem_xfer_patch(p)     =  ns%livestem_xfer_patch(p)  - nf%m_livestemn_xfer_to_litter_fire_patch(p)  * dt
         ns%deadstem_xfer_patch(p)     =  ns%deadstem_xfer_patch(p)  - nf%m_deadstemn_xfer_to_litter_fire_patch(p)  * dt
         ns%livecroot_xfer_patch(p)    =  ns%livecroot_xfer_patch(p) - nf%m_livecrootn_xfer_to_litter_fire_patch(p) * dt
         ns%deadcroot_xfer_patch(p)    =  ns%deadcroot_xfer_patch(p) - nf%m_deadcrootn_xfer_to_litter_fire_patch(p) * dt

         ! retranslocated N pool
         ns%retransn_patch(p)           =  ns%retransn_patch(p) - nf%m_retransn_to_fire_patch(p)        * dt
         ns%retransn_patch(p)           =  ns%retransn_patch(p) - nf%m_retransn_to_litter_fire_patch(p) * dt
         ns%pool_patch(p)              =  ns%pool_patch(p)    - nf%m_npool_to_fire_patch(p)           * dt
         ns%pool_patch(p)              =  ns%pool_patch(p)    - nf%m_npool_to_litter_fire_patch(p)    * dt
      end do

    end associate 

  end subroutine NStateUpdate3

end module CNNStateUpdate3Mod
