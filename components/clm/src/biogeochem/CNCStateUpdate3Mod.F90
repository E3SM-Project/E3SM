module CNCStateUpdate3Mod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for carbon state variable update, mortality fluxes.
  !
  ! !USES:
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use abortutils       , only : endrun
  use clm_time_manager , only : get_step_size
  use clm_varpar       , only : nlevdecomp, ndecomp_pools, i_cwd, i_met_lit, i_cel_lit, i_lig_lit
  use CNCarbonStateType, only : carbonstate_type
  use CNCarbonFluxType , only : carbonflux_type
  ! bgc interface & pflotran:
  use clm_varctl       , only : use_pflotran, pf_cmode
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: CStateUpdate3
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CStateUpdate3( num_soilc, filter_soilc, num_soilp, filter_soilp, &
       carbonflux_vars, carbonstate_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update all the prognostic carbon state
    ! variables affected by fire fluxes
    !
    use tracer_varcon       , only : is_active_betr_bgc
    use subgridAveMod       , only : p2c        
    ! !ARGUMENTS:
    integer                , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(carbonflux_type)  , intent(inout) :: carbonflux_vars
    type(carbonstate_type) , intent(inout) :: carbonstate_vars
    !
    ! !LOCAL VARIABLES:
    integer :: c,p,j,l,k ! indices
    integer :: fp,fc     ! lake filter indices
    real(r8):: dt        ! radiation time step (seconds)
    !-----------------------------------------------------------------------

    associate(                   & 
         cf => carbonflux_vars , &
         cs => carbonstate_vars  &
         )

      ! set time steps
      dt = real( get_step_size(), r8 )

      if ( .not.is_active_betr_bgc )then
         ! column level carbon fluxes from fire
         if (.not.(use_pflotran .and. pf_cmode)) then
             do j = 1, nlevdecomp
                do fc = 1,num_soilc
                   c = filter_soilc(fc)
                   ! pft-level wood to column-level CWD (uncombusted wood)
                   cs%decomp_pools_vr_col(c,j,i_cwd) = cs%decomp_pools_vr_col(c,j,i_cwd) &
                        + cf%fire_mortality_to_cwd_col(c,j) * dt

                   ! pft-level wood to column-level litter (uncombusted wood)
                   cs%decomp_pools_vr_col(c,j,i_met_lit) = cs%decomp_pools_vr_col(c,j,i_met_lit) &
                        + cf%m_to_litr_met_fire_col(c,j)* dt
                   cs%decomp_pools_vr_col(c,j,i_cel_lit) = cs%decomp_pools_vr_col(c,j,i_cel_lit) &
                        + cf%m_c_to_litr_cel_fire_col(c,j)* dt
                   cs%decomp_pools_vr_col(c,j,i_lig_lit) = cs%decomp_pools_vr_col(c,j,i_lig_lit) &
                        + cf%m_to_litr_lig_fire_col(c,j)* dt
                end do
             end do
         end if !(.not.(use_pflotran .and. pf_cmode))

         ! litter and CWD losses to fire
         do l = 1, ndecomp_pools
            do j = 1, nlevdecomp
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  cs%decomp_pools_vr_col(c,j,l) = cs%decomp_pools_vr_col(c,j,l) - cf%m_decomp_cpools_to_fire_vr_col(c,j,l) * dt
               end do
            end do
         end do
      endif !


      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)
         
         ! pft-level carbon fluxes from fire
         ! displayed pools
         cs%leaf_patch(p)              = cs%leaf_patch(p)               - cf%m_leaf_to_fire_patch(p)            * dt
         cs%leaf_patch(p)              = cs%leaf_patch(p)               - cf%m_leaf_to_litter_fire_patch(p)     * dt
         cs%froot_patch(p)             = cs%froot_patch(p)              - cf%m_froot_to_fire_patch(p)           * dt
         cs%froot_patch(p)             = cs%froot_patch(p)              - cf%m_froot_to_litter_fire_patch(p)    * dt
         cs%livestem_patch(p)          = cs%livestem_patch(p)           - cf%m_livestem_to_fire_patch(p)        * dt
         cs%livestem_patch(p)          = cs%livestem_patch(p)           - cf%m_livestem_to_litter_fire_patch(p) * dt
         cs%deadstem_patch(p)          = cs%deadstem_patch(p)           - cf%m_deadstem_to_fire_patch(p)        * dt
         cs%deadstem_patch(p)          = cs%deadstem_patch(p)           - cf%m_deadstem_to_litter_fire_patch(p) * dt
         cs%livecroot_patch(p)         = cs%livecroot_patch(p)          - cf%m_livecroot_to_fire_patch(p)       * dt
         cs%livecroot_patch(p)         = cs%livecroot_patch(p)          - cf%m_livecroot_to_litter_fire_patch(p)* dt
         cs%deadcroot_patch(p)         = cs%deadcroot_patch(p)          - cf%m_deadcroot_to_fire_patch(p)       * dt
         cs%deadcroot_patch(p)         = cs%deadcroot_patch(p)          - cf%m_deadcroot_to_litter_fire_patch(p)* dt

         ! storage pools
         cs%leaf_storage_patch(p)      = cs%leaf_storage_patch(p)       - cf%m_leaf_storage_to_fire_patch(p)            * dt
         cs%leaf_storage_patch(p)      = cs%leaf_storage_patch(p)       - cf%m_leaf_storage_to_litter_fire_patch(p)     * dt
         cs%froot_storage_patch(p)     = cs%froot_storage_patch(p)      - cf%m_froot_storage_to_fire_patch(p)           * dt
         cs%froot_storage_patch(p)     = cs%froot_storage_patch(p)      - cf%m_froot_storage_to_litter_fire_patch(p)    * dt
         cs%livestem_storage_patch(p)  = cs%livestem_storage_patch(p)   - cf%m_livestem_storage_to_fire_patch(p)        * dt
         cs%livestem_storage_patch(p)  = cs%livestem_storage_patch(p)   - cf%m_livestem_storage_to_litter_fire_patch(p) * dt
         cs%deadstem_storage_patch(p)  = cs%deadstem_storage_patch(p)   - cf%m_deadstem_storage_to_fire_patch(p)        * dt
         cs%deadstem_storage_patch(p)  = cs%deadstem_storage_patch(p)   - cf%m_deadstem_storage_to_litter_fire_patch(p) * dt
         cs%livecroot_storage_patch(p) = cs%livecroot_storage_patch(p)  - cf%m_livecroot_storage_to_fire_patch(p)       * dt
         cs%livecroot_storage_patch(p) = cs%livecroot_storage_patch(p)  - cf%m_livecroot_storage_to_litter_fire_patch(p)* dt
         cs%deadcroot_storage_patch(p) = cs%deadcroot_storage_patch(p)  - cf%m_deadcroot_storage_to_fire_patch(p)       * dt
         cs%deadcroot_storage_patch(p) = cs%deadcroot_storage_patch(p)  - cf%m_deadcroot_storage_to_litter_fire_patch(p)* dt
         cs%gresp_storage_patch(p)      = cs%gresp_storage_patch(p)       - cf%m_gresp_storage_to_fire_patch(p)            * dt
         cs%gresp_storage_patch(p)      = cs%gresp_storage_patch(p)       - cf%m_gresp_storage_to_litter_fire_patch(p)     * dt

         ! transfer pools
         cs%leaf_xfer_patch(p)         = cs%leaf_xfer_patch(p)          - cf%m_leaf_xfer_to_fire_patch(p)            * dt
         cs%leaf_xfer_patch(p)         = cs%leaf_xfer_patch(p)          - cf%m_leaf_xfer_to_litter_fire_patch(p)     * dt
         cs%froot_xfer_patch(p)        = cs%froot_xfer_patch(p)         - cf%m_froot_xfer_to_fire_patch(p)           * dt
         cs%froot_xfer_patch(p)        = cs%froot_xfer_patch(p)         - cf%m_froot_xfer_to_litter_fire_patch(p)    * dt
         cs%livestem_xfer_patch(p)     = cs%livestem_xfer_patch(p)      - cf%m_livestem_xfer_to_fire_patch(p)        * dt
         cs%livestem_xfer_patch(p)     = cs%livestem_xfer_patch(p)      - cf%m_livestem_xfer_to_litter_fire_patch(p) * dt
         cs%deadstem_xfer_patch(p)     = cs%deadstem_xfer_patch(p)      - cf%m_deadstem_xfer_to_fire_patch(p)        * dt
         cs%deadstem_xfer_patch(p)     = cs%deadstem_xfer_patch(p)      - cf%m_deadstem_xfer_to_litter_fire_patch(p) * dt
         cs%livecroot_xfer_patch(p)    = cs%livecroot_xfer_patch(p)     - cf%m_livecroot_xfer_to_fire_patch(p)       * dt
         cs%livecroot_xfer_patch(p)    = cs%livecroot_xfer_patch(p)     - cf%m_livecroot_xfer_to_litter_fire_patch(p)* dt
         cs%deadcroot_xfer_patch(p)    = cs%deadcroot_xfer_patch(p)     - cf%m_deadcroot_xfer_to_fire_patch(p)       * dt
         cs%deadcroot_xfer_patch(p)    = cs%deadcroot_xfer_patch(p)     - cf%m_deadcroot_xfer_to_litter_fire_patch(p)* dt
         cs%gresp_xfer_patch(p)         = cs%gresp_xfer_patch(p)          - cf%m_gresp_xfer_to_fire_patch(p)            * dt
         cs%gresp_xfer_patch(p)         = cs%gresp_xfer_patch(p)          - cf%m_gresp_xfer_to_litter_fire_patch(p)     * dt
         cs%pool_patch(p)              = cs%pool_patch(p)               - cf%m_pool_to_fire_patch(p)                 * dt
         cs%pool_patch(p)              = cs%pool_patch(p)               - cf%m_pool_to_litter_fire_patch(p)          * dt


      end do ! end of pft loop


      
    end associate

  end subroutine CStateUpdate3

end module CNCStateUpdate3Mod
