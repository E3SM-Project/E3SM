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
  !! bgc interface & pflotran:
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
         do j = 1, nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               ! pft-level wood to column-level CWD (uncombusted wood)
               cs%decomp_cpools_vr_col(c,j,i_cwd) = cs%decomp_cpools_vr_col(c,j,i_cwd) + cf%fire_mortality_c_to_cwdc_col(c,j) * dt

               ! pft-level wood to column-level litter (uncombusted wood)
               cs%decomp_cpools_vr_col(c,j,i_met_lit) = cs%decomp_cpools_vr_col(c,j,i_met_lit) + cf%m_c_to_litr_met_fire_col(c,j)* dt
               cs%decomp_cpools_vr_col(c,j,i_cel_lit) = cs%decomp_cpools_vr_col(c,j,i_cel_lit) + cf%m_c_to_litr_cel_fire_col(c,j)* dt
               cs%decomp_cpools_vr_col(c,j,i_lig_lit) = cs%decomp_cpools_vr_col(c,j,i_lig_lit) + cf%m_c_to_litr_lig_fire_col(c,j)* dt
            end do
         end do

         ! litter and CWD losses to fire
         do l = 1, ndecomp_pools
            do j = 1, nlevdecomp
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  cs%decomp_cpools_vr_col(c,j,l) = cs%decomp_cpools_vr_col(c,j,l) - cf%m_decomp_cpools_to_fire_vr_col(c,j,l) * dt
               end do
            end do
         end do

      else

         ! column level carbon fluxes from fire
         do j = 1, nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               ! pft-level wood to column-level CWD (uncombusted wood)
               cf%bgc_cpool_ext_inputs_vr_col(c,j,i_cwd) = cf%bgc_cpool_ext_inputs_vr_col(c,j,i_cwd) + cf%fire_mortality_c_to_cwdc_col(c,j) * dt
               
               ! pft-level wood to column-level litter (uncombusted wood)
               cf%bgc_cpool_ext_inputs_vr_col(c,j,i_met_lit) = cf%bgc_cpool_ext_inputs_vr_col(c,j,i_met_lit) + cf%m_c_to_litr_met_fire_col(c,j)* dt
               cf%bgc_cpool_ext_inputs_vr_col(c,j,i_cel_lit) = cf%bgc_cpool_ext_inputs_vr_col(c,j,i_cel_lit) + cf%m_c_to_litr_cel_fire_col(c,j)* dt
               cf%bgc_cpool_ext_inputs_vr_col(c,j,i_lig_lit) = cf%bgc_cpool_ext_inputs_vr_col(c,j,i_lig_lit) + cf%m_c_to_litr_lig_fire_col(c,j)* dt
            end do
         end do
         
         ! litter and CWD losses to fire
         do l = 1, ndecomp_pools
            do j = 1, nlevdecomp
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  cf%bgc_cpool_ext_loss_vr_col(c,j,l) = cf%bgc_cpool_ext_loss_vr_col(c,j,l) + cf%m_decomp_cpools_to_fire_vr_col(c,j,l) * dt
               end do
            end do
         end do

      endif !


      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)
         
         ! pft-level carbon fluxes from fire
         ! displayed pools
         cs%leafc_patch(p)              = cs%leafc_patch(p)               - cf%m_leafc_to_fire_patch(p)            * dt
         cs%leafc_patch(p)              = cs%leafc_patch(p)               - cf%m_leafc_to_litter_fire_patch(p)     * dt
         cs%frootc_patch(p)             = cs%frootc_patch(p)              - cf%m_frootc_to_fire_patch(p)           * dt
         cs%frootc_patch(p)             = cs%frootc_patch(p)              - cf%m_frootc_to_litter_fire_patch(p)    * dt
         cs%livestemc_patch(p)          = cs%livestemc_patch(p)           - cf%m_livestemc_to_fire_patch(p)        * dt
         cs%livestemc_patch(p)          = cs%livestemc_patch(p)           - cf%m_livestemc_to_litter_fire_patch(p) * dt
         cs%deadstemc_patch(p)          = cs%deadstemc_patch(p)           - cf%m_deadstemc_to_fire_patch(p)        * dt
         cs%deadstemc_patch(p)          = cs%deadstemc_patch(p)           - cf%m_deadstemc_to_litter_fire_patch(p) * dt
         cs%livecrootc_patch(p)         = cs%livecrootc_patch(p)          - cf%m_livecrootc_to_fire_patch(p)       * dt
         cs%livecrootc_patch(p)         = cs%livecrootc_patch(p)          - cf%m_livecrootc_to_litter_fire_patch(p)* dt
         cs%deadcrootc_patch(p)         = cs%deadcrootc_patch(p)          - cf%m_deadcrootc_to_fire_patch(p)       * dt
         cs%deadcrootc_patch(p)         = cs%deadcrootc_patch(p)          - cf%m_deadcrootc_to_litter_fire_patch(p)* dt

         ! storage pools
         cs%leafc_storage_patch(p)      = cs%leafc_storage_patch(p)       - cf%m_leafc_storage_to_fire_patch(p)            * dt
         cs%leafc_storage_patch(p)      = cs%leafc_storage_patch(p)       - cf%m_leafc_storage_to_litter_fire_patch(p)     * dt
         cs%frootc_storage_patch(p)     = cs%frootc_storage_patch(p)      - cf%m_frootc_storage_to_fire_patch(p)           * dt
         cs%frootc_storage_patch(p)     = cs%frootc_storage_patch(p)      - cf%m_frootc_storage_to_litter_fire_patch(p)    * dt
         cs%livestemc_storage_patch(p)  = cs%livestemc_storage_patch(p)   - cf%m_livestemc_storage_to_fire_patch(p)        * dt
         cs%livestemc_storage_patch(p)  = cs%livestemc_storage_patch(p)   - cf%m_livestemc_storage_to_litter_fire_patch(p) * dt
         cs%deadstemc_storage_patch(p)  = cs%deadstemc_storage_patch(p)   - cf%m_deadstemc_storage_to_fire_patch(p)        * dt
         cs%deadstemc_storage_patch(p)  = cs%deadstemc_storage_patch(p)   - cf%m_deadstemc_storage_to_litter_fire_patch(p) * dt
         cs%livecrootc_storage_patch(p) = cs%livecrootc_storage_patch(p)  - cf%m_livecrootc_storage_to_fire_patch(p)       * dt
         cs%livecrootc_storage_patch(p) = cs%livecrootc_storage_patch(p)  - cf%m_livecrootc_storage_to_litter_fire_patch(p)* dt
         cs%deadcrootc_storage_patch(p) = cs%deadcrootc_storage_patch(p)  - cf%m_deadcrootc_storage_to_fire_patch(p)       * dt
         cs%deadcrootc_storage_patch(p) = cs%deadcrootc_storage_patch(p)  - cf%m_deadcrootc_storage_to_litter_fire_patch(p)* dt
         cs%gresp_storage_patch(p)      = cs%gresp_storage_patch(p)       - cf%m_gresp_storage_to_fire_patch(p)            * dt
         cs%gresp_storage_patch(p)      = cs%gresp_storage_patch(p)       - cf%m_gresp_storage_to_litter_fire_patch(p)     * dt

         ! transfer pools
         cs%leafc_xfer_patch(p)         = cs%leafc_xfer_patch(p)          - cf%m_leafc_xfer_to_fire_patch(p)            * dt
         cs%leafc_xfer_patch(p)         = cs%leafc_xfer_patch(p)          - cf%m_leafc_xfer_to_litter_fire_patch(p)     * dt
         cs%frootc_xfer_patch(p)        = cs%frootc_xfer_patch(p)         - cf%m_frootc_xfer_to_fire_patch(p)           * dt
         cs%frootc_xfer_patch(p)        = cs%frootc_xfer_patch(p)         - cf%m_frootc_xfer_to_litter_fire_patch(p)    * dt
         cs%livestemc_xfer_patch(p)     = cs%livestemc_xfer_patch(p)      - cf%m_livestemc_xfer_to_fire_patch(p)        * dt
         cs%livestemc_xfer_patch(p)     = cs%livestemc_xfer_patch(p)      - cf%m_livestemc_xfer_to_litter_fire_patch(p) * dt
         cs%deadstemc_xfer_patch(p)     = cs%deadstemc_xfer_patch(p)      - cf%m_deadstemc_xfer_to_fire_patch(p)        * dt
         cs%deadstemc_xfer_patch(p)     = cs%deadstemc_xfer_patch(p)      - cf%m_deadstemc_xfer_to_litter_fire_patch(p) * dt
         cs%livecrootc_xfer_patch(p)    = cs%livecrootc_xfer_patch(p)     - cf%m_livecrootc_xfer_to_fire_patch(p)       * dt
         cs%livecrootc_xfer_patch(p)    = cs%livecrootc_xfer_patch(p)     - cf%m_livecrootc_xfer_to_litter_fire_patch(p)* dt
         cs%deadcrootc_xfer_patch(p)    = cs%deadcrootc_xfer_patch(p)     - cf%m_deadcrootc_xfer_to_fire_patch(p)       * dt
         cs%deadcrootc_xfer_patch(p)    = cs%deadcrootc_xfer_patch(p)     - cf%m_deadcrootc_xfer_to_litter_fire_patch(p)* dt
         cs%gresp_xfer_patch(p)         = cs%gresp_xfer_patch(p)          - cf%m_gresp_xfer_to_fire_patch(p)            * dt
         cs%gresp_xfer_patch(p)         = cs%gresp_xfer_patch(p)          - cf%m_gresp_xfer_to_litter_fire_patch(p)     * dt

      end do ! end of pft loop


      
    end associate

  end subroutine CStateUpdate3

end module CNCStateUpdate3Mod
