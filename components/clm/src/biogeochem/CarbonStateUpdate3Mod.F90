module CarbonStateUpdate3Mod

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
  use ColumnDataType   , only : column_carbon_state
  use VegetationDataType     , only : vegetation_carbon_state
  ! bgc interface & pflotran:
  use clm_varctl       , only : use_pflotran, pf_cmode
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: CarbonStateUpdate3
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CarbonStateUpdate3( num_soilc, filter_soilc, num_soilp, filter_soilp, &
       carbonflux_vars, carbonstate_vars, col_csv2, veg_csv2)
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
    type(column_carbon_state),intent(inout):: col_csv2
    type(vegetation_carbon_state),intent(inout) :: veg_csv2
    !
    ! !LOCAL VARIABLES:
    integer :: c,p,j,l,k ! indices
    integer :: fp,fc     ! lake filter indices
    real(r8):: dt        ! radiation time step (seconds)
    !-----------------------------------------------------------------------

    associate(                   & 
         cf => carbonflux_vars , &
         cs => carbonstate_vars , &
         csv2 => col_csv2       , &
         vcsv2 => veg_csv2        &
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
                   csv2%decomp_cpools_vr(c,j,i_cwd) = csv2%decomp_cpools_vr(c,j,i_cwd) &
                        + cf%fire_mortality_c_to_cwdc_col(c,j) * dt

                   ! pft-level wood to column-level litter (uncombusted wood)
                   csv2%decomp_cpools_vr(c,j,i_met_lit) = csv2%decomp_cpools_vr(c,j,i_met_lit) &
                        + cf%m_c_to_litr_met_fire_col(c,j)* dt
                   csv2%decomp_cpools_vr(c,j,i_cel_lit) = csv2%decomp_cpools_vr(c,j,i_cel_lit) &
                        + cf%m_c_to_litr_cel_fire_col(c,j)* dt
                   csv2%decomp_cpools_vr(c,j,i_lig_lit) = csv2%decomp_cpools_vr(c,j,i_lig_lit) &
                        + cf%m_c_to_litr_lig_fire_col(c,j)* dt
                end do
             end do
         end if !(.not.(use_pflotran .and. pf_cmode))

         ! litter and CWD losses to fire
         do l = 1, ndecomp_pools
            do j = 1, nlevdecomp
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  csv2%decomp_cpools_vr(c,j,l) = csv2%decomp_cpools_vr(c,j,l) - cf%m_decomp_cpools_to_fire_vr_col(c,j,l) * dt
               end do
            end do
         end do
      endif !


      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)
         
         ! pft-level carbon fluxes from fire
         ! displayed pools
         vcsv2%leafc(p)              = vcsv2%leafc(p)               - cf%m_leafc_to_fire_patch(p)            * dt
         vcsv2%leafc(p)              = vcsv2%leafc(p)               - cf%m_leafc_to_litter_fire_patch(p)     * dt
         vcsv2%frootc(p)             = vcsv2%frootc(p)              - cf%m_frootc_to_fire_patch(p)           * dt
         vcsv2%frootc(p)             = vcsv2%frootc(p)              - cf%m_frootc_to_litter_fire_patch(p)    * dt
         vcsv2%livestemc(p)          = vcsv2%livestemc(p)           - cf%m_livestemc_to_fire_patch(p)        * dt
         vcsv2%livestemc(p)          = vcsv2%livestemc(p)           - cf%m_livestemc_to_litter_fire_patch(p) * dt
         vcsv2%deadstemc(p)          = vcsv2%deadstemc(p)           - cf%m_deadstemc_to_fire_patch(p)        * dt
         vcsv2%deadstemc(p)          = vcsv2%deadstemc(p)           - cf%m_deadstemc_to_litter_fire_patch(p) * dt
         vcsv2%livecrootc(p)         = vcsv2%livecrootc(p)          - cf%m_livecrootc_to_fire_patch(p)       * dt
         vcsv2%livecrootc(p)         = vcsv2%livecrootc(p)          - cf%m_livecrootc_to_litter_fire_patch(p)* dt
         vcsv2%deadcrootc(p)         = vcsv2%deadcrootc(p)          - cf%m_deadcrootc_to_fire_patch(p)       * dt
         vcsv2%deadcrootc(p)         = vcsv2%deadcrootc(p)          - cf%m_deadcrootc_to_litter_fire_patch(p)* dt

         ! storage pools
         vcsv2%leafc_storage(p)      = vcsv2%leafc_storage(p)       - cf%m_leafc_storage_to_fire_patch(p)            * dt
         vcsv2%leafc_storage(p)      = vcsv2%leafc_storage(p)       - cf%m_leafc_storage_to_litter_fire_patch(p)     * dt
         vcsv2%frootc_storage(p)     = vcsv2%frootc_storage(p)      - cf%m_frootc_storage_to_fire_patch(p)           * dt
         vcsv2%frootc_storage(p)     = vcsv2%frootc_storage(p)      - cf%m_frootc_storage_to_litter_fire_patch(p)    * dt
         vcsv2%livestemc_storage(p)  = vcsv2%livestemc_storage(p)   - cf%m_livestemc_storage_to_fire_patch(p)        * dt
         vcsv2%livestemc_storage(p)  = vcsv2%livestemc_storage(p)   - cf%m_livestemc_storage_to_litter_fire_patch(p) * dt
         vcsv2%deadstemc_storage(p)  = vcsv2%deadstemc_storage(p)   - cf%m_deadstemc_storage_to_fire_patch(p)        * dt
         vcsv2%deadstemc_storage(p)  = vcsv2%deadstemc_storage(p)   - cf%m_deadstemc_storage_to_litter_fire_patch(p) * dt
         vcsv2%livecrootc_storage(p) = vcsv2%livecrootc_storage(p)  - cf%m_livecrootc_storage_to_fire_patch(p)       * dt
         vcsv2%livecrootc_storage(p) = vcsv2%livecrootc_storage(p)  - cf%m_livecrootc_storage_to_litter_fire_patch(p)* dt
         vcsv2%deadcrootc_storage(p) = vcsv2%deadcrootc_storage(p)  - cf%m_deadcrootc_storage_to_fire_patch(p)       * dt
         vcsv2%deadcrootc_storage(p) = vcsv2%deadcrootc_storage(p)  - cf%m_deadcrootc_storage_to_litter_fire_patch(p)* dt
         vcsv2%gresp_storage(p)      = vcsv2%gresp_storage(p)       - cf%m_gresp_storage_to_fire_patch(p)            * dt
         vcsv2%gresp_storage(p)      = vcsv2%gresp_storage(p)       - cf%m_gresp_storage_to_litter_fire_patch(p)     * dt

         ! transfer pools
         vcsv2%leafc_xfer(p)         = vcsv2%leafc_xfer(p)          - cf%m_leafc_xfer_to_fire_patch(p)            * dt
         vcsv2%leafc_xfer(p)         = vcsv2%leafc_xfer(p)          - cf%m_leafc_xfer_to_litter_fire_patch(p)     * dt
         vcsv2%frootc_xfer(p)        = vcsv2%frootc_xfer(p)         - cf%m_frootc_xfer_to_fire_patch(p)           * dt
         vcsv2%frootc_xfer(p)        = vcsv2%frootc_xfer(p)         - cf%m_frootc_xfer_to_litter_fire_patch(p)    * dt
         vcsv2%livestemc_xfer(p)     = vcsv2%livestemc_xfer(p)      - cf%m_livestemc_xfer_to_fire_patch(p)        * dt
         vcsv2%livestemc_xfer(p)     = vcsv2%livestemc_xfer(p)      - cf%m_livestemc_xfer_to_litter_fire_patch(p) * dt
         vcsv2%deadstemc_xfer(p)     = vcsv2%deadstemc_xfer(p)      - cf%m_deadstemc_xfer_to_fire_patch(p)        * dt
         vcsv2%deadstemc_xfer(p)     = vcsv2%deadstemc_xfer(p)      - cf%m_deadstemc_xfer_to_litter_fire_patch(p) * dt
         vcsv2%livecrootc_xfer(p)    = vcsv2%livecrootc_xfer(p)     - cf%m_livecrootc_xfer_to_fire_patch(p)       * dt
         vcsv2%livecrootc_xfer(p)    = vcsv2%livecrootc_xfer(p)     - cf%m_livecrootc_xfer_to_litter_fire_patch(p)* dt
         vcsv2%deadcrootc_xfer(p)    = vcsv2%deadcrootc_xfer(p)     - cf%m_deadcrootc_xfer_to_fire_patch(p)       * dt
         vcsv2%deadcrootc_xfer(p)    = vcsv2%deadcrootc_xfer(p)     - cf%m_deadcrootc_xfer_to_litter_fire_patch(p)* dt
         vcsv2%gresp_xfer(p)         = vcsv2%gresp_xfer(p)          - cf%m_gresp_xfer_to_fire_patch(p)            * dt
         vcsv2%gresp_xfer(p)         = vcsv2%gresp_xfer(p)          - cf%m_gresp_xfer_to_litter_fire_patch(p)     * dt
         vcsv2%cpool(p)              = vcsv2%cpool(p)               - cf%m_cpool_to_fire_patch(p)                 * dt
         vcsv2%cpool(p)              = vcsv2%cpool(p)               - cf%m_cpool_to_litter_fire_patch(p)          * dt


      end do ! end of pft loop


      
    end associate

  end subroutine CarbonStateUpdate3

end module CarbonStateUpdate3Mod
