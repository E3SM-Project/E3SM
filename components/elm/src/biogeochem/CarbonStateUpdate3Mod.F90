module CarbonStateUpdate3Mod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for carbon state variable update, mortality fluxes.
  !
  ! !USES:
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use abortutils       , only : endrun
  use elm_varpar       , only : nlevdecomp, ndecomp_pools, i_cwd, i_met_lit, i_cel_lit, i_lig_lit
  use elm_varctl       , only : ero_ccycle
  use CNCarbonStateType, only : carbonstate_type
  use CNCarbonFluxType , only : carbonflux_type
  use CNDecompCascadeConType , only : decomp_cascade_con
  use ColumnDataType         , only : column_carbon_state, column_carbon_flux
  use VegetationDataType     , only : vegetation_carbon_state, vegetation_carbon_flux
  ! bgc interface & pflotran:
  use elm_varctl       , only : use_pflotran, pf_cmode
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
        col_cs, veg_cs, col_cf, veg_cf, dt)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update all the prognostic carbon state
    ! variables affected by fire fluxes and also erosion flux
    !
      !$acc routine seq
    use tracer_varcon       , only : is_active_betr_bgc
    ! !ARGUMENTS:
    integer                , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(column_carbon_state),intent(inout):: col_cs
    type(vegetation_carbon_state),intent(inout) :: veg_cs
    type(column_carbon_flux)     ,intent(inout) :: col_cf
    type(vegetation_carbon_flux) ,intent(inout) :: veg_cf
    real(r8),   intent(in)    :: dt        ! radiation time step (seconds)

    !
    ! !LOCAL VARIABLES:
    integer :: c,p,j,l,k ! indices
    integer :: fp,fc     ! lake filter indices
    !-----------------------------------------------------------------------

      if ( .not.is_active_betr_bgc )then
         ! column level carbon fluxes from fire
         if (.not.(use_pflotran .and. pf_cmode)) then
             do j = 1, nlevdecomp
                do fc = 1,num_soilc
                   c = filter_soilc(fc)
                   ! pft-level wood to column-level CWD (uncombusted wood)
                   col_cs%decomp_cpools_vr(c,j,i_cwd) = col_cs%decomp_cpools_vr(c,j,i_cwd) &
                        + col_cf%fire_mortality_c_to_cwdc(c,j) * dt

                   ! pft-level wood to column-level litter (uncombusted wood)
                   col_cs%decomp_cpools_vr(c,j,i_met_lit) = col_cs%decomp_cpools_vr(c,j,i_met_lit) &
                        + col_cf%m_c_to_litr_met_fire(c,j)* dt
                   col_cs%decomp_cpools_vr(c,j,i_cel_lit) = col_cs%decomp_cpools_vr(c,j,i_cel_lit) &
                        + col_cf%m_c_to_litr_cel_fire(c,j)* dt
                   col_cs%decomp_cpools_vr(c,j,i_lig_lit) = col_cs%decomp_cpools_vr(c,j,i_lig_lit) &
                        + col_cf%m_c_to_litr_lig_fire(c,j)* dt
                end do
             end do
         end if !(.not.(use_pflotran .and. pf_cmode))

         ! litter and CWD losses to fire
         do l = 1, ndecomp_pools
            do j = 1, nlevdecomp
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  col_cs%decomp_cpools_vr(c,j,l) = col_cs%decomp_cpools_vr(c,j,l) - col_cf%m_decomp_cpools_to_fire_vr(c,j,l) * dt
               end do
            end do
         end do
      endif !

      ! SOM C losses due to erosion
      if ( ero_ccycle ) then
         do l = 1, ndecomp_pools
            if ( decomp_cascade_con%is_soil(l) ) then
               do j = 1, nlevdecomp
                  do fc = 1, num_soilc
                     c = filter_soilc(fc)
                     col_cs%decomp_cpools_vr(c,j,l) = col_cs%decomp_cpools_vr(c,j,l) - col_cf%decomp_cpools_yield_vr(c,j,l) * dt
                  end do
               end do
            end if
         end do
      end if


      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! pft-level carbon fluxes from fire
         ! displayed pools
         veg_cs%leafc(p)              = veg_cs%leafc(p)               - veg_cf%m_leafc_to_fire(p)            * dt
         veg_cs%leafc(p)              = veg_cs%leafc(p)               - veg_cf%m_leafc_to_litter_fire(p)     * dt
         veg_cs%frootc(p)             = veg_cs%frootc(p)              - veg_cf%m_frootc_to_fire(p)           * dt
         veg_cs%frootc(p)             = veg_cs%frootc(p)              - veg_cf%m_frootc_to_litter_fire(p)    * dt
         veg_cs%livestemc(p)          = veg_cs%livestemc(p)           - veg_cf%m_livestemc_to_fire(p)        * dt
         veg_cs%livestemc(p)          = veg_cs%livestemc(p)           - veg_cf%m_livestemc_to_litter_fire(p) * dt
         veg_cs%deadstemc(p)          = veg_cs%deadstemc(p)           - veg_cf%m_deadstemc_to_fire(p)        * dt
         veg_cs%deadstemc(p)          = veg_cs%deadstemc(p)           - veg_cf%m_deadstemc_to_litter_fire(p) * dt
         veg_cs%livecrootc(p)         = veg_cs%livecrootc(p)          - veg_cf%m_livecrootc_to_fire(p)       * dt
         veg_cs%livecrootc(p)         = veg_cs%livecrootc(p)          - veg_cf%m_livecrootc_to_litter_fire(p)* dt
         veg_cs%deadcrootc(p)         = veg_cs%deadcrootc(p)          - veg_cf%m_deadcrootc_to_fire(p)       * dt
         veg_cs%deadcrootc(p)         = veg_cs%deadcrootc(p)          - veg_cf%m_deadcrootc_to_litter_fire(p)* dt

         ! storage pools
         veg_cs%leafc_storage(p)      = veg_cs%leafc_storage(p)       - veg_cf%m_leafc_storage_to_fire(p)            * dt
         veg_cs%leafc_storage(p)      = veg_cs%leafc_storage(p)       - veg_cf%m_leafc_storage_to_litter_fire(p)     * dt
         veg_cs%frootc_storage(p)     = veg_cs%frootc_storage(p)      - veg_cf%m_frootc_storage_to_fire(p)           * dt
         veg_cs%frootc_storage(p)     = veg_cs%frootc_storage(p)      - veg_cf%m_frootc_storage_to_litter_fire(p)    * dt
         veg_cs%livestemc_storage(p)  = veg_cs%livestemc_storage(p)   - veg_cf%m_livestemc_storage_to_fire(p)        * dt
         veg_cs%livestemc_storage(p)  = veg_cs%livestemc_storage(p)   - veg_cf%m_livestemc_storage_to_litter_fire(p) * dt
         veg_cs%deadstemc_storage(p)  = veg_cs%deadstemc_storage(p)   - veg_cf%m_deadstemc_storage_to_fire(p)        * dt
         veg_cs%deadstemc_storage(p)  = veg_cs%deadstemc_storage(p)   - veg_cf%m_deadstemc_storage_to_litter_fire(p) * dt
         veg_cs%livecrootc_storage(p) = veg_cs%livecrootc_storage(p)  - veg_cf%m_livecrootc_storage_to_fire(p)       * dt
         veg_cs%livecrootc_storage(p) = veg_cs%livecrootc_storage(p)  - veg_cf%m_livecrootc_storage_to_litter_fire(p)* dt
         veg_cs%deadcrootc_storage(p) = veg_cs%deadcrootc_storage(p)  - veg_cf%m_deadcrootc_storage_to_fire(p)       * dt
         veg_cs%deadcrootc_storage(p) = veg_cs%deadcrootc_storage(p)  - veg_cf%m_deadcrootc_storage_to_litter_fire(p)* dt
         veg_cs%gresp_storage(p)      = veg_cs%gresp_storage(p)       - veg_cf%m_gresp_storage_to_fire(p)            * dt
         veg_cs%gresp_storage(p)      = veg_cs%gresp_storage(p)       - veg_cf%m_gresp_storage_to_litter_fire(p)     * dt

         ! transfer pools
         veg_cs%leafc_xfer(p)         = veg_cs%leafc_xfer(p)          - veg_cf%m_leafc_xfer_to_fire(p)            * dt
         veg_cs%leafc_xfer(p)         = veg_cs%leafc_xfer(p)          - veg_cf%m_leafc_xfer_to_litter_fire(p)     * dt
         veg_cs%frootc_xfer(p)        = veg_cs%frootc_xfer(p)         - veg_cf%m_frootc_xfer_to_fire(p)           * dt
         veg_cs%frootc_xfer(p)        = veg_cs%frootc_xfer(p)         - veg_cf%m_frootc_xfer_to_litter_fire(p)    * dt
         veg_cs%livestemc_xfer(p)     = veg_cs%livestemc_xfer(p)      - veg_cf%m_livestemc_xfer_to_fire(p)        * dt
         veg_cs%livestemc_xfer(p)     = veg_cs%livestemc_xfer(p)      - veg_cf%m_livestemc_xfer_to_litter_fire(p) * dt
         veg_cs%deadstemc_xfer(p)     = veg_cs%deadstemc_xfer(p)      - veg_cf%m_deadstemc_xfer_to_fire(p)        * dt
         veg_cs%deadstemc_xfer(p)     = veg_cs%deadstemc_xfer(p)      - veg_cf%m_deadstemc_xfer_to_litter_fire(p) * dt
         veg_cs%livecrootc_xfer(p)    = veg_cs%livecrootc_xfer(p)     - veg_cf%m_livecrootc_xfer_to_fire(p)       * dt
         veg_cs%livecrootc_xfer(p)    = veg_cs%livecrootc_xfer(p)     - veg_cf%m_livecrootc_xfer_to_litter_fire(p)* dt
         veg_cs%deadcrootc_xfer(p)    = veg_cs%deadcrootc_xfer(p)     - veg_cf%m_deadcrootc_xfer_to_fire(p)       * dt
         veg_cs%deadcrootc_xfer(p)    = veg_cs%deadcrootc_xfer(p)     - veg_cf%m_deadcrootc_xfer_to_litter_fire(p)* dt
         veg_cs%gresp_xfer(p)         = veg_cs%gresp_xfer(p)          - veg_cf%m_gresp_xfer_to_fire(p)            * dt
         veg_cs%gresp_xfer(p)         = veg_cs%gresp_xfer(p)          - veg_cf%m_gresp_xfer_to_litter_fire(p)     * dt
         veg_cs%cpool(p)              = veg_cs%cpool(p)               - veg_cf%m_cpool_to_fire(p)                 * dt
         veg_cs%cpool(p)              = veg_cs%cpool(p)               - veg_cf%m_cpool_to_litter_fire(p)          * dt


      end do ! end of pft loop

  end subroutine CarbonStateUpdate3

end module CarbonStateUpdate3Mod
