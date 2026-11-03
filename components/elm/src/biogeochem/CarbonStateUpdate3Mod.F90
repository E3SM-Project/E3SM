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
  use elm_varctl       , only : use_erosion, ero_ccycle
  use CNDecompCascadeConType , only : decomp_cascade_con
  use ColumnStateTypes, only : column_carbon_state
  use ColumnFluxTypes, only : column_carbon_flux
  use VegetationStateTypes, only : vegetation_carbon_state
  use VegetationFluxTypes, only : vegetation_carbon_flux
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

      ! column level carbon fluxes from fire
      if (.not.(use_pflotran .and. pf_cmode)) then
          do j = 1, nlevdecomp
             do fc = 1,num_soilc
                c = filter_soilc(fc)
                ! pft-level wood to column-level CWD (uncombusted wood)
                col_cs%core%decomp_pools_vr(c,j,i_cwd) = col_cs%core%decomp_pools_vr(c,j,i_cwd) &
                     + col_cf%core%fire_mortality_to_cwd(c,j) * dt

                ! pft-level wood to column-level litter (uncombusted wood)
                col_cs%core%decomp_pools_vr(c,j,i_met_lit) = col_cs%core%decomp_pools_vr(c,j,i_met_lit) &
                     + col_cf%core%m_to_litr_met_fire(c,j)* dt
                col_cs%core%decomp_pools_vr(c,j,i_cel_lit) = col_cs%core%decomp_pools_vr(c,j,i_cel_lit) &
                     + col_cf%core%m_to_litr_cel_fire(c,j)* dt
                col_cs%core%decomp_pools_vr(c,j,i_lig_lit) = col_cs%core%decomp_pools_vr(c,j,i_lig_lit) &
                     + col_cf%core%m_to_litr_lig_fire(c,j)* dt
             end do
          end do
      end if !(.not.(use_pflotran .and. pf_cmode))

      ! litter and CWD losses to fire
      do l = 1, ndecomp_pools
         do j = 1, nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               col_cs%core%decomp_pools_vr(c,j,l) = col_cs%core%decomp_pools_vr(c,j,l) - col_cf%core%m_decomp_pools_to_fire_vr(c,j,l) * dt
            end do
         end do
      end do

      ! SOM C losses due to erosion
      if ( ero_ccycle ) then
         do l = 1, ndecomp_pools
            if ( decomp_cascade_con%is_soil(l) ) then
               do j = 1, nlevdecomp
                  do fc = 1, num_soilc
                     c = filter_soilc(fc)
                     col_cs%core%decomp_pools_vr(c,j,l) = col_cs%core%decomp_pools_vr(c,j,l) - col_cf%core%decomp_pools_yield_vr(c,j,l) * dt
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
         veg_cs%core%leaf(p)              = veg_cs%core%leaf(p)               - veg_cf%core%m_leaf_to_fire(p)            * dt
         veg_cs%core%leaf(p)              = veg_cs%core%leaf(p)               - veg_cf%core%m_leaf_to_litter_fire(p)     * dt
         veg_cs%core%froot(p)             = veg_cs%core%froot(p)              - veg_cf%core%m_froot_to_fire(p)           * dt
         veg_cs%core%froot(p)             = veg_cs%core%froot(p)              - veg_cf%core%m_froot_to_litter_fire(p)    * dt
         veg_cs%core%livestem(p)          = veg_cs%core%livestem(p)           - veg_cf%core%m_livestem_to_fire(p)        * dt
         veg_cs%core%livestem(p)          = veg_cs%core%livestem(p)           - veg_cf%core%m_livestem_to_litter_fire(p) * dt
         veg_cs%core%deadstem(p)          = veg_cs%core%deadstem(p)           - veg_cf%core%m_deadstem_to_fire(p)        * dt
         veg_cs%core%deadstem(p)          = veg_cs%core%deadstem(p)           - veg_cf%core%m_deadstem_to_litter_fire(p) * dt
         veg_cs%core%livecroot(p)         = veg_cs%core%livecroot(p)          - veg_cf%core%m_livecroot_to_fire(p)       * dt
         veg_cs%core%livecroot(p)         = veg_cs%core%livecroot(p)          - veg_cf%core%m_livecroot_to_litter_fire(p)* dt
         veg_cs%core%deadcroot(p)         = veg_cs%core%deadcroot(p)          - veg_cf%core%m_deadcroot_to_fire(p)       * dt
         veg_cs%core%deadcroot(p)         = veg_cs%core%deadcroot(p)          - veg_cf%core%m_deadcroot_to_litter_fire(p)* dt

         ! storage pools
         veg_cs%core%leaf_storage(p)      = veg_cs%core%leaf_storage(p)       - veg_cf%core%m_leaf_storage_to_fire(p)            * dt
         veg_cs%core%leaf_storage(p)      = veg_cs%core%leaf_storage(p)       - veg_cf%core%m_leaf_storage_to_litter_fire(p)     * dt
         veg_cs%core%froot_storage(p)     = veg_cs%core%froot_storage(p)      - veg_cf%core%m_froot_storage_to_fire(p)           * dt
         veg_cs%core%froot_storage(p)     = veg_cs%core%froot_storage(p)      - veg_cf%core%m_froot_storage_to_litter_fire(p)    * dt
         veg_cs%core%livestem_storage(p)  = veg_cs%core%livestem_storage(p)   - veg_cf%core%m_livestem_storage_to_fire(p)        * dt
         veg_cs%core%livestem_storage(p)  = veg_cs%core%livestem_storage(p)   - veg_cf%core%m_livestem_storage_to_litter_fire(p) * dt
         veg_cs%core%deadstem_storage(p)  = veg_cs%core%deadstem_storage(p)   - veg_cf%core%m_deadstem_storage_to_fire(p)        * dt
         veg_cs%core%deadstem_storage(p)  = veg_cs%core%deadstem_storage(p)   - veg_cf%core%m_deadstem_storage_to_litter_fire(p) * dt
         veg_cs%core%livecroot_storage(p) = veg_cs%core%livecroot_storage(p)  - veg_cf%core%m_livecroot_storage_to_fire(p)       * dt
         veg_cs%core%livecroot_storage(p) = veg_cs%core%livecroot_storage(p)  - veg_cf%core%m_livecroot_storage_to_litter_fire(p)* dt
         veg_cs%core%deadcroot_storage(p) = veg_cs%core%deadcroot_storage(p)  - veg_cf%core%m_deadcroot_storage_to_fire(p)       * dt
         veg_cs%core%deadcroot_storage(p) = veg_cs%core%deadcroot_storage(p)  - veg_cf%core%m_deadcroot_storage_to_litter_fire(p)* dt
         veg_cs%gresp_storage(p)      = veg_cs%gresp_storage(p)       - veg_cf%m_gresp_storage_to_fire(p)            * dt
         veg_cs%gresp_storage(p)      = veg_cs%gresp_storage(p)       - veg_cf%m_gresp_storage_to_litter_fire(p)     * dt

         ! transfer pools
         veg_cs%core%leaf_xfer(p)         = veg_cs%core%leaf_xfer(p)          - veg_cf%core%m_leaf_xfer_to_fire(p)            * dt
         veg_cs%core%leaf_xfer(p)         = veg_cs%core%leaf_xfer(p)          - veg_cf%core%m_leaf_xfer_to_litter_fire(p)     * dt
         veg_cs%core%froot_xfer(p)        = veg_cs%core%froot_xfer(p)         - veg_cf%core%m_froot_xfer_to_fire(p)           * dt
         veg_cs%core%froot_xfer(p)        = veg_cs%core%froot_xfer(p)         - veg_cf%core%m_froot_xfer_to_litter_fire(p)    * dt
         veg_cs%core%livestem_xfer(p)     = veg_cs%core%livestem_xfer(p)      - veg_cf%core%m_livestem_xfer_to_fire(p)        * dt
         veg_cs%core%livestem_xfer(p)     = veg_cs%core%livestem_xfer(p)      - veg_cf%core%m_livestem_xfer_to_litter_fire(p) * dt
         veg_cs%core%deadstem_xfer(p)     = veg_cs%core%deadstem_xfer(p)      - veg_cf%core%m_deadstem_xfer_to_fire(p)        * dt
         veg_cs%core%deadstem_xfer(p)     = veg_cs%core%deadstem_xfer(p)      - veg_cf%core%m_deadstem_xfer_to_litter_fire(p) * dt
         veg_cs%core%livecroot_xfer(p)    = veg_cs%core%livecroot_xfer(p)     - veg_cf%core%m_livecroot_xfer_to_fire(p)       * dt
         veg_cs%core%livecroot_xfer(p)    = veg_cs%core%livecroot_xfer(p)     - veg_cf%core%m_livecroot_xfer_to_litter_fire(p)* dt
         veg_cs%core%deadcroot_xfer(p)    = veg_cs%core%deadcroot_xfer(p)     - veg_cf%core%m_deadcroot_xfer_to_fire(p)       * dt
         veg_cs%core%deadcroot_xfer(p)    = veg_cs%core%deadcroot_xfer(p)     - veg_cf%core%m_deadcroot_xfer_to_litter_fire(p)* dt
         veg_cs%gresp_xfer(p)         = veg_cs%gresp_xfer(p)          - veg_cf%m_gresp_xfer_to_fire(p)            * dt
         veg_cs%gresp_xfer(p)         = veg_cs%gresp_xfer(p)          - veg_cf%m_gresp_xfer_to_litter_fire(p)     * dt
         veg_cs%core%pool(p)              = veg_cs%core%pool(p)               - veg_cf%core%m_pool_to_fire(p)                 * dt
         veg_cs%core%pool(p)              = veg_cs%core%pool(p)               - veg_cf%core%m_pool_to_litter_fire(p)          * dt


      end do ! end of pft loop

  end subroutine CarbonStateUpdate3

end module CarbonStateUpdate3Mod
