module CarbonStateUpdate2Mod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for carbon state variable update, mortality fluxes.
  !
  ! !USES:
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use abortutils       , only : endrun
  use elm_varpar       , only : nlevdecomp, i_met_lit, i_cel_lit, i_lig_lit, i_cwd
  use pftvarcon        , only : iscft
  use elm_varctl       , only : use_pflotran, pf_cmode
  use VegetationType           , only : veg_pp
  use VegetationType        , only : veg_pp
  use ColumnStateTypes, only : column_carbon_state
  use ColumnFluxTypes, only : column_carbon_flux
  use VegetationStateTypes, only : vegetation_carbon_state
  use VegetationFluxTypes, only : vegetation_carbon_flux
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: CarbonStateUpdate2
  public:: CarbonStateUpdate2h
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CarbonStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
       col_cs, veg_cs, col_cf, veg_cf, dt)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update all the prognostic carbon state
    ! variables affected by gap-phase mortality fluxes
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
    real(r8), intent(in) :: dt     ! radiation time step (seconds)

    !
    ! !LOCAL VARIABLES:
    integer  :: c ,p,j ! indices
    integer  :: fp,fc  ! lake filter indices
    !-----------------------------------------------------------------------



     if ((.not.(use_pflotran .and. pf_cmode))) then
         ! column level carbon fluxes from gap-phase mortality
         do j = 1,nlevdecomp
            ! column loop
            do fc = 1,num_soilc
               c = filter_soilc(fc)

               ! column gap mortality fluxes
               col_cs%core%decomp_pools_vr(c,j,i_met_lit) = &
                    col_cs%core%decomp_pools_vr(c,j,i_met_lit) + col_cf%core%gap_mortality_to_litr_met(c,j) * dt
               col_cs%core%decomp_pools_vr(c,j,i_cel_lit) = &
                    col_cs%core%decomp_pools_vr(c,j,i_cel_lit) + col_cf%core%gap_mortality_to_litr_cel(c,j) * dt
               col_cs%core%decomp_pools_vr(c,j,i_lig_lit) = &
                    col_cs%core%decomp_pools_vr(c,j,i_lig_lit) + col_cf%core%gap_mortality_to_litr_lig(c,j) * dt
               col_cs%core%decomp_pools_vr(c,j,i_cwd) = &
                    col_cs%core%decomp_pools_vr(c,j,i_cwd) + col_cf%core%gap_mortality_to_cwd(c,j) * dt

            end do
         end do
      endif


      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! patch-level carbon fluxes from gap-phase mortality
         ! displayed pools
         veg_cs%core%leaf(p)               = veg_cs%core%leaf(p)              - veg_cf%core%m_leaf_to_litter(p)              * dt
         veg_cs%core%froot(p)              = veg_cs%core%froot(p)             - veg_cf%core%m_froot_to_litter(p)             * dt
         veg_cs%core%livestem(p)           = veg_cs%core%livestem(p)          - veg_cf%core%m_livestem_to_litter(p)          * dt
         veg_cs%core%deadstem(p)           = veg_cs%core%deadstem(p)          - veg_cf%core%m_deadstem_to_litter(p)          * dt
         veg_cs%core%livecroot(p)          = veg_cs%core%livecroot(p)         - veg_cf%core%m_livecroot_to_litter(p)         * dt
         veg_cs%core%deadcroot(p)          = veg_cs%core%deadcroot(p)         - veg_cf%core%m_deadcroot_to_litter(p)         * dt
         ! storage pools
         veg_cs%core%leaf_storage(p)       = veg_cs%core%leaf_storage(p)      - veg_cf%core%m_leaf_storage_to_litter(p)      * dt
         veg_cs%core%froot_storage(p)      = veg_cs%core%froot_storage(p)     - veg_cf%core%m_froot_storage_to_litter(p)     * dt
         veg_cs%core%livestem_storage(p)   = veg_cs%core%livestem_storage(p)  - veg_cf%core%m_livestem_storage_to_litter(p)  * dt
         veg_cs%core%deadstem_storage(p)   = veg_cs%core%deadstem_storage(p)  - veg_cf%core%m_deadstem_storage_to_litter(p)  * dt
         veg_cs%core%livecroot_storage(p)  = veg_cs%core%livecroot_storage(p) - veg_cf%core%m_livecroot_storage_to_litter(p) * dt
         veg_cs%core%deadcroot_storage(p)  = veg_cs%core%deadcroot_storage(p) - veg_cf%core%m_deadcroot_storage_to_litter(p) * dt
         veg_cs%gresp_storage(p)       = veg_cs%gresp_storage(p)      - veg_cf%m_gresp_storage_to_litter(p)      * dt
         veg_cs%core%pool(p)               = veg_cs%core%pool(p)              - veg_cf%core%m_pool_to_litter(p)              * dt

         ! transfer pools
         veg_cs%core%leaf_xfer(p)          = veg_cs%core%leaf_xfer(p)         - veg_cf%core%m_leaf_xfer_to_litter(p)         * dt
         veg_cs%core%froot_xfer(p)         = veg_cs%core%froot_xfer(p)        - veg_cf%core%m_froot_xfer_to_litter(p)        * dt
         veg_cs%core%livestem_xfer(p)      = veg_cs%core%livestem_xfer(p)     - veg_cf%core%m_livestem_xfer_to_litter(p)     * dt
         veg_cs%core%deadstem_xfer(p)      = veg_cs%core%deadstem_xfer(p)     - veg_cf%core%m_deadstem_xfer_to_litter(p)     * dt
         veg_cs%core%livecroot_xfer(p)     = veg_cs%core%livecroot_xfer(p)    - veg_cf%core%m_livecroot_xfer_to_litter(p)    * dt
         veg_cs%core%deadcroot_xfer(p)     = veg_cs%core%deadcroot_xfer(p)    - veg_cf%core%m_deadcroot_xfer_to_litter(p)    * dt
         veg_cs%gresp_xfer(p)          = veg_cs%gresp_xfer(p)         - veg_cf%m_gresp_xfer_to_litter(p)         * dt
      end do ! end of patch loop

  end subroutine CarbonStateUpdate2

  !-----------------------------------------------------------------------
  subroutine CarbonStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
       col_cs, veg_cs, col_cf, veg_cf, dt)
    !
    ! !DESCRIPTION:
    ! Update all the prognostic carbon state
    ! variables affected by harvest mortality fluxes
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
    real(r8),    intent(in)    :: dt        ! radiation time step (seconds)

    !
    ! !LOCAL VARIABLES:
    integer :: c,p,j,k,l ! indices
    integer :: fp,fc     ! lake filter indices
    !-----------------------------------------------------------------------

    associate(                   &
         ivt => veg_pp%itype     & ! Input:  [integer (:)]  pft vegetation type
        )

      ! set time steps

      if (.not.(use_pflotran .and. pf_cmode)) then
         ! column level carbon fluxes from harvest mortality
         do j = 1, nlevdecomp
            ! column loop
            do fc = 1,num_soilc
               c = filter_soilc(fc)

               ! column harvest fluxes
               col_cs%core%decomp_pools_vr(c,j,i_met_lit) = &
                    col_cs%core%decomp_pools_vr(c,j,i_met_lit) + col_cf%core%harvest_to_litr_met(c,j) * dt
               col_cs%core%decomp_pools_vr(c,j,i_cel_lit) = &
                    col_cs%core%decomp_pools_vr(c,j,i_cel_lit) + col_cf%core%harvest_to_litr_cel(c,j) * dt
               col_cs%core%decomp_pools_vr(c,j,i_lig_lit) = &
                    col_cs%core%decomp_pools_vr(c,j,i_lig_lit) + col_cf%core%harvest_to_litr_lig(c,j) * dt
               col_cs%core%decomp_pools_vr(c,j,i_cwd) = &
                    col_cs%core%decomp_pools_vr(c,j,i_cwd) + col_cf%core%harvest_to_cwd(c,j)  * dt

               ! wood to product pools - states updated in WoodProducts()
            end do
         end do

      endif

      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! patch-level carbon fluxes from harvest mortality
         ! displayed pools
         veg_cs%core%leaf(p)               = veg_cs%core%leaf(p)              - veg_cf%core%hrv_leaf_to_litter(p)              * dt
         veg_cs%core%froot(p)              = veg_cs%core%froot(p)             - veg_cf%core%hrv_froot_to_litter(p)             * dt
         veg_cs%core%livestem(p)           = veg_cs%core%livestem(p)          - veg_cf%core%hrv_livestem_to_litter(p)          * dt
         veg_cs%core%deadstem(p)           = veg_cs%core%deadstem(p)          - veg_cf%core%hrv_deadstem_to_prod10(p)         * dt
         veg_cs%core%deadstem(p)           = veg_cs%core%deadstem(p)          - veg_cf%core%hrv_deadstem_to_prod100(p)        * dt
         veg_cs%core%livecroot(p)          = veg_cs%core%livecroot(p)         - veg_cf%core%hrv_livecroot_to_litter(p)         * dt
         veg_cs%core%deadcroot(p)          = veg_cs%core%deadcroot(p)         - veg_cf%core%hrv_deadcroot_to_litter(p)         * dt

         ! crops
         if (iscft(ivt(p))) then ! skip 2 generic crops
             veg_cs%core%livestem(p)       = veg_cs%core%livestem(p)          - veg_cf%core%hrv_livestem_to_prod1(p)          *dt
             veg_cs%core%leaf(p)           = veg_cs%core%leaf(p)              - veg_cf%core%hrv_leaf_to_prod1(p)              *dt
             veg_cs%core%grain(p)          = veg_cs%core%grain(p)             - veg_cf%core%hrv_grain_to_prod1(p)             *dt
         end if

         ! xsmrpool
         veg_cs%xsmrpool(p)            = veg_cs%xsmrpool(p)           - veg_cf%hrv_xsmrpool_to_atm(p)              * dt
         ! storage pools
         veg_cs%core%leaf_storage(p)       = veg_cs%core%leaf_storage(p)      - veg_cf%core%hrv_leaf_storage_to_litter(p)      * dt
         veg_cs%core%froot_storage(p)      = veg_cs%core%froot_storage(p)     - veg_cf%core%hrv_froot_storage_to_litter(p)     * dt
         veg_cs%core%livestem_storage(p)   = veg_cs%core%livestem_storage(p)  - veg_cf%core%hrv_livestem_storage_to_litter(p)  * dt
         veg_cs%core%deadstem_storage(p)   = veg_cs%core%deadstem_storage(p)  - veg_cf%core%hrv_deadstem_storage_to_litter(p)  * dt
         veg_cs%core%livecroot_storage(p)  = veg_cs%core%livecroot_storage(p) - veg_cf%core%hrv_livecroot_storage_to_litter(p) * dt
         veg_cs%core%deadcroot_storage(p)  = veg_cs%core%deadcroot_storage(p) - veg_cf%core%hrv_deadcroot_storage_to_litter(p) * dt
         veg_cs%gresp_storage(p)       = veg_cs%gresp_storage(p)      - veg_cf%hrv_gresp_storage_to_litter(p)      * dt
         veg_cs%core%pool(p)               = veg_cs%core%pool(p)              - veg_cf%core%hrv_pool_to_litter(p)              * dt

         ! transfer pools
         veg_cs%core%leaf_xfer(p)          = veg_cs%core%leaf_xfer(p)         - veg_cf%core%hrv_leaf_xfer_to_litter(p)         * dt
         veg_cs%core%froot_xfer(p)         = veg_cs%core%froot_xfer(p)        - veg_cf%core%hrv_froot_xfer_to_litter(p)        * dt
         veg_cs%core%livestem_xfer(p)      = veg_cs%core%livestem_xfer(p)     - veg_cf%core%hrv_livestem_xfer_to_litter(p)     * dt
         veg_cs%core%deadstem_xfer(p)      = veg_cs%core%deadstem_xfer(p)     - veg_cf%core%hrv_deadstem_xfer_to_litter(p)     * dt
         veg_cs%core%livecroot_xfer(p)     = veg_cs%core%livecroot_xfer(p)    - veg_cf%core%hrv_livecroot_xfer_to_litter(p)    * dt
         veg_cs%core%deadcroot_xfer(p)     = veg_cs%core%deadcroot_xfer(p)    - veg_cf%core%hrv_deadcroot_xfer_to_litter(p)    * dt
         veg_cs%gresp_xfer(p)          = veg_cs%gresp_xfer(p)         - veg_cf%hrv_gresp_xfer_to_litter(p)         * dt

      end do ! end of patch loop

    end associate

  end subroutine CarbonStateUpdate2h

end module CarbonStateUpdate2Mod
