module PhosphorusStateUpdate2Mod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for phosphorus state variable update, mortality fluxes.
  ! X.YANG
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use elm_varpar          , only : nlevsoi, nlevdecomp
  use elm_varpar          , only : i_met_lit, i_cel_lit, i_lig_lit, i_cwd
  use elm_varctl          , only : iulog
  !use PhosphorusStateType , only : phosphorusstate_type
  !use PhosphorusFLuxType  , only : phosphorusflux_type
  use VegetationType           , only : veg_pp
  use pftvarcon           , only : iscft
  ! bgc interface & pflotran:
  use elm_varctl          , only : use_pflotran, pf_cmode
  use ColumnDataType      , only : col_ps, col_pf
  use VegetationDataType  , only : veg_ps, veg_pf
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: PhosphorusStateUpdate2
  public:: PhosphorusStateUpdate2h
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine PhosphorusStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, dt)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update all the prognostic phosporus state
    ! variables affected by gap-phase mortality fluxes
    ! NOTE - associate statements have been removed where there are
    ! no science equations. This increases readability and maintainability
    !
    ! !ARGUMENTS:
      !$acc routine seq
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                  , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:) ! filter for soil patches
    real(r8), intent(in) :: dt      ! radiation time step (seconds)

    !
    ! !LOCAL VARIABLES:
    integer  :: c,p,j,l ! indices
    integer  :: fp,fc   ! lake filter indices
    !-----------------------------------------------------------------------


      !------------------------------------------------------------------
      ! if coupled with pflotran, the following updates are NOT needed
!      if (.not.(use_pflotran .and. pf_cmode)) then
      !------------------------------------------------------------------

      ! column-level phosporus fluxes from gap-phase mortality
      do j = 1, nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)

            col_ps%core%decomp_pools_vr(c,j,i_met_lit) = &
                 col_ps%core%decomp_pools_vr(c,j,i_met_lit) + col_pf%core%gap_mortality_to_litr_met(c,j) * dt
            col_ps%core%decomp_pools_vr(c,j,i_cel_lit) = &
                 col_ps%core%decomp_pools_vr(c,j,i_cel_lit) + col_pf%core%gap_mortality_to_litr_cel(c,j) * dt
            col_ps%core%decomp_pools_vr(c,j,i_lig_lit) = &
                 col_ps%core%decomp_pools_vr(c,j,i_lig_lit) + col_pf%core%gap_mortality_to_litr_lig(c,j) * dt
            col_ps%core%decomp_pools_vr(c,j,i_cwd)     = &
                 col_ps%core%decomp_pools_vr(c,j,i_cwd)     + col_pf%core%gap_mortality_to_cwd(c,j)       * dt
         end do
      end do
      !------------------------------------------------------------------

      ! patch -level phosporus fluxes from gap-phase mortality

      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! displayed pools
         veg_ps%core%leaf(p)              =  veg_ps%core%leaf(p)      - veg_pf%core%m_leaf_to_litter(p)      * dt
         veg_ps%core%froot(p)             =  veg_ps%core%froot(p)     - veg_pf%core%m_froot_to_litter(p)     * dt
         veg_ps%core%livestem(p)          =  veg_ps%core%livestem(p)  - veg_pf%core%m_livestem_to_litter(p)  * dt
         veg_ps%core%deadstem(p)          =  veg_ps%core%deadstem(p)  - veg_pf%core%m_deadstem_to_litter(p)  * dt
         veg_ps%core%livecroot(p)         =  veg_ps%core%livecroot(p) - veg_pf%core%m_livecroot_to_litter(p) * dt
         veg_ps%core%deadcroot(p)         =  veg_ps%core%deadcroot(p) - veg_pf%core%m_deadcroot_to_litter(p) * dt
         veg_ps%retrans%pool(p)           =  veg_ps%retrans%pool(p)   - veg_pf%retrans%m_retrans_to_litter(p)   * dt
         veg_ps%core%pool(p)              =  veg_ps%core%pool(p)      - veg_pf%core%m_pool_to_litter(p)   * dt

         ! storage pools
         veg_ps%core%leaf_storage(p)      =  veg_ps%core%leaf_storage(p)      - veg_pf%core%m_leaf_storage_to_litter(p)      * dt
         veg_ps%core%froot_storage(p)     =  veg_ps%core%froot_storage(p)     - veg_pf%core%m_froot_storage_to_litter(p)     * dt
         veg_ps%core%livestem_storage(p)  =  veg_ps%core%livestem_storage(p)  - veg_pf%core%m_livestem_storage_to_litter(p)  * dt
         veg_ps%core%deadstem_storage(p)  =  veg_ps%core%deadstem_storage(p)  - veg_pf%core%m_deadstem_storage_to_litter(p)  * dt
         veg_ps%core%livecroot_storage(p) =  veg_ps%core%livecroot_storage(p) - veg_pf%core%m_livecroot_storage_to_litter(p) * dt
         veg_ps%core%deadcroot_storage(p) =  veg_ps%core%deadcroot_storage(p) - veg_pf%core%m_deadcroot_storage_to_litter(p) * dt

         ! transfer pools
         veg_ps%core%leaf_xfer(p)         =  veg_ps%core%leaf_xfer(p)      - veg_pf%core%m_leaf_xfer_to_litter(p)      * dt
         veg_ps%core%froot_xfer(p)        =  veg_ps%core%froot_xfer(p)     - veg_pf%core%m_froot_xfer_to_litter(p)     * dt
         veg_ps%core%livestem_xfer(p)     =  veg_ps%core%livestem_xfer(p)  - veg_pf%core%m_livestem_xfer_to_litter(p)  * dt
         veg_ps%core%deadstem_xfer(p)     =  veg_ps%core%deadstem_xfer(p)  - veg_pf%core%m_deadstem_xfer_to_litter(p)  * dt
         veg_ps%core%livecroot_xfer(p)    =  veg_ps%core%livecroot_xfer(p) - veg_pf%core%m_livecroot_xfer_to_litter(p) * dt
         veg_ps%core%deadcroot_xfer(p)    =  veg_ps%core%deadcroot_xfer(p) - veg_pf%core%m_deadcroot_xfer_to_litter(p) * dt

      end do

  end subroutine PhosphorusStateUpdate2

  !-----------------------------------------------------------------------
  subroutine PhosphorusStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, dt)
    !
    ! !DESCRIPTION:
    ! Update all the prognostic phosphorus state
    ! variables affected by harvest mortality fluxes
    ! NOTE - associate statements have been removed where there are
    ! no science equations. This increases readability and maintainability
    !
    ! !ARGUMENTS:
      !$acc routine seq
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                  , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:) ! filter for soil patches
    !type(phosphorusflux_type)  , intent(in)    :: phosphorusflux_vars
    !type(phosphorusstate_type) , intent(inout) :: phosphorusstate_vars
    real(r8),     intent(in)     :: dt      ! radiation time step (seconds)

    !
    ! !LOCAL VARIABLES:
    integer :: c,p,j,l ! indices
    integer :: fp,fc   ! lake filter indices
    !-----------------------------------------------------------------------

    associate(                      &
         ivt => veg_pp%itype            & ! Input:  [integer  (:) ]  pft vegetation type

         )

      !------------------------------------------------------------------
      ! if coupled with pflotran, the following updates are NOT needed
      if (.not.(use_pflotran .and. pf_cmode)) then
      !------------------------------------------------------------------

      ! column-level phosporus fluxes from harvest mortality

      do j = 1,nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            col_ps%core%decomp_pools_vr(c,j,i_met_lit) = &
                 col_ps%core%decomp_pools_vr(c,j,i_met_lit) + col_pf%core%harvest_to_litr_met(c,j) * dt
            col_ps%core%decomp_pools_vr(c,j,i_cel_lit) = &
                 col_ps%core%decomp_pools_vr(c,j,i_cel_lit) + col_pf%core%harvest_to_litr_cel(c,j) * dt
            col_ps%core%decomp_pools_vr(c,j,i_lig_lit) = &
                 col_ps%core%decomp_pools_vr(c,j,i_lig_lit) + col_pf%core%harvest_to_litr_lig(c,j) * dt
            col_ps%core%decomp_pools_vr(c,j,i_cwd)     = &
                 col_ps%core%decomp_pools_vr(c,j,i_cwd)     + col_pf%core%harvest_to_cwd(c,j)       * dt
         end do
      end do
      endif ! if (.not.(use_pflotran .and. pf_cmode))
      !------------------------------------------------------------------

      ! patch-level phosporus fluxes from harvest mortality

      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! displayed pools
         veg_ps%core%leaf(p)      = veg_ps%core%leaf(p)      - veg_pf%core%hrv_leaf_to_litter(p)      * dt
         veg_ps%core%froot(p)     = veg_ps%core%froot(p)     - veg_pf%core%hrv_froot_to_litter(p)     * dt
         veg_ps%core%livestem(p)  = veg_ps%core%livestem(p)  - veg_pf%core%hrv_livestem_to_litter(p)  * dt
         veg_ps%core%deadstem(p)  = veg_ps%core%deadstem(p)  - veg_pf%core%hrv_deadstem_to_prod10(p) * dt
         veg_ps%core%deadstem(p)  = veg_ps%core%deadstem(p)  - veg_pf%core%hrv_deadstem_to_prod100(p)* dt
         veg_ps%core%livecroot(p) = veg_ps%core%livecroot(p) - veg_pf%core%hrv_livecroot_to_litter(p) * dt
         veg_ps%core%deadcroot(p) = veg_ps%core%deadcroot(p) - veg_pf%core%hrv_deadcroot_to_litter(p) * dt
         veg_ps%retrans%pool(p)   = veg_ps%retrans%pool(p)   - veg_pf%retrans%hrv_retrans_to_litter(p)   * dt
         veg_ps%core%pool(p)      = veg_ps%core%pool(p)      - veg_pf%core%hrv_pool_to_litter(p)      * dt

       if (iscft(ivt(p))) then ! skip 2 generic crops
           veg_ps%core%livestem(p)= veg_ps%core%livestem(p)  - veg_pf%core%hrv_livestem_to_prod1(p)  * dt
           veg_ps%core%leaf(p)    = veg_ps%core%leaf(p)      - veg_pf%core%hrv_leaf_to_prod1(p)      * dt
           veg_ps%core%grain(p)   = veg_ps%core%grain(p)     - veg_pf%core%hrv_grain_to_prod1(p)     * dt
       end if

         ! storage pools
         veg_ps%core%leaf_storage(p)      = veg_ps%core%leaf_storage(p)      - veg_pf%core%hrv_leaf_storage_to_litter(p)      * dt
         veg_ps%core%froot_storage(p)     = veg_ps%core%froot_storage(p)     - veg_pf%core%hrv_froot_storage_to_litter(p)     * dt
         veg_ps%core%livestem_storage(p)  = veg_ps%core%livestem_storage(p)  - veg_pf%core%hrv_livestem_storage_to_litter(p)  * dt
         veg_ps%core%deadstem_storage(p)  = veg_ps%core%deadstem_storage(p)  - veg_pf%core%hrv_deadstem_storage_to_litter(p)  * dt
         veg_ps%core%livecroot_storage(p) = veg_ps%core%livecroot_storage(p) - veg_pf%core%hrv_livecroot_storage_to_litter(p) * dt
         veg_ps%core%deadcroot_storage(p) = veg_ps%core%deadcroot_storage(p) - veg_pf%core%hrv_deadcroot_storage_to_litter(p) * dt

         ! transfer pools
         veg_ps%core%leaf_xfer(p)      = veg_ps%core%leaf_xfer(p)      - veg_pf%core%hrv_leaf_xfer_to_litter(p)      *dt
         veg_ps%core%froot_xfer(p)     = veg_ps%core%froot_xfer(p)     - veg_pf%core%hrv_froot_xfer_to_litter(p)     *dt
         veg_ps%core%livestem_xfer(p)  = veg_ps%core%livestem_xfer(p)  - veg_pf%core%hrv_livestem_xfer_to_litter(p)  *dt
         veg_ps%core%deadstem_xfer(p)  = veg_ps%core%deadstem_xfer(p)  - veg_pf%core%hrv_deadstem_xfer_to_litter(p)  *dt
         veg_ps%core%livecroot_xfer(p) = veg_ps%core%livecroot_xfer(p) - veg_pf%core%hrv_livecroot_xfer_to_litter(p) *dt
         veg_ps%core%deadcroot_xfer(p) = veg_ps%core%deadcroot_xfer(p) - veg_pf%core%hrv_deadcroot_xfer_to_litter(p) *dt

      end do

    end associate

  end subroutine PhosphorusStateUpdate2h

end module PhosphorusStateUpdate2Mod
