module PStateUpdate2Mod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for phosphorus state variable update, mortality fluxes.
  ! X.YANG
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use clm_time_manager    , only : get_step_size
  use clm_varpar          , only : nlevsoi, nlevdecomp
  use clm_varpar          , only : i_met_lit, i_cel_lit, i_lig_lit, i_cwd
  use clm_varctl          , only : iulog
  use PhosphorusStateType , only : phosphorusstate_type
  use PhosphorusFLuxType  , only : phosphorusflux_type
  use VegetationType           , only : veg_pp
  use pftvarcon           , only : npcropmin
  use tracer_varcon       , only : is_active_betr_bgc
  ! bgc interface & pflotran:
  use clm_varctl          , only : use_pflotran, pf_cmode
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: PStateUpdate2
  public:: PStateUpdate2h
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine PStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
       phosphorusflux_vars, phosphorusstate_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update all the prognostic phosporus state
    ! variables affected by gap-phase mortality fluxes
    ! NOTE - associate statements have been removed where there are
    ! no science equations. This increases readability and maintainability
    !
    ! !ARGUMENTS:
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                  , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(phosphorusflux_type)  , intent(in)    :: phosphorusflux_vars
    type(phosphorusstate_type) , intent(inout) :: phosphorusstate_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: c,p,j,l ! indices
    integer  :: fp,fc   ! lake filter indices
    real(r8) :: dt      ! radiation time step (seconds)
    !-----------------------------------------------------------------------

    associate(                      & 
         pf => phosphorusflux_vars  , &
         ps => phosphorusstate_vars   &
         )

      ! set time steps
      dt = real( get_step_size(), r8 )

      !------------------------------------------------------------------
      ! if coupled with pflotran, the following updates are NOT needed
!      if (.not.(use_pflotran .and. pf_cmode)) then
      !------------------------------------------------------------------

      ! column-level phosporus fluxes from gap-phase mortality
      if (.not. is_active_betr_bgc) then
      do j = 1, nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)

            ps%decomp_pools_vr_col(c,j,i_met_lit) = &
                 ps%decomp_pools_vr_col(c,j,i_met_lit) + pf%gap_mortality_to_litr_met_col(c,j) * dt
            ps%decomp_pools_vr_col(c,j,i_cel_lit) = &
                 ps%decomp_pools_vr_col(c,j,i_cel_lit) + pf%gap_mortality_to_litr_cel_col(c,j) * dt
            ps%decomp_pools_vr_col(c,j,i_lig_lit) = &
                 ps%decomp_pools_vr_col(c,j,i_lig_lit) + pf%gap_mortality_to_litr_lig_col(c,j) * dt
            ps%decomp_pools_vr_col(c,j,i_cwd)     = &
                 ps%decomp_pools_vr_col(c,j,i_cwd)     + pf%gap_mortality_to_cwd_col(c,j)       * dt
         end do
      end do
      endif ! if (.not.is_active_betr_bgc))
      !------------------------------------------------------------------

      ! patch -level phosporus fluxes from gap-phase mortality

      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! displayed pools
         ps%leaf_patch(p)              =  ps%leaf_patch(p)      - pf%m_leaf_to_litter_patch(p)      * dt
         ps%froot_patch(p)             =  ps%froot_patch(p)     - pf%m_froot_to_litter_patch(p)     * dt
         ps%livestem_patch(p)          =  ps%livestem_patch(p)  - pf%m_livestem_to_litter_patch(p)  * dt
         ps%deadstem_patch(p)          =  ps%deadstem_patch(p)  - pf%m_deadstem_to_litter_patch(p)  * dt
         ps%livecroot_patch(p)         =  ps%livecroot_patch(p) - pf%m_livecroot_to_litter_patch(p) * dt
         ps%deadcroot_patch(p)         =  ps%deadcroot_patch(p) - pf%m_deadcroot_to_litter_patch(p) * dt
         ps%retransp_patch(p)           =  ps%retransp_patch(p)   - pf%m_retransp_to_litter_patch(p)   * dt
         ps%pool_patch(p)              =  ps%pool_patch(p)      - pf%m_pool_to_litter_patch(p)   * dt

         ! storage pools
         ps%leaf_storage_patch(p)      =  ps%leaf_storage_patch(p)      - pf%m_leaf_storage_to_litter_patch(p)      * dt
         ps%froot_storage_patch(p)     =  ps%froot_storage_patch(p)     - pf%m_froot_storage_to_litter_patch(p)     * dt
         ps%livestem_storage_patch(p)  =  ps%livestem_storage_patch(p)  - pf%m_livestem_storage_to_litter_patch(p)  * dt
         ps%deadstem_storage_patch(p)  =  ps%deadstem_storage_patch(p)  - pf%m_deadstem_storage_to_litter_patch(p)  * dt
         ps%livecroot_storage_patch(p) =  ps%livecroot_storage_patch(p) - pf%m_livecroot_storage_to_litter_patch(p) * dt
         ps%deadcroot_storage_patch(p) =  ps%deadcroot_storage_patch(p) - pf%m_deadcroot_storage_to_litter_patch(p) * dt

         ! transfer pools
         ps%leaf_xfer_patch(p)         =  ps%leaf_xfer_patch(p)      - pf%m_leaf_xfer_to_litter_patch(p)      * dt
         ps%froot_xfer_patch(p)        =  ps%froot_xfer_patch(p)     - pf%m_froot_xfer_to_litter_patch(p)     * dt
         ps%livestem_xfer_patch(p)     =  ps%livestem_xfer_patch(p)  - pf%m_livestem_xfer_to_litter_patch(p)  * dt
         ps%deadstem_xfer_patch(p)     =  ps%deadstem_xfer_patch(p)  - pf%m_deadstem_xfer_to_litter_patch(p)  * dt
         ps%livecroot_xfer_patch(p)    =  ps%livecroot_xfer_patch(p) - pf%m_livecroot_xfer_to_litter_patch(p) * dt
         ps%deadcroot_xfer_patch(p)    =  ps%deadcroot_xfer_patch(p) - pf%m_deadcroot_xfer_to_litter_patch(p) * dt

      end do

    end associate

  end subroutine PStateUpdate2

  !-----------------------------------------------------------------------
  subroutine PStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
       phosphorusflux_vars, phosphorusstate_vars)
    !
    ! !DESCRIPTION:
    ! Update all the prognostic phosphorus state
    ! variables affected by harvest mortality fluxes
    ! NOTE - associate statements have been removed where there are
    ! no science equations. This increases readability and maintainability
    !
    ! !ARGUMENTS:
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                  , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(phosphorusflux_type)  , intent(in)    :: phosphorusflux_vars
    type(phosphorusstate_type) , intent(inout) :: phosphorusstate_vars
    !
    ! !LOCAL VARIABLES:
    integer :: c,p,j,l ! indices
    integer :: fp,fc   ! lake filter indices
    real(r8):: dt      ! radiation time step (seconds)
    !-----------------------------------------------------------------------

    associate(                      & 
         ivt => veg_pp%itype           , & ! Input:  [integer  (:) ]  pft vegetation type
         pf => phosphorusflux_vars  , &
         ps => phosphorusstate_vars   &
         )

      ! set time steps
      dt = real( get_step_size(), r8 )

      !------------------------------------------------------------------
      ! if coupled with pflotran, the following updates are NOT needed
      if ((.not. is_active_betr_bgc) .and. .not.(use_pflotran .and. pf_cmode)) then
      !------------------------------------------------------------------

      ! column-level phosporus fluxes from harvest mortality

      do j = 1,nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            ps%decomp_pools_vr_col(c,j,i_met_lit) = &
                 ps%decomp_pools_vr_col(c,j,i_met_lit) + pf%harvest_to_litr_met_col(c,j) * dt
            ps%decomp_pools_vr_col(c,j,i_cel_lit) = &
                 ps%decomp_pools_vr_col(c,j,i_cel_lit) + pf%harvest_to_litr_cel_col(c,j) * dt
            ps%decomp_pools_vr_col(c,j,i_lig_lit) = &
                 ps%decomp_pools_vr_col(c,j,i_lig_lit) + pf%harvest_to_litr_lig_col(c,j) * dt
            ps%decomp_pools_vr_col(c,j,i_cwd)     = &
                 ps%decomp_pools_vr_col(c,j,i_cwd)     + pf%harvest_to_cwd_col(c,j)       * dt
         end do
      end do
      endif ! if (.not.(use_pflotran .and. pf_cmode))
      !------------------------------------------------------------------

      ! patch-level phosporus fluxes from harvest mortality

      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! displayed pools
         ps%leaf_patch(p)      = ps%leaf_patch(p)      - pf%hrv_leaf_to_litter_patch(p)      * dt
         ps%froot_patch(p)     = ps%froot_patch(p)     - pf%hrv_froot_to_litter_patch(p)     * dt
         ps%livestem_patch(p)  = ps%livestem_patch(p)  - pf%hrv_livestem_to_litter_patch(p)  * dt
         ps%deadstem_patch(p)  = ps%deadstem_patch(p)  - pf%hrv_deadstem_to_prod10_patch(p) * dt
         ps%deadstem_patch(p)  = ps%deadstem_patch(p)  - pf%hrv_deadstem_to_prod100_patch(p)* dt
         ps%livecroot_patch(p) = ps%livecroot_patch(p) - pf%hrv_livecroot_to_litter_patch(p) * dt
         ps%deadcroot_patch(p) = ps%deadcroot_patch(p) - pf%hrv_deadcroot_to_litter_patch(p) * dt
         ps%retransp_patch(p)   = ps%retransp_patch(p)   - pf%hrv_retransp_to_litter_patch(p)   * dt
         ps%pool_patch(p)      = ps%pool_patch(p)      - pf%hrv_pool_to_litter_patch(p)      * dt

       if (ivt(p) >= npcropmin) then ! skip 2 generic crops
           ps%livestem_patch(p)= ps%livestem_patch(p)  - pf%hrv_livestem_to_prod1_patch(p)  * dt
           ps%leaf_patch(p)    = ps%leaf_patch(p)      - pf%hrv_leaf_to_prod1_patch(p)      * dt
           ps%grain_patch(p)   = ps%grain_patch(p)     - pf%hrv_grain_to_prod1_patch(p)     * dt
       end if

         ! storage pools
         ps%leaf_storage_patch(p)      = ps%leaf_storage_patch(p)      - pf%hrv_leaf_storage_to_litter_patch(p)      * dt
         ps%froot_storage_patch(p)     = ps%froot_storage_patch(p)     - pf%hrv_froot_storage_to_litter_patch(p)     * dt
         ps%livestem_storage_patch(p)  = ps%livestem_storage_patch(p)  - pf%hrv_livestem_storage_to_litter_patch(p)  * dt
         ps%deadstem_storage_patch(p)  = ps%deadstem_storage_patch(p)  - pf%hrv_deadstem_storage_to_litter_patch(p)  * dt
         ps%livecroot_storage_patch(p) = ps%livecroot_storage_patch(p) - pf%hrv_livecroot_storage_to_litter_patch(p) * dt
         ps%deadcroot_storage_patch(p) = ps%deadcroot_storage_patch(p) - pf%hrv_deadcroot_storage_to_litter_patch(p) * dt

         ! transfer pools
         ps%leaf_xfer_patch(p)      = ps%leaf_xfer_patch(p)      - pf%hrv_leaf_xfer_to_litter_patch(p)      *dt
         ps%froot_xfer_patch(p)     = ps%froot_xfer_patch(p)     - pf%hrv_froot_xfer_to_litter_patch(p)     *dt
         ps%livestem_xfer_patch(p)  = ps%livestem_xfer_patch(p)  - pf%hrv_livestem_xfer_to_litter_patch(p)  *dt
         ps%deadstem_xfer_patch(p)  = ps%deadstem_xfer_patch(p)  - pf%hrv_deadstem_xfer_to_litter_patch(p)  *dt
         ps%livecroot_xfer_patch(p) = ps%livecroot_xfer_patch(p) - pf%hrv_livecroot_xfer_to_litter_patch(p) *dt
         ps%deadcroot_xfer_patch(p) = ps%deadcroot_xfer_patch(p) - pf%hrv_deadcroot_xfer_to_litter_patch(p) *dt

      end do

    end associate

  end subroutine PStateUpdate2h

end module PStateUpdate2Mod
