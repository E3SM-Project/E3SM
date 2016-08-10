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
  use PatchType           , only : pft
  use pftvarcon           , only : npcropmin
  !! bgc interface & pflotran:
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

      do j = 1, nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)

            ps%decomp_ppools_vr_col(c,j,i_met_lit) = &
                 ps%decomp_ppools_vr_col(c,j,i_met_lit) + pf%gap_mortality_p_to_litr_met_p_col(c,j) * dt
            ps%decomp_ppools_vr_col(c,j,i_cel_lit) = &
                 ps%decomp_ppools_vr_col(c,j,i_cel_lit) + pf%gap_mortality_p_to_litr_cel_p_col(c,j) * dt
            ps%decomp_ppools_vr_col(c,j,i_lig_lit) = &
                 ps%decomp_ppools_vr_col(c,j,i_lig_lit) + pf%gap_mortality_p_to_litr_lig_p_col(c,j) * dt
            ps%decomp_ppools_vr_col(c,j,i_cwd)     = &
                 ps%decomp_ppools_vr_col(c,j,i_cwd)     + pf%gap_mortality_p_to_cwdp_col(c,j)       * dt
         end do
      end do
!      endif ! if (.not.(use_pflotran .and. pf_cmode))
      !------------------------------------------------------------------

      ! patch -level phosporus fluxes from gap-phase mortality

      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! displayed pools
         ps%leafp_patch(p)              =  ps%leafp_patch(p)      - pf%m_leafp_to_litter_patch(p)      * dt
         ps%frootp_patch(p)             =  ps%frootp_patch(p)     - pf%m_frootp_to_litter_patch(p)     * dt
         ps%livestemp_patch(p)          =  ps%livestemp_patch(p)  - pf%m_livestemp_to_litter_patch(p)  * dt
         ps%deadstemp_patch(p)          =  ps%deadstemp_patch(p)  - pf%m_deadstemp_to_litter_patch(p)  * dt
         ps%livecrootp_patch(p)         =  ps%livecrootp_patch(p) - pf%m_livecrootp_to_litter_patch(p) * dt
         ps%deadcrootp_patch(p)         =  ps%deadcrootp_patch(p) - pf%m_deadcrootp_to_litter_patch(p) * dt
         ps%retransp_patch(p)           =  ps%retransp_patch(p)   - pf%m_retransp_to_litter_patch(p)   * dt

         ! storage pools
         ps%leafp_storage_patch(p)      =  ps%leafp_storage_patch(p)      - pf%m_leafp_storage_to_litter_patch(p)      * dt
         ps%frootp_storage_patch(p)     =  ps%frootp_storage_patch(p)     - pf%m_frootp_storage_to_litter_patch(p)     * dt
         ps%livestemp_storage_patch(p)  =  ps%livestemp_storage_patch(p)  - pf%m_livestemp_storage_to_litter_patch(p)  * dt
         ps%deadstemp_storage_patch(p)  =  ps%deadstemp_storage_patch(p)  - pf%m_deadstemp_storage_to_litter_patch(p)  * dt
         ps%livecrootp_storage_patch(p) =  ps%livecrootp_storage_patch(p) - pf%m_livecrootp_storage_to_litter_patch(p) * dt
         ps%deadcrootp_storage_patch(p) =  ps%deadcrootp_storage_patch(p) - pf%m_deadcrootp_storage_to_litter_patch(p) * dt

         ! transfer pools
         ps%leafp_xfer_patch(p)         =  ps%leafp_xfer_patch(p)      - pf%m_leafp_xfer_to_litter_patch(p)      * dt
         ps%frootp_xfer_patch(p)        =  ps%frootp_xfer_patch(p)     - pf%m_frootp_xfer_to_litter_patch(p)     * dt
         ps%livestemp_xfer_patch(p)     =  ps%livestemp_xfer_patch(p)  - pf%m_livestemp_xfer_to_litter_patch(p)  * dt
         ps%deadstemp_xfer_patch(p)     =  ps%deadstemp_xfer_patch(p)  - pf%m_deadstemp_xfer_to_litter_patch(p)  * dt
         ps%livecrootp_xfer_patch(p)    =  ps%livecrootp_xfer_patch(p) - pf%m_livecrootp_xfer_to_litter_patch(p) * dt
         ps%deadcrootp_xfer_patch(p)    =  ps%deadcrootp_xfer_patch(p) - pf%m_deadcrootp_xfer_to_litter_patch(p) * dt

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
         ivt => pft%itype           , & ! Input:  [integer  (:) ]  pft vegetation type
         pf => phosphorusflux_vars  , &
         ps => phosphorusstate_vars   &
         )

      ! set time steps
      dt = real( get_step_size(), r8 )

      !------------------------------------------------------------------
      ! if coupled with pflotran, the following updates are NOT needed
      if (.not.(use_pflotran .and. pf_cmode)) then
      !------------------------------------------------------------------

      ! column-level phosporus fluxes from harvest mortality

      do j = 1,nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            ps%decomp_ppools_vr_col(c,j,i_met_lit) = &
                 ps%decomp_ppools_vr_col(c,j,i_met_lit) + pf%harvest_p_to_litr_met_p_col(c,j) * dt
            ps%decomp_ppools_vr_col(c,j,i_cel_lit) = &
                 ps%decomp_ppools_vr_col(c,j,i_cel_lit) + pf%harvest_p_to_litr_cel_p_col(c,j) * dt
            ps%decomp_ppools_vr_col(c,j,i_lig_lit) = &
                 ps%decomp_ppools_vr_col(c,j,i_lig_lit) + pf%harvest_p_to_litr_lig_p_col(c,j) * dt
            ps%decomp_ppools_vr_col(c,j,i_cwd)     = &
                 ps%decomp_ppools_vr_col(c,j,i_cwd)     + pf%harvest_p_to_cwdp_col(c,j)       * dt
         end do
      end do
      endif ! if (.not.(use_pflotran .and. pf_cmode))
      !------------------------------------------------------------------

      ! patch-level phosporus fluxes from harvest mortality

      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! displayed pools
         ps%leafp_patch(p)      = ps%leafp_patch(p)      - pf%hrv_leafp_to_litter_patch(p)      * dt
         ps%frootp_patch(p)     = ps%frootp_patch(p)     - pf%hrv_frootp_to_litter_patch(p)     * dt
         ps%livestemp_patch(p)  = ps%livestemp_patch(p)  - pf%hrv_livestemp_to_litter_patch(p)  * dt
         ps%deadstemp_patch(p)  = ps%deadstemp_patch(p)  - pf%hrv_deadstemp_to_prod10p_patch(p) * dt
         ps%deadstemp_patch(p)  = ps%deadstemp_patch(p)  - pf%hrv_deadstemp_to_prod100p_patch(p)* dt
         ps%livecrootp_patch(p) = ps%livecrootp_patch(p) - pf%hrv_livecrootp_to_litter_patch(p) * dt
         ps%deadcrootp_patch(p) = ps%deadcrootp_patch(p) - pf%hrv_deadcrootp_to_litter_patch(p) * dt
         ps%retransp_patch(p)   = ps%retransp_patch(p)   - pf%hrv_retransp_to_litter_patch(p)   * dt

       if (ivt(p) >= npcropmin) then ! skip 2 generic crops
           ps%livestemp_patch(p)= ps%livestemp_patch(p)  - pf%hrv_livestemp_to_prod1p_patch(p)  * dt
           ps%leafp_patch(p)    = ps%leafp_patch(p)      - pf%hrv_leafp_to_prod1p_patch(p)      * dt
           ps%grainp_patch(p)   = ps%grainp_patch(p)     - pf%hrv_grainp_to_prod1p_patch(p)     * dt
       end if

         ! storage pools
         ps%leafp_storage_patch(p)      = ps%leafp_storage_patch(p)      - pf%hrv_leafp_storage_to_litter_patch(p)      * dt
         ps%frootp_storage_patch(p)     = ps%frootp_storage_patch(p)     - pf%hrv_frootp_storage_to_litter_patch(p)     * dt
         ps%livestemp_storage_patch(p)  = ps%livestemp_storage_patch(p)  - pf%hrv_livestemp_storage_to_litter_patch(p)  * dt
         ps%deadstemp_storage_patch(p)  = ps%deadstemp_storage_patch(p)  - pf%hrv_deadstemp_storage_to_litter_patch(p)  * dt
         ps%livecrootp_storage_patch(p) = ps%livecrootp_storage_patch(p) - pf%hrv_livecrootp_storage_to_litter_patch(p) * dt
         ps%deadcrootp_storage_patch(p) = ps%deadcrootp_storage_patch(p) - pf%hrv_deadcrootp_storage_to_litter_patch(p) * dt

         ! transfer pools
         ps%leafp_xfer_patch(p)      = ps%leafp_xfer_patch(p)      - pf%hrv_leafp_xfer_to_litter_patch(p)      *dt
         ps%frootp_xfer_patch(p)     = ps%frootp_xfer_patch(p)     - pf%hrv_frootp_xfer_to_litter_patch(p)     *dt
         ps%livestemp_xfer_patch(p)  = ps%livestemp_xfer_patch(p)  - pf%hrv_livestemp_xfer_to_litter_patch(p)  *dt
         ps%deadstemp_xfer_patch(p)  = ps%deadstemp_xfer_patch(p)  - pf%hrv_deadstemp_xfer_to_litter_patch(p)  *dt
         ps%livecrootp_xfer_patch(p) = ps%livecrootp_xfer_patch(p) - pf%hrv_livecrootp_xfer_to_litter_patch(p) *dt
         ps%deadcrootp_xfer_patch(p) = ps%deadcrootp_xfer_patch(p) - pf%hrv_deadcrootp_xfer_to_litter_patch(p) *dt

      end do

    end associate

  end subroutine PStateUpdate2h

end module PStateUpdate2Mod
