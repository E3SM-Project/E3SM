module PhosphorusStateUpdate2Mod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for phosphorus state variable update, mortality fluxes.
  ! X.YANG
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use clm_time_manager    , only : get_step_size
  use elm_varpar          , only : nlevsoi, nlevdecomp
  use elm_varpar          , only : i_met_lit, i_cel_lit, i_lig_lit, i_cwd
  use clm_varctl          , only : iulog
  use PhosphorusStateType , only : phosphorusstate_type
  use PhosphorusFLuxType  , only : phosphorusflux_type
  use VegetationType           , only : veg_pp
  use pftvarcon           , only : npcropmin
  use tracer_varcon       , only : is_active_betr_bgc
  ! bgc interface & pflotran:
  use clm_varctl          , only : use_pflotran, pf_cmode
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
  subroutine PhosphorusStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
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

            col_ps%decomp_ppools_vr(c,j,i_met_lit) = &
                 col_ps%decomp_ppools_vr(c,j,i_met_lit) + col_pf%gap_mortality_p_to_litr_met_p(c,j) * dt
            col_ps%decomp_ppools_vr(c,j,i_cel_lit) = &
                 col_ps%decomp_ppools_vr(c,j,i_cel_lit) + col_pf%gap_mortality_p_to_litr_cel_p(c,j) * dt
            col_ps%decomp_ppools_vr(c,j,i_lig_lit) = &
                 col_ps%decomp_ppools_vr(c,j,i_lig_lit) + col_pf%gap_mortality_p_to_litr_lig_p(c,j) * dt
            col_ps%decomp_ppools_vr(c,j,i_cwd)     = &
                 col_ps%decomp_ppools_vr(c,j,i_cwd)     + col_pf%gap_mortality_p_to_cwdp(c,j)       * dt
         end do
      end do
      endif ! if (.not.is_active_betr_bgc))
      !------------------------------------------------------------------

      ! patch -level phosporus fluxes from gap-phase mortality

      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! displayed pools
         veg_ps%leafp(p)              =  veg_ps%leafp(p)      - veg_pf%m_leafp_to_litter(p)      * dt
         veg_ps%frootp(p)             =  veg_ps%frootp(p)     - veg_pf%m_frootp_to_litter(p)     * dt
         veg_ps%livestemp(p)          =  veg_ps%livestemp(p)  - veg_pf%m_livestemp_to_litter(p)  * dt
         veg_ps%deadstemp(p)          =  veg_ps%deadstemp(p)  - veg_pf%m_deadstemp_to_litter(p)  * dt
         veg_ps%livecrootp(p)         =  veg_ps%livecrootp(p) - veg_pf%m_livecrootp_to_litter(p) * dt
         veg_ps%deadcrootp(p)         =  veg_ps%deadcrootp(p) - veg_pf%m_deadcrootp_to_litter(p) * dt
         veg_ps%retransp(p)           =  veg_ps%retransp(p)   - veg_pf%m_retransp_to_litter(p)   * dt
         veg_ps%ppool(p)              =  veg_ps%ppool(p)      - veg_pf%m_ppool_to_litter(p)   * dt

         ! storage pools
         veg_ps%leafp_storage(p)      =  veg_ps%leafp_storage(p)      - veg_pf%m_leafp_storage_to_litter(p)      * dt
         veg_ps%frootp_storage(p)     =  veg_ps%frootp_storage(p)     - veg_pf%m_frootp_storage_to_litter(p)     * dt
         veg_ps%livestemp_storage(p)  =  veg_ps%livestemp_storage(p)  - veg_pf%m_livestemp_storage_to_litter(p)  * dt
         veg_ps%deadstemp_storage(p)  =  veg_ps%deadstemp_storage(p)  - veg_pf%m_deadstemp_storage_to_litter(p)  * dt
         veg_ps%livecrootp_storage(p) =  veg_ps%livecrootp_storage(p) - veg_pf%m_livecrootp_storage_to_litter(p) * dt
         veg_ps%deadcrootp_storage(p) =  veg_ps%deadcrootp_storage(p) - veg_pf%m_deadcrootp_storage_to_litter(p) * dt

         ! transfer pools
         veg_ps%leafp_xfer(p)         =  veg_ps%leafp_xfer(p)      - veg_pf%m_leafp_xfer_to_litter(p)      * dt
         veg_ps%frootp_xfer(p)        =  veg_ps%frootp_xfer(p)     - veg_pf%m_frootp_xfer_to_litter(p)     * dt
         veg_ps%livestemp_xfer(p)     =  veg_ps%livestemp_xfer(p)  - veg_pf%m_livestemp_xfer_to_litter(p)  * dt
         veg_ps%deadstemp_xfer(p)     =  veg_ps%deadstemp_xfer(p)  - veg_pf%m_deadstemp_xfer_to_litter(p)  * dt
         veg_ps%livecrootp_xfer(p)    =  veg_ps%livecrootp_xfer(p) - veg_pf%m_livecrootp_xfer_to_litter(p) * dt
         veg_ps%deadcrootp_xfer(p)    =  veg_ps%deadcrootp_xfer(p) - veg_pf%m_deadcrootp_xfer_to_litter(p) * dt

      end do

    end associate

  end subroutine PhosphorusStateUpdate2

  !-----------------------------------------------------------------------
  subroutine PhosphorusStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
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
            col_ps%decomp_ppools_vr(c,j,i_met_lit) = &
                 col_ps%decomp_ppools_vr(c,j,i_met_lit) + col_pf%harvest_p_to_litr_met_p(c,j) * dt
            col_ps%decomp_ppools_vr(c,j,i_cel_lit) = &
                 col_ps%decomp_ppools_vr(c,j,i_cel_lit) + col_pf%harvest_p_to_litr_cel_p(c,j) * dt
            col_ps%decomp_ppools_vr(c,j,i_lig_lit) = &
                 col_ps%decomp_ppools_vr(c,j,i_lig_lit) + col_pf%harvest_p_to_litr_lig_p(c,j) * dt
            col_ps%decomp_ppools_vr(c,j,i_cwd)     = &
                 col_ps%decomp_ppools_vr(c,j,i_cwd)     + col_pf%harvest_p_to_cwdp(c,j)       * dt
         end do
      end do
      endif ! if (.not.(use_pflotran .and. pf_cmode))
      !------------------------------------------------------------------

      ! patch-level phosporus fluxes from harvest mortality

      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! displayed pools
         veg_ps%leafp(p)      = veg_ps%leafp(p)      - veg_pf%hrv_leafp_to_litter(p)      * dt
         veg_ps%frootp(p)     = veg_ps%frootp(p)     - veg_pf%hrv_frootp_to_litter(p)     * dt
         veg_ps%livestemp(p)  = veg_ps%livestemp(p)  - veg_pf%hrv_livestemp_to_litter(p)  * dt
         veg_ps%deadstemp(p)  = veg_ps%deadstemp(p)  - veg_pf%hrv_deadstemp_to_prod10p(p) * dt
         veg_ps%deadstemp(p)  = veg_ps%deadstemp(p)  - veg_pf%hrv_deadstemp_to_prod100p(p)* dt
         veg_ps%livecrootp(p) = veg_ps%livecrootp(p) - veg_pf%hrv_livecrootp_to_litter(p) * dt
         veg_ps%deadcrootp(p) = veg_ps%deadcrootp(p) - veg_pf%hrv_deadcrootp_to_litter(p) * dt
         veg_ps%retransp(p)   = veg_ps%retransp(p)   - veg_pf%hrv_retransp_to_litter(p)   * dt
         veg_ps%ppool(p)      = veg_ps%ppool(p)      - veg_pf%hrv_ppool_to_litter(p)      * dt

       if (ivt(p) >= npcropmin) then ! skip 2 generic crops
           veg_ps%livestemp(p)= veg_ps%livestemp(p)  - veg_pf%hrv_livestemp_to_prod1p(p)  * dt
           veg_ps%leafp(p)    = veg_ps%leafp(p)      - veg_pf%hrv_leafp_to_prod1p(p)      * dt
           veg_ps%grainp(p)   = veg_ps%grainp(p)     - veg_pf%hrv_grainp_to_prod1p(p)     * dt
       end if

         ! storage pools
         veg_ps%leafp_storage(p)      = veg_ps%leafp_storage(p)      - veg_pf%hrv_leafp_storage_to_litter(p)      * dt
         veg_ps%frootp_storage(p)     = veg_ps%frootp_storage(p)     - veg_pf%hrv_frootp_storage_to_litter(p)     * dt
         veg_ps%livestemp_storage(p)  = veg_ps%livestemp_storage(p)  - veg_pf%hrv_livestemp_storage_to_litter(p)  * dt
         veg_ps%deadstemp_storage(p)  = veg_ps%deadstemp_storage(p)  - veg_pf%hrv_deadstemp_storage_to_litter(p)  * dt
         veg_ps%livecrootp_storage(p) = veg_ps%livecrootp_storage(p) - veg_pf%hrv_livecrootp_storage_to_litter(p) * dt
         veg_ps%deadcrootp_storage(p) = veg_ps%deadcrootp_storage(p) - veg_pf%hrv_deadcrootp_storage_to_litter(p) * dt

         ! transfer pools
         veg_ps%leafp_xfer(p)      = veg_ps%leafp_xfer(p)      - veg_pf%hrv_leafp_xfer_to_litter(p)      *dt
         veg_ps%frootp_xfer(p)     = veg_ps%frootp_xfer(p)     - veg_pf%hrv_frootp_xfer_to_litter(p)     *dt
         veg_ps%livestemp_xfer(p)  = veg_ps%livestemp_xfer(p)  - veg_pf%hrv_livestemp_xfer_to_litter(p)  *dt
         veg_ps%deadstemp_xfer(p)  = veg_ps%deadstemp_xfer(p)  - veg_pf%hrv_deadstemp_xfer_to_litter(p)  *dt
         veg_ps%livecrootp_xfer(p) = veg_ps%livecrootp_xfer(p) - veg_pf%hrv_livecrootp_xfer_to_litter(p) *dt
         veg_ps%deadcrootp_xfer(p) = veg_ps%deadcrootp_xfer(p) - veg_pf%hrv_deadcrootp_xfer_to_litter(p) *dt

      end do

    end associate

  end subroutine PhosphorusStateUpdate2h

end module PhosphorusStateUpdate2Mod
