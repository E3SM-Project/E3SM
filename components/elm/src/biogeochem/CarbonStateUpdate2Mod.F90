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
  use tracer_varcon    , only : is_active_betr_bgc
  use VegetationType        , only : veg_pp
  use ColumnDataType         , only : column_carbon_state, column_carbon_flux
  use VegetationDataType     , only : vegetation_carbon_state, vegetation_carbon_flux
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
      !$acc routine seq
    use tracer_varcon, only : is_active_betr_bgc
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



     if (  .not. is_active_betr_bgc          .and. &
          (.not.(use_pflotran .and. pf_cmode))) then
         ! column level carbon fluxes from gap-phase mortality
         do j = 1,nlevdecomp
            ! column loop
            do fc = 1,num_soilc
               c = filter_soilc(fc)

               ! column gap mortality fluxes
               col_cs%decomp_cpools_vr(c,j,i_met_lit) = &
                    col_cs%decomp_cpools_vr(c,j,i_met_lit) + col_cf%gap_mortality_c_to_litr_met_c(c,j) * dt
               col_cs%decomp_cpools_vr(c,j,i_cel_lit) = &
                    col_cs%decomp_cpools_vr(c,j,i_cel_lit) + col_cf%gap_mortality_c_to_litr_cel_c(c,j) * dt
               col_cs%decomp_cpools_vr(c,j,i_lig_lit) = &
                    col_cs%decomp_cpools_vr(c,j,i_lig_lit) + col_cf%gap_mortality_c_to_litr_lig_c(c,j) * dt
               col_cs%decomp_cpools_vr(c,j,i_cwd) = &
                    col_cs%decomp_cpools_vr(c,j,i_cwd) + col_cf%gap_mortality_c_to_cwdc(c,j) * dt

            end do
         end do
      endif


      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! patch-level carbon fluxes from gap-phase mortality
         ! displayed pools
         veg_cs%leafc(p)               = veg_cs%leafc(p)              - veg_cf%m_leafc_to_litter(p)              * dt
         veg_cs%frootc(p)              = veg_cs%frootc(p)             - veg_cf%m_frootc_to_litter(p)             * dt
         veg_cs%livestemc(p)           = veg_cs%livestemc(p)          - veg_cf%m_livestemc_to_litter(p)          * dt
         veg_cs%deadstemc(p)           = veg_cs%deadstemc(p)          - veg_cf%m_deadstemc_to_litter(p)          * dt
         veg_cs%livecrootc(p)          = veg_cs%livecrootc(p)         - veg_cf%m_livecrootc_to_litter(p)         * dt
         veg_cs%deadcrootc(p)          = veg_cs%deadcrootc(p)         - veg_cf%m_deadcrootc_to_litter(p)         * dt
         ! storage pools
         veg_cs%leafc_storage(p)       = veg_cs%leafc_storage(p)      - veg_cf%m_leafc_storage_to_litter(p)      * dt
         veg_cs%frootc_storage(p)      = veg_cs%frootc_storage(p)     - veg_cf%m_frootc_storage_to_litter(p)     * dt
         veg_cs%livestemc_storage(p)   = veg_cs%livestemc_storage(p)  - veg_cf%m_livestemc_storage_to_litter(p)  * dt
         veg_cs%deadstemc_storage(p)   = veg_cs%deadstemc_storage(p)  - veg_cf%m_deadstemc_storage_to_litter(p)  * dt
         veg_cs%livecrootc_storage(p)  = veg_cs%livecrootc_storage(p) - veg_cf%m_livecrootc_storage_to_litter(p) * dt
         veg_cs%deadcrootc_storage(p)  = veg_cs%deadcrootc_storage(p) - veg_cf%m_deadcrootc_storage_to_litter(p) * dt
         veg_cs%gresp_storage(p)       = veg_cs%gresp_storage(p)      - veg_cf%m_gresp_storage_to_litter(p)      * dt
         veg_cs%cpool(p)               = veg_cs%cpool(p)              - veg_cf%m_cpool_to_litter(p)              * dt

         ! transfer pools
         veg_cs%leafc_xfer(p)          = veg_cs%leafc_xfer(p)         - veg_cf%m_leafc_xfer_to_litter(p)         * dt
         veg_cs%frootc_xfer(p)         = veg_cs%frootc_xfer(p)        - veg_cf%m_frootc_xfer_to_litter(p)        * dt
         veg_cs%livestemc_xfer(p)      = veg_cs%livestemc_xfer(p)     - veg_cf%m_livestemc_xfer_to_litter(p)     * dt
         veg_cs%deadstemc_xfer(p)      = veg_cs%deadstemc_xfer(p)     - veg_cf%m_deadstemc_xfer_to_litter(p)     * dt
         veg_cs%livecrootc_xfer(p)     = veg_cs%livecrootc_xfer(p)    - veg_cf%m_livecrootc_xfer_to_litter(p)    * dt
         veg_cs%deadcrootc_xfer(p)     = veg_cs%deadcrootc_xfer(p)    - veg_cf%m_deadcrootc_xfer_to_litter(p)    * dt
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
      !$acc routine seq
    use tracer_varcon,  only : is_active_betr_bgc
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

      if ( (.not. is_active_betr_bgc) .and. &
           .not.(use_pflotran .and. pf_cmode)) then
         ! column level carbon fluxes from harvest mortality
         do j = 1, nlevdecomp
            ! column loop
            do fc = 1,num_soilc
               c = filter_soilc(fc)

               ! column harvest fluxes
               col_cs%decomp_cpools_vr(c,j,i_met_lit) = &
                    col_cs%decomp_cpools_vr(c,j,i_met_lit) + col_cf%harvest_c_to_litr_met_c(c,j) * dt
               col_cs%decomp_cpools_vr(c,j,i_cel_lit) = &
                    col_cs%decomp_cpools_vr(c,j,i_cel_lit) + col_cf%harvest_c_to_litr_cel_c(c,j) * dt
               col_cs%decomp_cpools_vr(c,j,i_lig_lit) = &
                    col_cs%decomp_cpools_vr(c,j,i_lig_lit) + col_cf%harvest_c_to_litr_lig_c(c,j) * dt
               col_cs%decomp_cpools_vr(c,j,i_cwd) = &
                    col_cs%decomp_cpools_vr(c,j,i_cwd) + col_cf%harvest_c_to_cwdc(c,j)  * dt

               ! wood to product pools - states updated in WoodProducts()
            end do
         end do

      endif

      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! patch-level carbon fluxes from harvest mortality
         ! displayed pools
         veg_cs%leafc(p)               = veg_cs%leafc(p)              - veg_cf%hrv_leafc_to_litter(p)              * dt
         veg_cs%frootc(p)              = veg_cs%frootc(p)             - veg_cf%hrv_frootc_to_litter(p)             * dt
         veg_cs%livestemc(p)           = veg_cs%livestemc(p)          - veg_cf%hrv_livestemc_to_litter(p)          * dt
         veg_cs%deadstemc(p)           = veg_cs%deadstemc(p)          - veg_cf%hrv_deadstemc_to_prod10c(p)         * dt
         veg_cs%deadstemc(p)           = veg_cs%deadstemc(p)          - veg_cf%hrv_deadstemc_to_prod100c(p)        * dt
         veg_cs%livecrootc(p)          = veg_cs%livecrootc(p)         - veg_cf%hrv_livecrootc_to_litter(p)         * dt
         veg_cs%deadcrootc(p)          = veg_cs%deadcrootc(p)         - veg_cf%hrv_deadcrootc_to_litter(p)         * dt

         ! crops
         if (iscft(ivt(p)) >= 1) then ! skip 2 generic crops
             veg_cs%livestemc(p)       = veg_cs%livestemc(p)          - veg_cf%hrv_livestemc_to_prod1c(p)          *dt
             veg_cs%leafc(p)           = veg_cs%leafc(p)              - veg_cf%hrv_leafc_to_prod1c(p)              *dt
             veg_cs%grainc(p)          = veg_cs%grainc(p)             - veg_cf%hrv_grainc_to_prod1c(p)             *dt
         end if

         ! xsmrpool
         veg_cs%xsmrpool(p)            = veg_cs%xsmrpool(p)           - veg_cf%hrv_xsmrpool_to_atm(p)              * dt
         ! storage pools
         veg_cs%leafc_storage(p)       = veg_cs%leafc_storage(p)      - veg_cf%hrv_leafc_storage_to_litter(p)      * dt
         veg_cs%frootc_storage(p)      = veg_cs%frootc_storage(p)     - veg_cf%hrv_frootc_storage_to_litter(p)     * dt
         veg_cs%livestemc_storage(p)   = veg_cs%livestemc_storage(p)  - veg_cf%hrv_livestemc_storage_to_litter(p)  * dt
         veg_cs%deadstemc_storage(p)   = veg_cs%deadstemc_storage(p)  - veg_cf%hrv_deadstemc_storage_to_litter(p)  * dt
         veg_cs%livecrootc_storage(p)  = veg_cs%livecrootc_storage(p) - veg_cf%hrv_livecrootc_storage_to_litter(p) * dt
         veg_cs%deadcrootc_storage(p)  = veg_cs%deadcrootc_storage(p) - veg_cf%hrv_deadcrootc_storage_to_litter(p) * dt
         veg_cs%gresp_storage(p)       = veg_cs%gresp_storage(p)      - veg_cf%hrv_gresp_storage_to_litter(p)      * dt
         veg_cs%cpool(p)               = veg_cs%cpool(p)              - veg_cf%hrv_cpool_to_litter(p)              * dt

         ! transfer pools
         veg_cs%leafc_xfer(p)          = veg_cs%leafc_xfer(p)         - veg_cf%hrv_leafc_xfer_to_litter(p)         * dt
         veg_cs%frootc_xfer(p)         = veg_cs%frootc_xfer(p)        - veg_cf%hrv_frootc_xfer_to_litter(p)        * dt
         veg_cs%livestemc_xfer(p)      = veg_cs%livestemc_xfer(p)     - veg_cf%hrv_livestemc_xfer_to_litter(p)     * dt
         veg_cs%deadstemc_xfer(p)      = veg_cs%deadstemc_xfer(p)     - veg_cf%hrv_deadstemc_xfer_to_litter(p)     * dt
         veg_cs%livecrootc_xfer(p)     = veg_cs%livecrootc_xfer(p)    - veg_cf%hrv_livecrootc_xfer_to_litter(p)    * dt
         veg_cs%deadcrootc_xfer(p)     = veg_cs%deadcrootc_xfer(p)    - veg_cf%hrv_deadcrootc_xfer_to_litter(p)    * dt
         veg_cs%gresp_xfer(p)          = veg_cs%gresp_xfer(p)         - veg_cf%hrv_gresp_xfer_to_litter(p)         * dt

      end do ! end of patch loop

    end associate

  end subroutine CarbonStateUpdate2h

end module CarbonStateUpdate2Mod
