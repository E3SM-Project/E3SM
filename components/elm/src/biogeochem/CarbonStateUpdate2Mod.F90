module CarbonStateUpdate2Mod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for carbon state variable update, mortality fluxes.
  !
  ! !USES:
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use abortutils       , only : endrun
  use clm_time_manager , only : get_step_size
  use elm_varpar       , only : nlevdecomp, i_met_lit, i_cel_lit, i_lig_lit, i_cwd
  use CNCarbonStateType, only : carbonstate_type
  use CNCarbonFluxType , only : carbonflux_type
  use pftvarcon        , only : npcropmin
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
       carbonflux_vars, carbonstate_vars, col_csv2, veg_csv2, col_cfv2, veg_cfv2)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update all the prognostic carbon state
    ! variables affected by gap-phase mortality fluxes
    !
    use tracer_varcon, only : is_active_betr_bgc      
    ! !ARGUMENTS:
    integer                , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(carbonflux_type)  , intent(inout) :: carbonflux_vars
    type(carbonstate_type) , intent(inout) :: carbonstate_vars
    type(column_carbon_state),intent(inout):: col_csv2
    type(vegetation_carbon_state),intent(inout) :: veg_csv2
    type(column_carbon_flux)     ,intent(inout) :: col_cfv2
    type(vegetation_carbon_flux) ,intent(inout) :: veg_cfv2
    !
    ! !LOCAL VARIABLES:
    integer  :: c ,p,j ! indices
    integer  :: fp,fc  ! lake filter indices
    real(r8) :: dt     ! radiation time step (seconds)
    !-----------------------------------------------------------------------
    
    associate(                   & 
         cf => carbonflux_vars , &
         cs => carbonstate_vars, &
         csv2 => col_csv2      , &
         vcsv2             => veg_csv2                                , &
         ccfv2             => col_cfv2                                , &
         vcfv2             => veg_cfv2                                 &
         )

      ! set time steps
      dt = real( get_step_size(), r8 )

     if (  .not. is_active_betr_bgc          .and. &
          (.not.(use_pflotran .and. pf_cmode))) then
         ! column level carbon fluxes from gap-phase mortality
         do j = 1,nlevdecomp
            ! column loop
            do fc = 1,num_soilc
               c = filter_soilc(fc)               

               ! column gap mortality fluxes
               csv2%decomp_cpools_vr(c,j,i_met_lit) = &
                    csv2%decomp_cpools_vr(c,j,i_met_lit) + ccfv2%gap_mortality_c_to_litr_met_c(c,j) * dt
               csv2%decomp_cpools_vr(c,j,i_cel_lit) = &
                    csv2%decomp_cpools_vr(c,j,i_cel_lit) + ccfv2%gap_mortality_c_to_litr_cel_c(c,j) * dt
               csv2%decomp_cpools_vr(c,j,i_lig_lit) = &
                    csv2%decomp_cpools_vr(c,j,i_lig_lit) + ccfv2%gap_mortality_c_to_litr_lig_c(c,j) * dt
               csv2%decomp_cpools_vr(c,j,i_cwd) = &
                    csv2%decomp_cpools_vr(c,j,i_cwd) + ccfv2%gap_mortality_c_to_cwdc(c,j) * dt

            end do
         end do
      endif


      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! patch-level carbon fluxes from gap-phase mortality
         ! displayed pools
         vcsv2%leafc(p)               = vcsv2%leafc(p)              - vcfv2%m_leafc_to_litter(p)              * dt
         vcsv2%frootc(p)              = vcsv2%frootc(p)             - vcfv2%m_frootc_to_litter(p)             * dt
         vcsv2%livestemc(p)           = vcsv2%livestemc(p)          - vcfv2%m_livestemc_to_litter(p)          * dt
         vcsv2%deadstemc(p)           = vcsv2%deadstemc(p)          - vcfv2%m_deadstemc_to_litter(p)          * dt
         vcsv2%livecrootc(p)          = vcsv2%livecrootc(p)         - vcfv2%m_livecrootc_to_litter(p)         * dt
         vcsv2%deadcrootc(p)          = vcsv2%deadcrootc(p)         - vcfv2%m_deadcrootc_to_litter(p)         * dt
         ! storage pools
         vcsv2%leafc_storage(p)       = vcsv2%leafc_storage(p)      - vcfv2%m_leafc_storage_to_litter(p)      * dt
         vcsv2%frootc_storage(p)      = vcsv2%frootc_storage(p)     - vcfv2%m_frootc_storage_to_litter(p)     * dt
         vcsv2%livestemc_storage(p)   = vcsv2%livestemc_storage(p)  - vcfv2%m_livestemc_storage_to_litter(p)  * dt
         vcsv2%deadstemc_storage(p)   = vcsv2%deadstemc_storage(p)  - vcfv2%m_deadstemc_storage_to_litter(p)  * dt
         vcsv2%livecrootc_storage(p)  = vcsv2%livecrootc_storage(p) - vcfv2%m_livecrootc_storage_to_litter(p) * dt
         vcsv2%deadcrootc_storage(p)  = vcsv2%deadcrootc_storage(p) - vcfv2%m_deadcrootc_storage_to_litter(p) * dt
         vcsv2%gresp_storage(p)       = vcsv2%gresp_storage(p)      - vcfv2%m_gresp_storage_to_litter(p)      * dt
         vcsv2%cpool(p)               = vcsv2%cpool(p)              - vcfv2%m_cpool_to_litter(p)              * dt

         ! transfer pools
         vcsv2%leafc_xfer(p)          = vcsv2%leafc_xfer(p)         - vcfv2%m_leafc_xfer_to_litter(p)         * dt
         vcsv2%frootc_xfer(p)         = vcsv2%frootc_xfer(p)        - vcfv2%m_frootc_xfer_to_litter(p)        * dt
         vcsv2%livestemc_xfer(p)      = vcsv2%livestemc_xfer(p)     - vcfv2%m_livestemc_xfer_to_litter(p)     * dt
         vcsv2%deadstemc_xfer(p)      = vcsv2%deadstemc_xfer(p)     - vcfv2%m_deadstemc_xfer_to_litter(p)     * dt
         vcsv2%livecrootc_xfer(p)     = vcsv2%livecrootc_xfer(p)    - vcfv2%m_livecrootc_xfer_to_litter(p)    * dt
         vcsv2%deadcrootc_xfer(p)     = vcsv2%deadcrootc_xfer(p)    - vcfv2%m_deadcrootc_xfer_to_litter(p)    * dt
         vcsv2%gresp_xfer(p)          = vcsv2%gresp_xfer(p)         - vcfv2%m_gresp_xfer_to_litter(p)         * dt
      end do ! end of patch loop

    end associate
  end subroutine CarbonStateUpdate2

  !-----------------------------------------------------------------------
  subroutine CarbonStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
       carbonflux_vars, carbonstate_vars, col_csv2, veg_csv2, col_cfv2, veg_cfv2)
    !
    ! !DESCRIPTION:
    ! Update all the prognostic carbon state
    ! variables affected by harvest mortality fluxes
    !
    use tracer_varcon,  only : is_active_betr_bgc      
    ! !ARGUMENTS:
    integer                , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(carbonflux_type)  , intent(inout) :: carbonflux_vars
    type(carbonstate_type) , intent(inout) :: carbonstate_vars
    type(column_carbon_state),intent(inout):: col_csv2
    type(vegetation_carbon_state),intent(inout) :: veg_csv2
    type(column_carbon_flux)     ,intent(inout) :: col_cfv2
    type(vegetation_carbon_flux) ,intent(inout) :: veg_cfv2
    !
    ! !LOCAL VARIABLES:
    integer :: c,p,j,k,l ! indices
    integer :: fp,fc     ! lake filter indices
    real(r8):: dt        ! radiation time step (seconds)
    !-----------------------------------------------------------------------

    associate(                   & 
         ivt => veg_pp%itype      , & ! Input:  [integer (:)]  pft vegetation type
         cf => carbonflux_vars , &
         cs => carbonstate_vars, &
         csv2 => col_csv2      , &
         vcsv2             => veg_csv2                                , &
         ccfv2             => col_cfv2                                , &
         vcfv2             => veg_cfv2                                 &
         )
     
      ! set time steps
      dt = real( get_step_size(), r8 )

      if ( (.not. is_active_betr_bgc) .and. &
           .not.(use_pflotran .and. pf_cmode)) then
         ! column level carbon fluxes from harvest mortality
         do j = 1, nlevdecomp
            ! column loop
            do fc = 1,num_soilc
               c = filter_soilc(fc)

               ! column harvest fluxes
               csv2%decomp_cpools_vr(c,j,i_met_lit) = &
                    csv2%decomp_cpools_vr(c,j,i_met_lit) + ccfv2%harvest_c_to_litr_met_c(c,j) * dt
               csv2%decomp_cpools_vr(c,j,i_cel_lit) = &
                    csv2%decomp_cpools_vr(c,j,i_cel_lit) + ccfv2%harvest_c_to_litr_cel_c(c,j) * dt
               csv2%decomp_cpools_vr(c,j,i_lig_lit) = &
                    csv2%decomp_cpools_vr(c,j,i_lig_lit) + ccfv2%harvest_c_to_litr_lig_c(c,j) * dt
               csv2%decomp_cpools_vr(c,j,i_cwd) = &
                    csv2%decomp_cpools_vr(c,j,i_cwd) + ccfv2%harvest_c_to_cwdc(c,j)  * dt

               ! wood to product pools - states updated in WoodProducts()
            end do
         end do

      endif

      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! patch-level carbon fluxes from harvest mortality
         ! displayed pools
         vcsv2%leafc(p)               = vcsv2%leafc(p)              - vcfv2%hrv_leafc_to_litter(p)              * dt
         vcsv2%frootc(p)              = vcsv2%frootc(p)             - vcfv2%hrv_frootc_to_litter(p)             * dt
         vcsv2%livestemc(p)           = vcsv2%livestemc(p)          - vcfv2%hrv_livestemc_to_litter(p)          * dt
         vcsv2%deadstemc(p)           = vcsv2%deadstemc(p)          - vcfv2%hrv_deadstemc_to_prod10c(p)         * dt
         vcsv2%deadstemc(p)           = vcsv2%deadstemc(p)          - vcfv2%hrv_deadstemc_to_prod100c(p)        * dt
         vcsv2%livecrootc(p)          = vcsv2%livecrootc(p)         - vcfv2%hrv_livecrootc_to_litter(p)         * dt
         vcsv2%deadcrootc(p)          = vcsv2%deadcrootc(p)         - vcfv2%hrv_deadcrootc_to_litter(p)         * dt

         ! crops
         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
             vcsv2%livestemc(p)       = vcsv2%livestemc(p)          - vcfv2%hrv_livestemc_to_prod1c(p)          *dt
             vcsv2%leafc(p)           = vcsv2%leafc(p)              - vcfv2%hrv_leafc_to_prod1c(p)              *dt
             vcsv2%grainc(p)          = vcsv2%grainc(p)             - vcfv2%hrv_grainc_to_prod1c(p)             *dt
         end if

         ! xsmrpool
         vcsv2%xsmrpool(p)            = vcsv2%xsmrpool(p)           - vcfv2%hrv_xsmrpool_to_atm(p)              * dt
         ! storage pools
         vcsv2%leafc_storage(p)       = vcsv2%leafc_storage(p)      - vcfv2%hrv_leafc_storage_to_litter(p)      * dt
         vcsv2%frootc_storage(p)      = vcsv2%frootc_storage(p)     - vcfv2%hrv_frootc_storage_to_litter(p)     * dt
         vcsv2%livestemc_storage(p)   = vcsv2%livestemc_storage(p)  - vcfv2%hrv_livestemc_storage_to_litter(p)  * dt
         vcsv2%deadstemc_storage(p)   = vcsv2%deadstemc_storage(p)  - vcfv2%hrv_deadstemc_storage_to_litter(p)  * dt
         vcsv2%livecrootc_storage(p)  = vcsv2%livecrootc_storage(p) - vcfv2%hrv_livecrootc_storage_to_litter(p) * dt
         vcsv2%deadcrootc_storage(p)  = vcsv2%deadcrootc_storage(p) - vcfv2%hrv_deadcrootc_storage_to_litter(p) * dt
         vcsv2%gresp_storage(p)       = vcsv2%gresp_storage(p)      - vcfv2%hrv_gresp_storage_to_litter(p)      * dt
         vcsv2%cpool(p)               = vcsv2%cpool(p)              - vcfv2%hrv_cpool_to_litter(p)              * dt

         ! transfer pools
         vcsv2%leafc_xfer(p)          = vcsv2%leafc_xfer(p)         - vcfv2%hrv_leafc_xfer_to_litter(p)         * dt
         vcsv2%frootc_xfer(p)         = vcsv2%frootc_xfer(p)        - vcfv2%hrv_frootc_xfer_to_litter(p)        * dt
         vcsv2%livestemc_xfer(p)      = vcsv2%livestemc_xfer(p)     - vcfv2%hrv_livestemc_xfer_to_litter(p)     * dt
         vcsv2%deadstemc_xfer(p)      = vcsv2%deadstemc_xfer(p)     - vcfv2%hrv_deadstemc_xfer_to_litter(p)     * dt
         vcsv2%livecrootc_xfer(p)     = vcsv2%livecrootc_xfer(p)    - vcfv2%hrv_livecrootc_xfer_to_litter(p)    * dt
         vcsv2%deadcrootc_xfer(p)     = vcsv2%deadcrootc_xfer(p)    - vcfv2%hrv_deadcrootc_xfer_to_litter(p)    * dt
         vcsv2%gresp_xfer(p)          = vcsv2%gresp_xfer(p)         - vcfv2%hrv_gresp_xfer_to_litter(p)         * dt

      end do ! end of patch loop

    end associate

  end subroutine CarbonStateUpdate2h

end module CarbonStateUpdate2Mod
