module CNCStateUpdate2Mod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for carbon state variable update, mortality fluxes.
  !
  ! !USES:
  use shr_kind_mod                   , only : r8 => shr_kind_r8
  use shr_log_mod                    , only : errMsg => shr_log_errMsg
  use abortutils                     , only : endrun
  use clm_time_manager               , only : get_step_size
  use clm_varpar                     , only : nlevdecomp, i_met_lit, i_cel_lit, i_lig_lit, i_cwd
  use CNvegCarbonStateType           , only : cnveg_carbonstate_type
  use CNVegCarbonFluxType            , only : cnveg_carbonflux_type
  use SoilBiogeochemCarbonStatetype  , only : soilbiogeochem_carbonstate_type
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: CStateUpdate2
  public:: CStateUpdate2h
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
       cnveg_carbonflux_inst, cnveg_carbonstate_inst, soilbiogeochem_carbonstate_inst)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update all the prognostic carbon state
    ! variables affected by gap-phase mortality fluxes
    !
    ! !ARGUMENTS:
    integer                                , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                                , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                                , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnveg_carbonflux_type)            , intent(in)    :: cnveg_carbonflux_inst
    type(cnveg_carbonstate_type)           , intent(inout) :: cnveg_carbonstate_inst
    type(soilbiogeochem_carbonstate_type)  , intent(inout) :: soilbiogeochem_carbonstate_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: c ,p,j ! indices
    integer  :: fp,fc  ! lake filter indices
    real(r8) :: dt     ! radiation time step (seconds)
    !-----------------------------------------------------------------------
    
    associate(                                      & 
         cf_veg => cnveg_carbonflux_inst          , &
         cs_veg => cnveg_carbonstate_inst         , &

         cs_soil => soilbiogeochem_carbonstate_inst  &
         )

      ! set time steps
      dt = real( get_step_size(), r8 )

      ! column level carbon fluxes from gap-phase mortality
      do j = 1,nlevdecomp
         ! column loop
         do fc = 1,num_soilc
            c = filter_soilc(fc)

            ! column gap mortality fluxes
            cs_soil%decomp_cpools_vr_col(c,j,i_met_lit) = &
                 cs_soil%decomp_cpools_vr_col(c,j,i_met_lit) + cf_veg%gap_mortality_c_to_litr_met_c_col(c,j) * dt
            cs_soil%decomp_cpools_vr_col(c,j,i_cel_lit) = &
                 cs_soil%decomp_cpools_vr_col(c,j,i_cel_lit) + cf_veg%gap_mortality_c_to_litr_cel_c_col(c,j) * dt
            cs_soil%decomp_cpools_vr_col(c,j,i_lig_lit) = &
                 cs_soil%decomp_cpools_vr_col(c,j,i_lig_lit) + cf_veg%gap_mortality_c_to_litr_lig_c_col(c,j) * dt
            cs_soil%decomp_cpools_vr_col(c,j,i_cwd) = &
                 cs_soil%decomp_cpools_vr_col(c,j,i_cwd) + cf_veg%gap_mortality_c_to_cwdc_col(c,j) * dt

         end do
      end do

      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! patch-level carbon fluxes from gap-phase mortality
         ! displayed pools
         cs_veg%leafc_patch(p) = cs_veg%leafc_patch(p)                           &
              - cf_veg%m_leafc_to_litter_patch(p) * dt
         cs_veg%frootc_patch(p) = cs_veg%frootc_patch(p)                         &
              - cf_veg%m_frootc_to_litter_patch(p) * dt
         cs_veg%livestemc_patch(p) = cs_veg%livestemc_patch(p)                   &
              - cf_veg%m_livestemc_to_litter_patch(p) * dt
         cs_veg%deadstemc_patch(p) = cs_veg%deadstemc_patch(p)                   &
              - cf_veg%m_deadstemc_to_litter_patch(p) * dt
         cs_veg%livecrootc_patch(p) = cs_veg%livecrootc_patch(p)                 &
              - cf_veg%m_livecrootc_to_litter_patch(p) * dt
         cs_veg%deadcrootc_patch(p) = cs_veg%deadcrootc_patch(p)                 &
              - cf_veg%m_deadcrootc_to_litter_patch(p) * dt

         ! storage pools
         cs_veg%leafc_storage_patch(p) = cs_veg%leafc_storage_patch(p)           &
              - cf_veg%m_leafc_storage_to_litter_patch(p) * dt
         cs_veg%frootc_storage_patch(p) = cs_veg%frootc_storage_patch(p)         &
              - cf_veg%m_frootc_storage_to_litter_patch(p) * dt
         cs_veg%livestemc_storage_patch(p) = cs_veg%livestemc_storage_patch(p)   &
              - cf_veg%m_livestemc_storage_to_litter_patch(p) * dt
         cs_veg%deadstemc_storage_patch(p) = cs_veg%deadstemc_storage_patch(p)   &
              - cf_veg%m_deadstemc_storage_to_litter_patch(p) * dt
         cs_veg%livecrootc_storage_patch(p) = cs_veg%livecrootc_storage_patch(p) &
              - cf_veg%m_livecrootc_storage_to_litter_patch(p) * dt
         cs_veg%deadcrootc_storage_patch(p) = cs_veg%deadcrootc_storage_patch(p) &
              - cf_veg%m_deadcrootc_storage_to_litter_patch(p) * dt
         cs_veg%gresp_storage_patch(p) = cs_veg%gresp_storage_patch(p)           &
              - cf_veg%m_gresp_storage_to_litter_patch(p) * dt

         ! transfer pools
         cs_veg%leafc_xfer_patch(p) = cs_veg%leafc_xfer_patch(p)                 &
              - cf_veg%m_leafc_xfer_to_litter_patch(p) * dt
         cs_veg%frootc_xfer_patch(p) = cs_veg%frootc_xfer_patch(p)               &
              - cf_veg%m_frootc_xfer_to_litter_patch(p) * dt
         cs_veg%livestemc_xfer_patch(p) = cs_veg%livestemc_xfer_patch(p)         &
              - cf_veg%m_livestemc_xfer_to_litter_patch(p) * dt
         cs_veg%deadstemc_xfer_patch(p) = cs_veg%deadstemc_xfer_patch(p)         &
              - cf_veg%m_deadstemc_xfer_to_litter_patch(p) * dt
         cs_veg%livecrootc_xfer_patch(p) = cs_veg%livecrootc_xfer_patch(p)       &
              - cf_veg%m_livecrootc_xfer_to_litter_patch(p) * dt
         cs_veg%deadcrootc_xfer_patch(p) = cs_veg%deadcrootc_xfer_patch(p)       &
              - cf_veg%m_deadcrootc_xfer_to_litter_patch(p) * dt
         cs_veg%gresp_xfer_patch(p) = cs_veg%gresp_xfer_patch(p)                 &
              - cf_veg%m_gresp_xfer_to_litter_patch(p) * dt
      end do ! end of patch loop

    end associate

  end subroutine CStateUpdate2

  !-----------------------------------------------------------------------
  subroutine CStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
       cnveg_carbonflux_inst, cnveg_carbonstate_inst, soilbiogeochem_carbonstate_inst)
    !
    ! !DESCRIPTION:
    ! Update all the prognostic carbon state
    ! variables affected by harvest mortality fluxes
    !
    ! !ARGUMENTS:
    integer                                , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                                , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                                , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnveg_carbonflux_type)            , intent(in)    :: cnveg_carbonflux_inst
    type(cnveg_carbonstate_type)           , intent(inout) :: cnveg_carbonstate_inst
    type(soilbiogeochem_carbonstate_type)  , intent(inout) :: soilbiogeochem_carbonstate_inst
    !
    ! !LOCAL VARIABLES:
    integer :: c,p,j,k,l ! indices
    integer :: fp,fc     ! lake filter indices
    real(r8):: dt        ! radiation time step (seconds)
    !-----------------------------------------------------------------------

    associate(                                     & 
         cf_veg => cnveg_carbonflux_inst         , &
         cs_veg => cnveg_carbonstate_inst        , &
         cs_soil => soilbiogeochem_carbonstate_inst &
         )
     
      ! set time steps
      dt = real( get_step_size(), r8 )

      ! column level carbon fluxes from harvest mortality
      do j = 1, nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)

            ! column harvest fluxes
            cs_soil%decomp_cpools_vr_col(c,j,i_met_lit) = &
                 cs_soil%decomp_cpools_vr_col(c,j,i_met_lit) + cf_veg%harvest_c_to_litr_met_c_col(c,j) * dt
            cs_soil%decomp_cpools_vr_col(c,j,i_cel_lit) = &
                 cs_soil%decomp_cpools_vr_col(c,j,i_cel_lit) + cf_veg%harvest_c_to_litr_cel_c_col(c,j) * dt
            cs_soil%decomp_cpools_vr_col(c,j,i_lig_lit) = &
                 cs_soil%decomp_cpools_vr_col(c,j,i_lig_lit) + cf_veg%harvest_c_to_litr_lig_c_col(c,j) * dt
            cs_soil%decomp_cpools_vr_col(c,j,i_cwd) = &
                 cs_soil%decomp_cpools_vr_col(c,j,i_cwd) + cf_veg%harvest_c_to_cwdc_col(c,j)  * dt

            ! wood to product pools - states updated in CNWoodProducts()
         end do
      end do

      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! patch-level carbon fluxes from harvest mortality
         ! displayed pools
         cs_veg%leafc_patch(p) = cs_veg%leafc_patch(p)                           &
              - cf_veg%hrv_leafc_to_litter_patch(p) * dt
         cs_veg%frootc_patch(p) = cs_veg%frootc_patch(p)                         &
              - cf_veg%hrv_frootc_to_litter_patch(p) * dt
         cs_veg%livestemc_patch(p) = cs_veg%livestemc_patch(p)                   &
              - cf_veg%hrv_livestemc_to_litter_patch(p) * dt
         cs_veg%deadstemc_patch(p) = cs_veg%deadstemc_patch(p)                   &
              - cf_veg%hrv_deadstemc_to_prod10c_patch(p) * dt
         cs_veg%deadstemc_patch(p) = cs_veg%deadstemc_patch(p)                   &
              - cf_veg%hrv_deadstemc_to_prod100c_patch(p) * dt
         cs_veg%livecrootc_patch(p) = cs_veg%livecrootc_patch(p)                 &
              - cf_veg%hrv_livecrootc_to_litter_patch(p) * dt
         cs_veg%deadcrootc_patch(p) = cs_veg%deadcrootc_patch(p)                 &
              - cf_veg%hrv_deadcrootc_to_litter_patch(p) * dt

         ! xsmrpool
         cs_veg%xsmrpool_patch(p) = cs_veg%xsmrpool_patch(p)                     &
              - cf_veg%hrv_xsmrpool_to_atm_patch(p) * dt

         ! storage pools
         cs_veg%leafc_storage_patch(p) = cs_veg%leafc_storage_patch(p)           &
              - cf_veg%hrv_leafc_storage_to_litter_patch(p) * dt
         cs_veg%frootc_storage_patch(p) = cs_veg%frootc_storage_patch(p)         &
              - cf_veg%hrv_frootc_storage_to_litter_patch(p) * dt
         cs_veg%livestemc_storage_patch(p) = cs_veg%livestemc_storage_patch(p)   &
              - cf_veg%hrv_livestemc_storage_to_litter_patch(p) * dt
         cs_veg%deadstemc_storage_patch(p) = cs_veg%deadstemc_storage_patch(p)   &
              - cf_veg%hrv_deadstemc_storage_to_litter_patch(p) * dt
         cs_veg%livecrootc_storage_patch(p) = cs_veg%livecrootc_storage_patch(p) &
              - cf_veg%hrv_livecrootc_storage_to_litter_patch(p) * dt
         cs_veg%deadcrootc_storage_patch(p) = cs_veg%deadcrootc_storage_patch(p) &
              - cf_veg%hrv_deadcrootc_storage_to_litter_patch(p) * dt
         cs_veg%gresp_storage_patch(p) = cs_veg%gresp_storage_patch(p)           &
              - cf_veg%hrv_gresp_storage_to_litter_patch(p) * dt

         ! transfer pools
         cs_veg%leafc_xfer_patch(p) = cs_veg%leafc_xfer_patch(p)                 &
              - cf_veg%hrv_leafc_xfer_to_litter_patch(p) * dt
         cs_veg%frootc_xfer_patch(p) = cs_veg%frootc_xfer_patch(p)               &
              - cf_veg%hrv_frootc_xfer_to_litter_patch(p) * dt
         cs_veg%livestemc_xfer_patch(p) = cs_veg%livestemc_xfer_patch(p)         &
              - cf_veg%hrv_livestemc_xfer_to_litter_patch(p) * dt
         cs_veg%deadstemc_xfer_patch(p) = cs_veg%deadstemc_xfer_patch(p)         &
              - cf_veg%hrv_deadstemc_xfer_to_litter_patch(p) * dt
         cs_veg%livecrootc_xfer_patch(p) = cs_veg%livecrootc_xfer_patch(p)       &
              - cf_veg%hrv_livecrootc_xfer_to_litter_patch(p) * dt
         cs_veg%deadcrootc_xfer_patch(p) = cs_veg%deadcrootc_xfer_patch(p)       &
              - cf_veg%hrv_deadcrootc_xfer_to_litter_patch(p) * dt
         cs_veg%gresp_xfer_patch(p) = cs_veg%gresp_xfer_patch(p)                 &
              - cf_veg%hrv_gresp_xfer_to_litter_patch(p) * dt

      end do ! end of patch loop

    end associate

  end subroutine CStateUpdate2h

end module CNCStateUpdate2Mod
