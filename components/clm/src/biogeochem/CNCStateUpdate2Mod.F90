module CNCStateUpdate2Mod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for carbon state variable update, mortality fluxes.
  !
  ! !USES:
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use abortutils       , only : endrun
  use clm_time_manager , only : get_step_size
  use clm_varpar       , only : nlevdecomp, i_met_lit, i_cel_lit, i_lig_lit, i_cwd
  use CNCarbonStateType, only : carbonstate_type
  use CNCarbonFluxType , only : carbonflux_type
  use PatchType        , only : pft
  use pftvarcon        , only : npcropmin
  use clm_varctl       , only : use_pflotran, pf_cmode
  use PatchType           , only : pft   
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: CStateUpdate2
  public:: CStateUpdate2h
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
       carbonflux_vars, carbonstate_vars)
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
    !
    ! !LOCAL VARIABLES:
    integer  :: c ,p,j ! indices
    integer  :: fp,fc  ! lake filter indices
    real(r8) :: dt     ! radiation time step (seconds)
    !-----------------------------------------------------------------------
    
    associate(                   & 
         cf => carbonflux_vars , &
         cs => carbonstate_vars  &
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
               cs%decomp_cpools_vr_col(c,j,i_met_lit) = &
                    cs%decomp_cpools_vr_col(c,j,i_met_lit) + cf%gap_mortality_c_to_litr_met_c_col(c,j) * dt
               cs%decomp_cpools_vr_col(c,j,i_cel_lit) = &
                    cs%decomp_cpools_vr_col(c,j,i_cel_lit) + cf%gap_mortality_c_to_litr_cel_c_col(c,j) * dt
               cs%decomp_cpools_vr_col(c,j,i_lig_lit) = &
                    cs%decomp_cpools_vr_col(c,j,i_lig_lit) + cf%gap_mortality_c_to_litr_lig_c_col(c,j) * dt
               cs%decomp_cpools_vr_col(c,j,i_cwd) = &
                    cs%decomp_cpools_vr_col(c,j,i_cwd) + cf%gap_mortality_c_to_cwdc_col(c,j) * dt

            end do
         end do
      else if (is_active_betr_bgc) then

         do j = 1,nlevdecomp
            ! column loop
            do fc = 1,num_soilc
               c = filter_soilc(fc)

               ! column gap mortality fluxes
               cf%bgc_cpool_ext_inputs_vr_col(c,j,i_met_lit) = &
                    cf%bgc_cpool_ext_inputs_vr_col(c,j,i_met_lit) + cf%gap_mortality_c_to_litr_met_c_col(c,j) * dt
               cf%bgc_cpool_ext_inputs_vr_col(c,j,i_cel_lit) = &
                    cf%bgc_cpool_ext_inputs_vr_col(c,j,i_cel_lit) + cf%gap_mortality_c_to_litr_cel_c_col(c,j) * dt
               cf%bgc_cpool_ext_inputs_vr_col(c,j,i_lig_lit) = &
                    cf%bgc_cpool_ext_inputs_vr_col(c,j,i_lig_lit) + cf%gap_mortality_c_to_litr_lig_c_col(c,j) * dt
               cf%bgc_cpool_ext_inputs_vr_col(c,j,i_cwd) = &
                    cf%bgc_cpool_ext_inputs_vr_col(c,j,i_cwd) + cf%gap_mortality_c_to_cwdc_col(c,j) * dt

            end do
         end do
      endif


      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! patch-level carbon fluxes from gap-phase mortality
         ! displayed pools
         cs%leafc_patch(p)               = cs%leafc_patch(p)              - cf%m_leafc_to_litter_patch(p)              * dt
         cs%frootc_patch(p)              = cs%frootc_patch(p)             - cf%m_frootc_to_litter_patch(p)             * dt
         cs%livestemc_patch(p)           = cs%livestemc_patch(p)          - cf%m_livestemc_to_litter_patch(p)          * dt
         cs%deadstemc_patch(p)           = cs%deadstemc_patch(p)          - cf%m_deadstemc_to_litter_patch(p)          * dt
         cs%livecrootc_patch(p)          = cs%livecrootc_patch(p)         - cf%m_livecrootc_to_litter_patch(p)         * dt
         cs%deadcrootc_patch(p)          = cs%deadcrootc_patch(p)         - cf%m_deadcrootc_to_litter_patch(p)         * dt
         ! storage pools
         cs%leafc_storage_patch(p)       = cs%leafc_storage_patch(p)      - cf%m_leafc_storage_to_litter_patch(p)      * dt
         cs%frootc_storage_patch(p)      = cs%frootc_storage_patch(p)     - cf%m_frootc_storage_to_litter_patch(p)     * dt
         cs%livestemc_storage_patch(p)   = cs%livestemc_storage_patch(p)  - cf%m_livestemc_storage_to_litter_patch(p)  * dt
         cs%deadstemc_storage_patch(p)   = cs%deadstemc_storage_patch(p)  - cf%m_deadstemc_storage_to_litter_patch(p)  * dt
         cs%livecrootc_storage_patch(p)  = cs%livecrootc_storage_patch(p) - cf%m_livecrootc_storage_to_litter_patch(p) * dt
         cs%deadcrootc_storage_patch(p)  = cs%deadcrootc_storage_patch(p) - cf%m_deadcrootc_storage_to_litter_patch(p) * dt
         cs%gresp_storage_patch(p)       = cs%gresp_storage_patch(p)      - cf%m_gresp_storage_to_litter_patch(p)      * dt

         ! transfer pools
         cs%leafc_xfer_patch(p)          = cs%leafc_xfer_patch(p)         - cf%m_leafc_xfer_to_litter_patch(p)         * dt
         cs%frootc_xfer_patch(p)         = cs%frootc_xfer_patch(p)        - cf%m_frootc_xfer_to_litter_patch(p)        * dt
         cs%livestemc_xfer_patch(p)      = cs%livestemc_xfer_patch(p)     - cf%m_livestemc_xfer_to_litter_patch(p)     * dt
         cs%deadstemc_xfer_patch(p)      = cs%deadstemc_xfer_patch(p)     - cf%m_deadstemc_xfer_to_litter_patch(p)     * dt
         cs%livecrootc_xfer_patch(p)     = cs%livecrootc_xfer_patch(p)    - cf%m_livecrootc_xfer_to_litter_patch(p)    * dt
         cs%deadcrootc_xfer_patch(p)     = cs%deadcrootc_xfer_patch(p)    - cf%m_deadcrootc_xfer_to_litter_patch(p)    * dt
         cs%gresp_xfer_patch(p)          = cs%gresp_xfer_patch(p)         - cf%m_gresp_xfer_to_litter_patch(p)         * dt
      end do ! end of patch loop

    end associate
  end subroutine CStateUpdate2

  !-----------------------------------------------------------------------
  subroutine CStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
       carbonflux_vars, carbonstate_vars)
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
    !
    ! !LOCAL VARIABLES:
    integer :: c,p,j,k,l ! indices
    integer :: fp,fc     ! lake filter indices
    real(r8):: dt        ! radiation time step (seconds)
    !-----------------------------------------------------------------------

    associate(                   & 
         ivt => pft%itype      , & ! Input:  [integer (:)]  pft vegetation type
         cf => carbonflux_vars , &
         cs => carbonstate_vars  &
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
               cs%decomp_cpools_vr_col(c,j,i_met_lit) = &
                    cs%decomp_cpools_vr_col(c,j,i_met_lit) + cf%harvest_c_to_litr_met_c_col(c,j) * dt
               cs%decomp_cpools_vr_col(c,j,i_cel_lit) = &
                    cs%decomp_cpools_vr_col(c,j,i_cel_lit) + cf%harvest_c_to_litr_cel_c_col(c,j) * dt
               cs%decomp_cpools_vr_col(c,j,i_lig_lit) = &
                    cs%decomp_cpools_vr_col(c,j,i_lig_lit) + cf%harvest_c_to_litr_lig_c_col(c,j) * dt
               cs%decomp_cpools_vr_col(c,j,i_cwd) = &
                    cs%decomp_cpools_vr_col(c,j,i_cwd) + cf%harvest_c_to_cwdc_col(c,j)  * dt

               ! wood to product pools - states updated in CNWoodProducts()
            end do
         end do

      else if (is_active_betr_bgc) then
         do j = 1, nlevdecomp
            ! column loop
            do fc = 1,num_soilc
               c = filter_soilc(fc)          
               cf%bgc_cpool_ext_inputs_vr_col(c,j,i_met_lit) = &
                    cf%bgc_cpool_ext_inputs_vr_col(c,j,i_met_lit) + cf%harvest_c_to_litr_met_c_col(c,j) * dt
               cf%bgc_cpool_ext_inputs_vr_col(c,j,i_cel_lit) = &
                    cf%bgc_cpool_ext_inputs_vr_col(c,j,i_cel_lit) + cf%harvest_c_to_litr_cel_c_col(c,j) * dt
               cf%bgc_cpool_ext_inputs_vr_col(c,j,i_lig_lit) = &
                    cf%bgc_cpool_ext_inputs_vr_col(c,j,i_lig_lit) + cf%harvest_c_to_litr_lig_c_col(c,j) * dt
               cf%bgc_cpool_ext_inputs_vr_col(c,j,i_cwd) = &
                    cf%bgc_cpool_ext_inputs_vr_col(c,j,i_cwd) + cf%harvest_c_to_cwdc_col(c,j)  * dt
            end do
         end do
      endif

      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! patch-level carbon fluxes from harvest mortality
         ! displayed pools
         cs%leafc_patch(p)               = cs%leafc_patch(p)              - cf%hrv_leafc_to_litter_patch(p)              * dt
         cs%frootc_patch(p)              = cs%frootc_patch(p)             - cf%hrv_frootc_to_litter_patch(p)             * dt
         cs%livestemc_patch(p)           = cs%livestemc_patch(p)          - cf%hrv_livestemc_to_litter_patch(p)          * dt
         cs%deadstemc_patch(p)           = cs%deadstemc_patch(p)          - cf%hrv_deadstemc_to_prod10c_patch(p)         * dt
         cs%deadstemc_patch(p)           = cs%deadstemc_patch(p)          - cf%hrv_deadstemc_to_prod100c_patch(p)        * dt
         cs%livecrootc_patch(p)          = cs%livecrootc_patch(p)         - cf%hrv_livecrootc_to_litter_patch(p)         * dt
         cs%deadcrootc_patch(p)          = cs%deadcrootc_patch(p)         - cf%hrv_deadcrootc_to_litter_patch(p)         * dt

         ! crops
         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
             cs%livestemc_patch(p)       = cs%livestemc_patch(p)          - cf%hrv_livestemc_to_prod1c_patch(p)          *dt
             cs%leafc_patch(p)           = cs%leafc_patch(p)              - cf%hrv_leafc_to_prod1c_patch(p)              *dt
             cs%grainc_patch(p)          = cs%grainc_patch(p)             - cf%hrv_grainc_to_prod1c_patch(p)             *dt
         end if

         ! xsmrpool
         cs%xsmrpool_patch(p)            = cs%xsmrpool_patch(p)           - cf%hrv_xsmrpool_to_atm_patch(p)              * dt
         ! storage pools
         cs%leafc_storage_patch(p)       = cs%leafc_storage_patch(p)      - cf%hrv_leafc_storage_to_litter_patch(p)      * dt
         cs%frootc_storage_patch(p)      = cs%frootc_storage_patch(p)     - cf%hrv_frootc_storage_to_litter_patch(p)     * dt
         cs%livestemc_storage_patch(p)   = cs%livestemc_storage_patch(p)  - cf%hrv_livestemc_storage_to_litter_patch(p)  * dt
         cs%deadstemc_storage_patch(p)   = cs%deadstemc_storage_patch(p)  - cf%hrv_deadstemc_storage_to_litter_patch(p)  * dt
         cs%livecrootc_storage_patch(p)  = cs%livecrootc_storage_patch(p) - cf%hrv_livecrootc_storage_to_litter_patch(p) * dt
         cs%deadcrootc_storage_patch(p)  = cs%deadcrootc_storage_patch(p) - cf%hrv_deadcrootc_storage_to_litter_patch(p) * dt
         cs%gresp_storage_patch(p)       = cs%gresp_storage_patch(p)      - cf%hrv_gresp_storage_to_litter_patch(p)      * dt

         ! transfer pools
         cs%leafc_xfer_patch(p)          = cs%leafc_xfer_patch(p)         - cf%hrv_leafc_xfer_to_litter_patch(p)         * dt
         cs%frootc_xfer_patch(p)         = cs%frootc_xfer_patch(p)        - cf%hrv_frootc_xfer_to_litter_patch(p)        * dt
         cs%livestemc_xfer_patch(p)      = cs%livestemc_xfer_patch(p)     - cf%hrv_livestemc_xfer_to_litter_patch(p)     * dt
         cs%deadstemc_xfer_patch(p)      = cs%deadstemc_xfer_patch(p)     - cf%hrv_deadstemc_xfer_to_litter_patch(p)     * dt
         cs%livecrootc_xfer_patch(p)     = cs%livecrootc_xfer_patch(p)    - cf%hrv_livecrootc_xfer_to_litter_patch(p)    * dt
         cs%deadcrootc_xfer_patch(p)     = cs%deadcrootc_xfer_patch(p)    - cf%hrv_deadcrootc_xfer_to_litter_patch(p)    * dt
         cs%gresp_xfer_patch(p)          = cs%gresp_xfer_patch(p)         - cf%hrv_gresp_xfer_to_litter_patch(p)         * dt

      end do ! end of patch loop

    end associate

  end subroutine CStateUpdate2h

end module CNCStateUpdate2Mod
