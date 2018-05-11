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
  use VegetationType        , only : veg_pp
  use pftvarcon        , only : npcropmin
  use clm_varctl       , only : use_pflotran, pf_cmode
  use VegetationType           , only : veg_pp   
  use tracer_varcon    , only : is_active_betr_bgc
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
               cs%decomp_pools_vr_col(c,j,i_met_lit) = &
                    cs%decomp_pools_vr_col(c,j,i_met_lit) + cf%gap_mortality_to_litr_met_col(c,j) * dt
               cs%decomp_pools_vr_col(c,j,i_cel_lit) = &
                    cs%decomp_pools_vr_col(c,j,i_cel_lit) + cf%gap_mortality_to_litr_cel_col(c,j) * dt
               cs%decomp_pools_vr_col(c,j,i_lig_lit) = &
                    cs%decomp_pools_vr_col(c,j,i_lig_lit) + cf%gap_mortality_to_litr_lig_col(c,j) * dt
               cs%decomp_pools_vr_col(c,j,i_cwd) = &
                    cs%decomp_pools_vr_col(c,j,i_cwd) + cf%gap_mortality_to_cwd_col(c,j) * dt

            end do
         end do
      endif


      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! patch-level carbon fluxes from gap-phase mortality
         ! displayed pools
         cs%leaf_patch(p)               = cs%leaf_patch(p)              - cf%m_leaf_to_litter_patch(p)              * dt
         cs%froot_patch(p)              = cs%froot_patch(p)             - cf%m_froot_to_litter_patch(p)             * dt
         cs%livestem_patch(p)           = cs%livestem_patch(p)          - cf%m_livestem_to_litter_patch(p)          * dt
         cs%deadstem_patch(p)           = cs%deadstem_patch(p)          - cf%m_deadstem_to_litter_patch(p)          * dt
         cs%livecroot_patch(p)          = cs%livecroot_patch(p)         - cf%m_livecroot_to_litter_patch(p)         * dt
         cs%deadcroot_patch(p)          = cs%deadcroot_patch(p)         - cf%m_deadcroot_to_litter_patch(p)         * dt
         ! storage pools
         cs%leaf_storage_patch(p)       = cs%leaf_storage_patch(p)      - cf%m_leaf_storage_to_litter_patch(p)      * dt
         cs%froot_storage_patch(p)      = cs%froot_storage_patch(p)     - cf%m_froot_storage_to_litter_patch(p)     * dt
         cs%livestem_storage_patch(p)   = cs%livestem_storage_patch(p)  - cf%m_livestem_storage_to_litter_patch(p)  * dt
         cs%deadstem_storage_patch(p)   = cs%deadstem_storage_patch(p)  - cf%m_deadstem_storage_to_litter_patch(p)  * dt
         cs%livecroot_storage_patch(p)  = cs%livecroot_storage_patch(p) - cf%m_livecroot_storage_to_litter_patch(p) * dt
         cs%deadcroot_storage_patch(p)  = cs%deadcroot_storage_patch(p) - cf%m_deadcroot_storage_to_litter_patch(p) * dt
         cs%gresp_storage_patch(p)       = cs%gresp_storage_patch(p)      - cf%m_gresp_storage_to_litter_patch(p)      * dt
         cs%pool_patch(p)               = cs%pool_patch(p)              - cf%m_pool_to_litter_patch(p)              * dt

         ! transfer pools
         cs%leaf_xfer_patch(p)          = cs%leaf_xfer_patch(p)         - cf%m_leaf_xfer_to_litter_patch(p)         * dt
         cs%froot_xfer_patch(p)         = cs%froot_xfer_patch(p)        - cf%m_froot_xfer_to_litter_patch(p)        * dt
         cs%livestem_xfer_patch(p)      = cs%livestem_xfer_patch(p)     - cf%m_livestem_xfer_to_litter_patch(p)     * dt
         cs%deadstem_xfer_patch(p)      = cs%deadstem_xfer_patch(p)     - cf%m_deadstem_xfer_to_litter_patch(p)     * dt
         cs%livecroot_xfer_patch(p)     = cs%livecroot_xfer_patch(p)    - cf%m_livecroot_xfer_to_litter_patch(p)    * dt
         cs%deadcroot_xfer_patch(p)     = cs%deadcroot_xfer_patch(p)    - cf%m_deadcroot_xfer_to_litter_patch(p)    * dt
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
         ivt => veg_pp%itype      , & ! Input:  [integer (:)]  pft vegetation type
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
               cs%decomp_pools_vr_col(c,j,i_met_lit) = &
                    cs%decomp_pools_vr_col(c,j,i_met_lit) + cf%harvest_to_litr_met_col(c,j) * dt
               cs%decomp_pools_vr_col(c,j,i_cel_lit) = &
                    cs%decomp_pools_vr_col(c,j,i_cel_lit) + cf%harvest_to_litr_cel_col(c,j) * dt
               cs%decomp_pools_vr_col(c,j,i_lig_lit) = &
                    cs%decomp_pools_vr_col(c,j,i_lig_lit) + cf%harvest_to_litr_lig_col(c,j) * dt
               cs%decomp_pools_vr_col(c,j,i_cwd) = &
                    cs%decomp_pools_vr_col(c,j,i_cwd) + cf%harvest_to_cwd_col(c,j)  * dt

               ! wood to product pools - states updated in CNWoodProducts()
            end do
         end do

      endif

      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! patch-level carbon fluxes from harvest mortality
         ! displayed pools
         cs%leaf_patch(p)               = cs%leaf_patch(p)              - cf%hrv_leaf_to_litter_patch(p)              * dt
         cs%froot_patch(p)              = cs%froot_patch(p)             - cf%hrv_froot_to_litter_patch(p)             * dt
         cs%livestem_patch(p)           = cs%livestem_patch(p)          - cf%hrv_livestem_to_litter_patch(p)          * dt
         cs%deadstem_patch(p)           = cs%deadstem_patch(p)          - cf%hrv_deadstem_to_prod10_patch(p)         * dt
         cs%deadstem_patch(p)           = cs%deadstem_patch(p)          - cf%hrv_deadstem_to_prod100_patch(p)        * dt
         cs%livecroot_patch(p)          = cs%livecroot_patch(p)         - cf%hrv_livecroot_to_litter_patch(p)         * dt
         cs%deadcroot_patch(p)          = cs%deadcroot_patch(p)         - cf%hrv_deadcroot_to_litter_patch(p)         * dt

         ! crops
         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
             cs%livestem_patch(p)       = cs%livestem_patch(p)          - cf%hrv_livestem_to_prod1_patch(p)          *dt
             cs%leaf_patch(p)           = cs%leaf_patch(p)              - cf%hrv_leaf_to_prod1_patch(p)              *dt
             cs%grain_patch(p)          = cs%grain_patch(p)             - cf%hrv_grain_to_prod1_patch(p)             *dt
         end if

         ! xsmrpool
         cs%xsmrpool_patch(p)            = cs%xsmrpool_patch(p)           - cf%hrv_xsmrpool_to_atm_patch(p)              * dt
         ! storage pools
         cs%leaf_storage_patch(p)       = cs%leaf_storage_patch(p)      - cf%hrv_leaf_storage_to_litter_patch(p)      * dt
         cs%froot_storage_patch(p)      = cs%froot_storage_patch(p)     - cf%hrv_froot_storage_to_litter_patch(p)     * dt
         cs%livestem_storage_patch(p)   = cs%livestem_storage_patch(p)  - cf%hrv_livestem_storage_to_litter_patch(p)  * dt
         cs%deadstem_storage_patch(p)   = cs%deadstem_storage_patch(p)  - cf%hrv_deadstem_storage_to_litter_patch(p)  * dt
         cs%livecroot_storage_patch(p)  = cs%livecroot_storage_patch(p) - cf%hrv_livecroot_storage_to_litter_patch(p) * dt
         cs%deadcroot_storage_patch(p)  = cs%deadcroot_storage_patch(p) - cf%hrv_deadcroot_storage_to_litter_patch(p) * dt
         cs%gresp_storage_patch(p)       = cs%gresp_storage_patch(p)      - cf%hrv_gresp_storage_to_litter_patch(p)      * dt
         cs%pool_patch(p)               = cs%pool_patch(p)              - cf%hrv_pool_to_litter_patch(p)              * dt

         ! transfer pools
         cs%leaf_xfer_patch(p)          = cs%leaf_xfer_patch(p)         - cf%hrv_leaf_xfer_to_litter_patch(p)         * dt
         cs%froot_xfer_patch(p)         = cs%froot_xfer_patch(p)        - cf%hrv_froot_xfer_to_litter_patch(p)        * dt
         cs%livestem_xfer_patch(p)      = cs%livestem_xfer_patch(p)     - cf%hrv_livestem_xfer_to_litter_patch(p)     * dt
         cs%deadstem_xfer_patch(p)      = cs%deadstem_xfer_patch(p)     - cf%hrv_deadstem_xfer_to_litter_patch(p)     * dt
         cs%livecroot_xfer_patch(p)     = cs%livecroot_xfer_patch(p)    - cf%hrv_livecroot_xfer_to_litter_patch(p)    * dt
         cs%deadcroot_xfer_patch(p)     = cs%deadcroot_xfer_patch(p)    - cf%hrv_deadcroot_xfer_to_litter_patch(p)    * dt
         cs%gresp_xfer_patch(p)          = cs%gresp_xfer_patch(p)         - cf%hrv_gresp_xfer_to_litter_patch(p)         * dt

      end do ! end of patch loop

    end associate

  end subroutine CStateUpdate2h

end module CNCStateUpdate2Mod
