module CNNStateUpdate2Mod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for nitrogen state variable update, mortality fluxes.
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use clm_time_manager    , only : get_step_size
  use clm_varpar          , only : nlevsoi, nlevdecomp
  use clm_varpar          , only : i_met_lit, i_cel_lit, i_lig_lit, i_cwd
  use clm_varctl          , only : iulog
  use CNNitrogenStateType , only : nitrogenstate_type
  use CNNitrogenFLuxType  , only : nitrogenflux_type
  use VegetationType           , only : veg_pp
  use pftvarcon           , only : npcropmin
  ! bgc interface & pflotran:
  use clm_varctl          , only : use_pflotran, pf_cmode
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: NStateUpdate2
  public:: NStateUpdate2h
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine NStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
       nitrogenflux_vars, nitrogenstate_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update all the prognostic nitrogen state
    ! variables affected by gap-phase mortality fluxes
    ! NOTE - associate statements have been removed where there are
    ! no science equations. This increases readability and maintainability
    !
    use tracer_varcon, only : is_active_betr_bgc      
    ! !ARGUMENTS:
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                  , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: c,p,j,l ! indices
    integer  :: fp,fc   ! lake filter indices
    real(r8) :: dt      ! radiation time step (seconds)
    !-----------------------------------------------------------------------

    associate(                      & 
         nf => nitrogenflux_vars  , &
         ns => nitrogenstate_vars   &
         )

      ! set time steps
      dt = real( get_step_size(), r8 )


      ! column-level nitrogen fluxes from gap-phase mortality
      if ( .not. is_active_betr_bgc .and. &
           .not.(use_pflotran .and. pf_cmode)) then
         do j = 1, nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               
               ns%decomp_pools_vr_col(c,j,i_met_lit) = &
                    ns%decomp_pools_vr_col(c,j,i_met_lit) + nf%gap_mortality_to_litr_met_col(c,j) * dt
               ns%decomp_pools_vr_col(c,j,i_cel_lit) = &
                    ns%decomp_pools_vr_col(c,j,i_cel_lit) + nf%gap_mortality_to_litr_cel_col(c,j) * dt
               ns%decomp_pools_vr_col(c,j,i_lig_lit) = &
                    ns%decomp_pools_vr_col(c,j,i_lig_lit) + nf%gap_mortality_to_litr_lig_col(c,j) * dt
               ns%decomp_pools_vr_col(c,j,i_cwd)     = &
                    ns%decomp_pools_vr_col(c,j,i_cwd)     + nf%gap_mortality_to_cwd_col(c,j)       * dt
            end do
         end do

     endif

      ! patch -level nitrogen fluxes from gap-phase mortality

      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! displayed pools
         ns%leaf_patch(p)              =  ns%leaf_patch(p)      - nf%m_leaf_to_litter_patch(p)      * dt
         ns%froot_patch(p)             =  ns%froot_patch(p)     - nf%m_froot_to_litter_patch(p)     * dt
         ns%livestem_patch(p)          =  ns%livestem_patch(p)  - nf%m_livestem_to_litter_patch(p)  * dt
         ns%deadstem_patch(p)          =  ns%deadstem_patch(p)  - nf%m_deadstem_to_litter_patch(p)  * dt
         ns%livecroot_patch(p)         =  ns%livecroot_patch(p) - nf%m_livecroot_to_litter_patch(p) * dt
         ns%deadcroot_patch(p)         =  ns%deadcroot_patch(p) - nf%m_deadcroot_to_litter_patch(p) * dt
         ns%retransn_patch(p)           =  ns%retransn_patch(p)   - nf%m_retransn_to_litter_patch(p)   * dt
         ns%pool_patch(p)              =  ns%pool_patch(p)      - nf%m_pool_to_litter_patch(p)      * dt

         ! storage pools
         ns%leaf_storage_patch(p)      =  ns%leaf_storage_patch(p)      - nf%m_leaf_storage_to_litter_patch(p)      * dt
         ns%froot_storage_patch(p)     =  ns%froot_storage_patch(p)     - nf%m_froot_storage_to_litter_patch(p)     * dt
         ns%livestem_storage_patch(p)  =  ns%livestem_storage_patch(p)  - nf%m_livestem_storage_to_litter_patch(p)  * dt
         ns%deadstem_storage_patch(p)  =  ns%deadstem_storage_patch(p)  - nf%m_deadstem_storage_to_litter_patch(p)  * dt
         ns%livecroot_storage_patch(p) =  ns%livecroot_storage_patch(p) - nf%m_livecroot_storage_to_litter_patch(p) * dt
         ns%deadcroot_storage_patch(p) =  ns%deadcroot_storage_patch(p) - nf%m_deadcroot_storage_to_litter_patch(p) * dt

         ! transfer pools
         ns%leaf_xfer_patch(p)         =  ns%leaf_xfer_patch(p)      - nf%m_leaf_xfer_to_litter_patch(p)      * dt
         ns%froot_xfer_patch(p)        =  ns%froot_xfer_patch(p)     - nf%m_froot_xfer_to_litter_patch(p)     * dt
         ns%livestem_xfer_patch(p)     =  ns%livestem_xfer_patch(p)  - nf%m_livestem_xfer_to_litter_patch(p)  * dt
         ns%deadstem_xfer_patch(p)     =  ns%deadstem_xfer_patch(p)  - nf%m_deadstem_xfer_to_litter_patch(p)  * dt
         ns%livecroot_xfer_patch(p)    =  ns%livecroot_xfer_patch(p) - nf%m_livecroot_xfer_to_litter_patch(p) * dt
         ns%deadcroot_xfer_patch(p)    =  ns%deadcroot_xfer_patch(p) - nf%m_deadcroot_xfer_to_litter_patch(p) * dt

      end do

    end associate

  end subroutine NStateUpdate2

  !-----------------------------------------------------------------------
  subroutine NStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
       nitrogenflux_vars, nitrogenstate_vars)
    !
    ! !DESCRIPTION:
    ! Update all the prognostic nitrogen state
    ! variables affected by harvest mortality fluxes
    ! NOTE - associate statements have been removed where there are
    ! no science equations. This increases readability and maintainability
    !
    use tracer_varcon, only : is_active_betr_bgc      
    ! !ARGUMENTS:
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                  , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    !
    ! !LOCAL VARIABLES:
    integer :: c,p,j,l ! indices
    integer :: fp,fc   ! lake filter indices
    real(r8):: dt      ! radiation time step (seconds)
    !-----------------------------------------------------------------------

    associate(                      & 
         ivt => veg_pp%itype         , & ! Input:  [integer  (:) ]  pft vegetation type
         nf => nitrogenflux_vars  , &
         ns => nitrogenstate_vars   &
         )

      ! set time steps
      dt = real( get_step_size(), r8 )

      if (.not. is_active_betr_bgc .and. &
           .not.(use_pflotran .and. pf_cmode)) then
         ! column-level nitrogen fluxes from harvest mortality

         do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               ns%decomp_pools_vr_col(c,j,i_met_lit) = &
                    ns%decomp_pools_vr_col(c,j,i_met_lit) + nf%harvest_to_litr_met_col(c,j) * dt
               ns%decomp_pools_vr_col(c,j,i_cel_lit) = &
                    ns%decomp_pools_vr_col(c,j,i_cel_lit) + nf%harvest_to_litr_cel_col(c,j) * dt
               ns%decomp_pools_vr_col(c,j,i_lig_lit) = &
                    ns%decomp_pools_vr_col(c,j,i_lig_lit) + nf%harvest_to_litr_lig_col(c,j) * dt
               ns%decomp_pools_vr_col(c,j,i_cwd)     = &
                    ns%decomp_pools_vr_col(c,j,i_cwd)     + nf%harvest_to_cwd_col(c,j)       * dt
            end do
         end do

      endif

      ! patch-level nitrogen fluxes from harvest mortality

      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! displayed pools
         ns%leaf_patch(p)      = ns%leaf_patch(p)      - nf%hrv_leaf_to_litter_patch(p)      * dt
         ns%froot_patch(p)     = ns%froot_patch(p)     - nf%hrv_froot_to_litter_patch(p)     * dt
         ns%livestem_patch(p)  = ns%livestem_patch(p)  - nf%hrv_livestem_to_litter_patch(p)  * dt
         ns%deadstem_patch(p)  = ns%deadstem_patch(p)  - nf%hrv_deadstem_to_prod10_patch(p) * dt
         ns%deadstem_patch(p)  = ns%deadstem_patch(p)  - nf%hrv_deadstem_to_prod100_patch(p)* dt
         ns%livecroot_patch(p) = ns%livecroot_patch(p) - nf%hrv_livecroot_to_litter_patch(p) * dt
         ns%deadcroot_patch(p) = ns%deadcroot_patch(p) - nf%hrv_deadcroot_to_litter_patch(p) * dt
         ns%retransn_patch(p)   = ns%retransn_patch(p)   - nf%hrv_retransn_to_litter_patch(p)   * dt
         ns%pool_patch(p)      = ns%pool_patch(p)      - nf%hrv_pool_to_litter_patch(p)     * dt

       if (ivt(p) >= npcropmin) then ! skip 2 generic crops
           ns%livestem_patch(p)= ns%livestem_patch(p)  - nf%hrv_livestem_to_prod1_patch(p)  * dt
           ns%leaf_patch(p)    = ns%leaf_patch(p)      - nf%hrv_leaf_to_prod1_patch(p)      * dt
           ns%grain_patch(p)   = ns%grain_patch(p)     - nf%hrv_grain_to_prod1_patch(p)     * dt
       end if

         ! storage pools
         ns%leaf_storage_patch(p)      = ns%leaf_storage_patch(p)      - nf%hrv_leaf_storage_to_litter_patch(p)      * dt
         ns%froot_storage_patch(p)     = ns%froot_storage_patch(p)     - nf%hrv_froot_storage_to_litter_patch(p)     * dt
         ns%livestem_storage_patch(p)  = ns%livestem_storage_patch(p)  - nf%hrv_livestem_storage_to_litter_patch(p)  * dt
         ns%deadstem_storage_patch(p)  = ns%deadstem_storage_patch(p)  - nf%hrv_deadstem_storage_to_litter_patch(p)  * dt
         ns%livecroot_storage_patch(p) = ns%livecroot_storage_patch(p) - nf%hrv_livecroot_storage_to_litter_patch(p) * dt
         ns%deadcroot_storage_patch(p) = ns%deadcroot_storage_patch(p) - nf%hrv_deadcroot_storage_to_litter_patch(p) * dt

         ! transfer pools
         ns%leaf_xfer_patch(p)      = ns%leaf_xfer_patch(p)      - nf%hrv_leaf_xfer_to_litter_patch(p)      *dt
         ns%froot_xfer_patch(p)     = ns%froot_xfer_patch(p)     - nf%hrv_froot_xfer_to_litter_patch(p)     *dt
         ns%livestem_xfer_patch(p)  = ns%livestem_xfer_patch(p)  - nf%hrv_livestem_xfer_to_litter_patch(p)  *dt
         ns%deadstem_xfer_patch(p)  = ns%deadstem_xfer_patch(p)  - nf%hrv_deadstem_xfer_to_litter_patch(p)  *dt
         ns%livecroot_xfer_patch(p) = ns%livecroot_xfer_patch(p) - nf%hrv_livecroot_xfer_to_litter_patch(p) *dt
         ns%deadcroot_xfer_patch(p) = ns%deadcroot_xfer_patch(p) - nf%hrv_deadcroot_xfer_to_litter_patch(p) *dt

      end do

    end associate

  end subroutine NStateUpdate2h

end module CNNStateUpdate2Mod
