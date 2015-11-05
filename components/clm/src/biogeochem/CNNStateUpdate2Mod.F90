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
               
               ns%decomp_npools_vr_col(c,j,i_met_lit) = &
                    ns%decomp_npools_vr_col(c,j,i_met_lit) + nf%gap_mortality_n_to_litr_met_n_col(c,j) * dt
               ns%decomp_npools_vr_col(c,j,i_cel_lit) = &
                    ns%decomp_npools_vr_col(c,j,i_cel_lit) + nf%gap_mortality_n_to_litr_cel_n_col(c,j) * dt
               ns%decomp_npools_vr_col(c,j,i_lig_lit) = &
                    ns%decomp_npools_vr_col(c,j,i_lig_lit) + nf%gap_mortality_n_to_litr_lig_n_col(c,j) * dt
               ns%decomp_npools_vr_col(c,j,i_cwd)     = &
                    ns%decomp_npools_vr_col(c,j,i_cwd)     + nf%gap_mortality_n_to_cwdn_col(c,j)       * dt
            end do
         end do

      elseif (is_active_betr_bgc) then

         do j = 1, nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               
               nf%bgc_npool_ext_inputs_vr_col(c,j,i_met_lit) = &
                    nf%bgc_npool_ext_inputs_vr_col(c,j,i_met_lit) + nf%gap_mortality_n_to_litr_met_n_col(c,j) * dt
               nf%bgc_npool_ext_inputs_vr_col(c,j,i_cel_lit) = &
                    nf%bgc_npool_ext_inputs_vr_col(c,j,i_cel_lit) + nf%gap_mortality_n_to_litr_cel_n_col(c,j) * dt
               nf%bgc_npool_ext_inputs_vr_col(c,j,i_lig_lit) = &
                    nf%bgc_npool_ext_inputs_vr_col(c,j,i_lig_lit) + nf%gap_mortality_n_to_litr_lig_n_col(c,j) * dt
               nf%bgc_npool_ext_inputs_vr_col(c,j,i_cwd)     = &
                    nf%bgc_npool_ext_inputs_vr_col(c,j,i_cwd)     + nf%gap_mortality_n_to_cwdn_col(c,j)       * dt
            end do
         end do         
     endif

      ! patch -level nitrogen fluxes from gap-phase mortality

      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! displayed pools
         ns%leafn_patch(p)              =  ns%leafn_patch(p)      - nf%m_leafn_to_litter_patch(p)      * dt
         ns%frootn_patch(p)             =  ns%frootn_patch(p)     - nf%m_frootn_to_litter_patch(p)     * dt
         ns%livestemn_patch(p)          =  ns%livestemn_patch(p)  - nf%m_livestemn_to_litter_patch(p)  * dt
         ns%deadstemn_patch(p)          =  ns%deadstemn_patch(p)  - nf%m_deadstemn_to_litter_patch(p)  * dt
         ns%livecrootn_patch(p)         =  ns%livecrootn_patch(p) - nf%m_livecrootn_to_litter_patch(p) * dt
         ns%deadcrootn_patch(p)         =  ns%deadcrootn_patch(p) - nf%m_deadcrootn_to_litter_patch(p) * dt
         ns%retransn_patch(p)           =  ns%retransn_patch(p)   - nf%m_retransn_to_litter_patch(p)   * dt

         ! storage pools
         ns%leafn_storage_patch(p)      =  ns%leafn_storage_patch(p)      - nf%m_leafn_storage_to_litter_patch(p)      * dt
         ns%frootn_storage_patch(p)     =  ns%frootn_storage_patch(p)     - nf%m_frootn_storage_to_litter_patch(p)     * dt
         ns%livestemn_storage_patch(p)  =  ns%livestemn_storage_patch(p)  - nf%m_livestemn_storage_to_litter_patch(p)  * dt
         ns%deadstemn_storage_patch(p)  =  ns%deadstemn_storage_patch(p)  - nf%m_deadstemn_storage_to_litter_patch(p)  * dt
         ns%livecrootn_storage_patch(p) =  ns%livecrootn_storage_patch(p) - nf%m_livecrootn_storage_to_litter_patch(p) * dt
         ns%deadcrootn_storage_patch(p) =  ns%deadcrootn_storage_patch(p) - nf%m_deadcrootn_storage_to_litter_patch(p) * dt

         ! transfer pools
         ns%leafn_xfer_patch(p)         =  ns%leafn_xfer_patch(p)      - nf%m_leafn_xfer_to_litter_patch(p)      * dt
         ns%frootn_xfer_patch(p)        =  ns%frootn_xfer_patch(p)     - nf%m_frootn_xfer_to_litter_patch(p)     * dt
         ns%livestemn_xfer_patch(p)     =  ns%livestemn_xfer_patch(p)  - nf%m_livestemn_xfer_to_litter_patch(p)  * dt
         ns%deadstemn_xfer_patch(p)     =  ns%deadstemn_xfer_patch(p)  - nf%m_deadstemn_xfer_to_litter_patch(p)  * dt
         ns%livecrootn_xfer_patch(p)    =  ns%livecrootn_xfer_patch(p) - nf%m_livecrootn_xfer_to_litter_patch(p) * dt
         ns%deadcrootn_xfer_patch(p)    =  ns%deadcrootn_xfer_patch(p) - nf%m_deadcrootn_xfer_to_litter_patch(p) * dt

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
         ivt => pft%itype         , & ! Input:  [integer  (:) ]  pft vegetation type
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
               ns%decomp_npools_vr_col(c,j,i_met_lit) = &
                    ns%decomp_npools_vr_col(c,j,i_met_lit) + nf%harvest_n_to_litr_met_n_col(c,j) * dt
               ns%decomp_npools_vr_col(c,j,i_cel_lit) = &
                    ns%decomp_npools_vr_col(c,j,i_cel_lit) + nf%harvest_n_to_litr_cel_n_col(c,j) * dt
               ns%decomp_npools_vr_col(c,j,i_lig_lit) = &
                    ns%decomp_npools_vr_col(c,j,i_lig_lit) + nf%harvest_n_to_litr_lig_n_col(c,j) * dt
               ns%decomp_npools_vr_col(c,j,i_cwd)     = &
                    ns%decomp_npools_vr_col(c,j,i_cwd)     + nf%harvest_n_to_cwdn_col(c,j)       * dt
            end do
         end do

      elseif (is_active_betr_bgc) then

         do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               nf%bgc_npool_ext_inputs_vr_col(c,j,i_met_lit) = &
                    nf%bgc_npool_ext_inputs_vr_col(c,j,i_met_lit) + nf%harvest_n_to_litr_met_n_col(c,j) * dt
               nf%bgc_npool_ext_inputs_vr_col(c,j,i_cel_lit) = &
                    nf%bgc_npool_ext_inputs_vr_col(c,j,i_cel_lit) + nf%harvest_n_to_litr_cel_n_col(c,j) * dt
               nf%bgc_npool_ext_inputs_vr_col(c,j,i_lig_lit) = &
                    nf%bgc_npool_ext_inputs_vr_col(c,j,i_lig_lit) + nf%harvest_n_to_litr_lig_n_col(c,j) * dt
               nf%bgc_npool_ext_inputs_vr_col(c,j,i_cwd)     = &
                    nf%bgc_npool_ext_inputs_vr_col(c,j,i_cwd)     + nf%harvest_n_to_cwdn_col(c,j)       * dt
            end do
         end do         
      endif

      ! patch-level nitrogen fluxes from harvest mortality

      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! displayed pools
         ns%leafn_patch(p)      = ns%leafn_patch(p)      - nf%hrv_leafn_to_litter_patch(p)      * dt
         ns%frootn_patch(p)     = ns%frootn_patch(p)     - nf%hrv_frootn_to_litter_patch(p)     * dt
         ns%livestemn_patch(p)  = ns%livestemn_patch(p)  - nf%hrv_livestemn_to_litter_patch(p)  * dt
         ns%deadstemn_patch(p)  = ns%deadstemn_patch(p)  - nf%hrv_deadstemn_to_prod10n_patch(p) * dt
         ns%deadstemn_patch(p)  = ns%deadstemn_patch(p)  - nf%hrv_deadstemn_to_prod100n_patch(p)* dt
         ns%livecrootn_patch(p) = ns%livecrootn_patch(p) - nf%hrv_livecrootn_to_litter_patch(p) * dt
         ns%deadcrootn_patch(p) = ns%deadcrootn_patch(p) - nf%hrv_deadcrootn_to_litter_patch(p) * dt
         ns%retransn_patch(p)   = ns%retransn_patch(p)   - nf%hrv_retransn_to_litter_patch(p)   * dt

       if (ivt(p) >= npcropmin) then ! skip 2 generic crops
           ns%livestemn_patch(p)= ns%livestemn_patch(p)  - nf%hrv_livestemn_to_prod1n_patch(p)  * dt
           ns%leafn_patch(p)    = ns%leafn_patch(p)      - nf%hrv_leafn_to_prod1n_patch(p)      * dt
           ns%grainn_patch(p)   = ns%grainn_patch(p)     - nf%hrv_grainn_to_prod1n_patch(p)     * dt
       end if

         ! storage pools
         ns%leafn_storage_patch(p)      = ns%leafn_storage_patch(p)      - nf%hrv_leafn_storage_to_litter_patch(p)      * dt
         ns%frootn_storage_patch(p)     = ns%frootn_storage_patch(p)     - nf%hrv_frootn_storage_to_litter_patch(p)     * dt
         ns%livestemn_storage_patch(p)  = ns%livestemn_storage_patch(p)  - nf%hrv_livestemn_storage_to_litter_patch(p)  * dt
         ns%deadstemn_storage_patch(p)  = ns%deadstemn_storage_patch(p)  - nf%hrv_deadstemn_storage_to_litter_patch(p)  * dt
         ns%livecrootn_storage_patch(p) = ns%livecrootn_storage_patch(p) - nf%hrv_livecrootn_storage_to_litter_patch(p) * dt
         ns%deadcrootn_storage_patch(p) = ns%deadcrootn_storage_patch(p) - nf%hrv_deadcrootn_storage_to_litter_patch(p) * dt

         ! transfer pools
         ns%leafn_xfer_patch(p)      = ns%leafn_xfer_patch(p)      - nf%hrv_leafn_xfer_to_litter_patch(p)      *dt
         ns%frootn_xfer_patch(p)     = ns%frootn_xfer_patch(p)     - nf%hrv_frootn_xfer_to_litter_patch(p)     *dt
         ns%livestemn_xfer_patch(p)  = ns%livestemn_xfer_patch(p)  - nf%hrv_livestemn_xfer_to_litter_patch(p)  *dt
         ns%deadstemn_xfer_patch(p)  = ns%deadstemn_xfer_patch(p)  - nf%hrv_deadstemn_xfer_to_litter_patch(p)  *dt
         ns%livecrootn_xfer_patch(p) = ns%livecrootn_xfer_patch(p) - nf%hrv_livecrootn_xfer_to_litter_patch(p) *dt
         ns%deadcrootn_xfer_patch(p) = ns%deadcrootn_xfer_patch(p) - nf%hrv_deadcrootn_xfer_to_litter_patch(p) *dt

      end do

    end associate

  end subroutine NStateUpdate2h

end module CNNStateUpdate2Mod
