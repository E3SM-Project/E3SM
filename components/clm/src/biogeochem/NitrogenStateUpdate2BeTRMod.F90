module NitrogenStateUpdate2BeTRMod

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
  use ColumnDataType          , only : column_nitrogen_state, column_nitrogen_flux
  use VegetationDataType      , only : vegetation_nitrogen_state, vegetation_nitrogen_flux
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: NitrogenStateUpdate2
  public:: NitrogenStateUpdate2h
  public:: NitrogenStateUpdate2Soil
  public:: NitrogenStateUpdate2hSoil
  public:: NitrogenStateUpdate2Veg
  public:: NitrogenStateUpdate2hVeg


  !-----------------------------------------------------------------------

contains
  !-----------------------------------------------------------------------
  subroutine NitrogenStateUpdate2Soil(num_soilc, filter_soilc, num_soilp, filter_soilp, &
       col_ns, col_nf)
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
    type(column_nitrogen_state)    , intent(inout) :: col_ns
    type(column_nitrogen_flux)     , intent(inout) :: col_nf
    !
    ! !LOCAL VARIABLES:
    integer  :: c,p,j,l ! indices
    integer  :: fp,fc   ! lake filter indices
    real(r8) :: dt      ! radiation time step (seconds)
    !-----------------------------------------------------------------------


      ! set time steps
      dt = real( get_step_size(), r8 )

      ! column-level nitrogen fluxes from gap-phase mortality
         do j = 1, nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)

              col_ns%decomp_npools_vr(c,j,i_met_lit) = &
                   col_ns%decomp_npools_vr(c,j,i_met_lit) +col_nf%gap_mortality_n_to_litr_met_n(c,j) * dt
              col_ns%decomp_npools_vr(c,j,i_cel_lit) = &
                   col_ns%decomp_npools_vr(c,j,i_cel_lit) +col_nf%gap_mortality_n_to_litr_cel_n(c,j) * dt
              col_ns%decomp_npools_vr(c,j,i_lig_lit) = &
                   col_ns%decomp_npools_vr(c,j,i_lig_lit) +col_nf%gap_mortality_n_to_litr_lig_n(c,j) * dt
              col_ns%decomp_npools_vr(c,j,i_cwd)     = &
                   col_ns%decomp_npools_vr(c,j,i_cwd)     +col_nf%gap_mortality_n_to_cwdn(c,j)       * dt
            end do
         end do

  end subroutine NitrogenStateUpdate2Soil
  !-----------------------------------------------------------------------
  subroutine NitrogenStateUpdate2Veg(num_soilc, filter_soilc, num_soilp, filter_soilp, &
       veg_ns, veg_nf)
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
    type(vegetation_nitrogen_state),intent(inout) :: veg_ns
    type(vegetation_nitrogen_flux) ,intent(inout) :: veg_nf    !
    !
    ! !LOCAL VARIABLES:
    integer  :: c,p,j,l ! indices
    integer  :: fp,fc   ! lake filter indices
    real(r8) :: dt      ! radiation time step (seconds)
    !-----------------------------------------------------------------------


      ! set time steps
      dt = real( get_step_size(), r8 )

      ! patch -level nitrogen fluxes from gap-phase mortality

      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! displayed pools
        veg_ns%leafn(p)              = veg_ns%leafn(p)      -veg_nf%m_leafn_to_litter(p)      * dt
        veg_ns%frootn(p)             = veg_ns%frootn(p)     -veg_nf%m_frootn_to_litter(p)     * dt
        veg_ns%livestemn(p)          = veg_ns%livestemn(p)  -veg_nf%m_livestemn_to_litter(p)  * dt
        veg_ns%deadstemn(p)          = veg_ns%deadstemn(p)  -veg_nf%m_deadstemn_to_litter(p)  * dt
        veg_ns%livecrootn(p)         = veg_ns%livecrootn(p) -veg_nf%m_livecrootn_to_litter(p) * dt
        veg_ns%deadcrootn(p)         = veg_ns%deadcrootn(p) -veg_nf%m_deadcrootn_to_litter(p) * dt
        veg_ns%retransn(p)           = veg_ns%retransn(p)   -veg_nf%m_retransn_to_litter(p)   * dt
        veg_ns%npool(p)              = veg_ns%npool(p)      -veg_nf%m_npool_to_litter(p)      * dt

         ! storage pools
        veg_ns%leafn_storage(p)      = veg_ns%leafn_storage(p)      -veg_nf%m_leafn_storage_to_litter(p)      * dt
        veg_ns%frootn_storage(p)     = veg_ns%frootn_storage(p)     -veg_nf%m_frootn_storage_to_litter(p)     * dt
        veg_ns%livestemn_storage(p)  = veg_ns%livestemn_storage(p)  -veg_nf%m_livestemn_storage_to_litter(p)  * dt
        veg_ns%deadstemn_storage(p)  = veg_ns%deadstemn_storage(p)  -veg_nf%m_deadstemn_storage_to_litter(p)  * dt
        veg_ns%livecrootn_storage(p) = veg_ns%livecrootn_storage(p) -veg_nf%m_livecrootn_storage_to_litter(p) * dt
        veg_ns%deadcrootn_storage(p) = veg_ns%deadcrootn_storage(p) -veg_nf%m_deadcrootn_storage_to_litter(p) * dt

         ! transfer pools
        veg_ns%leafn_xfer(p)         = veg_ns%leafn_xfer(p)      -veg_nf%m_leafn_xfer_to_litter(p)      * dt
        veg_ns%frootn_xfer(p)        = veg_ns%frootn_xfer(p)     -veg_nf%m_frootn_xfer_to_litter(p)     * dt
        veg_ns%livestemn_xfer(p)     = veg_ns%livestemn_xfer(p)  -veg_nf%m_livestemn_xfer_to_litter(p)  * dt
        veg_ns%deadstemn_xfer(p)     = veg_ns%deadstemn_xfer(p)  -veg_nf%m_deadstemn_xfer_to_litter(p)  * dt
        veg_ns%livecrootn_xfer(p)    = veg_ns%livecrootn_xfer(p) -veg_nf%m_livecrootn_xfer_to_litter(p) * dt
        veg_ns%deadcrootn_xfer(p)    = veg_ns%deadcrootn_xfer(p) -veg_nf%m_deadcrootn_xfer_to_litter(p) * dt

      end do


  end subroutine NitrogenStateUpdate2Veg
  !-----------------------------------------------------------------------
  subroutine NitrogenStateUpdate2hSoil(num_soilc, filter_soilc, num_soilp, filter_soilp, &
       col_ns, col_nf)
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
    type(column_nitrogen_state)    , intent(inout) :: col_ns
    type(column_nitrogen_flux)     , intent(inout) :: col_nf
    !
    ! !LOCAL VARIABLES:
    integer :: c,p,j,l ! indices
    integer :: fp,fc   ! lake filter indices
    real(r8):: dt      ! radiation time step (seconds)
    !-----------------------------------------------------------------------

    associate(                         &
         ivt => veg_pp%itype           & ! Input:  [integer  (:) ]  pft vegetation type

         )

      ! set time steps
      dt = real( get_step_size(), r8 )

         do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
              col_ns%decomp_npools_vr(c,j,i_met_lit) = &
                   col_ns%decomp_npools_vr(c,j,i_met_lit) +col_nf%harvest_n_to_litr_met_n(c,j) * dt
              col_ns%decomp_npools_vr(c,j,i_cel_lit) = &
                   col_ns%decomp_npools_vr(c,j,i_cel_lit) +col_nf%harvest_n_to_litr_cel_n(c,j) * dt
              col_ns%decomp_npools_vr(c,j,i_lig_lit) = &
                   col_ns%decomp_npools_vr(c,j,i_lig_lit) +col_nf%harvest_n_to_litr_lig_n(c,j) * dt
              col_ns%decomp_npools_vr(c,j,i_cwd)     = &
                   col_ns%decomp_npools_vr(c,j,i_cwd)     +col_nf%harvest_n_to_cwdn(c,j)       * dt
            end do
         end do
    end associate
  end subroutine NitrogenStateUpdate2hSoil
  !-----------------------------------------------------------------------
  subroutine NitrogenStateUpdate2hVeg(num_soilc, filter_soilc, num_soilp, filter_soilp, &
       veg_ns, veg_nf)
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
    type(vegetation_nitrogen_state),intent(inout) :: veg_ns
    type(vegetation_nitrogen_flux) ,intent(inout) :: veg_nf    !
    !
    ! !LOCAL VARIABLES:
    integer :: c,p,j,l ! indices
    integer :: fp,fc   ! lake filter indices
    real(r8):: dt      ! radiation time step (seconds)
    !-----------------------------------------------------------------------

    associate(                        &
         ivt => veg_pp%itype          & ! Input:  [integer  (:) ]  pft vegetation type

         )

      ! set time steps
      dt = real( get_step_size(), r8 )



      ! patch-level nitrogen fluxes from harvest mortality

      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! displayed pools
        veg_ns%leafn(p)      =veg_ns%leafn(p)      -veg_nf%hrv_leafn_to_litter(p)      * dt
        veg_ns%frootn(p)     =veg_ns%frootn(p)     -veg_nf%hrv_frootn_to_litter(p)     * dt
        veg_ns%livestemn(p)  =veg_ns%livestemn(p)  -veg_nf%hrv_livestemn_to_litter(p)  * dt
        veg_ns%deadstemn(p)  =veg_ns%deadstemn(p)  -veg_nf%hrv_deadstemn_to_prod10n(p) * dt
        veg_ns%deadstemn(p)  =veg_ns%deadstemn(p)  -veg_nf%hrv_deadstemn_to_prod100n(p)* dt
        veg_ns%livecrootn(p) =veg_ns%livecrootn(p) -veg_nf%hrv_livecrootn_to_litter(p) * dt
        veg_ns%deadcrootn(p) =veg_ns%deadcrootn(p) -veg_nf%hrv_deadcrootn_to_litter(p) * dt
        veg_ns%retransn(p)   =veg_ns%retransn(p)   -veg_nf%hrv_retransn_to_litter(p)   * dt
        veg_ns%npool(p)      =veg_ns%npool(p)      -veg_nf%hrv_npool_to_litter(p)     * dt

       if (ivt(p) >= npcropmin) then ! skip 2 generic crops
          veg_ns%livestemn(p)=veg_ns%livestemn(p)  -veg_nf%hrv_livestemn_to_prod1n(p)  * dt
          veg_ns%leafn(p)    =veg_ns%leafn(p)      -veg_nf%hrv_leafn_to_prod1n(p)      * dt
          veg_ns%grainn(p)   =veg_ns%grainn(p)     -veg_nf%hrv_grainn_to_prod1n(p)     * dt
       end if

         ! storage pools
        veg_ns%leafn_storage(p)      =veg_ns%leafn_storage(p)      -veg_nf%hrv_leafn_storage_to_litter(p)      * dt
        veg_ns%frootn_storage(p)     =veg_ns%frootn_storage(p)     -veg_nf%hrv_frootn_storage_to_litter(p)     * dt
        veg_ns%livestemn_storage(p)  =veg_ns%livestemn_storage(p)  -veg_nf%hrv_livestemn_storage_to_litter(p)  * dt
        veg_ns%deadstemn_storage(p)  =veg_ns%deadstemn_storage(p)  -veg_nf%hrv_deadstemn_storage_to_litter(p)  * dt
        veg_ns%livecrootn_storage(p) =veg_ns%livecrootn_storage(p) -veg_nf%hrv_livecrootn_storage_to_litter(p) * dt
        veg_ns%deadcrootn_storage(p) =veg_ns%deadcrootn_storage(p) -veg_nf%hrv_deadcrootn_storage_to_litter(p) * dt

         ! transfer pools
        veg_ns%leafn_xfer(p)      =veg_ns%leafn_xfer(p)      -veg_nf%hrv_leafn_xfer_to_litter(p)      *dt
        veg_ns%frootn_xfer(p)     =veg_ns%frootn_xfer(p)     -veg_nf%hrv_frootn_xfer_to_litter(p)     *dt
        veg_ns%livestemn_xfer(p)  =veg_ns%livestemn_xfer(p)  -veg_nf%hrv_livestemn_xfer_to_litter(p)  *dt
        veg_ns%deadstemn_xfer(p)  =veg_ns%deadstemn_xfer(p)  -veg_nf%hrv_deadstemn_xfer_to_litter(p)  *dt
        veg_ns%livecrootn_xfer(p) =veg_ns%livecrootn_xfer(p) -veg_nf%hrv_livecrootn_xfer_to_litter(p) *dt
        veg_ns%deadcrootn_xfer(p) =veg_ns%deadcrootn_xfer(p) -veg_nf%hrv_deadcrootn_xfer_to_litter(p) *dt

      end do

    end associate

  end subroutine NitrogenStateUpdate2hVeg

  !-----------------------------------------------------------------------
  subroutine NitrogenStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
       col_ns, veg_ns, col_nf, veg_nf)
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
    type(column_nitrogen_state),intent(inout):: col_ns
    type(vegetation_nitrogen_state),intent(inout) :: veg_ns
    type(column_nitrogen_flux)     ,intent(inout) :: col_nf
    type(vegetation_nitrogen_flux) ,intent(inout) :: veg_nf


  call NitrogenStateUpdate2Soil(num_soilc, filter_soilc, num_soilp, filter_soilp, &
       col_ns, col_nf)

  call NitrogenStateUpdate2Veg(num_soilc, filter_soilc, num_soilp, filter_soilp, &
     veg_ns, veg_nf)

  end subroutine NitrogenStateUpdate2


  !-----------------------------------------------------------------------
  subroutine NitrogenStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
       col_ns, veg_ns, col_nf, veg_nf)
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
    type(column_nitrogen_state),intent(inout):: col_ns
    type(vegetation_nitrogen_state),intent(inout) :: veg_ns
    type(column_nitrogen_flux)     ,intent(inout) :: col_nf
    type(vegetation_nitrogen_flux) ,intent(inout) :: veg_nf

    call NitrogenStateUpdate2hSoil(num_soilc, filter_soilc, num_soilp, filter_soilp, &
       col_ns, col_nf)

    call NitrogenStateUpdate2hVeg(num_soilc, filter_soilc, num_soilp, filter_soilp, &
       veg_ns, veg_nf)

  end subroutine NitrogenStateUpdate2h
end module NitrogenStateUpdate2BeTRMod
