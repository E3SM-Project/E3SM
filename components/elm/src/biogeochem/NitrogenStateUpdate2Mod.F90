module NitrogenStateUpdate2Mod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for nitrogen state variable update, mortality fluxes.
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use elm_varpar          , only : nlevsoi, nlevdecomp
  use elm_varpar          , only : i_met_lit, i_cel_lit, i_lig_lit, i_cwd
  use elm_varctl          , only : iulog
  use ColumnDataType      , only : col_ns, col_nf
  use VegetationType      , only : veg_pp
  use VegetationDataType  , only : veg_ns, veg_nf
  use pftvarcon           , only : iscft
  ! bgc interface & pflotran:
  use elm_varctl          , only : use_pflotran, pf_cmode
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: NitrogenStateUpdate2
  public:: NitrogenStateUpdate2h
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine NitrogenStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp,dt)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update all the prognostic nitrogen state
    ! variables affected by gap-phase mortality fluxes
    ! NOTE - associate statements have been removed where there are
    ! no science equations. This increases readability and maintainability
    !
    ! !ARGUMENTS:
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                  , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:) ! filter for soil patches
    real(r8) , intent(in)   :: dt      ! radiation time step (seconds)

    !
    ! !LOCAL VARIABLES:
    integer  :: c,p,j,l ! indices
    integer  :: fp,fc   ! lake filter indices
    !-----------------------------------------------------------------------
      
    ! column-level nitrogen fluxes from gap-phase mortality
      if (.not.(use_pflotran .and. pf_cmode)) then
         do j = 1, nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)

               col_ns%core%decomp_pools_vr(c,j,i_met_lit) = &
                    col_ns%core%decomp_pools_vr(c,j,i_met_lit) + col_nf%core%gap_mortality_to_litr_met(c,j) * dt
               col_ns%core%decomp_pools_vr(c,j,i_cel_lit) = &
                    col_ns%core%decomp_pools_vr(c,j,i_cel_lit) + col_nf%core%gap_mortality_to_litr_cel(c,j) * dt
               col_ns%core%decomp_pools_vr(c,j,i_lig_lit) = &
                    col_ns%core%decomp_pools_vr(c,j,i_lig_lit) + col_nf%core%gap_mortality_to_litr_lig(c,j) * dt
               col_ns%core%decomp_pools_vr(c,j,i_cwd)     = &
                    col_ns%core%decomp_pools_vr(c,j,i_cwd)     + col_nf%core%gap_mortality_to_cwd(c,j)       * dt
            end do
         end do

     endif

      ! patch -level nitrogen fluxes from gap-phase mortality

      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! displayed pools
         veg_ns%core%leaf(p)              =  veg_ns%core%leaf(p)      - veg_nf%core%m_leaf_to_litter(p)      * dt
         veg_ns%core%froot(p)             =  veg_ns%core%froot(p)     - veg_nf%core%m_froot_to_litter(p)     * dt
         veg_ns%core%livestem(p)          =  veg_ns%core%livestem(p)  - veg_nf%core%m_livestem_to_litter(p)  * dt
         veg_ns%core%deadstem(p)          =  veg_ns%core%deadstem(p)  - veg_nf%core%m_deadstem_to_litter(p)  * dt
         veg_ns%core%livecroot(p)         =  veg_ns%core%livecroot(p) - veg_nf%core%m_livecroot_to_litter(p) * dt
         veg_ns%core%deadcroot(p)         =  veg_ns%core%deadcroot(p) - veg_nf%core%m_deadcroot_to_litter(p) * dt
         veg_ns%retrans%pool(p)           =  veg_ns%retrans%pool(p)   - veg_nf%retrans%m_retrans_to_litter(p)   * dt
         veg_ns%core%pool(p)              =  veg_ns%core%pool(p)      - veg_nf%core%m_pool_to_litter(p)      * dt

         ! storage pools
         veg_ns%core%leaf_storage(p)      =  veg_ns%core%leaf_storage(p)      - veg_nf%core%m_leaf_storage_to_litter(p)      * dt
         veg_ns%core%froot_storage(p)     =  veg_ns%core%froot_storage(p)     - veg_nf%core%m_froot_storage_to_litter(p)     * dt
         veg_ns%core%livestem_storage(p)  =  veg_ns%core%livestem_storage(p)  - veg_nf%core%m_livestem_storage_to_litter(p)  * dt
         veg_ns%core%deadstem_storage(p)  =  veg_ns%core%deadstem_storage(p)  - veg_nf%core%m_deadstem_storage_to_litter(p)  * dt
         veg_ns%core%livecroot_storage(p) =  veg_ns%core%livecroot_storage(p) - veg_nf%core%m_livecroot_storage_to_litter(p) * dt
         veg_ns%core%deadcroot_storage(p) =  veg_ns%core%deadcroot_storage(p) - veg_nf%core%m_deadcroot_storage_to_litter(p) * dt

         ! transfer pools
         veg_ns%core%leaf_xfer(p)         =  veg_ns%core%leaf_xfer(p)      - veg_nf%core%m_leaf_xfer_to_litter(p)      * dt
         veg_ns%core%froot_xfer(p)        =  veg_ns%core%froot_xfer(p)     - veg_nf%core%m_froot_xfer_to_litter(p)     * dt
         veg_ns%core%livestem_xfer(p)     =  veg_ns%core%livestem_xfer(p)  - veg_nf%core%m_livestem_xfer_to_litter(p)  * dt
         veg_ns%core%deadstem_xfer(p)     =  veg_ns%core%deadstem_xfer(p)  - veg_nf%core%m_deadstem_xfer_to_litter(p)  * dt
         veg_ns%core%livecroot_xfer(p)    =  veg_ns%core%livecroot_xfer(p) - veg_nf%core%m_livecroot_xfer_to_litter(p) * dt
         veg_ns%core%deadcroot_xfer(p)    =  veg_ns%core%deadcroot_xfer(p) - veg_nf%core%m_deadcroot_xfer_to_litter(p) * dt

      end do

  end subroutine NitrogenStateUpdate2

  !-----------------------------------------------------------------------
  subroutine NitrogenStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, dt)
    !
    ! !DESCRIPTION:
    ! Update all the prognostic nitrogen state
    ! variables affected by harvest mortality fluxes
    ! NOTE - associate statements have been removed where there are
    ! no science equations. This increases readability and maintainability
    !
    ! !ARGUMENTS:
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                  , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:) ! filter for soil patches
    real(r8)                 , intent(in)    :: dt      ! radiation time step (seconds)

    !
    ! !LOCAL VARIABLES:
    integer :: c,p,j,l ! indices
    integer :: fp,fc   ! lake filter indices
    !-----------------------------------------------------------------------

    associate(                      &
         ivt => veg_pp%itype        & ! Input:  [integer  (:) ]  pft vegetation type
         )

      if (.not.(use_pflotran .and. pf_cmode)) then
         ! column-level nitrogen fluxes from harvest mortality

         do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               col_ns%core%decomp_pools_vr(c,j,i_met_lit) = &
                    col_ns%core%decomp_pools_vr(c,j,i_met_lit) + col_nf%core%harvest_to_litr_met(c,j) * dt
               col_ns%core%decomp_pools_vr(c,j,i_cel_lit) = &
                    col_ns%core%decomp_pools_vr(c,j,i_cel_lit) + col_nf%core%harvest_to_litr_cel(c,j) * dt
               col_ns%core%decomp_pools_vr(c,j,i_lig_lit) = &
                    col_ns%core%decomp_pools_vr(c,j,i_lig_lit) + col_nf%core%harvest_to_litr_lig(c,j) * dt
               col_ns%core%decomp_pools_vr(c,j,i_cwd)     = &
                    col_ns%core%decomp_pools_vr(c,j,i_cwd)     + col_nf%core%harvest_to_cwd(c,j)       * dt
            end do
         end do

      endif

      ! patch-level nitrogen fluxes from harvest mortality

      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! displayed pools
         veg_ns%core%leaf(p)      = veg_ns%core%leaf(p)      - veg_nf%core%hrv_leaf_to_litter(p)      * dt
         veg_ns%core%froot(p)     = veg_ns%core%froot(p)     - veg_nf%core%hrv_froot_to_litter(p)     * dt
         veg_ns%core%livestem(p)  = veg_ns%core%livestem(p)  - veg_nf%core%hrv_livestem_to_litter(p)  * dt
         veg_ns%core%deadstem(p)  = veg_ns%core%deadstem(p)  - veg_nf%core%hrv_deadstem_to_prod10(p) * dt
         veg_ns%core%deadstem(p)  = veg_ns%core%deadstem(p)  - veg_nf%core%hrv_deadstem_to_prod100(p)* dt
         veg_ns%core%livecroot(p) = veg_ns%core%livecroot(p) - veg_nf%core%hrv_livecroot_to_litter(p) * dt
         veg_ns%core%deadcroot(p) = veg_ns%core%deadcroot(p) - veg_nf%core%hrv_deadcroot_to_litter(p) * dt
         veg_ns%retrans%pool(p)   = veg_ns%retrans%pool(p)   - veg_nf%retrans%hrv_retrans_to_litter(p)   * dt
         veg_ns%core%pool(p)      = veg_ns%core%pool(p)      - veg_nf%core%hrv_pool_to_litter(p)     * dt

       if (iscft(ivt(p))) then ! skip 2 generic crops
           veg_ns%core%livestem(p)= veg_ns%core%livestem(p)  - veg_nf%core%hrv_livestem_to_prod1(p)  * dt
           veg_ns%core%leaf(p)    = veg_ns%core%leaf(p)      - veg_nf%core%hrv_leaf_to_prod1(p)      * dt
           veg_ns%core%grain(p)   = veg_ns%core%grain(p)     - veg_nf%core%hrv_grain_to_prod1(p)     * dt
       end if

         ! storage pools
         veg_ns%core%leaf_storage(p)      = veg_ns%core%leaf_storage(p)      - veg_nf%core%hrv_leaf_storage_to_litter(p)      * dt
         veg_ns%core%froot_storage(p)     = veg_ns%core%froot_storage(p)     - veg_nf%core%hrv_froot_storage_to_litter(p)     * dt
         veg_ns%core%livestem_storage(p)  = veg_ns%core%livestem_storage(p)  - veg_nf%core%hrv_livestem_storage_to_litter(p)  * dt
         veg_ns%core%deadstem_storage(p)  = veg_ns%core%deadstem_storage(p)  - veg_nf%core%hrv_deadstem_storage_to_litter(p)  * dt
         veg_ns%core%livecroot_storage(p) = veg_ns%core%livecroot_storage(p) - veg_nf%core%hrv_livecroot_storage_to_litter(p) * dt
         veg_ns%core%deadcroot_storage(p) = veg_ns%core%deadcroot_storage(p) - veg_nf%core%hrv_deadcroot_storage_to_litter(p) * dt

         ! transfer pools
         veg_ns%core%leaf_xfer(p)      = veg_ns%core%leaf_xfer(p)      - veg_nf%core%hrv_leaf_xfer_to_litter(p)      *dt
         veg_ns%core%froot_xfer(p)     = veg_ns%core%froot_xfer(p)     - veg_nf%core%hrv_froot_xfer_to_litter(p)     *dt
         veg_ns%core%livestem_xfer(p)  = veg_ns%core%livestem_xfer(p)  - veg_nf%core%hrv_livestem_xfer_to_litter(p)  *dt
         veg_ns%core%deadstem_xfer(p)  = veg_ns%core%deadstem_xfer(p)  - veg_nf%core%hrv_deadstem_xfer_to_litter(p)  *dt
         veg_ns%core%livecroot_xfer(p) = veg_ns%core%livecroot_xfer(p) - veg_nf%core%hrv_livecroot_xfer_to_litter(p) *dt
         veg_ns%core%deadcroot_xfer(p) = veg_ns%core%deadcroot_xfer(p) - veg_nf%core%hrv_deadcroot_xfer_to_litter(p) *dt

      end do

    end associate

  end subroutine NitrogenStateUpdate2h

end module NitrogenStateUpdate2Mod
