module NitrogenStateUpdate3Mod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for nitrogen state variable update, mortality fluxes.
  ! Also, sminn leaching flux.
  !
  ! !USES:
  use shr_kind_mod        , only: r8 => shr_kind_r8
  use elm_varpar          , only: nlevdecomp, ndecomp_pools
  use elm_varctl          , only : iulog
  use elm_varpar          , only : i_cwd, i_met_lit, i_cel_lit, i_lig_lit
  use elm_varctl          , only : use_erosion, ero_ccycle
  use CNDecompCascadeConType , only : decomp_cascade_con
  use ColumnDataType      , only : col_ns, col_nf
  use VegetationDataType  , only : veg_ns, veg_nf
  ! bgc interface & pflotran:
  use elm_varctl          , only : use_pflotran, pf_cmode, use_fates
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: NitrogenStateUpdate3
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine NitrogenStateUpdate3(num_soilc, filter_soilc, num_soilp, filter_soilp, dt)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update all the prognostic nitrogen state
    ! variables affected by gap-phase mortality fluxes. Also the Sminn leaching flux.
    ! Also include the erosion flux.
    ! NOTE - associate statements have been removed where there are
    ! no science equations. This increases readability and maintainability.
    !
    ! !ARGUMENTS:
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                  , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:) ! filter for soil patches
    real(r8)                  , intent(in)   :: dt         ! radiation time step (seconds)

    !
    ! !LOCAL VARIABLES:
    integer :: c,p,j,l,k        ! indices
    integer :: fp,fc      ! lake filter indices
    !-----------------------------------------------------------------------
      do j = 1, nlevdecomp
         ! column loop
         do fc = 1,num_soilc
            c = filter_soilc(fc)
               
            ! mineral N loss due to leaching and runoff
            col_ns%smin_no3_vr(c,j) = max( col_ns%smin_no3_vr(c,j) - &
                 ( col_nf%smin_no3_leached_vr(c,j) + col_nf%smin_no3_runoff_vr(c,j) ) * dt, 0._r8)
               
            col_ns%mineral%smin_vr(c,j) = col_ns%smin_no3_vr(c,j) + col_ns%smin_nh4_vr(c,j)
            if (use_pflotran .and. pf_cmode) then 
               col_ns%mineral%smin_vr(c,j) = col_ns%mineral%smin_vr(c,j) + col_ns%smin_nh4sorb_vr(c,j)
            end if

            if (.not.(use_pflotran .and. pf_cmode)) then
                ! column level nitrogen fluxes from fire
                ! pft-level wood to column-level CWD (uncombusted wood)
                col_ns%core%decomp_pools_vr(c,j,i_cwd) = col_ns%core%decomp_pools_vr(c,j,i_cwd) &
                     + col_nf%core%fire_mortality_to_cwd(c,j) * dt

                ! pft-level wood to column-level litter (uncombusted wood)
                col_ns%core%decomp_pools_vr(c,j,i_met_lit) = col_ns%core%decomp_pools_vr(c,j,i_met_lit) &
                     + col_nf%core%m_to_litr_met_fire(c,j)* dt
                col_ns%core%decomp_pools_vr(c,j,i_cel_lit) = col_ns%core%decomp_pools_vr(c,j,i_cel_lit) &
                     + col_nf%core%m_to_litr_cel_fire(c,j)* dt
                col_ns%core%decomp_pools_vr(c,j,i_lig_lit) = col_ns%core%decomp_pools_vr(c,j,i_lig_lit) &
                     + col_nf%core%m_to_litr_lig_fire(c,j)* dt
            end if !(.not.(use_pflotran .and. pf_cmode))
         end do ! end of column loop
      end do

      ! litter and CWD losses to fire
      do l = 1, ndecomp_pools
         do j = 1, nlevdecomp
            ! column loop
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               col_ns%core%decomp_pools_vr(c,j,l) = col_ns%core%decomp_pools_vr(c,j,l) - col_nf%core%m_decomp_pools_to_fire_vr(c,j,l) * dt
            end do
         end do
      end do


      ! SOM N losses due to erosion
      if ( ero_ccycle ) then
         do l = 1, ndecomp_pools
            if ( decomp_cascade_con%is_soil(l) ) then
               do j = 1, nlevdecomp
                  do fc = 1, num_soilc
                     c = filter_soilc(fc)
                     col_ns%core%decomp_pools_vr(c,j,l) = col_ns%core%decomp_pools_vr(c,j,l) - col_nf%core%decomp_pools_yield_vr(c,j,l) * dt
                  end do
               end do
            end if
         end do
      end if

      ! patch-level nitrogen fluxes
      if(.not.use_fates) then

         do fp = 1,num_soilp
            p = filter_soilp(fp)

            !from fire displayed pools
            veg_ns%core%leaf(p)              =  veg_ns%core%leaf(p)      - veg_nf%core%m_leaf_to_fire(p)      * dt
            veg_ns%core%froot(p)             =  veg_ns%core%froot(p)     - veg_nf%core%m_froot_to_fire(p)     * dt
            veg_ns%core%livestem(p)          =  veg_ns%core%livestem(p)  - veg_nf%core%m_livestem_to_fire(p)  * dt
            veg_ns%core%deadstem(p)          =  veg_ns%core%deadstem(p)  - veg_nf%core%m_deadstem_to_fire(p)  * dt
            veg_ns%core%livecroot(p)         =  veg_ns%core%livecroot(p) - veg_nf%core%m_livecroot_to_fire(p) * dt
            veg_ns%core%deadcroot(p)         =  veg_ns%core%deadcroot(p) - veg_nf%core%m_deadcroot_to_fire(p) * dt

            veg_ns%core%leaf(p)              =  veg_ns%core%leaf(p)      - veg_nf%core%m_leaf_to_litter_fire(p)           * dt
            veg_ns%core%froot(p)             =  veg_ns%core%froot(p)     - veg_nf%core%m_froot_to_litter_fire(p)          * dt
            veg_ns%core%livestem(p)          =  veg_ns%core%livestem(p)  - veg_nf%core%m_livestem_to_litter_fire(p)       * dt
            veg_ns%core%deadstem(p)          =  veg_ns%core%deadstem(p)  - veg_nf%core%m_deadstem_to_litter_fire(p)       * dt
            veg_ns%core%livecroot(p)         =  veg_ns%core%livecroot(p) - veg_nf%core%m_livecroot_to_litter_fire(p)      * dt
            veg_ns%core%deadcroot(p)         =  veg_ns%core%deadcroot(p) - veg_nf%core%m_deadcroot_to_litter_fire(p)      * dt

            ! storage pools
            veg_ns%core%leaf_storage(p)      =  veg_ns%core%leaf_storage(p)      - veg_nf%core%m_leaf_storage_to_fire(p)      * dt
            veg_ns%core%froot_storage(p)     =  veg_ns%core%froot_storage(p)     - veg_nf%core%m_froot_storage_to_fire(p)     * dt
            veg_ns%core%livestem_storage(p)  =  veg_ns%core%livestem_storage(p)  - veg_nf%core%m_livestem_storage_to_fire(p)  * dt
            veg_ns%core%deadstem_storage(p)  =  veg_ns%core%deadstem_storage(p)  - veg_nf%core%m_deadstem_storage_to_fire(p)  * dt
            veg_ns%core%livecroot_storage(p) =  veg_ns%core%livecroot_storage(p) - veg_nf%core%m_livecroot_storage_to_fire(p) * dt
            veg_ns%core%deadcroot_storage(p) =  veg_ns%core%deadcroot_storage(p) - veg_nf%core%m_deadcroot_storage_to_fire(p) * dt

            veg_ns%core%leaf_storage(p)      =  veg_ns%core%leaf_storage(p)      - veg_nf%core%m_leaf_storage_to_litter_fire(p)      * dt
            veg_ns%core%froot_storage(p)     =  veg_ns%core%froot_storage(p)     - veg_nf%core%m_froot_storage_to_litter_fire(p)     * dt
            veg_ns%core%livestem_storage(p)  =  veg_ns%core%livestem_storage(p)  - veg_nf%core%m_livestem_storage_to_litter_fire(p)  * dt
            veg_ns%core%deadstem_storage(p)  =  veg_ns%core%deadstem_storage(p)  - veg_nf%core%m_deadstem_storage_to_litter_fire(p)  * dt
            veg_ns%core%livecroot_storage(p) =  veg_ns%core%livecroot_storage(p) - veg_nf%core%m_livecroot_storage_to_litter_fire(p) * dt
            veg_ns%core%deadcroot_storage(p) =  veg_ns%core%deadcroot_storage(p) - veg_nf%core%m_deadcroot_storage_to_litter_fire(p) * dt


            ! transfer pools
            veg_ns%core%leaf_xfer(p)         =  veg_ns%core%leaf_xfer(p)      - veg_nf%core%m_leaf_xfer_to_fire(p)      * dt
            veg_ns%core%froot_xfer(p)        =  veg_ns%core%froot_xfer(p)     - veg_nf%core%m_froot_xfer_to_fire(p)     * dt
            veg_ns%core%livestem_xfer(p)     =  veg_ns%core%livestem_xfer(p)  - veg_nf%core%m_livestem_xfer_to_fire(p)  * dt
            veg_ns%core%deadstem_xfer(p)     =  veg_ns%core%deadstem_xfer(p)  - veg_nf%core%m_deadstem_xfer_to_fire(p)  * dt
            veg_ns%core%livecroot_xfer(p)    =  veg_ns%core%livecroot_xfer(p) - veg_nf%core%m_livecroot_xfer_to_fire(p) * dt
            veg_ns%core%deadcroot_xfer(p)    =  veg_ns%core%deadcroot_xfer(p) - veg_nf%core%m_deadcroot_xfer_to_fire(p) * dt

            veg_ns%core%leaf_xfer(p)         =  veg_ns%core%leaf_xfer(p)      - veg_nf%core%m_leaf_xfer_to_litter_fire(p)      * dt
            veg_ns%core%froot_xfer(p)        =  veg_ns%core%froot_xfer(p)     - veg_nf%core%m_froot_xfer_to_litter_fire(p)     * dt
            veg_ns%core%livestem_xfer(p)     =  veg_ns%core%livestem_xfer(p)  - veg_nf%core%m_livestem_xfer_to_litter_fire(p)  * dt
            veg_ns%core%deadstem_xfer(p)     =  veg_ns%core%deadstem_xfer(p)  - veg_nf%core%m_deadstem_xfer_to_litter_fire(p)  * dt
            veg_ns%core%livecroot_xfer(p)    =  veg_ns%core%livecroot_xfer(p) - veg_nf%core%m_livecroot_xfer_to_litter_fire(p) * dt
            veg_ns%core%deadcroot_xfer(p)    =  veg_ns%core%deadcroot_xfer(p) - veg_nf%core%m_deadcroot_xfer_to_litter_fire(p) * dt

            ! retranslocated N pool
            veg_ns%retrans%pool(p)           =  veg_ns%retrans%pool(p) - veg_nf%retrans%m_retrans_to_fire(p)        * dt
            veg_ns%retrans%pool(p)           =  veg_ns%retrans%pool(p) - veg_nf%retrans%m_retrans_to_litter_fire(p) * dt
            veg_ns%core%pool(p)              =  veg_ns%core%pool(p)    - veg_nf%core%m_pool_to_fire(p)           * dt
            veg_ns%core%pool(p)              =  veg_ns%core%pool(p)    - veg_nf%core%m_pool_to_litter_fire(p)    * dt
         end do
      end if


  end subroutine NitrogenStateUpdate3

end module NitrogenStateUpdate3Mod
