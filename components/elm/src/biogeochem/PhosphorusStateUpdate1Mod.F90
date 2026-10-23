module PhosphorusStateUpdate1Mod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for phosphorus state variable updates, non-mortality fluxes.
  ! X.YANG
  ! !USES:
  use shr_kind_mod           , only: r8 => shr_kind_r8
  use elm_varpar             , only : nlevdecomp, ndecomp_pools, ndecomp_cascade_transitions
  use elm_varpar             , only : crop_prog, i_met_lit, i_cel_lit, i_lig_lit, i_cwd
  use elm_varctl             , only : iulog
  use pftvarcon              , only : iscft
  use soilorder_varcon       , only : smax,ks_sorption
  use VegetationPropertiesType         , only : veg_vp
  use CNDecompCascadeConType , only : decomp_cascade_con
  use CNStateType            , only : cnstate_type
  use VegetationType              , only : veg_pp
  ! bgc interface & pflotran:
  use elm_varctl             , only : use_pflotran, pf_cmode
  use elm_varctl             , only : nu_com
  ! forest fertilization experiment
  use elm_time_manager       , only : get_curr_date
  use CNStateType            , only : fert_type , fert_continue, fert_dose, fert_start, fert_end
  use elm_varctl             , only : forest_fert_exp
  use elm_varctl             , only : NFIX_PTASE_plant
  use decompMod              , only : bounds_type
  use elm_varcon             , only : dzsoi_decomp
  use elm_varctl             , only : use_fates
  use GridcellDataType       , only : grc_ps, grc_pf
  use ColumnDataType         , only : col_ps, col_pf
  use VegetationDataType     , only : veg_ps, veg_pf

  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: PhosphorusStateUpdateDynPatch
  public :: PhosphorusStateUpdate1
  !-----------------------------------------------------------------------

contains
  subroutine PhosphorusStateUpdateDynPatch(bounds, num_soilc_with_inactive,&
     filter_soilc_with_inactive, dt)
    !
    ! !DESCRIPTION:
    ! Update phosphorus states based on fluxes from dyn_cnbal_patch
    !$acc routine seq
    ! !ARGUMENTS:
    type(bounds_type)          , intent(in)    :: bounds
    integer                    , intent(in)    :: num_soilc_with_inactive       ! number of columns in soil filter
    integer                    , intent(in)    :: filter_soilc_with_inactive(:) ! soil column filter that includes inactive points
    !
    real(r8)                   , intent(in)   :: dt
    ! !LOCAL VARIABLES:
    integer                                    :: c                             ! column index
    integer                                    :: fc                            ! column filter index
    integer                                    :: g                             ! gridcell index
    integer                                    :: j                             ! level index

    !-----------------------------------------------------------------------

      if (.not.use_fates) then

         do g = bounds%begg, bounds%endg
            grc_ps%core%seed(g) = grc_ps%core%seed(g) &
                 - grc_pf%core%dwt_seed_to_leaf(g)     * dt &
                 - grc_pf%core%dwt_seed_to_deadstem(g) * dt &
                 - grc_pf%dwt_seedp_to_ppool(g)    * dt
         end do

         do fc = 1, num_soilc_with_inactive

            c = filter_soilc_with_inactive(fc)
            col_ps%core%prod10(c) = col_ps%core%prod10(c) + col_pf%core%dwt_prod10_gain(c)*dt
            col_ps%core%prod100(c) = col_ps%core%prod100(c) + col_pf%core%dwt_prod100_gain(c)*dt
            col_ps%core%prod1(c) = col_ps%core%prod1(c) + col_pf%core%dwt_crop_product_gain(c)*dt

            do j = 1,nlevdecomp

               col_ps%core%decomp_pools_vr(c,j,i_met_lit) = col_ps%core%decomp_pools_vr(c,j,i_met_lit) + &
                    col_pf%core%dwt_froot_to_litr_met(c,j) * dt
               col_ps%core%decomp_pools_vr(c,j,i_cel_lit) = col_ps%core%decomp_pools_vr(c,j,i_cel_lit) + &
                    col_pf%core%dwt_froot_to_litr_cel(c,j) * dt
               col_ps%core%decomp_pools_vr(c,j,i_lig_lit) = col_ps%core%decomp_pools_vr(c,j,i_lig_lit) + &
                    col_pf%core%dwt_froot_to_litr_lig(c,j) * dt
               col_ps%core%decomp_pools_vr(c,j,i_cwd) = col_ps%core%decomp_pools_vr(c,j,i_cwd) + &
                    ( col_pf%core%dwt_livecroot_to_cwd(c,j) + col_pf%core%dwt_deadcroot_to_cwd(c,j) ) * dt

            end do
         end do
      end if


  end subroutine PhosphorusStateUpdateDynPatch

  !-----------------------------------------------------------------------
  subroutine PhosphorusStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
       cnstate_vars, dt)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update all the prognostic phosphorus state
    ! variables (except for gap-phase mortality and fire fluxes)
    !
    ! !ARGUMENTS:
      !$acc routine seq
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                  , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnstate_type)       , intent(in)    :: cnstate_vars
    !
    real(r8)                 , intent(in)    :: dt !radiation time step
    ! !LOCAL VARIABLES:
    integer :: c,p,j,l,k ! indices
    integer :: fp,fc     ! lake filter indices

    integer:: kyr                     ! current year
    integer:: kmo                     ! month of year  (1, ..., 12)
    integer:: kda                     ! day of month   (1, ..., 31)
    integer:: mcsec                   ! seconds of day (0, ..., seconds/day)
    !-----------------------------------------------------------------------

    associate(                                                                                           &
         ivt                   => veg_pp%itype                                , & ! Input:  [integer  (:)     ]  pft vegetation type

         woody                 => veg_vp%woody                         , & ! Input:  [real(r8) (:)     ]  woody lifeform flag (0 = non-woody, 1 = tree, 2 = shrub)

         cascade_donor_pool    => decomp_cascade_con%cascade_donor_pool    , & ! Input:  [integer  (:)     ]  which pool is C taken from for a given decomposition step
         cascade_receiver_pool => decomp_cascade_con%cascade_receiver_pool , & ! Input:  [integer  (:)     ]  which pool is C added to for a given decomposition step

         !!! N deposition profile, will weathering profile be needed?  -X.YANG
         ndep_prof             => cnstate_vars%ndep_prof_col               & ! Input:  [real(r8) (:,:)   ]  profile over which N deposition is distributed through column (1/m)
         )


      ! column-level fluxes

      !------------------------------------------------------------------
      ! if coupled with pflotran, the following updates are NOT needed
      ! if (.not.(use_pflotran .and. pf_cmode)) then
      !------------------------------------------------------------------

         if(.not.use_fates)then
            do j = 1, nlevdecomp
               do fc = 1,num_soilc
                  c = filter_soilc(fc)

                  ! plant to litter fluxes
                  ! phenology and dynamic landcover fluxes
                  col_pf%core%decomp_pools_sourcesink(c,j,i_met_lit) = &
                       col_pf%core%phenology_to_litr_met(c,j) * dt

                  col_pf%core%decomp_pools_sourcesink(c,j,i_cel_lit) = &
                       col_pf%core%phenology_to_litr_cel(c,j) * dt

                  col_pf%core%decomp_pools_sourcesink(c,j,i_lig_lit) = &
                       col_pf%core%phenology_to_litr_lig(c,j) * dt

               end do
            end do
         end if

        ! P fertilization for crops
      if ( crop_prog )then
         do j = 1, nlevdecomp

            ! column loop
            do fc = 1,num_soilc
               c = filter_soilc(fc)
                  ! P fertilization
                  col_ps%solutionp_vr(c,j) = col_ps%solutionp_vr(c,j) + col_pf%fert_p_to_sminp(c)*dt * ndep_prof(c,j)
            end do
         end do
      end if

      ! decomposition fluxes
      do k = 1, ndecomp_cascade_transitions
         do j = 1, nlevdecomp
            ! column loop
            do fc = 1,num_soilc
               c = filter_soilc(fc)

               col_pf%core%decomp_pools_sourcesink(c,j,cascade_donor_pool(k)) = &
                    col_pf%core%decomp_pools_sourcesink(c,j,cascade_donor_pool(k)) - &
                    col_pf%core%decomp_cascade_transfer_vr(c,j,k) * dt
            end do
         end do
      end do
      do k = 1, ndecomp_cascade_transitions
         if ( cascade_receiver_pool(k) /= 0 ) then  ! skip terminal transitions
            do j = 1, nlevdecomp
               ! column loop
               do fc = 1,num_soilc
                  c = filter_soilc(fc)

                  col_pf%core%decomp_pools_sourcesink(c,j,cascade_receiver_pool(k)) = &
                       col_pf%core%decomp_pools_sourcesink(c,j,cascade_receiver_pool(k)) + &
                       (col_pf%core%decomp_cascade_transfer_vr(c,j,k) + col_pf%mineral%decomp_cascade_smin_flux_vr(c,j,k)) * dt
               end do
            end do
         else  ! terminal transitions
            do j = 1, nlevdecomp
               ! column loop
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  col_pf%core%decomp_pools_sourcesink(c,j,cascade_donor_pool(k)) = &
                       col_pf%core%decomp_pools_sourcesink(c,j,cascade_donor_pool(k)) - &
                       col_pf%mineral%decomp_cascade_smin_flux_vr(c,j,k) * dt
               end do
            end do
         end if
      end do
      !------------------------------------------------------------------

      ! forest fertilization
      call get_curr_date(kyr, kmo, kda, mcsec)
      if (forest_fert_exp) then
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            if ( ((fert_continue(c) == 1 .and. kyr > fert_start(c) .and. kyr <= fert_end(c)) .or.  kyr == fert_start(c)) &
               .and. fert_type(c) == 2 &
               .and. kda == 1  .and. mcsec == 1800) then ! fertilization assumed to occur at the begnining of each month
               do j = 1, nlevdecomp
                  col_ps%solutionp_vr(c,j) = col_ps%solutionp_vr(c,j) + fert_dose(c,kmo)*ndep_prof(c,j)
               end do
            end if
         end do
      end if

      if(.not.use_fates)then

          ! patch loop (veg)

          do fp = 1,num_soilp
              p = filter_soilp(fp)

              ! phenology: transfer growth fluxes
              veg_ps%core%leaf(p)       = veg_ps%core%leaf(p)       + veg_pf%core%leaf_xfer_to_leaf(p)*dt
              veg_ps%core%leaf_xfer(p)  = veg_ps%core%leaf_xfer(p)  - veg_pf%core%leaf_xfer_to_leaf(p)*dt
              veg_ps%core%froot(p)      = veg_ps%core%froot(p)      + veg_pf%core%froot_xfer_to_froot(p)*dt
              veg_ps%core%froot_xfer(p) = veg_ps%core%froot_xfer(p) - veg_pf%core%froot_xfer_to_froot(p)*dt

              if (woody(ivt(p)) >= 1.0_r8) then
                  veg_ps%core%livestem(p)       = veg_ps%core%livestem(p)       + veg_pf%core%livestem_xfer_to_livestem(p)*dt
                  veg_ps%core%livestem_xfer(p)  = veg_ps%core%livestem_xfer(p)  - veg_pf%core%livestem_xfer_to_livestem(p)*dt
                  veg_ps%core%deadstem(p)       = veg_ps%core%deadstem(p)       + veg_pf%core%deadstem_xfer_to_deadstem(p)*dt
                  veg_ps%core%deadstem_xfer(p)  = veg_ps%core%deadstem_xfer(p)  - veg_pf%core%deadstem_xfer_to_deadstem(p)*dt
                  veg_ps%core%livecroot(p)      = veg_ps%core%livecroot(p)      + veg_pf%core%livecroot_xfer_to_livecroot(p)*dt
                  veg_ps%core%livecroot_xfer(p) = veg_ps%core%livecroot_xfer(p) - veg_pf%core%livecroot_xfer_to_livecroot(p)*dt
                  veg_ps%core%deadcroot(p)      = veg_ps%core%deadcroot(p)      + veg_pf%core%deadcroot_xfer_to_deadcroot(p)*dt
                  veg_ps%core%deadcroot_xfer(p) = veg_ps%core%deadcroot_xfer(p) - veg_pf%core%deadcroot_xfer_to_deadcroot(p)*dt
              end if

              if (iscft(ivt(p))) then ! skip 2 generic crops
                  ! lines here for consistency; the transfer terms are zero
                  veg_ps%core%livestem(p)       = veg_ps%core%livestem(p)      + veg_pf%core%livestem_xfer_to_livestem(p)*dt
                  veg_ps%core%livestem_xfer(p)  = veg_ps%core%livestem_xfer(p) - veg_pf%core%livestem_xfer_to_livestem(p)*dt
                  veg_ps%core%grain(p)          = veg_ps%core%grain(p)         + veg_pf%core%grain_xfer_to_grain(p)*dt
                  veg_ps%core%grain_xfer(p)     = veg_ps%core%grain_xfer(p)    - veg_pf%core%grain_xfer_to_grain(p)*dt
              end if

              ! phenology: litterfall and retranslocation fluxes
              veg_ps%core%leaf(p)    = veg_ps%core%leaf(p)    - veg_pf%core%leaf_to_litter(p)*dt
              veg_ps%core%froot(p)   = veg_ps%core%froot(p)   - veg_pf%core%froot_to_litter(p)*dt
              veg_ps%core%leaf(p)    = veg_ps%core%leaf(p)    - veg_pf%retrans%leaf_to_retrans(p)*dt
              veg_ps%retrans%pool(p) = veg_ps%retrans%pool(p) + veg_pf%retrans%leaf_to_retrans(p)*dt

              ! live wood turnover and retranslocation fluxes
              if (woody(ivt(p)) >= 1.0_r8) then
                  veg_ps%core%livestem(p)  = veg_ps%core%livestem(p)  - veg_pf%core%livestem_to_deadstem(p)*dt
                  veg_ps%core%deadstem(p)  = veg_ps%core%deadstem(p)  + veg_pf%core%livestem_to_deadstem(p)*dt
                  veg_ps%core%livestem(p)  = veg_ps%core%livestem(p)  - veg_pf%retrans%livestem_to_retrans(p)*dt
                  veg_ps%retrans%pool(p)   = veg_ps%retrans%pool(p)   + veg_pf%retrans%livestem_to_retrans(p)*dt
                  veg_ps%core%livecroot(p) = veg_ps%core%livecroot(p) - veg_pf%core%livecroot_to_deadcroot(p)*dt
                  veg_ps%core%deadcroot(p) = veg_ps%core%deadcroot(p) + veg_pf%core%livecroot_to_deadcroot(p)*dt
                  veg_ps%core%livecroot(p) = veg_ps%core%livecroot(p) - veg_pf%retrans%livecroot_to_retrans(p)*dt
                  veg_ps%retrans%pool(p)   = veg_ps%retrans%pool(p)   + veg_pf%retrans%livecroot_to_retrans(p)*dt
              end if
              if (iscft(ivt(p))) then ! Beth adds retrans from froot
                  veg_ps%core%froot(p)     = veg_ps%core%froot(p)     - veg_pf%retrans%froot_to_retrans(p)*dt
                  veg_ps%retrans%pool(p)   = veg_ps%retrans%pool(p)   + veg_pf%retrans%froot_to_retrans(p)*dt
                  veg_ps%core%livestem(p)  = veg_ps%core%livestem(p)  - veg_pf%core%livestem_to_litter(p)*dt
                  veg_ps%core%livestem(p)  = veg_ps%core%livestem(p)  - veg_pf%retrans%livestem_to_retrans(p)*dt
                  veg_ps%retrans%pool(p)   = veg_ps%retrans%pool(p)   + veg_pf%retrans%livestem_to_retrans(p)*dt
                  veg_ps%core%grain(p)     = veg_ps%core%grain(p)     - veg_pf%core%grain_to_food(p)*dt

                  veg_ps%core%cropseed_deficit(p) = veg_ps%core%cropseed_deficit(p) &
                        - veg_pf%core%crop_seed_to_leaf(p) * dt
              end if

              ! uptake from soil mineral N pool
              veg_ps%core%pool(p) = &
                    veg_ps%core%pool(p) + veg_pf%uptake%smin_to_pool(p)*dt
              if (nu_com .ne. 'RD') veg_ps%core%pool(p) = veg_ps%core%pool(p) + veg_pf%uptake%supplement_to_plant(p)*dt
              if (NFIX_PTASE_plant) veg_ps%core%pool(p) = veg_ps%core%pool(p) + veg_pf%biochem_pmin_to_plant(p)*dt

              ! deployment from retranslocation pool
              veg_ps%core%pool(p)    = veg_ps%core%pool(p)    + veg_pf%retrans%retrans_to_pool(p)*dt
              veg_ps%retrans%pool(p) = veg_ps%retrans%pool(p) - veg_pf%retrans%retrans_to_pool(p)*dt

              ! allocation fluxes
              veg_ps%core%pool(p)           = veg_ps%core%pool(p)          - veg_pf%core%pool_to_leaf(p)*dt
              veg_ps%core%leaf(p)           = veg_ps%core%leaf(p)          + veg_pf%core%pool_to_leaf(p)*dt
              veg_ps%core%pool(p)           = veg_ps%core%pool(p)          - veg_pf%core%pool_to_leaf_storage(p)*dt
              veg_ps%core%leaf_storage(p)   = veg_ps%core%leaf_storage(p)  + veg_pf%core%pool_to_leaf_storage(p)*dt
              veg_ps%core%pool(p)           = veg_ps%core%pool(p)          - veg_pf%core%pool_to_froot(p)*dt
              veg_ps%core%froot(p)          = veg_ps%core%froot(p)         + veg_pf%core%pool_to_froot(p)*dt
              veg_ps%core%pool(p)           = veg_ps%core%pool(p)          - veg_pf%core%pool_to_froot_storage(p)*dt
              veg_ps%core%froot_storage(p)  = veg_ps%core%froot_storage(p) + veg_pf%core%pool_to_froot_storage(p)*dt

              if (woody(ivt(p)) >= 1.0_r8) then
                  veg_ps%core%pool(p)              = veg_ps%core%pool(p)              - veg_pf%core%pool_to_livestem(p)*dt
                  veg_ps%core%livestem(p)          = veg_ps%core%livestem(p)          + veg_pf%core%pool_to_livestem(p)*dt
                  veg_ps%core%pool(p)              = veg_ps%core%pool(p)              - veg_pf%core%pool_to_livestem_storage(p)*dt
                  veg_ps%core%livestem_storage(p)  = veg_ps%core%livestem_storage(p)  + veg_pf%core%pool_to_livestem_storage(p)*dt
                  veg_ps%core%pool(p)              = veg_ps%core%pool(p)              - veg_pf%core%pool_to_deadstem(p)*dt
                  veg_ps%core%deadstem(p)          = veg_ps%core%deadstem(p)          + veg_pf%core%pool_to_deadstem(p)*dt
                  veg_ps%core%pool(p)              = veg_ps%core%pool(p)              - veg_pf%core%pool_to_deadstem_storage(p)*dt
                  veg_ps%core%deadstem_storage(p)  = veg_ps%core%deadstem_storage(p)  + veg_pf%core%pool_to_deadstem_storage(p)*dt
                  veg_ps%core%pool(p)              = veg_ps%core%pool(p)              - veg_pf%core%pool_to_livecroot(p)*dt
                  veg_ps%core%livecroot(p)         = veg_ps%core%livecroot(p)         + veg_pf%core%pool_to_livecroot(p)*dt
                  veg_ps%core%pool(p)              = veg_ps%core%pool(p)              - veg_pf%core%pool_to_livecroot_storage(p)*dt
                  veg_ps%core%livecroot_storage(p) = veg_ps%core%livecroot_storage(p) + veg_pf%core%pool_to_livecroot_storage(p)*dt
                  veg_ps%core%pool(p)              = veg_ps%core%pool(p)              - veg_pf%core%pool_to_deadcroot(p)*dt
                  veg_ps%core%deadcroot(p)         = veg_ps%core%deadcroot(p)         + veg_pf%core%pool_to_deadcroot(p)*dt
                  veg_ps%core%pool(p)              = veg_ps%core%pool(p)              - veg_pf%core%pool_to_deadcroot_storage(p)*dt
                  veg_ps%core%deadcroot_storage(p) = veg_ps%core%deadcroot_storage(p) + veg_pf%core%pool_to_deadcroot_storage(p)*dt
              end if

              if (iscft(ivt(p))) then ! skip 2 generic crops
                  veg_ps%core%pool(p)              = veg_ps%core%pool(p)              - veg_pf%core%pool_to_livestem(p)*dt
                  veg_ps%core%livestem(p)          = veg_ps%core%livestem(p)          + veg_pf%core%pool_to_livestem(p)*dt
                  veg_ps%core%pool(p)              = veg_ps%core%pool(p)              - veg_pf%core%pool_to_livestem_storage(p)*dt
                  veg_ps%core%livestem_storage(p)  = veg_ps%core%livestem_storage(p)  + veg_pf%core%pool_to_livestem_storage(p)*dt
                  veg_ps%core%pool(p)              = veg_ps%core%pool(p)              - veg_pf%core%pool_to_grain(p)*dt
                  veg_ps%core%grain(p)             = veg_ps%core%grain(p)             + veg_pf%core%pool_to_grain(p)*dt
                  veg_ps%core%pool(p)              = veg_ps%core%pool(p)              - veg_pf%core%pool_to_grain_storage(p)*dt
                  veg_ps%core%grain_storage(p)     = veg_ps%core%grain_storage(p)     + veg_pf%core%pool_to_grain_storage(p)*dt
              end if

              ! move storage pools into transfer pools
              veg_ps%core%leaf_storage(p)  = veg_ps%core%leaf_storage(p)  - veg_pf%core%leaf_storage_to_xfer(p)*dt
              veg_ps%core%leaf_xfer(p)     = veg_ps%core%leaf_xfer(p)     + veg_pf%core%leaf_storage_to_xfer(p)*dt
              veg_ps%core%froot_storage(p) = veg_ps%core%froot_storage(p) - veg_pf%core%froot_storage_to_xfer(p)*dt
              veg_ps%core%froot_xfer(p)    = veg_ps%core%froot_xfer(p)    + veg_pf%core%froot_storage_to_xfer(p)*dt

              if (woody(ivt(p)) >= 1.0_r8) then
                  veg_ps%core%livestem_storage(p)  = veg_ps%core%livestem_storage(p)  - veg_pf%core%livestem_storage_to_xfer(p)*dt
                  veg_ps%core%livestem_xfer(p)     = veg_ps%core%livestem_xfer(p)     + veg_pf%core%livestem_storage_to_xfer(p)*dt
                  veg_ps%core%deadstem_storage(p)  = veg_ps%core%deadstem_storage(p)  - veg_pf%core%deadstem_storage_to_xfer(p)*dt
                  veg_ps%core%deadstem_xfer(p)     = veg_ps%core%deadstem_xfer(p)     + veg_pf%core%deadstem_storage_to_xfer(p)*dt
                  veg_ps%core%livecroot_storage(p) = veg_ps%core%livecroot_storage(p) - veg_pf%core%livecroot_storage_to_xfer(p)*dt
                  veg_ps%core%livecroot_xfer(p)    = veg_ps%core%livecroot_xfer(p)    + veg_pf%core%livecroot_storage_to_xfer(p)*dt
                  veg_ps%core%deadcroot_storage(p) = veg_ps%core%deadcroot_storage(p) - veg_pf%core%deadcroot_storage_to_xfer(p)*dt
                  veg_ps%core%deadcroot_xfer(p)    = veg_ps%core%deadcroot_xfer(p)    + veg_pf%core%deadcroot_storage_to_xfer(p)*dt
              end if

              if (iscft(ivt(p))) then ! skip 2 generic crops
                  ! lines here for consistency; the transfer terms are zero
                  veg_ps%core%livestem_storage(p)  = veg_ps%core%livestem_storage(p) - veg_pf%core%livestem_storage_to_xfer(p)*dt
                  veg_ps%core%livestem_xfer(p)     = veg_ps%core%livestem_xfer(p)    + veg_pf%core%livestem_storage_to_xfer(p)*dt
                  veg_ps%core%grain_storage(p)     = veg_ps%core%grain_storage(p)    - veg_pf%core%grain_storage_to_xfer(p)*dt
                  veg_ps%core%grain_xfer(p)        = veg_ps%core%grain_xfer(p)       + veg_pf%core%grain_storage_to_xfer(p)*dt
              end if

          end do ! num_soilp
      end if  !if(.not.use_fates)

    end associate

  end subroutine PhosphorusStateUpdate1

end module PhosphorusStateUpdate1Mod
