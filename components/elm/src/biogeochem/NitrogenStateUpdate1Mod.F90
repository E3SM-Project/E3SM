module NitrogenStateUpdate1Mod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for nitrogen state variable updates, non-mortality fluxes.
  !
  ! !USES:
  use shr_kind_mod           , only: r8 => shr_kind_r8
  use elm_varpar             , only : nlevdecomp, ndecomp_pools, ndecomp_cascade_transitions
  use elm_varpar             , only : crop_prog, i_met_lit, i_cel_lit, i_lig_lit, i_cwd
  use elm_varctl             , only : iulog
  use elm_varcon             , only : nitrif_n2o_loss_frac
  use pftvarcon              , only : iscft
  use VegetationPropertiesType         , only : veg_vp
  use CNDecompCascadeConType , only : decomp_cascade_con
  use CNStateType            , only : cnstate_type
  use GridcellDataType       , only : grc_ns, grc_nf
  use ColumnDataType         , only : col_ns, col_nf
  use VegetationType         , only : veg_pp
  use VegetationDataType     , only : veg_ns, veg_nf
  ! bgc interface & pflotran:
  use elm_varctl             , only : use_pflotran, pf_cmode
  ! forest fertilization experiment
  use elm_time_manager       , only : get_curr_date
  use CNStateType            , only : fert_type , fert_continue, fert_dose, fert_start, fert_end
  use elm_varctl             , only : forest_fert_exp
  use elm_varctl             , only : nu_com
  use elm_varctl             , only : NFIX_PTASE_plant
  use decompMod              , only : bounds_type
  use elm_varcon             , only : dzsoi_decomp
  use elm_varctl             , only : use_fates
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: NitrogenStateUpdateDynPatch
  public :: NitrogenStateUpdate1
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine NitrogenStateUpdateDynPatch(bounds, num_soilc_with_inactive,&
     filter_soilc_with_inactive, dt)
    !
    ! !DESCRIPTION:
    ! Update nitrogen states based on fluxes from dyn_cnbal_patch
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc_with_inactive       ! number of columns in soil filter
    integer                  , intent(in)    :: filter_soilc_with_inactive(:) ! soil column filter that includes inactive points
    real(r8), intent(in)                      :: dt                            ! time step (seconds)

    !
    ! !LOCAL VARIABLES:
    integer                                  :: c                             ! column index
    integer                                  :: fc                            ! column filter index
    integer                                  :: g                             ! gridcell index
    integer                                  :: j                             ! level index

    !character(len=*)         , parameter     :: subname = 'NitrogenStateUpdateDynPatch'
    !-----------------------------------------------------------------------

      if (.not.use_fates) then

         do g = bounds%begg, bounds%endg
            grc_ns%core%seed(g) = grc_ns%core%seed(g) &
                 - grc_nf%core%dwt_seed_to_leaf(g)     * dt &
                 - grc_nf%core%dwt_seed_to_deadstem(g) * dt &
                 - grc_nf%dwt_seedn_to_npool(g)    * dt
         end do

         do fc = 1, num_soilc_with_inactive

            c = filter_soilc_with_inactive(fc)
            col_ns%core%prod10(c) = col_ns%core%prod10(c) + col_nf%core%dwt_prod10_gain(c)*dt
            col_ns%core%prod100(c) = col_ns%core%prod100(c) + col_nf%core%dwt_prod100_gain(c)*dt
            col_ns%core%prod1(c) = col_ns%core%prod1(c) + col_nf%core%dwt_crop_product_gain(c)*dt

            do j = 1,nlevdecomp

               col_ns%core%decomp_pools_vr(c,j,i_met_lit) = col_ns%core%decomp_pools_vr(c,j,i_met_lit) + &
                    col_nf%core%dwt_froot_to_litr_met(c,j) * dt
               col_ns%core%decomp_pools_vr(c,j,i_cel_lit) = col_ns%core%decomp_pools_vr(c,j,i_cel_lit) + &
                    col_nf%core%dwt_froot_to_litr_cel(c,j) * dt
               col_ns%core%decomp_pools_vr(c,j,i_lig_lit) = col_ns%core%decomp_pools_vr(c,j,i_lig_lit) + &
                    col_nf%core%dwt_froot_to_litr_lig(c,j) * dt
               col_ns%core%decomp_pools_vr(c,j,i_cwd) = col_ns%core%decomp_pools_vr(c,j,i_cwd) + &
                    ( col_nf%core%dwt_livecroot_to_cwd(c,j) + col_nf%core%dwt_deadcroot_to_cwd(c,j) ) * dt

            end do
         end do
      end if

  end subroutine NitrogenStateUpdateDynPatch

  !-----------------------------------------------------------------------
  subroutine NitrogenStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
       cnstate_vars, dt)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update all the prognostic nitrogen state
    ! variables (except for gap-phase mortality and fire fluxes)
    !
    ! !ARGUMENTS:
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                  , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnstate_type)       , intent(in)    :: cnstate_vars
    real(r8)                  , intent(in)    :: dt        ! radiation time step (seconds)

    !
    ! !LOCAL VARIABLES:
    integer :: c,p,j,l,k ! indices
    integer :: fp,fc     ! lake filter indices
    real(r8), parameter :: frootc_nfix_thc = 10._r8  !threshold fine root carbon for nitrogen fixation gC/m2

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

         ndep_prof             => cnstate_vars%ndep_prof_col               , & ! Input:  [real(r8) (:,:)   ]  profile over which N deposition is distributed through column (1/m)
         nfixation_prof        => cnstate_vars%nfixation_prof_col           & ! Input:  [real(r8) (:,:)   ]  profile over which N fixation is distributed through column (1/m)

         )

      ! column-level fluxes

      if (.not.(use_pflotran .and. pf_cmode)) then

         do j = 1, nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)

               ! N deposition and fixation (put all into NH4 pool)
               col_ns%smin_nh4_vr(c,j) = col_ns%smin_nh4_vr(c,j) + col_nf%mineral%dep_to_smin(c)*dt * ndep_prof(c,j)
               col_ns%smin_nh4_vr(c,j) = col_ns%smin_nh4_vr(c,j) + col_nf%nfix_to_sminn(c)*dt * nfixation_prof(c,j)

               ! plant to litter fluxes
               ! phenology and dynamic landcover fluxes
               if(.not.use_fates) then
                  col_nf%core%decomp_pools_sourcesink(c,j,i_met_lit) = &
                       col_nf%core%phenology_to_litr_met(c,j) * dt

                  col_nf%core%decomp_pools_sourcesink(c,j,i_cel_lit) = &
                       col_nf%core%phenology_to_litr_cel(c,j) * dt

                  col_nf%core%decomp_pools_sourcesink(c,j,i_lig_lit) = &
                       col_nf%core%phenology_to_litr_lig(c,j) * dt
               end if
            end do
         end do

         ! repeating N dep and fixation for crops
         if ( crop_prog )then
            do j = 1, nlevdecomp

               ! column loop
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  ! N deposition and fixation (put all into NH4 pool)
                  col_ns%smin_nh4_vr(c,j) = col_ns%smin_nh4_vr(c,j) + col_nf%fert_to_sminn(c)*dt * ndep_prof(c,j)
                  col_ns%smin_nh4_vr(c,j) = col_ns%smin_nh4_vr(c,j) + col_nf%soyfixn_to_sminn(c)*dt * nfixation_prof(c,j)
               end do
            end do
         end if

         ! decomposition fluxes
         do k = 1, ndecomp_cascade_transitions
            do j = 1, nlevdecomp
               ! column loop
               do fc = 1,num_soilc
                  c = filter_soilc(fc)

                  col_nf%core%decomp_pools_sourcesink(c,j,cascade_donor_pool(k)) = &
                       col_nf%core%decomp_pools_sourcesink(c,j,cascade_donor_pool(k)) - &
                       col_nf%core%decomp_cascade_transfer_vr(c,j,k) * dt
               end do
            end do
         end do

         do k = 1, ndecomp_cascade_transitions
            if ( cascade_receiver_pool(k) /= 0 ) then  ! skip terminal transitions
               do j = 1, nlevdecomp
                  ! column loop
                  do fc = 1,num_soilc
                     c = filter_soilc(fc)

                     col_nf%core%decomp_pools_sourcesink(c,j,cascade_receiver_pool(k)) = &
                          col_nf%core%decomp_pools_sourcesink(c,j,cascade_receiver_pool(k)) + &
                          (col_nf%core%decomp_cascade_transfer_vr(c,j,k) + col_nf%mineral%decomp_cascade_smin_flux_vr(c,j,k)) * dt
                  end do
               end do
            else  ! terminal transitions
               do j = 1, nlevdecomp
                  ! column loop
                  do fc = 1,num_soilc
                     c = filter_soilc(fc)
                     col_nf%core%decomp_pools_sourcesink(c,j,cascade_donor_pool(k)) = &
                          col_nf%core%decomp_pools_sourcesink(c,j,cascade_donor_pool(k)) - &
                          col_nf%mineral%decomp_cascade_smin_flux_vr(c,j,k) * dt
                  end do
               end do
            end if
         end do

         do j = 1, nlevdecomp
            ! column loop
            do fc = 1,num_soilc
               c = filter_soilc(fc)

               ! mineralization fluxes (divert a fraction of this stream to nitrification flux, add the rest to NH4 pool)
               col_ns%smin_nh4_vr(c,j) = col_ns%smin_nh4_vr(c,j) + col_nf%mineral%gross_min_vr(c,j)*dt

               ! immobilization fluxes
               col_ns%smin_nh4_vr(c,j) = col_ns%smin_nh4_vr(c,j) - col_nf%actual_immob_nh4_vr(c,j)*dt

               col_ns%smin_no3_vr(c,j) = col_ns%smin_no3_vr(c,j) - col_nf%actual_immob_no3_vr(c,j)*dt

               ! plant uptake fluxes
               col_ns%smin_nh4_vr(c,j) = col_ns%smin_nh4_vr(c,j) - col_nf%smin_nh4_to_plant_vr(c,j)*dt

               col_ns%smin_no3_vr(c,j) = col_ns%smin_no3_vr(c,j) - col_nf%smin_no3_to_plant_vr(c,j)*dt

               ! Account for nitrification fluxes
               col_ns%smin_nh4_vr(c,j) = col_ns%smin_nh4_vr(c,j) - col_nf%f_nit_vr(c,j) * dt

               col_ns%smin_no3_vr(c,j) = col_ns%smin_no3_vr(c,j) + col_nf%f_nit_vr(c,j) * dt * (1._r8 - nitrif_n2o_loss_frac)

               ! Account for denitrification fluxes
               col_ns%smin_no3_vr(c,j) = col_ns%smin_no3_vr(c,j) - col_nf%f_denit_vr(c,j) * dt

               ! flux that prevents N limitation (when Carbon_only is set; put all into NH4)
               col_ns%smin_nh4_vr(c,j) = col_ns%smin_nh4_vr(c,j) + col_nf%mineral%supplement_to_smin_vr(c,j)*dt

               ! update diagnostic total
               col_ns%mineral%smin_vr(c,j) = col_ns%smin_nh4_vr(c,j) + col_ns%smin_no3_vr(c,j)

            end do ! end of column loop
         end do

      endif  

      ! forest fertilization
      call get_curr_date(kyr, kmo, kda, mcsec)
      if (forest_fert_exp) then
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            if ( ((fert_continue(c) == 1 .and. kyr > fert_start(c) .and. kyr <= fert_end(c)) .or.  kyr == fert_start(c)) &
               .and. fert_type(c) == 1 &
               .and. kda == 1  .and. mcsec == 1800) then ! fertilization assumed to occur at the begnining of each month
               do j = 1, nlevdecomp
                  col_ns%smin_nh4_vr(c,j) = col_ns%smin_nh4_vr(c,j) + fert_dose(c,kmo)/2._r8*ndep_prof(c,j)
                  col_ns%smin_no3_vr(c,j) = col_ns%smin_no3_vr(c,j) + fert_dose(c,kmo)/2._r8*ndep_prof(c,j)
                  col_ns%mineral%smin_vr(c,j) = col_ns%smin_nh4_vr(c,j) + col_ns%smin_no3_vr(c,j)
               end do
            end if
         end do
      end if


      ! patch loop (veg)
      if(.not.use_fates) then
          do fp = 1,num_soilp
              p = filter_soilp(fp)

              ! phenology: transfer growth fluxes
              veg_ns%core%leaf(p)       = veg_ns%core%leaf(p)       + veg_nf%core%leaf_xfer_to_leaf(p)*dt
              veg_ns%core%leaf_xfer(p)  = veg_ns%core%leaf_xfer(p)  - veg_nf%core%leaf_xfer_to_leaf(p)*dt
              veg_ns%core%froot(p)      = veg_ns%core%froot(p)      + veg_nf%core%froot_xfer_to_froot(p)*dt
              veg_ns%core%froot_xfer(p) = veg_ns%core%froot_xfer(p) - veg_nf%core%froot_xfer_to_froot(p)*dt

              if (woody(ivt(p)) >= 1.0_r8) then
                  veg_ns%core%livestem(p)       = veg_ns%core%livestem(p)       + veg_nf%core%livestem_xfer_to_livestem(p)*dt
                  veg_ns%core%livestem_xfer(p)  = veg_ns%core%livestem_xfer(p)  - veg_nf%core%livestem_xfer_to_livestem(p)*dt
                  veg_ns%core%deadstem(p)       = veg_ns%core%deadstem(p)       + veg_nf%core%deadstem_xfer_to_deadstem(p)*dt
                  veg_ns%core%deadstem_xfer(p)  = veg_ns%core%deadstem_xfer(p)  - veg_nf%core%deadstem_xfer_to_deadstem(p)*dt
                  veg_ns%core%livecroot(p)      = veg_ns%core%livecroot(p)      + veg_nf%core%livecroot_xfer_to_livecroot(p)*dt
                  veg_ns%core%livecroot_xfer(p) = veg_ns%core%livecroot_xfer(p) - veg_nf%core%livecroot_xfer_to_livecroot(p)*dt
                  veg_ns%core%deadcroot(p)      = veg_ns%core%deadcroot(p)      + veg_nf%core%deadcroot_xfer_to_deadcroot(p)*dt
                  veg_ns%core%deadcroot_xfer(p) = veg_ns%core%deadcroot_xfer(p) - veg_nf%core%deadcroot_xfer_to_deadcroot(p)*dt
              end if

              if (iscft(ivt(p))) then ! skip 2 generic crops
                  ! lines here for consistency; the transfer terms are zero
                  veg_ns%core%livestem(p)       = veg_ns%core%livestem(p)      + veg_nf%core%livestem_xfer_to_livestem(p)*dt
                  veg_ns%core%livestem_xfer(p)  = veg_ns%core%livestem_xfer(p) - veg_nf%core%livestem_xfer_to_livestem(p)*dt
                  veg_ns%core%grain(p)          = veg_ns%core%grain(p)         + veg_nf%core%grain_xfer_to_grain(p)*dt
                  veg_ns%core%grain_xfer(p)     = veg_ns%core%grain_xfer(p)    - veg_nf%core%grain_xfer_to_grain(p)*dt
              end if

              ! phenology: litterfall and retranslocation fluxes
              veg_ns%core%leaf(p)    = veg_ns%core%leaf(p)    - veg_nf%core%leaf_to_litter(p)*dt
              veg_ns%core%froot(p)   = veg_ns%core%froot(p)   - veg_nf%core%froot_to_litter(p)*dt
              veg_ns%core%leaf(p)    = veg_ns%core%leaf(p)    - veg_nf%retrans%leaf_to_retrans(p)*dt
              veg_ns%retrans%pool(p) = veg_ns%retrans%pool(p) + veg_nf%retrans%leaf_to_retrans(p)*dt

              ! live wood turnover and retranslocation fluxes
              if (woody(ivt(p)) >= 1.0_r8) then
                  veg_ns%core%livestem(p)  = veg_ns%core%livestem(p)  - veg_nf%core%livestem_to_deadstem(p)*dt
                  veg_ns%core%deadstem(p)  = veg_ns%core%deadstem(p)  + veg_nf%core%livestem_to_deadstem(p)*dt
                  veg_ns%core%livestem(p)  = veg_ns%core%livestem(p)  - veg_nf%retrans%livestem_to_retrans(p)*dt
                  veg_ns%retrans%pool(p)   = veg_ns%retrans%pool(p)   + veg_nf%retrans%livestem_to_retrans(p)*dt
                  veg_ns%core%livecroot(p) = veg_ns%core%livecroot(p) - veg_nf%core%livecroot_to_deadcroot(p)*dt
                  veg_ns%core%deadcroot(p) = veg_ns%core%deadcroot(p) + veg_nf%core%livecroot_to_deadcroot(p)*dt
                  veg_ns%core%livecroot(p) = veg_ns%core%livecroot(p) - veg_nf%retrans%livecroot_to_retrans(p)*dt
                  veg_ns%retrans%pool(p)   = veg_ns%retrans%pool(p)   + veg_nf%retrans%livecroot_to_retrans(p)*dt
              end if
              if (iscft(ivt(p))) then ! Beth adds retrans from froot
                  veg_ns%core%froot(p)     = veg_ns%core%froot(p)     - veg_nf%retrans%froot_to_retrans(p)*dt
                  veg_ns%retrans%pool(p)   = veg_ns%retrans%pool(p)   + veg_nf%retrans%froot_to_retrans(p)*dt
                  veg_ns%core%livestem(p)  = veg_ns%core%livestem(p)  - veg_nf%core%livestem_to_litter(p)*dt
                  veg_ns%core%livestem(p)  = veg_ns%core%livestem(p)  - veg_nf%retrans%livestem_to_retrans(p)*dt
                  veg_ns%retrans%pool(p)   = veg_ns%retrans%pool(p)   + veg_nf%retrans%livestem_to_retrans(p)*dt
                  veg_ns%core%grain(p)     = veg_ns%core%grain(p)     - veg_nf%core%grain_to_food(p)*dt

                  veg_ns%core%cropseed_deficit(p) = veg_ns%core%cropseed_deficit(p) &
                        - veg_nf%core%crop_seed_to_leaf(p) * dt
              end if

              ! uptake from soil mineral N pool
              veg_ns%core%pool(p) = &
                    veg_ns%core%pool(p) + veg_nf%uptake%smin_to_pool(p)*dt
              if (nu_com .ne. 'RD') veg_ns%core%pool(p) = veg_ns%core%pool(p) + veg_nf%uptake%supplement_to_plant(p)*dt
              if (NFIX_PTASE_plant) veg_ns%core%pool(p) = veg_ns%core%pool(p) + veg_nf%nfix_to_plantn(p)*dt

              ! deployment from retranslocation pool
              veg_ns%core%pool(p)    = veg_ns%core%pool(p)    + veg_nf%retrans%retrans_to_pool(p)*dt
              veg_ns%retrans%pool(p) = veg_ns%retrans%pool(p) - veg_nf%retrans%retrans_to_pool(p)*dt

              ! allocation fluxes
              veg_ns%core%pool(p)           = veg_ns%core%pool(p)          - veg_nf%core%pool_to_leaf(p)*dt
              veg_ns%core%leaf(p)           = veg_ns%core%leaf(p)          + veg_nf%core%pool_to_leaf(p)*dt
              veg_ns%core%pool(p)           = veg_ns%core%pool(p)          - veg_nf%core%pool_to_leaf_storage(p)*dt
              veg_ns%core%leaf_storage(p)   = veg_ns%core%leaf_storage(p)  + veg_nf%core%pool_to_leaf_storage(p)*dt
              veg_ns%core%pool(p)           = veg_ns%core%pool(p)          - veg_nf%core%pool_to_froot(p)*dt
              veg_ns%core%froot(p)          = veg_ns%core%froot(p)         + veg_nf%core%pool_to_froot(p)*dt
              veg_ns%core%pool(p)           = veg_ns%core%pool(p)          - veg_nf%core%pool_to_froot_storage(p)*dt
              veg_ns%core%froot_storage(p)  = veg_ns%core%froot_storage(p) + veg_nf%core%pool_to_froot_storage(p)*dt

              if (woody(ivt(p)) >= 1.0_r8) then
                  veg_ns%core%pool(p)              = veg_ns%core%pool(p)              - veg_nf%core%pool_to_livestem(p)*dt
                  veg_ns%core%livestem(p)          = veg_ns%core%livestem(p)          + veg_nf%core%pool_to_livestem(p)*dt
                  veg_ns%core%pool(p)              = veg_ns%core%pool(p)              - veg_nf%core%pool_to_livestem_storage(p)*dt
                  veg_ns%core%livestem_storage(p)  = veg_ns%core%livestem_storage(p)  + veg_nf%core%pool_to_livestem_storage(p)*dt
                  veg_ns%core%pool(p)              = veg_ns%core%pool(p)              - veg_nf%core%pool_to_deadstem(p)*dt
                  veg_ns%core%deadstem(p)          = veg_ns%core%deadstem(p)          + veg_nf%core%pool_to_deadstem(p)*dt
                  veg_ns%core%pool(p)              = veg_ns%core%pool(p)              - veg_nf%core%pool_to_deadstem_storage(p)*dt
                  veg_ns%core%deadstem_storage(p)  = veg_ns%core%deadstem_storage(p)  + veg_nf%core%pool_to_deadstem_storage(p)*dt
                  veg_ns%core%pool(p)              = veg_ns%core%pool(p)              - veg_nf%core%pool_to_livecroot(p)*dt
                  veg_ns%core%livecroot(p)         = veg_ns%core%livecroot(p)         + veg_nf%core%pool_to_livecroot(p)*dt
                  veg_ns%core%pool(p)              = veg_ns%core%pool(p)              - veg_nf%core%pool_to_livecroot_storage(p)*dt
                  veg_ns%core%livecroot_storage(p) = veg_ns%core%livecroot_storage(p) + veg_nf%core%pool_to_livecroot_storage(p)*dt
                  veg_ns%core%pool(p)              = veg_ns%core%pool(p)              - veg_nf%core%pool_to_deadcroot(p)*dt
                  veg_ns%core%deadcroot(p)         = veg_ns%core%deadcroot(p)         + veg_nf%core%pool_to_deadcroot(p)*dt
                  veg_ns%core%pool(p)              = veg_ns%core%pool(p)              - veg_nf%core%pool_to_deadcroot_storage(p)*dt
                  veg_ns%core%deadcroot_storage(p) = veg_ns%core%deadcroot_storage(p) + veg_nf%core%pool_to_deadcroot_storage(p)*dt
              end if

              if (iscft(ivt(p))) then ! skip 2 generic crops
                  veg_ns%core%pool(p)              = veg_ns%core%pool(p)              - veg_nf%core%pool_to_livestem(p)*dt
                  veg_ns%core%livestem(p)          = veg_ns%core%livestem(p)          + veg_nf%core%pool_to_livestem(p)*dt
                  veg_ns%core%pool(p)              = veg_ns%core%pool(p)              - veg_nf%core%pool_to_livestem_storage(p)*dt
                  veg_ns%core%livestem_storage(p)  = veg_ns%core%livestem_storage(p)  + veg_nf%core%pool_to_livestem_storage(p)*dt
                  veg_ns%core%pool(p)              = veg_ns%core%pool(p)              - veg_nf%core%pool_to_grain(p)*dt
                  veg_ns%core%grain(p)             = veg_ns%core%grain(p)             + veg_nf%core%pool_to_grain(p)*dt
                  veg_ns%core%pool(p)              = veg_ns%core%pool(p)              - veg_nf%core%pool_to_grain_storage(p)*dt
                  veg_ns%core%grain_storage(p)     = veg_ns%core%grain_storage(p)     + veg_nf%core%pool_to_grain_storage(p)*dt
              end if

              ! move storage pools into transfer pools
              veg_ns%core%leaf_storage(p)  = veg_ns%core%leaf_storage(p)  - veg_nf%core%leaf_storage_to_xfer(p)*dt
              veg_ns%core%leaf_xfer(p)     = veg_ns%core%leaf_xfer(p)     + veg_nf%core%leaf_storage_to_xfer(p)*dt
              veg_ns%core%froot_storage(p) = veg_ns%core%froot_storage(p) - veg_nf%core%froot_storage_to_xfer(p)*dt
              veg_ns%core%froot_xfer(p)    = veg_ns%core%froot_xfer(p)    + veg_nf%core%froot_storage_to_xfer(p)*dt

              if (woody(ivt(p)) >= 1.0_r8) then
                  veg_ns%core%livestem_storage(p)  = veg_ns%core%livestem_storage(p)  - veg_nf%core%livestem_storage_to_xfer(p)*dt
                  veg_ns%core%livestem_xfer(p)     = veg_ns%core%livestem_xfer(p)     + veg_nf%core%livestem_storage_to_xfer(p)*dt
                  veg_ns%core%deadstem_storage(p)  = veg_ns%core%deadstem_storage(p)  - veg_nf%core%deadstem_storage_to_xfer(p)*dt
                  veg_ns%core%deadstem_xfer(p)     = veg_ns%core%deadstem_xfer(p)     + veg_nf%core%deadstem_storage_to_xfer(p)*dt
                  veg_ns%core%livecroot_storage(p) = veg_ns%core%livecroot_storage(p) - veg_nf%core%livecroot_storage_to_xfer(p)*dt
                  veg_ns%core%livecroot_xfer(p)    = veg_ns%core%livecroot_xfer(p)    + veg_nf%core%livecroot_storage_to_xfer(p)*dt
                  veg_ns%core%deadcroot_storage(p) = veg_ns%core%deadcroot_storage(p) - veg_nf%core%deadcroot_storage_to_xfer(p)*dt
                  veg_ns%core%deadcroot_xfer(p)    = veg_ns%core%deadcroot_xfer(p)    + veg_nf%core%deadcroot_storage_to_xfer(p)*dt
              end if

              if (iscft(ivt(p))) then ! skip 2 generic crops
                  ! lines here for consistency; the transfer terms are zero
                  veg_ns%core%livestem_storage(p)  = veg_ns%core%livestem_storage(p) - veg_nf%core%livestem_storage_to_xfer(p)*dt
                  veg_ns%core%livestem_xfer(p)     = veg_ns%core%livestem_xfer(p)    + veg_nf%core%livestem_storage_to_xfer(p)*dt
                  veg_ns%core%grain_storage(p)     = veg_ns%core%grain_storage(p)    - veg_nf%core%grain_storage_to_xfer(p)*dt
                  veg_ns%core%grain_xfer(p)        = veg_ns%core%grain_xfer(p)       + veg_nf%core%grain_storage_to_xfer(p)*dt
              end if

          end do  ! num_soilp

      end if ! if(.not.use_fates)

    end associate

  end subroutine NitrogenStateUpdate1

end module NitrogenStateUpdate1Mod
