module PStateUpdate1Mod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for phosphorus state variable updates, non-mortality fluxes.
  ! X.YANG
  ! !USES:
  use shr_kind_mod           , only: r8 => shr_kind_r8
  use clm_time_manager       , only : get_step_size
  use clm_varpar             , only : nlevdecomp, ndecomp_pools, ndecomp_cascade_transitions
  use clm_varpar             , only : crop_prog, i_met_lit, i_cel_lit, i_lig_lit, i_cwd
  use clm_varctl             , only : iulog, use_nitrif_denitrif
  use clm_varcon             , only : nitrif_n2o_loss_frac
  use pftvarcon              , only : npcropmin, nc3crop
  use soilorder_varcon       , only : smax,ks_sorption
  use VegetationPropertiesType         , only : veg_vp
  use CNDecompCascadeConType , only : decomp_cascade_con
  use CNStateType            , only : cnstate_type
  use PhosphorusFluxType     , only : phosphorusflux_type
  use PhosphorusStateType    , only : phosphorusstate_type
  use VegetationType              , only : veg_pp
  use tracer_varcon          , only : is_active_betr_bgc
  ! bgc interface & pflotran:
  use clm_varctl             , only : use_pflotran, pf_cmode
  use clm_varctl             , only : nu_com
  ! forest fertilization experiment
  use clm_time_manager       , only : get_curr_date
  use CNStateType            , only : fert_type , fert_continue, fert_dose, fert_start, fert_end
  use clm_varctl             , only : forest_fert_exp
  use decompMod              , only : bounds_type
  use clm_varcon             , only : dzsoi_decomp
  use clm_varctl             , only : use_fates
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: PStateUpdateDynPatch
  public :: PStateUpdate1
  !-----------------------------------------------------------------------

contains
  subroutine PStateUpdateDynPatch(bounds, num_soilc_with_inactive, filter_soilc_with_inactive, &
       phosphorusflux_vars, phosphorusstate_vars)
    !
    ! !DESCRIPTION:
    ! Update phosphorus states based on fluxes from dyn_cnbal_patch
    !
    ! !ARGUMENTS:
    type(bounds_type)          , intent(in)    :: bounds
    integer                    , intent(in)    :: num_soilc_with_inactive       ! number of columns in soil filter
    integer                    , intent(in)    :: filter_soilc_with_inactive(:) ! soil column filter that includes inactive points
    type(phosphorusflux_type)  , intent(in)    :: phosphorusflux_vars
    type(phosphorusstate_type) , intent(inout) :: phosphorusstate_vars
    !
    ! !LOCAL VARIABLES:
    integer                                    :: c                             ! column index
    integer                                    :: fc                            ! column filter index
    integer                                    :: g                             ! gridcell index
    integer                                    :: j                             ! level index
    real(r8)                                   :: dt                            ! time step (seconds)

    character(len=*)           , parameter     :: subname = 'PStateUpdateDynPatch'
    !-----------------------------------------------------------------------

    associate( &
         pf => phosphorusflux_vars  , &
         ps => phosphorusstate_vars   &
         )

      dt = real( get_step_size(), r8 )

      if (.not.use_fates) then

         do g = bounds%begg, bounds%endg
            ps%seed_grc(g) = ps%seed_grc(g) &
                 - pf%dwt_seed_to_leaf_grc(g)     * dt &
                 - pf%dwt_seed_to_deadstem_grc(g) * dt &
                 - pf%dwt_seedp_to_ppool_grc(g)    * dt
         end do

         do j = 1,nlevdecomp
            do fc = 1, num_soilc_with_inactive
               c = filter_soilc_with_inactive(fc)

               ps%decomp_pools_vr_col(c,j,i_met_lit) = ps%decomp_pools_vr_col(c,j,i_met_lit) + &
                    pf%dwt_froot_to_litr_met_col(c,j) * dt
               ps%decomp_pools_vr_col(c,j,i_cel_lit) = ps%decomp_pools_vr_col(c,j,i_cel_lit) + &
                    pf%dwt_froot_to_litr_cel_col(c,j) * dt
               ps%decomp_pools_vr_col(c,j,i_lig_lit) = ps%decomp_pools_vr_col(c,j,i_lig_lit) + &
                    pf%dwt_froot_to_litr_lig_col(c,j) * dt
               ps%decomp_pools_vr_col(c,j,i_cwd) = ps%decomp_pools_vr_col(c,j,i_cwd) + &
                    ( pf%dwt_livecroot_to_cwd_col(c,j) + pf%dwt_deadcroot_to_cwd_col(c,j) ) * dt

            end do
         end do
      end if

    end associate

  end subroutine PStateUpdateDynPatch

  !-----------------------------------------------------------------------
  subroutine PStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
       cnstate_vars, phosphorusflux_vars, phosphorusstate_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update all the prognostic phosphorus state
    ! variables (except for gap-phase mortality and fire fluxes)
    !
    ! !ARGUMENTS:
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                  , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnstate_type)       , intent(in)    :: cnstate_vars
    type(phosphorusflux_type)  , intent(inout) :: phosphorusflux_vars
    type(phosphorusstate_type) , intent(inout) :: phosphorusstate_vars
    !
    ! !LOCAL VARIABLES:
    integer :: c,p,j,l,k ! indices
    integer :: fp,fc     ! lake filter indices
    real(r8):: dt        ! radiation time step (seconds)

    integer:: kyr                     ! current year 
    integer:: kmo                     ! month of year  (1, ..., 12)
    integer:: kda                     ! day of month   (1, ..., 31) 
    integer:: mcsec                   ! seconds of day (0, ..., seconds/day)
    !-----------------------------------------------------------------------

    associate(                                                                                           & 
         ivt                   => veg_pp%itype                                , & ! Input:  [integer  (:)     ]  pft vegetation type                                

         woody                 => veg_vp%woody                         , & ! Input:  [real(r8) (:)     ]  binary flag for woody lifeform (1=woody, 0=not woody)

         cascade_donor_pool    => decomp_cascade_con%cascade_donor_pool    , & ! Input:  [integer  (:)     ]  which pool is C taken from for a given decomposition step
         cascade_receiver_pool => decomp_cascade_con%cascade_receiver_pool , & ! Input:  [integer  (:)     ]  which pool is C added to for a given decomposition step

         !!! N deposition profile, will weathering profile be needed?  -X.YANG
         ndep_prof             => cnstate_vars%ndep_prof_col               , & ! Input:  [real(r8) (:,:)   ]  profile over which N deposition is distributed through column (1/m)
!         nfixation_prof        => cnstate_vars%nfixation_prof_col          , & ! Input:  [real(r8) (:,:)   ]  profile over which N fixation is distributed through column (1/m)
         
         pf                    => phosphorusflux_vars                        , &
         ps                    => phosphorusstate_vars &
         )

      ! set time steps
      dt = real( get_step_size(), r8 )

      ! column-level fluxes


      !------------------------------------------------------------------
      ! if coupled with pflotran, the following updates are NOT needed
      ! if (.not.(use_pflotran .and. pf_cmode)) then
      !------------------------------------------------------------------
      if(.not. is_active_betr_bgc)then
      do j = 1, nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)

            ! plant to litter fluxes
            ! phenology and dynamic landcover fluxes
            pf%decomp_pools_sourcesink_col(c,j,i_met_lit) = &
                 pf%phenology_to_litr_met_col(c,j) * dt

            pf%decomp_pools_sourcesink_col(c,j,i_cel_lit) = &
                 pf%phenology_to_litr_cel_col(c,j) * dt

            pf%decomp_pools_sourcesink_col(c,j,i_lig_lit) = &
                 pf%phenology_to_litr_lig_col(c,j) * dt

         end do
      end do


      ! decomposition fluxes
      do k = 1, ndecomp_cascade_transitions
         do j = 1, nlevdecomp
            ! column loop
            do fc = 1,num_soilc
               c = filter_soilc(fc)

               pf%decomp_pools_sourcesink_col(c,j,cascade_donor_pool(k)) = &
                    pf%decomp_pools_sourcesink_col(c,j,cascade_donor_pool(k)) - &
                    pf%decomp_cascade_ptransfer_vr_col(c,j,k) * dt
            end do
         end do
      end do
      do k = 1, ndecomp_cascade_transitions
         if ( cascade_receiver_pool(k) /= 0 ) then  ! skip terminal transitions
            do j = 1, nlevdecomp
               ! column loop
               do fc = 1,num_soilc
                  c = filter_soilc(fc)

                  pf%decomp_pools_sourcesink_col(c,j,cascade_receiver_pool(k)) = &
                       pf%decomp_pools_sourcesink_col(c,j,cascade_receiver_pool(k)) + &
                       (pf%decomp_cascade_ptransfer_vr_col(c,j,k) + pf%decomp_cascade_sminp_flux_vr_col(c,j,k)) * dt
               end do
            end do
         else  ! terminal transitions
            do j = 1, nlevdecomp
               ! column loop
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  pf%decomp_pools_sourcesink_col(c,j,cascade_donor_pool(k)) = &
                       pf%decomp_pools_sourcesink_col(c,j,cascade_donor_pool(k)) - &
                       pf%decomp_cascade_sminp_flux_vr_col(c,j,k) * dt
               end do
            end do
         end if
      end do
      endif ! if (.not. is_active_betr_bgc))
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
                  ps%solutionp_vr_col(c,j) = ps%solutionp_vr_col(c,j) + fert_dose(c,kmo)*ndep_prof(c,j)
               end do
            end if
         end do
      end if


      ! patch loop

      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! phenology: transfer growth fluxes
         ps%leaf_patch(p)       = ps%leaf_patch(p)       + pf%leaf_xfer_to_leaf_patch(p)*dt
         ps%leaf_xfer_patch(p)  = ps%leaf_xfer_patch(p)  - pf%leaf_xfer_to_leaf_patch(p)*dt
         ps%froot_patch(p)      = ps%froot_patch(p)      + pf%froot_xfer_to_froot_patch(p)*dt
         ps%froot_xfer_patch(p) = ps%froot_xfer_patch(p) - pf%froot_xfer_to_froot_patch(p)*dt

         if (woody(ivt(p)) == 1.0_r8) then
            ps%livestem_patch(p)       = ps%livestem_patch(p)       + pf%livestem_xfer_to_livestem_patch(p)*dt
            ps%livestem_xfer_patch(p)  = ps%livestem_xfer_patch(p)  - pf%livestem_xfer_to_livestem_patch(p)*dt
            ps%deadstem_patch(p)       = ps%deadstem_patch(p)       + pf%deadstem_xfer_to_deadstem_patch(p)*dt
            ps%deadstem_xfer_patch(p)  = ps%deadstem_xfer_patch(p)  - pf%deadstem_xfer_to_deadstem_patch(p)*dt
            ps%livecroot_patch(p)      = ps%livecroot_patch(p)      + pf%livecroot_xfer_to_livecroot_patch(p)*dt
            ps%livecroot_xfer_patch(p) = ps%livecroot_xfer_patch(p) - pf%livecroot_xfer_to_livecroot_patch(p)*dt
            ps%deadcroot_patch(p)      = ps%deadcroot_patch(p)      + pf%deadcroot_xfer_to_deadcroot_patch(p)*dt
            ps%deadcroot_xfer_patch(p) = ps%deadcroot_xfer_patch(p) - pf%deadcroot_xfer_to_deadcroot_patch(p)*dt
         end if

         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            ! lines here for consistency; the transfer terms are zero
            ps%livestem_patch(p)       = ps%livestem_patch(p)      + pf%livestem_xfer_to_livestem_patch(p)*dt
            ps%livestem_xfer_patch(p)  = ps%livestem_xfer_patch(p) - pf%livestem_xfer_to_livestem_patch(p)*dt
            ps%grain_patch(p)          = ps%grain_patch(p)         + pf%grain_xfer_to_grain_patch(p)*dt
            ps%grain_xfer_patch(p)     = ps%grain_xfer_patch(p)    - pf%grain_xfer_to_grain_patch(p)*dt
         end if

         ! phenology: litterfall and retranslocation fluxes
         ps%leaf_patch(p)    = ps%leaf_patch(p)    - pf%leaf_to_litter_patch(p)*dt
         ps%froot_patch(p)   = ps%froot_patch(p)   - pf%froot_to_litter_patch(p)*dt
         ps%leaf_patch(p)    = ps%leaf_patch(p)    - pf%leafp_to_retransp_patch(p)*dt
         ps%retransp_patch(p) = ps%retransp_patch(p) + pf%leafp_to_retransp_patch(p)*dt

         ! live wood turnover and retranslocation fluxes
         if (woody(ivt(p)) == 1._r8) then
            ps%livestem_patch(p)  = ps%livestem_patch(p)  - pf%livestem_to_deadstem_patch(p)*dt
            ps%deadstem_patch(p)  = ps%deadstem_patch(p)  + pf%livestem_to_deadstem_patch(p)*dt
            ps%livestem_patch(p)  = ps%livestem_patch(p)  - pf%livestemp_to_retransp_patch(p)*dt
            ps%retransp_patch(p)   = ps%retransp_patch(p)   + pf%livestemp_to_retransp_patch(p)*dt
            ps%livecroot_patch(p) = ps%livecroot_patch(p) - pf%livecroot_to_deadcroot_patch(p)*dt
            ps%deadcroot_patch(p) = ps%deadcroot_patch(p) + pf%livecroot_to_deadcroot_patch(p)*dt
            ps%livecroot_patch(p) = ps%livecroot_patch(p) - pf%livecrootp_to_retransp_patch(p)*dt
            ps%retransp_patch(p)   = ps%retransp_patch(p)   + pf%livecrootp_to_retransp_patch(p)*dt
         end if
         if (ivt(p) >= npcropmin) then ! Beth adds retrans from froot
            ps%froot_patch(p)     = ps%froot_patch(p)     - pf%frootp_to_retransp_patch(p)*dt
            ps%retransp_patch(p)   = ps%retransp_patch(p)   + pf%frootp_to_retransp_patch(p)*dt
            ps%livestem_patch(p)  = ps%livestem_patch(p)  - pf%livestem_to_litter_patch(p)*dt
            ps%livestem_patch(p)  = ps%livestem_patch(p)  - pf%livestemp_to_retransp_patch(p)*dt
            ps%retransp_patch(p)   = ps%retransp_patch(p)   + pf%livestemp_to_retransp_patch(p)*dt
            ps%grain_patch(p)     = ps%grain_patch(p)     - pf%grain_to_food_patch(p)*dt

            ps%cropseed_deficit_patch(p) = ps%cropseed_deficit_patch(p) &
                 - pf%crop_seedp_to_leaf_patch(p) * dt
         end if

         ! uptake from soil mineral N pool
         ps%pool_patch(p) = &
              ps%pool_patch(p) + pf%sminp_to_ppool_patch(p)*dt
         if (nu_com .ne. 'RD') ps%pool_patch(p) = ps%pool_patch(p) + pf%supplement_to_plantp(p)*dt

         ! deployment from retranslocation pool
         ps%pool_patch(p)    = ps%pool_patch(p)    + pf%retransp_to_ppool_patch(p)*dt
         ps%retransp_patch(p) = ps%retransp_patch(p) - pf%retransp_to_ppool_patch(p)*dt

         ! allocation fluxes
         ps%pool_patch(p)           = ps%pool_patch(p)          - pf%pool_to_leaf_patch(p)*dt
         ps%leaf_patch(p)           = ps%leaf_patch(p)          + pf%pool_to_leaf_patch(p)*dt
         ps%pool_patch(p)           = ps%pool_patch(p)          - pf%pool_to_leaf_storage_patch(p)*dt
         ps%leaf_storage_patch(p)   = ps%leaf_storage_patch(p)  + pf%pool_to_leaf_storage_patch(p)*dt
         ps%pool_patch(p)           = ps%pool_patch(p)          - pf%pool_to_froot_patch(p)*dt
         ps%froot_patch(p)          = ps%froot_patch(p)         + pf%pool_to_froot_patch(p)*dt
         ps%pool_patch(p)           = ps%pool_patch(p)          - pf%pool_to_froot_storage_patch(p)*dt
         ps%froot_storage_patch(p)  = ps%froot_storage_patch(p) + pf%pool_to_froot_storage_patch(p)*dt

         if (woody(ivt(p)) == 1._r8) then
            ps%pool_patch(p)              = ps%pool_patch(p)              - pf%pool_to_livestem_patch(p)*dt
            ps%livestem_patch(p)          = ps%livestem_patch(p)          + pf%pool_to_livestem_patch(p)*dt
            ps%pool_patch(p)              = ps%pool_patch(p)              - pf%pool_to_livestem_storage_patch(p)*dt
            ps%livestem_storage_patch(p)  = ps%livestem_storage_patch(p)  + pf%pool_to_livestem_storage_patch(p)*dt
            ps%pool_patch(p)              = ps%pool_patch(p)              - pf%pool_to_deadstem_patch(p)*dt
            ps%deadstem_patch(p)          = ps%deadstem_patch(p)          + pf%pool_to_deadstem_patch(p)*dt
            ps%pool_patch(p)              = ps%pool_patch(p)              - pf%pool_to_deadstem_storage_patch(p)*dt
            ps%deadstem_storage_patch(p)  = ps%deadstem_storage_patch(p)  + pf%pool_to_deadstem_storage_patch(p)*dt
            ps%pool_patch(p)              = ps%pool_patch(p)              - pf%pool_to_livecroot_patch(p)*dt
            ps%livecroot_patch(p)         = ps%livecroot_patch(p)         + pf%pool_to_livecroot_patch(p)*dt
            ps%pool_patch(p)              = ps%pool_patch(p)              - pf%pool_to_livecroot_storage_patch(p)*dt
            ps%livecroot_storage_patch(p) = ps%livecroot_storage_patch(p) + pf%pool_to_livecroot_storage_patch(p)*dt
            ps%pool_patch(p)              = ps%pool_patch(p)              - pf%pool_to_deadcroot_patch(p)*dt
            ps%deadcroot_patch(p)         = ps%deadcroot_patch(p)         + pf%pool_to_deadcroot_patch(p)*dt
            ps%pool_patch(p)              = ps%pool_patch(p)              - pf%pool_to_deadcroot_storage_patch(p)*dt
            ps%deadcroot_storage_patch(p) = ps%deadcroot_storage_patch(p) + pf%pool_to_deadcroot_storage_patch(p)*dt
         end if

         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            ps%pool_patch(p)              = ps%pool_patch(p)              - pf%pool_to_livestem_patch(p)*dt
            ps%livestem_patch(p)          = ps%livestem_patch(p)          + pf%pool_to_livestem_patch(p)*dt
            ps%pool_patch(p)              = ps%pool_patch(p)              - pf%pool_to_livestem_storage_patch(p)*dt
            ps%livestem_storage_patch(p)  = ps%livestem_storage_patch(p)  + pf%pool_to_livestem_storage_patch(p)*dt
            ps%pool_patch(p)              = ps%pool_patch(p)              - pf%pool_to_grain_patch(p)*dt
            ps%grain_patch(p)             = ps%grain_patch(p)             + pf%pool_to_grain_patch(p)*dt
            ps%pool_patch(p)              = ps%pool_patch(p)              - pf%pool_to_grain_storage_patch(p)*dt
            ps%grain_storage_patch(p)     = ps%grain_storage_patch(p)     + pf%pool_to_grain_storage_patch(p)*dt
         end if

         ! move storage pools into transfer pools
         ps%leaf_storage_patch(p)  = ps%leaf_storage_patch(p)  - pf%leaf_storage_to_xfer_patch(p)*dt
         ps%leaf_xfer_patch(p)     = ps%leaf_xfer_patch(p)     + pf%leaf_storage_to_xfer_patch(p)*dt
         ps%froot_storage_patch(p) = ps%froot_storage_patch(p) - pf%froot_storage_to_xfer_patch(p)*dt
         ps%froot_xfer_patch(p)    = ps%froot_xfer_patch(p)    + pf%froot_storage_to_xfer_patch(p)*dt

         if (woody(ivt(p)) == 1._r8) then
            ps%livestem_storage_patch(p)  = ps%livestem_storage_patch(p)  - pf%livestem_storage_to_xfer_patch(p)*dt
            ps%livestem_xfer_patch(p)     = ps%livestem_xfer_patch(p)     + pf%livestem_storage_to_xfer_patch(p)*dt
            ps%deadstem_storage_patch(p)  = ps%deadstem_storage_patch(p)  - pf%deadstem_storage_to_xfer_patch(p)*dt
            ps%deadstem_xfer_patch(p)     = ps%deadstem_xfer_patch(p)     + pf%deadstem_storage_to_xfer_patch(p)*dt
            ps%livecroot_storage_patch(p) = ps%livecroot_storage_patch(p) - pf%livecroot_storage_to_xfer_patch(p)*dt
            ps%livecroot_xfer_patch(p)    = ps%livecroot_xfer_patch(p)    + pf%livecroot_storage_to_xfer_patch(p)*dt
            ps%deadcroot_storage_patch(p) = ps%deadcroot_storage_patch(p) - pf%deadcroot_storage_to_xfer_patch(p)*dt
            ps%deadcroot_xfer_patch(p)    = ps%deadcroot_xfer_patch(p)    + pf%deadcroot_storage_to_xfer_patch(p)*dt
         end if

         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            ! lines here for consistency; the transfer terms are zero
            ps%livestem_storage_patch(p)  = ps%livestem_storage_patch(p) - pf%livestem_storage_to_xfer_patch(p)*dt
            ps%livestem_xfer_patch(p)     = ps%livestem_xfer_patch(p)    + pf%livestem_storage_to_xfer_patch(p)*dt
            ps%grain_storage_patch(p)     = ps%grain_storage_patch(p)    - pf%grain_storage_to_xfer_patch(p)*dt
            ps%grain_xfer_patch(p)        = ps%grain_xfer_patch(p)       + pf%grain_storage_to_xfer_patch(p)*dt
         end if

      end do

    end associate

  end subroutine PStateUpdate1

end module PStateUpdate1Mod
