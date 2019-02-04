module PhosphorusStateUpdate1Mod
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
  use clm_varctl             , only : NFIX_PTASE_plant
  use decompMod              , only : bounds_type
  use clm_varcon             , only : dzsoi_decomp
  use clm_varctl             , only : use_fates
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
  subroutine PhosphorusStateUpdateDynPatch(bounds, num_soilc_with_inactive, filter_soilc_with_inactive, &
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

    character(len=*)           , parameter     :: subname = 'PhosphorusStateUpdateDynPatch'
    !-----------------------------------------------------------------------

    associate( &
         pf => phosphorusflux_vars  , &
         ps => phosphorusstate_vars   &
         )

      dt = real( get_step_size(), r8 )

      if (.not.use_fates) then

         do g = bounds%begg, bounds%endg
            ps%seedp_grc(g) = ps%seedp_grc(g) &
                 - pf%dwt_seedp_to_leaf_grc(g)     * dt &
                 - pf%dwt_seedp_to_deadstem_grc(g) * dt &
                 - pf%dwt_seedp_to_ppool_grc(g)    * dt
         end do

         do j = 1,nlevdecomp
            do fc = 1, num_soilc_with_inactive
               c = filter_soilc_with_inactive(fc)

               ps%decomp_ppools_vr_col(c,j,i_met_lit) = ps%decomp_ppools_vr_col(c,j,i_met_lit) + &
                    pf%dwt_frootp_to_litr_met_p_col(c,j) * dt
               ps%decomp_ppools_vr_col(c,j,i_cel_lit) = ps%decomp_ppools_vr_col(c,j,i_cel_lit) + &
                    pf%dwt_frootp_to_litr_cel_p_col(c,j) * dt
               ps%decomp_ppools_vr_col(c,j,i_lig_lit) = ps%decomp_ppools_vr_col(c,j,i_lig_lit) + &
                    pf%dwt_frootp_to_litr_lig_p_col(c,j) * dt
               ps%decomp_ppools_vr_col(c,j,i_cwd) = ps%decomp_ppools_vr_col(c,j,i_cwd) + &
                    ( pf%dwt_livecrootp_to_cwdp_col(c,j) + pf%dwt_deadcrootp_to_cwdp_col(c,j) ) * dt

            end do
         end do
      end if

    end associate

  end subroutine PhosphorusStateUpdateDynPatch

  !-----------------------------------------------------------------------
  subroutine PhosphorusStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
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
            pf%decomp_ppools_sourcesink_col(c,j,i_met_lit) = &
                 pf%phenology_p_to_litr_met_p_col(c,j) * dt

            pf%decomp_ppools_sourcesink_col(c,j,i_cel_lit) = &
                 pf%phenology_p_to_litr_cel_p_col(c,j) * dt

            pf%decomp_ppools_sourcesink_col(c,j,i_lig_lit) = &
                 pf%phenology_p_to_litr_lig_p_col(c,j) * dt

         end do
      end do


      ! decomposition fluxes
      do k = 1, ndecomp_cascade_transitions
         do j = 1, nlevdecomp
            ! column loop
            do fc = 1,num_soilc
               c = filter_soilc(fc)

               pf%decomp_ppools_sourcesink_col(c,j,cascade_donor_pool(k)) = &
                    pf%decomp_ppools_sourcesink_col(c,j,cascade_donor_pool(k)) - &
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

                  pf%decomp_ppools_sourcesink_col(c,j,cascade_receiver_pool(k)) = &
                       pf%decomp_ppools_sourcesink_col(c,j,cascade_receiver_pool(k)) + &
                       (pf%decomp_cascade_ptransfer_vr_col(c,j,k) + pf%decomp_cascade_sminp_flux_vr_col(c,j,k)) * dt
               end do
            end do
         else  ! terminal transitions
            do j = 1, nlevdecomp
               ! column loop
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  pf%decomp_ppools_sourcesink_col(c,j,cascade_donor_pool(k)) = &
                       pf%decomp_ppools_sourcesink_col(c,j,cascade_donor_pool(k)) - &
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
         ps%leafp_patch(p)       = ps%leafp_patch(p)       + pf%leafp_xfer_to_leafp_patch(p)*dt
         ps%leafp_xfer_patch(p)  = ps%leafp_xfer_patch(p)  - pf%leafp_xfer_to_leafp_patch(p)*dt
         ps%frootp_patch(p)      = ps%frootp_patch(p)      + pf%frootp_xfer_to_frootp_patch(p)*dt
         ps%frootp_xfer_patch(p) = ps%frootp_xfer_patch(p) - pf%frootp_xfer_to_frootp_patch(p)*dt

         if (woody(ivt(p)) == 1.0_r8) then
            ps%livestemp_patch(p)       = ps%livestemp_patch(p)       + pf%livestemp_xfer_to_livestemp_patch(p)*dt
            ps%livestemp_xfer_patch(p)  = ps%livestemp_xfer_patch(p)  - pf%livestemp_xfer_to_livestemp_patch(p)*dt
            ps%deadstemp_patch(p)       = ps%deadstemp_patch(p)       + pf%deadstemp_xfer_to_deadstemp_patch(p)*dt
            ps%deadstemp_xfer_patch(p)  = ps%deadstemp_xfer_patch(p)  - pf%deadstemp_xfer_to_deadstemp_patch(p)*dt
            ps%livecrootp_patch(p)      = ps%livecrootp_patch(p)      + pf%livecrootp_xfer_to_livecrootp_patch(p)*dt
            ps%livecrootp_xfer_patch(p) = ps%livecrootp_xfer_patch(p) - pf%livecrootp_xfer_to_livecrootp_patch(p)*dt
            ps%deadcrootp_patch(p)      = ps%deadcrootp_patch(p)      + pf%deadcrootp_xfer_to_deadcrootp_patch(p)*dt
            ps%deadcrootp_xfer_patch(p) = ps%deadcrootp_xfer_patch(p) - pf%deadcrootp_xfer_to_deadcrootp_patch(p)*dt
         end if

         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            ! lines here for consistency; the transfer terms are zero
            ps%livestemp_patch(p)       = ps%livestemp_patch(p)      + pf%livestemp_xfer_to_livestemp_patch(p)*dt
            ps%livestemp_xfer_patch(p)  = ps%livestemp_xfer_patch(p) - pf%livestemp_xfer_to_livestemp_patch(p)*dt
            ps%grainp_patch(p)          = ps%grainp_patch(p)         + pf%grainp_xfer_to_grainp_patch(p)*dt
            ps%grainp_xfer_patch(p)     = ps%grainp_xfer_patch(p)    - pf%grainp_xfer_to_grainp_patch(p)*dt
         end if

         ! phenology: litterfall and retranslocation fluxes
         ps%leafp_patch(p)    = ps%leafp_patch(p)    - pf%leafp_to_litter_patch(p)*dt
         ps%frootp_patch(p)   = ps%frootp_patch(p)   - pf%frootp_to_litter_patch(p)*dt
         ps%leafp_patch(p)    = ps%leafp_patch(p)    - pf%leafp_to_retransp_patch(p)*dt
         ps%retransp_patch(p) = ps%retransp_patch(p) + pf%leafp_to_retransp_patch(p)*dt

         ! live wood turnover and retranslocation fluxes
         if (woody(ivt(p)) == 1._r8) then
            ps%livestemp_patch(p)  = ps%livestemp_patch(p)  - pf%livestemp_to_deadstemp_patch(p)*dt
            ps%deadstemp_patch(p)  = ps%deadstemp_patch(p)  + pf%livestemp_to_deadstemp_patch(p)*dt
            ps%livestemp_patch(p)  = ps%livestemp_patch(p)  - pf%livestemp_to_retransp_patch(p)*dt
            ps%retransp_patch(p)   = ps%retransp_patch(p)   + pf%livestemp_to_retransp_patch(p)*dt
            ps%livecrootp_patch(p) = ps%livecrootp_patch(p) - pf%livecrootp_to_deadcrootp_patch(p)*dt
            ps%deadcrootp_patch(p) = ps%deadcrootp_patch(p) + pf%livecrootp_to_deadcrootp_patch(p)*dt
            ps%livecrootp_patch(p) = ps%livecrootp_patch(p) - pf%livecrootp_to_retransp_patch(p)*dt
            ps%retransp_patch(p)   = ps%retransp_patch(p)   + pf%livecrootp_to_retransp_patch(p)*dt
         end if
         if (ivt(p) >= npcropmin) then ! Beth adds retrans from froot
            ps%frootp_patch(p)     = ps%frootp_patch(p)     - pf%frootp_to_retransp_patch(p)*dt
            ps%retransp_patch(p)   = ps%retransp_patch(p)   + pf%frootp_to_retransp_patch(p)*dt
            ps%livestemp_patch(p)  = ps%livestemp_patch(p)  - pf%livestemp_to_litter_patch(p)*dt
            ps%livestemp_patch(p)  = ps%livestemp_patch(p)  - pf%livestemp_to_retransp_patch(p)*dt
            ps%retransp_patch(p)   = ps%retransp_patch(p)   + pf%livestemp_to_retransp_patch(p)*dt
            ps%grainp_patch(p)     = ps%grainp_patch(p)     - pf%grainp_to_food_patch(p)*dt

            ps%cropseedp_deficit_patch(p) = ps%cropseedp_deficit_patch(p) &
                 - pf%crop_seedp_to_leaf_patch(p) * dt
         end if

         ! uptake from soil mineral N pool
         ps%ppool_patch(p) = &
              ps%ppool_patch(p) + pf%sminp_to_ppool_patch(p)*dt
         if (nu_com .ne. 'RD') ps%ppool_patch(p) = ps%ppool_patch(p) + pf%supplement_to_plantp(p)*dt
         if (NFIX_PTASE_plant) ps%ppool_patch(p) = ps%ppool_patch(p) + pf%biochem_pmin_to_plant_patch(p)*dt

         ! deployment from retranslocation pool
         ps%ppool_patch(p)    = ps%ppool_patch(p)    + pf%retransp_to_ppool_patch(p)*dt
         ps%retransp_patch(p) = ps%retransp_patch(p) - pf%retransp_to_ppool_patch(p)*dt

         ! allocation fluxes
         ps%ppool_patch(p)           = ps%ppool_patch(p)          - pf%ppool_to_leafp_patch(p)*dt
         ps%leafp_patch(p)           = ps%leafp_patch(p)          + pf%ppool_to_leafp_patch(p)*dt
         ps%ppool_patch(p)           = ps%ppool_patch(p)          - pf%ppool_to_leafp_storage_patch(p)*dt
         ps%leafp_storage_patch(p)   = ps%leafp_storage_patch(p)  + pf%ppool_to_leafp_storage_patch(p)*dt
         ps%ppool_patch(p)           = ps%ppool_patch(p)          - pf%ppool_to_frootp_patch(p)*dt
         ps%frootp_patch(p)          = ps%frootp_patch(p)         + pf%ppool_to_frootp_patch(p)*dt
         ps%ppool_patch(p)           = ps%ppool_patch(p)          - pf%ppool_to_frootp_storage_patch(p)*dt
         ps%frootp_storage_patch(p)  = ps%frootp_storage_patch(p) + pf%ppool_to_frootp_storage_patch(p)*dt

         if (woody(ivt(p)) == 1._r8) then
            ps%ppool_patch(p)              = ps%ppool_patch(p)              - pf%ppool_to_livestemp_patch(p)*dt
            ps%livestemp_patch(p)          = ps%livestemp_patch(p)          + pf%ppool_to_livestemp_patch(p)*dt
            ps%ppool_patch(p)              = ps%ppool_patch(p)              - pf%ppool_to_livestemp_storage_patch(p)*dt
            ps%livestemp_storage_patch(p)  = ps%livestemp_storage_patch(p)  + pf%ppool_to_livestemp_storage_patch(p)*dt
            ps%ppool_patch(p)              = ps%ppool_patch(p)              - pf%ppool_to_deadstemp_patch(p)*dt
            ps%deadstemp_patch(p)          = ps%deadstemp_patch(p)          + pf%ppool_to_deadstemp_patch(p)*dt
            ps%ppool_patch(p)              = ps%ppool_patch(p)              - pf%ppool_to_deadstemp_storage_patch(p)*dt
            ps%deadstemp_storage_patch(p)  = ps%deadstemp_storage_patch(p)  + pf%ppool_to_deadstemp_storage_patch(p)*dt
            ps%ppool_patch(p)              = ps%ppool_patch(p)              - pf%ppool_to_livecrootp_patch(p)*dt
            ps%livecrootp_patch(p)         = ps%livecrootp_patch(p)         + pf%ppool_to_livecrootp_patch(p)*dt
            ps%ppool_patch(p)              = ps%ppool_patch(p)              - pf%ppool_to_livecrootp_storage_patch(p)*dt
            ps%livecrootp_storage_patch(p) = ps%livecrootp_storage_patch(p) + pf%ppool_to_livecrootp_storage_patch(p)*dt
            ps%ppool_patch(p)              = ps%ppool_patch(p)              - pf%ppool_to_deadcrootp_patch(p)*dt
            ps%deadcrootp_patch(p)         = ps%deadcrootp_patch(p)         + pf%ppool_to_deadcrootp_patch(p)*dt
            ps%ppool_patch(p)              = ps%ppool_patch(p)              - pf%ppool_to_deadcrootp_storage_patch(p)*dt
            ps%deadcrootp_storage_patch(p) = ps%deadcrootp_storage_patch(p) + pf%ppool_to_deadcrootp_storage_patch(p)*dt
         end if

         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            ps%ppool_patch(p)              = ps%ppool_patch(p)              - pf%ppool_to_livestemp_patch(p)*dt
            ps%livestemp_patch(p)          = ps%livestemp_patch(p)          + pf%ppool_to_livestemp_patch(p)*dt
            ps%ppool_patch(p)              = ps%ppool_patch(p)              - pf%ppool_to_livestemp_storage_patch(p)*dt
            ps%livestemp_storage_patch(p)  = ps%livestemp_storage_patch(p)  + pf%ppool_to_livestemp_storage_patch(p)*dt
            ps%ppool_patch(p)              = ps%ppool_patch(p)              - pf%ppool_to_grainp_patch(p)*dt
            ps%grainp_patch(p)             = ps%grainp_patch(p)             + pf%ppool_to_grainp_patch(p)*dt
            ps%ppool_patch(p)              = ps%ppool_patch(p)              - pf%ppool_to_grainp_storage_patch(p)*dt
            ps%grainp_storage_patch(p)     = ps%grainp_storage_patch(p)     + pf%ppool_to_grainp_storage_patch(p)*dt
         end if

         ! move storage pools into transfer pools
         ps%leafp_storage_patch(p)  = ps%leafp_storage_patch(p)  - pf%leafp_storage_to_xfer_patch(p)*dt
         ps%leafp_xfer_patch(p)     = ps%leafp_xfer_patch(p)     + pf%leafp_storage_to_xfer_patch(p)*dt
         ps%frootp_storage_patch(p) = ps%frootp_storage_patch(p) - pf%frootp_storage_to_xfer_patch(p)*dt
         ps%frootp_xfer_patch(p)    = ps%frootp_xfer_patch(p)    + pf%frootp_storage_to_xfer_patch(p)*dt

         if (woody(ivt(p)) == 1._r8) then
            ps%livestemp_storage_patch(p)  = ps%livestemp_storage_patch(p)  - pf%livestemp_storage_to_xfer_patch(p)*dt
            ps%livestemp_xfer_patch(p)     = ps%livestemp_xfer_patch(p)     + pf%livestemp_storage_to_xfer_patch(p)*dt
            ps%deadstemp_storage_patch(p)  = ps%deadstemp_storage_patch(p)  - pf%deadstemp_storage_to_xfer_patch(p)*dt
            ps%deadstemp_xfer_patch(p)     = ps%deadstemp_xfer_patch(p)     + pf%deadstemp_storage_to_xfer_patch(p)*dt
            ps%livecrootp_storage_patch(p) = ps%livecrootp_storage_patch(p) - pf%livecrootp_storage_to_xfer_patch(p)*dt
            ps%livecrootp_xfer_patch(p)    = ps%livecrootp_xfer_patch(p)    + pf%livecrootp_storage_to_xfer_patch(p)*dt
            ps%deadcrootp_storage_patch(p) = ps%deadcrootp_storage_patch(p) - pf%deadcrootp_storage_to_xfer_patch(p)*dt
            ps%deadcrootp_xfer_patch(p)    = ps%deadcrootp_xfer_patch(p)    + pf%deadcrootp_storage_to_xfer_patch(p)*dt
         end if

         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            ! lines here for consistency; the transfer terms are zero
            ps%livestemp_storage_patch(p)  = ps%livestemp_storage_patch(p) - pf%livestemp_storage_to_xfer_patch(p)*dt
            ps%livestemp_xfer_patch(p)     = ps%livestemp_xfer_patch(p)    + pf%livestemp_storage_to_xfer_patch(p)*dt
            ps%grainp_storage_patch(p)     = ps%grainp_storage_patch(p)    - pf%grainp_storage_to_xfer_patch(p)*dt
            ps%grainp_xfer_patch(p)        = ps%grainp_xfer_patch(p)       + pf%grainp_storage_to_xfer_patch(p)*dt
         end if

      end do

    end associate

  end subroutine PhosphorusStateUpdate1

end module PhosphorusStateUpdate1Mod
