module CarbonStateUpdate1Mod

  !-----------------------------------------------------------------------
  ! Module for carbon state variable update, non-mortality fluxes.
  !
  ! !USES:
  use shr_kind_mod           , only : r8 => shr_kind_r8
  use shr_log_mod            , only : errMsg => shr_log_errMsg
  use clm_varpar             , only : ndecomp_cascade_transitions, nlevdecomp
  use clm_time_manager       , only : get_step_size
  use clm_varpar             , only : i_met_lit, i_cel_lit, i_lig_lit, i_cwd
  use pftvarcon              , only : npcropmin, nc3crop
  use abortutils             , only : endrun
  use CNDecompCascadeConType , only : decomp_cascade_type
  use CNCarbonStateType      , only : carbonstate_type
  use CNCarbonFluxType       , only : carbonflux_type
  use CNStateType            , only : cnstate_type
  use CNDecompCascadeConType , only : decomp_cascade_con
  use VegetationPropertiesType, only : veg_vp
  use clm_varctl             , only : nu_com
  use VegetationType         , only : veg_pp
  use CropType               , only : crop_type
  use decompMod              , only : bounds_type
  use clm_varcon             , only : dzsoi_decomp
  ! bgc interface & pflotran:
  use clm_varctl             , only : use_pflotran, pf_cmode, use_fates
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CarbonStateUpdateDynPatch
  public :: CarbonStateUpdate1
  public :: CarbonStateUpdate0
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CarbonStateUpdateDynPatch(bounds, num_soilc_with_inactive, filter_soilc_with_inactive, &
       carbonflux_vars, carbonstate_vars)
    !
    ! !DESCRIPTION:
    ! Update carbon states based on fluxes from dyn_cnbal_patch
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: num_soilc_with_inactive       ! number of columns in soil filter
    integer                , intent(in)    :: filter_soilc_with_inactive(:) ! soil column filter that includes inactive points
    type(carbonflux_type)  , intent(in)    :: carbonflux_vars
    type(carbonstate_type) , intent(inout) :: carbonstate_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: c   ! column index
    integer  :: fc  ! column filter index
    integer  :: g   ! gridcell index
    integer  :: j   ! level index
    real(r8) :: dt  ! time step (seconds)

    character(len=*), parameter :: subname = 'CarbonStateUpdateDynPatch'
    !-----------------------------------------------------------------------

    associate( &
         cf => carbonflux_vars  , &
         cs => carbonstate_vars   &
         )

      dt = real( get_step_size(), r8 )

      if (.not.use_fates) then

         do g = bounds%begg, bounds%endg
            cs%seedc_grc(g) = cs%seedc_grc(g) &
                 - cf%dwt_seedc_to_leaf_grc(g)     * dt &
                 - cf%dwt_seedc_to_deadstem_grc(g) * dt
         end do

         do j = 1,nlevdecomp
            do fc = 1, num_soilc_with_inactive
               c = filter_soilc_with_inactive(fc)

               cs%decomp_cpools_vr_col(c,j,i_met_lit) = cs%decomp_cpools_vr_col(c,j,i_met_lit) + &
                    cf%dwt_frootc_to_litr_met_c_col(c,j) * dt
               cs%decomp_cpools_vr_col(c,j,i_cel_lit) = cs%decomp_cpools_vr_col(c,j,i_cel_lit) + &
                    cf%dwt_frootc_to_litr_cel_c_col(c,j) * dt
               cs%decomp_cpools_vr_col(c,j,i_lig_lit) = cs%decomp_cpools_vr_col(c,j,i_lig_lit) + &
                    cf%dwt_frootc_to_litr_lig_c_col(c,j) * dt
               cs%decomp_cpools_vr_col(c,j,i_cwd) = cs%decomp_cpools_vr_col(c,j,i_cwd) + &
                    ( cf%dwt_livecrootc_to_cwdc_col(c,j) + cf%dwt_deadcrootc_to_cwdc_col(c,j) ) * dt

            end do
         end do
      end if

    end associate

  end subroutine CarbonStateUpdateDynPatch

  !-----------------------------------------------------------------------
  subroutine CarbonStateUpdate0(&
       num_soilp, filter_soilp, &
       carbonflux_vars, carbonstate_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update cpool carbon state
    !

    ! !ARGUMENTS:
    integer                , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(carbonflux_type)  , intent(in)    :: carbonflux_vars
    type(carbonstate_type) , intent(inout) :: carbonstate_vars
    !
    ! !LOCAL VARIABLES:
    integer :: p  ! indices
    integer :: fp ! lake filter indices
    real(r8):: dt ! radiation time step (seconds)
    !-----------------------------------------------------------------------

    associate(                                                          & 
         cf                => carbonflux_vars                         , &
         cs                => carbonstate_vars                          &
         )

      ! set time steps
      dt = real( get_step_size(), r8 )

      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)
         ! gross photosynthesis fluxes
         cs%cpool_patch(p) = cs%cpool_patch(p) + cf%psnsun_to_cpool_patch(p)*dt
         cs%cpool_patch(p) = cs%cpool_patch(p) + cf%psnshade_to_cpool_patch(p)*dt
      end do

    end associate

  end subroutine CarbonStateUpdate0

  !-----------------------------------------------------------------------
  subroutine CarbonStateUpdate1(bounds, &
       num_soilc, filter_soilc, &
       num_soilp, filter_soilp, &
       crop_vars, carbonflux_vars, carbonstate_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update all the prognostic carbon state
    ! variables (except for gap-phase mortality and fire fluxes)
    !
    use tracer_varcon       , only : is_active_betr_bgc
    use subgridAveMod       , only : p2c
    use decompMod           , only : bounds_type    
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds  
    integer                , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(crop_type)        , intent(inout) :: crop_vars
    type(carbonflux_type)  , intent(inout) :: carbonflux_vars
    type(carbonstate_type) , intent(inout) :: carbonstate_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: c,p,j,k,l ! indices
    integer  :: fp,fc     ! lake filter indices
    real(r8) :: dt        ! radiation time step (seconds)
    !-----------------------------------------------------------------------

    associate(                                                                                     & 
         ivt                           =>    veg_pp%itype                                           , & ! Input:  [integer  (:)     ]  pft vegetation type                                

         woody                         =>    veg_vp%woody                                    , & ! Input:  [real(r8) (:)     ]  binary flag for woody lifeform (1=woody, 0=not woody)
         cascade_donor_pool            =>    decomp_cascade_con%cascade_donor_pool               , & ! Input:  [integer  (:)     ]  which pool is C taken from for a given decomposition step
         cascade_receiver_pool         =>    decomp_cascade_con%cascade_receiver_pool            , & ! Input:  [integer  (:)     ]  which pool is C added to for a given decomposition step

         harvdate                      =>    crop_vars%harvdate_patch                         , & ! Input:  [integer  (:)     ]  harvest date                                       
         
         cf => carbonflux_vars  , &
         cs => carbonstate_vars   &

         )

      ! set time steps
      dt = real( get_step_size(), r8 )

      ! column level fluxes

      if(.not.use_fates) then
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            cs%decomp_som2c_vr_col(c,1:nlevdecomp) = cs%decomp_cpools_vr_col(c,1:nlevdecomp,6)
         end do
      end if


      
      if (.not. is_active_betr_bgc .and. .not.(use_pflotran .and. pf_cmode) .and. .not.use_fates ) then

         ! plant to litter fluxes

         do j = 1,nlevdecomp
            ! column loop
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               ! phenology and dynamic land cover fluxes
               cf%decomp_cpools_sourcesink_col(c,j,i_met_lit) = &
                    cf%phenology_c_to_litr_met_c_col(c,j) * dt
               cf%decomp_cpools_sourcesink_col(c,j,i_cel_lit) = &
                    cf%phenology_c_to_litr_cel_c_col(c,j) * dt
               cf%decomp_cpools_sourcesink_col(c,j,i_lig_lit) = &
                    cf%phenology_c_to_litr_lig_c_col(c,j) * dt
            end do
         end do

         ! litter and SOM HR fluxes
         do k = 1, ndecomp_cascade_transitions
            do j = 1,nlevdecomp
               ! column loop
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  cf%decomp_cpools_sourcesink_col(c,j,cascade_donor_pool(k)) = &
                       cf%decomp_cpools_sourcesink_col(c,j,cascade_donor_pool(k)) &
                       - ( cf%decomp_cascade_hr_vr_col(c,j,k) + cf%decomp_cascade_ctransfer_vr_col(c,j,k)) *dt
               end do
            end do
         end do
         do k = 1, ndecomp_cascade_transitions
            if ( cascade_receiver_pool(k) /= 0 ) then  ! skip terminal transitions
               do j = 1,nlevdecomp
                  ! column loop
                  do fc = 1,num_soilc
                     c = filter_soilc(fc)
                     cf%decomp_cpools_sourcesink_col(c,j,cascade_receiver_pool(k)) = &
                          cf%decomp_cpools_sourcesink_col(c,j,cascade_receiver_pool(k)) &
                          + cf%decomp_cascade_ctransfer_vr_col(c,j,k)*dt
                  end do
               end do
            end if
         end do

      elseif( use_fates ) then

         ! The following pools were updated via the FATES interface
         ! cf%decomp_cpools_sourcesink_col(c,j,i_met_lit)
         ! cf%decomp_cpools_sourcesink_col(c,j,i_cel_lit)
         ! cf%decomp_cpools_sourcesink_col(c,j,i_lig_lit)

         ! litter and SOM HR fluxes
         do k = 1, ndecomp_cascade_transitions
            do j = 1,nlevdecomp
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  cf%decomp_cpools_sourcesink_col(c,j,cascade_donor_pool(k)) = &
                        cf%decomp_cpools_sourcesink_col(c,j,cascade_donor_pool(k)) &
                        - ( cf%decomp_cascade_hr_vr_col(c,j,k) + cf%decomp_cascade_ctransfer_vr_col(c,j,k)) *dt
               end do
            end do
         end do
         do k = 1, ndecomp_cascade_transitions
            if ( cascade_receiver_pool(k) /= 0 ) then  ! skip terminal transitions
               do j = 1,nlevdecomp
                  do fc = 1,num_soilc
                     c = filter_soilc(fc)
                     cf%decomp_cpools_sourcesink_col(c,j,cascade_receiver_pool(k)) = &
                           cf%decomp_cpools_sourcesink_col(c,j,cascade_receiver_pool(k)) &
                           + cf%decomp_cascade_ctransfer_vr_col(c,j,k)*dt
                  end do
               end do
            end if
         end do

         
   endif   !end if is_active_betr_bgc()
    

   if (.not.use_fates) then
    
      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! phenology: transfer growth fluxes
         cs%leafc_patch(p)           = cs%leafc_patch(p)       + cf%leafc_xfer_to_leafc_patch(p)*dt
         cs%leafc_xfer_patch(p)      = cs%leafc_xfer_patch(p)  - cf%leafc_xfer_to_leafc_patch(p)*dt
         cs%frootc_patch(p)          = cs%frootc_patch(p)      + cf%frootc_xfer_to_frootc_patch(p)*dt
         cs%frootc_xfer_patch(p)     = cs%frootc_xfer_patch(p) - cf%frootc_xfer_to_frootc_patch(p)*dt
             if (woody(ivt(p)) == 1._r8) then
                cs%livestemc_patch(p)       = cs%livestemc_patch(p)       + cf%livestemc_xfer_to_livestemc_patch(p)*dt
                cs%livestemc_xfer_patch(p)  = cs%livestemc_xfer_patch(p)  - cf%livestemc_xfer_to_livestemc_patch(p)*dt
                cs%deadstemc_patch(p)       = cs%deadstemc_patch(p)       + cf%deadstemc_xfer_to_deadstemc_patch(p)*dt
                cs%deadstemc_xfer_patch(p)  = cs%deadstemc_xfer_patch(p)  - cf%deadstemc_xfer_to_deadstemc_patch(p)*dt
                cs%livecrootc_patch(p)      = cs%livecrootc_patch(p)      + cf%livecrootc_xfer_to_livecrootc_patch(p)*dt
                cs%livecrootc_xfer_patch(p) = cs%livecrootc_xfer_patch(p) - cf%livecrootc_xfer_to_livecrootc_patch(p)*dt
                cs%deadcrootc_patch(p)      = cs%deadcrootc_patch(p)      + cf%deadcrootc_xfer_to_deadcrootc_patch(p)*dt
                cs%deadcrootc_xfer_patch(p) = cs%deadcrootc_xfer_patch(p) - cf%deadcrootc_xfer_to_deadcrootc_patch(p)*dt
         end if
         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            ! lines here for consistency; the transfer terms are zero
            cs%livestemc_patch(p)       = cs%livestemc_patch(p)      + cf%livestemc_xfer_to_livestemc_patch(p)*dt
            cs%livestemc_xfer_patch(p)  = cs%livestemc_xfer_patch(p) - cf%livestemc_xfer_to_livestemc_patch(p)*dt
            cs%grainc_patch(p)          = cs%grainc_patch(p)         + cf%grainc_xfer_to_grainc_patch(p)*dt
            cs%grainc_xfer_patch(p)     = cs%grainc_xfer_patch(p)    - cf%grainc_xfer_to_grainc_patch(p)*dt
         end if

         ! phenology: litterfall fluxes
         cs%leafc_patch(p) = cs%leafc_patch(p) - cf%leafc_to_litter_patch(p)*dt
         cs%frootc_patch(p) = cs%frootc_patch(p) - cf%frootc_to_litter_patch(p)*dt

         ! livewood turnover fluxes
         if (woody(ivt(p)) == 1._r8) then
            cs%livestemc_patch(p)  = cs%livestemc_patch(p)  - cf%livestemc_to_deadstemc_patch(p)*dt
            cs%deadstemc_patch(p)  = cs%deadstemc_patch(p)  + cf%livestemc_to_deadstemc_patch(p)*dt
            cs%livecrootc_patch(p) = cs%livecrootc_patch(p) - cf%livecrootc_to_deadcrootc_patch(p)*dt
            cs%deadcrootc_patch(p) = cs%deadcrootc_patch(p) + cf%livecrootc_to_deadcrootc_patch(p)*dt
         end if
         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            cs%livestemc_patch(p)  = cs%livestemc_patch(p)  - cf%livestemc_to_litter_patch(p)*dt
            cs%grainc_patch(p)     = cs%grainc_patch(p)     - cf%grainc_to_food_patch(p)*dt

            cs%cropseedc_deficit_patch(p) = cs%cropseedc_deficit_patch(p) &
                 - cf%crop_seedc_to_leaf_patch(p) * dt
         end if

         ! maintenance respiration fluxes from cpool
         cs%cpool_patch(p) = cs%cpool_patch(p) - cf%cpool_to_xsmrpool_patch(p)*dt
         cs%cpool_patch(p) = cs%cpool_patch(p) - cf%leaf_curmr_patch(p)*dt
         cs%cpool_patch(p) = cs%cpool_patch(p) - cf%froot_curmr_patch(p)*dt
         if (woody(ivt(p)) == 1._r8) then
            cs%cpool_patch(p) = cs%cpool_patch(p) - cf%livestem_curmr_patch(p)*dt
            cs%cpool_patch(p) = cs%cpool_patch(p) - cf%livecroot_curmr_patch(p)*dt
         end if
         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            cs%cpool_patch(p) = cs%cpool_patch(p) - cf%livestem_curmr_patch(p)*dt
            cs%cpool_patch(p) = cs%cpool_patch(p) - cf%grain_curmr_patch(p)*dt
         end if
         ! excess respiration flux from cpool
         cs%cpool_patch(p) = cs%cpool_patch(p) - cf%xr_patch(p)*dt

         ! maintenance respiration fluxes from xsmrpool
         cs%xsmrpool_patch(p) = cs%xsmrpool_patch(p) + cf%cpool_to_xsmrpool_patch(p)*dt
         cs%xsmrpool_patch(p) = cs%xsmrpool_patch(p) - cf%leaf_xsmr_patch(p)*dt
         cs%xsmrpool_patch(p) = cs%xsmrpool_patch(p) - cf%froot_xsmr_patch(p)*dt
         if (nu_com .ne. 'RD') then
            cs%xsmrpool_patch(p) = cs%xsmrpool_patch(p) - cf%xsmrpool_turnover_patch(p)*dt
         end if
         if (woody(ivt(p)) == 1._r8) then
            cs%xsmrpool_patch(p) = cs%xsmrpool_patch(p) - cf%livestem_xsmr_patch(p)*dt
            cs%xsmrpool_patch(p) = cs%xsmrpool_patch(p) - cf%livecroot_xsmr_patch(p)*dt
         end if
         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            cs%xsmrpool_patch(p) = cs%xsmrpool_patch(p) - cf%livestem_xsmr_patch(p)*dt
            cs%xsmrpool_patch(p) = cs%xsmrpool_patch(p) - cf%grain_xsmr_patch(p)*dt
            if (harvdate(p) < 999) then ! beginning at harvest, send to atm
               cf%xsmrpool_to_atm_patch(p) = cf%xsmrpool_to_atm_patch(p) + cs%xsmrpool_patch(p)/dt
               cs%xsmrpool_patch(p)        = cs%xsmrpool_patch(p)        - cf%xsmrpool_to_atm_patch(p)*dt
            end if
         end if

         ! allocation fluxes
         cs%cpool_patch(p)           = cs%cpool_patch(p)          - cf%cpool_to_leafc_patch(p)*dt
         cs%leafc_patch(p)           = cs%leafc_patch(p)          + cf%cpool_to_leafc_patch(p)*dt
         cs%cpool_patch(p)           = cs%cpool_patch(p)          - cf%cpool_to_leafc_storage_patch(p)*dt
         cs%leafc_storage_patch(p)   = cs%leafc_storage_patch(p)  + cf%cpool_to_leafc_storage_patch(p)*dt
         cs%cpool_patch(p)           = cs%cpool_patch(p)          - cf%cpool_to_frootc_patch(p)*dt
         cs%frootc_patch(p)          = cs%frootc_patch(p)         + cf%cpool_to_frootc_patch(p)*dt
         cs%cpool_patch(p)           = cs%cpool_patch(p)          - cf%cpool_to_frootc_storage_patch(p)*dt
         cs%frootc_storage_patch(p)  = cs%frootc_storage_patch(p) + cf%cpool_to_frootc_storage_patch(p)*dt
         if (woody(ivt(p)) == 1._r8) then
            cs%cpool_patch(p)               = cs%cpool_patch(p)              - cf%cpool_to_livestemc_patch(p)*dt
            cs%livestemc_patch(p)           = cs%livestemc_patch(p)          + cf%cpool_to_livestemc_patch(p)*dt
            cs%cpool_patch(p)               = cs%cpool_patch(p)              - cf%cpool_to_livestemc_storage_patch(p)*dt
            cs%livestemc_storage_patch(p)   = cs%livestemc_storage_patch(p)  + cf%cpool_to_livestemc_storage_patch(p)*dt
            cs%cpool_patch(p)               = cs%cpool_patch(p)              - cf%cpool_to_deadstemc_patch(p)*dt
            cs%deadstemc_patch(p)           = cs%deadstemc_patch(p)          + cf%cpool_to_deadstemc_patch(p)*dt
            cs%cpool_patch(p)               = cs%cpool_patch(p)              - cf%cpool_to_deadstemc_storage_patch(p)*dt
            cs%deadstemc_storage_patch(p)   = cs%deadstemc_storage_patch(p)  + cf%cpool_to_deadstemc_storage_patch(p)*dt
            cs%cpool_patch(p)               = cs%cpool_patch(p)              - cf%cpool_to_livecrootc_patch(p)*dt
            cs%livecrootc_patch(p)          = cs%livecrootc_patch(p)         + cf%cpool_to_livecrootc_patch(p)*dt
            cs%cpool_patch(p)               = cs%cpool_patch(p)              - cf%cpool_to_livecrootc_storage_patch(p)*dt
            cs%livecrootc_storage_patch(p)  = cs%livecrootc_storage_patch(p) + cf%cpool_to_livecrootc_storage_patch(p)*dt
            cs%cpool_patch(p)               = cs%cpool_patch(p)              - cf%cpool_to_deadcrootc_patch(p)*dt
            cs%deadcrootc_patch(p)          = cs%deadcrootc_patch(p)         + cf%cpool_to_deadcrootc_patch(p)*dt
            cs%cpool_patch(p)               = cs%cpool_patch(p)              - cf%cpool_to_deadcrootc_storage_patch(p)*dt
            cs%deadcrootc_storage_patch(p)  = cs%deadcrootc_storage_patch(p) + cf%cpool_to_deadcrootc_storage_patch(p)*dt
         end if
         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            cs%cpool_patch(p)               = cs%cpool_patch(p)              - cf%cpool_to_livestemc_patch(p)*dt
            cs%livestemc_patch(p)           = cs%livestemc_patch(p)          + cf%cpool_to_livestemc_patch(p)*dt
            cs%cpool_patch(p)               = cs%cpool_patch(p)              - cf%cpool_to_livestemc_storage_patch(p)*dt
            cs%livestemc_storage_patch(p)   = cs%livestemc_storage_patch(p)  + cf%cpool_to_livestemc_storage_patch(p)*dt
            cs%cpool_patch(p)               = cs%cpool_patch(p)              - cf%cpool_to_grainc_patch(p)*dt
            cs%grainc_patch(p)              = cs%grainc_patch(p)             + cf%cpool_to_grainc_patch(p)*dt
            cs%cpool_patch(p)               = cs%cpool_patch(p)              - cf%cpool_to_grainc_storage_patch(p)*dt
            cs%grainc_storage_patch(p)      = cs%grainc_storage_patch(p)     + cf%cpool_to_grainc_storage_patch(p)*dt
         end if

         ! growth respiration fluxes for current growth
         cs%cpool_patch(p) = cs%cpool_patch(p) - cf%cpool_leaf_gr_patch(p)*dt
         cs%cpool_patch(p) = cs%cpool_patch(p) - cf%cpool_froot_gr_patch(p)*dt
         if (woody(ivt(p)) == 1._r8) then
            cs%cpool_patch(p) = cs%cpool_patch(p) - cf%cpool_livestem_gr_patch(p)*dt
            cs%cpool_patch(p) = cs%cpool_patch(p) - cf%cpool_deadstem_gr_patch(p)*dt
            cs%cpool_patch(p) = cs%cpool_patch(p) - cf%cpool_livecroot_gr_patch(p)*dt
            cs%cpool_patch(p) = cs%cpool_patch(p) - cf%cpool_deadcroot_gr_patch(p)*dt
         end if
         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            cs%cpool_patch(p) = cs%cpool_patch(p) - cf%cpool_livestem_gr_patch(p)*dt
            cs%cpool_patch(p) = cs%cpool_patch(p) - cf%cpool_grain_gr_patch(p)*dt
         end if

         ! growth respiration for transfer growth
         cs%gresp_xfer_patch(p) = cs%gresp_xfer_patch(p) - cf%transfer_leaf_gr_patch(p)*dt
         cs%gresp_xfer_patch(p) = cs%gresp_xfer_patch(p) - cf%transfer_froot_gr_patch(p)*dt
         if (woody(ivt(p)) == 1._r8) then
            cs%gresp_xfer_patch(p) = cs%gresp_xfer_patch(p) - cf%transfer_livestem_gr_patch(p)*dt
            cs%gresp_xfer_patch(p) = cs%gresp_xfer_patch(p) - cf%transfer_deadstem_gr_patch(p)*dt
            cs%gresp_xfer_patch(p) = cs%gresp_xfer_patch(p) - cf%transfer_livecroot_gr_patch(p)*dt
            cs%gresp_xfer_patch(p) = cs%gresp_xfer_patch(p) - cf%transfer_deadcroot_gr_patch(p)*dt
         end if
         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            cs%gresp_xfer_patch(p) = cs%gresp_xfer_patch(p) - cf%transfer_livestem_gr_patch(p)*dt
            cs%gresp_xfer_patch(p) = cs%gresp_xfer_patch(p) - cf%transfer_grain_gr_patch(p)*dt
         end if

         ! growth respiration at time of storage
         cs%cpool_patch(p) = cs%cpool_patch(p) - cf%cpool_leaf_storage_gr_patch(p)*dt
         cs%cpool_patch(p) = cs%cpool_patch(p) - cf%cpool_froot_storage_gr_patch(p)*dt
         if (woody(ivt(p)) == 1._r8) then
            cs%cpool_patch(p) = cs%cpool_patch(p) - cf%cpool_livestem_storage_gr_patch(p)*dt
            cs%cpool_patch(p) = cs%cpool_patch(p) - cf%cpool_deadstem_storage_gr_patch(p)*dt
            cs%cpool_patch(p) = cs%cpool_patch(p) - cf%cpool_livecroot_storage_gr_patch(p)*dt
            cs%cpool_patch(p) = cs%cpool_patch(p) - cf%cpool_deadcroot_storage_gr_patch(p)*dt
         end if
         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            cs%cpool_patch(p) = cs%cpool_patch(p) - cf%cpool_livestem_storage_gr_patch(p)*dt
            cs%cpool_patch(p) = cs%cpool_patch(p) - cf%cpool_grain_storage_gr_patch(p)*dt
         end if

         ! growth respiration stored for release during transfer growth
         cs%cpool_patch(p)         = cs%cpool_patch(p)         - cf%cpool_to_gresp_storage_patch(p)*dt
         cs%gresp_storage_patch(p) = cs%gresp_storage_patch(p) + cf%cpool_to_gresp_storage_patch(p)*dt

         ! move storage pools into transfer pools
         cs%leafc_storage_patch(p)  = cs%leafc_storage_patch(p)  - cf%leafc_storage_to_xfer_patch(p)*dt
         cs%leafc_xfer_patch(p)     = cs%leafc_xfer_patch(p)     + cf%leafc_storage_to_xfer_patch(p)*dt
         cs%frootc_storage_patch(p) = cs%frootc_storage_patch(p) - cf%frootc_storage_to_xfer_patch(p)*dt
         cs%frootc_xfer_patch(p)    = cs%frootc_xfer_patch(p)    + cf%frootc_storage_to_xfer_patch(p)*dt
         if (woody(ivt(p)) == 1._r8) then
            cs%livestemc_storage_patch(p)  = cs%livestemc_storage_patch(p) - cf%livestemc_storage_to_xfer_patch(p)*dt
            cs%livestemc_xfer_patch(p)     = cs%livestemc_xfer_patch(p)    + cf%livestemc_storage_to_xfer_patch(p)*dt
            cs%deadstemc_storage_patch(p)  = cs%deadstemc_storage_patch(p) - cf%deadstemc_storage_to_xfer_patch(p)*dt
            cs%deadstemc_xfer_patch(p)     = cs%deadstemc_xfer_patch(p)    + cf%deadstemc_storage_to_xfer_patch(p)*dt
            cs%livecrootc_storage_patch(p) = cs%livecrootc_storage_patch(p)- cf%livecrootc_storage_to_xfer_patch(p)*dt
            cs%livecrootc_xfer_patch(p)    = cs%livecrootc_xfer_patch(p)   + cf%livecrootc_storage_to_xfer_patch(p)*dt
            cs%deadcrootc_storage_patch(p) = cs%deadcrootc_storage_patch(p)- cf%deadcrootc_storage_to_xfer_patch(p)*dt
            cs%deadcrootc_xfer_patch(p)    = cs%deadcrootc_xfer_patch(p)   + cf%deadcrootc_storage_to_xfer_patch(p)*dt
            cs%gresp_storage_patch(p)      = cs%gresp_storage_patch(p)     - cf%gresp_storage_to_xfer_patch(p)*dt
            cs%gresp_xfer_patch(p)         = cs%gresp_xfer_patch(p)        + cf%gresp_storage_to_xfer_patch(p)*dt
         end if
         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            ! lines here for consistency; the transfer terms are zero
            cs%livestemc_storage_patch(p)  = cs%livestemc_storage_patch(p) - cf%livestemc_storage_to_xfer_patch(p)*dt
            cs%livestemc_xfer_patch(p)     = cs%livestemc_xfer_patch(p)    + cf%livestemc_storage_to_xfer_patch(p)*dt
            cs%grainc_storage_patch(p)     = cs%grainc_storage_patch(p)    - cf%grainc_storage_to_xfer_patch(p)*dt
            cs%grainc_xfer_patch(p)        = cs%grainc_xfer_patch(p)       + cf%grainc_storage_to_xfer_patch(p)*dt
         end if

      end do ! end of patch loop

   end if

 end associate

end subroutine CarbonStateUpdate1

end module CarbonStateUpdate1Mod
