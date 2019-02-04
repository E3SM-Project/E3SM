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
  use ColumnDataType         , only : column_carbon_state, column_carbon_flux
  use VegetationDataType     , only : vegetation_carbon_state, vegetation_carbon_flux
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
  subroutine CStateUpdateDynPatch(bounds, num_soilc_with_inactive, filter_soilc_with_inactive, &
       carbonflux_vars, col_cs, col_cfv2)
    !
    ! !DESCRIPTION:
    ! Update carbon states based on fluxes from dyn_cnbal_patch
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: num_soilc_with_inactive       ! number of columns in soil filter
    integer                , intent(in)    :: filter_soilc_with_inactive(:) ! soil column filter that includes inactive points
    type(carbonflux_type)  , intent(in)    :: carbonflux_vars
    type(column_carbon_state),intent(inout):: col_cs
    type(column_carbon_flux),intent(inout) :: col_cfv2
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
         cs => col_cs           , &
         ccfv2 => col_cfv2
         )

      dt = real( get_step_size(), r8 )

      if (.not.use_fates) then

         do g = bounds%begg, bounds%endg
            cs%seedc(g) = cs%seedc(g) &
                 - cf%dwt_seedc_to_leaf_grc(g)     * dt &
                 - cf%dwt_seedc_to_deadstem_grc(g) * dt
         end do

         do j = 1,nlevdecomp
            do fc = 1, num_soilc_with_inactive
               c = filter_soilc_with_inactive(fc)

               cs%decomp_cpools_vr(c,j,i_met_lit) = cs%decomp_cpools_vr(c,j,i_met_lit) + &
                    cf%dwt_frootc_to_litr_met_c_col(c,j) * dt
               cs%decomp_cpools_vr(c,j,i_cel_lit) = cs%decomp_cpools_vr(c,j,i_cel_lit) + &
                    cf%dwt_frootc_to_litr_cel_c_col(c,j) * dt
               cs%decomp_cpools_vr(c,j,i_lig_lit) = cs%decomp_cpools_vr(c,j,i_lig_lit) + &
                    cf%dwt_frootc_to_litr_lig_c_col(c,j) * dt
               cs%decomp_cpools_vr(c,j,i_cwd) = cs%decomp_cpools_vr(c,j,i_cwd) + &
                    ( cf%dwt_livecrootc_to_cwdc_col(c,j) + cf%dwt_deadcrootc_to_cwdc_col(c,j) ) * dt

            end do
         end do
      end if

    end associate

  end subroutine CarbonStateUpdateDynPatch

  !-----------------------------------------------------------------------
  subroutine CarbonStateUpdate0(&
       num_soilp, filter_soilp, &
       carbonflux_vars, carbonstate_vars, col_csv2, veg_csv2, col_cfv2, veg_cfv2)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update cpool carbon state
    !

    ! !ARGUMENTS:
    integer                , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(carbonflux_type)  , intent(in)    :: carbonflux_vars
    type(carbonstate_type) , intent(inout) :: carbonstate_vars
    type(column_carbon_state)    ,intent(inout) :: col_csv2
    type(vegetation_carbon_state),intent(inout) :: veg_csv2
    type(column_carbon_flux)     ,intent(inout) :: col_cfv2
    type(vegetation_carbon_flux) ,intent(inout) :: veg_cfv2
    !
    ! !LOCAL VARIABLES:
    integer :: p  ! indices
    integer :: fp ! lake filter indices
    real(r8):: dt ! radiation time step (seconds)
    !-----------------------------------------------------------------------

    associate(                                                          & 
         cf                => carbonflux_vars                         , &
         cs                => carbonstate_vars                        , &
         csv2              => col_csv2                                , &
         vcsv2             => veg_csv2                                , &
         ccfv2             => col_cfv2                                , &
         vcfv2             => veg_cfv2                                 &
         )

      ! set time steps
      dt = real( get_step_size(), r8 )

      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)
         ! gross photosynthesis fluxes
         vcsv2%cpool(p) = vcsv2%cpool(p) + vcfv2%psnsun_to_cpool(p)*dt
         vcsv2%cpool(p) = vcsv2%cpool(p) + vcfv2%psnshade_to_cpool(p)*dt
      end do

    end associate

  end subroutine CarbonStateUpdate0

  !-----------------------------------------------------------------------
  subroutine CarbonStateUpdate1(bounds, &
       num_soilc, filter_soilc, &
       num_soilp, filter_soilp, &
       crop_vars, carbonflux_vars, carbonstate_vars, col_csv2, veg_csv2, col_cfv2, veg_cfv2)
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
    type(column_carbon_state),intent(inout):: col_csv2
    type(vegetation_carbon_state),intent(inout) :: veg_csv2
    type(column_carbon_flux)     ,intent(inout) :: col_cfv2
    type(vegetation_carbon_flux) ,intent(inout) :: veg_cfv2
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
         cs => carbonstate_vars , &
         csv2 => col_csv2       , &
         vcsv2             => veg_csv2                                , &
         ccfv2             => col_cfv2                                , &
         vcfv2             => veg_cfv2                                 &

         )

      ! set time steps
      dt = real( get_step_size(), r8 )

      ! column level fluxes

      if(.not.use_fates) then
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            csv2%decomp_som2c_vr(c,1:nlevdecomp) = csv2%decomp_cpools_vr(c,1:nlevdecomp,6)
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
         vcsv2%leafc(p)           = vcsv2%leafc(p)       + vcfv2%leafc_xfer_to_leafc(p)*dt
         vcsv2%leafc_xfer(p)      = vcsv2%leafc_xfer(p)  - vcfv2%leafc_xfer_to_leafc(p)*dt
         vcsv2%frootc(p)          = vcsv2%frootc(p)      + vcfv2%frootc_xfer_to_frootc(p)*dt
         vcsv2%frootc_xfer(p)     = vcsv2%frootc_xfer(p) - vcfv2%frootc_xfer_to_frootc(p)*dt
             if (woody(ivt(p)) == 1._r8) then
                vcsv2%livestemc(p)       = vcsv2%livestemc(p)       + vcfv2%livestemc_xfer_to_livestemc(p)*dt
                vcsv2%livestemc_xfer(p)  = vcsv2%livestemc_xfer(p)  - vcfv2%livestemc_xfer_to_livestemc(p)*dt
                vcsv2%deadstemc(p)       = vcsv2%deadstemc(p)       + vcfv2%deadstemc_xfer_to_deadstemc(p)*dt
                vcsv2%deadstemc_xfer(p)  = vcsv2%deadstemc_xfer(p)  - vcfv2%deadstemc_xfer_to_deadstemc(p)*dt
                vcsv2%livecrootc(p)      = vcsv2%livecrootc(p)      + vcfv2%livecrootc_xfer_to_livecrootc(p)*dt
                vcsv2%livecrootc_xfer(p) = vcsv2%livecrootc_xfer(p) - vcfv2%livecrootc_xfer_to_livecrootc(p)*dt
                vcsv2%deadcrootc(p)      = vcsv2%deadcrootc(p)      + vcfv2%deadcrootc_xfer_to_deadcrootc(p)*dt
                vcsv2%deadcrootc_xfer(p) = vcsv2%deadcrootc_xfer(p) - vcfv2%deadcrootc_xfer_to_deadcrootc(p)*dt
         end if
         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            ! lines here for consistency; the transfer terms are zero
            vcsv2%livestemc(p)       = vcsv2%livestemc(p)      + vcfv2%livestemc_xfer_to_livestemc(p)*dt
            vcsv2%livestemc_xfer(p)  = vcsv2%livestemc_xfer(p) - vcfv2%livestemc_xfer_to_livestemc(p)*dt
            vcsv2%grainc(p)          = vcsv2%grainc(p)         + vcfv2%grainc_xfer_to_grainc(p)*dt
            vcsv2%grainc_xfer(p)     = vcsv2%grainc_xfer(p)    - vcfv2%grainc_xfer_to_grainc(p)*dt
         end if

         ! phenology: litterfall fluxes
         vcsv2%leafc(p) = vcsv2%leafc(p) - vcfv2%leafc_to_litter(p)*dt
         vcsv2%frootc(p) = vcsv2%frootc(p) - vcfv2%frootc_to_litter(p)*dt

         ! livewood turnover fluxes
         if (woody(ivt(p)) == 1._r8) then
            vcsv2%livestemc(p)  = vcsv2%livestemc(p)  - vcfv2%livestemc_to_deadstemc(p)*dt
            vcsv2%deadstemc(p)  = vcsv2%deadstemc(p)  + vcfv2%livestemc_to_deadstemc(p)*dt
            vcsv2%livecrootc(p) = vcsv2%livecrootc(p) - vcfv2%livecrootc_to_deadcrootc(p)*dt
            vcsv2%deadcrootc(p) = vcsv2%deadcrootc(p) + vcfv2%livecrootc_to_deadcrootc(p)*dt
         end if
         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            vcsv2%livestemc(p)  = vcsv2%livestemc(p)  - vcfv2%livestemc_to_litter(p)*dt
            vcsv2%grainc(p)     = vcsv2%grainc(p)     - vcfv2%grainc_to_food(p)*dt

            vcsv2%cropseedc_deficit(p) = vcsv2%cropseedc_deficit(p) &
                 - vcfv2%crop_seedc_to_leaf(p) * dt
         end if

         ! maintenance respiration fluxes from cpool
         vcsv2%cpool(p) = vcsv2%cpool(p) - vcfv2%cpool_to_xsmrpool(p)*dt
         vcsv2%cpool(p) = vcsv2%cpool(p) - vcfv2%leaf_curmr(p)*dt
         vcsv2%cpool(p) = vcsv2%cpool(p) - vcfv2%froot_curmr(p)*dt
         if (woody(ivt(p)) == 1._r8) then
            vcsv2%cpool(p) = vcsv2%cpool(p) - vcfv2%livestem_curmr(p)*dt
            vcsv2%cpool(p) = vcsv2%cpool(p) - vcfv2%livecroot_curmr(p)*dt
         end if
         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            vcsv2%cpool(p) = vcsv2%cpool(p) - vcfv2%livestem_curmr(p)*dt
            vcsv2%cpool(p) = vcsv2%cpool(p) - vcfv2%grain_curmr(p)*dt
         end if
         ! excess respiration flux from cpool
         vcsv2%cpool(p) = vcsv2%cpool(p) - vcfv2%xr(p)*dt

         ! maintenance respiration fluxes from xsmrpool
         vcsv2%xsmrpool(p) = vcsv2%xsmrpool(p) + vcfv2%cpool_to_xsmrpool(p)*dt
         vcsv2%xsmrpool(p) = vcsv2%xsmrpool(p) - vcfv2%leaf_xsmr(p)*dt
         vcsv2%xsmrpool(p) = vcsv2%xsmrpool(p) - vcfv2%froot_xsmr(p)*dt
         if (nu_com .ne. 'RD') then
            vcsv2%xsmrpool(p) = vcsv2%xsmrpool(p) - vcfv2%xsmrpool_turnover(p)*dt
         end if
         if (woody(ivt(p)) == 1._r8) then
            vcsv2%xsmrpool(p) = vcsv2%xsmrpool(p) - vcfv2%livestem_xsmr(p)*dt
            vcsv2%xsmrpool(p) = vcsv2%xsmrpool(p) - vcfv2%livecroot_xsmr(p)*dt
         end if
         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            vcsv2%xsmrpool(p) = vcsv2%xsmrpool(p) - vcfv2%livestem_xsmr(p)*dt
            vcsv2%xsmrpool(p) = vcsv2%xsmrpool(p) - vcfv2%grain_xsmr(p)*dt
            if (harvdate(p) < 999) then ! beginning at harvest, send to atm
               vcfv2%xsmrpool_to_atm(p) = vcfv2%xsmrpool_to_atm(p) + vcsv2%xsmrpool(p)/dt
               vcsv2%xsmrpool(p)        = vcsv2%xsmrpool(p)        - vcfv2%xsmrpool_to_atm(p)*dt
            end if
         end if

         ! allocation fluxes
         vcsv2%cpool(p)           = vcsv2%cpool(p)          - vcfv2%cpool_to_leafc(p)*dt
         vcsv2%leafc(p)           = vcsv2%leafc(p)          + vcfv2%cpool_to_leafc(p)*dt
         vcsv2%cpool(p)           = vcsv2%cpool(p)          - vcfv2%cpool_to_leafc_storage(p)*dt
         vcsv2%leafc_storage(p)   = vcsv2%leafc_storage(p)  + vcfv2%cpool_to_leafc_storage(p)*dt
         vcsv2%cpool(p)           = vcsv2%cpool(p)          - vcfv2%cpool_to_frootc(p)*dt
         vcsv2%frootc(p)          = vcsv2%frootc(p)         + vcfv2%cpool_to_frootc(p)*dt
         vcsv2%cpool(p)           = vcsv2%cpool(p)          - vcfv2%cpool_to_frootc_storage(p)*dt
         vcsv2%frootc_storage(p)  = vcsv2%frootc_storage(p) + vcfv2%cpool_to_frootc_storage(p)*dt
         if (woody(ivt(p)) == 1._r8) then
            vcsv2%cpool(p)               = vcsv2%cpool(p)              - vcfv2%cpool_to_livestemc(p)*dt
            vcsv2%livestemc(p)           = vcsv2%livestemc(p)          + vcfv2%cpool_to_livestemc(p)*dt
            vcsv2%cpool(p)               = vcsv2%cpool(p)              - vcfv2%cpool_to_livestemc_storage(p)*dt
            vcsv2%livestemc_storage(p)   = vcsv2%livestemc_storage(p)  + vcfv2%cpool_to_livestemc_storage(p)*dt
            vcsv2%cpool(p)               = vcsv2%cpool(p)              - vcfv2%cpool_to_deadstemc(p)*dt
            vcsv2%deadstemc(p)           = vcsv2%deadstemc(p)          + vcfv2%cpool_to_deadstemc(p)*dt
            vcsv2%cpool(p)               = vcsv2%cpool(p)              - vcfv2%cpool_to_deadstemc_storage(p)*dt
            vcsv2%deadstemc_storage(p)   = vcsv2%deadstemc_storage(p)  + vcfv2%cpool_to_deadstemc_storage(p)*dt
            vcsv2%cpool(p)               = vcsv2%cpool(p)              - vcfv2%cpool_to_livecrootc(p)*dt
            vcsv2%livecrootc(p)          = vcsv2%livecrootc(p)         + vcfv2%cpool_to_livecrootc(p)*dt
            vcsv2%cpool(p)               = vcsv2%cpool(p)              - vcfv2%cpool_to_livecrootc_storage(p)*dt
            vcsv2%livecrootc_storage(p)  = vcsv2%livecrootc_storage(p) + vcfv2%cpool_to_livecrootc_storage(p)*dt
            vcsv2%cpool(p)               = vcsv2%cpool(p)              - vcfv2%cpool_to_deadcrootc(p)*dt
            vcsv2%deadcrootc(p)          = vcsv2%deadcrootc(p)         + vcfv2%cpool_to_deadcrootc(p)*dt
            vcsv2%cpool(p)               = vcsv2%cpool(p)              - vcfv2%cpool_to_deadcrootc_storage(p)*dt
            vcsv2%deadcrootc_storage(p)  = vcsv2%deadcrootc_storage(p) + vcfv2%cpool_to_deadcrootc_storage(p)*dt
         end if
         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            vcsv2%cpool(p)               = vcsv2%cpool(p)              - vcfv2%cpool_to_livestemc(p)*dt
            vcsv2%livestemc(p)           = vcsv2%livestemc(p)          + vcfv2%cpool_to_livestemc(p)*dt
            vcsv2%cpool(p)               = vcsv2%cpool(p)              - vcfv2%cpool_to_livestemc_storage(p)*dt
            vcsv2%livestemc_storage(p)   = vcsv2%livestemc_storage(p)  + vcfv2%cpool_to_livestemc_storage(p)*dt
            vcsv2%cpool(p)               = vcsv2%cpool(p)              - vcfv2%cpool_to_grainc(p)*dt
            vcsv2%grainc(p)              = vcsv2%grainc(p)             + vcfv2%cpool_to_grainc(p)*dt
            vcsv2%cpool(p)               = vcsv2%cpool(p)              - vcfv2%cpool_to_grainc_storage(p)*dt
            vcsv2%grainc_storage(p)      = vcsv2%grainc_storage(p)     + vcfv2%cpool_to_grainc_storage(p)*dt
         end if

         ! growth respiration fluxes for current growth
         vcsv2%cpool(p) = vcsv2%cpool(p) - vcfv2%cpool_leaf_gr(p)*dt
         vcsv2%cpool(p) = vcsv2%cpool(p) - vcfv2%cpool_froot_gr(p)*dt
         if (woody(ivt(p)) == 1._r8) then
            vcsv2%cpool(p) = vcsv2%cpool(p) - vcfv2%cpool_livestem_gr(p)*dt
            vcsv2%cpool(p) = vcsv2%cpool(p) - vcfv2%cpool_deadstem_gr(p)*dt
            vcsv2%cpool(p) = vcsv2%cpool(p) - vcfv2%cpool_livecroot_gr(p)*dt
            vcsv2%cpool(p) = vcsv2%cpool(p) - vcfv2%cpool_deadcroot_gr(p)*dt
         end if
         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            vcsv2%cpool(p) = vcsv2%cpool(p) - vcfv2%cpool_livestem_gr(p)*dt
            vcsv2%cpool(p) = vcsv2%cpool(p) - vcfv2%cpool_grain_gr(p)*dt
         end if

         ! growth respiration for transfer growth
         vcsv2%gresp_xfer(p) = vcsv2%gresp_xfer(p) - vcfv2%transfer_leaf_gr(p)*dt
         vcsv2%gresp_xfer(p) = vcsv2%gresp_xfer(p) - vcfv2%transfer_froot_gr(p)*dt
         if (woody(ivt(p)) == 1._r8) then
            vcsv2%gresp_xfer(p) = vcsv2%gresp_xfer(p) - vcfv2%transfer_livestem_gr(p)*dt
            vcsv2%gresp_xfer(p) = vcsv2%gresp_xfer(p) - vcfv2%transfer_deadstem_gr(p)*dt
            vcsv2%gresp_xfer(p) = vcsv2%gresp_xfer(p) - vcfv2%transfer_livecroot_gr(p)*dt
            vcsv2%gresp_xfer(p) = vcsv2%gresp_xfer(p) - vcfv2%transfer_deadcroot_gr(p)*dt
         end if
         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            vcsv2%gresp_xfer(p) = vcsv2%gresp_xfer(p) - vcfv2%transfer_livestem_gr(p)*dt
            vcsv2%gresp_xfer(p) = vcsv2%gresp_xfer(p) - vcfv2%transfer_grain_gr(p)*dt
         end if

         ! growth respiration at time of storage
         vcsv2%cpool(p) = vcsv2%cpool(p) - vcfv2%cpool_leaf_storage_gr(p)*dt
         vcsv2%cpool(p) = vcsv2%cpool(p) - vcfv2%cpool_froot_storage_gr(p)*dt
         if (woody(ivt(p)) == 1._r8) then
            vcsv2%cpool(p) = vcsv2%cpool(p) - vcfv2%cpool_livestem_storage_gr(p)*dt
            vcsv2%cpool(p) = vcsv2%cpool(p) - vcfv2%cpool_deadstem_storage_gr(p)*dt
            vcsv2%cpool(p) = vcsv2%cpool(p) - vcfv2%cpool_livecroot_storage_gr(p)*dt
            vcsv2%cpool(p) = vcsv2%cpool(p) - vcfv2%cpool_deadcroot_storage_gr(p)*dt
         end if
         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            vcsv2%cpool(p) = vcsv2%cpool(p) - vcfv2%cpool_livestem_storage_gr(p)*dt
            vcsv2%cpool(p) = vcsv2%cpool(p) - vcfv2%cpool_grain_storage_gr(p)*dt
         end if

         ! growth respiration stored for release during transfer growth
         vcsv2%cpool(p)         = vcsv2%cpool(p)         - vcfv2%cpool_to_gresp_storage(p)*dt
         vcsv2%gresp_storage(p) = vcsv2%gresp_storage(p) + vcfv2%cpool_to_gresp_storage(p)*dt

         ! move storage pools into transfer pools
         vcsv2%leafc_storage(p)  = vcsv2%leafc_storage(p)  - vcfv2%leafc_storage_to_xfer(p)*dt
         vcsv2%leafc_xfer(p)     = vcsv2%leafc_xfer(p)     + vcfv2%leafc_storage_to_xfer(p)*dt
         vcsv2%frootc_storage(p) = vcsv2%frootc_storage(p) - vcfv2%frootc_storage_to_xfer(p)*dt
         vcsv2%frootc_xfer(p)    = vcsv2%frootc_xfer(p)    + vcfv2%frootc_storage_to_xfer(p)*dt
         if (woody(ivt(p)) == 1._r8) then
            vcsv2%livestemc_storage(p)  = vcsv2%livestemc_storage(p) - vcfv2%livestemc_storage_to_xfer(p)*dt
            vcsv2%livestemc_xfer(p)     = vcsv2%livestemc_xfer(p)    + vcfv2%livestemc_storage_to_xfer(p)*dt
            vcsv2%deadstemc_storage(p)  = vcsv2%deadstemc_storage(p) - vcfv2%deadstemc_storage_to_xfer(p)*dt
            vcsv2%deadstemc_xfer(p)     = vcsv2%deadstemc_xfer(p)    + vcfv2%deadstemc_storage_to_xfer(p)*dt
            vcsv2%livecrootc_storage(p) = vcsv2%livecrootc_storage(p)- vcfv2%livecrootc_storage_to_xfer(p)*dt
            vcsv2%livecrootc_xfer(p)    = vcsv2%livecrootc_xfer(p)   + vcfv2%livecrootc_storage_to_xfer(p)*dt
            vcsv2%deadcrootc_storage(p) = vcsv2%deadcrootc_storage(p)- vcfv2%deadcrootc_storage_to_xfer(p)*dt
            vcsv2%deadcrootc_xfer(p)    = vcsv2%deadcrootc_xfer(p)   + vcfv2%deadcrootc_storage_to_xfer(p)*dt
            vcsv2%gresp_storage(p)      = vcsv2%gresp_storage(p)     - vcfv2%gresp_storage_to_xfer(p)*dt
            vcsv2%gresp_xfer(p)         = vcsv2%gresp_xfer(p)        + vcfv2%gresp_storage_to_xfer(p)*dt
         end if
         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            ! lines here for consistency; the transfer terms are zero
            vcsv2%livestemc_storage(p)  = vcsv2%livestemc_storage(p) - vcfv2%livestemc_storage_to_xfer(p)*dt
            vcsv2%livestemc_xfer(p)     = vcsv2%livestemc_xfer(p)    + vcfv2%livestemc_storage_to_xfer(p)*dt
            vcsv2%grainc_storage(p)     = vcsv2%grainc_storage(p)    - vcfv2%grainc_storage_to_xfer(p)*dt
            vcsv2%grainc_xfer(p)        = vcsv2%grainc_xfer(p)       + vcfv2%grainc_storage_to_xfer(p)*dt
         end if

      end do ! end of patch loop

   end if

 end associate

end subroutine CarbonStateUpdate1

end module CarbonStateUpdate1Mod
