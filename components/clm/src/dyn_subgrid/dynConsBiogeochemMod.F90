module dynConsBiogeochemMod

#include "shr_assert.h"

  !---------------------------------------------------------------------------
  !
  ! !DESCRIPTION:
  ! Handle conservation of biogeochemical quantities (C & N) with dynamic land cover.
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use decompMod           , only : bounds_type
  use abortutils          , only : endrun
  use clm_varctl          , only : iulog, use_c13, use_c14
  use VegetationPropertiesType      , only : veg_vp
  use CanopyStateType     , only : canopystate_type
  use PhotosynthesisType  , only : photosyns_type
  use CNStateType         , only : cnstate_type
  use CNCarbonFluxType    , only : carbonflux_type
  use CNCarbonStateType   , only : carbonstate_type
  use CNNitrogenFluxType  , only : nitrogenflux_type
  use CNNitrogenStateType , only : nitrogenstate_type
  use PhosphorusFluxType  , only : phosphorusflux_type
  use PhosphorusStateType , only : phosphorusstate_type
  use LandunitType        , only : lun_pp                
  use ColumnType          , only : col_pp                
  use VegetationType           , only : veg_pp                
  use CNSpeciesMod        , only : CN_SPECIES_C12, CN_SPECIES_C13, CN_SPECIES_C14
  use CNSpeciesMod        , only : CN_SPECIES_N, CN_SPECIES_P
  use clm_varcon          , only : c3_r2, c4_r2, c14ratio
  use dynPatchStateUpdaterMod      , only : patch_state_updater_type
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  integer  , parameter :: COMPONENT_LEAF       = 1
  integer  , parameter :: COMPONENT_DEADWOOD   = 2
  integer  , parameter :: COMPONENT_SEED       = 3
  real(r8) , parameter :: leafc_seed_param     = 1._r8
  real(r8) , parameter :: npool_seed_param     = 0.1_r8
  real(r8) , parameter :: ppool_seed_param     = 0.01_r8
  real(r8) , parameter :: deadstemc_seed_param = 0.1_r8

  save
  
  public :: dyn_cnbal_patch
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine dyn_cnbal_patch(bounds, &
       num_soilp_with_inactive, filter_soilp_with_inactive, &
       num_soilc_with_inactive, filter_soilc_with_inactive, &
       prior_weights, &
       patch_state_updater, &
       canopystate_vars, photosyns_vars, cnstate_vars, &
       carbonstate_vars, c13_carbonstate_vars, c14_carbonstate_vars, &
       carbonflux_vars, c13_carbonflux_vars, c14_carbonflux_vars, &
       nitrogenstate_vars, nitrogenflux_vars, &
       phosphorusstate_vars,phosphorusflux_vars)
    !
    ! !DESCRIPTION:
    ! Modify pft-level state and flux variables to maintain carbon and nitrogen balance with
    ! dynamic pft-weights.
    !
    ! !USES:
    use shr_const_mod      , only : SHR_CONST_PDB
    use landunit_varcon    , only : istsoil, istcrop
    use clm_varpar         , only : numveg, nlevdecomp, max_patch_per_col
    use pftvarcon          , only : pconv, pprod10, pprod100
    use clm_varcon         , only : c13ratio, c14ratio
    use clm_time_manager   , only : get_step_size
    use dynPriorWeightsMod , only : prior_weights_type
    !
    ! !ARGUMENTS:
    type(bounds_type)              , intent(in)    :: bounds        
    integer                        , intent(in)    :: num_soilp_with_inactive ! number of points in filter
    integer                        , intent(in)    :: filter_soilp_with_inactive(:) ! soil patch filter that includes inactive points
    integer                        , intent(in)    :: num_soilc_with_inactive ! number of points in filter
    integer                        , intent(in)    :: filter_soilc_with_inactive(:) ! soil column filter that includes inactive points
    type(prior_weights_type)       , intent(in)    :: prior_weights ! weights prior to the subgrid weight updates
    type(patch_state_updater_type) , intent(in)    :: patch_state_updater
    type(canopystate_type)         , intent(inout) :: canopystate_vars
    type(photosyns_type)           , intent(inout) :: photosyns_vars
    type(cnstate_type)             , intent(inout) :: cnstate_vars
    type(carbonstate_type)         , intent(inout) :: carbonstate_vars
    type(carbonstate_type)         , intent(inout) :: c13_carbonstate_vars
    type(carbonstate_type)         , intent(inout) :: c14_carbonstate_vars
    type(carbonflux_type)          , intent(inout) :: carbonflux_vars
    type(carbonflux_type)          , intent(inout) :: c13_carbonflux_vars
    type(carbonflux_type)          , intent(inout) :: c14_carbonflux_vars
    type(nitrogenstate_type)       , intent(inout) :: nitrogenstate_vars
    type(nitrogenflux_type)        , intent(inout) :: nitrogenflux_vars

    type(phosphorusstate_type)     , intent(inout) :: phosphorusstate_vars
    type(phosphorusflux_type)      , intent(inout) :: phosphorusflux_vars
    !
    ! !LOCAL VARIABLES:
    integer                       :: pi,p,c,l,g,j                  ! indices
    integer                       :: ier                           ! error code
    real(r8)                      :: dwt                           ! change in pft weight (relative to column)
    real(r8)                      :: dt                            ! land model time step (sec)
    real(r8)                      :: init_h2ocan                   ! initial canopy water mass
    real(r8)                      :: new_h2ocan                    ! canopy water mass after weight shift
    real(r8), allocatable         :: dwt_leafc_seed(:)             ! pft-level mass gain due to seeding of new area
    real(r8), allocatable         :: dwt_leafn_seed(:)             ! pft-level mass gain due to seeding of new area
    real(r8), allocatable         :: dwt_deadstemc_seed(:)         ! pft-level mass gain due to seeding of new area
    real(r8), allocatable         :: dwt_deadstemn_seed(:)         ! pft-level mass gain due to seeding of new area
    real(r8), allocatable         :: dwt_npool_seed(:)             ! pft-level mass gain due to seeding of new area
    real(r8), allocatable         :: dwt_frootc_to_litter(:)       ! pft-level mass loss due to weight shift
    real(r8), allocatable         :: dwt_livecrootc_to_litter(:)   ! pft-level mass loss due to weight shift
    real(r8), allocatable         :: dwt_deadcrootc_to_litter(:)   ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_frootn_to_litter(:)       ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_livecrootn_to_litter(:)   ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_deadcrootn_to_litter(:)   ! pft-level mass loss due to weight shift
    real(r8), allocatable         :: conv_cflux(:)                 ! pft-level mass loss due to weight shift
    real(r8), allocatable         :: prod10_cflux(:)               ! pft-level mass loss due to weight shift
    real(r8), allocatable         :: prod100_cflux(:)              ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: conv_nflux(:)                 ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod10_nflux(:)               ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod100_nflux(:)              ! pft-level mass loss due to weight shift
    real(r8)                      :: t1,t2,wt_new,wt_old
    real(r8)                      :: init_state, change_state, new_state
    real(r8)                      :: tot_leaf, pleaf, pstor, pxfer
    real(r8)                      :: leafc_seed, leafn_seed
    real(r8)                      :: deadstemc_seed, deadstemn_seed, npool_seed
    real(r8), pointer             :: dwt_ptr0, dwt_ptr1, dwt_ptr2, dwt_ptr3, ptr
    character(len=32)             :: subname='dyn_cbal'            ! subroutine name

    ! ! add phosphorus local variables
    real(r8), allocatable         :: dwt_leafp_seed(:)             ! pft-level mass gain due to seeding of new area    
    real(r8), allocatable         :: dwt_deadstemp_seed(:)         ! pft-level mass gain due to seeding of new area
    real(r8), allocatable         :: dwt_ppool_seed(:)             ! pft-level mass gain due to seeding of new area
    real(r8), allocatable, target :: dwt_frootp_to_litter(:)       ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_livecrootp_to_litter(:)   ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_deadcrootp_to_litter(:)   ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: conv_pflux(:)                 ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod10_pflux(:)               ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod100_pflux(:)              ! pft-level mass loss due to weight shift
    real(r8)                      :: leafp_seed
    real(r8)                      :: deadstemp_seed, ppool_seed

    !! C13
    real(r8), allocatable         :: dwt_leafc13_seed(:)           ! pft-level mass gain due to seeding of new area
    real(r8), allocatable         :: dwt_deadstemc13_seed(:)       ! pft-level mass gain due to seeding of new area
    real(r8), allocatable, target :: dwt_frootc13_to_litter(:)     ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_livecrootc13_to_litter(:) ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_deadcrootc13_to_litter(:) ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: conv_c13flux(:)               ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod10_c13flux(:)             ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod100_c13flux(:)            ! pft-level mass loss due to weight shift
    real(r8)                      :: leafc13_seed, deadstemc13_seed
    !! C14
    real(r8), allocatable         :: dwt_leafc14_seed(:)           ! pft-level mass gain due to seeding of new area
    real(r8), allocatable         :: dwt_deadstemc14_seed(:)       ! pft-level mass gain due to seeding of new area
    real(r8), allocatable, target :: dwt_frootc14_to_litter(:)     ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_livecrootc14_to_litter(:) ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_deadcrootc14_to_litter(:) ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: conv_c14flux(:)               ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod10_c14flux(:)             ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod100_c14flux(:)            ! pft-level mass loss due to weight shift
    real(r8)                      :: leafc14_seed, deadstemc14_seed
    real(r8) :: froot, croot
    real(r8) :: fr_flab, fr_fcel, fr_flig
    !-----------------------------------------------------------------------
    
    associate(&
         cs      => carbonstate_vars,    &
         c13_cs => c13_carbonstate_vars, &
         c14_cs => c14_carbonstate_vars, &
         cf     => carbonflux_vars     , &
         c13_cf => c13_carbonflux_vars , &
         c14_cf => c14_carbonflux_vars , &
         ns     => nitrogenstate_vars  , &
         nf     => nitrogenflux_vars   ,  &
         ps     => phosphorusstate_vars  , &
         pf     => phosphorusflux_vars     &
         )


    ! Allocate pft-level mass loss arrays
    allocate(dwt_leafc_seed           (bounds%begp:bounds%endp), stat=ier)
    allocate(dwt_leafn_seed           (bounds%begp:bounds%endp), stat=ier)
    allocate(dwt_deadstemc_seed       (bounds%begp:bounds%endp), stat=ier)
    allocate(dwt_deadstemn_seed       (bounds%begp:bounds%endp), stat=ier)
    allocate(dwt_npool_seed           (bounds%begp:bounds%endp), stat=ier)
    allocate(dwt_frootc_to_litter     (bounds%begp:bounds%endp), stat=ier)
    allocate(dwt_livecrootc_to_litter (bounds%begp:bounds%endp), stat=ier)
    allocate(dwt_deadcrootc_to_litter (bounds%begp:bounds%endp), stat=ier)
    allocate(dwt_frootn_to_litter     (bounds%begp:bounds%endp), stat=ier)
    allocate(dwt_livecrootn_to_litter (bounds%begp:bounds%endp), stat=ier)
    allocate(dwt_deadcrootn_to_litter (bounds%begp:bounds%endp), stat=ier)
    allocate(conv_cflux               (bounds%begp:bounds%endp), stat=ier)
    allocate(prod10_cflux             (bounds%begp:bounds%endp), stat=ier)
    allocate(prod100_cflux            (bounds%begp:bounds%endp), stat=ier)
    allocate(conv_nflux               (bounds%begp:bounds%endp), stat=ier)
    allocate(prod10_nflux             (bounds%begp:bounds%endp), stat=ier)
    allocate(prod100_nflux            (bounds%begp:bounds%endp), stat=ier)

    ! Allocate P arrays
    allocate(dwt_leafp_seed           (bounds%begp:bounds%endp), stat=ier)
    allocate(dwt_deadstemp_seed       (bounds%begp:bounds%endp), stat=ier)
    allocate(dwt_ppool_seed           (bounds%begp:bounds%endp), stat=ier)
    allocate(dwt_frootp_to_litter     (bounds%begp:bounds%endp), stat=ier)
    allocate(dwt_livecrootp_to_litter (bounds%begp:bounds%endp), stat=ier)
    allocate(dwt_deadcrootp_to_litter (bounds%begp:bounds%endp), stat=ier)
    allocate(conv_pflux               (bounds%begp:bounds%endp), stat=ier)
    allocate(prod10_pflux             (bounds%begp:bounds%endp), stat=ier)
    allocate(prod100_pflux            (bounds%begp:bounds%endp), stat=ier)

    if ( use_c13 ) then
       allocate(dwt_leafc13_seed           (bounds%begp:bounds%endp), stat=ier)
       allocate(dwt_deadstemc13_seed       (bounds%begp:bounds%endp), stat=ier)
       allocate(dwt_frootc13_to_litter     (bounds%begp:bounds%endp), stat=ier)
       allocate(dwt_livecrootc13_to_litter (bounds%begp:bounds%endp), stat=ier)
       allocate(dwt_deadcrootc13_to_litter (bounds%begp:bounds%endp), stat=ier)
       allocate(conv_c13flux               (bounds%begp:bounds%endp), stat=ier)
       allocate(prod10_c13flux             (bounds%begp:bounds%endp), stat=ier)
       allocate(prod100_c13flux            (bounds%begp:bounds%endp), stat=ier)
    endif
    if ( use_c14 ) then
       allocate(dwt_leafc14_seed           (bounds%begp:bounds%endp), stat=ier)
       allocate(dwt_deadstemc14_seed       (bounds%begp:bounds%endp), stat=ier)
       allocate(dwt_frootc14_to_litter     (bounds%begp:bounds%endp), stat=ier)
       allocate(dwt_livecrootc14_to_litter (bounds%begp:bounds%endp), stat=ier)
       allocate(dwt_deadcrootc14_to_litter (bounds%begp:bounds%endp), stat=ier)
       allocate(conv_c14flux               (bounds%begp:bounds%endp), stat=ier)
       allocate(prod10_c14flux             (bounds%begp:bounds%endp), stat=ier)
       allocate(prod100_c14flux            (bounds%begp:bounds%endp), stat=ier)
    endif
    
    ! Get time step
    dt = real( get_step_size(), r8 )
    
    do p = bounds%begp,bounds%endp
       c = veg_pp%column(p)
       ! initialize all the pft-level local flux arrays
       dwt_leafc_seed(p)           = 0._r8
       dwt_deadstemc_seed(p)       = 0._r8
       dwt_frootc_to_litter(p)     = 0._r8
       dwt_livecrootc_to_litter(p) = 0._r8
       dwt_deadcrootc_to_litter(p) = 0._r8
       conv_cflux(p)               = 0._r8
       prod10_cflux(p)             = 0._r8
       prod100_cflux(p)            = 0._r8
       
       dwt_leafn_seed(p)           = 0._r8
       dwt_deadstemn_seed(p)       = 0._r8
       dwt_npool_seed(p)           = 0._r8
       dwt_frootn_to_litter(p)     = 0._r8
       dwt_livecrootn_to_litter(p) = 0._r8
       dwt_deadcrootn_to_litter(p) = 0._r8
       conv_nflux(p)               = 0._r8
       prod10_nflux(p)             = 0._r8
       prod100_nflux(p)            = 0._r8
       
       dwt_leafp_seed(p)           = 0._r8
       dwt_deadstemp_seed(p)       = 0._r8
       dwt_ppool_seed(p)           = 0._r8
       dwt_frootp_to_litter(p)     = 0._r8
       dwt_livecrootp_to_litter(p) = 0._r8
       dwt_deadcrootp_to_litter(p) = 0._r8
       conv_pflux(p)               = 0._r8
       prod10_pflux(p)             = 0._r8
       prod100_pflux(p)            = 0._r8

       if ( use_c13 ) then
          dwt_leafc13_seed(p)           = 0._r8
          dwt_deadstemc13_seed(p)       = 0._r8
          dwt_frootc13_to_litter(p)     = 0._r8
          dwt_livecrootc13_to_litter(p) = 0._r8
          dwt_deadcrootc13_to_litter(p) = 0._r8
          conv_c13flux(p)               = 0._r8
          prod10_c13flux(p)             = 0._r8
          prod100_c13flux(p)            = 0._r8
       endif
       
       if ( use_c14 ) then
          dwt_leafc14_seed(p)           = 0._r8
          dwt_deadstemc14_seed(p)       = 0._r8
          dwt_frootc14_to_litter(p)     = 0._r8
          dwt_livecrootc14_to_litter(p) = 0._r8
          dwt_deadcrootc14_to_litter(p) = 0._r8
          conv_c14flux(p)               = 0._r8
          prod10_c14flux(p)             = 0._r8
          prod100_c14flux(p)            = 0._r8
       endif
       
       l = veg_pp%landunit(p)
       if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then
          
          ! calculate the change in weight for the timestep
          dwt = veg_pp%wtcol(p)-prior_weights%pwtcol(p)
          cnstate_vars%lfpftd_patch(p) = -dwt

          ! Patches for which weight increases on this timestep
          if (dwt > 0._r8) then
             
             ! first identify Patches that are initiating on this timestep
             ! and set all the necessary state and flux variables
             if (prior_weights%pwtcol(p) == 0._r8) then
                
                ! set initial conditions for PFT that is being initiated
                ! in this time step.  Based on the settings in cnIniTimeVar.
                
                ! pft-level carbon state variables
                call CarbonStateVarsInit     (cs, p)
                call NitrogenStateVarsInit   (ns, p)
                call PhosphorusStateVarsInit (ps, p)
                call CanopyStateVarsInit     (canopystate_vars, p)
                call CNStateVarsInit         (cnstate_vars, p, c)
                call CabronFluxVarsInit      (cf, p)
                call NitrogenFluxVarsInit    (nf, p)
                call PhosphorusFluxVarsInit  (pf, p)
                
                if ( use_c13 ) then
                   call CarbonStateVarsInit(c13_cs, p)
                endif
                
                if ( use_c14 ) then
                   call CarbonStateVarsInit(c14_cs, p)
                endif
                                 
                ! add phosphorus related variables

                call photosyns_vars%NewPatchinit(p)
                
             end if  ! end initialization of new pft
                                                    
          end if       ! weight decreasing
       end if           ! is soil
    end do               ! patch loop
    
    call CStateDynamicPatchAdjustments( &
         bounds,                        &
         num_soilp_with_inactive,       &
         filter_soilp_with_inactive,    &
         prior_weights,                 &
         patch_state_updater,           &
         cs,                            &
         dwt_leafc_seed,                &
         dwt_deadstemc_seed,            &
         conv_cflux,                    &
         dwt_frootc_to_litter,          &
         dwt_livecrootc_to_litter,      &
         dwt_deadcrootc_to_litter,      &
         prod10_cflux,                  &
         prod100_cflux                  &
         )

    if (use_c13) then
       call CStateDynamicPatchAdjustments( &
            bounds,                        &
            num_soilp_with_inactive,       &
            filter_soilp_with_inactive,    &
            prior_weights,                 &
            patch_state_updater,           &
            c13_cs,                        &
            dwt_leafc13_seed,              &
            dwt_deadstemc13_seed,          &
            conv_c13flux,                  &
            dwt_frootc13_to_litter,        &
            dwt_livecrootc13_to_litter,    &
            dwt_deadcrootc13_to_litter,    &
            prod10_c13flux,                &
            prod100_c13flux                &
            )
    endif

    if (use_c14) then
       call CStateDynamicPatchAdjustments( &
            bounds,                        &
            num_soilp_with_inactive,       &
            filter_soilp_with_inactive,    &
            prior_weights,                 &
            patch_state_updater,           &
            c14_cs,                        &
            dwt_leafc14_seed,              &
            dwt_deadstemc14_seed,          &
            conv_c14flux,                  &
            dwt_frootc14_to_litter,        &
            dwt_livecrootc14_to_litter,    &
            dwt_deadcrootc14_to_litter,    &
            prod10_c14flux,                &
            prod100_c14flux                &
            )
    endif

    call NStateDynamicPatchAdjustments( &
         bounds,                        &
         num_soilp_with_inactive,       &
         filter_soilp_with_inactive,    &
         prior_weights,                 &
         patch_state_updater,           &
         ns,                            &
         dwt_leafn_seed,                &
         dwt_deadstemn_seed,            &
         dwt_npool_seed,                &
         conv_nflux,                    &
         dwt_frootn_to_litter,          &
         dwt_livecrootn_to_litter,      &
         dwt_deadcrootn_to_litter,      &
         prod10_nflux,                  &
         prod100_nflux                  &
         )

    call PStateDynamicPatchAdjustments( &
         bounds,                        &
         num_soilp_with_inactive,       &
         filter_soilp_with_inactive,    &
         prior_weights,                 &
         patch_state_updater,           &
         ps,                            &
         dwt_leafp_seed,                &
         dwt_deadstemp_seed,            &
         dwt_ppool_seed,                &
         conv_pflux,                    &
         dwt_frootp_to_litter,          &
         dwt_livecrootp_to_litter,      &
         dwt_deadcrootp_to_litter,      &
         prod10_pflux,                  &
         prod100_pflux                  &
         )

    ! calculate column-level seeding fluxes
    do pi = 1,max_patch_per_col
       do c = bounds%begc, bounds%endc
          if ( pi <=  col_pp%npfts(c) ) then
             p = col_pp%pfti(c) + pi - 1
             
             ! C fluxe
             cf%dwt_seedc_to_leaf_col(c) = cf%dwt_seedc_to_leaf_col(c) + dwt_leafc_seed(p)/dt
             cf%dwt_seedc_to_deadstem_col(c) = cf%dwt_seedc_to_deadstem_col(c) &
                  + dwt_deadstemc_seed(p)/dt

             if ( use_c13 ) then
                c13_cf%dwt_seedc_to_leaf_col(c) = c13_cf%dwt_seedc_to_leaf_col(c) + dwt_leafc13_seed(p)/dt
                c13_cf%dwt_seedc_to_deadstem_col(c) = c13_cf%dwt_seedc_to_deadstem_col(c) &
                     + dwt_deadstemc13_seed(p)/dt
             endif
             
             if ( use_c14 ) then	
                c14_cf%dwt_seedc_to_leaf_col(c) = c14_cf%dwt_seedc_to_leaf_col(c) + dwt_leafc14_seed(p)/dt
                c14_cf%dwt_seedc_to_deadstem_col(c) = c14_cf%dwt_seedc_to_deadstem_col(c) &
                     + dwt_deadstemc14_seed(p)/dt
             endif
             
             ! N fluxes
             nf%dwt_seedn_to_leaf_col(c) = nf%dwt_seedn_to_leaf_col(c) + dwt_leafn_seed(p)/dt
             nf%dwt_seedn_to_deadstem_col(c) = nf%dwt_seedn_to_deadstem_col(c) &
                  + dwt_deadstemn_seed(p)/dt
             nf%dwt_seedn_to_npool_col(c) = nf%dwt_seedn_to_npool_col(c) &
                  + dwt_npool_seed(p)/dt
             ! P fluxes
             pf%dwt_seedp_to_leaf_col(c) = pf%dwt_seedp_to_leaf_col(c) + dwt_leafp_seed(p)/dt
             pf%dwt_seedp_to_deadstem_col(c) = pf%dwt_seedp_to_deadstem_col(c) &
                  + dwt_deadstemp_seed(p)/dt
             pf%dwt_seedp_to_ppool_col(c) = pf%dwt_seedp_to_ppool_col(c) &
                  + dwt_ppool_seed(p)/dt
          end if
       end do
    end do
    
    
    ! calculate pft-to-column for fluxes into litter and CWD pools
    do j = 1, nlevdecomp
       do pi = 1,max_patch_per_col
          do c = bounds%begc, bounds%endc
             if ( pi <=  col_pp%npfts(c) ) then
                p = col_pp%pfti(c) + pi - 1

                froot   = cnstate_vars%froot_prof_patch(p,j)
                croot   = cnstate_vars%croot_prof_patch(p,j)
                fr_flab = veg_vp%fr_flab(veg_pp%itype(p))
                fr_fcel = veg_vp%fr_fcel(veg_pp%itype(p))
                fr_flig = veg_vp%fr_flig(veg_pp%itype(p))
                
                
                ! fine root litter carbon fluxes
                cf%dwt_frootc_to_litr_met_c_col(c,j) = &
                     cf%dwt_frootc_to_litr_met_c_col(c,j) + &
                     (dwt_frootc_to_litter(p)* fr_flab)/dt * froot

                cf%dwt_frootc_to_litr_cel_c_col(c,j) = &
                     cf%dwt_frootc_to_litr_cel_c_col(c,j) + &
                     (dwt_frootc_to_litter(p)* fr_fcel)/dt * froot

                cf%dwt_frootc_to_litr_lig_c_col(c,j) = &
                     cf%dwt_frootc_to_litr_lig_c_col(c,j) + &
                     (dwt_frootc_to_litter(p)* fr_flig)/dt * froot
                
                
                ! fine root litter nitrogen fluxes
                nf%dwt_frootn_to_litr_met_n_col(c,j) = &
                     nf%dwt_frootn_to_litr_met_n_col(c,j) + &
                     (dwt_frootn_to_litter(p)* fr_flab)/dt * froot
                nf%dwt_frootn_to_litr_cel_n_col(c,j) = &

                     nf%dwt_frootn_to_litr_cel_n_col(c,j) + &
                     (dwt_frootn_to_litter(p)* fr_fcel)/dt * froot

                nf%dwt_frootn_to_litr_lig_n_col(c,j) = &
                     nf%dwt_frootn_to_litr_lig_n_col(c,j) + &
                     (dwt_frootn_to_litter(p)* fr_flig)/dt * froot
                

                ! fine root litter phosphorus fluxes
                pf%dwt_frootp_to_litr_met_p_col(c,j) = &
                     pf%dwt_frootp_to_litr_met_p_col(c,j) + &
                     (dwt_frootp_to_litter(p)* fr_flab)/dt * froot
                pf%dwt_frootp_to_litr_cel_p_col(c,j) = &

                     pf%dwt_frootp_to_litr_cel_p_col(c,j) + &
                     (dwt_frootp_to_litter(p)* fr_fcel)/dt * froot

                pf%dwt_frootp_to_litr_lig_p_col(c,j) = &
                     pf%dwt_frootp_to_litr_lig_p_col(c,j) + &
                     (dwt_frootp_to_litter(p)* fr_flig)/dt * froot

                ! livecroot fluxes to cwd
                cf%dwt_livecrootc_to_cwdc_col(c,j) = &
                     cf%dwt_livecrootc_to_cwdc_col(c,j) + &
                     (dwt_livecrootc_to_litter(p))/dt * croot

                nf%dwt_livecrootn_to_cwdn_col(c,j) = &
                     nf%dwt_livecrootn_to_cwdn_col(c,j) + &
                     (dwt_livecrootn_to_litter(p))/dt * croot
                
                pf%dwt_livecrootp_to_cwdp_col(c,j) = &
                     pf%dwt_livecrootp_to_cwdp_col(c,j) + &
                     (dwt_livecrootp_to_litter(p))/dt * croot

                ! deadcroot fluxes to cwd
                cf%dwt_deadcrootc_to_cwdc_col(c,j) = &
                     cf%dwt_deadcrootc_to_cwdc_col(c,j) + &
                     (dwt_deadcrootc_to_litter(p))/dt * croot

                nf%dwt_deadcrootn_to_cwdn_col(c,j) = &
                     nf%dwt_deadcrootn_to_cwdn_col(c,j) + &
                     (dwt_deadcrootn_to_litter(p))/dt * croot
             
                pf%dwt_deadcrootp_to_cwdp_col(c,j) = &
                     pf%dwt_deadcrootp_to_cwdp_col(c,j) + &
                     (dwt_deadcrootp_to_litter(p))/dt * croot

                if ( use_c13 ) then
                   ! C13 fine root litter fluxes
                   c13_cf%dwt_frootc_to_litr_met_c_col(c,j) = &
                        c13_cf%dwt_frootc_to_litr_met_c_col(c,j) + &
                        (dwt_frootc13_to_litter(p)* fr_flab)/dt * froot

                   c13_cf%dwt_frootc_to_litr_cel_c_col(c,j) = &
                        c13_cf%dwt_frootc_to_litr_cel_c_col(c,j) + &
                        (dwt_frootc13_to_litter(p)* fr_fcel)/dt * froot

                   c13_cf%dwt_frootc_to_litr_lig_c_col(c,j) = &
                        c13_cf%dwt_frootc_to_litr_lig_c_col(c,j) + &
                        (dwt_frootc13_to_litter(p)* fr_flig)/dt * froot

                   ! livecroot fluxes to cwd
                   c13_cf%dwt_livecrootc_to_cwdc_col(c,j) = &
                        c13_cf%dwt_livecrootc_to_cwdc_col(c,j) + &
                        (dwt_livecrootc13_to_litter(p))/dt * croot

                   ! deadcroot fluxes to cwd
                   c13_cf%dwt_deadcrootc_to_cwdc_col(c,j) = &
                        c13_cf%dwt_deadcrootc_to_cwdc_col(c,j) + &
                        (dwt_deadcrootc13_to_litter(p))/dt * croot
                   
                endif
                
                if ( use_c14 ) then                   
                   ! C14 fine root litter fluxes
                   c14_cf%dwt_frootc_to_litr_met_c_col(c,j) = &
                        c14_cf%dwt_frootc_to_litr_met_c_col(c,j) + &
                        (dwt_frootc14_to_litter(p)* fr_flab)/dt * froot

                   c14_cf%dwt_frootc_to_litr_cel_c_col(c,j) = &
                        c14_cf%dwt_frootc_to_litr_cel_c_col(c,j) + &
                        (dwt_frootc14_to_litter(p)* fr_fcel)/dt * froot

                   c14_cf%dwt_frootc_to_litr_lig_c_col(c,j) = &
                        c14_cf%dwt_frootc_to_litr_lig_c_col(c,j) + &
                        (dwt_frootc14_to_litter(p)* fr_flig)/dt * froot

                   ! livecroot fluxes to cwd
                   c14_cf%dwt_livecrootc_to_cwdc_col(c,j) = &
                        c14_cf%dwt_livecrootc_to_cwdc_col(c,j) + &
                        (dwt_livecrootc14_to_litter(p))/dt * croot

                   ! deadcroot fluxes to cwd
                   c14_cf%dwt_deadcrootc_to_cwdc_col(c,j) = &
                        c14_cf%dwt_deadcrootc_to_cwdc_col(c,j) + &
                        (dwt_deadcrootc14_to_litter(p))/dt * croot
                endif
                
             end if
          end do
       end do
    end do
    ! calculate pft-to-column for fluxes into product pools and conversion flux
    do pi = 1,max_patch_per_col
       do c = bounds%begc,bounds%endc
          if (pi <= col_pp%npfts(c)) then
             p = col_pp%pfti(c) + pi - 1
             
             ! column-level fluxes are accumulated as positive fluxes.
             ! column-level C flux updates
             cf%dwt_conv_cflux_col(c) = cf%dwt_conv_cflux_col(c) - conv_cflux(p)/dt
             cf%dwt_prod10c_gain_col(c) = cf%dwt_prod10c_gain_col(c) - prod10_cflux(p)/dt
             cf%dwt_prod100c_gain_col(c) = cf%dwt_prod100c_gain_col(c) - prod100_cflux(p)/dt

             ! These magic constants should be replaced with: nbrdlf_evr_trp_tree and nbrdlf_dcd_trp_tree
             if(veg_pp%itype(p)==4.or.veg_pp%itype(p)==6)then
                cf%lf_conv_cflux_col(c) = cf%lf_conv_cflux_col(c) - conv_cflux(p)/dt
             end if
             
             if ( use_c13 ) then
                ! C13 column-level flux updates
                c13_cf%dwt_conv_cflux_col(c) = c13_cf%dwt_conv_cflux_col(c) - conv_c13flux(p)/dt
                c13_cf%dwt_prod10c_gain_col(c) = c13_cf%dwt_prod10c_gain_col(c) - prod10_c13flux(p)/dt
                c13_cf%dwt_prod100c_gain_col(c) = c13_cf%dwt_prod100c_gain_col(c) - prod100_c13flux(p)/dt
             endif
             
             if ( use_c14 ) then
                ! C14 column-level flux updates
                c14_cf%dwt_conv_cflux_col(c) = c14_cf%dwt_conv_cflux_col(c) - conv_c14flux(p)/dt
                c14_cf%dwt_prod10c_gain_col(c) = c14_cf%dwt_prod10c_gain_col(c) - prod10_c14flux(p)/dt
                c14_cf%dwt_prod100c_gain_col(c) = c14_cf%dwt_prod100c_gain_col(c) - prod100_c14flux(p)/dt
             endif
             
             ! column-level N flux updates
             nf%dwt_conv_nflux_col(c) = nf%dwt_conv_nflux_col(c) - conv_nflux(p)/dt
             nf%dwt_prod10n_gain_col(c) = nf%dwt_prod10n_gain_col(c) - prod10_nflux(p)/dt
             nf%dwt_prod100n_gain_col(c) = nf%dwt_prod100n_gain_col(c) - prod100_nflux(p)/dt
             
             ! column-level P flux updates

             pf%dwt_conv_pflux_col(c) = pf%dwt_conv_pflux_col(c) - conv_pflux(p)/dt
             pf%dwt_prod10p_gain_col(c) = pf%dwt_prod10p_gain_col(c) - prod10_pflux(p)/dt
             pf%dwt_prod100p_gain_col(c) = pf%dwt_prod100p_gain_col(c) - prod100_pflux(p)/dt

          end if
       end do
    end do
    
    ! Deallocate pft-level flux arrays
    deallocate(dwt_leafc_seed)
    deallocate(dwt_leafn_seed)
    deallocate(dwt_deadstemc_seed)
    deallocate(dwt_deadstemn_seed)
    deallocate(dwt_npool_seed)
    deallocate(dwt_frootc_to_litter)
    deallocate(dwt_livecrootc_to_litter)
    deallocate(dwt_deadcrootc_to_litter)
    deallocate(dwt_frootn_to_litter)
    deallocate(dwt_livecrootn_to_litter)
    deallocate(dwt_deadcrootn_to_litter)
    deallocate(conv_cflux)
    deallocate(prod10_cflux)
    deallocate(prod100_cflux)
    deallocate(conv_nflux)
    deallocate(prod10_nflux)
    deallocate(prod100_nflux)

    deallocate(dwt_leafp_seed)
    deallocate(dwt_deadstemp_seed)
    deallocate(dwt_ppool_seed)
    deallocate(dwt_frootp_to_litter)
    deallocate(dwt_livecrootp_to_litter)
    deallocate(dwt_deadcrootp_to_litter)
    deallocate(conv_pflux)
    deallocate(prod10_pflux)
    deallocate(prod100_pflux)
             
    if ( use_c13 ) then
       deallocate(dwt_leafc13_seed)
       deallocate(dwt_deadstemc13_seed)
       deallocate(dwt_frootc13_to_litter)
       deallocate(dwt_livecrootc13_to_litter)
       deallocate(dwt_deadcrootc13_to_litter)
       deallocate(conv_c13flux)
       deallocate(prod10_c13flux)
       deallocate(prod100_c13flux)
    endif
             
    if ( use_c14 ) then
       deallocate(dwt_leafc14_seed)
       deallocate(dwt_deadstemc14_seed)
       deallocate(dwt_frootc14_to_litter)
       deallocate(dwt_livecrootc14_to_litter)
       deallocate(dwt_deadcrootc14_to_litter)
       deallocate(conv_c14flux)
       deallocate(prod10_c14flux)
       deallocate(prod100_c14flux)
    endif
    
   end associate
 end subroutine dyn_cnbal_patch
 
 !-----------------------------------------------------------------------
 subroutine CarbonStateVarsInit(cs, p)
   !
   ! !DESCRIPTION:
   ! Initializes p-th patch of carbonstate_type
   !
   implicit none
   !
   ! !ARGUMENT
   type(carbonstate_type), intent(inout) :: cs
   integer               , intent(in)    :: p
   
   cs%leafc_patch(p)              = 0._r8
   cs%leafc_storage_patch(p)      = 0._r8
   cs%leafc_xfer_patch(p)         = 0._r8
   cs%frootc_patch(p)             = 0._r8
   cs%frootc_storage_patch(p)     = 0._r8
   cs%frootc_xfer_patch(p)        = 0._r8
   cs%livestemc_patch(p)          = 0._r8
   cs%livestemc_storage_patch(p)  = 0._r8
   cs%livestemc_xfer_patch(p)     = 0._r8
   cs%deadstemc_patch(p)          = 0._r8
   cs%deadstemc_storage_patch(p)  = 0._r8
   cs%deadstemc_xfer_patch(p)     = 0._r8
   cs%livecrootc_patch(p)         = 0._r8
   cs%livecrootc_storage_patch(p) = 0._r8
   cs%livecrootc_xfer_patch(p)    = 0._r8
   cs%deadcrootc_patch(p)         = 0._r8
   cs%deadcrootc_storage_patch(p) = 0._r8
   cs%deadcrootc_xfer_patch(p)    = 0._r8
   cs%gresp_storage_patch(p)      = 0._r8
   cs%gresp_xfer_patch(p)         = 0._r8
   cs%cpool_patch(p)              = 0._r8
   cs%xsmrpool_patch(p)           = 0._r8
   cs%ctrunc_patch(p)             = 0._r8
   cs%dispvegc_patch(p)           = 0._r8
   cs%storvegc_patch(p)           = 0._r8
   cs%totvegc_patch(p)            = 0._r8
   cs%totpftc_patch(p)            = 0._r8

 end subroutine CarbonStateVarsInit

 !-----------------------------------------------------------------------
 subroutine NitrogenStateVarsInit(ns, p)
   !
   ! !DESCRIPTION:
   ! Initializes p-th patch of nitrogenstate_type
   !
   implicit none
   !
   ! !ARGUMENT
   type(nitrogenstate_type), intent(inout) :: ns
   integer                 , intent(in)    :: p
   
   ns%leafn_patch(p)              = 0._r8
   ns%leafn_storage_patch(p)      = 0._r8
   ns%leafn_xfer_patch(p)         = 0._r8
   ns%frootn_patch(p)             = 0._r8
   ns%frootn_storage_patch(p)     = 0._r8
   ns%frootn_xfer_patch(p)        = 0._r8
   ns%livestemn_patch(p)          = 0._r8
   ns%livestemn_storage_patch(p)  = 0._r8
   ns%livestemn_xfer_patch(p)     = 0._r8
   ns%deadstemn_patch(p)          = 0._r8
   ns%deadstemn_storage_patch(p)  = 0._r8
   ns%deadstemn_xfer_patch(p)     = 0._r8
   ns%livecrootn_patch(p)         = 0._r8
   ns%livecrootn_storage_patch(p) = 0._r8
   ns%livecrootn_xfer_patch(p)    = 0._r8
   ns%deadcrootn_patch(p)         = 0._r8
   ns%deadcrootn_storage_patch(p) = 0._r8
   ns%deadcrootn_xfer_patch(p)    = 0._r8
   ns%retransn_patch(p)           = 0._r8
   ns%npool_patch(p)              = 0._r8
   ns%ntrunc_patch(p)             = 0._r8
   ns%dispvegn_patch(p)           = 0._r8
   ns%storvegn_patch(p)           = 0._r8
   ns%totvegn_patch(p)            = 0._r8
   ns%totpftn_patch (p)           = 0._r8

 end subroutine NitrogenStateVarsInit

 !-----------------------------------------------------------------------
 subroutine PhosphorusStateVarsInit(ps, p)
   !
   ! !DESCRIPTION:
   ! Initializes p-th patch of phosphorusstate_type
   !
   implicit none
   !
   ! !ARGUMENT
   type(phosphorusstate_type), intent(inout) :: ps
   integer                   , intent(in)    :: p

   ps%leafp_patch(p)              = 0._r8
   ps%leafp_storage_patch(p)      = 0._r8
   ps%leafp_xfer_patch(p)         = 0._r8
   ps%frootp_patch(p)             = 0._r8
   ps%frootp_storage_patch(p)     = 0._r8
   ps%frootp_xfer_patch(p)        = 0._r8
   ps%livestemp_patch(p)          = 0._r8
   ps%livestemp_storage_patch(p)  = 0._r8
   ps%livestemp_xfer_patch(p)     = 0._r8
   ps%deadstemp_patch(p)          = 0._r8
   ps%deadstemp_storage_patch(p)  = 0._r8
   ps%deadstemp_xfer_patch(p)     = 0._r8
   ps%livecrootp_patch(p)         = 0._r8
   ps%livecrootp_storage_patch(p) = 0._r8
   ps%livecrootp_xfer_patch(p)    = 0._r8
   ps%deadcrootp_patch(p)         = 0._r8
   ps%deadcrootp_storage_patch(p) = 0._r8
   ps%deadcrootp_xfer_patch(p)    = 0._r8
   ps%retransp_patch(p)           = 0._r8
   ps%ppool_patch(p)              = 0._r8
   ps%ptrunc_patch(p)             = 0._r8
   ps%dispvegp_patch(p)           = 0._r8
   ps%storvegp_patch(p)           = 0._r8
   ps%totvegp_patch(p)            = 0._r8
   ps%totpftp_patch (p)           = 0._r8

 end subroutine PhosphorusStateVarsInit

 !-----------------------------------------------------------------------
 subroutine CanopyStateVarsInit(canopystate_vars, p)
   !
   ! !DESCRIPTION:
   ! Initializes p-th patch of canopystate_type
   !
   implicit none
   !
   ! !ARGUMENT
   type(canopystate_type), intent(inout) :: canopystate_vars
   integer               , intent(in)    :: p

   canopystate_vars%laisun_patch(p) = 0._r8
   canopystate_vars%laisha_patch(p) = 0._r8

 end subroutine CanopyStateVarsInit

 !-----------------------------------------------------------------------
 subroutine CNStateVarsInit(cnstate_vars, p, c)
   !
   ! !DESCRIPTION:
   ! Initializes p-th patch of cnstate_type
   !
   use clm_varcon, only : c14ratio
   implicit none
   !
   ! !ARGUMENT
   type(cnstate_type), intent(inout) :: cnstate_vars
   integer           , intent(in)    :: p
   integer           , intent(in)    :: c

   cnstate_vars%dormant_flag_patch(p)          = 1._r8
   cnstate_vars%days_active_patch(p)           = 0._r8
   cnstate_vars%onset_flag_patch(p)            = 0._r8
   cnstate_vars%onset_counter_patch(p)         = 0._r8
   cnstate_vars%onset_gddflag_patch(p)         = 0._r8
   cnstate_vars%onset_fdd_patch(p)             = 0._r8
   cnstate_vars%onset_gdd_patch(p)             = 0._r8
   cnstate_vars%onset_swi_patch(p)             = 0._r8
   cnstate_vars%offset_flag_patch(p)           = 0._r8
   cnstate_vars%offset_counter_patch(p)        = 0._r8
   cnstate_vars%offset_fdd_patch(p)            = 0._r8
   cnstate_vars%offset_swi_patch(p)            = 0._r8
   cnstate_vars%lgsf_patch(p)                  = 0._r8
   cnstate_vars%bglfr_patch(p)                 = 0._r8
   cnstate_vars%bglfr_leaf_patch(p)            = 0._r8
   cnstate_vars%bglfr_froot_patch(p)           = 0._r8
   cnstate_vars%bgtr_patch(p)                  = 0._r8
   cnstate_vars%annavg_t2m_patch(p)            = cnstate_vars%annavg_t2m_col(c)
   cnstate_vars%tempavg_t2m_patch(p)           = 0._r8
   cnstate_vars%alloc_pnow_patch(p)            = 1._r8
   cnstate_vars%c_allometry_patch(p)           = 0._r8
   cnstate_vars%n_allometry_patch(p)           = 0._r8
   cnstate_vars%p_allometry_patch(p)           = 0._r8
   cnstate_vars%tempsum_potential_gpp_patch(p) = 0._r8
   cnstate_vars%annsum_potential_gpp_patch(p)  = 0._r8
   cnstate_vars%tempmax_retransn_patch(p)      = 0._r8
   cnstate_vars%annmax_retransn_patch(p)       = 0._r8
   cnstate_vars%downreg_patch(p)               = 0._r8

   cnstate_vars%tempmax_retransp_patch(p)      = 0._r8
   cnstate_vars%annmax_retransp_patch(p)       = 0._r8

   if ( use_c14 ) then
      cnstate_vars%rc14_atm_patch(p) = c14ratio
      cnstate_vars%rc14_atm_patch(p) = 0._r8
   endif

 end subroutine CNStateVarsInit

 !-----------------------------------------------------------------------
 subroutine CabronFluxVarsInit(cf, p)
   !
   ! !DESCRIPTION:
   ! Initializes p-th patch of carbonflux_type
   !
   use clm_varcon, only : c13ratio
   !
   implicit none
   !
   ! !ARGUMENT
   type(carbonflux_type), intent(inout) :: cf
   integer              , intent(in)    :: p

   cf%xsmrpool_recover_patch(p)      = 0._r8
   cf%plant_calloc_patch(p)          = 0._r8
   cf%excess_cflux_patch(p)          = 0._r8
   cf%prev_leafc_to_litter_patch(p)  = 0._r8
   cf%prev_frootc_to_litter_patch(p) = 0._r8
   cf%availc_patch(p)                = 0._r8
   cf%gpp_before_downreg_patch(p)    = 0._r8
   cf%tempsum_npp_patch(p)           = 0._r8
   cf%annsum_npp_patch(p)            = 0._r8

   if ( use_c13 ) then
      cf%xsmrpool_c13ratio_patch(p) = c13ratio
   end if

 end subroutine CabronFluxVarsInit

 !-----------------------------------------------------------------------
 subroutine NitrogenFluxVarsInit(nf, p)
   !
   ! !DESCRIPTION:
   ! Initializes p-th patch of nitrogenflux_type
   !
   implicit none
   !
   ! !ARGUMENT
   type(nitrogenflux_type), intent(inout) :: nf
   integer                , intent(in)    :: p

   nf%plant_ndemand_patch(p)         = 0._r8
   nf%avail_retransn_patch(p)        = 0._r8
   nf%plant_nalloc_patch(p)          = 0._r8

 end subroutine NitrogenFluxVarsInit

 !-----------------------------------------------------------------------
 subroutine PhosphorusFluxVarsInit(pf, p)
   !
   ! !DESCRIPTION:
   ! Initializes p-th patch of phosphorusflux_type
   !
   implicit none
   !
   ! !ARGUMENT
   type(phosphorusflux_type), intent(inout) :: pf
   integer                  , intent(in)    :: p

   pf%plant_pdemand_patch(p)         = 0._r8
   pf%avail_retransp_patch(p)        = 0._r8
   pf%plant_palloc_patch(p)          = 0._r8

 end subroutine PhosphorusFluxVarsInit

 !-----------------------------------------------------------------------
 subroutine ComputeMassLossDueToWtDec(wt_old, wt_new, mass, mass_loss, &
      factor)
   !
   implicit none
   !
   real(r8), intent (in)           :: wt_old
   real(r8), intent (in)           :: wt_new
   real(r8), intent (inout)        :: mass
   real(r8), intent (inout)        :: mass_loss
   real(r8), intent (in), optional :: factor
   !
   real(r8) :: dwt
   real(r8) :: init_state
   real(r8) :: change_state
   real(r8) :: new_state

   dwt = wt_new - wt_old

   init_state   = mass*wt_old
   change_state = mass*dwt
   new_state    = init_state + change_state

   if (.not.present(factor)) then
      if (wt_new /= 0._r8) then
         mass      = new_state/wt_new
         mass_loss = mass_loss + change_state
      else
         mass      = 0._r8
         mass_loss = mass_loss - init_state
      end if
   else
      if (wt_new /= 0._r8) then
         mass      = new_state/wt_new
         mass_loss = mass_loss + change_state*factor
      else
         mass      = 0._r8
         mass_loss = mass_loss - init_state*factor
      end if
   endif
  
 end subroutine ComputeMassLossDueToWtDec
 
  !-----------------------------------------------------------------------
  subroutine LeafProportions(pft_type, ignore_current_state, &
       leaf, leaf_storage, leaf_xfer, &
       pleaf, pstorage, pxfer)
    !
    ! !DESCRIPTION:
    ! Compute leaf proportions (leaf, storage and xfer)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer , intent(in)  :: pft_type
    logical , intent(in)  :: ignore_current_state
    real(r8), intent(in)  :: leaf         ! g/m2 leaf C or N
    real(r8), intent(in)  :: leaf_storage ! g/m2 leaf C or N storage
    real(r8), intent(in)  :: leaf_xfer    ! g/m2 leaf C or N transfer

    real(r8), intent(out) :: pleaf        ! proportion in leaf itself
    real(r8), intent(out) :: pstorage     ! proportion in leaf storage
    real(r8), intent(out) :: pxfer        ! proportion in leaf xfer
    !
    ! !LOCAL VARIABLES:
    real(r8) :: tot_leaf

    character(len=*), parameter :: subname = 'LeafProportions'
    !-----------------------------------------------------------------------

    tot_leaf = leaf + leaf_storage + leaf_xfer
    pleaf    = 0._r8
    pstorage = 0._r8
    pxfer    = 0._r8

    if (tot_leaf == 0._r8 .or. ignore_current_state) then
       if (veg_vp%evergreen(pft_type) == 1._r8) then
          pleaf    = 1._r8
       else
          pstorage = 1._r8
       end if
    else
       pleaf    = leaf        /tot_leaf
       pstorage = leaf_storage/tot_leaf
       pxfer    = leaf_xfer   /tot_leaf
    end if

  end subroutine LeafProportions

  !-----------------------------------------------------------------------
  function SpeciesTypeMultiplier(species, pft_type, component) result(multiplier)
    !
    ! !DESCRIPTION:
    ! Returns a multiplier based on the species type. This multiplier is
    ! meant to be applied to some state variable expressed in terms of g C, translating
    ! this value into an appropriate value for c13, c14, n or p.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    real(r8)            :: multiplier ! function result
    integer, intent(in) :: species    ! which C/N species we're operating on; should be one of the values in CNSpeciesMod
    integer, intent(in) :: pft_type
    integer, intent(in) :: component  ! which plant component; should be one of the COMPONENT_* parameters defined in this module
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'SpeciesTypeMultiplier'
    !-----------------------------------------------------------------------

    select case (species)
    case (CN_SPECIES_C12)
       multiplier = 1._r8

    case (CN_SPECIES_C13)
       if (veg_vp%c3psn(pft_type) == 1._r8) then
          multiplier = c3_r2
       else
          multiplier = c4_r2
       end if

    case (CN_SPECIES_C14)
       ! 14c state is initialized assuming initial "modern" 14C of 1.e-12
       multiplier = c14ratio

    case (CN_SPECIES_N)
       select case (component)
       case (COMPONENT_LEAF)
          multiplier = 1._r8 / veg_vp%leafcn(pft_type)
       case (COMPONENT_DEADWOOD)
          multiplier = 1._r8 / veg_vp%deadwdcn(pft_type)
       case (COMPONENT_SEED)
          if (pft_type /= 0) then
             multiplier = 1._r8
          else
             multiplier = 0._r8
          end if
       case default
          write(iulog,*) subname//' ERROR: unknown component: ', component
          call endrun(subname//': unknown component')
       end select

    case (CN_SPECIES_P)
       select case (component)
       case (COMPONENT_LEAF)
          multiplier = 1._r8 / veg_vp%leafcp(pft_type)
       case (COMPONENT_DEADWOOD)
          multiplier = 1._r8 / veg_vp%deadwdcp(pft_type)
       case (COMPONENT_SEED)
          if (pft_type /= 0) then
             multiplier = 1._r8
          else
             multiplier = 0._r8
          end if
       case default
          write(iulog,*) subname//' ERROR: unknown component: ', component
          call endrun(subname//': unknown component')
       end select

    case default
       write(iulog,*) subname//' ERROR: unknown species: ', species
       call endrun(subname//': unknown species')
    end select

  end function SpeciesTypeMultiplier

  !-----------------------------------------------------------------------
  subroutine ComputeSeedAmounts(bounds,                                &
       species,                                                        &
       leaf_patch, leaf_storage_patch, leaf_xfer_patch,                &
       compute_here_patch, ignore_current_state_patch,                 &
       seed_leaf_patch, seed_leaf_storage_patch, seed_leaf_xfer_patch, &
       seed_deadstem_patch, pool_seed_param, pool_seed_patch)
    !
    ! !DESCRIPTION:
    ! Compute seed amounts for patches that increase in area, for various variables, for
    ! the given species (c12, c13, c14 or n).
    !
    ! The output variables are only set for patches inside the filter, where
    ! compute_here_patch is true; for other patches, they remain at their original values.
    !
    ! Note that, regardless of the species, leafc_seed and deadstemc_seed are specified
    ! in terms of gC/m2; these amounts are converted to the amount of the given species
    ! here.
    !
    ! !USES:
    use pftvarcon       , only : noveg
    use landunit_varcon , only : istsoil, istcrop
    !
                                                                                  ! !ARGUMENTS:
    type(bounds_type) , intent(in)     :: bounds
    integer           , intent(in)     :: species                                 ! which C/N species we're operating on; should be one of the values in CNSpeciesMod
    real(r8)          , intent(in)     :: leaf_patch( bounds%begp: )              ! current leaf C or N content (g/m2)
    real(r8)          , intent(in)     :: leaf_storage_patch( bounds%begp: )      ! current leaf C or N storage content (g/m2)
    real(r8)          , intent(in)     :: leaf_xfer_patch( bounds%begp: )         ! current leaf C or N xfer content (g/m2)

                                                                                  ! whether to compute outputs for each patch
    logical           , intent(in)     :: compute_here_patch( bounds%begp: )

                                                                                  ! If ignore_current_state is true, then use default leaf proportions rather than
                                                                                  ! proportions based on current state.
    logical           , intent(in)     :: ignore_current_state_patch( bounds%begp: )

    real(r8)          , intent(inout)  :: seed_leaf_patch( bounds%begp: )         ! seed amount for leaf itself for this species (g/m2)
    real(r8)          , intent(inout)  :: seed_leaf_storage_patch( bounds%begp: ) ! seed amount for leaf storage for this species (g/m2)
    real(r8)          , intent(inout)  :: seed_leaf_xfer_patch( bounds%begp: )    ! seed amount for leaf xfer for this species (g/m2)
    real(r8)          , intent(inout)  :: seed_deadstem_patch( bounds%begp: )     ! seed amount for deadstem for this species (g/m2)
    real(r8), optional, intent(in)     :: pool_seed_param
    real(r8), optional, intent(inout)  :: pool_seed_patch( bounds%begp: )
    !
    ! !LOCAL VARIABLES:
    integer  :: fp, p, c, l
    integer  :: begp, endp
    real(r8) :: my_leaf_seed
    real(r8) :: my_deadstem_seed
    real(r8) :: my_pool_seed
    integer  :: pft_type
    real(r8) :: pleaf
    real(r8) :: pstor
    real(r8) :: pxfer

    character(len=*), parameter :: subname = 'ComputeSeedAmounts'
    !-----------------------------------------------------------------------

    begp = bounds%begp
    endp = bounds%endp

    SHR_ASSERT_ALL((ubound(leaf_patch                 ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(leaf_storage_patch         ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(leaf_xfer_patch            ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(compute_here_patch         ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(ignore_current_state_patch ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(seed_leaf_patch            ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(seed_leaf_storage_patch    ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(seed_leaf_xfer_patch       ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(seed_deadstem_patch        ) == (/endp/)), errMsg(__FILE__, __LINE__))
    
    if (present(pool_seed_patch)) then
       SHR_ASSERT_ALL((ubound(pool_seed_patch         ) == (/endp/)), errMsg(__FILE__, __LINE__))
       if (.not. present(pool_seed_param)) then
          call endrun(subname//': pool_seed_patch can only be provided with pool_seed_param')          
       end if
    end if

    
    do p = bounds%begp,bounds%endp
       c = veg_pp%column(p)

       l = veg_pp%landunit(p)

       if (compute_here_patch(p) .and. (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop)) then

          my_leaf_seed = 0._r8
          my_deadstem_seed = 0._r8

          pft_type = veg_pp%itype(p)

          call LeafProportions( &
               ignore_current_state = ignore_current_state_patch(p), &
               pft_type = pft_type, &
               leaf = leaf_patch(p), &
               leaf_storage = leaf_storage_patch(p), &
               leaf_xfer = leaf_xfer_patch(p), &
               pleaf = pleaf, &
               pstorage = pstor, &
               pxfer = pxfer)

          if (pft_type /= noveg) then
             my_leaf_seed = leafc_seed_param * &
                  SpeciesTypeMultiplier(species, pft_type, COMPONENT_LEAF)
             if (veg_vp%woody(pft_type) == 1._r8) then
                my_deadstem_seed = deadstemc_seed_param * &
                     SpeciesTypeMultiplier(species, pft_type, COMPONENT_DEADWOOD)
             end if
             if (present(pool_seed_param)) then
                my_pool_seed = pool_seed_param     * &
                     SpeciesTypeMultiplier(species, pft_type, COMPONENT_SEED)
             end if
          end if

          seed_leaf_patch(p)         = my_leaf_seed * pleaf
          seed_leaf_storage_patch(p) = my_leaf_seed * pstor
          seed_leaf_xfer_patch(p)    = my_leaf_seed * pxfer
          seed_deadstem_patch(p)     = my_deadstem_seed
          if (present(pool_seed_param)) then
             pool_seed_patch(p) = my_pool_seed
          end if
       end if

    end do

  end subroutine ComputeSeedAmounts

  !-----------------------------------------------------------------------
  subroutine CStateDynamicPatchAdjustments( &
       bounds,                              &
       num_filterp_with_inactive,           &
       filterp_with_inactive,               &
       prior_weights,                       &
       patch_state_updater,                 &
       cs,                                  &
       dwt_leafc_seed,                      &
       dwt_deadstemc_seed,                  &
       conv_cflux,                          &
       dwt_frootc_to_litter,                &
       dwt_livecrootc_to_litter,            &
       dwt_deadcrootc_to_litter,            &
       prod10_cflux,                        &
       prod100_cflux                        &
       )
    !
    ! !DESCRIPTION:
    ! Adjust state variables and compute associated fluxes when patch areas change due to
    ! dynamic landuse
    !
    ! !USES:
    use pftvarcon          , only : pconv, pprod10, pprod100
    use dynPriorWeightsMod , only : prior_weights_type
    use landunit_varcon    , only : istsoil, istcrop
    !
    ! !ARGUMENTS:
    type(bounds_type)              , intent(in)    :: bounds
    integer                        , intent(in)    :: num_filterp_with_inactive ! number of points in filterp_with_inactive
    integer                        , intent(in)    :: filterp_with_inactive(:) ! patch filter that includes inactive points
    type(prior_weights_type)       , intent(in)    :: prior_weights
    type(patch_state_updater_type) , intent(in)    :: patch_state_updater
    type(carbonstate_type)         , intent(inout) :: cs
    real(r8)                       , intent(inout) :: dwt_leafc_seed           (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_deadstemc_seed       (bounds%begp:)
    real(r8)                       , intent(inout) :: conv_cflux               (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_frootc_to_litter     (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_livecrootc_to_litter (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_deadcrootc_to_litter (bounds%begp:)
    real(r8)                       , intent(inout) :: prod10_cflux             (bounds%begp:)
    real(r8)                       , intent(inout) :: prod100_cflux            (bounds%begp:)
    !
    ! !LOCAL VARIABLES:
    integer                     :: begp, endp
    integer                     :: l, c, p
    logical                     :: old_weight_was_zero(bounds%begp:bounds%endp)
    logical                     :: patch_grew(bounds%begp:bounds%endp)

    ! The following are only set for growing patches:
    real(r8)                    :: seed_leafc_patch(bounds%begp:bounds%endp)
    real(r8)                    :: seed_leafc_storage_patch(bounds%begp:bounds%endp)
    real(r8)                    :: seed_leafc_xfer_patch(bounds%begp:bounds%endp)
    real(r8)                    :: seed_deadstemc_patch(bounds%begp:bounds%endp)

    real(r8)                    :: wood_product_cflux(bounds%begp:bounds%endp)
    real(r8)                    :: deadstemc_patch_temp(bounds%begp:bounds%endp)

    character(len=*), parameter :: subname = 'CStateDynamicPatchAdjustments'
    !-----------------------------------------------------------------------

    begp = bounds%begp
    endp = bounds%endp

    SHR_ASSERT_ALL((ubound(dwt_leafc_seed           ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_deadstemc_seed       ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(conv_cflux               ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_frootc_to_litter     ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_livecrootc_to_litter ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_deadcrootc_to_litter ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(prod10_cflux             ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(prod100_cflux            ) == (/endp/)), errMsg(__FILE__, __LINE__))
   
    old_weight_was_zero = patch_state_updater%old_weight_was_zero(bounds)
    patch_grew          = patch_state_updater%patch_grew(bounds)

    call ComputeSeedAmounts(bounds                                        , &
         species                    = cs%species                          , &
         leaf_patch                 = cs%leafc_patch(begp:endp)           , &
         leaf_storage_patch         = cs%leafc_storage_patch(begp:endp)   , &
         leaf_xfer_patch            = cs%leafc_xfer_patch(begp:endp)      , &

         ! Calculations only needed for patches that grew:
         compute_here_patch         = patch_grew(begp:endp)               , &

         ! For patches that previously had zero area, ignore the current state for the
         ! sake of computing leaf proportions:
         ignore_current_state_patch = old_weight_was_zero(begp:endp)      , &

         seed_leaf_patch            = seed_leafc_patch(begp:endp)         , &
         seed_leaf_storage_patch    = seed_leafc_storage_patch(begp:endp) , &
         seed_leaf_xfer_patch       = seed_leafc_xfer_patch(begp:endp)    , &
         seed_deadstem_patch        = seed_deadstemc_patch(begp:endp))

    ! 1) LEAFC_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = cs%leafc_patch   (begp:endp)           , &
         flux_out_grc_area = conv_cflux       (begp:endp)           , &
         seed              = seed_leafc_patch (begp:endp)           , &
         seed_addition     = dwt_leafc_seed   (begp:endp))

    ! 2) LEAFC_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = cs%leafc_storage_patch   (begp:endp)   , &
         flux_out_grc_area = conv_cflux               (begp:endp)   , &
         seed              = seed_leafc_storage_patch (begp:endp)   , &
         seed_addition     = dwt_leafc_seed           (begp:endp))

    ! 3) LEAF_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = cs%leafc_xfer_patch   (begp:endp)      , &
         flux_out_grc_area = conv_cflux            (begp:endp)      , &
         seed              = seed_leafc_xfer_patch (begp:endp)      , &
         seed_addition     = dwt_leafc_seed        (begp:endp))

    ! 4) FROOTC_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = cs%frootc_patch(begp:endp)             , &
         flux_out_grc_area = dwt_frootc_to_litter(begp:endp))

    ! 5) FROOTC_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = cs%frootc_storage_patch(begp:endp)     , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 6) FROOTC_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = cs%frootc_xfer_patch(begp:endp)        , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 7) LIVESTEMC_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = cs%livestemc_patch(begp:endp)          , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 8) LIVESTEMC_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = cs%livestemc_storage_patch(begp:endp)  , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 9) LIVESTEMC_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = cs%livestemc_xfer_patch(begp:endp)     , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 10) LIVECROOTC_PATCH 
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = cs%livecrootc_patch(begp:endp)         , &
         flux_out_grc_area = dwt_frootc_to_litter(begp:endp))

    ! 11) LIVECROOTC_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = cs%livecrootc_storage_patch(begp:endp) , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 12) LIVECROOTC_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = cs%livecrootc_xfer_patch(begp:endp)    , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 13) DEADCROOTC_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = cs%deadcrootc_patch(begp:endp)         , &
         flux_out_grc_area = dwt_frootc_to_litter(begp:endp))

    ! 14) DEADCROOTC_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = cs%deadcrootc_storage_patch(begp:endp) , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 15) DEADCROOT_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = cs%deadcrootc_xfer_patch(begp:endp)    , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 16) GRESP_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = cs%gresp_storage_patch(begp:endp)      , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 17) GRESP_XFER_STORAGE
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = cs%gresp_xfer_patch(begp:endp)         , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 18) CPOOL_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = cs%cpool_patch(begp:endp)              , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 19) XSMRPOOL_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = cs%xsmrpool_patch(begp:endp)           , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 20) CTRUNC_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = cs%ctrunc_patch(begp:endp)             , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 21) DISPVEGC_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = cs%dispvegc_patch(begp:endp))

    ! 22) STORVEGC_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = cs%storvegc_patch(begp:endp))

    ! 23) TOTVEGC_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = cs%totvegc_patch(begp:endp))

    ! 24) TOTPFTC_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = cs%totpftc_patch(begp:endp))

    ! 25) DEADSTEMC_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = cs%deadstemc_storage_patch(begp:endp)  , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 26) DEADSTEMC_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = cs%deadstemc_xfer_patch(begp:endp)     , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 27) PROD10_CFLUX
    wood_product_cflux(begp:endp)      = 0._r8
    deadstemc_patch_temp(begp:endp)    = cs%deadstemc_patch(begp:endp)
    call patch_state_updater%update_patch_state_partition_flux_by_type( &
         bounds                                                        , &
         num_filterp_with_inactive                                     , &
         filterp_with_inactive                                         , &
         flux1_fraction_by_pft_type = pprod10                          , &
         var                        = deadstemc_patch_temp    (begp:endp) , &
         flux1_out                  = prod10_cflux            (begp:endp) , &
         flux2_out                  = wood_product_cflux      (begp:endp) , &
         seed                       = seed_deadstemc_patch    (begp:endp) )

    ! 28) PROD100_CFLUX
    wood_product_cflux(begp:endp)      = 0._r8
    deadstemc_patch_temp(begp:endp)    = cs%deadstemc_patch(begp:endp)
    call patch_state_updater%update_patch_state_partition_flux_by_type( &
         bounds                                                        , &
         num_filterp_with_inactive                                     , &
         filterp_with_inactive                                         , &
         flux1_fraction_by_pft_type = pprod100                         , &
         var                        = deadstemc_patch_temp    (begp:endp) , &
         flux1_out                  = prod100_cflux            (begp:endp) , &
         flux2_out                  = wood_product_cflux      (begp:endp) , &
         seed                       = seed_deadstemc_patch    (begp:endp))

    ! 29) DEADSTEMC_PATCH
    wood_product_cflux(begp:endp)      = 0._r8
    call patch_state_updater%update_patch_state_partition_flux_by_type( &
         bounds                                                        , &
         num_filterp_with_inactive                                     , &
         filterp_with_inactive                                         , &
         flux1_fraction_by_pft_type = pconv                          , &
         var                        = cs%deadstemc_patch(begp:endp) , &
         flux1_out                  = conv_cflux              (begp:endp) , &
         flux2_out                  = wood_product_cflux       (begp:endp) , &
         seed                       = seed_deadstemc_patch    (begp:endp) , &
         seed_addition              = dwt_deadstemc_seed     (begp:endp))

    ! These fluxes are computed as negative quantities, but are expected to be positive,
    ! so flip the signs
    do p = begp, endp
       dwt_frootc_to_litter(p)     = -1._r8 * dwt_frootc_to_litter(p)
       dwt_livecrootc_to_litter(p) = -1._r8 * dwt_livecrootc_to_litter(p)
       dwt_deadcrootc_to_litter(p) = -1._r8 * dwt_deadcrootc_to_litter(p)
    end do

  end subroutine CStateDynamicPatchAdjustments
  

  !-----------------------------------------------------------------------
  subroutine NStateDynamicPatchAdjustments( &
       bounds,                              &
       num_filterp_with_inactive,           &
       filterp_with_inactive,               &
       prior_weights,                       &
       patch_state_updater,                 &
       ns,                                  &
       dwt_leafn_seed,                      &
       dwt_deadstemn_seed,                  &
       dwt_npool_seed,                      &
       conv_nflux,                          &
       dwt_frootn_to_litter,                &
       dwt_livecrootn_to_litter,            &
       dwt_deadcrootn_to_litter,            &
       prod10_nflux,                        &
       prod100_nflux                        &
       )
    !
    ! !DESCRIPTION:
    ! Adjust state variables and compute associated fluxes when patch areas change due to
    ! dynamic landuse
    !
    ! !USES:
    use pftvarcon          , only : pconv, pprod10, pprod100
    use dynPriorWeightsMod , only : prior_weights_type
    use landunit_varcon    , only : istsoil, istcrop
    !
    ! !ARGUMENTS:
    type(bounds_type)              , intent(in)    :: bounds
    integer                        , intent(in)    :: num_filterp_with_inactive
    integer                        , intent(in)    :: filterp_with_inactive(:)
    type(prior_weights_type)       , intent(in)    :: prior_weights
    type(patch_state_updater_type) , intent(in)    :: patch_state_updater
    type(nitrogenstate_type)       , intent(inout) :: ns
    real(r8)                       , intent(inout) :: dwt_leafn_seed           (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_deadstemn_seed       (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_npool_seed           (bounds%begp:)
    real(r8)                       , intent(inout) :: conv_nflux               (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_frootn_to_litter     (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_livecrootn_to_litter (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_deadcrootn_to_litter (bounds%begp:)
    real(r8)                       , intent(inout) :: prod10_nflux             (bounds%begp:)
    real(r8)                       , intent(inout) :: prod100_nflux            (bounds%begp:)
    !
    ! !LOCAL VARIABLES:
    integer                     :: begp, endp
    integer                     :: l, c, p
    logical                     :: old_weight_was_zero      (bounds%begp:bounds%endp)
    logical                     :: patch_grew               (bounds%begp:bounds%endp)

    ! The following are only set for growing patches:
    real(r8)                    :: seed_leafn_patch         (bounds%begp:bounds%endp)
    real(r8)                    :: seed_leafn_storage_patch (bounds%begp:bounds%endp)
    real(r8)                    :: seed_leafn_xfer_patch    (bounds%begp:bounds%endp)
    real(r8)                    :: seed_deadstemn_patch     (bounds%begp:bounds%endp)
    real(r8)                    :: seed_npool_patch         (bounds%begp:bounds%endp)

    real(r8)                    :: wood_product_nflux       (bounds%begp:bounds%endp)
    real(r8)                    :: deadstemn_patch_temp     (bounds%begp:bounds%endp)

    character(len=*), parameter :: subname = 'NStateDynamicPatchAdjustments'
    !-----------------------------------------------------------------------

    begp = bounds%begp
    endp = bounds%endp

    SHR_ASSERT_ALL((ubound(dwt_leafn_seed           ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_deadstemn_seed       ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_npool_seed           ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(conv_nflux               ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_frootn_to_litter     ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_livecrootn_to_litter ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_deadcrootn_to_litter ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(prod10_nflux             ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(prod100_nflux            ) == (/endp/)), errMsg(__FILE__, __LINE__))
   
    old_weight_was_zero = patch_state_updater%old_weight_was_zero(bounds)
    patch_grew          = patch_state_updater%patch_grew(bounds)

    call ComputeSeedAmounts(bounds                                        , &
         species                    = CN_SPECIES_N                        , &
         leaf_patch                 = ns%leafn_patch(begp:endp)           , &
         leaf_storage_patch         = ns%leafn_storage_patch(begp:endp)   , &
         leaf_xfer_patch            = ns%leafn_xfer_patch(begp:endp)      , &

         ! Calculations only needed for patches that grew:
         compute_here_patch         = patch_grew(begp:endp)               , &

         ! For patches that previously had zero area, ignore the current state for the
         ! sake of computing leaf proportions:
         ignore_current_state_patch = old_weight_was_zero(begp:endp)      , &

         seed_leaf_patch            = seed_leafn_patch(begp:endp)         , &
         seed_leaf_storage_patch    = seed_leafn_storage_patch(begp:endp) , &
         seed_leaf_xfer_patch       = seed_leafn_xfer_patch(begp:endp)    , &
         seed_deadstem_patch        = seed_deadstemn_patch(begp:endp)     , &
         pool_seed_param            = npool_seed_param                    , &
         pool_seed_patch            = seed_npool_patch(begp:endp))

    ! 1) LEAFN_PATCH
    call patch_state_updater%update_patch_state(            &
         bounds                                           , &
         num_filterp_with_inactive                        , &
         filterp_with_inactive                            , &
         var               = ns%leafn_patch   (begp:endp) , &
         flux_out_grc_area = conv_nflux       (begp:endp) , &
         seed              = seed_leafn_patch (begp:endp) , &
         seed_addition     = dwt_leafn_seed   (begp:endp))

    ! 2) LEAFN_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                    &
         bounds                                                   , &
         num_filterp_with_inactive                                , &
         filterp_with_inactive                                    , &
         var               = ns%leafn_storage_patch   (begp:endp) , &
         flux_out_grc_area = conv_nflux               (begp:endp) , &
         seed              = seed_leafn_storage_patch (begp:endp) , &
         seed_addition     = dwt_leafn_seed           (begp:endp))

    ! 3) LEAFN_XFER_PATCH
    call patch_state_updater%update_patch_state( &
         bounds                                                        , &
         num_filterp_with_inactive                                     , &
         filterp_with_inactive                                         , &
         var               = ns%leafn_xfer_patch   (begp:endp), &
         flux_out_grc_area = conv_nflux            (begp:endp), &
         seed              = seed_leafn_xfer_patch (begp:endp), &
         seed_addition     = dwt_leafn_seed        (begp:endp))

    ! 4) FROOTN_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ns%frootn_patch(begp:endp)             , &
         flux_out_grc_area = dwt_frootn_to_litter(begp:endp))

    ! 5) FROOTN_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ns%frootn_storage_patch(begp:endp)     , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 6) FROOTN_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ns%frootn_xfer_patch(begp:endp)        , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 7) LIVESTEMN_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ns%livestemn_patch(begp:endp)          , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 8) LIVESTEMN_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ns%livestemn_storage_patch(begp:endp)  , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 9) LIVESTEMN_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ns%livestemn_xfer_patch(begp:endp)     , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 10) LIVECROOTN_PATCH 
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ns%livecrootn_patch(begp:endp)         , &
         flux_out_grc_area = dwt_frootn_to_litter(begp:endp))

    ! 11) LIVECROOTN_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ns%livecrootn_storage_patch(begp:endp) , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 12) LIVECROOTN_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ns%livecrootn_xfer_patch(begp:endp)    , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 13) DEADCROOTN_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ns%deadcrootn_patch(begp:endp)         , &
         flux_out_grc_area = dwt_frootn_to_litter(begp:endp))

    ! 14) DEADCROOTN_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ns%deadcrootn_storage_patch(begp:endp) , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 15) DEADCROOT_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ns%deadcrootn_xfer_patch(begp:endp)    , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 16) RETRANSN_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ns%retransn_patch(begp:endp)           , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 17) NTRUNC_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ns%ntrunc_patch(begp:endp)             , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 18) CPOOL_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ns%npool_patch(begp:endp)              , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 19) DISPVEGN_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ns%dispvegn_patch(begp:endp))

    ! 20) STORVEGN_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ns%storvegn_patch(begp:endp))

    ! 21) TOTVEGN_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ns%totvegn_patch(begp:endp))

    ! 22) TOTPFTN_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ns%totpftn_patch(begp:endp))

    ! 23) DEADSTEMN_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ns%deadstemn_storage_patch(begp:endp)  , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 24) DEADSTEMN_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ns%deadstemn_xfer_patch(begp:endp)     , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 25) PROD10_NFLUX
    wood_product_nflux(begp:endp)      = 0._r8
    deadstemn_patch_temp(begp:endp)    = ns%deadstemn_patch(begp:endp)
    call patch_state_updater%update_patch_state_partition_flux_by_type(     &
         bounds                                                           , &
         num_filterp_with_inactive                                        , &
         filterp_with_inactive                                            , &
         flux1_fraction_by_pft_type = pprod10                             , &
         var                        = deadstemn_patch_temp    (begp:endp) , &
         flux1_out                  = prod10_nflux            (begp:endp) , &
         flux2_out                  = wood_product_nflux      (begp:endp) , &
         seed                       = seed_deadstemn_patch    (begp:endp) )

    ! 26) PROD100_NFLUX
    wood_product_nflux(begp:endp)      = 0._r8
    deadstemn_patch_temp(begp:endp)    = ns%deadstemn_patch(begp:endp)
    call patch_state_updater%update_patch_state_partition_flux_by_type(     &
         bounds                                                           , &
         num_filterp_with_inactive                                        , &
         filterp_with_inactive                                            , &
         flux1_fraction_by_pft_type = pprod100                            , &
         var                        = deadstemn_patch_temp    (begp:endp) , &
         flux1_out                  = prod100_nflux           (begp:endp) , &
         flux2_out                  = wood_product_nflux      (begp:endp) , &
         seed                       = seed_deadstemn_patch    (begp:endp))

    ! 27) DEADSTEMN_PATCH
    wood_product_nflux(begp:endp)      = 0._r8
    call patch_state_updater%update_patch_state_partition_flux_by_type(     &
         bounds                                                           , &
         num_filterp_with_inactive                                        , &
         filterp_with_inactive                                            , &
         flux1_fraction_by_pft_type = pconv                               , &
         var                        = ns%deadstemn_patch   (begp:endp)    , &
         flux1_out                  = conv_nflux           (begp:endp)    , &
         flux2_out                  = wood_product_nflux   (begp:endp)    , &
         seed                       = seed_deadstemn_patch (begp:endp)    , &
         seed_addition              = dwt_deadstemn_seed   (begp:endp))


    ! These fluxes are computed as negative quantities, but are expected to be positive,
    ! so flip the signs
    do p = begp,endp
       dwt_frootn_to_litter(p)     = -1._r8 * dwt_frootn_to_litter(p)
       dwt_livecrootn_to_litter(p) = -1._r8 * dwt_livecrootn_to_litter(p)
       dwt_deadcrootn_to_litter(p) = -1._r8 * dwt_deadcrootn_to_litter(p)
    end do

  end subroutine NStateDynamicPatchAdjustments
  
  !-----------------------------------------------------------------------
  subroutine PStateDynamicPatchAdjustments( &
       bounds,                              &
       num_filterp_with_inactive,           &
       filterp_with_inactive,               &
       prior_weights,                       &
       patch_state_updater,                 &
       ps,                                  &
       dwt_leafp_seed,                      &
       dwt_deadstemp_seed,                  &
       dwt_ppool_seed,                      &
       conv_pflux,                          &
       dwt_frootp_to_litter,                &
       dwt_livecrootp_to_litter,            &
       dwt_deadcrootp_to_litter,            &
       prod10_pflux,                        &
       prod100_pflux                        &
       )
    !
    ! !DESCRIPTION:
    ! Adjust state variables and compute associated fluxes when patch areas change due to
    ! dynamic landuse
    !
    ! !USES:
    use pftvarcon          , only : pconv, pprod10, pprod100
    use dynPriorWeightsMod , only : prior_weights_type
    use landunit_varcon    , only : istsoil, istcrop
    !
    ! !ARGUMENTS:
    type(bounds_type)              , intent(in)    :: bounds
    integer                        , intent(in)    :: num_filterp_with_inactive
    integer                        , intent(in)    :: filterp_with_inactive(:)
    type(prior_weights_type)       , intent(in)    :: prior_weights
    type(patch_state_updater_type) , intent(in)    :: patch_state_updater
    type(phosphorusstate_type)       , intent(inout) :: ps
    real(r8)                       , intent(inout) :: dwt_leafp_seed           (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_deadstemp_seed       (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_ppool_seed           (bounds%begp:)
    real(r8)                       , intent(inout) :: conv_pflux               (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_frootp_to_litter     (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_livecrootp_to_litter (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_deadcrootp_to_litter (bounds%begp:)
    real(r8)                       , intent(inout) :: prod10_pflux             (bounds%begp:)
    real(r8)                       , intent(inout) :: prod100_pflux            (bounds%begp:)
    !
    ! !LOCAL VARIABLES:
    integer                     :: begp, endp
    integer                     :: l, c, p
    logical                     :: old_weight_was_zero      (bounds%begp:bounds%endp)
    logical                     :: patch_grew               (bounds%begp:bounds%endp)

    ! The following are only set for growing patches:
    real(r8)                    :: seed_leafp_patch         (bounds%begp:bounds%endp)
    real(r8)                    :: seed_leafp_storage_patch (bounds%begp:bounds%endp)
    real(r8)                    :: seed_leafp_xfer_patch    (bounds%begp:bounds%endp)
    real(r8)                    :: seed_deadstemp_patch     (bounds%begp:bounds%endp)
    real(r8)                    :: seed_ppool_patch         (bounds%begp:bounds%endp)

    real(r8)                    :: wood_product_pflux       (bounds%begp:bounds%endp)
    real(r8)                    :: deadstemp_patch_temp     (bounds%begp:bounds%endp)

    character(len=*), parameter :: subname = 'PStateDynamicPatchAdjustments'
    !-----------------------------------------------------------------------

    begp = bounds%begp
    endp = bounds%endp

    SHR_ASSERT_ALL((ubound(dwt_leafp_seed           ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_deadstemp_seed       ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_ppool_seed           ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(conv_pflux               ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_frootp_to_litter     ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_livecrootp_to_litter ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_deadcrootp_to_litter ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(prod10_pflux             ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(prod100_pflux            ) == (/endp/)), errMsg(__FILE__, __LINE__))
   
    old_weight_was_zero = patch_state_updater%old_weight_was_zero(bounds)
    patch_grew          = patch_state_updater%patch_grew(bounds)

    call ComputeSeedAmounts(bounds                                        , &
         species                    = CN_SPECIES_P                        , &
         leaf_patch                 = ps%leafp_patch(begp:endp)           , &
         leaf_storage_patch         = ps%leafp_storage_patch(begp:endp)   , &
         leaf_xfer_patch            = ps%leafp_xfer_patch(begp:endp)      , &

         ! Calculations only needed for patches that grew:
         compute_here_patch         = patch_grew(begp:endp)               , &

         ! For patches that previously had zero area, ignore the current state for the
         ! sake of computing leaf proportions:
         ignore_current_state_patch = old_weight_was_zero(begp:endp)      , &

         seed_leaf_patch            = seed_leafp_patch(begp:endp)         , &
         seed_leaf_storage_patch    = seed_leafp_storage_patch(begp:endp) , &
         seed_leaf_xfer_patch       = seed_leafp_xfer_patch(begp:endp)    , &
         seed_deadstem_patch        = seed_deadstemp_patch(begp:endp)     , &
         pool_seed_param            = ppool_seed_param                    , &
         pool_seed_patch            = seed_ppool_patch(begp:endp))

    ! 1) LEAFP_PATCH
    call patch_state_updater%update_patch_state(            &
         bounds                                           , &
         num_filterp_with_inactive                        , &
         filterp_with_inactive                            , &
         var               = ps%leafp_patch   (begp:endp) , &
         flux_out_grc_area = conv_pflux       (begp:endp) , &
         seed              = seed_leafp_patch (begp:endp) , &
         seed_addition     = dwt_leafp_seed   (begp:endp))

    ! 2) LEAFP_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                    &
         bounds                                                   , &
         num_filterp_with_inactive                                , &
         filterp_with_inactive                                    , &
         var               = ps%leafp_storage_patch   (begp:endp) , &
         flux_out_grc_area = conv_pflux               (begp:endp) , &
         seed              = seed_leafp_storage_patch (begp:endp) , &
         seed_addition     = dwt_leafp_seed           (begp:endp))

    ! 3) LEAFP_XFER_PATCH
    call patch_state_updater%update_patch_state( &
         bounds                                                        , &
         num_filterp_with_inactive                                     , &
         filterp_with_inactive                                         , &
         var               = ps%leafp_xfer_patch   (begp:endp), &
         flux_out_grc_area = conv_pflux            (begp:endp), &
         seed              = seed_leafp_xfer_patch (begp:endp), &
         seed_addition     = dwt_leafp_seed        (begp:endp))

    ! 4) FROOTP_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ps%frootp_patch(begp:endp)             , &
         flux_out_grc_area = dwt_frootp_to_litter(begp:endp))

    ! 5) FROOTP_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ps%frootp_storage_patch(begp:endp)     , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 6) FROOTP_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ps%frootp_xfer_patch(begp:endp)        , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 7) LIVESTEMP_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ps%livestemp_patch(begp:endp)          , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 8) LIVESTEMP_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ps%livestemp_storage_patch(begp:endp)  , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 9) LIVESTEMP_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ps%livestemp_xfer_patch(begp:endp)     , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 10) LIVECROOTP_PATCH 
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ps%livecrootp_patch(begp:endp)         , &
         flux_out_grc_area = dwt_frootp_to_litter(begp:endp))

    ! 11) LIVECROOTP_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ps%livecrootp_storage_patch(begp:endp) , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 12) LIVECROOTP_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ps%livecrootp_xfer_patch(begp:endp)    , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 13) DEADCROOTP_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ps%deadcrootp_patch(begp:endp)         , &
         flux_out_grc_area = dwt_frootp_to_litter(begp:endp))

    ! 14) DEADCROOTP_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ps%deadcrootp_storage_patch(begp:endp) , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 15) DEADCROOT_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ps%deadcrootp_xfer_patch(begp:endp)    , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 16) RETRANSP_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ps%retransp_patch(begp:endp)           , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 17) PTRUNC_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ps%ptrunc_patch(begp:endp)             , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 18) PPOOL_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ps%ppool_patch(begp:endp)              , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 19) DISPVEGP_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ps%dispvegp_patch(begp:endp))

    ! 20) STORVEGP_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ps%storvegp_patch(begp:endp))

    ! 21) TOTVEGP_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ps%totvegp_patch(begp:endp))

    ! 22) TOTPFTP_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ps%totpftp_patch(begp:endp))

    ! 23) DEADSTEMP_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ps%deadstemp_storage_patch(begp:endp)  , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 24) DEADSTEMP_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = ps%deadstemp_xfer_patch(begp:endp)     , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 25) PROD10_PFLUX
    wood_product_pflux(begp:endp)      = 0._r8
    deadstemp_patch_temp(begp:endp)    = ps%deadstemp_patch(begp:endp)
    call patch_state_updater%update_patch_state_partition_flux_by_type(     &
         bounds                                                           , &
         num_filterp_with_inactive                                        , &
         filterp_with_inactive                                            , &
         flux1_fraction_by_pft_type = pprod10                             , &
         var                        = deadstemp_patch_temp    (begp:endp) , &
         flux1_out                  = prod10_pflux            (begp:endp) , &
         flux2_out                  = wood_product_pflux      (begp:endp) , &
         seed                       = seed_deadstemp_patch    (begp:endp) )

    ! 26) PROD100_PFLUX
    wood_product_pflux(begp:endp)      = 0._r8
    deadstemp_patch_temp(begp:endp)    = ps%deadstemp_patch(begp:endp)
    call patch_state_updater%update_patch_state_partition_flux_by_type(     &
         bounds                                                           , &
         num_filterp_with_inactive                                        , &
         filterp_with_inactive                                            , &
         flux1_fraction_by_pft_type = pprod100                            , &
         var                        = deadstemp_patch_temp    (begp:endp) , &
         flux1_out                  = prod100_pflux           (begp:endp) , &
         flux2_out                  = wood_product_pflux      (begp:endp) , &
         seed                       = seed_deadstemp_patch    (begp:endp))

    ! 27) DEADSTEMP_PATCH
    wood_product_pflux(begp:endp)      = 0._r8
    call patch_state_updater%update_patch_state_partition_flux_by_type(     &
         bounds                                                           , &
         num_filterp_with_inactive                                        , &
         filterp_with_inactive                                            , &
         flux1_fraction_by_pft_type = pconv                               , &
         var                        = ps%deadstemp_patch   (begp:endp)    , &
         flux1_out                  = conv_pflux           (begp:endp)    , &
         flux2_out                  = wood_product_pflux   (begp:endp)    , &
         seed                       = seed_deadstemp_patch (begp:endp)    , &
         seed_addition              = dwt_deadstemp_seed   (begp:endp))

    ! These fluxes are computed as negative quantities, but are expected to be positive,
    ! so flip the signs
    do p = begp,endp
       dwt_frootp_to_litter(p)     = -1._r8 * dwt_frootp_to_litter(p)
       dwt_livecrootp_to_litter(p) = -1._r8 * dwt_livecrootp_to_litter(p)
       dwt_deadcrootp_to_litter(p) = -1._r8 * dwt_deadcrootp_to_litter(p)
    end do

  end subroutine PStateDynamicPatchAdjustments

end module dynConsBiogeochemMod
