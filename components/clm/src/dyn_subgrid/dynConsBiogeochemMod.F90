module dynConsBiogeochemMod

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
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  save
  
  public :: dyn_cnbal_patch
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine dyn_cnbal_patch(bounds, prior_weights, &
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
    type(bounds_type)        , intent(in)    :: bounds        
    type(prior_weights_type) , intent(in)    :: prior_weights ! weights prior to the subgrid weight updates
    type(canopystate_type)   , intent(inout) :: canopystate_vars
    type(photosyns_type)     , intent(inout) :: photosyns_vars
    type(cnstate_type)       , intent(inout) :: cnstate_vars
    type(carbonstate_type)   , intent(inout) :: carbonstate_vars
    type(carbonstate_type)   , intent(inout) :: c13_carbonstate_vars
    type(carbonstate_type)   , intent(inout) :: c14_carbonstate_vars
    type(carbonflux_type)    , intent(inout) :: carbonflux_vars
    type(carbonflux_type)    , intent(inout) :: c13_carbonflux_vars
    type(carbonflux_type)    , intent(inout) :: c14_carbonflux_vars
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars

    type(phosphorusstate_type) , intent(inout) :: phosphorusstate_vars
    type(phosphorusflux_type)  , intent(inout) :: phosphorusflux_vars
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
       dwt_leafn_seed(p)           = 0._r8
       dwt_deadstemc_seed(p)       = 0._r8
       dwt_deadstemn_seed(p)       = 0._r8
       dwt_npool_seed(p)           = 0._r8
       dwt_frootc_to_litter(p)     = 0._r8
       dwt_livecrootc_to_litter(p) = 0._r8
       dwt_deadcrootc_to_litter(p) = 0._r8
       dwt_frootn_to_litter(p)     = 0._r8
       dwt_livecrootn_to_litter(p) = 0._r8
       dwt_deadcrootn_to_litter(p) = 0._r8
       conv_cflux(p)               = 0._r8
       prod10_cflux(p)             = 0._r8
       prod100_cflux(p)            = 0._r8
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
             
             ! set the seed sources for leaf and deadstem
             ! leaf source is split later between leaf, leaf_storage, leaf_xfer
             call ComputeLeafSeedPools(veg_pp%itype(p), &
                  leafc_seed, leafn_seed, leafp_seed, leafc13_seed, leafc14_seed, &
                  npool_seed, ppool_seed)
             
             call ComputeStemSeedPools(veg_pp%itype(p), &
                  deadstemc_seed, deadstemn_seed, deadstemp_seed, &
                  deadstemc13_seed, deadstemc14_seed)
             
             call UpdateCStateVarsDueToWtInc(cs, p, &
                  veg_pp%wtcol(p), prior_weights%pwtcol(p),     &
                  leafc_seed, deadstemc_seed)
             
             if ( use_c13 ) then
                call UpdateCStateVarsDueToWtInc(c13_cs, p, &
                     veg_pp%wtcol(p), prior_weights%pwtcol(p),         &
                  leafc13_seed, deadstemc13_seed)
             endif
             
             if ( use_c14 ) then
                call UpdateCStateVarsDueToWtInc(c14_cs, p, &
                     veg_pp%wtcol(p), prior_weights%pwtcol(p),         &
                     leafc14_seed, deadstemc14_seed)
             endif

             call UpdateNStateVarsDueToWtInc(ns, p, &
                  veg_pp%wtcol(p), prior_weights%pwtcol(p),       &
                  leafn_seed, deadstemn_seed, npool_seed)
             
             call UpdatePStateVarsDueToWtInc(ps, p, &
                  veg_pp%wtcol(p), prior_weights%pwtcol(p),         &
                  leafp_seed, deadstemp_seed, ppool_seed)

             ! update temporary seed source arrays
             ! These are calculated in terms of the required contributions from
             ! column-level seed source
             dwt_leafc_seed(p)     = leafc_seed     * dwt

             dwt_leafn_seed(p)     = leafn_seed     * dwt
             dwt_deadstemc_seed(p) = deadstemc_seed * dwt
             dwt_deadstemn_seed(p) = deadstemn_seed * dwt
             dwt_npool_seed(p)     = npool_seed     * dwt             

             dwt_leafp_seed(p)     = leafp_seed     * dwt
             dwt_deadstemp_seed(p) = deadstemp_seed * dwt
             dwt_ppool_seed(p)     = ppool_seed     * dwt

             if ( use_c13 ) then
                dwt_leafc13_seed(p)     = leafc13_seed     * dwt
                dwt_deadstemc13_seed(p) = deadstemc13_seed * dwt
             endif

             if ( use_c14 ) then
                dwt_leafc14_seed(p)     = leafc14_seed     * dwt
                dwt_deadstemc14_seed(p) = deadstemc14_seed * dwt
             endif

          else if (dwt < 0._r8) then
             
             ! if the pft lost weight on the timestep, then the carbon and nitrogen state
             ! variables are directed to litter, CWD, and wood product pools.
             
             ! N.B. : the conv_cflux, prod10_cflux, and prod100_cflux fluxes are accumulated
             ! as negative values, but the fluxes for pft-to-litter are accumulated as 
             ! positive values
             
             ! set local weight variables for this pft
             wt_new = veg_pp%wtcol(p)
             wt_old = prior_weights%pwtcol(p)
             
             call UpdateCStateDueToWtDec( &
                  cs, p, wt_old, wt_new,         &
                  conv_cflux(p),                 &
                  dwt_frootc_to_litter(p),       &
                  dwt_livecrootc_to_litter(p),   &
                  dwt_deadcrootc_to_litter(p),   &
                  prod10_cflux(p),               &
                  prod100_cflux(p),              &
                  pprod10(veg_pp%itype(p)),      &
                  pprod100(veg_pp%itype(p)))

             if ( use_c13 ) then
                call UpdateCStateDueToWtDec( &
                     c13_cs, p, wt_old, wt_new,     &
                     conv_c13flux(p),               &
                     dwt_frootc13_to_litter(p),     &
                     dwt_livecrootc13_to_litter(p), &
                     dwt_deadcrootc13_to_litter(p), &
                     prod10_c13flux(p),             &
                     prod100_c13flux(p),            &
                     pprod10(veg_pp%itype(p)),      &
                     pprod100(veg_pp%itype(p)))
             endif

             if ( use_c14 ) then
                call UpdateCStateDueToWtDec( &
                     c14_cs, p, wt_old, wt_new,     &
                     conv_c14flux(p),               &
                     dwt_frootc14_to_litter(p),     &
                     dwt_livecrootc14_to_litter(p), &
                     dwt_deadcrootc14_to_litter(p), &
                     prod10_c14flux(p),             &
                     prod100_c14flux(p),            &
                     pprod10(veg_pp%itype(p)),      &
                     pprod100(veg_pp%itype(p)))

             endif
             
             call UpdateNStateDueToWtDec( &
                  ns, p, wt_old, wt_new,           &
                  conv_nflux(p),                   &
                  dwt_frootn_to_litter(p),         &
                  dwt_livecrootn_to_litter(p),     &
                  dwt_deadcrootn_to_litter(p),     &
                  prod10_nflux(p),                 &
                  prod100_nflux(p),                &
                  pprod10(veg_pp%itype(p)),        &
                  pprod100(veg_pp%itype(p)))

             call UpdatePStateDueToWtDec( &
                  ps, p, wt_old, wt_new,             &
                  conv_pflux(p),                     &
                  dwt_frootn_to_litter(p),           &
                  dwt_livecrootp_to_litter(p),       &
                  dwt_deadcrootp_to_litter(p),       &
                  prod10_pflux(p),                   &
                  prod100_pflux(p),                  &
                  pprod10(veg_pp%itype(p)),          &
                  pprod100(veg_pp%itype(p)))
                                                    
          end if       ! weight decreasing
       end if           ! is soil
    end do               ! patch loop
    
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
 subroutine ComputeLeafSeedPools(veg_type, &
     leafc_seed, leafn_seed, leafp_seed, leafc13_seed, leafc14_seed, &
      npool_seed, ppool_seed)
   !
   ! !DESCRIPTION:
   ! Computes values of seed for various leaf pool
   !
   use shr_const_mod      , only : SHR_CONST_PDB
   use clm_varcon         , only : c14ratio, c3_r2, c4_r2
   !
   implicit none
   !
   ! !ARGUMENT
   integer               , intent(in)  :: veg_type
   real(r8)              , intent(out) :: leafc_seed
   real(r8)              , intent(out) :: leafn_seed
   real(r8)              , intent(out) :: leafp_seed
   real(r8)              , intent(out) :: leafc13_seed
   real(r8)              , intent(out) :: leafc14_seed
   real(r8)              , intent(out) :: npool_seed
   real(r8)              , intent(out) :: ppool_seed
   !
   ! !LOCAL
   real(r8)              , parameter :: leafc_seed_param = 1._r8
   real(r8)              , parameter :: npool_seed_param = 0.1_r8
   real(r8)              , parameter :: ppool_seed_param = 0.01_r8

   leafc_seed     = 0._r8
   leafn_seed     = 0._r8
   leafp_seed     = 0._r8
   npool_seed     = 0._r8
   ppool_seed     = 0._r8

   if (veg_type /= 0) then
      leafc_seed = leafc_seed_param
      leafn_seed = leafc_seed / veg_vp%leafcn(veg_type)
      leafp_seed = leafc_seed / veg_vp%leafcp(veg_type)

      if (veg_vp%nstor(veg_type) > 1e-6_r8) then
         npool_seed     = npool_seed_param
         ppool_seed     = ppool_seed_param
      end if

      if ( use_c13 ) then
         if (veg_vp%c3psn(veg_type) == 1._r8) then
            leafc13_seed     = leafc_seed     * c3_r2
         else
            leafc13_seed     = leafc_seed     * c4_r2
         end if
      endif

      if ( use_c14 ) then
         ! 14c state is initialized assuming initial "modern" 14C of 1.e-12
         if (veg_vp%c3psn(veg_type) == 1._r8) then
            leafc14_seed     = leafc_seed     * c14ratio
         else
            leafc14_seed     = leafc_seed     * c14ratio
         end if
      endif
   end if

 end subroutine ComputeLeafSeedPools

 !-----------------------------------------------------------------------
 subroutine ComputeStemSeedPools(veg_type, &
      deadstemc_seed, deadstemn_seed, deadstemp_seed, &
      deadstemc13_seed, deadstemc14_seed)
   !
   ! !DESCRIPTION:
   ! Computes values of seed for various stem pool
   !
   use shr_const_mod      , only : SHR_CONST_PDB
   use clm_varcon         , only : c14ratio, c3_r2, c4_r2
   !
   implicit none
   !
   ! !ARGUMENT
   integer               , intent(in)  :: veg_type
   real(r8)              , intent(out) :: deadstemc_seed
   real(r8)              , intent(out) :: deadstemn_seed
   real(r8)              , intent(out) :: deadstemp_seed
   real(r8)              , intent(out) :: deadstemc13_seed
   real(r8)              , intent(out) :: deadstemc14_seed
   !
   ! !LOCAL
   real(r8)              , parameter :: deadstemc_seed_param = 0.1_r8
  
   deadstemc_seed = 0._r8
   deadstemn_seed = 0._r8
   deadstemp_seed = 0._r8

   if (veg_type /= 0) then
      if (veg_vp%woody(veg_type) == 1._r8) then
         deadstemc_seed = deadstemc_seed_param
         deadstemn_seed = deadstemc_seed / veg_vp%deadwdcn(veg_type)
         deadstemp_seed = deadstemc_seed / veg_vp%deadwdcp(veg_type)
      end if

      if ( use_c13 ) then
         if (veg_vp%c3psn(veg_type) == 1._r8) then
            deadstemc13_seed = deadstemc_seed * c3_r2
         else
            deadstemc13_seed = deadstemc_seed * c4_r2
         end if
      endif

      if ( use_c14 ) then
         ! 14c state is initialized assuming initial "modern" 14C of 1.e-12
         if (veg_vp%c3psn(veg_type) == 1._r8) then
            deadstemc14_seed = deadstemc_seed * c14ratio
         else
            deadstemc14_seed = deadstemc_seed * c14ratio
         end if
      endif
   end if

 end subroutine ComputeStemSeedPools

 !-----------------------------------------------------------------------
 subroutine UpdateCStateVarsDueToWtInc(cs, p, wt_new, wt_old, &
      leafc_seed, deadstemc_seed)
   !
   ! !DESCRIPTION:
   ! Updates p-th patch of carbonstate_type
   ! When PFT area expands (dwt > 0), the pft-level mass density 
   ! is modified to conserve the original pft mass distributed
   ! over the new (larger) area, plus a term to account for the
   ! introduction of new seed source for leaf and deadstem
   !
   implicit none
   !
   ! !ARGUMENT
   type(carbonstate_type), intent(inout) :: cs
   integer               , intent(in)    :: p
   real(r8)              , intent(in)    :: wt_new
   real(r8)              , intent(in)    :: wt_old
   real(r8)              , intent(in)    :: leafc_seed
   real(r8)              , intent(in)    :: deadstemc_seed
   !
   ! !LOCAL
   real(r8)                              :: pleaf
   real(r8)                              :: pstor
   real(r8)                              :: pxfer
   real(r8)                              :: t1
   real(r8)                              :: t2

   t1 = wt_old           /wt_new
   t2 = (wt_new - wt_old)/wt_old
   
   call LeafProportions(veg_pp%itype(p), &
       cs%leafc_patch(p), cs%leafc_storage_patch(p), cs%leafc_xfer_patch(p), &
       pleaf, pstor, pxfer)

   cs%leafc_patch(p)              = cs%leafc_patch(p)              * t1 + leafc_seed*pleaf*t2
   cs%leafc_storage_patch(p)      = cs%leafc_storage_patch(p)      * t1 + leafc_seed*pstor*t2
   cs%leafc_xfer_patch(p)         = cs%leafc_xfer_patch(p)         * t1 + leafc_seed*pxfer*t2
   cs%frootc_patch(p)             = cs%frootc_patch(p)             * t1
   cs%frootc_storage_patch(p)     = cs%frootc_storage_patch(p)     * t1
   cs%frootc_xfer_patch(p)        = cs%frootc_xfer_patch(p)        * t1
   cs%livestemc_patch(p)          = cs%livestemc_patch(p)          * t1
   cs%livestemc_storage_patch(p)  = cs%livestemc_storage_patch(p)  * t1
   cs%livestemc_xfer_patch(p)     = cs%livestemc_xfer_patch(p)     * t1
   cs%deadstemc_patch(p)          = cs%deadstemc_patch(p)          * t1  + deadstemc_seed*t2
   cs%deadstemc_storage_patch(p)  = cs%deadstemc_storage_patch(p)  * t1
   cs%deadstemc_xfer_patch(p)     = cs%deadstemc_xfer_patch(p)     * t1
   cs%livecrootc_patch(p)         = cs%livecrootc_patch(p)         * t1
   cs%livecrootc_storage_patch(p) = cs%livecrootc_storage_patch(p) * t1
   cs%livecrootc_xfer_patch(p)    = cs%livecrootc_xfer_patch(p)    * t1
   cs%deadcrootc_patch(p)         = cs%deadcrootc_patch(p)         * t1
   cs%deadcrootc_storage_patch(p) = cs%deadcrootc_storage_patch(p) * t1
   cs%deadcrootc_xfer_patch(p)    = cs%deadcrootc_xfer_patch(p)    * t1
   cs%gresp_storage_patch(p)      = cs%gresp_storage_patch(p)      * t1
   cs%gresp_xfer_patch(p)         = cs%gresp_xfer_patch(p)         * t1
   cs%cpool_patch(p)              = cs%cpool_patch(p)              * t1
   cs%xsmrpool_patch(p)           = cs%xsmrpool_patch(p)           * t1
   cs%ctrunc_patch(p)             = cs%ctrunc_patch(p)             * t1
   cs%dispvegc_patch(p)           = cs%dispvegc_patch(p)           * t1
   cs%storvegc_patch(p)           = cs%storvegc_patch(p)           * t1
   cs%totvegc_patch(p)            = cs%totvegc_patch(p)            * t1
   cs%totpftc_patch(p)            = cs%totpftc_patch(p)            * t1

 end subroutine UpdateCStateVarsDueToWtInc

 !-----------------------------------------------------------------------
 subroutine UpdateNStateVarsDueToWtInc(ns, p, wt_new, wt_old, &
      leafn_seed, deadstemn_seed, npool_seed)
   !
   ! !DESCRIPTION:
   ! Updates p-th patch of nitrogenstate_type
   ! When PFT area expands (dwt > 0), the pft-level mass density
   ! is modified to conserve the original pft mass distributed
   ! over the new (larger) area, plus a term to account for the
   ! introduction of new seed source for leaf and deadstem
   !
   implicit none
   !
   ! !ARGUMENT
   type(nitrogenstate_type) , intent(inout) :: ns
   integer                  , intent(in)    :: p
   real(r8)                 , intent(in)    :: wt_new
   real(r8)                 , intent(in)    :: wt_old
   real(r8)                 , intent(in)    :: leafn_seed
   real(r8)                 , intent(in)    :: deadstemn_seed
   real(r8)                 , intent(in)    :: npool_seed
   !
   ! !LOCAL
   real(r8)                                 :: pleaf
   real(r8)                                 :: pstor
   real(r8)                                 :: pxfer
   real(r8)                                 :: t1
   real(r8)                                 :: t2

   t1 = wt_old           /wt_new
   t2 = (wt_new - wt_old)/wt_old
   
   call LeafProportions(veg_pp%itype(p), &
       ns%leafn_patch(p), ns%leafn_storage_patch(p), ns%leafn_xfer_patch(p), &
       pleaf, pstor, pxfer)

   ns%leafn_patch(p)              = ns%leafn_patch(p)              * t1 + leafn_seed*pleaf*t2
   ns%leafn_storage_patch(p)      = ns%leafn_storage_patch(p)      * t1 + leafn_seed*pstor*t2
   ns%leafn_xfer_patch(p)         = ns%leafn_xfer_patch(p)         * t1 + leafn_seed*pxfer*t2
   ns%frootn_patch(p)             = ns%frootn_patch(p)             * t1
   ns%frootn_storage_patch(p)     = ns%frootn_storage_patch(p)     * t1
   ns%frootn_xfer_patch(p)        = ns%frootn_xfer_patch(p)        * t1
   ns%livestemn_patch(p)	  = ns%livestemn_patch(p)          * t1
   ns%livestemn_storage_patch(p)  = ns%livestemn_storage_patch(p)  * t1
   ns%livestemn_xfer_patch(p)     = ns%livestemn_xfer_patch(p)     * t1
   ns%deadstemn_patch(p)          = ns%deadstemn_patch(p)          * t1 + deadstemn_seed*t2
   ns%deadstemn_storage_patch(p)  = ns%deadstemn_storage_patch(p)  * t1
   ns%deadstemn_xfer_patch(p)     = ns%deadstemn_xfer_patch(p)     * t1
   ns%npool_patch(p)              = ns%npool_patch(p)              * t1 + npool_seed*t2
   ns%livecrootn_patch(p)         = ns%livecrootn_patch(p)         * t1
   ns%livecrootn_storage_patch(p) = ns%livecrootn_storage_patch(p) * t1
   ns%livecrootn_xfer_patch(p)    = ns%livecrootn_xfer_patch(p)    * t1
   ns%deadcrootn_patch(p)         = ns%deadcrootn_patch(p)         * t1
   ns%deadcrootn_storage_patch(p) = ns%deadcrootn_storage_patch(p) * t1
   ns%deadcrootn_xfer_patch(p)    = ns%deadcrootn_xfer_patch(p)    * t1
   ns%retransn_patch(p)           = ns%retransn_patch(p)           * t1
   ns%ntrunc_patch(p)             = ns%ntrunc_patch(p)             * t1
   ns%dispvegn_patch(p)           = ns%dispvegn_patch(p)           * t1
   ns%storvegn_patch(p)           = ns%storvegn_patch(p)           * t1
   ns%totvegn_patch(p)            = ns%totvegn_patch(p)            * t1
   ns%totpftn_patch(p)            = ns%totpftn_patch(p)            * t1

 end subroutine UpdateNStateVarsDueToWtInc

 !-----------------------------------------------------------------------
 subroutine UpdatePStateVarsDueToWtInc(ps, p, wt_new, wt_old, &
      leafp_seed, deadstemp_seed, ppool_seed)
   !
   ! !DESCRIPTION:
   ! Updates p-th patch of phosphorusstate_type
   ! When PFT area expands (dwt > 0), the pft-level mass density
   ! is modified to conserve the original pft mass distributed
   ! over the new (larger) area, plus a term to account for the
   ! introduction of new seed source for leaf and deadstem
   !
   implicit none
   !
   ! !ARGUMENT
   type(phosphorusstate_type) , intent(inout) :: ps
   integer                    , intent(in)    :: p
   real(r8)                   , intent(in)    :: wt_new
   real(r8)                   , intent(in)    :: wt_old
   real(r8)                   , intent(in)    :: leafp_seed
   real(r8)                   , intent(in)    :: deadstemp_seed
   real(r8)                   , intent(in)    :: ppool_seed
   !
   ! !LOCAL
   real(r8)                                   :: pleaf
   real(r8)                                   :: pstor
   real(r8)                                   :: pxfer
   real(r8)                                   :: t1
   real(r8)                                   :: t2

   t1 = wt_old           /wt_new
   t2 = (wt_new - wt_old)/wt_old
   
   call LeafProportions(veg_pp%itype(p), &
       ps%leafp_patch(p), ps%leafp_storage_patch(p), ps%leafp_xfer_patch(p), &
       pleaf, pstor, pxfer)

   ps%leafp_patch(p)              = ps%leafp_patch(p)              * t1 + leafp_seed*pleaf*t2
   ps%leafp_storage_patch(p)      = ps%leafp_storage_patch(p)      * t1 + leafp_seed*pstor*t2
   ps%leafp_xfer_patch(p)         = ps%leafp_xfer_patch(p)         * t1 + leafp_seed*pxfer*t2
   ps%frootp_patch(p)             = ps%frootp_patch(p)             * t1
   ps%frootp_storage_patch(p)     = ps%frootp_storage_patch(p)     * t1
   ps%frootp_xfer_patch(p)        = ps%frootp_xfer_patch(p)        * t1
   ps%livestemp_patch(p)	  = ps%livestemp_patch(p)          * t1
   ps%livestemp_storage_patch(p)  = ps%livestemp_storage_patch(p)  * t1
   ps%livestemp_xfer_patch(p)     = ps%livestemp_xfer_patch(p)     * t1
   ps%deadstemp_patch(p)          = ps%deadstemp_patch(p)          * t1 + deadstemp_seed*t2
   ps%deadstemp_storage_patch(p)  = ps%deadstemp_storage_patch(p)  * t1
   ps%deadstemp_xfer_patch(p)     = ps%deadstemp_xfer_patch(p)     * t1
   ps%livecrootp_patch(p)         = ps%livecrootp_patch(p)         * t1
   ps%livecrootp_storage_patch(p) = ps%livecrootp_storage_patch(p) * t1
   ps%livecrootp_xfer_patch(p)    = ps%livecrootp_xfer_patch(p)    * t1
   ps%deadcrootp_patch(p)         = ps%deadcrootp_patch(p)         * t1
   ps%deadcrootp_storage_patch(p) = ps%deadcrootp_storage_patch(p) * t1
   ps%deadcrootp_xfer_patch(p)    = ps%deadcrootp_xfer_patch(p)    * t1
   ps%retransp_patch(p)           = ps%retransp_patch(p)           * t1
   ps%ppool_patch(p)              = ps%ppool_patch(p)              * t1 + ppool_seed*t2
   ps%ptrunc_patch(p)             = ps%ptrunc_patch(p)             * t1
   ps%dispvegp_patch(p)           = ps%dispvegp_patch(p)           * t1
   ps%storvegp_patch(p)           = ps%storvegp_patch(p)           * t1
   ps%totvegp_patch(p)            = ps%totvegp_patch(p)            * t1
   ps%totpftp_patch(p)            = ps%totpftp_patch(p)            * t1

 end subroutine UpdatePStateVarsDueToWtInc

  !-----------------------------------------------------------------------
 subroutine UpdateCStateDueToWtDec(cs, p, wt_old, wt_new, conv_flux, &
      dwt_froot_to_litter, dwt_livecroot_to_litter, dwt_deadcrootc_to_litter, &
      prod10_flux, prod100_flux, prod10, prod100)
   !
   ! !DESCRIPTION:
   ! Updates p-th patch of carbonstate_type due to decrease in fractional cover
   !
   implicit none
   !
   ! !ARGUMENT
   type(carbonstate_type), intent(inout)  :: cs
   integer               , intent(in)     :: p
   real(r8)              , intent(in)     :: wt_old
   real(r8)              , intent(in)     :: wt_new
   real(r8)              , intent (inout) :: conv_flux
   real(r8)              , intent (inout) :: dwt_froot_to_litter
   real(r8)              , intent (inout) :: dwt_livecroot_to_litter
   real(r8)              , intent (inout) :: dwt_deadcrootc_to_litter
   real(r8)              , intent (inout) :: prod10_flux
   real(r8)              , intent (inout) :: prod100_flux
   real(r8)              , intent (in)    :: prod10
   real(r8)              , intent (in)    :: prod100
   !
   ! !LOCAL
   real(r8)                               :: mass_tmp

   call ComputeMassLossDueToWtDec(wt_old, wt_new, cs%leafc_patch(p)              , conv_flux)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, cs%leafc_storage_patch(p)      , conv_flux)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, cs%leafc_xfer_patch(p)         , conv_flux)

   call ComputeMassLossDueToWtDec(wt_old, wt_new, cs%frootc_patch(p)             , dwt_froot_to_litter)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, cs%frootc_storage_patch(p)     , conv_flux)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, cs%frootc_xfer_patch(p)        , conv_flux)

   call ComputeMassLossDueToWtDec(wt_old, wt_new, cs%livestemc_patch(p)          , conv_flux)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, cs%livestemc_storage_patch(p)  , conv_flux)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, cs%livestemc_xfer_patch(p)     , conv_flux)

   call ComputeMassLossDueToWtDec(wt_old, wt_new, cs%deadstemc_storage_patch(p)  , conv_flux)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, cs%deadstemc_xfer_patch(p)     , conv_flux)

   call ComputeMassLossDueToWtDec(wt_old, wt_new, cs%livecrootc_patch(p)         , dwt_livecroot_to_litter)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, cs%livecrootc_storage_patch(p) , conv_flux)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, cs%livecrootc_xfer_patch(p)    , conv_flux)

   call ComputeMassLossDueToWtDec(wt_old, wt_new, cs%deadcrootc_patch(p)         , dwt_deadcrootc_to_litter)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, cs%deadcrootc_storage_patch(p) , conv_flux)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, cs%deadcrootc_xfer_patch(p)    , conv_flux)

   call ComputeMassLossDueToWtDec(wt_old, wt_new, cs%gresp_storage_patch(p)      , conv_flux)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, cs%gresp_xfer_patch(p)         , conv_flux)

   call ComputeMassLossDueToWtDec(wt_old, wt_new, cs%cpool_patch(p)              , conv_flux)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, cs%xsmrpool_patch(p)           , conv_flux)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, cs%ctrunc_patch(p)             , conv_flux)

   ! deadstemc
   mass_tmp = cs%deadstemc_patch(p)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, mass_tmp, prod10_flux, prod10)
   mass_tmp = cs%deadstemc_patch(p)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, mass_tmp, prod100_flux, prod100)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, cs%deadstemc_patch(p), conv_flux)

 end subroutine UpdateCStateDueToWtDec
 
   !-----------------------------------------------------------------------
 subroutine UpdateNStateDueToWtDec(ns, p, wt_old, wt_new, conv_flux, &
      dwt_froot_to_litter, dwt_livecroot_to_litter, dwt_deadcrootc_to_litter, &
      prod10_flux, prod100_flux, prod10, prod100)
   !
   ! !DESCRIPTION:
   ! Updates p-th patch of nitrogenstate_type due to decrease in fractional cover
   !
   implicit none
   !
   ! !ARGUMENT
   type(nitrogenstate_type) , intent (inout) :: ns
   integer                  , intent (in)    :: p
   real(r8)                 , intent (in)    :: wt_old
   real(r8)                 , intent (in)    :: wt_new
   real(r8)                 , intent (inout) :: conv_flux
   real(r8)                 , intent (inout) :: dwt_froot_to_litter
   real(r8)                 , intent (inout) :: dwt_livecroot_to_litter
   real(r8)                 , intent (inout) :: dwt_deadcrootc_to_litter
   real(r8)                 , intent (inout) :: prod10_flux
   real(r8)                 , intent (inout) :: prod100_flux
   real(r8)                 , intent (in)    :: prod10
   real(r8)                 , intent (in)    :: prod100
   !
   ! !LOCAL
   real(r8)                                  :: mass_tmp

   call ComputeMassLossDueToWtDec(wt_old, wt_new, ns%leafn_patch(p)              , conv_flux)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, ns%leafn_storage_patch(p)      , conv_flux)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, ns%leafn_xfer_patch(p)         , conv_flux)

   call ComputeMassLossDueToWtDec(wt_old, wt_new, ns%frootn_patch(p)             , dwt_froot_to_litter)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, ns%frootn_storage_patch(p)     , conv_flux)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, ns%frootn_xfer_patch(p)        , conv_flux)

   call ComputeMassLossDueToWtDec(wt_old, wt_new, ns%livestemn_patch(p)          , conv_flux)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, ns%livestemn_storage_patch(p)  , conv_flux)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, ns%livestemn_xfer_patch(p)     , conv_flux)

   call ComputeMassLossDueToWtDec(wt_old, wt_new, ns%deadstemn_storage_patch(p)  , conv_flux)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, ns%deadstemn_xfer_patch(p)     , conv_flux)

   call ComputeMassLossDueToWtDec(wt_old, wt_new, ns%livecrootn_patch(p)         , dwt_livecroot_to_litter)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, ns%livecrootn_storage_patch(p) , conv_flux)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, ns%livecrootn_xfer_patch(p)    , conv_flux)

   call ComputeMassLossDueToWtDec(wt_old, wt_new, ns%deadcrootn_patch(p)         , dwt_deadcrootc_to_litter)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, ns%deadcrootn_storage_patch(p) , conv_flux)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, ns%deadcrootn_xfer_patch(p)    , conv_flux)

   call ComputeMassLossDueToWtDec(wt_old, wt_new, ns%retransn_patch(p)           , conv_flux)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, ns%npool_patch(p)              , conv_flux)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, ns%ntrunc_patch(p)             , conv_flux)

   ! deadstemc
   mass_tmp = ns%deadstemn_patch(p)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, mass_tmp, prod10_flux, prod10)
   mass_tmp = ns%deadstemn_patch(p)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, mass_tmp, prod100_flux, prod100)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, ns%deadstemn_patch(p), conv_flux)

 end subroutine UpdateNStateDueToWtDec

   !-----------------------------------------------------------------------
 subroutine UpdatePStateDueToWtDec(ps, p, wt_old, wt_new, conv_flux, &
      dwt_froot_to_litter, dwt_livecroot_to_litter, dwt_deadcrootc_to_litter, &
      prod10_flux, prod100_flux, prod10, prod100)
   !
   ! !DESCRIPTION:
   ! Updates p-th patch of nitrogenstate_type due to decrease in fractional cover
   !
   implicit none
   !
   ! !ARGUMENT
   type(phosphorusstate_type) , intent (inout) :: ps
   integer                    , intent (in)    :: p
   real(r8)                   , intent (in)    :: wt_old
   real(r8)                   , intent (in)    :: wt_new
   real(r8)                   , intent (inout) :: conv_flux
   real(r8)                   , intent (inout) :: dwt_froot_to_litter
   real(r8)                   , intent (inout) :: dwt_livecroot_to_litter
   real(r8)                   , intent (inout) :: dwt_deadcrootc_to_litter
   real(r8)                   , intent (inout) :: prod10_flux
   real(r8)                   , intent (inout) :: prod100_flux
   real(r8)                   , intent (in)    :: prod10
   real(r8)                   , intent (in)    :: prod100
   !
   ! !LOCAL
   real(r8)                                    :: mass_tmp

   call ComputeMassLossDueToWtDec(wt_old, wt_new, ps%leafp_patch(p)              , conv_flux)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, ps%leafp_storage_patch(p)      , conv_flux)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, ps%leafp_xfer_patch(p)         , conv_flux)

   call ComputeMassLossDueToWtDec(wt_old, wt_new, ps%frootp_patch(p)             , dwt_froot_to_litter)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, ps%frootp_storage_patch(p)     , conv_flux)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, ps%frootp_xfer_patch(p)        , conv_flux)

   call ComputeMassLossDueToWtDec(wt_old, wt_new, ps%livestemp_patch(p)          , conv_flux)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, ps%livestemp_storage_patch(p)  , conv_flux)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, ps%livestemp_xfer_patch(p)     , conv_flux)

   call ComputeMassLossDueToWtDec(wt_old, wt_new, ps%deadstemp_storage_patch(p)  , conv_flux)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, ps%deadstemp_xfer_patch(p)     , conv_flux)

   call ComputeMassLossDueToWtDec(wt_old, wt_new, ps%livecrootp_patch(p)         , dwt_livecroot_to_litter)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, ps%livecrootp_storage_patch(p) , conv_flux)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, ps%livecrootp_xfer_patch(p)    , conv_flux)

   call ComputeMassLossDueToWtDec(wt_old, wt_new, ps%deadcrootp_patch(p)         , dwt_deadcrootc_to_litter)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, ps%deadcrootp_storage_patch(p) , conv_flux)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, ps%deadcrootp_xfer_patch(p)    , conv_flux)

   call ComputeMassLossDueToWtDec(wt_old, wt_new, ps%retransp_patch(p)           , conv_flux)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, ps%ppool_patch(p)              , conv_flux)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, ps%ptrunc_patch(p)             , conv_flux)

   ! deadstemc
   mass_tmp = ps%deadstemp_patch(p)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, mass_tmp, prod10_flux, prod10)
   mass_tmp = ps%deadstemp_patch(p)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, mass_tmp, prod100_flux, prod100)
   call ComputeMassLossDueToWtDec(wt_old, wt_new, ps%deadstemp_patch(p), conv_flux)

 end subroutine UpdatePStateDueToWtDec
 
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
  subroutine LeafProportions(pft_type, &
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

    if (tot_leaf == 0._r8) then
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

end module dynConsBiogeochemMod
