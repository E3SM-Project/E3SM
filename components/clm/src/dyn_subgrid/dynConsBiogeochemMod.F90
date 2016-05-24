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
  use EcophysConType      , only : ecophyscon
  use CanopyStateType     , only : canopystate_type
  use PhotosynthesisType  , only : photosyns_type
  use CNStateType         , only : cnstate_type
  use CNCarbonFluxType    , only : carbonflux_type
  use CNCarbonStateType   , only : carbonstate_type
  use CNNitrogenFluxType  , only : nitrogenflux_type
  use CNNitrogenStateType , only : nitrogenstate_type
  use PhosphorusFluxType  , only : phosphorusflux_type
  use PhosphorusStateType , only : phosphorusstate_type
  use LandunitType        , only : lun                
  use ColumnType          , only : col                
  use PatchType           , only : pft                
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
    real(r8)                      :: deadstemc_seed, deadstemn_seed
    real(r8), pointer             :: dwt_ptr0, dwt_ptr1, dwt_ptr2, dwt_ptr3, ptr
    character(len=32)             :: subname='dyn_cbal'            ! subroutine name

    ! ! add phosphorus local variables
    real(r8), allocatable         :: dwt_leafp_seed(:)             ! pft-level mass gain due to seeding of new area    
    real(r8), allocatable         :: dwt_deadstemp_seed(:)         ! pft-level mass gain due to seeding of new area
    real(r8), allocatable, target :: dwt_frootp_to_litter(:)       ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_livecrootp_to_litter(:)   ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_deadcrootp_to_litter(:)   ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: conv_pflux(:)                 ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod10_pflux(:)               ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod100_pflux(:)              ! pft-level mass loss due to weight shift
    real(r8)                      :: leafp_seed
    real(r8)                      :: deadstemp_seed



    !! C13
    real(r8), allocatable         :: dwt_leafc13_seed(:)           ! pft-level mass gain due to seeding of new area
    real(r8), allocatable         :: dwt_deadstemc13_seed(:)       ! pft-level mass gain due to seeding of new area
    real(r8), allocatable, target :: dwt_frootc13_to_litter(:)     ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_livecrootc13_to_litter(:) ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_deadcrootc13_to_litter(:) ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: conv_c13flux(:)               ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod10_c13flux(:)             ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod100_c13flux(:)            ! pft-level mass loss due to weight shift
    real(r8)                      :: c3_del13c                     ! typical del13C for C3 photosynthesis (permil, relative to PDB)
    real(r8)                      :: c4_del13c                     ! typical del13C for C4 photosynthesis (permil, relative to PDB)
    real(r8)                      :: c3_r1_c13                     ! isotope ratio (13c/12c) for C3 photosynthesis
    real(r8)                      :: c4_r1_c13                     ! isotope ratio (13c/12c) for C4 photosynthesis
    real(r8)                      :: c3_r2_c13                     ! isotope ratio (13c/[12c+13c]) for C3 photosynthesis
    real(r8)                      :: c4_r2_c13                     ! isotope ratio (13c/[12c+13c]) for C4 photosynthesis
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
    real(r8)                      :: c3_del14c                     ! typical del14C for C3 photosynthesis (permil, relative to PDB)
    real(r8)                      :: c4_del14c                     ! typical del14C for C4 photosynthesis (permil, relative to PDB)
    real(r8)                      :: c3_r1_c14                     ! isotope ratio (14c/12c) for C3 photosynthesis
    real(r8)                      :: c4_r1_c14                     ! isotope ratio (14c/12c) for C4 photosynthesis
    real(r8)                      :: c3_r2_c14                     ! isotope ratio (14c/[12c+14c]) for C3 photosynthesis
    real(r8)                      :: c4_r2_c14                     ! isotope ratio (14c/[12c+14c]) for C4 photosynthesis
    real(r8)                      :: leafc14_seed, deadstemc14_seed
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
    allocate(dwt_leafc_seed(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_leafc_seed'
          call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    allocate(dwt_leafn_seed(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_leafn_seed'
          call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    allocate(dwt_deadstemc_seed(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadstemc_seed'
          call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    allocate(dwt_deadstemn_seed(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadstemn_seed'
          call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    allocate(dwt_frootc_to_litter(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_frootc_to_litter'
          call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    allocate(dwt_livecrootc_to_litter(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_livecrootc_to_litter'
          call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    allocate(dwt_deadcrootc_to_litter(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadcrootc_to_litter'
          call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    allocate(dwt_frootn_to_litter(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_frootn_to_litter'
          call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    allocate(dwt_livecrootn_to_litter(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_livecrootn_to_litter'
          call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    allocate(dwt_deadcrootn_to_litter(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadcrootn_to_litter'
          call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    allocate(conv_cflux(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for conv_cflux'
          call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    allocate(prod10_cflux(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for prod10_cflux'
          call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    allocate(prod100_cflux(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for prod100_cflux'
          call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    allocate(conv_nflux(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for conv_nflux'
          call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    allocate(prod10_nflux(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for prod10_nflux'
          call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    allocate(prod100_nflux(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for prod100_nflux'
          call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    ! Allocate P arrays
    allocate(dwt_leafp_seed(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_leafp_seed'
          call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    allocate(dwt_deadstemp_seed(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadstemp_seed'
          call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    allocate(dwt_frootp_to_litter(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_frootp_to_litter'
          call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    allocate(dwt_livecrootp_to_litter(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_livecrootp_to_litter'
          call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    allocate(dwt_deadcrootp_to_litter(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadcrootp_to_litter'
          call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    allocate(conv_pflux(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for conv_pflux'
          call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    allocate(prod10_pflux(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for prod10_pflux'
          call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    allocate(prod100_pflux(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
          write(iulog,*)subname,' allocation error for prod100_pflux'
          call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    if ( use_c13 ) then
       allocate(dwt_leafc13_seed(bounds%begp:bounds%endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_leafc13_seed'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
       allocate(dwt_deadstemc13_seed(bounds%begp:bounds%endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadstemc13_seed'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
       allocate(dwt_frootc13_to_litter(bounds%begp:bounds%endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_frootc13_to_litter'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
       allocate(dwt_livecrootc13_to_litter(bounds%begp:bounds%endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_livecrootc13_to_litter'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
       allocate(dwt_deadcrootc13_to_litter(bounds%begp:bounds%endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadcrootc13_to_litter'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
       allocate(conv_c13flux(bounds%begp:bounds%endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for conv_c13flux'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
       allocate(prod10_c13flux(bounds%begp:bounds%endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for prod10_c13flux'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
       allocate(prod100_c13flux(bounds%begp:bounds%endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for prod100_c13flux'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
    endif
    if ( use_c14 ) then
       allocate(dwt_leafc14_seed(bounds%begp:bounds%endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_leafc14_seed'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
       allocate(dwt_deadstemc14_seed(bounds%begp:bounds%endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadstemc14_seed'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
       allocate(dwt_frootc14_to_litter(bounds%begp:bounds%endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_frootc14_to_litter'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
       allocate(dwt_livecrootc14_to_litter(bounds%begp:bounds%endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_livecrootc14_to_litter'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
       allocate(dwt_deadcrootc14_to_litter(bounds%begp:bounds%endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for dwt_deadcrootc14_to_litter'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
       allocate(conv_c14flux(bounds%begp:bounds%endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for conv_c14flux'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
       allocate(prod10_c14flux(bounds%begp:bounds%endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for prod10_c14flux'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
       allocate(prod100_c14flux(bounds%begp:bounds%endp), stat=ier)
       if (ier /= 0) then
          write(iulog,*)subname,' allocation error for prod100_c14flux'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
    endif
    
    ! Get time step
    dt = real( get_step_size(), r8 )
    
    do p = bounds%begp,bounds%endp
       c = pft%column(p)
       ! initialize all the pft-level local flux arrays
       dwt_leafc_seed(p) = 0._r8
       dwt_leafn_seed(p) = 0._r8
       dwt_deadstemc_seed(p) = 0._r8
       dwt_deadstemn_seed(p) = 0._r8
       dwt_frootc_to_litter(p) = 0._r8
       dwt_livecrootc_to_litter(p) = 0._r8
       dwt_deadcrootc_to_litter(p) = 0._r8
       dwt_frootn_to_litter(p) = 0._r8
       dwt_livecrootn_to_litter(p) = 0._r8
       dwt_deadcrootn_to_litter(p) = 0._r8
       conv_cflux(p) = 0._r8
       prod10_cflux(p) = 0._r8
       prod100_cflux(p) = 0._r8
       conv_nflux(p) = 0._r8
       prod10_nflux(p) = 0._r8
       prod100_nflux(p) = 0._r8
       
       dwt_leafp_seed(p) = 0._r8
       dwt_deadstemp_seed(p) = 0._r8
       dwt_frootp_to_litter(p) = 0._r8
       dwt_livecrootp_to_litter(p) = 0._r8
       dwt_deadcrootp_to_litter(p) = 0._r8
       conv_pflux(p) = 0._r8
       prod10_pflux(p) = 0._r8
       prod100_pflux(p) = 0._r8

       if ( use_c13 ) then
          dwt_leafc13_seed(p) = 0._r8
          dwt_deadstemc13_seed(p) = 0._r8
          dwt_frootc13_to_litter(p) = 0._r8
          dwt_livecrootc13_to_litter(p) = 0._r8
          dwt_deadcrootc13_to_litter(p) = 0._r8
          conv_c13flux(p) = 0._r8
          prod10_c13flux(p) = 0._r8
          prod100_c13flux(p) = 0._r8
       endif
       
       if ( use_c14 ) then
          dwt_leafc14_seed(p) = 0._r8
          dwt_deadstemc14_seed(p) = 0._r8
          dwt_frootc14_to_litter(p) = 0._r8
          dwt_livecrootc14_to_litter(p) = 0._r8
          dwt_deadcrootc14_to_litter(p) = 0._r8
          conv_c14flux(p) = 0._r8
          prod10_c14flux(p) = 0._r8
          prod100_c14flux(p) = 0._r8
       endif
       
       l = pft%landunit(p)
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          
          ! calculate the change in weight for the timestep
          dwt = pft%wtcol(p)-prior_weights%pwtcol(p)
          cnstate_vars%lfpftd_patch(p) = -dwt

          ! Patches for which weight increases on this timestep
          if (dwt > 0._r8) then
             
             ! first identify Patches that are initiating on this timestep
             ! and set all the necessary state and flux variables
             if (prior_weights%pwtcol(p) == 0._r8) then
                
                ! set initial conditions for PFT that is being initiated
                ! in this time step.  Based on the settings in cnIniTimeVar.
                
                ! pft-level carbon state variables
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
                
                if ( use_c13 ) then
                   ! pft-level carbon-13 state variables
                   c13_cs%leafc_patch(p)              = 0._r8
                   c13_cs%leafc_storage_patch(p)      = 0._r8
                   c13_cs%leafc_xfer_patch(p)         = 0._r8
                   c13_cs%frootc_patch(p)             = 0._r8
                   c13_cs%frootc_storage_patch(p)     = 0._r8
                   c13_cs%frootc_xfer_patch(p)        = 0._r8
                   c13_cs%livestemc_patch(p)          = 0._r8
                   c13_cs%livestemc_storage_patch(p)  = 0._r8
                   c13_cs%livestemc_xfer_patch(p)     = 0._r8
                   c13_cs%deadstemc_patch(p)          = 0._r8
                   c13_cs%deadstemc_storage_patch(p)  = 0._r8
                   c13_cs%deadstemc_xfer_patch(p)     = 0._r8
                   c13_cs%livecrootc_patch(p)         = 0._r8
                   c13_cs%livecrootc_storage_patch(p) = 0._r8
                   c13_cs%livecrootc_xfer_patch(p)    = 0._r8
                   c13_cs%deadcrootc_patch(p)         = 0._r8
                   c13_cs%deadcrootc_storage_patch(p) = 0._r8
                   c13_cs%deadcrootc_xfer_patch(p)    = 0._r8
                   c13_cs%gresp_storage_patch(p)      = 0._r8
                   c13_cs%gresp_xfer_patch(p)         = 0._r8
                   c13_cs%cpool_patch(p)              = 0._r8
                   c13_cs%xsmrpool_patch(p)           = 0._r8
                   c13_cs%ctrunc_patch(p)             = 0._r8
                   c13_cs%dispvegc_patch(p)           = 0._r8
                   c13_cs%storvegc_patch(p)           = 0._r8
                   c13_cs%totvegc_patch(p)            = 0._r8
                   c13_cs%totpftc_patch(p)            = 0._r8
                endif
                
                if ( use_c14 ) then
                   ! pft-level carbon-14 state variables
                   c14_cs%leafc_patch(p)              = 0._r8
                   c14_cs%leafc_storage_patch(p)      = 0._r8
                   c14_cs%leafc_xfer_patch(p)         = 0._r8
                   c14_cs%frootc_patch(p)             = 0._r8
                   c14_cs%frootc_storage_patch(p)     = 0._r8
                   c14_cs%frootc_xfer_patch(p)        = 0._r8
                   c14_cs%livestemc_patch(p)          = 0._r8
                   c14_cs%livestemc_storage_patch(p)  = 0._r8
                   c14_cs%livestemc_xfer_patch(p)     = 0._r8
                   c14_cs%deadstemc_patch(p)          = 0._r8
                   c14_cs%deadstemc_storage_patch(p)  = 0._r8
                   c14_cs%deadstemc_xfer_patch(p)     = 0._r8
                   c14_cs%livecrootc_patch(p)         = 0._r8
                   c14_cs%livecrootc_storage_patch(p) = 0._r8
                   c14_cs%livecrootc_xfer_patch(p)    = 0._r8
                   c14_cs%deadcrootc_patch(p)         = 0._r8
                   c14_cs%deadcrootc_storage_patch(p) = 0._r8
                   c14_cs%deadcrootc_xfer_patch(p)    = 0._r8
                   c14_cs%gresp_storage_patch(p)      = 0._r8
                   c14_cs%gresp_xfer_patch(p)         = 0._r8
                   c14_cs%cpool_patch(p)              = 0._r8
                   c14_cs%xsmrpool_patch(p)           = 0._r8
                   c14_cs%ctrunc_patch(p)             = 0._r8
                   c14_cs%dispvegc_patch(p)           = 0._r8
                   c14_cs%storvegc_patch(p)           = 0._r8
                   c14_cs%totvegc_patch(p)            = 0._r8
                   c14_cs%totpftc_patch(p)            = 0._r8
                endif
                
                ! pft-level nitrogen state variables
                ns%leafn_patch(p)	       = 0._r8
                ns%leafn_storage_patch(p)      = 0._r8
                ns%leafn_xfer_patch(p)         = 0._r8
                ns%frootn_patch(p)	       = 0._r8
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
                ns%npool_patch(p)	       = 0._r8
                ns%ntrunc_patch(p)             = 0._r8
                ns%dispvegn_patch(p)           = 0._r8
                ns%storvegn_patch(p)           = 0._r8
                ns%totvegn_patch(p)            = 0._r8
                ns%totpftn_patch (p)           = 0._r8
                
                ! pft-level phosphorus state variables
                ps%leafp_patch(p)	       = 0._r8
                ps%leafp_storage_patch(p)      = 0._r8
                ps%leafp_xfer_patch(p)         = 0._r8
                ps%frootp_patch(p)	       = 0._r8
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
                ps%ppool_patch(p)	       = 0._r8
                ps%ptrunc_patch(p)             = 0._r8
                ps%dispvegp_patch(p)           = 0._r8
                ps%storvegp_patch(p)           = 0._r8
                ps%totvegp_patch(p)            = 0._r8
                ps%totpftp_patch (p)           = 0._r8

                ! initialize same flux and epv variables that are set
                canopystate_vars%laisun_patch(p) = 0._r8
                canopystate_vars%laisha_patch(p) = 0._r8
                
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
 
                ! add phosphorus related variables
                cnstate_vars%tempmax_retransp_patch(p)      = 0._r8
                cnstate_vars%annmax_retransp_patch(p)       = 0._r8
           

                cf%xsmrpool_recover_patch(p)      = 0._r8
                cf%plant_calloc_patch(p)          = 0._r8
                cf%excess_cflux_patch(p)          = 0._r8
                cf%prev_leafc_to_litter_patch(p)  = 0._r8
                cf%prev_frootc_to_litter_patch(p) = 0._r8
                cf%availc_patch(p)                = 0._r8
                cf%gpp_before_downreg_patch(p)    = 0._r8
                cf%tempsum_npp_patch(p)           = 0._r8
                cf%annsum_npp_patch(p)            = 0._r8

                nf%plant_ndemand_patch(p)         = 0._r8
                nf%avail_retransn_patch(p)        = 0._r8
                nf%plant_nalloc_patch(p)          = 0._r8
                 
                ! add phosphorus related variables
                pf%plant_pdemand_patch(p)         = 0._r8
                pf%avail_retransp_patch(p)        = 0._r8
                pf%plant_palloc_patch(p)          = 0._r8

                if ( use_c13 ) then
                   cf%xsmrpool_c13ratio_patch(p) = c13ratio
                end if
                if ( use_c14 ) then
                   cnstate_vars%rc14_atm_patch(p) = c14ratio
                   cnstate_vars%rc14_atm_patch(p) = 0._r8
                endif

                call photosyns_vars%NewPatchinit(p)
                
             end if  ! end initialization of new pft
             
             ! (still in dwt > 0 block)
             
             ! set the seed sources for leaf and deadstem
             ! leaf source is split later between leaf, leaf_storage, leaf_xfer
             leafc_seed   = 0._r8
             leafn_seed   = 0._r8
             leafp_seed   = 0._r8
             deadstemc_seed   = 0._r8
             deadstemn_seed   = 0._r8
             deadstemp_seed   = 0._r8
             if ( use_c13 ) then
                leafc13_seed = 0._r8
                deadstemc13_seed = 0._r8
             endif
             if ( use_c14 ) then
                leafc14_seed = 0._r8
                deadstemc14_seed = 0._r8
             endif
             if (pft%itype(p) /= 0) then
                leafc_seed = 1._r8
                leafn_seed  = leafc_seed / ecophyscon%leafcn(pft%itype(p))
                leafp_seed  = leafc_seed / ecophyscon%leafcp(pft%itype(p))
                if (ecophyscon%woody(pft%itype(p)) == 1._r8) then
                   deadstemc_seed = 0.1_r8
                   deadstemn_seed = deadstemc_seed / ecophyscon%deadwdcn(pft%itype(p))
                   deadstemp_seed = deadstemc_seed / ecophyscon%deadwdcp(pft%itype(p))
                end if
                
                if ( use_c13 ) then
                   ! 13c state is initialized assuming del13c = -28 permil for C3, and -13 permil for C4.
                   ! That translates to ratios of (13c/(12c+13c)) of 0.01080455 for C3, and 0.01096945 for C4
                   ! based on the following formulae: 
                   ! r1 (13/12) = PDB + (del13c * PDB)/1000.0
                   ! r2 (13/(13+12)) = r1/(1+r1)
                   ! PDB = 0.0112372_R8  (ratio of 13C/12C in Pee Dee Belemnite, C isotope standard)
                   c3_del13c = -28._r8
                   c4_del13c = -13._r8
                   c3_r1_c13 = SHR_CONST_PDB + ((c3_del13c*SHR_CONST_PDB)/1000._r8)
                   c3_r2_c13 = c3_r1_c13/(1._r8 + c3_r1_c13)
                   c4_r1_c13 = SHR_CONST_PDB + ((c4_del13c*SHR_CONST_PDB)/1000._r8)
                   c4_r2_c13 = c4_r1_c13/(1._r8 + c4_r1_c13)
                   
                   if (ecophyscon%c3psn(pft%itype(p)) == 1._r8) then
                      leafc13_seed     = leafc_seed     * c3_r2_c13
                      deadstemc13_seed = deadstemc_seed * c3_r2_c13
                   else
                      leafc13_seed     = leafc_seed     * c4_r2_c13
                      deadstemc13_seed = deadstemc_seed * c4_r2_c13
                   end if
                endif
                
                if ( use_c14 ) then
                   ! 14c state is initialized assuming initial "modern" 14C of 1.e-12
                   if (ecophyscon%c3psn(pft%itype(p)) == 1._r8) then
                      leafc14_seed     = leafc_seed     * c14ratio
                      deadstemc14_seed = deadstemc_seed * c14ratio
                   else
                      leafc14_seed     = leafc_seed     * c14ratio
                      deadstemc14_seed = deadstemc_seed * c14ratio
                   end if
                endif
             end if
             
             ! When PFT area expands (dwt > 0), the pft-level mass density 
             ! is modified to conserve the original pft mass distributed
             ! over the new (larger) area, plus a term to account for the 
             ! introduction of new seed source for leaf and deadstem
             t1 = prior_weights%pwtcol(p)/pft%wtcol(p)
             t2 = dwt/pft%wtcol(p)
             
             tot_leaf = cs%leafc_patch(p) + cs%leafc_storage_patch(p) + cs%leafc_xfer_patch(p)
             pleaf = 0._r8
             pstor = 0._r8
             pxfer = 0._r8
             if (tot_leaf /= 0._r8) then
                ! when adding seed source to non-zero leaf state, use current proportions
                pleaf = cs%leafc_patch(p)/tot_leaf
                pstor = cs%leafc_storage_patch(p)/tot_leaf
                pxfer = cs%leafc_xfer_patch(p)/tot_leaf
             else
                ! when initiating from zero leaf state, use evergreen flag to set proportions
                if (ecophyscon%evergreen(pft%itype(p)) == 1._r8) then
                   pleaf = 1._r8
                else
                   pstor = 1._r8
                end if
             end if
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
             
             if ( use_c13 ) then
                ! pft-level carbon-13 state variables 
                tot_leaf = c13_cs%leafc_patch(p) + c13_cs%leafc_storage_patch(p) + c13_cs%leafc_xfer_patch(p)
                pleaf = 0._r8
                pstor = 0._r8
                pxfer = 0._r8
                if (tot_leaf /= 0._r8) then
                   pleaf = c13_cs%leafc_patch(p)/tot_leaf
                   pstor = c13_cs%leafc_storage_patch(p)/tot_leaf
                   pxfer = c13_cs%leafc_xfer_patch(p)/tot_leaf
                else
                   ! when initiating from zero leaf state, use evergreen flag to set proportions
                   if (ecophyscon%evergreen(pft%itype(p)) == 1._r8) then
                      pleaf = 1._r8
                   else
                      pstor = 1._r8
                   end if
                end if
                c13_cs%leafc_patch(p)              = c13_cs%leafc_patch(p)              * t1 + leafc13_seed*pleaf*t2
                c13_cs%leafc_storage_patch(p)      = c13_cs%leafc_storage_patch(p)      * t1 + leafc13_seed*pstor*t2
                c13_cs%leafc_xfer_patch(p)         = c13_cs%leafc_xfer_patch(p)         * t1 + leafc13_seed*pxfer*t2
                c13_cs%frootc_patch(p)             = c13_cs%frootc_patch(p)             * t1
                c13_cs%frootc_storage_patch(p)     = c13_cs%frootc_storage_patch(p)     * t1
                c13_cs%frootc_xfer_patch(p)        = c13_cs%frootc_xfer_patch(p)        * t1
                c13_cs%livestemc_patch(p)          = c13_cs%livestemc_patch(p)          * t1
                c13_cs%livestemc_storage_patch(p)  = c13_cs%livestemc_storage_patch(p)  * t1
                c13_cs%livestemc_xfer_patch(p)     = c13_cs%livestemc_xfer_patch(p)     * t1
                c13_cs%deadstemc_patch(p)          = c13_cs%deadstemc_patch(p)          * t1 + deadstemc13_seed*t2
                c13_cs%deadstemc_storage_patch(p)  = c13_cs%deadstemc_storage_patch(p)  * t1
                c13_cs%deadstemc_xfer_patch(p)     = c13_cs%deadstemc_xfer_patch(p)     * t1
                c13_cs%livecrootc_patch(p)         = c13_cs%livecrootc_patch(p)         * t1
                c13_cs%livecrootc_storage_patch(p) = c13_cs%livecrootc_storage_patch(p) * t1
                c13_cs%livecrootc_xfer_patch(p)    = c13_cs%livecrootc_xfer_patch(p)    * t1
                c13_cs%deadcrootc_patch(p)         = c13_cs%deadcrootc_patch(p)         * t1
                c13_cs%deadcrootc_storage_patch(p) = c13_cs%deadcrootc_storage_patch(p) * t1
                c13_cs%deadcrootc_xfer_patch(p)    = c13_cs%deadcrootc_xfer_patch(p)    * t1
                c13_cs%gresp_storage_patch(p)      = c13_cs%gresp_storage_patch(p)      * t1
                c13_cs%gresp_xfer_patch(p)         = c13_cs%gresp_xfer_patch(p)         * t1
                c13_cs%cpool_patch(p)              = c13_cs%cpool_patch(p)              * t1
                c13_cs%xsmrpool_patch(p)           = c13_cs%xsmrpool_patch(p)           * t1
                c13_cs%ctrunc_patch(p)             = c13_cs%ctrunc_patch(p)             * t1
                c13_cs%dispvegc_patch(p)           = c13_cs%dispvegc_patch(p)           * t1
                c13_cs%storvegc_patch(p)           = c13_cs%storvegc_patch(p)           * t1
                c13_cs%totvegc_patch(p)            = c13_cs%totvegc_patch(p)            * t1
                c13_cs%totpftc_patch(p)            = c13_cs%totpftc_patch(p)            * t1
                
             endif
             
             if ( use_c14 ) then
                ! pft-level carbon-14 state variables 
                tot_leaf = c14_cs%leafc_patch(p) + c14_cs%leafc_storage_patch(p) + c14_cs%leafc_xfer_patch(p)
                pleaf = 0._r8
                pstor = 0._r8
                pxfer = 0._r8
                if (tot_leaf /= 0._r8) then
                   pleaf = c14_cs%leafc_patch(p)/tot_leaf
                   pstor = c14_cs%leafc_storage_patch(p)/tot_leaf
                   pxfer = c14_cs%leafc_xfer_patch(p)/tot_leaf
                else
                   ! when initiating from zero leaf state, use evergreen flag to set proportions
                   if (ecophyscon%evergreen(pft%itype(p)) == 1._r8) then
                      pleaf = 1._r8
                   else
                      pstor = 1._r8
                   end if
                end if
                c14_cs%leafc_patch(p)              = c14_cs%leafc_patch(p)              * t1 + leafc14_seed*pleaf*t2
                c14_cs%leafc_storage_patch(p)      = c14_cs%leafc_storage_patch(p)      * t1 + leafc14_seed*pstor*t2
                c14_cs%leafc_xfer_patch(p)         = c14_cs%leafc_xfer_patch(p)         * t1 + leafc14_seed*pxfer*t2
                c14_cs%frootc_patch(p)             = c14_cs%frootc_patch(p)             * t1
                c14_cs%frootc_storage_patch(p)     = c14_cs%frootc_storage_patch(p)     * t1
                c14_cs%frootc_xfer_patch(p)        = c14_cs%frootc_xfer_patch(p)        * t1
                c14_cs%livestemc_patch(p)          = c14_cs%livestemc_patch(p)          * t1
                c14_cs%livestemc_storage_patch(p)  = c14_cs%livestemc_storage_patch(p)  * t1
                c14_cs%livestemc_xfer_patch(p)     = c14_cs%livestemc_xfer_patch(p)     * t1
                c14_cs%deadstemc_patch(p)          = c14_cs%deadstemc_patch(p)          * t1 + deadstemc14_seed*t2
                c14_cs%deadstemc_storage_patch(p)  = c14_cs%deadstemc_storage_patch(p)  * t1
                c14_cs%deadstemc_xfer_patch(p)     = c14_cs%deadstemc_xfer_patch(p)     * t1
                c14_cs%livecrootc_patch(p)         = c14_cs%livecrootc_patch(p)         * t1
                c14_cs%livecrootc_storage_patch(p) = c14_cs%livecrootc_storage_patch(p) * t1
                c14_cs%livecrootc_xfer_patch(p)    = c14_cs%livecrootc_xfer_patch(p)    * t1
                c14_cs%deadcrootc_patch(p)         = c14_cs%deadcrootc_patch(p)         * t1
                c14_cs%deadcrootc_storage_patch(p) = c14_cs%deadcrootc_storage_patch(p) * t1
                c14_cs%deadcrootc_xfer_patch(p)    = c14_cs%deadcrootc_xfer_patch(p)    * t1
                c14_cs%gresp_storage_patch(p)      = c14_cs%gresp_storage_patch(p)      * t1
                c14_cs%gresp_xfer_patch(p)         = c14_cs%gresp_xfer_patch(p)         * t1
                c14_cs%cpool_patch(p)              = c14_cs%cpool_patch(p)              * t1
                c14_cs%xsmrpool_patch(p)           = c14_cs%xsmrpool_patch(p)           * t1
                c14_cs%ctrunc_patch(p)             = c14_cs%ctrunc_patch(p)             * t1
                c14_cs%dispvegc_patch(p)           = c14_cs%dispvegc_patch(p)           * t1
                c14_cs%storvegc_patch(p)           = c14_cs%storvegc_patch(p)           * t1
                c14_cs%totvegc_patch(p)            = c14_cs%totvegc_patch(p)            * t1
                c14_cs%totpftc_patch(p)            = c14_cs%totpftc_patch(p)            * t1
             endif
             
             tot_leaf = ns%leafn_patch(p) + ns%leafn_storage_patch(p) + ns%leafn_xfer_patch(p)
             pleaf = 0._r8
             pstor = 0._r8
             pxfer = 0._r8
             if (tot_leaf /= 0._r8) then
                pleaf = ns%leafn_patch(p)/tot_leaf
                pstor = ns%leafn_storage_patch(p)/tot_leaf
                pxfer = ns%leafn_xfer_patch(p)/tot_leaf
             else
                ! when initiating from zero leaf state, use evergreen flag to set proportions
                if (ecophyscon%evergreen(pft%itype(p)) == 1._r8) then
                   pleaf = 1._r8
                else
                   pstor = 1._r8
                end if
             end if
             ! pft-level nitrogen state variables
             ns%leafn_patch(p)              = ns%leafn_patch(p)              * t1 + leafn_seed*pleaf*t2
             ns%leafn_storage_patch(p)      = ns%leafn_storage_patch(p)      * t1 + leafn_seed*pstor*t2
             ns%leafn_xfer_patch(p)         = ns%leafn_xfer_patch(p)         * t1 + leafn_seed*pxfer*t2
             ns%frootn_patch(p)  	    = ns%frootn_patch(p)             * t1
             ns%frootn_storage_patch(p)     = ns%frootn_storage_patch(p)     * t1
             ns%frootn_xfer_patch(p)        = ns%frootn_xfer_patch(p)        * t1
             ns%livestemn_patch(p)	    = ns%livestemn_patch(p)          * t1
             ns%livestemn_storage_patch(p)  = ns%livestemn_storage_patch(p)  * t1
             ns%livestemn_xfer_patch(p)     = ns%livestemn_xfer_patch(p)     * t1
             ns%deadstemn_patch(p)          = ns%deadstemn_patch(p)          * t1 + deadstemn_seed*t2
             ns%deadstemn_storage_patch(p)  = ns%deadstemn_storage_patch(p)  * t1
             ns%deadstemn_xfer_patch(p)     = ns%deadstemn_xfer_patch(p)     * t1
             ns%livecrootn_patch(p)         = ns%livecrootn_patch(p)         * t1
             ns%livecrootn_storage_patch(p) = ns%livecrootn_storage_patch(p) * t1
             ns%livecrootn_xfer_patch(p)    = ns%livecrootn_xfer_patch(p)    * t1
             ns%deadcrootn_patch(p)         = ns%deadcrootn_patch(p)         * t1
             ns%deadcrootn_storage_patch(p) = ns%deadcrootn_storage_patch(p) * t1
             ns%deadcrootn_xfer_patch(p)    = ns%deadcrootn_xfer_patch(p)    * t1
             ns%retransn_patch(p)	    = ns%retransn_patch(p)           * t1
             ns%npool_patch(p)              = ns%npool_patch(p)              * t1
             ns%ntrunc_patch(p)             = ns%ntrunc_patch(p)             * t1
             ns%dispvegn_patch(p)	    = ns%dispvegn_patch(p)           * t1
             ns%storvegn_patch(p)	    = ns%storvegn_patch(p)           * t1
             ns%totvegn_patch(p) 	    = ns%totvegn_patch(p)            * t1
             ns%totpftn_patch(p) 	    = ns%totpftn_patch(p)            * t1
             
             ! add phosphorus -X.YANG
             tot_leaf = ps%leafp_patch(p) + ps%leafp_storage_patch(p) + ps%leafp_xfer_patch(p)
             pleaf = 0._r8
             pstor = 0._r8
             pxfer = 0._r8
             if (tot_leaf /= 0._r8) then
                pleaf = ps%leafp_patch(p)/tot_leaf
                pstor = ps%leafp_storage_patch(p)/tot_leaf
                pxfer = ps%leafp_xfer_patch(p)/tot_leaf
             else
                ! when initiating from zero leaf state, use evergreen flag to set proportions
                if (ecophyscon%evergreen(pft%itype(p)) == 1._r8) then
                   pleaf = 1._r8
                else
                   pstor = 1._r8
                end if
             end if
             ! pft-level nitrogen state variables
             ps%leafp_patch(p)              = ps%leafp_patch(p)              * t1 + leafp_seed*pleaf*t2
             ps%leafp_storage_patch(p)      = ps%leafp_storage_patch(p)      * t1 + leafp_seed*pstor*t2
             ps%leafp_xfer_patch(p)         = ps%leafp_xfer_patch(p)         * t1 + leafp_seed*pxfer*t2
             ps%frootp_patch(p)  	    = ps%frootp_patch(p)             * t1
             ps%frootp_storage_patch(p)     = ps%frootp_storage_patch(p)     * t1
             ps%frootp_xfer_patch(p)        = ps%frootp_xfer_patch(p)        * t1
             ps%livestemp_patch(p)	    = ps%livestemp_patch(p)          * t1
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
             ps%retransp_patch(p)	    = ps%retransp_patch(p)           * t1
             ps%ppool_patch(p)              = ps%ppool_patch(p)              * t1
             ps%ptrunc_patch(p)             = ps%ptrunc_patch(p)             * t1
             ps%dispvegp_patch(p)	    = ps%dispvegp_patch(p)           * t1
             ps%storvegp_patch(p)	    = ps%storvegp_patch(p)           * t1
             ps%totvegp_patch(p) 	    = ps%totvegp_patch(p)            * t1
             ps%totpftp_patch(p) 	    = ps%totpftp_patch(p)            * t1


             ! update temporary seed source arrays
             ! These are calculated in terms of the required contributions from
             ! column-level seed source
             dwt_leafc_seed(p)   = leafc_seed   * dwt
             if ( use_c13 ) then
                dwt_leafc13_seed(p) = leafc13_seed * dwt
                dwt_deadstemc13_seed(p) = deadstemc13_seed * dwt
             endif
             if ( use_c14 ) then
                dwt_leafc14_seed(p) = leafc14_seed * dwt
                dwt_deadstemc14_seed(p) = deadstemc14_seed * dwt
             endif
             dwt_leafn_seed(p)   = leafn_seed   * dwt
             dwt_deadstemc_seed(p)   = deadstemc_seed   * dwt
             dwt_deadstemn_seed(p)   = deadstemn_seed   * dwt
             
             dwt_leafp_seed(p)   = leafp_seed   * dwt
             dwt_deadstemp_seed(p)   = deadstemp_seed   * dwt
          else if (dwt < 0._r8) then
             
             ! if the pft lost weight on the timestep, then the carbon and nitrogen state
             ! variables are directed to litter, CWD, and wood product pools.
             
             ! N.B. : the conv_cflux, prod10_cflux, and prod100_cflux fluxes are accumulated
             ! as negative values, but the fluxes for pft-to-litter are accumulated as 
             ! positive values
             
             ! set local weight variables for this pft
             wt_new = pft%wtcol(p)
             wt_old = prior_weights%pwtcol(p)
             
             !---------------
             ! C state update
             !---------------
             
             ! leafc 
             ptr => cs%leafc_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                conv_cflux(p) = conv_cflux(p) + change_state
             else
                ptr = 0._r8
                conv_cflux(p) = conv_cflux(p) - init_state
             end if
             
             ! leafc_storage 
             ptr => cs%leafc_storage_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                conv_cflux(p) = conv_cflux(p) + change_state
             else
                ptr = 0._r8
                conv_cflux(p) = conv_cflux(p) - init_state
             end if
             
             ! leafc_xfer 
             ptr => cs%leafc_xfer_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                conv_cflux(p) = conv_cflux(p) + change_state
             else
                ptr = 0._r8
                conv_cflux(p) = conv_cflux(p) - init_state
             end if
             
             ! frootc 
             ptr => cs%frootc_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_frootc_to_litter(p) = dwt_frootc_to_litter(p) - change_state
             else
                ptr = 0._r8
                dwt_frootc_to_litter(p) = dwt_frootc_to_litter(p) + init_state
             end if
             
             ! frootc_storage 
             ptr => cs%frootc_storage_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                conv_cflux(p) = conv_cflux(p) + change_state
             else
                ptr = 0._r8
                conv_cflux(p) = conv_cflux(p) - init_state
             end if
             
             ! frootc_xfer 
             ptr => cs%frootc_xfer_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                conv_cflux(p) = conv_cflux(p) + change_state
             else
                ptr = 0._r8
                conv_cflux(p) = conv_cflux(p) - init_state
             end if
             
             ! livestemc 
             ptr => cs%livestemc_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                conv_cflux(p) = conv_cflux(p) + change_state
             else
                ptr = 0._r8
                conv_cflux(p) = conv_cflux(p) - init_state
             end if
             
             ! livestemc_storage 
             ptr => cs%livestemc_storage_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                conv_cflux(p) = conv_cflux(p) + change_state
             else
                ptr = 0._r8
                conv_cflux(p) = conv_cflux(p) - init_state
             end if
             
             ! livestemc_xfer 
             ptr => cs%livestemc_xfer_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                conv_cflux(p) = conv_cflux(p) + change_state
             else
                ptr = 0._r8
                conv_cflux(p) = conv_cflux(p) - init_state
             end if
             
             ! deadstemc 
             ptr => cs%deadstemc_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                conv_cflux(p) = conv_cflux(p) + change_state*pconv(pft%itype(p))
                prod10_cflux(p) = prod10_cflux(p) + change_state*pprod10(pft%itype(p))
                prod100_cflux(p) = prod100_cflux(p) + change_state*pprod100(pft%itype(p))
             else
                ptr = 0._r8
                conv_cflux(p) = conv_cflux(p) - init_state*pconv(pft%itype(p))
                prod10_cflux(p) = prod10_cflux(p) - init_state*pprod10(pft%itype(p))
                prod100_cflux(p) = prod100_cflux(p) - init_state*pprod100(pft%itype(p))
             end if
             
             ! deadstemc_storage 
             ptr => cs%deadstemc_storage_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                conv_cflux(p) = conv_cflux(p) + change_state
             else
                ptr = 0._r8
                conv_cflux(p) = conv_cflux(p) - init_state
             end if
             
             ! deadstemc_xfer 
             ptr => cs%deadstemc_xfer_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                conv_cflux(p) = conv_cflux(p) + change_state
             else
                ptr = 0._r8
                conv_cflux(p) = conv_cflux(p) - init_state
             end if
             
             ! livecrootc 
             ptr => cs%livecrootc_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_livecrootc_to_litter(p) = dwt_livecrootc_to_litter(p) - change_state
             else
                ptr = 0._r8
                dwt_livecrootc_to_litter(p) = dwt_livecrootc_to_litter(p) + init_state
             end if
             
             ! livecrootc_storage 
             ptr => cs%livecrootc_storage_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                conv_cflux(p) = conv_cflux(p) + change_state
             else
                ptr = 0._r8
                conv_cflux(p) = conv_cflux(p) - init_state
             end if
             
             ! livecrootc_xfer 
             ptr => cs%livecrootc_xfer_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                conv_cflux(p) = conv_cflux(p) + change_state
             else
                ptr = 0._r8
                conv_cflux(p) = conv_cflux(p) - init_state
             end if
             
             ! deadcrootc 
             ptr => cs%deadcrootc_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_deadcrootc_to_litter(p) = dwt_deadcrootc_to_litter(p) - change_state
             else
                ptr = 0._r8
                dwt_deadcrootc_to_litter(p) = dwt_deadcrootc_to_litter(p) + init_state
             end if
             
             ! deadcrootc_storage 
             ptr => cs%deadcrootc_storage_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                conv_cflux(p) = conv_cflux(p) + change_state
             else
                ptr = 0._r8
                conv_cflux(p) = conv_cflux(p) - init_state
             end if
             
             ! deadcrootc_xfer 
             ptr => cs%deadcrootc_xfer_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                conv_cflux(p) = conv_cflux(p) + change_state
             else
                ptr = 0._r8
                conv_cflux(p) = conv_cflux(p) - init_state
             end if
             
             ! gresp_storage 
             ptr => cs%gresp_storage_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                conv_cflux(p) = conv_cflux(p) + change_state
             else
                ptr = 0._r8
                conv_cflux(p) = conv_cflux(p) - init_state
             end if
             
             ! gresp_xfer 
             ptr => cs%gresp_xfer_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                conv_cflux(p) = conv_cflux(p) + change_state
             else
                ptr = 0._r8
                conv_cflux(p) = conv_cflux(p) - init_state
             end if
             
             ! cpool 
             ptr => cs%cpool_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                conv_cflux(p) = conv_cflux(p) + change_state
             else
                ptr = 0._r8
                conv_cflux(p) = conv_cflux(p) - init_state
             end if
             
             ! xsmrpool 
             ptr => cs%xsmrpool_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                conv_cflux(p) = conv_cflux(p) + change_state
             else
                ptr = 0._r8
                conv_cflux(p) = conv_cflux(p) - init_state
             end if
             
             ! pft_ctrunc 
             ptr => cs%ctrunc_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                conv_cflux(p) = conv_cflux(p) + change_state
             else
                ptr = 0._r8
                conv_cflux(p) = conv_cflux(p) - init_state
             end if

             if ( use_c13 ) then
                !-------------------
                ! C13 state update
                !-------------------
                
                ! set pointers to the conversion and product pool fluxes for this pft
                ! dwt_ptr0 is reserved for local assignment to dwt_xxx_to_litter fluxes
                dwt_ptr1 => conv_c13flux(p)
                dwt_ptr2 => prod10_c13flux(p)
                dwt_ptr3 => prod100_c13flux(p)
                
                ! leafc 
                ptr => c13_cs%leafc_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state
                end if
                
                ! leafc_storage 
                ptr => c13_cs%leafc_storage_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state
                end if
                
                ! leafc_xfer 
                ptr => c13_cs%leafc_xfer_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state
                end if
                
                ! frootc 
                ptr => c13_cs%frootc_patch(p)
                dwt_ptr0 => dwt_frootc13_to_litter(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr0 = dwt_ptr0 - change_state
                else
                   ptr = 0._r8
                   dwt_ptr0 = dwt_ptr0 + init_state
                end if
                
                ! frootc_storage 
                ptr => c13_cs%frootc_storage_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state
                end if
                
                ! frootc_xfer 
                ptr => c13_cs%frootc_xfer_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state
                end if
                
                ! livestemc 
                ptr => c13_cs%livestemc_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state
                end if
                
                ! livestemc_storage 
                ptr => c13_cs%livestemc_storage_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state
                end if
                
                ! livestemc_xfer 
                ptr => c13_cs%livestemc_xfer_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state
                end if
                
                ! deadstemc 
                ptr => c13_cs%deadstemc_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state*pconv(pft%itype(p))
                   dwt_ptr2 = dwt_ptr2 + change_state*pprod10(pft%itype(p))
                   dwt_ptr3 = dwt_ptr3 + change_state*pprod100(pft%itype(p))
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state*pconv(pft%itype(p))
                   dwt_ptr2 = dwt_ptr2 - init_state*pprod10(pft%itype(p))
                   dwt_ptr3 = dwt_ptr3 - init_state*pprod100(pft%itype(p))
                end if
                
                ! deadstemc_storage 
                ptr => c13_cs%deadstemc_storage_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state
                end if
                
                ! deadstemc_xfer 
                ptr => c13_cs%deadstemc_xfer_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state
                end if
                
                ! livecrootc 
                ptr => c13_cs%livecrootc_patch(p)
                dwt_ptr0 => dwt_livecrootc13_to_litter(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr0 = dwt_ptr0 - change_state
                else
                   ptr = 0._r8
                   dwt_ptr0 = dwt_ptr0 + init_state
                end if
                
                ! livecrootc_storage 
                ptr => c13_cs%livecrootc_storage_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state
                end if
                
                ! livecrootc_xfer 
                ptr => c13_cs%livecrootc_xfer_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state
                end if
                
                ! deadcrootc 
                ptr => c13_cs%deadcrootc_patch(p)
                dwt_ptr0 => dwt_deadcrootc13_to_litter(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr0 = dwt_ptr0 - change_state
                else
                   ptr = 0._r8
                   dwt_ptr0 = dwt_ptr0 + init_state
                end if
                
                ! deadcrootc_storage 
                ptr => c13_cs%deadcrootc_storage_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state
                end if
                
                ! deadcrootc_xfer 
                ptr => c13_cs%deadcrootc_xfer_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state
                end if
                
                ! gresp_storage 
                ptr => c13_cs%gresp_storage_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state
                end if
                
                ! gresp_xfer 
                ptr => c13_cs%gresp_xfer_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state
                end if
                
                ! cpool 
                ptr => c13_cs%cpool_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state
                end if
                
                ! pft_ctrunc 
                ptr => c13_cs%ctrunc_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state
                end if
                
             endif

             if ( use_c14 ) then
                !-------------------
                ! C14 state update
                !-------------------
                
                ! set pointers to the conversion and product pool fluxes for this pft
                ! dwt_ptr0 is reserved for local assignment to dwt_xxx_to_litter fluxes
                dwt_ptr1 => conv_c14flux(p)
                dwt_ptr2 => prod10_c14flux(p)
                dwt_ptr3 => prod100_c14flux(p)
                
                ! leafc 
                ptr => c14_cs%leafc_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state
                end if
                
                ! leafc_storage 
                ptr => c14_cs%leafc_storage_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state
                end if
                
                ! leafc_xfer 
                ptr => c14_cs%leafc_xfer_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state
                end if
                
                ! frootc 
                ptr => c14_cs%frootc_patch(p)
                dwt_ptr0 => dwt_frootc14_to_litter(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr0 = dwt_ptr0 - change_state
                else
                   ptr = 0._r8
                   dwt_ptr0 = dwt_ptr0 + init_state
                end if
                
                ! frootc_storage 
                ptr => c14_cs%frootc_storage_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state
                end if
                
                ! frootc_xfer 
                ptr => c14_cs%frootc_xfer_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state
                end if
                
                ! livestemc 
                ptr => c14_cs%livestemc_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state
                end if
                
                ! livestemc_storage 
                ptr => c14_cs%livestemc_storage_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state
                end if
                
                ! livestemc_xfer 
                ptr => c14_cs%livestemc_xfer_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state
                end if
                
                ! deadstemc 
                ptr => c14_cs%deadstemc_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state*pconv(pft%itype(p))
                   dwt_ptr2 = dwt_ptr2 + change_state*pprod10(pft%itype(p))
                   dwt_ptr3 = dwt_ptr3 + change_state*pprod100(pft%itype(p))
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state*pconv(pft%itype(p))
                   dwt_ptr2 = dwt_ptr2 - init_state*pprod10(pft%itype(p))
                   dwt_ptr3 = dwt_ptr3 - init_state*pprod100(pft%itype(p))
                end if
                
                ! deadstemc_storage 
                ptr => c14_cs%deadstemc_storage_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state
                end if
                
                ! deadstemc_xfer 
                ptr => c14_cs%deadstemc_xfer_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state
                end if
                
                ! livecrootc 
                ptr => c14_cs%livecrootc_patch(p)
                dwt_ptr0 => dwt_livecrootc14_to_litter(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr0 = dwt_ptr0 - change_state
                else
                   ptr = 0._r8
                   dwt_ptr0 = dwt_ptr0 + init_state
                end if
                
                ! livecrootc_storage 
                ptr => c14_cs%livecrootc_storage_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state
                end if
                
                ! livecrootc_xfer 
                ptr => c14_cs%livecrootc_xfer_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state
                end if
                
                ! deadcrootc 
                ptr => c14_cs%deadcrootc_patch(p)
                dwt_ptr0 => dwt_deadcrootc14_to_litter(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr0 = dwt_ptr0 - change_state
                else
                   ptr = 0._r8
                   dwt_ptr0 = dwt_ptr0 + init_state
                end if
                
                ! deadcrootc_storage 
                ptr => c14_cs%deadcrootc_storage_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state
                end if
                
                ! deadcrootc_xfer 
                ptr => c14_cs%deadcrootc_xfer_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state
                end if
                
                ! gresp_storage 
                ptr => c14_cs%gresp_storage_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state
                end if
                
                ! gresp_xfer 
                ptr => c14_cs%gresp_xfer_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state
                end if
                
                ! cpool 
                ptr => c14_cs%cpool_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state
                end if
                
                ! pft_ctrunc 
                ptr => c14_cs%ctrunc_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state
                end if
             endif
             
             
             !---------------
             ! N state update
             !---------------
             
             ! set pointers to the conversion and product pool fluxes for this pft
             ! dwt_ptr0 is reserved for local assignment to dwt_xxx_to_litter fluxes
             dwt_ptr1 => conv_nflux(p)
             dwt_ptr2 => prod10_nflux(p)
             dwt_ptr3 => prod100_nflux(p)
             
             ! leafn 
             ptr => ns%leafn_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr1 = dwt_ptr1 + change_state
             else
                ptr = 0._r8
                dwt_ptr1 = dwt_ptr1 - init_state
             end if
             
             ! leafn_storage  
             ptr => ns%leafn_storage_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr1 = dwt_ptr1 + change_state
             else
                ptr = 0._r8
                dwt_ptr1 = dwt_ptr1 - init_state
             end if
             
             ! leafn_xfer  
             ptr => ns%leafn_xfer_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr1 = dwt_ptr1 + change_state
             else
                ptr = 0._r8
                dwt_ptr1 = dwt_ptr1 - init_state
             end if
             
             ! frootn 
             ptr => ns%frootn_patch(p)
             dwt_ptr0 => dwt_frootn_to_litter(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr0 = dwt_ptr0 - change_state
             else
                ptr = 0._r8
                dwt_ptr0 = dwt_ptr0 + init_state
             end if
             
             ! frootn_storage 
             ptr => ns%frootn_storage_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr1 = dwt_ptr1 + change_state
             else
                ptr = 0._r8
                dwt_ptr1 = dwt_ptr1 - init_state
             end if
             
             ! frootn_xfer  
             ptr => ns%frootn_xfer_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr1 = dwt_ptr1 + change_state
             else
                ptr = 0._r8
                dwt_ptr1 = dwt_ptr1 - init_state
             end if
             
             ! livestemn  
             ptr => ns%livestemn_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr1 = dwt_ptr1 + change_state
             else
                ptr = 0._r8
                dwt_ptr1 = dwt_ptr1 - init_state
             end if
             
             ! livestemn_storage 
             ptr => ns%livestemn_storage_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr1 = dwt_ptr1 + change_state
             else
                ptr = 0._r8
                dwt_ptr1 = dwt_ptr1 - init_state
             end if
             
             ! livestemn_xfer 
             ptr => ns%livestemn_xfer_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr1 = dwt_ptr1 + change_state
             else
                ptr = 0._r8
                dwt_ptr1 = dwt_ptr1 - init_state
             end if
             
             ! deadstemn 
             ptr => ns%deadstemn_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr1 = dwt_ptr1 + change_state*pconv(pft%itype(p))
                dwt_ptr2 = dwt_ptr2 + change_state*pprod10(pft%itype(p))
                dwt_ptr3 = dwt_ptr3 + change_state*pprod100(pft%itype(p))
             else
                ptr = 0._r8
                dwt_ptr1 = dwt_ptr1 - init_state*pconv(pft%itype(p))
                dwt_ptr2 = dwt_ptr2 - init_state*pprod10(pft%itype(p))
                dwt_ptr3 = dwt_ptr3 - init_state*pprod100(pft%itype(p))
             end if
             
             ! deadstemn_storage 
             ptr => ns%deadstemn_storage_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr1 = dwt_ptr1 + change_state
             else
                ptr = 0._r8
                dwt_ptr1 = dwt_ptr1 - init_state
             end if
             
             ! deadstemn_xfer 
             ptr => ns%deadstemn_xfer_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr1 = dwt_ptr1 + change_state
             else
                ptr = 0._r8
                dwt_ptr1 = dwt_ptr1 - init_state
             end if
             
             ! livecrootn 
             ptr => ns%livecrootn_patch(p)
             dwt_ptr0 => dwt_livecrootn_to_litter(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr0 = dwt_ptr0 - change_state
             else
                ptr = 0._r8
                dwt_ptr0 = dwt_ptr0 + init_state
             end if
             
             ! livecrootn_storage  
             ptr => ns%livecrootn_storage_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr1 = dwt_ptr1 + change_state
             else
                ptr = 0._r8
                dwt_ptr1 = dwt_ptr1 - init_state
             end if
             
             ! livecrootn_xfer  
             ptr => ns%livecrootn_xfer_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr1 = dwt_ptr1 + change_state
             else
                ptr = 0._r8
                dwt_ptr1 = dwt_ptr1 - init_state
             end if
             
             ! deadcrootn 
             ptr => ns%deadcrootn_patch(p)
             dwt_ptr0 => dwt_deadcrootn_to_litter(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr0 = dwt_ptr0 - change_state
             else
                ptr = 0._r8
                dwt_ptr0 = dwt_ptr0 + init_state
             end if
             
             ! deadcrootn_storage  
             ptr => ns%deadcrootn_storage_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr1 = dwt_ptr1 + change_state
             else
                ptr = 0._r8
                dwt_ptr1 = dwt_ptr1 - init_state
             end if
             
             ! deadcrootn_xfer  
             ptr => ns%deadcrootn_xfer_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr1 = dwt_ptr1 + change_state
             else
                ptr = 0._r8
                dwt_ptr1 = dwt_ptr1 - init_state
             end if
             
             ! retransn  
             ptr => ns%retransn_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr1 = dwt_ptr1 + change_state
             else
                ptr = 0._r8
                dwt_ptr1 = dwt_ptr1 - init_state
             end if
             
             ! npool  
             ptr => ns%npool_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr1 = dwt_ptr1 + change_state
             else
                ptr = 0._r8
                dwt_ptr1 = dwt_ptr1 - init_state
             end if
             
             ! pft_ntrunc  
             ptr => ns%ntrunc_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr1 = dwt_ptr1 + change_state
             else
                ptr = 0._r8
                dwt_ptr1 = dwt_ptr1 - init_state
             end if


             !---------------
             ! P state update  X.YANG
             !---------------
             
             ! set pointers to the conversion and product pool fluxes for this pft
             ! dwt_ptr0 is reserved for local assignment to dwt_xxx_to_litter fluxes
             dwt_ptr1 => conv_pflux(p)
             dwt_ptr2 => prod10_pflux(p)
             dwt_ptr3 => prod100_pflux(p)
             
             ! leafp 
             ptr => ps%leafp_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr1 = dwt_ptr1 + change_state
             else
                ptr = 0._r8
                dwt_ptr1 = dwt_ptr1 - init_state
             end if
             
             ! leafp_storage  
             ptr => ps%leafp_storage_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr1 = dwt_ptr1 + change_state
             else
                ptr = 0._r8
                dwt_ptr1 = dwt_ptr1 - init_state
             end if
             
             ! leafp_xfer  
             ptr => ps%leafp_xfer_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr1 = dwt_ptr1 + change_state
             else
                ptr = 0._r8
                dwt_ptr1 = dwt_ptr1 - init_state
             end if
             
             ! frootp 
             ptr => ps%frootp_patch(p)
             dwt_ptr0 => dwt_frootp_to_litter(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr0 = dwt_ptr0 - change_state
             else
                ptr = 0._r8
                dwt_ptr0 = dwt_ptr0 + init_state
             end if
             
             ! frootp_storage 
             ptr => ps%frootp_storage_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr1 = dwt_ptr1 + change_state
             else
                ptr = 0._r8
                dwt_ptr1 = dwt_ptr1 - init_state
             end if
             
             ! frootp_xfer  
             ptr => ps%frootp_xfer_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr1 = dwt_ptr1 + change_state
             else
                ptr = 0._r8
                dwt_ptr1 = dwt_ptr1 - init_state
             end if
             
             ! livestemp  
             ptr => ps%livestemp_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr1 = dwt_ptr1 + change_state
             else
                ptr = 0._r8
                dwt_ptr1 = dwt_ptr1 - init_state
             end if
             
             ! livestemp_storage 
             ptr => ps%livestemp_storage_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr1 = dwt_ptr1 + change_state
             else
                ptr = 0._r8
                dwt_ptr1 = dwt_ptr1 - init_state
             end if
             
             ! livestemp_xfer 
             ptr => ps%livestemp_xfer_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr1 = dwt_ptr1 + change_state
             else
                ptr = 0._r8
                dwt_ptr1 = dwt_ptr1 - init_state
             end if
             
             ! deadstemp 
             ptr => ps%deadstemp_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr1 = dwt_ptr1 + change_state*pconv(pft%itype(p))
                dwt_ptr2 = dwt_ptr2 + change_state*pprod10(pft%itype(p))
                dwt_ptr3 = dwt_ptr3 + change_state*pprod100(pft%itype(p))
             else
                ptr = 0._r8
                dwt_ptr1 = dwt_ptr1 - init_state*pconv(pft%itype(p))
                dwt_ptr2 = dwt_ptr2 - init_state*pprod10(pft%itype(p))
                dwt_ptr3 = dwt_ptr3 - init_state*pprod100(pft%itype(p))
             end if
             
             ! deadstemp_storage 
             ptr => ps%deadstemp_storage_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr1 = dwt_ptr1 + change_state
             else
                ptr = 0._r8
                dwt_ptr1 = dwt_ptr1 - init_state
             end if
             
             ! deadstemp_xfer 
             ptr => ps%deadstemp_xfer_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr1 = dwt_ptr1 + change_state
             else
                ptr = 0._r8
                dwt_ptr1 = dwt_ptr1 - init_state
             end if
             
             ! livecrootp 
             ptr => ps%livecrootp_patch(p)
             dwt_ptr0 => dwt_livecrootp_to_litter(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr0 = dwt_ptr0 - change_state
             else
                ptr = 0._r8
                dwt_ptr0 = dwt_ptr0 + init_state
             end if
             
             ! livecrootp_storage  
             ptr => ps%livecrootp_storage_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr1 = dwt_ptr1 + change_state
             else
                ptr = 0._r8
                dwt_ptr1 = dwt_ptr1 - init_state
             end if
             
             ! livecrootp_xfer  
             ptr => ps%livecrootp_xfer_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr1 = dwt_ptr1 + change_state
             else
                ptr = 0._r8
                dwt_ptr1 = dwt_ptr1 - init_state
             end if
             
             ! deadcrootp 
             ptr => ps%deadcrootp_patch(p)
             dwt_ptr0 => dwt_deadcrootp_to_litter(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr0 = dwt_ptr0 - change_state
             else
                ptr = 0._r8
                dwt_ptr0 = dwt_ptr0 + init_state
             end if
             
             ! deadcrootp_storage  
             ptr => ps%deadcrootp_storage_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr1 = dwt_ptr1 + change_state
             else
                ptr = 0._r8
                dwt_ptr1 = dwt_ptr1 - init_state
             end if
             
             ! deadcrootp_xfer  
             ptr => ps%deadcrootp_xfer_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr1 = dwt_ptr1 + change_state
             else
                ptr = 0._r8
                dwt_ptr1 = dwt_ptr1 - init_state
             end if
             
             ! retransp  
             ptr => ps%retransp_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr1 = dwt_ptr1 + change_state
             else
                ptr = 0._r8
                dwt_ptr1 = dwt_ptr1 - init_state
             end if
             
             ! ppool  
             ptr => ps%ppool_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr1 = dwt_ptr1 + change_state
             else
                ptr = 0._r8
                dwt_ptr1 = dwt_ptr1 - init_state
             end if
             
             ! pft_ptrunc  
             ptr => ps%ptrunc_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr1 = dwt_ptr1 + change_state
             else
                ptr = 0._r8
                dwt_ptr1 = dwt_ptr1 - init_state
             end if
             
          end if       ! weight decreasing
       end if           ! is soil
    end do               ! patch loop
    
    ! calculate column-level seeding fluxes
    do pi = 1,max_patch_per_col
       do c = bounds%begc, bounds%endc
          if ( pi <=  col%npfts(c) ) then
             p = col%pfti(c) + pi - 1
             
             ! C fluxes
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

             ! P fluxes
             pf%dwt_seedp_to_leaf_col(c) = pf%dwt_seedp_to_leaf_col(c) + dwt_leafp_seed(p)/dt
             pf%dwt_seedp_to_deadstem_col(c) = pf%dwt_seedp_to_deadstem_col(c) &
                  + dwt_deadstemp_seed(p)/dt
          end if
       end do
    end do
    
    
    ! calculate pft-to-column for fluxes into litter and CWD pools
    do j = 1, nlevdecomp
       do pi = 1,max_patch_per_col
          do c = bounds%begc, bounds%endc
             if ( pi <=  col%npfts(c) ) then
                p = col%pfti(c) + pi - 1
                
                ! fine root litter carbon fluxes
                cf%dwt_frootc_to_litr_met_c_col(c,j) = &
                     cf%dwt_frootc_to_litr_met_c_col(c,j) + &
                     (dwt_frootc_to_litter(p)*ecophyscon%fr_flab(pft%itype(p)))/dt * cnstate_vars%froot_prof_patch(p,j)

                cf%dwt_frootc_to_litr_cel_c_col(c,j) = &
                     cf%dwt_frootc_to_litr_cel_c_col(c,j) + &
                     (dwt_frootc_to_litter(p)*ecophyscon%fr_fcel(pft%itype(p)))/dt * cnstate_vars%froot_prof_patch(p,j)

                cf%dwt_frootc_to_litr_lig_c_col(c,j) = &
                     cf%dwt_frootc_to_litr_lig_c_col(c,j) + &
                     (dwt_frootc_to_litter(p)*ecophyscon%fr_flig(pft%itype(p)))/dt * cnstate_vars%froot_prof_patch(p,j)
                
                
                ! fine root litter nitrogen fluxes
                nf%dwt_frootn_to_litr_met_n_col(c,j) = &
                     nf%dwt_frootn_to_litr_met_n_col(c,j) + &
                     (dwt_frootn_to_litter(p)*ecophyscon%fr_flab(pft%itype(p)))/dt * cnstate_vars%froot_prof_patch(p,j)
                nf%dwt_frootn_to_litr_cel_n_col(c,j) = &

                     nf%dwt_frootn_to_litr_cel_n_col(c,j) + &
                     (dwt_frootn_to_litter(p)*ecophyscon%fr_fcel(pft%itype(p)))/dt * cnstate_vars%froot_prof_patch(p,j)

                nf%dwt_frootn_to_litr_lig_n_col(c,j) = &
                     nf%dwt_frootn_to_litr_lig_n_col(c,j) + &
                     (dwt_frootn_to_litter(p)*ecophyscon%fr_flig(pft%itype(p)))/dt * cnstate_vars%froot_prof_patch(p,j)
                

                ! fine root litter phosphorus fluxes
                pf%dwt_frootp_to_litr_met_p_col(c,j) = &
                     pf%dwt_frootp_to_litr_met_p_col(c,j) + &
                     (dwt_frootp_to_litter(p)*ecophyscon%fr_flab(pft%itype(p)))/dt * cnstate_vars%froot_prof_patch(p,j)
                pf%dwt_frootp_to_litr_cel_p_col(c,j) = &

                     pf%dwt_frootp_to_litr_cel_p_col(c,j) + &
                     (dwt_frootp_to_litter(p)*ecophyscon%fr_fcel(pft%itype(p)))/dt * cnstate_vars%froot_prof_patch(p,j)

                pf%dwt_frootp_to_litr_lig_p_col(c,j) = &
                     pf%dwt_frootp_to_litr_lig_p_col(c,j) + &
                     (dwt_frootp_to_litter(p)*ecophyscon%fr_flig(pft%itype(p)))/dt * cnstate_vars%froot_prof_patch(p,j)

                ! livecroot fluxes to cwd
                cf%dwt_livecrootc_to_cwdc_col(c,j) = &
                     cf%dwt_livecrootc_to_cwdc_col(c,j) + &
                     (dwt_livecrootc_to_litter(p))/dt * cnstate_vars%croot_prof_patch(p,j)

                nf%dwt_livecrootn_to_cwdn_col(c,j) = &
                     nf%dwt_livecrootn_to_cwdn_col(c,j) + &
                     (dwt_livecrootn_to_litter(p))/dt * cnstate_vars%croot_prof_patch(p,j)
                
                pf%dwt_livecrootp_to_cwdp_col(c,j) = &
                     pf%dwt_livecrootp_to_cwdp_col(c,j) + &
                     (dwt_livecrootp_to_litter(p))/dt * cnstate_vars%croot_prof_patch(p,j)

                ! deadcroot fluxes to cwd
                cf%dwt_deadcrootc_to_cwdc_col(c,j) = &
                     cf%dwt_deadcrootc_to_cwdc_col(c,j) + &
                     (dwt_deadcrootc_to_litter(p))/dt * cnstate_vars%croot_prof_patch(p,j)

                nf%dwt_deadcrootn_to_cwdn_col(c,j) = &
                     nf%dwt_deadcrootn_to_cwdn_col(c,j) + &
                     (dwt_deadcrootn_to_litter(p))/dt * cnstate_vars%croot_prof_patch(p,j)
             
                pf%dwt_deadcrootp_to_cwdp_col(c,j) = &
                     pf%dwt_deadcrootp_to_cwdp_col(c,j) + &
                     (dwt_deadcrootp_to_litter(p))/dt * cnstate_vars%croot_prof_patch(p,j)

                if ( use_c13 ) then
                   ! C13 fine root litter fluxes
                   c13_cf%dwt_frootc_to_litr_met_c_col(c,j) = &
                        c13_cf%dwt_frootc_to_litr_met_c_col(c,j) + &
                        (dwt_frootc13_to_litter(p)*ecophyscon%fr_flab(pft%itype(p)))/dt * cnstate_vars%froot_prof_patch(p,j)

                   c13_cf%dwt_frootc_to_litr_cel_c_col(c,j) = &
                        c13_cf%dwt_frootc_to_litr_cel_c_col(c,j) + &
                        (dwt_frootc13_to_litter(p)*ecophyscon%fr_fcel(pft%itype(p)))/dt * cnstate_vars%froot_prof_patch(p,j)

                   c13_cf%dwt_frootc_to_litr_lig_c_col(c,j) = &
                        c13_cf%dwt_frootc_to_litr_lig_c_col(c,j) + &
                        (dwt_frootc13_to_litter(p)*ecophyscon%fr_flig(pft%itype(p)))/dt * cnstate_vars%froot_prof_patch(p,j)

                   ! livecroot fluxes to cwd
                   c13_cf%dwt_livecrootc_to_cwdc_col(c,j) = &
                        c13_cf%dwt_livecrootc_to_cwdc_col(c,j) + &
                        (dwt_livecrootc13_to_litter(p))/dt * cnstate_vars%croot_prof_patch(p,j)

                   ! deadcroot fluxes to cwd
                   c13_cf%dwt_deadcrootc_to_cwdc_col(c,j) = &
                        c13_cf%dwt_deadcrootc_to_cwdc_col(c,j) + &
                        (dwt_deadcrootc13_to_litter(p))/dt * cnstate_vars%croot_prof_patch(p,j)
                   
                endif
                
                if ( use_c14 ) then                   
                   ! C14 fine root litter fluxes
                   c14_cf%dwt_frootc_to_litr_met_c_col(c,j) = &
                        c14_cf%dwt_frootc_to_litr_met_c_col(c,j) + &
                        (dwt_frootc14_to_litter(p)*ecophyscon%fr_flab(pft%itype(p)))/dt * cnstate_vars%froot_prof_patch(p,j)

                   c14_cf%dwt_frootc_to_litr_cel_c_col(c,j) = &
                        c14_cf%dwt_frootc_to_litr_cel_c_col(c,j) + &
                        (dwt_frootc14_to_litter(p)*ecophyscon%fr_fcel(pft%itype(p)))/dt * cnstate_vars%froot_prof_patch(p,j)

                   c14_cf%dwt_frootc_to_litr_lig_c_col(c,j) = &
                        c14_cf%dwt_frootc_to_litr_lig_c_col(c,j) + &
                        (dwt_frootc14_to_litter(p)*ecophyscon%fr_flig(pft%itype(p)))/dt * cnstate_vars%froot_prof_patch(p,j)

                   ! livecroot fluxes to cwd
                   c14_cf%dwt_livecrootc_to_cwdc_col(c,j) = &
                        c14_cf%dwt_livecrootc_to_cwdc_col(c,j) + &
                        (dwt_livecrootc14_to_litter(p))/dt * cnstate_vars%croot_prof_patch(p,j)

                   ! deadcroot fluxes to cwd
                   c14_cf%dwt_deadcrootc_to_cwdc_col(c,j) = &
                        c14_cf%dwt_deadcrootc_to_cwdc_col(c,j) + &
                        (dwt_deadcrootc14_to_litter(p))/dt * cnstate_vars%croot_prof_patch(p,j)
                endif
                
             end if
          end do
       end do
    end do
    ! calculate pft-to-column for fluxes into product pools and conversion flux
    do pi = 1,max_patch_per_col
       do c = bounds%begc,bounds%endc
          if (pi <= col%npfts(c)) then
             p = col%pfti(c) + pi - 1
             
             ! column-level fluxes are accumulated as positive fluxes.
             ! column-level C flux updates
             cf%dwt_conv_cflux_col(c) = cf%dwt_conv_cflux_col(c) - conv_cflux(p)/dt
             cf%dwt_prod10c_gain_col(c) = cf%dwt_prod10c_gain_col(c) - prod10_cflux(p)/dt
             cf%dwt_prod100c_gain_col(c) = cf%dwt_prod100c_gain_col(c) - prod100_cflux(p)/dt

             ! These magic constants should be replaced with: nbrdlf_evr_trp_tree and nbrdlf_dcd_trp_tree
             if(pft%itype(p)==4.or.pft%itype(p)==6)then
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

end module dynConsBiogeochemMod
