module dynConsBiogeochemMod

  !---------------------------------------------------------------------------
  !
  ! !DESCRIPTION:
  ! Handle conservation of biogeochemical quantities (C & N) with dynamic land cover.
  !
  ! !USES:
  use shr_kind_mod                 , only : r8 => shr_kind_r8
  use shr_log_mod                  , only : errMsg => shr_log_errMsg
  use decompMod                    , only : bounds_type
  use abortutils                   , only : endrun
  use clm_varctl                   , only : iulog, use_c13, use_c14
  use pftconMod                    , only : pftcon
  use CanopyStateType              , only : canopystate_type
  use PhotosynthesisMod            , only : photosyns_type
  use CNVegStateType               , only : cnveg_state_type
  use CNVegCarbonStateType         , only : cnveg_carbonstate_type
  use CNVegCarbonFluxType          , only : cnveg_carbonflux_type
  use CNVegNitrogenStateType       , only : cnveg_nitrogenstate_type
  use CNVegNitrogenFluxType        , only : cnveg_nitrogenflux_type
  use SoilBiogeochemStateType      , only : soilBiogeochem_state_type
  use SoilBiogeochemCarbonFluxType , only : soilBiogeochem_carbonflux_type
  use LandunitType                 , only : lun                
  use ColumnType                   , only : col                
  use PatchType                    , only : patch                
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  !
  public :: dyn_cnbal_patch
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine dyn_cnbal_patch(bounds, prior_weights,                                       &
       canopystate_inst, photosyns_inst, cnveg_state_inst,                                &
       cnveg_carbonstate_inst, c13_cnveg_carbonstate_inst, c14_cnveg_carbonstate_inst,    &
       cnveg_carbonflux_inst, c13_cnveg_carbonflux_inst, c14_cnveg_carbonflux_inst,       &
       cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst, soilbiogeochem_carbonflux_inst, &
       soilbiogeochem_state_inst)
    !
    ! !DESCRIPTION:
    ! Modify patch-level state and flux variables to maintain carbon and nitrogen balance with
    ! dynamic patch-weights.
    !
    ! !USES:
    use shr_const_mod      , only : SHR_CONST_PDB
    use landunit_varcon    , only : istsoil, istcrop
    use clm_varpar         , only : numveg, nlevdecomp, max_patch_per_col
    use clm_varcon         , only : c13ratio, c14ratio
    use clm_time_manager   , only : get_step_size
    use dynPriorWeightsMod , only : prior_weights_type
    !
    ! !ARGUMENTS:
    type(bounds_type)                    , intent(in)    :: bounds        
    type(prior_weights_type)             , intent(in)    :: prior_weights ! weights prior to the subgrid weight updates
    type(canopystate_type)               , intent(inout) :: canopystate_inst
    type(photosyns_type)                 , intent(inout) :: photosyns_inst
    type(cnveg_state_type)               , intent(inout) :: cnveg_state_inst
    type(cnveg_carbonstate_type)         , intent(inout) :: cnveg_carbonstate_inst
    type(cnveg_carbonstate_type)         , intent(inout) :: c13_cnveg_carbonstate_inst
    type(cnveg_carbonstate_type)         , intent(inout) :: c14_cnveg_carbonstate_inst
    type(cnveg_carbonflux_type)          , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_carbonflux_type)          , intent(inout) :: c13_cnveg_carbonflux_inst
    type(cnveg_carbonflux_type)          , intent(inout) :: c14_cnveg_carbonflux_inst
    type(cnveg_nitrogenstate_type)       , intent(inout) :: cnveg_nitrogenstate_inst
    type(cnveg_nitrogenflux_type)        , intent(inout) :: cnveg_nitrogenflux_inst
    type(soilbiogeochem_carbonflux_type) , intent(inout) :: soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_state_type)      , intent(in)    :: soilbiogeochem_state_inst
    !
    ! !LOCAL VARIABLES:
    integer                       :: pi,p,c,l,g,j                  ! indices
    integer                       :: ier                           ! error code
    real(r8)                      :: dwt                           ! change in patch weight (relative to column)
    real(r8)                      :: dt                            ! land model time step (sec)
    real(r8)                      :: init_h2ocan                   ! initial canopy water mass
    real(r8)                      :: new_h2ocan                    ! canopy water mass after weight shift
    real(r8), allocatable         :: dwt_leafc_seed(:)             ! patch-level mass gain due to seeding of new area
    real(r8), allocatable         :: dwt_leafn_seed(:)             ! patch-level mass gain due to seeding of new area
    real(r8), allocatable         :: dwt_deadstemc_seed(:)         ! patch-level mass gain due to seeding of new area
    real(r8), allocatable         :: dwt_deadstemn_seed(:)         ! patch-level mass gain due to seeding of new area
    real(r8), allocatable         :: dwt_frootc_to_litter(:)       ! patch-level mass loss due to weight shift
    real(r8), allocatable         :: dwt_livecrootc_to_litter(:)   ! patch-level mass loss due to weight shift
    real(r8), allocatable         :: dwt_deadcrootc_to_litter(:)   ! patch-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_frootn_to_litter(:)       ! patch-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_livecrootn_to_litter(:)   ! patch-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_deadcrootn_to_litter(:)   ! patch-level mass loss due to weight shift
    real(r8), allocatable         :: conv_cflux(:)                 ! patch-level mass loss due to weight shift
    real(r8), allocatable         :: prod10_cflux(:)               ! patch-level mass loss due to weight shift
    real(r8), allocatable         :: prod100_cflux(:)              ! patch-level mass loss due to weight shift
    real(r8), allocatable, target :: conv_nflux(:)                 ! patch-level mass loss due to weight shift
    real(r8), allocatable, target :: prod10_nflux(:)               ! patch-level mass loss due to weight shift
    real(r8), allocatable, target :: prod100_nflux(:)              ! patch-level mass loss due to weight shift
    real(r8)                      :: t1,t2,wt_new,wt_old
    real(r8)                      :: init_state, change_state, new_state
    real(r8)                      :: tot_leaf, pleaf, pstor, pxfer
    real(r8)                      :: leafc_seed, leafn_seed
    real(r8)                      :: deadstemc_seed, deadstemn_seed
    real(r8), pointer             :: dwt_ptr0, dwt_ptr1, dwt_ptr2, dwt_ptr3, ptr
    character(len=32)             :: subname='dyn_cbal'            ! subroutine name
    !! C13
    real(r8), allocatable         :: dwt_leafc13_seed(:)           ! patch-level mass gain due to seeding of new area
    real(r8), allocatable         :: dwt_deadstemc13_seed(:)       ! patch-level mass gain due to seeding of new area
    real(r8), allocatable, target :: dwt_frootc13_to_litter(:)     ! patch-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_livecrootc13_to_litter(:) ! patch-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_deadcrootc13_to_litter(:) ! patch-level mass loss due to weight shift
    real(r8), allocatable, target :: conv_c13flux(:)               ! patch-level mass loss due to weight shift
    real(r8), allocatable, target :: prod10_c13flux(:)             ! patch-level mass loss due to weight shift
    real(r8), allocatable, target :: prod100_c13flux(:)            ! patch-level mass loss due to weight shift
    real(r8)                      :: c3_del13c                     ! typical del13C for C3 photosynthesis (permil, relative to PDB)
    real(r8)                      :: c4_del13c                     ! typical del13C for C4 photosynthesis (permil, relative to PDB)
    real(r8)                      :: c3_r1_c13                     ! isotope ratio (13c/12c) for C3 photosynthesis
    real(r8)                      :: c4_r1_c13                     ! isotope ratio (13c/12c) for C4 photosynthesis
    real(r8)                      :: c3_r2_c13                     ! isotope ratio (13c/[12c+13c]) for C3 photosynthesis
    real(r8)                      :: c4_r2_c13                     ! isotope ratio (13c/[12c+13c]) for C4 photosynthesis
    real(r8)                      :: leafc13_seed, deadstemc13_seed
    !! C14
    real(r8), allocatable         :: dwt_leafc14_seed(:)           ! patch-level mass gain due to seeding of new area
    real(r8), allocatable         :: dwt_deadstemc14_seed(:)       ! patch-level mass gain due to seeding of new area
    real(r8), allocatable, target :: dwt_frootc14_to_litter(:)     ! patch-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_livecrootc14_to_litter(:) ! patch-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_deadcrootc14_to_litter(:) ! patch-level mass loss due to weight shift
    real(r8), allocatable, target :: conv_c14flux(:)               ! patch-level mass loss due to weight shift
    real(r8), allocatable, target :: prod10_c14flux(:)             ! patch-level mass loss due to weight shift
    real(r8), allocatable, target :: prod100_c14flux(:)            ! patch-level mass loss due to weight shift
    real(r8)                      :: c3_del14c                     ! typical del14C for C3 photosynthesis (permil, relative to PDB)
    real(r8)                      :: c4_del14c                     ! typical del14C for C4 photosynthesis (permil, relative to PDB)
    real(r8)                      :: c3_r1_c14                     ! isotope ratio (14c/12c) for C3 photosynthesis
    real(r8)                      :: c4_r1_c14                     ! isotope ratio (14c/12c) for C4 photosynthesis
    real(r8)                      :: c3_r2_c14                     ! isotope ratio (14c/[12c+14c]) for C3 photosynthesis
    real(r8)                      :: c4_r2_c14                     ! isotope ratio (14c/[12c+14c]) for C4 photosynthesis
    real(r8)                      :: leafc14_seed, deadstemc14_seed
    !-----------------------------------------------------------------------
    
    ! Allocate patch-level mass loss arrays
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
       c = patch%column(p)
       ! initialize all the patch-level local flux arrays
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
       
       l = patch%landunit(p)
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          
          ! calculate the change in weight for the timestep
          dwt = patch%wtcol(p)-prior_weights%pwtcol(p)
          CNveg_state_inst%lfpftd_patch(p) = -dwt

          ! Patches for which weight increases on this timestep
          if (dwt > 0._r8) then
             
             ! first identify Patches that are initiating on this timestep
             ! and set all the necessary state and flux variables
             if (prior_weights%pwtcol(p) == 0._r8) then
                
                ! set initial conditions for PFT that is being initiated
                ! in this time step.  Based on the settings in cnIniTimeVar.
                
                ! patch-level carbon state variables
                cnveg_carbonstate_inst%leafc_patch(p)              = 0._r8
                cnveg_carbonstate_inst%leafc_storage_patch(p)      = 0._r8
                cnveg_carbonstate_inst%leafc_xfer_patch(p)         = 0._r8
                cnveg_carbonstate_inst%frootc_patch(p)             = 0._r8
                cnveg_carbonstate_inst%frootc_storage_patch(p)     = 0._r8
                cnveg_carbonstate_inst%frootc_xfer_patch(p)        = 0._r8
                cnveg_carbonstate_inst%livestemc_patch(p)          = 0._r8
                cnveg_carbonstate_inst%livestemc_storage_patch(p)  = 0._r8
                cnveg_carbonstate_inst%livestemc_xfer_patch(p)     = 0._r8
                cnveg_carbonstate_inst%deadstemc_patch(p)          = 0._r8
                cnveg_carbonstate_inst%deadstemc_storage_patch(p)  = 0._r8
                cnveg_carbonstate_inst%deadstemc_xfer_patch(p)     = 0._r8
                cnveg_carbonstate_inst%livecrootc_patch(p)         = 0._r8
                cnveg_carbonstate_inst%livecrootc_storage_patch(p) = 0._r8
                cnveg_carbonstate_inst%livecrootc_xfer_patch(p)    = 0._r8
                cnveg_carbonstate_inst%deadcrootc_patch(p)         = 0._r8
                cnveg_carbonstate_inst%deadcrootc_storage_patch(p) = 0._r8
                cnveg_carbonstate_inst%deadcrootc_xfer_patch(p)    = 0._r8
                cnveg_carbonstate_inst%gresp_storage_patch(p)      = 0._r8
                cnveg_carbonstate_inst%gresp_xfer_patch(p)         = 0._r8
                cnveg_carbonstate_inst%cpool_patch(p)              = 0._r8
                cnveg_carbonstate_inst%xsmrpool_patch(p)           = 0._r8
                cnveg_carbonstate_inst%ctrunc_patch(p)             = 0._r8
                cnveg_carbonstate_inst%dispvegc_patch(p)           = 0._r8
                cnveg_carbonstate_inst%storvegc_patch(p)           = 0._r8
                cnveg_carbonstate_inst%totc_patch(p)               = 0._r8
                cnveg_carbonstate_inst%totvegc_patch(p)            = 0._r8
                
                if ( use_c13 ) then
                   ! patch-level carbon-13 state variables
                   c13_cnveg_carbonstate_inst%leafc_patch(p)              = 0._r8
                   c13_cnveg_carbonstate_inst%leafc_storage_patch(p)      = 0._r8
                   c13_cnveg_carbonstate_inst%leafc_xfer_patch(p)         = 0._r8
                   c13_cnveg_carbonstate_inst%frootc_patch(p)             = 0._r8
                   c13_cnveg_carbonstate_inst%frootc_storage_patch(p)     = 0._r8
                   c13_cnveg_carbonstate_inst%frootc_xfer_patch(p)        = 0._r8
                   c13_cnveg_carbonstate_inst%livestemc_patch(p)          = 0._r8
                   c13_cnveg_carbonstate_inst%livestemc_storage_patch(p)  = 0._r8
                   c13_cnveg_carbonstate_inst%livestemc_xfer_patch(p)     = 0._r8
                   c13_cnveg_carbonstate_inst%deadstemc_patch(p)          = 0._r8
                   c13_cnveg_carbonstate_inst%deadstemc_storage_patch(p)  = 0._r8
                   c13_cnveg_carbonstate_inst%deadstemc_xfer_patch(p)     = 0._r8
                   c13_cnveg_carbonstate_inst%livecrootc_patch(p)         = 0._r8
                   c13_cnveg_carbonstate_inst%livecrootc_storage_patch(p) = 0._r8
                   c13_cnveg_carbonstate_inst%livecrootc_xfer_patch(p)    = 0._r8
                   c13_cnveg_carbonstate_inst%deadcrootc_patch(p)         = 0._r8
                   c13_cnveg_carbonstate_inst%deadcrootc_storage_patch(p) = 0._r8
                   c13_cnveg_carbonstate_inst%deadcrootc_xfer_patch(p)    = 0._r8
                   c13_cnveg_carbonstate_inst%gresp_storage_patch(p)      = 0._r8
                   c13_cnveg_carbonstate_inst%gresp_xfer_patch(p)         = 0._r8
                   c13_cnveg_carbonstate_inst%cpool_patch(p)              = 0._r8
                   c13_cnveg_carbonstate_inst%xsmrpool_patch(p)           = 0._r8
                   c13_cnveg_carbonstate_inst%ctrunc_patch(p)             = 0._r8
                   c13_cnveg_carbonstate_inst%dispvegc_patch(p)           = 0._r8
                   c13_cnveg_carbonstate_inst%storvegc_patch(p)           = 0._r8
                   c13_cnveg_carbonstate_inst%totc_patch(p)               = 0._r8
                   c13_cnveg_carbonstate_inst%totvegc_patch(p)            = 0._r8
                endif
                
                if ( use_c14 ) then
                   ! patch-level carbon-14 state variables
                   c14_cnveg_carbonstate_inst%leafc_patch(p)              = 0._r8
                   c14_cnveg_carbonstate_inst%leafc_storage_patch(p)      = 0._r8
                   c14_cnveg_carbonstate_inst%leafc_xfer_patch(p)         = 0._r8
                   c14_cnveg_carbonstate_inst%frootc_patch(p)             = 0._r8
                   c14_cnveg_carbonstate_inst%frootc_storage_patch(p)     = 0._r8
                   c14_cnveg_carbonstate_inst%frootc_xfer_patch(p)        = 0._r8
                   c14_cnveg_carbonstate_inst%livestemc_patch(p)          = 0._r8
                   c14_cnveg_carbonstate_inst%livestemc_storage_patch(p)  = 0._r8
                   c14_cnveg_carbonstate_inst%livestemc_xfer_patch(p)     = 0._r8
                   c14_cnveg_carbonstate_inst%deadstemc_patch(p)          = 0._r8
                   c14_cnveg_carbonstate_inst%deadstemc_storage_patch(p)  = 0._r8
                   c14_cnveg_carbonstate_inst%deadstemc_xfer_patch(p)     = 0._r8
                   c14_cnveg_carbonstate_inst%livecrootc_patch(p)         = 0._r8
                   c14_cnveg_carbonstate_inst%livecrootc_storage_patch(p) = 0._r8
                   c14_cnveg_carbonstate_inst%livecrootc_xfer_patch(p)    = 0._r8
                   c14_cnveg_carbonstate_inst%deadcrootc_patch(p)         = 0._r8
                   c14_cnveg_carbonstate_inst%deadcrootc_storage_patch(p) = 0._r8
                   c14_cnveg_carbonstate_inst%deadcrootc_xfer_patch(p)    = 0._r8
                   c14_cnveg_carbonstate_inst%gresp_storage_patch(p)      = 0._r8
                   c14_cnveg_carbonstate_inst%gresp_xfer_patch(p)         = 0._r8
                   c14_cnveg_carbonstate_inst%cpool_patch(p)              = 0._r8
                   c14_cnveg_carbonstate_inst%xsmrpool_patch(p)           = 0._r8
                   c14_cnveg_carbonstate_inst%ctrunc_patch(p)             = 0._r8
                   c14_cnveg_carbonstate_inst%dispvegc_patch(p)           = 0._r8
                   c14_cnveg_carbonstate_inst%storvegc_patch(p)           = 0._r8
                   c14_cnveg_carbonstate_inst%totc_patch(p)               = 0._r8
                   c14_cnveg_carbonstate_inst%totvegc_patch(p)            = 0._r8
                endif
                
                ! patch-level nitrogen state variables
                cnveg_nitrogenstate_inst%leafn_patch(p)	             = 0._r8
                cnveg_nitrogenstate_inst%leafn_storage_patch(p)      = 0._r8
                cnveg_nitrogenstate_inst%leafn_xfer_patch(p)         = 0._r8
                cnveg_nitrogenstate_inst%frootn_patch(p)	     = 0._r8
                cnveg_nitrogenstate_inst%frootn_storage_patch(p)     = 0._r8
                cnveg_nitrogenstate_inst%frootn_xfer_patch(p)        = 0._r8
                cnveg_nitrogenstate_inst%livestemn_patch(p)          = 0._r8
                cnveg_nitrogenstate_inst%livestemn_storage_patch(p)  = 0._r8
                cnveg_nitrogenstate_inst%livestemn_xfer_patch(p)     = 0._r8
                cnveg_nitrogenstate_inst%deadstemn_patch(p)          = 0._r8
                cnveg_nitrogenstate_inst%deadstemn_storage_patch(p)  = 0._r8
                cnveg_nitrogenstate_inst%deadstemn_xfer_patch(p)     = 0._r8
                cnveg_nitrogenstate_inst%livecrootn_patch(p)         = 0._r8
                cnveg_nitrogenstate_inst%livecrootn_storage_patch(p) = 0._r8
                cnveg_nitrogenstate_inst%livecrootn_xfer_patch(p)    = 0._r8
                cnveg_nitrogenstate_inst%deadcrootn_patch(p)         = 0._r8
                cnveg_nitrogenstate_inst%deadcrootn_storage_patch(p) = 0._r8
                cnveg_nitrogenstate_inst%deadcrootn_xfer_patch(p)    = 0._r8
                cnveg_nitrogenstate_inst%retransn_patch(p)           = 0._r8
                cnveg_nitrogenstate_inst%npool_patch(p)	             = 0._r8
                cnveg_nitrogenstate_inst%ntrunc_patch(p)             = 0._r8
                cnveg_nitrogenstate_inst%dispvegn_patch(p)           = 0._r8
                cnveg_nitrogenstate_inst%storvegn_patch(p)           = 0._r8
                cnveg_nitrogenstate_inst%totvegn_patch(p)            = 0._r8
                cnveg_nitrogenstate_inst%totn_patch (p)              = 0._r8
                
                ! initialize same flux and epv variables that are set
                canopystate_inst%laisun_patch(p) = 0._r8
                canopystate_inst%laisha_patch(p) = 0._r8
                
                cnveg_state_inst%dormant_flag_patch(p)          = 1._r8
                cnveg_state_inst%days_active_patch(p)           = 0._r8
                cnveg_state_inst%onset_flag_patch(p)            = 0._r8
                cnveg_state_inst%onset_counter_patch(p)         = 0._r8
                cnveg_state_inst%onset_gddflag_patch(p)         = 0._r8
                cnveg_state_inst%onset_fdd_patch(p)             = 0._r8
                cnveg_state_inst%onset_gdd_patch(p)             = 0._r8
                cnveg_state_inst%onset_swi_patch(p)             = 0._r8
                cnveg_state_inst%offset_flag_patch(p)           = 0._r8
                cnveg_state_inst%offset_counter_patch(p)        = 0._r8
                cnveg_state_inst%offset_fdd_patch(p)            = 0._r8
                cnveg_state_inst%offset_swi_patch(p)            = 0._r8
                cnveg_state_inst%lgsf_patch(p)                  = 0._r8
                cnveg_state_inst%bglfr_patch(p)                 = 0._r8
                cnveg_state_inst%bgtr_patch(p)                  = 0._r8
                cnveg_state_inst%annavg_t2m_patch(p)            = cnveg_state_inst%annavg_t2m_col(c)
                cnveg_state_inst%tempavg_t2m_patch(p)           = 0._r8
                cnveg_state_inst%c_allometry_patch(p)           = 0._r8
                cnveg_state_inst%n_allometry_patch(p)           = 0._r8
                cnveg_state_inst%tempsum_potential_gpp_patch(p) = 0._r8
                cnveg_state_inst%annsum_potential_gpp_patch(p)  = 0._r8
                cnveg_state_inst%tempmax_retransn_patch(p)      = 0._r8
                cnveg_state_inst%annmax_retransn_patch(p)       = 0._r8
                cnveg_state_inst%downreg_patch(p)               = 0._r8

                cnveg_carbonflux_inst%xsmrpool_recover_patch(p)      = 0._r8
                cnveg_carbonflux_inst%plant_calloc_patch(p)          = 0._r8
                cnveg_carbonflux_inst%excess_cflux_patch(p)          = 0._r8
                cnveg_carbonflux_inst%prev_leafc_to_litter_patch(p)  = 0._r8
                cnveg_carbonflux_inst%prev_frootc_to_litter_patch(p) = 0._r8
                cnveg_carbonflux_inst%availc_patch(p)                = 0._r8
                cnveg_carbonflux_inst%gpp_before_downreg_patch(p)    = 0._r8

                cnveg_carbonflux_inst%tempsum_npp_patch(p)       = 0._r8
                cnveg_carbonflux_inst%annsum_npp_patch(p)        = 0._r8

                cnveg_nitrogenflux_inst%plant_ndemand_patch(p)       = 0._r8
                cnveg_nitrogenflux_inst%avail_retransn_patch(p)      = 0._r8
                cnveg_nitrogenflux_inst%plant_nalloc_patch(p)        = 0._r8
                
                if ( use_c13 ) then
                   c13_cnveg_carbonflux_inst%xsmrpool_c13ratio_patch(p) = c13ratio
                end if

                call photosyns_inst%NewPatchinit(p)
                
             end if  ! end initialization of new patch
             
             ! (still in dwt > 0 block)
             
             ! set the seed sources for leaf and deadstem
             ! leaf source is split later between leaf, leaf_storage, leaf_xfer
             leafc_seed   = 0._r8
             leafn_seed   = 0._r8
             deadstemc_seed   = 0._r8
             deadstemn_seed   = 0._r8
             if ( use_c13 ) then
                leafc13_seed = 0._r8
                deadstemc13_seed = 0._r8
             endif
             if ( use_c14 ) then
                leafc14_seed = 0._r8
                deadstemc14_seed = 0._r8
             endif
             if (patch%itype(p) /= 0) then
                leafc_seed = 1._r8
                leafn_seed  = leafc_seed / pftcon%leafcn(patch%itype(p))
                if (pftcon%woody(patch%itype(p)) == 1._r8) then
                   deadstemc_seed = 0.1_r8
                   deadstemn_seed = deadstemc_seed / pftcon%deadwdcn(patch%itype(p))
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
                   
                   if (pftcon%c3psn(patch%itype(p)) == 1._r8) then
                      leafc13_seed     = leafc_seed     * c3_r2_c13
                      deadstemc13_seed = deadstemc_seed * c3_r2_c13
                   else
                      leafc13_seed     = leafc_seed     * c4_r2_c13
                      deadstemc13_seed = deadstemc_seed * c4_r2_c13
                   end if
                endif
                
                if ( use_c14 ) then
                   ! 14c state is initialized assuming initial "modern" 14C of 1.e-12
                   if (pftcon%c3psn(patch%itype(p)) == 1._r8) then
                      leafc14_seed     = leafc_seed     * c14ratio
                      deadstemc14_seed = deadstemc_seed * c14ratio
                   else
                      leafc14_seed     = leafc_seed     * c14ratio
                      deadstemc14_seed = deadstemc_seed * c14ratio
                   end if
                endif
             end if
             
             ! When PATCH area expands (dwt > 0), the patch-level mass density 
             ! is modified to conserve the original patch mass distributed
             ! over the new (larger) area, plus a term to account for the 
             ! introduction of new seed source for leaf and deadstem
             t1 = prior_weights%pwtcol(p)/patch%wtcol(p)
             t2 = dwt/patch%wtcol(p)
             
             tot_leaf = cnveg_carbonstate_inst%leafc_patch(p) + &
                  cnveg_carbonstate_inst%leafc_storage_patch(p) + &
                  cnveg_carbonstate_inst%leafc_xfer_patch(p)
             pleaf = 0._r8
             pstor = 0._r8
             pxfer = 0._r8
             if (tot_leaf /= 0._r8) then
                ! when adding seed source to non-zero leaf state, use current proportions
                pleaf = cnveg_carbonstate_inst%leafc_patch(p)/tot_leaf
                pstor = cnveg_carbonstate_inst%leafc_storage_patch(p)/tot_leaf
                pxfer = cnveg_carbonstate_inst%leafc_xfer_patch(p)/tot_leaf
             else
                ! when initiating from zero leaf state, use evergreen flag to set proportions
                if (pftcon%evergreen(patch%itype(p)) == 1._r8) then
                   pleaf = 1._r8
                else
                   pstor = 1._r8
                end if
             end if
             cnveg_carbonstate_inst%leafc_patch(p)              = cnveg_carbonstate_inst%leafc_patch(p)              *t1 &
                  + leafc_seed*pleaf*t2
             cnveg_carbonstate_inst%leafc_storage_patch(p)      = cnveg_carbonstate_inst%leafc_storage_patch(p)      *t1 &
                  + leafc_seed*pstor*t2
             cnveg_carbonstate_inst%leafc_xfer_patch(p)         = cnveg_carbonstate_inst%leafc_xfer_patch(p)         *t1 &
                  + leafc_seed*pxfer*t2
             cnveg_carbonstate_inst%frootc_patch(p)             = cnveg_carbonstate_inst%frootc_patch(p)             *t1
             cnveg_carbonstate_inst%frootc_storage_patch(p)     = cnveg_carbonstate_inst%frootc_storage_patch(p)     *t1
             cnveg_carbonstate_inst%frootc_xfer_patch(p)        = cnveg_carbonstate_inst%frootc_xfer_patch(p)        *t1
             cnveg_carbonstate_inst%livestemc_patch(p)          = cnveg_carbonstate_inst%livestemc_patch(p)          *t1
             cnveg_carbonstate_inst%livestemc_storage_patch(p)  = cnveg_carbonstate_inst%livestemc_storage_patch(p)  *t1
             cnveg_carbonstate_inst%livestemc_xfer_patch(p)     = cnveg_carbonstate_inst%livestemc_xfer_patch(p)     *t1
             cnveg_carbonstate_inst%deadstemc_patch(p)          = cnveg_carbonstate_inst%deadstemc_patch(p)          *t1  &
                  + deadstemc_seed*t2
             cnveg_carbonstate_inst%deadstemc_storage_patch(p)  = cnveg_carbonstate_inst%deadstemc_storage_patch(p)  *t1
             cnveg_carbonstate_inst%deadstemc_xfer_patch(p)     = cnveg_carbonstate_inst%deadstemc_xfer_patch(p)     *t1
             cnveg_carbonstate_inst%livecrootc_patch(p)         = cnveg_carbonstate_inst%livecrootc_patch(p)         *t1
             cnveg_carbonstate_inst%livecrootc_storage_patch(p) = cnveg_carbonstate_inst%livecrootc_storage_patch(p) *t1
             cnveg_carbonstate_inst%livecrootc_xfer_patch(p)    = cnveg_carbonstate_inst%livecrootc_xfer_patch(p)    *t1
             cnveg_carbonstate_inst%deadcrootc_patch(p)         = cnveg_carbonstate_inst%deadcrootc_patch(p)         *t1
             cnveg_carbonstate_inst%deadcrootc_storage_patch(p) = cnveg_carbonstate_inst%deadcrootc_storage_patch(p) *t1
             cnveg_carbonstate_inst%deadcrootc_xfer_patch(p)    = cnveg_carbonstate_inst%deadcrootc_xfer_patch(p)    *t1
             cnveg_carbonstate_inst%gresp_storage_patch(p)      = cnveg_carbonstate_inst%gresp_storage_patch(p)      *t1
             cnveg_carbonstate_inst%gresp_xfer_patch(p)         = cnveg_carbonstate_inst%gresp_xfer_patch(p)         *t1
             cnveg_carbonstate_inst%cpool_patch(p)              = cnveg_carbonstate_inst%cpool_patch(p)              *t1
             cnveg_carbonstate_inst%xsmrpool_patch(p)           = cnveg_carbonstate_inst%xsmrpool_patch(p)           *t1
             cnveg_carbonstate_inst%ctrunc_patch(p)             = cnveg_carbonstate_inst%ctrunc_patch(p)             *t1
             cnveg_carbonstate_inst%dispvegc_patch(p)           = cnveg_carbonstate_inst%dispvegc_patch(p)           *t1
             cnveg_carbonstate_inst%storvegc_patch(p)           = cnveg_carbonstate_inst%storvegc_patch(p)           *t1
             cnveg_carbonstate_inst%totc_patch(p)               = cnveg_carbonstate_inst%totc_patch(p)               *t1
             cnveg_carbonstate_inst%totvegc_patch(p)            = cnveg_carbonstate_inst%totvegc_patch(p)            *t1
             
             if ( use_c13 ) then
                ! patch-level carbon-13 state variables 
                tot_leaf = &
                     c13_cnveg_carbonstate_inst%leafc_patch(p) + &
                     c13_cnveg_carbonstate_inst%leafc_storage_patch(p) + &
                     c13_cnveg_carbonstate_inst%leafc_xfer_patch(p)
                pleaf = 0._r8
                pstor = 0._r8
                pxfer = 0._r8
                if (tot_leaf /= 0._r8) then
                   pleaf = c13_cnveg_carbonstate_inst%leafc_patch(p)/tot_leaf
                   pstor = c13_cnveg_carbonstate_inst%leafc_storage_patch(p)/tot_leaf
                   pxfer = c13_cnveg_carbonstate_inst%leafc_xfer_patch(p)/tot_leaf
                else
                   ! when initiating from zero leaf state, use evergreen flag to set proportions
                   if (pftcon%evergreen(patch%itype(p)) == 1._r8) then
                      pleaf = 1._r8
                   else
                      pstor = 1._r8
                   end if
                end if
                c13_cnveg_carbonstate_inst%leafc_patch(p)              = c13_cnveg_carbonstate_inst%leafc_patch(p)              *t1 &
                     + leafc13_seed*pleaf*t2
                c13_cnveg_carbonstate_inst%leafc_storage_patch(p)      = c13_cnveg_carbonstate_inst%leafc_storage_patch(p)      *t1 &
                     + leafc13_seed*pstor*t2
                c13_cnveg_carbonstate_inst%leafc_xfer_patch(p)         = c13_cnveg_carbonstate_inst%leafc_xfer_patch(p)         *t1 &
                     + leafc13_seed*pxfer*t2
                c13_cnveg_carbonstate_inst%frootc_patch(p)             = c13_cnveg_carbonstate_inst%frootc_patch(p)             *t1
                c13_cnveg_carbonstate_inst%frootc_storage_patch(p)     = c13_cnveg_carbonstate_inst%frootc_storage_patch(p)     *t1
                c13_cnveg_carbonstate_inst%frootc_xfer_patch(p)        = c13_cnveg_carbonstate_inst%frootc_xfer_patch(p)        *t1
                c13_cnveg_carbonstate_inst%livestemc_patch(p)          = c13_cnveg_carbonstate_inst%livestemc_patch(p)          *t1
                c13_cnveg_carbonstate_inst%livestemc_storage_patch(p)  = c13_cnveg_carbonstate_inst%livestemc_storage_patch(p)  *t1
                c13_cnveg_carbonstate_inst%livestemc_xfer_patch(p)     = c13_cnveg_carbonstate_inst%livestemc_xfer_patch(p)     *t1
                c13_cnveg_carbonstate_inst%deadstemc_patch(p)          = c13_cnveg_carbonstate_inst%deadstemc_patch(p)          *t1 &
                     + deadstemc13_seed*t2
                c13_cnveg_carbonstate_inst%deadstemc_storage_patch(p)  = c13_cnveg_carbonstate_inst%deadstemc_storage_patch(p)  *t1
                c13_cnveg_carbonstate_inst%deadstemc_xfer_patch(p)     = c13_cnveg_carbonstate_inst%deadstemc_xfer_patch(p)     *t1
                c13_cnveg_carbonstate_inst%livecrootc_patch(p)         = c13_cnveg_carbonstate_inst%livecrootc_patch(p)         *t1
                c13_cnveg_carbonstate_inst%livecrootc_storage_patch(p) = c13_cnveg_carbonstate_inst%livecrootc_storage_patch(p) *t1
                c13_cnveg_carbonstate_inst%livecrootc_xfer_patch(p)    = c13_cnveg_carbonstate_inst%livecrootc_xfer_patch(p)    *t1
                c13_cnveg_carbonstate_inst%deadcrootc_patch(p)         = c13_cnveg_carbonstate_inst%deadcrootc_patch(p)         *t1
                c13_cnveg_carbonstate_inst%deadcrootc_storage_patch(p) = c13_cnveg_carbonstate_inst%deadcrootc_storage_patch(p) *t1
                c13_cnveg_carbonstate_inst%deadcrootc_xfer_patch(p)    = c13_cnveg_carbonstate_inst%deadcrootc_xfer_patch(p)    *t1
                c13_cnveg_carbonstate_inst%gresp_storage_patch(p)      = c13_cnveg_carbonstate_inst%gresp_storage_patch(p)      *t1
                c13_cnveg_carbonstate_inst%gresp_xfer_patch(p)         = c13_cnveg_carbonstate_inst%gresp_xfer_patch(p)         *t1
                c13_cnveg_carbonstate_inst%cpool_patch(p)              = c13_cnveg_carbonstate_inst%cpool_patch(p)              *t1
                c13_cnveg_carbonstate_inst%xsmrpool_patch(p)           = c13_cnveg_carbonstate_inst%xsmrpool_patch(p)           *t1
                c13_cnveg_carbonstate_inst%ctrunc_patch(p)             = c13_cnveg_carbonstate_inst%ctrunc_patch(p)             *t1
                c13_cnveg_carbonstate_inst%dispvegc_patch(p)           = c13_cnveg_carbonstate_inst%dispvegc_patch(p)           *t1
                c13_cnveg_carbonstate_inst%storvegc_patch(p)           = c13_cnveg_carbonstate_inst%storvegc_patch(p)           *t1
                c13_cnveg_carbonstate_inst%totc_patch(p)               = c13_cnveg_carbonstate_inst%totc_patch(p)               *t1
                c13_cnveg_carbonstate_inst%totvegc_patch(p)            = c13_cnveg_carbonstate_inst%totvegc_patch(p)            *t1
                
             endif
             
             if ( use_c14 ) then
                ! patch-level carbon-14 state variables 
                tot_leaf = &
                     c14_cnveg_carbonstate_inst%leafc_patch(p) + &
                     c14_cnveg_carbonstate_inst%leafc_storage_patch(p) + &
                     c14_cnveg_carbonstate_inst%leafc_xfer_patch(p)
                pleaf = 0._r8
                pstor = 0._r8
                pxfer = 0._r8
                if (tot_leaf /= 0._r8) then
                   pleaf = c14_cnveg_carbonstate_inst%leafc_patch(p)/tot_leaf
                   pstor = c14_cnveg_carbonstate_inst%leafc_storage_patch(p)/tot_leaf
                   pxfer = c14_cnveg_carbonstate_inst%leafc_xfer_patch(p)/tot_leaf
                else
                   ! when initiating from zero leaf state, use evergreen flag to set proportions
                   if (pftcon%evergreen(patch%itype(p)) == 1._r8) then
                      pleaf = 1._r8
                   else
                      pstor = 1._r8
                   end if
                end if
                c14_cnveg_carbonstate_inst%leafc_patch(p)              = c14_cnveg_carbonstate_inst%leafc_patch(p)              *t1 &
                     + leafc14_seed*pleaf*t2
                c14_cnveg_carbonstate_inst%leafc_storage_patch(p)      = c14_cnveg_carbonstate_inst%leafc_storage_patch(p)      *t1 &
                     + leafc14_seed*pstor*t2
                c14_cnveg_carbonstate_inst%leafc_xfer_patch(p)         = c14_cnveg_carbonstate_inst%leafc_xfer_patch(p)         *t1 &
                     + leafc14_seed*pxfer*t2
                c14_cnveg_carbonstate_inst%frootc_patch(p)             = c14_cnveg_carbonstate_inst%frootc_patch(p)             *t1
                c14_cnveg_carbonstate_inst%frootc_storage_patch(p)     = c14_cnveg_carbonstate_inst%frootc_storage_patch(p)     *t1
                c14_cnveg_carbonstate_inst%frootc_xfer_patch(p)        = c14_cnveg_carbonstate_inst%frootc_xfer_patch(p)        *t1
                c14_cnveg_carbonstate_inst%livestemc_patch(p)          = c14_cnveg_carbonstate_inst%livestemc_patch(p)          *t1
                c14_cnveg_carbonstate_inst%livestemc_storage_patch(p)  = c14_cnveg_carbonstate_inst%livestemc_storage_patch(p)  *t1
                c14_cnveg_carbonstate_inst%livestemc_xfer_patch(p)     = c14_cnveg_carbonstate_inst%livestemc_xfer_patch(p)     *t1
                c14_cnveg_carbonstate_inst%deadstemc_patch(p)          = c14_cnveg_carbonstate_inst%deadstemc_patch(p)          *t1 &
                     + deadstemc14_seed*t2
                c14_cnveg_carbonstate_inst%deadstemc_storage_patch(p)  = c14_cnveg_carbonstate_inst%deadstemc_storage_patch(p)  *t1
                c14_cnveg_carbonstate_inst%deadstemc_xfer_patch(p)     = c14_cnveg_carbonstate_inst%deadstemc_xfer_patch(p)     *t1
                c14_cnveg_carbonstate_inst%livecrootc_patch(p)         = c14_cnveg_carbonstate_inst%livecrootc_patch(p)         *t1
                c14_cnveg_carbonstate_inst%livecrootc_storage_patch(p) = c14_cnveg_carbonstate_inst%livecrootc_storage_patch(p) *t1
                c14_cnveg_carbonstate_inst%livecrootc_xfer_patch(p)    = c14_cnveg_carbonstate_inst%livecrootc_xfer_patch(p)    *t1
                c14_cnveg_carbonstate_inst%deadcrootc_patch(p)         = c14_cnveg_carbonstate_inst%deadcrootc_patch(p)         *t1
                c14_cnveg_carbonstate_inst%deadcrootc_storage_patch(p) = c14_cnveg_carbonstate_inst%deadcrootc_storage_patch(p) *t1
                c14_cnveg_carbonstate_inst%deadcrootc_xfer_patch(p)    = c14_cnveg_carbonstate_inst%deadcrootc_xfer_patch(p)    *t1
                c14_cnveg_carbonstate_inst%gresp_storage_patch(p)      = c14_cnveg_carbonstate_inst%gresp_storage_patch(p)      *t1
                c14_cnveg_carbonstate_inst%gresp_xfer_patch(p)         = c14_cnveg_carbonstate_inst%gresp_xfer_patch(p)         *t1
                c14_cnveg_carbonstate_inst%cpool_patch(p)              = c14_cnveg_carbonstate_inst%cpool_patch(p)              *t1
                c14_cnveg_carbonstate_inst%xsmrpool_patch(p)           = c14_cnveg_carbonstate_inst%xsmrpool_patch(p)           *t1
                c14_cnveg_carbonstate_inst%ctrunc_patch(p)             = c14_cnveg_carbonstate_inst%ctrunc_patch(p)             *t1
                c14_cnveg_carbonstate_inst%dispvegc_patch(p)           = c14_cnveg_carbonstate_inst%dispvegc_patch(p)           *t1
                c14_cnveg_carbonstate_inst%storvegc_patch(p)           = c14_cnveg_carbonstate_inst%storvegc_patch(p)           *t1
                c14_cnveg_carbonstate_inst%totc_patch(p)               = c14_cnveg_carbonstate_inst%totc_patch(p)               *t1
                c14_cnveg_carbonstate_inst%totvegc_patch(p)            = c14_cnveg_carbonstate_inst%totvegc_patch(p)            *t1
             endif
             
             tot_leaf = cnveg_nitrogenstate_inst%leafn_patch(p) + &
                  cnveg_nitrogenstate_inst%leafn_storage_patch(p) + &
                  cnveg_nitrogenstate_inst%leafn_xfer_patch(p)
             pleaf = 0._r8
             pstor = 0._r8
             pxfer = 0._r8
             if (tot_leaf /= 0._r8) then
                pleaf = cnveg_nitrogenstate_inst%leafn_patch(p)/tot_leaf
                pstor = cnveg_nitrogenstate_inst%leafn_storage_patch(p)/tot_leaf
                pxfer = cnveg_nitrogenstate_inst%leafn_xfer_patch(p)/tot_leaf
             else
                ! when initiating from zero leaf state, use evergreen flag to set proportions
                if (pftcon%evergreen(patch%itype(p)) == 1._r8) then
                   pleaf = 1._r8
                else
                   pstor = 1._r8
                end if
             end if
             ! patch-level nitrogen state variables
             cnveg_nitrogenstate_inst%leafn_patch(p)              = cnveg_nitrogenstate_inst%leafn_patch(p)              *t1 &
                  + leafn_seed*pleaf*t2
             cnveg_nitrogenstate_inst%leafn_storage_patch(p)      = cnveg_nitrogenstate_inst%leafn_storage_patch(p)      *t1 &
                  + leafn_seed*pstor*t2
             cnveg_nitrogenstate_inst%leafn_xfer_patch(p)         = cnveg_nitrogenstate_inst%leafn_xfer_patch(p)         *t1 &
                  + leafn_seed*pxfer*t2
             cnveg_nitrogenstate_inst%frootn_patch(p)  	          = cnveg_nitrogenstate_inst%frootn_patch(p)             *t1
             cnveg_nitrogenstate_inst%frootn_storage_patch(p)     = cnveg_nitrogenstate_inst%frootn_storage_patch(p)     *t1
             cnveg_nitrogenstate_inst%frootn_xfer_patch(p)        = cnveg_nitrogenstate_inst%frootn_xfer_patch(p)        *t1
             cnveg_nitrogenstate_inst%livestemn_patch(p)	  = cnveg_nitrogenstate_inst%livestemn_patch(p)          *t1
             cnveg_nitrogenstate_inst%livestemn_storage_patch(p)  = cnveg_nitrogenstate_inst%livestemn_storage_patch(p)  *t1
             cnveg_nitrogenstate_inst%livestemn_xfer_patch(p)     = cnveg_nitrogenstate_inst%livestemn_xfer_patch(p)     *t1
             cnveg_nitrogenstate_inst%deadstemn_patch(p)          = cnveg_nitrogenstate_inst%deadstemn_patch(p)          *t1 &
                  + deadstemn_seed*t2
             cnveg_nitrogenstate_inst%deadstemn_storage_patch(p)  = cnveg_nitrogenstate_inst%deadstemn_storage_patch(p)  *t1
             cnveg_nitrogenstate_inst%deadstemn_xfer_patch(p)     = cnveg_nitrogenstate_inst%deadstemn_xfer_patch(p)     *t1
             cnveg_nitrogenstate_inst%livecrootn_patch(p)         = cnveg_nitrogenstate_inst%livecrootn_patch(p)         *t1
             cnveg_nitrogenstate_inst%livecrootn_storage_patch(p) = cnveg_nitrogenstate_inst%livecrootn_storage_patch(p) *t1
             cnveg_nitrogenstate_inst%livecrootn_xfer_patch(p)    = cnveg_nitrogenstate_inst%livecrootn_xfer_patch(p)    *t1
             cnveg_nitrogenstate_inst%deadcrootn_patch(p)         = cnveg_nitrogenstate_inst%deadcrootn_patch(p)         *t1
             cnveg_nitrogenstate_inst%deadcrootn_storage_patch(p) = cnveg_nitrogenstate_inst%deadcrootn_storage_patch(p) *t1
             cnveg_nitrogenstate_inst%deadcrootn_xfer_patch(p)    = cnveg_nitrogenstate_inst%deadcrootn_xfer_patch(p)    *t1
             cnveg_nitrogenstate_inst%retransn_patch(p)	          = cnveg_nitrogenstate_inst%retransn_patch(p)           *t1
             cnveg_nitrogenstate_inst%npool_patch(p)              = cnveg_nitrogenstate_inst%npool_patch(p)              *t1
             cnveg_nitrogenstate_inst%ntrunc_patch(p)             = cnveg_nitrogenstate_inst%ntrunc_patch(p)             *t1
             cnveg_nitrogenstate_inst%dispvegn_patch(p)	          = cnveg_nitrogenstate_inst%dispvegn_patch(p)           *t1
             cnveg_nitrogenstate_inst%storvegn_patch(p)	          = cnveg_nitrogenstate_inst%storvegn_patch(p)           *t1
             cnveg_nitrogenstate_inst%totvegn_patch(p) 	          = cnveg_nitrogenstate_inst%totvegn_patch(p)            *t1
             cnveg_nitrogenstate_inst%totn_patch(p) 	          = cnveg_nitrogenstate_inst%totn_patch(p)               *t1
             
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
             
          else if (dwt < 0._r8) then
             
             ! if the pft lost weight on the timestep, then the carbon and nitrogen state
             ! variables are directed to litter, CWD, and wood product pools.
             
             ! N.B. : the conv_cflux, prod10_cflux, and prod100_cflux fluxes are accumulated
             ! as negative values, but the fluxes for pft-to-litter are accumulated as 
             ! positive values
             
             ! set local weight variables for this pft
             wt_new = patch%wtcol(p)
             wt_old = prior_weights%pwtcol(p)
             
             !---------------
             ! C state update
             !---------------
             
             ! leafc 
             ptr => cnveg_carbonstate_inst%leafc_patch(p)
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
             ptr => cnveg_carbonstate_inst%leafc_storage_patch(p)
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
             ptr => cnveg_carbonstate_inst%leafc_xfer_patch(p)
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
             ptr => cnveg_carbonstate_inst%frootc_patch(p)
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
             ptr => cnveg_carbonstate_inst%frootc_storage_patch(p)
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
             ptr => cnveg_carbonstate_inst%frootc_xfer_patch(p)
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
             ptr => cnveg_carbonstate_inst%livestemc_patch(p)
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
             ptr => cnveg_carbonstate_inst%livestemc_storage_patch(p)
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
             ptr => cnveg_carbonstate_inst%livestemc_xfer_patch(p)
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
             ptr => cnveg_carbonstate_inst%deadstemc_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                conv_cflux(p) = conv_cflux(p) + change_state * pftcon%pconv(patch%itype(p))
                prod10_cflux(p) = prod10_cflux(p) + change_state * pftcon%pprod10(patch%itype(p))
                prod100_cflux(p) = prod100_cflux(p) + change_state * pftcon%pprod100(patch%itype(p))
             else
                ptr = 0._r8
                conv_cflux(p) = conv_cflux(p) - init_state * pftcon%pconv(patch%itype(p))
                prod10_cflux(p) = prod10_cflux(p) - init_state * pftcon%pprod10(patch%itype(p))
                prod100_cflux(p) = prod100_cflux(p) - init_state * pftcon%pprod100(patch%itype(p))
             end if
             
             ! deadstemc_storage 
             ptr => cnveg_carbonstate_inst%deadstemc_storage_patch(p)
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
             ptr => cnveg_carbonstate_inst%deadstemc_xfer_patch(p)
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
             ptr => cnveg_carbonstate_inst%livecrootc_patch(p)
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
             ptr => cnveg_carbonstate_inst%livecrootc_storage_patch(p)
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
             ptr => cnveg_carbonstate_inst%livecrootc_xfer_patch(p)
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
             ptr => cnveg_carbonstate_inst%deadcrootc_patch(p)
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
             ptr => cnveg_carbonstate_inst%deadcrootc_storage_patch(p)
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
             ptr => cnveg_carbonstate_inst%deadcrootc_xfer_patch(p)
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
             ptr => cnveg_carbonstate_inst%gresp_storage_patch(p)
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
             ptr => cnveg_carbonstate_inst%gresp_xfer_patch(p)
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
             ptr => cnveg_carbonstate_inst%cpool_patch(p)
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
             ptr => cnveg_carbonstate_inst%xsmrpool_patch(p)
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
             ptr => cnveg_carbonstate_inst%ctrunc_patch(p)
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
                ptr => cnveg_carbonstate_inst%leafc_patch(p)
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
                ptr => cnveg_carbonstate_inst%leafc_storage_patch(p)
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
                ptr => cnveg_carbonstate_inst%leafc_xfer_patch(p)
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
                ptr => cnveg_carbonstate_inst%frootc_patch(p)
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
                ptr => cnveg_carbonstate_inst%frootc_storage_patch(p)
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
                ptr => cnveg_carbonstate_inst%frootc_xfer_patch(p)
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
                ptr => cnveg_carbonstate_inst%livestemc_patch(p)
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
                ptr => cnveg_carbonstate_inst%livestemc_storage_patch(p)
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
                ptr => cnveg_carbonstate_inst%livestemc_xfer_patch(p)
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
                ptr => cnveg_carbonstate_inst%deadstemc_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state * pftcon%pconv(patch%itype(p))
                   dwt_ptr2 = dwt_ptr2 + change_state * pftcon%pprod10(patch%itype(p))
                   dwt_ptr3 = dwt_ptr3 + change_state * pftcon%pprod100(patch%itype(p))
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state * pftcon%pconv(patch%itype(p))
                   dwt_ptr2 = dwt_ptr2 - init_state * pftcon%pprod10(patch%itype(p))
                   dwt_ptr3 = dwt_ptr3 - init_state * pftcon%pprod100(patch%itype(p))
                end if
                
                ! deadstemc_storage 
                ptr => cnveg_carbonstate_inst%deadstemc_storage_patch(p)
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
                ptr => cnveg_carbonstate_inst%deadstemc_xfer_patch(p)
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
                ptr => cnveg_carbonstate_inst%livecrootc_patch(p)
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
                ptr => cnveg_carbonstate_inst%livecrootc_storage_patch(p)
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
                ptr => cnveg_carbonstate_inst%livecrootc_xfer_patch(p)
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
                ptr => cnveg_carbonstate_inst%deadcrootc_patch(p)
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
                ptr => cnveg_carbonstate_inst%deadcrootc_storage_patch(p)
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
                ptr => cnveg_carbonstate_inst%deadcrootc_xfer_patch(p)
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
                ptr => cnveg_carbonstate_inst%gresp_storage_patch(p)
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
                ptr => cnveg_carbonstate_inst%gresp_xfer_patch(p)
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
                ptr => cnveg_carbonstate_inst%cpool_patch(p)
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
                ptr => cnveg_carbonstate_inst%ctrunc_patch(p)
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
                
                ! set pointers to the conversion and product pool fluxes for this patch
                ! dwt_ptr0 is reserved for local assignment to dwt_xxx_to_litter fluxes
                dwt_ptr1 => conv_c14flux(p)
                dwt_ptr2 => prod10_c14flux(p)
                dwt_ptr3 => prod100_c14flux(p)
                
                ! leafc 
                ptr => cnveg_carbonstate_inst%leafc_patch(p)
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
                ptr => cnveg_carbonstate_inst%leafc_storage_patch(p)
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
                ptr => cnveg_carbonstate_inst%leafc_xfer_patch(p)
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
                ptr => cnveg_carbonstate_inst%frootc_patch(p)
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
                ptr => cnveg_carbonstate_inst%frootc_storage_patch(p)
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
                ptr => cnveg_carbonstate_inst%frootc_xfer_patch(p)
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
                ptr => cnveg_carbonstate_inst%livestemc_patch(p)
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
                ptr => cnveg_carbonstate_inst%livestemc_storage_patch(p)
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
                ptr => cnveg_carbonstate_inst%livestemc_xfer_patch(p)
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
                ptr => cnveg_carbonstate_inst%deadstemc_patch(p)
                init_state = ptr*wt_old
                change_state = ptr*dwt
                new_state = init_state+change_state
                if (wt_new /= 0._r8) then
                   ptr = new_state/wt_new
                   dwt_ptr1 = dwt_ptr1 + change_state * pftcon%pconv(patch%itype(p))
                   dwt_ptr2 = dwt_ptr2 + change_state * pftcon%pprod10(patch%itype(p))
                   dwt_ptr3 = dwt_ptr3 + change_state * pftcon%pprod100(patch%itype(p))
                else
                   ptr = 0._r8
                   dwt_ptr1 = dwt_ptr1 - init_state * pftcon%pconv(patch%itype(p))
                   dwt_ptr2 = dwt_ptr2 - init_state * pftcon%pprod10(patch%itype(p))
                   dwt_ptr3 = dwt_ptr3 - init_state * pftcon%pprod100(patch%itype(p))
                end if
                
                ! deadstemc_storage 
                ptr => cnveg_carbonstate_inst%deadstemc_storage_patch(p)
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
                ptr => cnveg_carbonstate_inst%deadstemc_xfer_patch(p)
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
                ptr => cnveg_carbonstate_inst%livecrootc_patch(p)
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
                ptr => cnveg_carbonstate_inst%livecrootc_storage_patch(p)
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
                ptr => cnveg_carbonstate_inst%livecrootc_xfer_patch(p)
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
                ptr => cnveg_carbonstate_inst%deadcrootc_patch(p)
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
                ptr => cnveg_carbonstate_inst%deadcrootc_storage_patch(p)
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
                ptr => cnveg_carbonstate_inst%deadcrootc_xfer_patch(p)
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
                ptr => cnveg_carbonstate_inst%gresp_storage_patch(p)
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
                ptr => cnveg_carbonstate_inst%gresp_xfer_patch(p)
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
                ptr => cnveg_carbonstate_inst%cpool_patch(p)
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
                ptr => cnveg_carbonstate_inst%ctrunc_patch(p)
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
             
             ! set pointers to the conversion and product pool fluxes for this patch
             ! dwt_ptr0 is reserved for local assignment to dwt_xxx_to_litter fluxes
             dwt_ptr1 => conv_nflux(p)
             dwt_ptr2 => prod10_nflux(p)
             dwt_ptr3 => prod100_nflux(p)
             
             ! leafn 
             ptr => cnveg_nitrogenstate_inst%leafn_patch(p)
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
             ptr => cnveg_nitrogenstate_inst%leafn_storage_patch(p)
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
             ptr => cnveg_nitrogenstate_inst%leafn_xfer_patch(p)
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
             ptr => cnveg_nitrogenstate_inst%frootn_patch(p)
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
             ptr => cnveg_nitrogenstate_inst%frootn_storage_patch(p)
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
             ptr => cnveg_nitrogenstate_inst%frootn_xfer_patch(p)
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
             ptr => cnveg_nitrogenstate_inst%livestemn_patch(p)
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
             ptr => cnveg_nitrogenstate_inst%livestemn_storage_patch(p)
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
             ptr => cnveg_nitrogenstate_inst%livestemn_xfer_patch(p)
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
             ptr => cnveg_nitrogenstate_inst%deadstemn_patch(p)
             init_state = ptr*wt_old
             change_state = ptr*dwt
             new_state = init_state+change_state
             if (wt_new /= 0._r8) then
                ptr = new_state/wt_new
                dwt_ptr1 = dwt_ptr1 + change_state * pftcon%pconv(patch%itype(p))
                dwt_ptr2 = dwt_ptr2 + change_state * pftcon%pprod10(patch%itype(p))
                dwt_ptr3 = dwt_ptr3 + change_state * pftcon%pprod100(patch%itype(p))
             else
                ptr = 0._r8
                dwt_ptr1 = dwt_ptr1 - init_state * pftcon%pconv(patch%itype(p))
                dwt_ptr2 = dwt_ptr2 - init_state * pftcon%pprod10(patch%itype(p))
                dwt_ptr3 = dwt_ptr3 - init_state * pftcon%pprod100(patch%itype(p))
             end if
             
             ! deadstemn_storage 
             ptr => cnveg_nitrogenstate_inst%deadstemn_storage_patch(p)
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
             ptr => cnveg_nitrogenstate_inst%deadstemn_xfer_patch(p)
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
             ptr => cnveg_nitrogenstate_inst%livecrootn_patch(p)
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
             ptr => cnveg_nitrogenstate_inst%livecrootn_storage_patch(p)
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
             ptr => cnveg_nitrogenstate_inst%livecrootn_xfer_patch(p)
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
             ptr => cnveg_nitrogenstate_inst%deadcrootn_patch(p)
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
             ptr => cnveg_nitrogenstate_inst%deadcrootn_storage_patch(p)
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
             ptr => cnveg_nitrogenstate_inst%deadcrootn_xfer_patch(p)
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
             ptr => cnveg_nitrogenstate_inst%retransn_patch(p)
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
             ptr => cnveg_nitrogenstate_inst%npool_patch(p)
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
             ptr => cnveg_nitrogenstate_inst%ntrunc_patch(p)
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
       end if          ! is soil
    end do             ! patch loop
    
    ! calculate column-level seeding fluxes
    do pi = 1,max_patch_per_col
       do c = bounds%begc, bounds%endc
          if ( pi <=  col%npatches(c) ) then
             p = col%patchi(c) + pi - 1
             
             ! C fluxes
             cnveg_carbonflux_inst%dwt_seedc_to_leaf_col(c) = &
                  cnveg_carbonflux_inst%dwt_seedc_to_leaf_col(c) + dwt_leafc_seed(p)/dt
             cnveg_carbonflux_inst%dwt_seedc_to_deadstem_col(c) = &
                  cnveg_carbonflux_inst%dwt_seedc_to_deadstem_col(c) + dwt_deadstemc_seed(p)/dt
             
             if ( use_c13 ) then
                c13_cnveg_carbonflux_inst%dwt_seedc_to_leaf_col(c) = &
                     c13_cnveg_carbonflux_inst%dwt_seedc_to_leaf_col(c) + dwt_leafc13_seed(p)/dt
                c13_cnveg_carbonflux_inst%dwt_seedc_to_deadstem_col(c) = &
                     c13_cnveg_carbonflux_inst%dwt_seedc_to_deadstem_col(c) + dwt_deadstemc13_seed(p)/dt
             endif
             
             if ( use_c14 ) then	
                c14_cnveg_carbonflux_inst%dwt_seedc_to_leaf_col(c) = &
                     c14_cnveg_carbonflux_inst%dwt_seedc_to_leaf_col(c) + dwt_leafc14_seed(p)/dt
                c14_cnveg_carbonflux_inst%dwt_seedc_to_deadstem_col(c) = &
                     c14_cnveg_carbonflux_inst%dwt_seedc_to_deadstem_col(c) + dwt_deadstemc14_seed(p)/dt
             endif
             
             ! N fluxes
             cnveg_nitrogenflux_inst%dwt_seedn_to_leaf_col(c) = &
                  cnveg_nitrogenflux_inst%dwt_seedn_to_leaf_col(c) + dwt_leafn_seed(p)/dt
             cnveg_nitrogenflux_inst%dwt_seedn_to_deadstem_col(c) = &
                  cnveg_nitrogenflux_inst%dwt_seedn_to_deadstem_col(c) + dwt_deadstemn_seed(p)/dt
          end if
       end do
    end do
    
    
    ! calculate patch-to-column for fluxes into litter and CWD pools
    do j = 1, nlevdecomp
       do pi = 1,max_patch_per_col
          do c = bounds%begc, bounds%endc
             if ( pi <=  col%npatches(c) ) then
                p = col%patchi(c) + pi - 1
                
                ! fine root litter carbon fluxes
                cnveg_carbonflux_inst%dwt_frootc_to_litr_met_c_col(c,j) = &
                     cnveg_carbonflux_inst%dwt_frootc_to_litr_met_c_col(c,j) + &
                     (dwt_frootc_to_litter(p)*pftcon%fr_flab(patch%itype(p)))/dt &
                     * soilbiogeochem_state_inst%froot_prof_patch(p,j)

                cnveg_carbonflux_inst%dwt_frootc_to_litr_cel_c_col(c,j) = &
                     cnveg_carbonflux_inst%dwt_frootc_to_litr_cel_c_col(c,j) + &
                     (dwt_frootc_to_litter(p)*pftcon%fr_fcel(patch%itype(p)))/dt &
                     * soilbiogeochem_state_inst%froot_prof_patch(p,j)

                cnveg_carbonflux_inst%dwt_frootc_to_litr_lig_c_col(c,j) = &
                     cnveg_carbonflux_inst%dwt_frootc_to_litr_lig_c_col(c,j) + &
                     (dwt_frootc_to_litter(p)*pftcon%fr_flig(patch%itype(p)))/dt &
                     * soilbiogeochem_state_inst%froot_prof_patch(p,j)
                
                
                ! fine root litter nitrogen fluxes
                cnveg_nitrogenflux_inst%dwt_frootn_to_litr_met_n_col(c,j) = &
                     cnveg_nitrogenflux_inst%dwt_frootn_to_litr_met_n_col(c,j) + &
                     (dwt_frootn_to_litter(p)*pftcon%fr_flab(patch%itype(p)))/dt &
                     * soilbiogeochem_state_inst%froot_prof_patch(p,j)
                cnveg_nitrogenflux_inst%dwt_frootn_to_litr_cel_n_col(c,j) = &
                     cnveg_nitrogenflux_inst%dwt_frootn_to_litr_cel_n_col(c,j) + &
                     (dwt_frootn_to_litter(p)*pftcon%fr_fcel(patch%itype(p)))/dt &
                     * soilbiogeochem_state_inst%froot_prof_patch(p,j)

                cnveg_nitrogenflux_inst%dwt_frootn_to_litr_lig_n_col(c,j) = &
                     cnveg_nitrogenflux_inst%dwt_frootn_to_litr_lig_n_col(c,j) + &
                     (dwt_frootn_to_litter(p)*pftcon%fr_flig(patch%itype(p)))/dt &
                     * soilbiogeochem_state_inst%froot_prof_patch(p,j)
                
                ! livecroot fluxes to cwd
                cnveg_carbonflux_inst%dwt_livecrootc_to_cwdc_col(c,j) = &
                     cnveg_carbonflux_inst%dwt_livecrootc_to_cwdc_col(c,j) + &
                     (dwt_livecrootc_to_litter(p))/dt * soilbiogeochem_state_inst%croot_prof_patch(p,j)

                cnveg_nitrogenflux_inst%dwt_livecrootn_to_cwdn_col(c,j) = &
                     cnveg_nitrogenflux_inst%dwt_livecrootn_to_cwdn_col(c,j) + &
                     (dwt_livecrootn_to_litter(p))/dt * soilbiogeochem_state_inst%croot_prof_patch(p,j)
                
                ! deadcroot fluxes to cwd
                cnveg_carbonflux_inst%dwt_deadcrootc_to_cwdc_col(c,j) = &
                     cnveg_carbonflux_inst%dwt_deadcrootc_to_cwdc_col(c,j) + &
                     (dwt_deadcrootc_to_litter(p))/dt * soilbiogeochem_state_inst%croot_prof_patch(p,j)

                cnveg_nitrogenflux_inst%dwt_deadcrootn_to_cwdn_col(c,j) = &
                     cnveg_nitrogenflux_inst%dwt_deadcrootn_to_cwdn_col(c,j) + &
                     (dwt_deadcrootn_to_litter(p))/dt * soilbiogeochem_state_inst%croot_prof_patch(p,j)
             
                if ( use_c13 ) then
                   ! C13 fine root litter fluxes
                   c13_cnveg_carbonflux_inst%dwt_frootc_to_litr_met_c_col(c,j) = &
                        c13_cnveg_carbonflux_inst%dwt_frootc_to_litr_met_c_col(c,j) + &
                        (dwt_frootc13_to_litter(p)*pftcon%fr_flab(patch%itype(p)))/dt &
                        * soilbiogeochem_state_inst%froot_prof_patch(p,j)

                   c13_cnveg_carbonflux_inst%dwt_frootc_to_litr_cel_c_col(c,j) = &
                        c13_cnveg_carbonflux_inst%dwt_frootc_to_litr_cel_c_col(c,j) + &
                        (dwt_frootc13_to_litter(p)*pftcon%fr_fcel(patch%itype(p)))/dt &
                        * soilbiogeochem_state_inst%froot_prof_patch(p,j)
                   
                   c13_cnveg_carbonflux_inst%dwt_frootc_to_litr_lig_c_col(c,j) = &
                        c13_cnveg_carbonflux_inst%dwt_frootc_to_litr_lig_c_col(c,j) + &
                        (dwt_frootc13_to_litter(p)*pftcon%fr_flig(patch%itype(p)))/dt &
                        * soilbiogeochem_state_inst%froot_prof_patch(p,j)

                   ! livecroot fluxes to cwd
                   c13_cnveg_carbonflux_inst%dwt_livecrootc_to_cwdc_col(c,j) = &
                        c13_cnveg_carbonflux_inst%dwt_livecrootc_to_cwdc_col(c,j) + &
                        (dwt_livecrootc13_to_litter(p))/dt * soilbiogeochem_state_inst%croot_prof_patch(p,j)

                   ! deadcroot fluxes to cwd
                   c13_cnveg_carbonflux_inst%dwt_deadcrootc_to_cwdc_col(c,j) = &
                        c13_cnveg_carbonflux_inst%dwt_deadcrootc_to_cwdc_col(c,j) + &
                        (dwt_deadcrootc13_to_litter(p))/dt * soilbiogeochem_state_inst%croot_prof_patch(p,j)
                   
                endif
                
                if ( use_c14 ) then                   
                   ! C14 fine root litter fluxes
                   c14_cnveg_carbonflux_inst%dwt_frootc_to_litr_met_c_col(c,j) = &
                        c14_cnveg_carbonflux_inst%dwt_frootc_to_litr_met_c_col(c,j) + &
                        (dwt_frootc14_to_litter(p)*pftcon%fr_flab(patch%itype(p)))/dt &
                        * soilbiogeochem_state_inst%froot_prof_patch(p,j)

                   c14_cnveg_carbonflux_inst%dwt_frootc_to_litr_cel_c_col(c,j) = &
                        c14_cnveg_carbonflux_inst%dwt_frootc_to_litr_cel_c_col(c,j) + &
                        (dwt_frootc14_to_litter(p)*pftcon%fr_fcel(patch%itype(p)))/dt &
                        * soilbiogeochem_state_inst%froot_prof_patch(p,j)

                   c14_cnveg_carbonflux_inst%dwt_frootc_to_litr_lig_c_col(c,j) = &
                        c14_cnveg_carbonflux_inst%dwt_frootc_to_litr_lig_c_col(c,j) + &
                        (dwt_frootc14_to_litter(p)*pftcon%fr_flig(patch%itype(p)))/dt &
                        * soilbiogeochem_state_inst%froot_prof_patch(p,j)

                   ! livecroot fluxes to cwd
                   c14_cnveg_carbonflux_inst%dwt_livecrootc_to_cwdc_col(c,j) = &
                        c14_cnveg_carbonflux_inst%dwt_livecrootc_to_cwdc_col(c,j) + &
                        (dwt_livecrootc14_to_litter(p))/dt * soilbiogeochem_state_inst%croot_prof_patch(p,j)

                   ! deadcroot fluxes to cwd
                   c14_cnveg_carbonflux_inst%dwt_deadcrootc_to_cwdc_col(c,j) = &
                        c14_cnveg_carbonflux_inst%dwt_deadcrootc_to_cwdc_col(c,j) + &
                        (dwt_deadcrootc14_to_litter(p))/dt * soilbiogeochem_state_inst%croot_prof_patch(p,j)
                endif
                
             end if
          end do
       end do
    end do

    ! calculate patch-to-column for fluxes into product pools and conversion flux
    do pi = 1,max_patch_per_col
       do c = bounds%begc,bounds%endc
          if (pi <= col%npatches(c)) then
             p = col%patchi(c) + pi - 1
             
             ! column-level fluxes are accumulated as positive fluxes.
             ! column-level C flux updates
             cnveg_carbonflux_inst%dwt_conv_cflux_col(c) = &
                  cnveg_carbonflux_inst%dwt_conv_cflux_col(c) - conv_cflux(p)/dt
             cnveg_carbonflux_inst%dwt_prod10c_gain_col(c) = &
                  cnveg_carbonflux_inst%dwt_prod10c_gain_col(c) - prod10_cflux(p)/dt
             cnveg_carbonflux_inst%dwt_prod100c_gain_col(c) = &
                  cnveg_carbonflux_inst%dwt_prod100c_gain_col(c) - prod100_cflux(p)/dt

             ! These magic constants should be replaced with: nbrdlf_evr_trp_tree and nbrdlf_dcd_trp_tree
             if(patch%itype(p)==4.or.patch%itype(p)==6)then
                soilbiogeochem_carbonflux_inst%lf_conv_cflux_col(c) = &
                     soilbiogeochem_carbonflux_inst%lf_conv_cflux_col(c) - conv_cflux(p)/dt
             end if
             
             if ( use_c13 ) then
                ! C13 column-level flux updates
                c13_cnveg_carbonflux_inst%dwt_conv_cflux_col(c) = &
                     c13_cnveg_carbonflux_inst%dwt_conv_cflux_col(c) - conv_c13flux(p)/dt
                c13_cnveg_carbonflux_inst%dwt_prod10c_gain_col(c) = &
                     c13_cnveg_carbonflux_inst%dwt_prod10c_gain_col(c) - prod10_c13flux(p)/dt
                c13_cnveg_carbonflux_inst%dwt_prod100c_gain_col(c) = &
                     c13_cnveg_carbonflux_inst%dwt_prod100c_gain_col(c) - prod100_c13flux(p)/dt
             endif
             
             if ( use_c14 ) then
                ! C14 column-level flux updates
                c14_cnveg_carbonflux_inst%dwt_conv_cflux_col(c) = &
                     c14_cnveg_carbonflux_inst%dwt_conv_cflux_col(c) - conv_c14flux(p)/dt
                c14_cnveg_carbonflux_inst%dwt_prod10c_gain_col(c) = &
                     c14_cnveg_carbonflux_inst%dwt_prod10c_gain_col(c) - prod10_c14flux(p)/dt
                c14_cnveg_carbonflux_inst%dwt_prod100c_gain_col(c) = &
                     c14_cnveg_carbonflux_inst%dwt_prod100c_gain_col(c) - prod100_c14flux(p)/dt
             endif
             
             ! column-level N flux updates
             cnveg_nitrogenflux_inst%dwt_conv_nflux_col(c) = &
                  cnveg_nitrogenflux_inst%dwt_conv_nflux_col(c) - conv_nflux(p)/dt
             cnveg_nitrogenflux_inst%dwt_prod10n_gain_col(c) = &
                  cnveg_nitrogenflux_inst%dwt_prod10n_gain_col(c) - prod10_nflux(p)/dt
             cnveg_nitrogenflux_inst%dwt_prod100n_gain_col(c) = &
                  cnveg_nitrogenflux_inst%dwt_prod100n_gain_col(c) - prod100_nflux(p)/dt
             
          end if
       end do
    end do
    
    ! Deallocate patch-level flux arrays
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
    
   end subroutine dyn_cnbal_patch

end module dynConsBiogeochemMod
