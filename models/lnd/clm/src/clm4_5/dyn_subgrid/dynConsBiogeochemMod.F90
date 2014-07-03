module dynConsBiogeochemMod

  !---------------------------------------------------------------------------
  !
  ! !DESCRIPTION:
  ! Handle conservation of biogeochemical quantities (C & N) with dynamic land cover.
  !
  ! !USES:
  use clmtype
  use decompMod  , only : bounds_type
  use abortutils , only : endrun
  use clm_varctl , only : iulog, use_c13, use_c14
  use shr_log_mod, only : errMsg => shr_log_errMsg
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  save

  public :: dyn_cnbal_pft

contains

  !-----------------------------------------------------------------------
  subroutine dyn_cnbal_pft(bounds, prior_weights)
    !
    ! !DESCRIPTION:
    ! Modify pft-level state and flux variables to maintain carbon and nitrogen balance with
    ! dynamic pft-weights.
    !
    ! !USES:
    use shr_const_mod      , only : SHR_CONST_PDB
    use clm_varcon         , only : istsoil
    use clm_varpar         , only : numveg, nlevdecomp, max_pft_per_col
    use clm_varcon         , only : istcrop
    use pftvarcon          , only : pconv, pprod10, pprod100
    use clm_varcon         , only : c13ratio, c14ratio
    use clm_time_manager   , only : get_step_size
    use dynPriorWeightsMod , only : prior_weights_type
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in) :: bounds        ! bounds
    type(prior_weights_type) , intent(in) :: prior_weights ! weights prior to the subgrid weight updates
    !
    ! !LOCAL VARIABLES:
    integer  :: pi,p,c,l,g,j    ! indices
    integer  :: ier           ! error code
    real(r8) :: dwt           ! change in pft weight (relative to column)
    real(r8) :: dt            ! land model time step (sec)
    real(r8) :: init_h2ocan   ! initial canopy water mass
    real(r8) :: new_h2ocan    ! canopy water mass after weight shift
    real(r8), allocatable :: dwt_leafc_seed(:)       ! pft-level mass gain due to seeding of new area
    real(r8), allocatable :: dwt_leafn_seed(:)       ! pft-level mass gain due to seeding of new area
    real(r8), allocatable :: dwt_deadstemc_seed(:)       ! pft-level mass gain due to seeding of new area
    real(r8), allocatable :: dwt_deadstemn_seed(:)       ! pft-level mass gain due to seeding of new area
    real(r8), allocatable :: dwt_frootc_to_litter(:)       ! pft-level mass loss due to weight shift
    real(r8), allocatable :: dwt_livecrootc_to_litter(:)   ! pft-level mass loss due to weight shift
    real(r8), allocatable :: dwt_deadcrootc_to_litter(:)   ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_frootn_to_litter(:)       ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_livecrootn_to_litter(:)   ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_deadcrootn_to_litter(:)   ! pft-level mass loss due to weight shift
    real(r8), allocatable :: conv_cflux(:)         ! pft-level mass loss due to weight shift
    real(r8), allocatable :: prod10_cflux(:)       ! pft-level mass loss due to weight shift
    real(r8), allocatable :: prod100_cflux(:)      ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: conv_nflux(:)         ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod10_nflux(:)       ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod100_nflux(:)      ! pft-level mass loss due to weight shift
    real(r8) :: t1,t2,wt_new,wt_old
    real(r8) :: init_state, change_state, new_state
    real(r8) :: tot_leaf, pleaf, pstor, pxfer
    real(r8) :: leafc_seed, leafn_seed
    real(r8) :: deadstemc_seed, deadstemn_seed
    real(r8), pointer :: dwt_ptr0, dwt_ptr1, dwt_ptr2, dwt_ptr3, ptr
    character(len=32) :: subname='dyn_cbal' ! subroutine name
    !! C13
    real(r8), allocatable :: dwt_leafc13_seed(:)     ! pft-level mass gain due to seeding of new area
    real(r8), allocatable :: dwt_deadstemc13_seed(:)     ! pft-level mass gain due to seeding of new area
    real(r8), allocatable, target :: dwt_frootc13_to_litter(:)     ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_livecrootc13_to_litter(:) ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_deadcrootc13_to_litter(:) ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: conv_c13flux(:)       ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod10_c13flux(:)     ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod100_c13flux(:)    ! pft-level mass loss due to weight shift
    real(r8) :: c3_del13c     ! typical del13C for C3 photosynthesis (permil, relative to PDB)
    real(r8) :: c4_del13c     ! typical del13C for C4 photosynthesis (permil, relative to PDB)
    real(r8) :: c3_r1_c13         ! isotope ratio (13c/12c) for C3 photosynthesis
    real(r8) :: c4_r1_c13         ! isotope ratio (13c/12c) for C4 photosynthesis
    real(r8) :: c3_r2_c13         ! isotope ratio (13c/[12c+13c]) for C3 photosynthesis
    real(r8) :: c4_r2_c13         ! isotope ratio (13c/[12c+13c]) for C4 photosynthesis
    real(r8) :: leafc13_seed, deadstemc13_seed
    !! C14
    real(r8), allocatable :: dwt_leafc14_seed(:)     ! pft-level mass gain due to seeding of new area
    real(r8), allocatable :: dwt_deadstemc14_seed(:)     ! pft-level mass gain due to seeding of new area
    real(r8), allocatable, target :: dwt_frootc14_to_litter(:)     ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_livecrootc14_to_litter(:) ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: dwt_deadcrootc14_to_litter(:) ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: conv_c14flux(:)       ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod10_c14flux(:)     ! pft-level mass loss due to weight shift
    real(r8), allocatable, target :: prod100_c14flux(:)    ! pft-level mass loss due to weight shift
    real(r8) :: c3_del14c     ! typical del14C for C3 photosynthesis (permil, relative to PDB)
    real(r8) :: c4_del14c     ! typical del14C for C4 photosynthesis (permil, relative to PDB)
    real(r8) :: c3_r1_c14         ! isotope ratio (14c/12c) for C3 photosynthesis
    real(r8) :: c4_r1_c14         ! isotope ratio (14c/12c) for C4 photosynthesis
    real(r8) :: c3_r2_c14         ! isotope ratio (14c/[12c+14c]) for C3 photosynthesis
    real(r8) :: c4_r2_c14         ! isotope ratio (14c/[12c+14c]) for C4 photosynthesis
    real(r8) :: leafc14_seed, deadstemc14_seed
    !-----------------------------------------------------------------------
    
   associate(& 
   lfpftd  =>  pps%lfpftd  & ! Output:  [real(r8) (:)] F. Li and S. Levis                                        
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
            lfpftd(p)=-dwt
          ! PFTs for which weight increases on this timestep
          if (dwt > 0._r8) then
             
             ! first identify PFTs that are initiating on this timestep
             ! and set all the necessary state and flux variables
             if (prior_weights%pwtcol(p) == 0._r8) then
                
                ! set initial conditions for PFT that is being initiated
                ! in this time step.  Based on the settings in cnIniTimeVar.
                
                ! pft-level carbon state variables
                pcs%leafc(p)              = 0._r8
                pcs%leafc_storage(p)      = 0._r8
                pcs%leafc_xfer(p)         = 0._r8
                pcs%frootc(p)             = 0._r8
                pcs%frootc_storage(p)     = 0._r8
                pcs%frootc_xfer(p)        = 0._r8
                pcs%livestemc(p)          = 0._r8
                pcs%livestemc_storage(p)  = 0._r8
                pcs%livestemc_xfer(p)     = 0._r8
                pcs%deadstemc(p)          = 0._r8
                pcs%deadstemc_storage(p)  = 0._r8
                pcs%deadstemc_xfer(p)     = 0._r8
                pcs%livecrootc(p)         = 0._r8
                pcs%livecrootc_storage(p) = 0._r8
                pcs%livecrootc_xfer(p)    = 0._r8
                pcs%deadcrootc(p)         = 0._r8
                pcs%deadcrootc_storage(p) = 0._r8
                pcs%deadcrootc_xfer(p)    = 0._r8
                pcs%gresp_storage(p)      = 0._r8
                pcs%gresp_xfer(p)         = 0._r8
                pcs%cpool(p)              = 0._r8
                pcs%xsmrpool(p)           = 0._r8
                pcs%pft_ctrunc(p)         = 0._r8
                pcs%dispvegc(p)           = 0._r8
                pcs%storvegc(p)           = 0._r8
                pcs%totvegc(p)            = 0._r8
                pcs%totpftc(p)            = 0._r8
                
                if ( use_c13 ) then
                   ! pft-level carbon-13 state variables
                   pc13s%leafc(p)              = 0._r8
                   pc13s%leafc_storage(p)      = 0._r8
                   pc13s%leafc_xfer(p)         = 0._r8
                   pc13s%frootc(p)             = 0._r8
                   pc13s%frootc_storage(p)     = 0._r8
                   pc13s%frootc_xfer(p)        = 0._r8
                   pc13s%livestemc(p)          = 0._r8
                   pc13s%livestemc_storage(p)  = 0._r8
                   pc13s%livestemc_xfer(p)     = 0._r8
                   pc13s%deadstemc(p)          = 0._r8
                   pc13s%deadstemc_storage(p)  = 0._r8
                   pc13s%deadstemc_xfer(p)     = 0._r8
                   pc13s%livecrootc(p)         = 0._r8
                   pc13s%livecrootc_storage(p) = 0._r8
                   pc13s%livecrootc_xfer(p)    = 0._r8
                   pc13s%deadcrootc(p)         = 0._r8
                   pc13s%deadcrootc_storage(p) = 0._r8
                   pc13s%deadcrootc_xfer(p)    = 0._r8
                   pc13s%gresp_storage(p)      = 0._r8
                   pc13s%gresp_xfer(p)         = 0._r8
                   pc13s%cpool(p)              = 0._r8
                   pc13s%xsmrpool(p)           = 0._r8
                   pc13s%pft_ctrunc(p)         = 0._r8
                   pc13s%dispvegc(p)           = 0._r8
                   pc13s%storvegc(p)           = 0._r8
                   pc13s%totvegc(p)            = 0._r8
                   pc13s%totpftc(p)            = 0._r8
                endif
                
                if ( use_c14 ) then
                   ! pft-level carbon-14 state variables
                   pc14s%leafc(p)              = 0._r8
                   pc14s%leafc_storage(p)      = 0._r8
                   pc14s%leafc_xfer(p)         = 0._r8
                   pc14s%frootc(p)             = 0._r8
                   pc14s%frootc_storage(p)     = 0._r8
                   pc14s%frootc_xfer(p)        = 0._r8
                   pc14s%livestemc(p)          = 0._r8
                   pc14s%livestemc_storage(p)  = 0._r8
                   pc14s%livestemc_xfer(p)     = 0._r8
                   pc14s%deadstemc(p)          = 0._r8
                   pc14s%deadstemc_storage(p)  = 0._r8
                   pc14s%deadstemc_xfer(p)     = 0._r8
                   pc14s%livecrootc(p)         = 0._r8
                   pc14s%livecrootc_storage(p) = 0._r8
                   pc14s%livecrootc_xfer(p)    = 0._r8
                   pc14s%deadcrootc(p)         = 0._r8
                   pc14s%deadcrootc_storage(p) = 0._r8
                   pc14s%deadcrootc_xfer(p)    = 0._r8
                   pc14s%gresp_storage(p)      = 0._r8
                   pc14s%gresp_xfer(p)         = 0._r8
                   pc14s%cpool(p)              = 0._r8
                   pc14s%xsmrpool(p)           = 0._r8
                   pc14s%pft_ctrunc(p)         = 0._r8
                   pc14s%dispvegc(p)           = 0._r8
                   pc14s%storvegc(p)           = 0._r8
                   pc14s%totvegc(p)            = 0._r8
                   pc14s%totpftc(p)            = 0._r8
                endif
                
                ! pft-level nitrogen state variables
                pns%leafn(p)	           = 0._r8
                pns%leafn_storage(p)      = 0._r8
                pns%leafn_xfer(p)         = 0._r8
                pns%frootn(p)	           = 0._r8
                pns%frootn_storage(p)     = 0._r8
                pns%frootn_xfer(p)        = 0._r8
                pns%livestemn(p)	       = 0._r8
                pns%livestemn_storage(p)  = 0._r8
                pns%livestemn_xfer(p)     = 0._r8
                pns%deadstemn(p)	       = 0._r8
                pns%deadstemn_storage(p)  = 0._r8
                pns%deadstemn_xfer(p)     = 0._r8
                pns%livecrootn(p)         = 0._r8
                pns%livecrootn_storage(p) = 0._r8
                pns%livecrootn_xfer(p)    = 0._r8
                pns%deadcrootn(p)         = 0._r8
                pns%deadcrootn_storage(p) = 0._r8
                pns%deadcrootn_xfer(p)    = 0._r8
                pns%retransn(p)	       = 0._r8
                pns%npool(p)	           = 0._r8
                pns%pft_ntrunc(p)         = 0._r8
                pns%dispvegn(p)           = 0._r8
                pns%storvegn(p)           = 0._r8
                pns%totvegn(p)            = 0._r8
                pns%totpftn (p)           = 0._r8
                
                ! initialize same flux and epv variables that are set
                ! in CNiniTimeVar
                pcf%psnsun(p) = 0._r8
                pcf%psnsha(p) = 0._r8
                pps%laisun(p) = 0._r8
                pps%laisha(p) = 0._r8
                
                pepv%dormant_flag(p) = 1._r8
                pepv%days_active(p) = 0._r8
                pepv%onset_flag(p) = 0._r8
                pepv%onset_counter(p) = 0._r8
                pepv%onset_gddflag(p) = 0._r8
                pepv%onset_fdd(p) = 0._r8
                pepv%onset_gdd(p) = 0._r8
                pepv%onset_swi(p) = 0.0_r8
                pepv%offset_flag(p) = 0._r8
                pepv%offset_counter(p) = 0._r8
                pepv%offset_fdd(p) = 0._r8
                pepv%offset_swi(p) = 0._r8
                pepv%lgsf(p) = 0._r8
                pepv%bglfr(p) = 0._r8
                pepv%bgtr(p) = 0._r8
                ! difference from CNiniTimeVar: using column-level
                ! information to initialize annavg_t2m.
                pepv%annavg_t2m(p) = cps%cannavg_t2m(c)
                pepv%tempavg_t2m(p) = 0._r8
                pepv%gpp(p) = 0._r8
                pepv%availc(p) = 0._r8
                pepv%xsmrpool_recover(p) = 0._r8
                pepv%alloc_pnow(p) = 1._r8
                pepv%c_allometry(p) = 0._r8
                pepv%n_allometry(p) = 0._r8
                pepv%plant_ndemand(p) = 0._r8
                pepv%tempsum_potential_gpp(p) = 0._r8
                pepv%annsum_potential_gpp(p) = 0._r8
                pepv%tempmax_retransn(p) = 0._r8
                pepv%annmax_retransn(p) = 0._r8
                pepv%avail_retransn(p) = 0._r8
                pepv%plant_nalloc(p) = 0._r8
                pepv%plant_calloc(p) = 0._r8
                pepv%excess_cflux(p) = 0._r8
                pepv%downreg(p) = 0._r8
                pepv%prev_leafc_to_litter(p) = 0._r8
                pepv%prev_frootc_to_litter(p) = 0._r8
                pepv%tempsum_npp(p) = 0._r8
                pepv%annsum_npp(p) = 0._r8
                
                if ( use_c13 ) then
                   pc13f%psnsun(p) = 0._r8
                   pc13f%psnsha(p) = 0._r8
                   
                   pps%alphapsnsun(p) = 0._r8
                   pps%alphapsnsha(p) = 0._r8
                   
                   pepv%xsmrpool_c13ratio(p) = c13ratio
                   
                   pepv%rc13_canair(p) = 0._r8
                   pepv%rc13_psnsun(p) = 0._r8
                   pepv%rc13_psnsha(p) = 0._r8
                endif
                
                if ( use_c14 ) then
                   pc14f%psnsun(p) = 0._r8
                   pc14f%psnsha(p) = 0._r8
                   pepv%rc14_atm(p) = c14ratio
                   pepv%rc14_atm(p) = 0._r8
                endif
                
             end if  ! end initialization of new pft
             
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
             if (pft%itype(p) /= 0) then
                leafc_seed = 1._r8
                leafn_seed  = leafc_seed / pftcon%leafcn(pft%itype(p))
                if (pftcon%woody(pft%itype(p)) == 1._r8) then
                   deadstemc_seed = 0.1_r8
                   deadstemn_seed = deadstemc_seed / pftcon%deadwdcn(pft%itype(p))
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
                   
                   if (pftcon%c3psn(pft%itype(p)) == 1._r8) then
                      leafc13_seed     = leafc_seed     * c3_r2_c13
                      deadstemc13_seed = deadstemc_seed * c3_r2_c13
                   else
                      leafc13_seed     = leafc_seed     * c4_r2_c13
                      deadstemc13_seed = deadstemc_seed * c4_r2_c13
                   end if
                endif
                
                if ( use_c14 ) then
                   ! 14c state is initialized assuming initial "modern" 14C of 1.e-12
                   if (pftcon%c3psn(pft%itype(p)) == 1._r8) then
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
             
             tot_leaf = pcs%leafc(p) + pcs%leafc_storage(p) + pcs%leafc_xfer(p)
             pleaf = 0._r8
             pstor = 0._r8
             pxfer = 0._r8
             if (tot_leaf /= 0._r8) then
                ! when adding seed source to non-zero leaf state, use current proportions
                pleaf = pcs%leafc(p)/tot_leaf
                pstor = pcs%leafc_storage(p)/tot_leaf
                pxfer = pcs%leafc_xfer(p)/tot_leaf
             else
                ! when initiating from zero leaf state, use evergreen flag to set proportions
                if (pftcon%evergreen(pft%itype(p)) == 1._r8) then
                   pleaf = 1._r8
                else
                   pstor = 1._r8
                end if
             end if
             pcs%leafc(p)         = pcs%leafc(p)*t1         + leafc_seed*pleaf*t2
             pcs%leafc_storage(p) = pcs%leafc_storage(p)*t1 + leafc_seed*pstor*t2
             pcs%leafc_xfer(p)    = pcs%leafc_xfer(p)*t1    + leafc_seed*pxfer*t2
             pcs%frootc(p)  		   = pcs%frootc(p) 			* t1
             pcs%frootc_storage(p)     = pcs%frootc_storage(p) 	* t1
             pcs%frootc_xfer(p) 	   = pcs%frootc_xfer(p)		* t1
             pcs%livestemc(p)		   = pcs%livestemc(p)  		* t1
             pcs%livestemc_storage(p)  = pcs%livestemc_storage(p)  * t1
             pcs%livestemc_xfer(p)     = pcs%livestemc_xfer(p) 	* t1
             pcs%deadstemc(p)     = pcs%deadstemc(p)*t1     + deadstemc_seed*t2
             pcs%deadstemc_storage(p)  = pcs%deadstemc_storage(p)  * t1
             pcs%deadstemc_xfer(p)     = pcs%deadstemc_xfer(p) 	* t1
             pcs%livecrootc(p)  	   = pcs%livecrootc(p) 		* t1
             pcs%livecrootc_storage(p) = pcs%livecrootc_storage(p) * t1
             pcs%livecrootc_xfer(p)    = pcs%livecrootc_xfer(p)	* t1
             pcs%deadcrootc(p)  	   = pcs%deadcrootc(p) 		* t1
             pcs%deadcrootc_storage(p) = pcs%deadcrootc_storage(p) * t1
             pcs%deadcrootc_xfer(p)    = pcs%deadcrootc_xfer(p)	* t1
             pcs%gresp_storage(p)	   = pcs%gresp_storage(p)  	* t1
             pcs%gresp_xfer(p)  	   = pcs%gresp_xfer(p) 		* t1
             pcs%cpool(p)			   = pcs%cpool(p)  			* t1
             pcs%xsmrpool(p)		   = pcs%xsmrpool(p)			* t1
             pcs%pft_ctrunc(p)  	   = pcs%pft_ctrunc(p) 		* t1
             pcs%dispvegc(p)		   = pcs%dispvegc(p)			* t1
             pcs%storvegc(p)		   = pcs%storvegc(p)			* t1
             pcs%totvegc(p) 		   = pcs%totvegc(p)			* t1
             pcs%totpftc(p) 		   = pcs%totpftc(p)			* t1
             
             if ( use_c13 ) then
                ! pft-level carbon-13 state variables 
                tot_leaf = pc13s%leafc(p) + pc13s%leafc_storage(p) + pc13s%leafc_xfer(p)
                pleaf = 0._r8
                pstor = 0._r8
                pxfer = 0._r8
                if (tot_leaf /= 0._r8) then
                   pleaf = pc13s%leafc(p)/tot_leaf
                   pstor = pc13s%leafc_storage(p)/tot_leaf
                   pxfer = pc13s%leafc_xfer(p)/tot_leaf
                else
                   ! when initiating from zero leaf state, use evergreen flag to set proportions
                   if (pftcon%evergreen(pft%itype(p)) == 1._r8) then
                      pleaf = 1._r8
                   else
                      pstor = 1._r8
                   end if
                end if
                pc13s%leafc(p)         = pc13s%leafc(p)*t1         + leafc13_seed*pleaf*t2
                pc13s%leafc_storage(p) = pc13s%leafc_storage(p)*t1 + leafc13_seed*pstor*t2
                pc13s%leafc_xfer(p)    = pc13s%leafc_xfer(p)*t1    + leafc13_seed*pxfer*t2
                pc13s%frootc(p)			 = pc13s%frootc(p) 		* t1
                pc13s%frootc_storage(p)	         = pc13s%frootc_storage(p) 	* t1
                pc13s%frootc_xfer(p)		 = pc13s%frootc_xfer(p)		* t1
                pc13s%livestemc(p) 		 = pc13s%livestemc(p)  		* t1
                pc13s%livestemc_storage(p)          = pc13s%livestemc_storage(p)      * t1
                pc13s%livestemc_xfer(p)	         = pc13s%livestemc_xfer(p) 	* t1
                pc13s%deadstemc(p)                  = pc13s%deadstemc(p)*t1     + deadstemc13_seed*t2
                pc13s%deadstemc_storage(p)          = pc13s%deadstemc_storage(p)      * t1
                pc13s%deadstemc_xfer(p)	         = pc13s%deadstemc_xfer(p) 	* t1
                pc13s%livecrootc(p)		 = pc13s%livecrootc(p) 		* t1
                pc13s%livecrootc_storage(p)         = pc13s%livecrootc_storage(p)     * t1
                pc13s%livecrootc_xfer(p)	         = pc13s%livecrootc_xfer(p)	* t1
                pc13s%deadcrootc(p)		 = pc13s%deadcrootc(p) 		* t1
                pc13s%deadcrootc_storage(p)         = pc13s%deadcrootc_storage(p)     * t1
                pc13s%deadcrootc_xfer(p)	         = pc13s%deadcrootc_xfer(p)	* t1
                pc13s%gresp_storage(p) 	         = pc13s%gresp_storage(p)  	* t1
                pc13s%gresp_xfer(p)		 = pc13s%gresp_xfer(p) 		* t1
                pc13s%cpool(p) 			 = pc13s%cpool(p)  		* t1
                pc13s%xsmrpool(p)  		 = pc13s%xsmrpool(p)		* t1
                pc13s%pft_ctrunc(p)		 = pc13s%pft_ctrunc(p) 		* t1
                pc13s%dispvegc(p)  		 = pc13s%dispvegc(p)		* t1
                pc13s%storvegc(p)  		 = pc13s%storvegc(p)		* t1
                pc13s%totvegc(p)			 = pc13s%totvegc(p)		* t1
                pc13s%totpftc(p)			 = pc13s%totpftc(p)		* t1
                
             endif
             
             if ( use_c14 ) then
                ! pft-level carbon-14 state variables 
                tot_leaf = pc14s%leafc(p) + pc14s%leafc_storage(p) + pc14s%leafc_xfer(p)
                pleaf = 0._r8
                pstor = 0._r8
                pxfer = 0._r8
                if (tot_leaf /= 0._r8) then
                   pleaf = pc14s%leafc(p)/tot_leaf
                   pstor = pc14s%leafc_storage(p)/tot_leaf
                   pxfer = pc14s%leafc_xfer(p)/tot_leaf
                else
                   ! when initiating from zero leaf state, use evergreen flag to set proportions
                   if (pftcon%evergreen(pft%itype(p)) == 1._r8) then
                      pleaf = 1._r8
                   else
                      pstor = 1._r8
                   end if
                end if
                pc14s%leafc(p)         = pc14s%leafc(p)*t1         + leafc14_seed*pleaf*t2
                pc14s%leafc_storage(p) = pc14s%leafc_storage(p)*t1 + leafc14_seed*pstor*t2
                pc14s%leafc_xfer(p)    = pc14s%leafc_xfer(p)*t1    + leafc14_seed*pxfer*t2
                pc14s%frootc(p)			 = pc14s%frootc(p) 		* t1
                pc14s%frootc_storage(p)	         = pc14s%frootc_storage(p) 	* t1
                pc14s%frootc_xfer(p)		 = pc14s%frootc_xfer(p)		* t1
                pc14s%livestemc(p) 		 = pc14s%livestemc(p)  		* t1
                pc14s%livestemc_storage(p)          = pc14s%livestemc_storage(p)      * t1
                pc14s%livestemc_xfer(p)	         = pc14s%livestemc_xfer(p) 	* t1
                pc14s%deadstemc(p)                  = pc14s%deadstemc(p)*t1     + deadstemc14_seed*t2
                pc14s%deadstemc_storage(p)          = pc14s%deadstemc_storage(p)      * t1
                pc14s%deadstemc_xfer(p)	         = pc14s%deadstemc_xfer(p) 	* t1
                pc14s%livecrootc(p)		 = pc14s%livecrootc(p) 		* t1
                pc14s%livecrootc_storage(p)         = pc14s%livecrootc_storage(p)     * t1
                pc14s%livecrootc_xfer(p)	         = pc14s%livecrootc_xfer(p)	* t1
                pc14s%deadcrootc(p)		 = pc14s%deadcrootc(p) 		* t1
                pc14s%deadcrootc_storage(p)         = pc14s%deadcrootc_storage(p)     * t1
                pc14s%deadcrootc_xfer(p)	         = pc14s%deadcrootc_xfer(p)	* t1
                pc14s%gresp_storage(p) 	         = pc14s%gresp_storage(p)  	* t1
                pc14s%gresp_xfer(p)		 = pc14s%gresp_xfer(p) 		* t1
                pc14s%cpool(p) 			 = pc14s%cpool(p)  		* t1
                pc14s%xsmrpool(p)  		 = pc14s%xsmrpool(p)		* t1
                pc14s%pft_ctrunc(p)		 = pc14s%pft_ctrunc(p) 		* t1
                pc14s%dispvegc(p)  		 = pc14s%dispvegc(p)		* t1
                pc14s%storvegc(p)  		 = pc14s%storvegc(p)		* t1
                pc14s%totvegc(p)			 = pc14s%totvegc(p)		* t1
                pc14s%totpftc(p)			 = pc14s%totpftc(p)		* t1
             endif
             
             
             tot_leaf = pns%leafn(p) + pns%leafn_storage(p) + pns%leafn_xfer(p)
             pleaf = 0._r8
             pstor = 0._r8
             pxfer = 0._r8
             if (tot_leaf /= 0._r8) then
                pleaf = pns%leafn(p)/tot_leaf
                pstor = pns%leafn_storage(p)/tot_leaf
                pxfer = pns%leafn_xfer(p)/tot_leaf
             else
                ! when initiating from zero leaf state, use evergreen flag to set proportions
                if (pftcon%evergreen(pft%itype(p)) == 1._r8) then
                   pleaf = 1._r8
                else
                   pstor = 1._r8
                end if
             end if
             ! pft-level nitrogen state variables
             pns%leafn(p)         = pns%leafn(p)*t1         + leafn_seed*pleaf*t2
             pns%leafn_storage(p) = pns%leafn_storage(p)*t1 + leafn_seed*pstor*t2
             pns%leafn_xfer(p)    = pns%leafn_xfer(p)*t1    + leafn_seed*pxfer*t2
             pns%frootn(p)  		   = pns%frootn(p) 		* t1
             pns%frootn_storage(p)         = pns%frootn_storage(p) 	* t1
             pns%frootn_xfer(p) 	   = pns%frootn_xfer(p)		* t1
             pns%livestemn(p)		   = pns%livestemn(p)  		* t1
             pns%livestemn_storage(p)      = pns%livestemn_storage(p)      * t1
             pns%livestemn_xfer(p)         = pns%livestemn_xfer(p) 	* t1
             pns%deadstemn(p)              = pns%deadstemn(p)*t1     + deadstemn_seed*t2
             pns%deadstemn_storage(p)      = pns%deadstemn_storage(p)      * t1
             pns%deadstemn_xfer(p)         = pns%deadstemn_xfer(p) 	* t1
             pns%livecrootn(p)  	   = pns%livecrootn(p) 		* t1
             pns%livecrootn_storage(p)     = pns%livecrootn_storage(p)     * t1
             pns%livecrootn_xfer(p)        = pns%livecrootn_xfer(p)	* t1
             pns%deadcrootn(p)  	   = pns%deadcrootn(p) 		* t1
             pns%deadcrootn_storage(p)     = pns%deadcrootn_storage(p)     * t1
             pns%deadcrootn_xfer(p)        = pns%deadcrootn_xfer(p)        * t1
             pns%retransn(p)		   = pns%retransn(p)		* t1
             pns%npool(p)		   = pns%npool(p)  		* t1
             pns%pft_ntrunc(p)  	   = pns%pft_ntrunc(p)        	* t1
             pns%dispvegn(p)		   = pns%dispvegn(p)		* t1
             pns%storvegn(p)		   = pns%storvegn(p)		* t1
             pns%totvegn(p) 		   = pns%totvegn(p)		* t1
             pns%totpftn(p) 		   = pns%totpftn(p)		* t1
             
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
             wt_new = pft%wtcol(p)
             wt_old = prior_weights%pwtcol(p)
             
             !---------------
             ! C state update
             !---------------
             
             ! leafc 
             ptr => pcs%leafc(p)
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
             ptr => pcs%leafc_storage(p)
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
             ptr => pcs%leafc_xfer(p)
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
             ptr => pcs%frootc(p)
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
             ptr => pcs%frootc_storage(p)
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
             ptr => pcs%frootc_xfer(p)
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
             ptr => pcs%livestemc(p)
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
             ptr => pcs%livestemc_storage(p)
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
             ptr => pcs%livestemc_xfer(p)
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
             ptr => pcs%deadstemc(p)
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
             ptr => pcs%deadstemc_storage(p)
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
             ptr => pcs%deadstemc_xfer(p)
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
             ptr => pcs%livecrootc(p)
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
             ptr => pcs%livecrootc_storage(p)
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
             ptr => pcs%livecrootc_xfer(p)
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
             ptr => pcs%deadcrootc(p)
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
             ptr => pcs%deadcrootc_storage(p)
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
             ptr => pcs%deadcrootc_xfer(p)
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
             ptr => pcs%gresp_storage(p)
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
             ptr => pcs%gresp_xfer(p)
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
             ptr => pcs%cpool(p)
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
             ptr => pcs%xsmrpool(p)
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
             ptr => pcs%pft_ctrunc(p)
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
                ptr => pc13s%leafc(p)
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
                ptr => pc13s%leafc_storage(p)
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
                ptr => pc13s%leafc_xfer(p)
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
                ptr => pc13s%frootc(p)
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
                ptr => pc13s%frootc_storage(p)
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
                ptr => pc13s%frootc_xfer(p)
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
                ptr => pc13s%livestemc(p)
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
                ptr => pc13s%livestemc_storage(p)
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
                ptr => pc13s%livestemc_xfer(p)
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
                ptr => pc13s%deadstemc(p)
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
                ptr => pc13s%deadstemc_storage(p)
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
                ptr => pc13s%deadstemc_xfer(p)
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
                ptr => pc13s%livecrootc(p)
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
                ptr => pc13s%livecrootc_storage(p)
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
                ptr => pc13s%livecrootc_xfer(p)
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
                ptr => pc13s%deadcrootc(p)
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
                ptr => pc13s%deadcrootc_storage(p)
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
                ptr => pc13s%deadcrootc_xfer(p)
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
                ptr => pc13s%gresp_storage(p)
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
                ptr => pc13s%gresp_xfer(p)
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
                ptr => pc13s%cpool(p)
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
                ptr => pc13s%pft_ctrunc(p)
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
                ptr => pc14s%leafc(p)
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
                ptr => pc14s%leafc_storage(p)
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
                ptr => pc14s%leafc_xfer(p)
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
                ptr => pc14s%frootc(p)
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
                ptr => pc14s%frootc_storage(p)
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
                ptr => pc14s%frootc_xfer(p)
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
                ptr => pc14s%livestemc(p)
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
                ptr => pc14s%livestemc_storage(p)
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
                ptr => pc14s%livestemc_xfer(p)
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
                ptr => pc14s%deadstemc(p)
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
                ptr => pc14s%deadstemc_storage(p)
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
                ptr => pc14s%deadstemc_xfer(p)
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
                ptr => pc14s%livecrootc(p)
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
                ptr => pc14s%livecrootc_storage(p)
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
                ptr => pc14s%livecrootc_xfer(p)
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
                ptr => pc14s%deadcrootc(p)
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
                ptr => pc14s%deadcrootc_storage(p)
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
                ptr => pc14s%deadcrootc_xfer(p)
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
                ptr => pc14s%gresp_storage(p)
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
                ptr => pc14s%gresp_xfer(p)
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
                ptr => pc14s%cpool(p)
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
                ptr => pc14s%pft_ctrunc(p)
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
             ptr => pns%leafn(p)
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
             ptr => pns%leafn_storage(p)
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
             ptr => pns%leafn_xfer(p)
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
             ptr => pns%frootn(p)
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
             ptr => pns%frootn_storage(p)
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
             ptr => pns%frootn_xfer(p)
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
             ptr => pns%livestemn(p)
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
             ptr => pns%livestemn_storage(p)
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
             ptr => pns%livestemn_xfer(p)
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
             ptr => pns%deadstemn(p)
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
             ptr => pns%deadstemn_storage(p)
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
             ptr => pns%deadstemn_xfer(p)
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
             ptr => pns%livecrootn(p)
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
             ptr => pns%livecrootn_storage(p)
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
             ptr => pns%livecrootn_xfer(p)
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
             ptr => pns%deadcrootn(p)
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
             ptr => pns%deadcrootn_storage(p)
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
             ptr => pns%deadcrootn_xfer(p)
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
             ptr => pns%retransn(p)
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
             ptr => pns%npool(p)
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
             ptr => pns%pft_ntrunc(p)
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
    end do               ! pft loop
    
    ! calculate column-level seeding fluxes
    do pi = 1,max_pft_per_col
       do c = bounds%begc, bounds%endc
          if ( pi <=  col%npfts(c) ) then
             p = col%pfti(c) + pi - 1
             
             ! C fluxes
             ccf%dwt_seedc_to_leaf(c) = ccf%dwt_seedc_to_leaf(c) + dwt_leafc_seed(p)/dt
             ccf%dwt_seedc_to_deadstem(c) = ccf%dwt_seedc_to_deadstem(c) &
                  + dwt_deadstemc_seed(p)/dt
             
             if ( use_c13 ) then
                cc13f%dwt_seedc_to_leaf(c) = cc13f%dwt_seedc_to_leaf(c) + dwt_leafc13_seed(p)/dt
                cc13f%dwt_seedc_to_deadstem(c) = cc13f%dwt_seedc_to_deadstem(c) &
                     + dwt_deadstemc13_seed(p)/dt
             endif
             
             if ( use_c14 ) then	
                cc14f%dwt_seedc_to_leaf(c) = cc14f%dwt_seedc_to_leaf(c) + dwt_leafc14_seed(p)/dt
                cc14f%dwt_seedc_to_deadstem(c) = cc14f%dwt_seedc_to_deadstem(c) &
                     + dwt_deadstemc14_seed(p)/dt
             endif
             
             ! N fluxes
             cnf%dwt_seedn_to_leaf(c) = cnf%dwt_seedn_to_leaf(c) + dwt_leafn_seed(p)/dt
             cnf%dwt_seedn_to_deadstem(c) = cnf%dwt_seedn_to_deadstem(c) &
                  + dwt_deadstemn_seed(p)/dt
          end if
       end do
    end do
    
    
    ! calculate pft-to-column for fluxes into litter and CWD pools
    do j = 1, nlevdecomp
       do pi = 1,max_pft_per_col
          do c = bounds%begc, bounds%endc
             if ( pi <=  col%npfts(c) ) then
                p = col%pfti(c) + pi - 1
                
                ! fine root litter carbon fluxes
                ccf%dwt_frootc_to_litr_met_c(c,j) = ccf%dwt_frootc_to_litr_met_c(c,j) + &
                     (dwt_frootc_to_litter(p)*pftcon%fr_flab(pft%itype(p)))/dt * pps%froot_prof(p,j)
                ccf%dwt_frootc_to_litr_cel_c(c,j) = ccf%dwt_frootc_to_litr_cel_c(c,j) + &
                     (dwt_frootc_to_litter(p)*pftcon%fr_fcel(pft%itype(p)))/dt * pps%froot_prof(p,j)
                ccf%dwt_frootc_to_litr_lig_c(c,j) = ccf%dwt_frootc_to_litr_lig_c(c,j) + &
                     (dwt_frootc_to_litter(p)*pftcon%fr_flig(pft%itype(p)))/dt * pps%froot_prof(p,j)
                
                
                ! fine root litter nitrogen fluxes
                cnf%dwt_frootn_to_litr_met_n(c,j) = cnf%dwt_frootn_to_litr_met_n(c,j) + &
                     (dwt_frootn_to_litter(p)*pftcon%fr_flab(pft%itype(p)))/dt * pps%froot_prof(p,j)
                cnf%dwt_frootn_to_litr_cel_n(c,j) = cnf%dwt_frootn_to_litr_cel_n(c,j) + &
                     (dwt_frootn_to_litter(p)*pftcon%fr_fcel(pft%itype(p)))/dt * pps%froot_prof(p,j)
                cnf%dwt_frootn_to_litr_lig_n(c,j) = cnf%dwt_frootn_to_litr_lig_n(c,j) + &
                     (dwt_frootn_to_litter(p)*pftcon%fr_flig(pft%itype(p)))/dt * pps%froot_prof(p,j)
                
                ! livecroot fluxes to cwd
                ccf%dwt_livecrootc_to_cwdc(c,j) = ccf%dwt_livecrootc_to_cwdc(c,j) + &
                     (dwt_livecrootc_to_litter(p))/dt * pps%croot_prof(p,j)
                cnf%dwt_livecrootn_to_cwdn(c,j) = cnf%dwt_livecrootn_to_cwdn(c,j) + &
                     (dwt_livecrootn_to_litter(p))/dt * pps%croot_prof(p,j)
                
                ! deadcroot fluxes to cwd
                ccf%dwt_deadcrootc_to_cwdc(c,j) = ccf%dwt_deadcrootc_to_cwdc(c,j) + &
                     (dwt_deadcrootc_to_litter(p))/dt * pps%croot_prof(p,j)
                cnf%dwt_deadcrootn_to_cwdn(c,j) = cnf%dwt_deadcrootn_to_cwdn(c,j) + &
                     (dwt_deadcrootn_to_litter(p))/dt * pps%croot_prof(p,j)
             
                if ( use_c13 ) then
                   ! C13 fine root litter fluxes
                   cc13f%dwt_frootc_to_litr_met_c(c,j) = cc13f%dwt_frootc_to_litr_met_c(c,j) + &
                        (dwt_frootc13_to_litter(p)*pftcon%fr_flab(pft%itype(p)))/dt * pps%froot_prof(p,j)
                   cc13f%dwt_frootc_to_litr_cel_c(c,j) = cc13f%dwt_frootc_to_litr_cel_c(c,j) + &
                        (dwt_frootc13_to_litter(p)*pftcon%fr_fcel(pft%itype(p)))/dt * pps%froot_prof(p,j)
                   cc13f%dwt_frootc_to_litr_lig_c(c,j) = cc13f%dwt_frootc_to_litr_lig_c(c,j) + &
                        (dwt_frootc13_to_litter(p)*pftcon%fr_flig(pft%itype(p)))/dt * pps%froot_prof(p,j)
                   ! livecroot fluxes to cwd
                   cc13f%dwt_livecrootc_to_cwdc(c,j) = cc13f%dwt_livecrootc_to_cwdc(c,j) + &
                        (dwt_livecrootc13_to_litter(p))/dt * pps%croot_prof(p,j)
                   ! deadcroot fluxes to cwd
                   cc13f%dwt_deadcrootc_to_cwdc(c,j) = cc13f%dwt_deadcrootc_to_cwdc(c,j) + &
                        (dwt_deadcrootc13_to_litter(p))/dt * pps%croot_prof(p,j)
                   
                endif
                
                if ( use_c14 ) then                   
                   ! C14 fine root litter fluxes
                   cc14f%dwt_frootc_to_litr_met_c(c,j) = cc14f%dwt_frootc_to_litr_met_c(c,j) + &
                        (dwt_frootc14_to_litter(p)*pftcon%fr_flab(pft%itype(p)))/dt * pps%froot_prof(p,j)
                   cc14f%dwt_frootc_to_litr_cel_c(c,j) = cc14f%dwt_frootc_to_litr_cel_c(c,j) + &
                        (dwt_frootc14_to_litter(p)*pftcon%fr_fcel(pft%itype(p)))/dt * pps%froot_prof(p,j)
                   cc14f%dwt_frootc_to_litr_lig_c(c,j) = cc14f%dwt_frootc_to_litr_lig_c(c,j) + &
                        (dwt_frootc14_to_litter(p)*pftcon%fr_flig(pft%itype(p)))/dt * pps%froot_prof(p,j)
                   ! livecroot fluxes to cwd
                   cc14f%dwt_livecrootc_to_cwdc(c,j) = cc14f%dwt_livecrootc_to_cwdc(c,j) + &
                        (dwt_livecrootc14_to_litter(p))/dt * pps%croot_prof(p,j)
                   ! deadcroot fluxes to cwd
                   cc14f%dwt_deadcrootc_to_cwdc(c,j) = cc14f%dwt_deadcrootc_to_cwdc(c,j) + &
                        (dwt_deadcrootc14_to_litter(p))/dt * pps%croot_prof(p,j)
                endif
                
             end if
          end do
       end do
    end do
    ! calculate pft-to-column for fluxes into product pools and conversion flux
    do pi = 1,max_pft_per_col
       do c = bounds%begc,bounds%endc
          if (pi <= col%npfts(c)) then
             p = col%pfti(c) + pi - 1
             
             ! column-level fluxes are accumulated as positive fluxes.
             ! column-level C flux updates
             ccf%dwt_conv_cflux(c) = ccf%dwt_conv_cflux(c) - conv_cflux(p)/dt
             ccf%dwt_prod10c_gain(c) = ccf%dwt_prod10c_gain(c) - prod10_cflux(p)/dt
             ccf%dwt_prod100c_gain(c) = ccf%dwt_prod100c_gain(c) - prod100_cflux(p)/dt

             ! These magic constants should be replaced with: nbrdlf_evr_trp_tree and nbrdlf_dcd_trp_tree
             if(pft%itype(p)==4.or.pft%itype(p)==6)then
                ccf%lf_conv_cflux(c) = ccf%lf_conv_cflux(c) - conv_cflux(p)/dt
             end if
             
             if ( use_c13 ) then
                ! C13 column-level flux updates
                cc13f%dwt_conv_cflux(c) = cc13f%dwt_conv_cflux(c) - conv_c13flux(p)/dt
                cc13f%dwt_prod10c_gain(c) = cc13f%dwt_prod10c_gain(c) - prod10_c13flux(p)/dt
                cc13f%dwt_prod100c_gain(c) = cc13f%dwt_prod100c_gain(c) - prod100_c13flux(p)/dt
             endif
             
             if ( use_c14 ) then
                ! C14 column-level flux updates
                cc14f%dwt_conv_cflux(c) = cc14f%dwt_conv_cflux(c) - conv_c14flux(p)/dt
                cc14f%dwt_prod10c_gain(c) = cc14f%dwt_prod10c_gain(c) - prod10_c14flux(p)/dt
                cc14f%dwt_prod100c_gain(c) = cc14f%dwt_prod100c_gain(c) - prod100_c14flux(p)/dt
             endif
             
             ! column-level N flux updates
             cnf%dwt_conv_nflux(c) = cnf%dwt_conv_nflux(c) - conv_nflux(p)/dt
             cnf%dwt_prod10n_gain(c) = cnf%dwt_prod10n_gain(c) - prod10_nflux(p)/dt
             cnf%dwt_prod100n_gain(c) = cnf%dwt_prod100n_gain(c) - prod100_nflux(p)/dt
             
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
   end subroutine dyn_cnbal_pft


end module dynConsBiogeochemMod
