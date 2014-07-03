module eul_control_mod
!----------------------------------------------------------------------- 
! 
! Purpose: control info for eul dynamics core
!
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!- use statements ------------------------------------------------------
!-----------------------------------------------------------------------
  use shr_kind_mod, only : r8=>shr_kind_r8
  use pmgrid,       only: plat, plon, plev
  use spmd_utils,   only: masterproc
  use dycore,       only: get_resolution
  use pspect,       only: pnmax
  use infnan,       only: posinf, assignment(=)
!-----------------------------------------------------------------------
!- module boilerplate --------------------------------------------------
!-----------------------------------------------------------------------
   implicit none
   private                   ! Make the default access private
   save

!---------------------------------------------------------------------------------
! Public interfaces
!---------------------------------------------------------------------------------
  public dyn_eul_readnl                    ! read dynamics namelist 

!-----------------------------------------------------------------------
! Public data ----------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Purpose: Global integrals for moisture and mass conservation and geopotential height
!-----------------------------------------------------------------------

   real(r8) ,public ::  tmass(plat)  ! Mass integral for each latitude pair
   real(r8) ,public ::  tmass0       ! Specified dry mass of atmosphere
   real(r8) ,public ::  tmassf       ! Global mass integral
   real(r8) ,public ::  qmassf       ! Global moisture integral
   real(r8) ,public ::  fixmas       ! Proportionality factor for ps in dry mass fixer
   real(r8) ,public ::  qmass1       ! Contribution to global moisture integral (mass
                                     !  weighting is based upon the "A" part of the hybrid grid)
   real(r8) ,public ::  qmass2       ! Contribution to global moisture integral (mass
                                     !  weighting is based upon the "B" part of the hybrid grid)
   real(r8) ,public ::  pdela(plon,plev)! pressure difference between interfaces (pressure
                                     !  defined using the "A" part of hybrid grid only)
   real(r8) ,public ::  zgsint       ! global integral of geopotential height

! from comfft.h

   integer  ,public :: pcray                   ! length of vector register (words)
   parameter (pcray=64)
!
   real(r8) ,public :: trig (3*plon/2+1,plat)  ! trigonometric funct values used by fft
   integer  ,public :: ifax(19,plat)           ! fft factorization of plon/2
   real(r8), public :: cnfac                   ! Courant num factor(multiply by max |V|)
   real(r8), public :: cnlim                   ! Maximum allowable courant number
   real(r8), public :: hdfsd2(pnmax)   	       ! Del^2 mult. for each wave (vort-div)
   real(r8), public :: hdfst2(pnmax)   	       ! Del^2 multiplier for each wave (t-q)
   real(r8), public :: hdfsd4(pnmax)   	       ! Del^4 mult. for each wave (vort-div)
   real(r8), public :: hdfst4(pnmax)   	       ! Del^4 multiplier for each wave (t-q)
   real(r8), public :: hdiftq(pnmax,plev)      ! Temperature-tracer diffusion factors
   real(r8), public :: hdifzd(pnmax,plev)      ! Vorticity-divergence diffusion factors
   integer,  public :: kmnhd4                  ! Top level for del^4 diffusion
   integer,  public :: kmxhd2                  ! Bottom level for increased del^2 diffusion
   integer,  public :: nindex(plev)            ! Starting index for spectral truncation
   integer,  public :: nmaxhd                  ! Maximum two dimensional wave number

!-----------------------------------------------------------------------
!
! Purpose: Spectral Namelist variables
!
!-----------------------------------------------------------------------

   real(r8), public :: dif2 = 2.5e5_r8  !del2 horizontal diffusion coeff.
   real(r8), public :: dif4             !del4 horizontal diffusion coeff.
   real(r8) ,public :: divdampn = 0._r8    !Number of days (from nstep 0) to run divergence
   real(r8) ,public :: eps = 0.06_r8    !time filter coefficient. Defaults to 0.06.
   integer,  public :: kmxhdc = 5       !number of levels (starting from model top) to apply Courant limiter.
   integer  ,public :: eul_nsplit = 1   !Intended number of dynamics timesteps per physics timestep
   
 contains
   
   subroutine dyn_eul_readnl(nlfile)
     
     ! Read dynamics namelist group.

!-----------------------------------------------------------------------
!- use statements ------------------------------------------------------
!-----------------------------------------------------------------------
     use abortutils,      only: endrun
     use namelist_utils,  only: find_group_name
     use units,           only: getunit, freeunit
     use mpishorthand
     use cam_logfile,     only: iulog

! args
     
     character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input
    
! local vars
     integer :: unitn, ierr
     
    namelist /dyn_spectral_inparm/ eul_nsplit,dif2,dif4,divdampn,eps,kmxhdc    
! 
! Read namelist 
!
    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'dyn_spectral_inparm', status=ierr)
       if (ierr == 0) then
          read(unitn, dyn_spectral_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun('eul_readnl: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    call mpibcast (divdampn,1,mpir8,0,mpicom)
    call mpibcast (eps     ,1,mpir8,0,mpicom)
    call mpibcast (dif2    ,1,mpir8,0,mpicom)
    call mpibcast (dif4    ,1,mpir8,0,mpicom)
    call mpibcast (kmxhdc  ,1,mpiint,0,mpicom)
    call mpibcast (eul_nsplit  ,1,mpiint,0,mpicom)
#endif

    if (divdampn > 0._r8) then
       write(iulog,*) 'Divergence damper for spectral dycore invoked for days 0. to ',divdampn,' of this case'
    elseif (divdampn < 0._r8) then
       call endrun ('READ_NAMELIST: divdampn must be a positive number')
    else
       write(iulog,*) 'Divergence damper for spectral dycore NOT invoked'
    endif

    ! Ensure dif4 is not used unset
    dif4 = posinf
    
    select case ( get_resolution() )
    case ( 'T5' )
       dif2 = 2.5e7_r8
       dif4 = 1.0e18_r8
    case ( 'T21' )
       dif2 = 2.5e5_r8
       dif4 = 2.0e16_r8
    case ( 'T31' )
       dif2 = 2.5e5_r8
       dif4 = 2.0e16_r8
    case ( 'T42' )
       dif2 = 2.5e5_r8
       dif4 = 1.0e16_r8
    case ( 'T85' )
       dif2 = 2.5e5_r8
       dif4 = 1.0e15_r8
    case ( 'T170' )
       dif2 = 2.5e5_r8
       dif4 = 1.5e14_r8
    end select

    if (kmxhdc >= plev .or. kmxhdc < 0) then
       call endrun ('DYN_EUL_READNL:  ERROR:  KMXHDC must be between 0 and plev-1')
    end if

!
! Write physics variables from namelist to std. output
!
    write(iulog,9108) eps,dif2,dif4,kmxhdc,eul_nsplit

9108 format(' Time filter coefficient (EPS)                 ',f10.3,/,&
            ' DEL2 Horizontal diffusion coefficient (DIF2)  ',e10.3/, &
            ' DEL4 Horizontal diffusion coefficient (DIF4)  ',e10.3/, &
            ' Number of levels Courant limiter applied      ',i10/,   &
            ' Dynamics Subcycling                           ',i10)

end subroutine dyn_eul_readnl

end module eul_control_mod
