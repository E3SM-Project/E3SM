
module iop
!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: iop
! 
! !DESCRIPTION: 
! iop specific routines
!
! !USES:
!
  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_scam_mod, only: shr_scam_GetCloseLatLon
  use constituents, only: readtrace, cnst_get_ind, pcnst, cnst_name
  use string_utils, only: to_lower
  use pmgrid
  use prognostics
  use time_manager, only: timemgr_init, get_curr_date, get_curr_calday,&
                          get_nstep,get_start_date,timemgr_time_inc
  use cam_abortutils,   only: endrun
  use scamMod
  use wrap_nf
  use cam_logfile,  only: iulog
  use phys_control, only: phys_getopts
  use eul_control_mod,only: eul_nsplit
!
! !PUBLIC TYPES:
  implicit none

  private

  real(r8), allocatable, target :: dqfx3sav(:,:,:,:)       
  real(r8), allocatable, target :: t2sav(:,:,:)       
  real(r8), allocatable, target :: divq3dsav(:,:,:,:)
  real(r8), allocatable, target :: divt3dsav(:,:,:)       
  real(r8), allocatable, target :: betasav(:)
  integer :: closelatidx,closelonidx,latid,lonid,levid,timeid

  real(r8):: closelat,closelon

!
! !PUBLIC MEMBER FUNCTIONS:
  public :: init_iop_fields
!  public :: scam_use_iop_srf
! !PUBLIC DATA:
  public betasav, &
         dqfx3sav, divq3dsav, divt3dsav,t2sav

!
! !REVISION HISTORY:
! Created by John Truesdale
!
!EOP
!
! !PRIVATE MEMBER FUNCTIONS:
!----------------------------------------------------------------------- 

contains
   subroutine init_iop_fields()
!------------------------------------------------------------------------------
! Coupler for converting dynamics output variables into physics input variables
! also writes dynamics variables (on physics grid) to history file
!------------------------------------------------------------------------------
   implicit none

!-----------------------------------------------------------------------
   if (eul_nsplit>1) then
      call endrun('iop module cannot be used with eul_nsplit>1')
   endif
	        
   if(.not.allocated(betasav)) then
      allocate (betasav(beglat:endlat))
      betasav(:)=0._r8
   endif

   if(.not.allocated(dqfx3sav)) then
      allocate (dqfx3sav(plon,plev,pcnst,beglat:endlat))
      dqfx3sav(:,:,:,:)=0._r8
   endif
   if(.not.allocated(divq3dsav)) then
      allocate (divq3dsav(plon,plev,pcnst,beglat:endlat))
      divq3dsav(:,:,:,:)=0._r8
   endif
   if(.not.allocated(divt3dsav)) then
      allocate (divt3dsav(plon,plev,beglat:endlat))
      divt3dsav(:,:,:)=0._r8
   endif
   if(.not.allocated(t2sav)) then
      allocate (t2sav(plon,plev,beglat:endlat))  ! temp tendency
      t2sav(:,:,:)=0._r8
   endif
  end subroutine init_iop_fields

end module iop

