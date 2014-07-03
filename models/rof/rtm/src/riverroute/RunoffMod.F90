module RunoffMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: RunoffMod
!
! !DESCRIPTION:
! Module containing utilities for history file and coupler runoff data
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use mct_mod
  use RtmVar         , only : iulog, spval
  use rtm_cpl_indices, only : nt_rtm

! !PUBLIC TYPES:
  implicit none
  private

  public :: runoff_flow
  type runoff_flow
     !    - local initialization
     real(r8), pointer :: lonc(:)          ! lon of cell
     real(r8), pointer :: latc(:)          ! lat of cell
     real(r8), pointer :: area(:)          ! area of cell
     integer , pointer :: gindex(:)        ! global index
     integer , pointer :: dsi(:)           ! downstream index

     !    - local runtime
     real(r8), pointer :: runoff(:,:)      ! RTM flow (m**3 H2O/s)
     real(r8), pointer :: runofflnd(:,:)   ! runoff masked for land (m**3 H2O/s)
     real(r8), pointer :: runoffocn(:,:)   ! runoff masked for ocn  (m**3 H2O/s)
     real(r8), pointer :: dvolrdt(:,:)     ! RTM change in storage (mm/s)
     real(r8), pointer :: dvolrdtlnd(:,:)  ! dvolrdt masked for land (mm/s)
     real(r8), pointer :: dvolrdtocn(:,:)  ! dvolrdt masked for ocn  (mm/s)
     real(r8), pointer :: volr(:,:)        ! RTM storage (m**3)
     real(r8), pointer :: volrlnd(:,:)     ! RTM storage masked for land (m**3)
     real(r8), pointer :: fluxout(:,:)     ! RTM cell tracer outlflux (m3/s)
     real(r8), pointer :: fthresh(:)       ! RTM water flood threshold
     real(r8), pointer :: flood(:)         ! RTM water (flood) sent back to clm (mm/s)

     !    - global 
     integer , pointer :: mask(:)          ! mask of cell 0=none, 1=lnd, 2=ocn
     real(r8), pointer :: rlon(:)          ! rtm longitude list, 1d
     real(r8), pointer :: rlat(:)          ! rtm latitude list, 1d
     real(r8)          :: totarea          ! global area
     integer           :: numr             ! rtm gdc global number of cells
     integer           :: numrl            ! rtm gdc global number of lnd cells
     integer           :: numro            ! rtm gdc global number of ocn cells

     !    - local
     integer           :: begr,endr        ! local start/stop indices
     integer           :: lnumr            ! local number of cells

     !    - 1d field pointers for history files (currently needed)
     real(r8), pointer :: runofflnd_nt1(:)
     real(r8), pointer :: runofflnd_nt2(:)
     real(r8), pointer :: runoffocn_nt1(:)
     real(r8), pointer :: runoffocn_nt2(:)
     real(r8), pointer :: dvolrdtlnd_nt1(:)
     real(r8), pointer :: dvolrdtlnd_nt2(:)
     real(r8), pointer :: dvolrdtocn_nt1(:)
     real(r8), pointer :: dvolrdtocn_nt2(:)
     real(r8), pointer :: volr_nt1(:)
     real(r8), pointer :: volr_nt2(:)
  end type runoff_flow
  !
  type (runoff_flow), public :: runoff

  public :: RunoffInit

contains

  subroutine RunoffInit(begr, endr, numr)

    integer, intent(in) :: begr, endr, numr

    integer :: ier

    allocate(runoff%runoff(begr:endr,nt_rtm),     &
             runoff%dvolrdt(begr:endr,nt_rtm),    &
             runoff%runofflnd(begr:endr,nt_rtm),  &
             runoff%dvolrdtlnd(begr:endr,nt_rtm), &
             runoff%runoffocn(begr:endr,nt_rtm),  &
             runoff%dvolrdtocn(begr:endr,nt_rtm), &
             runoff%area(begr:endr),              &
             runoff%volr(begr:endr,nt_rtm),       &
             runoff%volrlnd(begr:endr,nt_rtm),    &
             runoff%fluxout(begr:endr,nt_rtm),    &
             runoff%lonc(begr:endr),              &
             runoff%latc(begr:endr),              &
             runoff%dsi(begr:endr),               &
             runoff%runofflnd_nt1(begr:endr),     &
             runoff%runofflnd_nt2(begr:endr),     &
             runoff%runoffocn_nt1(begr:endr),     &
             runoff%runoffocn_nt2(begr:endr),     &
             runoff%volr_nt1(begr:endr),          &
             runoff%volr_nt2(begr:endr),          &
             runoff%dvolrdtlnd_nt1(begr:endr),    &
             runoff%dvolrdtlnd_nt2(begr:endr),    &
             runoff%dvolrdtocn_nt1(begr:endr),    &
             runoff%dvolrdtocn_nt2(begr:endr),    &
             runoff%mask(numr),                   &
             runoff%gindex(begr:endr),            &
             runoff%fthresh(begr:endr),           &
             runoff%flood(begr:endr),             &
             stat=ier)
    if (ier /= 0) then
       write(iulog,*)'Rtmini ERROR allocation of runoff local arrays'
       call shr_sys_abort
    end if

    runoff%runoff(:,:)     = 0._r8
    runoff%runofflnd(:,:)  = spval
    runoff%runoffocn(:,:)  = spval
    runoff%dvolrdt(:,:)    = 0._r8
    runoff%dvolrdtlnd(:,:) = spval
    runoff%dvolrdtocn(:,:) = spval
    runoff%volr(:,:)       = 0._r8
    runoff%volrlnd(:,:)    = 0._r8
    runoff%volr_nt1(:)     = 0._r8
    runoff%volr_nt2(:)     = 0._r8
    runoff%flood(:)        = 0._r8

  end subroutine RunoffInit

end module RunoffMod
