module mmf_surface_mod
  
  use shr_kind_mod,     only: r8 => shr_kind_r8
  ! use physics_buffer,   only: physics_buffer_desc
  use physics_types,    only: physics_state, physics_ptend, physics_ptend_init
  use camsrfexch,       only: cam_in_t
  use constituents,     only: pcnst
  use ppgrid,           only: pver, pcols
  use physconst,        only: gravit

  implicit none
  ! private      
  ! save

  public mmf_surface_ac
  public mmf_surface_bc

contains
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine mmf_surface_ac(state, cam_in, ptend)
  !-----------------------------------------------------------------------------
  ! Interface arguments
  type(physics_state), intent(in)  :: state
  type(cam_in_t),      intent(in)  :: cam_in
  type(physics_ptend), intent(out) :: ptend
  !-----------------------------------------------------------------------------
  ! Local variables
  logical  :: lq(pcnst)
  integer  :: ncol, i, m
  real(r8) :: g_dp
  !-----------------------------------------------------------------------------
  ncol   = state%ncol

  lq(2:) = .TRUE.
  call physics_ptend_init(ptend, state%psetcols, "mmf_surface", ls=.true., lu=.true., lv=.true., lq=lq)
  ! call physics_ptend_init(ptend, state%psetcols, "mmf_surface", lq=lq)

! #ifdef MMF_SFC1
  do i = 1,ncol
    g_dp = gravit * state%rpdel(i,pver)
    ! ptend%s(i,pver) = g_dp * cam_in%shf(i)
    ptend%u(i,pver) = g_dp * cam_in%wsx(i)
    ptend%v(i,pver) = g_dp * cam_in%wsy(i)
    do m = 2,pcnst
      ptend%q(i,pver,m) = g_dp * cam_in%cflx(i,m)
    end do
  end do
! #endif

! #ifdef MMF_SFC2
!   do i = 1,ncol
!     g_dp = gravit * state%rpdel(i,pver)
!     ! ptend%s(i,pver) = g_dp * cam_in%shf(i)
!     ! ptend%u(i,pver) = g_dp * cam_in%wsx(i)
!     ! ptend%v(i,pver) = g_dp * cam_in%wsy(i)
!     do m = 2,pcnst
!       ptend%q(i,pver,m) = g_dp * cam_in%cflx(i,m)
!     end do
!   end do
! #endif

! #ifdef MMF_SFC3
!   do i = 1,ncol
!     g_dp = gravit * state%rpdel(i,pver)
!     ! ptend%s(i,pver) = g_dp * cam_in%shf(i)
!     ptend%u(i,pver) = g_dp * cam_in%wsx(i) * 0.5
!     ptend%v(i,pver) = g_dp * cam_in%wsy(i) * 0.5
!     do m = 2,pcnst
!       ptend%q(i,pver,m) = g_dp * cam_in%cflx(i,m)
!     end do
!   end do
! #endif

end subroutine mmf_surface_ac
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine mmf_surface_bc(state, cam_in, ptend)
  !-----------------------------------------------------------------------------
  ! Interface arguments
  type(physics_state), intent(in)  :: state
  type(cam_in_t),      intent(in)  :: cam_in
  type(physics_ptend), intent(out) :: ptend
  !-----------------------------------------------------------------------------
  ! Local variables
  logical  :: lq(pcnst)
  integer  :: ncol, i, m
  real(r8) :: g_dp
  !-----------------------------------------------------------------------------
  ncol   = state%ncol

  lq(1) = .TRUE.
  call physics_ptend_init(ptend, state%psetcols, "mmf_surface", ls=.true., lu=.true., lv=.true., lq=lq)

! #ifdef MMF_SFC1
  do i = 1,ncol
    g_dp = gravit * state%rpdel(i,pver)
    ptend%s(i,pver) = g_dp * cam_in%shf(i)
    ! ptend%u(i,pver) = g_dp * cam_in%wsx(i)
    ! ptend%v(i,pver) = g_dp * cam_in%wsy(i)
    ptend%q(i,pver,1) = g_dp * cam_in%cflx(i,1)
    ! do m = 1,pcnst
    !   ptend%q(i,pver,m) = g_dp * cam_in%cflx(i,m)
    ! end do
  end do
! #endif

! #ifdef MMF_SFC2
!   do i = 1,ncol
!     g_dp = gravit * state%rpdel(i,pver)
!     ptend%s(i,pver) = g_dp * cam_in%shf(i)
!     ptend%u(i,pver) = g_dp * cam_in%wsx(i)
!     ptend%v(i,pver) = g_dp * cam_in%wsy(i)
!     ptend%q(i,pver,1) = g_dp * cam_in%cflx(i,1)
!     ! do m = 1,pcnst
!     !   ptend%q(i,pver,m) = g_dp * cam_in%cflx(i,m)
!     ! end do
!   end do
! #endif

! #ifdef MMF_SFC3
!   do i = 1,ncol
!     g_dp = gravit * state%rpdel(i,pver)
!     ptend%s(i,pver) = g_dp * cam_in%shf(i)
!     ptend%u(i,pver) = g_dp * cam_in%wsx(i) * 0.5
!     ptend%v(i,pver) = g_dp * cam_in%wsy(i) * 0.5
!     ptend%q(i,pver,1) = g_dp * cam_in%cflx(i,1)
!     ! do m = 1,pcnst
!     !   ptend%q(i,pver,m) = g_dp * cam_in%cflx(i,m)
!     ! end do
!   end do
! #endif

end subroutine mmf_surface_bc
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

end module mmf_surface_mod