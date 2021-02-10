module mmf_surface_mod
  
  use shr_kind_mod,     only: r8 => shr_kind_r8
  ! use physics_buffer,   only: physics_buffer_desc
  use physics_types,    only: physics_state, physics_ptend, physics_ptend_init
  use camsrfexch,       only: cam_in_t
  use constituents,     only: pcnst, cnst_type
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
  integer  :: ncol, icol, m
  real(r8) :: g_dp
  !-----------------------------------------------------------------------------
  ncol   = state%ncol

  lq(:) = .TRUE.
  call physics_ptend_init(ptend, state%psetcols, "mmf_surface_ac", ls=.true., lu=.true., lv=.true., lq=lq)

  do icol = 1,ncol

    g_dp = gravit * state%rpdel(icol,pver)

#ifdef MMF_SFC1
    ! all fluxes added in ac
    ptend%s(icol,pver)   = g_dp * cam_in%shf(icol)
    ptend%q(icol,pver,1) = g_dp * cam_in%cflx(icol,1)
    ptend%u(icol,pver)   = g_dp * cam_in%wsx(icol)
    ptend%v(icol,pver)   = g_dp * cam_in%wsy(icol)
#endif

#ifdef MMF_SFC2
    ptend%u(icol,pver)   = g_dp * cam_in%wsx(icol)
    ptend%v(icol,pver)   = g_dp * cam_in%wsy(icol)
#endif

#ifdef MMF_SFC3
    ! momentum fluxes in bc for this option
#endif

#ifdef MMF_SFC4
    ptend%u(icol,pver)   = g_dp * cam_in%wsx(icol) * 0.5
    ptend%v(icol,pver)   = g_dp * cam_in%wsy(icol) * 0.5
#endif

#ifdef MMF_SFC5
    ptend%q(icol,pver,1) = g_dp * cam_in%cflx(icol,1)
    ptend%u(icol,pver)   = g_dp * cam_in%wsx(icol)
    ptend%v(icol,pver)   = g_dp * cam_in%wsy(icol)
#endif

#ifdef MMF_SFC6
    ptend%q(icol,pver,1) = g_dp * cam_in%cflx(icol,1)
    ptend%u(icol,pver)   = g_dp * cam_in%wsx(icol) * 0.5
    ptend%v(icol,pver)   = g_dp * cam_in%wsy(icol) * 0.5
#endif

#ifdef MMF_SFC7
    ! all fluxes added with parallel split
    ptend%s(icol,pver)   = g_dp * cam_in%shf(icol)    * 0.5
    ptend%q(icol,pver,1) = g_dp * cam_in%cflx(icol,1) * 0.5
    ptend%u(icol,pver)   = g_dp * cam_in%wsx(icol)    * 0.5
    ptend%v(icol,pver)   = g_dp * cam_in%wsy(icol)    * 0.5
#endif

    do m = 2,pcnst
      ptend%q(icol,pver,m) = g_dp * cam_in%cflx(icol,m)
      ! Convert tendencies of dry constituents to dry basis
      if (cnst_type(m).eq.'dry') then
        ptend%q(icol,pver,m) = ptend%q(icol,pver,m) * state%pdel(icol,pver) / state%pdeldry(icol,pver)
      endif
    end do ! m

  end do ! icol

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
  integer  :: ncol, icol, m
  real(r8) :: g_dp
  !-----------------------------------------------------------------------------
  ncol   = state%ncol

  lq(:) = .TRUE.
  call physics_ptend_init(ptend, state%psetcols, "mmf_surface_bc", ls=.true., lu=.true., lv=.true., lq=lq)

  do icol = 1,ncol

    g_dp = gravit * state%rpdel(icol,pver)

#ifdef MMF_SFC1
    ! all fluxes added in ac
#endif

#ifdef MMF_SFC2
    ptend%s(icol,pver)   = g_dp * cam_in%shf(icol)
    ptend%q(icol,pver,1) = g_dp * cam_in%cflx(icol,1)
    ! momentum fluxes in ac for this option
#endif

#ifdef MMF_SFC3
    ptend%s(icol,pver)   = g_dp * cam_in%shf(icol)
    ptend%q(icol,pver,1) = g_dp * cam_in%cflx(icol,1)
    ptend%u(icol,pver)   = g_dp * cam_in%wsx(icol)
    ptend%v(icol,pver)   = g_dp * cam_in%wsy(icol)
#endif

#ifdef MMF_SFC4
    ptend%s(icol,pver)   = g_dp * cam_in%shf(icol)
    ptend%q(icol,pver,1) = g_dp * cam_in%cflx(icol,1)
    ptend%u(icol,pver)   = g_dp * cam_in%wsx(icol) * 0.5
    ptend%v(icol,pver)   = g_dp * cam_in%wsy(icol) * 0.5
#endif

#ifdef MMF_SFC5
    ! only SHF added in bc - will not work with MMF_CRM_SFC_FLX
    ptend%s(icol,pver)   = g_dp * cam_in%shf(icol)
#endif

#ifdef MMF_SFC6
    ! combination of #4-5 - will not work with MMF_CRM_SFC_FLX
    ptend%s(icol,pver)   = g_dp * cam_in%shf(icol)
    ptend%u(icol,pver)   = g_dp * cam_in%wsx(icol) * 0.5
    ptend%v(icol,pver)   = g_dp * cam_in%wsy(icol) * 0.5
#endif

#ifdef MMF_SFC7
    ! all fluxes added with parallel split
    ptend%s(icol,pver)   = g_dp * cam_in%shf(icol)    * 0.5
    ptend%q(icol,pver,1) = g_dp * cam_in%cflx(icol,1) * 0.5
    ptend%u(icol,pver)   = g_dp * cam_in%wsx(icol)    * 0.5
    ptend%v(icol,pver)   = g_dp * cam_in%wsy(icol)    * 0.5
#endif

  end do ! icol


#ifdef MMF_CRM_SFC_FLX
  ! fluxes will be injected within the CRM, so undo flux addition
  do icol = 1,ncol
    ptend%s(icol,pver)   = 0
    ptend%q(icol,pver,1) = 0
  end do ! icol
#endif

end subroutine mmf_surface_bc
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

end module mmf_surface_mod