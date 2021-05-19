subroutine forecast(chunk, ps_in, u_in, v_in, &         ! In
                    t_in, q_in, scm_dt, t_phys_frc,&    ! In
                    t_fcst, q_fcst, u_fcst, v_fcst)     ! Out

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Eularian forecast of t, u, and v.   Advection terms are also converted
! to flux form and integrated to check conservation
! 
! Author: 
! Original version: Adopted from CAM3.5/CAM5
! Updated version for E3SM: Peter Bogenschutz (bogenschutz1@llnl.gov)
!
!-----------------------------------------------------------------------

  use shr_kind_mod,   only: r8 => shr_kind_r8, i8 => shr_kind_i8
  use pmgrid
  use constituents,   only: pcnst, cnst_get_ind
  use pspect
  use physconst,      only: rair,cpair
  use cam_logfile,    only: iulog
  use scamMod
  use cam_history,    only: outfld
!-----------------------------------------------------------------------
  implicit none
!-----------------------------------------------------------------------
!
! Input arguments
!
  integer, intent(in) :: chunk              ! begchunk identifier for output
  real(r8), intent(in) :: ps_in             ! surface pressure [Pa]
  real(r8), intent(in) :: u_in(plev)        ! zonal wind [m/s]
  real(r8), intent(in) :: v_in(plev)        ! meridional wind [m/s]
  real(r8), intent(in) :: t_in(plev)        ! temperature [K]
  real(r8), intent(in) :: q_in(plev,pcnst)  ! q tracer array [units vary]
  real(r8), intent(in) :: t_phys_frc(plev)  ! temperature forcing from physics [K/s]
  real(r8), intent(in) :: scm_dt            ! model time step [s]

! Output arguments
  real(r8), intent(out) :: t_fcst(plev)      ! updated temperature [K]
  real(r8), intent(out) :: q_fcst(plev,pcnst)! updated q tracer array [units vary]
  real(r8), intent(out) :: u_fcst(plev)      ! updated zonal wind [m/s]
  real(r8), intent(out) :: v_fcst(plev)      ! updated meridional wind [m/s]

!
!---------------------------Local variables-----------------------------
!
  real(r8) pmidm1(plev)  ! pressure at model levels
  real(r8) pintm1(plevp) ! pressure at model interfaces
  real(r8) pdelm1(plev)  ! pdel(k)   = pint  (k+1)-pint  (k)
  real(r8) t_lsf(plev)       ! storage for temperature large scale forcing
  real(r8) q_lsf(plev,pcnst) ! storage for moisture large scale forcing
  real(r8) weight,fac
  real(r8) wfldint(plevp)
  real(r8) t_expan

  integer i,k,m           ! longitude, level, constituent indices
  integer nlon
!
!  variables for relaxation addition
!
  real(r8) dist
  real(r8) denom
  real(r8) rtau(plev)
  real(r8) relaxt(plev)
  real(r8) relaxq(plev)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Perform weight averaged in pressure of the omega profile to get from
!   the midpoint grid to the interface grid

  nlon = 1 ! number of columns for plevs0 routine
  call plevs0(nlon    ,plon   ,plev    ,ps_in   ,pintm1  ,pmidm1 ,pdelm1)

  wfldint(1) = 0.0_r8

  do k=2,plev
     weight = (pintm1(k) - pmidm1(k-1))/(pmidm1(k) - pmidm1(k-1))
     wfldint(k) = (1.0_r8 - weight)*wfld(k-1) + weight*wfld(k)
  end do

  wfldint(plevp) = 0.0_r8

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!  Advance T and Q due to large scale forcing

  if (use_3dfrc) then
    t_lsf(:plev) = divt3d(:plev)
    q_lsf(:plev,:pcnst) = divq3d(:plev,:pcnst)
  else
    t_lsf(:plev) = divt(:plev)
    q_lsf(:plev,:pcnst) = divq(:plev,:pcnst)
  endif

  do k=1,plev
    ! Initialize thermal expansion term to zero.  This term is only
    !  considered if using the preq-x dycore and if three dimensional
    !  forcing is not provided by IOP forcing file.
    t_expan = 0._r8
#ifndef MODEL_THETA_L
    ! this term is already taken into account through
    !  LS vertical advection in theta-l dycore  
    if (.not. use_3dfrc) then 
      t_expan = scm_dt*wfld(k)*t_in(k)*rair/(cpair*pmidm1(k))
    endif
#endif  
    t_fcst(k) = t_in(k) + t_expan + scm_dt*(t_phys_frc(k) + t_lsf(k))
    do m=1,pcnst
      q_fcst(k,m) = q_in(k,m) + scm_dt*q_lsf(k,m)
    end do
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!  Set U and V fields

  if (have_v .and. have_u) then
    do k=1,plev
      u_fcst(k) = uobs(k)
      v_fcst(k) = vobs(k)
    enddo
  else
    do k=1,plev
      u_fcst(k) = u_in(k)
      v_fcst(k) = v_in(k)
    enddo
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!  Set U and V fields

  if (scm_relaxation) then
!
!    THIS IS WHERE WE RELAX THE SOLUTION IF REQUESTED
!    The relaxation can be thought of as a part of the "adjustment" physics
!
    do k=1,plev
      relaxt(k) = 0.0_r8
      relaxq(k) = 0.0_r8
    end do

    do k=1,plev

      if (pmidm1(k) .le. scm_relaxation_low*100._r8 .and. &
        pmidm1(k) .ge. scm_relaxation_high*100._r8) then

        rtau(k)   = 10800._r8          ! 3-hr adj. time scale
        rtau(k)   = max(scm_dt,rtau(k))
        relaxt(k) = -(t_fcst(k)   - tobs(k))/rtau(k)
        relaxq(k) = -(q_fcst(k,1) - qobs(k))/rtau(k)

        t_fcst(k)   = t_fcst(k)   + relaxt(k)*scm_dt
        q_fcst(k,1) = q_fcst(k,1) + relaxq(k)*scm_dt

      endif

    end do

    call outfld('TRELAX',relaxt,plon,chunk)
    call outfld('QRELAX',relaxq,plon,chunk)
    call outfld('TAURELAX',rtau,plon,chunk)
    
  end if ! end relaxation logical
     
!  evaluate the difference in state information from observed
  do k = 1, plev
    tdiff(k) = t_fcst(k)   - tobs(k)
    qdiff(k) = q_fcst(k,1) - qobs(k)
    udiff(k) = u_fcst(k)   - uobs(k)
    vdiff(k) = v_fcst(k)   - vobs(k)
  end do
!
!===============================================================
!
!  outfld calls related to SCM

  call outfld('TOBS',tobs,plon,chunk)
  call outfld('QOBS',qobs,plon,chunk)
  call outfld('TDIFF',tdiff,plon,chunk)
  call outfld('QDIFF',qdiff,plon,chunk)
  call outfld('DIVQ',divq,plon,chunk)
  call outfld('DIVT',divt,plon,chunk)
  call outfld('DIVQ3D',divq3d,plon,chunk)
  call outfld('DIVT3D',divt3d,plon,chunk)
  call outfld('PRECOBS',precobs,plon,chunk)
  call outfld('LHFLXOBS',lhflxobs,plon,chunk)
  call outfld('SHFLXOBS',shflxobs,plon,chunk)

end subroutine forecast


!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------


