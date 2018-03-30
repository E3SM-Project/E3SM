module simple_cloud_fraction

  use shr_kind_mod, only: r8=>shr_kind_r8
  use cam_abortutils, only: endrun
  use cam_logfile,    only: iulog

  implicit none
  private

  public :: smpl_frc

contains


  subroutine smpl_frc( q, ql, qsat, ast, rhu00, gbmrh, dastdrh, dlnastdrh, &
                       smpl_frc_schm, smpl_frc_cld, smpl_frc_clr, pcols, pver, ncol )

  integer, intent(in) :: pcols, pver, ncol
  integer, intent(in) :: smpl_frc_schm
  real(r8),intent(in) :: smpl_frc_cld
  real(r8),intent(in) :: smpl_frc_clr
 
  real(r8),intent(in)    :: q(pcols,pver)
  real(r8),intent(in)    :: ql(pcols,pver)
  real(r8),intent(in)    :: qsat(pcols,pver)
  real(r8),intent(inout) :: ast(pcols,pver)
  real(r8),intent(inout) :: gbmrh(pcols,pver)
  real(r8),intent(out)   :: dastdrh(pcols,pver)
  real(r8),intent(out)   :: rhu00 
  real(r8),intent(out)   :: dlnastdrh(pcols,pver)

  integer :: i,k

  real(r8) :: ztmp (pcols,pver)

  real(r8),parameter :: cnst_qlmin = 1.e-18_r8  ! for the simple constant/binary cloud fraction scheme
  real(r8) :: rhpert, rhlim, rhdif(pcols,pver)  ! for the Slingo formula of cloud fraction
  real(r8) :: rhmin, cldrh, dv                  ! for Park's pdf scheme of cloud fraction
  real(r8),parameter :: pi = 3.141592653589793
 !real(r8),parameter :: fmax= 0.999_r8 ! upper limit for the cloud fraction
  real(r8),parameter :: fmax= 1._r8 ! upper limit for the cloud fraction

  select case (smpl_frc_schm)
  case (0) ! Constant or binary

    where ( ql(:ncol,:pver) > cnst_qlmin )
      ast(:ncol,:pver) = smpl_frc_cld
    elsewhere
      ast(:ncol,:pver) = smpl_frc_clr
    end where

    rhu00 = 0._r8
 
  case (1) ! Slingo formula

    rhpert = 0._r8
    gbmrh(:ncol,:pver) = q(:ncol,:pver) /qsat(:ncol,:pver) *(1.0_r8+real(0,r8)*rhpert)

    rhlim = 0.8_r8
    rhdif(:ncol,:pver) = (gbmrh(:ncol,:pver) - rhlim)/(1.0_r8-rhlim)
    rhdif(:ncol,:pver) = max( min(rhdif(:ncol,:pver),1.0_r8), 0._r8)

    ! Cloud fraction f
    ast(:ncol,:pver) = rhdif(:ncol,:pver)**2
    ast(:ncol,:pver) = min(fmax,ast(:ncol,:pver))

    ! df/dRH
    dastdrh(:ncol,:pver) = 0._r8
    where( (ast(:ncol,:pver).gt. 0._r8) .and. (ast(:ncol,:pver).lt.fmax) )
      dastdrh(:ncol,:pver) = rhdif(:ncol,:pver)/(1.0_r8-rhlim)*2._r8
    end where

    ! dln(f)/dRH
    dlnastdrh(:ncol,:pver) = 0._r8
    where( (ast(:ncol,:pver).gt. 0._r8) .and. (ast(:ncol,:pver).lt.fmax) )
      dlnastdrh(:ncol,:pver) = 2._r8/rhdif(:ncol,:pver)/(1.0_r8-rhlim)
    end where

    rhu00 = rhlim

  case (2) ! Park et al. (2014) statistical scheme based on triangular PDF

    rhmin = 0.8_r8              ! critical gribox-mean RH (gbmrh) for st clouds to form
    cldrh = 1.0_r8              ! in-cloud RH
    dV    = cldrh - rhmin       ! half-width of triangular PDF

    gbmrh(:ncol,:pver) = ( q(:ncol,:pver) + ql(:ncol,:pver) ) /qsat(:ncol,:pver)

    do k=1,pver
       do i=1,ncol

       if( gbmrh(i,k) .ge. 1._r8 ) then

           ast(i,k) = 1._r8

       elseif( gbmrh(i,k) .gt. (cldrh-dV/6._r8) .and. gbmrh(i,k) .lt. 1._r8 ) then

           ast(i,k) = 1._r8 - (-3._r8/sqrt(2._r8)*(gbmrh(i,k)-cldrh)/dV)**(2._r8/3._r8)

       elseif( gbmrh(i,k) .gt. (cldrh-dV) .and. gbmrh(i,k) .le. (cldrh-dV/6._r8) ) then

           ast(i,k) = 4._r8*(cos((1._r8/3._r8)*(acos((3._r8/2._r8/sqrt(2._r8))* &
                      (1._r8+(gbmrh(i,k)-cldrh)/dV))-2._r8*3.141592_r8)))**2._r8

       elseif( gbmrh(i,k) .le. (cldrh-dV) ) then
           ast(i,k) = 0._r8
       endif

       end do
    end do
    rhu00 = rhmin

  case (3) ! like the slingo formula, but with a different functional form that is C0-continuous and C1-continuous

    rhpert = 0._r8
    gbmrh(:ncol,:pver) = q(:ncol,:pver) /qsat(:ncol,:pver) *(1.0_r8+real(0,r8)*rhpert)

    rhlim = 0.8_r8
    rhdif(:ncol,:pver) = (gbmrh(:ncol,:pver) - rhlim)/(1.0_r8-rhlim)
    rhdif(:ncol,:pver) =  max( min(rhdif(:ncol,:pver),1.0_r8), 0._r8)

    ztmp(:ncol,:pver) = ( rhdif(:ncol,:pver) - 0.5_r8 )*pi 

    ! Cloud fraction f
    ast (:ncol,:pver) = 0.5_r8*( sin(ztmp(:ncol,:pver)) + 1._r8 )
    ast (:ncol,:pver) = min(fmax, ast (:ncol,:pver) )

    ! df/dRH
    dastdrh(:ncol,:pver) = 0._r8
    where( (ast(:ncol,:pver).gt. 0._r8) .and. (ast(:ncol,:pver).lt.fmax) )
      dastdrh(:ncol,:pver) =  0.5_r8* cos(ztmp(:ncol,:pver)) *pi/(1.0_r8-rhlim)
    end where

    ! dln(f)/dRH
    dlnastdrh(:ncol,:pver) = 0._r8
    where( (ast(:ncol,:pver).gt. 0._r8) .and. (ast(:ncol,:pver).lt.fmax) )
      dlnastdrh(:ncol,:pver) = 1._r8/ast(:ncol,:pver) * dastdrh(:ncol,:pver) 
    end where

    rhu00 = rhlim

  case default
    write(iulog,*) "Unrecognized value of smpl_frc_schm:",smpl_frc_schm,". Abort."
    call endrun
  end select

  end subroutine smpl_frc

end module simple_cloud_fraction
