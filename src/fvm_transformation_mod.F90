!***********************************************************************************!
! Author: Christoph Erath                                                           !
! Date: 18.March 2011                                                               !
! Collect all transformation using in home                                          !
! - vortex_rotatedsphere                                                            !
! - vortex_rotatedsphereback                                                        !
!***********************************************************************************!

module fvm_transformation_mod
  use kinds, only : real_kind
  use physical_constants, only : DD_PI
  
implicit none
private

public :: vortex_rotatedsphere
public :: vortex_rotatedsphereback
contains
!-----------------------------------------------------------------------------------!
! INPUT:  Transformation between a rotated system with north pole (lap,thp) and     !
!         the unrotated system (lau,thu)                                            !
! OUTPUT: coordinate system (lar,thr) of the rotated system,                        !
!         lar\in [-pi,pi]; thr\in [-pi/2,pi/2]                                      !
! Literature: *W. T. M. Verkley. The construction of barotropic modons on a         !
!              sphere. J. Atmos. Sci., 41, 2492-2504, 1984                          !
!             *R. D. Nair and C. Jablonowski. Moving Vortices on the Sphere: A Test !
!              Case for Horizontal Advection Problems. Mon. Wea. Rev, 136, 699-711, !
!              2008                                                                 !  
! Remark: if thr=pi/2, lar can be arbitrary (is pole, singularity)                  !
!-----------------------------------------------------------------------------------!

subroutine vortex_rotatedsphere(lap,thp,lau,thu,lar,thr)
  ! Rotate to new North Pole at (lap,thp)
  ! (lar,thr) are the rotated coordinates coorsponding to (lau,thu) in 
  ! the regular sphere
  implicit none
  real (kind=real_kind), intent(in)  :: lap,thp,lau,thu
  real (kind=real_kind), intent(out) :: lar,thr

  real (kind=real_kind) :: cost,sint,sinp,cosp
  real (kind=real_kind) ::  trm, trm1,trm2,trm3
  real (kind=real_kind) ::  pi2

  pi2=2.0D0*DD_PI

  sinp = sin(thp)
  cosp = cos(thp)
  cost = cos(thu)
  sint = sin(thu)

  trm  = cost * cos(lau- lap)
  trm1 = cost * sin(lau- lap)
  trm2 = sinp * trm  - cosp * sint
  trm3 = sinp * sint + cosp * trm

  lar = atan2(trm1,trm2)
  if (lar < 0.0D0 ) lar = lar + pi2
  if (lar > pi2)    lar = lar - pi2
  thr = asin(trm3)      

end subroutine vortex_rotatedsphere


!-----------------------------------------------------------------------------------!
! INPUT:  (BACK) Transformation between a rotated system with north pole (lap,thp)  !
!         with the coordinate system (lar, thr) to the unrotated system             !
! OUTPUT: coordinate system (lau,thu) of the unrotated system,                      !
!         rla\in [-pi,pi]; rth\in [-pi/2,pi/2]                                      !
! Literature: *W. T. M. Verkley. The construction of barotropic modons on a         !
!              sphere. J. Atmos. Sci., 41, 2492-2504, 1984                          !
!             *R. D. Nair and C. Jablonowski. Moving Vortices on the Sphere: A Test !
!              Case for Horizontal Advection Problems. Mon. Wea. Rev, 136, 699-711, !
!              2008                                                                 !
! Remark: if thr=pi/2, lar can be arbitrary (is pole, singularity)                  !
!-----------------------------------------------------------------------------------!

subroutine vortex_rotatedsphereback(lap,thp,lar,thr,lau,thu)

  implicit none
  real (kind=real_kind), intent(in)  :: lap,thp,lar,thr
  real (kind=real_kind), intent(out) :: lau,thu
  !
  real (kind=real_kind) :: cost,sint,cosp,sinp,clam,slam 
  real (kind=real_kind) :: trm, t1,t2,t3
  real (kind=real_kind) ::  pi2

  pi2=2.0D0*DD_PI

  !
  !* Back to unrotated system

  cost = cos(thr)
  sint = sin(thr)
  clam = cos(lar)
  slam = sin(lar)
  cosp = cos(thp)
  sinp = sin(thp)

  t1 = slam * cost
  t2 = sint*cosp + cost*clam*sinp
  t3 = sint*sinp - cost*clam*cosp
  lau =  lap+ atan2(t1,t2)
  if (lau < 0.0D0 )  lau = lau + pi2
  if (lau > pi2)     lau = lau - pi2
  thu =  asin(t3)

end subroutine vortex_rotatedsphereback

end module fvm_transformation_mod