  subroutine cldefrint(m,n,landfrac,tlayer,rel,rei,psurface, player,landm,icefrac,snowh)
!-----------------------------------------------------------------------
!
! interface for cldefr to work with isccp simulator calls and CAM3 radiation
!
!-----------------------------------------------------------------------
    use rrtm_grid, only: nzm
    use cam_rad_parameterizations, only : &
        computeRe_Liquid, computeRe_Ice
    implicit none
!------------------------------Parameters-------------------------------

! Input arguments
!
  integer m,n, k
  real landfrac(1)
  real icefrac(1)       ! Ice fraction
  real psurface(1)      ! Surface pressure
  real tlayer(1,nzm)   ! Temperature
  real player(1,nzm)   !  Midpoint pressures
  real landm(1)         ! Land fraction
  real snowh(1)      ! snow depth, water equivalent (meters)

!
! Output arguments
!
  real rel(nzm)   ! Liquid effective drop size (microns)
  real rei(nzm)   ! Ice effective drop size (microns)

  do k = 1,nzm
    rel(k) = computeRe_Liquid(tlayer(1,k), landm(1), icefrac(1), snowh(1))
    rei(k) = computeRe_Ice(tlayer(1,k))
  end do

  return
  end
 
