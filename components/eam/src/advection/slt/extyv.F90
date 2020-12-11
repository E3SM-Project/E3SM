
subroutine extyv(pkcnst  ,pkdim   ,coslam  ,sinlam  ,ub      , &
                 vb      )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Fill latitude extensions of a vector component extended array.
! 
! Method: 
! This is done in 2 steps:
! 1) interpolate to the pole points;
!    use coefficients for zonal wave number 1 on the Gaussian
!    latitude closest to the pole.
! 2) add latitude lines beyond the poles.
! 
! Author: J. Olson
! 
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid,       only: plat
  use scanslt,      only: nxpt, platd, nlonex, beglatex, endlatex, plond, &
                          jintmx
  implicit none

!------------------------------Parameters-------------------------------
  integer, parameter :: istart = nxpt+1           ! index to start computation
  integer, parameter :: js = 1    + nxpt + jintmx ! index of southernmost model lat
  integer, parameter :: jn = plat + nxpt + jintmx ! index of northernmost model lat
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
!
! Input arguments
!
  integer , intent(in) :: pkcnst              ! dimensioning construct for 3-D arrays
  integer , intent(in) :: pkdim               ! vertical dimension
  real(r8), intent(in) :: coslam(plond,platd) ! Cos of long at x-grid points (global grid)
  real(r8), intent(in) :: sinlam(plond,platd) ! Sin of long at x-grid points (global grid)
  real(r8), intent(inout):: ub(plond,pkdim*pkcnst,beglatex:endlatex) ! U-wind with extents
  real(r8), intent(inout):: vb(plond,pkdim*pkcnst,beglatex:endlatex) ! V-wind with extents
!
!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  integer i                 ! index
  integer ig                ! index
  integer j                 ! index
  integer k                 ! index
  integer nlon2             ! half the number of real longitudes
  integer istop             ! index to stop  computation
  real(r8) zavecv           ! accumulator for wavenumber 1 of v
  real(r8) zavesv           ! accumulator for wavenumber 1 of v
  real(r8) zavecu           ! accumulator for wavenumber 1 of u
  real(r8) zavesu           ! accumulator for wavenumber 1 of u
  real(r8) zaucvs           ! used to couple u and v (wavenumber 1)
  real(r8) zavcus           ! used to couple u and v (wavenumber 1)
!-----------------------------------------------------------------------
!
! Fill north pole line.
!
#if ( defined SPMD )
  if (jn+1<=endlatex) then   ! north pole is on-processor
#endif
     do k = 1,pkdim
        zavecv = 0.0_r8
        zavesv = 0.0_r8
        zavecu = 0.0_r8
        zavesu = 0.0_r8
        ig     = 0
        istop  = nxpt + nlonex(jn)
        do i = istart,istop
           ig     = ig + 1
           zavecv = zavecv + vb(i,k,jn  )*coslam(ig,jn)
           zavesv = zavesv + vb(i,k,jn  )*sinlam(ig,jn)
           zavecu = zavecu + ub(i,k,jn  )*coslam(ig,jn)
           zavesu = zavesu + ub(i,k,jn  )*sinlam(ig,jn)
        end do
        zavcus = (zavecv + zavesu)/nlonex(jn)
        zaucvs = (zavecu - zavesv)/nlonex(jn)
        ig     = 0
        istop  = nxpt + nlonex(jn+1)
        do i = istart,istop
           ig           = ig + 1
           vb(i,k,jn+1) = zavcus*coslam(ig,jn+1) - zaucvs*sinlam(ig,jn+1)
           ub(i,k,jn+1) = zaucvs*coslam(ig,jn+1) + zavcus*sinlam(ig,jn+1)
        end do
     end do
#if ( defined SPMD )
  end if
#endif
!
! Fill northern lines beyond pole line.
!
  if( jn+2 <= platd )then
     do j = jn+2,platd
#if ( defined SPMD )
        if (j<=endlatex) then
#endif
           nlon2 = nlonex(j)/2
           do k = 1,pkdim
              do i = istart,istart+nlon2-1
                 vb(      i,k,j) = -vb(nlon2+i,k,2*jn+2-j)
                 vb(nlon2+i,k,j) = -vb(      i,k,2*jn+2-j)
                 ub(      i,k,j) = -ub(nlon2+i,k,2*jn+2-j)
                 ub(nlon2+i,k,j) = -ub(      i,k,2*jn+2-j)
              end do
           end do
#if ( defined SPMD )
        end if
#endif
     end do
  end if
! 
! Fill south pole line.
!
#if ( defined SPMD )
  if (js-1>=beglatex) then   ! south pole is on-processor
#endif
     do k = 1,pkdim
        zavecv = 0.0_r8
        zavesv = 0.0_r8
        zavecu = 0.0_r8
        zavesu = 0.0_r8
        ig     = 0
        istop  = nxpt + nlonex(js)
        do i = istart,istop
           ig     = ig + 1
           zavecv = zavecv + vb(i,k,js  )*coslam(ig,js)
           zavesv = zavesv + vb(i,k,js  )*sinlam(ig,js)
           zavecu = zavecu + ub(i,k,js  )*coslam(ig,js)
           zavesu = zavesu + ub(i,k,js  )*sinlam(ig,js)
        end do
        zavcus = (zavecv - zavesu)/nlonex(js)
        zaucvs = (zavecu + zavesv)/nlonex(js)
        ig     = 0
        istop  = nxpt + nlonex(js-1)
        do i = istart,istop
           ig           = ig + 1
           vb(i,k,js-1) = zavcus*coslam(ig,js-1) + zaucvs*sinlam(ig,js-1)
           ub(i,k,js-1) = zaucvs*coslam(ig,js-1) - zavcus*sinlam(ig,js-1)
        end do
     end do
#if ( defined SPMD )
  end if
#endif
!
! Fill southern lines beyond pole line.
!
  if( js-2 >= 1 )then
     do j = 1,js-2
#if ( defined SPMD )
        if (j>=beglatex) then
#endif
           nlon2 = nlonex(j)/2
           do k = 1,pkdim
              do i = istart,istart+nlon2-1
                 vb(      i,k,j) = -vb(nlon2+i,k,2*js-2-j)
                 vb(nlon2+i,k,j) = -vb(      i,k,2*js-2-j)
                 ub(      i,k,j) = -ub(nlon2+i,k,2*js-2-j)
                 ub(nlon2+i,k,j) = -ub(      i,k,2*js-2-j)
              end do
           end do
#if ( defined SPMD )
        end if
#endif
     end do
  end if
!
  return
end subroutine extyv
