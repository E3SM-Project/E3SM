
subroutine kdpfnd(pkdim   ,pmap    ,sig     ,sigdp   ,kdpmap  , &
                  kdp     ,nlon    )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Determine vertical departure point indices that point into a grid
! containing the full or half sigma levels.  Use an artificial evenly 
! spaced vertical grid to map into the true model levels.
! 
! Method: 
! Indices are computed assuming the the sigdp values have
! been constrained so that sig(1) .le. sigdp(i,j) .lt. sig(pkdim).
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
  use pmgrid,       only: plon, plev
  implicit none

!------------------------------Arguments--------------------------------
  integer , intent(in) :: pkdim             ! dimension of "sig"
  integer , intent(in) :: pmap              ! dimension of "kdpmap"
  real(r8), intent(in) :: sig  (pkdim)      ! vertical grid coordinates
  integer , intent(in) :: kdpmap(pmap)      ! array of model grid indices which
  real(r8), intent(in) :: sigdp(plon,plev)  ! vertical coords. of departure points
  integer , intent(out):: kdp(plon,plev)    ! vertical index for each dep. pt.
  integer , intent(in) :: nlon              ! longitude dimensio
!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  integer i,k,ii            ! indices
  real(r8) rdel             ! recip. of interval in artificial grid
  real(r8) sig1ln           ! ln (sig(1))
!-----------------------------------------------------------------------
!
  rdel   = real(pmap,r8)/( log(sig(pkdim)) - log(sig(1)) )
  sig1ln = log( sig(1) )
!
!$OMP PARALLEL DO PRIVATE (K, I, II)
  do k=1,plev
     do i=1,nlon
!
! First guess of the departure point's location in the model grid
!
        ii = max0(1,min0(pmap,int((log(sigdp(i,k))-sig1ln)*rdel+1._r8)))
        kdp(i,k) = kdpmap(ii)
!
! Determine if location is in next interval
!
        if(sigdp(i,k) >= sig(kdp(i,k)+1)) then
           kdp(i,k) = kdp(i,k) + 1
        end if
     end do
  end do

  return
end subroutine kdpfnd
