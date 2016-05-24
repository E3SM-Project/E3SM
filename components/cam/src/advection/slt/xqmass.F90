
subroutine xqmass(cwava   ,etamid  ,w       ,qo      ,qn      , &
                  xo      ,xn      ,pdela   ,pdelb   ,hwxal   , &
                  hwxbl   ,nlon    )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute comtribution of current latitude to global integrals necessary
! to compute the fixer for the non-water constituents.
! 
! Method: 
! 
! Author: J. Olson, March 1994
! 
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid,       only: plev, plon
  use constituents, only: pcnst,  cnst_get_type_byind

  implicit none

!---------------------------Arguments-----------------------------------
  real(r8), intent(in)  :: cwava                ! normalization factor
  real(r8), intent(in)  :: etamid(plev)         ! vertical coords at midpoints 
  real(r8), intent(in)  :: w                    ! gaussian weight this latitude
  real(r8), intent(in)  :: qo(plon,plev      )  ! q old            (pre -SLT)
  real(r8), intent(in)  :: qn(plon,plev      )  ! q new            (post-SLT)
  real(r8), intent(in)  :: xo(plon,plev,pcnst)  ! old constituents (pre -SLT)
  real(r8), intent(in)  :: xn(plon,plev,pcnst)  ! new constituents (post-SLT)
  real(r8), intent(in)  :: pdela(plon,plev)     ! pressure diff between interfaces
  integer , intent(in) :: nlon                  ! number of longitudes
                                                ! based pure pressure part of hybrid grid
  real(r8), intent(in)  :: pdelb(plon,plev)     ! pressure diff between interfaces
                                                ! based sigma part of hybrid grid
  real(r8), intent(inout) :: hwxal(pcnst,4)     ! partial integrals (weighted by pure
                                                ! pressure part of hybrid pressures)
  real(r8), intent(inout) :: hwxbl(pcnst,4)     ! partial integrals (weighted by sigma
                                                ! part of hybrid pressures)
!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  integer i                     ! longitude index
  integer k                     ! level index
  integer m                     ! constituent index
  integer n                     ! index for partial integral
  real(r8) a                    ! integral constant
  real(r8) xdx,xq1,xqdq,xdxq1   ! work elements
  real(r8) xdxqdq               ! work elements
  real(r8) hwak(4),hwbk(4)      ! work arrays
  real(r8) q1 (plon,plev)       ! work array
  real(r8) qdq(plon,plev)       ! work array
  real(r8) hwalat(4)            ! partial integrals (weighted by pure
!                               ! pressure part of hybrid pressures)
  real(r8) hwblat(4)            ! partial integrals (weighted by sigma
!                               ! part of hybrid pressures)
  real(r8) etamsq(plev)         ! etamid*etamid
  real(r8) xnt(plon)            ! temp version of xn
  character*3 cnst_type         ! 'dry' or 'wet' mixing ratio
!-----------------------------------------------------------------------
!
  a = cwava*w*0.5_r8
  do k = 1,plev
     etamsq(k) = etamid(k)*etamid(k)
  end do
!
! Compute terms involving water vapor mixing ratio
!
!$OMP PARALLEL DO PRIVATE (K, I)
  do k = 1,plev
     do i = 1,nlon
        q1 (i,k) = 1._r8 - qn(i,k)
        qdq(i,k) = qn(i,k)*abs(qn(i,k) - qo(i,k))
     end do
  end do
!
! Compute partial integrals for non-water constituents
!
!$OMP PARALLEL DO PRIVATE (M, CNST_TYPE, N, HWALAT, HWBLAT, K, HWAK, HWBK, &
!$OMP                      I, XNT, XDX, XQ1, XQDQ, XDXQ1, XDXQDQ)
  do m = 2,pcnst
     cnst_type = cnst_get_type_byind(m)
     do n = 1,4
        hwalat(n) = 0._r8
        hwblat(n) = 0._r8
     end do
     do k = 1,plev
        do n = 1,4
           hwak(n) = 0._r8
           hwbk(n) = 0._r8
        end do

        if (cnst_type.eq.'dry' ) then
           do i = 1, nlon
              if (abs(xn(i,k,m) - xo(i,k,m)) &
                 .lt.1.0e-13_r8 * max(abs(xn(i,k,m)), abs(xo(i,k,m)))) then
                 xnt(i) = xo(i,k,m)
              else
                 xnt(i) = xn(i,k,m)
              end if
           end do
        else
           do i = 1, nlon
              xnt(i) = xn(i,k,m)
           end do
        end if

        do i = 1,nlon
           xdx    = xnt(i)*abs(xn(i,k,m) - xo(i,k,m))
           xq1    = xnt(i)*q1 (i,k)
           xqdq   = xnt(i)*qdq(i,k)
           xdxq1  = xdx      *q1 (i,k)
           xdxqdq = xdx      *qdq(i,k)

           hwak(1) = hwak(1) + xq1   *pdela(i,k)
           hwbk(1) = hwbk(1) + xq1   *pdelb(i,k)
           hwak(2) = hwak(2) + xqdq  *pdela(i,k)
           hwbk(2) = hwbk(2) + xqdq  *pdelb(i,k)
           hwak(3) = hwak(3) + xdxq1 *pdela(i,k)
           hwbk(3) = hwbk(3) + xdxq1 *pdelb(i,k)
           hwak(4) = hwak(4) + xdxqdq*pdela(i,k)
           hwbk(4) = hwbk(4) + xdxqdq*pdelb(i,k)
        end do

        hwalat(1) = hwalat(1) + hwak(1)
        hwblat(1) = hwblat(1) + hwbk(1)
        hwalat(2) = hwalat(2) + hwak(2)*etamid(k)
        hwblat(2) = hwblat(2) + hwbk(2)*etamid(k)
        hwalat(3) = hwalat(3) + hwak(3)*etamid(k)
        hwblat(3) = hwblat(3) + hwbk(3)*etamid(k)
        hwalat(4) = hwalat(4) + hwak(4)*etamsq(k)
        hwblat(4) = hwblat(4) + hwbk(4)*etamsq(k)
     end do
!
! The 0.5 factor arises because Gaussian weights sum to 2
!
     do n = 1,4
        hwxal(m,n) = hwxal(m,n) + hwalat(n)*a
        hwxbl(m,n) = hwxbl(m,n) + hwblat(n)*a
     end do
  end do

  return
end subroutine xqmass
