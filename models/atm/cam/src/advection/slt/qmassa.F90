module qmassa


contains

subroutine qmassarun(cwava   ,w       ,q3      ,pdel    ,hw1lat  , &
                  nlon    ,q0  ,lat     ,pdeld   )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Calculate contribution of current latitude to mass of constituents
! being advected by slt.
! 
! Method: 
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
  use pmgrid,       only: plev, plon
  use constituents, only: pcnst, cnst_get_type_byind
  use dycore, only: dycore_is
  use abortutils, only: endrun

  implicit none

!------------------------------Arguments--------------------------------
  integer , intent(in)  :: nlon                 ! longitude dimension  
  real(r8), intent(in)  :: cwava                ! normalization factor    l/(g*plon)
  real(r8), intent(in)  :: w                    ! gaussian weight this latitude
  real(r8), intent(in)  :: q3(plon,plev,pcnst)  ! constituents
  real(r8), intent(in)  :: q0(plon,plev,pcnst)  ! constituents at begining of time step
  real(r8), intent(in)  :: pdel(plon,plev)      ! pressure diff between interfaces
  real(r8), intent(out) :: hw1lat(pcnst)        ! accumulator
  real(r8), intent(in),optional  :: pdeld(:,:)  ! dry pressure difference for dry-type constituents
                                                ! only used when called from eularian dynamics

  
  integer lat
!-----------------------------------------------------------------------
!
!---------------------------Local variables-----------------------------
  integer i,k,m             ! longitude, level, constituent indices
  real(r8) const            ! temporary constant
!-----------------------------------------------------------------------
!
! Integration factor (the 0.5 factor arises because gaussian weights sum to 2)
!
  const = cwava*w*0.5_r8
  do m=1,pcnst
     hw1lat(m) = 0._r8
  end do

!$OMP PARALLEL DO PRIVATE (M, K, I)
  do m=1,pcnst
     if (m == 1) then
!
! Compute mass integral for water
!
        do k=1,plev
           do i=1,nlon
              hw1lat(1) = hw1lat(1) + q3(i,k,1)*pdel(i,k)
           end do
        end do
!
! Compute mass integral for non-water constituents (on either WET or DRY basis)
!
     elseif (cnst_get_type_byind(m).eq.'dry' ) then  ! dry type constituents
        if (  dycore_is ('EUL') ) then           ! EUL dycore computes pdeld in time filter
           if ( .not. present(pdeld) ) &
                call endrun('for dry type cnst with eul dycore, qmassa requires pdeld argument')
           do k=1,plev
              do i=1,nlon
                 hw1lat(m) = hw1lat(m) + q3(i,k,m)*pdeld(i,k)
              end do
           end do
        else !dycore SLD
           do k=1,plev
              do i=1,nlon
                 hw1lat(m) = hw1lat(m) + q3(i,k,m)*(1._r8 - q0(i,k,1))*pdel(i,k)
              end do
           end do
        endif ! dycore
     else  !wet type constituents
        do k=1,plev
           do i=1,nlon
              hw1lat(m) = hw1lat(m) + q3(i,k,m)*(1._r8 - q3(i,k,1))*pdel(i,k)
           end do
        end do
     end if !dry or wet
  end do  

  do m = 1,pcnst
     hw1lat(m) = hw1lat(m)*const
  end do

  return
end subroutine qmassarun

end module qmassa




