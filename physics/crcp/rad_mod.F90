module rad_mod
  use grid_init_mod, only : height
  use prof_init_mod, only : th_e, tm_e
  real :: hmx, hmn
contains  
  subroutine radcool(ft,radav, nx, nz)
    implicit none
    integer, intent(in) :: nx, nz
    real, intent(inout) :: ft(nx,nz)
    real, intent(out)   :: radav(nz) 
!    
    integer :: i, k
    real :: coe
!    real, parameter :: hmx=15.e3, hmn=12.e3, day=24.*3600.
    real, parameter :: day=24.*3600.

    real :: fun23, a, b, c

    fun23(a,b,c)=max(0.,min(1.,(b-a)/(b-c)))

!       day=24.*3600.
!       hmx=15.e3
!       hmn=12.e3
       do k=1,nz
       coe=-1.5 * fun23(height(k),hmx,hmn) * th_e(k)/tm_e(k)
       do i=1,nx
       ft(i,k)=ft(i,k)+2.*coe/day
       enddo
       radav(k)=coe/day
       enddo

      return
    end subroutine radcool
  end module rad_mod
