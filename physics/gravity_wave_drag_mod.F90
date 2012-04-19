#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module gravity_wave_drag_mod
  use kinds,          only : real_kind
  use dimensions_mod, only : nlev, np
  use element_mod,    only : element_t
  use time_mod,       only : TimeLevel_t
  use control_mod,    only : integration
  implicit none
  
!
! Ue, Ve, Tme, The, Pre, Qve initial profiles.
! These are currently initialized in aquaplanet.
! They can be used also to other than G-W purposes.
! If this happen they should be moved from here
! to aquaplanet initialization module.
!
  real (kind=real_kind) :: Rayleigh(nlev)    ! G-W absorber 
  real (kind=real_kind) :: Ue(nlev)          ! G-W absorber   
  real (kind=real_kind) :: Ve(nlev)          ! G-W absorber   
  real (kind=real_kind) :: Tme(nlev)         ! G-W absorber   
  real (kind=real_kind) :: The(nlev)         ! G-W absorber   
  real (kind=real_kind) :: Pre(nlev)         ! G-W absorber   
  real (kind=real_kind) :: Qve(nlev)         ! G-W absorber   

  contains
  subroutine gravity_wave_drag_forcing(dt,elemin,tl)
    type(element_t), intent(inout) :: elemin
    type (TimeLevel_t)   :: tl
    real(kind=real_kind), intent(in) :: dt

    real(kind=real_kind) :: v1,v2,U,V,UGWdrag,VGWdrag
    integer :: i, j, k, nm1, m

    nm1=tl%nm1

    if(integration == "explicit")then
    do k=1,nlev
      do m=1,2
        do j=1,np
          do i=1,np

	    elemin%derived%FM(i,j,m,k,nm1) = elemin%derived%FM(i,j,m,k,nm1) - &
               elemin%state%v(i,j,m,k,nm1)*Rayleigh(k)

          enddo
        enddo
      enddo
    enddo
    else
    do k=1,nlev
      do j=1,np
        do i=1,np
          v1 = elemin%state%v(i,j,1,k,nm1)
          v2 = elemin%state%v(i,j,2,k,nm1)

          U  = v1*elemin%D(1,1,i,j) + v2*elemin%D(1,2,i,j)
          V  = v1*elemin%D(2,1,i,j) + v2*elemin%D(2,2,i,j)

          UGWdrag = U*2.0_real_kind*dt*Rayleigh(k)
          VGWdrag = V*2.0_real_kind*dt*Rayleigh(k)

          elemin%derived%FM(i,j,1,k,nm1) = &
            elemin%derived%FM(i,j,1,k,nm1) - &
            (UGWdrag*elemin%Dinv(1,1,i,j) + VGWdrag*elemin%Dinv(1,2,i,j))

          elemin%derived%FM(i,j,2,k,nm1) = &
            elemin%derived%FM(i,j,2,k,nm1) - &
            (UGWdrag*elemin%Dinv(2,1,i,j) + VGWdrag*elemin%Dinv(2,2,i,j))
        enddo
      enddo
    enddo
    endif
    
  end subroutine gravity_wave_drag_forcing
    	

end module gravity_wave_drag_mod

