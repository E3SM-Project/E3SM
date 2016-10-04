#ifdef HAVE_CONFIG_H
#include "config.h"
#endif



module rk_mod

  ! =======================
  use dimensions_mod, only : nlev, np
  ! =======================
  use element_mod, only : element_t
  ! =======================
  use hybrid_mod, only : hybrid_t
  ! =======================
  use kinds, only : real_kind
  ! =======================
  use types_mod, only : rk_t
  ! =======================
  use hybrid_mod,  only : hybrid_t
  ! =======================
  use reduction_mod, only : ParallelMax
  ! =======================
  use physical_constants, only : dd_pi
  implicit none

  private  

  public :: RkInit, ComputeCFL  

contains

  subroutine ComputeCFL(dt_grv,pmean,mindx,n0,elem,nets,nete,hybrid,cfl,dt_adv)

    real (kind=real_kind),intent(in) :: pmean,mindx
    integer, intent(in)              :: n0
    type (element_t), intent(in)     :: elem(:)
    integer, intent(in)              :: nets
    integer, intent(in)              :: nete
    type (hybrid_t), intent(in)      :: hybrid   ! distributed parallel structure (shared)
    
    ! Result

    integer,               intent(out) :: cfl
    real (kind=real_kind), intent(out) :: dt_adv
    real (kind=real_kind)              :: dt_grv

    ! local

    real (kind=real_kind) :: max_vel,v1,v2,dt_grvin
    integer :: ie,i,j,k

    max_vel = 30.0D0 ! lower bound
    
    do ie=nets,nete
       do k=1,nlev
          do j=1,np
             do i=1,np
                v1     = elem(ie)%D(i,j,1,1)*elem(ie)%state%v(i,j,1,k,n0)+elem(ie)%D(i,j,1,2)*elem(ie)%state%v(i,j,2,k,n0)
                v2     = elem(ie)%D(i,j,2,1)*elem(ie)%state%v(i,j,1,k,n0)+elem(ie)%D(i,j,2,2)*elem(ie)%state%v(i,j,2,k,n0)
                max_vel = max(max_vel, sqrt(v1*v1+v2*v2))
             enddo
          enddo
       enddo
    enddo

    max_vel=ParallelMax(max_vel,hybrid)
    
    dt_adv = (mindx/dd_pi)/max_vel

    do ie=nets,nete
       do k=1,nlev
          do j=1,np
             do i=1,np
                v1     = elem(ie)%D(i,j,1,1)*elem(ie)%state%v(i,j,1,k,n0)+elem(ie)%D(i,j,1,2)*elem(ie)%state%v(i,j,2,k,n0)
                v2     = elem(ie)%D(i,j,2,1)*elem(ie)%state%v(i,j,1,k,n0)+elem(ie)%D(i,j,2,2)*elem(ie)%state%v(i,j,2,k,n0)
                max_vel = max(max_vel, sqrt(v1*v1+v2*v2) + sqrt(pmean + elem(ie)%state%p(i,j,k,n0)))
             enddo
          enddo
       enddo
    enddo

    max_vel=ParallelMax(max_vel,hybrid)
    
    dt_grvin = (mindx/dd_pi)/max_vel

    
    cfl = floor(0.75D0*dt_adv/dt_grv)+1
    
    dt_adv = min(dt_adv,dt_grv*real(cfl))
    cfl=cfl+1

  end subroutine ComputeCFL


  subroutine RkInit_(cfl,MyRk)
    
    implicit none
    type (Rk_t)                     :: MyRk
    integer                         :: cfl,stages,i
    logical :: Debug =.FALSE.
     
    ! No Butcher tableau we use
    ! the Kraaijevanger form
    
    stages = 3

    MyRk%beta(1)   = 1.0D0
    MyRk%alpha(1)  = 1.0D0
    MyRk%alpha0(1) = 0.0D0

    MyRk%beta(2)   = 0.250D0
    MyRk%alpha(2)  = 0.250D0
    MyRk%alpha0(2) = 0.750D0

    MyRk%beta(3)   = 2.0D0/3.0D0
    MyRk%alpha(3)  = 2.0D0/3.0D0
    MyRk%alpha0(3) = 1.0D0/3.0D0

    MyRk%RKCFL  = real(cfl,kind=real_kind)

    MyRk%Stages = 3

#if 0
    print *,"RK STAGES = ", MyRk%Stages
    do i=1,stages
       print *,i," b,a,a0 = ",MyRk%beta(i),",",MyRk%alpha(i),",",MyRk%alpha0(i)
    end do
#endif

  end subroutine RkInit_


  subroutine RkInit(cfl,MyRk)
    
    implicit none
    type (Rk_t)                     :: MyRk
    integer                         :: cfl,stages,i
    logical :: Debug =.FALSE.
     
    ! No Butcher tableau we use
    ! the Kraaijevanger form
    
    stages = cfl+1

    do i=1,stages
       MyRk%beta(i)   = 1.0D0/real(stages-1,kind=real_kind)
       MyRk%alpha(i)  = 1.0D0
       MyRk%alpha0(i) = 0.0D0
    end do

    MyRk%beta(stages)   = 1.0D0/real(stages,kind=real_kind)
    MyRk%alpha(stages)  = real(stages-1,kind=real_kind)/real(stages,kind=real_kind)
    MyRk%alpha0(stages) = 1.0D0/real(stages,kind=real_kind)

    MyRk%RKCFL  = real(cfl,kind=real_kind)

    MyRk%Stages = Stages

#if 0
    print *,"RK STAGES = ", MyRk%Stages
    do i=1,stages
       print *,i," b,a,a0 = ",MyRk%beta(i),",",MyRk%alpha(i),",",MyRk%alpha0(i)
    end do
#endif

  end subroutine RkInit

end module rk_mod
