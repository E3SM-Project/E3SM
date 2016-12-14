#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module forcing_mod
  use kinds,          only : real_kind

  implicit none

  private

  integer,public,parameter :: INTERP_BDF2=1,EXP_EULER=2,EXP_EULER_NP1=3

  ! Used for the BDF scheme
  real (kind=real_kind),parameter :: an   = 4.0D0/3.0D0
  real (kind=real_kind),parameter :: anm1 =-1.0D0/3.0D0
  real (kind=real_kind),parameter :: bn   = 4.0D0/3.0D0
  real (kind=real_kind),parameter :: bnm1 =-2.0D0/3.0D0

  public :: Apply_Forcing, Zero_Forcing

contains

  !
  ! Called after timestep, before Robert filter
  ! 
  ! EXP_EULER:  uses F^n-1
  ! EXTRAPOLATED BDF-2:  uses F^n-1 and F^n
  subroutine Apply_Forcing(METHOD,elemin,hvcoord,tld,dt,tlp)
    use hybvcoord_mod, only : hvcoord_t
    use element_mod,    only : element_t
    use dimensions_mod, only : np, nlev, qsize
    use time_mod,       only : TimeLevel_t

    !input
    type (hvcoord_t), intent(inout)      :: hvcoord         ! hybrid vertical coordinate struct
    type (element_t), intent(inout)  :: elemin
    integer,intent(in)            :: METHOD
    type (TimeLevel_t),intent(in) :: tld
    type (TimeLevel_t),intent(in),optional :: tlp
    real (kind=real_kind)         :: dt

    !local
    real (kind=real_kind) :: dt2,dp, dt2_q
    integer :: i,j,k,ie, m,q
    integer :: np1_d,n0_d,nm1_d
    integer :: np1_p,n0_p,nm1_p

    np1_d = tld%np1
    n0_d  = tld%n0
    nm1_d = tld%nm1
    if(present(tlp)) then
       np1_p = tlp%np1	
       n0_p  = tlp%n0
       nm1_p = tlp%nm1
    else
       np1_p = tld%np1	
       n0_p  = tld%n0
       nm1_p = tld%nm1
    end if
    dt2   = 2_real_kind*dt
    dt2_q = dt2

    select case(METHOD)
    case(INTERP_BDF2)

       do k=1,nlev
          do j=1,np
             do i=1,np
                elemin%state%v(i,j,1,k,np1_d) = &
                     an*elemin%state%v(i,j,1,k,n0_d) + anm1*elemin%state%v(i,j,1,k,nm1_d) + &
                     bn*dt*elemin%derived%FM(i,j,1,k,n0_p) + bnm1*dt*elemin%derived%FM(i,j,1,k,nm1_p)
                elemin%state%v(i,j,2,k,np1_d) = &
                     an*elemin%state%v(i,j,2,k,n0_d) + anm1*elemin%state%v(i,j,2,k,nm1_d) + &
                     bn*dt*elemin%derived%FM(i,j,2,k,n0_p) + bnm1*dt*elemin%derived%FM(i,j,2,k,nm1_p)
!THETA
!                elemin%state%T(i,j,k,np1_d) = &
!                     an*elemin%state%T(i,j,k,n0_d) + anm1*elemin%state%T(i,j,k,nm1_d) + &
!                     bn*dt*elemin%derived%FT(i,j,k,n0_p) + bnm1*dt*elemin%derived%FT(i,j,k,nm1_p)                 
             enddo
          enddo
       enddo
       do m=1,qsize
          do k=1,nlev
             do j=1,np
                do i=1,np
                   elemin%state%Q(i,j,k,m) = elemin%state%Q(i,j,k,m) + dt2_q*elemin%derived%FQ(i,j,k,m,nm1_p)
                end do
             end do
          end do
       end do

    case(EXP_EULER)
       do k=1,nlev
          do j=1,np
             do i=1,np
                elemin%state%v(i,j,1,k,np1_d) = elemin%state%v(i,j,1,k,np1_d) + dt2*elemin%derived%FM(i,j,1,k,nm1_p)
                elemin%state%v(i,j,2,k,np1_d) = elemin%state%v(i,j,2,k,np1_d) + dt2*elemin%derived%FM(i,j,2,k,nm1_p)
!THETA
!                elemin%state%T(i,j,k,np1_d)   = elemin%state%T(i,j,k,np1_d)   + dt2*elemin%derived%FT(i,j,k,nm1_p)
             enddo
          enddo
       enddo
       do m=1,qsize
          elemin%state%Q(:,:,:,m) = elemin%state%Q(:,:,:,m)&
               + dt2_q*elemin%derived%FQ(:,:,:,m,nm1_p) 
       end do

    case default
       ! Selects no forcing       
    end select

    ! forcing always applied to Q, update Qdp:
    do q=1,qsize
       do k=1,nlev
          do j=1,np
             do i=1,np
                dp = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                     ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elemin%state%ps_v(i,j,np1_d)
                elemin%state%Qdp(i,j,k,q,np1_d) = elemin%state%Q(i,j,k,q)*dp
             enddo
          enddo
       enddo
    enddo
  end subroutine Apply_Forcing


  subroutine Zero_Forcing(elemin,tl)
    use element_mod, only : element_t
    use time_mod, only : timelevel_t
    use dimensions_mod, only : np, nlev, qsize
    type(element_t), intent(inout) :: elemin
    type(timelevel_t), intent(in) :: tl

    integer :: i, j, k, nm1, m

    nm1=tl%nm1
    do k=1,nlev
       do j=1,np
          do i=1,np
             elemin%derived%FM(i,j,1,k,nm1) = 0.0_real_kind
             elemin%derived%FM(i,j,2,k,nm1) = 0.0_real_kind
             elemin%derived%FT(i,j,k,nm1) = 0.0_real_kind
          enddo
       enddo
    enddo
    do m=1,qsize
       do k=1,nlev
          do j=1,np
             do i=1,np
                elemin%derived%FQ(i,j,k,m,nm1) = 0.0_real_kind
             end do
          end do
       end do
    end do
    
   end subroutine Zero_Forcing


end module forcing_mod
