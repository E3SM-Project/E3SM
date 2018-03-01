#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
! hevi mod contains subroutines used by the implicit solver and some element and statesaved
! loop operations used in the HEVI-IMEX methods
module hevi_mod

  use dimensions_mod, only: np, nlev
  use element_mod,    only: element_t
  use kinds,          only: real_kind
  implicit none

contains

! state_save is a subroutine that saves the value of the state variables at a 
! give time-step to a storage vector
  subroutine state_save(elem,state,n,nets,nete)

    type (element_t),                   intent(inout), target :: elem(:)!
    real (kind=real_kind), intent(inout) :: state(nets:nete,np,np,nlev,6)
    integer                             :: n,nets,nete

    integer :: ie
    do ie=nets,nete
      state(ie,:,:,:,1) = elem(ie)%state%v(:,:,1,:,n)
      state(ie,:,:,:,2) = elem(ie)%state%v(:,:,2,:,n)
      state(ie,:,:,:,3) = elem(ie)%state%w(:,:,:,n)
      state(ie,:,:,:,4) = elem(ie)%state%phinh(:,:,:,n)
      state(ie,:,:,:,5) = elem(ie)%state%theta_dp_cp(:,:,:,n)
      state(ie,:,:,:,6) = elem(ie)%state%dp3d(:,:,:,n)
    end do

  end subroutine state_save

! state_read is a subroutine that reads from a storage vector to 
! overwrite the value of the state variables at a given time-step
  subroutine state_read(elem,state,n,nets,nete)

    type (element_t),                   intent(inout), target :: elem(:)!
    real (kind=real_kind), intent(inout) :: state(nets:nete,np,np,nlev,6)
    integer                            :: n,nets,nete

    integer :: ie
    do ie=nets,nete
      elem(ie)%state%v(:,:,1,:,n)         = state(ie,:,:,:,1)
      elem(ie)%state%v(:,:,2,:,n)         = state(ie,:,:,:,2)
      elem(ie)%state%w(:,:,:,n)           = state(ie,:,:,:,3)
      elem(ie)%state%phinh(:,:,:,n)       = state(ie,:,:,:,4)
      elem(ie)%state%theta_dp_cp(:,:,:,n) = state(ie,:,:,:,5)
      elem(ie)%state%dp3d(:,:,:,n)        = state(ie,:,:,:,6)
    end do

  end subroutine state_read

! Subroutine mgs implements the modified (stable) Gram-Schmidt algorithm
! and forms a QR factorization of a matrix A = Q*R where the QR factorization
! is over the vertical components 
  subroutine mgs(A,Q,R)

    real (kind=real_kind), intent(inout) :: A(np,np,nlev,nlev)
    real (kind=real_kind), intent(inout) :: Q(np,np,nlev,nlev)
    real (kind=real_kind), intent(inout) :: R(np,np,nlev,nlev)

    ! local variables
    integer :: i,j,k,l
    real (kind=real_kind) :: Atemp(nlev,nlev)
    real (kind=real_kind) :: Qtemp(nlev,nlev)
    real (kind=real_kind) :: Rtemp(nlev,nlev)

    R=0.0
    do i=1,np
      do j=1,np
        Atemp(:,:)=A(i,j,:,:)
        Rtemp=0.0
#if (defined COLUMN_OPENMP)
  !$omp parallel do default(shared), private(k)
#endif
        do k=1,nlev
          Rtemp(k,k)=sqrt(dot_product(Atemp(:,k),Atemp(:,k))) !norm2(Atemp(:,k))
          Qtemp(1:nlev,k)=Atemp(1:nlev,k)/Rtemp(k,k)
          do l=k+1,nlev
            Rtemp(k,l)=dot_product(Qtemp(:,k),Atemp(:,l))
            Atemp(:,l)=Atemp(:,l)-Rtemp(k,l)*Qtemp(:,k)
          end do
        end do
        A(i,j,:,:)=Atemp(:,:)
        Q(i,j,:,:)=Qtemp(:,:)
        R(i,j,:,:)=Rtemp(:,:)
      end do
    end do
  end subroutine mgs

! Backsubsitution computes the solution of R x = b where R is 
! assumed upper triangular and the solution is done over the 
! vertical hevi part of the arrays
  subroutine backsubstitution(R,b,x)

    real (kind=real_kind), intent(in) :: R(np,np,nlev,nlev)
    real (kind=real_kind), intent(in) :: b(np,np,nlev)
    real (kind=real_kind), intent(inout) :: x(np,np,nlev)

    integer :: i,j
    real (kind=real_kind) :: sum(np,np)
    real (kind=real_kind) :: Rtemp(nlev,nlev)
    real (kind=real_kind) :: btemp(nlev,1)
    real (kind=real_kind) :: xtemp(nlev,1)


    do i=nlev,1,-1
      sum(:,:) = b(:,:,i)
      do j=i+1,nlev
        sum(:,:)=sum(:,:)-R(:,:,i,j)*x(:,:,j)
      end do
      x(:,:,i)=sum(:,:)/R(:,:,i,i)
    end do
    do i=1,np
      do j=1,np
        Rtemp(:,:)=R(i,j,:,:)
        xtemp(:,1)=x(i,j,:)
        btemp(:,1)=b(i,j,:)
     end do
   end do

  end subroutine backsubstitution

! elemstate_add is an assortment of linear combinations of state varables
! and storage vectors at specified time-levels
!
! alpha==1, elem(n1) = a1*elem(n2) + a2*elem(n3)
! alpha==2, elem(n1) = a1*elem(n2) + a3*state
! alpha==3, elem(n1) = a1*elem(n2) + a2*elem(n3) + a3*state
! alpha==4, elem(n1) = elem(n1) + a1*elem(n2)+a2*elem(n3)
 subroutine elemstate_add(elem,state,nets,nete,alpha,n1,n2,n3,a1,a2,a3)

    type (element_t)     , intent(inout), target :: elem(:)
    real (kind=real_kind), intent(in)  :: state(nets:nete,np,np,nlev,6)
    real (kind=real_kind), intent(in)  :: a1,a2,a3
    integer                            :: n1,n2,n3,nets,nete,alpha

    ! local variable
    integer :: ie

    if (alpha==1) then ! add two elements
      do ie=nets,nete
        elem(ie)%state%v(:,:,1,:,n1)         = a1*elem(ie)%state%v(:,:,1,:,n2) + &
          a2*elem(ie)%state%v(:,:,1,:,n3)
        elem(ie)%state%v(:,:,2,:,n1)         = a1*elem(ie)%state%v(:,:,2,:,n2) + &
          a2*elem(ie)%state%v(:,:,2,:,n3)
        elem(ie)%state%w(:,:,:,n1)           = a1*elem(ie)%state%w(:,:,:,n2) + &
          a2*elem(ie)%state%w(:,:,:,n3)
        elem(ie)%state%phinh(:,:,:,n1)         = a1*elem(ie)%state%phinh(:,:,:,n2) + &
          a2*elem(ie)%state%phinh(:,:,:,n3)
        elem(ie)%state%theta_dp_cp(:,:,:,n1) = a1*elem(ie)%state%theta_dp_cp(:,:,:,n2) + &
          a2*elem(ie)%state%theta_dp_cp(:,:,:,n3)
        elem(ie)%state%dp3d(:,:,:,n1)        = a1*elem(ie)%state%dp3d(:,:,:,n2) + &
          a2*elem(ie)%state%dp3d(:,:,:,n3)
      end do
    elseif (alpha==2) then ! add element with state
     do ie=nets,nete
       elem(ie)%state%v(:,:,1,:,n1)         = a3*state(ie,:,:,:,1) + &
         a1*elem(ie)%state%v(:,:,1,:,n2)
       elem(ie)%state%v(:,:,2,:,n1)         = a3*state(ie,:,:,:,2) + &
         a1*elem(ie)%state%v(:,:,2,:,n2)
       elem(ie)%state%w(:,:,:,n1)           = a3*state(ie,:,:,:,3) + &
         a1*elem(ie)%state%w(:,:,:,n2)
       elem(ie)%state%phinh(:,:,:,n1)         = a3*state(ie,:,:,:,4) + &
         a1*elem(ie)%state%phinh(:,:,:,n2)
       elem(ie)%state%theta_dp_cp(:,:,:,n1) = a3*state(ie,:,:,:,5) + &
         a1*elem(ie)%state%theta_dp_cp(:,:,:,n2)
       elem(ie)%state%dp3d(:,:,:,n1)        = a3*state(ie,:,:,:,6) + &
         a1*elem(ie)%state%dp3d(:,:,:,n2)
     end do
    elseif (alpha==3) then ! add 2 elements with state
     do ie=nets,nete
        elem(ie)%state%v(:,:,1,:,n1)         = a1*elem(ie)%state%v(:,:,1,:,n2) + &
          a2*elem(ie)%state%v(:,:,1,:,n3) + a3*state(ie,:,:,:,1)
        elem(ie)%state%v(:,:,2,:,n1)         = a1*elem(ie)%state%v(:,:,2,:,n2) + &
          a2*elem(ie)%state%v(:,:,2,:,n3) + a3*state(ie,:,:,:,2)
        elem(ie)%state%w(:,:,:,n1)           = a1*elem(ie)%state%w(:,:,:,n2) + &
          a2*elem(ie)%state%w(:,:,:,n3) + a3*state(ie,:,:,:,3)
        elem(ie)%state%phinh(:,:,:,n1)         = a1*elem(ie)%state%phinh(:,:,:,n2) + &
          a2*elem(ie)%state%phinh(:,:,:,n3) + a3*state(ie,:,:,:,4)
        elem(ie)%state%theta_dp_cp(:,:,:,n1) = a1*elem(ie)%state%theta_dp_cp(:,:,:,n2) + &
          a2*elem(ie)%state%theta_dp_cp(:,:,:,n3) + a3*state(ie,:,:,:,5)
        elem(ie)%state%dp3d(:,:,:,n1)        = a1*elem(ie)%state%dp3d(:,:,:,n2) + &
          a2*elem(ie)%state%dp3d(:,:,:,n3) + a3*state(ie,:,:,:,6)
      end do
    else ! if alpha ==4 then += an element with 2 other elements
      do ie=nets,nete
        elem(ie)%state%v(:,:,1,:,n1)         = a1*elem(ie)%state%v(:,:,1,:,n2) + &
          a2*elem(ie)%state%v(:,:,1,:,n3) + elem(ie)%state%v(:,:,1,:,n1)
        elem(ie)%state%v(:,:,2,:,n1)         = a1*elem(ie)%state%v(:,:,2,:,n2) + &
          a2*elem(ie)%state%v(:,:,2,:,n3) + elem(ie)%state%v(:,:,2,:,n1)
        elem(ie)%state%w(:,:,:,n1)           = a1*elem(ie)%state%w(:,:,:,n2) + &
          a2*elem(ie)%state%w(:,:,:,n3) + elem(ie)%state%w(:,:,:,n1)
        elem(ie)%state%phinh(:,:,:,n1)         = a1*elem(ie)%state%phinh(:,:,:,n2) + &
          a2*elem(ie)%state%phinh(:,:,:,n3) + elem(ie)%state%phinh(:,:,:,n1)
        elem(ie)%state%theta_dp_cp(:,:,:,n1) = a1*elem(ie)%state%theta_dp_cp(:,:,:,n2) + &
          a2*elem(ie)%state%theta_dp_cp(:,:,:,n3) + elem(ie)%state%theta_dp_cp(:,:,:,n1)
        elem(ie)%state%dp3d(:,:,:,n1)        = a1*elem(ie)%state%dp3d(:,:,:,n2) + &
          a2*elem(ie)%state%dp3d(:,:,:,n3) + elem(ie)%state%dp3d(:,:,:,n1)
      end do
    end if
  end subroutine elemstate_add


end module

