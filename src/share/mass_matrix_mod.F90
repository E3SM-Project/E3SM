#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module mass_matrix_mod
  use kinds, only : real_kind
  use dimensions_mod, only : np, nelemd
  use quadrature_mod, only : quadrature_t, gauss ,gausslobatto
  use element_mod, only : element_t
  use parallel_mod, only : parallel_t
  use edge_mod, only : edgebuffer_t,edgevpack,edgevunpack, &
       freeedgebuffer,initedgebuffer  
  use bndry_mod, only : bndry_exchangev
implicit none
private

  public :: mass_matrix

contains

! ===========================================
! mass_matrix:
!
! Compute the mass matrix for each element...
! ===========================================

  subroutine mass_matrix(par,elem)

    type (parallel_t),intent(in) :: par
    type (element_t) :: elem(:)

    type (EdgeBuffer_t)    :: edge

    real(kind=real_kind)  da                     ! area element

    type (quadrature_t) :: gp

    integer ii
    integer i,j
    integer kptr
    integer iptr

    ! ===================
    ! begin code
    ! ===================

    call initEdgeBuffer(edge,1)

    ! =================================================
    ! mass matrix on the velocity grid
    ! =================================================    

    gp=gausslobatto(np)
 
    do ii=1,nelemd
       do j=1,np
          do i=1,np
              ! MNL: metric term for map to reference element is now in metdet!
             elem(ii)%mp(i,j)=gp%weights(i)*gp%weights(j)
             elem(ii)%rmp(i,j)=elem(ii)%mp(i,j)
          end do
       end do

       kptr=0
       call edgeVpack(edge,elem(ii)%rmp,1,kptr,elem(ii)%desc)

    end do

    ! ==============================
    ! Insert boundary exchange here
    ! ==============================

    call bndry_exchangeV(par,edge)

    do ii=1,nelemd

       kptr=0
       call edgeVunpack(edge,elem(ii)%rmp,1,kptr,elem(ii)%desc)

       do j=1,np
          do i=1,np
             elem(ii)%rmp(i,j)=1.0D0/elem(ii)%rmp(i,j)
          end do
       end do

    end do
#if (! defined VERT_OPENMP)
!$OMP BARRIER
#endif

    deallocate(gp%points)
    deallocate(gp%weights)

    ! =============================================
    ! compute spherical element mass matrix
    ! =============================================
    do ii=1,nelemd
       do j=1,np
          do i=1,np
             elem(ii)%spheremp(i,j)=elem(ii)%mp(i,j)*elem(ii)%metdet(i,j)
             elem(ii)%rspheremp(i,j)=elem(ii)%spheremp(i,j)
          end do
       end do
       kptr=0
       call edgeVpack(edge,elem(ii)%rspheremp,1,kptr,elem(ii)%desc)
    end do
    call bndry_exchangeV(par,edge)
    do ii=1,nelemd
       kptr=0
       call edgeVunpack(edge,elem(ii)%rspheremp,1,kptr,elem(ii)%desc)
       do j=1,np
          do i=1,np
             elem(ii)%rspheremp(i,j)=1.0D0/elem(ii)%rspheremp(i,j)
          end do
       end do
    end do
#if (! defined VERT_OPENMP)
!$OMP BARRIER
#endif

    ! =============================================
    ! compute the mass matrix 
    ! =============================================
    ! Jose Garcia: Not sure but I think this code is just dead code
    !do ii=1,nelemd
    !   iptr=1
    !   do j=1,np
    !      do i=1,np
    !         elem(ii)%mp(i,j)=elem(ii)%mp(i,j)
    !         iptr=iptr+1
    !      end do
    !   end do
    !end do
   
    call FreeEdgeBuffer(edge)
       
  end subroutine mass_matrix

end module mass_matrix_mod

