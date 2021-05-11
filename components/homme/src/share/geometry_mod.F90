#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module geometry_mod
  use kinds, only : real_kind, long_kind, longdouble_kind
  use coordinate_systems_mod, only : sphere_tri_area
  use dimensions_mod, only : np
  use planar_mod, only : plane_init_atomic
  use cube_mod, only : cube_init_atomic
  use parallel_mod, only : abortmp

implicit none

public  :: set_area_correction_map0, set_area_correction_map2

contains

subroutine set_area_correction_map0(elem, nelemd, par, gp)
! Numerical area of the domain (sphere) is sum of integration weights. The sum
! is not exactly equal geometric area (4\pi R for cube). It is required that
! numerical area = geometric area.
! This correction 'butters'
! the difference between numerical and geometrical areas evenly among DOFs. Then
! geometric areas of individual elements still do not equal their numerical areas,
! only whole domain's areas coinside.


#ifndef CAM
  use repro_sum_mod,      only: repro_sum
#else
  use shr_reprosum_mod,   only: repro_sum => shr_reprosum_calc
#endif

  use quadrature_mod,   only: quadrature_t
  use element_mod, only: element_t
  use parallel_mod, only: parallel_t
  use kinds, only: iulog
  use control_mod, only: geometry
  use physical_constants, only: domain_size
  implicit none

  type (element_t),   pointer, intent(inout)     :: elem(:)
  type (parallel_t),  intent(in)  :: par
  integer, intent(in) :: nelemd
  type (quadrature_t), intent(in)   :: gp

  real(kind=real_kind) :: aratio(nelemd,1)
  real(kind=real_kind) :: area(1)
  integer :: ie, i, j

  area = 0.0d0
  do ie=1,nelemd
     aratio(ie,1) = sum(elem(ie)%mp(:,:)*elem(ie)%metdet(:,:))
  enddo

  call repro_sum(aratio, area, nelemd, nelemd, 1, commid=par%comm)

  area(1) = domain_size/area(1)  ! ratio correction


  if (par%masterproc) &
     write(iulog,'(a,f20.17)') " re-initializing elements: alpha area correction=",&
     area(1)

  if ( geometry == "sphere" ) then
     do ie=1,nelemd
        call cube_init_atomic(elem(ie),gp%points,area(1))
     enddo
  else if ( geometry == "plane" ) then
    do ie=1,nelemd
       call plane_init_atomic(elem(ie),gp%points,area(1))
    enddo
  end if

end subroutine set_area_correction_map0


subroutine set_area_correction_map2(elem, nelemd, par, gp)
! Numerical area of the domain (sphere) is sum of integration weights. The sum
! is not exactly equal geometric area (4\pi R for cube). It is required that
! numerical area = geometric area.
! The 'epsilon bubble' approach modifies inner weights in each element so that
! geometic and numerical areas of each element match.
#ifndef CAM
  use repro_sum_mod,      only: repro_sum
#else
  use shr_reprosum_mod,   only: repro_sum => shr_reprosum_calc
#endif

  use quadrature_mod,   only: quadrature_t
  use element_mod, only: element_t
  use parallel_mod, only: parallel_t
  use kinds, only: iulog
  use control_mod, only: geometry
  use physical_constants, only: domain_size, dx, dy

  implicit none

  type (element_t),   pointer, intent(inout)     :: elem(:)
  type (parallel_t),  intent(in)  :: par
  integer, intent(in) :: nelemd
  type (quadrature_t), intent(in)   :: gp

  real(kind=real_kind) :: aratio(nelemd,1)
  real(kind=real_kind) :: area(1), area_elem, area_num, area_dummy,sum_w, delta
  integer :: ie, i, j
  real(kind=real_kind) :: tol_local = 1e-15


  if ( geometry == "sphere" ) then
     do ie=1,nelemd
        ! Obtain area of element = sum of areas of 2 triangles.
        call sphere_tri_area(elem(ie)%corners3D(1), elem(ie)%corners3D(2),&
                             elem(ie)%corners3D(3), area_elem)
        call sphere_tri_area(elem(ie)%corners3D(1), elem(ie)%corners3D(3),&
                             elem(ie)%corners3D(4), area_dummy)
        ! Store element's area in area_elem.
        area_elem = area_elem + area_dummy

        ! Compute 'numerical area' of the element as sum of integration
        ! weights.
        area_num = 0.0d0
        do j = 1,np
           do i = 1,np
              area_num = area_num + gp%weights(i)*gp%weights(j)*elem(ie)%metdet(i,j)
           enddo
        enddo

        ! Compute sum of inner integration weights for correction.
        sum_w = 0.0d0 ! or sum_w = sum(elem(ie)%mp(2:np-1,2:np-1)*elem(ie)%metdet(2:np-1,2:np-1))
        do j = 2, np-1
           do i = 2, np-1
              sum_w = sum_w + gp%weights(i)*gp%weights(j)*elem(ie)%metdet(i,j)
           enddo
        enddo
        ! Which tol is to use here?
        if ( sum_w > tol_local ) then
           delta = (area_elem - area_num)/sum_w
           call cube_init_atomic(elem(ie),gp%points,1.0d0 + delta)
           else
              ! Abort since the denominator in correction is too small.
              call abortmp('Cube_mod,set_area_correction_map2(): sum_w is too small.')
           endif
        enddo ! loop over elements
      else if ( geometry == "plane" ) then
           do ie=1,nelemd

              ! Store element's area in area_elem.
              area_elem = dx * dy

              ! Compute 'numerical area' of the element as sum of integration
              ! weights.
              area_num = 0.0d0
              do j = 1,np
                 do i = 1,np
                    area_num = area_num + gp%weights(i)*gp%weights(j)*elem(ie)%metdet(i,j)
                 enddo
              enddo

              ! Compute sum of inner integration weights for correction.
              sum_w = 0.0d0 ! or sum_w = sum(elem(ie)%mp(2:np-1,2:np-1)*elem(ie)%metdet(2:np-1,2:np-1))
              do j = 2, np-1
                 do i = 2, np-1
                    sum_w = sum_w + gp%weights(i)*gp%weights(j)*elem(ie)%metdet(i,j)
                 enddo
              enddo
              ! Which tol is to use here?
              if ( sum_w > tol_local ) then
                 delta = (area_elem - area_num)/sum_w
                 call plane_init_atomic(elem(ie),gp%points,1.0d0 + delta)
                 else
                    ! Abort since the denominator in correction is too small.
                    call abortmp('Cube_mod,set_area_correction_map2(): sum_w is too small.')
                 endif
              enddo ! loop over elements

      endif



        ! code for verification.
        area = 0.0d0
        do ie = 1,nelemd
           aratio(ie,1) = 0.0
           do j = 1,np
              do i = 1,np
                 aratio(ie,1) = aratio(ie,1) + gp%weights(i)*gp%weights(j)*elem(ie)%metdet(i,j)
              enddo
           enddo
        enddo

        call repro_sum(aratio, area, nelemd, nelemd, 1, commid=par%comm)
        if (par%masterproc) &
           write(iulog,'(a,f20.17)') "Epsilon bubble correction: Corrected area - domain size ",area(1) - domain_size


end subroutine set_area_correction_map2

end module geometry_mod
