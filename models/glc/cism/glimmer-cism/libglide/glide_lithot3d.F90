!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glide_lithot3d.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2013
!   Glimmer-CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of Glimmer-CISM.
!
!   Glimmer-CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   Glimmer-CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with Glimmer-CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

#include "glide_mask.inc"

! module for 3D temperature calculations in the upper lithosphere
! (serial only)

!TODO - Is this something we want to support in parallel?

module glide_lithot3d

  implicit none

  private
  public :: init_lithot3d, calc_lithot3d, finalise_lithot3d


contains  

  subroutine init_lithot3d(model)

    use glide_types
    use glimmer_paramets, only: len0,tim0
    implicit none
    type(glide_global_type),intent(inout) :: model       !*FD model instance

    ! local variables
    integer i,j,k,r,icount,jcount,ewn,nsn

    ! allocate memory for 3D code
    ewn=model%general%ewn
    nsn=model%general%nsn

    call new_sparse_matrix(ewn*nsn*model%lithot%nlayer, &
                          (model%lithot%nlayer-1)*ewn*nsn*7+ewn*nsn+1,model%lithot%fd_coeff)
    call new_sparse_matrix(ewn*nsn*model%lithot%nlayer, &
                           (model%lithot%nlayer-1)*ewn*nsn*7+ewn*nsn+1,model%lithot%fd_coeff_slap)
    allocate(model%lithot%rhs(model%lithot%nlayer*ewn*nsn))
    allocate(model%lithot%answer(model%lithot%nlayer*ewn*nsn))
    model%lithot%mxnelt = 20 * model%lithot%nlayer*ewn*nsn
    allocate(model%lithot%rwork(model%lithot%mxnelt))
    allocate(model%lithot%iwork(model%lithot%mxnelt))

    ! set up factors for horizontal finite differences
    model%lithot%xfactor = 0.5*model%lithot%diffu*tim0*model%numerics%dt / (model%numerics%dew*len0)**2
    model%lithot%yfactor = 0.5*model%lithot%diffu*tim0*model%numerics%dt / (model%numerics%dns*len0)**2


    ! calculate finite difference coefficient matrix
    ! top face
    ! simply match air temperature where no ice and basal temperature where ice
    k = 1
    do j=1,model%general%nsn
       do i=1,model%general%ewn
          r = linearise(model,i,j,k)
          call sparse_insert_val(model%lithot%fd_coeff,r,r, 1.d0)
       end do
    end do
    do k=2, model%lithot%nlayer-1
       do j=1,model%general%nsn
          do i=1,model%general%ewn
             icount = 0
             jcount = 0
             r = linearise(model,i,j,k)
             ! i-1,j,k
             if (i /= 1) then
                call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i-1,j,k), -model%lithot%xfactor)
                icount = icount + 1
             end if
             ! i+1, j, k
             if (i /= model%general%ewn) then
                call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i+1,j,k), -model%lithot%xfactor)
                icount = icount + 1
             end if
             ! i,j-1,k
             if (j /= 1) then
                call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j-1,k), -model%lithot%yfactor)
                jcount = jcount + 1
             end if
             ! i,j+1,k
             if (j /= model%general%nsn) then
                call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j+1,k), -model%lithot%yfactor)
                jcount = jcount + 1
             end if
             ! i,j,k-1
             call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j,k-1), -model%lithot%zfactors(1,k))
             ! i,j,k+1
             call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j,k+1), -model%lithot%zfactors(3,k))
             ! i,j,k
             call sparse_insert_val(model%lithot%fd_coeff,r,r, &
                  icount*model%lithot%xfactor + jcount*model%lithot%yfactor + model%lithot%zfactors(2,k) + 1.)
          end do
       end do
    end do
    
    ! bottom face
    ! keep constant
    k = model%lithot%nlayer
    do j=1,model%general%nsn
       do i=1,model%general%ewn
          r = linearise(model,i,j,k)
          call sparse_insert_val(model%lithot%fd_coeff,r,r, 1.d0)
       end do
    end do

    ! convert from SLAP Triad to SLAP Column format
    call copy_sparse_matrix(model%lithot%fd_coeff,model%lithot%fd_coeff_slap)
    call ds2y(model%general%nsn*model%general%ewn*model%lithot%nlayer,model%lithot%fd_coeff_slap%nonzeros, &
         model%lithot%fd_coeff_slap%col,model%lithot%fd_coeff_slap%row,model%lithot%fd_coeff_slap%val, 0)

    ! initialise result vector
    do k=1,model%lithot%nlayer
       do j=1,model%general%nsn
          do i=1,model%general%ewn
             model%lithot%answer(linearise(model,i,j,k)) = model%lithot%temp(i,j,k)
          end do
       end do
    end do
  end subroutine init_lithot3d

  subroutine calc_lithot3d(model)
    use glide_types
    use glide_stop
    use glimmer_log
    implicit none
    type(glide_global_type),intent(inout) :: model       !*FD model instance

    integer i,j,k,r
    integer iter
    real(dp) err
    real(dp), parameter :: tol = 1.0d-12
    integer, parameter :: isym = 0, itol = 2, itmax = 101
    integer :: ierr

    ! calculate RHS
    call sparse_matrix_vec_prod(model%lithot%fd_coeff,model%lithot%answer,model%lithot%rhs)
    model%lithot%rhs = -model%lithot%rhs + 2. * model%lithot%answer
    ! calc RHS on upper boundary
    k = 1
    do j=1,model%general%nsn
       do i=1,model%general%ewn
          r = linearise(model,i,j,k)
          if (GLIDE_IS_GROUND(model%geometry%thkmask(i,j)) .and. .not. GLIDE_IS_THIN(model%geometry%thkmask(i,j)) ) then
             model%lithot%rhs(r) = model%temper%temp(model%general%upn,i,j) ! ice basal temperature
             model%lithot%mask(i,j) = .true.
          else
             if (model%lithot%mask(i,j)) then
                if (GLIDE_IS_OCEAN(model%geometry%thkmask(i,j))) then
                   model%lithot%rhs(r) = model%lithot%mart
                else if (GLIDE_IS_LAND(model%geometry%thkmask(i,j))) then
                   model%lithot%rhs(r) = model%climate%artm(i,j) ! air temperature outside ice sheet
                end if
             end if
          end if
       end do
    end do

    ! solve matrix equation
    call dslucs(model%general%nsn*model%general%ewn*model%lithot%nlayer, model%lithot%rhs, model%lithot%answer, &
         model%lithot%fd_coeff_slap%nonzeros, model%lithot%fd_coeff_slap%col,model%lithot%fd_coeff_slap%row, &
         model%lithot%fd_coeff_slap%val, isym,itol,tol,itmax,iter,err,ierr,0, &
         model%lithot%rwork, model%lithot%mxnelt, model%lithot%iwork, model%lithot%mxnelt)

    if (ierr /= 0) then
      print *, 'pcg error ', ierr, itmax, iter
      write(*,*) model%numerics%time
      call glide_finalise(model,.true.)
      call close_log
      stop
    end if

    ! de-linearise results
    do k=1, model%lithot%nlayer
       do j=1,model%general%nsn
          do i=1,model%general%ewn
             model%lithot%temp(i,j,k) = model%lithot%answer(linearise(model,i,j,k))
          end do
       end do
    end do
      
  end subroutine calc_lithot3d

  subroutine finalise_lithot3d(model)
    use glide_types
    implicit none
    type(glide_global_type),intent(inout) :: model       !*FD model instance

    call del_sparse_matrix(model%lithot%fd_coeff)
    call del_sparse_matrix(model%lithot%fd_coeff_slap)
    deallocate(model%lithot%rhs)
    deallocate(model%lithot%answer)
  end subroutine finalise_lithot3d

  function linearise(model,i,j,k)
    use glide_types
    implicit none
    type(glide_global_type),intent(in) :: model   
    integer, intent(in) :: i,j,k
    integer :: linearise
    
    linearise = i + (j-1)*model%general%ewn + (k-1)*model%general%ewn*model%general%nsn
  end function linearise


end module glide_lithot3d
