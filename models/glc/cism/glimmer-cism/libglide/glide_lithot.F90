!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glide_lithot.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

! module for temperature calculations in the upper lithosphere

!TODO - Any changes needed for parallel code?

module glide_lithot

    implicit none

contains  

  subroutine init_lithot(model)
    use glide_types
    use glide_setup
    use glimmer_paramets, only: tim0
    use glimmer_log
    use glide_lithot1d
    use glide_lithot3d
    implicit none
    type(glide_global_type),intent(inout) :: model       !*FD model instance

    ! local variables
    integer k
    real(dp) :: factor

    ! allocate memory for common arrays
    allocate(model%lithot%deltaz(model%lithot%nlayer)); model%lithot%deltaz = 0.0
    allocate(model%lithot%zfactors(3,model%lithot%nlayer)); model%lithot%zfactors = 0.0    

    ! set up vertical grid
    do k=1,model%lithot%nlayer
       model%lithot%deltaz(k) = (1-glide_calc_sigma(1D0*real(model%lithot%nlayer-k)/real(model%lithot%nlayer-1),2D0)) &
            *model%lithot%rock_base
    end do

    ! calculate diffusion coefficient
    model%lithot%diffu = model%lithot%con_r/(model%lithot%rho_r*model%lithot%shc_r)

    ! set up factors for vertical finite differences
    do k=2,model%lithot%nlayer-1
       model%lithot%zfactors(1,k) =  model%lithot%diffu*tim0*model%numerics%dt / &
            ((model%lithot%deltaz(k)-model%lithot%deltaz(k-1)) * (model%lithot%deltaz(k+1)-model%lithot%deltaz(k-1)))
       model%lithot%zfactors(2,k) = model%lithot%diffu*tim0*model%numerics%dt / &
            ((model%lithot%deltaz(k+1)-model%lithot%deltaz(k)) * (model%lithot%deltaz(k)-model%lithot%deltaz(k-1)))
       model%lithot%zfactors(3,k) = model%lithot%diffu*tim0*model%numerics%dt / &
            ((model%lithot%deltaz(k+1)-model%lithot%deltaz(k)) * (model%lithot%deltaz(k+1)-model%lithot%deltaz(k-1)))
    end do
    k = model%lithot%nlayer
    model%lithot%zfactors(:,k) = 0.5*model%lithot%diffu*tim0*model%numerics%dt / &
         (model%lithot%deltaz(k)-model%lithot%deltaz(k-1))**2

    !TODO - Make sure the sign is correct here.
    !NOTE: Glimmer-CISM convention is that geot is positive down, so geot < 0 for upward geothermal flux

    if (model%options%is_restart == RESTART_FALSE) then
       ! set initial temp distribution to thermal gradient
       factor = model%paramets%geot / model%lithot%con_r
       do k=1,model%lithot%nlayer
          model%lithot%temp(:,:,k) = model%lithot%surft + model%lithot%deltaz(k)*factor
       end do
    end if

    if (model%lithot%num_dim==1) then
       call init_lithot1d(model)
    else if (model%lithot%num_dim==3) then
       call init_lithot3d(model)
    else
       call write_log('Error, init_lithot: Wrong number of dimensions',GM_FATAL,__FILE__,__LINE__)
    end if

  end subroutine init_lithot    

  subroutine spinup_lithot(model)
    use parallel
    use glide_types
    use glimmer_log
    use glide_mask
    implicit none
    type(glide_global_type),intent(inout) :: model       !*FD model instance

    integer t

    if (model%options%is_restart == RESTART_FALSE .and. model%lithot%numt > 0) then
       call write_log('Spinning up GTHF calculations',type=GM_INFO)
       call not_parallel(__FILE__,__LINE__)
       do t=1,model%lithot%numt
          call calc_lithot(model)
       end do

    end if
  end subroutine spinup_lithot

  !TODO - Pretty sure that calc_lithot3d is not parallel.  What about 1d?

  subroutine calc_lithot(model)
    use glide_types
    use glimmer_log
    use glide_lithot1d
    use glide_lithot3d
    implicit none
    type(glide_global_type),intent(inout) :: model       !*FD model instance

    if (model%lithot%num_dim==1) then
       call calc_lithot1d(model)
    else if (model%lithot%num_dim==3) then
       call calc_lithot3d(model)
    else
       call write_log('Wrong number of dimensions.',GM_FATAL,__FILE__,__LINE__)
    end if
      
    call calc_geoth(model)

  end subroutine calc_lithot

  subroutine calc_geoth(model)
    !*FD calculate geothermal heat flux
    use glide_types
    implicit none
    type(glide_global_type),intent(inout) :: model       !*FD model instance

    real(dp) factor

    factor = model%lithot%con_r/(model%lithot%deltaz(2)-model%lithot%deltaz(1))
    model%temper%bheatflx(:,:) = factor*(model%lithot%temp(:,:,2)-model%lithot%temp(:,:,1))

  end subroutine calc_geoth

  subroutine finalise_lithot(model)
    use glide_types
    use glide_lithot1d
    use glimmer_log
    use glide_lithot3d
    implicit none
    type(glide_global_type),intent(inout) :: model       !*FD model instance

    deallocate(model%lithot%deltaz)
    deallocate(model%lithot%zfactors)

    if (model%lithot%num_dim==1) then
       call finalise_lithot1d(model)
    else if (model%lithot%num_dim==3) then
       call finalise_lithot3d(model)
    else
       call write_log('Wrong number of dimensions.',GM_FATAL,__FILE__,__LINE__)
    end if
  end subroutine finalise_lithot

end module glide_lithot
