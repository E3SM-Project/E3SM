!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   isostasy.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

module isostasy

  !TODO - Test the isostasy for parallel simulations.
  !       Elastic lithosphere take some work to parallelize, but local calculations
  !        should be OK.
  
  !*FD calculate isostatic adjustment due to changing surface loads

  use glimmer_global, only : dp
  use isostasy_elastic

  implicit none

  private :: relaxing_mantle

!-------------------------------------------------------------------------

contains

!-------------------------------------------------------------------------

  subroutine init_isostasy(model)

    !*FD initialise isostasy calculations
    use parallel
    use glide_types
    use glimmer_physcon,  only: scyr
    use glimmer_paramets, only: tim0
    implicit none

    type(glide_global_type) :: model

    if (model%isostasy%lithosphere == LITHOSPHERE_ELASTIC) then
       call not_parallel(__FILE__,__LINE__)
       call init_elastic(model%isostasy%rbel,model%numerics%dew)
    end if

    model%isostasy%next_calc = model%numerics%tstart
    model%isostasy%relaxed_tau = model%isostasy%relaxed_tau * scyr / tim0

  end subroutine init_isostasy

!-------------------------------------------------------------------------
  
  subroutine isos_icewaterload(model)

    !*FD calculate surface load factors due to water and ice distribution

    use glimmer_physcon
    use glide_types
    implicit none

    type(glide_global_type) :: model

    real(dp) :: ice_mass, water_depth, water_mass
    integer :: ew,ns
  
     do ns=1,model%general%nsn
       do ew=1,model%general%ewn
          ice_mass = rhoi * model%geometry%thck(ew,ns)

          if (model%geometry%topg(ew,ns) - model%climate%eus < 0.d0) then   ! check if we are below sea level

             water_depth = model%climate%eus - model%geometry%topg(ew,ns)
             water_mass = rhoo * water_depth

             ! Just the water load due to changes in sea-level
             model%isostasy%load_factors(ew,ns) = rhoo* model%climate%eus/rhom

             ! Check if ice is not floating
             if ( ice_mass > water_mass ) then
                model%isostasy%load_factors(ew,ns) = model%isostasy%load_factors(ew,ns) + (ice_mass - water_mass)/rhom
             end if

          else                                       ! bedrock is above sea level

             model%isostasy%load_factors(ew,ns) = ice_mass/rhom

          end if

       end do
    end do

  end subroutine isos_icewaterload

!-------------------------------------------------------------------------

  subroutine isos_compute(model)

    !*FD calculate isostatic adjustment due to changing surface loads

    use glide_types
    implicit none

    type(glide_global_type) :: model

    ! update load if necessary
    if (model%isostasy%new_load) then
       call isos_lithosphere(model, model%isostasy%load, model%isostasy%load_factors)
       ! update bedrock with (non-viscous) fluid mantle
       if (model%isostasy%asthenosphere == ASTHENOSPHERE_FLUID) then
          model%geometry%topg = model%isostasy%relx - model%isostasy%load
       end if
       model%isostasy%new_load = .false.
    end if

    ! update bedrock with relaxing mantle
    if (model%isostasy%asthenosphere == ASTHENOSPHERE_RELAXING) then
       call relaxing_mantle(model)
    end if

  end subroutine isos_compute

!-------------------------------------------------------------------------

  subroutine isos_lithosphere(model,load,load_factors)

    use glide_types
    implicit none
    type(glide_global_type) :: model
    real(dp), dimension(:,:), intent(out) :: load !*FD loading effect due to load_factors
    real(dp), dimension(:,:), intent(in)  :: load_factors !*FD load mass divided by mantle density

    if (model%isostasy%lithosphere == LITHOSPHERE_LOCAL) then
       load = load_factors
    else if (model%isostasy%lithosphere == LITHOSPHERE_ELASTIC) then
       call calc_elastic(model%isostasy%rbel, load, load_factors)
    end if

  end subroutine isos_lithosphere

!-------------------------------------------------------------------------

  subroutine isos_relaxed(model)

    !*FD Calculate the relaxed topography, assuming the isostatic depression
    !*FD is the equilibrium state for the current topography.

    use glide_types
    implicit none
    type(glide_global_type) :: model

    ! Calculate the load
    call isos_icewaterload(model)

    ! Apply lithosphere model
    call isos_lithosphere(model, model%isostasy%load, model%isostasy%load_factors)

    ! Add to present topography to get relaxed topography
    model%isostasy%relx = model%geometry%topg + model%isostasy%load

  end subroutine isos_relaxed

!-------------------------------------------------------------------------
! private subroutines
!-------------------------------------------------------------------------

  subroutine relaxing_mantle(model)

    !*FD approximate mantle with a relaxing half-space: dh/dt=-1/tau*(w-h)
    use glide_types
    implicit none
    type(glide_global_type) :: model
    
    integer :: ew,ns
    real(dp) :: ft1, ft2

    ft1 = exp(-model%numerics%dt/model%isostasy%relaxed_tau)
    ft2 = 1. - ft1  !TODO - 1.d0

    do ns=1,model%general%nsn
       do ew=1,model%general%ewn
          model%geometry%topg(ew,ns) = ft2 * (model%isostasy%relx(ew,ns) - model%isostasy%load(ew,ns)) &
                                     + ft1 *  model%geometry%topg(ew,ns)
       end do
    end do

  end subroutine relaxing_mantle

!-------------------------------------------------------------------------

end module isostasy

!-------------------------------------------------------------------------
