!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   erosion_types.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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
#include <config.inc>
#endif

module erosion_types
  !*FD type definition for erosion calcluations

  use glimmer_global, only : dp  
  use glimmer_sparse
  use glimmer_coordinates
  use erosion_transport_type

  type er_prof_type
     integer :: erate
     integer :: calc_lag
     integer :: trans_sed
     integer :: sed_eros
     integer :: sed_dep
  end type er_prof_type

  type er_sed_type
     real(kind=dp) :: effective_pressure = 50              !*FD effective pressure in [kPa]
     real(kind=dp) :: eff_press_grad     = -10.            !*FD effective pressure gradient [kPa/m]
     real(kind=dp) :: phi                = 30              !*FD angle of internal friction [degree]
     real(kind=dp) :: c                  = 15              !*FD cohesion [kPa]
     real(kind=dp) :: alpha, beta                          !*FD parameters calculated
     real(kind=dp) :: a = 107.11                           !*FD factor for flow law
     real(kind=dp) :: m = 1.35                             !*FD exponent of effective pressure
     real(kind=dp) :: n = 0.77                             !*FD exponent of shear stress
     logical :: calc_btrc = .False.                        !*FD set to .True. if sliding velos should be calculated by
                                                           !*FD sediment module

     real(kind=dp),dimension(11) :: params                 !*FD parameters for sediment flow law

  end type er_sed_type
     
  type erosion_type
     logical :: doerosion = .False.                        !*FD set to true when erosion should be included
     integer :: ndt = 1                                    !*FD erosion time step (multiplier of main time step)
     real(kind=dp) :: dt                                   !*FD erosion time step
     real(kind=dp) :: hb_erosion_factor =   1.d-10         !*FD constant of proportionality for erosion rate calcs
     real :: density = 3000.                               !*FD density of hard bedrock (kg m$^{-3}$)
     ! sediment transport stuff
     logical :: dotransport = .False.                      !*FD set to true to move sediments about
     logical :: simple_seds = .false.                      !*FD toggle between sediment models
     logical :: update_topo = .True.                       !*FD whether erosion/sediment transport changes bedrock topo
     real(kind=dp) :: transport_fac = 0.2                  !*FD multiplier for velos in deformable beds
     real(kind=dp) :: dirty_ice_max = 0.1                  !*FD maximum thickness of dirty basal ice layer
     real(kind=dp) :: soft_b = 1.d-5                       !*FD param A for max def thick calculations
     real(kind=dp) :: soft_a = 0.                          !*FD param B for max def thick calculations
     ! internal fields, etc
     type(er_transport_type) :: trans                      !*FD type holding transport stuff
     type(er_sed_type) :: sediment                         !*FD deforming sediment layer stuff
     real(kind=dp), dimension(:,:), pointer :: tau_mag => null()    !*FD magnitude of basal shear stress
     real(kind=dp), dimension(:,:), pointer :: tau_dir => null()    !*FD direction of basal shear stress
     real(kind=dp),dimension(:,:),pointer :: erosion_rate => null() !*FD hard bedrock erosion rate
     real(kind=dp),dimension(:,:),pointer :: erosion => null()      !*FD total hard bedrock erosion
     real(kind=dp),dimension(:,:),pointer :: er_accu => null()      !*FD accumulated erosion during one erosion time step
     real(kind=dp),dimension(:,:),pointer :: er_isos => null()      !*FD accumulated erosion for isostasy calcs
     real(kind=dp),dimension(:,:),pointer :: er_load => null()      !*FD load due to erosion
     real(kind=dp),dimension(:,:),pointer :: seds_init => null()    !*FD initial sediment distribution
     real(kind=dp),dimension(:,:),pointer :: seds1 => null()        !*FD thickness of dirty basal ice layer
     real(kind=dp),dimension(:,:),pointer :: seds2 => null()        !*FD thickness of deforming sediment layer
     real(kind=dp),dimension(:,:),pointer :: seds2_vx => null()     !*FD x-component of sediment velocity
     real(kind=dp),dimension(:,:),pointer :: seds2_vy => null()     !*FD y-component of sediment velocity
     real(kind=dp),dimension(:,:),pointer :: seds2_max => null()    !*FD maximum thickness of deforming sediment layer
     real(kind=dp),dimension(:,:),pointer :: seds3 => null()        !*FD thickness of non-deforming sediment layer
     ! profiling
     type(er_prof_type) :: er_prof
  end type erosion_type

contains
  subroutine er_allocate(erosion,model)
    !*FD allocate erosion data
    use glide_types
    implicit none
    type(erosion_type) :: erosion     !*FD data structure holding erosion stuff
    type(glide_global_type) :: model  !*FD model instance

    call coordsystem_allocate(model%general%velo_grid, erosion%tau_mag)
    call coordsystem_allocate(model%general%velo_grid, erosion%tau_dir)
    call coordsystem_allocate(model%general%ice_grid,  erosion%er_isos)
    call coordsystem_allocate(model%general%ice_grid,  erosion%er_load)
    call coordsystem_allocate(model%general%velo_grid, erosion%erosion_rate)
    call coordsystem_allocate(model%general%velo_grid, erosion%erosion)
    call coordsystem_allocate(model%general%velo_grid, erosion%er_accu)
    call coordsystem_allocate(model%general%velo_grid, erosion%seds_init)
    call coordsystem_allocate(model%general%velo_grid, erosion%seds1)
    call coordsystem_allocate(model%general%velo_grid, erosion%seds2)
    call coordsystem_allocate(model%general%velo_grid, erosion%seds2_vx)
    call coordsystem_allocate(model%general%velo_grid, erosion%seds2_vy)
    call coordsystem_allocate(model%general%velo_grid, erosion%seds2_max)
    call coordsystem_allocate(model%general%velo_grid, erosion%seds3)
  end subroutine er_allocate
    
  subroutine er_deallocate(erosion)
    !*FD free memory used by erosion data structure
    implicit none
    type(erosion_type) :: erosion     !*FD data structure holding erosion stuff

    deallocate(erosion%erosion_rate)
    deallocate(erosion%erosion)
    deallocate(erosion%er_accu)
    deallocate(erosion%er_isos)
    deallocate(erosion%er_load)
    deallocate(erosion%seds_init)
    deallocate(erosion%seds1)
    deallocate(erosion%seds2)
    deallocate(erosion%seds2_max)
    deallocate(erosion%seds3)
  end subroutine er_deallocate
end module erosion_types
