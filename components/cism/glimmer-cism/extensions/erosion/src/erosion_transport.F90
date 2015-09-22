!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   erosion_transport.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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
! this module transports some scalar quantity c through a 2D velo field using
! the conservation of 2nd order moments algorithm

#ifdef HAVE_CONFIG_H
#include <config.inc>
#endif

module erosion_transport
  
contains
  subroutine init_transport(trans,model,erosion)
    use glide_types
    use erosion_types
    implicit none
    type(er_transport_type) :: trans       ! structure holding transport stuff
    type(glide_global_type) :: model       ! model instance
    type(erosion_type) :: erosion          !*FD structure holding erosion data

    call init_advect(trans%mo_seds1,erosion%seds1,model%numerics%dew,model%numerics%dns)
    call init_advect(trans%mo_seds2,erosion%seds2,model%numerics%dew,model%numerics%dns)

  end subroutine init_transport

  subroutine transport_sediments(erosion,model)
    use erosion_types
    use glide_types
    use erosion_sediment
    
    implicit none
    type(erosion_type) :: erosion          !*FD structure holding erosion data
    type(glide_global_type) :: model

    integer ns,ew


    ! transform basal shear from cartesian to radial components
    call er_trans_tau(erosion,model)

    ! transport in ice base    
    call set_species(erosion%trans%mo_seds1,erosion%seds1)
    call advect(erosion%trans%mo_seds1,model%velocity%ubas,model%velocity%vbas,erosion%dt)
    call get_species(erosion%trans%mo_seds1,erosion%seds1)

    ! transport in deformable sediment layer
    if (erosion%simple_seds) then
       call er_calc_dthick(erosion,model)
       erosion%seds2_vx = erosion%transport_fac*model%velocity%ubas
       erosion%seds2_vy = erosion%transport_fac*model%velocity%vbas
    else
       call er_sediment_tstep(erosion,model)
    end if
    call set_species(erosion%trans%mo_seds2,erosion%seds2)
    call advect(erosion%trans%mo_seds2,erosion%seds2_vx,erosion%seds2_vy,erosion%dt)
    call get_species(erosion%trans%mo_seds2,erosion%seds2)
    
  end subroutine transport_sediments

  subroutine er_calc_dthick(erosion,model)
    !*FD calculate thickness of deformable sediment bed
    use glide_types
    use glide_velo
    use erosion_types
    implicit none
    type(erosion_type) :: erosion          !*FD structure holding erosion data
    type(glide_global_type) :: model       !*FD model instance

    integer ew,ns

    erosion%seds2_max = 0.
    do ns=1,model%general%nsn-1
       do ew=1,model%general%ewn-1
          if (abs(model%velocity%ubas(ew,ns))+abs(model%velocity%vbas(ew,ns)) .gt. 0.) then
             erosion%seds2_max(ew,ns) = erosion%soft_a + erosion%soft_b * erosion%tau_mag(ew,ns)
          end if
       end do
    end do    
  end subroutine er_calc_dthick

end module erosion_transport
