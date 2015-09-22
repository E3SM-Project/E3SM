!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   eis_cony.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

module eis_cony
  !*FD climate forcing similar to the old Edinburgh Ice Sheet model
  !*FD Magnus Hagdorn, June 2004
  
  !*FD This module reproduces the continentality forcing used to drive the 
  !*FD Edinburgh Ice Sheet model. 

  use glimmer_searchcircle
  use glimmer_global, only: dp

  type eis_cony_type
     real(dp) :: period = 0.d0                        !*FD how often cony field should be updated, set to 0 to switch off
     real(dp) :: cony_radius = 600000                 !*FD continentality radius
     integer :: update_from_file = 0                  !*FD load cony from file
     real(dp) :: next_update                          !*FD when the next update occurs
     real(dp),dimension(:,:),pointer :: cony => null()!*FD cony field
     real(dp),dimension(:,:),pointer :: topo => null()!*FD topography mask
     type(searchdata) :: sdata                        !*FD search circle structure
  end type eis_cony_type

contains  
  subroutine eis_cony_config(config,cony)
    !*FD get ELA configuration from config file
    use glimmer_config
    implicit none
    type(eis_cony_type)          :: cony    !*FD cony data
    type(ConfigSection), pointer :: config  !*FD structure holding sections of configuration file   
    ! local variables
    type(ConfigSection), pointer :: section

    call GetSection(config,section,'EIS CONY')
    if (associated(section)) then
       cony%period = 500.d0
       call GetValue(section,'period',cony%period)
       call GetValue(section,'radius',cony%cony_radius)
       call GetValue(section,'file',cony%update_from_file)
       if (cony%update_from_file.eq.1) then
          cony%period = 0.d0
       end if
    end if
  end subroutine eis_cony_config

  subroutine eis_cony_printconfig(cony)
    !*FD print configuration to log
    use glimmer_log
    implicit none
    type(eis_cony_type)      :: cony   !*FD cony data
    ! local variables
    character(len=100) :: message
    call write_log('EIS Continentality')
    call write_log('------------------')
    write(message,*) 'Update Period        : ',cony%period
    call write_log(message)
    write(message,*) 'Continentality radius: ',cony%cony_radius
    call write_log(message)
    if (cony%update_from_file.eq.1) then
       call write_log(' Load continentality from input netCDF file')
    end if
    call write_log('')
  end subroutine eis_cony_printconfig

  subroutine eis_init_cony(cony,model)
    !*FD initialise cony forcing
    use glide_types
    use glimmer_paramets, only: len0
    implicit none
    type(eis_cony_type)     :: cony  !*FD cony data
    type(glide_global_type) :: model !*FD model instance

    integer radius

    ! scale parameters
    cony%cony_radius = cony%cony_radius / len0

    cony%next_update = model%numerics%tstart

    ! allocate data
    allocate(cony%cony(model%general%ewn,model%general%nsn))
    cony%cony = 0.d0

    if (cony%period .gt. 0.d0) then
       allocate(cony%topo(model%general%ewn,model%general%nsn))
       ! initialise search circles
       radius = int(cony%cony_radius/sqrt(model%numerics%dew**2 + model%numerics%dns**2))
       cony%sdata = sc_initdata(radius,1,1,model%general%ewn,model%general%nsn)
    end if
  end subroutine eis_init_cony
  
  subroutine eis_continentality(cony,model,time)
    !*FD calculate continentality field
    use glide_types
    use glimmer_global, only : dp
    implicit none
    type(eis_cony_type)       :: cony  !*FD cony data
    type(glide_global_type)   :: model !*FD model instance
    real(dp), intent(in) :: time  !*FD current time
    
    if (cony%period.gt.0.d0 .and. cony%next_update.ge.time) then
       
       where (.not. GLIDE_IS_OCEAN(model%geometry%thkmask) .and. .not. GLIDE_IS_FLOAT(model%geometry%thkmask))
          cony%topo = 1.d0
       elsewhere
          cony%topo = 0.d0
       end where

       call sc_search(cony%sdata,cony%topo,cony%cony)

       where(cony%topo.eq.0.d0)
          cony%cony = 0.d0
       elsewhere
          cony%cony =2.d0 * max(cony%cony/cony%sdata%total_area-0.5d0, 0.d0) 
       end where

       cony%next_update = cony%next_update+cony%period
    end if
  end subroutine eis_continentality
end module eis_cony
