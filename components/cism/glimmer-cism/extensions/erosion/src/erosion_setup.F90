!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   erosion_setup.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

module erosion_setup
contains
  subroutine er_readconfig(erosion,config)
    !*FD read erosion configuration
    use erosion_types
    use glimmer_config
    implicit none
    type(erosion_type) :: erosion           !*FD structure holding erosion data
    type(ConfigSection), pointer :: config  !*FD structure holding sections of configuration file   
    ! local variables
    type(ConfigSection), pointer :: section
    integer :: dummy

    call GetSection(config,section,'Erosion')
    if (associated(section)) then
       erosion%doerosion = .True.
       call GetValue(section,'hb_erosion',erosion%hb_erosion_factor)
       call GetValue(section,'ntime',erosion%ndt)
       call GetValue(section,'update_topo',erosion%update_topo)
    end if
    call GetSection(config,section,'Basic_Transport')
    if (associated(section)) then
       erosion%dotransport = .True.
       erosion%simple_seds = .True.       
       call GetValue(section,'deformable_velo',erosion%transport_fac)
       call GetValue(section,'dirty_ice_thick',erosion%dirty_ice_max)
       call GetValue(section,'soft_a',erosion%soft_a)
       call GetValue(section,'soft_b',erosion%soft_b)
    end if
    call GetSection(config,section,'Transport')
    if (associated(section)) then
       erosion%dotransport = .True.
       erosion%simple_seds = .False.
       call GetValue(section,'dirty_ice_thick',erosion%dirty_ice_max)
       dummy = 0
       call GetValue(section,'calc_btrc',dummy)
       erosion%sediment%calc_btrc=(dummy.eq.1)
       call GetValue(section,'effective_pressure',erosion%sediment%effective_pressure)
       call GetValue(section,'pressure_gradient',erosion%sediment%eff_press_grad)
       call GetValue(section,'phi',erosion%sediment%phi)
       call GetValue(section,'cohesion',erosion%sediment%c)
       call GetValue(section,'a',erosion%sediment%a)
       call GetValue(section,'m',erosion%sediment%m)
       call GetValue(section,'n',erosion%sediment%n)
    end if
  end subroutine er_readconfig

  subroutine er_printconfig(erosion)
    !*FD print erosion config to log
    use glimmer_log
    use erosion_types
    implicit none
    type(erosion_type) :: erosion           !*FD structure holding erosion data
    ! local variables
    character(len=100) :: message

    if (erosion%doerosion) then
       call write_log('Erosion')
       call write_log('-------')
       write(message,*) 'Updating erosion every ',erosion%ndt,' time steps'
       call write_log(message)
       write(message,*) 'hard bedrock erosion constant : ',erosion%hb_erosion_factor
       call write_log(message)
       if (erosion%update_topo) then
          call write_log('erosion modifies bedrock topography')
       else
          call write_log('erosion does not modifies bedrock topography')
       end if
       call write_log('')
       if (erosion%dotransport) then
          if (erosion%simple_seds) then
             call write_log('Basic Sediment Transport')
             call write_log('------------------------')
             write(message,*) 'deformable sediment velo factor: ',erosion%transport_fac
             call write_log(message)
             write(message,*) 'max thickness of dirty basal ice layer: ',erosion%dirty_ice_max
             call write_log(message)
             write(message,*) 'deformable sediment param a: ',erosion%soft_a
             call write_log(message)
             write(message,*) 'deformable sediment param b: ',erosion%soft_b
             call write_log(message)
          else
             call write_log('Full Sediment Transport')
             call write_log('-----------------------')
             if (erosion%sediment%calc_btrc) then
                call write_log('Calculate basal sliding velocities')
             end if
             write(message,*) 'max thickness of dirty basal ice layer: ',erosion%dirty_ice_max
             call write_log(message)
             write(message,*) 'effective pressure : ',erosion%sediment%effective_pressure, ' kPa'
             call write_log(message)
             write(message,*) 'effective pressure gradient : ',erosion%sediment%eff_press_grad, 'kPa/m'
             call write_log(message)
             write(message,*) 'angle of internal friction : ',erosion%sediment%phi, ' degree'
             call write_log(message)
             write(message,*) 'cohesion : ',erosion%sediment%c, ' kPa'
             call write_log(message)
              write(message,*) 'factor for flow law (a): ',erosion%sediment%a
             call write_log(message)
             write(message,*) 'exponent of effective pressure (m) : ',erosion%sediment%m
             call write_log(message)
             write(message,*) 'exponent of shear stress (n) : ',erosion%sediment%n
             call write_log(message)
          end if
          call write_log('')
      end if
    end if
  end subroutine er_printconfig

  subroutine erosion_prof_init(model,erosion)
    !*FD initialise profiling
    use profile
    use glide_types
    use erosion_types
    implicit none
    type(glide_global_type) :: model        !*FD model instance
    type(erosion_type) :: erosion           !*FD structure holding erosion data
    
    if (model%profile%profile_unit .eq. 0) then
       call profile_init(model%profile,'erosion.profile')
       write(model%profile%profile_unit,*) '# take a profile every ',model%numerics%profile_period,' time steps'
    end if

    erosion%er_prof%erate     = profile_register(model%profile,'erosion rate')
    erosion%er_prof%calc_lag  = profile_register(model%profile,'calc lagrangian')
    erosion%er_prof%trans_sed = profile_register(model%profile,'trans sed')
    erosion%er_prof%sed_eros  = profile_register(model%profile,'erode sediments')
    erosion%er_prof%sed_dep   = profile_register(model%profile,'deposit sediments')
  end subroutine erosion_prof_init

end module erosion_setup
