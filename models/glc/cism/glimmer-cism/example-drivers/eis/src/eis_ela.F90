!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   eis_ela.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

module eis_ela
  !*FD climate forcing similar to the old Edinburgh Ice Sheet model
  !*FD Magnus Hagdorn, June 2004
  
  !*FD This module reproduces the ELA forcing used to drive the 
  !*FD Edinburgh Ice Sheet model. 

  use glimmer_ts
  use glimmer_global, only : fname_length

  type eis_ela_type
     !*FD Parameters for the EIS climate forcing
     real :: ela_a = 10821.                      !*FD ELA paramters
     real :: ela_b = -238.                       !*FD ELA paramters
     real :: ela_c = 1.312                       !*FD ELA paramters
     real :: zmax_mar = 1200.                    !*FD parameters describing how the MB
     real :: bmax_mar = 1.5                      !*FD varies around the ELA
     real :: zmax_cont = 500.                    !*FD parameters describing how the MB
     real :: bmax_cont = 0.3                     !*FD varies around the ELA
     real :: shelf_ablation = -1.0               !*FD ablation over ice shelf
     character(len=fname_length) :: fname=''     !*FD name of file containing ELA ts
     character(len=fname_length) :: ew_fname=''  !*FD name of file containing longitude dependance of ELA field
     type(glimmer_tseries) :: ela_ts             !*FD ELA time series 
     real,dimension(:,:),pointer :: ela => null()!*FD ELA field
  end type eis_ela_type

  private :: ela_lat, calc_mb

contains
  subroutine eis_ela_config(config,ela)
    !*FD get ELA configuration from config file
    use glimmer_config
    use glimmer_filenames, only : filenames_inputname
    implicit none
    type(eis_ela_type)           :: ela     !*FD ela data
    type(ConfigSection), pointer :: config  !*FD structure holding sections of configuration file   
    ! local variables
    type(ConfigSection), pointer :: section

    ela%ew_fname = ''
    ela%fname = ''

    call GetSection(config,section,'EIS ELA')
    if (associated(section)) then
       call GetValue(section,'ela_file',ela%fname)
       call GetValue(section,'ela_ew',ela%ew_fname)
       call GetValue(section,'zmax_mar',ela%zmax_mar)
       call GetValue(section,'bmax_mar',ela%bmax_mar)
       call GetValue(section,'zmax_cont',ela%zmax_cont)
       call GetValue(section,'bmax_cont',ela%bmax_cont)
       call GetValue(section,'ela_a',ela%ela_a)
       call GetValue(section,'ela_b',ela%ela_b)
       call GetValue(section,'ela_c',ela%ela_c)
       if (trim(ela%fname).ne.'') then
          ela%fname = trim(filenames_inputname(ela%fname))
       end if
       if (trim(ela%ew_fname).ne.'') then
          ela%ew_fname = trim(filenames_inputname(ela%ew_fname))
       end if
    end if
  end subroutine eis_ela_config
    
  subroutine eis_ela_printconfig(ela)
    !*FD print configuration to log
    use glimmer_log
    use glimmer_global, only : msg_length
    implicit none
    type(eis_ela_type)      :: ela   !*FD ela data
    ! local variables
    character(len=msg_length) :: message
    call write_log('EIS ELA')
    call write_log('-------')
    write(message,*) 'ELA file     : ',trim(ela%fname)
    call write_log(message)
    if (len(trim(ela%ew_fname)).ne.0) then
       write(message,*) 'ELA depends on longitude: ',trim(ela%ew_fname)
       call write_log(message)
    end if
    write(message,*) 'maritime zmax: ',ela%zmax_mar
    call write_log(message)
    write(message,*) 'maritime bmax: ',ela%bmax_mar
    call write_log(message)
    write(message,*) 'cont zmax    : ',ela%zmax_cont
    call write_log(message)
    write(message,*) 'cont bmax    : ',ela%bmax_cont
    call write_log(message)
    write(message,*) 'ela A        : ',ela%ela_a
    call write_log(message)
    write(message,*) 'ela B        : ',ela%ela_B
    call write_log(message)
    write(message,*) 'ela C        : ',ela%ela_C
    call write_log(message)
    call write_log('')
  end subroutine eis_ela_printconfig

  subroutine eis_init_ela(ela,model)
    !*FD initialise ELA forcing
    use glide_types
    use glimmer_paramets, only: thk0, acc0, scyr
    implicit none
    type(eis_ela_type)      :: ela   !*FD ela data
    type(glide_global_type) :: model !*FD model instance
    
    call glimmer_read_ts(ela%ela_ts,ela%fname)

    ! scale parameters
    ela%ela_ts%values = ela%ela_ts%values/thk0
    ela%zmax_mar = ela%zmax_mar/thk0
    ela%bmax_mar = ela%bmax_mar/ (acc0 * scyr)
    ela%zmax_cont = ela%zmax_cont/thk0
    ela%bmax_cont = ela%bmax_cont/ (acc0 * scyr)
    ela%shelf_ablation = ela%shelf_ablation / (acc0 * scyr)
    ela%ela_a = ela%ela_a/thk0
    ela%ela_b = ela%ela_b/thk0
    ela%ela_c = ela%ela_c/thk0

    ! calculate shape of mass balance field
    call coordsystem_allocate(model%general%ice_grid, ela%ela)
    ela%ela = 0.
  end subroutine eis_init_ela
    
  subroutine eis_calc_ela(ela,model)
    !*FD calculate latitude dependence of ELA field
    use glide_types
    use glimmer_paramets, only: thk0
    implicit none
    type(eis_ela_type)      :: ela   !*FD ela data
    type(glide_global_type) :: model !*FD model instance
    
    ! local variables
    integer ew,ns
    real :: ela_m
    type(glimmer_tseries) :: ela_ew

    ela%ela = ela%ela + ela_lat(ela%ela_a,ela%ela_b,ela%ela_c,model%climate%lati)
    if (len(trim(ela%ew_fname)).ne.0) then
       call glimmer_read_ts(ela_ew,ela%ew_fname)
       ! loop over grid
       do ns = 1,model%general%nsn
          do ew = 1,model%general%ewn
             call glimmer_ts_linear(ela_ew,model%climate%loni(ew,ns),ela_m)
             ela%ela(ew,ns) = ela%ela(ew,ns) + ela_m/thk0
          end do
       end do
    end if
  end subroutine eis_calc_ela

  subroutine eis_massbalance(ela,cony,model,time)
    !*FD calculate mass balance
    use eis_cony
    use glide_types
    use glimmer_global, only : rk
    implicit none
    type(eis_ela_type)        :: ela   !*FD ela data
    type(eis_cony_type)       :: cony  !*FD cony data
    type(glide_global_type)   :: model !*FD model instance
    real(kind=rk), intent(in) :: time  !*FD current time

    ! local variables
    real :: ela_time = 0.

    call glimmer_ts_step(ela%ela_ts,real(time),ela_time)

    model%climate%acab = calc_mb(ela%ela+ela_time, &
         model%geometry%topg, &
         model%geometry%thck, &
         cony%cony, &
         model%climate%eus, &
         ela%zmax_mar,ela%bmax_mar, &
         ela%zmax_cont,ela%bmax_cont, &
         ela%shelf_ablation)
  end subroutine eis_massbalance

  !*****************************************************************************
  ! private procedures
  !*****************************************************************************
  elemental function calc_mb(ela,topo,thick,cony,eus,mzmax,mbmax,czmax,cbmax,shelf_ablation)
    !*FD calculate mass balance
    use glimmer_global, only : dp
    implicit none
    real, intent(in) :: ela       !*FD equilibrium line altitude
    real(kind=dp), intent(in) :: topo      !*FD topography
    real(kind=dp), intent(in) :: thick     !*FD ice thickness
    real,intent(in)           :: cony      !*FD continentality
    real, intent(in) :: eus       !*FD eustatic sea level
    real, intent(in) :: mzmax,mbmax !*FD parameters describing MB variation around ELA
    real, intent(in) :: czmax,cbmax !*FD parameters describing MB variation around ELA
    real, intent(in) :: shelf_ablation !*FD ablation over ice shelf
    real calc_mb

    ! local variables
    real z
    real zmax,bmax

    if (topo.ge.eus .or. thick.gt.0) then
       z = topo+thick-eus
       if (z.lt.0.) then
          calc_mb = shelf_ablation   ! ablation on ice shelf
          return
       end if
    else
       calc_mb = 0
       return
    end if
    z = z - ela
    zmax = mzmax - (mzmax - czmax)*cony
    bmax = mbmax - (mbmax - cbmax)*cony
    if ((zmax.le.0.01).or.(z .ge. zmax)) then
       ! first condition allows forcing of mb=bmax.
       calc_mb = bmax
       return
    else
       calc_mb = 2.*z*bmax/zmax - z*z*bmax/(zmax*zmax)
       return
    end if
   end function calc_mb

  elemental function ela_lat(a,b,c,lat)
    !*FD calculate ELA variation with latitude
    implicit none
    real, intent(in) :: a,b,c !*FD shape of ELA field
    real,intent(in)  :: lat   !*FD latitude
    real ela_lat

    ela_lat = a +  b * lat +  c * lat * lat
  end function ela_lat
  
end module eis_ela
