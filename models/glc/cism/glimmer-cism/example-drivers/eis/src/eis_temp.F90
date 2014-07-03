!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   eis_temp.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

module eis_temp
  !*FD climate forcing similar to the old Edinburgh Ice Sheet model
  !*FD Magnus Hagdorn, June 2004
  !*FD temperature forcing

  use glimmer_ts
  use glimmer_global, only : fname_length

  type eis_temp_type
     real :: lapse_rate  = -0.008                !*FD lapse rate
     character(len=fname_length) :: fname=''     !*FD name of file containing temperature ts
     integer :: lat_type = 0                     !*FD type of lat dependend function (0 polynomial,
                                                 !*FD 1 exponential: a+b*exp(c*(lat-lat0))
     real :: lat0 = 44.95                        !*FD parameter used for exponential function
     integer :: torder = 2                       !*FD order of temperature polynomial
     type(glimmer_tseries) :: temp_ts            !*FD temperature time series 
     real, dimension(:),pointer :: tvalue        !*FD temperature value
  end type eis_temp_type

contains
   subroutine eis_temp_config(config,temp)
    !*FD get temperature configuration from config file
    use glimmer_config
    use glimmer_filenames, only : filenames_inputname
    implicit none
    type(eis_temp_type)           :: temp     !*FD ela data
    type(ConfigSection), pointer :: config  !*FD structure holding sections of configuration file   
    ! local variables
    type(ConfigSection), pointer :: section

    temp%fname=''
    
    call GetSection(config,section,'EIS Temperature')
    if (associated(section)) then
       call GetValue(section,'temp_file',temp%fname)
       call GetValue(section,'type',temp%lat_type)
       if (temp%lat_type.eq.1) then
          call GetValue(section,'lat0',temp%lat0)
          temp%torder = 3
       else
          call GetValue(section,'order',temp%torder)
       end if
       call GetValue(section,'lapse_rate',temp%lapse_rate)
       if (trim(temp%fname).ne.'') then
          temp%fname = trim(filenames_inputname(temp%fname))
       end if
    end if
  end subroutine eis_temp_config

  subroutine eis_temp_printconfig(temp)
    !*FD print configuration to log
    use glimmer_log
    use glimmer_global, only : msg_length
    implicit none
    type(eis_temp_type)      :: temp   !*FD temperature data
    ! local variables
    character(len=msg_length) :: message
    call write_log('EIS Temperature')
    call write_log('---------------')
    write(message,*) 'temperature file  : ',trim(temp%fname)
    call write_log(message)
    if (temp%lat_type.eq.1) then
       call write_log('lat dependance is exponential')
       write(message,*) 'lat0               : ',temp%lat0
       call write_log(message)
    else
       call write_log('lat dependance is polynomial')
       write(message,*) 'order of temp poly: ',temp%torder
       call write_log(message)
    end if
    write(message,*) 'lapse rate        : ',temp%lapse_rate
    call write_log(message)
    call write_log('')
  end subroutine eis_temp_printconfig

  subroutine eis_init_temp(temp)
    !*FD initialise temperature forcing
    use glide_types
    use glimmer_paramets, only: thk0
    implicit none
    type(eis_temp_type)     :: temp  !*FD temp data

    call glimmer_read_ts(temp%temp_ts,temp%fname,temp%torder)
    allocate(temp%tvalue(temp%torder))
    ! scale parameters
    temp%lapse_rate = temp%lapse_rate*thk0
  end subroutine eis_init_temp

  subroutine eis_surftemp(temp,model,time)
    !*FD calculate surface temperature
    use glide_types
    use glimmer_global, only : rk
    implicit none
    type(eis_temp_type)       :: temp  !*FD temperature data
    type(glide_global_type)   :: model !*FD model instance
    real(kind=rk), intent(in) :: time  !*FD current time

    integer i

    call glimmer_ts_step(temp%temp_ts,real(time),temp%tvalue)
    ! spatial temp distrib
    model%climate%artm = temp%tvalue(1)
    if (temp%lat_type.eq.1) then
       model%climate%artm(:,:) = model%climate%artm(:,:) + &
            temp%tvalue(2)*exp(temp%tvalue(3)*(model%climate%lati(:,:)-temp%lat0))
    else
       do i=2,temp%torder
          model%climate%artm(:,:) = model%climate%artm(:,:) + temp%tvalue(i)*model%climate%lati(:,:)**(i-1)
       end do
    end if
    ! vertical temp distrib
    model%climate%artm(:,:) = model%climate%artm(:,:) + temp%lapse_rate * & 
         (model%geometry%usrf(:,:) - model%climate%eus)
  end subroutine eis_surftemp
end module eis_temp
