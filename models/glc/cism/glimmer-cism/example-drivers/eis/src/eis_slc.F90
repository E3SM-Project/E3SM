!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   eis_slc.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

module eis_slc
  !*FD climate forcing similar to the old Edinburgh Ice Sheet model
  !*FD Magnus Hagdorn, June 2004
  !*FD this modules handles the sea level component

  use glimmer_ts
  use glimmer_global, only : fname_length

  type eis_slc_type
     !*FD parameters for the EIS sea level forcing
     character(len=fname_length) :: fname=''     !*FD name of file containing ELA ts
     type(glimmer_tseries) :: slc_ts             !*FD ELA time series 
  end type eis_slc_type
  
contains  
  subroutine eis_slc_config(config,slc)
    !*FD get SLC configuration from config file
    use glimmer_config
    use glimmer_filenames, only : filenames_inputname
    implicit none
    type(eis_slc_type)           :: slc     !*FD slc data
    type(ConfigSection), pointer :: config  !*FD structure holding sections of configuration file   
    ! local variables
    type(ConfigSection), pointer :: section

    slc%fname=''

    call GetSection(config,section,'EIS SLC')
    if (associated(section)) then
       call GetValue(section,'slc_file',slc%fname)
       if (trim(slc%fname).ne.'') then
          slc%fname = trim(filenames_inputname(slc%fname))
       end if
    end if
  end subroutine eis_slc_config

  subroutine eis_slc_printconfig(slc)
    !*FD print configuration to log
    use glimmer_log
    use glimmer_global, only : msg_length
    implicit none
    type(eis_slc_type)           :: slc     !*FD slc data
    ! local variables
    character(len=msg_length) :: message

    call write_log('EIS SLC')
    call write_log('-------')
    write(message,*) 'SLC file: ',trim(slc%fname)
    call write_log(message)
    call write_log('')
  end subroutine eis_slc_printconfig

  subroutine eis_init_slc(slc)
    !*FD initialise SLC forcing
    use glimmer_paramets, only: thk0
    implicit none
    type(eis_slc_type)           :: slc     !*FD slc data
    
    call glimmer_read_ts(slc%slc_ts,slc%fname)
    ! scale parameters
    slc%slc_ts%values = slc%slc_ts%values/thk0
  end subroutine eis_init_slc

  subroutine eis_eus(slc,model,time)
    !*FD calculate mass balance
    use glide_types
    use glimmer_global, only : dp
    implicit none
    type(eis_slc_type)        :: slc   !*FD slc data
    type(glide_global_type)   :: model !*FD model instance
    real(dp), intent(in)      :: time  !*FD current time

    call glimmer_ts_linear(slc%slc_ts, time, model%climate%eus)
  end subroutine eis_eus
end module eis_slc
