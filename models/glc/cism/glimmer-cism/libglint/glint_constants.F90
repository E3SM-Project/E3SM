!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glint_constants.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

module glint_constants

  use glimmer_global
  use glimmer_physcon, only: pi   ! make this constant available to glint modules that use pi

  implicit none

  ! ------------------------------------------------------------
  ! global parameters
  ! ------------------------------------------------------------

  !TODO - Change to dp
!!  real(rk),parameter :: pi=3.141592654          !*FD The value of pi (defined in glimmer_paramets)
  real(rk),parameter :: days2hours = 24.0
  real(rk),parameter :: hours2seconds = 3600.0    !*FD Hours to seconds conversion factor

  integer, parameter :: default_diy = 360                    !*FD Default number of days in year
  integer, parameter :: default_y2h = days2hours*default_diy !*FD Default years to hours conversion

  ! Constants set at run-time

  integer  :: days_in_year = default_diy        !*FD The number of days in a year  
  real(rk) :: years2hours = default_y2h        !*FD Years to hours conversion factor
  real(rk) :: hours2years = 1.0_rk/default_y2h !*FD Hours to years conversion factor

  private :: default_diy, default_y2h

contains

  subroutine glint_set_year_length(daysinyear)

    integer, intent(in) :: daysinyear

    days_in_year = daysinyear
    years2hours = days2hours*days_in_year 
    hours2years = 1.0_rk/years2hours      

  end subroutine glint_set_year_length

end module glint_constants
