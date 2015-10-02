!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   eismint3_types.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

module eismint3_types

  use glimmer_global, only: dp
  use glint_pdd

  implicit none

  type eismint3_climate
     real(dp),dimension(:,:),pointer :: prcp !*FD Annual accumulation
     real(dp),dimension(:,:),pointer :: acab !*FD Mass-balance
     real(dp),dimension(:,:),pointer :: artm !*FD Surface temp
     real(dp),dimension(:,:),pointer :: arng !*FD Surface temp half-range
     real(dp),dimension(:,:),pointer :: usrf !*FD Surface elevation
     real(dp),dimension(:,:),pointer :: ablt !*FD Ablation
     real(dp),dimension(:,:),pointer :: presusurf !*FD Present-day upper-surface elevation
     real(dp),dimension(:,:),pointer :: presartm  !*FD Present-day surface temperature
     real(dp),dimension(:,:),pointer :: presprcp  !*FD Present-day precip (water-equivalent)
     logical,dimension(:,:),pointer :: landsea !*FD Land-sea mask
     type(glint_pdd_params) :: pdd_scheme
     integer :: loadthk = 0          !*FD Load thickness from file or start from scratch
     real(dp) :: pfac=1.0533d0       !*FD Precip enhancement factor (default is supposed EISMINT value)
     real(dp) :: temp_perturb = 0.d0 !*FD Climate temperature perturbation
  end type eismint3_climate

end module eismint3_types
