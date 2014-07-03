!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glex_ebm.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

program glint_ebm_ex

!*FD This program demonstrates the use of GLIMMER. It loads in
!*FD some example global fields and associated grid data,
!*FD Initialises the model, and then runs it for 1000 years.

  use glint_main
  use glimmer_log
  use glint_global_interp
  use glint_ebm_ex_clim
  use glint_commandline
  use glimmer_writestats
  implicit none

  ! Program variables -------------------------------------------------------------------

  type(glint_params)      :: ice_sheet   ! This is the derived type variable that holds all 
                                         ! domains of the ice model
  type(glex_ebm_climate)      :: climate     ! Climate parameters and fields

  ! Arrays which hold the global fields used as input to GLIMMER ------------------------

  real(rk),dimension(:,:),allocatable :: temp      ! Temperature     (degC)
  real(rk),dimension(:,:),allocatable :: precip    ! Precipitation   (mm/s)
  real(rk),dimension(:,:),allocatable :: zonwind   ! Zonal wind      (m/s)
  real(rk),dimension(:,:),allocatable :: merwind   ! Meridional wind (m/s)
  real(rk),dimension(:,:),allocatable :: orog      ! Orography       (m)
  real(rk),dimension(:,:),allocatable :: humid     ! Surface Humidity (%)
  real(rk),dimension(:,:),allocatable :: lwdown    ! Downwelling LW  (W/m^2)
  real(rk),dimension(:,:),allocatable :: swdown    ! Downwelling SW  (W/m^2)
  real(rk),dimension(:,:),allocatable :: airpress  ! Surface air pressure (Pa)

  ! Arrays which hold information about the ice model instances -------------------------

  real(rk),dimension(:,:),allocatable :: coverage ! Coverage map for normal global grid
  real(rk),dimension(:,:),allocatable :: cov_orog ! Coverage map for orography grid

  ! Arrays which hold output from the model ---------------------------------------------
  ! These are all on the normal global grid, except for the orography

  real(rk),dimension(:,:),allocatable :: albedo   ! Fractional albedo
  real(rk),dimension(:,:),allocatable :: orog_out ! Output orography (m)
  real(rk),dimension(:,:),allocatable :: ice_frac ! Ice coverage fraction
  real(rk),dimension(:,:),allocatable :: fw       ! Freshwater output flux (mm/s)
  real(rk),dimension(:,:),allocatable :: fw_in    ! Freshwater input flux (mm/s)

  ! Arrays which hold information about the global grid ---------------------------------

  real(rk),dimension(:),  allocatable :: lats_orog      ! Latitudes of global orography gridpoints
  real(rk),dimension(:),  allocatable :: lons_orog      ! Longitudes of global oropraphy gridpoints

  ! Scalars which hold information about the global grid --------------------------------

  integer :: nx,ny   ! Size of normal global grid
  integer :: nxo,nyo ! Size of global orography grid

  ! Scalar model outputs ----------------------------------------------------------------

  real(rk) :: twin     ! Timestep-integrated input water flux (kg)
  real(rk) :: twout    ! Timestep-integrated output water flux (kg)
  real(rk) :: ice_vol  ! Total ice volume (m^3)
  
  ! Other variables ---------------------------------------------------------------------

  logical :: out    ! Outputs set flag
  integer :: i,j    ! Array index counters
  integer :: time   ! Current time (hours)
  real(kind=dp) t1,t2
  integer clock,clock_rate
  integer :: ierr

  ! -------------------------------------------------------------------------------------
  ! Executable code starts here - Basic initialisation
  ! -------------------------------------------------------------------------------------

  call glint_GetCommandline()

  ! start timing
  call system_clock(clock,clock_rate)
  t1 = real(clock,kind=dp)/real(clock_rate,kind=dp)

  ! Initialise climate

  call glex_ebm_clim_init(climate,commandline_climatename)

  ! Set dimensions of global grids

  call get_grid_dims(climate%grid,nx,ny) ! Normal global grid
  nxo=200 ; nyo=100                      ! Example grid used for orographic output

  ! start logging
  call open_log(unit=101)  

  ! Allocate arrays appropriately

  allocate(temp(nx,ny),precip(nx,ny),zonwind(nx,ny))
  allocate(merwind(nx,ny),orog(nx,ny))
  allocate(coverage(nx,ny),orog_out(nxo,nyo),albedo(nx,ny),ice_frac(nx,ny),fw(nx,ny))
  allocate(lats_orog(nyo),lons_orog(nxo),cov_orog(nxo,nyo),fw_in(nx,ny))
  allocate(humid(nx,ny),lwdown(nx,ny),swdown(nx,ny),airpress(nx,ny))

  ! Initialise array contents

  temp=0.0
  precip=0.0
  zonwind=0.0
  merwind=0.0
  humid=0.0
  lwdown=0.0
  swdown=0.0
  airpress=0.0
  albedo=0.0
  orog_out=0.0
  orog=real(climate%orog_clim)                    ! Put orography where it belongs

  ! Set up global grids ----------------------------------------------------------------

  ! Calculate example orographic latitudes

  do j=1,nyo
    lats_orog(j)=-(180.0/nyo)*j+90.0+(90.0/nyo)
  enddo

  ! Calculate example orographic longitudes

  do i=1,nxo
    lons_orog(i)=(360.0/nxo)*i-(180.0/nxo)
  enddo

  ! Initialise the ice model

  call initialise_glint(ice_sheet, &
       climate%grid%lats, &
       climate%grid%lons, &
       climate%climate_tstep, &
       (/commandline_configname/), &
       orog=orog_out, &
       ice_frac=ice_frac, &
       albedo=albedo, &
       orog_longs=lons_orog, &
       orog_lats=lats_orog, &
       daysinyear=365)

  ! Set the message level (1 is the default - only fatal errors)
  ! N.B. Must do this after initialisation

  call glimmer_set_msg_level(6)

  ! Get coverage maps for the ice model instances

  if (glint_coverage_map(ice_sheet,coverage,cov_orog).ne.0) then
    call write_log('Unable to get coverage maps',GM_FATAL,__FILE__,__LINE__)
    stop
  endif

  ! Do timesteps ---------------------------------------------------------------------------

  time=climate%climate_tstep

  do
     call ebm_ex_climate(climate,precip,temp,zonwind,merwind,humid,lwdown,swdown,airpress,real(time,rk))
     call glint(ice_sheet,time,temp,precip,orog, &
          zonwind=zonwind,     merwind=merwind,       humid=humid, &
          lwdown=lwdown,       swdown=swdown,         airpress=airpress, &
          orog_out=orog_out,   albedo=albedo,         output_flag=out, &
          ice_frac=ice_frac,   water_out=fw,          water_in=fw_in, &
          total_water_in=twin, total_water_out=twout, ice_volume=ice_vol) 
     time=time+climate%climate_tstep
     if (time>climate%total_years*climate%hours_in_year) exit
  end do

  ! Finalise/tidy up everything ------------------------------------------------------------

  call end_glint(ice_sheet)
  call system_clock(clock,clock_rate)
  t2 = real(clock,kind=dp)/real(clock_rate,kind=dp)
  call glimmer_write_stats(commandline_resultsname,commandline_configname,t2-t1)

100 format(f9.5)
101 format(e12.5)

end program glint_ebm_ex
