!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glint_example.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

program glint_example

  !*FD This program demonstrates the use of GLIMMER. It loads in
  !*FD some example global fields and associated grid data,
  !*FD Initialises the model, and then runs it for 1000 years.

  use glint_main
  use glimmer_log
  use glint_global_interp
  use glint_example_clim
  use glint_commandline
  use glimmer_writestats
  use glimmer_paramets, only: GLC_DEBUG

!WHL - debug
  use glimmer_physcon, only: scyr
  use glint_type, only: iglint_global, jglint_global

  implicit none

  ! Program variables -------------------------------------------------------------------

  type(glint_params)      :: ice_sheet   ! This is the derived type variable that holds all 
                                         ! domains of the ice model
  type(glex_climate)      :: climate     ! Climate parameters and fields

  ! Arrays which hold the global fields used as input to GLIMMER ------------------------

  real(rk),dimension(:,:),allocatable :: temp      ! Temperature     (degC)
  real(rk),dimension(:,:),allocatable :: precip    ! Precipitation   (mm/s)
  real(rk),dimension(:,:),allocatable :: orog      ! Orography       (m)

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

  ! fields passed to and from a GCM 
  ! (useful for testing the GCM subroutines in standalone mode)

  real(rk),dimension(:,:,:), allocatable :: qsmb     ! surface mass balance (kg/m^2/s)
  real(rk),dimension(:,:,:), allocatable :: tsfc     ! surface temperature (degC) 
  real(rk),dimension(:,:,:), allocatable :: topo     ! surface elevation (m)

  real(rk),dimension(:,:,:), allocatable :: gfrac    ! fractional glacier area [0,1] 
  real(rk),dimension(:,:,:), allocatable :: gtopo    ! glacier surface elevation (m) 
  real(rk),dimension(:,:,:), allocatable :: grofi    ! ice runoff (calving) flux (kg/m^2/s)
  real(rk),dimension(:,:,:), allocatable :: grofl    ! ice runoff (liquid) flux (kg/m^2/s)
  real(rk),dimension(:,:,:), allocatable :: ghflx    ! heat flux from glacier interior, positive down (W/m^2)

  integer, parameter :: glc_nec = 10               ! number of elevation classes

  real(rk), dimension(0:glc_nec) ::   &
      glc_topomax = (/ 0.d0,  200.d0,  400.d0,  700.d0, 1000.d0,  1300.d0,   &
                             1600.d0, 2000.d0, 2500.d0, 3000.d0, 10000.d0 /)  ! upper limit of each class (m)

  logical :: ice_tstep                      ! true if ice timestep was done
  logical :: output_flag                    ! true if outputs have been set

  integer :: ig, jg, k

  ! -------------------------------------------------------------------------------------
  ! Executable code starts here - Basic initialisation
  ! -------------------------------------------------------------------------------------

  call glint_GetCommandline()

  ! start timing
  call system_clock(clock,clock_rate)
  t1 = real(clock,kind=dp)/real(clock_rate,kind=dp)

  ! Initialise climate

  call glex_clim_init(climate,commandline_climatename)

  ! Set dimensions of global grids

  call get_grid_dims(climate%clim_grid,nx,ny) ! Normal global grid
  nxo=200 ; nyo=100                          ! Example grid used for orographic output

  ! start logging
  call open_log(unit=101, fname=logname(commandline_configname))  

  print*, 'Starting program glint_example'
  print*, 'climatename = ', trim(commandline_climatename)
  print*, 'configname = ', trim(commandline_configname)
  print*, 'climate%gcm_smb =', climate%gcm_smb
  print*, ' '

  ! Allocate global arrays

  allocate(temp(nx,ny),precip(nx,ny),orog(nx,ny))
  allocate(coverage(nx,ny),orog_out(nxo,nyo),albedo(nx,ny),ice_frac(nx,ny),fw(nx,ny))
  allocate(lats_orog(nyo),lons_orog(nxo),cov_orog(nxo,nyo),fw_in(nx,ny))

  !TODO - Change to dp
  ! Initialize global arrays

  temp     = 0.0
  precip   = 0.0
  albedo   = 0.0
  orog_out = 0.0
  orog     = real(climate%orog_clim)    ! Put orography where it belongs

  ! Allocate and initialize GCM arrays
  
  if (climate%gcm_smb) then

     ! input from GCM
     allocate(tsfc(nx,ny,glc_nec))
     allocate(qsmb(nx,ny,glc_nec))
     allocate(topo(nx,ny,glc_nec))

     tsfc(:,:,:)   = 0.d0
     qsmb(:,:,:)   = 0.d0
     topo(:,:,:)   = 0.d0

     ! output to GCM
     allocate(gfrac(nx,ny,glc_nec))
     allocate(gtopo(nx,ny,glc_nec))
     allocate(grofi(nx,ny,glc_nec))
     allocate(grofl(nx,ny,glc_nec))
     allocate(ghflx(nx,ny,glc_nec))

     gfrac(:,:,:) = 0.d0
     gtopo(:,:,:) = 0.d0
     grofi(:,:,:) = 0.d0
     grofl(:,:,:) = 0.d0
     ghflx(:,:,:) = 0.d0

  endif

  ! Set up global grids ----------------------------------------------------------------

  ! Calculate example orographic latitudes

  do j=1,nyo
     lats_orog(j) = -(180.0/nyo)*j + 90.0 + (90.0/nyo)
  enddo

  ! Calculate example orographic longitudes

  do i=1,nxo
     lons_orog(i) = (360.0/nxo)*i - (180.0/nxo)
  enddo

  ! Initialise the ice model

  if (climate%gcm_smb) then   ! act as if we are receiving the SMB from a GCM

     call initialise_glint_gcm(ice_sheet,                       &
                               climate%clim_grid%lats,          &
                               climate%clim_grid%lons,          &
                               climate%climate_tstep,           &
                               (/commandline_configname/),      &
                               daysinyear=climate%days_in_year, &
                               glc_nec = glc_nec,               &
                               gfrac = gfrac,                   &
                               gtopo = gtopo,                   &
                               grofi = grofi,                   &
                               grofl = grofl,                   &
                               ghflx = ghflx)

  else   ! standard Glint initialization

     call initialise_glint(ice_sheet, &
                           climate%clim_grid%lats, &
                           climate%clim_grid%lons, &
                           climate%climate_tstep, &
                           (/commandline_configname/), &
                           orog=orog_out, &
                           albedo=albedo, &
                           ice_frac=ice_frac, &
                           orog_longs=lons_orog, &
                           orog_lats=lats_orog, &
                           daysinyear=climate%days_in_year)

  endif  ! gcm_smb

  ! Set the message level (1 is the default - only fatal errors)
  ! N.B. Must do this after initialisation

  call glimmer_set_msg_level(6)

  ! Get coverage maps for the ice model instances

  if (climate%gcm_smb) then     ! not using cov_orog
     if (glint_coverage_map(ice_sheet, coverage) .ne. 0) then
        call write_log('Unable to get coverage maps',GM_FATAL,__FILE__,__LINE__)
        stop
     endif
  else
     if (glint_coverage_map(ice_sheet, coverage, cov_orog) .ne. 0) then
        call write_log('Unable to get coverage maps',GM_FATAL,__FILE__,__LINE__)
        stop
     endif
  endif

  ! Do timesteps ---------------------------------------------------------------------------

  !TODO - Timestepping as in simple_glide?  Initialize with time = 0, then update time right after 'do'
  !       This would require changing some time logic inside the Glint subroutines.

  time = climate%climate_tstep     ! time in integer hours

  do

    !      The SMB is computed crudely for now, just to test the GCM interfaces.
    !      At some point we could read in a realistic SMB as in CESM TG runs.

     ! get current temp and precip fields

     call example_climate(climate, precip, temp, real(time,rk))

     if (climate%gcm_smb) then   ! act as if we are receiving the SMB from a GCM
                                  
        ! call a simple subroutine to estimate qsmb and tsfc in different elevation classes

        call compute_gcm_smb(temp,        precip,   &
                             orog,                  &
                             qsmb,        tsfc,     &
                             topo,                  &
                             glc_nec,     glc_topomax)

        !WHL - debug
        !     ig = iglint_global  ! in glint_type
        !     jg = jglint_global 
        !     print*, ' '
        !     print*, 'Global i, j, time (days):', ig, jg, real(time)/24.d0 
        !     print*, ' '
        !     print*, 'orog (m), temp (C), prcp (m/yr):'
        !     print*, orog(ig,jg), temp(ig,jg), precip(ig,jg)*scyr
        !     print*, ' '
        !     print*, 'topo (m), tsfc (C), qsmb (m/yr):'
        !     do k = 1, glc_nec
        !        print*, topo(ig,jg,k), tsfc(ig,jg,k), qsmb(ig,jg,k)*scyr
        !     enddo
        !

        ! call glint

        call glint_gcm (ice_sheet,        time,            &
                        qsmb,             tsfc,            &
                        topo,                              &
                        output_flag = output_flag,         &
                        ice_tstep = ice_tstep,             & 
                        gfrac = gfrac,    gtopo = gtopo,   &
                        grofi = grofi,    grofl = grofl,   &
                        ghflx = ghflx)

     else    ! standard Glint timestepping 

        call glint(ice_sheet,   time,   temp,      precip,     orog,            &
                   orog_out=orog_out,   albedo=albedo,         output_flag=out, &
                   ice_frac=ice_frac,   water_out=fw,          water_in=fw_in,  &
                   total_water_in=twin, total_water_out=twout, ice_volume=ice_vol) 

     endif   ! gcm_smb

     time = time + climate%climate_tstep
     if (time > climate%total_years*climate%hours_in_year) exit

  end do  ! main timestep loop

  if (GLC_DEBUG) then
     ! Print time so as to have something to watch while the code runs
     if (mod(real(time),8760.) < 0.01)   &
         print*, 'time (yr) =', real(time)/8760.
  end if

  ! Finalise/tidy up everything ------------------------------------------------------------

  call end_glint(ice_sheet)
  call system_clock(clock,clock_rate)
  t2 = real(clock,kind=dp)/real(clock_rate,kind=dp)
  call glimmer_write_stats(commandline_resultsname,commandline_configname,t2-t1)

100 format(f9.5)
101 format(e12.5)

!---------------------------------------------------------------------------------

end program glint_example

!---------------------------------------------------------------------------------
