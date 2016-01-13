!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   gcm_to_cism_glint.F90 - part of the Community Ice Sheet Model (CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2014
!   CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of CISM.
!
!   CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module gcm_to_cism_glint

  ! This module demonstrates the use of the Glint interface. 
  ! It loads in some example global fields and associated grid data,
  !  initialises the model, and then runs it for a user-prescribed period.
  ! The surface mass balance can be computed either with a PDD scheme
  !  (as in the original Glimmer code), or with a crude scheme 
  !  that imitates SMB input from a climate model.

  use glimmer_global, only: dp
  use glint_main
  use glimmer_log
  use glint_global_interp
  use glint_example_clim
  use glint_commandline
  use glimmer_writestats
!  use glimmer_commandline
  use glimmer_paramets, only: GLC_DEBUG
  use parallel, only: main_task

type gcm_to_cism_type

  ! Program variables -------------------------------------------------------------------

  integer :: which_gcm = 0  ! type of global climate model being used, 0=minimal model, 1=data model, 2=CESM model
  type(glide_global_type) :: glide_model ! ice sheet model used for glide

  type(glint_params)      :: ice_sheet   ! This is the derived type variable that holds all 
                                         ! domains of the ice model
  type(glex_climate)      :: climate     ! Climate parameters and fields

  ! Arrays which hold the global fields used as input to Glint ------------------------

  real(dp),dimension(:,:),pointer :: temp => null()     ! Temperature     (degC)
  real(dp),dimension(:,:),allocatable :: precip    ! Precipitation   (mm/s)
  real(dp),dimension(:,:),allocatable :: orog      ! Orography       (m)

  ! Arrays which hold information about the ice model instances -------------------------

  real(dp),dimension(:,:),allocatable :: coverage ! Coverage map for normal global grid
  real(dp),dimension(:,:),allocatable :: cov_orog ! Coverage map for orography grid

  ! Arrays which hold output from the model ---------------------------------------------
  ! These are all on the normal global grid, except for the orography

  real(dp),dimension(:,:),allocatable :: albedo   ! Fractional albedo
  real(dp),dimension(:,:),allocatable :: orog_out ! Output orography (m)
  real(dp),dimension(:,:),allocatable :: ice_frac ! Ice coverage fraction
  real(dp),dimension(:,:),allocatable :: fw       ! Freshwater output flux (mm/s)
  real(dp),dimension(:,:),allocatable :: fw_in    ! Freshwater input flux (mm/s)

  ! Arrays which hold information about the global grid ---------------------------------

  real(dp),dimension(:),  allocatable :: lats_orog      ! Latitudes of global orography gridpoints
  real(dp),dimension(:),  allocatable :: lons_orog      ! Longitudes of global oropraphy gridpoints

  ! Scalars which hold information about the global grid --------------------------------

  integer :: nx,ny   ! Size of normal global grid
  integer :: nxo,nyo ! Size of global orography grid

  ! Scalar model outputs ----------------------------------------------------------------

  real(dp) :: twin     ! Timestep-integrated input water flux (kg)
  real(dp) :: twout    ! Timestep-integrated output water flux (kg)
  real(dp) :: ice_vol  ! Total ice volume (m^3)

  ! Other variables ---------------------------------------------------------------------

  logical :: out    ! Outputs set flag
  integer :: i,j    ! Array index counters
  integer :: time   ! Current time (hours)
  real(dp):: t1,t2
  integer :: clock,clock_rate

  ! fields passed to and from a GCM 
  ! (useful for testing the GCM subroutines in standalone mode)
  !
  ! Note that, for fields that possess a third dimension, this dimension is the elevation
  ! class. Elevation class goes from 0 to glc_nec, where class 0 represents the bare land
  ! "elevation class".

  real(dp),dimension(:,:,:), allocatable :: qsmb     ! surface mass balance (kg/m^2/s)
  real(dp),dimension(:,:,:), allocatable :: tsfc     ! surface temperature (degC) 
  real(dp),dimension(:,:,:), allocatable :: topo     ! surface elevation (m)

  real(dp),dimension(:,:,:), allocatable :: gfrac    ! fractional glacier area [0,1] 
  real(dp),dimension(:,:,:), allocatable :: gtopo    ! glacier surface elevation (m) 
  real(dp),dimension(:,:,:), allocatable :: ghflx    ! heat flux from glacier interior, positive down (W/m^2)
  real(dp),dimension(:,:),   allocatable :: grofi    ! ice runoff (calving) flux (kg/m^2/s)
  real(dp),dimension(:,:),   allocatable :: grofl    ! ice runoff (liquid) flux (kg/m^2/s)
  real(dp),dimension(:,:),   allocatable :: ice_sheet_grid_mask    ! mask of ice sheet grid coverage
  real(dp),dimension(:,:),   allocatable :: icemask_coupled_fluxes ! mask of ice sheet grid coverage where we are potentially sending non-zero fluxes

  integer :: glc_nec ! , parameter :: glc_nec = 10               ! number of elevation classes

  real(dp),dimension(:), allocatable :: glc_topomax
!          dimension(0:integer(glc_nec)) ::   &
!      glc_topomax = (/ 0.d0,  200.d0,  400.d0,  700.d0, 1000.d0,  1300.d0,   &
!                             1600.d0, 2000.d0, 2500.d0, 3000.d0, 10000.d0 /)  ! upper limit of each class (m)

  logical :: ice_tstep                      ! true if ice timestep was done
  logical :: output_flag                    ! true if outputs have been set

  ! from glint_commandline.F90:
  character(len=5000)         :: commandline_history     !< complete command line
  character(len=fname_length) :: commandline_configname  !< name of the configuration file
  character(len=fname_length) :: commandline_results_fname !< name of results file
  character(len=fname_length) :: commandline_climate_fname !< name of climate configur 

end type gcm_to_cism_type

   logical, parameter :: verbose_glint = .true.  ! set to true for debugging

contains

subroutine g2c_glint_init(g2c)

  ! Initialise glint

  implicit none

  type(gcm_to_cism_type) :: g2c

  integer :: i,j    ! Array index counters

  ! -------------------------------------------------------------------------------------
  ! Executable code starts here - Basic initialisation
  ! -------------------------------------------------------------------------------------
  integer, parameter :: glc_nec = 10               ! number of elevation classes

  real(dp), dimension(0:glc_nec) ::   &
      glc_topomax = (/ 0.d0,  200.d0,  400.d0,  700.d0, 1000.d0,  1300.d0,   &
                             1600.d0, 2000.d0, 2500.d0, 3000.d0, 10000.d0 /)  ! upper limit of each class (m)  

  g2c%glc_nec = glc_nec
  g2c%glc_topomax = glc_topomax
  
  g2c%which_gcm = 1

  call glint_GetCommandline()

  g2c%commandline_history =  commandline_history         !< complete command line
  g2c%commandline_configname =  commandline_configname   !< name of the configuration file
  g2c%commandline_results_fname = commandline_resultsname  !< name of results file
  g2c%commandline_climate_fname = commandline_climatename  !< name of climate configuration 

  call system_clock(g2c%clock,g2c%clock_rate)  
  g2c%t1 = real(g2c%clock,kind=dp)/real(g2c%clock_rate,kind=dp)
  !print *,"g2c%clock, g2c%clock_rate, t1",g2c%clock,g2c%clock_rate,g2c%t1

  if (verbose_glint .and. main_task) print*, 'call glex_clim_init'

  ! Initialise climate

  call glex_clim_init(g2c%climate,g2c%commandline_climate_fname)

  ! Set dimensions of global grids

  call get_grid_dims(g2c%climate%clim_grid,g2c%nx,g2c%ny) ! Normal global grid
  g2c%nxo=200 ; g2c%nyo=100                          ! Example grid used for orographic output

!print *,"g2c% nxo, nyo, nx, ny: ",g2c%nxo,g2c%nyo,g2c%nx,g2c%ny,nxo,nyo

  ! start logging
!  call open_log(unit=101, fname=logname(g2c%commandline_configname))  

  if (verbose_glint .and. main_task) then
     print*, ' '
     print*, 'Starting glint_example:'
     print*, 'climatename = ', trim(g2c%commandline_climate_fname)
     print*, 'configname = ', trim(g2c%commandline_configname)
     print*, 'climate%gcm_smb:', g2c%climate%gcm_smb
     print*, ' '
  endif

  ! Allocate global arrays

  allocate(g2c%temp(g2c%nx,g2c%ny),g2c%precip(g2c%nx,g2c%ny),g2c%orog(g2c%nx,g2c%ny))
  allocate(g2c%coverage(g2c%nx,g2c%ny),g2c%orog_out(g2c%nxo,g2c%nyo),g2c%albedo(g2c%nx,g2c%ny))
!!Check this:
  allocate(g2c%ice_frac(g2c%nx,g2c%ny),g2c%fw(g2c%nx,g2c%ny))
  allocate(g2c%lats_orog(g2c%nyo),g2c%lons_orog(g2c%nxo),g2c%cov_orog(g2c%nxo,g2c%nyo),g2c%fw_in(g2c%nx,g2c%ny))

  ! Initialize global arrays

  g2c%temp     = 0.d0
  g2c%precip   = 0.d0
  g2c%albedo   = 0.d0
  g2c%orog_out = 0.d0
  g2c%orog     = real(g2c%climate%orog_clim,dp)    ! Put orography where it belongs

  ! Allocate and initialize GCM arrays
  
  if (g2c%climate%gcm_smb) then

     ! input from GCM
     allocate(g2c%tsfc(g2c%nx,g2c%ny, 0:g2c%glc_nec))
     allocate(g2c%qsmb(g2c%nx,g2c%ny, 0:g2c%glc_nec))
     allocate(g2c%topo(g2c%nx,g2c%ny, 0:g2c%glc_nec))

     g2c%tsfc(:,:,:)   = 0.d0
     g2c%qsmb(:,:,:)   = 0.d0
     g2c%topo(:,:,:)   = 0.d0

     ! output to GCM
     allocate(g2c%gfrac(g2c%nx,g2c%ny, 0:g2c%glc_nec))
     allocate(g2c%gtopo(g2c%nx,g2c%ny, 0:g2c%glc_nec))
     allocate(g2c%ghflx(g2c%nx,g2c%ny, 0:g2c%glc_nec))
     allocate(g2c%grofi(g2c%nx,g2c%ny))
     allocate(g2c%grofl(g2c%nx,g2c%ny))
     allocate(g2c%ice_sheet_grid_mask(g2c%nx,g2c%ny))
     allocate(g2c%icemask_coupled_fluxes(g2c%nx,g2c%ny))

     g2c%gfrac(:,:,:) = 0.d0
     g2c%gtopo(:,:,:) = 0.d0
     g2c%ghflx(:,:,:) = 0.d0
     g2c%grofi(:,:)   = 0.d0
     g2c%grofl(:,:)   = 0.d0
     g2c%ice_sheet_grid_mask(:,:)    = 0.d0
     g2c%icemask_coupled_fluxes(:,:) = 0.d0

  endif

  ! Set up global grids ----------------------------------------------------------------

  ! Calculate example orographic latitudes

  do j=1,g2c%nyo
     g2c%lats_orog(j) = -(180.d0/g2c%nyo)*j + 90.d0 + (90.d0/g2c%nyo)
  enddo

  ! Calculate example orographic longitudes

  do i=1,g2c%nxo
     g2c%lons_orog(i) = (360.d0/g2c%nxo)*i - (180.d0/g2c%nxo)
  enddo

  ! Initialise the ice model

  if (g2c%climate%gcm_smb) then   ! act as if we are receiving the SMB from a GCM

     call initialise_glint_gcm(g2c%ice_sheet,                       &
                               g2c%climate%clim_grid%lats,          &
                               g2c%climate%clim_grid%lons,          &
                               g2c%climate%climate_tstep,           &
                               (/g2c%commandline_configname/),      &
                               daysinyear=g2c%climate%days_in_year, &
                               glc_nec = g2c%glc_nec,               &
                               gfrac = g2c%gfrac,                   &
                               gtopo = g2c%gtopo,                   &
                               grofi = g2c%grofi,                   &
                               grofl = g2c%grofl,                   &
                               ghflx = g2c%ghflx,                   &
                               ice_sheet_grid_mask    = g2c%ice_sheet_grid_mask, &
                               icemask_coupled_fluxes = g2c%icemask_coupled_fluxes)

  else   ! standard Glint initialization

     call initialise_glint(g2c%ice_sheet, &
                           g2c%climate%clim_grid%lats, &
                           g2c%climate%clim_grid%lons, &
                           g2c%climate%climate_tstep, &
                           (/g2c%commandline_configname/), &
                           orog=g2c%orog_out, &
                           albedo=g2c%albedo, &
                           ice_frac=g2c%ice_frac, &
                           orog_longs=g2c%lons_orog, &
                           orog_lats=g2c%lats_orog, &
                           daysinyear=g2c%climate%days_in_year)

  endif  ! gcm_smb

  ! Set the message level (1 is the default - only fatal errors)
  ! N.B. Must do this after initialisation

  call glimmer_set_msg_level(6)

  ! Get coverage maps for the ice model instances

  if (g2c%climate%gcm_smb) then     ! not using cov_orog
     if (glint_coverage_map(g2c%ice_sheet, g2c%coverage) .ne. 0) then
        call write_log('Unable to get coverage maps',GM_FATAL,__FILE__,__LINE__)
        stop
     endif
  else
     if (glint_coverage_map(g2c%ice_sheet, g2c%coverage, g2c%cov_orog) .ne. 0) then
        call write_log('Unable to get coverage maps',GM_FATAL,__FILE__,__LINE__)
        stop
     endif
  endif

  g2c%time = g2c%climate%climate_tstep     ! time in integer hours
                                 
!  if (main_task) print*, 'Done in g2c_glint_init'

end subroutine g2c_glint_init


subroutine g2c_glint_run(g2c)

  type(gcm_to_cism_type) :: g2c

  ! Do timesteps ---------------------------------------------------------------------------

  !TODO - Timestepping as in simple_glide?  Initialize with time = 0, then update time right after 'do'
  !       This would require changing some time logic inside the Glint subroutines.

!  g2c%time = g2c%climate%climate_tstep     ! time in integer hours

!  do

    !      The SMB is computed crudely for now, just to test the GCM interfaces.
    !      At some point we could read in a realistic SMB as in CESM TG runs.

     ! get current temp and precip fields

     call example_climate(g2c%climate, g2c%precip, g2c%temp, real(g2c%time,dp))

     if (g2c%climate%gcm_smb) then   ! act as if we are receiving the SMB from a GCM

        !TODO - For some reason, the gcm code is much slower than the pdd code.
        !       Figure out why.

        ! call a simple subroutine to estimate qsmb and tsfc in different elevation classes

        call compute_gcm_smb(g2c%temp,        g2c%precip,   &
                             g2c%orog,                  &
                             g2c%qsmb,        g2c%tsfc,     &
                             g2c%topo,                  &
                             g2c%glc_nec,     g2c%glc_topomax)

        call glint_gcm (g2c%ice_sheet,        g2c%time,            &
                        g2c%qsmb,             g2c%tsfc,            &
                        g2c%topo,                              &
                        output_flag = g2c%output_flag,         &
                        ice_tstep = g2c%ice_tstep,             & 
                        gfrac = g2c%gfrac,    &
                        gtopo = g2c%gtopo,    &
                        grofi = g2c%grofi,    &
                        grofl = g2c%grofl,    &
                        ghflx = g2c%ghflx,    &
                        ice_sheet_grid_mask    = g2c%ice_sheet_grid_mask, &
                        icemask_coupled_fluxes = g2c%icemask_coupled_fluxes)

     else    ! standard Glint timestepping 

        call glint(g2c%ice_sheet,  g2c%time,   g2c%temp,      g2c%precip,     g2c%orog,            &
                   orog_out=g2c%orog_out,   albedo=g2c%albedo,         output_flag=g2c%out, &
                   ice_frac=g2c%ice_frac,   water_out=g2c%fw,          water_in=g2c%fw_in,  &
                   total_water_in=g2c%twin, total_water_out=g2c%twout, ice_volume=g2c%ice_vol) 

     endif   ! gcm_smb

     !g2c%time = g2c%time + g2c%climate%climate_tstep
     !  if (g2c%time > g2c%climate%total_years*g2c%climate%hours_in_year) exit

!  end do  ! main timestep loop

  if (GLC_DEBUG) then
     ! Print time so as to have something to watch while the code runs
     if (mod(real(g2c%time,dp),8760.d0) < 0.01) print*, 'time (yr) =', real(g2c%time,dp)/8760.d0
  end if
end subroutine g2c_glint_run


subroutine g2c_glint_climate_time_step(g2c)
  type(gcm_to_cism_type) :: g2c

  g2c%time = g2c%time + g2c%climate%climate_tstep
end subroutine g2c_glint_climate_time_step

subroutine g2c_glint_check_finished(g2c,finished)
   type(gcm_to_cism_type) :: g2c
   logical :: finished

   if (g2c%time > g2c%climate%total_years*g2c%climate%hours_in_year) then
      finished = .true.
   else
      finished = .false.
   endif

end subroutine g2c_glint_check_finished


subroutine g2c_glint_end(g2c)
  type(gcm_to_cism_type) :: g2c

  ! Finalise/tidy up everything -----------------------------------------------------------

  call end_glint(g2c%ice_sheet)
  call system_clock(g2c%clock,g2c%clock_rate)
  t2 = real(g2c%clock,kind=dp)/real(g2c%clock_rate,kind=dp)
  call glimmer_write_stats(g2c%commandline_results_fname,g2c%commandline_configname,g2c%t2-g2c%t1)

  ! 101 format(e12.5)

end subroutine g2c_glint_end

!---------------------------------------------------------------------------------


end module gcm_to_cism_glint
