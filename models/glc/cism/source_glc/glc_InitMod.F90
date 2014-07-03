!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module glc_InitMod

!BOP
! !MODULE: glc_InitMod
! !DESCRIPTION:
!  This module contains the glc initialization method and initializes
!  everything needed by a glc simulation.  Primarily it is a driver
!  that calls individual initialization routines for each glc module.
!
! !USERDOC:
!
! !REFDOC:
!
! !REVISION HISTORY:
!  SVN:$Id: POP_InitMod.F90 808 2006-04-28 17:06:38Z njn01 $
!  Adapted by William Lipscomb from POP_InitMod.F90

! !USES:

   use glc_kinds_mod
   use glc_ErrorMod
   use glc_communicate, only: my_task, master_task
   use glc_broadcast,   only: broadcast_scalar, broadcast_array
   use glc_time_management, only: iyear0, imonth0, iday0, elapsed_days0,  &
                                  iyear,  imonth,  iday,  elapsed_days,   &
                                  ihour,  iminute, isecond, nsteps_total, &
                                  ymd2eday, eday2ymd, runtype
   use glc_constants, only: nml_in, stdout, glc_nec
   use glc_io,        only: glc_io_read_restart_time
   use glc_files,     only: nml_filename
   use glc_exit_mod
   use shr_kind_mod,  only: CL=>SHR_KIND_CL
   use shr_sys_mod, only: shr_sys_flush

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: glc_initialize

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     module variables
!
!-----------------------------------------------------------------------

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: glc_initialize
! !INTERFACE:

 subroutine glc_initialize(errorCode)

! !DESCRIPTION:
!  This routine is the initialization driver that initializes a glc run 
!  by calling individual module initialization routines.
!
! !USERDOC:
!
! !REFDOC:
!
! !REVISION HISTORY:
!  same as module

! !USES:
   use glint_main
   use glint_example_clim

   use glc_global_fields, only: glc_allocate_global, climate, ice_sheet,   &
                                temp,    precip,  orog,    ice_frac,       &
                                tsfc,    qsmb,    topo,                    &
                                gfrac,   gtopo,   grofi,   grofl,  ghflx,  &
				ice_sheet_grid_mask

   use glc_global_grid, only: init_glc_grid, glc_grid
   use glc_override_frac, only: init_glc_frac_overrides
   use glc_constants
   use glc_communicate, only: init_communicate
   use glc_io, only: history_vars
   use glc_time_management, only: init_time1, init_time2, dtt, ihour
   use glimmer_log
   use glc_global_grid, only: glc_landmask
   use glc_route_ice_runoff, only: set_routing
   use shr_file_mod, only : shr_file_getunit, shr_file_freeunit

     !TODO - probably not needed; commented out for now
!!   use glc_global_fields, only: albedo, lats_orog, lons_orog, orog_out
!!   use glc_global_fields, only: time, coverage, cov_orog  ! to be removed?

! !INPUT/OUTPUT PARAMETERS:

   integer (i4), intent(inout) :: &
      errorCode              ! Returns an error code if any init fails

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

  character(fname_length) ::  &
      paramfile        ! Name of the top-level configuration file
 
  character(fname_length) ::  &
      cesm_restart_file  ! Name of the hotstart file to be used for a restart
 
  character(CL) ::  &
      cesm_history_vars  ! Name of the CISM variables to be output in cesm
                         ! history files

  character(CL) :: &
       ice_flux_routing  ! Code for how solid ice should be routed to ocean or sea ice

  ! Scalars which hold information about the global grid --------------
 
  integer (i4) ::  &
      nx,ny,          &! Size of global glc grid
      nxo,nyo          ! Size of global orography grid
 
  integer (i4) ::  &
      i,j              ! Array index counters

  integer (i4) :: &
      nml_error        ! namelist i/o error flag

  integer (i4) :: &
      nhour_glint      ! number of hours since start of complete glint/glimmer run

  logical :: &
      cesm_restart = .false. ! Logical flag to pass to glimmer, telling it to hotstart
                             ! from a CESM restart

  logical :: &
      cism_debug   = .false. ! Logical flag to pass to glimmer, telling it to output extra
                             ! debug diagnostics

  real(dp), dimension(:), allocatable ::  &    
      glint_lats     ,&! lats on glint grid (N to S indexing, instead of S to N as on glc_grid)  
      glint_lons     ,&! lons on glint grid
      glint_latb       ! lat_bound on glint grid

  integer, dimension(:,:), allocatable ::  &    
      glint_landmask   ! landmask on glint grid (N to S indexing)

  integer :: unit      ! fileunit passed to Glint 

  namelist /cism_params/  paramfile, cism_debug, cesm_history_vars, ice_flux_routing
 
!-----------------------------------------------------------------------
!  initialize return flag
!-----------------------------------------------------------------------

   ErrorCode = glc_Success

! TODO - Write version info?
!-----------------------------------------------------------------------
!  write version information to output log after output redirection
!-----------------------------------------------------------------------
!!   if (my_task == master_task) then
!!      write(stdout,blank_fmt)
!!      write(stdout,ndelim_fmt)
!!      write(stdout,blank_fmt)
!!      write(stdout,'(a)') ' GLC version xxx '
!!      write(stdout,blank_fmt)
!!      call shr_sys_flush(stdout)
!!   endif

!-----------------------------------------------------------------------
!
!  compute time step and initialize time-related quantities
!
!-----------------------------------------------------------------------
 
   call init_time1
 
!-----------------------------------------------------------------------
!
!  output delimiter to log file
!
!-----------------------------------------------------------------------
 
   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      call shr_sys_flush (stdout)
   endif
 
!--------------------------------------------------------------------
! Initialize ice sheet model, grid, and coupling.
! The following code is largely based on GLIMMER.
!-----------------------------------------------------------------------

   paramfile  = 'unknown_paramfile'

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=cism_params,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif
   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call exit_glc(sigAbort,'ERROR reading cism_params nml')
   endif

   call broadcast_scalar(paramfile,         master_task)
   call broadcast_scalar(cism_debug,        master_task)
   call broadcast_scalar(cesm_history_vars, master_task)
   history_vars = trim(cesm_history_vars)
   call broadcast_scalar(ice_flux_routing,  master_task)
   call set_routing(ice_flux_routing)

   if (verbose .and. my_task==master_task) then
      write (stdout,*) 'paramfile =   ', paramfile
      write (stdout,*) 'dtt =', dtt
      write (stdout,*) 'Initialize glc_grid'
      call shr_sys_flush(stdout)
   endif

 ! Initialize glc_grid (used for coupling)

   call init_glc_grid   ! glc_landmask is computed here

   ! switch latitude indices for sending info to Glint
   ! latitude is S to N on glc_grid, N to S on glint grid

   nx = glc_grid%nx
   ny = glc_grid%ny

   allocate(glint_lons(nx))    !JW necessary to initialize glint w/o glc_grid
   allocate(glint_lats(ny))

   do j = 1, nx
      glint_lons(j) = glc_grid%lons(j)
   enddo
   do j = 1, ny
      glint_lats(j) = glc_grid%lats(ny-j+1)
   enddo

   allocate(glint_landmask(nx,ny))

   ! Reverse j index for glint_landmask, which assumes increasing index from N to S
   do j = 1, ny
   do i = 1, nx
      glint_landmask(i,j) = glc_landmask(i,ny-j+1)
   enddo
   enddo

   ! set values of climate derived type

   climate%climate_tstep = nint(dtt/3600._r8)   ! convert from sec to integer hours

   if (verbose .and. my_task==master_task) then
      write (stdout,*) 'climate_tstep (hr) =', climate%climate_tstep
      write (stdout,*) 'Set glimmer_unit =', stdout
      write (stdout,*) 'Initialize glint'
   endif

  ! Set glimmer_unit for diagnostic output from Glimmer. (Log file is already open)
!  call open_log(unit=101)

  call set_glimmer_unit(stdout)
 
  ! Allocate global arrays
  call glc_allocate_global(nx, ny, glc_nec)

  ! Initialize global arrays
 
  ! The following are passed to Glint if running if SMB mode.
  tsfc(:,:,:)   = 0._r8
  qsmb(:,:,:)   = 0._r8
  topo(:,:,:)   = 0._r8

  ! These are passed to Glint only if running in PDD mode.
  temp(:,:)     = 0._r8
  precip(:,:)   = 0._r8
  orog(:,:)     = 0._r8

! This array probably is not needed; comment out for now.
!!  orog = real(climate%orog_clim)     

  ! Initialize the ice sheet model

  nhour_glint = 0     ! number of hours glint has run since start of complete simulation
                      ! must be set to correct value if reading from a restart file
 
  call init_glc_frac_overrides()

  ! if this is a continuation run, then set up to read restart file and get the restart time
  if (runtype == 'continue') then
    cesm_restart = .true.
    call glc_io_read_restart_time(nhour_glint, cesm_restart_file)
    call ymd2eday (iyear0, imonth0, iday0, elapsed_days0)
    elapsed_days = elapsed_days0 + nhour_glint/24     
    call eday2ymd(elapsed_days, iyear, imonth, iday)
    ihour = 0
    iminute = 0
    isecond = 0
    nsteps_total = nhour_glint / climate%climate_tstep     
    if (verbose .and. my_task==master_task) then
       write(stdout,*) 'Successfully read restart, nhour_glint =', nhour_glint
       write(stdout,*) 'Initial eday/y/m/d:', elapsed_days0, iyear0, imonth0, iday0
       write(stdout,*) 'eday/y/m/d after restart:', elapsed_days, iyear, imonth, iday
       write(stdout,*) 'nsteps_total =', nsteps_total
       write(stdout,*) 'Initialize glint:'
    endif
  endif

  if (verbose .and. my_task==master_task) then
     write(stdout,*) 'Initialize glint, nhour_glint =', nhour_glint
  endif

  unit = shr_file_getUnit()

  call initialise_glint_gcm(ice_sheet,                            &
                            glint_lats,                           &   ! indexing is N to S for Glint
                            glint_lons,                           &
                            climate%climate_tstep,                &
                            (/paramfile/),                        &
                            daysinyear = climate%days_in_year,    &
                            start_time = nhour_glint,             &
                            glc_nec = glc_nec,                    &
                            gfrac = gfrac,                        &
                            gtopo = gtopo,                        &
                            grofi = grofi,                        &
                            grofl = grofl,                        &
                            ghflx = ghflx,                        &
			    ice_sheet_grid_mask=ice_sheet_grid_mask,&
                            gmask = glint_landmask,               &
                            gcm_restart = cesm_restart,           &
                            gcm_restart_file = cesm_restart_file, &
                            gcm_debug = cism_debug,               &
                            gcm_fileunit = unit)

!TODO - Implement PDD option
 
   call shr_file_freeunit(unit)

! Do the following:
! For each instance, convert ice_sheet%instances(i)%glide_time to hours and compare to nhour_glint.
! If different: Reset params%instances(i)%next_time, params%start_time, params%next_av_start
! Do this here or in initialise_glint?

  ! If restarting (nhour_glint > 0), recompute the year, month, and day
  ! By default, iyear0 = imonth0 = iday0 = 1 (set in namelist file)
  ! Assume that ihour0 = iminute0 = isecond0 = 0
  ! Note that glint does not handle leap years

!jw  Moved this code up to before glint_initialization 
!jw! if this is done before glint initialization, we can simply pass start_time = nint(thour)
!jw! Then set elapsed_days = thour/24 (without call to ymd2eday)
!jw  if (cesm_restart) then
!jw     call ymd2eday (iyear0, imonth0, iday0, elapsed_days0)
!jw     elapsed_days = elapsed_days0 + nhour_glint/24     
!jw     call eday2ymd(elapsed_days, iyear, imonth, iday)
!jw     ihour = 0
!jw     iminute = 0
!jw     isecond = 0
!jw     nsteps_total = nhour_glint / climate%climate_tstep     
!jw  endif

  ! Set the message level (1 is the default - only fatal errors)
  ! N.B. Must do this after initialization
 
  call glimmer_set_msg_level(6)
 
!-----------------------------------------------------------------------
!
!  finish computing time-related quantities after restart info
!  available (including iyear, imonth, and iday)
!
!-----------------------------------------------------------------------
 
   call init_time2
 
!-----------------------------------------------------------------------
!
!  output delimiter to log file
!
!-----------------------------------------------------------------------
 
   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,'(" End of GLC initialization")')
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      call shr_sys_flush (stdout)
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine glc_initialize

!***********************************************************************

 end module glc_InitMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
