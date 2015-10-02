!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module restart

!BOP
! !MODULE: restart
! !DESCRIPTION:
!  This module contains routins for reading and writing data necessary
!  for restarting a POP simulation.
!
! !REVISION HISTORY:
!  SVN:$Id: restart.F90 59359 2014-04-20 22:10:41Z mlevy@ucar.edu $
!
! !USES:

   use POP_KindsMod
   use POP_ErrorMod
   use POP_IOUnitsMod
   use POP_FieldMod
   use POP_GridHorzMod
   use POP_HaloMod

   use domain_size
   use domain
   use constants, only: char_blank, field_loc_NEcorner, field_type_vector,  &
       field_loc_center, field_type_scalar, blank_fmt, c0, grav
   use blocks, only: nx_block, ny_block, block, get_block
   use prognostic, only: UBTROP, VBTROP, PSURF, GRADPX, GRADPY, UVEL, VVEL, &
       PGUESS, TRACER, nt, nx_global, ny_global, km, curtime, newtime, oldtime,      &
       tracer_d
   use broadcast, only: broadcast_scalar 
   use communicate, only: my_task, master_task
   use operators, only: div
   use grid, only: sfc_layer_type, sfc_layer_varthick, CALCU, CALCT, KMU,   &
       KMT, HU, TAREA_R
   use io, only: data_set
   use io_types, only: io_field_desc, datafile, io_dim, luse_pointer_files, &
       pointer_filename, stdout, construct_io_field, construct_file,        &
       rec_type_dbl, construct_io_dim, nml_in, nml_filename, get_unit,      &
       release_unit, destroy_file, add_attrib_file, destroy_io_field,       &
       extract_attrib_file
   use time_management
   use ice, only: tlast_ice, liceform, AQICE, FW_FREEZE, QFLUX
   use forcing_fields, only: FW_OLD
   use forcing_ap, only: ap_interp_last
   use forcing_ws, only: ws_interp_last
   use forcing_shf, only: shf_interp_last
   use forcing_sfwf, only: sfwf_interp_last, sum_precip, precip_fact,       &
       ssh_initial, sal_initial
   use forcing_pt_interior, only: pt_interior_interp_last
   use forcing_s_interior, only: s_interior_interp_last
   use exit_mod, only: sigAbort, exit_pop, flushm
   use registry
   use passive_tracers, only: write_restart_passive_tracers
   use overflows
   use overflow_type

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_restart,  &
             write_restart, &
             read_restart

! !PUBLIC DATA MEMBERS:
   public :: restart_fmt,   &
             read_restart_filename, &
             lrestart_write

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  module variables
!
!-----------------------------------------------------------------------

   character (POP_charLength) :: &
      restart_outfile       ! restart output filename root

   character (POP_charLength) :: &
      restart_fmt           ! format (bin or nc) of output restart

   character (POP_charLength) ::  &
      read_restart_filename = 'undefined' ! file name for restart file

   character (POP_charLength) ::  &
      exit_string           = 'undefined' ! error-exit string

   logical (POP_logical) ::   &
      pressure_correction, &! fix pressure for exact restart
      lrestart_on,         &! flag to turn restarts on/off
      leven_odd_on,        &! flag to turn even_odd restarts on/off
      lrestart_write        ! flag to determine whether restart is written


   integer (POP_i4) ::  &
      even_odd_freq,      &! even/odd restart files every freq steps
      last_even_odd,      &! last even/odd dump
      restart_flag,       &! time flag id for restarts
      evenodd_flag,       &! time flag id for even-odd restarts
      out_stop_now,       &! time flag id for stop_now flag
      restart_cpl_ts,     &! time flag id for coupled_ts time flag
      restart_freq_iopt,  &! restart frequency option
      restart_freq,       &! restart frequency
      restart_start_iopt, &! start after option 
      restart_start        ! start regular restart writes after restart_start

   integer (POP_i4), parameter :: &
      even = 0,           &! integer for which even/odd dump
      odd  = 1

   !*** field descriptors for all output fields

   type (io_field_desc) ::    &
      UBTROP_CUR, UBTROP_OLD, &! barotropic U at current, old times
      VBTROP_CUR, VBTROP_OLD, &! barotropic U at current, old times
      PSURF_CUR,  PSURF_OLD,  &! surface press at current, old times
      GRADPX_CUR, GRADPX_OLD, &! sfc press gradient in x at both times
      GRADPY_CUR, GRADPY_OLD, &! sfc press gradient in y at both times
      PGUESSd,                &! guess for next surface pressure
      FW_OLDd,                &! freshwater input at old time
      FW_FREEZEd,             &! water flux at T points due to frazil ice formation
      AQICEd,                 &! accumulated ice melt/freeze
      QFLUXd,                 &! internal ocn heat flux due to ice formation
      UVEL_CUR, UVEL_OLD,     &! U at current, old times
      VVEL_CUR, VVEL_OLD       ! V at current, old times

   type (io_field_desc), dimension(nt) :: &
      TRACER_CUR, TRACER_OLD    ! tracers at current, old times

!-----------------------------------------------------------------------
!   ccsm coupling variable
!-----------------------------------------------------------------------
      integer (POP_i4) ::  &
         cpl_write_restart     ! flag id for restart-file signal from cpl


!-----------------------------------------------------------------------
!
!  scalar data to be written/read from restart file
!
!     runid,
!     iyear,  imonth,  iday,  ihour,  iminute,  isecond
!     iyear0, imonth0, iday0, ihour0, iminute0, isecond0
!     dtt,    iday_of_year,   iday_of_year_last
!     elapsed_days, elapsed_months, elapsed_years
!     elapsed_days_this_year
!     seconds_this_day,  seconds_this_day_next
!     seconds_this_year, seconds_this_year_next
!     nsteps_total
!     eod, eod_last, eom, eom_last, eom_next, eoy, eoy_last
!     midnight_last, adjust_year_next, newday, newhour
!     leapyear, days_in_year, days_in_prior_year
!     seconds_in_year, hours_in_year
!     tlast_ice
!     lcoupled_ts
!     shf_interp_last, sfwf_interp_last, ws_interp_last
!     ap_interp_last, pt_interior_interp_last
!     s_interior_interp_last
!     sal_initial, sum_precip, precip_fact, ssh_initial
!
!-----------------------------------------------------------------------

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: read_restart
! !INTERFACE:

 subroutine read_restart(in_filename,lccsm_branch,lccsm_hybrid, &
                         in_restart_fmt, errorCode)

! !DESCRIPTION:
!  This routine reads restart data from a file.
!
!  Prognostic fields read are:
!     UBTROP,VBTROP  : barotropic velocities
!     PSURF          : surface pressure
!     GRADPX,GRADPY  : surface pressure gradient
!     PGUESS         : next guess for pressure
!     UVEL,VVEL      : 3d velocities
!     TRACER         : 3d tracers
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      in_filename,              &! filename of restart file
      in_restart_fmt             ! format of restart file (bin,nc)

   logical (POP_logical), intent(in) :: &
      lccsm_branch      ,&! flag if ccsm branch initialization
      lccsm_hybrid        ! flag if ccsm hybrid initialization

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode           ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      n, k,              &! dummy counters
      nu,                &! i/o unit for pointer file reads
      iblock,            &! local block index
      cindx,cindx2        ! indices into character strings

   real (POP_r8), dimension(nx_block,ny_block) :: &
      WORK1,WORK2        ! work space for pressure correction

   character (POP_charLength) ::  &
      restart_pointer_file, &! file name for restart pointer file
      short_name, long_name  ! tracer name temporaries

   logical (POP_logical) ::   &
      lcoupled_ts           ! flag to check whether coupled time step

   type (block) ::         &
      this_block            ! block information for current block

   type (datafile) :: &
      restart_file    ! io file descriptor

   type (io_dim) :: &
      i_dim, j_dim, &! dimension descriptors for horiz dims
      k_dim          ! dimension descriptor  for vertical levels

!-----------------------------------------------------------------------
!
!  if pointer files are used, pointer file must be read to get
!  actual filenames - skip this for branch initialization
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   read_restart_filename = char_blank
   restart_pointer_file = char_blank

   if (luse_pointer_files) then
      call get_unit(nu)
      if (my_task == master_task) then
         restart_pointer_file = pointer_filename
         cindx = len_trim(pointer_filename) + 1
         cindx2= cindx + 7
         restart_pointer_file(cindx:cindx2) = '.restart'
         write(stdout,*) 'Reading pointer file: ', &
                          trim(restart_pointer_file)
         call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
         open(nu, file=trim(restart_pointer_file), form='formatted', &
                  status='old')
         read(nu,'(a)') read_restart_filename
         close(nu)
      endif
      call release_unit(nu)

      call broadcast_scalar(read_restart_filename, master_task)

!-----------------------------------------------------------------------
!
!  otherwise use input filename
!
!-----------------------------------------------------------------------

   else
      cindx2 = len_trim(in_filename)
      read_restart_filename(1:cindx2) = trim(in_filename)
   endif

!-----------------------------------------------------------------------
!
!  create input file and define scalars with default values
!
!-----------------------------------------------------------------------

   restart_file =  construct_file(in_restart_fmt,                   &
                                  full_name=trim(read_restart_filename), &
                                  record_length=rec_type_dbl,       &
                                  recl_words=nx_global*ny_global)

   !*** set some defaults for namelist variables not initialized
   !*** under some options

   tlast_ice = c0
   lcoupled_ts = .false.

   !*** add defaults as file attributes

   call add_attrib_file(restart_file, 'runid',    runid   )
   call add_attrib_file(restart_file, 'iyear',    iyear   )
   call add_attrib_file(restart_file, 'imonth',   imonth  )
   call add_attrib_file(restart_file, 'iday',     iday    )
   call add_attrib_file(restart_file, 'ihour',    ihour   )
   call add_attrib_file(restart_file, 'iminute',  iminute )
   call add_attrib_file(restart_file, 'isecond',  isecond )
   call add_attrib_file(restart_file, 'iyear0',   iyear0  )
   call add_attrib_file(restart_file, 'imonth0',  imonth0 )
   call add_attrib_file(restart_file, 'iday0',    iday0   )
   call add_attrib_file(restart_file, 'ihour0',   ihour0  )
   call add_attrib_file(restart_file, 'iminute0', iminute0)
   call add_attrib_file(restart_file, 'isecond0', isecond0)
   call add_attrib_file(restart_file, 'dtt',      dtt     )
   call add_attrib_file(restart_file, 'iday_of_year', iday_of_year)
   call add_attrib_file(restart_file, 'iday_of_year_last',       &
                                       iday_of_year_last)
   call add_attrib_file(restart_file, 'elapsed_days',    elapsed_days)
   call add_attrib_file(restart_file, 'elapsed_months',  elapsed_months)
   call add_attrib_file(restart_file, 'elapsed_years',   elapsed_years)
   call add_attrib_file(restart_file, 'elapsed_days_this_year',  &
                                       elapsed_days_this_year)
   call add_attrib_file(restart_file, 'seconds_this_day',        &
                                       seconds_this_day)
   call add_attrib_file(restart_file, 'seconds_this_day_next',   &
                                       seconds_this_day_next)
   call add_attrib_file(restart_file, 'seconds_this_year',       &
                                       seconds_this_year)
   call add_attrib_file(restart_file, 'seconds_this_year_next',  &
                                       seconds_this_year_next)
   call add_attrib_file(restart_file, 'nsteps_total'  , nsteps_total)
   call add_attrib_file(restart_file, 'eod'     , eod     )
   call add_attrib_file(restart_file, 'eod_last', eod_last)
   call add_attrib_file(restart_file, 'eom'     , eom     )
   call add_attrib_file(restart_file, 'eom_last', eom_last)
   call add_attrib_file(restart_file, 'eom_next', eom_next)
   call add_attrib_file(restart_file, 'eoy'     , eoy     )
   call add_attrib_file(restart_file, 'eoy_last', eoy_last)
   call add_attrib_file(restart_file, 'midnight_last', midnight_last)
   call add_attrib_file(restart_file, 'adjust_year_next',        &
                                       adjust_year_next)
   call add_attrib_file(restart_file, 'newday' , newday  )
   call add_attrib_file(restart_file, 'newhour', newhour )
   call add_attrib_file(restart_file, 'leapyear',leapyear)
   call add_attrib_file(restart_file, 'days_in_year', days_in_year)
   call add_attrib_file(restart_file, 'days_in_prior_year',      &
                                       days_in_prior_year)
   call add_attrib_file(restart_file, 'seconds_in_year',         &
                                       seconds_in_year)
   call add_attrib_file(restart_file, 'hours_in_year', hours_in_year)
   call add_attrib_file(restart_file, 'tlast_ice',     tlast_ice  )
   call add_attrib_file(restart_file, 'lcoupled_ts',   lcoupled_ts)
   call add_attrib_file(restart_file, 'shf_interp_last',         &
                                       shf_interp_last)
   call add_attrib_file(restart_file, 'sfwf_interp_last',        &
                                       sfwf_interp_last)
   call add_attrib_file(restart_file, 'ws_interp_last',          &
                                       ws_interp_last)
   call add_attrib_file(restart_file, 'ap_interp_last',          &
                                       ap_interp_last)
   call add_attrib_file(restart_file, 'pt_interior_interp_last', &
                                       pt_interior_interp_last)
   call add_attrib_file(restart_file, 's_interior_interp_last',  &
                                       s_interior_interp_last)
   call add_attrib_file(restart_file, 'sum_precip' , sum_precip )
   call add_attrib_file(restart_file, 'precip_fact', precip_fact)
   call add_attrib_file(restart_file, 'ssh_initial', ssh_initial)

   short_name = char_blank
   do k=1,km
      write(short_name,'(a11,i3.3)') 'sal_initial',k
      call add_attrib_file(restart_file,trim(short_name),sal_initial(k))
   end do

!-----------------------------------------------------------------------
!
!  open a file and extract scalars as file attributes
!  do not extract if this is a ccsm branch initialization - values are set elsewhere
!
!-----------------------------------------------------------------------

   !*** open file and read attributes

   call data_set(restart_file, 'open_read')

   !*** extract scalars if not a ccsm branch initialization

   if (.not. lccsm_branch .and. .not. lccsm_hybrid) then
      call extract_attrib_file(restart_file, 'runid',    runid   )
   endif


   if (.not. lccsm_hybrid) then
      call extract_attrib_file(restart_file, 'iyear',    iyear   )
      call extract_attrib_file(restart_file, 'imonth',   imonth  )
      call extract_attrib_file(restart_file, 'iday',     iday    )
      call extract_attrib_file(restart_file, 'ihour',    ihour   )
      call extract_attrib_file(restart_file, 'iminute',  iminute )
      call extract_attrib_file(restart_file, 'isecond',  isecond )
      call extract_attrib_file(restart_file, 'iyear0',   iyear0  )
      call extract_attrib_file(restart_file, 'imonth0',  imonth0 )
      call extract_attrib_file(restart_file, 'iday0',    iday0   )
      call extract_attrib_file(restart_file, 'ihour0',   ihour0  )
      call extract_attrib_file(restart_file, 'iminute0', iminute0)
      call extract_attrib_file(restart_file, 'isecond0', isecond0)
      call extract_attrib_file(restart_file, 'dtt',      dtt     )
      call extract_attrib_file(restart_file, 'iday_of_year',     &
                                              iday_of_year)
      call extract_attrib_file(restart_file, 'iday_of_year_last',&
                                              iday_of_year_last)
      call extract_attrib_file(restart_file, 'elapsed_days',     &
                                              elapsed_days)
      call extract_attrib_file(restart_file, 'elapsed_months',   &
                                              elapsed_months)
      call extract_attrib_file(restart_file, 'elapsed_years',    &
                                              elapsed_years)
      call extract_attrib_file(restart_file, 'elapsed_days_this_year', &
                                              elapsed_days_this_year)
      call extract_attrib_file(restart_file, 'seconds_this_day',       &
                                              seconds_this_day)
      call extract_attrib_file(restart_file, 'seconds_this_day_next',  &
                                              seconds_this_day_next)
      call extract_attrib_file(restart_file, 'seconds_this_year',      &
                                              seconds_this_year)
      call extract_attrib_file(restart_file, 'seconds_this_year_next',  &
                                              seconds_this_year_next)
      call extract_attrib_file(restart_file, 'nsteps_total',     &
                                                    nsteps_total)
      call extract_attrib_file(restart_file, 'eod'     , eod     )
      call extract_attrib_file(restart_file, 'eod_last', eod_last)
      call extract_attrib_file(restart_file, 'eom'     , eom     )
      call extract_attrib_file(restart_file, 'eom_last', eom_last)
      call extract_attrib_file(restart_file, 'eom_next', eom_next)
      call extract_attrib_file(restart_file, 'eoy'     , eoy     )
      call extract_attrib_file(restart_file, 'eoy_last', eoy_last)
      call extract_attrib_file(restart_file, 'midnight_last',     &
                                              midnight_last)
      call extract_attrib_file(restart_file, 'adjust_year_next',  &
                                              adjust_year_next)
      call extract_attrib_file(restart_file, 'newday' , newday  )
      call extract_attrib_file(restart_file, 'newhour', newhour )
      call extract_attrib_file(restart_file, 'leapyear',leapyear)
      call extract_attrib_file(restart_file, 'days_in_year',      &
                                              days_in_year)
      call extract_attrib_file(restart_file, 'days_in_prior_year',&
                                              days_in_prior_year)
      call extract_attrib_file(restart_file, 'seconds_in_year',   &
                                              seconds_in_year)
      call extract_attrib_file(restart_file, 'hours_in_year',     &
                                              hours_in_year)
      call extract_attrib_file(restart_file, 'tlast_ice', tlast_ice  )
   endif ! .not. lccsm_hybrid


      call extract_attrib_file(restart_file, 'lcoupled_ts',       &
                                              lcoupled_ts)

   if (.not. lccsm_hybrid) then
      call extract_attrib_file(restart_file, 'shf_interp_last',   &
                                              shf_interp_last)
      call extract_attrib_file(restart_file, 'sfwf_interp_last',  &
                                              sfwf_interp_last)
      call extract_attrib_file(restart_file, 'ws_interp_last',    &
                                              ws_interp_last)
      call extract_attrib_file(restart_file, 'ap_interp_last',    &
                                              ap_interp_last)
      call extract_attrib_file(restart_file, 'pt_interior_interp_last',&
                                              pt_interior_interp_last)
      call extract_attrib_file(restart_file, 's_interior_interp_last', &
                                              s_interior_interp_last)
      call extract_attrib_file(restart_file, 'sum_precip' , sum_precip )
      call extract_attrib_file(restart_file, 'precip_fact', precip_fact)
      call extract_attrib_file(restart_file, 'ssh_initial', ssh_initial)


      short_name = char_blank
      do k=1,km
         write(short_name,'(a11,i3.3)') 'sal_initial',k
         call extract_attrib_file(restart_file, trim(short_name), &
                                                sal_initial(k))
      end do

      call int_to_char(4, iyear,cyear)
      call int_to_char(2, imonth,cmonth)
      call int_to_char(2, iday,cday)
      cmonth3 = month3_all(imonth)

      !*** set old value for the time flag 'coupled_ts'

      if (lcoupled_ts) then
         call register_string('coupled_ts_last_true')
         ! coupled_ts will be initialized accordingly in pop_init_coupled
      endif

   endif ! .not. lccsm_hybrid

!-----------------------------------------------------------------------
!
!  define all fields to be read
!
!-----------------------------------------------------------------------

   !*** define dimensions

   i_dim = construct_io_dim('i', nx_global)
   j_dim = construct_io_dim('j', ny_global)
   k_dim = construct_io_dim('k', km)

   UBTROP_CUR = construct_io_field('UBTROP_CUR', dim1=i_dim, dim2=j_dim,        &
                   long_name='U barotropic velocity at current time', &
                   units    ='cm/s', grid_loc ='2220',                &
                   field_loc = field_loc_NEcorner,                    &
                   field_type = field_type_vector,                    &
                   d2d_array =UBTROP(:,:,curtime,:))
   call data_set (restart_file, 'define', UBTROP_CUR)

   UBTROP_OLD = construct_io_field('UBTROP_OLD', dim1=i_dim, dim2=j_dim,        &
                   long_name='U barotropic velocity at old time',     &
                   units    ='cm/s', grid_loc ='2220',                &
                   field_loc = field_loc_NEcorner,                    &
                   field_type = field_type_vector,                    &
                   d2d_array =UBTROP(:,:,oldtime,:))
   call data_set (restart_file, 'define', UBTROP_OLD)

   VBTROP_CUR = construct_io_field('VBTROP_CUR', dim1=i_dim, dim2=j_dim,        &
                   long_name='V barotropic velocity at current time', &
                   units    ='cm/s', grid_loc ='2220',                &
                   field_loc = field_loc_NEcorner,                    &
                   field_type = field_type_vector,                    &
                   d2d_array =VBTROP(:,:,curtime,:))
   call data_set (restart_file, 'define', VBTROP_CUR)

   VBTROP_OLD = construct_io_field('VBTROP_OLD', dim1=i_dim, dim2=j_dim,        &
                   long_name='V barotropic velocity at old time',     &
                   units    ='cm/s', grid_loc ='2220',                &
                   field_loc = field_loc_NEcorner,                    &
                   field_type = field_type_vector,                    &
                   d2d_array =VBTROP(:,:,oldtime,:))
   call data_set (restart_file, 'define', VBTROP_OLD)

   PSURF_CUR = construct_io_field('PSURF_CUR', dim1=i_dim, dim2=j_dim,          &
                   long_name='surface pressure at current time',      &
                   units    ='dyne/cm2', grid_loc ='2110',            &
                   field_loc = field_loc_center,                      &
                   field_type = field_type_scalar,                    &
                   d2d_array =PSURF(:,:,curtime,:))
   call data_set (restart_file, 'define', PSURF_CUR)

   PSURF_OLD = construct_io_field('PSURF_OLD', dim1=i_dim, dim2=j_dim,          &
                   long_name='surface pressure at old time',          &
                   units    ='dyne/cm2', grid_loc ='2110',            &
                   field_loc = field_loc_center,                      &
                   field_type = field_type_scalar,                    &
                   d2d_array =PSURF(:,:,oldtime,:))
   call data_set (restart_file, 'define', PSURF_OLD)

   GRADPX_CUR = construct_io_field('GRADPX_CUR', dim1=i_dim, dim2=j_dim,        &
                   long_name='sfc press gradient in x at current time',&
                   units    ='dyne/cm3', grid_loc ='2220',            &
                   field_loc = field_loc_NEcorner,                    &
                   field_type = field_type_vector,                    &
                   d2d_array =GRADPX(:,:,curtime,:))
   call data_set (restart_file, 'define', GRADPX_CUR)

   GRADPX_OLD = construct_io_field('GRADPX_OLD', dim1=i_dim, dim2=j_dim,        &
                   long_name='sfc press gradient in x at old time',   &
                   units    ='dyne/cm3', grid_loc ='2220',            &
                   field_loc = field_loc_NEcorner,                    &
                   field_type = field_type_vector,                    &
                   d2d_array =GRADPX(:,:,oldtime,:))
   call data_set (restart_file, 'define', GRADPX_OLD)

   GRADPY_CUR = construct_io_field('GRADPY_CUR', dim1=i_dim, dim2=j_dim,        &
                   long_name='sfc press gradient in y at current time',&
                   units    ='dyne/cm3', grid_loc ='2220',            &
                   field_loc = field_loc_NEcorner,                    &
                   field_type = field_type_vector,                    &
                   d2d_array =GRADPY(:,:,curtime,:))
   call data_set (restart_file, 'define', GRADPY_CUR)

   GRADPY_OLD = construct_io_field('GRADPY_OLD', dim1=i_dim, dim2=j_dim,        &
                   long_name='sfc press gradient in y at old time',   &
                   units    ='dyne/cm3', grid_loc ='2220',            &
                   field_loc = field_loc_NEcorner,                    &
                   field_type = field_type_vector,                    &
                   d2d_array =GRADPY(:,:,oldtime,:))
   call data_set (restart_file, 'define', GRADPY_OLD)

   PGUESSd = construct_io_field('PGUESS', dim1=i_dim, dim2=j_dim,               &
                   long_name='guess for sfc pressure at new time',    &
                   units    ='dyne/cm2', grid_loc ='2110',            &
                   field_loc = field_loc_center,                      &
                   field_type = field_type_scalar,                    &
                   d2d_array =PGUESS)
   call data_set (restart_file, 'define', PGUESSd)

   if (sfc_layer_type == sfc_layer_varthick) then
      FW_OLDd = construct_io_field('FW_OLD', dim1=i_dim, dim2=j_dim,            &
                   long_name='fresh water input at old time',         &
                   grid_loc ='2110',                                  &
                   field_loc = field_loc_center,                      &
                   field_type = field_type_scalar,                    &
                   d2d_array = FW_OLD)
      call data_set (restart_file, 'define', FW_OLDd)
      FW_FREEZEd = construct_io_field('FW_FREEZE', dim1=i_dim, dim2=j_dim,      &
                   long_name='water flux due to frazil ice formation',&
                   grid_loc ='2110',                                  &
                   field_loc = field_loc_center,                      &
                   field_type = field_type_scalar,                    &
                   d2d_array = FW_FREEZE)
      call data_set (restart_file, 'define', FW_FREEZEd)
   endif

   if (liceform) then
      if (lcoupled_ts) then
        QFLUXd = construct_io_field('QFLUX', dim1=i_dim, dim2=j_dim,            &
                     long_name='Internal Ocean Heat Flux Due to Ice Formation',&
                     grid_loc ='2110',                                &
                     field_loc = field_loc_center,                    &
                     field_type = field_type_scalar,                  &
                     d2d_array = QFLUX)
        call data_set (restart_file, 'define', QFLUXd)
      else
        AQICEd = construct_io_field('AQICE', dim1=i_dim, dim2=j_dim,            &
                     long_name='accumulated ice melt/heat',           &
                     grid_loc ='2110',                                &
                     field_loc = field_loc_center,                    &
                     field_type = field_type_scalar,                  &
                     d2d_array = AQICE)
        call data_set (restart_file, 'define', AQICEd)
      endif
   endif

   UVEL_CUR = construct_io_field('UVEL_CUR', dim1=i_dim, dim2=j_dim, dim3=k_dim,&
                   long_name='U velocity at current time',            &
                   units    ='cm/s',                                  &
                   grid_loc ='3221',                                  &
                   field_loc = field_loc_NEcorner,                    &
                   field_type = field_type_vector,                    &
                   d3d_array = UVEL(:,:,:,curtime,:))
   call data_set (restart_file, 'define', UVEL_CUR)


   UVEL_OLD = construct_io_field('UVEL_OLD', dim1=i_dim, dim2=j_dim, dim3=k_dim,&
                   long_name='U velocity at old time',                &
                   units    ='cm/s',                                  &
                   grid_loc ='3221',                                  &
                   field_loc = field_loc_NEcorner,                    &
                   field_type = field_type_vector,                    &
                   d3d_array = UVEL(:,:,:,oldtime,:))
   call data_set (restart_file, 'define', UVEL_OLD)

   VVEL_CUR = construct_io_field('VVEL_CUR', dim1=i_dim, dim2=j_dim, dim3=k_dim,&
                   long_name='V velocity at current time',            &
                   units    ='cm/s',                                  &
                   grid_loc ='3221',                                  &
                   field_loc = field_loc_NEcorner,                    &
                   field_type = field_type_vector,                    &
                   d3d_array = VVEL(:,:,:,curtime,:))
   call data_set (restart_file, 'define', VVEL_CUR)

   VVEL_OLD = construct_io_field('VVEL_OLD', dim1=i_dim, dim2=j_dim, dim3=k_dim,&
                   long_name='V velocity at old time',                &
                   units    ='cm/s',                                  &
                   grid_loc ='3221',                                  &
                   field_loc = field_loc_NEcorner,                    &
                   field_type = field_type_vector,                    &
                   d3d_array = VVEL(:,:,:,oldtime,:))
   call data_set (restart_file, 'define', VVEL_OLD)

   do n=1,2
      short_name = char_blank
      short_name = trim(tracer_d(n)%short_name)/&
                                                &/'_CUR'
      long_name = char_blank
      long_name = trim(tracer_d(n)%long_name)/&
                                              &/' at current time'
      TRACER_CUR(n) = construct_io_field(trim(short_name),            &
                   dim1=i_dim, dim2=j_dim, dim3=k_dim,                          &
                   long_name=trim(long_name),                         &
                   units    =trim(tracer_d(n)%units),                 &
                   grid_loc ='3111',                                  &
                   field_loc = field_loc_center,                      &
                   field_type = field_type_scalar,                    &
                   d3d_array = TRACER(:,:,:,n,curtime,:))
      call data_set (restart_file, 'define', TRACER_CUR(n))
   end do


     do n=1,2

      short_name = char_blank
      short_name = trim(tracer_d(n)%short_name)/&
                                                &/'_OLD'
      long_name = char_blank
      long_name = trim(tracer_d(n)%long_name)/&
                                              &/' at old time'

      TRACER_OLD(n) = construct_io_field(trim(short_name),            &
                      dim1=i_dim, dim2=j_dim, dim3=k_dim,                       &
                      long_name=trim(long_name),                      &
                      units    =trim(tracer_d(n)%units),              &
                      grid_loc ='3111',                               &
                      field_loc = field_loc_center,                   &
                      field_type = field_type_scalar,                 &
                      d3d_array = TRACER(:,:,:,n,oldtime,:))

      call data_set (restart_file, 'define', TRACER_OLD(n))
   end do

!-----------------------------------------------------------------------
!
!  now we actually read each field
!  after reading, get rid of io field descriptors and close file
!
!-----------------------------------------------------------------------

   call data_set (restart_file, 'read', UBTROP_CUR)
   call data_set (restart_file, 'read', UBTROP_OLD)
   call data_set (restart_file, 'read', VBTROP_CUR)
   call data_set (restart_file, 'read', VBTROP_OLD)
   call data_set (restart_file, 'read', PSURF_CUR)
   call data_set (restart_file, 'read', PSURF_OLD)
   call data_set (restart_file, 'read', GRADPX_CUR)
   call data_set (restart_file, 'read', GRADPX_OLD)
   call data_set (restart_file, 'read', GRADPY_CUR)
   call data_set (restart_file, 'read', GRADPY_OLD)
   call data_set (restart_file, 'read', PGUESSd)

   if (sfc_layer_type == sfc_layer_varthick) then
      call data_set (restart_file, 'read', FW_OLDd)
      call data_set (restart_file, 'read', FW_FREEZEd)
   endif
   if (liceform) then
      if (lcoupled_ts) then
        call data_set (restart_file, 'read', QFLUXd)
      else
        call data_set (restart_file, 'read', AQICEd)
      endif
   endif

   call data_set (restart_file, 'read', UVEL_CUR)
   call data_set (restart_file, 'read', UVEL_OLD)
   call data_set (restart_file, 'read', VVEL_CUR)
   call data_set (restart_file, 'read', VVEL_OLD)

   do n=1,2
      call data_set (restart_file, 'read', TRACER_CUR(n))
      call data_set (restart_file, 'read', TRACER_OLD(n))
   end do

   call destroy_io_field (UBTROP_CUR)
   call destroy_io_field (UBTROP_OLD)
   call destroy_io_field (VBTROP_CUR)
   call destroy_io_field (VBTROP_OLD)
   call destroy_io_field (PSURF_CUR)
   call destroy_io_field (PSURF_OLD)
   call destroy_io_field (GRADPX_CUR)
   call destroy_io_field (GRADPX_OLD)
   call destroy_io_field (GRADPY_CUR)
   call destroy_io_field (GRADPY_OLD)
   call destroy_io_field (PGUESSd)

   if (sfc_layer_type == sfc_layer_varthick) then
      call destroy_io_field (FW_OLDd)
      call destroy_io_field (FW_FREEZEd)
   endif
   if (liceform) then
      if (lcoupled_ts) then
        call destroy_io_field (QFLUXd)
      else
        call destroy_io_field (AQICEd)
      endif
   endif

   call destroy_io_field (UVEL_CUR)
   call destroy_io_field (UVEL_OLD)
   call destroy_io_field (VVEL_CUR)
   call destroy_io_field (VVEL_OLD)
   do n=1,2
      call destroy_io_field (TRACER_CUR(n))
      call destroy_io_field (TRACER_OLD(n))
   end do

   call data_set (restart_file, 'close')

   if (my_task == master_task) then
     write(stdout,blank_fmt)
     write(stdout,*) ' file read: ', trim(read_restart_filename)
   endif

   call destroy_file(restart_file)

!-----------------------------------------------------------------------
!
!  zero prognostic variables at land points
!
!-----------------------------------------------------------------------

   do iblock = 1,nblocks_clinic

      this_block = get_block(blocks_clinic(iblock),iblock)

      where (.not. CALCU(:,:,iblock))
         UBTROP(:,:,curtime,iblock) = c0
         VBTROP(:,:,curtime,iblock) = c0
         GRADPX(:,:,curtime,iblock) = c0
         GRADPY(:,:,curtime,iblock) = c0
         UBTROP(:,:,oldtime,iblock) = c0
         VBTROP(:,:,oldtime,iblock) = c0
         GRADPX(:,:,oldtime,iblock) = c0
         GRADPY(:,:,oldtime,iblock) = c0
      endwhere

      where (.not. CALCT(:,:,iblock))
         PSURF(:,:,curtime,iblock) = c0
         PSURF(:,:,oldtime,iblock) = c0
         PGUESS(:,:,iblock) = c0
      endwhere

      if (liceform) then
         if (lcoupled_ts) then
         where (.not. CALCT(:,:,iblock))
            QFLUX(:,:,iblock) = c0
         endwhere
         else
         where (.not. CALCT(:,:,iblock))
            AQICE(:,:,iblock) = c0
         endwhere
         endif
      endif

      if (sfc_layer_type == sfc_layer_varthick) then
         where (.not. CALCT(:,:,iblock)) 
            FW_OLD   (:,:,iblock) = c0
            FW_FREEZE(:,:,iblock) = c0
         endwhere
      endif

      if( overflows_on .and. overflows_interactive &
          .and. overflows_restart_type /= 'ccsm_startup' ) then
              ! Do not set sidewall velocities to zero when overflows
              ! on and interactive; otherwise, valid overflow velocities
              ! will be lost
      else
      do k = 1,km
         where (k > KMU(:,:,iblock))
            UVEL(:,:,k,curtime,iblock) = c0
            VVEL(:,:,k,curtime,iblock) = c0
         endwhere
      enddo
      endif

      do n = 1,2
         do k = 1,km
            where (k > KMT(:,:,iblock))
               TRACER(:,:,k,n,curtime,iblock) = c0
               TRACER(:,:,k,n,oldtime,iblock) = c0
            endwhere
         enddo
      enddo

!-----------------------------------------------------------------------
!
!     reset PSURF(oldtime) to eliminate error in barotropic continuity
!     eqn due to (possible) use of different timestep after restart
!
!     NOTE: use pressure_correction = .false. for exact restart
!
!-----------------------------------------------------------------------

      if (pressure_correction) then

         WORK1 = HU(:,:,iblock)*UBTROP(:,:,curtime,iblock)
         WORK2 = HU(:,:,iblock)*VBTROP(:,:,curtime,iblock)

         !*** use PSURF(oldtime) as temp
         call div(1, PSURF(:,:,oldtime,iblock),WORK1,WORK2,this_block)

         PSURF(:,:,oldtime,iblock) = PSURF(:,:,curtime,iblock) +  &
                            grav*dtp*PSURF(:,:,oldtime,iblock)*   &
                            TAREA_R(:,:,iblock)

      endif
   end do !block loop

   if (pressure_correction) then
      call POP_HaloUpdate(PSURF(:,:,oldtime,:), POP_haloClinic,       &
                          POP_gridHorzLocCenter, POP_fieldKindScalar, &
                          errorCode, fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'read_restart: error updating sfc pressure halo')
         return
      endif
   endif

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

   call POP_HaloUpdate(UBTROP,                         &
                       POP_haloClinic,                 &
                       POP_gridHorzLocNECorner,        &
                       POP_fieldKindVector, errorCode, &
                       fillValue = 0.0_POP_r8)

   call POP_HaloUpdate(VBTROP,                         &
                       POP_haloClinic,                 &
                       POP_gridHorzLocNECorner,        &
                       POP_fieldKindVector, errorCode, &
                       fillValue = 0.0_POP_r8)

   call POP_HaloUpdate(UVEL,                           &
                       POP_haloClinic,                 &
                       POP_gridHorzLocNECorner,        &
                       POP_fieldKindVector, errorCode, &
                       fillValue = 0.0_POP_r8)

   call POP_HaloUpdate(VVEL,                           &
                       POP_haloClinic,                 &
                       POP_gridHorzLocNECorner,        &
                       POP_fieldKindVector, errorCode, &
                       fillValue = 0.0_POP_r8)

   call POP_HaloUpdate(TRACER(:,:,:,:,curtime,:),      &
                       POP_haloClinic,                 &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

   call POP_HaloUpdate(TRACER(:,:,:,:,newtime,:),      &
                       POP_haloClinic,                 &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

   call POP_HaloUpdate(GRADPX,                         &
                       POP_haloClinic,                 &
                       POP_gridHorzLocNECorner,        &
                       POP_fieldKindVector, errorCode, &
                       fillValue = 0.0_POP_r8)

   call POP_HaloUpdate(GRADPY,                         &
                       POP_haloClinic,                 &
                       POP_gridHorzLocNECorner,        &
                       POP_fieldKindVector, errorCode, &
                       fillValue = 0.0_POP_r8)

   call POP_HaloUpdate(PSURF,                          &
                       POP_haloClinic,                 &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

   call POP_HaloUpdate(PGUESS,                         &
                       POP_haloClinic,                 &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

   if (sfc_layer_type == sfc_layer_varthick) then
   call POP_HaloUpdate(FW_OLD,                         &
                       POP_haloClinic,                 &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

   call POP_HaloUpdate(FW_FREEZE,                      &
                       POP_haloClinic,                 &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)
   if (liceform) then
      if (lcoupled_ts) then
   call POP_HaloUpdate(QFLUX,                      &
                       POP_haloClinic,                 &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)
      else
   call POP_HaloUpdate(AQICE,                      &
                       POP_haloClinic,                 &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)
      endif
   endif
   endif


!-----------------------------------------------------------------------
!EOC

 end subroutine read_restart

!***********************************************************************
!BOP
! !IROUTINE: write_restart
! !INTERFACE:

 subroutine write_restart(restart_type)

! !DESCRIPTION:
!  This routine writes all the data necessary for restarting a POP
!  simulation if it is determined that the time has come to write
!  the data.  It also returns the type of restart that was written
!  so that the tavg module can determine whether it need to write
!  a restart file for tavg fields.
!
! !REVISION HISTORY:
!  same as module
!
! !OUTPUT PARAMETERS:

   character(POP_charLength), intent(out) :: &
      restart_type  ! type of restart file written if any
                    ! possible values are: none,restart,even,odd,end

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   character (POP_charLength) :: &
      file_suffix          ! suffix to append to root filename

   integer (POP_i4) :: &
      k, n,              &! dummy counters
      nu                  ! i/o unit for pointer file writes

   character (POP_charLength) ::  &
      write_restart_filename, &! modified file name for restart file
      restart_pointer_file, &! file name for restart pointer file
      short_name,           &! temporary for short name for io fields
      long_name              ! temporary for long  name for io fields

   logical (POP_logical) ::   &
      lcoupled_ts           ! flag to check whether coupled time step

   type (datafile) :: &
      restart_file    ! io file descriptor

   type (io_dim) :: &
      i_dim, j_dim, &! dimension descriptors for horiz dims
      k_dim          ! dimension descriptor  for vertical levels

!-----------------------------------------------------------------------
!
!  always set restart_type, because it is used in write_tavg
!
!-----------------------------------------------------------------------

   restart_type = char_blank
   restart_type = 'none'

!-----------------------------------------------------------------------
!
!  check to see if it is time to begin regularly writing restart files
!
!-----------------------------------------------------------------------


   if (check_time_flag(out_stop_now) ) then
    ! procede regardless of time_to_start option
   else
    if ( .not. time_to_start(restart_start_iopt,restart_start)) return
   endif

!-----------------------------------------------------------------------
!
!  check to see whether it is time to write files
!
!-----------------------------------------------------------------------

   lrestart_write = .false.

   !*** write restart files if code is stopping for any reason

   if (check_time_flag(out_stop_now) .and. &
       (lrestart_on .or. leven_odd_on)) then

      lrestart_write = .true.
      restart_type = char_blank
      restart_type = 'end'
   endif

   !*** check if it is time for even/odd output
   !*** (but not if an end file is written)

   if (.not. lrestart_write .and. check_time_flag(evenodd_flag)  &
                      .and. .not. check_time_flag(out_stop_now)) then

      lrestart_write = .true.
      restart_type = char_blank

      if (last_even_odd == even) then
         restart_type = 'odd'
         last_even_odd = odd
      else
         restart_type = 'even'
         last_even_odd = even
      endif

   endif

   !*** check if it is time for regularly-scheduled restart output
   !*** note that this option over-rides others

   if (check_time_flag(restart_flag) .or.   &
       (check_time_flag(cpl_write_restart) .and. &
        (nsteps_this_interval.eq.nsteps_per_interval)) ) then

      lrestart_write = .true.
      restart_type = char_blank
      restart_type = 'restart'

   endif

   !*** turn off cpl_write_restart if necessary

   if (check_time_flag(cpl_write_restart) .and. & 
       (nsteps_this_interval.eq.nsteps_per_interval)) &
     call override_time_flag(cpl_write_restart, value=.false.)


!-----------------------------------------------------------------------
!
!  the rest of this routine is only executed if it is time to write a
!  restart file
!
!-----------------------------------------------------------------------

   if (lrestart_write) then


!-----------------------------------------------------------------------
!
!  create filename for user-supplied root and append date
!
!-----------------------------------------------------------------------

   write_restart_filename = char_blank
   file_suffix = char_blank

   if (registry_match('lccsm')) then
     call create_restart_suffix_ccsm(file_suffix, restart_type,freq_opt_nsecond)
   else
     call create_restart_suffix(file_suffix, restart_type)
   endif

   !*** must split concatenation operator to avoid preprocessor mangling

   write_restart_filename = trim(restart_outfile)/&
                                                  &/'.'/&
                                                  &/trim(file_suffix)

!-----------------------------------------------------------------------
!
!  create output file
!
!-----------------------------------------------------------------------

   restart_file =  construct_file(restart_fmt,                      &
                                  full_name=trim(write_restart_filename), &
                                  record_length=rec_type_dbl,       &
                                  recl_words=nx_global*ny_global)

!-----------------------------------------------------------------------
!
!  scalar variables in restart file are output as file attributes
!  so define them here
!
!-----------------------------------------------------------------------

   !*** query time_flag routine for present value of lcoupled_ts

   lcoupled_ts = check_time_flag (restart_cpl_ts)

   !*** add defaults as file attributes

   call add_attrib_file(restart_file, 'runid',    runid   )
   call add_attrib_file(restart_file, 'iyear',    iyear   )
   call add_attrib_file(restart_file, 'imonth',   imonth  )
   call add_attrib_file(restart_file, 'iday',     iday    )
   call add_attrib_file(restart_file, 'ihour',    ihour   )
   call add_attrib_file(restart_file, 'iminute',  iminute )
   call add_attrib_file(restart_file, 'isecond',  isecond )
   call add_attrib_file(restart_file, 'iyear0',   iyear0  )
   call add_attrib_file(restart_file, 'imonth0',  imonth0 )
   call add_attrib_file(restart_file, 'iday0',    iday0   )
   call add_attrib_file(restart_file, 'ihour0',   ihour0  )
   call add_attrib_file(restart_file, 'iminute0', iminute0)
   call add_attrib_file(restart_file, 'isecond0', isecond0)
   call add_attrib_file(restart_file, 'dtt',      dtt     )
   call add_attrib_file(restart_file, 'iday_of_year', iday_of_year)
   call add_attrib_file(restart_file, 'iday_of_year_last',       &
                                       iday_of_year_last)
   call add_attrib_file(restart_file, 'elapsed_days',   elapsed_days)
   call add_attrib_file(restart_file, 'elapsed_months', elapsed_months)
   call add_attrib_file(restart_file, 'elapsed_years',  elapsed_years)
   call add_attrib_file(restart_file, 'elapsed_days_this_year',  &
                                       elapsed_days_this_year)
   call add_attrib_file(restart_file, 'seconds_this_day',        &
                                       seconds_this_day)
   call add_attrib_file(restart_file, 'seconds_this_day_next',   &
                                       seconds_this_day_next)
   call add_attrib_file(restart_file, 'seconds_this_year',       &
                                       seconds_this_year)
   call add_attrib_file(restart_file, 'seconds_this_year_next',  &
                                       seconds_this_year_next)
   call add_attrib_file(restart_file, 'nsteps_total',            &
                                             nsteps_total)
   call add_attrib_file(restart_file, 'eod'     , eod     )
   call add_attrib_file(restart_file, 'eod_last', eod_last)
   call add_attrib_file(restart_file, 'eom'     , eom     )
   call add_attrib_file(restart_file, 'eom_last', eom_last)
   call add_attrib_file(restart_file, 'eom_next', eom_next)
   call add_attrib_file(restart_file, 'eoy'     , eoy     )
   call add_attrib_file(restart_file, 'eoy_last', eoy_last)
   call add_attrib_file(restart_file, 'midnight_last', midnight_last)
   call add_attrib_file(restart_file, 'adjust_year_next',        &
                                       adjust_year_next)
   call add_attrib_file(restart_file, 'newday' , newday  )
   call add_attrib_file(restart_file, 'newhour', newhour )
   call add_attrib_file(restart_file, 'leapyear',leapyear)
   call add_attrib_file(restart_file, 'days_in_year', days_in_year)
   call add_attrib_file(restart_file, 'days_in_prior_year',      &
                                       days_in_prior_year)
   call add_attrib_file(restart_file, 'seconds_in_year', &
                                       seconds_in_year)
   call add_attrib_file(restart_file, 'hours_in_year', hours_in_year)
   call add_attrib_file(restart_file, 'tlast_ice',     tlast_ice  )
   call add_attrib_file(restart_file, 'lcoupled_ts',   lcoupled_ts)
   call add_attrib_file(restart_file, 'shf_interp_last',         &
                                       shf_interp_last)
   call add_attrib_file(restart_file, 'sfwf_interp_last',        &
                                       sfwf_interp_last)
   call add_attrib_file(restart_file, 'ws_interp_last',          &
                                       ws_interp_last)
   call add_attrib_file(restart_file, 'ap_interp_last',          &
                                       ap_interp_last)
   call add_attrib_file(restart_file, 'pt_interior_interp_last', &
                                       pt_interior_interp_last)
   call add_attrib_file(restart_file, 's_interior_interp_last',  &
                                       s_interior_interp_last)
   call add_attrib_file(restart_file, 'sum_precip' , sum_precip )
   call add_attrib_file(restart_file, 'precip_fact', precip_fact)
   call add_attrib_file(restart_file, 'ssh_initial', ssh_initial)

   short_name = char_blank
   do k=1,km
      write(short_name,'(a11,i3.3)') 'sal_initial',k
      call add_attrib_file(restart_file,trim(short_name),sal_initial(k))
   end do

   if (nt > 2) call write_restart_passive_tracers(restart_file,'add_attrib_file')

!-----------------------------------------------------------------------
!
!  open a file (also writes scalars as attributes to file)
!
!-----------------------------------------------------------------------

   call data_set(restart_file, 'open')

!-----------------------------------------------------------------------
!
!  construct all fields to be written
!
!-----------------------------------------------------------------------

   !*** define dimensions

   i_dim = construct_io_dim('i', nx_global)
   j_dim = construct_io_dim('j', ny_global)
   k_dim = construct_io_dim('k', km)

   UBTROP_CUR = construct_io_field('UBTROP_CUR', dim1=i_dim, dim2=j_dim,        &
                   long_name='U barotropic velocity at current time', &
                   units    ='cm/s',                                  &
                   grid_loc ='2220',                                  &
                   d2d_array =UBTROP(:,:,curtime,:))

   UBTROP_OLD = construct_io_field('UBTROP_OLD', dim1=i_dim, dim2=j_dim,        &
                   long_name='U barotropic velocity at old time',     &
                   units    ='cm/s',                                  &
                   grid_loc ='2220',                                  &
                   d2d_array =UBTROP(:,:,oldtime,:))

   VBTROP_CUR = construct_io_field('VBTROP_CUR', dim1=i_dim, dim2=j_dim,        &
                   long_name='V barotropic velocity at current time', &
                   units    ='cm/s',                                  &
                   grid_loc ='2220',                                  &
                   d2d_array =VBTROP(:,:,curtime,:))

   VBTROP_OLD = construct_io_field('VBTROP_OLD', dim1=i_dim, dim2=j_dim,        &
                   long_name='V barotropic velocity at old time',     &
                   units    ='cm/s',                                  &
                   grid_loc ='2220',                                  &
                   d2d_array =VBTROP(:,:,oldtime,:))

   PSURF_CUR = construct_io_field('PSURF_CUR', dim1=i_dim, dim2=j_dim,          &
                   long_name='surface pressure at current time',      &
                   units    ='dyne/cm2',                              &
                   grid_loc ='2110',                                  &
                   d2d_array =PSURF(:,:,curtime,:))

   PSURF_OLD = construct_io_field('PSURF_OLD', dim1=i_dim, dim2=j_dim,          &
                   long_name='surface pressure at old time',          &
                   units    ='dyne/cm2',                              &
                   grid_loc ='2110',                                  &
                   d2d_array =PSURF(:,:,oldtime,:))

   GRADPX_CUR = construct_io_field('GRADPX_CUR', dim1=i_dim, dim2=j_dim,        &
                   long_name='sfc press gradient in x at current time',&
                   units    ='dyne/cm3',                              &
                   grid_loc ='2220',                                  &
                   d2d_array =GRADPX(:,:,curtime,:))

   GRADPX_OLD = construct_io_field('GRADPX_OLD', dim1=i_dim, dim2=j_dim,        &
                   long_name='sfc press gradient in x at old time',   &
                   units    ='dyne/cm3',                              &
                   grid_loc ='2220',                                  &
                   d2d_array =GRADPX(:,:,oldtime,:))

   GRADPY_CUR = construct_io_field('GRADPY_CUR', dim1=i_dim, dim2=j_dim,        &
                   long_name='sfc press gradient in y at current time',&
                   units    ='dyne/cm3',                              &
                   grid_loc ='2220',                                  &
                   d2d_array =GRADPY(:,:,curtime,:))

   GRADPY_OLD = construct_io_field('GRADPY_OLD', dim1=i_dim, dim2=j_dim,        &
                   long_name='sfc press gradient in y at old time',   &
                   units    ='dyne/cm3',                              &
                   grid_loc ='2220',                                  &
                   d2d_array =GRADPY(:,:,oldtime,:))

   PGUESSd = construct_io_field('PGUESS', dim1=i_dim, dim2=j_dim,               &
                   long_name='guess for sfc pressure at new time',    &
                   units    ='dyne/cm2',                              &
                   grid_loc ='2110',                                  &
                   d2d_array =PGUESS)

   if (sfc_layer_type == sfc_layer_varthick) then
      FW_OLDd = construct_io_field('FW_OLD', dim1=i_dim, dim2=j_dim,            &
                   long_name='fresh water input at old time',         &
                   grid_loc ='2110',                                  &
                   d2d_array = FW_OLD)
      FW_FREEZEd = construct_io_field('FW_FREEZE', dim1=i_dim, dim2=j_dim,      &
                   long_name='water flux due to frazil ice formation',&
                   grid_loc ='2110',                                  &
                   d2d_array = FW_FREEZE)
   endif

   if (liceform) then
      if (lcoupled_ts) then
        QFLUXd = construct_io_field('QFLUX', dim1=i_dim, dim2=j_dim,            &
                     long_name='Internal Ocean Heat Flux Due to Ice Formation',&
                     grid_loc ='2110',                                &
                     d2d_array = QFLUX)
      else
        AQICEd = construct_io_field('AQICE', dim1=i_dim, dim2=j_dim,            &
                     long_name='accumulated ice melt/heat',           &
                     grid_loc ='2110',                                &
                     d2d_array = AQICE)
      endif
   endif

   UVEL_CUR = construct_io_field('UVEL_CUR', dim1=i_dim, dim2=j_dim, dim3=k_dim,&
                   long_name='U velocity at current time',            &
                   units    ='cm/s',                                  &
                   grid_loc ='3221',                                  &
                   d3d_array = UVEL(:,:,:,curtime,:))

   UVEL_OLD = construct_io_field('UVEL_OLD', dim1=i_dim, dim2=j_dim, dim3=k_dim,&
                   long_name='U velocity at old time',                &
                   units    ='cm/s',                                  &
                   grid_loc ='3221',                                  &
                   d3d_array = UVEL(:,:,:,oldtime,:))

   VVEL_CUR = construct_io_field('VVEL_CUR', dim1=i_dim, dim2=j_dim, dim3=k_dim,&
                   long_name='V velocity at current time',            &
                   units    ='cm/s',                                  &
                   grid_loc ='3221',                                  &
                   d3d_array = VVEL(:,:,:,curtime,:))

   VVEL_OLD = construct_io_field('VVEL_OLD', dim1=i_dim, dim2=j_dim, dim3=k_dim,&
                   long_name='V velocity at old time',                &
                   units    ='cm/s',                                  &
                   grid_loc ='3221',                                  &
                   d3d_array = VVEL(:,:,:,oldtime,:))

   do n=1,nt
      short_name = char_blank
      short_name = trim(tracer_d(n)%short_name)/&
                                                &/'_CUR'
      long_name = char_blank
      long_name = trim(tracer_d(n)%long_name)/&
                                              &/' at current time'

      TRACER_CUR(n) = construct_io_field(trim(short_name),            &
                   dim1=i_dim, dim2=j_dim, dim3=k_dim,                          &
                   long_name=trim(long_name),                         &
                   units    =trim(tracer_d(n)%units),                 &
                   grid_loc ='3111',                                  &
                   d3d_array = TRACER(:,:,:,n,curtime,:))
   end do

   do n=1,nt
      short_name = char_blank
      short_name = trim(tracer_d(n)%short_name)/&
                                                &/'_OLD'
      long_name = char_blank
      long_name = trim(tracer_d(n)%long_name)/&
                                              &/' at old time'

      TRACER_OLD(n) = construct_io_field(trim(short_name),            &
                   dim1=i_dim, dim2=j_dim, dim3=k_dim,                          &
                   long_name=trim(long_name),                         &
                   units    =trim(tracer_d(n)%units),                 &
                   grid_loc ='3111',                                  &
                   d3d_array = TRACER(:,:,:,n,oldtime,:))
   end do

!-----------------------------------------------------------------------
!
!  define all fields to be read
!
!-----------------------------------------------------------------------

   !*** must call in this order for back compatibility

   call data_set (restart_file, 'define', UBTROP_CUR)
   call data_set (restart_file, 'define', UBTROP_OLD)
   call data_set (restart_file, 'define', VBTROP_CUR)
   call data_set (restart_file, 'define', VBTROP_OLD)
   call data_set (restart_file, 'define', PSURF_CUR)
   call data_set (restart_file, 'define', PSURF_OLD)
   call data_set (restart_file, 'define', GRADPX_CUR)
   call data_set (restart_file, 'define', GRADPX_OLD)
   call data_set (restart_file, 'define', GRADPY_CUR)
   call data_set (restart_file, 'define', GRADPY_OLD)
   call data_set (restart_file, 'define', PGUESSd)

   if (sfc_layer_type == sfc_layer_varthick) then
      call data_set (restart_file, 'define', FW_OLDd)
      call data_set (restart_file, 'define', FW_FREEZEd)
   endif
   if (liceform) then
      if (lcoupled_ts) then
        call data_set (restart_file, 'define', QFLUXd)
      else
        call data_set (restart_file, 'define', AQICEd)
      endif
   endif

   call data_set (restart_file, 'define', UVEL_CUR)
   call data_set (restart_file, 'define', UVEL_OLD)
   call data_set (restart_file, 'define', VVEL_CUR)
   call data_set (restart_file, 'define', VVEL_OLD)

   do n=1,nt
      call data_set (restart_file, 'define', TRACER_CUR(n))
      call data_set (restart_file, 'define', TRACER_OLD(n))
   end do

   if (nt > 2) call write_restart_passive_tracers(restart_file,'define')

!-----------------------------------------------------------------------
!
!  now we actually write each field
!
!-----------------------------------------------------------------------

   call data_set (restart_file, 'write', UBTROP_CUR)
   call data_set (restart_file, 'write', UBTROP_OLD)
   call data_set (restart_file, 'write', VBTROP_CUR)
   call data_set (restart_file, 'write', VBTROP_OLD)
   call data_set (restart_file, 'write', PSURF_CUR)
   call data_set (restart_file, 'write', PSURF_OLD)
   call data_set (restart_file, 'write', GRADPX_CUR)
   call data_set (restart_file, 'write', GRADPX_OLD)
   call data_set (restart_file, 'write', GRADPY_CUR)
   call data_set (restart_file, 'write', GRADPY_OLD)
   call data_set (restart_file, 'write', PGUESSd)

   if (sfc_layer_type == sfc_layer_varthick) then
      call data_set (restart_file, 'write', FW_OLDd)
      call data_set (restart_file, 'write', FW_FREEZEd)
   endif
   if (liceform) then
      if (lcoupled_ts) then
        call data_set (restart_file, 'write', QFLUXd)
      else
        call data_set (restart_file, 'write', AQICEd)
      endif
   endif

   call data_set (restart_file, 'write', UVEL_CUR)
   call data_set (restart_file, 'write', UVEL_OLD)
   call data_set (restart_file, 'write', VVEL_CUR)
   call data_set (restart_file, 'write', VVEL_OLD)
   do n=1,nt
      call data_set (restart_file, 'write', TRACER_CUR(n))
      call data_set (restart_file, 'write', TRACER_OLD(n))
   end do

   if (nt > 2) call write_restart_passive_tracers(restart_file,'write')

!-----------------------------------------------------------------------
!
!  close and destroy file
!
!-----------------------------------------------------------------------

   call data_set (restart_file, 'close')

   if (my_task == master_task) then
     write(stdout,blank_fmt)
     write(stdout,*) ' restart file written: ', trim(write_restart_filename)
   endif

   call destroy_file(restart_file)

!-----------------------------------------------------------------------
!
!  if pointer files are used, write filename to pointer file
!
!-----------------------------------------------------------------------

   if (luse_pointer_files) then
      call get_unit(nu)
      if (my_task == master_task) then
       restart_pointer_file = trim(pointer_filename)/&
                                                       &/'.restart'

       open(nu, file=restart_pointer_file, form='formatted', &
                status='unknown')
       write(nu,'(a)') trim(write_restart_filename)
       write(nu,'(a,a)') 'RESTART_FMT=',trim(restart_fmt)
       close(nu)
       write(stdout,blank_fmt)
       write(stdout,*) ' restart pointer file written: ',trim(restart_pointer_file)
       call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
     endif
     call release_unit(nu)
   endif


!-----------------------------------------------------------------------
!
!  finished writing file
!
!-----------------------------------------------------------------------

   endif ! lrestart_write

!-----------------------------------------------------------------------
!EOC

 end subroutine write_restart


!***********************************************************************
!BOP
! !IROUTINE: init_restart
! !INTERFACE:

 subroutine init_restart

! !DESCRIPTION:
!  Initializes quantities associated with writing all the data
!  necessary for restarting a simulation.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables and input namelist
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      n,                  &! tracer loop index
      nml_error            ! namelist i/o error flag

   character (POP_charLength) :: &
      restart_freq_opt,    &! input option for freq of restart dumps
      restart_start_opt     ! choice for starting regular restart writes 

   character (POP_charLength), parameter :: &
      start_fmt = "('regular restart writes will start at ',a,i8)"

   namelist /restart_nml/ restart_freq_opt, restart_freq, &
                          restart_outfile, restart_fmt,   &
                          leven_odd_on, even_odd_freq,    &
                          pressure_correction,            &
                          restart_start_opt, restart_start

!-----------------------------------------------------------------------
!
!     register init_restart
!
!-----------------------------------------------------------------------
      call register_string('init_restart')

!-----------------------------------------------------------------------
!
!  read namelist input and broadcast variables
!
!-----------------------------------------------------------------------

   restart_outfile   = 'd'
   restart_fmt       = 'bin'
   restart_freq_iopt = freq_opt_never
   restart_start_iopt= start_opt_nstep 
   restart_start     = 0
   restart_freq      = 100000
   leven_odd_on      = .false.
   even_odd_freq     = 100000
   pressure_correction = .false.

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=restart_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      exit_string = 'FATAL ERROR: reading restart_nml'
      call document ('init_restart', exit_string)
      call exit_POP (sigAbort, exit_string, out_unit=stdout)
   endif

   if (my_task == master_task) then
      select case (restart_freq_opt)
      case ('never')
         restart_freq_iopt = freq_opt_never
      case ('nyear')
         restart_freq_iopt = freq_opt_nyear
      case ('nmonth')
         restart_freq_iopt = freq_opt_nmonth
      case ('nday')
         restart_freq_iopt = freq_opt_nday
      case ('nhour')
         restart_freq_iopt = freq_opt_nhour
      case ('nsecond')
         restart_freq_iopt = freq_opt_nsecond
      case ('nstep')
         restart_freq_iopt = freq_opt_nstep
      case default
         restart_freq_iopt = -1000
      end select

      if (restart_freq_iopt /= freq_opt_never) then
         select case (restart_start_opt)
         case ('nstep')
            restart_start_iopt = start_opt_nstep
            write(stdout,start_fmt) 'step ', restart_start
         case ('nday')
            restart_start_iopt = start_opt_nday
            write(stdout,start_fmt) 'day  ', restart_start
         case ('nyear')
            restart_start_iopt = start_opt_nyear
            write(stdout,start_fmt) 'year ', restart_start
         case ('date')
            restart_start_iopt = start_opt_date
            write(stdout,start_fmt) '     ', restart_start
         case default
            restart_start_iopt = -1000
         end select
      endif

   endif

   call broadcast_scalar (restart_outfile,      master_task)
   call broadcast_scalar (restart_freq_iopt,    master_task)
   call broadcast_scalar (restart_freq,         master_task)
   call broadcast_scalar (restart_start_iopt,   master_task)
   call broadcast_scalar (restart_start,        master_task)
   call broadcast_scalar (restart_fmt,          master_task)
   call broadcast_scalar (leven_odd_on,         master_task)
   call broadcast_scalar (even_odd_freq,        master_task)
   call broadcast_scalar (pressure_correction,  master_task)

   if (restart_freq_iopt == -1000) then
      exit_string = 'FATAL ERROR: unknown restart frequency option'
      call document ('init_restart', exit_string)
      call exit_POP (sigAbort, exit_string, out_unit=stdout)
   else if (restart_start_iopt == -1000) then
      exit_string = 'FATAL ERROR: unknown restart start option'
      call document ('init_restart', exit_string)
      call exit_POP (sigAbort, exit_string, out_unit=stdout)
   else if (restart_freq_iopt == freq_opt_never) then
      lrestart_on = .false.
   else
      lrestart_on = .true.
   endif

!-----------------------------------------------------------------------
!
!  create some time flags
!
!-----------------------------------------------------------------------

   call init_time_flag('restart',restart_flag,  default=.false.,    &
                        owner    = 'init_restart',                  &
                        freq_opt = restart_freq_iopt,               &
                        freq     = restart_freq)

   if (leven_odd_on) then
      last_even_odd = even
      call init_time_flag('evenodd', evenodd_flag, default=.false., &
                           owner    = 'init_restart',               &
                           freq_opt = freq_opt_nstep,               &
                           freq     = even_odd_freq)
   else
      call init_time_flag('evenodd',evenodd_flag, default=.false.,  &
                           owner    = 'init_restart',               &
                           freq_opt = freq_opt_never,               &
                           freq     = even_odd_freq)
   endif

!-----------------------------------------------------------------------
!
!  get handle for time flags defined in other modules
!
!-----------------------------------------------------------------------

   call access_time_flag('cpl_write_restart',cpl_write_restart)
   call access_time_flag('coupled_ts',restart_cpl_ts)
   call access_time_flag('stop_now',out_stop_now)

!-----------------------------------------------------------------------
!EOC

 end subroutine init_restart

!***********************************************************************
!BOP
! !IROUTINE: create_restart_suffix
! !INTERFACE:

 subroutine create_restart_suffix(file_suffix, restart_type)

! !DESCRIPTION:
!  Determines suffix to append to restart files based on restart type.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      restart_type           ! type of restart file to be written
                             ! (restart,even,odd,end)

! !OUTPUT PARAMETERS:

   character (POP_charLength), intent(out) :: &
      file_suffix            ! suffix to append to root filename

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     local variable
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      cindx, cindx2,     &! indices into character strings
      len_date            ! length of date string

   character (POP_charLength) :: &
      char_temp            ! temp character space

   character (10) :: &
      cdate          ! date string

!-----------------------------------------------------------------------
!
!  clear character strings
!
!-----------------------------------------------------------------------

   file_suffix = char_blank
   char_temp   = char_blank

!-----------------------------------------------------------------------
!
!  for even, odd, or end, simply add the appropriate string
!
!-----------------------------------------------------------------------

   select case (trim(restart_type))
   case('end')
      file_suffix = trim(runid)/&
                                &/'.end'
   case('even')
      file_suffix = trim(runid)/&
                                &/'.even'
   case('odd')
      file_suffix = trim(runid)/&
                                &/'.odd'

!-----------------------------------------------------------------------
!
!  for a regular restart file, append a date/time string
!
!-----------------------------------------------------------------------

   case('restart')

      if (date_separator == ' ') then
         len_date = 8
         cdate(1:4) = cyear
         cdate(5:6) = cmonth
         cdate(7:8) = cday
         cdate(9:10)= '  '
      else
         len_date = 10
         cdate(1:4) = cyear
         cdate(5:5) = date_separator
         cdate(6:7) = cmonth
         cdate(8:8) = date_separator
         cdate(9:10) = cday
      endif

      select case (restart_freq_iopt)
      case (freq_opt_nyear, freq_opt_nmonth, freq_opt_nday)

         !*** append the date after the runid

         file_suffix = trim(runid)/&
                                   &/'.'/&
                                         &/trim(cdate)

      case (freq_opt_nhour)

         !*** append the date to runid and add hour

         write(file_suffix,'(i2)') ihour
         char_temp = adjustl(file_suffix)

         file_suffix = trim(runid)/&
                                   &/'.'/&
                                   &/trim(cdate)/&
                                   &/'.h'/&
                                   &/trim(char_temp)

      case (freq_opt_nsecond)

         !*** append the date to runid and the elapsed seconds in day

         write (file_suffix,'(i6)') isecond
         char_temp = adjustl(file_suffix)

         file_suffix = trim(runid)/&
                                   &/'.'/&
                                   &/trim(cdate)/&
                                   &/'.s'/&
                                   &/trim(char_temp)

      case (freq_opt_nstep)

         !*** append the step number

         write (file_suffix,'(i10)') nsteps_total
         char_temp = adjustl(file_suffix)

         file_suffix = trim(runid)/&
                                   &/'.'/&
                                         &/trim(char_temp)

      case default
         file_suffix = trim(runid)
      end select

   end select

!-----------------------------------------------------------------------
!EOC

 end subroutine create_restart_suffix
!***********************************************************************

!BOP
! !IROUTINE: create_restart_suffix_ccsm
! !INTERFACE:

 subroutine create_restart_suffix_ccsm(file_suffix, restart_type,in_freq_opt)

! !DESCRIPTION:
!  Determines suffix to append to CCSM restart files based on restart type.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in)      :: &
      restart_type                  ! type of restart file to be written
                                    ! (restart,even,odd,end)
   integer (POP_i4), intent(in) :: &
      in_freq_opt                   ! type of ccsm date string
                                    ! (annual, monthly, daily, or instantaneous)

! !OUTPUT PARAMETERS:

   character (POP_charLength), intent(out) :: &
      file_suffix            ! suffix to append to root filename

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     local variable
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      cindx, cindx2,     &! indices into character strings
      len_date            ! length of date string

   character (POP_charLength) :: &
      char_temp,           &! temp character space
      ccsm_date_string

   character (10) :: &
      cdate          ! date string

!-----------------------------------------------------------------------
!
!  clear character strings
!
!-----------------------------------------------------------------------

   file_suffix = char_blank
   char_temp   = char_blank

!-----------------------------------------------------------------------
!
!  for even, odd, or end, simply add the appropriate string
!
!-----------------------------------------------------------------------

   select case (trim(restart_type))
   case('end')
      file_suffix = trim(runid)/&
                                &/'.end'
   case('even')
      file_suffix = trim(runid)/&
                                &/'.even'
   case('odd')
      file_suffix = trim(runid)/&
                                &/'.odd'

!-----------------------------------------------------------------------
!
!  for a regular restart file, append a date/time string
!
!-----------------------------------------------------------------------

   case('restart')

      char_temp   = char_blank
      file_suffix = char_blank

      select case (in_freq_opt)
      case (freq_opt_nyear)
        char_temp = 'y'

      case (freq_opt_nmonth)
        char_temp = 'ym'

      case (freq_opt_nday)
        char_temp = 'ymd'

      case (freq_opt_nhour)
        char_temp = 'ymds'

      case (freq_opt_nsecond)
        char_temp = 'ymds'

      case (freq_opt_nstep)
        char_temp = 'ymds'
 
      case default
        char_temp = 'ymds'
      end select


      call ccsm_date_stamp (ccsm_date_string, char_temp)
 
      file_suffix = trim(ccsm_date_string)

   end select


!-----------------------------------------------------------------------
!
!  for a restart file in netCDF format, append the suffix '.nc'
!
!-----------------------------------------------------------------------

   select case (trim(restart_fmt))
     case('nc')
       file_suffix = trim(file_suffix)/&
                                 &/'.'/&
                                 &/'nc'
   end select

!-----------------------------------------------------------------------
!EOC

 end subroutine create_restart_suffix_ccsm



!***********************************************************************

 end module restart

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
