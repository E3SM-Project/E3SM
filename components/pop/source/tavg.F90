!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module tavg

!BOP
! !MODULE: tavg
! !DESCRIPTION:
!  This module contains data types and routines for computing running 
!  time-averages of selected fields and writing this data to files.
!
! !REVISION HISTORY:
!  SVN:$Id: tavg.F90 56176 2013-12-20 18:35:46Z mlevy@ucar.edu $
!  

! !USES:

   use POP_KindsMod
   use POP_IOUnitsMod
   use POP_ErrorMod

   use kinds_mod
   use blocks
   use distribution
   use domain
   use constants
   use prognostic
   use grid
   use time_management
   use global_reductions
   use broadcast
   use io
   use io_types
   use exit_mod
 
   !*** ccsm
   use gather_scatter
   use operators, only: zcurl
   use io_ccsm
   use io_tools
   use diag_bsf, only: pcg_diag_bsf_solver, init_diag_bsf
   use diags_on_lat_aux_grid
   use registry
   use timers

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_tavg,              &
             define_tavg_field,      &
             tavg_increment_sum_qflux,&
             accumulate_tavg_field,  &
             accumulate_tavg_now,    &
             set_in_tavg_contents,   &
             tavg_in_this_stream,    &
             tavg_in_which_stream,   &
             write_tavg,             &
             read_tavg,              &
             tavg_set_flag,          &
             final_tavg
 
   !*** ccsm
   public :: tavg_id,                &
             tavg_global_sum_2D


! !PUBLIC DATA MEMBERS:

   logical (log_kind), public :: &
      ltavg_restart = .false.    ! run started from restart

   integer (i4), parameter, public ::  &
      tavg_method_unknown  = 0,         &
      tavg_method_avg      = 1,         &
      tavg_method_min      = 2,         &
      tavg_method_max      = 3,         &
      tavg_method_qflux    = 4,         &
      tavg_method_constant = 5

   integer (int_kind), parameter, public :: &
      max_avail_tavg_streams = 9     ! limit on number of streams; coding limitations restrict this to <= 9 

   real (r8),public,dimension(max_avail_tavg_streams) ::  &
      tavg_sum           ! accumulated time (in seconds)

   !*** ccsm
   real (r8),dimension(max_avail_tavg_streams), public ::  &
      tavg_sum_qflux


!-----------------------------------------------------------------------
!
!  tavg field descriptor data type and array of such types
!
!-----------------------------------------------------------------------

! !PUBLIC TYPES:
   !*** ccsm
   type,public :: tavg_field_desc_ccsm
      character(char_len)     :: short_name     ! short name for field
      character(char_len)     :: long_name      ! long descriptive name
      character(char_len)     :: units          ! units
      character(char_len)     :: coordinates    ! coordinates
      character(char_len)     :: nftype         ! indicates data type 
      character(4)            :: grid_loc       ! location in grid
      real (rtavg)            :: fill_value     ! _FillValue
      real (rtavg)            :: scale_factor   ! r4 scale factor
      real (r4), dimension(2) :: valid_range    ! min/max
      integer (i4)            :: ndims          ! num dims (2 or 3)
      integer (i4)            :: buf_loc        ! location in buffer
      integer (i4)            :: method         ! method for averaging
      integer (i4)            :: field_loc      ! grid location and field
      integer (i4)            :: field_type     ! type for io, ghost cells
      integer (i4)            :: stream_number  ! identifies tavg_contents "stream"
   end type

!EOP
!BOC

   integer (int_kind), parameter :: &
      max_avail_tavg_fields = 500+21*nt   ! limit on available fields - can
                                          !   be pushed as high as necessary practical
                                          !   (total of all fields in all streams)

   !*** ccsm
   type (tavg_field_desc_ccsm), dimension(max_avail_tavg_fields) :: &
      avail_tavg_fields

   integer (int_kind) ::                &
      num_avail_tavg_fields      = 0,   &! current number of defined fields
      tavg_num_requested_fields,        &! number of fields requested
      tavg_num_contents_lines            ! number of lines in tavg_contents file

!-----------------------------------------------------------------------
!
!  tavg stream information (support for separate tavg output "streams")
!
!-----------------------------------------------------------------------

   type, public :: tavg_stream
      character (char_len) :: infile
      character (char_len) :: outfile
      character (char_len) :: outfile_orig
      character (char_len) :: stream_filestring
      character (char_len) :: fmt_in
      character (char_len) :: fmt_out
      integer (int_kind)   :: freq_iopt
      integer (int_kind)   :: start_iopt
      integer (int_kind)   :: num_requested_fields
      integer (int_kind)   :: field_flag           ! time flag id for writing tavg fields
      integer (int_kind)   :: file_flag            ! time flag id for writing tavg FILES
      integer (int_kind)   :: tavg_offset_year
      integer (int_kind)   :: tavg_offset_month
      integer (int_kind)   :: tavg_offset_day
      integer (int_kind)   :: tavg_num_time_slices
      logical (log_kind)   :: ltavg_on
      logical (log_kind)   :: ltavg_file_is_open
      logical (log_kind)   :: ltavg_fmt_in_nc
      logical (log_kind)   :: ltavg_fmt_out_nc
      logical (log_kind)   :: ltavg_has_offset_date
      logical (log_kind)   :: ltavg_one_time_header
      logical (log_kind)   :: ltavg_first_header
      logical (log_kind)   :: ltavg_qflux_method_on
      real (r8)            :: lower_time_bound
      real (r8)            :: upper_time_bound
      type (io_field_desc), dimension(:), allocatable :: tavg_fields
   end type

   type (tavg_stream), dimension(max_avail_tavg_streams), public :: tavg_streams

   integer (int_kind), public ::  &
      n_tavg_streams               ! actual number of tavg output "streams" requested

   integer (int_kind) ::  &
      nstreams,           &! shorthand name for n_tavg_streams
      qflux_stream         ! stream in which qflux is requested

   logical (log_kind), dimension (max_avail_tavg_fields), public :: &
      ltavg_on      = .false.     ! tavg file output wanted

   character (char_len),dimension(max_avail_tavg_streams) ::    &
      tavg_fmt_out,        &! format (nc or bin) for writing    
      tavg_stream_filestrings     ! output filename string (eg, h2.nday1)

   type (io_dim) ::   &
      i_dim, j_dim,   &! dimension descriptors for horiz dims
      k_dim,          &! dimension descriptor for vert levels (z_t, z_w_top, or z_w_bot grid)
      time_dim         ! dimension descriptor for (unlimited) time dim
 

!-----------------------------------------------------------------------
!
!  buffers for holding running tavg variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      tavg_bufsize_2d,   &    ! size of buffer for 2d fields
      tavg_bufsize_3d         ! size of buffer for 3d fields

   real (rtavg), dimension(:,:,:,:), allocatable :: &
      TAVG_BUF_2D         ! buffer for holding accumulated sums

   real (rtavg), dimension(:,:,:,:,:), allocatable :: &
      TAVG_BUF_3D         ! buffer for holding accumulated sums

   integer (i4), dimension(:), allocatable :: &
      TAVG_BUF_2D_METHOD,  &! method for each requested 2d field
      TAVG_BUF_3D_METHOD    ! method for each requested 3d field

   real (rtavg), dimension (:,:,:), allocatable ::  &
         TAVG_TEMP          ! work array in write_restart
!-----------------------------------------------------------------------
!
!  variables for writing data
!
!-----------------------------------------------------------------------

   !*** field frequency support
   character (char_len), dimension(max_avail_tavg_streams) ::  &
      tavg_freq_opt,    &! choice for frequency of tavg output
      tavg_start_opt,   &! choice for starting averaging
      tavg_fmt_in        ! format (nc or bin) for reading

   integer (i4), dimension(max_avail_tavg_streams) ::  &
      tavg_freq_iopt,   &! frequency option for writing tavg
      tavg_freq,        &! frequency of tavg output
      tavg_start_iopt,  &! start after option
      tavg_start         ! start tavg after tavg_start

  !*** FILE frequency (timeseries) support
   character (char_len), dimension(max_avail_tavg_streams) ::  &
      tavg_file_freq_opt  ! choice for frequency of tavg FILE output 

   integer (i4), dimension(max_avail_tavg_streams) ::  &
      tavg_file_freq_iopt,   &! frequency option for writing tavg FILE
      tavg_file_freq          ! frequency of tavg output FILE creation

   character (char_len) ::  & 
      tavg_contents,    &! filename for choosing fields for output
      tavg_infile,      &! filename for restart input
      tavg_outfile,     &! root filename for tavg output
      tavg_outfile_orig  ! root filename for tavg output (original)

   type (datafile), dimension(max_avail_tavg_streams) :: &
      tavg_file_desc    ! IO file descriptor
 

!-----------------------------------------------------------------------
!
!  scalars
!
!-----------------------------------------------------------------------

   real (r8) ::        &
      dtavg             ! current time step

   character (10) :: &
      beg_date       ! date on which the current accumulated sum
                     ! was started (not the tavg_start date)

!-----------------------------------------------------------------------
!
!  coupled code -- ccsm-related
!
!-----------------------------------------------------------------------
   logical (log_kind) ::          &
      lccsm,                      &
      ldiag_bsf,                  &
      ldiag_gm_bolus,             &! logical for diag_gm_bolus
      lsubmeso,                   &! logical for submesoscale_mixing
      lactive_time_dim,           &
      ltavg_fmt_out_nc,           &! true if netCDF output format
      ltavg_streams_index_present  ! true if using new streams tavg_contents;
                                   ! false if using standard tavg_contents

   logical (log_kind), dimension (:,:,:,:), allocatable ::  &
      MASK_22
 
   integer (int_kind), parameter ::      &
      max_avail_tavg_nstd_fields =  20,  &! limit on available fields
      max_avail_tavg_labels      =  10,  &
      max_num_ccsm_coordinates   =  10,  &
      max_num_ccsm_time_invar    =  50,  &
      max_num_ccsm_scalars       = 100
 
   integer (int_kind) ::     &
      num_avail_tavg_nstd_fields   = 0, &! current number of defined nonstandard fields
      num_avail_tavg_labels        = 0, &! current number of ccsm labels
      num_ccsm_coordinates         = 0, &
      zt_150m_levs

   integer (int_kind), dimension (max_avail_tavg_streams) ::     &
      num_ccsm_time_invar          = 0, &
      num_ccsm_scalars             = 0

   real (rtavg), dimension(:), allocatable, target :: &
      ZT_150m_R   ! single/double precision array
 
   integer (int_kind) ::  &
      tavg_BSF,           &
      tavg_MOC,           &
      tavg_N_HEAT,        &
      tavg_N_SALT

   integer (int_kind) ::  & !indices needed for tavg diagnostics
      tavg_id_WVEL,       &
      tavg_id_VVEL,       &
      tavg_id_WISOP,      &
      tavg_id_VISOP,      &
      tavg_id_WSUBM,      &
      tavg_id_VSUBM,      &
      tavg_id_TEMP,       &
      tavg_loc_WVEL,      &
      tavg_loc_VVEL,      &
      tavg_loc_WISOP,     &
      tavg_loc_VISOP,     &
      tavg_loc_WSUBM,     &
      tavg_loc_VSUBM,     &
      tavg_loc_TEMP       

   integer (int_kind) ::    &
      tavg_id_ADVT,         &
      tavg_id_ADVS ,        &
      tavg_id_VNT,          &
      tavg_id_VNS,          &
      tavg_id_HDIFT,        &
      tavg_id_HDIFS,        &
      tavg_id_ADVT_ISOP,    &
      tavg_id_ADVS_ISOP,    &
      tavg_id_VNT_ISOP,     &
      tavg_id_VNS_ISOP,     &
      tavg_id_ADVT_SUBM,    &
      tavg_id_ADVS_SUBM,    &
      tavg_id_VNT_SUBM,     &
      tavg_id_VNS_SUBM,     &
      tavg_loc_ADVT,        &
      tavg_loc_ADVS ,       &
      tavg_loc_VNT,         &
      tavg_loc_VNS,         &
      tavg_loc_HDIFT,       &
      tavg_loc_HDIFS,       &
      tavg_loc_ADVT_ISOP,   &
      tavg_loc_ADVS_ISOP,   &
      tavg_loc_VNT_ISOP,    &
      tavg_loc_VNS_ISOP,    &
      tavg_loc_ADVT_SUBM,   &
      tavg_loc_ADVS_SUBM,   &
      tavg_loc_VNT_SUBM,    &
      tavg_loc_VNS_SUBM    

   integer (int_kind), dimension(:,:), allocatable    ::  &
      KMU_G            ! k index of deepest grid cell on global U grid

   integer (i4) ::  &
      time_bound_id,&
      moc_id,       &
      n_heat_id,    &
      n_salt_id

   integer (i4), dimension(max_avail_tavg_labels) ::  &
      avail_tavg_labels_id

   type (tavg_field_desc_ccsm) ::  &
      avail_tavg_nstd_fields (max_avail_tavg_nstd_fields),  &
      avail_tavg_labels      (max_avail_tavg_labels)


   type (io_field_desc) ::                          &
      ccsm_coordinates (max_num_ccsm_coordinates,max_avail_tavg_streams),  &
      ccsm_time_invar  (max_num_ccsm_time_invar, max_avail_tavg_streams),  & 
      ccsm_scalars     (max_num_ccsm_scalars,    max_avail_tavg_streams),  &
      time_coordinate  (1,                       max_avail_tavg_streams)

   type (io_dim) ::       &
      z_dim,              &! dimension descriptor for vert (z_t, z_w_top, or z_w_bot grid)
      zt_dim,             &! dimension descriptor for vert (z_t grid)
      zt_150m_dim,        &! dimension descriptor for near-surf vert (z_t grid)
      zw_dim,             &! dimension descriptor for vert (same as z_w_top; keep for backwards compatability)
      zw_dim_top,         &! dimension descriptor for vert (z_w_top grid)
      zw_dim_bot,         &! dimension descriptor for vert (z_w_bot grid)
      tr_dim,             &! dimension descriptor 
      nchar_dim,          &! dimension descriptor for character arrays
      d2_dim,             &! dimension descriptor  
      lat_aux_grid_dim,   &! dimension descriptor  
      moc_z_dim,          &! dimension descriptor  
      moc_comp_dim,       &! dimension descriptor  
      transport_comp_dim, &! dimension descriptor  
      transport_reg_dim    ! dimension descriptor  
 
   type(io_dim)                                          ::  &
      time_bound_dims   (2),                                 &
      io_dims_nstd_ccsm (5,max_avail_tavg_nstd_fields),      &
      io_dims_labels    (2,max_avail_tavg_labels)

   integer, dimension (max_avail_tavg_nstd_fields) ::  &
      ndims_nstd_ccsm
 

   integer (int_kind) ::  &
      tavg_debug  = 0      ! debug level [0,1]  1 ==> messages
 
!-----------------------------------------------------------------------
!
!     variables for local spatial averaging of some time-averaged fields
!
!-----------------------------------------------------------------------

   logical (log_kind), dimension(max_avail_tavg_streams) ::   &
      ltavg_nino_diags,  &! true if ltavg_nino_diags_requested is true,
                          ! TEMP is requested and is in the stream
      ltavg_moc_diags,   &! true if moc_requested is true, all necessary
                          ! fields are requested and are in the same stream
      ltavg_n_heat_trans,&! true if n_heat_trans_requested is true, all necessary
                          ! fields are requested and are in the same stream
      ltavg_n_salt_trans  ! true if n_salt_trans_requested is true, all necessary
                          ! fields are requested and are in the same stream

   logical (log_kind)  ::   &
      ltavg_nino_diags_requested      ! namelist 

   integer (int_kind), parameter ::  &
      n_reg_0D = 4                    ! number of regions

   real (r8), dimension(:), allocatable :: &
      SAVG_0D,                             &! local- and time-mean value 
      SAVG_0D_AREA                          ! area of the region

   real (r8), dimension(:,:,:,:), allocatable ::  &
      SAVG_0D_MASK                          ! mask for the region, i.e. 0 and 1 mean
                                            ! outside and inside the region, respectively

   character (char_len), dimension(:), allocatable ::  &
      SAVG_0D_NAME                          ! name of the region

   integer (int_kind) :: &
      timer_write_std,   &
      timer_write_nstd,  &
      timer_tavg_ccsm_diags_bsf, &
      timer_tavg_ccsm_diags_moc, &
      timer_tavg_ccsm_diags_trans

   character (char_len) :: exit_string

!EOC
!***********************************************************************

 contains

!***********************************************************************
!EOP
! !IROUTINE: init_tavg
! !INTERFACE:

 subroutine init_tavg

! !DESCRIPTION:
!  This routine initializes tavg options and reads in contents file to
!  determine which fields for which the user wants time-averaged data.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   save

   integer (POP_i4) ::     &
      errorCode

   integer (i4) ::         &
      n,                   &! dummy index
      nn,                  &! dummy index
      ns,                  &! dummy index for tavg streams
      i,ip1,j,k,           &! dummy indices
      iblock,              &! local block index
      loc,                 &! location of field in buffer
      nu,                  &! unit for contents input file
      cindex,              &! character index for manipulating strings
      nml_error,           &! namelist i/o error flag
      contents_error        ! error flag for contents file read

   logical (log_kind), dimension(max_avail_tavg_streams) ::    &
      ltavg_one_time_header,&! T if tavg file contains full header info only on first file of run
      ltavg_first_header     ! T if this is the first header of the run

   integer (int_kind) ::   &
      max_days,            &! maximum number of days per month in a year
      duplicate,           &! number of duplicate tavg_contents entries
      ignored,             &! number of ignored tavg_contents entries
      id_temp

   logical (log_kind) ::    &
      reject,               &! true if duplicate tavg_contents entry
      skip,                 &! true if tavg_contents entry is to be skipped (stream = 0)
      ltavg_ignore_extra_streams !ignores tavg_contents streams > n_tavg_streams
                                 ! and allows model to continue running

   logical (log_kind), dimension(max_avail_tavg_streams) ::    &
      ltavg_has_offset_date  ! T if tavg time-flag has an offset date

   integer (int_kind), dimension(max_avail_tavg_streams) ::    &
      tavg_offset_years,     &! tavg-flag offset year 
      tavg_offset_months,    &! tavg-flag offset month
      tavg_offset_days        ! tavg-flag offset day

   character (char_len) ::  &
      char_temp              ! temporary for manipulating fields

   character (3) :: char_ns

   character (64), parameter :: &
      freq_fmt = "('stream #',i3, ': tavg diagnostics every ',i6,a8)"

   character (64), parameter :: &
      start_fmt = "('stream #',i3, ': tavg sums accumulated starting at ',a5,i8)"

   character (char_len),dimension(:),allocatable ::  &
      tavg_contents_request   ! list of all fields requested in tavg_contents file

   type (block) ::        &
      this_block          ! block information for current block

   namelist /tavg_nml/ tavg_freq_opt, tavg_freq, tavg_infile,                  &
                       tavg_outfile, tavg_contents, tavg_start_opt,            &
                       tavg_start, tavg_fmt_in, tavg_fmt_out,                  &
                       tavg_stream_filestrings,                                &
                       ltavg_nino_diags_requested, n_tavg_streams,             &
                       ltavg_ignore_extra_streams,                             &
                       ltavg_streams_index_present, ltavg_has_offset_date,     &
                       tavg_offset_years, tavg_offset_months, tavg_offset_days,&
                       ltavg_one_time_header,                                  &
                       tavg_file_freq_opt, tavg_file_freq

!-----------------------------------------------------------------------
!
!  determine if this is a ccsm coupled run
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   lccsm = registry_match('lcoupled')

!-----------------------------------------------------------------------
!
!  read tavg file output frequency and filenames from namelist
!
!-----------------------------------------------------------------------

   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      write(stdout,blank_fmt)
      write(stdout,'(a12)') ' Tavg:'
      write(stdout,blank_fmt)
      write(stdout,delim_fmt)
      call POP_IOUnitsFlush(POP_stdout); call POP_IOUnitsFlush(stdout)
   endif

   qflux_stream                = 1
   tavg_freq_iopt              = freq_opt_never     ! field frequency
   tavg_freq_opt               = 'never'                
   tavg_freq                   = 100000
   tavg_file_freq_iopt         = freq_opt_never     ! file  frequency
   tavg_file_freq_opt          = 'never'
   tavg_file_freq              = 100000
   tavg_start_iopt             = start_opt_nstep
   tavg_start_opt              = 'nstep'
   tavg_start                  = 0
   tavg_infile                 = 'unknown_tavg_infile'
   tavg_fmt_in                 = 'nc'
   tavg_outfile                = 't'
   tavg_fmt_out                = 'nc'
   tavg_stream_filestrings     = char_blank
   tavg_contents               = 'unknown_tavg_contents'
   ltavg_one_time_header       = .false.
   ltavg_first_header          = .true.
   ltavg_streams_index_present = .true.
   ltavg_has_offset_date       = .false.
   ltavg_ignore_extra_streams  = .false.

   avail_tavg_fields(:)%stream_number    = -999
   tavg_streams(:)%infile                = 'unknown_tavg_infile'
   tavg_streams(:)%outfile               = 't'
   tavg_streams(:)%outfile_orig          = 't'
   tavg_streams(:)%stream_filestring     = char_blank
   tavg_streams(:)%fmt_in                = 'nc'
   tavg_streams(:)%fmt_out               = 'nc'
   tavg_streams(:)%freq_iopt             = freq_opt_never
   tavg_streams(:)%start_iopt            = start_opt_nstep
   tavg_streams(:)%num_requested_fields  = 0
   tavg_streams(:)%tavg_num_time_slices  = 0
   tavg_streams(:)%ltavg_on              = .false.
   tavg_streams(:)%ltavg_file_is_open    = .false.
   tavg_streams(:)%ltavg_has_offset_date = .false.
   tavg_streams(:)%ltavg_one_time_header = .false.
   tavg_streams(:)%ltavg_first_header    = .true.
   tavg_streams(:)%ltavg_qflux_method_on = .false.


   if (registry_match ('init_time1')) then
      tavg_offset_years  = iyear0
      tavg_offset_months = imonth0
      tavg_offset_days   = iday0
   else
      tavg_offset_years  = 1
      tavg_offset_months = 1
      tavg_offset_days   = 2
   endif    

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=tavg_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      exit_string = 'FATAL ERROR: reading tavg_nml'
      call document ('init_tavg', exit_string)
      call exit_POP (sigAbort, exit_string, out_unit=stdout)
   endif

   call broadcast_scalar(n_tavg_streams, master_task)
   if (n_tavg_streams > max_avail_tavg_streams) then
      exit_string = 'FATAL ERROR: reading tavg_nml -- too many tavg streams'
      call document ('init_tavg', exit_string)
      call exit_POP (sigAbort,exit_string,out_unit=stdout)
   endif
  
   if (n_tavg_streams > 9) then
      !*** filename creation will fail; contents read will fail; do you really need more than 9 streams?
      exit_string = 'FATAL ERROR: reading tavg_nml -- tavg streams must be <= 9'
      call document ('init_tavg', exit_string)
      call exit_POP (sigAbort,exit_string,out_unit=stdout)
   endif
 
   nstreams = n_tavg_streams

   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,'(a28)') ' tavg_nml namelist settings:'
      write(stdout,blank_fmt)
      write(stdout,tavg_nml)
      write(stdout,blank_fmt)
      write(stdout,*) ' There will be ', nstreams, ' tavg output streams created during this run'
      write(stdout,blank_fmt)
      write(stdout,delim_fmt)
      call POP_IOUnitsFlush(POP_stdout); call POP_IOUnitsFlush(stdout)
   endif

   if (my_task == master_task) then
     do ns=1,nstreams
      
      !*** field frequency
      select case (tavg_freq_opt(ns))
      case ('never')
         tavg_freq_iopt(ns) = freq_opt_never
         write(stdout,'(a,i3)') 'tavg diagnostics disabled for stream ', ns
      case ('once')
         tavg_freq_iopt(ns) = freq_opt_once
         write(stdout,'(a,i3,a)') 'stream #', ns, 'tavg fields written once (time-invariant fields)'
      case ('nyear')
         tavg_freq_iopt(ns) = freq_opt_nyear
         write(stdout,freq_fmt) ns, tavg_freq(ns),' years  '
      case ('nmonth')
         tavg_freq_iopt(ns) = freq_opt_nmonth
         write(stdout,freq_fmt) ns, tavg_freq(ns),' months '
      case ('nday')
         tavg_freq_iopt(ns) = freq_opt_nday
         write(stdout,freq_fmt) ns, tavg_freq(ns),' days   '
      case ('nhour')
         tavg_freq_iopt(ns) = freq_opt_nhour
         write(stdout,freq_fmt) ns, tavg_freq(ns),' hours  '
      case ('nsecond')
         tavg_freq_iopt(ns) = freq_opt_nsecond
         write(stdout,freq_fmt) ns, tavg_freq(ns),' seconds'
      case ('nstep')
         tavg_freq_iopt(ns) = freq_opt_nstep
         write(stdout,freq_fmt) ns, tavg_freq(ns),' steps  '
      case default
         tavg_freq_iopt(ns) = -1000
      end select

      if (tavg_freq_iopt(ns) /= freq_opt_never) then
         select case (tavg_start_opt(ns))
         case ('nstep')
            tavg_start_iopt(ns) = start_opt_nstep
            write(stdout,start_fmt) ns, 'step ', tavg_start(ns)
         case ('nday')
            tavg_start_iopt(ns) = start_opt_nday
            write(stdout,start_fmt) ns, 'day  ', tavg_start(ns)
         case ('nyear')
            tavg_start_iopt(ns) = start_opt_nyear
            write(stdout,start_fmt) ns, 'year ', tavg_start(ns)
         case ('date')
            tavg_start_iopt(ns) = start_opt_date
            write(stdout,start_fmt) ns, '     ', tavg_start(ns)
         case default
            tavg_start_iopt(ns) = -1000
         end select
      endif

      if (tavg_freq_iopt(ns) /= freq_opt_never) then
        !*** FILE frequency
        select case (tavg_file_freq_opt(ns))
        case ('never')
           tavg_file_freq_iopt(ns) = freq_opt_never
           write(stdout,'(a,i3)') 'tavg diagnostics disabled for stream ', ns
        case ('once')
           tavg_file_freq_iopt(ns) = freq_opt_once
           write(stdout,'(a,i3,a)') 'stream #', ns, 'tavg file written once (time-invariant fields)'
        case ('nyear')
           tavg_file_freq_iopt(ns) = freq_opt_nyear
           write(stdout,freq_fmt) ns, tavg_file_freq(ns),' years  '
        case ('nmonth')
           tavg_file_freq_iopt(ns) = freq_opt_nmonth
           write(stdout,freq_fmt) ns, tavg_file_freq(ns),' months '
        case ('nday')
           tavg_file_freq_iopt(ns) = freq_opt_nday
           write(stdout,freq_fmt) ns, tavg_file_freq(ns),' days   '
        case ('nhour')
           tavg_file_freq_iopt(ns) = freq_opt_nhour
           write(stdout,freq_fmt) ns, tavg_file_freq(ns),' hours  '
        case ('nsecond')
           tavg_file_freq_iopt(ns) = freq_opt_nsecond
           write(stdout,freq_fmt) ns, tavg_file_freq(ns),' seconds'
        case ('nstep')
           tavg_file_freq_iopt(ns) = freq_opt_nstep
           write(stdout,freq_fmt) ns, tavg_file_freq(ns),' steps  '
        case default
           tavg_file_freq_iopt(ns) = -1000
        end select
      endif

      call POP_IOUnitsFlush(POP_stdout); call POP_IOUnitsFlush(stdout)
     enddo ! ns

   endif

   call broadcast_array(tavg_freq_iopt, master_task)

   if (ANY(tavg_freq_iopt == -1000)) then
      exit_string = 'FATAL ERROR: unknown option for tavg FIELD frequency'
      call document ('init_tavg', exit_string)
      call exit_POP (sigAbort,exit_string,out_unit=stdout)
   endif

   if (ANY(tavg_start_iopt == -1000)) then
      exit_string = 'FATAL ERROR: unknown option for tavg start option'
      call document ('init_tavg', exit_string)
      call exit_POP (sigAbort,exit_string,out_unit=stdout)
   endif

   call broadcast_array(tavg_file_freq_iopt, master_task)
   if (ANY(tavg_file_freq_iopt == -1000)) then
      exit_string = 'FATAL ERROR: unknown option for tavg FILE frequency'
      call document ('init_tavg', exit_string)
      call exit_POP (sigAbort,exit_string,out_unit=stdout)
   endif

   call broadcast_scalar(tavg_contents,                  master_task)
   call broadcast_scalar(tavg_infile,                    master_task)
   call broadcast_scalar(tavg_outfile,                   master_task)
   call broadcast_scalar(ltavg_ignore_extra_streams,     master_task)
   call broadcast_scalar(ltavg_nino_diags_requested,     master_task)
   call broadcast_scalar(ltavg_streams_index_present,    master_task)

   do ns=1,nstreams
      call broadcast_scalar(tavg_freq(ns),               master_task)
      call broadcast_scalar(tavg_file_freq(ns),          master_task)
      call broadcast_scalar(tavg_start_iopt(ns),         master_task)
      call broadcast_scalar(tavg_start(ns),              master_task)
      call broadcast_scalar(tavg_fmt_in(ns),             master_task)
      call broadcast_scalar(tavg_fmt_out(ns),            master_task)
      call broadcast_scalar(tavg_stream_filestrings(ns), master_task)
      call broadcast_scalar(ltavg_has_offset_date(ns),   master_task) 
      call broadcast_scalar(tavg_offset_years(ns),       master_task) 
      call broadcast_scalar(tavg_offset_months(ns),      master_task)
      call broadcast_scalar(tavg_offset_days(ns),        master_task)
      call broadcast_scalar(ltavg_one_time_header(ns),   master_task)
    enddo ! ns

!-----------------------------------------------------------------------
!
!  initialize tavg_stream
!
!-----------------------------------------------------------------------


   do ns=1,nstreams
      if (ns == 1) then
        tavg_streams(ns)%infile  = trim(tavg_infile)
        tavg_streams(ns)%outfile = trim(tavg_outfile)
      else
        char_temp = trim(tavg_stream_filestrings(ns))
        if (char_temp .eq. ' ') then
          write(tavg_streams(ns)%infile, *)trim(tavg_infile),'.',ns
          write(tavg_streams(ns)%outfile,'(a,i1)')trim(tavg_outfile),ns
        else
          write(tavg_streams(ns)%infile, '(a,a,a)')trim(tavg_infile), '.',trim(char_temp)
          write(tavg_streams(ns)%outfile,'(a,a,a)')trim(tavg_outfile),'.',trim(char_temp)
        endif
      endif
      tavg_streams(ns)%outfile_orig           = tavg_streams(ns)%outfile
      tavg_streams(ns)%stream_filestring      = tavg_stream_filestrings(ns)
      tavg_streams(ns)%fmt_in                 = tavg_fmt_in (ns)
      tavg_streams(ns)%fmt_out                = tavg_fmt_out(ns)
      tavg_streams(ns)%freq_iopt              = tavg_freq_iopt(ns)
      tavg_streams(ns)%start_iopt             = tavg_start_iopt(ns)
      tavg_streams(ns)%tavg_offset_year       = tavg_offset_years(ns)
      tavg_streams(ns)%tavg_offset_month      = tavg_offset_months(ns)
      tavg_streams(ns)%tavg_offset_day        = tavg_offset_days(ns)
      tavg_streams(ns)%ltavg_has_offset_date  = ltavg_has_offset_date(ns)
      tavg_streams(ns)%ltavg_one_time_header  = ltavg_one_time_header(ns)

!-----------------------------------------------------------------------
!
!  initialize time flags for writing tavg FIELDS
!                                         
!-----------------------------------------------------------------------

      write(char_temp,1100) 'tavg',ns ;  1100 format (a,i1)

      if (ltavg_has_offset_date(ns)) then
        call init_time_flag (trim(char_temp),                       &
                             tavg_streams(ns)%field_flag,           &
                             owner        = 'init_tavg',            &
                             default      =.false.,                 &
                             freq_opt     = tavg_freq_iopt(ns),     &
                             freq         = tavg_freq(ns),          &
                             offset_year  = tavg_offset_years(ns),  &
                             offset_month = tavg_offset_months(ns), & 
                             offset_day   = tavg_offset_days(ns)    ) 
      else
        call init_time_flag(trim(char_temp),                    &
                             tavg_streams(ns)%field_flag,       &
                             owner        = 'init_tavg',        &
                             default      =.false.,             &
                             freq_opt     = tavg_freq_iopt(ns), &
                             freq         = tavg_freq(ns)       )
      endif

!-----------------------------------------------------------------------
!
!  initialize time flags for writing tavg FILES
!                                         
!-----------------------------------------------------------------------

      write(char_temp,1100) 'tavg_file',ns 

      if (ltavg_has_offset_date(ns)) then
        call init_time_flag (trim(char_temp),                       &
                             tavg_streams(ns)%file_flag,            &
                             owner        = 'init_tavg',            &
                             default      =.false.,                 &
                             freq_opt     = tavg_freq_iopt(ns),     &
                             freq         = tavg_freq(ns),          &
                             offset_year  = tavg_offset_years(ns),  &
                             offset_month = tavg_offset_months(ns), & 
                             offset_day   = tavg_offset_days(ns)    ) 
      else
        call init_time_flag(trim(char_temp),                        &
                            tavg_streams(ns)%file_flag,             &
                            owner        = 'init_tavg',             &
                            default      =.false.,                  &
                            freq_opt     = tavg_file_freq_iopt(ns), &
                            freq         = tavg_file_freq(ns)       )
      endif

     if (trim(tavg_fmt_in(ns)) == 'nc') then
        tavg_streams(ns)%ltavg_fmt_in_nc = .true.
     else
        tavg_streams(ns)%ltavg_fmt_in_nc = .false.
     endif

     if (trim(tavg_fmt_out(ns)) == 'nc') then
        tavg_streams(ns)%ltavg_fmt_out_nc = .true.
     else
        tavg_streams(ns)%ltavg_fmt_out_nc = .false.
     endif

   enddo ! nstreams


!-----------------------------------------------------------------------
!
!  if all tavg frequencies are set to 'never', return
!  if not, confirm that the tavg_contents file exists
!
!-----------------------------------------------------------------------
   if (ALL(tavg_freq_iopt == freq_opt_never)) then
      exit_string = 'NOTE: no tavg output has been requested; return control to calling routine'
      call document ('init_tavg', exit_string)
      RETURN
   else
     if (tavg_contents(1:8) == 'unknown_') then
      exit_string = 'FATAL ERROR: tavg option is active, but tavg_contents file = '//trim(tavg_contents)
      call document ('init_tavg', exit_string)
      call exit_POP (sigAbort,exit_string,out_unit=stdout)
     endif
   endif


!-----------------------------------------------------------------------
!
!  define time-averaged BSF field; it may or may not be requested in 
!  the contents file
!
!  All calls to define_tavg_field must be completed prior to calling
!    request_tavg_field
!-----------------------------------------------------------------------

   call define_tavg_field(tavg_BSF,'BSF',2,                                 &
                          long_name='Diagnostic barotropic streamfunction', &
                          units='Sv', grid_loc='2220',                      &
                          coordinates='ULONG ULAT time')


   tavg_bufsize_2d = 0
   tavg_bufsize_3d = 0

   call get_unit(nu)

!-----------------------------------------------------------------------
!
!  count the number of lines in tavg_contents file
!
!-----------------------------------------------------------------------

   call tavg_count_contents(tavg_num_contents_lines,tavg_contents)

   allocate (tavg_contents_request(tavg_num_contents_lines))

   if (my_task == master_task) then
      write(stdout,*) '(init_tavg) total number of lines in tavg_contents file = ', tavg_num_contents_lines
      call POP_IOUnitsFlush(POP_stdout); call POP_IOUnitsFlush(stdout)
      
      open(nu, file=tavg_contents, status='old')
      write(stdout,'(a38)') 'tavg diagnostics requested for fields:'
      call POP_IOUnitsFlush(POP_stdout); call POP_IOUnitsFlush(stdout)
   endif

!-----------------------------------------------------------------------
!
!  read contents files to determine which fields to dump; eliminate 
!    duplicate requests and any field with a leading "0"
!
!-----------------------------------------------------------------------

   contents_error = 0
   duplicate = 0
   ignored = 0

   read_tavg_contents_loop: do n=1,tavg_num_contents_lines

      if (my_task == master_task) then
         char_temp = char_blank

         if (ltavg_streams_index_present) then
           read(nu,'(a3,a)',iostat=contents_error) char_ns,char_temp
           if (char_ns(1:1) == '#' .or. char_ns(1:1) == '!' ) then
             ns = 0
           else
             read(char_ns(1:1), '(i1)') ns
           endif 
         else
           read(nu,'(a)',iostat=contents_error) char_temp
           ns = 1
         endif
      endif ! master_task

      call broadcast_scalar(ns,             master_task)
      call broadcast_scalar(char_temp,      master_task)
      call broadcast_scalar(contents_error, master_task)


      !*** error trapping
      if (contents_error /= 0) then
        exit_string = 'FATAL ERROR: reading tavg contents'
        call document ('init_tavg', exit_string)
        call exit_POP (sigAbort,exit_string,out_unit=stdout)
      endif

      if (ns < 0 .or. ns > max_avail_tavg_streams) then
        exit_string = 'FATAL ERROR: invalid stream number in tavg_contents file'
        call document ('init_tavg', exit_string)
        call exit_POP (sigAbort,exit_string,out_unit=stdout)
      endif

      if (ns > n_tavg_streams) then
        if (ltavg_ignore_extra_streams) then
           exit_string = 'WARNING: you have requested a stream number > n_tavg_streams in tavg_contents file'
           call document ('init_tavg', exit_string)
        else
           write(exit_string,'(a,i2,2x,a,2x,a)') ' stream requested = ', ns, 'field requested = ', char_temp
           call document ('init_tavg', exit_string)
           exit_string = 'FATAL ERROR: you have requested a stream number > n_tavg_streams in tavg_contents file'
           call document ('init_tavg', exit_string)
           call exit_POP (sigAbort,exit_string,out_unit=stdout)
        endif
      endif

    
      !*** look for inactive or duplicate field requests

      if (ns == 0) then
        !*** inactive field; skip the rest of this block
        ignored = ignored + 1
        call document ('init_tavg', 'inactive tavg_contents field: ',  trim(char_temp))
      else
        char_temp = adjustl(char_temp)
        cindex = index(char_temp,' ')
        char_temp(cindex:) = ' '

        if (ltavg_streams_index_present) then
          if (trim(char_temp) == 'QFLUX') qflux_stream = ns
        else
          if (trim(char_temp) == 'QFLUX') qflux_stream = 1
          ns = 1
        endif

        !*** reject any duplicate requests
        reject = .false.
        if (n == 1) then
          tavg_contents_request(n) = trim(char_temp)
        else
          dup_loop: do nn=1,n-1
            if (trim(char_temp) == trim(tavg_contents_request(nn))) then
               duplicate = duplicate+1
               reject = .true.
               exit dup_loop
            endif
          enddo dup_loop

          if (.not. reject) then
              tavg_contents_request(n) = trim(char_temp)
              if (my_task == master_task) then
                write(stdout,*)  ns, '  ',trim(char_temp)
                call POP_IOUnitsFlush(POP_stdout);call POP_IOUnitsFlush(stdout)
              endif
          else
              tavg_contents_request(n) = 'duplicate'
          endif ! reject
        endif !  n == 1

        if (reject) then
          call document ('init_tavg', 'duplicate tavg_contents field: ',  trim(char_temp))
        else
          !*** activate requested tavg fields
          !    determines tavg_bufsize_2d, tavg_bufsize_3d

          call request_tavg_field(trim(char_temp),ns)  
          !*** now that requested tavg field is successfully requested,
          !    store stream information

          id_temp = tavg_id(trim(char_temp))
          avail_tavg_fields(id_temp)%stream_number = ns

          tavg_streams(ns)%num_requested_fields = tavg_streams(ns)%num_requested_fields + 1

        endif ! reject

     endif  ! ignored field ("0")
   enddo read_tavg_contents_loop

   call release_unit(nu)

   !*** adjust tavg_num_requested_fields to account for rejected duplicate requests:
   tavg_num_requested_fields = tavg_num_contents_lines - duplicate - ignored

   call document ('init_tavg', 'Total number of tavg fields requested ',  tavg_num_requested_fields)

   !*** ensure that there are no empty streams
   do ns=1,nstreams
     if (tavg_streams(ns)%num_requested_fields == 0) then
       call document ('init_tavg', 'ERROR: Empty stream number', ns)
       exit_string = 'FATAL ERROR: Empty stream'
       call document ('init_tavg', exit_string)
       call exit_POP (sigAbort,exit_string,out_unit=stdout)
     endif
   enddo

   if (tavg_num_requested_fields > 0) then
   !*** allocate and initialize running tavg buffers

     if (tavg_bufsize_2d > 0) then
      allocate(TAVG_BUF_2D(nx_block,ny_block,nblocks_clinic,tavg_bufsize_2d))
      allocate(TAVG_BUF_2D_METHOD(tavg_bufsize_2d))
      TAVG_BUF_2D = c0
     endif

     if (tavg_bufsize_3d > 0) then
      allocate(TAVG_BUF_3D(nx_block,ny_block,km,nblocks_clinic,tavg_bufsize_3d))
      allocate(TAVG_BUF_3D_METHOD(tavg_bufsize_3d))
      TAVG_BUF_3D = c0
     endif

     allocate(TAVG_TEMP(nx_block,ny_block,nblocks_clinic))

     tavg_sum = c0

     call time_stamp('now','ymd',date_string=beg_date)
     do ns=1,nstreams
       if (tavg_freq_iopt(ns) == freq_opt_nstep) &
          write(beg_date,'(i10)') nsteps_total
     enddo ! ns


     do n = 1,num_avail_tavg_fields  ! check all available fields
       loc = abs(avail_tavg_fields(n)%buf_loc)
       if (loc /= 0) then  ! field is actually requested and in buffer
          if (avail_tavg_fields(n)%ndims == 2) then
             TAVG_BUF_2D_METHOD(loc) = avail_tavg_fields(n)%method
          else if (avail_tavg_fields(n)%ndims == 3) then
             TAVG_BUF_3D_METHOD(loc) = avail_tavg_fields(n)%method
          endif

          !*** determine which streams use tavg_method_qflux
          if (avail_tavg_fields(n)%method == tavg_method_qflux) then
             tavg_streams(avail_tavg_fields(n)%stream_number)%ltavg_qflux_method_on = .true.
          endif
       endif
     end do

     !*** initialize buffers based on requested method
     call tavg_reset_field_all

   endif !tavg_num_requested_fields

  !*** document which streams are using tavg_method_qflux
   if (my_task == master_task) then
      do ns=1,nstreams
         write(stdout,*) '(init_tavg)  tavg_streams(',ns,  &
               ')%ltavg_qflux_method_on = ', tavg_streams(ns)%ltavg_qflux_method_on
      enddo ! ns
   endif


!-----------------------------------------------------------------------
!
!  define dimensions for tavg output files
!
!-----------------------------------------------------------------------

   !*** define dimensions for tavg output files
   !    redefined below if ccsm
   i_dim     = construct_io_dim('i',nx_global)
   j_dim     = construct_io_dim('j',ny_global)
   k_dim     = construct_io_dim('k',km)
   time_dim  = construct_io_dim('time',0,start=1,stop=1) ! used only to set %active=.false.
 
!-----------------------------------------------------------------------
!
!  initialize ccsm NINO and transport diagnostics
!
!-----------------------------------------------------------------------
 
   !*** external diagnostics initialization routines
   call init_lat_aux_grid         
   call init_moc_ts_transport_arrays   
   call init_diag_bsf(set_in_tavg_contents(tavg_BSF), errorCode)

   !*** important: the following logical variables must be set after calls to 
   !      external diagnostics initialization routines
   ldiag_bsf      = registry_match('ldiag_bsf')
   ldiag_gm_bolus = registry_match('diag_gm_bolus')
   lsubmeso       = registry_match('init_submeso')

   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,'(a)') ' Internal tavg Diagnostics Control Variables:'
      write(stdout,blank_fmt)
      write(stdout,*) '   (init_tavg) ldiag_bsf      = ', ldiag_bsf
      write(stdout,*) '   (init_tavg) ldiag_gm_bolus = ', ldiag_gm_bolus
      write(stdout,*) '   (init_tavg) lsubmeso       = ', lsubmeso
      write(stdout,blank_fmt)
      write(stdout,delim_fmt)
      call POP_IOUnitsFlush(POP_stdout); call POP_IOUnitsFlush(stdout)
   endif

   !*** internal diagnostics initialization routines
   call tavg_init_local_spatial_avg
   call tavg_init_moc_diags
   call tavg_init_transport_diags

!-----------------------------------------------------------------------
!
!  ccsm-specific initializations
!
!-----------------------------------------------------------------------
 
   if (lccsm) then

     !*** how many levels have their midpoint shallower than 150m
     zt_150m_levs = count(zt < 150.0e2_r8)
     if (.not. allocated(ZT_150m_R)) &
       allocate(ZT_150m_R(zt_150m_levs))

     !*** define dimensions for tavg output files
     i_dim      = construct_io_dim('nlon',nx_global)
     j_dim      = construct_io_dim('nlat',ny_global)
     zt_dim     = construct_io_dim('z_t',km)
     zt_150m_dim= construct_io_dim('z_t_150m',zt_150m_levs)
     zw_dim     = construct_io_dim('z_w',km)  ! same as zw_dim_top
     zw_dim_top = construct_io_dim('z_w_top',km)
     zw_dim_bot = construct_io_dim('z_w_bot',km)
     time_dim   = construct_io_dim('time',0)    ! "0" ==> unlimited dimension
     tr_dim     = construct_io_dim('tracers',nt)
     nchar_dim  = construct_io_dim('nchar',char_len)
     d2_dim     = construct_io_dim('d2',2)

     if (moc_requested .or. n_heat_trans_requested .or. n_salt_trans_requested) then
       lat_aux_grid_dim  = construct_io_dim('lat_aux_grid',n_lat_aux_grid+1)
       transport_reg_dim = construct_io_dim('transport_reg',n_transport_reg)
       if (moc_requested) then
         moc_z_dim    = construct_io_dim('moc_z',km+1)
         moc_comp_dim = construct_io_dim('moc_comp',n_moc_comp)
       endif
       if (n_heat_trans_requested .or. n_salt_trans_requested) &
         transport_comp_dim = construct_io_dim('transport_comp',n_transport_comp)
     endif ! moc_requested

!-----------------------------------------------------------------------
!
!  set up masking array
!
!-----------------------------------------------------------------------

     allocate (MASK_22(nx_block,ny_block,km,nblocks_clinic))

 
     !--------------------------------------------------------
     !    create MASK_22, layer by layer
     !--------------------------------------------------------
 
 
     do iblock = 1,nblocks_clinic 
     this_block = get_block(blocks_clinic(iblock),iblock)  
       do k=1,km
         MASK_22(:,:,k,iblock) = .false.
         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            if (k > KMT(i  ,j,   iblock)   .and.  &
                k > KMT(i+1,j,   iblock)   .and.  &
                k > KMT(i  ,j+1, iblock)   .and.  &
                k > KMT(i+1,j+1, iblock))   then
                  MASK_22(i,j,k,iblock) = .true.
            endif
         enddo ! i
         enddo ! j
       enddo ! k
     enddo ! iblock
 
   endif ! lccsm
 

!-----------------------------------------------------------------------
!
!  check: time to begin accumulating time averages?
!
!-----------------------------------------------------------------------

   call tavg_set_flag(update_time=.false.)


!-----------------------------------------------------------------------
!
!  finally, read restart files if necessary
!  must do this after i_dim, j_dim, etc are defined
!
!-----------------------------------------------------------------------

   do ns=1,nstreams

     !*** make sure tavg flag is set correctly
     !    evaluate time_flag(tavg_streams(ns)%field_flag)%value via time_to_do
     call eval_time_flag(tavg_streams(ns)%field_flag) 

     if (ltavg_on(ns) .and. ltavg_restart) then
        !*** do not read restart if last restart was at a regular tavg write
        !*** interval (should start new tavg sums in this case)

        if (.not. check_time_flag(tavg_streams(ns)%field_flag)) then
           call read_tavg(ns)
        endif
     endif
   enddo ! ns

!-----------------------------------------------------------------------
!
!  error checking
!
!-----------------------------------------------------------------------

   if (moc_requested) then
     if (.not. set_in_tavg_contents(tavg_id('WVEL')) .or.  &
         .not. set_in_tavg_contents(tavg_id('VVEL')) )     then
        exit_string = 'FATAL ERROR: for moc diagnostics, WVEL and VVEL must be requested in tavg_contents file'
        call document ('init_tavg', exit_string)
        call exit_POP (sigAbort,exit_string,out_unit=stdout)
     endif
   endif

   if ( n_heat_trans_requested .or. n_salt_trans_requested) then
     if (.not. set_in_tavg_contents(tavg_id('ADVT'))   .or. &
         .not. set_in_tavg_contents(tavg_id('ADVS'))   .or. &
         .not. set_in_tavg_contents(tavg_id('VNT'))    .or. &
         .not. set_in_tavg_contents(tavg_id('VNS'))    .or. &
         .not. set_in_tavg_contents(tavg_id('HDIFT'))  .or. &
         .not. set_in_tavg_contents(tavg_id('HDIFS')) ) then
        exit_string = &
        'FATAL ERROR: diag_gm_bolus must have ADVT ADVS VNT VNS HDIFT HDIFS in tavg_contents file'
        call document ('init_tavg', exit_string)
        call exit_POP (sigAbort,exit_string,out_unit=stdout)
     endif
   endif

   if (n_heat_trans_requested .or. n_salt_trans_requested) then
   if (registry_match('diag_gm_bolus')) then
     if (.not. set_in_tavg_contents(tavg_id('ADVT_ISOP'))   .or. &
         .not. set_in_tavg_contents(tavg_id('ADVS_ISOP'))   .or. &
         .not. set_in_tavg_contents(tavg_id('VNT_ISOP'))    .or. &
         .not. set_in_tavg_contents(tavg_id('VNS_ISOP'))  ) then
        exit_string = &
        'FATAL ERROR: diag_gm_bolus must have ADVT_ISOP ADVS_ISOP VNT_ISOP VNS_ISOP in tavg_contents file'
        call document ('init_tavg', exit_string)
        call exit_POP (sigAbort,exit_string,out_unit=stdout)
     endif
   endif
   endif

   deallocate (tavg_contents_request)

!-----------------------------------------------------------------------
!
!  initialize timers
!
!-----------------------------------------------------------------------

   call get_timer(timer_write_std,'TAVG_WRITE_STD', nblocks_clinic, distrb_clinic%nprocs)
   call get_timer(timer_write_nstd,'TAVG_WRITE_NONSTD', nblocks_clinic, distrb_clinic%nprocs)
   if (ldiag_bsf)  &
   call get_timer(timer_tavg_ccsm_diags_bsf,'TAVG_CCSM_DIAGS_BSF', nblocks_clinic, distrb_clinic%nprocs)
   if (moc_requested)  &
   call get_timer(timer_tavg_ccsm_diags_moc,'TAVG_CCSM_DIAGS_MOC', nblocks_clinic, distrb_clinic%nprocs)
   if (n_heat_trans_requested .or. n_salt_trans_requested)  &
   call get_timer(timer_tavg_ccsm_diags_trans,'TAVG_CCSM_DIAGS_TRANS', nblocks_clinic, distrb_clinic%nprocs)

!-----------------------------------------------------------------------
!EOC

 end subroutine init_tavg

!***********************************************************************
!BOP
! !IROUTINE: tavg_set_flag
! !INTERFACE:

 subroutine tavg_set_flag(update_time)

! !DESCRIPTION:
!  This routine checks the time avg option and tavg start condition
!  to see whether tavg sums should be accumulated.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   logical (log_kind), intent(in) :: &
      update_time        ! if F, only sets ltavg_on without advancing time interval

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

 integer (int_kind) ::   &
    n,                   &! loop index
    n1, n2,              &! loop bounds
    ns                    ! loop index

!-----------------------------------------------------------------------
!
!  if tavg requested and tavg not already turned on, check to see
!  if it is time to start time averaging
!
!-----------------------------------------------------------------------

   do ns=1,nstreams

   if (.not. ltavg_on(ns) .and. tavg_freq_iopt(ns) /= freq_opt_never) then

      ltavg_on(ns) = time_to_start(tavg_start_iopt(ns), tavg_start(ns))
      call time_stamp('now','ymd',date_string=beg_date)
      if (tavg_freq_iopt(ns) == freq_opt_nstep) write(beg_date,'(i10)') nsteps_total

      !*** if it is time to start, make sure requested fields
      !*** get triggered by the requested function

      if (ltavg_on(ns)) then
         do n=1,num_avail_tavg_fields
            if (avail_tavg_fields(n)%buf_loc < 0) &
                avail_tavg_fields(n)%buf_loc =    &
                abs(avail_tavg_fields(n)%buf_loc)
         end do
         tavg_streams(ns)%lower_time_bound = tday00
      endif

   endif

   tavg_streams(ns)%ltavg_on = ltavg_on(ns)

!-----------------------------------------------------------------------
!
!  setup time step and total integrated time for time average
!  adjust for averaging timesteps: if this is an averaging timestep,
!  the past values only contribute for 1/4 of a step and the
!  values for the step just before an averaging timestep contribute
!  for 1 1/4 steps.
!
!-----------------------------------------------------------------------

   if (tavg_streams(ns)%ltavg_on .and. update_time) then
      if (avg_ts .or. back_to_back) then
         dtavg = p5*dtt
      else
         dtavg = dtt
      endif

      tavg_sum(ns) = tavg_sum(ns) + dtavg
   endif
   enddo ! ns

!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_set_flag

!***********************************************************************
!BOP
! !IROUTINE: write_tavg
! !INTERFACE:

 subroutine write_tavg(restart_type)

! !DESCRIPTION:
!  This routine writes requested tavg fields to a file.  The fields are
!  normalized by the time interval before writing.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) ::  &
      restart_type           ! tells tavg whether to write restart

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
                       
   real (r8), dimension(1) ::  &
      TIME1D  

   integer (int_kind) ::  &
      ndims

   integer (i4) ::  &
      nu,           &! i/o unit for output file
      iblock,       &! dummy block index
      nfield,       &! dummy field index
      nindex,       &! dummy field index
      loc,          &! buffer location for field
      io_phase,     &!'define' or 'write'
      k,n,          &! indices
      nn,           &! index
      ns,           &! index
      i,j,          &! indices
      nstd_field_id,&
      field_id,     &!temporary id 
      field_counter,&
      method_integer

   character (char_len) ::  &
      string,               &! dummy character string
      file_suffix,          &! suffix to append to tavg file name
      hist_string,          &! string containing file history
      tavg_filename,        &! filename for tavg data
      tavg_pointer_file,    &! filename for pointer file containing
                             !   location/name of last restart file
      tavg_fmt_in,          &! format (nc or bin) for reading
      ns_temp,              &
      method_string          ! cell_methods string

   character (char_len),save ::  &
      cf_conventions ='CF-1.0; http://www.cgd.ucar.edu/cms/eaton/netcdf/CF-current.htm'

   character (8) ::  &
      date_created   ! string with (real) date this file created

   character (10) ::  &
      time_created   ! string with (real) date this file created

   logical (log_kind) ::     &
      ltavg_write_reg,       &! time to write regular tavg file
                              !   and reset time averages
      ltavg_write_rest,      &! time to write restart tavg file
      time_to_close           ! time to close file

   type (io_field_desc)  ::  &
      qflux_field

   type (block) ::        &
      this_block          ! block information for current block

   save

!-----------------------------------------------------------------------
!
!  loop over each tavg stream
!
!-----------------------------------------------------------------------

   do ns=1,nstreams


   ltavg_fmt_out_nc = tavg_streams(ns)%ltavg_fmt_out_nc

!-----------------------------------------------------------------------
!
!  is it time to write a file? (regular or restart)
!
!    The variable restart_type is set every timestep in write_restart.
!    If it is time to write a *restart* file, the restart_type
!    will have the value {restart, even, odd, end}
!
!-----------------------------------------------------------------------

   ltavg_write_reg  = .false.
   ltavg_write_rest = .false.


   if (tavg_streams(ns)%ltavg_on) then

     !*** time to write a regular tavg file?

     ltavg_write_reg = check_time_flag(tavg_streams(ns)%field_flag)

     !*** time to write a restart tavg file?
     !    if it is time to write a restart file, but not time to write a
     !      *regular* tavg file, then write a *restart* tavg file

     if (trim(restart_type) /= 'none' .and. .not. ltavg_write_reg) then
       ltavg_write_rest = .true.
       ! do not write to previously open file -- write to restart instead
       tavg_streams(ns)%ltavg_file_is_open = .false.
     endif

     !*** do not write a restart tavg file for a one-time stream

     if (tavg_streams(ns)%freq_iopt == freq_opt_once) then
       ltavg_write_rest = .false.
     endif

   endif !tavg_streams(ns)%ltavg_on

   if (ltavg_write_reg .and. ltavg_write_rest) then
     exit_string = 'FATAL ERROR: cannot have both regular and restart write'
     call document ('write_tavg', exit_string)
     call exit_POP (sigAbort,exit_string,out_unit=stdout)
   endif

!-----------------------------------------------------------------------
!
!  do the rest only if it is time to write a regular,restart, or one-time tavg file
!
!-----------------------------------------------------------------------

   if (ltavg_write_reg .or. ltavg_write_rest) then

!-----------------------------------------------------------------------
!
!  is time dimension active?
!     yes if ltavg_write_reg, netCDF, and using ccsm conventions
!     no if ltavg_write_rest
!     no if binary
!
!-----------------------------------------------------------------------

   if (lccsm .and. ltavg_fmt_out_nc .and. ltavg_write_reg) then
     time_dim%active = .true.
    !** support for timeseries output -- update time_dim 
    tavg_streams(ns)%tavg_num_time_slices = tavg_streams(ns)%tavg_num_time_slices + 1
    time_dim%start = tavg_streams(ns)%tavg_num_time_slices
    time_dim%stop  = tavg_streams(ns)%tavg_num_time_slices
   else
     time_dim%active = .false.
   endif


!-----------------------------------------------------------------------
!
!     compute and print global averages of tavg fields prior to normalizing
!
!-----------------------------------------------------------------------

      call tavg_global (ns)

!-----------------------------------------------------------------------
!
!     normalize time averages
!
!-----------------------------------------------------------------------

      call tavg_norm_field_all ('normalize',ns)
 
!-----------------------------------------------------------------------
!
!  compute ccsm diagnostcs from tavg quantities (after normalizing)
!    (on/off controls are internal to each subroutine)
!
!-----------------------------------------------------------------------

      if (ltavg_write_reg) then
       call tavg_bsf_diags (ns)         ! barotropic stream function
       call tavg_moc_diags (ns)         ! MOC diagnostics
       call tavg_transport_diags (ns)   ! northward heat/salt transport diagnostics
       call tavg_local_spatial_avg (ns) ! compute NINO diagnostics
      endif

!-----------------------------------------------------------------------
!
!  construct io_field descriptors for each 2D and 3D field
!  if writing a regular tavg file in netcdf format, mask the fields
!
!-----------------------------------------------------------------------


     if (.not. allocated(tavg_streams(ns)%tavg_fields)) then
       allocate(tavg_streams(ns)%tavg_fields(tavg_streams(ns)%num_requested_fields))
     endif

     field_counter = 0

     do nfield = 1,num_avail_tavg_fields  ! check all available fields

        loc = avail_tavg_fields(nfield)%buf_loc
        
        !*** continue only if field is requested, in buffer, and in this stream
        if (tavg_requested(nfield)) then
          if (tavg_in_this_stream(nfield,ns)) then 

          field_counter = field_counter + 1

          if (field_counter > tavg_streams(ns)%num_requested_fields ) then 
            exit_string = 'FATAL ERROR: too many fields requested'
            call document ('write_tavg', exit_string)
            call exit_POP (sigAbort,exit_string,out_unit=stdout)
          endif

          if (tavg_streams(ns)%ltavg_file_is_open) then
            field_id = tavg_streams(ns)%tavg_fields(field_counter)%id
          else
            field_id = 0
          endif

          if (avail_tavg_fields(nfield)%ndims == 2) then

            if (ltavg_fmt_out_nc .and. ltavg_write_reg) then
              call tavg_mask(TAVG_BUF_2D(:,:,:,loc),avail_tavg_fields(nfield),1)
            endif

            tavg_streams(ns)%tavg_fields(field_counter) = construct_io_field(&
                            avail_tavg_fields(nfield)%short_name,           &
                dim1=i_dim, dim2=j_dim,time_dim=time_dim,field_id=field_id, &
                  long_name=avail_tavg_fields(nfield)%long_name,            &
                      units=avail_tavg_fields(nfield)%units    ,            &
                   grid_loc=avail_tavg_fields(nfield)%grid_loc ,            &
                  field_loc=avail_tavg_fields(nfield)%field_loc,            &
                 field_type=avail_tavg_fields(nfield)%field_type,           &
                coordinates=avail_tavg_fields(nfield)%coordinates,          &
                valid_range=avail_tavg_fields(nfield)%valid_range,          &
#ifdef TAVG_R8
                 d2d_array=TAVG_BUF_2D(:,:,:,loc) ) 
#else
                r2d_array=TAVG_BUF_2D(:,:,:,loc) )
#endif

          else if (avail_tavg_fields(nfield)%ndims == 3) then

            if (ltavg_fmt_out_nc) then
              select case (trim(avail_tavg_fields(nfield)%grid_loc(4:4)))
                case('1')
                  z_dim = zt_dim
                case('2')
                  z_dim = zw_dim_top
                case('3')
                  z_dim = zw_dim_bot
                case('4')
                  z_dim = zt_150m_dim
              end select
            else
              z_dim = k_dim
            endif ! lccsm

            if (ltavg_fmt_out_nc .and. ltavg_write_reg) then
              do k=1,km
                 TAVG_TEMP(:,:,:)=TAVG_BUF_3D(:,:,k,:,loc)
                 call tavg_mask(TAVG_TEMP(:,:,:),avail_tavg_fields(nfield),k)
                 TAVG_BUF_3D(:,:,k,:,loc)=TAVG_TEMP(:,:,:)
               enddo 
            endif

            tavg_streams(ns)%tavg_fields(field_counter) = construct_io_field(&
                           avail_tavg_fields(nfield)%short_name,                        &
               dim1=i_dim, dim2=j_dim, dim3=z_dim, time_dim=time_dim,field_id=field_id, &
                 long_name=avail_tavg_fields(nfield)%long_name,                         &
                     units=avail_tavg_fields(nfield)%units    ,                         &
                  grid_loc=avail_tavg_fields(nfield)%grid_loc ,                         &
                 field_loc=avail_tavg_fields(nfield)%field_loc,                         &
                field_type=avail_tavg_fields(nfield)%field_type,                        &
               coordinates=avail_tavg_fields(nfield)%coordinates,                       &
               valid_range=avail_tavg_fields(nfield)%valid_range,                       &
#ifdef TAVG_R8
                 d3d_array=TAVG_BUF_3D(:,:,:,:,loc) )
#else
                 r3d_array=TAVG_BUF_3D(:,:,:,:,loc) ) 
#endif
          endif ! 2D/3D test

          if (lccsm .and. ltavg_write_reg) &
             call tavg_add_attrib_io_field_ccsm (tavg_streams(ns)%tavg_fields(field_counter),nfield)
       endif ! tavg_in_this_stream 
       endif ! nfield and stream
    end do !nfield 

!-----------------------------------------------------------------------
!
!  construct time
!
!-----------------------------------------------------------------------
 
      if (tavg_streams(ns)%ltavg_file_is_open) then
        field_id = time_coordinate(1,ns)%id
      else
        field_id = 0
      endif

      call tavg_construct_ccsm_time (field_id,ns)  !time_coordinate

!-----------------------------------------------------------------------
!
!  create regular or restart tavg output filenames 
!
!-----------------------------------------------------------------------

    if (.not. tavg_streams(ns)%ltavg_file_is_open) then
!   ==========================================================================

    tavg_outfile_orig = trim(tavg_streams(ns)%outfile_orig)

    if (ltavg_write_reg) then
      tavg_outfile = trim(tavg_outfile_orig)
      if (lccsm) then
          call tavg_create_suffix_ccsm(file_suffix,ns)
      else
        call tavg_create_suffix(file_suffix,ns)
      endif
    endif

    if (ltavg_write_rest) then
      if (lccsm) then
        call tavg_create_outfile_ccsm(tavg_outfile_orig,tavg_outfile,ns)
           call tavg_create_suffix_ccsm(file_suffix,ns,date_string='ymds')
        string = trim(file_suffix)
      else
        string = trim(runid)
      endif

      select case (trim(restart_type))
        case('even','odd','end')
          file_suffix = trim(string)/&
                                     &/'.'/&
                                     &/trim(restart_type)
        case default
          if (.not. lccsm) then
             call tavg_create_suffix(file_suffix,ns)
             file_suffix = trim(file_suffix)/&
                                             &/'.restart'
          endif
      end select
         
    endif ! ltavg_write_rest

!-----------------------------------------------------------------------
!
!   create data file descriptor
!
!-----------------------------------------------------------------------

    if (my_task.eq.master_task) then
      call date_and_time(date=date_created, time=time_created)
    end if
    call broadcast_scalar(date_created, master_task)
    call broadcast_scalar(time_created, master_task)
    hist_string = char_blank
    write(hist_string,'(a23,a8,1x,a10)') & 
    'POP TAVG file created: ',date_created,time_created

    if (ltavg_fmt_out_nc) then
      tavg_file_desc(ns)  = construct_file(tavg_fmt_out(ns),  &
                          root_name  = trim(tavg_outfile),    &
                          file_suffix= trim(file_suffix),     &
                          title      = trim(runid),           &
                          conventions=trim(cf_conventions),   &
                          record_length = rec_type_real,      &
                          recl_words=nx_global*ny_global)
    else
      tavg_file_desc(ns) = construct_file(tavg_fmt_out(ns),   &
                          root_name  = trim(tavg_outfile),    &
                          file_suffix= trim(file_suffix),     &
                          title      ='POP TAVG file',        &
                          conventions='POP TAVG conventions', &
                          history    = trim(hist_string),     &
                          record_length = rec_type_real,      &
                          recl_words=nx_global*ny_global)
    endif !ltavg_fmt_out_nc

!-----------------------------------------------------------------------
!
!   add scalar fields to file as file attributes
!
!-----------------------------------------------------------------------

    call add_attrib_file(tavg_file_desc(ns), 'tavg_sum'    , tavg_sum(ns))
    call add_attrib_file(tavg_file_desc(ns), 'nsteps_total', nsteps_total)

    if (tavg_streams(ns)%ltavg_qflux_method_on) &
    call add_attrib_file(tavg_file_desc(ns), 'tavg_sum_qflux'  , tavg_sum_qflux(ns))

    if (ltavg_fmt_out_nc .and. ltavg_write_reg) then
      call tavg_add_attrib_file_ccsm (tavg_file_desc(ns)) 
    else
      call add_attrib_file(tavg_file_desc(ns), 'tday'      , tday)
      call add_attrib_file(tavg_file_desc(ns), 'iyear'     , iyear)
      call add_attrib_file(tavg_file_desc(ns), 'imonth'    , imonth)
      call add_attrib_file(tavg_file_desc(ns), 'iday'      , iday)
      call add_attrib_file(tavg_file_desc(ns), 'beg_date'  , beg_date)
    endif

!-----------------------------------------------------------------------
!
!   open output file
!
!-----------------------------------------------------------------------

    call data_set (tavg_file_desc(ns), 'open')

    if (ltavg_fmt_out_nc .and. ltavg_write_reg) then
   
      call tavg_define_time_bounds     (tavg_file_desc(ns))
      call tavg_define_labels_ccsm     (tavg_file_desc(ns),ns)

      !*** construct_io_fields
      call tavg_construct_ccsm_coordinates  (tavg_file_desc(ns),ns)  !ccsm_coordinates
      call tavg_construct_ccsm_time_invar   (tavg_file_desc(ns),ns)  !ccsm_time_invar
      call tavg_construct_ccsm_scalars      (tavg_file_desc(ns),ns)  !ccsm_scalars

      !*** define fields 
      call data_set (tavg_file_desc(ns), 'define', time_coordinate(1,ns))

      do n=1,num_ccsm_coordinates
       call data_set (tavg_file_desc(ns), 'define', ccsm_coordinates(n,ns))
      enddo


      do n=1,num_ccsm_time_invar(ns)
       call data_set (tavg_file_desc(ns), 'define', ccsm_time_invar(n,ns))
      enddo

      do n=1,num_ccsm_scalars(ns)
       call data_set (tavg_file_desc(ns), 'define', ccsm_scalars(n,ns))
      enddo
    endif  !ltavg_fmt_out_nc .and. ltavg_write_reg

    field_counter = 0

    do nfield = 1,num_avail_tavg_fields
    !*** define only if field is requested, in buffer, and in this stream
    if (tavg_requested(nfield)) then
      if (tavg_in_this_stream(nfield,ns)) then
        field_counter = field_counter + 1
        call data_set (tavg_file_desc(ns), 'define', tavg_streams(ns)%tavg_fields(field_counter))
      endif
    endif
    enddo ! nfield

!-----------------------------------------------------------------------
!
!   define nonstandard fields
!
!-----------------------------------------------------------------------
  
    if (ltavg_fmt_out_nc) then

      do nn = 1, num_avail_tavg_nstd_fields
        if ( (nn == tavg_MOC    .and. ltavg_moc_diags(ns)   ) .or.  &
             (nn == tavg_N_HEAT .and. ltavg_n_heat_trans(ns)) .or.  &
             (nn == tavg_N_SALT .and. ltavg_n_salt_trans(ns)))      then
          
          !*** get info on cell_methods
          if (nn == tavg_MOC) then
             method_integer = TAVG_BUF_3D_METHOD(tavg_loc_WVEL)
             call tavg_get_cell_method_string (method_integer,method_string)
          else if (nn == tavg_N_HEAT) then
             method_integer = TAVG_BUF_2D_METHOD(tavg_loc_ADVT)
             call tavg_get_cell_method_string (method_integer,method_string)
          else if (nn == tavg_N_SALT) then
             method_integer = TAVG_BUF_2D_METHOD(tavg_loc_ADVS)
             call tavg_get_cell_method_string (method_integer,method_string)
          endif
          call data_set_nstd_ccsm (                               &
              tavg_file_desc(ns), 'define', nstd_field_id,        &
              ndims_nstd_ccsm(nn), io_dims_nstd_ccsm(:,nn),       &
              short_name=avail_tavg_nstd_fields(nn)%short_name,   &
               long_name=avail_tavg_nstd_fields(nn)%long_name,    &
                   units=avail_tavg_nstd_fields(nn)%units,        &
                time_dim=time_dim,                                &
             coordinates=avail_tavg_nstd_fields(nn)%coordinates,  &
              fill_value=undefined_nf_r4,                         & ! kludge for r4 MOC,etc
              method_string=trim(method_string),                  &
                  nftype=avail_tavg_nstd_fields(nn)%nftype        )

          if (nn == tavg_MOC) then
              moc_id = nstd_field_id
          elseif (nn == tavg_N_HEAT) then
              n_heat_id = nstd_field_id
          elseif (nn == tavg_N_SALT) then
              n_salt_id = nstd_field_id
          endif
         endif ! streams check
        enddo !nn 
    endif !ltavg_fmt_out_nc

!-----------------------------------------------------------------------
!
!   write nonstandard fields
!
!-----------------------------------------------------------------------

    if (ltavg_fmt_out_nc .and. ltavg_write_reg) then
        call timer_start(timer_write_nstd)
        call tavg_write_vars_ccsm (tavg_file_desc(ns),num_ccsm_coordinates,   ccsm_coordinates(:,ns))
        call tavg_write_vars_ccsm (tavg_file_desc(ns),num_ccsm_time_invar(ns),ccsm_time_invar(:,ns))
        call tavg_write_vars_ccsm (tavg_file_desc(ns),num_ccsm_scalars(ns),   ccsm_scalars(:,ns))
        call timer_stop(timer_write_nstd)
    endif !ltavg_fmt_out_nc .and. ltavg_write_reg

   endif !.not. tavg_streams(ns)%ltavg_file_is_open
!  ================================================

!-----------------------------------------------------------------------
!
!   write fields to file
!   in this second phase, we actually write the data for all the fields.
!   after writing all fields, the field descriptors are destroyed and the
!   file can be closed
!
!-----------------------------------------------------------------------

   if (ltavg_write_reg) then
     time_to_close = check_time_flag(tavg_streams(ns)%file_flag)
   else
     time_to_close = .true.
   endif

   if (ltavg_fmt_out_nc .and. ltavg_write_reg) then

       call tavg_write_vars_ccsm (tavg_file_desc(ns),1,time_coordinate(1,ns))

       call tavg_write_time_bounds(tavg_file_desc(ns), ns)

       if (time_to_close) then
         call destroy_io_field(time_coordinate (1,ns))

         do n=1,num_ccsm_coordinates
          call destroy_io_field(ccsm_coordinates(n,ns))
         enddo

         do n=1,num_ccsm_time_invar(ns)
          call destroy_io_field(ccsm_time_invar (n,ns))
         enddo

         do n=1,num_ccsm_scalars(ns)
          call destroy_io_field(ccsm_scalars    (n,ns))
         enddo
       endif! time_to_close

!-----------------------------------------------------------------------
!
!  use data_set_nstd_ccsm to write labels and transport diags
!
!-----------------------------------------------------------------------

!==================> can this be moved up???
     call timer_start(timer_write_nstd)
     call tavg_write_vars_nstd_ccsm (tavg_file_desc(ns),ns) 
     call timer_stop(timer_write_nstd)

    endif !ltavg_fmt_out_nc .and. ltavg_write_reg


   !*** standard 2D and 3D fields
   call timer_start(timer_write_std)

   field_counter = 0

   do nfield = 1,num_avail_tavg_fields
      !*** write only if field is requested, in buffer, and in this stream
      if (tavg_requested(nfield)) then
        if (tavg_in_this_stream(nfield,ns)) then
         field_counter = field_counter + 1
         call data_set (tavg_file_desc(ns), 'write', tavg_streams(ns)%tavg_fields(field_counter))
         if (time_to_close) call destroy_io_field(tavg_streams(ns)%tavg_fields(field_counter))
        endif ! tavg_in_this_stream
      endif ! tavg_requested
   end do ! nfield

   call timer_stop(timer_write_std)

!-----------------------------------------------------------------------
!
!  after writing all fields, determine if the file can be closed
!
!-----------------------------------------------------------------------

   if (time_to_close) then

     tavg_streams(ns)%tavg_num_time_slices = 0

     deallocate(tavg_streams(ns)%tavg_fields)
     call data_set (tavg_file_desc(ns), 'close')

     if (my_task == master_task) then
       write(stdout,blank_fmt)
       write(stdout,*) 'tavg file written: ', trim(tavg_file_desc(ns)%full_name)
       call POP_IOUnitsFlush(POP_stdout);call POP_IOUnitsFlush(stdout)
     endif
   else
     call data_set (tavg_file_desc(ns), 'flush')
     if (my_task == master_task) then
       write(stdout,blank_fmt)
       write(stdout,*) 'data appended to tavg file: ', trim(tavg_file_desc(ns)%full_name)
       call POP_IOUnitsFlush(POP_stdout);call POP_IOUnitsFlush(stdout)
     endif
   endif ! time_to_close

!-----------------------------------------------------------------------
!
!  if pointer files are used, write tavg filenames to pointer file
!  do this only for restart tavg files, not regular tavg files
!
!-----------------------------------------------------------------------

   if (luse_pointer_files .and. .not. ltavg_write_reg) then
     write(ns_temp,'(i1)') ns
     call get_unit(nu)
     if (my_task == master_task) then
       if (ns == 1) then
         tavg_pointer_file = trim(pointer_filename)/&
                                                    &/'.tavg'
       else
         tavg_pointer_file = trim(pointer_filename)/&
                                                    &/'.tavg.'/&
                                                    &/trim(ns_temp)
       endif

       open(nu,file=tavg_pointer_file,form='formatted', status='unknown')
       write(nu,'(a)') trim(tavg_file_desc(ns)%full_name)
       close(nu)
     endif !master_task
     call release_unit(nu)
   endif !luse_pointer_files


!-----------------------------------------------------------------------
!
!  if this is a regular tavg write, reset time averages
!  if this is a restart tavg write, denormalize
!   in case a normal restart dump is written on the same timestep
!
!-----------------------------------------------------------------------

   if (ltavg_write_reg) then
     tavg_sum(ns) = c0
     call time_stamp('now', 'ymd',date_string=beg_date)
     if (tavg_freq_iopt(ns) == freq_opt_nstep) write(beg_date,'(i10)') nsteps_total
     call tavg_reset_field_stream(ns)
   else
     call tavg_norm_field_all ('denormalize',ns)
   endif !ltavg_write_reg

!-----------------------------------------------------------------------
!
!  get rid of file descriptor
!
!-----------------------------------------------------------------------

   if (time_to_close) call destroy_file(tavg_file_desc(ns))

   tavg_streams(ns)%ltavg_first_header = .false.

   if (tavg_streams(ns)%freq_iopt == freq_opt_once) then
     call override_time_flag(tavg_streams(ns)%field_flag, done=.true.)
   endif

   tavg_streams(ns)%ltavg_file_is_open = .not. time_to_close
 
  endif !ltavg_write_reg .or. ltavg_write_rest

  enddo ! nstreams

!-----------------------------------------------------------------------
!EOC

 end subroutine write_tavg


!***********************************************************************
!BOP
! !IROUTINE: read_tavg
! !INTERFACE:

 subroutine read_tavg (ns)

! !DESCRIPTION:
!  This routine reads a time average restart dump to continue
!  running time averages of requested fields.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) ::  &
      ns

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (i4) ::      &
     nu,                &! i/o unit
     iblock,            &! dummy block index
     n,                 &! dummy for indexing character string
     in_fields,         &! num of fields in restart file
     nfield,            &! dummy field counter
     hdr_error,         &! error file for reading restart hdr
     in_nsteps_total,   &! nsteps_total according to tavg file
     in_iyear,          &! iyear according to tavg file
     in_imonth,         &! imonth according to tavg file
     in_iday,           &! iday according to tavg file
     loc,               &! buffer location
     errVal              ! internal error flag

   real (r8) ::         &
     in_tday             ! tday according to tavg file

   logical (log_kind  ) ::  &
     file_exists

   character (char_len) ::  &
     header_filename,   &! filename for restart contents
     char_temp,         &! for string manipulation
     tavg_pointer_file, &! filename for pointer file containing
                         !   location/name of last restart file
     ns_temp

   type (io_field_desc), dimension(:), allocatable :: &
      tavg_fields_in       ! io field description for each field in file

   type (datafile) ::  &
      tavg_file_desc_in    ! IO file descriptor


!-----------------------------------------------------------------------
!
!  if pointer files are used, pointer file must be read to get 
!  actual filenames
!
!-----------------------------------------------------------------------

   errVal = 0

   call get_unit(nu)

   if (luse_pointer_files) then

      if (my_task == master_task) then
         tavg_pointer_file = char_blank
         write(ns_temp,'(i1)') ns
         if (ns > 1) then
            tavg_pointer_file = trim(pointer_filename)/&
                                                      &/'.tavg.'/&
                                                      &/trim(ns_temp)
         else
            tavg_pointer_file = trim(pointer_filename)/&
                                                      &/'.tavg'
         endif

         write(stdout,*) 'Reading pointer file: ', trim(tavg_pointer_file)
         call POP_IOUnitsFlush(POP_stdout); call POP_IOUnitsFlush(stdout)

         inquire(file=trim(tavg_pointer_file),exist=file_exists)
         if (file_exists) then
            open(nu, file=trim(tavg_pointer_file), form='formatted', &
                     status='old')
            read(nu,'(a)') tavg_infile
            close(nu)
         else
            errVal = -1
         endif

      endif
      call broadcast_scalar(errVal, master_task)

      if (errVal .ne. 0) then
        exit_string = 'FATAL ERROR: tavg_rpointer_file does not exist. pop2 model will exit '
        call document ('read_tavg', exit_string)
        call exit_POP (sigAbort,exit_string,out_unit=stdout)
      endif

      call broadcast_scalar(tavg_infile, master_task)

   endif

   call release_unit(nu)

!-----------------------------------------------------------------------
!
!  define input file
!
!-----------------------------------------------------------------------

   tavg_file_desc_in = construct_file (tavg_streams(ns)%fmt_in,       &
                                    full_name=trim(tavg_infile),   &
                                    record_length = rec_type_real, &
                                    recl_words=nx_global*ny_global)

!-----------------------------------------------------------------------
!
!  define scalar fields in file as file attributes to be read during
!  open
!
!-----------------------------------------------------------------------

   if (lccsm) then
   else
      call add_attrib_file(tavg_file_desc_in, 'tday'        , tday)
      call add_attrib_file(tavg_file_desc_in, 'iyear'       , iyear)
      call add_attrib_file(tavg_file_desc_in, 'imonth'      , imonth)
      call add_attrib_file(tavg_file_desc_in, 'iday'        , iday)
      call add_attrib_file(tavg_file_desc_in, 'beg_date'    , beg_date)
   endif
   call add_attrib_file(tavg_file_desc_in, 'nsteps_total', nsteps_total)
   call add_attrib_file(tavg_file_desc_in, 'tavg_sum'    , tavg_sum(ns))
   if (tavg_streams(ns)%ltavg_qflux_method_on) &
   call add_attrib_file(tavg_file_desc_in, 'tavg_sum_qflux'  , tavg_sum_qflux(ns))


!-----------------------------------------------------------------------
!
!  open input file
!  this will also extract scalar variables which are defined as
!  file attributes
!
!-----------------------------------------------------------------------

   call data_set (tavg_file_desc_in, 'open_read')

   if (lccsm) then
   else
      call extract_attrib_file(tavg_file_desc_in, 'beg_date', beg_date)
      !call extract_attrib_file(tavg_file_desc_in, 'tday', in_tday)
      !call extract_attrib_file(tavg_file_desc_in, 'iyear', in_iyear)
      !call extract_attrib_file(tavg_file_desc_in, 'imonth', in_imonth)
      !call extract_attrib_file(tavg_file_desc_in, 'iday', in_iday)
   endif
   call extract_attrib_file(tavg_file_desc_in, 'nsteps_total', &
                                          in_nsteps_total)
   call extract_attrib_file(tavg_file_desc_in, 'tavg_sum', tavg_sum(ns))
   if (tavg_streams(ns)%ltavg_qflux_method_on) &
   call extract_attrib_file(tavg_file_desc_in, 'tavg_sum_qflux', tavg_sum_qflux(ns))

   !*** report nsteps total and tavg_sum
   if (my_task == master_task) then
      write(stdout,'(i6,a29,i6,a35)') &
      in_nsteps_total,' nsteps_total in tavg restart', &
      nsteps_total,   ' nsteps_total in current simulation'
      write(stdout,*) ' tavg_sum = ', tavg_sum(ns), ' in tavg restart'
      call POP_IOUnitsFlush(POP_stdout); call POP_IOUnitsFlush(stdout)
   endif

   !*** check nsteps total for validity
   if (in_nsteps_total /= nsteps_total) then
      write(stdout,'(a,i6,a29,i6,a35)') 'ERROR: ', &
        in_nsteps_total,' nsteps_total in tavg restart', &
        nsteps_total,   ' nsteps_total in current simulation'
      exit_string = 'FATAL ERROR: TAVG restart file has wrong time step?'
      call document ('read_tavg', exit_string)
      call exit_POP (sigAbort,exit_string,out_unit=stdout)
   endif

!-----------------------------------------------------------------------
!
!  define requested fields to read in from file
!  NOTE: This requires that the tavg_contents file is consistent
!  with the tavg restart file.  There are currently no checks on this.
!
!-----------------------------------------------------------------------

   !*** define dimensions

   allocate(tavg_fields_in(num_avail_tavg_fields))

   do nfield = 1,num_avail_tavg_fields
      loc = avail_tavg_fields(nfield)%buf_loc
      if (loc > 0) then
        if (tavg_in_this_stream(nfield,ns)) then
         if (avail_tavg_fields(nfield)%ndims == 2) then

            tavg_fields_in(nfield) = construct_io_field(                &
                            avail_tavg_fields(nfield)%short_name,    &
                            dim1=i_dim, dim2=j_dim,                  &
                  long_name=avail_tavg_fields(nfield)%long_name,     &
                  units    =avail_tavg_fields(nfield)%units    ,     &
                  grid_loc =avail_tavg_fields(nfield)%grid_loc ,     &
                 field_loc =avail_tavg_fields(nfield)%field_loc,     &
                field_type =avail_tavg_fields(nfield)%field_type,    &
                valid_range=avail_tavg_fields(nfield)%valid_range,   &
#ifdef TAVG_R8
                 d2d_array =TAVG_BUF_2D(:,:,:,loc) ) !nonstandard, debugging only
#else
                 r2d_array =TAVG_BUF_2D(:,:,:,loc) )
#endif
         else if (avail_tavg_fields(nfield)%ndims == 3) then

            if (tavg_streams(ns)%ltavg_fmt_in_nc) then
              select case (trim(avail_tavg_fields(nfield)%grid_loc(4:4)))
                case('1')
                  z_dim = zt_dim
                case('2')
                  z_dim = zw_dim_top
                case('3')
                  z_dim = zw_dim_bot
                case('4')
                  z_dim = zt_150m_dim
              end select
            else
              z_dim = k_dim
            endif ! lccsm

            tavg_fields_in(nfield) = construct_io_field(                &
                            avail_tavg_fields(nfield)%short_name,    &
                            dim1=i_dim, dim2=j_dim, dim3=z_dim,      &
                  long_name=avail_tavg_fields(nfield)%long_name,     &
                  units    =avail_tavg_fields(nfield)%units    ,     &
                  grid_loc =avail_tavg_fields(nfield)%grid_loc ,     &
                 field_loc =avail_tavg_fields(nfield)%field_loc,     &
                field_type =avail_tavg_fields(nfield)%field_type,    &
                valid_range=avail_tavg_fields(nfield)%valid_range,   &
#ifdef TAVG_R8
                 d3d_array =TAVG_BUF_3D(:,:,:,:,loc) ) !nonstandard, debugging only
#else
                 r3d_array =TAVG_BUF_3D(:,:,:,:,loc) )
#endif
         endif

         call data_set (tavg_file_desc_in, 'define', tavg_fields_in(nfield))
        endif ! tavg_in_this_stream(nfield,ns) 
      endif
   end do

!-----------------------------------------------------------------------
!
!  now we actually read each field
!  after reading, get rid of io field descriptors and close file
!
!-----------------------------------------------------------------------

   do nfield = 1,num_avail_tavg_fields
      loc = avail_tavg_fields(nfield)%buf_loc
      if (loc > 0) then
        if (tavg_in_this_stream(nfield,ns)) then
         call data_set (tavg_file_desc_in, 'read', tavg_fields_in(nfield))
         call destroy_io_field(tavg_fields_in(nfield))
        endif !tavg_in_this_stream(nfield,ns)
      endif
   end do

   deallocate(tavg_fields_in)
   call data_set (tavg_file_desc_in, 'close')

   if (my_task == master_task) then
     write(stdout,blank_fmt)
     write(stdout,*) ' file read: ', tavg_infile
   endif

   call destroy_file(tavg_file_desc_in)
   call release_unit(nu)

!-----------------------------------------------------------------------
!
!  de-normalize sums
!
!-----------------------------------------------------------------------

   call tavg_norm_field_all ('denormalize',ns)

!-----------------------------------------------------------------------

   call tavg_global (ns)   ! print global sums of time averages

!-----------------------------------------------------------------------
!EOC

 end subroutine read_tavg

!***********************************************************************

! !IROUTINE: tavg_global
! !INTERFACE:

 subroutine tavg_global(ns)

! !DESCRIPTION:
!  Calculates and print global integrals of time average fields
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) ::  &
      ns

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (i4) ::     &
      k,               &   ! vertical level index
      ifield,          &   ! field identifier
      iblock,          &   ! block index
      nfield,          &   ! dummy field index
      field_loc,       &   ! field location (center,Nface,Eface,NEcorner)
      field_type           ! field type (scalar, vector, angle)

   real (r8) ::        &
      tavg_field_sum,  &   ! sum of tavg field
      tavg_norm            ! normalization for average

   real (r8), dimension (:,:,:), allocatable ::  &
      WORK               ! temp for holding area_weighted field

   real (r8), dimension (:,:), allocatable ::  &
      RMASK              ! topography mask for global sum

   character (char_len) ::  &
      time_string,          &
      date_string

!-----------------------------------------------------------------------
!
!  calculate globally-integrated time average of each chosen 2d field
!
!-----------------------------------------------------------------------

   allocate (RMASK(nx_block,ny_block), &
             WORK (nx_block,ny_block,nblocks_clinic))

   call time_stamp ('now','mdy',date_string=date_string,time_string=time_string)

   if (my_task == master_task) then
     write (stdout,blank_fmt)
     write (stdout,*) 'Global Time Averages: ' // trim(date_string) // ' ' // trim(time_string)
   endif

   fields_loop: do nfield=1,num_avail_tavg_fields
      ifield = avail_tavg_fields(nfield)%buf_loc

      if (ifield > 0 ) then
        if ( tavg_in_this_stream(nfield,ns)) then

         field_loc  = avail_tavg_fields(nfield)%field_loc
         field_type = avail_tavg_fields(nfield)%field_type

         if (avail_tavg_fields(nfield)%method == tavg_method_avg) then
            tavg_norm = tavg_sum(ns)
         else if (avail_tavg_fields(nfield)%method == tavg_method_qflux) then
            tavg_norm = tavg_sum_qflux(ns)
         else
            tavg_norm = c1
         endif

         if (tavg_norm == c0) then
            if (my_task == master_task) then
              write(stdout,*) 'Cannot compute global integral of ', &
                               trim (avail_tavg_fields(nfield)%short_name)
              call POP_IOUnitsFlush(POP_stdout);call POP_IOUnitsFlush(stdout)
            endif
            cycle fields_loop 
            !*** call exit_POP (SigAbort,'(tavg_gloal) ERROR: tavg_norm = 0 in tavg_global',out_unit=stdout)
         endif

         !*** 2-d fields

         if (avail_tavg_fields(nfield)%ndims == 2) then

!njn01      !$OMP PARALLEL DO PRIVATE(iblock)
            do iblock = 1,nblocks_clinic
               select case(field_loc)
               case(field_loc_center)
                  WORK(:,:,iblock)  = TAVG_BUF_2D(:,:,iblock,ifield)* &
                                    TAREA(:,:,iblock)*RCALCT(:,:,iblock)
               case(field_loc_NEcorner)
                  WORK(:,:,iblock)  = TAVG_BUF_2D(:,:,iblock,ifield)* &
                                    UAREA(:,:,iblock)*RCALCU(:,:,iblock)
               case default ! make U cell the default for all other cases
                  WORK(:,:,iblock)  = TAVG_BUF_2D(:,:,iblock,ifield)* &
                                    UAREA(:,:,iblock)*RCALCU(:,:,iblock)
               end select
            end do
!njn01      !$OMP END PARALLEL DO

            tavg_field_sum = global_sum(WORK, distrb_clinic, field_loc)

            select case(field_loc)
            case(field_loc_center)
               tavg_field_sum = tavg_field_sum/(tavg_norm*area_t)
            case(field_loc_NEcorner)
               tavg_field_sum = tavg_field_sum/(tavg_norm*area_u)
            case default ! make U cell the default for all other cases
               tavg_field_sum = tavg_field_sum/(tavg_norm*area_u)
            end select

         !*** 3-d fields

         else
      
!njn01      !$OMP PARALLEL DO PRIVATE(iblock,k,RMASK)
            do iblock = 1,nblocks_clinic
               WORK(:,:,iblock) = c0

               select case(field_loc)

               case(field_loc_center)
                  do k=1,km
                     RMASK(:,:) = merge(c1, c0, k <= KMT(:,:,iblock)) 
                     WORK(:,:,iblock) = WORK(:,:,iblock) + dz(k)* &
                                        TAVG_BUF_3D(:,:,k,iblock,ifield)* &
                                        TAREA(:,:,iblock)*RMASK
                  end do

               case(field_loc_NEcorner)
                  do k=1,km
                     RMASK(:,:) = merge(c1, c0, k <= KMU(:,:,iblock)) 
                     WORK(:,:,iblock) = WORK(:,:,iblock) + dz(k)* &
                                        TAVG_BUF_3D(:,:,k,iblock,ifield)* &
                                        UAREA(:,:,iblock)*RMASK
                  end do

               case default ! make U cell the default for all other cases
                  do k=1,km
                     RMASK(:,:) = merge(c1, c0, k <= KMU(:,:,iblock)) 
                     WORK(:,:,iblock) = WORK(:,:,iblock) + dz(k)* &
                                        TAVG_BUF_3D(:,:,k,iblock,ifield)* &
                                        UAREA(:,:,iblock)*RMASK
                  end do

               end select
            end do
!njn01      !$OMP END PARALLEL DO

            tavg_field_sum = global_sum(WORK, distrb_clinic, field_loc)

            select case(field_loc)
            case(field_loc_center)
               tavg_field_sum = tavg_field_sum/(tavg_norm*volume_t)
            case(field_loc_NEcorner)
               tavg_field_sum = tavg_field_sum/(tavg_norm*volume_u)
            case default ! make U cell the default for all other cases
               tavg_field_sum = tavg_field_sum/(tavg_norm*volume_u)
            end select

         endif

         if (my_task == master_task) then
            write (stdout,*) trim(avail_tavg_fields(nfield)%short_name), &
                             ': ', tavg_field_sum
         endif
       endif ! tavg_in_this_stream
      endif ! ifield > 0
   end do fields_loop

   deallocate (RMASK, WORK)

!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_global
!***********************************************************************
!BOP
! !IROUTINE: tavg_global_sum_2D
! !INTERFACE:

 function tavg_global_sum_2D (id)

! !DESCRIPTION:
!  Calculates the global sum of a requested 2D time-averaged field
!  Presently, this is a special-purpose routine used only by the
!  budget_diagnostics module.  It could be extended and generalized,
!  and used with a revised version of subroutine tavg_global, if time
!  and interest allow.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) ::  &
      id                               ! identifier of time-averaged field 

! !OUTPUT PARAMETERS:

   real (r8) :: &
      tavg_global_sum_2D     ! result of this function: the global sum of
                             !   a requested 2D time-averaged field


!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (i4) ::     &
      loc,             &! field buffer location
      ns,              &! stream index
      nstream,         &! stream identifier
      iblock,          &! block index
      nfield,          &! dummy field index
      field_loc,       &! field location (center,Nface,Eface,NEcorner)
      field_type        ! field type (scalar, vector, angle)

   real (r8) ::        &
      tavg_norm         ! normalization for average

   real (r8), dimension (nx_block,ny_block,max_blocks_clinic) ::  &
      WORK              ! temp for holding area_weighted field

   real (r8), dimension (nx_block, ny_block) ::  &
      RMASK             ! topography mask for global sum


!-----------------------------------------------------------------------
!
!  does this field exist?
!
!-----------------------------------------------------------------------

   if (.not. tavg_requested(id)) then
      exit_string = 'FATAL ERROR: invalid field request'
      call document ('tavg_global_sum_2D', 'id  ', id)
      call document ('tavg_global_sum_2D', 'tavg_requested(id) ', tavg_requested(id))
      call document ('tavg_global_sum_2D', exit_string)
      call exit_POP (sigAbort,exit_string,out_unit=stdout)
   endif

!-----------------------------------------------------------------------
!
!  is this a 2D field?
!
!-----------------------------------------------------------------------

   if (avail_tavg_fields(id)%ndims /= 2) then
      exit_string = 'FATAL ERROR: invalid dimension'
      call document ('tavg_global_sum_2D', exit_string)
      call exit_POP (sigAbort,exit_string,out_unit=stdout)
   endif

!-----------------------------------------------------------------------
!
!  identify the requested field
!
!-----------------------------------------------------------------------

   field_loc  = avail_tavg_fields(id)%field_loc
   field_type = avail_tavg_fields(id)%field_type

!-----------------------------------------------------------------------
!
!  identify the buffer location in TAVG_BUF_2D
!
!-----------------------------------------------------------------------

   loc = avail_tavg_fields(id)%buf_loc
   if (loc <= 0) then
      write(exit_string,*)  'FATAL ERROR -- invalid loc; loc = ', loc
      call document ('tavg_global_sum_2D', exit_string)
      call exit_POP (sigAbort,exit_string,out_unit=stdout)
   endif

!-----------------------------------------------------------------------
!
!  locate the stream in which the field belongs
!
!-----------------------------------------------------------------------

    nstream = tavg_in_which_stream(id)

!-----------------------------------------------------------------------
!
!  compute the global average of the 2D field, using tavg_sum
!   from the fields' stream
!
!-----------------------------------------------------------------------

   if (avail_tavg_fields(id)%method == tavg_method_avg) then
      tavg_norm = tavg_sum(nstream)
   else if (avail_tavg_fields(id)%method == tavg_method_qflux) then
      tavg_norm = tavg_sum_qflux(nstream)
   else
      tavg_norm = c1
   endif

   if (tavg_norm == c0) then
     exit_string = 'FATAL ERROR: tavg_norm = 0'
     call document ('tavg_global_sum_2D', exit_string)
     call exit_POP (sigAbort,exit_string,out_unit=stdout)
   endif

!njn01 !$OMP PARALLEL DO PRIVATE(iblock)
   do iblock = 1,nblocks_clinic
      select case(field_loc)
        case(field_loc_center)
           WORK(:,:,iblock)  = TAVG_BUF_2D(:,:,iblock,loc)* &
                               TAREA(:,:,iblock)*RCALCT(:,:,iblock)
        case(field_loc_NEcorner)
           WORK(:,:,iblock)  = TAVG_BUF_2D(:,:,iblock,loc)* &
                               UAREA(:,:,iblock)*RCALCU(:,:,iblock)
        case default ! make U cell the default for all other cases
           WORK(:,:,iblock)  = TAVG_BUF_2D(:,:,iblock,loc)* &
                               UAREA(:,:,iblock)*RCALCU(:,:,iblock)
      end select
   end do
!njn01 !$OMP END PARALLEL DO

   tavg_global_sum_2D = global_sum(WORK, distrb_clinic, field_loc)

   select case(field_loc)
     case(field_loc_center)
        tavg_global_sum_2D = tavg_global_sum_2D/(tavg_norm*area_t)
     case(field_loc_NEcorner)
        tavg_global_sum_2D = tavg_global_sum_2D/(tavg_norm*area_u)
     case default ! make U cell the default for all other cases
        tavg_global_sum_2D = tavg_global_sum_2D/(tavg_norm*area_u)
   end select



!-----------------------------------------------------------------------
!EOC

 end function tavg_global_sum_2D

!***********************************************************************
!BOP
! !IROUTINE: tavg_increment_sum_qflux
! !INTERFACE:

 subroutine tavg_increment_sum_qflux(const)

! !DESCRIPTION:
!  Increment the scalar tavg_sum_qflux with the weight for this timestep.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), intent(in) ::  &
      const
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::  &
      n

   do n=1,nstreams
     tavg_sum_qflux(n) = tavg_sum_qflux(n) + const ! const = tlast_ice
   enddo

!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_increment_sum_qflux

!***********************************************************************
!BOP
! !IROUTINE: accumulate_tavg_field
! !INTERFACE:

 subroutine accumulate_tavg_field(ARRAY,field_id,block,k,const)

! !DESCRIPTION:
!  This routine updates a tavg field.  If the time average of the
!  field is requested, it accumulates a time sum of a field by 
!  multiplying by the time step and accumulating the sum into the 
!  tavg buffer array.  If the min or max of a field is requested, it
!  checks the current value and replaces the min, max if the current
!  value is less than or greater than the stored value.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      block,           &! local block address (in baroclinic distribution)
      k,               &! vertical level
      field_id          ! index into available fields for tavg field info

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      ARRAY             ! array of data for this block to add to 
                        !  accumulated sum in tavg buffer
   real (r8), optional, intent(in) ::  &
      const
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      bufloc,            &! location of field in tavg buffer
      ndims               ! rank of field (2=2d,3=3d)

!-----------------------------------------------------------------------
! 
!  test: mix_pass, ltavg_on, and tavg_requested.                
!                                                              
!-----------------------------------------------------------------------
                                                                
   if (.not. accumulate_tavg_now(field_id)) return             
                                                              

!-----------------------------------------------------------------------
!
!  get buffer location and field info from avail_tavg_field array
!
!-----------------------------------------------------------------------

   bufloc = avail_tavg_fields(field_id)%buf_loc

   if (bufloc <= 0) then
     exit_string = 'FATAL ERROR: attempt to accumulate bad tavg field'
     call document ('accumulate_tavg_field', exit_string)
     call exit_POP (sigAbort,exit_string,out_unit=stdout)
   endif

   ndims = avail_tavg_fields(field_id)%ndims

   if ((ndims == 3) .and. &
       (avail_tavg_fields(field_id)%grid_loc(4:4) == '4') .and. &
       (k > zt_150m_levs)) return

!-----------------------------------------------------------------------
!
!  update the field into the tavg buffer
!
!-----------------------------------------------------------------------

   select case (avail_tavg_fields(field_id)%method)

   case (tavg_method_avg)  ! accumulate running time sum for time avg
      if (ndims == 2) then
         TAVG_BUF_2D(:,:,block,bufloc) = &
         TAVG_BUF_2D(:,:,block,bufloc) + dtavg*ARRAY
      else
         TAVG_BUF_3D(:,:,k,block,bufloc) = &
         TAVG_BUF_3D(:,:,k,block,bufloc) + dtavg*ARRAY
      endif
   case (tavg_method_qflux)  
         TAVG_BUF_2D(:,:,block,bufloc) =  &
         TAVG_BUF_2D(:,:,block,bufloc) + const*max (c0,ARRAY)
   case (tavg_method_min)  ! replace with current minimum value
      if (ndims == 2) then
         where (ARRAY < TAVG_BUF_2D(:,:,block,bufloc))
            TAVG_BUF_2D(:,:,block,bufloc) = ARRAY
         end where
      else
         where (ARRAY < TAVG_BUF_3D(:,:,k,block,bufloc))
            TAVG_BUF_3D(:,:,k,block,bufloc) = ARRAY
         end where
      endif
   case (tavg_method_max)  ! replace with current minimum value
      if (ndims == 2) then
         where (ARRAY > TAVG_BUF_2D(:,:,block,bufloc))
            TAVG_BUF_2D(:,:,block,bufloc) = ARRAY
         end where
      else
         where (ARRAY > TAVG_BUF_3D(:,:,k,block,bufloc))
            TAVG_BUF_3D(:,:,k,block,bufloc) = ARRAY
         end where
      endif
   case (tavg_method_constant)  ! overwrite with current value; intended for time-invariant fields
      if (ndims == 2) then
         TAVG_BUF_2D(:,:,block,bufloc) = ARRAY
      else
         TAVG_BUF_3D(:,:,k,block,bufloc) = ARRAY
      endif
   case default
   end select

!-----------------------------------------------------------------------
!EOC

 end subroutine accumulate_tavg_field

!***********************************************************************
!BOP
! !IROUTINE: accumulate_tavg_now
! !INTERFACE:

 function accumulate_tavg_now(field_id)

! !DESCRIPTION:
!  This function determines whether an available (defined) tavg field
!  can be accumulated at this time.
!
!  If the following are true:
!    1) mix_pass .ne. 1
!    2) the field is requested
!    3) tavg is on for the stream which contains the field
!  then accumulate_tavg_now is true 
!
! !REVISION HISTORY:
!  Nancy J. Norton 7 April 2010

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      field_id            
                         

! !OUTPUT PARAMETERS:

   logical (log_kind) :: &
      accumulate_tavg_now 
                         

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      tavg_stream_num            


   accumulate_tavg_now   = .false.

   if (mix_pass /= 1) then
         
     if (tavg_requested(field_id)) then
       tavg_stream_num  = tavg_in_which_stream(field_id)
       accumulate_tavg_now = ltavg_on(tavg_stream_num)
     endif

   endif


!-----------------------------------------------------------------------
!EOC

 end function accumulate_tavg_now

!***********************************************************************
!BOP
! !IROUTINE: tavg_reset_field_all
! !INTERFACE:

 subroutine tavg_reset_field_all
 
! !DESCRIPTION:
!  Resets all TAVG_BUF_2D and TAVG_BUF_3D arrays 
!
! !REVISION HISTORY:
!  same as module
 
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
 
   integer (int_kind) ::  &
      iblock,             &
      nfield


 !$OMP PARALLEL DO PRIVATE(iblock,nfield)
  do iblock=1,nblocks_clinic
     do nfield=1,tavg_bufsize_2d
       select case (TAVG_BUF_2D_METHOD(nfield))
        case (tavg_method_avg)
           TAVG_BUF_2D(:,:,  iblock,nfield) = c0
        case (tavg_method_qflux)
           TAVG_BUF_2D(:,:,  iblock,nfield) = c0
           tavg_sum_qflux = c0
        case (tavg_method_min)
           TAVG_BUF_2D(:,:,  iblock,nfield) = bignum
        case (tavg_method_max)
           TAVG_BUF_2D(:,:,  iblock,nfield) = -bignum
        case default
           TAVG_BUF_2D(:,:,  iblock,nfield) = c0
        end select
     end do
 
     do nfield=1,tavg_bufsize_3d
        select case (TAVG_BUF_3D_METHOD(nfield))
        case (tavg_method_avg)
           TAVG_BUF_3D(:,:,:,iblock,nfield) = c0
        case (tavg_method_min)
           TAVG_BUF_3D(:,:,:,iblock,nfield) = bignum
        case (tavg_method_max)
           TAVG_BUF_3D(:,:,:,iblock,nfield) = -bignum
        case default
           TAVG_BUF_3D(:,:,:,iblock,nfield) = c0
        end select
     end do
  end do
 !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_reset_field_all


!***********************************************************************
!BOP
! !IROUTINE: tavg_reset_field_stream
! !INTERFACE:

 subroutine tavg_reset_field_stream(ns)

! !DESCRIPTION:
!  Resets all TAVG_BUF_2D and TAVG_BUF_3D arrays
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) ::  &
      ns

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::  &
      iblock,             &
      nfield,             &
      loc


   do nfield = 1,num_avail_tavg_fields  ! check all available fields
     loc = avail_tavg_fields(nfield)%buf_loc ! locate field in buffer

     !*** continue only if field is in buffer and in this stream
     if (loc > 0) then
     if (tavg_in_this_stream(nfield,ns)) then

       if (avail_tavg_fields(nfield)%ndims == 2) then

         select case (TAVG_BUF_2D_METHOD(loc))
            case (tavg_method_avg)
               TAVG_BUF_2D(:,:,:,loc) = c0
            case (tavg_method_qflux)
               TAVG_BUF_2D(:,:,:,loc) = c0
              tavg_sum_qflux(ns) = c0
            case (tavg_method_min)
               TAVG_BUF_2D(:,:,:,loc) = bignum
            case (tavg_method_max)
               TAVG_BUF_2D(:,:,:,loc) = -bignum
            case default
               TAVG_BUF_2D(:,:,:,loc) = c0
            end select

       elseif (avail_tavg_fields(nfield)%ndims == 3) then

         select case (TAVG_BUF_3D_METHOD(loc))
         case (tavg_method_avg)
            TAVG_BUF_3D(:,:,:,:,loc) = c0
         case (tavg_method_min)
            TAVG_BUF_3D(:,:,:,:,loc) = bignum
         case (tavg_method_max)
            TAVG_BUF_3D(:,:,:,:,loc) = -bignum
         case default
            TAVG_BUF_3D(:,:,:,:,loc) = c0
         end select

       endif ! 2D/3D
     endif ! in stream
     endif ! loc
   enddo ! nfield

!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_reset_field_stream


!***********************************************************************
!BOP
! !IROUTINE: tavg_norm_field_all
! !INTERFACE:

 subroutine tavg_norm_field_all (norm_flag,ns)
 
! !DESCRIPTION:
!  Normalizes or de-normalizes all TAVG_BUF_2D and TAVG_BUF_3D arrays 
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:
 
   character (*), intent(in) ::  &
      norm_flag

   integer (int_kind), intent(in) ::  &
      ns

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
 
   real (r8) ::  &
      factor,    &
      factorq

   integer (int_kind) ::  &
      iblock,             &
      nfield,             &
      loc

   
   select case (trim(norm_flag))
 
      case ('normalize')
        if (tavg_sum(ns) /= 0) then
          factor  = c1/tavg_sum(ns)
        else
          exit_string = 'FATAL ERROR: attempt to divide by zero tavg_sum'
          call document ('tavg_norm_field_all', exit_string)
          call exit_POP (sigAbort,exit_string,out_unit=stdout)
        endif
        if (tavg_streams(ns)%ltavg_qflux_method_on) then
          if (tavg_sum_qflux(ns) /= 0) then
            factorq = c1/tavg_sum_qflux(ns)
          elseif (tavg_freq_iopt(ns) == freq_opt_nhour   .or.  &
                  tavg_freq_iopt(ns) == freq_opt_nsecond .or.  &
                  tavg_freq_iopt(ns) == freq_opt_nstep         ) then
            ! do nothing; these frequencies are not compatable with time-averaged qflux
            factorq = c1
          else
            call document ('tavg_norm_field_all', ' ns ', ns)
            exit_string = 'FATAL ERROR: attempt to divide by zero tavg_sum_qflux(ns)'
            call document ('tavg_norm_field_all', exit_string)
            call exit_POP (sigAbort,exit_string,out_unit=stdout)
          endif
        endif

      case ('denormalize')
          factor  = tavg_sum(ns)
          factorq = tavg_sum_qflux(ns)
      case default
          exit_string = 'FATAL ERROR: unknown option'
          call document ('tavg_norm_field_all', exit_string)
          call exit_POP (sigAbort,exit_string,out_unit=stdout)
   end select
 
 
  !$OMP PARALLEL DO PRIVATE(iblock,nfield,loc)
   do nfield = 1,num_avail_tavg_fields  ! check all available fields
      loc = avail_tavg_fields(nfield)%buf_loc ! locate field in buffer

      !*** continue only if field is requested, in buffer, and in this stream
      if (loc > 0 ) then
        if (tavg_in_this_stream(nfield,ns)) then

        if (avail_tavg_fields(nfield)%ndims == 2) then

           select case (TAVG_BUF_2D_METHOD(loc))
             case (tavg_method_avg)
               TAVG_BUF_2D(:,:,  :,loc) = TAVG_BUF_2D(:,:,  :,loc)*factor
             case (tavg_method_qflux)
               if (tavg_streams(ns)%ltavg_qflux_method_on) then
                 TAVG_BUF_2D(:,:,  :,loc) = TAVG_BUF_2D(:,:,  :,loc)*factorq
               endif
           end select

        else if (avail_tavg_fields(nfield)%ndims == 3) then

           if (TAVG_BUF_3D_METHOD(loc) == tavg_method_avg) then
              TAVG_BUF_3D(:,:,:,:,loc) = TAVG_BUF_3D(:,:,:,:,loc)*factor
           endif

        endif ! avail_tavg_fields
      endif ! stream
      endif ! loc
   end do

  !$OMP END PARALLEL DO
 

!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_norm_field_all


!***********************************************************************
!BOP
! !IROUTINE: define_tavg_field
! !INTERFACE:

 subroutine define_tavg_field(id, short_name, ndims, tavg_method,       &
                                  long_name, units,                     &
                                  grid_loc, valid_range,                &
                                  field_loc, field_type,coordinates,    &
                                  scale_factor,                         &
                                  nftype,                               &
                                  nstd_fields,                          &
                                  num_nstd_fields, max_nstd_fields      )

! !DESCRIPTION:
!  Initializes description of an available field and returns location
!  in the available fields array for use in later tavg calls.
!
! !REVISION HISTORY:
!  same as module

! !OUTPUT PARAMETERS:

   integer (int_kind), intent(out) :: &
      id                ! location in avail_fields array for use in
                        ! later tavg routines 

! !INPUT PARAMETERS:

   character(*), intent(in) :: &
      short_name                ! short name for field

   integer (i4), intent(in) :: &
      ndims                     ! number of dims of the field

   integer (i4), intent(in), optional :: &
      field_loc,              &! location in grid 
      field_type,             &! type of field (scalar, vector, angle)
      tavg_method              ! id for method of averaging
                               ! default is tavg_method_avg

   character(*), intent(in), optional :: &
      long_name,              &! long descriptive name for field
      units,                  &! physical units for field
      coordinates              ! CF coordinates

   character(4), intent(in), optional :: &
      grid_loc                 ! location in grid (in 4-digit code)

   real (rtavg), intent(in), optional :: &
      scale_factor             ! scale factor

   real (r4), dimension(2), intent(in), optional :: &
      valid_range              ! min/max

   character(*), intent(in), optional :: &
      nftype                    ! type string 

   integer (int_kind), intent(in), optional ::  &
      max_nstd_fields

! !INPUT/OUTPUT PARAMETERS:

   type (tavg_field_desc_ccsm), dimension(:), intent(inout), optional ::  &
      nstd_fields

   integer (int_kind), intent(inout), optional ::  &
      num_nstd_fields

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   type (tavg_field_desc_ccsm) ::  &
      tavg_field

   integer (int_kind) ::  &
      num_fields

   logical (log_kind) ::  &
      error,              &
      nonstandard_fields

!-----------------------------------------------------------------------
!
!  increment the number of defined fields and make sure it does not
!  exceed the maximum
!  return the id as the current number
!
!-----------------------------------------------------------------------

   error = .false.
   nonstandard_fields = .false.

   if (present (nstd_fields) ) then
      nonstandard_fields = .true.
      num_nstd_fields = num_nstd_fields + 1
      num_fields = num_nstd_fields
      if (num_nstd_fields > max_nstd_fields) error = .true.
   else
      num_avail_tavg_fields = num_avail_tavg_fields + 1
      num_fields = num_avail_tavg_fields
      if (num_avail_tavg_fields > max_avail_tavg_fields) error = .true.
   endif

   if (error) then
      exit_string = 'FATAL ERROR: defined tavg fields > max allowed'
      call document ('define_tavg_field', exit_string)
      call exit_POP (sigAbort,exit_string,out_unit=stdout)
   endif
 
   id = num_fields

   if (my_task == master_task .and. nsteps_run <= 1) then
     write(stdout,*) 'define_tavg_field: id = ', id, ' ', short_name
     call POP_IOUnitsFlush(POP_stdout);call POP_IOUnitsFlush(stdout)
   endif

!-----------------------------------------------------------------------
!
!  now fill the local field descriptor
!
!-----------------------------------------------------------------------

   tavg_field%ndims      = ndims
   tavg_field%short_name = short_name
   tavg_field%buf_loc    = 0  ! will be reset later

   if (present(long_name)) then
      tavg_field%long_name = long_name
   else
      tavg_field%long_name = char_blank
   endif

   if (present(units)) then
      tavg_field%units = units
   else
      tavg_field%units = char_blank
   endif

   if (present(grid_loc)) then
      tavg_field%grid_loc = grid_loc
   else
      tavg_field%grid_loc = '    '
   endif

   if (present(tavg_method)) then
      tavg_field%method = tavg_method
   else
      tavg_field%method = tavg_method_avg
   endif


!  kludge for  MOC and transport diagnostics -- always r4
   if (.not. nonstandard_fields) then
   if (present(scale_factor)) then
      if (scale_factor /= 1.0_rtavg) then
        tavg_field%scale_factor = scale_factor
        if (scale_factor /= 0.0_rtavg) then
          tavg_field%fill_value    = undefined_nf/scale_factor
        endif
      else
        tavg_field%scale_factor  = undefined_nf
        tavg_field%fill_value    = undefined_nf
      endif
   else
      tavg_field%scale_factor  = undefined_nf
      tavg_field%fill_value    = undefined_nf
   endif
   endif

   if (present(valid_range)) then
      tavg_field%valid_range = valid_range
   else
      tavg_field%valid_range = undefined
   endif

   if (present(coordinates)) then
      tavg_field%coordinates = coordinates
   else
      tavg_field%coordinates = char_blank
   endif

   !*** set field location, field type used by i/o, ghost cell update
   !*** and other communication routines.  because ghost cells for tavg
   !*** fields are not typically used, the default is field_xxx_noupdate

   if (present(field_loc)) then
      tavg_field%field_loc = field_loc
   else
      tavg_field%field_loc = field_loc_noupdate
      if (present(grid_loc)) then
         !*** try to decode field location from grid_loc
         if (grid_loc(2:2) == '1' .and. grid_loc(3:3) == '1') then
            tavg_field%field_loc = field_loc_center
         else if (grid_loc(2:2) == '2' .and. grid_loc(3:3) == '2') then
            tavg_field%field_loc = field_loc_NEcorner
         else if (grid_loc(2:2) == '1' .and. grid_loc(3:3) == '2') then
            tavg_field%field_loc = field_loc_Nface
         else if (grid_loc(2:2) == '2' .and. grid_loc(3:3) == '1') then
            tavg_field%field_loc = field_loc_Eface
         endif
      endif
   endif

   if (present(field_type)) then
      tavg_field%field_type = field_type
   else
      tavg_field%field_type = field_type_noupdate
   endif

   if (present(nftype)) then
      tavg_field%nftype = trim(nftype)
   else
      tavg_field%nftype = 'float'
   endif


   if (nonstandard_fields) then
     nstd_fields(id) = tavg_field
   else
     avail_tavg_fields(id) = tavg_field
   endif

   if (my_task == master_task .and. tavg_debug > 0) then
     call document ('define_tavg_field',  trim(tavg_field%short_name))
     call document ('define_tavg_field', 'buffer id number ', id)
     call document ('define_tavg_field',  trim(tavg_field%long_name))
     call document ('define_tavg_field',  trim(tavg_field%units))
     call document ('define_tavg_field',  trim(tavg_field%grid_loc))
     call document ('define_tavg_field', '_FillValue',   tavg_field%fill_value)
     call document ('define_tavg_field', 'scale_factor', tavg_field%scale_factor)
     call document ('define_tavg_field',  trim(tavg_field%coordinates))
     call document ('define_tavg_field', 'field_type',   tavg_field%field_type)
     call document ('define_tavg_field',  trim(tavg_field%nftype))
     call document ('define_tavg_field', 'end subroutine ')
   endif
!-----------------------------------------------------------------------
!EOC

 end subroutine define_tavg_field

!***********************************************************************
!BOP
! !IROUTINE: request_tavg_field
! !INTERFACE:

 subroutine request_tavg_field(short_name,ns)

! !DESCRIPTION:
!  This field marks an available field as requested and computes
!  the location in the tavg buffer array.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      short_name                ! the short name of the field

   integer (int_kind), intent(in) :: &
      ns                              ! stream in which field is located
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind)  ::  &
      n,                   &! loop index
      id                    ! location of field in avail_fields array

!-----------------------------------------------------------------------
!
!  search for field with same name
!
!-----------------------------------------------------------------------

   id = 0

   srch_loop: do n=1,num_avail_tavg_fields
      if (trim(avail_tavg_fields(n)%short_name) == short_name) then
         id = n
         exit srch_loop
      endif
   end do srch_loop

   if (id == 0) then
      write(stdout,*) 'Requested ', trim(short_name)
      exit_string = 'FATAL ERROR: requested field unknown'
      call document ('request_tavg_field', exit_string)
      call exit_POP (sigAbort,exit_string,out_unit=stdout)
   endif

!-----------------------------------------------------------------------
!
!  set the position in the buffer and advance the buffer position
!  for the next field
!
!-----------------------------------------------------------------------

   if (avail_tavg_fields(id)%ndims == 2) then
      tavg_bufsize_2d = tavg_bufsize_2d + 1
      avail_tavg_fields(id)%buf_loc = tavg_bufsize_2d
   else if (avail_tavg_fields(id)%ndims == 3) then
      tavg_bufsize_3d = tavg_bufsize_3d + 1
      avail_tavg_fields(id)%buf_loc = tavg_bufsize_3d
   endif

!-----------------------------------------------------------------------
!
!  if tavg is on, but not started yet, set the buf_loc to -buf_loc
!  to show that it is requested, but will not return true for
!  requested_tavg_field
!
!-----------------------------------------------------------------------

   if (.not. tavg_streams(ns)%ltavg_on) then
      avail_tavg_fields(id)%buf_loc = -avail_tavg_fields(id)%buf_loc
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine request_tavg_field

!***********************************************************************
!BOP
! !IROUTINE: tavg_requested
! !INTERFACE:

 function tavg_requested(id)

! !DESCRIPTION:
!  This function determines whether an available (defined) tavg field
!  has been requested by a user (through the input contents file) and 
!  returns true if it has.  Note that if tavg has been turned off, 
!  the function will always return false.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      id                   ! id returned by the define function which
                           !   gives the location of the field

! !OUTPUT PARAMETERS:

   logical (log_kind) :: &
      tavg_requested     ! result of checking whether the field has
                         !   been requested

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  check the buffer location - if zero, the field has not been
!  requested
!
!-----------------------------------------------------------------------

   if (id < 1 .or. id > num_avail_tavg_fields) then
      write(stdout,*) '(tavg_requested) id = ', id; call POP_IOUnitsFlush(stdout)
      exit_string = 'FATAL ERROR: invalid tavg id'
      call document ('tavg_requested', exit_string)
      call exit_POP (sigAbort,exit_string,out_unit=stdout)
   endif

   if (avail_tavg_fields(id)%buf_loc > 0) then
      tavg_requested = .true.
   else
      tavg_requested = .false.
   endif

!-----------------------------------------------------------------------
!EOC

 end function tavg_requested

!***********************************************************************
!BOP
! !IROUTINE: tavg_in_this_stream(id,stream_number)
! !INTERFACE:

 function tavg_in_this_stream(id,stream_number)

! !DESCRIPTION:
!  This function determines whether an available (defined) tavg field
!  has been selected in a given tavg_contents "stream."
!  returns true if it has.  
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      stream_number,     &! stream index
      id

! !OUTPUT PARAMETERS:

   logical (log_kind) :: &
      tavg_in_this_stream     ! result of checking whether the field has
                              !   been requested in this stream

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   integer (int_kind) ::  &
      ntest,              &
      buf_loc

!-----------------------------------------------------------------------
!
!  check to see if this field is in this stream
!
!-----------------------------------------------------------------------

  if(avail_tavg_fields(id)%stream_number == stream_number) then
    tavg_in_this_stream = .true.
  else
    tavg_in_this_stream = .false. 
  endif


!-----------------------------------------------------------------------
!EOC

 end function tavg_in_this_stream

!***********************************************************************
!BOP
! !IROUTINE: tavg_in_which_stream(id)
! !INTERFACE:

 function tavg_in_which_stream(id)

! !DESCRIPTION:
!  This subroutine locates the stream in which a given field 
!  resides. 
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      id   ! field buffer location ("buf_id")

! !OUTPUT PARAMETERS:

   integer (int_kind) ::  &
      tavg_in_which_stream ! stream in which this field resides

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   integer (int_kind) ::  &
      ns,                 &
      stream_number

   stream_number = avail_tavg_fields(id)%stream_number

   if (stream_number <= 0 .or. stream_number > max_avail_tavg_streams) then
     call document ('tavg_in_which_stream', 'id', id)
     write(stdout,*) '(tavg_in_which_stream) stream_number = ', stream_number
     exit_string = 'FATAL ERROR: not in any stream'
     call document ('tavg_in_which_stream', exit_string)
     call exit_POP (sigAbort,exit_string,out_unit=stdout)
   endif

   tavg_in_which_stream = stream_number

!-----------------------------------------------------------------------
!EOC

 end function tavg_in_which_stream

!***********************************************************************
!BOP
! !IROUTINE: set_in_tavg_contents
! !INTERFACE:

 function set_in_tavg_contents(id)

! !DESCRIPTION:
!  This function determines whether a tavg field has been set in
!  the input contents file and returns true if it has.  This function is
!  different from tavg_requested in that ltavg_on status is irrelevent.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      id                   ! id returned by the define function which
                           !   gives the location of the field

! !OUTPUT PARAMETERS:

   logical (log_kind) ::      &
      set_in_tavg_contents    ! result of checking whether the field has
                               !   been requested

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  check the buffer location - if not zero, then the field is in the
!  tavg contents file
!
!-----------------------------------------------------------------------

   if (id < 1 .or. id > num_avail_tavg_fields) then
     set_in_tavg_contents = .false.
   elseif (avail_tavg_fields(id)%buf_loc /= 0) then
     set_in_tavg_contents = .true.
   else
     set_in_tavg_contents = .false.
   endif

!-----------------------------------------------------------------------
!EOC

 end function set_in_tavg_contents

!***********************************************************************
!BOP
! !IROUTINE: tavg_id
! !INTERFACE:

 function tavg_id(short_name,quiet)

! !DESCRIPTION:
!  This function determines whether a tavg field has been defined
!  by some module.  If so, then the id of that field is returned.
!  This function does not cause a model exit if the field has not
!  been defined; error-checking must be done by the calling routine.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in)  :: &
      short_name                  ! the short name of the tavg field

   logical (log_kind), intent(in), optional :: &
      quiet                        ! do not print error message

! !OUTPUT PARAMETERS:

   integer (int_kind) :: &
      tavg_id               ! id of the tavg field, if it exists

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   integer (int_kind) ::  &
     id,                  &
     n

   logical (log_kind) ::  &
     msg

   id = 0

   srch_loop: do n=1,num_avail_tavg_fields
      if (trim(avail_tavg_fields(n)%short_name) == trim(short_name)) then
         id = n
         exit srch_loop
      endif
   end do srch_loop

   msg = .true.
   if (present(quiet)) then
     if (quiet) msg = .false.
   endif
     
   if (id == 0 .and. msg) then
      if (my_task == master_task)  & 
          write(stdout,*) 'Field ', trim(short_name), ' has not been defined.'
   endif

   tavg_id = id

!-----------------------------------------------------------------------
!EOC

 end function tavg_id

!***********************************************************************
!BOP
! !IROUTINE: tavg_create_suffix
! !INTERFACE:

 subroutine tavg_create_suffix(file_suffix,ns)

! !DESCRIPTION:
!  Creates a suffix to append to output filename based on frequency 
!  option and averaging interval.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      ns           ! stream index


! !OUTPUT PARAMETERS:

   character (char_len), intent(out) :: &
      file_suffix           ! suffix to append to root filename

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variable
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      cindx1, cindx2,    &! indices into character strings
      len_date            ! length of date string

   character (char_len) :: &
      char_temp            ! temp character space (for removing spaces)

   character (10) :: &
      cstep_beg,     &! beginning step  of this particular average
      cstep_end,     &! ending    step  of this particular average
      cdate           ! character string with yyyymmdd and optional 
                      ! separator (eg yyyy-mm-dd)

   character (4) :: &
      cyear_beg,    &! beginning year  of this particular average
      cyear_end      ! end       year  of this particular average

   character (2) :: &
      cmonth_beg,   &! beginning month of this particular average
      cmonth_end,   &! end       month of this particular average
      cday_beg,     &! beginning day   of this particular average
      cday_end       ! end       day   of this particular average

!-----------------------------------------------------------------------
!
!  start suffix with runid
!
!-----------------------------------------------------------------------

   file_suffix = char_blank
   cindx2 = len_trim(runid) + 1
   file_suffix(1:cindx2) = trim(runid)/&
                                       &/'.'
   cindx1 = cindx2 + 1
   
!-----------------------------------------------------------------------
!
!  extract beginning year, month, day or time step from beg_date
!  and determine end date
!
!-----------------------------------------------------------------------

   cdate = adjustl(beg_date)

   !***
   !*** use step numbers if tavg freq option is nstep
   !***

   cstep_beg  = trim(cdate) ! in case beg_date actually step number
   write(cstep_end,'(i10)') nsteps_total - 1
   cdate  = adjustl(cstep_end)
   cstep_end = trim(cdate)

   call time_stamp('last','ymd',date_string=cdate)  ! last date

   if (date_separator == ' ') then  ! no date separator
      cyear_beg  = beg_date(1:4)
      cmonth_beg = beg_date(5:6)
      cday_beg   = beg_date(7:8)
      cyear_end  = cdate(1:4)
      cmonth_end = cdate(5:6)
      cday_end   = cdate(7:8)
   else
      cyear_beg  = beg_date(1:4)
      cmonth_beg = beg_date(6:7)
      cday_beg   = beg_date(9:10)
      cyear_end  = cdate(1:4)
      cmonth_end = cdate(6:7)
      cday_end   = cdate(9:10)
   endif

!-----------------------------------------------------------------------
!
!  create time portion of suffix based on frequency option
!  note that concatenation operator split across lines to avoid
!   problems with some cpp preprocessors
!
!-----------------------------------------------------------------------

   select case (tavg_freq_iopt(ns))
   case (freq_opt_nyear)
      if (tavg_freq(ns) == 1) then
         cindx2 = cindx1 + 3
         file_suffix(cindx1:cindx2) = cyear_end
      else
         cindx2 = cindx1 + 8
         file_suffix(cindx1:cindx2) = cyear_beg/&
                                                &/'-'/&
                                                      &/cyear_end
      endif

   case (freq_opt_nmonth)
      if (tavg_freq(ns) == 1) then
         cindx2 = cindx1 + 5
         file_suffix(cindx1:cindx2) = cyear_end/&
                                    &/cmonth_end
      else
         cindx2 = cindx1 + 12
         file_suffix(cindx1:cindx2) = cyear_beg/&
                                    &/cmonth_beg/&
                                    &/'-'/&
                                    &/cyear_end/&
                                    &/cmonth_end
      endif

   case (freq_opt_nday)
      if (tavg_freq(ns) == 1) then
         cindx2 = cindx1 + 7
         file_suffix(cindx1:cindx2) = cyear_end/&
                                    &/cmonth_end/&
                                    &/cday_end
      else
         cindx2 = cindx1 + 16
         file_suffix(cindx1:cindx2) = cyear_beg/&
                                    &/cmonth_beg/&
                                    &/cday_beg/&
                                    &/'-'/&
                                    &/cyear_end/&
                                    &/cmonth_end/&
                                    &/cday_end
      endif

   case (freq_opt_nstep)
      cindx2 = cindx1 + len_trim(cstep_beg) + len_trim(cstep_end)
      file_suffix(cindx1:cindx2) = trim(cstep_beg)/&
                                 &/'-'/&
                                 &/trim(cstep_end)

   case default  ! use nstep for other options
      cindx2 = len_trim(cstep_beg) + len_trim(cstep_end) + 1
      file_suffix(cindx1:cindx2) = trim(cstep_beg)/&
                                 &/'-'/&
                                 &/trim(cstep_end)

   end select
 
!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_create_suffix

!***********************************************************************
!BOP
! !IROUTINE: tavg_create_suffix_ccsm
! !INTERFACE:

 subroutine tavg_create_suffix_ccsm(file_suffix,ns,date_string) 

! !DESCRIPTION:
!  Creates a suffix to append to output filename based on frequency 
!  option and averaging interval.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) ::  &
       ns

   character (*), intent(in), optional :: &
      date_string

! !OUTPUT PARAMETERS:

   character (char_len), intent(out) :: &
      file_suffix           ! suffix to append to root filename

!EOP
!BOC

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

   character (char_len) :: &
      char_temp,           &! temp character space
      ccsm_date_string


!-----------------------------------------------------------------------
!
!  clear character strings
!
!-----------------------------------------------------------------------

   file_suffix = char_blank
   char_temp   = char_blank


!-----------------------------------------------------------------------
!
!  for a ccsm tavg file, append a date/time string to the root name
!
!-----------------------------------------------------------------------

      select case (tavg_freq_iopt(ns))
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
 
      case (freq_opt_once)
        file_suffix = ''
        return 

      case default
        char_temp = 'ymds'
      end select

!-----------------------------------------------------------------------
!     override tavg_freq_iopt if this is a tavg restart file
!-----------------------------------------------------------------------

      if (present (date_string) ) then
         char_temp = trim(date_string)
      endif

      call ccsm_date_stamp (ccsm_date_string, char_temp)
 
      file_suffix = trim(ccsm_date_string)

 
!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_create_suffix_ccsm


!***********************************************************************
!BOP
! !IROUTINE: tavg_create_outfile_ccsm
! !INTERFACE:

 subroutine tavg_create_outfile_ccsm(tavg_outfile_orig,tavg_outfile,ns)

! !DESCRIPTION:
!    Forms the output filename for tavg RESTART history files.
!
!    Modifies the filename tavg_output such that it conforms to the CCSM4 requirements.
!    CCSM4 pop2 requires tavg RESTART history files to be of the form 
!    pathname/casename.rh.datestring
!
! !REVISION HISTORY:
!  same as module


! !INPUT PARAMETERS:

   character (char_len), intent(in) :: &
      tavg_outfile_orig           

   integer (int_kind), intent(in) ::  &
      ns                               ! stream number

! !INPUT/OUTPUT PARAMETERS:

   character (char_len), intent(inout) :: &
      tavg_outfile           

!-----------------------------------------------------------------------
!   NOTE: Assumptions for filename patterns:
!     
!     tavg_outfile ends in:
!        1) '.h' , if ns = 1
!        2) '.hN', N=1,2,...,9  OR '.h.string', if ns > 1 and base model
!        3) '.h.string' if extra-tracer modules
!     If your filename does not adhere to these CCSM conventions, the
!     model will abort.
!
!-----------------------------------------------------------------------
!EOP
!BOC

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------
   integer (i4) :: &
      iii,jjj       ! local indices

   character (char_len) :: &
      tavg_outfile_temp,   &! temp tavg output filename
      temp_string,         &
      string

!-----------------------------------------------------------------------
!    store original filename in tavg_outfile_temp
!-----------------------------------------------------------------------
     tavg_outfile_temp = trim(tavg_outfile_orig)
     string            = trim(tavg_outfile_orig)  ! shorter name
     iii = len_trim(tavg_outfile_temp)

!-----------------------------------------------------------------------
!    error checking
!-----------------------------------------------------------------------
     if (iii .lt. 1) then
       exit_string = 'FATAL ERROR: string length < 1 '
       call document ('tavg_create_outfile_ccsm', exit_string)
       call exit_POP (sigAbort,exit_string,out_unit=stdout)
     endif


!-----------------------------------------------------------------------
!    accounting for the possibility of multiple tavg streams,
!    replace the character ('h') with the string 'rh'
!-----------------------------------------------------------------------

     if (ns == 1) then 
         !-----------------------------------------------------
         !   assumes string ends in '.h' --  replace with '.rh'
         !-----------------------------------------------------
         tavg_outfile_temp(iii:iii+1) = 'rh'
     endif

     if (ns > 1) then
       if (string(iii-2:iii-1) == '.h') then  
         !-------------------------------------------------------
         !   assumes string ends in '.hN' -- replace with '.rhN'
         !-------------------------------------------------------
         tavg_outfile_temp(iii-1:iii) = 'rh'
         tavg_outfile_temp(iii+1:iii+1) = string(iii:iii)
       else
         !------------------------------------------------------------------
         !   assumes string ends in '.h.string' -- replace with '.rh.string'
         !------------------------------------------------------------------
         hist_search: do jjj=iii,3,-1
            if (string(jjj-2:jjj) == '.h.') then
              tavg_outfile_temp(jjj-1:jjj) = 'rh'
              tavg_outfile_temp(jjj+1:) = string(jjj:iii)
              exit hist_search  
            endif
         enddo hist_search
         if (trim(tavg_outfile_temp) == trim(string)) then
           exit_string = 'FATAL ERROR: string was not modified'
           call document ('tavg_create_outfile_ccsm', exit_string)
           call exit_POP (sigAbort,exit_string,out_unit=stdout)
         endif
       endif
     endif

!-----------------------------------------------------------------------
!    Fill tavg_outfile with modified filename
!-----------------------------------------------------------------------
     tavg_outfile = trim(tavg_outfile_temp)

!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_create_outfile_ccsm


!***********************************************************************
!BOP
! !IROUTINE: tavg_count_contents
! !INTERFACE:

 subroutine tavg_count_contents (num_fields,tavg_contents)
 
! !DESCRIPTION:
!  This routine counts the number of lines in the tavg_contents file
!
! !REVISION HISTORY:
!  same as module


! !INPUT PARAMETERS:
 
  character (char_len),intent(in)   :: tavg_contents
 
! !OUTPUT PARAMETERS:
 
  integer   (int_kind), intent(out) :: num_fields
 
!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

  integer   (int_kind)  :: nu
  integer   (int_kind)  :: ios
  integer   (int_kind)  :: mt
  integer   (int_kind)  :: grid_code_int
  character (char_len)  :: file_line
 

!---------------------------------------------------------------------
! get unit and open contents file
!---------------------------------------------------------------------
  call get_unit(nu)
  if (my_task == master_task) then
     open(nu, file=tavg_contents, status='old')
  endif

  num_fields = 0

  count_loop: do
 
!---------------------------------------------------------------------
!  read line from file, checking for end-of-file
!---------------------------------------------------------------------
   if (my_task == master_task) then
      read (nu,'(A)',iostat=ios) file_line
   endif

   call broadcast_scalar(ios, master_task)
   if (ios < 0) exit count_loop
   call broadcast_scalar(file_line, master_task)
 
!---------------------------------------------------------------------
!  increment 
!---------------------------------------------------------------------
   num_fields = num_fields + 1

  enddo count_loop
 
  if (my_task == master_task) then
     write(stdout,*) 'There are ',num_fields,  &
                    ' tavg fields requested via tavg_contents.'
     call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
  endif

!---------------------------------------------------------------------
! close file and release unit
!---------------------------------------------------------------------
 close(nu)
 call release_unit(nu)

!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_count_contents

!***********************************************************************
!BOP
! !IROUTINE: tavg_add_attrib_file_ccsm
! !INTERFACE:

 subroutine tavg_add_attrib_file_ccsm (tavg_file_desc)
 
! !DESCRIPTION:
!  This routine adds attributes to the ccsm tavg history file
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:
 type (datafile), intent(inout) :: tavg_file_desc    ! IO file descriptor
 
!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
 
 character (char_len) ::   &
    start_time,    & 
    current_date,  &
    current_time,  &
    cell_methods,  &
    calendar


 call add_attrib_file(tavg_file_desc, 'contents', 'Diagnostic and Prognostic Variables')
 call add_attrib_file(tavg_file_desc, 'source', 'CCSM POP2, the CCSM Ocean Component')
 call add_attrib_file(tavg_file_desc, 'revision', &
   '$Id: tavg.F90 56176 2013-12-20 18:35:46Z mlevy@ucar.edu $')

 if (allow_leapyear) then
    write(calendar,'(a,i5,a,i5,a)') &
       ' Leap years allowed. Normal years have',days_in_norm_year, &
       ' days. Leap years have ' ,days_in_leap_year, ' days.'
 else
    write(calendar,'(a,i5,a)')'All years have exactly', days_in_norm_year, ' days.'
 endif
 call add_attrib_file(tavg_file_desc, 'calendar', trim(calendar))

 
 if (my_task.eq.master_task) then
   call date_and_time(date=current_date, time=current_time)
 end if
 call broadcast_scalar(current_date, master_task)
 call broadcast_scalar(current_time, master_task)
 start_time = char_blank
 write(start_time,1000) current_date(1:4), current_date(5:6),  &
                        current_date(7:8), current_time(1:2),  &
                        current_time(3:4), current_time(5:8)
 call add_attrib_file(tavg_file_desc, 'start_time', trim(start_time))


 cell_methods = char_blank
 cell_methods = 'cell_methods = time: mean ==> the variable values ' /&
     &/  'are averaged over the time interval between the previous '      /&
     &/  'time coordinate and the current one.          '                 /&
     &/  'cell_methods  absent  ==> the variable values '                 /&
     &/  'are at the time given by the current time coordinate. '
 call add_attrib_file(tavg_file_desc, 'cell_methods', trim(cell_methods))

1000  format('This dataset was created on ',a,'-',a,'-',a,' at ',a,':',a,':',a)

!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_add_attrib_file_ccsm 
 

!***********************************************************************
!BOP
! !IROUTINE: tavg_add_attrib_io_field_ccsm
! !INTERFACE:

 subroutine tavg_add_attrib_io_field_ccsm (tavg_field,nfield)
 
! !DESCRIPTION:
!  This routine adds ccsm-required attributes to ccsm tavg fields
!  cell_methods, scale_factor, and _FillValue
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:
 integer (i4), intent(in) ::  &
    nfield

! !INPUT/OUTPUT PARAMETERS:
 type (io_field_desc), intent(inout) ::  &
    tavg_field    ! IO file descriptor
 
!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

 integer (i4)         ::  &
    method_integer
 character (char_len) ::  &
    method_string

   method_integer = avail_tavg_fields(nfield)%method
   call tavg_get_cell_method_string (method_integer,method_string)
   call add_attrib_io_field(tavg_field, 'cell_methods', trim(method_string))

   if (avail_tavg_fields(nfield)%scale_factor /= undefined_nf) then
     call add_attrib_io_field(tavg_field, 'scale_factor',avail_tavg_fields(nfield)%scale_factor)
   endif

!*** note: missing_value is a deprecated feature in CF1.4, and hence nco 4 versions, but
!          is added here because other software packages may require it
   call add_attrib_io_field(tavg_field,'_FillValue',   avail_tavg_fields(nfield)%fill_value )
   call add_attrib_io_field(tavg_field,'missing_value',avail_tavg_fields(nfield)%fill_value )

!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_add_attrib_io_field_ccsm 

!***********************************************************************
!BOP
! !IROUTINE: tavg_get_cell_method_string
! !INTERFACE:

 subroutine tavg_get_cell_method_string (method_integer,method_string)
 
! !DESCRIPTION:
!  This routine determines the string associated with the cell_method
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:
 integer (i4), intent(in) ::  &
    method_integer

! !OUTPUT PARAMETERS:
 character (char_len), intent(out) ::  &
    method_string

!EOP
!BOC

   method_string = char_blank

   select case (method_integer)

     case(tavg_method_avg)
        method_string='time: mean'

     case(tavg_method_min)
        method_string='time: minimum'

     case(tavg_method_max)
        method_string='time: maximum'

     case(tavg_method_qflux)
        method_string='time: mean'

     case default
        write(exit_string,'(a,1x,i4)') 'FATAL ERROR:  unknown method = ', method_integer
        call document ('tavg_add_attrib_io_field_ccsm', exit_string)
        call exit_POP (sigAbort,exit_string,out_unit=stdout)

   end select

!EOC

 end subroutine tavg_get_cell_method_string 

!***********************************************************************
!BOP
! !IROUTINE: tavg_construct_ccsm_coordinates
! !INTERFACE:

 subroutine tavg_construct_ccsm_coordinates(tavg_file_desc,ns)

! !DESCRIPTION:
!  This routine defines the netCDF coordinates that are used in the
!  ccsm netCDF output tavg files
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:
   integer (int_kind), intent(in) ::  &
      ns                ! tavg streams index

! !INPUT/OUTPUT PARAMETERS:
   type (datafile), intent(inout) :: tavg_file_desc    ! IO file descriptor

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   real (rtavg), dimension(km) ::  &
      ZT_R,      &! single/double precision array
      ZW_R,      &! single/double precision array
      ZW_BOT_R    ! single/double precision array

   real (rtavg), dimension(0:km) ::  &
      MOC_Z_R

   real (rtavg), dimension(1000) ::  &
      LAT_AUX_GRID_R

   integer (int_kind) ::  &
      ii, n, ndims   

   save  



   ii=0

   !*** z_t
   ii=ii+1
   ZT_R = zt 
   ccsm_coordinates(ii,ns) = construct_io_field('z_t',zt_dim,                 &
                         long_name='depth from surface to midpoint of layer', &
                         units    ='centimeters',                             &
#ifdef TAVG_R8
                         d1d_array =ZT_R)
#else
                         r1d_array =ZT_R)
#endif

   call add_attrib_io_field(ccsm_coordinates(ii,ns), 'positive', 'down')
   call add_attrib_io_field(ccsm_coordinates(ii,ns), 'valid_min', ZT_R(1))
   call add_attrib_io_field(ccsm_coordinates(ii,ns), 'valid_max', ZT_R(km))

   !*** z_t
   ii=ii+1
   ZT_150m_R = zt(1:zt_150m_levs)
   ccsm_coordinates(ii,ns) = construct_io_field('z_t_150m',zt_150m_dim,       &
                         long_name='depth from surface to midpoint of layer', &
                         units    ='centimeters',                             &
#ifdef TAVG_R8
                         d1d_array =ZT_150m_R)
#else
                         r1d_array =ZT_150m_R)
#endif

   call add_attrib_io_field(ccsm_coordinates(ii,ns), 'positive', 'down')
   call add_attrib_io_field(ccsm_coordinates(ii,ns), 'valid_min', ZT_150m_R(1))
   call add_attrib_io_field(ccsm_coordinates(ii,ns), 'valid_max', ZT_150m_R(zt_150m_levs))

   !*** z_w
   ii=ii+1
   ZW_R(1) = c0 
   ZW_R(2:km) = zw(1:km-1)
   ccsm_coordinates(ii,ns) = construct_io_field('z_w',zw_dim,                 &
                         long_name='depth from surface to top of layer',      &
                         units    ='centimeters',                             &
#ifdef TAVG_R8
                         d1d_array =ZW_R)
#else
                         r1d_array =ZW_R)
#endif

   call add_attrib_io_field(ccsm_coordinates(ii,ns), 'positive', 'down')
   call add_attrib_io_field(ccsm_coordinates(ii,ns), 'valid_min', ZW_R(1 ))
   call add_attrib_io_field(ccsm_coordinates(ii,ns), 'valid_max', ZW_R(km))

   !*** z_w_top
   ii=ii+1
   ZW_R(1) = c0  ! same as z_w
   ZW_R(2:km) = zw(1:km-1)
   ccsm_coordinates(ii,ns) = construct_io_field('z_w_top',zw_dim_top,         &
                         long_name='depth from surface to top of layer',      &
                         units    ='centimeters',                             &
#ifdef TAVG_R8
                         d1d_array =ZW_R)
#else
                         r1d_array =ZW_R)
#endif

   call add_attrib_io_field(ccsm_coordinates(ii,ns), 'positive', 'down')
   call add_attrib_io_field(ccsm_coordinates(ii,ns), 'valid_min', ZW_R(1 ))
   call add_attrib_io_field(ccsm_coordinates(ii,ns), 'valid_max', ZW_R(km))

   !*** z_w_bot
   ii=ii+1
   ZW_BOT_R(1:km) = zw(1:km)
   ccsm_coordinates(ii,ns) = construct_io_field('z_w_bot',zw_dim_bot,         &
                         long_name='depth from surface to bottom of layer',   &
                         units    ='centimeters',                             &
#ifdef TAVG_R8
                         d1d_array =ZW_BOT_R)
#else
                         r1d_array =ZW_BOT_R)
#endif

   call add_attrib_io_field(ccsm_coordinates(ii,ns), 'positive', 'down')
   call add_attrib_io_field(ccsm_coordinates(ii,ns), 'valid_min', ZW_BOT_R(1 ))
   call add_attrib_io_field(ccsm_coordinates(ii,ns), 'valid_max', ZW_BOT_R(km))


   if (ltavg_moc_diags(ns) .or. ltavg_n_heat_trans(ns)  .or. ltavg_n_salt_trans(ns)) then
     !*** lat_aux_grid
     ii=ii+1

     if (n_lat_aux_grid+1 > 1000) then
        exit_string = 'FATAL ERROR:  must increase dimension of LAT_AUX_GRID_R'
        call document ('tavg_construct_ccsm_coordinates', exit_string)
        call exit_POP (sigAbort,exit_string,out_unit=stdout)
     endif

     LAT_AUX_GRID_R = 0
     LAT_AUX_GRID_R(1:n_lat_aux_grid+1) = lat_aux_edge(1:n_lat_aux_grid+1)
     ccsm_coordinates(ii,ns) = construct_io_field('lat_aux_grid',lat_aux_grid_dim, &
                           long_name='latitude grid for transport diagnostics', &
                           units    ='degrees_north',                           &
#ifdef TAVG_R8
                           d1d_array =LAT_AUX_GRID_R(1:n_lat_aux_grid+1))
#else
                           r1d_array =LAT_AUX_GRID_R(1:n_lat_aux_grid+1))
#endif
     call add_attrib_io_field(ccsm_coordinates(ii,ns), 'valid_min', LAT_AUX_GRID_R(1 ))
     call add_attrib_io_field(ccsm_coordinates(ii,ns), 'valid_max', LAT_AUX_GRID_R(n_lat_aux_grid+1))
   endif

   if (ltavg_moc_diags(ns)) then
     !*** moc_z
     ii=ii+1

     MOC_Z_R(0) = c0 
     MOC_Z_R(1:km) = zw(1:km)
     ccsm_coordinates(ii,ns) = construct_io_field('moc_z',moc_z_dim,               &
                           long_name='depth from surface to top of layer',      &
                           units    ='centimeters',                             &
#ifdef TAVG_R8
                           d1d_array =MOC_Z_R)
#else
                           r1d_array =MOC_Z_R)
#endif

     call add_attrib_io_field(ccsm_coordinates(ii,ns), 'positive', 'down')
     call add_attrib_io_field(ccsm_coordinates(ii,ns), 'valid_min', MOC_Z_R(0 ))
     call add_attrib_io_field(ccsm_coordinates(ii,ns), 'valid_max', MOC_Z_R(km))

   endif ! ltavg_moc_diags


   ! after all coordinates are defined, define the total number of coordinates
   num_ccsm_coordinates = ii
 
   if (num_ccsm_coordinates > max_num_ccsm_coordinates) then
     exit_string = 'FATAL ERROR:  reset max_num_ccsm_coordinates'
     call document ('tavg_construct_ccsm_coordinates', exit_string)
     call exit_POP (sigAbort,exit_string,out_unit=stdout)
   endif



!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_construct_ccsm_coordinates
 
!***********************************************************************
!BOP
! !IROUTINE: tavg_construct_ccsm_time_invar
! !INTERFACE:

 subroutine tavg_construct_ccsm_time_invar(tavg_file_desc,ns)

! !DESCRIPTION:
!  This routine defines the netCDF time-invariant variables that are 
!  used in the ccsm netCDF output tavg files
!
! !REVISION HISTORY:
!  same as module


! !INPUT/OUTPUT PARAMETERS:
   type (datafile), intent(inout) ::  &
      tavg_file_desc    ! IO file descriptor

! !INPUT PARAMETERS:
   integer (int_kind), intent(in) ::  &
      ns                ! tavg streams index

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: ii, num

   real (rtavg), dimension(km)   ::  &
      DZ_R

   real (rtavg), dimension(0:km-1) ::  &
      DZW_R

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) ::  &
      ULON_DEG, ULAT_DEG

   integer (i4)       :: fill_value_i
   real (r8)          :: fill_value_d

   logical (log_kind) ::  &
      full_header

   save

   fill_value_i = undefined_nf_int
   fill_value_d = undefined_nf_r8

   ii=0

   !*** dz
   ii=ii+1
   DZ_R = dz 
   ccsm_time_invar(ii,ns) = construct_io_field('dz',zt_dim,                 &
                                long_name='thickness of layer k',           &
                                units    ='centimeters',                    &
#ifdef TAVG_R8
                                d1d_array =DZ_R)
#else
                                r1d_array =DZ_R)
#endif

   !*** dzw
   ii=ii+1
   DZW_R(0:) = dzw(0:km-1)
   ccsm_time_invar(ii,ns) = construct_io_field('dzw',zw_dim,                       &
                                long_name='midpoint of k to midpoint of k+1',      &
                                units    ='centimeters',                           &
#ifdef TAVG_R8
                                d1d_array =DZW_R)
#else
                                r1d_array =DZW_R)
#endif

!-----------------------------------------------------------------------
!
!  only write the following time-invariant fields if requested
!
!-----------------------------------------------------------------------

   if (.not. tavg_streams(ns)%ltavg_one_time_header) then
     full_header = .true.
   else 
     if (tavg_streams(ns)%ltavg_first_header) then
       full_header = .true.
     else
       full_header = .false.
     endif
   endif
 
   if (full_header) then

   !*** ULONG
   ii=ii+1

    ULON_DEG = ULON*radian
    ccsm_time_invar(ii,ns) = construct_io_field(  &
        'ULONG', dim1=i_dim, dim2=j_dim,              &
         long_name='array of u-grid longitudes',      &
         units    ='degrees_east',                    &
         d2d_array =ULON_DEG(:,:,:) )

   !*** ULAT
   ii=ii+1

    ULAT_DEG = ULAT*radian
    ccsm_time_invar(ii,ns) = construct_io_field(  &
        'ULAT', dim1=i_dim, dim2=j_dim,               &
         long_name='array of u-grid latitudes',       &
         units    ='degrees_north',                   &
         d2d_array =ULAT_DEG(:,:,:) )

   !*** TLONG (degrees)
   ii=ii+1

    ccsm_time_invar(ii,ns) = construct_io_field(  &
        'TLONG', dim1=i_dim, dim2=j_dim,              &
         long_name='array of t-grid longitudes',      &
         units    ='degrees_east',                    &
         d2d_array =TLOND(:,:,:) )

   !*** TLAT (degrees)
   ii=ii+1

    ccsm_time_invar(ii,ns) = construct_io_field(  &
        'TLAT', dim1=i_dim, dim2=j_dim,               &
         long_name='array of t-grid latitudes',       &
         units    ='degrees_north',                   &
         d2d_array =TLATD(:,:,:) )

   !*** KMT
   ii=ii+1

    ccsm_time_invar(ii,ns) = construct_io_field(          &
        'KMT', dim1=i_dim, dim2=j_dim,                        &
         long_name='k Index of Deepest Grid Cell on T Grid',  &
         coordinates = "TLONG TLAT",                          &
         i2d_array =KMT(:,:,:) )

   !*** KMU
   ii=ii+1

    ccsm_time_invar(ii,ns) = construct_io_field(          &
        'KMU', dim1=i_dim, dim2=j_dim,                        &
         long_name='k Index of Deepest Grid Cell on U Grid',  &
         coordinates = "ULONG ULAT",                          &
         i2d_array =KMU(:,:,:) )


   !*** REGION_MASK
   ii=ii+1

    ccsm_time_invar(ii,ns) = construct_io_field(          &
        'REGION_MASK', dim1=i_dim, dim2=j_dim,                &
         long_name='basin index number (signed integers)',    &
         coordinates = "TLONG TLAT",                          &
         i2d_array =REGION_MASK(:,:,:) )

   !*** UAREA
   ii=ii+1

    ccsm_time_invar(ii,ns) = construct_io_field(          &
        'UAREA', dim1=i_dim, dim2=j_dim,                      &
         long_name='area of U cells',                         &
         units    ='centimeter^2',                            &
         coordinates = "ULONG ULAT",                          &
         d2d_array =UAREA(:,:,:) )

   !*** TAREA
   ii=ii+1

    ccsm_time_invar(ii,ns) = construct_io_field(          &
        'TAREA', dim1=i_dim, dim2=j_dim,                      &
         long_name='area of T cells',                         &
         units    ='centimeter^2',                            &
         coordinates = "TLONG TLAT",                          &
         d2d_array =TAREA(:,:,:) )

   !*** HU
   ii=ii+1

    ccsm_time_invar(ii,ns) = construct_io_field(          &
        'HU', dim1=i_dim, dim2=j_dim,                         &
         long_name='ocean depth at U points',                 &
         units='centimeter',                                  &
         coordinates = "ULONG ULAT",                          &
         d2d_array =HU(:,:,:) )

   !*** HT
   ii=ii+1

    ccsm_time_invar(ii,ns) = construct_io_field(          &
        'HT', dim1=i_dim, dim2=j_dim,                         &
         long_name='ocean depth at T points',                 &
         units='centimeter',                                  &
         coordinates = "TLONG TLAT",                          &
         d2d_array =HT(:,:,:) )

   !*** DXU
   ii=ii+1

    ccsm_time_invar(ii,ns) = construct_io_field(          &
        'DXU', dim1=i_dim, dim2=j_dim,                        &
         long_name='x-spacing centered at U points',          &
         units='centimeters',                                 &
         coordinates = "ULONG ULAT",                          &
         d2d_array =DXU(:,:,:) )

   !*** DYU
   ii=ii+1

    ccsm_time_invar(ii,ns) = construct_io_field(          &
        'DYU', dim1=i_dim, dim2=j_dim,                        &
         long_name='y-spacing centered at U points',          &
         units='centimeters',                                 &
         coordinates = "ULONG ULAT",                          &
         d2d_array =DYU(:,:,:) )

   !*** DXT
   ii=ii+1

    ccsm_time_invar(ii,ns) = construct_io_field(          &
        'DXT', dim1=i_dim, dim2=j_dim,                        &
         long_name='x-spacing centered at T points',          &
         units='centimeters',                                 &
         coordinates = "TLONG TLAT",                          &
         d2d_array =DXT(:,:,:) )

   !*** DYT
   ii=ii+1

    ccsm_time_invar(ii,ns) = construct_io_field(          &
        'DYT', dim1=i_dim, dim2=j_dim,                        &
         long_name='y-spacing centered at T points',          &
         units='centimeters',                                 &
         coordinates = "TLONG TLAT",                          &
         d2d_array =DYT(:,:,:) )

   !*** HTN
   ii=ii+1

    ccsm_time_invar(ii,ns) = construct_io_field(          &
        'HTN', dim1=i_dim, dim2=j_dim,                        &
         long_name='cell widths on North sides of T cell',    &
         units='centimeters',                                 &
         coordinates = "TLONG TLAT",                          &
         d2d_array =HTN(:,:,:) )

   !*** HTE
   ii=ii+1

    ccsm_time_invar(ii,ns) = construct_io_field(          &
        'HTE', dim1=i_dim, dim2=j_dim,                        &
         long_name='cell widths on East sides of T cell',     &
         units='centimeters',                                 &
         coordinates = "TLONG TLAT",                          &
         d2d_array =HTE(:,:,:) )

   !*** HUS
   ii=ii+1

    ccsm_time_invar(ii,ns) = construct_io_field(          &
        'HUS', dim1=i_dim, dim2=j_dim,                        &
         long_name='cell widths on South sides of U cell',    &
         units='centimeters',                                 &
         coordinates = "ULONG ULAT",                          &
         d2d_array =HUS(:,:,:) )

   !*** HUW
   ii=ii+1

    ccsm_time_invar(ii,ns) = construct_io_field(          &
        'HUW', dim1=i_dim, dim2=j_dim,                        &
         long_name='cell widths on West sides of U cell',     &
         units='centimeters',                                 &
         coordinates = "ULONG ULAT",                          &
         d2d_array =HUW(:,:,:) )

   !*** ANGLE
   ii=ii+1

    ccsm_time_invar(ii,ns) = construct_io_field(          &
        'ANGLE', dim1=i_dim, dim2=j_dim,                      &
         long_name='angle grid makes with latitude line',     &
         units='radians',                                     &
         coordinates = "ULONG ULAT",                          &
         d2d_array =ANGLE(:,:,:) )

   !*** ANGLET
   ii=ii+1

    ccsm_time_invar(ii,ns) = construct_io_field(               &
        'ANGLET', dim1=i_dim, dim2=j_dim,                          &
         long_name='angle grid makes with latitude line on T grid',&
         units='radians',                                          &
         coordinates = "TLONG TLAT",                               &
         d2d_array =ANGLET(:,:,:) )

   endif ! full_header


!-----------------------------------------------------------------------
!
!   after all time-invariant arrays are defined, define the total number
!
!-----------------------------------------------------------------------
   num_ccsm_time_invar(ns) = ii

   if (num_ccsm_time_invar(ns)  > max_num_ccsm_time_invar) then
     exit_string = 'FATAL ERROR:  reset max_num_ccsm_time_invar -- too small'
     call document ('tavg_construct_ccsm_time_invar', exit_string)
     call exit_POP (sigAbort,exit_string,out_unit=stdout)
   endif


!-----------------------------------------------------------------------
!
!   Finally, add _FillValue and missing_value attributes to all fields
!   NOTE: missing_value is identical to _FillValue
!
!-----------------------------------------------------------------------

   do num = 1, num_ccsm_time_invar(ns)
     if (trim(ccsm_time_invar(num,ns)%short_name) == 'REGION_MASK'  .or.   &
         trim(ccsm_time_invar(num,ns)%short_name) == 'KMT'          .or.   &
         trim(ccsm_time_invar(num,ns)%short_name) == 'KMU'                 ) then
         call add_attrib_io_field(ccsm_time_invar(num,ns),'_FillValue'   ,fill_value_i )
         call add_attrib_io_field(ccsm_time_invar(num,ns),'missing_value',fill_value_i )
     elseif(trim(ccsm_time_invar(num,ns)%short_name) == 'dz'  .or.   &
         trim(ccsm_time_invar(num,ns)%short_name) == 'dzw'                 ) then
         call add_attrib_io_field(ccsm_time_invar(num,ns),'_FillValue'   ,undefined_nf )
         call add_attrib_io_field(ccsm_time_invar(num,ns),'missing_value',undefined_nf )
     else
         call add_attrib_io_field(ccsm_time_invar(num,ns),'_FillValue'   ,fill_value_d )
         call add_attrib_io_field(ccsm_time_invar(num,ns),'missing_value',fill_value_d )
     endif
   enddo

!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_construct_ccsm_time_invar

!***********************************************************************
!BOP
! !IROUTINE: tavg_construct_ccsm_scalars
! !INTERFACE:

 subroutine tavg_construct_ccsm_scalars(tavg_file_desc,ns)

! !DESCRIPTION:
!  This routine defines the netCDF scalars that are used in the
!  ccsm netCDF output tavg files
!
! !REVISION HISTORY:
!  same as module


! !INPUT/OUTPUT PARAMETERS:
   type (datafile), intent(inout) :: tavg_file_desc    ! IO file descriptor

! !INPUT PARAMETERS:
   integer (int_kind), intent(in) ::  &
      ns                ! tavg streams index

!EOP
!BOC

   real (r8) ::  &
      d0d_days_in_norm_year,  &
      d0d_days_in_leap_year,  &
      d0d_nsurface_t,         &
      d0d_nsurface_u

   integer (int_kind) ::  &
      ii,n

   logical (log_kind) ::  &
      full_header

   save  


   ii=0

!-----------------------------------------------------------------------
!
!  only write the following scalar fields if requested
!
!-----------------------------------------------------------------------

   if (.not. tavg_streams(ns)%ltavg_one_time_header) then
     full_header = .true.
   else 
     if (tavg_streams(ns)%ltavg_first_header) then
       full_header = .true.
     else
       full_header = .false.
     endif
   endif
 
   if (full_header) then

   !*** days_in_norm_year
   d0d_days_in_norm_year = days_in_norm_year
   ii=ii+1
   ccsm_scalars(ii,ns) = construct_io_field('days_in_norm_year',      &
                         long_name='Calendar Length',              &
                         units    ='days',                         &
                         d0d_array =d0d_days_in_norm_year)

   !*** grav
   ii=ii+1
   ccsm_scalars(ii,ns) = construct_io_field('grav',                   &
                         long_name='Acceleration Due to Gravity',  &
                         units    ='centimeter/s^2',               &
                         d0d_array =grav)

   !*** omega
   ii=ii+1
   ccsm_scalars(ii,ns) = construct_io_field('omega',                   &
                         long_name='Earths Angular Velocity',       &
                         units    ='1/second',                      &
                         d0d_array =omega)

   !*** radius
   ii=ii+1
   ccsm_scalars(ii,ns) = construct_io_field('radius',   &
                         long_name='Earths Radius',  &
                         units    ='centimeters',    &
                         d0d_array =radius)


   !*** cp_sw
   ii=ii+1
   ccsm_scalars(ii,ns) = construct_io_field('cp_sw',                 &
                         long_name='Specific Heat of Sea Water',  &
                         units    ='erg/g/K',                     &
                         d0d_array =cp_sw)

   !*** sound
   ii=ii+1
   ccsm_scalars(ii,ns) = construct_io_field('sound',                 &
                         long_name='Speed of Sound',              &
                         units    ='centimeter/s',                &
                         d0d_array =sound)

   !*** vonkar
   ii=ii+1
   ccsm_scalars(ii,ns) = construct_io_field('vonkar',                &
                         long_name='von Karman Constant',         &
                         d0d_array =vonkar)

   !*** cp_air
   ii=ii+1
   ccsm_scalars(ii,ns) = construct_io_field('cp_air',                &
                         long_name='Heat Capacity of Air',        &
                         units    ='joule/kg/degK',               &
                         d0d_array =cp_air)

   !*** rho_air
   ii=ii+1
   ccsm_scalars(ii,ns) = construct_io_field('rho_air',               &
                         long_name='Ambient Air Density',         &
                         units    ='kg/m^3',                      &
                         d0d_array =rho_air)

   !*** rho_sw
   ii=ii+1
   ccsm_scalars(ii,ns) = construct_io_field('rho_sw',                &
                         long_name='Density of Sea Water',        &
                         units    ='gram/centimeter^3',           &
                         d0d_array =rho_sw)

   !*** rho_fw
   ii=ii+1
   ccsm_scalars(ii,ns) = construct_io_field('rho_fw',                &
                         long_name='Density of Fresh Water',      &
                         units    ='gram/centimeter^3',           &
                         d0d_array =rho_fw)

   !*** stefan_boltzmann
   ii=ii+1
   ccsm_scalars(ii,ns) = construct_io_field('stefan_boltzmann',      &
                         long_name='Stefan-Boltzmann Constant',    &
                         units    ='watt/m^2/degK^4',             &
                         d0d_array =stefan_boltzmann)

   !*** latent_heat_vapor_mks
   ii=ii+1
   ccsm_scalars(ii,ns) = construct_io_field('latent_heat_vapor',  &
                         long_name='Latent Heat of Vaporization', &
                         units    ='J/kg',                        &
                         d0d_array =latent_heat_vapor_mks)

   !*** latent_heat_fusion
   ii=ii+1
   ccsm_scalars(ii,ns) = construct_io_field('latent_heat_fusion', &
                         long_name='Latent Heat of Fusion',       &
                         units    ='erg/g',                       &
                         d0d_array =latent_heat_fusion)

   !*** ocn_ref_salinity
   ii=ii+1
   ccsm_scalars(ii,ns) = construct_io_field('ocn_ref_salinity',      &
                         long_name='Ocean Reference Salinity',    &
                         units    ='g/kg',                        &
                         d0d_array =ocn_ref_salinity)

   !*** sea_ice_salinity
   ii=ii+1
   ccsm_scalars(ii,ns) = construct_io_field('sea_ice_salinity',      &
                         long_name='Salinity of Sea Ice',         &
                         units    ='g/kg',                        &
                         d0d_array =sea_ice_salinity)

   !*** T0_Kelvin
   ii=ii+1
   ccsm_scalars(ii,ns) = construct_io_field('T0_Kelvin',             &
                         long_name='Zero Point for Celsius',      &
                         units    ='degK',                        &
                         d0d_array =T0_Kelvin)

   !*** salt_to_ppt
   ii=ii+1
   ccsm_scalars(ii,ns) = construct_io_field('salt_to_ppt',           &
                         long_name='Convert Salt in gram/gram to g/kg',&
                         d0d_array =salt_to_ppt)

   !*** ppt_to_salt
   ii=ii+1
   ccsm_scalars(ii,ns) = construct_io_field('ppt_to_salt',           &
                         long_name='Convert Salt in g/kg to gram/gram',&
                         d0d_array =ppt_to_salt)
   !*** mass_to_Sv
   ii=ii+1
   ccsm_scalars(ii,ns) = construct_io_field('mass_to_Sv',           &
                         long_name='Convert Mass Flux to Sverdrups',&
                         d0d_array =mass_to_Sv)

   !*** heat_to_PW
   ii=ii+1
   ccsm_scalars(ii,ns) = construct_io_field('heat_to_PW',           &
                         long_name='Convert Heat Flux to Petawatts', &
                         d0d_array =heat_to_PW)

   !*** salt_to_Svppt
   ii=ii+1
   ccsm_scalars(ii,ns) = construct_io_field('salt_to_Svppt',           &
                         long_name='Convert Salt Flux to Sverdrups*g/kg', & 
                         d0d_array =salt_to_Svppt)
   !*** salt_to_mmday
   ii=ii+1
   ccsm_scalars(ii,ns) = construct_io_field('salt_to_mmday',           &
                         long_name='Convert Salt to Water (millimeters/day)', &
                         d0d_array =salt_to_mmday)

   !*** momentum_factor
   ii=ii+1
   ccsm_scalars(ii,ns) = construct_io_field('momentum_factor',           &
                         long_name='Convert Windstress to Velocity Flux', &
                         d0d_array =momentum_factor)

   !*** hflux_factor
   ii=ii+1
   ccsm_scalars(ii,ns) = construct_io_field('hflux_factor',           &
        long_name='Convert Heat and Solar Flux to Temperature Flux',&
                         d0d_array =hflux_factor)

   !*** fwflux_factor
   ii=ii+1
   ccsm_scalars(ii,ns) = construct_io_field('fwflux_factor',           &
  long_name='Convert Net Fresh Water Flux to Salt Flux (in model units)', &
                         d0d_array =fwflux_factor)
 
   !*** salinity_factor
   ii=ii+1
   ccsm_scalars(ii,ns) = construct_io_field('salinity_factor',           &
  long_name=' ', &
                         d0d_array =salinity_factor)
 
   !*** sflux_factor
   ii=ii+1
   ccsm_scalars(ii,ns) = construct_io_field('sflux_factor',           &
  long_name='Convert Salt Flux to Salt Flux (in model units)', &
                         d0d_array =sflux_factor)
 
   !*** nsurface_t
   ii=ii+1
   d0d_nsurface_t = nsurface_t
   ccsm_scalars(ii,ns) = construct_io_field('nsurface_t',           &
  long_name='Number of Ocean T Points at Surface', &
                         d0d_array =d0d_nsurface_t)
 
   !*** nsurface_u
   ii=ii+1
   d0d_nsurface_u = nsurface_u
   ccsm_scalars(ii,ns) = construct_io_field('nsurface_u',           &
  long_name='Number of Ocean U Points at Surface', &
                         d0d_array =d0d_nsurface_u)
 
   endif ! full_header
 
   ! after all scalars are defined, define the total number of scalars
   num_ccsm_scalars(ns) = ii

   if (num_ccsm_scalars(ns) > max_num_ccsm_scalars) then
     exit_string = 'FATAL ERROR:  reset num_ccsm_scalars -- too small'
     call document ('tavg_construct_ccsm_scalars', exit_string)
     call exit_POP (sigAbort,exit_string,out_unit=stdout)
   endif


!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_construct_ccsm_scalars

!***********************************************************************
!BOP
! !IROUTINE: tavg_define_labels_ccsm
! !INTERFACE:

 subroutine tavg_define_labels_ccsm(tavg_file_desc,ns)

! !DESCRIPTION:
!  This routine defines the netCDF labels that are used in the
!  ccsm netCDF output tavg files
!
! !REVISION HISTORY:
!  same as module


! !INPUT/OUTPUT PARAMETERS:
   type (datafile), intent(inout) :: tavg_file_desc    ! IO file descriptor

! !INPUT PARAMETERS:
   integer (int_kind), intent(in) ::  &
      ns                ! tavg streams index
!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::  &
      id, ii, n, ndims
 
   integer (i4) ::  &
      nstd_field_id
         

   save  

!-----------------------------------------------------------------------
!
! Note that avail_tavg_labels is type tavg_field_desc_ccsm, not io_field_desc
!
!-----------------------------------------------------------------------
   ii = 0
   num_avail_tavg_labels = 0

   if (ltavg_moc_diags(ns)) then
     ii = ii + 1
     ndims = 2
     io_dims_labels(1,ii)=nchar_dim
     io_dims_labels(2,ii)=moc_comp_dim
     call define_tavg_field(id, 'moc_components',ndims,            &
                            long_name='MOC component names',       &
                            nftype='char',                         &
                            nstd_fields=avail_tavg_labels,         &
                            num_nstd_fields=num_avail_tavg_labels, &
                            max_nstd_fields=max_avail_tavg_labels  )
   endif

   if (ltavg_n_heat_trans(ns) .or. ltavg_n_salt_trans(ns)) then
     ii = ii + 1
     ndims = 2
     io_dims_labels(1,ii)=nchar_dim
     io_dims_labels(2,ii)=transport_comp_dim
     call define_tavg_field(id, 'transport_components',ndims,      &
                            long_name='T,S transport components',  &
                            nftype='char',                         &
                            nstd_fields=avail_tavg_labels,         &
                            num_nstd_fields=num_avail_tavg_labels, &
                            max_nstd_fields=max_avail_tavg_labels  )
   endif

   if (ltavg_moc_diags(ns) .or. ltavg_n_heat_trans(ns) .or. ltavg_n_salt_trans(ns)) then
     ii = ii + 1
     ndims = 2
     io_dims_labels(1,ii)=nchar_dim
     io_dims_labels(2,ii)=transport_reg_dim
     call define_tavg_field(id, 'transport_regions',ndims,         &
                            long_name='regions for all transport diagnostics',  &
                            nftype='char',                         &
                            nstd_fields=avail_tavg_labels,         &
                            num_nstd_fields=num_avail_tavg_labels, &
                            max_nstd_fields=max_avail_tavg_labels  )
   endif

   if (ii /= num_avail_tavg_labels) then
      write (stdout,*) 'ii = ', ii
      write (stdout,*) 'num_avail_tavg_labels = ', num_avail_tavg_labels
      exit_string = 'FATAL ERROR:  mismatch in number of labels'
      call document ('tavg_define_labels_ccsm', exit_string)
      call exit_POP (sigAbort,exit_string,out_unit=stdout)
   endif

   do n=1,num_avail_tavg_labels
    call data_set_nstd_ccsm (                          &
                     tavg_file_desc, 'define',         &
                     nstd_field_id,                    &
                     avail_tavg_labels(n)%ndims,       &
                     io_dims_labels(:,n),              &
            short_name=avail_tavg_labels(n)%short_name,&
             long_name=avail_tavg_labels(n)%long_name, &
                 units=avail_tavg_labels(n)%units,     &
                nftype=avail_tavg_labels(n)%nftype     )

    avail_tavg_labels_id(n) = nstd_field_id

   enddo


!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_define_labels_ccsm
 

!***********************************************************************
!BOP
! !IROUTINE: tavg_construct_ccsm_time
! !INTERFACE:

 subroutine tavg_construct_ccsm_time (field_id,ns)

! !DESCRIPTION:
!  This routine defines the netCDF time variables that are used in the
!  ccsm netCDF output tavg files
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:
   integer (int_kind), intent(in) ::  &
      field_id,        &! time id
      ns                ! tavg streams index

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real (r8), dimension(1) ::  &
      TIME1D  

   character(char_len) ::  &
      nftype              

   integer (int_kind) ::  &
      ndims

   save  

   !*** time
   TIME1D(1)=tday00
   time_coordinate(1,ns) = construct_io_field('time',time_dim,        &
                     field_id  = field_id,                            &
                     long_name ='time',                               &
                     units     ='days since 0000-01-01 00:00:00',     &
                     d1d_array = TIME1D                               )

   if (field_id == 0) then
     call add_attrib_io_field(time_coordinate(1,ns), 'bounds', 'time_bound')
     call add_attrib_io_field(time_coordinate(1,ns), 'calendar', 'noleap')
   endif


!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_construct_ccsm_time
 
!***********************************************************************
!BOP
! !IROUTINE: tavg_define_time_bounds
! !INTERFACE:

 subroutine tavg_define_time_bounds(tavg_file_desc)

! !DESCRIPTION:
!  This routine defines the netCDF time variables that are used in the
!  ccsm netCDF output tavg files
!
! !REVISION HISTORY:
!  same as module


! !INPUT/OUTPUT PARAMETERS:
   type (datafile), intent(inout) :: tavg_file_desc    ! IO file descriptor

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real (r8), dimension(1) ::  &
      TIME1D  

   character(char_len) ::  &
      nftype              

   integer (int_kind) ::  &
      ndims

   save  


   !*** time_bound
   time_bound_dims(1) = d2_dim
   time_bound_dims(2) = time_dim
   ndims = 2
   nftype = 'double'

   call data_set_nstd_ccsm (tavg_file_desc, 'define',               &
                            time_bound_id,                          & 
                            ndims, time_bound_dims,                 &
                short_name='time_bound',                            &
                 long_name='boundaries for time-averaging interval',&
                     units='days since 0000-01-01 00:00:00',        &
                    nftype=nftype                                   )

!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_define_time_bounds
 

!***********************************************************************
!BOP
! !IROUTINE: tavg_write_vars_ccsm
! !INTERFACE:

 subroutine tavg_write_vars_ccsm (tavg_file_desc, nvars, ccsm_vars) 

! !DESCRIPTION:
!  This routine writes the ccsm variables (coordinates, scalars, and 
!  time-independent variables) to the ccsm version of the netCDF 
!  tavg output files
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:
   integer (int_kind), intent(in) ::  &
      nvars

! !INPUT/OUTPUT PARAMETERS:

   type (datafile), intent(inout) ::  &
      tavg_file_desc    ! IO file descriptor

   type (io_field_desc), dimension(nvars), intent(inout)  :: &
      ccsm_vars        

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::  &
      n, ndims             ! local index

   do n=1,nvars
      call data_set (tavg_file_desc, 'write', ccsm_vars(n))
   enddo

!-----------------------------------------------------------------------
!EOC
 end subroutine tavg_write_vars_ccsm


!***********************************************************************
!BOP
! !IROUTINE: tavg_write_vars_nstd_ccsm
! !INTERFACE:

 subroutine tavg_write_vars_nstd_ccsm(tavg_file_desc,ns) 

! !DESCRIPTION:
!  This routine writes the nonstandard ccsm variables to the
!  ccsm netCDF output tavg file.  These variables exist in
!  a variable of type (tavg_field_desc_ccsm), but do not have the
!  corresponding data bundled together via the function
!  construct_io_field. Therefore, this routine associates the
!  variable's data with the varible's definition, and sends
!  the information in separate pieces to be written.
!  Calls data_set_nstd_ccsm
!  
! 
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:
   integer (int_kind), intent(in) ::  &
      ns                ! tavg streams index

! !INPUT/OUTPUT PARAMETERS:
   type (datafile), intent(inout) :: tavg_file_desc    ! IO file descriptor

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   character (char_len) ::  &
      nftype

   integer (int_kind), parameter ::  &
      max_writes = 1

   integer (int_kind) ::  &
      n, ndims, nnn,      &! local index
      ii, indx

   integer (i4) ::  &
      id

   real (r8), dimension(2,max_writes) ::  &
      data_2d_r8

   character (char_len), dimension(30) :: &
      data_1d_ch

   type(io_dim), dimension(2) ::  & ! local array
       io_dims




   !*** moc-related variables
   if (ltavg_moc_diags(ns)) then
 
      !*** determine index of moc_components

      indx=0
      do ii=1,num_avail_tavg_labels
      if (trim(avail_tavg_labels(ii)%short_name) == 'moc_components')   &
          indx = ii
      enddo

     
      if (indx /= 0) then
        ndims=2
        data_1d_ch(1) = 'Eulerian Mean'
        if (n_moc_comp >= 2) &
          data_1d_ch(2) = 'Eddy-Induced (bolus)'
        if (n_moc_comp >= 3) & 
          data_1d_ch(3) = 'Submeso'
        io_dims(:) = io_dims_labels(:,indx)
        id = avail_tavg_labels_id(indx)
        nftype = 'char'
        lactive_time_dim = .false.
 

      call data_set_nstd_ccsm (tavg_file_desc, 'write',    &
                               id, ndims, io_dims, nftype, &
                               data_1d_ch = data_1d_ch     )

       endif

      !*** determine index of transport_regions

      indx=0
      do ii=1,num_avail_tavg_labels
      if (trim(avail_tavg_labels(ii)%short_name) == 'transport_regions')   &
          indx = ii
      enddo

     
      if (indx /= 0) then
        ndims=2
        data_1d_ch(1) = 'Global Ocean - Marginal Seas'
        if (n_transport_reg > 1 .and. nreg2_transport >=  1) then
          data_1d_ch(2) = trim(transport_region_info(1)%name)
          do nnn=2,nreg2_transport
             data_1d_ch(2)=trim(data_1d_ch(2)) /&
                                                 &/ ' + ' /&
                                                 &/trim(transport_region_info(nnn)%name)
          enddo
        endif
        io_dims(:) = io_dims_labels(:,indx)
        id = avail_tavg_labels_id(indx)
        nftype = 'char'
        lactive_time_dim = .false.
 
         call data_set_nstd_ccsm (tavg_file_desc, 'write',    &
                                  id, ndims, io_dims, nftype, &
                                  data_1d_ch = data_1d_ch     )
       endif

       lactive_time_dim = .true.
       call write_nstd_netcdf(tavg_file_desc,moc_id,5,  &
           io_dims_nstd_ccsm(:,1),'float',              &  ! <-- generalize later
           lactive_time_dim,                            &
           indata_4d_r4=TAVG_MOC_G)

   endif

   !*** transport variables
   if (ltavg_n_heat_trans(ns) .or. ltavg_n_salt_trans(ns)) then
 
      !*** determine index of transport_components

      indx=0
      do ii=1,num_avail_tavg_labels
      if (trim(avail_tavg_labels(ii)%short_name) == 'transport_components')   &
          indx = ii
      enddo

     
      if (indx /= 0) then
        ndims=2
        data_1d_ch(1) = 'Total'
        data_1d_ch(2) = 'Eulerian-Mean Advection'
        if (registry_match('init_gm')) then 
          data_1d_ch(3) = 'Eddy-Induced Advection (bolus) + Diffusion'
        else
          data_1d_ch(3) = 'Diffusion'
        endif
        if ( n_transport_comp >= 4 ) & 
          data_1d_ch(4) = 'Eddy-Induced (bolus) Advection'
        if ( n_transport_comp >= 5 ) & 
          data_1d_ch(5) = 'Submeso Advection'
        io_dims(:) = io_dims_labels(:,indx)
        id = avail_tavg_labels_id(indx)
        nftype = 'char'
        lactive_time_dim = .false.
 

      call data_set_nstd_ccsm (tavg_file_desc, 'write',    &
                               id, ndims, io_dims, nftype, &
                               data_1d_ch = data_1d_ch     )
       endif


   endif

   if (ltavg_n_heat_trans(ns)) then
       !*** N_HEAT
       lactive_time_dim = .true.
       call write_nstd_netcdf(                              &
          tavg_file_desc,n_heat_id,4,io_dims_nstd_ccsm(:,2),& !<-- generalize later
          'float',                                          &  ! <-- generalize later
           lactive_time_dim,                                &
          indata_3d_r4=TAVG_N_HEAT_TRANS_G)
   endif

   if (ltavg_n_salt_trans(ns)) then
       !*** N_SALT
       lactive_time_dim = .true.
       call write_nstd_netcdf(                             &
          tavg_file_desc,n_salt_id,4,io_dims_nstd_ccsm(:,3),& !<-- generalize later
          'float',                                         &  ! <-- generalize later
           lactive_time_dim,                               &
          indata_3d_r4=TAVG_N_SALT_TRANS_G)
   endif


!-----------------------------------------------------------------------
!EOC
 
 end subroutine tavg_write_vars_nstd_ccsm

!***********************************************************************
!BOP
! !IROUTINE: tavg_write_time_bounds
! !INTERFACE:

 subroutine tavg_write_time_bounds(tavg_file_desc,ns) 

! !DESCRIPTION:
!  This routine writes the ccsm time_bounds values to the
!  ccsm netCDF output tavg file. 
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: ns                ! tavg streams index

! !INPUT/OUTPUT PARAMETERS:

   type (datafile), intent(inout) :: tavg_file_desc    ! IO file descriptor

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   character (char_len) ::  &
      nftype

   integer (int_kind) ::  &
      ndims

   real (r8), dimension(2,1) ::  &   ! (2nd dimension is the unlimited dimension)
      data_2d_r8


   tavg_streams(ns)%upper_time_bound = tday00
   ndims=2
   data_2d_r8(1,1) = tavg_streams(ns)%lower_time_bound
   data_2d_r8(2,1) = tavg_streams(ns)%upper_time_bound

   time_bound_dims(2)%start = time_dim%start
   time_bound_dims(2)%stop  = time_dim%stop 

   call write_time_bounds (tavg_file_desc,time_bound_id,time_bound_dims, data_2d_r8 )

   tavg_streams(ns)%lower_time_bound = tday00 


!-----------------------------------------------------------------------
!EOC
 
 end subroutine tavg_write_time_bounds


!***********************************************************************
!BOP
! !IROUTINE:  tavg_mask
! !INTERFACE:

 subroutine tavg_mask(ARRAY,tavg_field,k)

! !DESCRIPTION:
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (rtavg), dimension (:,:,:), intent(inout) ::  &
      ARRAY

   type (tavg_field_desc_ccsm), intent(in) ::  &
      tavg_field

   integer (int_kind), intent(in) ::  &
      k

   type (block) ::        &
      this_block          ! block information for current block

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   integer (int_kind) :: iblock, i, j


   select case (trim(tavg_field%grid_loc(2:3)))

     case('11')
          
       do iblock=1,nblocks_clinic
        ARRAY(:,:,iblock) = merge (tavg_field%fill_value,ARRAY(:,:,iblock),  &
                              k > KMT(:,:,iblock))
       enddo ! iblock
 
     case('22')

       do iblock=1,nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
             if (MASK_22(i,j,k,iblock)) ARRAY(i,j,iblock) = tavg_field%fill_value
         enddo ! i
         enddo ! j
       enddo ! iblock
 
     case default
         !*** if field_desc is not on standard tracer or 
         !    velocity grid, never mind

   end select


!EOC
 
 end subroutine tavg_mask


!***********************************************************************
 
!BOP
! !IROUTINE: tavg_bsf_diags
! !INTERFACE:

 subroutine tavg_bsf_diags(ns)

! !DESCRIPTION:
!
!  Compute barotropic stream function diagnostically from other tavg
!  quantities
!
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) ::  &
      ns


!EOP
!BOC

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      errorCode

   real (r8), dimension (nx_block, ny_block,max_blocks_clinic) ::  &
      WORK1,WORK2,WORK3,  &
      PSI_T, PSI_U

   integer (int_kind) ::   &
      tavg_id_SU,          &
      tavg_id_SV,          &
      tavg_id_BSF,         &
      tavg_loc_SU,         &
      tavg_loc_SV,         &
      tavg_loc_BSF

   integer (int_kind) ::   &
      iblock,              &
      i,ii,                &
      j,jj,                &
      nstream 

   type (block) ::        &
      this_block          ! block information for current block


   errorCode = POP_Success

!-----------------------------------------------------------------------
!
!  test return conditions
!
!-----------------------------------------------------------------------

   if (.not. ldiag_bsf) return
 
   if (tavg_freq_iopt(ns) == freq_opt_nstep) then
     if (nsteps_run <= 1 .and. my_task == master_task) then
        write(stdout,*)  &
        'WARNING: BSF diagnostic is not computed if tavg_freq_iopt == freq_opt_nstep'
        call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout) ! temporary
     endif
     return
   endif

!-----------------------------------------------------------------------
!
!  begin computations
!
!-----------------------------------------------------------------------

   !*** zero out location identifiers
   tavg_loc_SU = 0; tavg_loc_SV  = 0; tavg_loc_BSF = 0

   !*** determine the tavg ids for tavg fields required by this module
   tavg_id_SU  = tavg_id('SU')
   tavg_id_SV  = tavg_id('SV')
   tavg_id_BSF = tavg_id('BSF')

   if (tavg_id_SU == 0 .or. tavg_id_SV == 0 .or. tavg_id_BSF == 0) then
     !*** do not abort; write warning message and proceed
     call document ('write_tavg', 'WARNING: cannot compute BSF diagnostic')
   endif
 
   !*** only compute BSF diagnostic if SU,SV, and BSF fields are in this stream
   nstream = tavg_in_which_stream(tavg_id_SU)

   if (ns /= nstream) then   
     !*** write warning message?
     return
   endif

   if (nstream /= tavg_in_which_stream(tavg_id_SV) .or.  &
       nstream /= tavg_in_which_stream(tavg_id_BSF)) then
     !*** write warning message?
     return
   endif

   call timer_start(timer_tavg_ccsm_diags_bsf)

   tavg_loc_SU  = avail_tavg_fields(tavg_id_SU)%buf_loc
   tavg_loc_SV  = avail_tavg_fields(tavg_id_SV)%buf_loc
   tavg_loc_BSF = avail_tavg_fields(tavg_id_BSF)%buf_loc

  !$OMP PARALLEL DO PRIVATE (this_block,iblock)
   do iblock=1,nblocks_clinic
     this_block = get_block(blocks_clinic(iblock),iblock)

     WORK2(:,:,iblock) = TAVG_BUF_2D(:,:,iblock,tavg_loc_SU)
     WORK3(:,:,iblock) = TAVG_BUF_2D(:,:,iblock,tavg_loc_SV)
     call zcurl (1,WORK1(:,:,iblock),WORK2(:,:,iblock),  &
                   WORK3(:,:,iblock),this_block)
   end do
  !$OMP END PARALLEL DO
      
   WORK2 = c0
   call pcg_diag_bsf_solver (WORK2,WORK1, errorCode)

   if (errorCode /= POP_Success) then
      exit_string = 'FATAL ERROR  in pcg_diag_bsf'
      call POP_ErrorSet(errorCode, exit_string)
      call document ('tavg_bsf_diags', exit_string)
      call exit_POP (sigAbort,exit_string,out_unit=stdout)
      return
   endif
 
   do iblock=1,nblocks_clinic
   TAVG_BUF_2D(:,:,iblock,tavg_loc_BSF) =  &
       merge(c0,WORK2(:,:,iblock), .not. CALCT(:,:,iblock))
   enddo ! iblock

   !*** convert to Sv
   TAVG_BUF_2D(:,:,:,tavg_loc_BSF) =  &
     TAVG_BUF_2D(:,:,:,tavg_loc_BSF)*mass_to_Sv

   do iblock=1,nblocks_clinic
     PSI_T(:,:,iblock)= TAVG_BUF_2D(:,:,iblock,tavg_loc_BSF)
   enddo

! !$OMP PARALLEL DO PRIVATE (iblock, i,j,ii,jj)  ! cannot exit inside a threaded loop
   do iblock=1,nblocks_clinic

      PSI_T(:,:,iblock) = TAVG_BUF_2D(:,:,iblock,tavg_loc_BSF)
      do j=2,ny_block-1
      do i=2,nx_block-1
        if (KMT(i,j,iblock) == 0) then
          ii_loop: do ii=i-1,i+1
          jj_loop: do jj=j-1,j+1
            if (KMT(ii,jj,iblock) /= 0) then
              PSI_T(i,j,iblock) = PSI_T(ii,jj,iblock)
              exit ii_loop
            endif
          enddo jj_loop
          enddo ii_loop
         endif
      enddo ! i
      enddo ! j
 
      call tgrid_to_ugrid (PSI_U(:,:,iblock), PSI_T(:,:,iblock),iblock)
 
      TAVG_BUF_2D(:,:,iblock,tavg_loc_BSF) = PSI_U(:,:,iblock)

   end do
! !$OMP END PARALLEL DO

 
   !*** stop bsf timer
   call timer_stop(timer_tavg_ccsm_diags_bsf)

!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_bsf_diags

!***********************************************************************
!BOP
! !IROUTINE: tavg_init_moc_diags
! !INTERFACE:

 subroutine tavg_init_moc_diags

! !DESCRIPTION:
!
!  Initialization routine for tavg_moc_diags, which computes meridional overturning 
!  circulation diagnostically from  other tavg quantities. This routine mainly
!  performs error checking: all the required tavg fields have been requested,
!  have been activated, and are in the same stream
!
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables:
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      moc_stream    ! stream in which WVEL is defined

!-----------------------------------------------------------------------
!
!  set default value regardless of moc_requested status
!   
!-----------------------------------------------------------------------

   ltavg_moc_diags = .false.

!-----------------------------------------------------------------------
!
!  is MOC diagnostic requested?
!
!-----------------------------------------------------------------------

   if (.not. moc_requested) return

!-----------------------------------------------------------------------
!
!  begin initialization
!
!-----------------------------------------------------------------------

   tavg_loc_WVEL      = 0 ; tavg_loc_VVEL      = 0
   tavg_loc_WISOP     = 0 ; tavg_loc_VISOP     = 0
   tavg_loc_WSUBM     = 0 ; tavg_loc_VSUBM     = 0
 
   tavg_id_WVEL   = tavg_id('WVEL')
   tavg_id_VVEL   = tavg_id('VVEL')

   if (tavg_id_WVEL  == 0 .or. tavg_id_VVEL  == 0 ) then
     call document ('tavg_moc', &
         'Error in moc diagnostics computations: none of the following should be zero')
     call document ('tavg_moc', 'tavg_id_WVEL', tavg_id_WVEL)
     call document ('tavg_moc', 'tavg_id_VVEL', tavg_id_VVEL)
     exit_string = 'FATAL ERROR  in tavg_moc'
     call document ('tavg_moc', exit_string)
     call exit_POP (sigAbort,exit_string,out_unit=stdout)
   endif
 
   !*** may not yet be activated, so use abs
   tavg_loc_WVEL = abs(avail_tavg_fields(tavg_id_WVEL)%buf_loc) 
   tavg_loc_VVEL = abs(avail_tavg_fields(tavg_id_VVEL)%buf_loc)

   if (tavg_loc_WVEL  == 0 .or. tavg_loc_VVEL  == 0 ) then
     call document ('tavg_moc', &
      'Error in moc diagnostics computations: none of the following should be zero')
     call document ('tavg_moc', 'tavg_loc_WVEL', tavg_loc_WVEL)
     call document ('tavg_moc', 'tavg_loc_VVEL', tavg_loc_VVEL)
     exit_string = 'FATAL ERROR  in tavg_moc'
     call document ('tavg_moc', exit_string)
     call exit_POP (sigAbort,exit_string,out_unit=stdout)
   endif

!-----------------------------------------------------------------------
!
!  checks for gm_bolus option terms
!
!-----------------------------------------------------------------------

   if (ldiag_gm_bolus) then
 
     tavg_id_WISOP  = tavg_id('WISOP')
     tavg_id_VISOP  = tavg_id('VISOP')
 
     if (tavg_id_WISOP  == 0 .or. tavg_id_VISOP  == 0 ) then
       call document ('tavg_moc', &
           'Error in moc diagnostics computations: none of the following should be zero')
       call document ('tavg_moc', 'tavg_id_WISOP',  tavg_id_WISOP)
       call document ('tavg_moc', 'tavg_id_VISOP',  tavg_id_VISOP)
       exit_string = 'FATAL ERROR  in tavg_moc'
       call document ('tavg_moc', exit_string)
       call exit_POP (sigAbort,exit_string,out_unit=stdout)
     endif

     !*** may not yet be activated, so use abs
     tavg_loc_WISOP = abs(avail_tavg_fields(tavg_id_WISOP)%buf_loc) 
     tavg_loc_VISOP = abs(avail_tavg_fields(tavg_id_VISOP)%buf_loc)

     if (tavg_loc_WISOP  == 0 .or. tavg_loc_VISOP  == 0 ) then
       call document ('tavg_moc', &
         'Error in moc diagnostics computations: none of the following should be zero')
       call document ('tavg_moc', 'tavg_loc_WISOP',  tavg_loc_WISOP)
       call document ('tavg_moc', 'tavg_loc_VISOP',  tavg_loc_VISOP)
       exit_string = 'FATAL ERROR  in tavg_moc'
       call document ('tavg_moc', exit_string)
       call exit_POP (sigAbort,exit_string,out_unit=stdout)
     endif

   endif ! ldiag_gm_bolus

!-----------------------------------------------------------------------
!
!  checks for submeso option terms
!
!-----------------------------------------------------------------------

   if (lsubmeso) then

     tavg_id_WSUBM  = tavg_id('WSUBM')
     tavg_id_VSUBM  = tavg_id('VSUBM')

     if (tavg_id_WSUBM  == 0 .or. tavg_id_VSUBM  == 0 ) then
       call document ('tavg_moc', &
         'Error in moc diagnostics computations: none of the following should be zero')
       call document ('tavg_moc', 'tavg_id_WSUBM',  tavg_id_WSUBM)
       call document ('tavg_moc', 'tavg_id_VSUBM',  tavg_id_VSUBM)
       exit_string = 'FATAL ERROR  in tavg_moc'
       call document ('tavg_moc', exit_string)
       call exit_POP (sigAbort,exit_string,out_unit=stdout)
     endif

     !*** may not yet be activated, so use abs
     tavg_loc_WSUBM = abs(avail_tavg_fields(tavg_id_WSUBM)%buf_loc)
     tavg_loc_VSUBM = abs(avail_tavg_fields(tavg_id_VSUBM)%buf_loc)

     if (tavg_loc_WSUBM  == 0 .or. tavg_loc_VSUBM  == 0 ) then
       call document ('tavg_moc', &
         'Error in moc diagnostics computations: none of the following should be zero')
       call document ('tavg_moc', 'tavg_loc_WSUBM',  tavg_loc_WSUBM)
       call document ('tavg_moc', 'tavg_loc_VSUBM',  tavg_loc_VSUBM)
       exit_string = 'FATAL ERROR  in tavg_moc'
       call document ('tavg_moc', exit_string)
       call exit_POP (sigAbort,exit_string,out_unit=stdout)
     endif

   endif ! lsubmeso

!-----------------------------------------------------------------------
!
!  Now confirm that all required fields are in the same stream
!
!-----------------------------------------------------------------------

   moc_stream =  tavg_in_which_stream(tavg_id_WVEL)


   if (lsubmeso) then
     if (tavg_in_this_stream(tavg_id_VVEL , moc_stream)  .and.   &
         tavg_in_this_stream(tavg_id_WISOP, moc_stream)  .and.   &
         tavg_in_this_stream(tavg_id_VISOP, moc_stream)  .and.   &
         tavg_in_this_stream(tavg_id_WSUBM, moc_stream)  .and.   &
         tavg_in_this_stream(tavg_id_VSUBM, moc_stream)          ) then
         ltavg_moc_diags(moc_stream) = .true.
     endif
   else if (ldiag_gm_bolus) then
     if (tavg_in_this_stream(tavg_id_VVEL , moc_stream)  .and.   &
         tavg_in_this_stream(tavg_id_WISOP, moc_stream)  .and.   &
         tavg_in_this_stream(tavg_id_VISOP, moc_stream)          ) then
         ltavg_moc_diags(moc_stream) = .true.
     endif
   else 
     if (tavg_in_this_stream(tavg_id_VVEL , moc_stream)          ) then
         ltavg_moc_diags(moc_stream) = .true.
     endif
   endif
   
!-----------------------------------------------------------------------
!    define MOC tavg variable and dimensions
!    note that this call fills num_avail_tavg_nstd_fields
!-----------------------------------------------------------------------

   if (ltavg_moc_diags(moc_stream)) then
     call define_tavg_field(                             &
        tavg_MOC, 'MOC', 5,                              &
        long_name='Meridional Overturning Circulation',  &
        units='Sverdrups',                               &
        coordinates='lat_aux_grid moc_z moc_components transport_region time',&
        nstd_fields=avail_tavg_nstd_fields,              &
        num_nstd_fields=num_avail_tavg_nstd_fields,      &
        max_nstd_fields=max_avail_tavg_nstd_fields       )

    endif

   call POP_IOUnitsFlush(POP_stdout); call POP_IOUnitsFlush(stdout)

 
!-----------------------------------------------------------------------
!EOC

  end subroutine tavg_init_moc_diags
 
!***********************************************************************
!BOP
! !IROUTINE: tavg_moc_diags
! !INTERFACE:

 subroutine tavg_moc_diags(ns)

! !DESCRIPTION:
!
!  Compute meridional overturning circulation diagnostically from 
!  other tavg quantities
!
!
! !REVISION HISTORY:
!  same as module


! !INPUT PARAMETERS:

   integer (int_kind), intent(in) ::  &
      ns


!EOP
!BOC

!-----------------------------------------------------------------------
!
!  is MOC diagnostic requested?
!
!-----------------------------------------------------------------------

   if (.not. ltavg_moc_diags(ns)) return
 
   if (tavg_freq_iopt(ns) == freq_opt_nstep) then
     if (nsteps_run <= 1 .and. my_task == master_task) then
        write(stdout,*)  &
        'WARNING: MOC diagnostic is not computed if tavg_freq_iopt == freq_opt_nstep'
        call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
     endif
     return
   endif

!-----------------------------------------------------------------------
!
!  begin computations
!
!-----------------------------------------------------------------------

   call timer_start(timer_tavg_ccsm_diags_moc)

   io_dims_nstd_ccsm(1,tavg_MOC) = lat_aux_grid_dim
   io_dims_nstd_ccsm(2,tavg_MOC) = moc_z_dim
   io_dims_nstd_ccsm(3,tavg_MOC) = moc_comp_dim
   io_dims_nstd_ccsm(4,tavg_MOC) = transport_reg_dim
   io_dims_nstd_ccsm(5,tavg_MOC) = time_dim
   ndims_nstd_ccsm  (  tavg_MOC) = 5

   if ( ldiag_gm_bolus  .and.  lsubmeso ) then
     call compute_moc ( TAVG_BUF_3D(:,:,:,:,tavg_loc_WVEL ),  &
                        TAVG_BUF_3D(:,:,:,:,tavg_loc_VVEL ),  &
                  W_I = TAVG_BUF_3D(:,:,:,:,tavg_loc_WISOP),  &
                  V_I = TAVG_BUF_3D(:,:,:,:,tavg_loc_VISOP),  &
                 W_SM = TAVG_BUF_3D(:,:,:,:,tavg_loc_WSUBM),  &
                 V_SM = TAVG_BUF_3D(:,:,:,:,tavg_loc_VSUBM))
  elseif ( ldiag_gm_bolus  .and.  .not.lsubmeso ) then
     call compute_moc ( TAVG_BUF_3D(:,:,:,:,tavg_loc_WVEL ),  &
                        TAVG_BUF_3D(:,:,:,:,tavg_loc_VVEL ),  &
                  W_I = TAVG_BUF_3D(:,:,:,:,tavg_loc_WISOP),  &
                  V_I = TAVG_BUF_3D(:,:,:,:,tavg_loc_VISOP))
  elseif ( .not.ldiag_gm_bolus  .and.  lsubmeso ) then
     call compute_moc ( TAVG_BUF_3D(:,:,:,:,tavg_loc_WVEL ),  &
                        TAVG_BUF_3D(:,:,:,:,tavg_loc_VVEL ),  &
                 W_SM = TAVG_BUF_3D(:,:,:,:,tavg_loc_WSUBM),  &
                 V_SM = TAVG_BUF_3D(:,:,:,:,tavg_loc_VSUBM))
  else
     call compute_moc ( TAVG_BUF_3D(:,:,:,:,tavg_loc_WVEL ),  &
                        TAVG_BUF_3D(:,:,:,:,tavg_loc_VVEL ))
  endif

  call timer_stop(timer_tavg_ccsm_diags_moc)
 
!-----------------------------------------------------------------------
!EOC

  end subroutine tavg_moc_diags

!***********************************************************************
!BOP
! !IROUTINE: tavg_init_transport_diags
! !INTERFACE:

 subroutine tavg_init_transport_diags

! !DESCRIPTION:
!
!  Compute northward transport of heat/salt diagnostically from 
!  other tavg quantities
!
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables:
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     trans_stream    ! stream in which ADVT is defined

!-----------------------------------------------------------------------
!
!  set default values regardless of n_heat_trans_requested/
!      n_salt_trans_requested status
!-----------------------------------------------------------------------

   ltavg_n_heat_trans = .false.
   ltavg_n_salt_trans = .false.

!-----------------------------------------------------------------------
!
!  are transport diagnostics requested?
!
!-----------------------------------------------------------------------

   if (.not. (n_heat_trans_requested .or. n_salt_trans_requested)) return

!-----------------------------------------------------------------------
!
!  begin initialization
!
!-----------------------------------------------------------------------
   tavg_loc_ADVT      = 0 ; tavg_loc_ADVS      = 0
   tavg_loc_VNT       = 0 ; tavg_loc_VNS       = 0
   tavg_loc_HDIFT     = 0 ; tavg_loc_HDIFS     = 0
   tavg_loc_ADVT_ISOP = 0 ; tavg_loc_ADVS_ISOP = 0
   tavg_loc_VNT_ISOP  = 0 ; tavg_loc_VNS_ISOP  = 0
   tavg_loc_ADVT_SUBM = 0 ; tavg_loc_ADVS_SUBM = 0
   tavg_loc_VNT_SUBM  = 0 ; tavg_loc_VNS_SUBM  = 0

   tavg_id_ADVT      = tavg_id('ADVT')
   tavg_id_ADVS      = tavg_id('ADVS')
   tavg_id_VNT       = tavg_id('VNT')
   tavg_id_VNS       = tavg_id('VNS')
   tavg_id_HDIFT     = tavg_id('HDIFT')
   tavg_id_HDIFS     = tavg_id('HDIFS')


   !*** may not yet be activated, so use abs
   tavg_loc_ADVT   = abs(avail_tavg_fields(tavg_id_ADVT )%buf_loc)
   tavg_loc_ADVS   = abs(avail_tavg_fields(tavg_id_ADVS )%buf_loc)
   tavg_loc_VNT    = abs(avail_tavg_fields(tavg_id_VNT  )%buf_loc)
   tavg_loc_VNS    = abs(avail_tavg_fields(tavg_id_VNS  )%buf_loc)
   tavg_loc_HDIFT  = abs(avail_tavg_fields(tavg_id_HDIFT)%buf_loc)
   tavg_loc_HDIFS  = abs(avail_tavg_fields(tavg_id_HDIFS)%buf_loc)


   if (tavg_loc_ADVT  == 0 .or. tavg_loc_ADVS  == 0 .or.  &
       tavg_loc_VNT   == 0 .or. tavg_loc_VNS   == 0 .or.  &
       tavg_loc_HDIFT == 0 .or. tavg_loc_HDIFS == 0) then
       call document ('tavg_transport', &
         'Error in heat/salt transport diags: none of the following should be zero')
       call document ('tavg_transport', 'tavg_loc_ADVT',  tavg_loc_ADVT)
       call document ('tavg_transport', 'tavg_loc_ADVS',  tavg_loc_ADVS)
       call document ('tavg_transport', 'tavg_loc_VNT ',  tavg_loc_VNT )
       call document ('tavg_transport', 'tavg_loc_VNS ',  tavg_loc_VNS )
       call document ('tavg_transport', 'tavg_loc_HDIFT', tavg_loc_HDIFT)
       call document ('tavg_transport', 'tavg_loc_HDIFS', tavg_loc_HDIFS)
       call exit_POP (SigAbort, '(tavg_transport) ERROR',out_unit=stdout)
       exit_string = 'FATAL ERROR  in tavg_transport'
       call document ('tavg_transport', exit_string)
       call exit_POP (sigAbort,exit_string,out_unit=stdout)
   endif

!-----------------------------------------------------------------------
!
!  checks for gm_bolus option terms
!
!-----------------------------------------------------------------------

   if (ldiag_gm_bolus) then
       tavg_id_ADVT_ISOP = tavg_id('ADVT_ISOP')
       tavg_id_ADVS_ISOP = tavg_id('ADVS_ISOP')
       tavg_id_VNT_ISOP  = tavg_id('VNT_ISOP')
       tavg_id_VNS_ISOP  = tavg_id('VNS_ISOP')

       !*** note: error checking for tavg_id is done in POP_checks
       !*** may not yet be activated, so use abs
       tavg_loc_ADVT_ISOP = abs(avail_tavg_fields(tavg_id_ADVT_ISOP)%buf_loc)
       tavg_loc_ADVS_ISOP = abs(avail_tavg_fields(tavg_id_ADVS_ISOP)%buf_loc)
       tavg_loc_VNT_ISOP  = abs(avail_tavg_fields(tavg_id_VNT_ISOP )%buf_loc)
       tavg_loc_VNS_ISOP  = abs(avail_tavg_fields(tavg_id_VNS_ISOP )%buf_loc)

       if (tavg_loc_ADVT_ISOP  == 0 .or. tavg_loc_ADVS_ISOP  == 0 .or.  &
           tavg_loc_VNT_ISOP   == 0 .or. tavg_loc_VNS_ISOP   == 0 ) then
          call document ('tavg_transport', &
            'Error in heat/salt transport diags: none of the following should be zero')
          call document ('tavg_transport', 'tavg_loc_ADVT_ISOP',  tavg_loc_ADVT)
          call document ('tavg_transport', 'tavg_loc_ADVS_ISOP',  tavg_loc_ADVS)
          call document ('tavg_transport', 'tavg_loc_VNT_ISOP ',  tavg_loc_VNT )
          call document ('tavg_transport', 'tavg_loc_VNS_ISOP ',  tavg_loc_VNS )
          exit_string = 'FATAL ERROR  in tavg_transport'
          call document ('tavg_transport', exit_string)
          call exit_POP (sigAbort,exit_string,out_unit=stdout)
       endif ! tavg_id testing
   endif ! ldiag_gm_bolus

!-----------------------------------------------------------------------
!
!  checks for submeso option terms
!
!-----------------------------------------------------------------------

   if (lsubmeso) then
       tavg_id_ADVT_SUBM = tavg_id('ADVT_SUBM')
       tavg_id_ADVS_SUBM = tavg_id('ADVS_SUBM')
       tavg_id_VNT_SUBM  = tavg_id('VNT_SUBM')
       tavg_id_VNS_SUBM  = tavg_id('VNS_SUBM')

       !*** may not yet be activated, so use abs
       tavg_loc_ADVT_SUBM = abs(avail_tavg_fields(tavg_id_ADVT_SUBM)%buf_loc)
       tavg_loc_ADVS_SUBM = abs(avail_tavg_fields(tavg_id_ADVS_SUBM)%buf_loc)
       tavg_loc_VNT_SUBM  = abs(avail_tavg_fields(tavg_id_VNT_SUBM )%buf_loc)
       tavg_loc_VNS_SUBM  = abs(avail_tavg_fields(tavg_id_VNS_SUBM )%buf_loc)

       if (tavg_loc_ADVT_SUBM  == 0 .or. tavg_loc_ADVS_SUBM  == 0 .or.  &
           tavg_loc_VNT_SUBM   == 0 .or. tavg_loc_VNS_SUBM   == 0 ) then
          call document ('tavg_transport', &
            'Error in heat/salt transport diags: none of the following should be zero')
          call document ('tavg_transport', 'tavg_loc_ADVT_SUBM',  tavg_loc_ADVT_SUBM)
          call document ('tavg_transport', 'tavg_loc_ADVS_SUBM',  tavg_loc_ADVS_SUBM)
          call document ('tavg_transport', 'tavg_loc_VNT_SUBM ',  tavg_loc_VNT_SUBM )
          call document ('tavg_transport', 'tavg_loc_VNS_SUBM ',  tavg_loc_VNS_SUBM )
          exit_string = 'FATAL ERROR  in tavg_transport'
          call document ('tavg_transport', exit_string)
          call exit_POP (sigAbort,exit_string,out_unit=stdout)
       endif ! tavg_id testing
   endif ! lsubmeso

!-----------------------------------------------------------------------
!
!  now confirm that all required fields are in the same stream
!
!-----------------------------------------------------------------------

   trans_stream =  tavg_in_which_stream(tavg_id_ADVT)

   if (lsubmeso) then
     if (tavg_in_this_stream(tavg_id_ADVS,     trans_stream)  .and.   &
         tavg_in_this_stream(tavg_id_VNT,      trans_stream)  .and.   &
         tavg_in_this_stream(tavg_id_VNS,      trans_stream)  .and.   &
         tavg_in_this_stream(tavg_id_ADVT_ISOP,trans_stream)  .and.   &
         tavg_in_this_stream(tavg_id_ADVS_ISOP,trans_stream)  .and.   &
         tavg_in_this_stream(tavg_id_VNT_ISOP, trans_stream)  .and.   &
         tavg_in_this_stream(tavg_id_VNS_ISOP, trans_stream)  .and.   &
         tavg_in_this_stream(tavg_id_ADVT_SUBM,trans_stream)  .and.   &
         tavg_in_this_stream(tavg_id_ADVS_SUBM,trans_stream)  .and.   &
         tavg_in_this_stream(tavg_id_VNT_SUBM, trans_stream)  .and.   &
         tavg_in_this_stream(tavg_id_VNS_SUBM, trans_stream)  .and.   &
         tavg_in_this_stream(tavg_id_HDIFT,    trans_stream)  .and.   &
         tavg_in_this_stream(tavg_id_HDIFS,    trans_stream)          ) then
         if (n_heat_trans_requested) ltavg_n_heat_trans(trans_stream) = .true.
         if (n_salt_trans_requested) ltavg_n_salt_trans(trans_stream) = .true.
     endif
   else if (ldiag_gm_bolus) then
     if (tavg_in_this_stream(tavg_id_ADVS,     trans_stream)  .and.   &
         tavg_in_this_stream(tavg_id_VNT,      trans_stream)  .and.   &
         tavg_in_this_stream(tavg_id_VNS,      trans_stream)  .and.   &
         tavg_in_this_stream(tavg_id_ADVT_ISOP,trans_stream)  .and.   &
         tavg_in_this_stream(tavg_id_ADVS_ISOP,trans_stream)  .and.   &
         tavg_in_this_stream(tavg_id_VNT_ISOP, trans_stream)  .and.   &
         tavg_in_this_stream(tavg_id_VNS_ISOP, trans_stream)  .and.   &
         tavg_in_this_stream(tavg_id_HDIFT,    trans_stream)  .and.   &
         tavg_in_this_stream(tavg_id_HDIFS,    trans_stream)          ) then
         if (n_heat_trans_requested) ltavg_n_heat_trans(trans_stream) = .true.
         if (n_salt_trans_requested) ltavg_n_salt_trans(trans_stream) = .true.
     endif
   else 
     if (tavg_in_this_stream(tavg_id_ADVS,     trans_stream)  .and.   &
         tavg_in_this_stream(tavg_id_VNT,      trans_stream)  .and.   &
         tavg_in_this_stream(tavg_id_VNS,      trans_stream)  .and.   &
         tavg_in_this_stream(tavg_id_HDIFT,    trans_stream)  .and.   &
         tavg_in_this_stream(tavg_id_HDIFS,    trans_stream)          ) then
         if (n_heat_trans_requested) ltavg_n_heat_trans(trans_stream) = .true.
         if (n_salt_trans_requested) ltavg_n_salt_trans(trans_stream) = .true.
     endif
   endif
 
!-----------------------------------------------------------------------
!
!  define heat and salt transport diagnostics fields
!    (note: these calls increment avail_tavg_nstd_fields)
!
!-----------------------------------------------------------------------
    if (ltavg_n_heat_trans(trans_stream)) then
      call define_tavg_field     (                  &
        tavg_N_HEAT, 'N_HEAT', 4,                   &
        long_name='Northward Heat Transport',       &
        units='Pwatt',                              &
        coordinates='lat_aux_grid transport_components transport_regions time',&
        nftype='float',                             &
        nstd_fields=avail_tavg_nstd_fields,         &
        num_nstd_fields=num_avail_tavg_nstd_fields, &
        max_nstd_fields=max_avail_tavg_nstd_fields  )
    endif

    if (ltavg_n_salt_trans(trans_stream)) then
     call define_tavg_field     (                       &
        tavg_N_SALT, 'N_SALT', 4,                       &
        long_name='Northward Salt Transport',           &
        units='gram centimeter^3/kg/s',                 &
        coordinates='lat_aux_grid transport_components transport_regions time',&
        nftype='float',                                 &
        nstd_fields=avail_tavg_nstd_fields,             &
        num_nstd_fields=num_avail_tavg_nstd_fields,     &
        max_nstd_fields=max_avail_tavg_nstd_fields  )
    endif

   call POP_IOUnitsFlush(POP_stdout); call POP_IOUnitsFlush(stdout)
!-----------------------------------------------------------------------
!EOC

  end subroutine tavg_init_transport_diags


!***********************************************************************
!BOP
! !IROUTINE: tavg_transport_diags
! !INTERFACE:

 subroutine tavg_transport_diags(ns)

! !DESCRIPTION:
!
!  Compute northward transport of heat/salt diagnostically from 
!  other tavg quantities
!
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) ::  &
      ns


!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables:
!
!-----------------------------------------------------------------------

   integer (int_kind) ::    &
      indx1,                &! index
      indx2,                &! index
      indx3,                &! index
      indx4,                &! index
      indx5,                &! index
      indx6,                &
      indx7,                &
      n

!-----------------------------------------------------------------------
!
!  are transport diagnostics requested?
!  is frequency compatible?
!-----------------------------------------------------------------------

   if (.not. (ltavg_n_heat_trans(ns) .or. ltavg_n_salt_trans(ns))) return

   if (tavg_freq_iopt(ns) == freq_opt_nstep) then
     if (nsteps_run <= 1 .and. my_task == master_task) then
        write(stdout,*)  &
        'WARNING: transport diagnostics are not computed if tavg_freq_iopt == freq_opt_nstep'
       call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
     endif
     return
   endif

!-----------------------------------------------------------------------
!
!  begin computations
!
!-----------------------------------------------------------------------

   call timer_start(timer_tavg_ccsm_diags_trans)

   io_dims_nstd_ccsm(1,tavg_N_HEAT) = lat_aux_grid_dim
   io_dims_nstd_ccsm(2,tavg_N_HEAT) = transport_comp_dim
   io_dims_nstd_ccsm(3,tavg_N_HEAT) = transport_reg_dim
   io_dims_nstd_ccsm(4,tavg_N_HEAT) = time_dim
   ndims_nstd_ccsm  (  tavg_N_HEAT) = 4

   io_dims_nstd_ccsm(1,tavg_N_SALT) = lat_aux_grid_dim
   io_dims_nstd_ccsm(2,tavg_N_SALT) = transport_comp_dim
   io_dims_nstd_ccsm(3,tavg_N_SALT) = transport_reg_dim
   io_dims_nstd_ccsm(4,tavg_N_SALT) = time_dim
   ndims_nstd_ccsm  (  tavg_N_SALT) = 4
 

   do n=1,2

     select case (n)
       case(1)
         indx1 = tavg_loc_ADVT
         indx2 = tavg_loc_HDIFT
         indx3 = tavg_loc_VNT
         if (ldiag_gm_bolus) then
           indx4 = tavg_loc_ADVT_ISOP
           indx5 = tavg_loc_VNT_ISOP
         endif
         if ( lsubmeso ) then
           indx6 = tavg_loc_ADVT_SUBM
           indx7 = tavg_loc_VNT_SUBM
         endif
       case(2)
         indx1 = tavg_loc_ADVS
         indx2 = tavg_loc_HDIFS
         indx3 = tavg_loc_VNS
         if (ldiag_gm_bolus) then
           indx4 = tavg_loc_ADVS_ISOP
           indx5 = tavg_loc_VNS_ISOP
         endif
         if ( lsubmeso ) then
           indx6 = tavg_loc_ADVS_SUBM
           indx7 = tavg_loc_VNS_SUBM
         endif
     end select
 
     if ( ldiag_gm_bolus  .and.  lsubmeso ) then
       call compute_tracer_transports (n,       &
                   TAVG_BUF_2D(:,:,:,  indx1),  &
                   TAVG_BUF_2D(:,:,:,  indx2),  &
                   TAVG_BUF_3D(:,:,:,:,indx3),  &
           ADV_I = TAVG_BUF_2D(:,:,:,  indx4),  &
           FN_I  = TAVG_BUF_3D(:,:,:,:,indx5),  &
          ADV_SM = TAVG_BUF_2D(:,:,:,  indx6),  &
           FN_SM = TAVG_BUF_3D(:,:,:,:,indx7))
     elseif ( ldiag_gm_bolus  .and.  .not.lsubmeso ) then
       call compute_tracer_transports (n,       &
                   TAVG_BUF_2D(:,:,:,  indx1),  &
                   TAVG_BUF_2D(:,:,:,  indx2),  &
                   TAVG_BUF_3D(:,:,:,:,indx3),  &
           ADV_I = TAVG_BUF_2D(:,:,:,  indx4),  &
           FN_I  = TAVG_BUF_3D(:,:,:,:,indx5))
     elseif ( .not.ldiag_gm_bolus  .and.  lsubmeso ) then
       call compute_tracer_transports (n,       &
                   TAVG_BUF_2D(:,:,:,  indx1),  &
                   TAVG_BUF_2D(:,:,:,  indx2),  &
                   TAVG_BUF_3D(:,:,:,:,indx3),  &
          ADV_SM = TAVG_BUF_2D(:,:,:,  indx6),  &
           FN_SM = TAVG_BUF_3D(:,:,:,:,indx7))
     else
       call compute_tracer_transports (n,       &
                   TAVG_BUF_2D(:,:,:,  indx1),  &
                   TAVG_BUF_2D(:,:,:,  indx2),  &
                   TAVG_BUF_3D(:,:,:,:,indx3))
     endif
   enddo
 
  !*** stop timer
   call timer_stop(timer_tavg_ccsm_diags_trans)
!-----------------------------------------------------------------------
!EOC

  end subroutine tavg_transport_diags


!***********************************************************************
!BOP
! !IROUTINE: tavg_init_local_spatial_avg
! !INTERFACE:

 subroutine tavg_init_local_spatial_avg

! !DESCRIPTION:
!
!  Initialize geographical masks and arrays for local spatial
!     averaging of some time-averaged fields
!
!       region key:
!        n_reg_0D = 1 --> Nino 1+2 region
!                 = 2 --> Nino  3  region
!                 = 3 --> Nino 3.4 region
!                 = 4 --> Nino  4  region
!
!       TLON_MIN_0D, TLON_MAX_0D, TLAT_MIN_0D, and TLAT_MAX_0D are all
!        in degrees. Also, TLON_MIN_0D and TLON_MAX_0D are in degrees
!        east.
!
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

!EOP
!BOC

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------
   type (block) ::        &
      this_block          ! block information for current block

   real (r8), dimension(n_reg_0D) ::  &
     TLON_MIN_0D = (/ 270.0_r8, 210.0_r8,      &
                      190.0_r8, 160.0_r8 /),   &
     TLON_MAX_0D = (/ 280.0_r8, 270.0_r8,      &
                      240.0_r8, 210.0_r8 /),   &
     TLAT_MIN_0D = (/ -10.0_r8,  -5.0_r8,      &
                       -5.0_r8,  -5.0_r8 /),   &
     TLAT_MAX_0D = (/   0.0_r8,   5.0_r8,      &
                        5.0_r8,   5.0_r8 /)     

   integer (int_kind) ::  &
      ns,                 &! local streams index
      iblock,             &! local block number
      n_reg,              &! loop index
      max_days             ! maximum number of days per month in a year





!-----------------------------------------------------------------------
!
!  set default value regardless of ltavg_nino_diags_requested status
! 
!-----------------------------------------------------------------------

   ltavg_nino_diags = .false.

!-----------------------------------------------------------------------
!
!  are NINO diagnostics requested?
!
!-----------------------------------------------------------------------

   if (.not. ltavg_nino_diags_requested) return

!-----------------------------------------------------------------------
!
!  begin initialization
!
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!    check if local spatial averaging is possible, based on tavg options
!    and namelist settings
!-----------------------------------------------------------------------

   !------------------------------------
   !    is TEMP available in any stream?
   !------------------------------------

   if (tavg_id('TEMP') == 0) then
     !*** do not abort; write warning message and proceed
     call document ('tavg_init_local_spatial_avg', &
       'WARNING: TEMP is not requested; cannot compute spatial averages from time-averaged TEMP')
     return
   endif

   max_days = maxval(days_in_month) 

   !------------------------------------
   !    is TEMP available in THIS stream?
   !------------------------------------

   do ns=1,nstreams

     if (tavg_in_this_stream(tavg_id('TEMP'),ns)) then
    
         !-------------------------------------------------------
         !    is tavg frequency consistent with nino diagnostics?
         !-------------------------------------------------------

      if ((tavg_freq_iopt(ns) == freq_opt_nmonth  .and. tavg_freq(ns) == 1)           .or. &
          (tavg_freq_iopt(ns) == freq_opt_nday    .and. tavg_freq(ns) <= max_days)    .or. &
          (tavg_freq_iopt(ns) == freq_opt_nhour   .and. tavg_freq(ns) <= max_days*24) .or. &
          (tavg_freq_iopt(ns) == freq_opt_nsecond .and. tavg_freq(ns) <= max_days*24*seconds_in_hour) .or. &
          (tavg_freq_iopt(ns) == freq_opt_nstep   .and. tavg_freq(ns)<= max_days*nsteps_per_day) ) then
          !*** ok to have ltavg_nino_diags enabled in this stream
          ltavg_nino_diags(ns) = .true.
      else
          !*** not ok to have ltavg_nino_diags enabled; disable it and proceed
         call document ('tavg_init_local_spatial_avg', &
                         'WARNING: cannot compute spatial averages from time-averaged TEMP')
         return
      endif
        
     endif ! tavg_in_this_stream

   enddo ! nstreams


!-----------------------------------------------------------------------
!
!     allocate and initialize arrays
!
!-----------------------------------------------------------------------

   allocate ( SAVG_0D(n_reg_0D),       &
              SAVG_0D_AREA(n_reg_0D),  &
              SAVG_0D_NAME(n_reg_0D),  &
              SAVG_0D_MASK(nx_block,ny_block,nblocks_clinic,n_reg_0D) )

   SAVG_0D_NAME = (/'NINO_1_PLUS_2 ',  &
                    'NINO_3        ',  &
                    'NINO_3_POINT_4',  &
                    'NINO_4        '/)

   SAVG_0D      = c0
   SAVG_0D_AREA = c0
   SAVG_0D_MASK = c0

   !*** determine masks for each region



   do iblock = 1,nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock) 
      do n_reg=1,n_reg_0D
        where (   CALCT(:,:,iblock)              .and.   &
                ( TLOND(:,:,iblock) >= TLON_MIN_0D(n_reg) )  .and.   &
                ( TLOND(:,:,iblock) <= TLON_MAX_0D(n_reg) )  .and.   &
                ( TLATD(:,:,iblock) >= TLAT_MIN_0D(n_reg) )  .and.   &
                ( TLATD(:,:,iblock) <= TLAT_MAX_0D(n_reg) ) )
          SAVG_0D_MASK(:,:,iblock,n_reg) = c1
        endwhere
      enddo ! n_reg
   enddo ! iblock

   !*** compute areas for each region to be used later for normalization
   do n_reg=1,n_reg_0D
     SAVG_0D_AREA(n_reg) = global_sum(TAREA(:,:,:),distrb_clinic,  &
                              field_loc_center,SAVG_0D_MASK(:,:,:,n_reg) )
     if ( SAVG_0D_AREA(n_reg) == c0 ) then
       exit_string = 'FATAL ERROR: SAVG_0D_AREA is zero.'
       call document ('tavg_init_local_spatial_avg', exit_string)
       call exit_POP (sigAbort,exit_string,out_unit=stdout)
     endif
   enddo

   call POP_IOUnitsFlush(POP_stdout); call POP_IOUnitsFlush(stdout)
!-----------------------------------------------------------------------
!EOC

  end subroutine tavg_init_local_spatial_avg

!***********************************************************************
!BOP
! !IROUTINE: tavg_local_spatial_avg
! !INTERFACE:

 subroutine tavg_local_spatial_avg(ns)

! !DESCRIPTION:
!
!  Compute local spatial averages from time-averaged TEMP
!
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) ::  &
      ns

!EOP
!BOC

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

   type (block) ::        &
      this_block          ! block information for current block

   character (char_len) ::  &
      time_string,          &
      date_string

   integer (int_kind) ::  &
      iblock,             &! local block number
      n_reg                ! loop index

   real (r8), dimension (:,:,:), allocatable ::  &
      WORK

!-----------------------------------------------------------------------
!
!  are NINO diagnotics requested?
!
!-----------------------------------------------------------------------

   if (.not. ltavg_nino_diags(ns)) return

!-----------------------------------------------------------------------
!
!  begin computations
!
!-----------------------------------------------------------------------

   allocate (WORK(nx_block,ny_block,nblocks_clinic))
   WORK = c0

   tavg_id_TEMP  = tavg_id('TEMP')

   if (tavg_id_TEMP == 0) then
     !*** do not abort; write warning message and proceed
     call document ('tavg_local_spatial_avg', &
       'WARNING: cannot compute spatial averages from time-averaged TEMP')
     return
   else
     tavg_loc_TEMP = avail_tavg_fields(tavg_id_TEMP)%buf_loc
   endif

   do iblock = 1,nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock) 
      WORK(:,:,iblock)=  &
          TAVG_BUF_3D(:,:,1,iblock,tavg_loc_TEMP)*TAREA(:,:,iblock)
   enddo ! iblock

   do n_reg=1,n_reg_0D
     SAVG_0D(n_reg) = global_sum ( WORK(:,:,:),distrb_clinic,  &
                         field_loc_center,SAVG_0D_MASK(:,:,:,n_reg))  &
                         / SAVG_0D_AREA(n_reg)
   enddo ! n_reg


   if ( my_task == master_task ) then

     call time_stamp('now','ymd',time_string=time_string)
     call time_stamp('now','ymd',date_string=date_string)

     write (stdout,*) ' '
     write (stdout,*)  &
     ' Local Time- and Space-Averages for Nino Regions: '  &
      // trim(time_string) // ' ' // trim(date_string)

     do n_reg=1,n_reg_0D
       write (stdout,*) trim(SAVG_0D_NAME(n_reg)),': ',  &
                        SAVG_0D(n_reg)
     enddo
     call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)

   endif
 
   deallocate (WORK)

!-----------------------------------------------------------------------
!EOC

  end subroutine tavg_local_spatial_avg
 
!***********************************************************************
!BOP
! !IROUTINE: final_tavg
! !INTERFACE:

 subroutine final_tavg

! !DESCRIPTION:
!  This routine closes all tavg files that are still open
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

    integer (int_kind) :: n, nfield, ns   ! loop indices

    do ns=1,nstreams
      if (tavg_streams(ns)%ltavg_file_is_open) then
        if (my_task.eq.master_task) then
          write(stdout,*) "tavg file ", trim(tavg_file_desc(ns)%full_name), &
                          " is still open... closing."
        end if
        ! Lots of variable clean up (necessary?)
        call destroy_io_field(time_coordinate(1,ns))
        do n=1,num_ccsm_coordinates
           call destroy_io_field(ccsm_coordinates(n,ns))
        end do
        do n=1,num_ccsm_time_invar(ns)
           call destroy_io_field(ccsm_time_invar(n,ns))
        end do
        do n=1,num_ccsm_scalars(ns)
           call destroy_io_field(ccsm_scalars(n,ns))
        end do
        n=0
        do nfield = 1, num_avail_tavg_fields
           if (tavg_in_this_stream(nfield, ns)) then
              n = n+1
              call destroy_io_field(tavg_streams(ns)%tavg_fields(n))
           end if
        end do

        tavg_streams(ns)%tavg_num_time_slices = 0
        deallocate(tavg_streams(ns)%tavg_fields)

        ! Close the file 
        call data_set(tavg_file_desc(ns), 'close')

        ! More variable clean up
        call destroy_file(tavg_file_desc(ns))
      end if
    end do

!-----------------------------------------------------------------------
!EOC

  end subroutine final_tavg
 
 end module tavg

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
