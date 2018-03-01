!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module movie

!BOP
! !MODULE: movie
! !DESCRIPTION:
!  This module contains data types and routines for computing running 
!  time-averages of selected fields and writing this data to files.
!
! !REVISION HISTORY:
!  CVS:$Id: movie.F90 41886 2012-11-13 16:56:30Z mlevy@ucar.edu $
!  CVS:$Name:  $

! !USES:

   use POP_IOUnitsMod
   use kinds_mod
   use blocks
   use distribution
   use domain
   use constants
   use prognostic
   use grid
   use time_management
   use registry
   use global_reductions
   use broadcast
   use io
   use exit_mod

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_movie,             &
             define_movie_field,     &
             update_movie_field,     &
             movie_requested,        &
             write_movie

! !PUBLIC DATA MEMBERS:

   logical (log_kind), public :: &
      lmovie_on      = .false. ! movie file output wanted

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  movie field descriptor data type and array of such types
!
!-----------------------------------------------------------------------

   type :: movie_field_desc
      character(char_len)     :: short_name     ! short name for field
      character(char_len)     :: long_name      ! long descriptive name
      character(char_len)     :: units          ! units
      character(4)            :: grid_loc       ! location in grid
      real (r4), dimension(2) :: valid_range    ! min/max
      real (r4)               :: fill_value     ! _FillValue
      integer (int_kind)      :: buf_loc        ! location in buffer
      integer (int_kind)      :: field_loc      ! grid location and field
      integer (int_kind)      :: field_type     ! type for io, ghost cells
      integer (r4)            :: field_depth_index  ! depth index of 2d slice
   end type

   integer (int_kind), parameter :: &
      max_avail_movie_fields = (4+nt)*km+50 ! limit on available fields - can
                                            !   be pushed as high as necessary

   integer (int_kind) ::           &
      num_avail_movie_fields = 0,   &! current number of defined fields
      num_requested_movie_fields,   &! number of fields requested
      movie_flag                     ! time flag for writing movie files

   type (movie_field_desc), dimension(max_avail_movie_fields) :: &
      avail_movie_fields

!-----------------------------------------------------------------------
!
!  buffers for holding running movie variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      movie_bufsize_2d      ! size of buffer for 2d fields

   real (r4), dimension(:,:,:,:), allocatable :: &
      MOVIE_BUF_2D         ! buffer for holding movie fields

!-----------------------------------------------------------------------
!
!  variables for writing data
!
!-----------------------------------------------------------------------

   integer (int_kind) ::     &
      movie_freq_iopt,  &! frequency option for writing movie
      movie_freq         ! frequency of movie output

   character (char_len) ::    &
      movie_outfile,           & ! root filename for movie output
      movie_fmt                  ! format (nc or bin) for writing

   type (datafile) :: movie_file_desc    ! IO file descriptor

   type (io_field_desc), target :: &
      MOVIE_iodesc                  ! io descriptor for movie fields

!-----------------------------------------------------------------------
!
!  ccsm variables
!
!-----------------------------------------------------------------------

   logical (log_kind) ::  &
      lccsm

!EOC
!***********************************************************************

 contains

!***********************************************************************
!EOP
! !IROUTINE: init_movie
! !INTERFACE:

 subroutine init_movie

! !DESCRIPTION:
!  This routine initializes movie options and reads in contents file to
!  determine which fields for which the user wants 2D snapshot data.
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

   integer (int_kind) ::         &
      n,                   &! dummy index
      k,                   &! depth index
      iblock,              &! local block index
      loc,                 &! location of field in buffer
      nu,                  &! unit for contents input file
      cindex,              &! character index for manipulating strings
      nml_error,           &! namelist i/o error flag
      contents_error        ! error flag for contents file read

   character (char_len) :: &
      movie_freq_opt,       &! choice for frequency of movie output
      movie_contents,       &! filename for choosing fields for output
      char_temp             ! temporary for manipulating fields

   character (34), parameter :: &
      freq_fmt = "('movie diagnostics every ',i6,a8)"

   namelist /movie_nml/ movie_freq_opt, movie_freq,        &
                       movie_outfile, movie_contents, movie_fmt

!-----------------------------------------------------------------------
!
!  read movie file output frequency and filenames from namelist
!
!-----------------------------------------------------------------------

   if (my_task == master_task) then
      write(stdout,delim_fmt)
      write(stdout,blank_fmt)
      write(stdout,'(a12)') 'Movie options'
      write(stdout,blank_fmt)
      write(stdout,delim_fmt)
   endif

   movie_freq_iopt = freq_opt_never
   movie_freq      = 100000
   movie_outfile   = 't'
   movie_contents  = 'unknown_movie_contents'

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=movie_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call exit_POP(sigAbort,'ERROR reading movie_nml')
   endif

   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      write(stdout,blank_fmt)
      write(stdout,*) ' Movie:'
      write(stdout,blank_fmt)
      write(stdout,*) ' movie_nml namelist settings:'
      write(stdout,blank_fmt)
      write(stdout,movie_nml)
      write(stdout,blank_fmt)
      call POP_IOUnitsFlush(stdout)
   endif

   if (my_task == master_task) then
      select case (movie_freq_opt)
      case ('never')
         movie_freq_iopt = freq_opt_never
         write(stdout,'(a21)') 'movie diagnostics off'
      case ('nyear')
         movie_freq_iopt = freq_opt_nyear
         write(stdout,freq_fmt) movie_freq,' years  '
      case ('nmonth')
         movie_freq_iopt = freq_opt_nmonth
         write(stdout,freq_fmt) movie_freq,' months '
      case ('nday')
         movie_freq_iopt = freq_opt_nday
         write(stdout,freq_fmt) movie_freq,' days   '
      case ('nhour')
         movie_freq_iopt = freq_opt_nhour
         write(stdout,freq_fmt) movie_freq,' hours  '
      case ('nsecond')
         movie_freq_iopt = freq_opt_nsecond
         write(stdout,freq_fmt) movie_freq,' seconds'
      case ('nstep')
         movie_freq_iopt = freq_opt_nstep
         write(stdout,freq_fmt) movie_freq,' steps  '
      case default
         movie_freq_iopt = -1000
      end select

   endif

   call POP_IOUnitsFlush(stdout)

   call broadcast_scalar(movie_freq_iopt, master_task)

   if (movie_freq_iopt == -1000) then
      call exit_POP(sigAbort,'unknown option for movie file frequency')
   else if (movie_freq_iopt /= freq_opt_never) then
      call broadcast_scalar(movie_freq,         master_task)
      call broadcast_scalar(movie_outfile,      master_task)
      call broadcast_scalar(movie_contents,     master_task)
      call broadcast_scalar(movie_fmt    ,      master_task)
   endif

!-----------------------------------------------------------------------
!
!  initialize time flag for writing movie files
!
!-----------------------------------------------------------------------

   call init_time_flag('movie',movie_flag, default=.false.,  &
                        owner    = 'init_movie',             &
                        freq_opt = movie_freq_iopt,          &
                        freq     = movie_freq)

!-----------------------------------------------------------------------
!
!  read contents file to determine which fields to dump
!
!-----------------------------------------------------------------------

   if (movie_freq_iopt /= freq_opt_never) then

      movie_bufsize_2d = 0

      call get_unit(nu)

      if (my_task == master_task) then
         open(nu, file=movie_contents, status='old')
         read(nu,*) num_requested_movie_fields
         write(stdout,'(a38)') 'movie diagnostics requested for fields:'
      endif

      call broadcast_scalar(num_requested_movie_fields, master_task)

      contents_error = 0

      do n=1,num_requested_movie_fields
         if (my_task == master_task) then
            read(nu,'(i3,a80)',iostat=contents_error) k, char_temp
            char_temp = adjustl(char_temp)
            cindex = index(char_temp,' ')
            char_temp(cindex:) = ' '
            write(stdout,*) '  ',trim(char_temp),' at level ',k
         endif

         call broadcast_scalar(contents_error, master_task)
         if (contents_error /= 0) then
            call exit_POP(sigAbort,'error reading movie contents')
         endif

         call broadcast_scalar(char_temp, master_task)
         call broadcast_scalar(k        , master_task)
         call request_movie_field(trim(char_temp), k)
      end do

      call release_unit(nu)

      !*** allocate and initialize running movie buffers

      allocate(                                                            &
         MOVIE_BUF_2D(nx_block,ny_block,   nblocks_clinic,movie_bufsize_2d) )

      lmovie_on = .true.

   endif

!-----------------------------------------------------------------------
!
!  determine if this is a ccsm coupled run
!
!-----------------------------------------------------------------------

   lccsm = registry_match('lcoupled')


!-----------------------------------------------------------------------
!EOC

 end subroutine init_movie

!***********************************************************************
!BOP
! !IROUTINE: write_movie
! !INTERFACE:

 subroutine write_movie

! !DESCRIPTION:
!  This routine writes requested movie fields to a file.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      nu,          &! i/o unit for output file
      iblock,      &! dummy block index
      nfield,      &! dummy field index
      loc,         &! buffer location for field
      io_phase      !'define' or 'write'

   character (char_len) ::  &
      file_suffix,          &! suffix to append to movie file name
      hist_string,          &! string containing file history
      movie_filename         ! filename for movie data

   character (8) :: &
      date_created   ! string with (real) date this file created

   character (10) :: &
      time_created   ! string with (real) date this file created

   type (io_field_desc), dimension(:), allocatable :: &
      movie_fields

   type (io_dim) :: &
      i_dim, j_dim ! dimension descriptors for horiz dims

   logical (log_kind) :: &
      lmovie_write       ! time to write a file

!-----------------------------------------------------------------------
!
!  is it time to write a file - if yes, create a file suffix
!
!-----------------------------------------------------------------------

   lmovie_write = .false.

   if (lmovie_on) then
      lmovie_write = check_time_flag(movie_flag)
   endif

   if (lmovie_write) then
      file_suffix = char_blank
      if (lccsm) then
         call create_suffix_movie_ccsm(file_suffix)
      else
         call create_suffix_movie(file_suffix)
      endif

!-----------------------------------------------------------------------
!
!     create data file descriptor
!
!-----------------------------------------------------------------------

      if (my_task.eq.master_task) then
        call date_and_time(date=date_created, time=time_created)
      end if
      call broadcast_scalar(date_created, master_task)
      call broadcast_scalar(time_created, master_task)
      hist_string = char_blank
      write(hist_string,'(a24,a8,1x,a10)') & 
         'POP MOVIE file created: ',date_created,time_created

      movie_file_desc = construct_file(movie_fmt,                    &
                                   root_name  = trim(movie_outfile),    &
                                   file_suffix= trim(file_suffix),     &
                                   title      ='POP MOVIE file',        &
                                   conventions='POP MOVIE conventions', &
                                   history    = trim(hist_string),     &
                                   record_length = rec_type_real,      &
                                   recl_words=nx_global*ny_global)

!-----------------------------------------------------------------------
!
!     add scalar fields to file as file attributes
!
!-----------------------------------------------------------------------

      call add_attrib_file(movie_file_desc, 'nsteps_total', nsteps_total)
      call add_attrib_file(movie_file_desc, 'tday'        , tday)
      call add_attrib_file(movie_file_desc, 'iyear'       , iyear)
      call add_attrib_file(movie_file_desc, 'imonth'      , imonth)
      call add_attrib_file(movie_file_desc, 'iday'        , iday)

!-----------------------------------------------------------------------
!
!     open output file
!
!-----------------------------------------------------------------------

      call data_set (movie_file_desc, 'open')

!-----------------------------------------------------------------------
!
!     write fields to file - this requires two phases
!     in this first phase, we define all the fields to be written
!
!-----------------------------------------------------------------------
 
      !*** define dimensions

      i_dim = construct_io_dim('i',nx_global)
      j_dim = construct_io_dim('j',ny_global)

      allocate(movie_fields(num_avail_movie_fields))

      do nfield = 1,num_avail_movie_fields  ! check all available fields

         loc = avail_movie_fields(nfield)%buf_loc ! locate field in buffer

         if (loc > 0) then  ! field is actually requested and in buffer

            !*** construct io_field descriptors for each field

               movie_fields(nfield) = construct_io_field(               &
                              avail_movie_fields(nfield)%short_name,    &
                              i_dim, j_dim,                            &
                    long_name=avail_movie_fields(nfield)%long_name,     &
                    units    =avail_movie_fields(nfield)%units    ,     &
                    grid_loc =avail_movie_fields(nfield)%grid_loc ,     &
                   field_loc =avail_movie_fields(nfield)%field_loc,     &
                  field_type =avail_movie_fields(nfield)%field_type,    &
                  valid_range=avail_movie_fields(nfield)%valid_range,   &
                   r2d_array =MOVIE_BUF_2D(:,:,:,loc) )

!-----------------------------------------------------------------------
!
!    missing_value is a deprecated feature in CF1.4, and hence nco 4 versions,
!    but it is added here because other software packages may require it
!-----------------------------------------------------------------------

           call add_attrib_io_field(movie_fields(nfield),'_FillValue',   &
                                    avail_movie_fields(nfield)%fill_value )
           call add_attrib_io_field(movie_fields(nfield),'missing_value',&
                                    avail_movie_fields(nfield)%fill_value )

            call data_set (movie_file_desc, 'define', movie_fields(nfield))
         endif
      end do

!-----------------------------------------------------------------------
!
!     write fields to file
!     in this second phase, we actually write the data for all the fields
!     after writing a field, the field descriptor is destroyed and the
!     file can be closed
!
!-----------------------------------------------------------------------
 
      do nfield = 1,num_avail_movie_fields  ! check all available fields

         loc = avail_movie_fields(nfield)%buf_loc ! locate field in buffer

         if (loc > 0) then  ! field is actually requested and in buffer
            call data_set (movie_file_desc, 'write', movie_fields(nfield))
            call destroy_io_field(movie_fields(nfield))
         endif
      end do

      deallocate(movie_fields)
      call data_set (movie_file_desc, 'close')

      if (my_task == master_task) then
         write(stdout,blank_fmt)
         write(stdout,*) 'Wrote file: ', trim(movie_file_desc%full_name)
      endif

!-----------------------------------------------------------------------
!
!     get rid of file descriptor
!
!-----------------------------------------------------------------------

      call destroy_file(movie_file_desc)
   endif ! lwrite_movie

!-----------------------------------------------------------------------
!EOC

 end subroutine write_movie

!***********************************************************************
!BOP
! !IROUTINE: movie_global
! !INTERFACE:

 subroutine movie_global

! !DESCRIPTION:
!  Calculates and print global integrals of time average fields
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

   integer (int_kind) ::     &
      k,               &   ! vertical level index
      ifield,          &   ! field identifier
      iblock,          &   ! block index
      nfield,          &   ! dummy field index
      field_loc,       &   ! field location (center,Nface,Eface,NEcorner)
      field_type           ! field type (scalar, vector, angle)

   real (r8) ::        &
      movie_field_sum,  &   ! sum of movie field
      movie_norm            ! normalization for average

   real (r8), dimension (:,:,:), allocatable ::  &
      WORK               ! temp for holding area_weighted field

   real (r8), dimension (:,:), allocatable ::  &
      RMASK              ! topography mask for global sum

!-----------------------------------------------------------------------
!
!  calculate globally-integrated time average of each chosen 2d field
!
!-----------------------------------------------------------------------

   allocate (RMASK(nx_block,ny_block), &
             WORK (nx_block,ny_block,nblocks_clinic))

   if (my_task == master_task) then
     write (stdout,blank_fmt)
     write (stdout,'(a22)') 'Global Time Averages: '
   endif

   do nfield=1,num_avail_movie_fields
      ifield = avail_movie_fields(nfield)%buf_loc
      if (ifield > 0) then

         field_loc  = avail_movie_fields(nfield)%field_loc
         field_type = avail_movie_fields(nfield)%field_type

            !$OMP PARALLEL DO PRIVATE(iblock)
            do iblock = 1,nblocks_clinic
               select case(field_loc)
               case(field_loc_center)
                  WORK(:,:,iblock)  = MOVIE_BUF_2D(:,:,iblock,ifield)* &
                                    TAREA(:,:,iblock)*RCALCT(:,:,iblock)
               case(field_loc_NEcorner)
                  WORK(:,:,iblock)  = MOVIE_BUF_2D(:,:,iblock,ifield)* &
                                    UAREA(:,:,iblock)*RCALCU(:,:,iblock)
               case default ! make U cell the default for all other cases
                  WORK(:,:,iblock)  = MOVIE_BUF_2D(:,:,iblock,ifield)* &
                                    UAREA(:,:,iblock)*RCALCU(:,:,iblock)
               end select
            end do
            !$OMP END PARALLEL DO

            movie_field_sum = global_sum(WORK, distrb_clinic, field_loc)

            select case(field_loc)
            case(field_loc_center)
               movie_field_sum = movie_field_sum/(area_t)
            case(field_loc_NEcorner)
               movie_field_sum = movie_field_sum/(area_u)
            case default ! make U cell the default for all other cases
               movie_field_sum = movie_field_sum/(area_u)
            end select

         if (my_task == master_task) then
            write (stdout,*) trim(avail_movie_fields(nfield)%short_name), &
                             ': ', movie_field_sum
         endif
      endif
   end do

   deallocate (RMASK, WORK)

!-----------------------------------------------------------------------
!EOC

 end subroutine movie_global

!***********************************************************************
!BOP
! !IROUTINE: update_movie_field
! !INTERFACE:

 subroutine update_movie_field(ARRAY,field_id,block,k)

! !DESCRIPTION:
!  This routine updates a movie field to the current value.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      block,           &! local block address (in baroclinic distribution)
      k,               &! vertical level
      field_id          ! index into available fields for movie field info

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      ARRAY             ! array of data for this block update movie buffer

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      bufloc            ! location of field in movie buffer

!-----------------------------------------------------------------------
!
!  get buffer location and field info from avail_movie_field array
!
!-----------------------------------------------------------------------

   bufloc = avail_movie_fields(field_id)%buf_loc
   if (bufloc <= 0) &
     call exit_POP(sigAbort, &
                    'movie: attempt to update bad movie field')

!-----------------------------------------------------------------------
!
!  update the field into the movie buffer
!
!-----------------------------------------------------------------------

   MOVIE_BUF_2D(:,:,block,bufloc) = ARRAY

!-----------------------------------------------------------------------
!EOC

 end subroutine update_movie_field

!***********************************************************************
!BOP
! !IROUTINE: define_movie_field
! !INTERFACE:

 subroutine define_movie_field(id, short_name, depth_index,  &
                                  long_name, units, &
                                  grid_loc,  valid_range, &
                                  field_loc, field_type)

! !DESCRIPTION:
!  Initializes description of an available field and returns location
!  in the available fields array for use in later movie calls.
!
! !REVISION HISTORY:
!  same as module

! !OUTPUT PARAMETERS:

   integer (int_kind), intent(out) :: &
      id                ! location in avail_fields array for use in
                        ! later movie routines

! !INPUT PARAMETERS:

   character(*), intent(in) :: &
      short_name               ! short name for field

   integer (int_kind), intent(in), optional :: &
      field_loc,              &! location in grid 
      field_type,             &! type of field (scalar, vector, angle)
      depth_index              ! depth index of 2d slice

   character(*), intent(in), optional :: &
      long_name,              &! long descriptive name for field
      units                    ! physical units for field

   character(4), intent(in), optional :: &
      grid_loc                 ! location in grid (in 4-digit code)

   real (r4), dimension(2), intent(in), optional :: &
      valid_range              ! min/max

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   character (char_len) :: &
      appended_long_name,  &  !  long name with depth appended
      appended_short_name     !  short name with depth appended

   character (len = 5) :: char_depth  ! character version of the depth of a 2d slice

   integer (int_kind) :: &
      cbegin, clen, cindx, & !  character indices
      nearest_integer_depth    !  integer version of the depth of a 2d slice

   real (r4) ::  &
      field_depth      ! floating point version of the depth of a 2d slice
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  increment the number of defined fields and make sure it does not
!  exceed the maximum
!  return the id as the current number
!
!-----------------------------------------------------------------------

   num_avail_movie_fields = num_avail_movie_fields + 1
   if (num_avail_movie_fields > max_avail_movie_fields) then
      call exit_POP(sigAbort, &
                    'movie: defined movie fields > max allowed')
   endif
 
   id = num_avail_movie_fields

!-----------------------------------------------------------------------
!
!  now fill the field descriptor
!
!-----------------------------------------------------------------------

   avail_movie_fields(id)%buf_loc    = 0  ! will be reset later

!-----------------------------------------------------------------------
!
!  now check the depth index since we will modify short and long field
!    names if not at surface (depth_index = 0)
!
!  for example, if we want TEMP for k=5 which corresponds to, say, 63 meters
!    depth, then the short_name becomes TEMP_63m
!  if long_name is defined, do a similar thing except replace _ with at
!    and m with meters
!
!-----------------------------------------------------------------------

   if (present(depth_index)) then
      avail_movie_fields(id)%field_depth_index = depth_index
   else
      avail_movie_fields(id)%field_depth_index = 0
   endif

   if (present(long_name) .and. avail_movie_fields(id)%field_depth_index <= 0) then
      avail_movie_fields(id)%long_name = long_name
   else
      avail_movie_fields(id)%long_name = char_blank
   endif

   if (avail_movie_fields(id)%field_depth_index > 0) then
      field_depth = zt(depth_index)*mpercm  !  assume mid-cell and convert to meters
      nearest_integer_depth = nint(field_depth)
      write(char_depth,'(i5)') 10000 + nearest_integer_depth
      if (nearest_integer_depth < 10) then
         cbegin = 5
      else if (nearest_integer_depth >= 10 .and. nearest_integer_depth < 100) then
         cbegin = 4
      else if (nearest_integer_depth >= 100 .and. nearest_integer_depth < 1000) then
         cbegin = 3
      else if (nearest_integer_depth >= 1000 .and. nearest_integer_depth < 10000) then
         cbegin = 2
      else
         cbegin = 1
      endif
      clen = 6 - cbegin
      appended_short_name = short_name
      cindx = len_trim(appended_short_name)
      cindx = cindx + 1
      appended_short_name(cindx:cindx) = '_'
      cindx = cindx + 1
      appended_short_name(cindx:cindx+clen-1) =   &
         char_depth(cbegin:cbegin+clen-1)
      cindx = cindx + clen
      appended_short_name(cindx:cindx) = 'm'  !  meters
      avail_movie_fields(id)%short_name = appended_short_name
      if (present(long_name)) then
         appended_long_name = long_name
         cindx = len_trim(appended_long_name)
         cindx = cindx + 1
         appended_long_name(cindx:cindx+3) = ' at '
         cindx = cindx + 4
         appended_long_name(cindx:cindx+clen-1) =   &
            char_depth(cbegin:cbegin+clen-1)
         cindx = cindx + clen
         appended_long_name(cindx:cindx+6) = ' meters'  !  meters
         avail_movie_fields(id)%long_name = appended_long_name
      endif   !  long_name is present
   else
      avail_movie_fields(id)%short_name = short_name
   endif

   if (present(units)) then
      avail_movie_fields(id)%units = units
   else
      avail_movie_fields(id)%units = char_blank
   endif

   if (present(grid_loc)) then
      avail_movie_fields(id)%grid_loc = grid_loc
   else
      avail_movie_fields(id)%grid_loc = '    '
   endif

   avail_movie_fields(id)%fill_value = undefined_nf_r4

   if (present(valid_range)) then
      avail_movie_fields(id)%valid_range = valid_range
   else
      avail_movie_fields(id)%valid_range = undefined
   endif

   !*** set field location, field type used by i/o, ghost cell update
   !*** and other communication routines.  because ghost cells for movie
   !*** fields are not typically used, the default is field_xxx_noupdate

   if (present(field_loc)) then
      avail_movie_fields(id)%field_loc = field_loc
   else
      !*** try to decode field location from grid_loc
      if (grid_loc(2:2) == '1' .and. grid_loc(3:3) == '1') then
         avail_movie_fields(id)%field_loc = field_loc_center
      else if (grid_loc(2:2) == '2' .and. grid_loc(3:3) == '2') then
         avail_movie_fields(id)%field_loc = field_loc_NEcorner
      else if (grid_loc(2:2) == '1' .and. grid_loc(3:3) == '2') then
         avail_movie_fields(id)%field_loc = field_loc_Nface
      else if (grid_loc(2:2) == '2' .and. grid_loc(3:3) == '1') then
         avail_movie_fields(id)%field_loc = field_loc_Eface
      else
         avail_movie_fields(id)%field_loc = field_loc_noupdate
      endif
   endif

   if (present(field_type)) then
      avail_movie_fields(id)%field_type = field_type
   else
      avail_movie_fields(id)%field_type = field_type_noupdate
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine define_movie_field

!***********************************************************************
!BOP
! !IROUTINE: request_movie_field
! !INTERFACE:

 subroutine request_movie_field(short_name,k)

! !DESCRIPTION:
!  This field marks an available field as requested and computes
!  the location in the movie buffer array.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      short_name                ! the short name of the field

   integer (int_kind), intent(in) :: &
      k                 ! depth index

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      n,                 &! loop index
      id                  ! location of field in avail_fields array

   character (char_len) :: &
      appended_short_name     !  short name with depth appended

   character (len = 5) :: char_depth  ! character version of the depth of a 2d slice

   integer (int_kind) :: &
      cbegin, clen, cindx, & !  character indices
      nearest_integer_depth    !  integer version of the depth of a 2d slice

   real (r4) ::  &
      field_depth      ! floating point version of the depth of a 2d slice

!-----------------------------------------------------------------------
!
!  search for field with same name
!
!-----------------------------------------------------------------------

   if ( k > 0 .and. k <= km) then
     field_depth = zt(k)*mpercm  !  assume mid-cell and convert to meters
     nearest_integer_depth = nint(field_depth)
     write(char_depth,'(i5)') 10000 + nearest_integer_depth
     if (nearest_integer_depth < 10) then
      cbegin = 5
     else if (nearest_integer_depth >= 10 .and. nearest_integer_depth < 100) then
      cbegin = 4
     else if (nearest_integer_depth >= 100 .and. nearest_integer_depth < 1000) then
      cbegin = 3
     else if (nearest_integer_depth >= 1000 .and. nearest_integer_depth < 10000) then
      cbegin = 2
     else
      cbegin = 1
     endif
     clen = 6 - cbegin
     appended_short_name = short_name
     cindx = len_trim(appended_short_name)
     cindx = cindx + 1
     appended_short_name(cindx:cindx) = '_'
     cindx = cindx + 1
     appended_short_name(cindx:cindx+clen-1) =   &
      char_depth(cbegin:cbegin+clen-1)
     cindx = cindx + clen
     appended_short_name(cindx:cindx) = 'm'  !  meters

   else

     appended_short_name = short_name

   endif

   id = 0
   srch_loop: do n=1,num_avail_movie_fields
      if (trim(avail_movie_fields(n)%short_name) == trim(appended_short_name)) then
         id = n
         exit srch_loop
      endif
   end do srch_loop

   if (id == 0) then
      if (my_task == master_task) write(stdout,*) 'Requested ', &
                                                  trim(appended_short_name)
      call exit_POP(sigAbort,'movie: requested field unknown')
   endif

!-----------------------------------------------------------------------
!
!  set the position in the buffer and advance the buffer position
!  for the next field
!
!-----------------------------------------------------------------------

   movie_bufsize_2d = movie_bufsize_2d + 1
   avail_movie_fields(id)%buf_loc = movie_bufsize_2d

!-----------------------------------------------------------------------
!EOC

 end subroutine request_movie_field

!***********************************************************************
!BOP
! !IROUTINE: movie_requested
! !INTERFACE:

 function movie_requested(id)

! !DESCRIPTION:
!  This function determines whether an available (defined) movie field
!  has been requested by a user (through the input contents file) and 
!  returns true if it has.  Note that if movie has been turned off, 
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
      movie_requested     ! result of checking whether the field has
                         !   been requested

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  check the buffer location - if zero, the field has not been
!  requested
!
!-----------------------------------------------------------------------

   if (id < 1 .or. id > num_avail_movie_fields) then
      call exit_POP(sigAbort,'movie_requested: invalid movie id')
   endif

   if (avail_movie_fields(id)%buf_loc > 0) then
      movie_requested = .true.
   else
      movie_requested = .false.
   endif

!-----------------------------------------------------------------------
!EOC

 end function movie_requested

!***********************************************************************
!BOP
! !IROUTINE: create_suffix_movie
! !INTERFACE:

 subroutine create_suffix_movie(file_suffix)

! !DESCRIPTION:
!  Creates a suffix to append to output filename based on frequency 
!  option and averaging interval.
!
! !REVISION HISTORY:
!  same as module

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

   !***
   !*** use step numbers if movie freq option is nstep
   !***

   write(cstep_end,'(i10)') nsteps_total - 1
   cdate  = adjustl(cstep_end)
   cstep_end = trim(cdate)

   call time_stamp('last', 'ymd', date_string=cdate)  ! last date

   if (date_separator == ' ') then  ! no date separator
      cyear_end  = cdate(1:4)
      cmonth_end = cdate(5:6)
      cday_end   = cdate(7:8)
   else
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

   select case (movie_freq_iopt)
   case (freq_opt_nyear, freq_opt_nmonth, freq_opt_nday)
      cindx2 = cindx1 + 7
      file_suffix(cindx1:cindx2) = cyear_end/&
                                 &/cmonth_end/&
                                 &/cday_end

   case (freq_opt_nstep)
      cindx2 = cindx1 + len_trim(cstep_end) - 1
      file_suffix(cindx1:cindx2) = trim(cstep_end)

   case default  ! use nstep for other options
      cindx2 = cindx1 + len_trim(cstep_end) - 1
      file_suffix(cindx1:cindx2) = trim(cstep_end)

   end select
 
!-----------------------------------------------------------------------
!EOC

 end subroutine create_suffix_movie


!***********************************************************************
!BOP
! !IROUTINE: create_suffix_movie_ccsm
! !INTERFACE:

 subroutine create_suffix_movie_ccsm(file_suffix)

! !DESCRIPTION:
!  Creates a suffix to append to output filename based on frequency 
!  option and averaging interval. Suffix conforms to CCSM output
!  file file-naming conventions.
!
! !REVISION HISTORY:
!  same as module

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
!  for a ccsm movie file, append a date/time string to the root name
!
!-----------------------------------------------------------------------

      select case (movie_freq_iopt)
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

 
!-----------------------------------------------------------------------
!EOC

 end subroutine create_suffix_movie_ccsm

!***********************************************************************

 end module movie

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
