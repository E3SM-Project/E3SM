! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
module mpas_io

   use mpas_derived_types
   use mpas_attlist
   use mpas_dmpar
   use mpas_log

   use pio
   use pio_kinds
   use pio_types

#ifdef SINGLE_PRECISION
   integer, parameter :: PIO_REALKIND = PIO_REAL
#else
   integer, parameter :: PIO_REALKIND = PIO_DOUBLE
#endif

#ifdef USE_PIO2
   integer, parameter :: MPAS_INT_FILLVAL = PIO_FILL_INT
   character, parameter :: MPAS_CHAR_FILLVAL = achar(0)  ! TODO: To be replaced with PIO_FILL_CHAR once PIO2 provides this variable
#else
   integer, parameter :: MPAS_INT_FILLVAL = NF_FILL_INT
   character, parameter :: MPAS_CHAR_FILLVAL = achar(NF_FILL_CHAR)
#endif

#ifdef USE_PIO2
#ifdef SINGLE_PRECISION
   real (kind=RKIND), parameter :: MPAS_REAL_FILLVAL = PIO_FILL_FLOAT
#else
   real (kind=RKIND), parameter :: MPAS_REAL_FILLVAL = PIO_FILL_DOUBLE
#endif
#else
#ifdef SINGLE_PRECISION
   real (kind=RKIND), parameter :: MPAS_REAL_FILLVAL = NF_FILL_FLOAT
#else
   real (kind=RKIND), parameter :: MPAS_REAL_FILLVAL = NF_FILL_DOUBLE
#endif
#endif

   interface MPAS_io_get_var
      module procedure MPAS_io_get_var_int0d
      module procedure MPAS_io_get_var_int1d
      module procedure MPAS_io_get_var_int2d
      module procedure MPAS_io_get_var_int3d
      module procedure MPAS_io_get_var_int4d
      module procedure MPAS_io_get_var_real0d
      module procedure MPAS_io_get_var_real1d
      module procedure MPAS_io_get_var_real2d
      module procedure MPAS_io_get_var_real3d
      module procedure MPAS_io_get_var_real4d
      module procedure MPAS_io_get_var_real5d
      module procedure MPAS_io_get_var_char0d
      module procedure MPAS_io_get_var_char1d
   end interface MPAS_io_get_var

   interface MPAS_io_put_var
      module procedure MPAS_io_put_var_int0d
      module procedure MPAS_io_put_var_int1d
      module procedure MPAS_io_put_var_int2d
      module procedure MPAS_io_put_var_int3d
      module procedure MPAS_io_put_var_int4d
      module procedure MPAS_io_put_var_real0d
      module procedure MPAS_io_put_var_real1d
      module procedure MPAS_io_put_var_real2d
      module procedure MPAS_io_put_var_real3d
      module procedure MPAS_io_put_var_real4d
      module procedure MPAS_io_put_var_real5d
      module procedure MPAS_io_put_var_char0d
      module procedure MPAS_io_put_var_char1d
   end interface MPAS_io_put_var

   interface MPAS_io_get_att
      module procedure MPAS_io_get_att_int0d
      module procedure MPAS_io_get_att_int1d
      module procedure MPAS_io_get_att_real0d
      module procedure MPAS_io_get_att_real1d
      module procedure MPAS_io_get_att_text
   end interface MPAS_io_get_att

   interface MPAS_io_put_att
      module procedure MPAS_io_put_att_int0d
      module procedure MPAS_io_put_att_int1d
      module procedure MPAS_io_put_att_real0d
      module procedure MPAS_io_put_att_real1d
      module procedure MPAS_io_put_att_text
   end interface MPAS_io_put_att

   contains


   subroutine MPAS_io_init(ioContext, io_task_count, io_task_stride, io_system, ierr)

      implicit none

      type (mpas_io_context_type), intent(inout) :: ioContext
      integer, intent(in) :: io_task_count
      integer, intent(in) :: io_task_stride
      type (iosystem_desc_t), optional, pointer :: io_system
      integer, intent(out), optional :: ierr

      integer :: local_ierr

      local_ierr = 0

!      call mpas_log_write('Called MPAS_io_init()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      if (present(io_system)) then
         ioContext % pio_iosystem => io_system
      else
!call mpas_log_write('MGD PIO_init')
         if ( io_task_count < 0 .or. io_task_count > ioContext % dminfo % nprocs ) then
            call mpas_log_write('PIO io_task_count has a value of $i.', MPAS_LOG_ERR, intArgs=(/io_task_count/))
            call mpas_log_write('        It must be between 1 and $i.', MPAS_LOG_ERR, intArgs=(/ioContext % dminfo % nprocs/) )
            local_ierr = 1
         end if

         if ( io_task_stride <= 0 .or. io_task_stride > ioContext % dminfo % nprocs ) then
            call mpas_log_write('PIO io_task_stride has a value of $i.', MPAS_LOG_ERR, intArgs=(/io_task_stride/) )
            call mpas_log_write('        It must be between 1 and $i.', MPAS_LOG_ERR, intArgs=(/ioContext % dminfo % nprocs/) )
            local_ierr = 1
         end if

         if ( io_task_stride * io_task_count > ioContext % dminfo % nprocs ) then
            call mpas_log_write('PIO io_task_stride * io_task_count has a value of $i.', MPAS_LOG_ERR, &
                    intArgs=(/io_task_stride * io_task_count/) )
            call mpas_log_write('        It must be between 1 and $i.', MPAS_LOG_ERR, intArgs=(/ioContext % dminfo % nprocs/) )
            local_ierr = 1
         end if

         if ( local_ierr /= 0 ) then
            call mpas_log_write('Invalid PIO configuration.', MPAS_LOG_CRIT)
         end if

         allocate(ioContext % pio_iosystem)
         call PIO_init(ioContext % dminfo % my_proc_id, &  ! comp_rank
                       ioContext % dminfo % comm,       &  ! comp_comm
                       io_task_count,             &     ! num_iotasks
                       0,                         &     ! num_aggregator
                       io_task_stride,            &     ! stride
                       PIO_rearr_box,             &     ! rearr
                       ioContext % pio_iosystem)           ! iosystem
  
      end if

      call pio_seterrorhandling(ioContext % pio_iosystem, PIO_BCAST_ERROR)

   end subroutine MPAS_io_init


!***********************************************************************
!
!  routine MPAS_io_set_iotype
!
!> \brief   Set master PIO io type
!> \author  Doug Jacobsen
!> \date    10/18/2013
!> \details 
!>  This routine sets the master io type for use with PIO.
!
!-----------------------------------------------------------------------
   subroutine MPAS_io_set_iotype(ioContext, io_type_in, ierr)

      implicit none

      type (mpas_io_context_type), intent(inout) :: ioContext
      integer, intent(in) :: io_type_in
      integer, intent(out), optional :: ierr

      if (present(ierr)) then
         ierr = MPAS_IO_NOERR
      end if

      ioContext % master_pio_iotype = io_type_in

   end subroutine MPAS_io_set_iotype


!***********************************************************************
!
!  routine MPAS_io_unset_iotype
!
!> \brief   Unset master PIO io type
!> \author  Doug Jacobsen
!> \date    10/18/2013
!> \details 
!>  This routine sets the master io type for use with PIO to it's default
!>  "unset" value.
!
!-----------------------------------------------------------------------
   subroutine MPAS_io_unset_iotype(ioContext, ierr)

      implicit none

      type (mpas_io_context_type), intent(inout) :: ioContext
      integer, intent(out), optional :: ierr

      if (present(ierr)) then
         ierr = MPAS_IO_NOERR
      end if

      ioContext % master_pio_iotype = -999

   end subroutine MPAS_io_unset_iotype

   
   type (MPAS_IO_Handle_type) function MPAS_io_open(filename, mode, ioformat, ioContext, clobber_file, truncate_file, ierr)

      implicit none

      character (len=*), intent(in) :: filename
      integer, intent(in) :: mode
      integer, intent(in) :: ioformat
      type (mpas_io_context_type), pointer :: ioContext
      logical, intent(in), optional :: clobber_file
      logical, intent(in), optional :: truncate_file
      integer, intent(out), optional :: ierr

      integer :: pio_ierr, pio_iotype, pio_mode
      logical :: local_clobber, local_truncate
      logical :: exists

!      call mpas_log_write('Called MPAS_io_open()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      if (present(clobber_file)) then
         local_clobber = clobber_file
      else
         local_clobber = .false.
      end if

      if (present(truncate_file)) then
         local_truncate = truncate_file
      else
         local_truncate = .false.
      end if

      ! Sanity checks
      if (mode /= MPAS_IO_READ .and. &
          mode /= MPAS_IO_WRITE) then
         if (present(ierr)) ierr = MPAS_IO_ERR_INVALID_MODE
         return 
      end if
      if (ioformat /= MPAS_IO_NETCDF .and. &
          ioformat /= MPAS_IO_NETCDF4 .and. &
          ioformat /= MPAS_IO_PNETCDF .and. &
          ioformat /= MPAS_IO_PNETCDF5) then
         if (present(ierr)) ierr = MPAS_IO_ERR_INVALID_FORMAT
         return 
      end if
      if (len(filename) > 1024) then
         if (present(ierr)) ierr = MPAS_IO_ERR_LONG_FILENAME
         return 
      end if

      MPAS_io_open % filename = filename
      MPAS_io_open % iomode   = mode
      MPAS_io_open % ioformat = ioformat
      MPAS_io_open % ioContext => ioContext

      if (ioContext % master_pio_iotype /= -999) then
         pio_iotype = ioContext % master_pio_iotype
         pio_mode = PIO_64BIT_OFFSET
      else
         if (ioformat == MPAS_IO_PNETCDF) then
            pio_iotype = PIO_iotype_pnetcdf
            pio_mode = PIO_64BIT_OFFSET
         else if (ioformat == MPAS_IO_PNETCDF5) then
            pio_iotype = PIO_iotype_pnetcdf
            pio_mode = PIO_64BIT_DATA
         else if (ioformat == MPAS_IO_NETCDF) then
            pio_iotype = PIO_iotype_netcdf
            pio_mode = PIO_64BIT_OFFSET
         else if (ioformat == MPAS_IO_NETCDF4) then
            pio_iotype = PIO_iotype_netcdf4p
#ifdef USE_PIO2
            pio_mode = 0
#else
            pio_mode = PIO_64BIT_OFFSET
#endif
         end if
      end if

      if (mode == MPAS_IO_WRITE) then
!call mpas_log_write('MGD PIO_createfile')
         if (ioContext % dminfo % my_proc_id == 0) then
            inquire(file=trim(filename), exist=exists)
         end if
         call mpas_dmpar_bcast_logical(ioContext % dminfo, exists)

         ! If the file exists and we are not allowed to clobber it, return an
         !    appropriate error code
         if (exists .and. (.not. local_clobber)) then
            if (present(ierr)) ierr = MPAS_IO_ERR_WOULD_CLOBBER
            return
         end if
 
         if (exists .and. (.not. local_truncate)) then
            pio_ierr = PIO_openfile(ioContext % pio_iosystem, MPAS_io_open % pio_file, pio_iotype, trim(filename), PIO_write)
            MPAS_io_open % preexisting_file = .true.
         else
            pio_ierr = PIO_createfile(ioContext % pio_iosystem, MPAS_io_open % pio_file, pio_iotype, trim(filename), pio_mode)
#ifdef MPAS_DEBUG
            if (exists) then
               call mpas_log_write('MPAS I/O: Truncating existing data in output file '//trim(filename), MPAS_LOG_WARN)
            end if
#endif
         end if
      else
         inquire(file=trim(filename), exist=exists)

         if (.not. exists) then
            if (present(ierr)) ierr = MPAS_IO_ERR_NOEXIST_READ
            return
         end if
!call mpas_log_write('MGD PIO_openfile')
         pio_ierr = PIO_openfile(ioContext % pio_iosystem, MPAS_io_open % pio_file, pio_iotype, trim(filename), PIO_nowrite)
      endif
      if (pio_ierr /= PIO_noerr) then
         if (present(ierr)) ierr = MPAS_IO_ERR_PIO
         return
      end if

      if (mode == MPAS_IO_READ .or. MPAS_io_open % preexisting_file) then
!MPAS_io_open % pio_unlimited_dimid = 44
         pio_ierr = PIO_inquire(MPAS_io_open % pio_file, unlimitedDimID=MPAS_io_open % pio_unlimited_dimid)
!call mpas_log_write('Found unlimited dim $i', intArgs=(/MPAS_io_open % pio_unlimited_dimid/) )
         if (pio_ierr /= PIO_noerr) then
            if (present(ierr)) ierr = MPAS_IO_ERR_PIO
            return
         end if

         ! Here we're depending on the undocumented behavior of PIO to return a
         ! negative dimension ID when an unlimited dimension is not found.  This
         ! might change in the future, causing this code to break, though it
         ! shouldn't break for files with unlimited dimensions.
         if ( MPAS_io_open % pio_unlimited_dimid >= 0 ) then
            pio_ierr = PIO_inq_dimlen(MPAS_io_open % pio_file, MPAS_io_open % pio_unlimited_dimid, MPAS_io_open % preexisting_records)
            if (pio_ierr /= PIO_noerr) then
               if (present(ierr)) ierr = MPAS_IO_ERR_PIO
               return
            end if
         else
            MPAS_io_open % preexisting_records = -1
         end if
      end if

      MPAS_io_open % initialized = .true.

      return

   end function MPAS_io_open


   subroutine MPAS_io_inq_unlimited_dim(handle, dimname, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(out) :: dimname
      integer, intent(out), optional :: ierr

      integer :: pio_ierr

!      call mpas_log_write('Called MPAS_io_inq_unlimited_dim()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      ! Sanity checks
      if (.not. handle % initialized) then
         if (present(ierr)) ierr = MPAS_IO_ERR_UNINIT_HANDLE
         return 
      end if
      if (handle % iomode /= MPAS_IO_READ) then       ! We could eventually handle this for write mode, too...
         if (present(ierr)) ierr = MPAS_IO_ERR_WRONG_MODE
         return
      end if

      pio_ierr = PIO_inq_dimname(handle % pio_file, handle % pio_unlimited_dimid, dimname) 
      if (pio_ierr /= PIO_noerr) then
         if (present(ierr)) ierr = MPAS_IO_ERR_NO_UNLIMITED_DIM
         dimname = ' '
         return
      end if

   end subroutine MPAS_io_inq_unlimited_dim


   subroutine MPAS_io_inq_dim(handle, dimname, dimsize, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: dimname
      integer, intent(out) :: dimsize
      integer, intent(out), optional :: ierr

      type (dimlist_type), pointer :: new_dimlist_node
      type (dimlist_type), pointer :: dim_cursor
      integer :: pio_ierr

!      call mpas_log_write('Called MPAS_io_inq_dim()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      ! Sanity checks
      if (.not. handle % initialized) then
         if (present(ierr)) ierr = MPAS_IO_ERR_UNINIT_HANDLE
         return 
      end if


      !
      ! First see if we already have this dimension in our list
      !
      dim_cursor => handle % dimlist_head
      do while (associated(dim_cursor))
         if (trim(dimname) == trim(dim_cursor % dimhandle % dimname)) then
            dimsize = dim_cursor % dimhandle % dimsize
            return
         end if
         dim_cursor => dim_cursor % next
      end do


      !
      ! Otherwise, query the file-level API for information about the dim
      !
      allocate(new_dimlist_node)
      nullify(new_dimlist_node % next)
      allocate(new_dimlist_node % dimhandle)

      new_dimlist_node % dimhandle % dimname = dimname

      pio_ierr = PIO_inq_dimid(handle % pio_file, trim(dimname), new_dimlist_node % dimhandle % dimid)
      if (pio_ierr /= PIO_noerr) then
         if (present(ierr)) ierr = MPAS_IO_ERR_MISSING_DIM
         deallocate(new_dimlist_node % dimhandle)
         deallocate(new_dimlist_node)
         dimsize = -1
         return
      end if

      if (new_dimlist_node % dimhandle % dimid == handle % pio_unlimited_dimid) new_dimlist_node % dimhandle % is_unlimited_dim = .true.

      pio_ierr = PIO_inq_dimlen(handle % pio_file, new_dimlist_node % dimhandle % dimid, new_dimlist_node % dimhandle % dimsize)
      if (pio_ierr /= PIO_noerr) then
         if (present(ierr)) ierr = MPAS_IO_ERR_PIO
         deallocate(new_dimlist_node % dimhandle)
         deallocate(new_dimlist_node)
         dimsize = -1
         return
      end if
   
      ! Keep dimension information for future reference
      if (.not. associated(handle % dimlist_head)) then
         handle % dimlist_head => new_dimlist_node
      end if
      if (.not. associated(handle % dimlist_tail)) then
         handle % dimlist_tail => new_dimlist_node
      else
         handle % dimlist_tail % next => new_dimlist_node
         handle % dimlist_tail => handle % dimlist_tail % next
      end if

      dimsize = new_dimlist_node % dimhandle % dimsize

   end subroutine MPAS_io_inq_dim


   subroutine MPAS_io_def_dim(handle, dimname, dimsize, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: dimname
      integer, intent(in) :: dimsize
      integer, intent(out), optional :: ierr

      integer :: pio_ierr
      integer :: inq_dimsize
      type (dimlist_type), pointer :: new_dimlist_node
      type (dimlist_type), pointer :: dim_cursor

!      call mpas_log_write('Called MPAS_io_def_dim()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      ! Sanity checks
      if (.not. handle % initialized) then
         if (present(ierr)) ierr = MPAS_IO_ERR_UNINIT_HANDLE
         return 
      end if
      if (handle % data_mode) then
         if (present(ierr)) ierr = MPAS_IO_ERR_DATA_MODE
         return 
      end if
      if (handle % iomode /= MPAS_IO_WRITE) then
         if (present(ierr)) ierr = MPAS_IO_ERR_NOWRITE
         return 
      end if


      !
      ! Check that this dimension hasn't already been defined
      !
      dim_cursor => handle % dimlist_head
      do while (associated(dim_cursor))
         if (trim(dimname) == trim(dim_cursor % dimhandle % dimname)) then

            ! The second half of the test below avoids raising errors in the case where 
            !    we are writing to an already existing file, in which case the dimlen for the
            !    unlimited dimension in the file will generally not be MPAS_IO_UNLIMITED_DIM
            if ((dimsize /= dim_cursor % dimhandle % dimsize) .and. &
               .not. (dimsize == MPAS_IO_UNLIMITED_DIM .and. dim_cursor % dimhandle % is_unlimited_dim)) then
               if (present(ierr)) ierr = MPAS_IO_ERR_REDEF_DIM
            end if
            return
         end if

         ! Also, check that the user is not trying to define more than one record dimension
         if (dimsize == MPAS_IO_UNLIMITED_DIM .and. dim_cursor % dimhandle % is_unlimited_dim) then
            if (present(ierr)) ierr = MPAS_IO_ERR_TWO_UNLIMITED_DIMS
            return
         end if
         dim_cursor => dim_cursor % next
      end do


      !
      ! If we're working with a preexisting file, see whether this dimension isn't already
      !    defined in the file; if it is, the dimension should match the dimsize specified in
      !    this call to define the dimension, at which point we can add it to our local cache
      !
      if (handle % preexisting_file) then
         call MPAS_io_inq_dim(handle, dimname, inq_dimsize, ierr=pio_ierr)
         if (pio_ierr /= MPAS_IO_ERR_PIO) then

            ! Verify that the dimsize matches...
            if (dimsize /= inq_dimsize .and. dimsize /= MPAS_IO_UNLIMITED_DIM) then
               if (present(ierr)) ierr = MPAS_IO_ERR_REDEF_DIM
               return
            end if

         end if
         return
      end if


      !
      ! Otherwise, define it
      !
      allocate(new_dimlist_node)
      nullify(new_dimlist_node % next)
      allocate(new_dimlist_node % dimhandle)

      new_dimlist_node % dimhandle % dimname = dimname
      new_dimlist_node % dimhandle % dimsize = dimsize
      if (dimsize == MPAS_IO_UNLIMITED_DIM) then
         new_dimlist_node % dimhandle % is_unlimited_dim = .true.
         pio_ierr = PIO_def_dim(handle % pio_file, trim(dimname), PIO_unlimited, new_dimlist_node % dimhandle % dimid)
      else
         pio_ierr = PIO_def_dim(handle % pio_file, trim(dimname), dimsize, new_dimlist_node % dimhandle % dimid)
      end if
      if (pio_ierr /= PIO_noerr) then
         if (present(ierr)) ierr = MPAS_IO_ERR_PIO
         deallocate(new_dimlist_node % dimhandle)
         deallocate(new_dimlist_node)
         return
      end if

      ! Keep dimension information
      if (.not. associated(handle % dimlist_head)) then
         handle % dimlist_head => new_dimlist_node
!call mpas_log_write('Assigning head for '//trim(dimname))
      end if
      if (.not. associated(handle % dimlist_tail)) then
         handle % dimlist_tail => new_dimlist_node
!call mpas_log_write('Assigning tail for '//trim(dimname))
      else
         handle % dimlist_tail % next => new_dimlist_node
         handle % dimlist_tail => handle % dimlist_tail % next
!call mpas_log_write('Extending tail for '//trim(dimname))
      end if

   end subroutine MPAS_io_def_dim


   subroutine MPAS_io_inq_var(handle, fieldname, fieldtype, ndims, dimnames, dimsizes, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      integer, intent(out), optional :: fieldtype
      integer, intent(out), optional :: ndims
      character (len=StrKIND), dimension(:), pointer, optional :: dimnames
      integer, dimension(:), pointer, optional :: dimsizes
      integer, intent(out), optional :: ierr

      integer :: i
      type (fieldlist_type), pointer :: new_fieldlist_node
      type (fieldlist_type), pointer :: field_cursor
      type (dimlist_type), pointer :: new_dimlist_node
      type (dimlist_type), pointer :: dim_cursor
      integer, dimension(:), pointer :: dimids
      logical :: found
      integer :: pio_ierr


!      call mpas_log_write('Called MPAS_io_inq_var()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      ! Sanity checks
      if (.not. handle % initialized) then
         if (present(ierr)) ierr = MPAS_IO_ERR_UNINIT_HANDLE
         return 
      end if


      !
      ! See if we already have this variable in our list
      !
      found = .false.
      field_cursor => handle % fieldlist_head
      do while (associated(field_cursor))
         if (trim(fieldname) == trim(field_cursor % fieldhandle % fieldname)) then
!call mpas_log_write('Already found variable in fieldlist')
            found = .true.
            exit
         end if
         field_cursor => field_cursor % next
      end do

      !
      ! Otherwise, inquire through the file-level API and add it to the list
      !
      if (.not. found) then

         allocate(new_fieldlist_node)
         nullify(new_fieldlist_node % next)
         allocate(new_fieldlist_node % fieldhandle)
      
         new_fieldlist_node % fieldhandle % fieldname = fieldname

         ! Get variable ID
         pio_ierr = PIO_inq_varid(handle % pio_file, trim(fieldname), new_fieldlist_node % fieldhandle % fieldid)
         pio_ierr = PIO_inq_varid(handle % pio_file, trim(fieldname), new_fieldlist_node % fieldhandle % field_desc)
         if (pio_ierr /= PIO_noerr) then
            if (present(ierr)) ierr = MPAS_IO_ERR_PIO
            deallocate(new_fieldlist_node % fieldhandle)
            deallocate(new_fieldlist_node)
            call mpas_log_write('Variable ' // trim(fieldname) // ' not in input file.', MPAS_LOG_WARN)
            return
         end if
!call mpas_log_write('Inquired about variable ID $i.', intArgs=(/new_fieldlist_node % fieldhandle % fieldid/) )

         ! Get field type
         pio_ierr = PIO_inq_vartype(handle % pio_file, new_fieldlist_node % fieldhandle % fieldid, new_fieldlist_node % fieldhandle % field_type)
         if (pio_ierr /= PIO_noerr) then
            if (present(ierr)) ierr = MPAS_IO_ERR_PIO
            deallocate(new_fieldlist_node % fieldhandle)
            deallocate(new_fieldlist_node)
            return
         end if
!call mpas_log_write('Inquired about variable type $i.', intArgs=(/new_fieldlist_node % fieldhandle % field_type/) )

         ! Convert to MPAS type
         new_fieldlist_node % fieldhandle % precision = MPAS_IO_NATIVE_PRECISION
         if (new_fieldlist_node % fieldhandle % field_type == PIO_double) then
            new_fieldlist_node % fieldhandle % field_type = MPAS_IO_DOUBLE
            new_fieldlist_node % fieldhandle % precision = MPAS_IO_DOUBLE_PRECISION
         else if (new_fieldlist_node % fieldhandle % field_type == PIO_real) then
            new_fieldlist_node % fieldhandle % field_type = MPAS_IO_REAL
            new_fieldlist_node % fieldhandle % precision = MPAS_IO_SINGLE_PRECISION
         else if (new_fieldlist_node % fieldhandle % field_type == PIO_int) then
            new_fieldlist_node % fieldhandle % field_type = MPAS_IO_INT
         else if (new_fieldlist_node % fieldhandle % field_type == PIO_char) then
            new_fieldlist_node % fieldhandle % field_type = MPAS_IO_CHAR
!!!!!!!! PIO DOES NOT SUPPORT LOGICAL !!!!!!!!
         end if

         ! Get number of dimensions
         pio_ierr = PIO_inq_varndims(handle % pio_file, new_fieldlist_node % fieldhandle % fieldid, new_fieldlist_node % fieldhandle % ndims)
         if (pio_ierr /= PIO_noerr) then
            if (present(ierr)) ierr = MPAS_IO_ERR_PIO
            deallocate(new_fieldlist_node % fieldhandle)
            deallocate(new_fieldlist_node)
            return
         end if
!call mpas_log_write('Inquired about number of dimensions $i', intArgs=(/new_fieldlist_node % fieldhandle % ndims/) )

         allocate(dimids(new_fieldlist_node % fieldhandle % ndims))

         ! Get dimension IDs
         if (new_fieldlist_node % fieldhandle % ndims > 0) then
            pio_ierr = PIO_inq_vardimid(handle % pio_file, new_fieldlist_node % fieldhandle % fieldid, dimids)
            if (pio_ierr /= PIO_noerr) then
               if (present(ierr)) ierr = MPAS_IO_ERR_PIO
               deallocate(new_fieldlist_node % fieldhandle)
               deallocate(new_fieldlist_node)
               deallocate(dimids)
               return
            end if
!call mpas_log_write('Inquired about dimension IDs $i.', intArgs=(/dimids/) )
         end if

         allocate(new_fieldlist_node % fieldhandle % dims(new_fieldlist_node % fieldhandle % ndims))

         ! Get information about dimensions
         do i=1,new_fieldlist_node % fieldhandle % ndims
            new_fieldlist_node % fieldhandle % dims(i) % dimid = dimids(i)
            if (dimids(i) == handle % pio_unlimited_dimid) then
               new_fieldlist_node % fieldhandle % dims(i) % is_unlimited_dim = .true.
               new_fieldlist_node % fieldhandle % has_unlimited_dim = .true.
            end if

            pio_ierr = PIO_inq_dimlen(handle % pio_file, dimids(i), new_fieldlist_node % fieldhandle % dims(i) % dimsize)
            if (pio_ierr /= PIO_noerr) then
               if (present(ierr)) ierr = MPAS_IO_ERR_PIO
               deallocate(new_fieldlist_node % fieldhandle)
               deallocate(new_fieldlist_node)
               deallocate(dimids)
               return
            end if
!call mpas_log_write('Inquired about dimension size $i', intArgs=(/new_fieldlist_node % fieldhandle % dims(i) % dimsize/)

            pio_ierr = PIO_inq_dimname(handle % pio_file, dimids(i), new_fieldlist_node % fieldhandle % dims(i) % dimname)
            if (pio_ierr /= PIO_noerr) then
               if (present(ierr)) ierr = MPAS_IO_ERR_PIO
               deallocate(new_fieldlist_node % fieldhandle)
               deallocate(new_fieldlist_node)
               deallocate(dimids)
               return
            end if
!call mpas_log_write('Inquired about dimension name ' // trim(new_fieldlist_node % fieldhandle % dims(i) % dimname))

         end do

         deallocate(dimids)

         ! Keep variable information for future reference
         if (.not. associated(handle % fieldlist_head)) then
            handle % fieldlist_head => new_fieldlist_node
!call mpas_log_write('Assigning head for '//trim(fieldname))
         end if
         if (.not. associated(handle % fieldlist_tail)) then
            handle % fieldlist_tail => new_fieldlist_node
!call mpas_log_write('Assigning tail for '//trim(fieldname))
         else
            handle % fieldlist_tail % next => new_fieldlist_node
            handle % fieldlist_tail => handle % fieldlist_tail % next
!call mpas_log_write('Extending tail for '//trim(fieldname))
         end if

         ! Keep dimension information for any new dimensions that were encountered
         do i=1,new_fieldlist_node % fieldhandle % ndims
            found = .false.
            dim_cursor => handle % dimlist_head
            do while (associated(dim_cursor))
               if (trim(dim_cursor % dimhandle % dimname) == trim(new_fieldlist_node % fieldhandle % dims(i) % dimname)) then
!call mpas_log_write('Already have dimension '//trim(new_fieldlist_node % fieldhandle % dims(i) % dimname)//' in our list...')
                  found = .true.
                  exit
               end if
               dim_cursor => dim_cursor % next
            end do

            if (.not. found) then
               allocate(new_dimlist_node)
               nullify(new_dimlist_node % next)
               allocate(new_dimlist_node % dimhandle)
               new_dimlist_node % dimhandle = new_fieldlist_node % fieldhandle % dims(i)
               if (.not. associated(handle % dimlist_head)) then
                  handle % dimlist_head => new_dimlist_node
!call mpas_log_write('Assigning head for '//trim(new_dimlist_node % dimhandle % dimname))
               end if
               if (.not. associated(handle % dimlist_tail)) then
                  handle % dimlist_tail => new_dimlist_node
!call mpas_log_write('Assigning tail for '//trim(new_dimlist_node % dimhandle % dimname))
               else
                  handle % dimlist_tail % next => new_dimlist_node
                  handle % dimlist_tail => handle % dimlist_tail % next
!call mpas_log_write('Extending tail for '//trim(new_dimlist_node % dimhandle % dimname))
               end if
            end if
         end do
         field_cursor => new_fieldlist_node
      end if


      !
      ! Set output arguments
      !
      if (present(fieldtype)) fieldtype = field_cursor % fieldhandle % field_type
      if (present(ndims)) ndims = field_cursor % fieldhandle % ndims
      if (present(dimnames)) then
         allocate(dimnames(field_cursor % fieldhandle % ndims))
         do i=1,field_cursor % fieldhandle % ndims
            dimnames(i) = field_cursor % fieldhandle % dims(i) % dimname
         end do
      end if
      if (present(dimsizes)) then
         allocate(dimsizes(field_cursor % fieldhandle % ndims))
         do i=1,field_cursor % fieldhandle % ndims
            dimsizes(i) = field_cursor % fieldhandle % dims(i) % dimsize
         end do
      end if

   end subroutine MPAS_io_inq_var


   subroutine MPAS_io_def_var(handle, fieldname, fieldtype, dimnames, precision, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      integer, intent(in) :: fieldtype
      character (len=StrKIND), dimension(:), intent(in) :: dimnames
      integer, intent(in), optional :: precision
      integer, intent(out), optional :: ierr

      integer :: i
      integer :: pio_ierr
      integer :: pio_type
      integer :: ndims
      integer :: inq_fieldtype
      integer :: inq_ndims
      character (len=StrKIND), dimension(:), pointer :: inq_dimnames
      type (fieldlist_type), pointer :: new_fieldlist_node
      type (fieldlist_type), pointer :: field_cursor
      type (dimlist_type), pointer :: dim_cursor
      integer, dimension(:), pointer :: dimids
      integer :: local_precision


!      call mpas_log_write('Called MPAS_io_def_var()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      ! Sanity checks
      if (.not. handle % initialized) then
         if (present(ierr)) ierr = MPAS_IO_ERR_UNINIT_HANDLE
         return 
      end if
      if (handle % data_mode) then
         if (present(ierr)) ierr = MPAS_IO_ERR_DATA_MODE
         return 
      end if
      if (handle % iomode /= MPAS_IO_WRITE) then
         if (present(ierr)) ierr = MPAS_IO_ERR_NOWRITE
         return 
      end if


      !
      ! Set precision to be use for reading/writing
      !
      if (present(precision)) then
         local_precision = precision
      else
         local_precision = MPAS_IO_NATIVE_PRECISION
      end if


      !
      ! Check whether this field has already been defined
      !
      ndims = size(dimnames)
!call mpas_log_write('Defining variable with $i dimensions', intArgs=(/ndims/))
      field_cursor => handle % fieldlist_head
      do while (associated(field_cursor))
         if (trim(fieldname) == trim(field_cursor % fieldhandle % fieldname)) then
            if (ndims /= field_cursor % fieldhandle % ndims) then
               if (present(ierr)) ierr = MPAS_IO_ERR_REDEF_VAR
!               call mpas_log_write('Field '//trim(fieldname)//' previously defined with conflicting number of dimensions: $i $i',&
!                           MPAS_LOG_ERR, intArgs=(/ndims, field_cursor % fieldhandle % ndims/) )
            end if
            if (fieldtype /= field_cursor % fieldhandle % field_type) then
               if (present(ierr)) ierr = MPAS_IO_ERR_REDEF_VAR
!               call mpas_log_write('Field '//trim(fieldname)//' previously defined with conflicting type: ', &
!                           MPAS_LOG_ERR, intArgs=(/fieldtype, field_cursor % fieldhandle % field_type/) )
            end if
            return
         end if
         field_cursor => field_cursor % next
      end do


      !
      ! If we're working with a preexisting file, see whether this variable isn't already
      !    defined in the file; if it is, the type and dimensions should match the type and dimensions
      !     specified in this call to define the variable, at which point we can add it to our local cache
      !
      if (handle % preexisting_file) then
         call MPAS_io_inq_var(handle, fieldname, inq_fieldtype, inq_ndims, inq_dimnames, ierr=pio_ierr)
         if (pio_ierr /= MPAS_IO_ERR_PIO) then

            ! Verify that the type and dimensions match...
            if (fieldtype == MPAS_IO_DOUBLE) then
               if (local_precision == MPAS_IO_SINGLE_PRECISION) then
                  pio_type = MPAS_IO_REAL
               else
                  pio_type = MPAS_IO_DOUBLE
               end if
            else
               pio_type = fieldtype
            end if
            if (pio_type /= inq_fieldtype) then
               if (present(ierr)) ierr = MPAS_IO_ERR_REDEF_DIM
               deallocate(inq_dimnames)
               return
            end if

            if (ndims /= inq_ndims) then
               if (present(ierr)) ierr = MPAS_IO_ERR_REDEF_DIM
               deallocate(inq_dimnames)
               return
            end if

            do i=1,ndims
               if (trim(dimnames(i)) /= trim(inq_dimnames(i))) then
                  if (present(ierr)) ierr = MPAS_IO_ERR_REDEF_DIM
                  deallocate(inq_dimnames)
                  return
               end if
            end do
   
            ! TODO: Can we get the dimension sizes to see whether they match those from the file?
   
         end if

         return
      end if


      !
      ! Otherwise, define it
      !
      allocate(new_fieldlist_node)
      nullify(new_fieldlist_node % next)
      allocate(new_fieldlist_node % fieldhandle)

      new_fieldlist_node % fieldhandle % fieldname = fieldname
      new_fieldlist_node % fieldhandle % field_type = fieldtype
      new_fieldlist_node % fieldhandle % ndims = ndims
      new_fieldlist_node % fieldhandle % precision = local_precision

      allocate(dimids(ndims))
      allocate(new_fieldlist_node % fieldhandle % dims(ndims))
      do i = 1, ndims
         dim_cursor => handle % dimlist_head
         do while (associated(dim_cursor))
            if (trim(dimnames(i)) == trim(dim_cursor % dimhandle % dimname)) then
               exit
            end if
            dim_cursor => dim_cursor % next
         end do
         if (associated(dim_cursor)) then
            dimids(i) = dim_cursor % dimhandle % dimid
            if (dim_cursor % dimhandle % is_unlimited_dim) new_fieldlist_node % fieldhandle % has_unlimited_dim = .true.
            new_fieldlist_node % fieldhandle % dims(i) = dim_cursor % dimhandle
!call mpas_log_write('Found dimension '//trim(new_fieldlist_node % fieldhandle % dims(i) % dimname)//' for field '//trim(fieldname))
         else
            if (present(ierr)) ierr = MPAS_IO_ERR_UNDEFINED_DIM
            deallocate(new_fieldlist_node % fieldhandle % dims)
            deallocate(new_fieldlist_node % fieldhandle)
            deallocate(new_fieldlist_node)
            deallocate(dimids)
            return
!            call mpas_log_write('Error finding dimension '//trim(dimnames(i))//' for field '//trim(fieldname), MPAS_LOG_ERR)
         end if
      end do

      ! Convert from MPAS type
      if (new_fieldlist_node % fieldhandle % field_type == MPAS_IO_DOUBLE) then
         if (local_precision == MPAS_IO_SINGLE_PRECISION) then
            pio_type = PIO_real
            new_fieldlist_node % fieldhandle % field_type = MPAS_IO_REAL
         else
            pio_type = PIO_double
         end if
      else if (new_fieldlist_node % fieldhandle % field_type == MPAS_IO_REAL) then
         if (local_precision == MPAS_IO_DOUBLE_PRECISION) then
            pio_type = PIO_double
            new_fieldlist_node % fieldhandle % field_type = MPAS_IO_DOUBLE
         else
            pio_type = PIO_real
         end if
      else if (new_fieldlist_node % fieldhandle % field_type == MPAS_IO_INT) then
         pio_type = PIO_int
      else if (new_fieldlist_node % fieldhandle % field_type == MPAS_IO_CHAR) then
         pio_type = PIO_char
!!!!!!!! PIO DOES NOT SUPPORT LOGICAL !!!!!!!!
      end if

      if (ndims == 0) then
         pio_ierr = PIO_def_var(handle % pio_file, trim(fieldname), pio_type, new_fieldlist_node % fieldhandle % field_desc)
      else
         pio_ierr = PIO_def_var(handle % pio_file, trim(fieldname), pio_type, dimids, new_fieldlist_node % fieldhandle % field_desc)
      end if
      if (pio_ierr /= PIO_noerr) then
         if (present(ierr)) ierr = MPAS_IO_ERR_PIO
         return
      end if

      ! Get the varid for use by put_att routines
      pio_ierr = PIO_inq_varid(handle % pio_file, trim(fieldname), new_fieldlist_node % fieldhandle % fieldid)
      if (pio_ierr /= PIO_noerr) then
         if (present(ierr)) ierr = MPAS_IO_ERR_PIO
         return
      end if

      deallocate(dimids)

      ! Keep variable information for future use
      if (.not. associated(handle % fieldlist_head)) then
         handle % fieldlist_head => new_fieldlist_node
!call mpas_log_write('Assigning head for '//trim(fieldname))
      end if
      if (.not. associated(handle % fieldlist_tail)) then
         handle % fieldlist_tail => new_fieldlist_node
!call mpas_log_write('Assigning tail for '//trim(fieldname))
      else
         handle % fieldlist_tail % next => new_fieldlist_node
         handle % fieldlist_tail => handle % fieldlist_tail % next
!call mpas_log_write('Extending tail for '//trim(fieldname))
      end if

   end subroutine MPAS_io_def_var


   subroutine MPAS_io_get_var_indices(handle, fieldname, indices, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(in) :: handle
      character (len=*), intent(in) :: fieldname
      integer, dimension(:), pointer :: indices
      integer, intent(out), optional :: ierr

      type (fieldlist_type), pointer :: field_cursor

!      call mpas_log_write('Called MPAS_io_get_var_indices()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      ! Sanity checks
      if (.not. handle % initialized) then
         if (present(ierr)) ierr = MPAS_IO_ERR_UNINIT_HANDLE
         return 
      end if


      !  
      ! Check whether the field has been defined
      !
      field_cursor => handle % fieldlist_head
      do while (associated(field_cursor))
         if (trim(fieldname) == trim(field_cursor % fieldhandle % fieldname)) then
            exit
         end if
         field_cursor => field_cursor % next
      end do
      if (.not. associated(field_cursor)) then
         if (present(ierr)) ierr = MPAS_IO_ERR_UNDEFINED_VAR
         return
      end if
!call mpas_log_write(trim(fieldname)//' has been defined')

      if (.not. associated(field_cursor % fieldhandle % decomp)) then
         if (present(ierr)) ierr = MPAS_IO_ERR_NO_DECOMP
         return
      end if

      allocate(indices(size(field_cursor % fieldhandle % decomp % indices)))
      indices(:) = field_cursor % fieldhandle % decomp % indices(:)

   end subroutine MPAS_io_get_var_indices


   subroutine MPAS_io_set_var_indices(handle, fieldname, indices, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      integer, dimension(:), intent(in) :: indices
      integer, intent(out), optional :: ierr

      type (fieldlist_type), pointer :: field_cursor
      integer :: pio_type
      integer :: ndims
      integer (kind=MPAS_IO_OFFSET_KIND) :: pd, indx
      integer :: i 
      integer :: early_return, early_return_global
      integer (kind=MPAS_IO_OFFSET_KIND) :: i1, i2, i3, i4, i5
      integer, dimension(:), pointer :: dimlist
      integer (kind=MPAS_IO_OFFSET_KIND), dimension(:), pointer :: compdof
      type (decomplist_type), pointer :: decomp_cursor, new_decomp

!      call mpas_log_write('Called MPAS_io_set_var_indices()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      ! Sanity checks
      if (.not. handle % initialized) then
         if (present(ierr)) ierr = MPAS_IO_ERR_UNINIT_HANDLE
         return 
      end if

!      call mpas_log_write('Assigning $i indices for '//trim(fieldname), intArgs=(/size(indices)/) )
      !  
      ! Check whether the field has been defined
      !
      field_cursor => handle % fieldlist_head
      do while (associated(field_cursor))
         if (trim(fieldname) == trim(field_cursor % fieldhandle % fieldname)) then
            exit
         end if
         field_cursor => field_cursor % next
      end do
      if (.not. associated(field_cursor)) then
         if (present(ierr)) ierr = MPAS_IO_ERR_UNDEFINED_VAR
         return
      end if
!call mpas_log_write(trim(fieldname)//' has been defined')

      !
      ! If this is a scalar field, just return
      !
      if (field_cursor % fieldhandle % ndims == 0 .or. &
          (field_cursor % fieldhandle % ndims == 1 .and. field_cursor % fieldhandle % has_unlimited_dim) .or. &
          field_cursor % fieldhandle % field_type == MPAS_IO_CHAR) then
!call mpas_log_write('No need to create a decomposition for a 0d field...')
         return
      end if


      !
      ! Check whether a suitable decomposition already exists
      !
      decomp_cursor => handle % ioContext % decomp_list
!if (.not. associated(decomp_cursor)) call mpas_log_write('No existing decompositions to check...')
      early_return = 0
      DECOMP_LOOP: do while (associated(decomp_cursor))
         if (decomp_cursor % decomphandle % field_type == field_cursor % fieldhandle % field_type) then
         if (size(decomp_cursor % decomphandle % dims) == field_cursor % fieldhandle % ndims) then
!call mpas_log_write('Number of dimensions matches...')
            do i=1,field_cursor % fieldhandle % ndims
!call mpas_log_write('Checking dimension ', intArgs=(/decomp_cursor % decomphandle % dims(i), field_cursor % fieldhandle % dims(i) % dimsize/) )
               if (decomp_cursor % decomphandle % dims(i) /= field_cursor % fieldhandle % dims(i) % dimsize) then
                  decomp_cursor => decomp_cursor % next
                  cycle DECOMP_LOOP
               end if
            end do

            if (size(decomp_cursor % decomphandle % indices) /= size(indices)) then
!call mpas_log_write('We do not have the same number of indices in this decomposition...')
               decomp_cursor => decomp_cursor % next
               cycle DECOMP_LOOP
            end if

            do i=1,size(decomp_cursor % decomphandle % indices)
               if (indices(i) /= decomp_cursor % decomphandle % indices(i)) then
!call mpas_log_write('One of the indices does not match... $i', intArgs=(/i/) )
                  decomp_cursor => decomp_cursor % next
                  cycle DECOMP_LOOP
               end if
            end do 
            
            ! OK, we have a match... just use this decomposition for the field and return
            field_cursor % fieldhandle % decomp => decomp_cursor % decomphandle 
!call mpas_log_write('Found a matching decomposition that we can use')
            early_return = 1
            exit DECOMP_LOOP
         else if ((size(decomp_cursor % decomphandle % dims) == field_cursor % fieldhandle % ndims - 1)  &
                  .and. field_cursor % fieldhandle % has_unlimited_dim  &
                 ) then
!call mpas_log_write('Number of non-record dimensions matches...')
            do i=1,field_cursor % fieldhandle % ndims
               if (field_cursor % fieldhandle % dims(i) % is_unlimited_dim) cycle
!call mpas_log_write('Checking dimension ', intArgs=(/decomp_cursor % decomphandle % dims(i), field_cursor % fieldhandle % dims(i) % dimsize/) )
               if (decomp_cursor % decomphandle % dims(i) /= field_cursor % fieldhandle % dims(i) % dimsize) then
                  decomp_cursor => decomp_cursor % next
                  cycle DECOMP_LOOP
               end if
            end do

            ! Type and dimensions match... what about indices?
            if (size(decomp_cursor % decomphandle % indices) /= size(indices)) then
!call mpas_log_write('We do not have the same number of indices in this decomposition...')
               decomp_cursor => decomp_cursor % next
               cycle DECOMP_LOOP
            end if

            do i=1,size(decomp_cursor % decomphandle % indices)
               if (indices(i) /= decomp_cursor % decomphandle % indices(i)) then
!call mpas_log_write('One of the indices does not match... $i', intArgs=(/i/) )
                  decomp_cursor => decomp_cursor % next
                  cycle DECOMP_LOOP
               end if
            end do 
            
            ! OK, we have a match... just use this decomposition for the field and return
            field_cursor % fieldhandle % decomp => decomp_cursor % decomphandle 
!call mpas_log_write('Found a matching decomposition that we can use (aside from record dimension)')
            early_return = 1
            exit DECOMP_LOOP
         end if
         end if
         decomp_cursor => decomp_cursor % next
      end do DECOMP_LOOP

      ! 
      ! If all tasks have set early_return to 1, then we have a usable decomposition and can return
      ! 
      call mpas_dmpar_min_int(handle % ioContext % dminfo, early_return, early_return_global)
      if (early_return_global == 1) then
         return
      end if


!call mpas_log_write('Creating a new decomposition')


      !
      ! Otherwise, we need to create a new decomposition
      !
      ndims = field_cursor % fieldhandle % ndims
      if (field_cursor % fieldhandle % has_unlimited_dim) ndims = ndims - 1
      

      allocate(new_decomp)
      nullify(new_decomp % next)
      allocate(new_decomp % decomphandle)
      allocate(new_decomp % decomphandle % dims(ndims))
      allocate(new_decomp % decomphandle % indices(size(indices)))

      new_decomp % decomphandle % field_type = field_cursor % fieldhandle % field_type
      new_decomp % decomphandle % indices(:) = indices(:)

      ! Convert from MPAS type
      if (field_cursor % fieldhandle % field_type == MPAS_IO_DOUBLE) then
         pio_type = PIO_double
      else if (field_cursor % fieldhandle % field_type == MPAS_IO_REAL) then
         pio_type = PIO_real
      else if (field_cursor % fieldhandle % field_type == MPAS_IO_INT) then
         pio_type = PIO_int
      else if (field_cursor % fieldhandle % field_type == MPAS_IO_CHAR) then
         pio_type = PIO_char
 !!!!!!! PIO DOES NOT SUPPORT LOGICAL !!!!!!!!
      end if

      allocate(dimlist(ndims))

      pd = 1
      do i=1,ndims-1
         dimlist(i) = field_cursor % fieldhandle % dims(i) % dimsize
         new_decomp % decomphandle % dims(i) = dimlist(i)
         pd = pd * int(dimlist(i),MPAS_IO_OFFSET_KIND)
      end do
      new_decomp % decomphandle % dims(ndims) = field_cursor % fieldhandle % dims(ndims) % dimsize
      dimlist(ndims) = size(indices)
      pd = pd * int(dimlist(ndims),MPAS_IO_OFFSET_KIND)

      allocate(compdof(pd)) 

      indx = 1
      if (ndims == 5) then
         do i5=1,dimlist(5)
         do i4=1,dimlist(4)
         do i3=1,dimlist(3)
         do i2=1,dimlist(2)
         do i1=1,dimlist(1)
            compdof(indx) = i1 &
                          + (i2-1)*int(dimlist(1),MPAS_IO_OFFSET_KIND) &
                          + (i3-1)*int(dimlist(2),MPAS_IO_OFFSET_KIND)*int(dimlist(1),MPAS_IO_OFFSET_KIND) &
                          + (i4-1)*int(dimlist(3),MPAS_IO_OFFSET_KIND)*int(dimlist(2),MPAS_IO_OFFSET_KIND)*int(dimlist(1),MPAS_IO_OFFSET_KIND) &
                          + int(indices(i5)-1,MPAS_IO_OFFSET_KIND)*int(dimlist(4),MPAS_IO_OFFSET_KIND)*int(dimlist(3),MPAS_IO_OFFSET_KIND)*int(dimlist(2),MPAS_IO_OFFSET_KIND)*int(dimlist(1),MPAS_IO_OFFSET_KIND)
            indx = indx + 1
         end do
         end do
         end do
         end do
         end do
      else if (ndims == 4) then
         do i4=1,dimlist(4)
         do i3=1,dimlist(3)
         do i2=1,dimlist(2)
         do i1=1,dimlist(1)
            compdof(indx) = i1 &
                          + (i2-1)*int(dimlist(1),MPAS_IO_OFFSET_KIND) &
                          + (i3-1)*int(dimlist(2),MPAS_IO_OFFSET_KIND)*int(dimlist(1),MPAS_IO_OFFSET_KIND) &
                          + int(indices(i4)-1,MPAS_IO_OFFSET_KIND)*int(dimlist(3),MPAS_IO_OFFSET_KIND)*int(dimlist(2),MPAS_IO_OFFSET_KIND)*int(dimlist(1),MPAS_IO_OFFSET_KIND)
            indx = indx + 1
         end do
         end do
         end do
         end do
      else if (ndims == 3) then
         do i3=1,dimlist(3)
         do i2=1,dimlist(2)
         do i1=1,dimlist(1)
            compdof(indx) = i1 + (i2-1)*int(dimlist(1),MPAS_IO_OFFSET_KIND) + int(indices(i3)-1,MPAS_IO_OFFSET_KIND)*int(dimlist(2),MPAS_IO_OFFSET_KIND)*int(dimlist(1),MPAS_IO_OFFSET_KIND)
            indx = indx + 1
         end do
         end do
         end do
      else if (ndims == 2) then
         do i2=1,dimlist(2)
         do i1=1,dimlist(1)
            compdof(indx) = i1 + int(indices(i2)-1,MPAS_IO_OFFSET_KIND)*int(dimlist(1),MPAS_IO_OFFSET_KIND)
            indx = indx + 1
         end do
         end do
      else if (ndims == 1) then
         do i1=1,dimlist(1)
            compdof(indx) = int(indices(i1),MPAS_IO_OFFSET_KIND)
            indx = indx + 1
         end do
      end if

      dimlist(ndims) = field_cursor % fieldhandle % dims(ndims) % dimsize
      call PIO_initdecomp(handle % ioContext % pio_iosystem, pio_type, dimlist, compdof, new_decomp % decomphandle % pio_iodesc)

      ! Add new decomposition to the list
      if (.not. associated(handle % ioContext % decomp_list)) then
         handle % ioContext % decomp_list => new_decomp
!call mpas_log_write('Adding first item to the decomp_list')
      else
         new_decomp % next => handle % ioContext % decomp_list
         handle % ioContext % decomp_list => new_decomp
!call mpas_log_write('Adding new decomp to the head of the list')
      end if

!call mpas_log_write('Setting decomp in fieldhandle')
      field_cursor % fieldhandle % decomp => new_decomp % decomphandle

      deallocate(compdof)
      deallocate(dimlist)
!call mpas_log_write('All finished.')

   end subroutine MPAS_io_set_var_indices


   subroutine MPAS_io_get_var_generic(handle, fieldname, intVal, intArray1d, intArray2d, intArray3d, intArray4d, &
                                                        realVal, realArray1d, realArray2d, realArray3d, realArray4d, realArray5d, &
                                                        charVal, charArray1d, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      integer, intent(out), optional :: intVal
      integer, dimension(:), intent(out), optional :: intArray1d
      integer, dimension(:,:), intent(out), optional :: intArray2d
      integer, dimension(:,:,:), intent(out), optional :: intArray3d
      integer, dimension(:,:,:,:), intent(out), optional :: intArray4d
      real (kind=RKIND), intent(out), optional :: realVal
      real (kind=RKIND), dimension(:), intent(out), optional :: realArray1d
      real (kind=RKIND), dimension(:,:), intent(out), optional :: realArray2d
      real (kind=RKIND), dimension(:,:,:), intent(out), optional :: realArray3d
      real (kind=RKIND), dimension(:,:,:,:), intent(out), optional :: realArray4d
      real (kind=RKIND), dimension(:,:,:,:,:), intent(out), optional :: realArray5d
      character (len=*), intent(out), optional :: charVal
      character (len=*), dimension(:), intent(out), optional :: charArray1d
      integer, intent(out), optional :: ierr

      integer :: pio_ierr
      integer, dimension(1) :: start1
      integer, dimension(1) :: count1
      integer, dimension(2) :: start2
      integer, dimension(2) :: count2
      integer, dimension(3) :: start3
      integer, dimension(3) :: count3
      integer, dimension(4) :: start4
      integer, dimension(4) :: count4
      integer, dimension(5) :: start5
      integer, dimension(5) :: count5
      integer, dimension(6) :: start6
      integer, dimension(6) :: count6
      character (len=StrKIND), dimension(1) :: tempchar
      type (fieldlist_type), pointer :: field_cursor

      integer i, j

      real (kind=R4KIND) :: singleVal
      real (kind=R4KIND), dimension(:), allocatable :: singleArray1d
      real (kind=R4KIND), dimension(:,:), allocatable :: singleArray2d
      real (kind=R4KIND), dimension(:,:,:), allocatable :: singleArray3d
      real (kind=R4KIND), dimension(:,:,:,:), allocatable :: singleArray4d
      real (kind=R4KIND), dimension(:,:,:,:,:), allocatable :: singleArray5d

      real (kind=R8KIND) :: doubleVal
      real (kind=R8KIND), dimension(:), allocatable :: doubleArray1d
      real (kind=R8KIND), dimension(:,:), allocatable :: doubleArray2d
      real (kind=R8KIND), dimension(:,:,:), allocatable :: doubleArray3d
      real (kind=R8KIND), dimension(:,:,:,:), allocatable :: doubleArray4d
      real (kind=R8KIND), dimension(:,:,:,:,:), allocatable :: doubleArray5d

      ! Sanity checks
      if (.not. handle % initialized) then
         if (present(ierr)) ierr = MPAS_IO_ERR_UNINIT_HANDLE
         return 
      end if

!      call mpas_log_write('Reading '// trim(fieldname))

      !
      ! Check whether the field has been defined
      !
!      call mpas_log_write('Checking if field is define')
      field_cursor => handle % fieldlist_head
      do while (associated(field_cursor))
         if (trim(fieldname) == trim(field_cursor % fieldhandle % fieldname)) then
            exit
         end if
         field_cursor => field_cursor % next
      end do
      if (.not. associated(field_cursor)) then
         if (present(ierr)) ierr = MPAS_IO_ERR_UNDEFINED_VAR
         return
      end if

!!!! Assume array was already allocated by the user

!      call mpas_log_write('Checking for unlimited dim')
      if (field_cursor % fieldhandle % has_unlimited_dim) then
#ifdef USE_PIO2
         call PIO_setframe(handle % pio_file, field_cursor % fieldhandle % field_desc, handle % frame_number)
#else
         call PIO_setframe(field_cursor % fieldhandle % field_desc, handle % frame_number)
#endif
         start1(1) = handle % frame_number
         count1(1) = 1
     
         start2(1) = 1
         start2(2) = handle % frame_number
         count2(2) = 1
      end if

!      call mpas_log_write('Checking for real, int, char, etc')
      if (present(realVal)) then
!         call mpas_log_write('  value is real')
         if ((field_cursor % fieldhandle % precision == MPAS_IO_SINGLE_PRECISION) .and. &
             (MPAS_IO_NATIVE_PRECISION /= MPAS_IO_SINGLE_PRECISION)) then
            if (field_cursor % fieldhandle % has_unlimited_dim) then
               pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start1, singleVal)
            else
               pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % field_desc, singleVal)
            end if
            realVal = real(singleVal,RKIND)
         else if ((field_cursor % fieldhandle % precision == MPAS_IO_DOUBLE_PRECISION) .and. &
             (MPAS_IO_NATIVE_PRECISION /= MPAS_IO_DOUBLE_PRECISION)) then
            if (field_cursor % fieldhandle % has_unlimited_dim) then
               pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start1, doubleVal)
            else
               pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % field_desc, doubleVal)
            end if
            realVal = real(doubleVal,RKIND)
         else
            if (field_cursor % fieldhandle % has_unlimited_dim) then
               pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start1, realVal)
            else
               pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % field_desc, realVal)
            end if
         end if
      else if (present(intVal)) then
!         call mpas_log_write('  value is int')
         if (field_cursor % fieldhandle % has_unlimited_dim) then
            pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start1, intVal)
         else
            pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % field_desc, intVal)
         end if
      else if (present(charVal)) then
!         call mpas_log_write('  value is char')
         if (field_cursor % fieldhandle % has_unlimited_dim) then
            count2(1) = field_cursor % fieldhandle % dims(1) % dimsize
            pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start2, count2, tempchar)
            charVal(1:count2(1)) = tempchar(1)(1:count2(1))
         else
            start1(1) = 1
            count1(1) = field_cursor % fieldhandle % dims(1) % dimsize
            pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start1, count1, tempchar)
            charVal(1:count1(1)) = tempchar(1)(1:count1(1))
         end if
      else if (present(charArray1d)) then
!         call mpas_log_write('  value is char1')
         if (field_cursor % fieldhandle % has_unlimited_dim) then
            ! Can only read one string at a time, since the sizes differ so much (i.e. StrLen != StrKIND)
            do i = 1, field_cursor % fieldhandle % dims(2) % dimsize
               start3(1) = 1
               start3(2) = i
               start3(3) = handle % frame_number
               count3(1) = field_cursor % fieldhandle % dims(1) % dimsize ! Should be StrLen
               count3(2) = 1
               count3(3) = 1
               pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start3, count3, tempchar)

               ! Copy all characters up to the first C null char or the end of
               ! the string, whichever comes first
               charArray1d(i)(:) = ' '
               do j = 1, len(tempchar(1))
                  if ( tempchar(1)(j:j) /= CHAR(0)) then
                     charArray1d(i)(j:j) = tempchar(1)(j:j)
                  else
                     exit
                  end if
               end do
            end do
         else
            ! Can only read one string at a time, since the sizes differ so much (i.e. StrLen != StrKIND)
            do i = 1, field_cursor % fieldhandle % dims(2) % dimsize
               start2(1) = 1
               start2(2) = i
               count2(1) = field_cursor % fieldhandle % dims(1) % dimsize ! Should be StrLen
               count2(2) = 1
               pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start2, count2, tempchar)

               ! Copy all characters up to the first C null char or the end of
               ! the string, whichever comes first
               charArray1d(i)(:) = ' '
               do j = 1, len(tempchar(1))
                  if ( tempchar(1)(j:j) /= CHAR(0)) then
                     charArray1d(i)(j:j) = tempchar(1)(j:j)
                  else
                     exit
                  end if
               end do
            end do
         end if
      else if (present(realArray1d)) then
!         call mpas_log_write('  value is real1')
         if ((field_cursor % fieldhandle % precision == MPAS_IO_SINGLE_PRECISION) .and. &
             (MPAS_IO_NATIVE_PRECISION /= MPAS_IO_SINGLE_PRECISION)) then
            allocate(singleArray1d(size(realArray1d,1)))
            if (associated(field_cursor % fieldhandle % decomp)) then
               call PIO_read_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                    singleArray1d, pio_ierr)
            else
               if (field_cursor % fieldhandle % has_unlimited_dim) then
                  start2(1) = 1
                  start2(2) = handle % frame_number
                  count2(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count2(2) = 1
                  pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start2, count2, singleArray1d)
               else
                  start1(:) = 1
                  count1(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start1, count1, singleArray1d)
               end if
            end if
            realArray1d(:) = real(singleArray1d(:),RKIND)
            deallocate(singleArray1d)
         else if ((field_cursor % fieldhandle % precision == MPAS_IO_DOUBLE_PRECISION) .and. &
             (MPAS_IO_NATIVE_PRECISION /= MPAS_IO_DOUBLE_PRECISION)) then
            allocate(doubleArray1d(size(realArray1d,1)))
            if (associated(field_cursor % fieldhandle % decomp)) then
               call PIO_read_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                    doubleArray1d, pio_ierr)
            else
               if (field_cursor % fieldhandle % has_unlimited_dim) then
                  start2(1) = 1
                  start2(2) = handle % frame_number
                  count2(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count2(2) = 1
                  pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start2, count2, doubleArray1d)
               else
                  start1(:) = 1
                  count1(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start1, count1, doubleArray1d)
               end if
            end if
            realArray1d(:) = real(doubleArray1d(:),RKIND)
            deallocate(doubleArray1d)
         else
            if (associated(field_cursor % fieldhandle % decomp)) then
               call PIO_read_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                    realArray1d, pio_ierr)
            else
               if (field_cursor % fieldhandle % has_unlimited_dim) then
                  start2(1) = 1
                  start2(2) = handle % frame_number
                  count2(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count2(2) = 1
                  pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start2, count2, realArray1d)
               else
                  start1(:) = 1
                  count1(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start1, count1, realArray1d)
               end if
            end if
         end if
      else if (present(realArray2d)) then
!         call mpas_log_write('  value is real2')
         if ((field_cursor % fieldhandle % precision == MPAS_IO_SINGLE_PRECISION) .and. &
             (MPAS_IO_NATIVE_PRECISION /= MPAS_IO_SINGLE_PRECISION)) then
            allocate(singleArray2d(size(realArray2d,1), size(realArray2d,2)))
            if (associated(field_cursor % fieldhandle % decomp)) then
               call PIO_read_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                    singleArray2d, pio_ierr)
            else
               if (field_cursor % fieldhandle % has_unlimited_dim) then
                  start3(:) = 1
                  start3(3) = handle % frame_number
                  count3(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count3(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count3(3) = 1
                  pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start3, count3, singleArray2d)
               else
                  start2(:) = 1
                  count2(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count2(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start2, count2, singleArray2d)
               end if
            end if
            realArray2d(:,:) = real(singleArray2d(:,:),RKIND)
            deallocate(singleArray2d)
         else if ((field_cursor % fieldhandle % precision == MPAS_IO_DOUBLE_PRECISION) .and. &
             (MPAS_IO_NATIVE_PRECISION /= MPAS_IO_DOUBLE_PRECISION)) then
            allocate(doubleArray2d(size(realArray2d,1), size(realArray2d,2)))
            if (associated(field_cursor % fieldhandle % decomp)) then
               call PIO_read_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                    doubleArray2d, pio_ierr)
            else
               if (field_cursor % fieldhandle % has_unlimited_dim) then
                  start3(:) = 1
                  start3(3) = handle % frame_number
                  count3(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count3(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count3(3) = 1
                  pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start3, count3, doubleArray2d)
               else
                  start2(:) = 1
                  count2(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count2(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start2, count2, doubleArray2d)
               end if
            end if
            realArray2d(:,:) = real(doubleArray2d(:,:),RKIND)
            deallocate(doubleArray2d)
         else
            if (associated(field_cursor % fieldhandle % decomp)) then
               call PIO_read_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                    realArray2d, pio_ierr)
            else
               if (field_cursor % fieldhandle % has_unlimited_dim) then
                  start3(:) = 1
                  start3(3) = handle % frame_number
                  count3(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count3(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count3(3) = 1
                  pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start3, count3, realArray2d)
               else
                  start2(:) = 1
                  count2(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count2(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start2, count2, realArray2d)
               end if
            end if
         end if
      else if (present(realArray3d)) then
!         call mpas_log_write('  value is real3')
         if ((field_cursor % fieldhandle % precision == MPAS_IO_SINGLE_PRECISION) .and. &
             (MPAS_IO_NATIVE_PRECISION /= MPAS_IO_SINGLE_PRECISION)) then
            allocate(singleArray3d(size(realArray3d,1),size(realArray3d,2),size(realArray3d,3)))
            if (associated(field_cursor % fieldhandle % decomp)) then
               call PIO_read_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                    singleArray3d, pio_ierr)
            else
               if (field_cursor % fieldhandle % has_unlimited_dim) then
                  start4(:) = 1
                  start4(4) = handle % frame_number
                  count4(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count4(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count4(3) = field_cursor % fieldhandle % dims(3) % dimsize
                  count4(4) = 1
                  pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start4, count4, singleArray3d)
               else
                  start3(:) = 1
                  count3(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count3(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count3(3) = field_cursor % fieldhandle % dims(3) % dimsize
                  pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start3, count3, singleArray3d)
               end if
            end if
            realArray3d(:,:,:) = real(singleArray3d(:,:,:),RKIND)
            deallocate(singleArray3d)
         else if ((field_cursor % fieldhandle % precision == MPAS_IO_DOUBLE_PRECISION) .and. &
             (MPAS_IO_NATIVE_PRECISION /= MPAS_IO_DOUBLE_PRECISION)) then
            allocate(doubleArray3d(size(realArray3d,1),size(realArray3d,2),size(realArray3d,3)))
            if (associated(field_cursor % fieldhandle % decomp)) then
               call PIO_read_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                    doubleArray3d, pio_ierr)
            else
               if (field_cursor % fieldhandle % has_unlimited_dim) then
                  start4(:) = 1
                  start4(4) = handle % frame_number
                  count4(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count4(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count4(3) = field_cursor % fieldhandle % dims(3) % dimsize
                  count4(4) = 1
                  pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start4, count4, doubleArray3d)
               else
                  start3(:) = 1
                  count3(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count3(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count3(3) = field_cursor % fieldhandle % dims(3) % dimsize
                  pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start3, count3, doubleArray3d)
               end if
            end if
            realArray3d(:,:,:) = real(doubleArray3d(:,:,:),RKIND)
            deallocate(doubleArray3d)
         else
            if (associated(field_cursor % fieldhandle % decomp)) then
               call PIO_read_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                    realArray3d, pio_ierr)
            else
               if (field_cursor % fieldhandle % has_unlimited_dim) then
                  start4(:) = 1
                  start4(4) = handle % frame_number
                  count4(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count4(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count4(3) = field_cursor % fieldhandle % dims(3) % dimsize
                  count4(4) = 1
                  pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start4, count4, realArray3d)
               else
                  start3(:) = 1
                  count3(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count3(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count3(3) = field_cursor % fieldhandle % dims(3) % dimsize
                  pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start3, count3, realArray3d)
               end if
            end if
         end if
      else if (present(realArray4d)) then
!         call mpas_log_write('  value is real4')
         if ((field_cursor % fieldhandle % precision == MPAS_IO_SINGLE_PRECISION) .and. &
             (MPAS_IO_NATIVE_PRECISION /= MPAS_IO_SINGLE_PRECISION)) then
            allocate(singleArray4d(size(realArray4d,1),size(realArray4d,2),size(realArray4d,3),size(realArray4d,4)))
            if (associated(field_cursor % fieldhandle % decomp)) then
               call PIO_read_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                    singleArray4d, pio_ierr)
            else
               if (field_cursor % fieldhandle % has_unlimited_dim) then
                  start5(:) = 1
                  start5(5) = handle % frame_number
                  count5(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count5(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count5(3) = field_cursor % fieldhandle % dims(3) % dimsize
                  count5(4) = field_cursor % fieldhandle % dims(4) % dimsize
                  count5(5) = 1
                  pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start5, count5, singleArray4d)
               else
                  start4(:) = 1
                  count4(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count4(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count4(3) = field_cursor % fieldhandle % dims(3) % dimsize
                  count4(4) = field_cursor % fieldhandle % dims(4) % dimsize
                  pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start4, count4, singleArray4d)
               end if
            end if
            realArray4d(:,:,:,:) = real(singleArray4d(:,:,:,:),RKIND)
            deallocate(singleArray4d)
         else if ((field_cursor % fieldhandle % precision == MPAS_IO_DOUBLE_PRECISION) .and. &
             (MPAS_IO_NATIVE_PRECISION /= MPAS_IO_DOUBLE_PRECISION)) then
            allocate(doubleArray4d(size(realArray4d,1),size(realArray4d,2),size(realArray4d,3),size(realArray4d,4)))
            if (associated(field_cursor % fieldhandle % decomp)) then
               call PIO_read_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                    doubleArray4d, pio_ierr)
            else
               if (field_cursor % fieldhandle % has_unlimited_dim) then
                  start5(:) = 1
                  start5(5) = handle % frame_number
                  count5(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count5(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count5(3) = field_cursor % fieldhandle % dims(3) % dimsize
                  count5(4) = field_cursor % fieldhandle % dims(4) % dimsize
                  count5(5) = 1
                  pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start5, count5, doubleArray4d)
               else
                  start4(:) = 1
                  count4(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count4(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count4(3) = field_cursor % fieldhandle % dims(3) % dimsize
                  count4(4) = field_cursor % fieldhandle % dims(4) % dimsize
                  pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start4, count4, doubleArray4d)
               end if
            end if
            realArray4d(:,:,:,:) = real(doubleArray4d(:,:,:,:),RKIND)
            deallocate(doubleArray4d)
         else
            if (associated(field_cursor % fieldhandle % decomp)) then
               call PIO_read_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                    realArray4d, pio_ierr)
            else
               if (field_cursor % fieldhandle % has_unlimited_dim) then
                  start5(:) = 1
                  start5(5) = handle % frame_number
                  count5(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count5(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count5(3) = field_cursor % fieldhandle % dims(3) % dimsize
                  count5(4) = field_cursor % fieldhandle % dims(4) % dimsize
                  count5(5) = 1
                  pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start5, count5, realArray4d)
               else
                  start4(:) = 1
                  count4(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count4(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count4(3) = field_cursor % fieldhandle % dims(3) % dimsize
                  count4(4) = field_cursor % fieldhandle % dims(4) % dimsize
                  pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start4, count4, realArray4d)
               end if
            end if
         end if
      else if (present(realArray5d)) then
!         call mpas_log_write('  value is real5')
         if ((field_cursor % fieldhandle % precision == MPAS_IO_SINGLE_PRECISION) .and. &
             (MPAS_IO_NATIVE_PRECISION /= MPAS_IO_SINGLE_PRECISION)) then
            allocate(singleArray5d(size(realArray5d,1),size(realArray5d,2),size(realArray5d,3),size(realArray5d,4),size(realArray5d,5)))
            if (associated(field_cursor % fieldhandle % decomp)) then
               call PIO_read_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                   singleArray5d, pio_ierr)
            else
               if (field_cursor % fieldhandle % has_unlimited_dim) then
                  start6(:) = 1
                  start6(6) = handle % frame_number
                  count6(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count6(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count6(3) = field_cursor % fieldhandle % dims(3) % dimsize
                  count6(4) = field_cursor % fieldhandle % dims(4) % dimsize
                  count6(5) = field_cursor % fieldhandle % dims(5) % dimsize
                  count6(6) = 1
                  pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start6, count6, singleArray5d)
               else
                  start5(:) = 1
                  count5(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count5(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count5(3) = field_cursor % fieldhandle % dims(3) % dimsize
                  count5(4) = field_cursor % fieldhandle % dims(4) % dimsize
                  count5(5) = field_cursor % fieldhandle % dims(5) % dimsize
                  pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start5, count5, singleArray5d)
               end if
            end if
            realArray5d(:,:,:,:,:) = real(singleArray5d(:,:,:,:,:),RKIND)
            deallocate(singleArray5d)
         else if ((field_cursor % fieldhandle % precision == MPAS_IO_DOUBLE_PRECISION) .and. &
             (MPAS_IO_NATIVE_PRECISION /= MPAS_IO_DOUBLE_PRECISION)) then
            allocate(doubleArray5d(size(realArray5d,1),size(realArray5d,2),size(realArray5d,3),size(realArray5d,4),size(realArray5d,5)))
            if (associated(field_cursor % fieldhandle % decomp)) then
               call PIO_read_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                   doubleArray5d, pio_ierr)
            else
               if (field_cursor % fieldhandle % has_unlimited_dim) then
                  start6(:) = 1
                  start6(6) = handle % frame_number
                  count6(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count6(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count6(3) = field_cursor % fieldhandle % dims(3) % dimsize
                  count6(4) = field_cursor % fieldhandle % dims(4) % dimsize
                  count6(5) = field_cursor % fieldhandle % dims(5) % dimsize
                  count6(6) = 1
                  pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start6, count6, doubleArray5d)
               else
                  start5(:) = 1
                  count5(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count5(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count5(3) = field_cursor % fieldhandle % dims(3) % dimsize
                  count5(4) = field_cursor % fieldhandle % dims(4) % dimsize
                  count5(5) = field_cursor % fieldhandle % dims(5) % dimsize
                  pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start5, count5, doubleArray5d)
               end if
            end if
            realArray5d(:,:,:,:,:) = real(doubleArray5d(:,:,:,:,:),RKIND)
            deallocate(doubleArray5d)
         else
            if (associated(field_cursor % fieldhandle % decomp)) then
               call PIO_read_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                   realArray5d, pio_ierr)
            else
               if (field_cursor % fieldhandle % has_unlimited_dim) then
                  start6(:) = 1
                  start6(6) = handle % frame_number
                  count6(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count6(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count6(3) = field_cursor % fieldhandle % dims(3) % dimsize
                  count6(4) = field_cursor % fieldhandle % dims(4) % dimsize
                  count6(5) = field_cursor % fieldhandle % dims(5) % dimsize
                  count6(6) = 1
                  pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start6, count6, realArray5d)
               else
                  start5(:) = 1
                  count5(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count5(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count5(3) = field_cursor % fieldhandle % dims(3) % dimsize
                  count5(4) = field_cursor % fieldhandle % dims(4) % dimsize
                  count5(5) = field_cursor % fieldhandle % dims(5) % dimsize
                  pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start5, count5, realArray5d)
               end if
            end if
         end if
      else if (present(intArray1d)) then
!         call mpas_log_write('  value is int1')
         if (associated(field_cursor % fieldhandle % decomp)) then
            call PIO_read_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                 intArray1d, pio_ierr)
         else
            if (field_cursor % fieldhandle % has_unlimited_dim) then
               start2(1) = 1
               start2(2) = handle % frame_number
               count2(1) = field_cursor % fieldhandle % dims(1) % dimsize
               count2(2) = 1
               pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start2, count2, intArray1d)
            else
               start1(:) = 1
               count1(1) = field_cursor % fieldhandle % dims(1) % dimsize
               pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start1, count1, intArray1d)
            end if
         end if
      else if (present(intArray2d)) then
!         call mpas_log_write('  value is int2')
         if (associated(field_cursor % fieldhandle % decomp)) then
            call PIO_read_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                 intArray2d, pio_ierr)
         else
            if (field_cursor % fieldhandle % has_unlimited_dim) then
               start3(:) = 1
               start3(3) = handle % frame_number
               count3(1) = field_cursor % fieldhandle % dims(1) % dimsize
               count3(2) = field_cursor % fieldhandle % dims(2) % dimsize
               count3(3) = 1
               pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start3, count3, intArray2d)
            else
               start2(:) = 1
               count2(1) = field_cursor % fieldhandle % dims(1) % dimsize
               count2(2) = field_cursor % fieldhandle % dims(2) % dimsize
               pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start2, count2, intArray2d)
            end if
         end if
      else if (present(intArray3d)) then
!         call mpas_log_write('  value is int3')
         if (associated(field_cursor % fieldhandle % decomp)) then
            call PIO_read_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                 intArray3d, pio_ierr)
         else
            if (field_cursor % fieldhandle % has_unlimited_dim) then
               start4(:) = 1
               start4(4) = handle % frame_number
               count4(1) = field_cursor % fieldhandle % dims(1) % dimsize
               count4(2) = field_cursor % fieldhandle % dims(2) % dimsize
               count4(3) = field_cursor % fieldhandle % dims(3) % dimsize
               count4(4) = 1
               pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start4, count4, intArray3d)
            else
               start3(:) = 1
               count3(1) = field_cursor % fieldhandle % dims(1) % dimsize
               count3(2) = field_cursor % fieldhandle % dims(2) % dimsize
               count3(3) = field_cursor % fieldhandle % dims(3) % dimsize
               pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start3, count3, intArray3d)
            end if
         end if
      else if (present(intArray4d)) then
!         call mpas_log_write('  value is int4')
         if (associated(field_cursor % fieldhandle % decomp)) then
            call PIO_read_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                 intArray4d, pio_ierr)
         else
            if (field_cursor % fieldhandle % has_unlimited_dim) then
               start5(:) = 1
               start5(5) = handle % frame_number
               count5(1) = field_cursor % fieldhandle % dims(1) % dimsize
               count5(2) = field_cursor % fieldhandle % dims(2) % dimsize
               count5(3) = field_cursor % fieldhandle % dims(3) % dimsize
               count5(4) = field_cursor % fieldhandle % dims(4) % dimsize
               count5(5) = 1
               pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start5, count5, intArray4d)
            else
               start4(:) = 1
               count4(1) = field_cursor % fieldhandle % dims(1) % dimsize
               count4(2) = field_cursor % fieldhandle % dims(2) % dimsize
               count4(3) = field_cursor % fieldhandle % dims(3) % dimsize
               count4(4) = field_cursor % fieldhandle % dims(4) % dimsize
               pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start4, count4, intArray4d)
            end if
         end if
      end if

!      call mpas_log_write('Checking for error')
      if (pio_ierr /= PIO_noerr) then
         if (present(ierr)) ierr = MPAS_IO_ERR_PIO
         return
      end if

   end subroutine MPAS_io_get_var_generic


   subroutine MPAS_io_get_var_int0d(handle, fieldname, val, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      integer, intent(out) :: val
      integer, intent(out), optional :: ierr

!      call mpas_log_write('Called MPAS_io_get_var_int0d()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_get_var_generic(handle, fieldname, intVal=val, ierr=ierr)

   end subroutine MPAS_io_get_var_int0d


   subroutine MPAS_io_get_var_int1d(handle, fieldname, array, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      integer, dimension(:), intent(out) :: array
      integer, intent(out), optional :: ierr

!      call mpas_log_write('Called MPAS_io_get_var_int1d()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_get_var_generic(handle, fieldname, intArray1d=array, ierr=ierr)

   end subroutine MPAS_io_get_var_int1d


   subroutine MPAS_io_get_var_int2d(handle, fieldname, array, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      integer, dimension(:,:), intent(out) :: array
      integer, intent(out), optional :: ierr

!      call mpas_log_write('Called MPAS_io_get_var_int2d()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_get_var_generic(handle, fieldname, intArray2d=array, ierr=ierr)

   end subroutine MPAS_io_get_var_int2d


   subroutine MPAS_io_get_var_int3d(handle, fieldname, array, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      integer, dimension(:,:,:), intent(out) :: array
      integer, intent(out), optional :: ierr

!      call mpas_log_write('Called MPAS_io_get_var_int3d()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_get_var_generic(handle, fieldname, intArray3d=array, ierr=ierr)

   end subroutine MPAS_io_get_var_int3d


   subroutine MPAS_io_get_var_int4d(handle, fieldname, array, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      integer, dimension(:,:,:,:), intent(out) :: array
      integer, intent(out), optional :: ierr

!      call mpas_log_write('Called MPAS_io_get_var_int4d()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_get_var_generic(handle, fieldname, intArray4d=array, ierr=ierr)

   end subroutine MPAS_io_get_var_int4d


   subroutine MPAS_io_get_var_real0d(handle, fieldname, val, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      real (kind=RKIND), intent(out) :: val
      integer, intent(out), optional :: ierr

!      call mpas_log_write('Called MPAS_io_get_var_real0d()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_get_var_generic(handle, fieldname, realVal=val, ierr=ierr)

   end subroutine MPAS_io_get_var_real0d


   subroutine MPAS_io_get_var_real1d(handle, fieldname, array, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      real (kind=RKIND), dimension(:), intent(out) :: array
      integer, intent(out), optional :: ierr

!      call mpas_log_write('Called MPAS_io_get_var_real1d()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_get_var_generic(handle, fieldname, realArray1d=array, ierr=ierr)

   end subroutine MPAS_io_get_var_real1d


   subroutine MPAS_io_get_var_real2d(handle, fieldname, array, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      real (kind=RKIND), dimension(:,:), intent(out) :: array
      integer, intent(out), optional :: ierr

!      call mpas_log_write('Called MPAS_io_get_var_real2d()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_get_var_generic(handle, fieldname, realArray2d=array, ierr=ierr)

   end subroutine MPAS_io_get_var_real2d


   subroutine MPAS_io_get_var_real3d(handle, fieldname, array, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      real (kind=RKIND), dimension(:,:,:), intent(out) :: array
      integer, intent(out), optional :: ierr

!      call mpas_log_write('Called MPAS_io_get_var_real3d()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_get_var_generic(handle, fieldname, realArray3d=array, ierr=ierr)

   end subroutine MPAS_io_get_var_real3d


   subroutine MPAS_io_get_var_real4d(handle, fieldname, array, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      real (kind=RKIND), dimension(:,:,:,:), intent(out) :: array
      integer, intent(out), optional :: ierr

!      call mpas_log_write('Called MPAS_io_get_var_real4d()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_get_var_generic(handle, fieldname, realArray4d=array, ierr=ierr)

   end subroutine MPAS_io_get_var_real4d


   subroutine MPAS_io_get_var_real5d(handle, fieldname, array, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      real (kind=RKIND), dimension(:,:,:,:,:), intent(out) :: array
      integer, intent(out), optional :: ierr

!      call mpas_log_write('Called MPAS_io_get_var_real5d()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_get_var_generic(handle, fieldname, realArray5d=array, ierr=ierr)

   end subroutine MPAS_io_get_var_real5d


   subroutine MPAS_io_get_var_char0d(handle, fieldname, val, ierr)

      use mpas_c_interfacing, only : MPAS_sanitize_string

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      character (len=*), intent(out) :: val
      integer, intent(out), optional :: ierr

!      call mpas_log_write('Called MPAS_io_get_var_char0d()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_get_var_generic(handle, fieldname, charVal=val, ierr=ierr)
      call MPAS_sanitize_string(val)

   end subroutine MPAS_io_get_var_char0d


   subroutine MPAS_io_get_var_char1d(handle, fieldname, val, ierr)

      use mpas_c_interfacing, only : MPAS_sanitize_string

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      character (len=*), dimension(:), intent(out) :: val
      integer, intent(out), optional :: ierr

      integer :: i

!     call mpas_log_write('Called MPAS_io_get_var_char1d()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_get_var_generic(handle, fieldname, charArray1d=val, ierr=ierr)
      do i=1,size(val)
         call MPAS_sanitize_string(val(i))
      end do

   end subroutine MPAS_io_get_var_char1d


   logical function MPAS_io_would_clobber_records(handle, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      integer, intent(out), optional :: ierr


      if (present(ierr)) ierr = MPAS_IO_NOERR

      ! Sanity checks
      if (.not. handle % initialized) then
         if (present(ierr)) ierr = MPAS_IO_ERR_UNINIT_HANDLE
         MPAS_io_would_clobber_records = .false.
         return 
      end if

      if (handle % frame_number <= handle % preexisting_records) then
         MPAS_io_would_clobber_records = .true.
      else
         MPAS_io_would_clobber_records = .false.
      end if

   end function MPAS_io_would_clobber_records


   subroutine MPAS_io_put_var_generic(handle, fieldname, intVal, intArray1d, intArray2d, intArray3d, intArray4d, &
                                                        realVal, realArray1d, realArray2d, realArray3d, realArray4d, realArray5d, &
                                                        charVal, charArray1d, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      integer, intent(in), optional :: intVal
      integer, dimension(:), intent(in), optional :: intArray1d
      integer, dimension(:,:), intent(in), optional :: intArray2d
      integer, dimension(:,:,:), intent(in), optional :: intArray3d
      integer, dimension(:,:,:,:), intent(in), optional :: intArray4d
      real (kind=RKIND), intent(in), optional :: realVal
      real (kind=RKIND), dimension(:), intent(in), optional :: realArray1d
      real (kind=RKIND), dimension(:,:), intent(in), optional :: realArray2d
      real (kind=RKIND), dimension(:,:,:), intent(in), optional :: realArray3d
      real (kind=RKIND), dimension(:,:,:,:), intent(in), optional :: realArray4d
      real (kind=RKIND), dimension(:,:,:,:,:), intent(in), optional :: realArray5d
      character (len=*), intent(in), optional :: charVal
      character (len=*), dimension(:), intent(in), optional :: charArray1d
      integer, intent(out), optional :: ierr

      integer :: pio_ierr
      integer, dimension(1) :: start1
      integer, dimension(1) :: count1
      integer, dimension(2) :: start2
      integer, dimension(2) :: count2
      integer, dimension(3) :: start3
      integer, dimension(3) :: count3
      integer, dimension(4) :: start4
      integer, dimension(4) :: count4
      integer, dimension(5) :: start5
      integer, dimension(5) :: count5
      integer, dimension(6) :: start6
      integer, dimension(6) :: count6
      type (fieldlist_type), pointer :: field_cursor

      integer :: i

      real (kind=R4KIND) :: singleVal
      real (kind=R4KIND), dimension(:), allocatable :: singleArray1d
      real (kind=R4KIND), dimension(:,:), allocatable :: singleArray2d
      real (kind=R4KIND), dimension(:,:,:), allocatable :: singleArray3d
      real (kind=R4KIND), dimension(:,:,:,:), allocatable :: singleArray4d
      real (kind=R4KIND), dimension(:,:,:,:,:), allocatable :: singleArray5d

      real (kind=R8KIND) :: doubleVal
      real (kind=R8KIND), dimension(:), allocatable :: doubleArray1d
      real (kind=R8KIND), dimension(:,:), allocatable :: doubleArray2d
      real (kind=R8KIND), dimension(:,:,:), allocatable :: doubleArray3d
      real (kind=R8KIND), dimension(:,:,:,:), allocatable :: doubleArray4d
      real (kind=R8KIND), dimension(:,:,:,:,:), allocatable :: doubleArray5d

      ! Sanity checks
      if (.not. handle % initialized) then
         if (present(ierr)) ierr = MPAS_IO_ERR_UNINIT_HANDLE
         return 
      end if

      if (.not. handle % data_mode) then
         handle % data_mode = .true.

         pio_ierr = PIO_enddef(handle % pio_file)
         ! If we are working with a preexisting file, we likely didn't define
         ! new dimensions or variables, in which case PIO_enddef() will return
         ! an error under harmless circumstances; so, don't return only for
         ! pre-existing files.
         if (pio_ierr /= PIO_noerr .and. (.not. handle % preexisting_file)) then
            if (present(ierr)) ierr = MPAS_IO_ERR_PIO
            return
         end if
      end if

!     call mpas_log_write('Writing '//trim(fieldname))

      !
      ! Check whether the field has been defined
      !
      field_cursor => handle % fieldlist_head
      do while (associated(field_cursor))
         if (trim(fieldname) == trim(field_cursor % fieldhandle % fieldname)) then
            exit
         end if
         field_cursor => field_cursor % next
      end do
      if (.not. associated(field_cursor)) then
         if (present(ierr)) ierr = MPAS_IO_ERR_UNDEFINED_VAR
         return
      end if

      if (field_cursor % fieldhandle % has_unlimited_dim) then

#ifdef USE_PIO2
         call PIO_setframe(handle % pio_file, field_cursor % fieldhandle % field_desc, handle % frame_number)
#else
         call PIO_setframe(field_cursor % fieldhandle % field_desc, handle % frame_number)
#endif
         start1(1) = handle % frame_number
         count1(1) = 1
     
         start2(1) = 1
         start2(2) = handle % frame_number
         count2(2) = 1
      else if (handle % frame_number > 1) then
         if (present(ierr)) ierr = MPAS_IO_NOERR
         return
      end if

      if (present(realVal)) then
         if ((field_cursor % fieldhandle % precision == MPAS_IO_SINGLE_PRECISION) .and. &
             (MPAS_IO_NATIVE_PRECISION /= MPAS_IO_SINGLE_PRECISION)) then
            singleVal = real(realVal,R4KIND)
            if (field_cursor % fieldhandle % has_unlimited_dim) then
               pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start1, singleVal)
            else
               pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, singleVal)
            end if
         else if ((field_cursor % fieldhandle % precision == MPAS_IO_DOUBLE_PRECISION) .and. &
             (MPAS_IO_NATIVE_PRECISION /= MPAS_IO_DOUBLE_PRECISION)) then
            doubleVal = real(realVal,R8KIND)
            if (field_cursor % fieldhandle % has_unlimited_dim) then
               pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start1, doubleVal)
            else
               pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, doubleVal)
            end if
         else
            if (field_cursor % fieldhandle % has_unlimited_dim) then
               pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start1, realVal)
            else
               pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, realVal)
            end if
         end if
      else if (present(intVal)) then
         if (field_cursor % fieldhandle % has_unlimited_dim) then
            pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start1, intVal)
         else
            pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, intVal)
         end if
      else if (present(charVal)) then
         if (field_cursor % fieldhandle % has_unlimited_dim) then
            count2(1) = field_cursor % fieldhandle % dims(1) % dimsize
            pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start2, charVal(1:count2(1)))
         else
            start1(1) = 1
            count1(1) = field_cursor % fieldhandle % dims(1) % dimsize
            pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start1, charVal(1:count1(1)))
         end if
      else if (present(charArray1d)) then
         if (field_cursor % fieldhandle % has_unlimited_dim) then
            ! Write one string at a time because the sizes differ so much (i.e. StrLen != StrKIND)
            do i = 1, field_cursor % fieldhandle % dims(2) % dimsize
               start3(1) = 1
               start3(2) = i
               start3(3) = handle % frame_number
               count3(1) = field_cursor % fieldhandle % dims(1) % dimsize ! Should be StrLen
               count3(2) = 1
               count3(3) = 1
               pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start3, charArray1d(i)(1:count3(1)))
            end do
         else
            ! Write one string at a time because the sizes differ so much (i.e. StrLen != StrKIND)
            do i = 1, field_cursor % fieldhandle % dims(2) % dimsize
               start2(1) = 1
               start2(2) = i
               count2(1) = field_cursor % fieldhandle % dims(1) % dimsize ! Should be StrLen
               count2(2) = 1
               pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start2, charArray1d(i)(1:count2(1)))
            end do
         end if
      else if (present(realArray1d)) then
         if ((field_cursor % fieldhandle % precision == MPAS_IO_SINGLE_PRECISION) .and. &
             (MPAS_IO_NATIVE_PRECISION /= MPAS_IO_SINGLE_PRECISION)) then
            allocate(singleArray1d(size(realArray1d)))
            singleArray1d(:) = real(realArray1d(:),R4KIND)
            if (associated(field_cursor % fieldhandle % decomp)) then
               call PIO_write_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                     singleArray1d, pio_ierr)
            else
               if (field_cursor % fieldhandle % has_unlimited_dim) then
                  start2(1) = 1
                  start2(2) = handle % frame_number
                  count2(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count2(2) = 1
                  pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start2, count2, singleArray1d)
               else
                  start1(1) = 1
                  count1(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start1, count1, singleArray1d)
               end if
            end if
            deallocate(singleArray1d)
         else if ((field_cursor % fieldhandle % precision == MPAS_IO_DOUBLE_PRECISION) .and. &
             (MPAS_IO_NATIVE_PRECISION /= MPAS_IO_DOUBLE_PRECISION)) then
            allocate(doubleArray1d(size(realArray1d)))
            doubleArray1d(:) = real(realArray1d(:),R8KIND)
            if (associated(field_cursor % fieldhandle % decomp)) then
               call PIO_write_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                     doubleArray1d, pio_ierr)
            else
               if (field_cursor % fieldhandle % has_unlimited_dim) then
                  start2(1) = 1
                  start2(2) = handle % frame_number
                  count2(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count2(2) = 1
                  pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start2, count2, doubleArray1d)
               else
                  start1(1) = 1
                  count1(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start1, count1, doubleArray1d)
               end if
            end if
            deallocate(doubleArray1d)
         else
            if (associated(field_cursor % fieldhandle % decomp)) then
               call PIO_write_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                     realArray1d, pio_ierr)
            else
               if (field_cursor % fieldhandle % has_unlimited_dim) then
                  start2(1) = 1
                  start2(2) = handle % frame_number
                  count2(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count2(2) = 1
                  pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start2, count2, realArray1d)
               else
                  start1(1) = 1
                  count1(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start1, count1, realArray1d)
               end if
            end if
         end if
      else if (present(realArray2d)) then
         if ((field_cursor % fieldhandle % precision == MPAS_IO_SINGLE_PRECISION) .and. &
             (MPAS_IO_NATIVE_PRECISION /= MPAS_IO_SINGLE_PRECISION)) then
            allocate(singleArray2d(size(realArray2d,1), size(realArray2d,2)))
            singleArray2d(:,:) = real(realArray2d(:,:),R4KIND)
            if (associated(field_cursor % fieldhandle % decomp)) then
               call PIO_write_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                     singleArray2d, pio_ierr)
            else
               if (field_cursor % fieldhandle % has_unlimited_dim) then
                  start3(1) = 1
                  start3(2) = 1
                  start3(3) = handle % frame_number
                  count3(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count3(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count3(3) = 1
                  pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start3, count3, singleArray2d)
               else
                  start2(:) = 1
                  count2(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count2(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start2, count2, singleArray2d)
               end if
            end if
            deallocate(singleArray2d)
         else if ((field_cursor % fieldhandle % precision == MPAS_IO_DOUBLE_PRECISION) .and. &
             (MPAS_IO_NATIVE_PRECISION /= MPAS_IO_DOUBLE_PRECISION)) then
            allocate(doubleArray2d(size(realArray2d,1), size(realArray2d,2)))
            doubleArray2d(:,:) = real(realArray2d(:,:),R8KIND)
            if (associated(field_cursor % fieldhandle % decomp)) then
               call PIO_write_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                     doubleArray2d, pio_ierr)
            else
               if (field_cursor % fieldhandle % has_unlimited_dim) then
                  start3(1) = 1
                  start3(2) = 1
                  start3(3) = handle % frame_number
                  count3(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count3(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count3(3) = 1
                  pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start3, count3, doubleArray2d)
               else
                  start2(:) = 1
                  count2(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count2(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start2, count2, doubleArray2d)
               end if
            end if
            deallocate(doubleArray2d)
         else
            if (associated(field_cursor % fieldhandle % decomp)) then
               call PIO_write_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                     realArray2d, pio_ierr)
            else
               if (field_cursor % fieldhandle % has_unlimited_dim) then
                  start3(1) = 1
                  start3(2) = 1
                  start3(3) = handle % frame_number
                  count3(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count3(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count3(3) = 1
                  pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start3, count3, realArray2d)
               else
                  start2(:) = 1
                  count2(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count2(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start2, count2, realArray2d)
               end if
            end if
         end if
      else if (present(realArray3d)) then
         if ((field_cursor % fieldhandle % precision == MPAS_IO_SINGLE_PRECISION) .and. &
             (MPAS_IO_NATIVE_PRECISION /= MPAS_IO_SINGLE_PRECISION)) then
            allocate(singleArray3d(size(realArray3d,1), size(realArray3d,2), size(realArray3d,3)))
            singleArray3d(:,:,:) = real(realArray3d(:,:,:),R4KIND)
            if (associated(field_cursor % fieldhandle % decomp)) then
               call PIO_write_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                     singleArray3d, pio_ierr)
            else
               if (field_cursor % fieldhandle % has_unlimited_dim) then
                  start4(1) = 1
                  start4(2) = 1
                  start4(3) = 1
                  start4(4) = handle % frame_number
                  count4(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count4(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count4(3) = field_cursor % fieldhandle % dims(3) % dimsize
                  count4(4) = 1
                  pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start4, count4, singleArray3d)
               else
                  start3(:) = 1
                  count3(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count3(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count3(3) = field_cursor % fieldhandle % dims(3) % dimsize
                  pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start3, count3, singleArray3d)
               end if
            end if
            deallocate(singleArray3d)
         else if ((field_cursor % fieldhandle % precision == MPAS_IO_DOUBLE_PRECISION) .and. &
             (MPAS_IO_NATIVE_PRECISION /= MPAS_IO_DOUBLE_PRECISION)) then
            allocate(doubleArray3d(size(realArray3d,1), size(realArray3d,2), size(realArray3d,3)))
            doubleArray3d(:,:,:) = real(realArray3d(:,:,:),R8KIND)
            if (associated(field_cursor % fieldhandle % decomp)) then
               call PIO_write_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                     doubleArray3d, pio_ierr)
            else
               if (field_cursor % fieldhandle % has_unlimited_dim) then
                  start4(1) = 1
                  start4(2) = 1
                  start4(3) = 1
                  start4(4) = handle % frame_number
                  count4(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count4(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count4(3) = field_cursor % fieldhandle % dims(3) % dimsize
                  count4(4) = 1
                  pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start4, count4, doubleArray3d)
               else
                  start3(:) = 1
                  count3(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count3(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count3(3) = field_cursor % fieldhandle % dims(3) % dimsize
                  pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start3, count3, doubleArray3d)
               end if
            end if
            deallocate(doubleArray3d)
         else
            if (associated(field_cursor % fieldhandle % decomp)) then
               call PIO_write_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                     realArray3d, pio_ierr)
            else
               if (field_cursor % fieldhandle % has_unlimited_dim) then
                  start4(1) = 1
                  start4(2) = 1
                  start4(3) = 1
                  start4(4) = handle % frame_number
                  count4(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count4(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count4(3) = field_cursor % fieldhandle % dims(3) % dimsize
                  count4(4) = 1
                  pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start4, count4, realArray3d)
               else
                  start3(:) = 1
                  count3(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count3(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count3(3) = field_cursor % fieldhandle % dims(3) % dimsize
                  pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start3, count3, realArray3d)
               end if
            end if
         end if
      else if (present(realArray4d)) then
         if ((field_cursor % fieldhandle % precision == MPAS_IO_SINGLE_PRECISION) .and. &
             (MPAS_IO_NATIVE_PRECISION /= MPAS_IO_SINGLE_PRECISION)) then
            allocate(singleArray4d(size(realArray4d,1), size(realArray4d,2), size(realArray4d,3), size(realArray4d,4)))
            singleArray4d(:,:,:,:) = real(realArray4d(:,:,:,:),R4KIND)
            if (associated(field_cursor % fieldhandle % decomp)) then
               call PIO_write_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                     singleArray4d, pio_ierr)
            else
               if (field_cursor % fieldhandle % has_unlimited_dim) then
                  start5(1) = 1
                  start5(2) = 1
                  start5(3) = 1
                  start5(4) = 1
                  start5(5) = handle % frame_number
                  count5(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count5(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count5(3) = field_cursor % fieldhandle % dims(3) % dimsize
                  count5(4) = field_cursor % fieldhandle % dims(4) % dimsize
                  count5(5) = 1
                  pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start5, count5, singleArray4d)
               else
                  start4(:) = 1
                  count4(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count4(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count4(3) = field_cursor % fieldhandle % dims(3) % dimsize
                  count4(4) = field_cursor % fieldhandle % dims(4) % dimsize
                  pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start4, count4, singleArray4d)
               end if
            end if
            deallocate(singleArray4d)
         else if ((field_cursor % fieldhandle % precision == MPAS_IO_DOUBLE_PRECISION) .and. &
             (MPAS_IO_NATIVE_PRECISION /= MPAS_IO_DOUBLE_PRECISION)) then
            allocate(doubleArray4d(size(realArray4d,1), size(realArray4d,2), size(realArray4d,3), size(realArray4d,4)))
            doubleArray4d(:,:,:,:) = real(realArray4d(:,:,:,:),R8KIND)
            if (associated(field_cursor % fieldhandle % decomp)) then
               call PIO_write_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                     doubleArray4d, pio_ierr)
            else
               if (field_cursor % fieldhandle % has_unlimited_dim) then
                  start5(1) = 1
                  start5(2) = 1
                  start5(3) = 1
                  start5(4) = 1
                  start5(5) = handle % frame_number
                  count5(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count5(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count5(3) = field_cursor % fieldhandle % dims(3) % dimsize
                  count5(4) = field_cursor % fieldhandle % dims(4) % dimsize
                  count5(5) = 1
                  pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start5, count5, doubleArray4d)
               else
                  start4(:) = 1
                  count4(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count4(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count4(3) = field_cursor % fieldhandle % dims(3) % dimsize
                  count4(4) = field_cursor % fieldhandle % dims(4) % dimsize
                  pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start4, count4, doubleArray4d)
               end if
            end if
            deallocate(doubleArray4d)
         else
            if (associated(field_cursor % fieldhandle % decomp)) then
               call PIO_write_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                     realArray4d, pio_ierr)
            else
               if (field_cursor % fieldhandle % has_unlimited_dim) then
                  start5(1) = 1
                  start5(2) = 1
                  start5(3) = 1
                  start5(4) = 1
                  start5(5) = handle % frame_number
                  count5(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count5(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count5(3) = field_cursor % fieldhandle % dims(3) % dimsize
                  count5(4) = field_cursor % fieldhandle % dims(4) % dimsize
                  count5(5) = 1
                  pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start5, count5, realArray4d)
               else
                  start4(:) = 1
                  count4(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count4(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count4(3) = field_cursor % fieldhandle % dims(3) % dimsize
                  count4(4) = field_cursor % fieldhandle % dims(4) % dimsize
                  pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start4, count4, realArray4d)
               end if
            end if
         end if
      else if (present(realArray5d)) then
         if ((field_cursor % fieldhandle % precision == MPAS_IO_SINGLE_PRECISION) .and. &
             (MPAS_IO_NATIVE_PRECISION /= MPAS_IO_SINGLE_PRECISION)) then
            allocate(singleArray5d(size(realArray5d,1), size(realArray5d,2), size(realArray5d,3), size(realArray5d,4), size(realArray5d,5)))
            singleArray5d(:,:,:,:,:) = real(realArray5d(:,:,:,:,:),R4KIND)
            if (associated(field_cursor % fieldhandle % decomp)) then
               call PIO_write_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                     singleArray5d, pio_ierr)
            else
               if (field_cursor % fieldhandle % has_unlimited_dim) then
                  start6(1) = 1
                  start6(2) = 1
                  start6(3) = 1
                  start6(4) = 1
                  start6(5) = 1
                  start6(6) = handle % frame_number
                  count6(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count6(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count6(3) = field_cursor % fieldhandle % dims(3) % dimsize
                  count6(4) = field_cursor % fieldhandle % dims(4) % dimsize
                  count6(5) = field_cursor % fieldhandle % dims(5) % dimsize
                  count6(6) = 1
                  pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start6, count6, singleArray5d)
               else
                  start5(:) = 1
                  count5(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count5(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count5(3) = field_cursor % fieldhandle % dims(3) % dimsize
                  count5(4) = field_cursor % fieldhandle % dims(4) % dimsize
                  count5(5) = field_cursor % fieldhandle % dims(5) % dimsize
                  pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start5, count5, singleArray5d)
               end if
            end if
            deallocate(singleArray5d)
         else if ((field_cursor % fieldhandle % precision == MPAS_IO_DOUBLE_PRECISION) .and. &
             (MPAS_IO_NATIVE_PRECISION /= MPAS_IO_DOUBLE_PRECISION)) then
            allocate(doubleArray5d(size(realArray5d,1), size(realArray5d,2), size(realArray5d,3), size(realArray5d,4), size(realArray5d,5)))
            doubleArray5d(:,:,:,:,:) = real(realArray5d(:,:,:,:,:),R8KIND)
            if (associated(field_cursor % fieldhandle % decomp)) then
               call PIO_write_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                     doubleArray5d, pio_ierr)
            else
               if (field_cursor % fieldhandle % has_unlimited_dim) then
                  start6(1) = 1
                  start6(2) = 1
                  start6(3) = 1
                  start6(4) = 1
                  start6(5) = 1
                  start6(6) = handle % frame_number
                  count6(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count6(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count6(3) = field_cursor % fieldhandle % dims(3) % dimsize
                  count6(4) = field_cursor % fieldhandle % dims(4) % dimsize
                  count6(5) = field_cursor % fieldhandle % dims(5) % dimsize
                  count6(6) = 1
                  pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start6, count6, doubleArray5d)
               else
                  start5(:) = 1
                  count5(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count5(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count5(3) = field_cursor % fieldhandle % dims(3) % dimsize
                  count5(4) = field_cursor % fieldhandle % dims(4) % dimsize
                  count5(5) = field_cursor % fieldhandle % dims(5) % dimsize
                  pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start5, count5, doubleArray5d)
               end if
            end if
            deallocate(doubleArray5d)
         else
            if (associated(field_cursor % fieldhandle % decomp)) then
               call PIO_write_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                     realArray5d, pio_ierr)
            else
               if (field_cursor % fieldhandle % has_unlimited_dim) then
                  start6(1) = 1
                  start6(2) = 1
                  start6(3) = 1
                  start6(4) = 1
                  start6(5) = 1
                  start6(6) = handle % frame_number
                  count6(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count6(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count6(3) = field_cursor % fieldhandle % dims(3) % dimsize
                  count6(4) = field_cursor % fieldhandle % dims(4) % dimsize
                  count6(5) = field_cursor % fieldhandle % dims(5) % dimsize
                  count6(6) = 1
                  pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start6, count6, realArray5d)
               else
                  start5(:) = 1
                  count5(1) = field_cursor % fieldhandle % dims(1) % dimsize
                  count5(2) = field_cursor % fieldhandle % dims(2) % dimsize
                  count5(3) = field_cursor % fieldhandle % dims(3) % dimsize
                  count5(4) = field_cursor % fieldhandle % dims(4) % dimsize
                  count5(5) = field_cursor % fieldhandle % dims(5) % dimsize
                  pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start5, count5, realArray5d)
               end if
            end if
         end if
      else if (present(intArray1d)) then
         if (associated(field_cursor % fieldhandle % decomp)) then
            call PIO_write_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                  intArray1d, pio_ierr)
         else
            if (field_cursor % fieldhandle % has_unlimited_dim) then
               start2(1) = 1
               start2(2) = handle % frame_number
               count2(1) = field_cursor % fieldhandle % dims(1) % dimsize
               count2(2) = 1
               pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start2, count2, intArray1d)
            else
               start1(1) = 1
               count1(1) = field_cursor % fieldhandle % dims(1) % dimsize
               pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start1, count1, intArray1d)
            end if
         end if
      else if (present(intArray2d)) then
         if (associated(field_cursor % fieldhandle % decomp)) then
            call PIO_write_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                  intArray2d, pio_ierr)
         else
            if (field_cursor % fieldhandle % has_unlimited_dim) then
               start3(1) = 1
               start3(2) = 1
               start3(3) = handle % frame_number
               count3(1) = field_cursor % fieldhandle % dims(1) % dimsize
               count3(2) = field_cursor % fieldhandle % dims(2) % dimsize
               count3(3) = 1
               pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start3, count3, intArray2d)
            else
               start2(:) = 1
               count2(1) = field_cursor % fieldhandle % dims(1) % dimsize
               count2(2) = field_cursor % fieldhandle % dims(2) % dimsize
               pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start2, count2, intArray2d)
            end if
         end if
      else if (present(intArray3d)) then
         if (associated(field_cursor % fieldhandle % decomp)) then
            call PIO_write_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                  intArray3d, pio_ierr)
         else
            if (field_cursor % fieldhandle % has_unlimited_dim) then
               start4(1) = 1
               start4(2) = 1
               start4(3) = 1
               start4(4) = handle % frame_number
               count4(1) = field_cursor % fieldhandle % dims(1) % dimsize
               count4(2) = field_cursor % fieldhandle % dims(2) % dimsize
               count4(3) = field_cursor % fieldhandle % dims(3) % dimsize
               count4(4) = 1
               pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start4, count4, intArray3d)
            else
               start3(:) = 1
               count3(1) = field_cursor % fieldhandle % dims(1) % dimsize
               count3(2) = field_cursor % fieldhandle % dims(2) % dimsize
               count3(3) = field_cursor % fieldhandle % dims(3) % dimsize
               pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start3, count3, intArray3d)
            end if
         end if
      else if (present(intArray4d)) then
         if (associated(field_cursor % fieldhandle % decomp)) then
            call PIO_write_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                  intArray4d, pio_ierr)
         else
            if (field_cursor % fieldhandle % has_unlimited_dim) then
               start5(1) = 1
               start5(2) = 1
               start5(3) = 1
               start5(4) = 1
               start5(5) = handle % frame_number
               count5(1) = field_cursor % fieldhandle % dims(1) % dimsize
               count5(2) = field_cursor % fieldhandle % dims(2) % dimsize
               count5(3) = field_cursor % fieldhandle % dims(3) % dimsize
               count5(4) = field_cursor % fieldhandle % dims(4) % dimsize
               count5(5) = 1
               pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start5, count5, intArray4d)
            else
               start4(:) = 1
               count4(1) = field_cursor % fieldhandle % dims(1) % dimsize
               count4(2) = field_cursor % fieldhandle % dims(2) % dimsize
               count4(3) = field_cursor % fieldhandle % dims(3) % dimsize
               count4(4) = field_cursor % fieldhandle % dims(4) % dimsize
               pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start4, count4, intArray4d)
            end if
         end if
      end if
      if (pio_ierr /= PIO_noerr) then
         if (present(ierr)) ierr = MPAS_IO_ERR_PIO
         return
      end if

   end subroutine MPAS_io_put_var_generic


   subroutine MPAS_io_put_var_int0d(handle, fieldname, val, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      integer, intent(in) :: val
      integer, intent(out), optional :: ierr

!      call mpas_log_write('Called MPAS_io_put_var_int0d()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_put_var_generic(handle, fieldname, intVal=val, ierr=ierr)

   end subroutine MPAS_io_put_var_int0d


   subroutine MPAS_io_put_var_int1d(handle, fieldname, array, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      integer, dimension(:), intent(in) :: array
      integer, intent(out), optional :: ierr

!      call mpas_log_write('Called MPAS_io_put_var_int1d()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_put_var_generic(handle, fieldname, intArray1d=array, ierr=ierr)

   end subroutine MPAS_io_put_var_int1d


   subroutine MPAS_io_put_var_int2d(handle, fieldname, array, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      integer, dimension(:,:), intent(in) :: array
      integer, intent(out), optional :: ierr

!      call mpas_log_write('Called MPAS_io_put_var_int2d()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_put_var_generic(handle, fieldname, intArray2d=array, ierr=ierr)

   end subroutine MPAS_io_put_var_int2d


   subroutine MPAS_io_put_var_int3d(handle, fieldname, array, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      integer, dimension(:,:,:), intent(in) :: array
      integer, intent(out), optional :: ierr

!      call mpas_log_write('Called MPAS_io_put_var_int3d()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_put_var_generic(handle, fieldname, intArray3d=array, ierr=ierr)

   end subroutine MPAS_io_put_var_int3d


   subroutine MPAS_io_put_var_int4d(handle, fieldname, array, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      integer, dimension(:,:,:,:), intent(in) :: array
      integer, intent(out), optional :: ierr

!      call mpas_log_write('Called MPAS_io_put_var_int4d()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_put_var_generic(handle, fieldname, intArray4d=array, ierr=ierr)

   end subroutine MPAS_io_put_var_int4d


   subroutine MPAS_io_put_var_real0d(handle, fieldname, val, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      real (kind=RKIND), intent(in) :: val
      integer, intent(out), optional :: ierr

!      call mpas_log_write('Called MPAS_io_put_var_real0d()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_put_var_generic(handle, fieldname, realVal=val, ierr=ierr)

   end subroutine MPAS_io_put_var_real0d


   subroutine MPAS_io_put_var_real1d(handle, fieldname, array, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      real (kind=RKIND), dimension(:), intent(in) :: array
      integer, intent(out), optional :: ierr

!      call mpas_log_write('Called MPAS_io_put_var_real1d()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_put_var_generic(handle, fieldname, realArray1d=array, ierr=ierr)

   end subroutine MPAS_io_put_var_real1d


   subroutine MPAS_io_put_var_real2d(handle, fieldname, array, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      real (kind=RKIND), dimension(:,:), intent(in) :: array
      integer, intent(out), optional :: ierr

!      call mpas_log_write('Called MPAS_io_put_var_real2d()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_put_var_generic(handle, fieldname, realArray2d=array, ierr=ierr)

   end subroutine MPAS_io_put_var_real2d


   subroutine MPAS_io_put_var_real3d(handle, fieldname, array, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      real (kind=RKIND), dimension(:,:,:), intent(in) :: array
      integer, intent(out), optional :: ierr

!      call mpas_log_write('Called MPAS_io_put_var_real3d()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_put_var_generic(handle, fieldname, realArray3d=array, ierr=ierr)

   end subroutine MPAS_io_put_var_real3d


   subroutine MPAS_io_put_var_real4d(handle, fieldname, array, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      real (kind=RKIND), dimension(:,:,:,:), intent(in) :: array
      integer, intent(out), optional :: ierr

!      call mpas_log_write('Called MPAS_io_put_var_real4d()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_put_var_generic(handle, fieldname, realArray4d=array, ierr=ierr)

   end subroutine MPAS_io_put_var_real4d


   subroutine MPAS_io_put_var_real5d(handle, fieldname, array, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      real (kind=RKIND), dimension(:,:,:,:,:), intent(in) :: array
      integer, intent(out), optional :: ierr

!      call mpas_log_write('Called MPAS_io_put_var_real5d()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_put_var_generic(handle, fieldname, realArray5d=array, ierr=ierr)

   end subroutine MPAS_io_put_var_real5d


   subroutine MPAS_io_put_var_char0d(handle, fieldname, val, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      character (len=*), intent(in) :: val
      integer, intent(out), optional :: ierr

!      call mpas_log_write('Called MPAS_io_put_var_char0d()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_put_var_generic(handle, fieldname, charVal=val, ierr=ierr)

   end subroutine MPAS_io_put_var_char0d


   subroutine MPAS_io_put_var_char1d(handle, fieldname, val, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      character (len=*), dimension(:), intent(in) :: val
      integer, intent(out), optional :: ierr

!      call mpas_log_write('Called MPAS_io_put_var_char1d()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_put_var_generic(handle, fieldname, charArray1d=val, ierr=ierr)

   end subroutine MPAS_io_put_var_char1d


   subroutine MPAS_io_get_att_int0d(handle, attName, attValue, fieldname, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: attName
      integer, intent(out) :: attValue
      character (len=*), intent(in), optional :: fieldname
      integer, intent(out), optional :: ierr

      integer :: pio_ierr
      integer :: varid
      integer :: xtype
#ifdef USE_PIO2
      integer (kind=MPAS_IO_OFFSET_KIND) :: len
#else
      integer :: len
#endif
      type (fieldlist_type), pointer :: field_cursor
      type (attlist_type), pointer :: att_cursor, new_att_node

!      call mpas_log_write('Called MPAS_io_get_att_int0d()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      ! Sanity checks
      if (.not. handle % initialized) then
         if (present(ierr)) ierr = MPAS_IO_ERR_UNINIT_HANDLE
         return 
      end if


      !
      ! For variable attributes, find the structure for fieldname
      !
      if (present(fieldname)) then
         field_cursor => handle % fieldlist_head
         do while (associated(field_cursor))
            if (trim(fieldname) == trim(field_cursor % fieldhandle % fieldname)) then
               varid = field_cursor % fieldhandle % fieldid
               exit
            end if
            field_cursor => field_cursor % next
         end do
         if (.not. associated(field_cursor)) then
            if (present(ierr)) ierr = MPAS_IO_ERR_UNDEFINED_VAR
            return
         end if

         ! Check whether we have this attribute cached
         att_cursor => field_cursor % fieldhandle % attlist_head
         do while (associated(att_cursor))
            if (trim(att_cursor % atthandle % attName) == trim(attName)) then
               if (att_cursor % atthandle % attType == MPAS_ATT_INT) then
!call mpas_log_write('Using cached attribute')
                  attValue = att_cursor % atthandle % attValueInt
               else
                  if (present(ierr)) ierr=MPAS_IO_ERR_WRONG_ATT_TYPE
               end if
               return
            end if
            att_cursor => att_cursor % next
         end do

      else

         ! Check whether we have this attribute cached
         att_cursor => handle % attlist_head
         do while (associated(att_cursor))
            if (trim(att_cursor % atthandle % attName) == trim(attName)) then
               if (att_cursor % atthandle % attType == MPAS_ATT_INT) then
!call mpas_log_write('Using cached attribute')
                  attValue = att_cursor % atthandle % attValueInt
               else
                  if (present(ierr)) ierr=MPAS_IO_ERR_WRONG_ATT_TYPE
               end if
               return
            end if
            att_cursor => att_cursor % next
         end do

         varid = PIO_global
      end if

      ! Query attribute value
      pio_ierr = PIO_inq_att(handle % pio_file, varid, attName, xtype, len)
      if (pio_ierr /= PIO_noerr) then
         if (present(ierr)) ierr = MPAS_IO_ERR_PIO
         return
      end if
      if (xtype /= PIO_int) then
         if (present(ierr)) ierr=MPAS_IO_ERR_WRONG_ATT_TYPE
         return
      end if

      pio_ierr = PIO_get_att(handle % pio_file, varid, attName, attValue)
      if (pio_ierr /= PIO_noerr) then
         if (present(ierr)) ierr = MPAS_IO_ERR_PIO
         return
      end if

      ! Keep attribute for future reference
      allocate(new_att_node)
      nullify(new_att_node % next)
      allocate(new_att_node % atthandle)
      new_att_node % atthandle % attName = attName
      new_att_node % atthandle % attType = MPAS_ATT_INT
      new_att_node % atthandle % attValueInt = attValue

      if (present(fieldname)) then
         if (.not. associated(field_cursor % fieldhandle % attlist_head)) then
            field_cursor % fieldhandle % attlist_head => new_att_node
!call mpas_log_write('Assigning att head for '//trim(attName))
         end if
         if (.not. associated(field_cursor % fieldhandle % attlist_tail)) then
            field_cursor % fieldhandle % attlist_tail => new_att_node
!call mpas_log_write('Assigning att tail for '//trim(attName))
         else
            field_cursor % fieldhandle % attlist_tail % next => new_att_node
            field_cursor % fieldhandle % attlist_tail => field_cursor % fieldhandle % attlist_tail % next
!call mpas_log_write('Extending att tail for '//trim(attName))
         end if
      else
         if (.not. associated(handle % attlist_head)) then
            handle % attlist_head => new_att_node
!call mpas_log_write('Assigning att head for '//trim(attName))
         end if
         if (.not. associated(handle % attlist_tail)) then
            handle % attlist_tail => new_att_node
!call mpas_log_write('Assigning att tail for '//trim(attName))
         else
            handle % attlist_tail % next => new_att_node
            handle % attlist_tail => handle % attlist_tail % next
!call mpas_log_write('Extending att tail for '//trim(attName))
         end if
      end if

   end subroutine MPAS_io_get_att_int0d


   subroutine MPAS_io_get_att_int1d(handle, attName, attValue, fieldname, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: attName
      integer, dimension(:), pointer :: attValue
      character (len=*), intent(in), optional :: fieldname
      integer, intent(out), optional :: ierr

      integer :: pio_ierr
      integer :: varid
      integer :: xtype
#ifdef USE_PIO2
      integer (kind=MPAS_IO_OFFSET_KIND) :: len, attlen
#else
      integer :: len, attlen
#endif
      type (fieldlist_type), pointer :: field_cursor
      type (attlist_type), pointer :: att_cursor, new_att_node

!      call mpas_log_write('Called MPAS_io_get_att_int1d()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      ! Sanity checks
      if (.not. handle % initialized) then
         if (present(ierr)) ierr = MPAS_IO_ERR_UNINIT_HANDLE
         return 
      end if


      !
      ! For variable attributes, find the structure for fieldname
      !
      if (present(fieldname)) then
         field_cursor => handle % fieldlist_head
         do while (associated(field_cursor))
            if (trim(fieldname) == trim(field_cursor % fieldhandle % fieldname)) then
               varid = field_cursor % fieldhandle % fieldid
               exit
            end if
            field_cursor => field_cursor % next
         end do
         if (.not. associated(field_cursor)) then
            if (present(ierr)) ierr = MPAS_IO_ERR_UNDEFINED_VAR
            return
         end if

         ! Check whether we have this attribute cached
         att_cursor => field_cursor % fieldhandle % attlist_head
         do while (associated(att_cursor))
            if (trim(att_cursor % atthandle % attName) == trim(attName)) then
               if (att_cursor % atthandle % attType == MPAS_ATT_INTA) then
!call mpas_log_write('Using cached attribute')
                  allocate(attValue(size(att_cursor % atthandle % attValueIntA)))
                  attValue = att_cursor % atthandle % attValueIntA
               else
                  if (present(ierr)) ierr=MPAS_IO_ERR_WRONG_ATT_TYPE
               end if
               return
            end if
            att_cursor => att_cursor % next
         end do

      else

         ! Check whether we have this attribute cached
         att_cursor => handle % attlist_head
         do while (associated(att_cursor))
            if (trim(att_cursor % atthandle % attName) == trim(attName)) then
               if (att_cursor % atthandle % attType == MPAS_ATT_INTA) then
!call mpas_log_write('Using cached attribute')
                  allocate(attValue(size(att_cursor % atthandle % attValueIntA)))
                  attValue = att_cursor % atthandle % attValueIntA
               else
                  if (present(ierr)) ierr=MPAS_IO_ERR_WRONG_ATT_TYPE
               end if
               return
            end if
            att_cursor => att_cursor % next
         end do

         varid = PIO_global
      end if

      ! Query attribute value
      pio_ierr = PIO_inq_att(handle % pio_file, varid, attName, xtype, len)
      if (pio_ierr /= PIO_noerr) then
         if (present(ierr)) ierr = MPAS_IO_ERR_PIO
         return
      end if

      if (xtype /= PIO_int) then
         if (present(ierr)) ierr=MPAS_IO_ERR_WRONG_ATT_TYPE
         return
      end if

      pio_ierr = PIO_inq_attlen(handle % pio_file, varid, attName, attlen)
      if (pio_ierr /= PIO_noerr) then
         if (present(ierr)) ierr = MPAS_IO_ERR_PIO
         return
      end if

      allocate(attValue(attlen))
      pio_ierr = PIO_get_att(handle % pio_file, varid, attName, attValue)
      if (pio_ierr /= PIO_noerr) then
         if (present(ierr)) ierr = MPAS_IO_ERR_PIO
         return
      end if

      ! Keep attribute for future reference
      allocate(new_att_node)
      nullify(new_att_node % next)
      allocate(new_att_node % atthandle)
      new_att_node % atthandle % attName = attName
      new_att_node % atthandle % attType = MPAS_ATT_INTA
      allocate(new_att_node % atthandle % attValueIntA(attlen))
      new_att_node % atthandle % attValueIntA = attValue

      if (present(fieldname)) then
         if (.not. associated(field_cursor % fieldhandle % attlist_head)) then
            field_cursor % fieldhandle % attlist_head => new_att_node
!call mpas_log_write('Assigning att head for '//trim(attName))
         end if
         if (.not. associated(field_cursor % fieldhandle % attlist_tail)) then
            field_cursor % fieldhandle % attlist_tail => new_att_node
!call mpas_log_write('Assigning att tail for '//trim(attName))
         else
            field_cursor % fieldhandle % attlist_tail % next => new_att_node
            field_cursor % fieldhandle % attlist_tail => field_cursor % fieldhandle % attlist_tail % next
!call mpas_log_write('Extending att tail for '//trim(attName))
         end if
      else
         if (.not. associated(handle % attlist_head)) then
            handle % attlist_head => new_att_node
!call mpas_log_write('Assigning att head for '//trim(attName))
         end if
         if (.not. associated(handle % attlist_tail)) then
            handle % attlist_tail => new_att_node
!call mpas_log_write('Assigning att tail for '//trim(attName))
         else
            handle % attlist_tail % next => new_att_node
            handle % attlist_tail => handle % attlist_tail % next
!call mpas_log_write('Extending att tail for '//trim(attName))
         end if
      end if

   end subroutine MPAS_io_get_att_int1d


   subroutine MPAS_io_get_att_real0d(handle, attName, attValue, fieldname, precision, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: attName
      real (kind=RKIND), intent(out) :: attValue
      character (len=*), intent(in), optional :: fieldname
      integer, intent(in), optional :: precision
      integer, intent(out), optional :: ierr

      integer :: pio_ierr
      integer :: varid
      integer :: local_precision
      real (kind=R4KIND) :: singleVal
      real (kind=R8KIND) :: doubleVal
      integer :: xtype
#ifdef USE_PIO2
      integer (kind=MPAS_IO_OFFSET_KIND) :: len
#else
      integer :: len
#endif
      type (fieldlist_type), pointer :: field_cursor
      type (attlist_type), pointer :: att_cursor, new_att_node

!      call mpas_log_write('Called MPAS_io_get_att_real0d()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      ! Sanity checks
      if (.not. handle % initialized) then
         if (present(ierr)) ierr = MPAS_IO_ERR_UNINIT_HANDLE
         return 
      end if


      !
      ! For variable attributes, find the structure for fieldname
      !
      if (present(fieldname)) then
         field_cursor => handle % fieldlist_head
         do while (associated(field_cursor))
            if (trim(fieldname) == trim(field_cursor % fieldhandle % fieldname)) then
               varid = field_cursor % fieldhandle % fieldid
               exit
            end if
            field_cursor => field_cursor % next
         end do
         if (.not. associated(field_cursor)) then
            if (present(ierr)) ierr = MPAS_IO_ERR_UNDEFINED_VAR
            return
         end if

         ! Check whether we have this attribute cached
         att_cursor => field_cursor % fieldhandle % attlist_head
         do while (associated(att_cursor))
            if (trim(att_cursor % atthandle % attName) == trim(attName)) then
               if (att_cursor % atthandle % attType == MPAS_ATT_REAL) then
!call mpas_log_write('Using cached attribute')
                  attValue = att_cursor % atthandle % attValueReal
               else
                  if (present(ierr)) ierr=MPAS_IO_ERR_WRONG_ATT_TYPE
               end if
               return
            end if
            att_cursor => att_cursor % next
         end do

      else

         ! Check whether we have this attribute cached
         att_cursor => handle % attlist_head
         do while (associated(att_cursor))
            if (trim(att_cursor % atthandle % attName) == trim(attName)) then
               if (att_cursor % atthandle % attType == MPAS_ATT_REAL) then
!call mpas_log_write('Using cached attribute')
                  attValue = att_cursor % atthandle % attValueReal
               else
                  if (present(ierr)) ierr=MPAS_IO_ERR_WRONG_ATT_TYPE
               end if
               return
            end if
            att_cursor => att_cursor % next
         end do

         varid = PIO_global
      end if

      if (present(precision)) then
         local_precision = precision
      else
         local_precision = MPAS_IO_NATIVE_PRECISION
      end if

      ! Query attribute value
      pio_ierr = PIO_inq_att(handle % pio_file, varid, attName, xtype, len)
      if (pio_ierr /= PIO_noerr) then
         if (present(ierr)) ierr = MPAS_IO_ERR_PIO
         return
      end if

      if ((local_precision == MPAS_IO_SINGLE_PRECISION) .and. &
          (MPAS_IO_NATIVE_PRECISION /= MPAS_IO_SINGLE_PRECISION)) then

         if (xtype /= PIO_REAL) then
            if (present(ierr)) ierr=MPAS_IO_ERR_WRONG_ATT_TYPE
            return
         end if
         pio_ierr = PIO_get_att(handle % pio_file, varid, attName, singleVal)
         attValue = real(singleVal,RKIND)

      else if ((local_precision == MPAS_IO_DOUBLE_PRECISION) .and. &
          (MPAS_IO_NATIVE_PRECISION /= MPAS_IO_DOUBLE_PRECISION)) then

         if (xtype /= PIO_DOUBLE) then
            if (present(ierr)) ierr=MPAS_IO_ERR_WRONG_ATT_TYPE
            return
         end if
         pio_ierr = PIO_get_att(handle % pio_file, varid, attName, doubleVal)
         attValue = real(doubleVal,RKIND)

      else

         if (xtype /= PIO_DOUBLE .and. xtype /= PIO_REAL) then
            if (present(ierr)) ierr=MPAS_IO_ERR_WRONG_ATT_TYPE
            return
         end if
         pio_ierr = PIO_get_att(handle % pio_file, varid, attName, attValue)

      end if
      if (pio_ierr /= PIO_noerr) then
         if (present(ierr)) ierr = MPAS_IO_ERR_PIO
         return
      end if

      ! Keep attribute for future reference
      allocate(new_att_node)
      nullify(new_att_node % next)
      allocate(new_att_node % atthandle)
      new_att_node % atthandle % attName = attName
      new_att_node % atthandle % attType = MPAS_ATT_REAL
      new_att_node % atthandle % attValueReal = attValue
      new_att_node % atthandle % precision = local_precision

      if (present(fieldname)) then
         if (.not. associated(field_cursor % fieldhandle % attlist_head)) then
            field_cursor % fieldhandle % attlist_head => new_att_node
!call mpas_log_write('Assigning att head for '//trim(attName))
         end if
         if (.not. associated(field_cursor % fieldhandle % attlist_tail)) then
            field_cursor % fieldhandle % attlist_tail => new_att_node
!call mpas_log_write('Assigning att tail for '//trim(attName))
         else
            field_cursor % fieldhandle % attlist_tail % next => new_att_node
            field_cursor % fieldhandle % attlist_tail => field_cursor % fieldhandle % attlist_tail % next
!call mpas_log_write('Extending att tail for '//trim(attName))
         end if
      else
         if (.not. associated(handle % attlist_head)) then
            handle % attlist_head => new_att_node
!call mpas_log_write('Assigning att head for '//trim(attName))
         end if
         if (.not. associated(handle % attlist_tail)) then
            handle % attlist_tail => new_att_node
!call mpas_log_write('Assigning att tail for '//trim(attName))
         else
            handle % attlist_tail % next => new_att_node
            handle % attlist_tail => handle % attlist_tail % next
!call mpas_log_write('Extending att tail for '//trim(attName))
         end if
      end if

   end subroutine MPAS_io_get_att_real0d


   subroutine MPAS_io_get_att_real1d(handle, attName, attValue, fieldname, precision, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: attName
      real (kind=RKIND), dimension(:), pointer :: attValue
      character (len=*), intent(in), optional :: fieldname
      integer, intent(in), optional :: precision
      integer, intent(out), optional :: ierr

      integer :: pio_ierr
      integer :: varid
      integer :: local_precision
      real (kind=R4KIND), dimension(:), allocatable :: singleVal
      real (kind=R8KIND), dimension(:), allocatable :: doubleVal
      integer :: xtype
#ifdef USE_PIO2
      integer (kind=MPAS_IO_OFFSET_KIND) :: len, attlen
#else
      integer :: len, attlen
#endif
      type (fieldlist_type), pointer :: field_cursor
      type (attlist_type), pointer :: att_cursor, new_att_node

!      call mpas_log_write('Called MPAS_io_get_att_real1d()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      ! Sanity checks
      if (.not. handle % initialized) then
         if (present(ierr)) ierr = MPAS_IO_ERR_UNINIT_HANDLE
         return 
      end if


      !
      ! For variable attributes, find the structure for fieldname
      !
      if (present(fieldname)) then
         field_cursor => handle % fieldlist_head
         do while (associated(field_cursor))
            if (trim(fieldname) == trim(field_cursor % fieldhandle % fieldname)) then
               varid = field_cursor % fieldhandle % fieldid
               exit
            end if
            field_cursor => field_cursor % next
         end do
         if (.not. associated(field_cursor)) then
            if (present(ierr)) ierr = MPAS_IO_ERR_UNDEFINED_VAR
            return
         end if

         ! Check whether we have this attribute cached
         att_cursor => field_cursor % fieldhandle % attlist_head
         do while (associated(att_cursor))
            if (trim(att_cursor % atthandle % attName) == trim(attName)) then
               if (att_cursor % atthandle % attType == MPAS_ATT_REALA) then
!call mpas_log_write('Using cached attribute')
                  allocate(attValue(size(att_cursor % atthandle % attValueRealA)))
                  attValue = att_cursor % atthandle % attValueRealA
               else
                  if (present(ierr)) ierr=MPAS_IO_ERR_WRONG_ATT_TYPE
               end if
               return
            end if
            att_cursor => att_cursor % next
         end do

      else

         ! Check whether we have this attribute cached
         att_cursor => handle % attlist_head
         do while (associated(att_cursor))
            if (trim(att_cursor % atthandle % attName) == trim(attName)) then
               if (att_cursor % atthandle % attType == MPAS_ATT_REALA) then
!call mpas_log_write('Using cached attribute')
                  allocate(attValue(size(att_cursor % atthandle % attValueRealA)))
                  attValue = att_cursor % atthandle % attValueRealA
               else
                  if (present(ierr)) ierr=MPAS_IO_ERR_WRONG_ATT_TYPE
               end if
               return
            end if
            att_cursor => att_cursor % next
         end do

         varid = PIO_global
      end if

      if (present(precision)) then
         local_precision = precision
      else
         local_precision = MPAS_IO_NATIVE_PRECISION
      end if

      ! Query attribute value
      pio_ierr = PIO_inq_att(handle % pio_file, varid, attName, xtype, len)
      if (pio_ierr /= PIO_noerr) then
         if (present(ierr)) ierr = MPAS_IO_ERR_PIO
         return
      end if

      pio_ierr = PIO_inq_attlen(handle % pio_file, varid, attName, attlen)
      if (pio_ierr /= PIO_noerr) then
         if (present(ierr)) ierr = MPAS_IO_ERR_PIO
         return
      end if

      if ((local_precision == MPAS_IO_SINGLE_PRECISION) .and. &
          (MPAS_IO_NATIVE_PRECISION /= MPAS_IO_SINGLE_PRECISION)) then

         if (xtype /= PIO_REAL) then
            if (present(ierr)) ierr=MPAS_IO_ERR_WRONG_ATT_TYPE
            return
         end if
         allocate(attValue(attlen))
         allocate(singleVal(attlen))
         pio_ierr = PIO_get_att(handle % pio_file, varid, attName, singleVal)
         attValue(:) = real(singleVal(:),RKIND)
         deallocate(singleVal)

      else if ((local_precision == MPAS_IO_DOUBLE_PRECISION) .and. &
          (MPAS_IO_NATIVE_PRECISION /= MPAS_IO_DOUBLE_PRECISION)) then

         if (xtype /= PIO_DOUBLE) then
            if (present(ierr)) ierr=MPAS_IO_ERR_WRONG_ATT_TYPE
            return
         end if
         allocate(attValue(attlen))
         allocate(doubleVal(attlen))
         pio_ierr = PIO_get_att(handle % pio_file, varid, attName, doubleVal)
         attValue(:) = real(doubleVal(:),RKIND)
         deallocate(doubleVal)

      else

         if (xtype /= PIO_DOUBLE .and. xtype /= PIO_REAL) then
            if (present(ierr)) ierr=MPAS_IO_ERR_WRONG_ATT_TYPE
            return
         end if
         allocate(attValue(attlen))
         pio_ierr = PIO_get_att(handle % pio_file, varid, attName, attValue)

      end if
      if (pio_ierr /= PIO_noerr) then
         if (present(ierr)) ierr = MPAS_IO_ERR_PIO
         return
      end if

      ! Keep attribute for future reference
      allocate(new_att_node)
      nullify(new_att_node % next)
      allocate(new_att_node % atthandle)
      new_att_node % atthandle % attName = attName
      new_att_node % atthandle % attType = MPAS_ATT_REALA
      allocate(new_att_node % atthandle % attValueRealA(attlen))
      new_att_node % atthandle % attValueRealA = attValue
      new_att_node % atthandle % precision = local_precision

      if (present(fieldname)) then
         if (.not. associated(field_cursor % fieldhandle % attlist_head)) then
            field_cursor % fieldhandle % attlist_head => new_att_node
!call mpas_log_write('Assigning att head for '//trim(attName))
         end if
         if (.not. associated(field_cursor % fieldhandle % attlist_tail)) then
            field_cursor % fieldhandle % attlist_tail => new_att_node
!call mpas_log_write('Assigning att tail for '//trim(attName))
         else
            field_cursor % fieldhandle % attlist_tail % next => new_att_node
            field_cursor % fieldhandle % attlist_tail => field_cursor % fieldhandle % attlist_tail % next
!call mpas_log_write('Extending att tail for '//trim(attName))
         end if
      else
         if (.not. associated(handle % attlist_head)) then
            handle % attlist_head => new_att_node
!call mpas_log_write('Assigning att head for '//trim(attName))
         end if
         if (.not. associated(handle % attlist_tail)) then
            handle % attlist_tail => new_att_node
!call mpas_log_write('Assigning att tail for '//trim(attName))
         else
            handle % attlist_tail % next => new_att_node
            handle % attlist_tail => handle % attlist_tail % next
!call mpas_log_write('Extending att tail for '//trim(attName))
         end if
      end if

   end subroutine MPAS_io_get_att_real1d


   subroutine MPAS_io_get_att_text(handle, attName, attValue, fieldname, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: attName
      character (len=*), intent(out) :: attValue
      character (len=*), intent(in), optional :: fieldname
      integer, intent(out), optional :: ierr

      integer :: pio_ierr
      integer :: varid
      integer :: xtype
#ifdef USE_PIO2
      integer (kind=MPAS_IO_OFFSET_KIND) :: len
#else
      integer :: len
#endif
      type (fieldlist_type), pointer :: field_cursor
      type (attlist_type), pointer :: att_cursor, new_att_node

!      call mpas_log_write('Called MPAS_io_get_att_text()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      ! Sanity checks
      if (.not. handle % initialized) then
         if (present(ierr)) ierr = MPAS_IO_ERR_UNINIT_HANDLE
         return 
      end if


      !
      ! For variable attributes, find the structure for fieldname
      !
      if (present(fieldname)) then
         field_cursor => handle % fieldlist_head
         do while (associated(field_cursor))
            if (trim(fieldname) == trim(field_cursor % fieldhandle % fieldname)) then
               varid = field_cursor % fieldhandle % fieldid
               exit
            end if
            field_cursor => field_cursor % next
         end do
         if (.not. associated(field_cursor)) then
            if (present(ierr)) ierr = MPAS_IO_ERR_UNDEFINED_VAR
            return
         end if

         ! Check whether we have this attribute cached
         att_cursor => field_cursor % fieldhandle % attlist_head
         do while (associated(att_cursor))
            if (trim(att_cursor % atthandle % attName) == trim(attName)) then
               if (att_cursor % atthandle % attType == MPAS_ATT_TEXT) then
!call mpas_log_write('Using cached attribute')
                  attValue = att_cursor % atthandle % attValueText
               else
                  if (present(ierr)) ierr=MPAS_IO_ERR_WRONG_ATT_TYPE
               end if
               return
            end if
            att_cursor => att_cursor % next
         end do

      else

         ! Check whether we have this attribute cached
         att_cursor => handle % attlist_head
         do while (associated(att_cursor))
            if (trim(att_cursor % atthandle % attName) == trim(attName)) then
               if (att_cursor % atthandle % attType == MPAS_ATT_TEXT) then
!call mpas_log_write('Using cached attribute')
                  attValue = att_cursor % atthandle % attValueText
               else
                  if (present(ierr)) ierr=MPAS_IO_ERR_WRONG_ATT_TYPE
               end if
               return
            end if
            att_cursor => att_cursor % next
         end do

         varid = PIO_global
      end if

      ! Query attribute value
      pio_ierr = PIO_inq_att(handle % pio_file, varid, attName, xtype, len)
      if (pio_ierr /= PIO_noerr) then
         if (present(ierr)) ierr = MPAS_IO_ERR_PIO
         return
      end if
      if (xtype /= PIO_char) then
         if (present(ierr)) ierr=MPAS_IO_ERR_WRONG_ATT_TYPE
         return
      end if

      pio_ierr = PIO_get_att(handle % pio_file, varid, attName, attValue)
      if (pio_ierr /= PIO_noerr) then
         if (present(ierr)) ierr = MPAS_IO_ERR_PIO
         return
      end if

      ! Keep attribute for future reference
      allocate(new_att_node)
      nullify(new_att_node % next)
      allocate(new_att_node % atthandle)
      new_att_node % atthandle % attName = attName
      new_att_node % atthandle % attType = MPAS_ATT_TEXT
      new_att_node % atthandle % attValueText = attValue

      if (present(fieldname)) then
         if (.not. associated(field_cursor % fieldhandle % attlist_head)) then
            field_cursor % fieldhandle % attlist_head => new_att_node
!call mpas_log_write('Assigning att head for '//trim(attName))
         end if
         if (.not. associated(field_cursor % fieldhandle % attlist_tail)) then
            field_cursor % fieldhandle % attlist_tail => new_att_node
!call mpas_log_write('Assigning att tail for '//trim(attName))
         else
            field_cursor % fieldhandle % attlist_tail % next => new_att_node
            field_cursor % fieldhandle % attlist_tail => field_cursor % fieldhandle % attlist_tail % next
!call mpas_log_write('Extending att tail for '//trim(attName))
         end if
      else
         if (.not. associated(handle % attlist_head)) then
            handle % attlist_head => new_att_node
!call mpas_log_write('Assigning att head for '//trim(attName))
         end if
         if (.not. associated(handle % attlist_tail)) then
            handle % attlist_tail => new_att_node
!call mpas_log_write('Assigning att tail for '//trim(attName))
         else
            handle % attlist_tail % next => new_att_node
            handle % attlist_tail => handle % attlist_tail % next
!call mpas_log_write('Extending att tail for '//trim(attName))
         end if
      end if

   end subroutine MPAS_io_get_att_text


   subroutine MPAS_io_put_att_int0d(handle, attName, attValue, fieldname, syncVal, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: attName
      integer, intent(in) :: attValue
      character (len=*), intent(in), optional :: fieldname
      logical, intent(in), optional :: syncVal
      integer, intent(out), optional :: ierr

      integer :: pio_ierr
      integer :: varid
      integer :: attValueLocal
      type (fieldlist_type), pointer :: field_cursor
      type (attlist_type), pointer :: attlist_cursor, new_attlist_node

      attValueLocal = attValue

!      call mpas_log_write('Called MPAS_io_put_att_int0d()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      ! Sanity checks
      if (.not. handle % initialized) then
         if (present(ierr)) ierr = MPAS_IO_ERR_UNINIT_HANDLE
         return 
      end if
      if (handle % data_mode) then
         if (present(ierr)) ierr = MPAS_IO_ERR_DATA_MODE
         return 
      end if
      if (handle % iomode /= MPAS_IO_WRITE) then
         if (present(ierr)) ierr = MPAS_IO_ERR_NOWRITE
         return 
      end if


      allocate(new_attlist_node) 
      nullify(new_attlist_node % next)
      allocate(new_attlist_node % attHandle)
      new_attlist_node % attHandle % attName = attName
      new_attlist_node % attHandle % attType = MPAS_ATT_INT
      new_attlist_node % attHandle % attValueInt = attValue


      !
      ! For variable attributes, find the structure for fieldname
      !
      if (present(fieldname)) then
         field_cursor => handle % fieldlist_head
         do while (associated(field_cursor))
            if (trim(fieldname) == trim(field_cursor % fieldhandle % fieldname)) then

               ! Check whether attribute was already defined
               attlist_cursor => field_cursor % fieldhandle % attlist_head
               do while (associated(attlist_cursor))
                  if (trim(attName) == trim(attlist_cursor % atthandle % attName)) then
!call mpas_log_write('Attribute already defined')
                     if (attlist_cursor % atthandle % attType /= MPAS_ATT_INT .or. &
                         attlist_cursor % atthandle % attValueInt /= attValue) then
                        if (present(ierr)) ierr = MPAS_IO_ERR_REDEF_ATT
                        deallocate(new_attlist_node % attHandle)
                        deallocate(new_attlist_node) 
                     end if
                     return
                  end if
                  attlist_cursor => attlist_cursor % next
               end do

               varid = field_cursor % fieldhandle % fieldid
               exit
            end if
            field_cursor => field_cursor % next
         end do
         if (.not. associated(field_cursor)) then
            if (present(ierr)) ierr = MPAS_IO_ERR_UNDEFINED_VAR
            deallocate(new_attlist_node % attHandle)
            deallocate(new_attlist_node) 
            return
         end if

         ! Add attribute to field attribute list
         if (.not. associated(field_cursor % fieldhandle % attlist_head)) then
            field_cursor % fieldhandle % attlist_head => new_attlist_node
!call mpas_log_write('Assigning att head for '//trim(attname))
         end if
         if (.not. associated(field_cursor % fieldhandle % attlist_tail)) then
            field_cursor % fieldhandle % attlist_tail => new_attlist_node
!call mpas_log_write('Assigning att tail for '//trim(attname))
         else
            field_cursor % fieldhandle % attlist_tail % next => new_attlist_node
            field_cursor % fieldhandle % attlist_tail => field_cursor % fieldhandle % attlist_tail % next
!call mpas_log_write('Extending att tail for '//trim(attname))
         end if

      else

         ! Check whether attribute was already defined
         attlist_cursor => handle % attlist_head
         do while (associated(attlist_cursor))
            if (trim(attName) == trim(attlist_cursor % atthandle % attName)) then
!call mpas_log_write('Attribute already defined')
               if (attlist_cursor % atthandle % attType /= MPAS_ATT_INT .or. &
                   attlist_cursor % atthandle % attValueInt /= attValue) then
                  if (present(ierr)) ierr = MPAS_IO_ERR_REDEF_ATT
                  deallocate(new_attlist_node % attHandle)
                  deallocate(new_attlist_node) 
               end if
               return
            end if
            attlist_cursor => attlist_cursor % next
         end do

         varid = PIO_global

         ! Add attribute to global attribute list
         if (.not. associated(handle % attlist_head)) then
            handle % attlist_head => new_attlist_node
!call mpas_log_write('Assigning att head for '//trim(attname))
         end if
         if (.not. associated(handle % attlist_tail)) then
            handle % attlist_tail => new_attlist_node
!call mpas_log_write('Assigning att tail for '//trim(attname))
         else
            handle % attlist_tail % next => new_attlist_node
            handle % attlist_tail => handle % attlist_tail % next
!call mpas_log_write('Extending att tail for '//trim(attname))
         end if
      end if

      if ( present(syncVal) ) then
         if ( syncVal ) then
            call mpas_dmpar_bcast_int(handle % iocontext % dminfo, attValueLocal, IO_NODE)
         end if
      end if

      pio_ierr = PIO_put_att(handle % pio_file, varid, attName, attValueLocal) 
      if (pio_ierr /= PIO_noerr) then
         if (present(ierr)) ierr = MPAS_IO_ERR_PIO
         return
      end if

      ! Maybe we should add attribute to list only after a successfull call to PIO?

   end subroutine MPAS_io_put_att_int0d


   subroutine MPAS_io_put_att_int1d(handle, attName, attValue, fieldname, syncVal, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: attName
      integer, dimension(:), intent(in) :: attValue
      character (len=*), intent(in), optional :: fieldname
      logical, intent(in), optional :: syncVal
      integer, intent(out), optional :: ierr

      integer :: pio_ierr
      integer :: varid
      integer, dimension(:), pointer :: attValueLocal
      type (fieldlist_type), pointer :: field_cursor
      type (attlist_type), pointer :: attlist_cursor, new_attlist_node

      allocate(attValueLocal(size(attValue, dim=1)))
      attValueLocal(:) = attValue(:)

!      call mpas_log_write('Called MPAS_io_put_att_int1d()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      ! Sanity checks
      if (.not. handle % initialized) then
         if (present(ierr)) ierr = MPAS_IO_ERR_UNINIT_HANDLE
         return 
      end if
      if (handle % data_mode) then
         if (present(ierr)) ierr = MPAS_IO_ERR_DATA_MODE
         return 
      end if
      if (handle % iomode /= MPAS_IO_WRITE) then
         if (present(ierr)) ierr = MPAS_IO_ERR_NOWRITE
         return 
      end if


      allocate(new_attlist_node) 
      nullify(new_attlist_node % next)
      allocate(new_attlist_node % attHandle)
      new_attlist_node % attHandle % attName = attName
      new_attlist_node % attHandle % attType = MPAS_ATT_INTA
      allocate(new_attlist_node % attHandle % attValueIntA(size(attValue)))
      new_attlist_node % attHandle % attValueIntA = attValue


      !
      ! For variable attributes, find the structure for fieldname
      !
      if (present(fieldname)) then
         field_cursor => handle % fieldlist_head
         do while (associated(field_cursor))
            if (trim(fieldname) == trim(field_cursor % fieldhandle % fieldname)) then

               ! Check whether attribute was already defined
               attlist_cursor => field_cursor % fieldhandle % attlist_head
               do while (associated(attlist_cursor))
                  if (trim(attName) == trim(attlist_cursor % atthandle % attName)) then
!call mpas_log_write('Attribute already defined')
                     if (attlist_cursor % atthandle % attType /= MPAS_ATT_INTA .or. &
                         size(attlist_cursor % atthandle % attValueIntA) /= size(attValue)) then
                        if (present(ierr)) ierr = MPAS_IO_ERR_REDEF_ATT
                        deallocate(new_attlist_node % attHandle)
                        deallocate(new_attlist_node) 
!                     else if (attlist_cursor % atthandle % attValueIntA(:) /= attValue(:)) then   ! array sizes should match based on previous if-test
!                        if (present(ierr)) ierr = MPAS_IO_ERR_REDEF_ATT
!                        deallocate(new_attlist_node % attHandle)
!                        deallocate(new_attlist_node) 
                     end if
                     return
                  end if
                  attlist_cursor => attlist_cursor % next
               end do

               varid = field_cursor % fieldhandle % fieldid
               exit
            end if
            field_cursor => field_cursor % next
         end do
         if (.not. associated(field_cursor)) then
            if (present(ierr)) ierr = MPAS_IO_ERR_UNDEFINED_VAR
            deallocate(new_attlist_node % attHandle)
            deallocate(new_attlist_node) 
            return
         end if

         ! Add attribute to field attribute list
         if (.not. associated(field_cursor % fieldhandle % attlist_head)) then
            field_cursor % fieldhandle % attlist_head => new_attlist_node
!call mpas_log_write('Assigning att head for '//trim(attname))
         end if
         if (.not. associated(field_cursor % fieldhandle % attlist_tail)) then
            field_cursor % fieldhandle % attlist_tail => new_attlist_node
!call mpas_log_write('Assigning att tail for '//trim(attname))
         else
            field_cursor % fieldhandle % attlist_tail % next => new_attlist_node
            field_cursor % fieldhandle % attlist_tail => field_cursor % fieldhandle % attlist_tail % next
!call mpas_log_write('Extending att tail for '//trim(attname))
         end if

      else

         ! Check whether attribute was already defined
         attlist_cursor => handle % attlist_head
         do while (associated(attlist_cursor))
            if (trim(attName) == trim(attlist_cursor % atthandle % attName)) then
!call mpas_log_write('Attribute already defined')
               if (attlist_cursor % atthandle % attType /= MPAS_ATT_INTA .or. &
                   size(attlist_cursor % atthandle % attValueIntA) /= size(attValue)) then
                  if (present(ierr)) ierr = MPAS_IO_ERR_REDEF_ATT
                  deallocate(new_attlist_node % attHandle)
                  deallocate(new_attlist_node) 
!               else if (attlist_cursor % atthandle % attValueIntA /= attValue) then
!               else if (attlist_cursor % atthandle % attValueIntA /= attValue) then
!               else if (attlist_cursor % atthandle % attValueIntA /= attValue) then
!               else if (attlist_cursor % atthandle % attValueIntA /= attValue) then
               end if
               return
            end if
            attlist_cursor => attlist_cursor % next
         end do

         varid = PIO_global

         ! Add attribute to global attribute list
         if (.not. associated(handle % attlist_head)) then
            handle % attlist_head => new_attlist_node
!call mpas_log_write('Assigning att head for '//trim(attname))
         end if
         if (.not. associated(handle % attlist_tail)) then
            handle % attlist_tail => new_attlist_node
!call mpas_log_write('Assigning att tail for '//trim(attname))
         else
            handle % attlist_tail % next => new_attlist_node
            handle % attlist_tail => handle % attlist_tail % next
!call mpas_log_write('Extending att tail for '//trim(attname))
         end if
      end if

      if ( present(syncVal) ) then
         if ( syncVal ) then
            call mpas_dmpar_bcast_ints(handle % iocontext % dminfo, size(attValueLocal, dim=1), attValueLocal, IO_NODE)
         end if
      end if

      pio_ierr = PIO_put_att(handle % pio_file, varid, attName, attValueLocal) 
      if (pio_ierr /= PIO_noerr) then
         if (present(ierr)) ierr = MPAS_IO_ERR_PIO
         return
      end if

      deallocate(attValueLocal)

      ! Maybe we should add attribute to list only after a successfull call to PIO?

   end subroutine MPAS_io_put_att_int1d


   subroutine MPAS_io_put_att_real0d(handle, attName, attValue, fieldname, syncVal, precision, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: attName
      real (kind=RKIND), intent(in) :: attValue
      character (len=*), intent(in), optional :: fieldname
      logical, intent(in), optional :: syncVal
      integer, intent(in), optional :: precision
      integer, intent(out), optional :: ierr

      integer :: pio_ierr
      integer :: varid
      real (kind=RKIND) :: attValueLocal
      real (kind=R4KIND) :: singleVal
      real (kind=R8KIND) :: doubleVal
      type (fieldlist_type), pointer :: field_cursor
      type (attlist_type), pointer :: attlist_cursor, new_attlist_node

      attValueLocal = attValue

!      call mpas_log_write('Called MPAS_io_put_att_real0d()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      ! Sanity checks
      if (.not. handle % initialized) then
         if (present(ierr)) ierr = MPAS_IO_ERR_UNINIT_HANDLE
         return 
      end if
      if (handle % data_mode) then
         if (present(ierr)) ierr = MPAS_IO_ERR_DATA_MODE
         return 
      end if
      if (handle % iomode /= MPAS_IO_WRITE) then
         if (present(ierr)) ierr = MPAS_IO_ERR_NOWRITE
         return 
      end if


      allocate(new_attlist_node) 
      nullify(new_attlist_node % next)
      allocate(new_attlist_node % attHandle)
      new_attlist_node % attHandle % attName = attName
      new_attlist_node % attHandle % attType = MPAS_ATT_REAL
      new_attlist_node % attHandle % attValueReal = attValue

      if (present(precision)) then
         new_attlist_node % attHandle % precision = precision
      else
         new_attlist_node % attHandle % precision = MPAS_IO_NATIVE_PRECISION
      end if

      !
      ! For variable attributes, find the structure for fieldname
      !
      if (present(fieldname)) then
         field_cursor => handle % fieldlist_head
         do while (associated(field_cursor))
            if (trim(fieldname) == trim(field_cursor % fieldhandle % fieldname)) then

               ! Check whether attribute was already defined
               attlist_cursor => field_cursor % fieldhandle % attlist_head
               do while (associated(attlist_cursor))
                  if (trim(attName) == trim(attlist_cursor % atthandle % attName)) then
!call mpas_log_write('Attribute already defined')
                     if (attlist_cursor % atthandle % attType /= MPAS_ATT_REAL .or. &
                         attlist_cursor % atthandle % attValueReal /= attValue) then
                        if (present(ierr)) ierr = MPAS_IO_ERR_REDEF_ATT
                        deallocate(new_attlist_node % attHandle)
                        deallocate(new_attlist_node) 
                     end if
                     return
                  end if
                  attlist_cursor => attlist_cursor % next
               end do

               varid = field_cursor % fieldhandle % fieldid
               exit
            end if
            field_cursor => field_cursor % next
         end do
         if (.not. associated(field_cursor)) then
            if (present(ierr)) ierr = MPAS_IO_ERR_UNDEFINED_VAR
            deallocate(new_attlist_node % attHandle)
            deallocate(new_attlist_node) 
            return
         end if

         ! Add attribute to field attribute list
         if (.not. associated(field_cursor % fieldhandle % attlist_head)) then
            field_cursor % fieldhandle % attlist_head => new_attlist_node
!call mpas_log_write('Assigning att head for '//trim(attname))
         end if
         if (.not. associated(field_cursor % fieldhandle % attlist_tail)) then
            field_cursor % fieldhandle % attlist_tail => new_attlist_node
!call mpas_log_write('Assigning att tail for '//trim(attname))
         else
            field_cursor % fieldhandle % attlist_tail % next => new_attlist_node
            field_cursor % fieldhandle % attlist_tail => field_cursor % fieldhandle % attlist_tail % next
!call mpas_log_write('Extending att tail for '//trim(attname))
         end if

      else

         ! Check whether attribute was already defined
         attlist_cursor => handle % attlist_head
         do while (associated(attlist_cursor))
            if (trim(attName) == trim(attlist_cursor % atthandle % attName)) then
!call mpas_log_write('Attribute already defined')
               if (attlist_cursor % atthandle % attType /= MPAS_ATT_REAL .or. &
                   attlist_cursor % atthandle % attValueReal /= attValue) then
                  if (present(ierr)) ierr = MPAS_IO_ERR_REDEF_ATT
                  deallocate(new_attlist_node % attHandle)
                  deallocate(new_attlist_node) 
               end if
               return
            end if
            attlist_cursor => attlist_cursor % next
         end do

         varid = PIO_global

         ! Add attribute to global attribute list
         if (.not. associated(handle % attlist_head)) then
            handle % attlist_head => new_attlist_node
!call mpas_log_write('Assigning att head for '//trim(attname))
         end if
         if (.not. associated(handle % attlist_tail)) then
            handle % attlist_tail => new_attlist_node
!call mpas_log_write('Assigning att tail for '//trim(attname))
         else
            handle % attlist_tail % next => new_attlist_node
            handle % attlist_tail => handle % attlist_tail % next
!call mpas_log_write('Extending att tail for '//trim(attname))
         end if
      end if

      if ( present(syncVal) ) then
         if ( syncVal ) then
            call mpas_dmpar_bcast_real(handle % iocontext % dminfo, attValueLocal, IO_NODE)
         end if
      end if

      if ((new_attlist_node % attHandle % precision == MPAS_IO_SINGLE_PRECISION) .and. &
          (MPAS_IO_NATIVE_PRECISION /= MPAS_IO_SINGLE_PRECISION)) then
         singleVal = real(attValueLocal,R4KIND)
         pio_ierr = PIO_put_att(handle % pio_file, varid, attName, singleVal) 
      else if ((new_attlist_node % attHandle % precision == MPAS_IO_DOUBLE_PRECISION) .and. &
          (MPAS_IO_NATIVE_PRECISION /= MPAS_IO_DOUBLE_PRECISION)) then
         doubleVal = real(attValueLocal,R8KIND)
         pio_ierr = PIO_put_att(handle % pio_file, varid, attName, doubleVal) 
      else
         pio_ierr = PIO_put_att(handle % pio_file, varid, attName, attValueLocal) 
      end if

      if (pio_ierr /= PIO_noerr) then
         if (present(ierr)) ierr = MPAS_IO_ERR_PIO
         return
      end if

      ! Maybe we should add attribute to list only after a successfull call to PIO?

   end subroutine MPAS_io_put_att_real0d


   subroutine MPAS_io_put_att_real1d(handle, attName, attValue, fieldname, syncVal, precision, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: attName
      real (kind=RKIND), dimension(:), intent(in) :: attValue
      character (len=*), intent(in), optional :: fieldname
      logical, intent(in), optional :: syncVal
      integer, intent(in), optional :: precision
      integer, intent(out), optional :: ierr

      integer :: pio_ierr
      integer :: varid
      real (kind=RKIND), dimension(:), allocatable :: attValueLocal
      real (kind=R4KIND), dimension(:), allocatable :: singleVal
      real (kind=R8KIND), dimension(:), allocatable :: doubleVal
      type (fieldlist_type), pointer :: field_cursor
      type (attlist_type), pointer :: attlist_cursor, new_attlist_node

      allocate(attValueLocal( size(attValue, dim=1) ) )

      attValueLocal(:) = attValue(:)

!      call mpas_log_write('Called MPAS_io_put_att_real1d()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      ! Sanity checks
      if (.not. handle % initialized) then
         if (present(ierr)) ierr = MPAS_IO_ERR_UNINIT_HANDLE
         return 
      end if
      if (handle % data_mode) then
         if (present(ierr)) ierr = MPAS_IO_ERR_DATA_MODE
         return 
      end if
      if (handle % iomode /= MPAS_IO_WRITE) then
         if (present(ierr)) ierr = MPAS_IO_ERR_NOWRITE
         return 
      end if


      allocate(new_attlist_node) 
      nullify(new_attlist_node % next)
      allocate(new_attlist_node % attHandle)
      new_attlist_node % attHandle % attName = attName
      new_attlist_node % attHandle % attType = MPAS_ATT_REALA
      allocate(new_attlist_node % attHandle % attValueRealA(size(attValue)))
      new_attlist_node % attHandle % attValueRealA = attValue

      if (present(precision)) then
         new_attlist_node % attHandle % precision = precision
      else
         new_attlist_node % attHandle % precision = MPAS_IO_NATIVE_PRECISION
      end if

      !
      ! For variable attributes, find the structure for fieldname
      !
      if (present(fieldname)) then
         field_cursor => handle % fieldlist_head
         do while (associated(field_cursor))
            if (trim(fieldname) == trim(field_cursor % fieldhandle % fieldname)) then

               ! Check whether attribute was already defined
               attlist_cursor => field_cursor % fieldhandle % attlist_head
               do while (associated(attlist_cursor))
                  if (trim(attName) == trim(attlist_cursor % atthandle % attName)) then
!call mpas_log_write('Attribute already defined')
                     if (attlist_cursor % atthandle % attType /= MPAS_ATT_REALA .or. &
                         size(attlist_cursor % atthandle % attValueRealA) /= size(attValue)) then
                        if (present(ierr)) ierr = MPAS_IO_ERR_REDEF_ATT
                        deallocate(new_attlist_node % attHandle)
                        deallocate(new_attlist_node) 
!                     else if (attlist_cursor % atthandle % attValueIntA(:) /= attValue(:)) then   ! array sizes should match based on previous if-test
!                        if (present(ierr)) ierr = MPAS_IO_ERR_REDEF_ATT
!                        deallocate(new_attlist_node % attHandle)
!                        deallocate(new_attlist_node) 
                     end if
                     return
                  end if
                  attlist_cursor => attlist_cursor % next
               end do

               varid = field_cursor % fieldhandle % fieldid
               exit
            end if
            field_cursor => field_cursor % next
         end do
         if (.not. associated(field_cursor)) then
            if (present(ierr)) ierr = MPAS_IO_ERR_UNDEFINED_VAR
            deallocate(new_attlist_node % attHandle)
            deallocate(new_attlist_node) 
            return
         end if

         ! Add attribute to field attribute list
         if (.not. associated(field_cursor % fieldhandle % attlist_head)) then
            field_cursor % fieldhandle % attlist_head => new_attlist_node
!call mpas_log_write('Assigning att head for '//trim(attname))
         end if
         if (.not. associated(field_cursor % fieldhandle % attlist_tail)) then
            field_cursor % fieldhandle % attlist_tail => new_attlist_node
!call mpas_log_write('Assigning att tail for '//trim(attname))
         else
            field_cursor % fieldhandle % attlist_tail % next => new_attlist_node
            field_cursor % fieldhandle % attlist_tail => field_cursor % fieldhandle % attlist_tail % next
!call mpas_log_write('Extending att tail for '//trim(attname))
         end if

      else

         ! Check whether attribute was already defined
         attlist_cursor => handle % attlist_head
         do while (associated(attlist_cursor))
            if (trim(attName) == trim(attlist_cursor % atthandle % attName)) then
!call mpas_log_write('Attribute already defined')
               if (attlist_cursor % atthandle % attType /= MPAS_ATT_REALA .or. &
                   size(attlist_cursor % atthandle % attValueRealA) /= size(attValue)) then
                  if (present(ierr)) ierr = MPAS_IO_ERR_REDEF_ATT
                  deallocate(new_attlist_node % attHandle)
                  deallocate(new_attlist_node) 
!               else if (attlist_cursor % atthandle % attValueIntA /= attValue) then
!               else if (attlist_cursor % atthandle % attValueIntA /= attValue) then
!               else if (attlist_cursor % atthandle % attValueIntA /= attValue) then
!               else if (attlist_cursor % atthandle % attValueIntA /= attValue) then
               end if
               return
            end if
            attlist_cursor => attlist_cursor % next
         end do

         varid = PIO_global

         ! Add attribute to global attribute list
         if (.not. associated(handle % attlist_head)) then
            handle % attlist_head => new_attlist_node
!call mpas_log_write('Assigning att head for '//trim(attname))
         end if
         if (.not. associated(handle % attlist_tail)) then
            handle % attlist_tail => new_attlist_node
!call mpas_log_write('Assigning att tail for '//trim(attname))
         else
            handle % attlist_tail % next => new_attlist_node
            handle % attlist_tail => handle % attlist_tail % next
!call mpas_log_write('Extending att tail for '//trim(attname))
         end if
      end if

      if ( present(syncVal) ) then
         if ( syncVal ) then
            call mpas_dmpar_bcast_reals(handle % iocontext % dminfo, size(attValueLocal, dim=1), attValueLocal, IO_NODE)
         end if
      end if

      if ((new_attlist_node % attHandle % precision == MPAS_IO_SINGLE_PRECISION) .and. &
          (MPAS_IO_NATIVE_PRECISION /= MPAS_IO_SINGLE_PRECISION)) then
         allocate(singleVal(size(attValueLocal)))
         singleVal(:) = real(attValueLocal(:),R4KIND)
         pio_ierr = PIO_put_att(handle % pio_file, varid, attName, singleVal) 
         deallocate(singleVal)
      else if ((new_attlist_node % attHandle % precision == MPAS_IO_DOUBLE_PRECISION) .and. &
          (MPAS_IO_NATIVE_PRECISION /= MPAS_IO_DOUBLE_PRECISION)) then
         allocate(doubleVal(size(attValueLocal)))
         doubleVal(:) = real(attValueLocal(:),R8KIND)
         pio_ierr = PIO_put_att(handle % pio_file, varid, attName, doubleVal) 
         deallocate(doubleVal)
      else
         pio_ierr = PIO_put_att(handle % pio_file, varid, attName, attValueLocal) 
      end if
      if (pio_ierr /= PIO_noerr) then
         if (present(ierr)) ierr = MPAS_IO_ERR_PIO
         return
      end if

      deallocate(attValueLocal)

      ! Maybe we should add attribute to list only after a successfull call to PIO?

   end subroutine MPAS_io_put_att_real1d


   subroutine MPAS_io_put_att_text(handle, attName, attValue, fieldname, syncVal, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: attName
      character (len=*), intent(in) :: attValue
      character (len=*), intent(in), optional :: fieldname
      logical, intent(in), optional :: syncVal
      integer, intent(out), optional :: ierr

      integer :: pio_ierr
      integer :: varid
      integer :: valLen
      character (len=StrKind) :: attValueLocal, trimmedVal
      type (fieldlist_type), pointer :: field_cursor
      type (attlist_type), pointer :: attlist_cursor, new_attlist_node

      valLen = len_trim(attValue)
      trimmedVal = trim(attValue)

      if ( valLen > StrKind ) then
         call mpas_log_write('Attribute ''' // trim(attName) // ''' has a value longer than StrKIND. ' // &
             'It will be cut to a length of $i', MPAS_LOG_WARN, intArgs=(/StrKIND/) )
         attValueLocal = trimmedVal(1:StrKIND)
      else
         attValueLocal = trimmedVal
      end if

!      call mpas_log_write('Called MPAS_io_put_att_text()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      ! Sanity checks
      if (.not. handle % initialized) then
         if (present(ierr)) ierr = MPAS_IO_ERR_UNINIT_HANDLE
         return 
      end if
      if (handle % data_mode) then
         if (present(ierr)) ierr = MPAS_IO_ERR_DATA_MODE
         return 
      end if
      if (handle % iomode /= MPAS_IO_WRITE) then
         if (present(ierr)) ierr = MPAS_IO_ERR_NOWRITE
         return 
      end if


      allocate(new_attlist_node) 
      nullify(new_attlist_node % next)
      allocate(new_attlist_node % attHandle)
      new_attlist_node % attHandle % attName = attName
      new_attlist_node % attHandle % attType = MPAS_ATT_TEXT
      new_attlist_node % attHandle % attValueText = attValue


      !
      ! For variable attributes, find the structure for fieldname
      !
      if (present(fieldname)) then
         field_cursor => handle % fieldlist_head
         do while (associated(field_cursor))
            if (trim(fieldname) == trim(field_cursor % fieldhandle % fieldname)) then

               ! Check whether attribute was already defined
               attlist_cursor => field_cursor % fieldhandle % attlist_head
               do while (associated(attlist_cursor))
                  if (trim(attName) == trim(attlist_cursor % atthandle % attName)) then
!call mpas_log_write('Attribute already defined')
                     if (attlist_cursor % atthandle % attType /= MPAS_ATT_TEXT .or. &
                         trim(attlist_cursor % atthandle % attValueText) /= trim(attValue)) then
                        if (present(ierr)) ierr = MPAS_IO_ERR_REDEF_ATT
                        deallocate(new_attlist_node % attHandle)
                        deallocate(new_attlist_node) 
                     end if
                     return
                  end if
                  attlist_cursor => attlist_cursor % next
               end do

               varid = field_cursor % fieldhandle % fieldid
               exit
            end if
            field_cursor => field_cursor % next
         end do
         if (.not. associated(field_cursor)) then
            if (present(ierr)) ierr = MPAS_IO_ERR_UNDEFINED_VAR
            deallocate(new_attlist_node % attHandle)
            deallocate(new_attlist_node) 
            return
         end if

         ! Add attribute to field attribute list
         if (.not. associated(field_cursor % fieldhandle % attlist_head)) then
            field_cursor % fieldhandle % attlist_head => new_attlist_node
!call mpas_log_write('Assigning att head for '//trim(attname))
         end if
         if (.not. associated(field_cursor % fieldhandle % attlist_tail)) then
            field_cursor % fieldhandle % attlist_tail => new_attlist_node
!call mpas_log_write('Assigning att tail for '//trim(attname))
         else
            field_cursor % fieldhandle % attlist_tail % next => new_attlist_node
            field_cursor % fieldhandle % attlist_tail => field_cursor % fieldhandle % attlist_tail % next
!call mpas_log_write('Extending att tail for '//trim(attname))
         end if

      else

         ! Check whether attribute was already defined
         attlist_cursor => handle % attlist_head
         do while (associated(attlist_cursor))
            if (trim(attName) == trim(attlist_cursor % atthandle % attName)) then
!call mpas_log_write('Attribute already defined')
               if (attlist_cursor % atthandle % attType /= MPAS_ATT_TEXT .or. &
                   trim(attlist_cursor % atthandle % attValueText) /= trim(attValue)) then
                  if (present(ierr)) ierr = MPAS_IO_ERR_REDEF_ATT
                  deallocate(new_attlist_node % attHandle)
                  deallocate(new_attlist_node) 
               end if
               return
            end if
            attlist_cursor => attlist_cursor % next
         end do

         varid = PIO_global

         ! Add attribute to global attribute list
         if (.not. associated(handle % attlist_head)) then
            handle % attlist_head => new_attlist_node
!call mpas_log_write('Assigning att head for '//trim(attname))
         end if
         if (.not. associated(handle % attlist_tail)) then
            handle % attlist_tail => new_attlist_node
!call mpas_log_write('Assigning att tail for '//trim(attname))
         else
            handle % attlist_tail % next => new_attlist_node
            handle % attlist_tail => handle % attlist_tail % next
!call mpas_log_write('Extending att tail for '//trim(attname))
         end if
      end if

      if ( present(syncVal) ) then
         if ( syncVal ) then
            call mpas_dmpar_bcast_char(handle % iocontext % dminfo, attValueLocal, IO_NODE)
         end if
      end if

      pio_ierr = PIO_put_att(handle % pio_file, varid, attName, trim(attValueLocal)) 
      if (pio_ierr /= PIO_noerr) then

         if (present(ierr)) ierr = MPAS_IO_ERR_PIO

         !
         ! If we are working with a pre-existing file and the text attribute is larger than in the file, we need 
         ! to enter define mode before writing the attribute. Note the PIO_redef documentation:
         !     'Entering and leaving netcdf define mode causes a file sync operation to occur, 
         !      these operations can be very expensive in parallel systems.'
         !
         if (handle % preexisting_file .and. .not. handle % data_mode) then
            pio_ierr = PIO_redef(handle % pio_file)
            if (pio_ierr /= PIO_noerr) then
               return
            end if

            pio_ierr = PIO_put_att(handle % pio_file, varid, attName, trim(attValueLocal)) 
            if (pio_ierr /= PIO_noerr) then
               return
            end if

            pio_ierr = PIO_enddef(handle % pio_file)
            if (pio_ierr /= PIO_noerr) then
               return
            end if

            if (present(ierr)) ierr = MPAS_IO_NOERR

         end if

         return
      end if

      ! Maybe we should add attribute to list only after a successfull call to PIO?

   end subroutine MPAS_io_put_att_text


   subroutine MPAS_io_set_frame(handle, frame, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      integer, intent(in) :: frame
      integer, intent(out), optional :: ierr

!      call mpas_log_write('Called MPAS_io_set_frame()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      handle % frame_number = frame

   end subroutine MPAS_io_set_frame


   subroutine MPAS_io_advance_frame(handle, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      integer, intent(out), optional :: ierr

!      call mpas_log_write('Called MPAS_io_advance_frame()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      handle % frame_number = handle % frame_number + 1

   end subroutine MPAS_io_advance_frame


   subroutine MPAS_io_sync(handle, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      integer, intent(out), optional :: ierr

!      call mpas_log_write('Called MPAS_io_sync()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      ! Sanity checks
      if (.not. handle % initialized) then
         if (present(ierr)) ierr = MPAS_IO_ERR_UNINIT_HANDLE
         return 
      end if

      call PIO_syncfile(handle % pio_file)

   end subroutine MPAS_io_sync


   subroutine MPAS_io_close(handle, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      integer, intent(out), optional :: ierr

      type (dimlist_type), pointer :: dimlist_ptr, dimlist_del
      type (fieldlist_type), pointer :: fieldlist_ptr, fieldlist_del
      type (attlist_type), pointer :: attlist_ptr, attlist_del

!      call mpas_log_write('Called MPAS_io_close()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      ! Sanity checks
      if (.not. handle % initialized) then
         if (present(ierr)) ierr = MPAS_IO_ERR_UNINIT_HANDLE
         return 
      end if

      ! Deallocate memory associated with the file
      fieldlist_ptr => handle % fieldlist_head
      do while (associated(fieldlist_ptr))
         fieldlist_del => fieldlist_ptr 
         fieldlist_ptr => fieldlist_ptr % next

         attlist_ptr => fieldlist_del % fieldhandle % attlist_head
         do while (associated(attlist_ptr))
            attlist_del => attlist_ptr 
            attlist_ptr => attlist_ptr % next
            if (attlist_del % atthandle % attType == MPAS_ATT_INTA) deallocate(attlist_del % atthandle % attValueIntA)
            if (attlist_del % atthandle % attType == MPAS_ATT_REALA) deallocate(attlist_del % atthandle % attValueRealA)
            deallocate(attlist_del % atthandle)
         end do
         nullify(fieldlist_del % fieldhandle % attlist_head)
         nullify(fieldlist_del % fieldhandle % attlist_tail)

         deallocate(fieldlist_del % fieldhandle % dims)

         deallocate(fieldlist_del % fieldhandle)
      end do
      nullify(handle % fieldlist_head)
      nullify(handle % fieldlist_tail)

      dimlist_ptr => handle % dimlist_head
      do while (associated(dimlist_ptr))
         dimlist_del => dimlist_ptr 
         dimlist_ptr => dimlist_ptr % next
         deallocate(dimlist_del % dimhandle)
      end do
      nullify(handle % dimlist_head)
      nullify(handle % dimlist_tail)

      attlist_ptr => handle % attlist_head
      do while (associated(attlist_ptr))
         attlist_del => attlist_ptr 
         attlist_ptr => attlist_ptr % next
         if (attlist_del % atthandle % attType == MPAS_ATT_INTA) deallocate(attlist_del % atthandle % attValueIntA)
         if (attlist_del % atthandle % attType == MPAS_ATT_REALA) deallocate(attlist_del % atthandle % attValueRealA)
         deallocate(attlist_del % atthandle)
      end do
      nullify(handle % attlist_head)
      nullify(handle % attlist_tail)

      handle % initialized = .false.

!call mpas_log_write('MGD PIO_closefile')
      call PIO_closefile(handle % pio_file)

   end subroutine MPAS_io_close


   subroutine MPAS_io_finalize(ioContext, finalize_iosystem, ierr)

      implicit none

      type (mpas_io_context_type), pointer :: ioContext
      logical, optional :: finalize_iosystem
      integer, intent(out), optional :: ierr

      integer :: pio_ierr
      type (decomplist_type), pointer :: decomp_cursor, decomp_del

!      call mpas_log_write('Called MPAS_io_finalize()')
      if (present(ierr)) ierr = MPAS_IO_NOERR

      decomp_cursor => ioContext % decomp_list
      do while (associated(decomp_cursor))
         decomp_del => decomp_cursor
         decomp_cursor => decomp_cursor % next
!call mpas_log_write('Deallocating a decomposition...')
!if (.not. associated(decomp_del % decomphandle)) call mpas_log_write('OOPS... do not have decomphandle')
         deallocate(decomp_del % decomphandle % dims)
         deallocate(decomp_del % decomphandle % indices)
         call PIO_freedecomp(ioContext % pio_iosystem, decomp_del % decomphandle % pio_iodesc)
         deallocate(decomp_del % decomphandle)
         deallocate(decomp_del)
      end do

!call mpas_log_write('MGD PIO_finalize')
      if (present(finalize_iosystem)) then
         if ( finalize_iosystem ) then
           call PIO_finalize(ioContext % pio_iosystem, pio_ierr)
           if (pio_ierr /= PIO_noerr) then
              if (present(ierr)) ierr = MPAS_IO_ERR_PIO
              return
           end if
           deallocate(ioContext % pio_iosystem)
         end if
      end if

   end subroutine MPAS_io_finalize


   type (dm_info) function MPAS_io_handle_dminfo(handle)

      implicit none

      type (MPAS_IO_Handle_type), intent(in) :: handle

      MPAS_io_handle_dminfo = handle % ioContext % dminfo

   end function MPAS_io_handle_dminfo


   subroutine MPAS_io_err_mesg(ierr, fatal)

      implicit none

      integer, intent(in) :: ierr
      logical, intent(in) :: fatal

      select case (ierr)
         case (MPAS_IO_NOERR)
            ! ... do nothing ...
         case (MPAS_IO_ERR_INVALID_MODE)
            call mpas_log_write('MPAS IO Error: Invalid file access mode', MPAS_LOG_ERR)
         case (MPAS_IO_ERR_INVALID_FORMAT)
            call mpas_log_write('MPAS IO Error: Invalid I/O format', MPAS_LOG_ERR)
         case (MPAS_IO_ERR_LONG_FILENAME)
            call mpas_log_write('MPAS IO Error: Filename too long', MPAS_LOG_ERR)
         case (MPAS_IO_ERR_UNINIT_HANDLE)
            call mpas_log_write('MPAS IO Error: Uninitialized I/O handle', MPAS_LOG_ERR)
         case (MPAS_IO_ERR_PIO)
            call mpas_log_write('MPAS IO Error: Bad return value from PIO', MPAS_LOG_ERR)
         case (MPAS_IO_ERR_DATA_MODE)
            call mpas_log_write('MPAS IO Error: Cannot define in data mode', MPAS_LOG_ERR)
         case (MPAS_IO_ERR_NOWRITE)
            call mpas_log_write('MPAS IO Error: File not opened for writing', MPAS_LOG_ERR)
         case (MPAS_IO_ERR_REDEF_DIM)
            call mpas_log_write('MPAS IO Error: Inconsistent redefinition of dimension', MPAS_LOG_ERR)
         case (MPAS_IO_ERR_REDEF_VAR)
            call mpas_log_write('MPAS IO Error: Inconsistent redefinition of field', MPAS_LOG_ERR)
         case (MPAS_IO_ERR_UNDEFINED_DIM)
            call mpas_log_write('MPAS IO Error: Field uses undefined dimension', MPAS_LOG_ERR)
         case (MPAS_IO_ERR_UNDEFINED_VAR)
            call mpas_log_write('MPAS IO Error: Undefined field', MPAS_LOG_ERR)
         case (MPAS_IO_ERR_REDEF_ATT)
            call mpas_log_write('MPAS IO Error: Inconsistent redefinition of attribute', MPAS_LOG_ERR)
         case (MPAS_IO_ERR_WRONG_ATT_TYPE)
            call mpas_log_write('MPAS IO Error: Wrong type for requested attribute', MPAS_LOG_ERR)
         case (MPAS_IO_ERR_NO_DECOMP)
            call mpas_log_write('MPAS IO Error: Decomposition indices not set for field', MPAS_LOG_ERR)
         case (MPAS_IO_ERR_TWO_UNLIMITED_DIMS)
            call mpas_log_write('MPAS IO Error: Defining more than one unlimited dimension', MPAS_LOG_ERR)
         case (MPAS_IO_ERR_WRONG_MODE)
            call mpas_log_write('MPAS IO Error: Operation not permitted in this file mode', MPAS_LOG_ERR)
         case (MPAS_IO_ERR_NO_UNLIMITED_DIM)
            call mpas_log_write('MPAS IO Error: No unlimited dimension found in dataset', MPAS_LOG_ERR)
         case (MPAS_IO_ERR_UNIMPLEMENTED)
            call mpas_log_write('MPAS IO Error: Unimplemented functionality', MPAS_LOG_ERR)
         case (MPAS_IO_ERR_WOULD_CLOBBER)
            call mpas_log_write('MPAS IO Error: Would clobber existing file', MPAS_LOG_ERR)
         case (MPAS_IO_ERR_NOEXIST_READ)
            call mpas_log_write('MPAS IO Error: Attempting to read a file which does not exist.', MPAS_LOG_ERR)
         case default
            call mpas_log_write('MPAS IO Error: Unrecognized error code...', MPAS_LOG_ERR)
      end select

      if (fatal .and. (ierr /= MPAS_IO_NOERR)) call mpas_log_write('ERROR In MPAS_IO', MPAS_LOG_CRIT)

   end subroutine MPAS_io_err_mesg
 
end module mpas_io
