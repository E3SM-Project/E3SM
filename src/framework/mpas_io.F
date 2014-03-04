! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module mpas_io

   use mpas_dmpar_types
   use mpas_dmpar
   use mpas_io_units

   use pio
   use piolib_mod
   use pionfatt_mod
   use pio_types

#ifdef SINGLE_PRECISION
   integer, parameter :: PIO_REALKIND = PIO_REAL
#else
   integer, parameter :: PIO_REALKIND = PIO_DOUBLE
#endif

   ! File access modes
   integer, parameter :: MPAS_IO_READ  = 1, &
                         MPAS_IO_WRITE = 2

   ! I/O formats
   integer, parameter :: MPAS_IO_NETCDF  = 3, &
                         MPAS_IO_PNETCDF = 4

   ! Field and attribute types
   integer, parameter :: MPAS_IO_REAL     = 5,  &
                         MPAS_IO_DOUBLE   = 6,  &
                         MPAS_IO_INT      = 7,  &
                         MPAS_IO_LOGICAL  = 8,  &
                         MPAS_IO_CHAR     = 9

   ! Unlimited / record dimension
   integer, parameter :: MPAS_IO_UNLIMITED_DIM = -123456

   ! Error codes
   integer, parameter :: MPAS_IO_NOERR              =  0, &
                         MPAS_IO_ERR_INVALID_MODE   = -1, &
                         MPAS_IO_ERR_INVALID_FORMAT = -2, &
                         MPAS_IO_ERR_LONG_FILENAME  = -3, &
                         MPAS_IO_ERR_UNINIT_HANDLE  = -4, &
                         MPAS_IO_ERR_PIO            = -5, &
                         MPAS_IO_ERR_DATA_MODE      = -6, &
                         MPAS_IO_ERR_NOWRITE        = -7, &
                         MPAS_IO_ERR_REDEF_DIM      = -8, &
                         MPAS_IO_ERR_REDEF_VAR      = -9, &
                         MPAS_IO_ERR_UNDEFINED_DIM  = -10, &
                         MPAS_IO_ERR_UNDEFINED_VAR  = -11, &
                         MPAS_IO_ERR_REDEF_ATT      = -12, &
                         MPAS_IO_ERR_WRONG_ATT_TYPE = -13, &
                         MPAS_IO_ERR_NO_DECOMP      = -14, &
                         MPAS_IO_ERR_TWO_UNLIMITED_DIMS = -15, &
                         MPAS_IO_ERR_WRONG_MODE         = -16, &
                         MPAS_IO_ERR_NO_UNLIMITED_DIM   = -17, &
                         MPAS_IO_ERR_UNIMPLEMENTED      = -18


   type MPAS_IO_Handle_type
      logical :: initialized = .false.
      logical :: data_mode = .false.
      type (file_desc_t) :: pio_file
      character (len=StrKIND) :: filename
      integer :: iomode
      integer :: ioformat
      integer :: pio_unlimited_dimid
      integer (kind=PIO_offset) :: frame_number = 1
      type (dimlist_type), pointer :: dimlist_head => null()
      type (dimlist_type), pointer :: dimlist_tail => null()
      type (fieldlist_type), pointer :: fieldlist_head => null()
      type (fieldlist_type), pointer :: fieldlist_tail => null()
      type (attlist_type), pointer :: attlist_head => null()
      type (attlist_type), pointer :: attlist_tail => null()
   end type MPAS_IO_Handle_type

   
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


   !!!!!!!! PRIVATE !!!!!!!!

!   integer, parameter :: ATT_INT   = 1
!   integer, parameter :: ATT_INTA  = 2
!   integer, parameter :: ATT_REAL  = 3
!   integer, parameter :: ATT_REALA = 4
!   integer, parameter :: ATT_TEXT  = 5

   type decomphandle_type
      integer :: field_type
      integer, dimension(:), pointer :: dims
      integer, dimension(:), pointer :: indices
      type (io_desc_t) :: pio_iodesc
   end type decomphandle_type

   type atthandle_type
      character (len=StrKIND) :: attName
      integer :: attType
      integer :: attValueInt
      integer, dimension(:), pointer :: attValueIntA => null()
      real (kind=RKIND) :: attValueReal
      real (kind=RKIND), dimension(:), pointer :: attValueRealA => null()
      character (len=StrKIND) :: attValueText
   end type atthandle_type

   type dimhandle_type
      character (len=StrKIND) :: dimname
      logical :: is_unlimited_dim = .false.
      integer :: dimsize
      integer :: dimid
   end type dimhandle_type

   type fieldhandle_type
      character (len=StrKIND) :: fieldname
      integer :: fieldid
      type (Var_desc_t) :: field_desc
      integer :: field_type
      logical :: has_unlimited_dim = .false.
      integer :: ndims
      type (dimhandle_type), pointer, dimension(:) :: dims
      type (attlist_type), pointer :: attlist_head => null()
      type (attlist_type), pointer :: attlist_tail => null()
      type (decomphandle_type), pointer :: decomp => null()
   end type fieldhandle_type

   type decomplist_type
      type (decomphandle_type), pointer :: decomphandle
      type (decomplist_type), pointer :: next => null()
   end type decomplist_type

   type attlist_type
      type (atthandle_type), pointer :: atthandle
      type (attlist_type), pointer :: next => null()
   end type attlist_type

   type dimlist_type
      type (dimhandle_type), pointer :: dimhandle
      type (dimlist_type), pointer :: next => null()
   end type dimlist_type

   type fieldlist_type
      type (fieldhandle_type), pointer :: fieldhandle
      type (fieldlist_type), pointer :: next => null()
   end type fieldlist_type

   type (iosystem_desc_t), pointer, private, save :: pio_iosystem
   type (decomplist_type), pointer, private :: decomp_list => null()
   type (dm_info), private :: local_dminfo
   integer, private:: master_pio_iotype = -999
   

   contains

   subroutine MPAS_io_init(dminfo, io_task_count, io_task_stride, io_system, ierr)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: io_task_count
      integer, intent(in) :: io_task_stride
      type (iosystem_desc_t), optional, pointer :: io_system
      integer, intent(out), optional :: ierr

!      write(stderrUnit,*) 'Called MPAS_io_init()'
      if (present(ierr)) ierr = MPAS_IO_NOERR

      local_dminfo = dminfo

      if(present(io_system)) then
        pio_iosystem => io_system
      else
!write(stderrUnit,*) 'MGD PIO_init'
        allocate(pio_iosystem)
        call PIO_init(local_dminfo % my_proc_id, &     ! comp_rank
                      local_dminfo % comm,       &     ! comp_comm
                      io_task_count,             &     ! num_iotasks
                      0,                         &     ! num_aggregator
                      io_task_stride,            &     ! stride
                      PIO_rearr_box,             &     ! rearr
                      pio_iosystem)                    ! iosystem
  
      end if

      call pio_seterrorhandling(pio_iosystem, PIO_BCAST_ERROR)

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
   subroutine MPAS_io_set_iotype(io_type_in, ierr)

      implicit none

      integer, intent(in) :: io_type_in
      integer, intent(out), optional :: ierr

      if (present(ierr)) then
         ierr = MPAS_IO_NOERR
      end if

      master_pio_iotype = io_type_in
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
   subroutine MPAS_io_unset_iotype(ierr)

      implicit none

      integer, intent(out), optional :: ierr

      if (present(ierr)) then
         ierr = MPAS_IO_NOERR
      end if

      master_pio_iotype = -999
   end subroutine MPAS_io_unset_iotype
   
   type (MPAS_IO_Handle_type) function MPAS_io_open(filename, mode, ioformat, ierr)

      implicit none

      character (len=*), intent(in) :: filename
      integer, intent(in) :: mode
      integer, intent(in) :: ioformat
      integer, intent(out), optional :: ierr

      integer :: pio_ierr, pio_iotype

!      write(stderrUnit,*) 'Called MPAS_io_open()'
      if (present(ierr)) ierr = MPAS_IO_NOERR


      ! Sanity checks
      if (mode /= MPAS_IO_READ .and. &
          mode /= MPAS_IO_WRITE) then
         if (present(ierr)) ierr = MPAS_IO_ERR_INVALID_MODE
         return 
      end if
      if (ioformat /= MPAS_IO_NETCDF .and. &
          ioformat /= MPAS_IO_PNETCDF) then
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

      if (master_pio_iotype /= -999) then
         pio_iotype = master_pio_iotype
      else
         if (ioformat == MPAS_IO_PNETCDF) then
            pio_iotype = PIO_iotype_pnetcdf
         else
            pio_iotype = PIO_iotype_netcdf
         end if
      end if

      if (mode == MPAS_IO_WRITE) then
!write(stderrUnit,*) 'MGD PIO_createfile'
         pio_ierr = PIO_createfile(pio_iosystem, MPAS_io_open % pio_file, pio_iotype, trim(filename), PIO_64BIT_OFFSET)
      else
!write(stderrUnit,*) 'MGD PIO_openfile'
         pio_ierr = PIO_openfile(pio_iosystem, MPAS_io_open % pio_file, pio_iotype, trim(filename), PIO_nowrite)
      endif
      if (pio_ierr /= PIO_noerr) then
         if (present(ierr)) ierr = MPAS_IO_ERR_PIO
         return
      end if

      if (mode == MPAS_IO_READ) then
!MPAS_io_open % pio_unlimited_dimid = 44
         pio_ierr = PIO_inquire(MPAS_io_open % pio_file, unlimitedDimID=MPAS_io_open % pio_unlimited_dimid)
!write(stderrUnit,*) 'Found unlimited dim ', MPAS_io_open % pio_unlimited_dimid
         if (pio_ierr /= PIO_noerr) then
            if (present(ierr)) ierr = MPAS_IO_ERR_PIO
            return
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

!      write(stderrUnit,*) 'Called MPAS_io_inq_unlimited_dim()'
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

!      write(stderrUnit,*) 'Called MPAS_io_inq_dim()'
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
         if (present(ierr)) ierr = MPAS_IO_ERR_PIO
         deallocate(new_dimlist_node % dimhandle)
         deallocate(new_dimlist_node)
         write(stderrUnit,*) 'WARNING: Dimension ', trim(dimname), ' not in input file.'
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
      type (dimlist_type), pointer :: new_dimlist_node
      type (dimlist_type), pointer :: dim_cursor

!      write(stderrUnit,*) 'Called MPAS_io_def_dim()'
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
            if (dimsize /= dim_cursor % dimhandle % dimsize) then
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
!write(stderrUnit,*) 'Assigning head for '//trim(dimname)
      end if
      if (.not. associated(handle % dimlist_tail)) then
         handle % dimlist_tail => new_dimlist_node
!write(stderrUnit,*) 'Assigning tail for '//trim(dimname)
      else
         handle % dimlist_tail % next => new_dimlist_node
         handle % dimlist_tail => handle % dimlist_tail % next
!write(stderrUnit,*) 'Extending tail for '//trim(dimname)
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

!      write(stderrUnit,*) 'Called MPAS_io_inq_var()'
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
!write(stderrUnit,*) 'Already found variable in fieldlist'
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
            write(stderrUnit,*) 'WARNING: Variable ', trim(fieldname), ' not in input file.'
            return
         end if
!write(stderrUnit,*) 'Inquired about variable ID', new_fieldlist_node % fieldhandle % fieldid

         ! Get field type
         pio_ierr = PIO_inq_vartype(handle % pio_file, new_fieldlist_node % fieldhandle % fieldid, new_fieldlist_node % fieldhandle % field_type)
         if (pio_ierr /= PIO_noerr) then
            if (present(ierr)) ierr = MPAS_IO_ERR_PIO
            deallocate(new_fieldlist_node % fieldhandle)
            deallocate(new_fieldlist_node)
            return
         end if
!write(stderrUnit,*) 'Inquired about variable type', new_fieldlist_node % fieldhandle % field_type

         ! Convert to MPAS type
         if (new_fieldlist_node % fieldhandle % field_type == PIO_double) then
            new_fieldlist_node % fieldhandle % field_type = MPAS_IO_DOUBLE
         else if (new_fieldlist_node % fieldhandle % field_type == PIO_real) then
            new_fieldlist_node % fieldhandle % field_type = MPAS_IO_REAL
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
!write(stderrUnit,*) 'Inquired about number of dimensions ', new_fieldlist_node % fieldhandle % ndims

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
!write(stderrUnit,*) 'Inquired about dimension IDs ', dimids
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
!write(stderrUnit,*) 'Inquired about dimension size ', new_fieldlist_node % fieldhandle % dims(i) % dimsize

            pio_ierr = PIO_inq_dimname(handle % pio_file, dimids(i), new_fieldlist_node % fieldhandle % dims(i) % dimname)
            if (pio_ierr /= PIO_noerr) then
               if (present(ierr)) ierr = MPAS_IO_ERR_PIO
               deallocate(new_fieldlist_node % fieldhandle)
               deallocate(new_fieldlist_node)
               deallocate(dimids)
               return
            end if
!write(stderrUnit,*) 'Inquired about dimension name ', trim(new_fieldlist_node % fieldhandle % dims(i) % dimname)

         end do

         deallocate(dimids)

         ! Keep variable information for future reference
         if (.not. associated(handle % fieldlist_head)) then
            handle % fieldlist_head => new_fieldlist_node
!write(stderrUnit,*) 'Assigning head for '//trim(fieldname)
         end if
         if (.not. associated(handle % fieldlist_tail)) then
            handle % fieldlist_tail => new_fieldlist_node
!write(stderrUnit,*) 'Assigning tail for '//trim(fieldname)
         else
            handle % fieldlist_tail % next => new_fieldlist_node
            handle % fieldlist_tail => handle % fieldlist_tail % next
!write(stderrUnit,*) 'Extending tail for '//trim(fieldname)
         end if

         ! Keep dimension information for any new dimensions that were encountered
         do i=1,new_fieldlist_node % fieldhandle % ndims
            found = .false.
            dim_cursor => handle % dimlist_head
            do while (associated(dim_cursor))
               if (trim(dim_cursor % dimhandle % dimname) == trim(new_fieldlist_node % fieldhandle % dims(i) % dimname)) then
!write(stderrUnit,*) 'Already have dimension '//trim(new_fieldlist_node % fieldhandle % dims(i) % dimname)//' in our list...'
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
!write(stderrUnit,*) 'Assigning head for '//trim(new_dimlist_node % dimhandle % dimname)
               end if
               if (.not. associated(handle % dimlist_tail)) then
                  handle % dimlist_tail => new_dimlist_node
!write(stderrUnit,*) 'Assigning tail for '//trim(new_dimlist_node % dimhandle % dimname)
               else
                  handle % dimlist_tail % next => new_dimlist_node
                  handle % dimlist_tail => handle % dimlist_tail % next
!write(stderrUnit,*) 'Extending tail for '//trim(new_dimlist_node % dimhandle % dimname)
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


   subroutine MPAS_io_def_var(handle, fieldname, fieldtype, dimnames, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      integer, intent(in) :: fieldtype
      character (len=StrKIND), dimension(:), intent(in) :: dimnames
      integer, intent(out), optional :: ierr

      integer :: i
      integer :: pio_ierr
      integer :: pio_type
      integer :: ndims
      type (fieldlist_type), pointer :: new_fieldlist_node
      type (fieldlist_type), pointer :: field_cursor
      type (dimlist_type), pointer :: dim_cursor
      integer, dimension(:), pointer :: dimids

!      write(stderrUnit,*) 'Called MPAS_io_def_var()'
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
      ! Check whether this field has already been defined
      !
      ndims = size(dimnames)
!write(stderrUnit,*) 'Defining variable with ',ndims,' dimensions'
      field_cursor => handle % fieldlist_head
      do while (associated(field_cursor))
         if (trim(fieldname) == trim(field_cursor % fieldhandle % fieldname)) then
            if (ndims /= field_cursor % fieldhandle % ndims) then
               if (present(ierr)) ierr = MPAS_IO_ERR_REDEF_VAR
!               write(stderrUnit,*) 'Error: Field '//trim(fieldname)//' previously defined with conflicting number of dimensions: ', &
!                           ndims, field_cursor % fieldhandle % ndims
            end if
            if (fieldtype /= field_cursor % fieldhandle % field_type) then
               if (present(ierr)) ierr = MPAS_IO_ERR_REDEF_VAR
!               write(stderrUnit,*) 'Error: Field '//trim(fieldname)//' previously defined with conflicting type: ', &
!                           fieldtype, field_cursor % fieldhandle % field_type
            end if
            return
         end if
         field_cursor => field_cursor % next
      end do

      !
      ! Otherwise, define it
      !
      allocate(new_fieldlist_node)
      nullify(new_fieldlist_node % next)
      allocate(new_fieldlist_node % fieldhandle)

      new_fieldlist_node % fieldhandle % fieldname = fieldname
      new_fieldlist_node % fieldhandle % field_type = fieldtype
      new_fieldlist_node % fieldhandle % ndims = ndims

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
!write(stderrUnit,*) 'Found dimension '//trim(new_fieldlist_node % fieldhandle % dims(i) % dimname)//' for field '//trim(fieldname)
         else
            if (present(ierr)) ierr = MPAS_IO_ERR_UNDEFINED_DIM
            deallocate(new_fieldlist_node % fieldhandle % dims)
            deallocate(new_fieldlist_node % fieldhandle)
            deallocate(new_fieldlist_node)
            deallocate(dimids)
            return
!            write(stderrUnit,*) 'Error finding dimension '//trim(dimnames(i))//' for field '//trim(fieldname)
         end if
      end do

      ! Convert from MPAS type
      if (new_fieldlist_node % fieldhandle % field_type == MPAS_IO_DOUBLE) then
         pio_type = PIO_double
      else if (new_fieldlist_node % fieldhandle % field_type == MPAS_IO_REAL) then
         pio_type = PIO_real
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
!write(stderrUnit,*) 'Assigning head for '//trim(fieldname)
      end if
      if (.not. associated(handle % fieldlist_tail)) then
         handle % fieldlist_tail => new_fieldlist_node
!write(stderrUnit,*) 'Assigning tail for '//trim(fieldname)
      else
         handle % fieldlist_tail % next => new_fieldlist_node
         handle % fieldlist_tail => handle % fieldlist_tail % next
!write(stderrUnit,*) 'Extending tail for '//trim(fieldname)
      end if

   end subroutine MPAS_io_def_var


   subroutine MPAS_io_get_var_indices(handle, fieldname, indices, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(in) :: handle
      character (len=*), intent(in) :: fieldname
      integer, dimension(:), pointer :: indices
      integer, intent(out), optional :: ierr

      type (fieldlist_type), pointer :: field_cursor

!      write(stderrUnit,*) 'Called MPAS_io_get_var_indices()'
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
!write(stderrUnit,*) trim(fieldname), ' has been defined'

      if (.not. associated(field_cursor % fieldhandle % decomp)) then
         if (present(ierr)) ierr = MPAS_IO_ERR_NO_DECOMP
         return
      end if

      allocate(indices(size(field_cursor % fieldhandle % decomp % indices)))
      indices(:) = field_cursor % fieldhandle % decomp % indices(:)

   end subroutine MPAS_io_get_var_indices


   subroutine MPAS_io_set_var_indices(handle, fieldname, indices, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(in) :: handle
      character (len=*), intent(in) :: fieldname
      integer, dimension(:), intent(in) :: indices
      integer, intent(out), optional :: ierr

      type (fieldlist_type), pointer :: field_cursor
      integer :: pio_type
      integer :: ndims
      integer (kind=PIO_OFFSET) :: pd, indx
      integer :: i 
      integer (kind=PIO_OFFSET) :: i1, i2, i3, i4, i5
      integer, dimension(:), pointer :: dimlist
      integer (kind=PIO_OFFSET), dimension(:), pointer :: compdof
      type (decomplist_type), pointer :: decomp_cursor, new_decomp

!      write(stderrUnit,*) 'Called MPAS_io_set_var_indices()'
      if (present(ierr)) ierr = MPAS_IO_NOERR

      ! Sanity checks
      if (.not. handle % initialized) then
         if (present(ierr)) ierr = MPAS_IO_ERR_UNINIT_HANDLE
         return 
      end if

!      write(stderrUnit,*) 'Assigning ', size(indices), ' indices for ', trim(fieldname)
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
!write(stderrUnit,*) trim(fieldname), ' has been defined'

      !
      ! If this is a scalar field, just return
      !
      if (field_cursor % fieldhandle % ndims == 0 .or. &
          (field_cursor % fieldhandle % ndims == 1 .and. field_cursor % fieldhandle % has_unlimited_dim) .or. &
          field_cursor % fieldhandle % field_type == MPAS_IO_CHAR) then
!write(stderrUnit,*) 'No need to create a decomposition for a 0d field...'
         return
      end if


      !
      ! Check whether a suitable decomposition already exists
      !
      decomp_cursor => decomp_list
!if (.not. associated(decomp_cursor)) write(stderrUnit,*) 'No existing decompositions to check...'
      DECOMP_LOOP: do while (associated(decomp_cursor))
         if (decomp_cursor % decomphandle % field_type == field_cursor % fieldhandle % field_type) then
         if (size(decomp_cursor % decomphandle % dims) == field_cursor % fieldhandle % ndims) then
!write(stderrUnit,*) 'Number of dimensions matches...'
            do i=1,field_cursor % fieldhandle % ndims
!write(stderrUnit,*) 'Checking dimension ', decomp_cursor % decomphandle % dims(i), field_cursor % fieldhandle % dims(i) % dimsize
               if (decomp_cursor % decomphandle % dims(i) /= field_cursor % fieldhandle % dims(i) % dimsize) then
                  decomp_cursor => decomp_cursor % next
                  cycle DECOMP_LOOP
               end if
            end do

            if (size(decomp_cursor % decomphandle % indices) /= size(indices)) then
!write(stderrUnit,*) 'We do not have the same number of indices in this decomposition...'
               decomp_cursor => decomp_cursor % next
               cycle DECOMP_LOOP
            end if

            do i=1,size(decomp_cursor % decomphandle % indices)
               if (indices(i) /= decomp_cursor % decomphandle % indices(i)) then
!write(stderrUnit,*) 'One of the indices does not match... ', i
                  decomp_cursor => decomp_cursor % next
                  cycle DECOMP_LOOP
               end if
            end do 
            
            ! OK, we have a match... just use this decomposition for the field and return
            field_cursor % fieldhandle % decomp => decomp_cursor % decomphandle 
!write(stderrUnit,*) 'Found a matching decomposition that we can use'
            return
         else if ((size(decomp_cursor % decomphandle % dims) == field_cursor % fieldhandle % ndims - 1)  &
                  .and. field_cursor % fieldhandle % has_unlimited_dim  &
                 ) then
!write(stderrUnit,*) 'Number of non-record dimensions matches...'
            do i=1,field_cursor % fieldhandle % ndims
               if (field_cursor % fieldhandle % dims(i) % is_unlimited_dim) cycle
!write(stderrUnit,*) 'Checking dimension ', decomp_cursor % decomphandle % dims(i), field_cursor % fieldhandle % dims(i) % dimsize
               if (decomp_cursor % decomphandle % dims(i) /= field_cursor % fieldhandle % dims(i) % dimsize) then
                  decomp_cursor => decomp_cursor % next
                  cycle DECOMP_LOOP
               end if
            end do

            ! Type and dimensions match... what about indices?
            
            ! OK, we have a match... just use this decomposition for the field and return
            field_cursor % fieldhandle % decomp => decomp_cursor % decomphandle 
!write(stderrUnit,*) 'Found a matching decomposition that we can use (aside from record dimension)'
            return
         end if
         end if
         decomp_cursor => decomp_cursor % next
      end do DECOMP_LOOP

!write(stderrUnit,*) 'Creating a new decomposition'


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
         pd = pd * int(dimlist(i),PIO_OFFSET)
      end do
      new_decomp % decomphandle % dims(ndims) = field_cursor % fieldhandle % dims(ndims) % dimsize
      dimlist(ndims) = size(indices)
      pd = pd * int(dimlist(ndims),PIO_OFFSET)

      allocate(compdof(pd)) 

      indx = 1
      if (ndims == 5) then
         do i5=1,dimlist(5)
         do i4=1,dimlist(4)
         do i3=1,dimlist(3)
         do i2=1,dimlist(2)
         do i1=1,dimlist(1)
            compdof(indx) = i1 &
                          + (i2-1)*int(dimlist(1),PIO_OFFSET) &
                          + (i3-1)*int(dimlist(2),PIO_OFFSET)*int(dimlist(1),PIO_OFFSET) &
                          + (i4-1)*int(dimlist(3),PIO_OFFSET)*int(dimlist(2),PIO_OFFSET)*int(dimlist(1),PIO_OFFSET) &
                          + int(indices(i5)-1,PIO_OFFSET)*int(dimlist(4),PIO_OFFSET)*int(dimlist(3),PIO_OFFSET)*int(dimlist(2),PIO_OFFSET)*int(dimlist(1),PIO_OFFSET)
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
                          + (i2-1)*int(dimlist(1),PIO_OFFSET) &
                          + (i3-1)*int(dimlist(2),PIO_OFFSET)*int(dimlist(1),PIO_OFFSET) &
                          + int(indices(i4)-1,PIO_OFFSET)*int(dimlist(3),PIO_OFFSET)*int(dimlist(2),PIO_OFFSET)*int(dimlist(1),PIO_OFFSET)
            indx = indx + 1
         end do
         end do
         end do
         end do
      else if (ndims == 3) then
         do i3=1,dimlist(3)
         do i2=1,dimlist(2)
         do i1=1,dimlist(1)
            compdof(indx) = i1 + (i2-1)*int(dimlist(1),PIO_OFFSET) + int(indices(i3)-1,PIO_OFFSET)*int(dimlist(2),PIO_OFFSET)*int(dimlist(1),PIO_OFFSET)
            indx = indx + 1
         end do
         end do
         end do
      else if (ndims == 2) then
         do i2=1,dimlist(2)
         do i1=1,dimlist(1)
            compdof(indx) = i1 + int(indices(i2)-1,PIO_OFFSET)*int(dimlist(1),PIO_OFFSET)
            indx = indx + 1
         end do
         end do
      else if (ndims == 1) then
         do i1=1,dimlist(1)
            compdof(indx) = int(indices(i1),PIO_OFFSET)
            indx = indx + 1
         end do
      end if

      dimlist(ndims) = field_cursor % fieldhandle % dims(ndims) % dimsize
      call PIO_initdecomp(pio_iosystem, pio_type, dimlist, compdof, new_decomp % decomphandle % pio_iodesc)

      ! Add new decomposition to the list
      if (.not. associated(decomp_list)) then
         decomp_list => new_decomp
!write(stderrUnit,*) 'Adding first item to the decomp_list'
      else
         new_decomp % next => decomp_list
         decomp_list => new_decomp
!write(stderrUnit,*) 'Adding new decomp to the head of the list'
      end if

!write(stderrUnit,*) 'Setting decomp in fieldhandle'
      field_cursor % fieldhandle % decomp => new_decomp % decomphandle

      deallocate(compdof)
      deallocate(dimlist)
!write(stderrUnit,*) 'All finished.'

   end subroutine MPAS_io_set_var_indices


   subroutine MPAS_io_get_var_generic(handle, fieldname, intVal, intArray1d, intArray2d, intArray3d, intArray4d, &
                                                        realVal, realArray1d, realArray2d, realArray3d, realArray4d, realArray5d, &
                                                        charVal, ierr)

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
      character (len=StrKIND), dimension(1) :: tempchar
      type (fieldlist_type), pointer :: field_cursor

      ! Sanity checks
      if (.not. handle % initialized) then
         if (present(ierr)) ierr = MPAS_IO_ERR_UNINIT_HANDLE
         return 
      end if

!      write(stderrUnit,*) 'Reading ', trim(fieldname)

      !
      ! Check whether the field has been defined
      !
!      write(stderrUnit,*) 'Checking if field is define'
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

!      write(stderrUnit,*) 'Checking for unlimited dim'
      if (field_cursor % fieldhandle % has_unlimited_dim) then
         call PIO_setframe(field_cursor % fieldhandle % field_desc, handle % frame_number)
         start1(1) = handle % frame_number
         count1(1) = 1
     
         start2(1) = 1
         start2(2) = handle % frame_number
         count2(2) = 1
      end if

!      write(stderrUnit,*) 'Checking for real, int, char, etc'
      if (present(realVal)) then
!         write (0,*) '  value is real'
         if (field_cursor % fieldhandle % has_unlimited_dim) then
            pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start1, realVal)
         else
            pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % field_desc, realVal)
         end if
      else if (present(intVal)) then
!         write (0,*) '  value is int'
         if (field_cursor % fieldhandle % has_unlimited_dim) then
            pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start1, intVal)
         else
            pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % field_desc, intVal)
         end if
      else if (present(charVal)) then
!         write (0,*) '  value is char'
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
      else if (present(realArray1d)) then
!         write (0,*) '  value is real1'
         if (associated(field_cursor % fieldhandle % decomp)) then
            call PIO_read_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                 realArray1d, pio_ierr)
         else
            start1(:) = 1
            count1(1) = field_cursor % fieldhandle % dims(1) % dimsize
            pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start1, count1, realArray1d)
         end if
      else if (present(realArray2d)) then
!         write (0,*) '  value is real2'
         if (associated(field_cursor % fieldhandle % decomp)) then
            call PIO_read_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                 realArray2d, pio_ierr)
         else
            start2(:) = 1
            count2(1) = field_cursor % fieldhandle % dims(1) % dimsize
            count2(2) = field_cursor % fieldhandle % dims(2) % dimsize
            pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start2, count2, realArray2d)
         end if
      else if (present(realArray3d)) then
!         write (0,*) '  value is real3'
         if (associated(field_cursor % fieldhandle % decomp)) then
            call PIO_read_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                 realArray3d, pio_ierr)
         else
            start3(:) = 1
            count3(1) = field_cursor % fieldhandle % dims(1) % dimsize
            count3(2) = field_cursor % fieldhandle % dims(2) % dimsize
            count3(3) = field_cursor % fieldhandle % dims(3) % dimsize
            pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start3, count3, realArray3d)
         end if
      else if (present(realArray4d)) then
!         write (0,*) '  value is real4'
         if (associated(field_cursor % fieldhandle % decomp)) then
            call PIO_read_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                 realArray4d, pio_ierr)
         else
            start4(:) = 1
            count4(1) = field_cursor % fieldhandle % dims(1) % dimsize
            count4(2) = field_cursor % fieldhandle % dims(2) % dimsize
            count4(3) = field_cursor % fieldhandle % dims(3) % dimsize
            count4(4) = field_cursor % fieldhandle % dims(4) % dimsize
            pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start4, count4, realArray4d)
         end if
      else if (present(realArray5d)) then
!         write (0,*) '  value is real5'
         if (associated(field_cursor % fieldhandle % decomp)) then
            call PIO_read_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                realArray5d, pio_ierr)
         else
            start5(:) = 1
            count5(1) = field_cursor % fieldhandle % dims(1) % dimsize
            count5(2) = field_cursor % fieldhandle % dims(2) % dimsize
            count5(3) = field_cursor % fieldhandle % dims(3) % dimsize
            count5(4) = field_cursor % fieldhandle % dims(4) % dimsize
            count5(5) = field_cursor % fieldhandle % dims(5) % dimsize
            pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start5, count5, realArray5d)
         end if
      else if (present(intArray1d)) then
!         write (0,*) '  value is int1'
         if (associated(field_cursor % fieldhandle % decomp)) then
            call PIO_read_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                 intArray1d, pio_ierr)
         else
            start1(:) = 1
            count1(1) = field_cursor % fieldhandle % dims(1) % dimsize
            pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start1, count1, intArray1d)
         end if
      else if (present(intArray2d)) then
!         write (0,*) '  value is int2'
         if (associated(field_cursor % fieldhandle % decomp)) then
            call PIO_read_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                 intArray2d, pio_ierr)
         else
            start2(:) = 1
            count2(1) = field_cursor % fieldhandle % dims(1) % dimsize
            count2(2) = field_cursor % fieldhandle % dims(2) % dimsize
            pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start2, count2, intArray2d)
         end if
      else if (present(intArray3d)) then
!         write (0,*) '  value is int3'
         if (associated(field_cursor % fieldhandle % decomp)) then
            call PIO_read_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                 intArray3d, pio_ierr)
         else
            start3(:) = 1
            count3(1) = field_cursor % fieldhandle % dims(1) % dimsize
            count3(2) = field_cursor % fieldhandle % dims(2) % dimsize
            count3(3) = field_cursor % fieldhandle % dims(3) % dimsize
            pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start3, count3, intArray3d)
         end if
      else if (present(intArray4d)) then
!         write (0,*) '  value is int4'
         if (associated(field_cursor % fieldhandle % decomp)) then
            call PIO_read_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                 intArray4d, pio_ierr)
         else
            start4(:) = 1
            count4(1) = field_cursor % fieldhandle % dims(1) % dimsize
            count4(2) = field_cursor % fieldhandle % dims(2) % dimsize
            count4(3) = field_cursor % fieldhandle % dims(3) % dimsize
            count4(4) = field_cursor % fieldhandle % dims(4) % dimsize
            pio_ierr = PIO_get_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start4, count4, intArray4d)
         end if
      end if

!      write (0,*) 'Checking for error'
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

!      write(stderrUnit,*) 'Called MPAS_io_get_var_int0d()'
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_get_var_generic(handle, fieldname, intVal=val, ierr=ierr)

   end subroutine MPAS_io_get_var_int0d


   subroutine MPAS_io_get_var_int1d(handle, fieldname, array, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      integer, dimension(:), intent(out) :: array
      integer, intent(out), optional :: ierr

!      write(stderrUnit,*) 'Called MPAS_io_get_var_int1d()'
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_get_var_generic(handle, fieldname, intArray1d=array, ierr=ierr)

   end subroutine MPAS_io_get_var_int1d


   subroutine MPAS_io_get_var_int2d(handle, fieldname, array, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      integer, dimension(:,:), intent(out) :: array
      integer, intent(out), optional :: ierr

!      write(stderrUnit,*) 'Called MPAS_io_get_var_int2d()'
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_get_var_generic(handle, fieldname, intArray2d=array, ierr=ierr)

   end subroutine MPAS_io_get_var_int2d


   subroutine MPAS_io_get_var_int3d(handle, fieldname, array, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      integer, dimension(:,:,:), intent(out) :: array
      integer, intent(out), optional :: ierr

!      write(stderrUnit,*) 'Called MPAS_io_get_var_int3d()'
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_get_var_generic(handle, fieldname, intArray3d=array, ierr=ierr)

   end subroutine MPAS_io_get_var_int3d


   subroutine MPAS_io_get_var_int4d(handle, fieldname, array, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      integer, dimension(:,:,:,:), intent(out) :: array
      integer, intent(out), optional :: ierr

!      write(stderrUnit,*) 'Called MPAS_io_get_var_int4d()'
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_get_var_generic(handle, fieldname, intArray4d=array, ierr=ierr)

   end subroutine MPAS_io_get_var_int4d


   subroutine MPAS_io_get_var_real0d(handle, fieldname, val, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      real (kind=RKIND), intent(out) :: val
      integer, intent(out), optional :: ierr

!      write(stderrUnit,*) 'Called MPAS_io_get_var_real0d()'
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_get_var_generic(handle, fieldname, realVal=val, ierr=ierr)

   end subroutine MPAS_io_get_var_real0d


   subroutine MPAS_io_get_var_real1d(handle, fieldname, array, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      real (kind=RKIND), dimension(:), intent(out) :: array
      integer, intent(out), optional :: ierr

!      write(stderrUnit,*) 'Called MPAS_io_get_var_real1d()'
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_get_var_generic(handle, fieldname, realArray1d=array, ierr=ierr)

   end subroutine MPAS_io_get_var_real1d


   subroutine MPAS_io_get_var_real2d(handle, fieldname, array, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      real (kind=RKIND), dimension(:,:), intent(out) :: array
      integer, intent(out), optional :: ierr

!      write(stderrUnit,*) 'Called MPAS_io_get_var_real2d()'
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_get_var_generic(handle, fieldname, realArray2d=array, ierr=ierr)

   end subroutine MPAS_io_get_var_real2d


   subroutine MPAS_io_get_var_real3d(handle, fieldname, array, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      real (kind=RKIND), dimension(:,:,:), intent(out) :: array
      integer, intent(out), optional :: ierr

!      write(stderrUnit,*) 'Called MPAS_io_get_var_real3d()'
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_get_var_generic(handle, fieldname, realArray3d=array, ierr=ierr)

   end subroutine MPAS_io_get_var_real3d


   subroutine MPAS_io_get_var_real4d(handle, fieldname, array, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      real (kind=RKIND), dimension(:,:,:,:), intent(out) :: array
      integer, intent(out), optional :: ierr

!      write(stderrUnit,*) 'Called MPAS_io_get_var_real4d()'
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_get_var_generic(handle, fieldname, realArray4d=array, ierr=ierr)

   end subroutine MPAS_io_get_var_real4d


   subroutine MPAS_io_get_var_real5d(handle, fieldname, array, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      real (kind=RKIND), dimension(:,:,:,:,:), intent(out) :: array
      integer, intent(out), optional :: ierr

!      write(stderrUnit,*) 'Called MPAS_io_get_var_real5d()'
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_get_var_generic(handle, fieldname, realArray5d=array, ierr=ierr)

   end subroutine MPAS_io_get_var_real5d


   subroutine MPAS_io_get_var_char0d(handle, fieldname, val, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      character (len=*), intent(out) :: val
      integer, intent(out), optional :: ierr

!      write(stderrUnit,*) 'Called MPAS_io_get_var_char0d()'
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_get_var_generic(handle, fieldname, charVal=val, ierr=ierr)

   end subroutine MPAS_io_get_var_char0d


   subroutine MPAS_io_put_var_generic(handle, fieldname, intVal, intArray1d, intArray2d, intArray3d, intArray4d, &
                                                        realVal, realArray1d, realArray2d, realArray3d, realArray4d, realArray5d, &
                                                        charVal, ierr)

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
      type (fieldlist_type), pointer :: field_cursor

      ! Sanity checks
      if (.not. handle % initialized) then
         if (present(ierr)) ierr = MPAS_IO_ERR_UNINIT_HANDLE
         return 
      end if

      if (.not. handle % data_mode) then
         handle % data_mode = .true.

         pio_ierr = PIO_enddef(handle % pio_file)
         if (pio_ierr /= PIO_noerr) then
            if (present(ierr)) ierr = MPAS_IO_ERR_PIO
            return
         end if
      end if

!     write(stderrUnit,*) 'Writing ', trim(fieldname)


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


      !
      ! Check that we have a decomposition for this field
      !
!     if (.not.present(intVal) .and. .not.present(realVal) .and. .not.present(charVal)) then
!        if (.not. associated(field_cursor % fieldhandle % decomp)) then
!           if (present(ierr)) ierr = MPAS_IO_ERR_NO_DECOMP
!           return
!        end if
!     end if

      if (field_cursor % fieldhandle % has_unlimited_dim) then
         call PIO_setframe(field_cursor % fieldhandle % field_desc, handle % frame_number)
         start1(1) = handle % frame_number
         count1(1) = 1
     
         start2(1) = 1
         start2(2) = handle % frame_number
         count2(2) = 1
      else if (handle % frame_number > 1) then
         if(present(ierr)) ierr = MPAS_IO_NOERR
         return
      end if

      if (present(realVal)) then
         if (field_cursor % fieldhandle % has_unlimited_dim) then
            pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, start1, realVal)
         else
            pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % field_desc, realVal)
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
            pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start2, count2, (/charVal/))
         else
            start1(1) = 1
            count1(1) = field_cursor % fieldhandle % dims(1) % dimsize
            pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start1, count1, (/charVal/))
         end if
      else if (present(realArray1d)) then
         if (associated(field_cursor % fieldhandle % decomp)) then
            call PIO_write_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                  realArray1d, pio_ierr)
         else
            if (field_cursor % fieldhandle % has_unlimited_dim) then
               start1(1) = handle % frame_number
               count1(1) = 1
            else
               start1(1) = 1
               count1(1) = field_cursor % fieldhandle % dims(1) % dimsize
            end if
            pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start1, count1, realArray1d)
         end if
      else if (present(realArray2d)) then
         if (associated(field_cursor % fieldhandle % decomp)) then
            call PIO_write_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                  realArray2d, pio_ierr)
         else
            if (field_cursor % fieldhandle % has_unlimited_dim) then
               start2(1) = 1
               start2(2) = handle % frame_number
               count2(1) = field_cursor % fieldhandle % dims(1) % dimsize
               count2(2) = 1
            else
               start2(:) = 1
               count2(1) = field_cursor % fieldhandle % dims(1) % dimsize
               count2(2) = field_cursor % fieldhandle % dims(2) % dimsize
            end if
            pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start2, count2, realArray2d)
         end if
      else if (present(realArray3d)) then
         if (associated(field_cursor % fieldhandle % decomp)) then
            call PIO_write_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                  realArray3d, pio_ierr)
         else
            if (field_cursor % fieldhandle % has_unlimited_dim) then
               start3(1) = 1
               start3(2) = 1
               start3(3) = handle % frame_number
               count3(1) = field_cursor % fieldhandle % dims(1) % dimsize
               count3(2) = field_cursor % fieldhandle % dims(2) % dimsize
               count3(3) = 1
            else
               start3(:) = 1
               count3(1) = field_cursor % fieldhandle % dims(1) % dimsize
               count3(2) = field_cursor % fieldhandle % dims(2) % dimsize
               count3(3) = field_cursor % fieldhandle % dims(3) % dimsize
            end if
            pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start3, count3, realArray3d)
         end if
      else if (present(realArray4d)) then
         if (associated(field_cursor % fieldhandle % decomp)) then
            call PIO_write_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                  realArray4d, pio_ierr)
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
            else
               start4(:) = 1
               count4(1) = field_cursor % fieldhandle % dims(1) % dimsize
               count4(2) = field_cursor % fieldhandle % dims(2) % dimsize
               count4(3) = field_cursor % fieldhandle % dims(3) % dimsize
               count4(4) = field_cursor % fieldhandle % dims(4) % dimsize
            end if
            pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start4, count4, realArray4d)
         end if
      else if (present(realArray5d)) then
         if (associated(field_cursor % fieldhandle % decomp)) then
            call PIO_write_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                  realArray5d, pio_ierr)
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
            else
               start5(:) = 1
               count5(1) = field_cursor % fieldhandle % dims(1) % dimsize
               count5(2) = field_cursor % fieldhandle % dims(2) % dimsize
               count5(3) = field_cursor % fieldhandle % dims(3) % dimsize
               count5(4) = field_cursor % fieldhandle % dims(4) % dimsize
               count5(5) = field_cursor % fieldhandle % dims(5) % dimsize
            end if
            pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start5, count5, realArray5d)
         end if
      else if (present(intArray1d)) then
         if (associated(field_cursor % fieldhandle % decomp)) then
            call PIO_write_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                  intArray1d, pio_ierr)
         else
            if (field_cursor % fieldhandle % has_unlimited_dim) then
               start1(1) = handle % frame_number
               count1(1) = 1
            else
               start1(1) = 1
               count1(1) = field_cursor % fieldhandle % dims(1) % dimsize
            end if
            pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start1, count1, intArray1d)
         end if
      else if (present(intArray2d)) then
         if (associated(field_cursor % fieldhandle % decomp)) then
            call PIO_write_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                  intArray2d, pio_ierr)
         else
            if (field_cursor % fieldhandle % has_unlimited_dim) then
               start2(1) = 1
               start2(2) = handle % frame_number
               count2(1) = field_cursor % fieldhandle % dims(1) % dimsize
               count2(2) = 1
            else
               start2(:) = 1
               count2(1) = field_cursor % fieldhandle % dims(1) % dimsize
               count2(2) = field_cursor % fieldhandle % dims(2) % dimsize
            end if
            pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start2, count2, intArray2d)
         end if
      else if (present(intArray3d)) then
         if (associated(field_cursor % fieldhandle % decomp)) then
            call PIO_write_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                  intArray3d, pio_ierr)
         else
            if (field_cursor % fieldhandle % has_unlimited_dim) then
               start3(1) = 1
               start3(2) = 1
               start3(3) = handle % frame_number
               count3(1) = field_cursor % fieldhandle % dims(1) % dimsize
               count3(2) = field_cursor % fieldhandle % dims(2) % dimsize
               count3(3) = 1
            else
               start3(:) = 1
               count3(1) = field_cursor % fieldhandle % dims(1) % dimsize
               count3(2) = field_cursor % fieldhandle % dims(2) % dimsize
               count3(3) = field_cursor % fieldhandle % dims(3) % dimsize
            end if
            pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start3, count3, intArray3d)
         end if
      else if (present(intArray4d)) then
         if (associated(field_cursor % fieldhandle % decomp)) then
            call PIO_write_darray(handle % pio_file, field_cursor % fieldhandle % field_desc, field_cursor % fieldhandle % decomp % pio_iodesc, &
                                  intArray4d, pio_ierr)
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
            else
               start4(:) = 1
               count4(1) = field_cursor % fieldhandle % dims(1) % dimsize
               count4(2) = field_cursor % fieldhandle % dims(2) % dimsize
               count4(3) = field_cursor % fieldhandle % dims(3) % dimsize
               count4(4) = field_cursor % fieldhandle % dims(4) % dimsize
            end if
            pio_ierr = PIO_put_var(handle % pio_file, field_cursor % fieldhandle % fieldid, start4, count4, intArray4d)
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

!      write(stderrUnit,*) 'Called MPAS_io_put_var_int0d()'
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_put_var_generic(handle, fieldname, intVal=val, ierr=ierr)

   end subroutine MPAS_io_put_var_int0d


   subroutine MPAS_io_put_var_int1d(handle, fieldname, array, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      integer, dimension(:), intent(in) :: array
      integer, intent(out), optional :: ierr

!      write(stderrUnit,*) 'Called MPAS_io_put_var_int1d()'
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_put_var_generic(handle, fieldname, intArray1d=array, ierr=ierr)

   end subroutine MPAS_io_put_var_int1d


   subroutine MPAS_io_put_var_int2d(handle, fieldname, array, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      integer, dimension(:,:), intent(in) :: array
      integer, intent(out), optional :: ierr

!      write(stderrUnit,*) 'Called MPAS_io_put_var_int2d()'
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_put_var_generic(handle, fieldname, intArray2d=array, ierr=ierr)

   end subroutine MPAS_io_put_var_int2d


   subroutine MPAS_io_put_var_int3d(handle, fieldname, array, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      integer, dimension(:,:,:), intent(in) :: array
      integer, intent(out), optional :: ierr

!      write(stderrUnit,*) 'Called MPAS_io_put_var_int3d()'
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_put_var_generic(handle, fieldname, intArray3d=array, ierr=ierr)

   end subroutine MPAS_io_put_var_int3d


   subroutine MPAS_io_put_var_int4d(handle, fieldname, array, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      integer, dimension(:,:,:,:), intent(in) :: array
      integer, intent(out), optional :: ierr

!      write(stderrUnit,*) 'Called MPAS_io_put_var_int4d()'
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_put_var_generic(handle, fieldname, intArray4d=array, ierr=ierr)

   end subroutine MPAS_io_put_var_int4d


   subroutine MPAS_io_put_var_real0d(handle, fieldname, val, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      real (kind=RKIND), intent(in) :: val
      integer, intent(out), optional :: ierr

!      write(stderrUnit,*) 'Called MPAS_io_put_var_real0d()'
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_put_var_generic(handle, fieldname, realVal=val, ierr=ierr)

   end subroutine MPAS_io_put_var_real0d


   subroutine MPAS_io_put_var_real1d(handle, fieldname, array, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      real (kind=RKIND), dimension(:), intent(in) :: array
      integer, intent(out), optional :: ierr

!      write(stderrUnit,*) 'Called MPAS_io_put_var_real1d()'
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_put_var_generic(handle, fieldname, realArray1d=array, ierr=ierr)

   end subroutine MPAS_io_put_var_real1d


   subroutine MPAS_io_put_var_real2d(handle, fieldname, array, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      real (kind=RKIND), dimension(:,:), intent(in) :: array
      integer, intent(out), optional :: ierr

!      write(stderrUnit,*) 'Called MPAS_io_put_var_real2d()'
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_put_var_generic(handle, fieldname, realArray2d=array, ierr=ierr)

   end subroutine MPAS_io_put_var_real2d


   subroutine MPAS_io_put_var_real3d(handle, fieldname, array, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      real (kind=RKIND), dimension(:,:,:), intent(in) :: array
      integer, intent(out), optional :: ierr

!      write(stderrUnit,*) 'Called MPAS_io_put_var_real3d()'
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_put_var_generic(handle, fieldname, realArray3d=array, ierr=ierr)

   end subroutine MPAS_io_put_var_real3d


   subroutine MPAS_io_put_var_real4d(handle, fieldname, array, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      real (kind=RKIND), dimension(:,:,:,:), intent(in) :: array
      integer, intent(out), optional :: ierr

!      write(stderrUnit,*) 'Called MPAS_io_put_var_real4d()'
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_put_var_generic(handle, fieldname, realArray4d=array, ierr=ierr)

   end subroutine MPAS_io_put_var_real4d


   subroutine MPAS_io_put_var_real5d(handle, fieldname, array, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      real (kind=RKIND), dimension(:,:,:,:,:), intent(in) :: array
      integer, intent(out), optional :: ierr

!      write(stderrUnit,*) 'Called MPAS_io_put_var_real5d()'
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_put_var_generic(handle, fieldname, realArray5d=array, ierr=ierr)

   end subroutine MPAS_io_put_var_real5d


   subroutine MPAS_io_put_var_char0d(handle, fieldname, val, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: fieldname
      character (len=*), intent(in) :: val
      integer, intent(out), optional :: ierr

!      write(stderrUnit,*) 'Called MPAS_io_put_var_char0d()'
      if (present(ierr)) ierr = MPAS_IO_NOERR

      call MPAS_io_put_var_generic(handle, fieldname, charVal=val, ierr=ierr)

   end subroutine MPAS_io_put_var_char0d


   subroutine MPAS_io_get_att_int0d(handle, attName, attValue, fieldname, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: attName
      integer, intent(out) :: attValue
      character (len=*), intent(in), optional :: fieldname
      integer, intent(out), optional :: ierr

      integer :: pio_ierr
      integer :: varid
      integer :: xtype, len
      type (fieldlist_type), pointer :: field_cursor
      type (attlist_type), pointer :: att_cursor, new_att_node

!      write(stderrUnit,*) 'Called MPAS_io_get_att_int0d()'
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
               if (att_cursor % atthandle % attType == ATT_INT) then
!write(stderrUnit,*) 'Using cached attribute'
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
               if (att_cursor % atthandle % attType == ATT_INT) then
!write(stderrUnit,*) 'Using cached attribute'
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
      new_att_node % atthandle % attType = ATT_INT
      new_att_node % atthandle % attValueInt = attValue

      if (present(fieldname)) then
         if (.not. associated(field_cursor % fieldhandle % attlist_head)) then
            field_cursor % fieldhandle % attlist_head => new_att_node
!write(stderrUnit,*) 'Assigning att head for '//trim(attName)
         end if
         if (.not. associated(field_cursor % fieldhandle % attlist_tail)) then
            field_cursor % fieldhandle % attlist_tail => new_att_node
!write(stderrUnit,*) 'Assigning att tail for '//trim(attName)
         else
            field_cursor % fieldhandle % attlist_tail % next => new_att_node
            field_cursor % fieldhandle % attlist_tail => field_cursor % fieldhandle % attlist_tail % next
!write(stderrUnit,*) 'Extending att tail for '//trim(attName)
         end if
      else
         if (.not. associated(handle % attlist_head)) then
            handle % attlist_head => new_att_node
!write(stderrUnit,*) 'Assigning att head for '//trim(attName)
         end if
         if (.not. associated(handle % attlist_tail)) then
            handle % attlist_tail => new_att_node
!write(stderrUnit,*) 'Assigning att tail for '//trim(attName)
         else
            handle % attlist_tail % next => new_att_node
            handle % attlist_tail => handle % attlist_tail % next
!write(stderrUnit,*) 'Extending att tail for '//trim(attName)
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
      integer :: xtype, len, attlen
      type (fieldlist_type), pointer :: field_cursor
      type (attlist_type), pointer :: att_cursor, new_att_node

!      write(stderrUnit,*) 'Called MPAS_io_get_att_int1d()'
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
               if (att_cursor % atthandle % attType == ATT_INTA) then
!write(stderrUnit,*) 'Using cached attribute'
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
               if (att_cursor % atthandle % attType == ATT_INTA) then
!write(stderrUnit,*) 'Using cached attribute'
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
      new_att_node % atthandle % attType = ATT_INTA
      allocate(new_att_node % atthandle % attValueIntA(attlen))
      new_att_node % atthandle % attValueIntA = attValue

      if (present(fieldname)) then
         if (.not. associated(field_cursor % fieldhandle % attlist_head)) then
            field_cursor % fieldhandle % attlist_head => new_att_node
!write(stderrUnit,*) 'Assigning att head for '//trim(attName)
         end if
         if (.not. associated(field_cursor % fieldhandle % attlist_tail)) then
            field_cursor % fieldhandle % attlist_tail => new_att_node
!write(stderrUnit,*) 'Assigning att tail for '//trim(attName)
         else
            field_cursor % fieldhandle % attlist_tail % next => new_att_node
            field_cursor % fieldhandle % attlist_tail => field_cursor % fieldhandle % attlist_tail % next
!write(stderrUnit,*) 'Extending att tail for '//trim(attName)
         end if
      else
         if (.not. associated(handle % attlist_head)) then
            handle % attlist_head => new_att_node
!write(stderrUnit,*) 'Assigning att head for '//trim(attName)
         end if
         if (.not. associated(handle % attlist_tail)) then
            handle % attlist_tail => new_att_node
!write(stderrUnit,*) 'Assigning att tail for '//trim(attName)
         else
            handle % attlist_tail % next => new_att_node
            handle % attlist_tail => handle % attlist_tail % next
!write(stderrUnit,*) 'Extending att tail for '//trim(attName)
         end if
      end if

   end subroutine MPAS_io_get_att_int1d


   subroutine MPAS_io_get_att_real0d(handle, attName, attValue, fieldname, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: attName
      real (kind=RKIND), intent(out) :: attValue
      character (len=*), intent(in), optional :: fieldname
      integer, intent(out), optional :: ierr

      integer :: pio_ierr
      integer :: varid
      integer :: xtype, len
      type (fieldlist_type), pointer :: field_cursor
      type (attlist_type), pointer :: att_cursor, new_att_node

!      write(stderrUnit,*) 'Called MPAS_io_get_att_real0d()'
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
               if (att_cursor % atthandle % attType == ATT_REAL) then
!write(stderrUnit,*) 'Using cached attribute'
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
               if (att_cursor % atthandle % attType == ATT_REAL) then
!write(stderrUnit,*) 'Using cached attribute'
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

      ! Query attribute value
      pio_ierr = PIO_inq_att(handle % pio_file, varid, attName, xtype, len)
      if (pio_ierr /= PIO_noerr) then
         if (present(ierr)) ierr = MPAS_IO_ERR_PIO
         return
      end if
      if (xtype /= PIO_REALKIND) then
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
      new_att_node % atthandle % attType = ATT_REAL
      new_att_node % atthandle % attValueReal = attValue

      if (present(fieldname)) then
         if (.not. associated(field_cursor % fieldhandle % attlist_head)) then
            field_cursor % fieldhandle % attlist_head => new_att_node
!write(stderrUnit,*) 'Assigning att head for '//trim(attName)
         end if
         if (.not. associated(field_cursor % fieldhandle % attlist_tail)) then
            field_cursor % fieldhandle % attlist_tail => new_att_node
!write(stderrUnit,*) 'Assigning att tail for '//trim(attName)
         else
            field_cursor % fieldhandle % attlist_tail % next => new_att_node
            field_cursor % fieldhandle % attlist_tail => field_cursor % fieldhandle % attlist_tail % next
!write(stderrUnit,*) 'Extending att tail for '//trim(attName)
         end if
      else
         if (.not. associated(handle % attlist_head)) then
            handle % attlist_head => new_att_node
!write(stderrUnit,*) 'Assigning att head for '//trim(attName)
         end if
         if (.not. associated(handle % attlist_tail)) then
            handle % attlist_tail => new_att_node
!write(stderrUnit,*) 'Assigning att tail for '//trim(attName)
         else
            handle % attlist_tail % next => new_att_node
            handle % attlist_tail => handle % attlist_tail % next
!write(stderrUnit,*) 'Extending att tail for '//trim(attName)
         end if
      end if

   end subroutine MPAS_io_get_att_real0d


   subroutine MPAS_io_get_att_real1d(handle, attName, attValue, fieldname, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: attName
      real (kind=RKIND), dimension(:), pointer :: attValue
      character (len=*), intent(in), optional :: fieldname
      integer, intent(out), optional :: ierr

      integer :: pio_ierr
      integer :: varid
      integer :: xtype, len, attlen
      type (fieldlist_type), pointer :: field_cursor
      type (attlist_type), pointer :: att_cursor, new_att_node

!      write(stderrUnit,*) 'Called MPAS_io_get_att_real1d()'
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
               if (att_cursor % atthandle % attType == ATT_REALA) then
!write(stderrUnit,*) 'Using cached attribute'
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
               if (att_cursor % atthandle % attType == ATT_REALA) then
!write(stderrUnit,*) 'Using cached attribute'
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

      ! Query attribute value
      pio_ierr = PIO_inq_att(handle % pio_file, varid, attName, xtype, len)
      if (pio_ierr /= PIO_noerr) then
         if (present(ierr)) ierr = MPAS_IO_ERR_PIO
         return
      end if

      if (xtype /= PIO_REALKIND) then
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
      new_att_node % atthandle % attType = ATT_REALA
      allocate(new_att_node % atthandle % attValueRealA(attlen))
      new_att_node % atthandle % attValueRealA = attValue

      if (present(fieldname)) then
         if (.not. associated(field_cursor % fieldhandle % attlist_head)) then
            field_cursor % fieldhandle % attlist_head => new_att_node
!write(stderrUnit,*) 'Assigning att head for '//trim(attName)
         end if
         if (.not. associated(field_cursor % fieldhandle % attlist_tail)) then
            field_cursor % fieldhandle % attlist_tail => new_att_node
!write(stderrUnit,*) 'Assigning att tail for '//trim(attName)
         else
            field_cursor % fieldhandle % attlist_tail % next => new_att_node
            field_cursor % fieldhandle % attlist_tail => field_cursor % fieldhandle % attlist_tail % next
!write(stderrUnit,*) 'Extending att tail for '//trim(attName)
         end if
      else
         if (.not. associated(handle % attlist_head)) then
            handle % attlist_head => new_att_node
!write(stderrUnit,*) 'Assigning att head for '//trim(attName)
         end if
         if (.not. associated(handle % attlist_tail)) then
            handle % attlist_tail => new_att_node
!write(stderrUnit,*) 'Assigning att tail for '//trim(attName)
         else
            handle % attlist_tail % next => new_att_node
            handle % attlist_tail => handle % attlist_tail % next
!write(stderrUnit,*) 'Extending att tail for '//trim(attName)
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
      integer :: xtype, len
      type (fieldlist_type), pointer :: field_cursor
      type (attlist_type), pointer :: att_cursor, new_att_node

!      write(stderrUnit,*) 'Called MPAS_io_get_att_text()'
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
               if (att_cursor % atthandle % attType == ATT_TEXT) then
!write(stderrUnit,*) 'Using cached attribute'
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
               if (att_cursor % atthandle % attType == ATT_TEXT) then
!write(stderrUnit,*) 'Using cached attribute'
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
      new_att_node % atthandle % attType = ATT_TEXT
      new_att_node % atthandle % attValueText = attValue

      if (present(fieldname)) then
         if (.not. associated(field_cursor % fieldhandle % attlist_head)) then
            field_cursor % fieldhandle % attlist_head => new_att_node
!write(stderrUnit,*) 'Assigning att head for '//trim(attName)
         end if
         if (.not. associated(field_cursor % fieldhandle % attlist_tail)) then
            field_cursor % fieldhandle % attlist_tail => new_att_node
!write(stderrUnit,*) 'Assigning att tail for '//trim(attName)
         else
            field_cursor % fieldhandle % attlist_tail % next => new_att_node
            field_cursor % fieldhandle % attlist_tail => field_cursor % fieldhandle % attlist_tail % next
!write(stderrUnit,*) 'Extending att tail for '//trim(attName)
         end if
      else
         if (.not. associated(handle % attlist_head)) then
            handle % attlist_head => new_att_node
!write(stderrUnit,*) 'Assigning att head for '//trim(attName)
         end if
         if (.not. associated(handle % attlist_tail)) then
            handle % attlist_tail => new_att_node
!write(stderrUnit,*) 'Assigning att tail for '//trim(attName)
         else
            handle % attlist_tail % next => new_att_node
            handle % attlist_tail => handle % attlist_tail % next
!write(stderrUnit,*) 'Extending att tail for '//trim(attName)
         end if
      end if

   end subroutine MPAS_io_get_att_text


   subroutine MPAS_io_put_att_int0d(handle, attName, attValue, fieldname, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: attName
      integer, intent(in) :: attValue
      character (len=*), intent(in), optional :: fieldname
      integer, intent(out), optional :: ierr

      integer :: pio_ierr
      integer :: varid
      type (fieldlist_type), pointer :: field_cursor
      type (attlist_type), pointer :: attlist_cursor, new_attlist_node

!      write(stderrUnit,*) 'Called MPAS_io_put_att_int0d()'
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
      new_attlist_node % attHandle % attType = ATT_INT
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
!write(stderrUnit,*) 'Attribute already defined'
                     if (attlist_cursor % atthandle % attType /= ATT_INT .or. &
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
!write(stderrUnit,*) 'Assigning att head for '//trim(attname)
         end if
         if (.not. associated(field_cursor % fieldhandle % attlist_tail)) then
            field_cursor % fieldhandle % attlist_tail => new_attlist_node
!write(stderrUnit,*) 'Assigning att tail for '//trim(attname)
         else
            field_cursor % fieldhandle % attlist_tail % next => new_attlist_node
            field_cursor % fieldhandle % attlist_tail => field_cursor % fieldhandle % attlist_tail % next
!write(stderrUnit,*) 'Extending att tail for '//trim(attname)
         end if

      else

         ! Check whether attribute was already defined
         attlist_cursor => handle % attlist_head
         do while (associated(attlist_cursor))
            if (trim(attName) == trim(attlist_cursor % atthandle % attName)) then
!write(stderrUnit,*) 'Attribute already defined'
               if (attlist_cursor % atthandle % attType /= ATT_INT .or. &
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
!write(stderrUnit,*) 'Assigning att head for '//trim(attname)
         end if
         if (.not. associated(handle % attlist_tail)) then
            handle % attlist_tail => new_attlist_node
!write(stderrUnit,*) 'Assigning att tail for '//trim(attname)
         else
            handle % attlist_tail % next => new_attlist_node
            handle % attlist_tail => handle % attlist_tail % next
!write(stderrUnit,*) 'Extending att tail for '//trim(attname)
         end if
      end if

      pio_ierr = PIO_put_att(handle % pio_file, varid, attName, attValue) 
      if (pio_ierr /= PIO_noerr) then
         if (present(ierr)) ierr = MPAS_IO_ERR_PIO
         return
      end if

      ! Maybe we should add attribute to list only after a successfull call to PIO?

   end subroutine MPAS_io_put_att_int0d


   subroutine MPAS_io_put_att_int1d(handle, attName, attValue, fieldname, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: attName
      integer, dimension(:), intent(in) :: attValue
      character (len=*), intent(in), optional :: fieldname
      integer, intent(out), optional :: ierr

      integer :: pio_ierr
      integer :: varid
      type (fieldlist_type), pointer :: field_cursor
      type (attlist_type), pointer :: attlist_cursor, new_attlist_node

!      write(stderrUnit,*) 'Called MPAS_io_put_att_int1d()'
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
      new_attlist_node % attHandle % attType = ATT_INTA
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
!write(stderrUnit,*) 'Attribute already defined'
                     if (attlist_cursor % atthandle % attType /= ATT_INTA .or. &
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
!write(stderrUnit,*) 'Assigning att head for '//trim(attname)
         end if
         if (.not. associated(field_cursor % fieldhandle % attlist_tail)) then
            field_cursor % fieldhandle % attlist_tail => new_attlist_node
!write(stderrUnit,*) 'Assigning att tail for '//trim(attname)
         else
            field_cursor % fieldhandle % attlist_tail % next => new_attlist_node
            field_cursor % fieldhandle % attlist_tail => field_cursor % fieldhandle % attlist_tail % next
!write(stderrUnit,*) 'Extending att tail for '//trim(attname)
         end if

      else

         ! Check whether attribute was already defined
         attlist_cursor => handle % attlist_head
         do while (associated(attlist_cursor))
            if (trim(attName) == trim(attlist_cursor % atthandle % attName)) then
!write(stderrUnit,*) 'Attribute already defined'
               if (attlist_cursor % atthandle % attType /= ATT_INTA .or. &
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
!write(stderrUnit,*) 'Assigning att head for '//trim(attname)
         end if
         if (.not. associated(handle % attlist_tail)) then
            handle % attlist_tail => new_attlist_node
!write(stderrUnit,*) 'Assigning att tail for '//trim(attname)
         else
            handle % attlist_tail % next => new_attlist_node
            handle % attlist_tail => handle % attlist_tail % next
!write(stderrUnit,*) 'Extending att tail for '//trim(attname)
         end if
      end if

      pio_ierr = PIO_put_att(handle % pio_file, varid, attName, attValue) 
      if (pio_ierr /= PIO_noerr) then
         if (present(ierr)) ierr = MPAS_IO_ERR_PIO
         return
      end if

      ! Maybe we should add attribute to list only after a successfull call to PIO?

   end subroutine MPAS_io_put_att_int1d


   subroutine MPAS_io_put_att_real0d(handle, attName, attValue, fieldname, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: attName
      real (kind=RKIND), intent(in) :: attValue
      character (len=*), intent(in), optional :: fieldname
      integer, intent(out), optional :: ierr

      integer :: pio_ierr
      integer :: varid
      type (fieldlist_type), pointer :: field_cursor
      type (attlist_type), pointer :: attlist_cursor, new_attlist_node

!      write(stderrUnit,*) 'Called MPAS_io_put_att_real0d()'
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
      new_attlist_node % attHandle % attType = ATT_REAL
      new_attlist_node % attHandle % attValueReal = attValue


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
!write(stderrUnit,*) 'Attribute already defined'
                     if (attlist_cursor % atthandle % attType /= ATT_REAL .or. &
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
!write(stderrUnit,*) 'Assigning att head for '//trim(attname)
         end if
         if (.not. associated(field_cursor % fieldhandle % attlist_tail)) then
            field_cursor % fieldhandle % attlist_tail => new_attlist_node
!write(stderrUnit,*) 'Assigning att tail for '//trim(attname)
         else
            field_cursor % fieldhandle % attlist_tail % next => new_attlist_node
            field_cursor % fieldhandle % attlist_tail => field_cursor % fieldhandle % attlist_tail % next
!write(stderrUnit,*) 'Extending att tail for '//trim(attname)
         end if

      else

         ! Check whether attribute was already defined
         attlist_cursor => handle % attlist_head
         do while (associated(attlist_cursor))
            if (trim(attName) == trim(attlist_cursor % atthandle % attName)) then
!write(stderrUnit,*) 'Attribute already defined'
               if (attlist_cursor % atthandle % attType /= ATT_REAL .or. &
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
!write(stderrUnit,*) 'Assigning att head for '//trim(attname)
         end if
         if (.not. associated(handle % attlist_tail)) then
            handle % attlist_tail => new_attlist_node
!write(stderrUnit,*) 'Assigning att tail for '//trim(attname)
         else
            handle % attlist_tail % next => new_attlist_node
            handle % attlist_tail => handle % attlist_tail % next
!write(stderrUnit,*) 'Extending att tail for '//trim(attname)
         end if
      end if

      pio_ierr = PIO_put_att(handle % pio_file, varid, attName, attValue) 
      if (pio_ierr /= PIO_noerr) then
         if (present(ierr)) ierr = MPAS_IO_ERR_PIO
         return
      end if

      ! Maybe we should add attribute to list only after a successfull call to PIO?

   end subroutine MPAS_io_put_att_real0d


   subroutine MPAS_io_put_att_real1d(handle, attName, attValue, fieldname, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: attName
      real (kind=RKIND), dimension(:), intent(in) :: attValue
      character (len=*), intent(in), optional :: fieldname
      integer, intent(out), optional :: ierr

      integer :: pio_ierr
      integer :: varid
      type (fieldlist_type), pointer :: field_cursor
      type (attlist_type), pointer :: attlist_cursor, new_attlist_node

!      write(stderrUnit,*) 'Called MPAS_io_put_att_real1d()'
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
      new_attlist_node % attHandle % attType = ATT_REALA
      allocate(new_attlist_node % attHandle % attValueRealA(size(attValue)))
      new_attlist_node % attHandle % attValueRealA = attValue


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
!write(stderrUnit,*) 'Attribute already defined'
                     if (attlist_cursor % atthandle % attType /= ATT_REALA .or. &
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
!write(stderrUnit,*) 'Assigning att head for '//trim(attname)
         end if
         if (.not. associated(field_cursor % fieldhandle % attlist_tail)) then
            field_cursor % fieldhandle % attlist_tail => new_attlist_node
!write(stderrUnit,*) 'Assigning att tail for '//trim(attname)
         else
            field_cursor % fieldhandle % attlist_tail % next => new_attlist_node
            field_cursor % fieldhandle % attlist_tail => field_cursor % fieldhandle % attlist_tail % next
!write(stderrUnit,*) 'Extending att tail for '//trim(attname)
         end if

      else

         ! Check whether attribute was already defined
         attlist_cursor => handle % attlist_head
         do while (associated(attlist_cursor))
            if (trim(attName) == trim(attlist_cursor % atthandle % attName)) then
!write(stderrUnit,*) 'Attribute already defined'
               if (attlist_cursor % atthandle % attType /= ATT_REALA .or. &
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
!write(stderrUnit,*) 'Assigning att head for '//trim(attname)
         end if
         if (.not. associated(handle % attlist_tail)) then
            handle % attlist_tail => new_attlist_node
!write(stderrUnit,*) 'Assigning att tail for '//trim(attname)
         else
            handle % attlist_tail % next => new_attlist_node
            handle % attlist_tail => handle % attlist_tail % next
!write(stderrUnit,*) 'Extending att tail for '//trim(attname)
         end if
      end if

      pio_ierr = PIO_put_att(handle % pio_file, varid, attName, attValue) 
      if (pio_ierr /= PIO_noerr) then
         if (present(ierr)) ierr = MPAS_IO_ERR_PIO
         return
      end if

      ! Maybe we should add attribute to list only after a successfull call to PIO?

   end subroutine MPAS_io_put_att_real1d


   subroutine MPAS_io_put_att_text(handle, attName, attValue, fieldname, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      character (len=*), intent(in) :: attName
      character (len=*), intent(in) :: attValue
      character (len=*), intent(in), optional :: fieldname
      integer, intent(out), optional :: ierr

      integer :: pio_ierr
      integer :: varid
      type (fieldlist_type), pointer :: field_cursor
      type (attlist_type), pointer :: attlist_cursor, new_attlist_node

!      write(stderrUnit,*) 'Called MPAS_io_put_att_text()'
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
      new_attlist_node % attHandle % attType = ATT_TEXT
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
!write(stderrUnit,*) 'Attribute already defined'
                     if (attlist_cursor % atthandle % attType /= ATT_TEXT .or. &
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
!write(stderrUnit,*) 'Assigning att head for '//trim(attname)
         end if
         if (.not. associated(field_cursor % fieldhandle % attlist_tail)) then
            field_cursor % fieldhandle % attlist_tail => new_attlist_node
!write(stderrUnit,*) 'Assigning att tail for '//trim(attname)
         else
            field_cursor % fieldhandle % attlist_tail % next => new_attlist_node
            field_cursor % fieldhandle % attlist_tail => field_cursor % fieldhandle % attlist_tail % next
!write(stderrUnit,*) 'Extending att tail for '//trim(attname)
         end if

      else

         ! Check whether attribute was already defined
         attlist_cursor => handle % attlist_head
         do while (associated(attlist_cursor))
            if (trim(attName) == trim(attlist_cursor % atthandle % attName)) then
!write(stderrUnit,*) 'Attribute already defined'
               if (attlist_cursor % atthandle % attType /= ATT_TEXT .or. &
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
!write(stderrUnit,*) 'Assigning att head for '//trim(attname)
         end if
         if (.not. associated(handle % attlist_tail)) then
            handle % attlist_tail => new_attlist_node
!write(stderrUnit,*) 'Assigning att tail for '//trim(attname)
         else
            handle % attlist_tail % next => new_attlist_node
            handle % attlist_tail => handle % attlist_tail % next
!write(stderrUnit,*) 'Extending att tail for '//trim(attname)
         end if
      end if

      pio_ierr = PIO_put_att(handle % pio_file, varid, attName, trim(attValue)) 
      if (pio_ierr /= PIO_noerr) then
         if (present(ierr)) ierr = MPAS_IO_ERR_PIO
         return
      end if

      ! Maybe we should add attribute to list only after a successfull call to PIO?

   end subroutine MPAS_io_put_att_text


   subroutine MPAS_io_set_frame(handle, frame, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      integer, intent(in) :: frame
      integer, intent(out), optional :: ierr

!      write(stderrUnit,*) 'Called MPAS_io_set_frame()'
      if (present(ierr)) ierr = MPAS_IO_NOERR

      handle % frame_number = frame

   end subroutine MPAS_io_set_frame


   subroutine MPAS_io_advance_frame(handle, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      integer, intent(out), optional :: ierr

!      write(stderrUnit,*) 'Called MPAS_io_advance_frame()'
      if (present(ierr)) ierr = MPAS_IO_NOERR

      handle % frame_number = handle % frame_number + 1

   end subroutine MPAS_io_advance_frame


   subroutine MPAS_io_sync(handle, ierr)

      implicit none

      type (MPAS_IO_Handle_type), intent(inout) :: handle
      integer, intent(out), optional :: ierr

!      write(stderrUnit,*) 'Called MPAS_io_sync()'
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

!      write(stderrUnit,*) 'Called MPAS_io_close()'
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
            if (attlist_del % atthandle % attType == ATT_INTA) deallocate(attlist_del % atthandle % attValueIntA)
            if (attlist_del % atthandle % attType == ATT_REALA) deallocate(attlist_del % atthandle % attValueRealA)
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
         if (attlist_del % atthandle % attType == ATT_INTA) deallocate(attlist_del % atthandle % attValueIntA)
         if (attlist_del % atthandle % attType == ATT_REALA) deallocate(attlist_del % atthandle % attValueRealA)
         deallocate(attlist_del % atthandle)
      end do
      nullify(handle % attlist_head)
      nullify(handle % attlist_tail)

      handle % initialized = .false.

!write(stderrUnit,*) 'MGD PIO_closefile'
      call PIO_closefile(handle % pio_file)

   end subroutine MPAS_io_close


   subroutine MPAS_io_finalize(io_system, ierr)

      implicit none

      type (iosystem_desc_t), optional, pointer :: io_system
      integer, intent(out), optional :: ierr

      integer :: pio_ierr
      type (decomplist_type), pointer :: decomp_cursor, decomp_del

!      write(stderrUnit,*) 'Called MPAS_io_finalize()'
      if (present(ierr)) ierr = MPAS_IO_NOERR

      decomp_cursor => decomp_list
      do while (associated(decomp_cursor))
         decomp_del => decomp_cursor
         decomp_cursor => decomp_cursor % next
!write(stderrUnit,*) 'Deallocating a decomposition...'
!if (.not. associated(decomp_del % decomphandle)) write(stderrUnit,*) 'OOPS... do not have decomphandle'
         deallocate(decomp_del % decomphandle % dims)
         deallocate(decomp_del % decomphandle % indices)
         call PIO_freedecomp(pio_iosystem, decomp_del % decomphandle % pio_iodesc)
         deallocate(decomp_del % decomphandle)
         deallocate(decomp_del)
      end do

!write(stderrUnit,*) 'MGD PIO_finalize'
      if(.not.present(io_system)) then
        call PIO_finalize(pio_iosystem, pio_ierr)
        if (pio_ierr /= PIO_noerr) then
           if (present(ierr)) ierr = MPAS_IO_ERR_PIO
           return
        end if
        deallocate(pio_iosystem)
      end if

   end subroutine MPAS_io_finalize


   subroutine MPAS_io_err_mesg(ierr, fatal)

      implicit none

      integer, intent(in) :: ierr
      logical, intent(in) :: fatal

      select case (ierr)
         case (MPAS_IO_NOERR)
            ! ... do nothing ...
         case (MPAS_IO_ERR_INVALID_MODE)
            write(stderrUnit,*) 'MPAS IO Error: Invalid file access mode'
         case (MPAS_IO_ERR_INVALID_FORMAT)
            write(stderrUnit,*) 'MPAS IO Error: Invalid I/O format'
         case (MPAS_IO_ERR_LONG_FILENAME)
            write(stderrUnit,*) 'MPAS IO Error: Filename too long'
         case (MPAS_IO_ERR_UNINIT_HANDLE)
            write(stderrUnit,*) 'MPAS IO Error: Uninitialized I/O handle'
         case (MPAS_IO_ERR_PIO)
            write(stderrUnit,*) 'MPAS IO Error: Bad return value from PIO'
         case (MPAS_IO_ERR_DATA_MODE)
            write(stderrUnit,*) 'MPAS IO Error: Cannot define in data mode'
         case (MPAS_IO_ERR_NOWRITE)
            write(stderrUnit,*) 'MPAS IO Error: File not opened for writing'
         case (MPAS_IO_ERR_REDEF_DIM)
            write(stderrUnit,*) 'MPAS IO Error: Inconsistent redefinition of dimension'
         case (MPAS_IO_ERR_REDEF_VAR)
            write(stderrUnit,*) 'MPAS IO Error: Inconsistent redefinition of field'
         case (MPAS_IO_ERR_UNDEFINED_DIM)
            write(stderrUnit,*) 'MPAS IO Error: Field uses undefined dimension'
         case (MPAS_IO_ERR_UNDEFINED_VAR)
            write(stderrUnit,*) 'MPAS IO Error: Undefined field'
         case (MPAS_IO_ERR_REDEF_ATT)
            write(stderrUnit,*) 'MPAS IO Error: Inconsistent redefinition of attribute'
         case (MPAS_IO_ERR_WRONG_ATT_TYPE)
            write(stderrUnit,*) 'MPAS IO Error: Wrong type for requested attribute'
         case (MPAS_IO_ERR_NO_DECOMP)
            write(stderrUnit,*) 'MPAS IO Error: Decomposition indices not set for field'
         case (MPAS_IO_ERR_TWO_UNLIMITED_DIMS)
            write(stderrUnit,*) 'MPAS IO Error: Defining more than one unlimited dimension'
         case (MPAS_IO_ERR_WRONG_MODE)
            write(stderrUnit,*) 'MPAS IO Error: Operation not permitted in this file mode'
         case (MPAS_IO_ERR_NO_UNLIMITED_DIM)
            write(stderrUnit,*) 'MPAS IO Error: No unlimited dimension found in dataset'
         case (MPAS_IO_ERR_UNIMPLEMENTED)
            write(stderrUnit,*) 'MPAS IO Error: Unimplemented functionality'
         case default
            write(stderrUnit,*) 'MPAS IO Error: Unrecognized error code...'
      end select

      if (fatal .and. (ierr /= MPAS_IO_NOERR)) call mpas_dmpar_abort(local_dminfo)

   end subroutine MPAS_io_err_mesg
 
end module mpas_io
