!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module io_netcdf

!BOP
! !MODULE: io_netcdf
! !DESCRIPTION:
!  This module provides a generic input/output interface
!  for writing arrays in netCDF format using pio.
!
! !REVISION HISTORY:
!  SVN:$Id: io_netcdf.F90 45958 2013-04-12 23:10:42Z mlevy@ucar.edu $

! !USES:

   use POP_KindsMod
   use POP_IOUnitsMod
   use POP_ErrorMod
   use kinds_mod
   use domain_size
   use domain
   use constants
   use communicate
   use broadcast
   use gather_scatter
   use exit_mod
   use io_types
   use io_tools
   use io_pio
   use pio
   use shr_sys_mod

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: open_read_netcdf,    &
             open_netcdf,         &
             close_netcdf,        &
             sync_netcdf,         &
             define_field_netcdf, &
             read_field_netcdf,   &
             write_field_netcdf,  &
             define_nstd_netcdf,  &
             write_nstd_netcdf,   &
             write_time_bounds

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  module variables
!
!-----------------------------------------------------------------------

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: open_read_netcdf
! !INTERFACE:

 subroutine open_read_netcdf(data_file)

! !INPUT/OUTPUT PARAMETERS:

   type (datafile), target, intent (inout)  :: data_file

! !DESCRIPTION:
!  This routine opens a netcdf data file and extracts global file
!  attributes.
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


   character (char_len) :: &
      path            ! filename to read

   character (80) :: &
      work_line,     &! temporary to use for parsing file lines
      att_name        ! temporary to use for attribute names


   integer (i4) ::  &
      iostat,       &! status flag
      nsize,        &! size parameter returned by inquire function
      n,            &! loop index
      itype,        &! netCDF data type
      att_ival,     &! netCDF data type
      num_atts,     &! number of global attributes
      xtype

   logical (log_kind) :: &
      att_lval           ! temp space for logical attribute

   real (r4) ::     &
      att_rval       ! temp space for real attribute

   real (r8) ::     &
      att_dval       ! temp space for double attribute

   logical (log_kind) :: &
      attrib_error        ! error flag for reading attributes

!-----------------------------------------------------------------------
!
!  set the readonly flag in the data file descriptor
!
!-----------------------------------------------------------------------

   data_file%readonly = .true.

!-----------------------------------------------------------------------
!
!  open the netCDF file
!
!-----------------------------------------------------------------------

   path = trim(data_file%full_name)

   call io_pio_init('read', path, data_file%File)
   
!-----------------------------------------------------------------------
!
!  determine number of global file attributes
!
!-----------------------------------------------------------------------

   iostat = pio_inquire(data_file%File, nAttributes = num_atts)

!-----------------------------------------------------------------------
!
!  now read each attribute and set attribute values
!
!-----------------------------------------------------------------------

   do n=1,num_atts

      !***
      !*** get attribute name
      !***

      att_name = char_blank
      iostat = pio_inq_attname(data_file%File, PIO_GLOBAL, n, att_name)
      !***
      !*** check to see if name matches any of the standard file
      !*** attributes
      !***

      select case(trim(att_name))

      case('title')

         data_file%title = char_blank

         iostat = pio_inq_att(data_file%File, PIO_GLOBAL, 'title', &
                               xtype, nsize)

         if (iostat == pio_noerr) then
            if (nsize <= len(data_file%title)) then
               iostat = pio_get_att(data_file%File, PIO_GLOBAL, 'title', &
                                     data_file%title(1:nsize))
            else
               if (my_task == master_task) then
                  call document('open_read_netcdf', 'nsize', nsize)
                  call document('open_read_netcdf', 'len(data_file%title)', &
                                len(data_file%title))
                  write(stdout,*) 'string too short; not enough room to read title from ' /&
                                  &/ trim(path)
               endif
            endif
         endif

      case('history')
         
         data_file%history = char_blank
         iostat = pio_inq_attlen(data_file%File, PIO_GLOBAL, 'history', nsize)

         if (iostat == pio_noerr) then
            if (nsize <= len(data_file%history)) then
               iostat = pio_get_att(data_file%File, PIO_GLOBAL, 'history', &
                                     data_file%history(1:nsize))
            else
               if (my_task == master_task) then
                  call document('open_read_netcdf', 'nsize', nsize)
                  call document('open_read_netcdf', 'len(data_file%history)', &
                                len(data_file%history))
                  write(stdout,*) 'string too short; not enough room to read history attribute from ' /&
                                  &/ trim(path)
               endif
            endif
         endif

      case('conventions','Conventions','CONVENTIONS')

         data_file%conventions = char_blank
         iostat = pio_inq_att(data_file%File, PIO_GLOBAL,  trim(att_name), &
                                xtype, nsize)
         if (iostat == pio_noerr) then
            if (nsize <= len(data_file%conventions)) then
               iostat = pio_get_att(data_file%File, PIO_GLOBAL, trim(att_name), &
                                      data_file%conventions(1:nsize))
            else
               if (my_task == master_task) then
                  call document('open_read_netcdf', 'nsize', nsize)
                  call document('open_read_netcdf', 'len(data_file%conventions)', &
                                len(data_file%conventions))
                  write(stdout,*) 'string too short; not enough room to read conventions from ' /&
                                  &/ trim(path)
               endif
            endif
         endif

      case default

         !***
         !*** if does not match any of the standard file attributes
         !*** add the attribute to the datafile
         !***

         iostat = pio_inq_att(data_file%File, PIO_GLOBAL, trim(att_name), &
              itype,  nsize) 

         select case (itype)

         case (PIO_CHAR)
            work_line = char_blank
            if (nsize <= len(work_line)) then
               iostat = pio_get_att(data_file%File, PIO_GLOBAL, trim(att_name), &
                                     work_line(1:nsize))
            else
               if (my_task == master_task) then
                  call document('open_read_netcdf', 'nsize', nsize)
                  call document('open_read_netcdf', 'len(work_line)', &
                                len(work_line))
                  write(stdout,*) 'string too short; not enough room to read ' /&
                                  &/ trim(att_name) /&
                                  &/ ' from ' /&
                                  &/ trim(path)
               endif
            endif
            call add_attrib_file(data_file, trim(att_name), trim(work_line))

         case (PIO_INT)
            iostat = pio_get_att(data_file%File, PIO_GLOBAL, trim(att_name), &
                                  att_ival)

            if (att_name(1:4) == 'LOG_') then !*** attribute logical
               work_line = att_name
               work_line(1:4) = '    '
               att_name = adjustl(work_line)

               if (att_ival == 1) then
                  att_lval = .true.
               else
                  att_lval = .false.
               endif
               call add_attrib_file(data_file, trim(att_name), att_lval)
            else
               call add_attrib_file(data_file, trim(att_name), att_ival)
            endif

         case (PIO_REAL) 
            iostat = pio_get_att(data_file%File, PIO_GLOBAL, trim(att_name), &
                                  att_rval)
            call add_attrib_file(data_file, trim(att_name), att_rval)

         case (PIO_DOUBLE) 
            iostat = pio_get_att(data_file%File, PIO_GLOBAL, trim(att_name), &
                                  att_dval)
            call add_attrib_file(data_file, trim(att_name), att_dval)

         end select

      end select

   end do ! num_atts

!-----------------------------------------------------------------------
!EOC

 end subroutine open_read_netcdf

!***********************************************************************
!BOP
! !IROUTINE: open_netcdf
! !INTERFACE:

 subroutine open_netcdf(data_file)

! !INPUT/OUTPUT PARAMETERS:

   type (datafile), target, intent (inout)  :: data_file

! !DESCRIPTION:
!  This routine opens a data file for writing and
!  writes global file attributes.
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

   character (char_len) :: &
      path             ! temp to use for filename

   character (255) :: &
      work_line        ! temp to use for character manipulation

   integer (i4) ::  &
      iostat,       &! status flag for netCDF function calls
      itmp,         &! integer temp for equivalent logical attribute
      n,            &! loop index
      ncvals,       &! counter for number of character attributes
      nlvals,       &! counter for number of logical   attributes
      nivals,       &! counter for number of integer   attributes
      nrvals,       &! counter for number of real      attributes
      ndvals         ! counter for number of double    attributes

   logical (log_kind) :: &
      attrib_error       ! error flag for reading attributes

!-----------------------------------------------------------------------
!
!  open the netCDF file
!
!-----------------------------------------------------------------------

   path = trim(data_file%full_name)

   call io_pio_init(mode='write', filename=path, File=data_file%File, &
        clobber=.true., cdf64=luse_nf_64bit_offset)
   
   data_file%ldefine = .true.  ! file in netCDF define mode

!-----------------------------------------------------------------------
!
!  define global file attributes
!
!-----------------------------------------------------------------------

   attrib_error = .false.

   !*** standard attributes


   iostat = pio_put_att(data_file%File, PIO_GLOBAL, 'title', &
                         trim(data_file%title))

   iostat = pio_put_att(data_file%File, PIO_GLOBAL, 'history', &
                         trim(data_file%history))

   iostat = pio_put_att(data_file%File, PIO_GLOBAL, 'Conventions', &
                         trim(data_file%conventions))

   !*** additional attributes

   if (associated(data_file%add_attrib_cval)) then
      ncvals = size(data_file%add_attrib_cval)
   else
      ncvals = 0
   endif
   if (associated(data_file%add_attrib_lval)) then
      nlvals = size(data_file%add_attrib_lval)
   else
      nlvals = 0
   endif
   if (associated(data_file%add_attrib_ival)) then
      nivals = size(data_file%add_attrib_ival)
   else
      nivals = 0
   endif
   if (associated(data_file%add_attrib_rval)) then
      nrvals = size(data_file%add_attrib_rval)
   else
      nrvals = 0
   endif
   if (associated(data_file%add_attrib_dval)) then
      ndvals = size(data_file%add_attrib_dval)
   else
      ndvals = 0
   endif

   do n=1,ncvals
      work_line = data_file%add_attrib_cname(n)
      iostat = pio_put_att(data_file%File, PIO_GLOBAL, trim(work_line), &
                            trim(data_file%add_attrib_cval(n)))
   end do

   do n=1,nlvals
      work_line = 'LOG_'/&
                         &/data_file%add_attrib_lname(n)
      if (data_file%add_attrib_lval(n)) then
         itmp = 1
      else
         itmp = 0
      endif
      iostat = pio_put_att(data_file%File, PIO_GLOBAL, trim(work_line), &
                            itmp)
   end do

   do n=1,nivals
      work_line = data_file%add_attrib_iname(n)

      iostat = pio_put_att(data_file%File, PIO_GLOBAL, trim(work_line), &
                            data_file%add_attrib_ival(n))
   end do

   do n=1,nrvals
      work_line = data_file%add_attrib_rname(n)

      iostat = pio_put_att(data_file%File, PIO_GLOBAL, trim(work_line), &
                            data_file%add_attrib_rval(n))
   end do

   do n=1,ndvals
      work_line = data_file%add_attrib_dname(n)
      iostat = pio_put_att(data_file%File, PIO_GLOBAL, trim(work_line), &
                            data_file%add_attrib_dval(n))
   end do

   if (attrib_error) call exit_POP(sigAbort, &
                                   'Error writing file attributes')

!-----------------------------------------------------------------------
!EOC

 end subroutine open_netcdf

!***********************************************************************
!BOP
! !IROUTINE: close_netcdf
! !INTERFACE:

 subroutine close_netcdf(data_file)

! !INPUT/OUTPUT PARAMETERS:

   type (datafile), intent (inout)  :: data_file
   integer :: iostat
! !DESCRIPTION:
!  This routine closes an open netcdf data file.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  close a data file
!
!-----------------------------------------------------------------------

   call pio_closefile(data_file%File)

!-----------------------------------------------------------------------
!EOC

 end subroutine close_netcdf

!***********************************************************************
!BOP
! !IROUTINE: sync_netcdf
! !INTERFACE:

 subroutine sync_netcdf(data_file)

! !INPUT/OUTPUT PARAMETERS:

   type (datafile), intent (inout)  :: data_file

! !DESCRIPTION:
!  This routine uses pio_syncfile to flush an open netcdf data file.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  close a data file
!
!-----------------------------------------------------------------------

   call pio_syncfile(data_file%File)

!-----------------------------------------------------------------------
!EOC

 end subroutine sync_netcdf

!***********************************************************************
!BOP
! !IROUTINE: define_field_netcdf
! !INTERFACE:

 subroutine define_field_netcdf(data_file, io_field)

! !DESCRIPTION:
!  This routine defines an io field for a netCDF file.
!  When reading a file, the define routine will attempt to fill an 
!  io field structure with meta-data information from the netCDF file.
!  When writing a file, it calls the appropriate netCDF routines
!  to define all the field attributes and assign a field id.
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

   type (datafile), target, intent (inout)  :: &
      data_file       ! data file in which field contained

   type (io_field_desc), intent (inout) :: &
      io_field        ! field descriptor for this field

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   character (80) :: &
      work_line,     &! workspace for manipulating input string
      comp_line,     &! comparison string
      att_name        ! attribute name

   integer (i4) :: &
      iostat,      &! status flag for netCDF calls
      varid,       &! variable id for field
      ndims,       &! number of dimensions
      dimid,       &! dimension id
      n,           &! loop index
      ncount,      &! num additional attributes
      nsize,       &! length of character strings
      itype,       &! netCDF data type
      num_atts,    &! number of variable attributes
      att_ival,    &! temp for integer attribute
      ncvals,      &! counter for number of character attributes
      nlvals,      &! counter for number of logical   attributes
      nivals,      &! counter for number of integer   attributes
      nrvals,      &! counter for number of real      attributes
      ndvals        ! counter for number of double    attributes

   logical (log_kind) ::    &
      att_lval      ! temp for logical attribute

   real (r4) ::            &
      att_rval              ! temp for real attribute

   real (r8) ::    &
      att_dval      ! temp for double attribute

   logical (log_kind) :: &
      define_error       ! error flag

   integer (i4) :: xtype

   define_error = .false.
   
   data_file%ldefine = .true.  ! file in netCDF define mode

!-----------------------------------------------------------------------
!
!  for input files, get the variable id and determine number of field
!  attributes
!
!-----------------------------------------------------------------------

   call pio_seterrorhandling(data_file%File, PIO_BCAST_ERROR)

   if (data_file%readonly) then

      ! Note that currently a lot of pio inquire functions need a 
      ! netcdf varid and not a pio vardesc. Currently pio_inq_varnatts 
      ! can only be accessed through a pio vardesc.  

      iostat = pio_inq_varid(data_file%File, io_field%short_name, io_field%id) 
      if (iostat /= pio_noerr) &
           call exit_POP(sigAbort,'Error in getting varid for netCDF field')

      iostat = pio_inq_varid(data_file%File, io_field%short_name, io_field%varDesc) 
      if (iostat /= pio_noerr) &
           call exit_POP(sigAbort,'Error in getting varDesc for netCDF field')

      iostat = pio_inq_varnatts(data_file%File, io_field%varDesc, nAtts=num_atts)
      if (iostat /= pio_noerr) &
           call exit_POP(sigAbort,'Error getting attrib count for netCDF field')

      !***
      !*** for each attribute, define standard attributes or add
      !*** attribute to io_field
      !***

      do n=1,num_atts

         !***
         !*** get attribute name
         !***

         att_name = char_blank
         iostat = pio_inq_attname(data_file%File, io_field%id, n, att_name)
         if (iostat /= pio_noerr) &
            call exit_POP(sigAbort,'Error getting netCDF field attribute name')
   
         !***
         !*** check to see if name matches any of the standard field
         !*** attributes
         !***

         select case(trim(att_name))

         case('long_name')

            io_field%long_name = char_blank

            iostat = pio_inq_att(data_file%File, io_field%id, 'long_name', &
	                          xtype, nsize)
            if (iostat == pio_noerr) then
               if (nsize <= len(io_field%long_name)) then
                  iostat = pio_get_att(data_file%File, io_field%id, 'long_name', &
                                        io_field%long_name(1:nsize))
               else
                  if (my_task == master_task) then
                     call document('define_field_netcdf', 'nsize', nsize)
                     call document('define_field_netcdf', 'len(io_field%long_name)', &
                                   len(io_field%long_name))
                     write(stdout,*) 'string too short; not enough room to read long_name of ' /&
                                     &/ trim(io_field%short_name) /&
                                     &/ ' from ' /&
                                     &/ trim(data_file%full_name)
                  end if
               endif
            endif
            if (iostat /= pio_noerr) then
               call exit_POP(sigAbort, &
                   'Error reading long_name from netCDF file')
            endif

         case('units')

            io_field%units = char_blank

            iostat = pio_inq_att(data_file%File, io_field%id, 'units', &
	                          xtype, nsize)

            if (iostat == pio_noerr) then
               if (nsize <= len(io_field%units)) then
                  iostat = pio_get_att(data_file%File, io_field%id, 'units', &
	                                io_field%units(1:nsize))
               else
                  if (my_task == master_task) then
                     call document('define_field_netcdf', 'nsize', nsize)
                     call document('define_field_netcdf', 'len(io_field%units)', &
                                   len(io_field%units))
                     write(stdout,*) 'string too short; not enough room to read units of ' /&
                                     &/ trim(io_field%short_name) /&
                                     &/ ' from ' /&
                                     &/ trim(data_file%full_name)
                  end if
               endif
            endif

         case('coordinates')

            io_field%coordinates = char_blank

            iostat = pio_inq_att(data_file%File, io_field%id, 'coordinates', &
	                          xtype, nsize)

            if (iostat == pio_noerr) then
               if (nsize <= len(io_field%coordinates)) then
                  iostat = pio_get_att(data_file%File, io_field%id, 'coordinates', &
                                        io_field%coordinates(1:nsize))
                  if (iostat /= pio_noerr) then
                     call exit_POP(sigAbort, &
                          'Error reading coordinates from netCDF file')
            endif
               else
                  if (my_task == master_task) then
                     call document('define_field_netcdf', 'nsize', nsize)
                     call document('define_field_netcdf', 'len(io_field%coordinates)', &
                                   len(io_field%coordinates))
                     write(stdout,*) 'string too short; not enough room to read coordinates of ' /&
                                     &/ trim(io_field%short_name) /&
                                     &/ ' from ' /&
                                     &/ trim(data_file%full_name)
                  endif
               endif
            endif


         case('grid_loc')

            io_field%grid_loc = '    '

            iostat = pio_inq_att(data_file%File, io_field%id, 'grid_loc', &
                                  xtype, nsize)

            if (iostat == pio_noerr) then
               if (nsize <= len(io_field%grid_loc)) then
                  iostat = pio_get_att(data_file%File, io_field%id, 'grid_loc', &
                                        io_field%grid_loc(1:nsize))
               else
                  call document('define_field_netcdf', 'nsize', nsize)
                  call document('define_field_netcdf', 'len(io_field%grid_loc)', &
                                len(io_field%grid_loc))
                  write(stdout,*) 'string too short; not enough room to read grid_loc of ' /&
                                  &/ trim(io_field%short_name) /&
                                  &/ ' from ' /&
                                  &/ trim(data_file%full_name)
               endif
            endif

            if (iostat /= pio_noerr) then
               call exit_POP(sigAbort, &
                   'Error reading grid_loc from netCDF file')
            endif

         case('valid_range')

            iostat = pio_get_att(data_file%File, io_field%id, &
                                  'valid_range',   &
                                  io_field%valid_range)
            if (iostat /= pio_noerr) then
               call exit_POP(sigAbort, &
                   'Error reading valid_range from netCDF file')
            endif

         case default

            !***
            !*** if does not match any of the standard file attributes
            !*** add the attribute to the datafile
            !***

            iostat = pio_inq_att(data_file%File, io_field%id, trim(att_name), &
                                 itype,  nsize) 

            if (iostat /= pio_noerr) then
               call exit_POP(sigAbort, &
                   'Error reading netCDF file attribute')
            endif
   
            select case (itype)

            case (PIO_CHAR)
               work_line = char_blank
               if (nsize <= len(work_line)) then
                  iostat = pio_get_att(data_file%File, io_field%id, trim(att_name), &
                                        work_line(1:nsize))
               else
                  if (my_task == master_task) then
                     call document('define_field_netcdf', 'nsize', nsize)
                     call document('define_field_netcdf', 'len(work_line)', &
                                len(work_line))
                     write(stdout,*) 'string too short; not enough room to read ' /&
                                     &/ trim(att_name) /&
                                     &/ ' of ' /&
                                     &/ trim(io_field%short_name) /&
                                     &/ ' from ' /&
                                     &/ trim(data_file%full_name)
                  endif
               endif
               if (iostat /= pio_noerr) then
                  call exit_POP(sigAbort, &
                                'Error reading netCDF file attribute')
               endif

               call add_attrib_io_field(io_field, trim(att_name), &
                                                  trim(work_line))

            case (PIO_INT) !*** both integer and logical attributes
               iostat = pio_get_att(data_file%File, io_field%id, &
                                     trim(att_name), att_ival)
               if (iostat /= pio_noerr) then
                  call exit_POP(sigAbort, &
                                'Error reading netCDF file attribute')
               endif
   
               if (att_name(1:4) == 'LOG_') then !*** attribute logical
                  work_line = att_name
                  work_line(1:4) = '    '
                  att_name = adjustl(work_line)

                  if (att_ival == 1) then
                     att_lval = .true.
                  else
                     att_lval = .false.
                  endif
                  call add_attrib_file(data_file, trim(att_name), &
                                                  att_lval)

               else
                  call add_attrib_file(data_file, trim(att_name), &
                                                  att_ival)
               endif

            case (PIO_REAL)
               iostat = pio_get_att(data_file%File, io_field%id, &
                                     trim(att_name), att_rval)
               if (iostat /= pio_noerr) then
                  call exit_POP(sigAbort, &
                                'Error reading netCDF file attribute')
               endif

               call add_attrib_io_field(io_field, trim(att_name), &
                                                  att_rval)

            case (PIO_DOUBLE)
               iostat = pio_get_att(data_file%File, io_field%id, &
                                     trim(att_name), att_dval)
               if (iostat /= pio_noerr) then
                  call exit_POP(sigAbort, &
                                'Error reading netCDF file attribute')
               endif
   
               call add_attrib_io_field(io_field, trim(att_name), &
                                                  att_dval)
   
            end select

         end select

      end do ! num_atts

!-----------------------------------------------------------------------
!
!  for output files, need to define everything
!  make sure file is in define mode
!
!-----------------------------------------------------------------------

   else ! output file

      if (.not. data_file%ldefine) &
        call exit_POP(sigAbort, &
                      'attempt to define field but not in define mode')

!-----------------------------------------------------------------------
!
!     define the dimensions
!
!-----------------------------------------------------------------------

      ndims = io_field%nfield_dims

      do n = 1,ndims
         dimid = 0

         !*** check to see whether already defined
         
         iostat = pio_inq_dimid(data_file%file,                                 &
                                name=trim(io_field%field_dim(n)%name),&
                                dimid=dimid)

         if (iostat /= PIO_NOERR) then ! dimension not yet defined
            iostat = pio_def_dim (data_file%File,                          &
                             name=trim(io_field%field_dim(n)%name), &
                             len=io_field%field_dim(n)%length,      &
                             dimid=io_field%field_dim(n)%id)
         else
            io_field%field_dim(n)%id = dimid
         end if
      end do

!-----------------------------------------------------------------------
!
!        now define the field
!
!-----------------------------------------------------------------------

      !*** check to see whether field of this name already defined.
      
      iostat = pio_inq_varid(data_file%File, trim(io_field%short_name), varid)
      
      if (iostat /= PIO_NOERR) then ! variable was not yet defined

         if (associated (io_field%field_r_1d).or. &
             associated (io_field%field_r_2d).or. &
             associated (io_field%field_r_3d)) then
            iostat = pio_def_var (data_file%File,                            &
                                  name=trim(io_field%short_name),  &
                                  type=PIO_REAL,                 &
                dimids=(/ (io_field%field_dim(n)%id, n=1,ndims) /),&
                                  varDesc=io_field%varDesc)

         else if (            io_field%nfield_dims == c0) then 
            ! do not supply optional dimids for scalars
            iostat = pio_def_var (data_file%File,                           &
                                  name=trim(io_field%short_name), &
                                  type=PIO_DOUBLE,               &
                                  varDesc=io_field%varDesc)
         else if (associated (io_field%field_d_1d).or. &
                  associated (io_field%field_d_2d).or. &
                  associated (io_field%field_d_3d)) then
            iostat = pio_def_var (data_file%File,                           &
                                  name=trim(io_field%short_name), &
                                  type=PIO_DOUBLE,               &
               dimids=(/ (io_field%field_dim(n)%id, n=1,ndims) /),&
                                  varDesc=io_field%varDesc)
         else if (associated (io_field%field_i_1d).or. &
                  associated (io_field%field_i_2d).or. &
                  associated (io_field%field_i_3d)) then
            iostat = pio_def_var (data_file%File,                           &
                                  name=trim(io_field%short_name), &
                                  type=PIO_INT,                  &
               dimids=(/ (io_field%field_dim(n)%id, n=1,ndims) /),&
                                  varDesc=io_field%varDesc)
         else
            define_error = .true.
         end if
         if (iostat /= pio_noerr) define_error = .true.

      end if

      ! Now get a valid netcdf varid for the variable and fill in 
      ! the io_field%id setting 

      iostat = pio_inq_varid(data_file%File, trim(io_field%short_name), varid)
      io_field%id = varid
      if (iostat /= PIO_NOERR) define_error = .true.

      iostat = pio_inq_varid(data_file%File, trim(io_field%short_name), io_field%vardesc)
      if (iostat /= pio_noerr) define_error = .true.

      if (define_error) then
         write(stdout,*) '(define_field_netcdf) ', trim(io_field%short_name)
         call exit_POP(sigAbort, 'Error defining netCDF field')
      endif

!-----------------------------------------------------------------------
!
!     Now define the field attributes
!
!-----------------------------------------------------------------------

      !*** long_name

      if (io_field%long_name /= char_blank) then
         iostat = pio_inq_att(data_file%File, varid, 'long_name', &
                               xtype, nsize)
         if (iostat /= PIO_NOERR) then ! attrib probably not defined
            iostat = pio_put_att(data_file%File, varid=varid, &
                                  name='long_name',       &
                                  value=trim(io_field%long_name))
            if (iostat /= PIO_NOERR) define_error = .true.
         end if
      endif

      !*** units

      if (io_field%units /= char_blank) then
         iostat = pio_inq_att(data_file%File, varid, 'units', &
                               xtype, nsize)
         if (iostat /= PIO_NOERR) then ! attrib probably not defined
            iostat = pio_put_att(data_file%File, varid=varid,       &
                                  name='units',           &
                                  value=trim(io_field%units))
            if (iostat /= PIO_NOERR) define_error = .true.
         end if
      endif

      !*** coordinates

      if (io_field%coordinates /= char_blank) then
         iostat = pio_inq_att(data_file%File, varid, 'coordinates', &
                               xtype, nsize)
         if (iostat /= PIO_NOERR) then ! attrib probably not defined
            iostat = pio_put_att(data_file%File, varid=varid, &
                                  name='coordinates',           &
                                  value=trim(io_field%coordinates))
            if (iostat /= PIO_NOERR) define_error = .true.
         end if
      endif

      !*** grid_loc

      if (io_field%grid_loc /= '    ') then
         iostat = pio_inq_att(data_file%File, varid, 'grid_loc', &
                               xtype, nsize)
         if (iostat /= PIO_NOERR) then ! attrib probably not defined
            iostat = pio_put_att(data_file%File, varid=varid, &
                                  name='grid_loc',        &
                                  value=io_field%grid_loc)
            if (iostat /= PIO_NOERR) define_error = .true.
         end if
      endif


      !*** valid_range(1:2)

      if (any(io_field%valid_range /= undefined)) then
         iostat = pio_inq_att(data_file%File, varid, 'valid_range', &
                               xtype, nsize)
         if (iostat /= PIO_NOERR) then ! attrib probably not yet defined
            iostat = pio_put_att(data_file%File, varid=varid, &
                                  name='valid_range',       &
                                  value=io_field%valid_range(:))
            if (iostat /= PIO_NOERR) define_error = .true.
         end if
      endif

      !*** additional attributes if defined

      ncvals = 0
      nlvals = 0
      nivals = 0
      nrvals = 0
      ndvals = 0
      if (associated(io_field%add_attrib_cval)) &
         ncvals = size(io_field%add_attrib_cval)
      if (associated(io_field%add_attrib_lval)) &
         nlvals = size(io_field%add_attrib_lval)
      if (associated(io_field%add_attrib_ival)) &
         nivals = size(io_field%add_attrib_ival)
      if (associated(io_field%add_attrib_rval)) &
         nrvals = size(io_field%add_attrib_rval)
      if (associated(io_field%add_attrib_dval)) &
         ndvals = size(io_field%add_attrib_dval)

      do n=1,ncvals
         iostat = pio_put_att(data_file%File, varid=varid,             &
                      name=trim(io_field%add_attrib_cname(n)), &
                      value=trim(io_field%add_attrib_cval(n)))
         if (iostat /= PIO_NOERR) define_error = .true.
      end do

      do n=1,nlvals
         work_line = 'LOG_'/&
                            &/trim(io_field%add_attrib_lname(n))
         iostat = pio_put_att(data_file%File, varid=varid,             &
                      name=trim(work_line),                        &
                      value=io_field%add_attrib_ival(n))
         if (iostat /= PIO_NOERR) define_error = .true.
      end do

      do n=1,nivals
         iostat = pio_put_att(data_file%File, varid=varid,             &
                      name=trim(io_field%add_attrib_iname(n)),     &
                      value=io_field%add_attrib_ival(n))
         if (iostat /= PIO_NOERR) define_error = .true.
      end do

      do n=1,nrvals
         iostat = pio_put_att(data_file%file, varid=varid,             &
                      name=trim(io_field%add_attrib_rname(n)),     &
                      value=io_field%add_attrib_rval(n))
         if (iostat /= PIO_NOERR) define_error = .true.
      end do

      do n=1,ndvals
         iostat = pio_put_att(data_file%File, varid=varid,             &
                      name=trim(io_field%add_attrib_dname(n)),     &
                      value=io_field%add_attrib_dval(n))
         if (iostat /= PIO_NOERR) define_error = .true.
      end do

      if (define_error) call exit_POP(sigAbort, &
                        'Error adding attributes to field')

   endif ! input/output file

   call pio_seterrorhandling(data_file%File, PIO_INTERNAL_ERROR)

!-----------------------------------------------------------------------
!EOC

 end subroutine define_field_netcdf

!***********************************************************************
!BOP
! !IROUTINE: write_field_netcdf
! !INTERFACE:

 subroutine write_field_netcdf(data_file, io_field)

! !DESCRIPTION:
!  This routine writes a field to a netCDF data file.
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

   type (datafile), target, intent (inout)  :: &
      data_file             ! file to which field will be written

   type (io_field_desc), intent (inout) :: &
      io_field              ! field to write to file

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (i4) :: &
      iostat,       &! netCDF status flag
      ndims,        &! dimension index 	
      k,n            ! loop counters

   logical (log_kind) :: &
      write_error         ! error flag

   integer (i4), dimension(1) ::  &
      start,count        ! dimension quantities for netCDF

!-----------------------------------------------------------------------
!
!  exit define mode if necessary
!
!-----------------------------------------------------------------------

   write_error = .false.
   

   if (data_file%ldefine) then
      iostat = pio_enddef(data_file%File)
      data_file%ldefine = .false.
   endif

!-----------------------------------------------------------------------
!
!  make sure field has been defined
!
!-----------------------------------------------------------------------

   if (io_field%id == 0) then
      call exit_POP(sigAbort,'Attempt to write undefined field in netCDF write')
   end if

!-----------------------------------------------------------------------
!
!  write data based on type
!
!-----------------------------------------------------------------------

   if (trim(io_field%short_name) == 'time') then
      ndims    = io_field%nfield_dims
      start(1) = io_field%field_dim(ndims)%start
      count(1) = 1	
      iostat = pio_put_var(data_file%File, varid=io_field%id, start=start(:), count=count(:), &
                           ival=io_field%field_d_1d) 
      if (iostat /= pio_noerr) then
         call document('write_field_netcdf', 'short_name', io_field%short_name)
         call exit_POP(sigAbort,'Error writing field time to netCDF file')
      end if
      RETURN
   end if

   ! Set the unlimited dimension pointer for the variable

   if (io_field%set_iodesc) then
      if (associated(io_field%field_r_3d)) then
         call io_pio_initdecomp(PIO_REAL,   ndim3=io_field%field_dim(3)%length, &
              kdim3=size(io_field%field_r_3d,3), iodesc=io_field%ioDesc)
      else if (associated(io_field%field_d_3d)) then
         call io_pio_initdecomp(PIO_DOUBLE, ndim3=io_field%field_dim(3)%length, &
              kdim3=size(io_field%field_d_3d,3), iodesc=io_field%ioDesc)
      else if (associated(io_field%field_i_3d)) then
         call io_pio_initdecomp(PIO_INT,    ndim3=io_field%field_dim(3)%length, &
              kdim3=size(io_field%field_i_3d,3), iodesc=io_field%ioDesc)
      else if (associated(io_field%field_r_2d)) then
         call io_pio_initdecomp(PIO_REAL,   ndim3=0, kdim3=0, iodesc=io_field%ioDesc)
      else if (associated(io_field%field_d_2d)) then
         call io_pio_initdecomp(PIO_DOUBLE, ndim3=0, kdim3=0, iodesc=io_field%ioDesc)
      else if (associated(io_field%field_i_2d)) then
         call io_pio_initdecomp(PIO_INT,    ndim3=0, kdim3=0, iodesc=io_field%ioDesc)
      end if
      io_field%set_iodesc = .false.
   end if
      
   if (io_field%set_ioFrame) then	
      ndims = io_field%nfield_dims
      call pio_setframe(io_field%vardesc, int(io_field%field_dim(ndims)%start,kind=PIO_OFFSET))
   end if

   if (associated(io_field%field_r_3d)) then

      call pio_write_darray(data_file%File, io_field%vardesc, io_field%iodesc, &
                            io_field%field_r_3d, iostat)

   else if (associated(io_field%field_r_2d)) then

      call pio_write_darray(data_file%File, io_field%vardesc, io_field%iodesc, &
                            io_field%field_r_2d, iostat)

   else if (associated(io_field%field_r_1d)) then

      ! 1d vectors are not distributed to blocks
      iostat = pio_put_var(data_file%File, io_field%vardesc, io_field%field_r_1d)

   else if (associated(io_field%field_d_3d)) then

      call pio_write_darray(data_file%File, io_field%vardesc, io_field%iodesc, &
                            io_field%field_d_3d, iostat)

   else if (associated(io_field%field_d_2d)) then

      call pio_write_darray(data_file%File, io_field%vardesc, io_field%iodesc, &
                            io_field%field_d_2d, iostat)

   else if (associated(io_field%field_d_1d)) then

      ! 1d vectors are not distributed to blocks; no need for gather_global
      iostat = pio_put_var(data_file%File, io_field%vardesc, io_field%field_d_1d)
 
   else if (io_field%nfield_dims == c0) then

      ! scalars are not distributed to blocks; no need for gather_global
      ! for now, all scalars are r8   and are not pointers or targets 
      iostat = pio_put_var(data_file%File, io_field%vardesc, io_field%field_d_0d)

   else if (associated(io_field%field_i_3d)) then

      call pio_write_darray(data_file%File, io_field%vardesc, io_field%iodesc, &
                            io_field%field_i_3d, iostat)

   else if (associated(io_field%field_i_2d)) then

      call pio_write_darray(data_file%File, io_field%vardesc, io_field%iodesc, &
                            io_field%field_i_2d, iostat)

   else if (associated(io_field%field_i_1d)) then

      ! 1d vectors are not distributed to blocks; no need for gather_global
      iostat = pio_put_var(data_file%File, io_field%vardesc, io_field%field_i_1d)

   else
      call exit_POP(sigAbort, &
                    'No field associated for writing to netCDF')
   end if

   if (iostat /= pio_noerr) then
      write_error = .true.
   endif

   if (write_error) then
      call document('write_field_netcdf', 'short_name', io_field%short_name)
      call exit_POP(sigAbort, &
                    'Error writing field to netCDF file')
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine write_field_netcdf

!***********************************************************************
!BOP
! !IROUTINE: read_field_netcdf
! !INTERFACE:

 subroutine read_field_netcdf(data_file, io_field)

! !DESCRIPTION:
!  This routine reads a field from a netcdf input file.
!
! !REVISION HISTORY:
!  same as module
!
! !USES

   use POP_FieldMod	           
   use POP_GridHorzMod
   use Pop_HaloMod

! !INPUT/OUTPUT PARAMETERS:

   type (datafile), target, intent (inout)  :: &
      data_file              ! file from which to read field

   type (io_field_desc), intent (inout) :: &
      io_field               ! field to be read

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (i4) :: &
      iostat,      &! netCDF status flag
      k,n           ! loop counters
 
   character(len=8) :: fieldtype, fieldloc

   integer (POP_i4) ::  errorCode           ! returned error code

   logical (log_kind) :: lhalo_update

!-----------------------------------------------------------------------
!
!  make sure field has been defined
!
!-----------------------------------------------------------------------

   
   iostat = pio_inq_varid(data_file%File, trim(io_field%short_name), io_field%varDesc)

!-----------------------------------------------------------------------
!
!  if no boundary update type defined, assume center location scalar
!
!-----------------------------------------------------------------------

   if (io_field%field_loc == field_loc_unknown) then
      io_field%field_loc  = field_loc_center
      io_field%field_type = field_type_scalar
   endif

!-----------------------------------------------------------------------
!
!  read data based on type
!
!-----------------------------------------------------------------------

   if (io_field%set_iodesc) then
      call pio_setframe(io_field%vardesc, int(1,kind=PIO_OFFSET))
      if (associated(io_field%field_r_3d)) then
         call io_pio_initdecomp(PIO_REAL, ndim3=io_field%field_dim(3)%length, &
              kdim3=io_field%field_dim(3)%length, iodesc=io_field%ioDesc)
      else if (associated(io_field%field_d_3d)) then
         call io_pio_initdecomp(PIO_DOUBLE, ndim3=io_field%field_dim(3)%length, &
              kdim3=io_field%field_dim(3)%length, iodesc=io_field%ioDesc)
      else if (associated(io_field%field_i_3d)) then
         call io_pio_initdecomp(PIO_INT, ndim3=io_field%field_dim(3)%length, &
              kdim3=io_field%field_dim(3)%length, iodesc=io_field%ioDesc)
      else if (associated(io_field%field_r_2d)) then
         call io_pio_initdecomp(PIO_REAL, ndim3=0, kdim3=0, iodesc=io_field%ioDesc)
      else if (associated(io_field%field_d_2d)) then
         call io_pio_initdecomp(PIO_DOUBLE, ndim3=0, kdim3=0, iodesc=io_field%ioDesc)
      else if (associated(io_field%field_i_2d)) then
         call io_pio_initdecomp(PIO_INT, ndim3=0, kdim3=0, iodesc=io_field%ioDesc)
      end if
      io_field%set_iodesc = .false.
   end if

   ! Set values for halo updates if needed
   if (io_field%field_loc == field_loc_center) then
      fieldLoc = POP_gridHorzLocCenter
   else if (io_field%field_loc == field_loc_NEcorner) then
      fieldLoc = POP_gridHorzLocNECorner
   else if (io_field%field_loc == field_loc_Nface) then
      fieldLoc = POP_gridHorzLocNface
   else if (io_field%field_loc == field_loc_Eface) then
      fieldLoc = POP_gridHorzLocEface
   end if

   if (io_field%field_type == field_type_vector) then
      fieldType = POP_fieldKindVector
   else if (io_field%field_type == field_type_scalar) then
      fieldType = POP_fieldKindScalar
   else if (io_field%field_type == field_type_angle) then
      fieldType = POP_fieldKindAngle
   else if (io_field%field_type == field_type_noupdate) then
      fieldType = POP_fieldKindNoUpdate
   else
      call exit_POP(sigAbort, 'read_field_netcdf field_type is not supported')
   end if

   ! Currently halo update is not supported for tripole grid
   if (ltripole_grid) then
      if (io_field%field_type == field_type_noupdate .or. &
	  io_field%field_loc  == field_loc_noupdate) then
         lhalo_update = .false.
      else
         lhalo_update = .true.
      end if
   else 
      if (io_field%field_loc == field_loc_noupdate .or. &
	  io_field%field_loc == field_loc_unknown) then
    	 lhalo_update = .false.
      else
         lhalo_update = .true.
      end if
   end if

   if (associated(io_field%field_r_3d)) then

      call pio_read_darray(data_file%File, io_field%varDesc, io_field%iodesc, &
                           io_field%field_r_3d(:,:,:,:), iostat)

      if (lhalo_update) then
         call POP_HaloUpdate(array=io_field%field_r_3d(:,:,:,:),       &
                             halo=POP_haloClinic,                      &
                             fieldLoc=FieldLoc,                        &
                             fieldKind=FieldType, errorCode=errorCode, &
                             fillValue=0.0_POP_r4)
         if (errorCode /= POP_Success) then
            call exit_POP(sigAbort, &
               'read_field_netcdf: error updating halo for field_r_3d')
         endif
      end if

   else if (associated(io_field%field_r_2d)) then

      call pio_read_darray(data_file%File, io_field%varDesc, io_field%iodesc, &
                           io_field%field_r_2d, iostat)

      if (lhalo_update) then
         call POP_HaloUpdate(array=io_field%field_r_2d(:,:,:),         &
                             halo=POP_haloClinic,                      &
                             fieldLoc=FieldLoc,                        &
                             fieldKind=FieldType, errorCode=errorCode, &
                             fillValue=0.0_POP_r4)
         if (errorCode /= POP_Success) then
            call exit_POP(sigAbort, &
               'read_field_netcdf: error updating halo for field_r_2d')
         endif
      end if

   else if (associated(io_field%field_r_1d)) then

      ! 1d vectors are not distributed to blocks; therefore, no scatter_global needed
      iostat = pio_get_var (data_file%File,io_field%varDesc,&
                            io_field%field_r_1d)

   else if (associated(io_field%field_r_1d)) then

      ! scalars are not distributed to blocks; therefore, no scatter_global needed
      iostat = pio_get_var (data_file%File,io_field%varDesc, &
                            io_field%field_r_0d)

   else if (associated(io_field%field_d_3d)) then

      call pio_read_darray(data_file%File, io_field%varDesc, io_field%ioDesc, &
                           io_field%field_d_3d, iostat)

      if (lhalo_update) then
         call POP_HaloUpdate(array=io_field%field_d_3d(:,:,:,:),       &
                             halo=POP_haloClinic,                      &
                             fieldLoc=FieldLoc,                        &
                             fieldKind=FieldType, errorCode=errorCode, &
                             fillValue=0.0_POP_r8)
         if (errorCode /= POP_Success) then
            call exit_POP(sigAbort, &
               'read_field_netcdf: error updating halo for field_d_3d')
         endif
      end if

   else if (associated(io_field%field_d_2d)) then

      call pio_read_darray(data_file%File, io_field%varDesc, io_field%ioDesc, &
                           io_field%field_d_2d, iostat)

      if (lhalo_update) then
         call POP_HaloUpdate(array=io_field%field_d_2d(:,:,:),         &
                             halo=POP_haloClinic,                      &
                             fieldLoc=FieldLoc,                        &
                             fieldKind=FieldType, errorCode=errorCode, &
                             fillValue=0.0_POP_r8)
         if (errorCode /= POP_Success) then
            call exit_POP(sigAbort, &
               'read_field_netcdf: error updating halo for field_d_2d')
         endif
      end if

   else if (associated(io_field%field_d_1d)) then

      ! 1d vectors are not distributed to blocks; therefore, no scatter_global needed
      iostat = pio_get_var (data_file%File,io_field%varDesc, &
                            io_field%field_d_1d)

   else if (associated(io_field%field_d_1d)) then

      ! scalars are not distributed to blocks; therefore, no scatter_global needed
      iostat = pio_get_var (data_file%File, io_field%varDesc, &
                            io_field%field_d_0d)

   else if (associated(io_field%field_i_3d)) then

      call pio_read_darray(data_file%File, io_field%varDesc, io_field%ioDesc, &
                           io_field%field_i_3d, iostat)

      if (lhalo_update) then
         call POP_HaloUpdate(array=io_field%field_i_3d(:,:,:,:),       &
                             halo=POP_haloClinic,                      &
                             fieldLoc=FieldLoc,                        &
                             fieldKind=FieldType, errorCode=errorCode, &
                             fillValue=0_POP_i4)
         if (errorCode /= POP_Success) then
            call exit_POP(sigAbort, &
               'read_field_netcdf: error updating halo for field_i_3d')
         endif
      end if

   else if (associated(io_field%field_i_2d)) then

      call pio_read_darray(data_file%File, io_field%varDesc, io_field%ioDesc, &
                           io_field%field_i_2d, iostat)

      if (lhalo_update) then
         call POP_HaloUpdate(array=io_field%field_i_2d(:,:,:),         &
                             halo=POP_haloClinic,                      &
                             fieldLoc=FieldLoc,                        &
                             fieldKind=FieldType, errorCode=errorCode, &
                             fillValue=0_POP_i4)
         if (errorCode /= POP_Success) then
            call exit_POP(sigAbort, &
               'read_field_netcdf: error updating halo for field_i_2d')
         endif
      end if

   else if (associated(io_field%field_i_1d)) then

      ! 1d vectors are not distributed to blocks; therefore, no scatter_global needed
      iostat = pio_get_var (data_file%File,io_field%varDesc, &
                            io_field%field_i_1d)

   else if (associated(io_field%field_i_1d)) then

      ! scalars are not distributed to blocks; therefore, no scatter_global needed
      iostat = pio_get_var (data_file%File, io_field%varDesc, &
                            io_field%field_i_0d)
   else
      call exit_POP(sigAbort, &
                    'No field associated for reading from netCDF')
   end if

!-----------------------------------------------------------------------
!EOC

 end subroutine read_field_netcdf

!***********************************************************************
!BOP
! !IROUTINE: define_nstd_netcdf
! !INTERFACE:

 subroutine define_nstd_netcdf(data_file,ndims,io_dims,field_id,         &
                                 short_name,long_name,units,coordinates, &
                                 fill_value,method_string,nftype)

! !DESCRIPTION:
!  This routine defines the nonstandard CCSM time-averaged diagnostic fields
!  on nonstandard grids: MOC, N_HEAT, and N_SALT
!  This routine is totally CCSM-specific
!
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (datafile), target, intent (inout)  :: &
      data_file       ! data file in which field contained

   real (r4), intent (in)  ::  &
      fill_value

   integer (int_kind), intent(in) ::  &
      ndims                ! number of dimensions for nonstandard field

   character (*), intent (in)  ::  &
      short_name,                  &
      long_name,                   &
      units,                       &
      coordinates,                 &
      nftype,                      &
      method_string
    
! !INPUT/OUTPUT PARAMETERS:

   type (io_dim), dimension(:), intent (inout) :: &
      io_dims         

   integer (i4), intent (inout) :: &
      field_id                      ! variable id 

   optional :: coordinates,fill_value,nftype,method_string

!EOP
!BOP
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (i4) ::    &
      iostat,         &! status flag for netCDF calls
      n,              &! loop index
      dimid,          &! dimension id
      xtype,          &
      len	


   logical (log_kind) :: &
      define_error        ! error flag

   ! Note this is a local variable since it is only used to satisfy 
   ! pio_def_var interface needs. Only pio_put_var is used for 
   ! non-standard variablesand all output is from the master processor.
   ! So vardesc is never directly.

   type(var_desc_t) :: &  
      vardesc	

   define_error = .false.
   

!-----------------------------------------------------------------------
!
!  make sure file has been opened and is in define mode
!
!-----------------------------------------------------------------------

   call check_definemode (data_file, 'define_nstd_netcdf')

!-----------------------------------------------------------------------
!
!  define the dimensions
!
!-----------------------------------------------------------------------
!
!  Set pio to return errors to this subroutine instead of handling them internally
!
   call pio_seterrorhandling(data_file%File, PIO_BCAST_ERROR)
   do n = 1,ndims
      dimid = 0

      !*** check to see whether dimension is already defined
      iostat = PIO_INQ_DIMID(data_file%File, name=trim(io_dims(n)%name),&
                              dimid=dimid)
      if (iostat /= PIO_NOERR) then ! dimension not yet defined
         iostat = PIO_DEF_DIM (data_file%File, name=trim(io_dims(n)%name), &
                                len=io_dims(n)%length, dimid=io_dims(n)%id)
      else
         io_dims(n)%id = dimid
      end if
   end do

!-----------------------------------------------------------------------
!
!  define the field
!
!-----------------------------------------------------------------------

   if (present(nftype)) then
      select case (trim(nftype))
        case ('float','FLOAT')
          xtype = PIO_REAL
        case ('double','DOUBLE')
          xtype = PIO_DOUBLE
        case ('integer','INTEGER')
          xtype = PIO_INT
        case ('char','CHAR','character', 'CHARACTER')
          xtype = PIO_CHAR
        case default
          call exit_POP(sigAbort,'unknown nftype') 
      end select
   else
    xtype = PIO_REAL
   endif

   !*** check to see whether field of this name already defined.

   iostat = PIO_INQ_VARID(data_file%File, trim(short_name), field_id)
   if (iostat /= PIO_NOERR) then ! variable was not yet defined
      ! Note currently must use vardesc to define var 
      iostat = PIO_DEF_VAR (data_file%File,name=trim(short_name), type=xtype,&
                             dimids=(/ (io_dims(n)%id, n=1,ndims) /),&
                             vardesc=vardesc)
      if (iostat /= pio_noerr) define_error = .true.

      iostat = PIO_INQ_VARID(data_file%File, trim(short_name), field_id)
      if (iostat /= pio_noerr) define_error = .true.
   end if

   if (define_error) then
       write(stdout,*) '(define_nstd_netcdf) Error for field = ', trim(short_name)
       call exit_POP(sigAbort, 'Error defining nonstandard CCSM netCDF field')
   endif

!-----------------------------------------------------------------------
!
!     Now define the field attributes
!
!-----------------------------------------------------------------------

    !*** long_name
    iostat = pio_inq_att(data_file%File, field_id, 'long_name', &
	                  xtype, len)  
    if (iostat /= PIO_NOERR) then ! attrib probably not defined
       iostat = pio_put_att(data_file%File, varid=field_id, &
                             name='long_name',    &
                             value=trim(long_name))
       if (iostat /= PIO_NOERR) define_error = .true.
    end if

    !*** units
    iostat = pio_inq_att(data_file%File, field_id, 'units', &
	                  xtype, len)  
    if (iostat /= PIO_NOERR) then ! attrib probably not defined
       iostat = pio_put_att(data_file%File, varid=field_id, &
                             name='units',        &
                             value=trim(units))
       if (iostat /= PIO_NOERR) define_error = .true.
    end if

    !*** coordinates
    if (present(coordinates)) then
       iostat = pio_inq_att(data_file%File, field_id, 'coordinates', &
                             xtype, len)  
       if (iostat /= PIO_NOERR) then ! attrib probably not defined
          iostat = pio_put_att(data_file%File, varid=field_id,    &
                                name='coordinates',     &
                                value=trim(coordinates))
          if (iostat /= PIO_NOERR) define_error = .true.
       end if
    endif

    !*** cell_methods
    if (present(method_string)) then
       iostat = pio_inq_att(data_file%File, field_id, 'cell_methods', &
                             xtype, len)  
       if (iostat /= PIO_NOERR) then ! attrib probably not defined
          iostat = pio_put_att(data_file%File, varid=field_id,    &
                                name='cell_methods',     &
                                value=trim(method_string))
          if (iostat /= PIO_NOERR) define_error = .true.
       end if
    endif

    !*** fill_value -- and missing_value, for now
    if (present(fill_value)) then
       iostat = pio_inq_att(data_file%File, varid=field_id, name='_FillValue', &
                             xtype=xtype, len=len)
       if (iostat /= PIO_NOERR) then ! attrib probably not defined
          iostat = pio_put_att(data_file%File, varid=field_id,        &
                                name='_FillValue',       &
                                value=fill_value)
          if (iostat /= PIO_NOERR) define_error = .true.
       end if
       iostat = pio_inq_att(data_file%File, varid=field_id, name='missing_value', &
                             xtype=xtype, len=len)
       if (iostat /= PIO_NOERR) then ! attrib probably not defined
          iostat = pio_put_att(data_file%File, varid=field_id,        &
                                name='missing_value',       &
                                value=fill_value)
          if (iostat /= PIO_NOERR) define_error = .true.
       end if
    endif


   if (define_error) call exit_POP(sigAbort, &
                     '(define_nstd_netcdf) Error adding attributes to field')

!
!  Reset PIO to handle errors internally
! 
   call pio_seterrorhandling(data_file%File, PIO_INTERNAL_ERROR)

!-----------------------------------------------------------------------
!EOC

 end subroutine define_nstd_netcdf

!***********************************************************************
!BOP
! !IROUTINE: write_time_bounds
! !INTERFACE:

 subroutine  write_time_bounds (data_file, time_bound_id, &
	                        time_bound_dims, time_bound_data)

! !INPUT PARAMETERS:
   integer (i4), intent (in)                ::  time_bound_id  
   type (io_dim), dimension(:), intent (in) ::  time_bound_dims
   real (r8), dimension(2,1),intent (in)    ::  time_bound_data

! !INPUT/OUTPUT PARAMETERS:
   type (datafile), target, intent (inout)  :: &
      data_file                         ! file to which field will be written

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (i4), dimension(2) ::  &
      start,length,count          ! dimension quantities for netCDF

   integer  :: &
      iostat,  &! netCDF status flag
      n         ! index

   integer  :: ncid, nout(5)

   logical (log_kind) :: &
      write_error         ! error flag

!-----------------------------------------------------------------------
!
!  exit define mode if necessary
!
!-----------------------------------------------------------------------

   write_error = .false.
   

   if (data_file%ldefine) then
      iostat = pio_enddef(data_file%File)
      data_file%ldefine = .false.
   endif

!-----------------------------------------------------------------------
!
!  make sure field has been defined
!
!-----------------------------------------------------------------------

   if (time_bound_id == 0) write_error = .true.

   if (write_error) then
      write(stdout,*) '(write_time_bounds) ERROR: undefined field -- time_bound'
      call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
      call exit_POP(sigAbort,' Attempt to write undefined time_bound in netCDF write')
   endif

!-----------------------------------------------------------------------
!
!  allocate dimension start,stop quantities
!
!-----------------------------------------------------------------------

   do n=1,2
      start (n) = time_bound_dims(n)%start
      length(n) = time_bound_dims(n)%stop - start(n) + 1
   end do
   
   iostat = pio_put_var(data_file%File, varid=time_bound_id, start=start(:), count=length(:), &
                        ival=time_bound_data) 

 end subroutine write_time_bounds

!***********************************************************************
!BOP
! !IROUTINE: write_nstd_netcdf
! !INTERFACE:

 subroutine write_nstd_netcdf(data_file,field_id,            &
                              ndims, io_dims,                &
                              nftype,                        &
                              lactive_time_dim,              &
                              indata_1d_r8,                  &
                              indata_2d_r8,                  &
                              indata_2d_r4,                  &
                              indata_3d_r4 ,                 &
                              indata_4d_r4,                  &
                              indata_1d_ch,                  &
                              indata_2d_ch                   )

! !DESCRIPTION:
!  This is a specialized, CCSM-speicific routine to write any desired
!  output field that cannot presently be defined through construct_io_field
!  to the CCSM version of the netCDF time-averaged history output files 
!
! !REVISION HISTORY:
!  same as module
!  USES
   use shr_pio_mod, only : shr_pio_getioroot
! !INPUT PARAMETERS:
   
   character (*), intent (in) ::  &
       nftype

   integer (i4), intent (in) :: &
      field_id                   ! netCDF id for the nonstandard variables

   integer (int_kind), intent (in) :: &
      ndims

   type (io_dim), dimension(:), intent (in) ::  &
      io_dims

   real (r8), dimension(:,:),intent (in) ::  &
      indata_2d_r8
   real (r8), dimension(:),  intent (in) ::  &
      indata_1d_r8
         
   real (r4), dimension(:,:,:,:), intent (in) ::  &
      indata_4d_r4
   real (r4), dimension(:,:,:),   intent (in) ::  &
      indata_3d_r4
   real (r4), dimension(:,:),     intent (in) ::  &
      indata_2d_r4

   character (*), dimension(:,:), intent (in) ::  &
      indata_2d_ch
   character (*), dimension(:),   intent (in) ::  &
      indata_1d_ch

   logical (log_kind), intent(in) ::  &
      lactive_time_dim


! !INPUT/OUTPUT PARAMETERS:

   type (datafile), target, intent (inout)  :: &
      data_file                         ! file to which field will be written

   optional ::           &
     indata_1d_r8,       &
     indata_2d_r8,       &
     indata_2d_r4,       &
     indata_3d_r4,       &
     indata_4d_r4,       &
     indata_1d_ch,       &
     indata_2d_ch

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer, parameter ::  &
      max_dims = 20

   integer , dimension(max_dims) :: &
      start,count               ! dimension quantities for netCDF

   integer  :: &
      iostat,  &! netCDF status flag
      n,m,     &! general indices
      tb        ! time indices

   integer  :: nout(5)

   logical (log_kind) :: &
      write_error,       &! error flag
      supported

   real (r4), allocatable, dimension (:,:,:,:,:)  ::  &
      outdata_5d_r4
   
   real (r4), allocatable, dimension (:,:,:,:)    ::  &
      outdata_4d_r4

   real (r4), allocatable, dimension (:,:,:)      ::  &
      outdata_3d_r4
   
   real (r4), allocatable, dimension (:,:)        ::  &
      outdata_2d_r4
   
   real (r8), allocatable, dimension (:)          ::  &
      outdata_1d_r8
   real (r8), allocatable, dimension (:,:)        ::  &
      outdata_2d_r8

   character(char_len), allocatable, dimension (:,:)  ::  &
      outdata_2d_ch

   character(1), dimension(char_len) :: &
      tmpString                          ! temp for manipulating output string

   integer :: ioroot

!-----------------------------------------------------------------------
!
!  exit define mode if necessary
!
!-----------------------------------------------------------------------


   write_error = .false.

   if (data_file%ldefine) then
      iostat = pio_enddef(data_file%File)
      data_file%ldefine = .false.
   endif
   
!-----------------------------------------------------------------------
!
!  make sure field has been defined
!
!-----------------------------------------------------------------------
   if (field_id == 0) write_error = .true.

   if (write_error) &
      call exit_POP(sigAbort, &
          '(write_nstd_netcdf) Attempt to write undefined field in netCDF write')

   supported = .true.
   
   ioroot = shr_pio_getioroot(inst_name)

!-----------------------------------------------------------------------
!
!  define start, count for all dimensions; do not allow out-of-bounds
!
!-----------------------------------------------------------------------
   if (ndims > max_dims) &
      call exit_POP(sigAbort, &
          '(write_nstd_netcdf) ndims > max_dims -- increase max_dims')

   do n=1,ndims
     start (n) = io_dims(n)%start
     count(n) = io_dims(n)%stop - start(n) + 1
   end do

   select case (trim(nftype))

      case('double','DOUBLE')
         select case (lactive_time_dim)
           case (.true.)
              select case (ndims)
                  case(2)
                     if (my_task == ioroot) then
                       nout(1) = size(indata_1d_r8,DIM=1)
                       allocate (outdata_2d_r8(nout(1),1))
                       outdata_2d_r8(:,1) = indata_1d_r8(:)
                     else
                       allocate (outdata_2d_r8(1,1))
                     endif
                     iostat = pio_put_var (data_file%File, field_id, outdata_2d_r8 )
                     deallocate (outdata_2d_r8)
                  case default
                   supported = .false. 
              end select ! ndims
           case (.false.)
              select case (ndims)
                  case(1)
                     iostat = pio_put_var (data_file%File, field_id, indata_1d_r8 )
                  case(2)
                     iostat = pio_put_var (data_file%File, field_id, indata_2d_r8 )
                  case default
                   supported = .false. 
              End select ! ndims
           end select ! lactive_time_dim

      case('float','FLOAT') 
         select case (lactive_time_dim)
           case (.true.)
              select case (ndims)
                  case(1)
                     supported = .false. 
                  case(2)
                     supported = .false. 
                  case(3)
                     if (my_task == ioroot) then
                       do n=1,ndims-1
                         nout(n) = size(indata_2d_r4,DIM=n)
                       enddo
                       tb = io_dims(ndims)%start
                       allocate (outdata_3d_r4(nout(1),nout(2),tb:tb))
                       outdata_3d_r4(:,:,tb) = indata_2d_r4(:,:)
                     else
                       allocate (outdata_3d_r4(1,1,1))
                     endif
                     iostat = pio_put_var (data_file%File, field_id, ival=outdata_3d_r4,  &
                                           start=start(:), count=count(:))
                     deallocate (outdata_3d_r4)
                  case(4)
                     if (my_task == ioroot) then
                       do n=1,ndims-1
                         nout(n) = size(indata_3d_r4,DIM=n)
                       enddo
                       tb = io_dims(ndims)%start
                       allocate (outdata_4d_r4(nout(1),nout(2),nout(3),tb:tb))
                       outdata_4d_r4(:,:,:,tb) = indata_3d_r4(:,:,:)
                     else
                       allocate (outdata_4d_r4(1,1,1,1))
                     endif
                     iostat = pio_put_var (data_file%File, field_id, ival=outdata_4d_r4,  &
                                           start=start(:), count=count(:))
                     deallocate (outdata_4d_r4)
                  case(5)
                     if (my_task == ioroot) then
                       do n=1,ndims-1
                         nout(n) = size(indata_4d_r4,DIM=n)
                       enddo
                       tb = io_dims(ndims)%start
                       allocate (outdata_5d_r4(nout(1),nout(2),nout(3),nout(4),tb:tb))
                       outdata_5d_r4(:,:,:,:,tb) = indata_4d_r4(:,:,:,:)
                     else
                       allocate (outdata_5d_r4(1,1,1,1,1))
                     endif
                     iostat = pio_put_var (data_file%File, field_id, ival=outdata_5d_r4,  &
                                           start=start(:), count=count(:))
                     deallocate (outdata_5d_r4)
                  case default
                   supported = .false. 
              end select ! ndims
           case (.false.)
              select case (ndims)
                  case(1)
                     supported = .false. 
                  case(2)
                     iostat = pio_put_var (data_file%File, field_id, indata_2d_r4 )
                  case(3)
                     iostat = pio_put_var (data_file%File, field_id, indata_3d_r4 )
                  case(4)
                     iostat = pio_put_var (data_file%File, field_id, indata_4d_r4 )
                  case default
                   supported = .false. 
              end select ! ndims
         end select ! lactive_time_dim

      case('char','character','CHAR','CHARACTER') 
         select case (lactive_time_dim)
           case (.true.)
              select case (ndims)
                  case default
                   supported = .false. 
              end select ! ndims
           case (.false.)
              select case (ndims)
                  case(2)
                     do n=1,io_dims(2)%length
                       start(1) = 1
                       start(2) = n
                       count(1) = len_trim(indata_1d_ch(n))
                       count(2) = 1
                       do m = 1,count(1)
                          tmpString(m:m) = indata_1d_ch(n)(m:m)
                       end do
                       iostat = pio_put_var (data_file%File, field_id, &
                               ival=tmpString(1:count(1)), start=start, count=count)
                     enddo

                  case default
                   supported = .false. 
              end select ! ndims
         end select ! lactive_time_dim

      case default

    end select ! nftype

!-----------------------------------------------------------------------
!
!  check for errors
!
!-----------------------------------------------------------------------

   if (.not. supported) then
     call document('write_nstd_netcdf', 'ndims', ndims)
     call document('write_nstd_netcdf', 'nftype', trim(nftype))
     call document('write_nstd_netcdf', 'lactive_time_dim', lactive_time_dim)
     call exit_POP(sigAbort, '(write_nstd_netcdf) option not supported')
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine write_nstd_netcdf


!***********************************************************************
!BOP
! !IROUTINE: check_definemode 
! !INTERFACE:

 subroutine check_definemode (data_file, name)

! !DESCRIPTION:
!  This utility routine checks if the data file is in define mode
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

   type (datafile), target, intent (inout)  :: &
      data_file                        

   character(*),intent (in):: name

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer  :: &
      iostat        ! netCDF status flag

   logical (log_kind) :: &
      write_error         ! error flag

   character (char_len) :: string


!-----------------------------------------------------------------------
!
!  make sure file is in define mode
!
!-----------------------------------------------------------------------


   if (.not. data_file%ldefine) &
     call exit_POP(sigAbort, &
                   '('//trim(name)//') attempt to define field but not in define mode')

!-----------------------------------------------------------------------
!EOC

 end subroutine check_definemode


!***********************************************************************
end module io_netcdf

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
