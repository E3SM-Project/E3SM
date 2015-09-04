!=======================================================================
!BOP
!
! !MODULE: ice_read_write
!
! !DESCRIPTION:
!
! Routines for opening, reading and writing external files
!
! !REVISION HISTORY:
!  SVN:$Id$
!
! author: Tony Craig, NCAR
!
! 2004: Block structure added by William Lipscomb, LANL
! 2006: Converted to free source form (F90) by Elizabeth Hunke
! 2007: netcdf versions added by Alison McLaren & Ann Keen, Met Office
!
! !INTERFACE:
!
      module ice_read_write
!
! !USES:
!
      use ice_kinds_mod
      use ice_constants
      use ice_communicate, only: my_task, master_task
      use ice_broadcast
      use ice_domain_size
      use ice_blocks
      use ice_fileunits
#ifdef ncdf
      use netcdf      
#endif
!
!EOP

      implicit none

       public :: ice_read_global_nc

       interface ice_read_global_nc
          module procedure ice_read_global_nc_dbl, &
                           ice_read_global_nc_r4
       end interface

!=======================================================================

      contains

!=======================================================================
!
!BOP
!
! !IROUTINE: ice_open - opens an unformatted file for reading
!
! !INTERFACE:
!
      subroutine ice_open(nu, filename, nbits)
!
! !DESCRIPTION:
!
! Opens an unformatted file for reading \\
! nbits indicates whether the file is sequential or direct access
!
! !REVISION HISTORY:
!
! author: Tony Craig, NCAR
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
           nu        , & ! unit number
           nbits         ! no. of bits per variable (0 for sequential access)

      character (*) :: filename
!
!EOP
!
      if (my_task == master_task) then

         if (nbits == 0) then   ! sequential access

            open(nu,file=filename,form='unformatted')

         else                   ! direct access
            open(nu,file=filename,recl=nx_global*ny_global*nbits/8, &
                  form='unformatted',access='direct')
         endif                   ! nbits = 0

      endif                      ! my_task = master_task

      end subroutine ice_open

!=======================================================================
!BOP
!
! !IROUTINE: ice_read - read and scatter an unformatted file
!
! !INTERFACE:
!
      subroutine ice_read(nu,  nrec,  work, atype, diag, &
                          field_loc, field_type, &
                          ignore_eof, hit_eof)
!
! !DESCRIPTION:
!
! Read an unformatted file and scatter to processors\\
! work is a real array, atype indicates the format of the data\\
! If the optional variables field_loc and field_type are present \\
! the ghost cells are filled using values from the global array.\\
! This prevents them from being filled with zeroes in land cells \\
! (subroutine ice_HaloUpdate need not be called).
!
! !REVISION HISTORY:
!
! author: Tony Craig, NCAR
!
! !USES:
!
      use ice_domain
      use ice_gather_scatter
      use ice_work, only: work_g1, work_gr, work_gi4, work_gi8
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
           nu            , & ! unit number
           nrec              ! record number (0 for sequential access)

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks), &
           intent(out) :: &
           work              ! output array (real, 8-byte)

      character (len=4), intent(in) :: &
           atype             ! format for input array
                             ! (real/integer, 4-byte/8-byte)

      logical (kind=log_kind), intent(in) :: &
           diag              ! if true, write diagnostic output

      integer (kind=int_kind), optional, intent(in) :: &
           field_loc, &      ! location of field on staggered grid
           field_type        ! type of field (scalar, vector, angle)

      logical (kind=log_kind), optional, intent(in)  :: ignore_eof
      logical (kind=log_kind), optional, intent(out) :: hit_eof
!
!EOP
!
      integer (kind=int_kind) :: i, j, ios

      real (kind=dbl_kind) :: &
         amin, amax         ! min and max values of input array

      logical (kind=log_kind) :: ignore_eof_use

      if (my_task == master_task) then
         allocate(work_g1(nx_global,ny_global))
      else
         allocate(work_g1(1,1))   ! to save memory
      endif

      if (my_task == master_task) then

    !-------------------------------------------------------------------
    ! Read global array according to format atype
    !-------------------------------------------------------------------
         if (present(hit_eof)) hit_eof = .false.

         if (atype == 'ida4') then
            allocate(work_gi4(nx_global,ny_global))
            read(nu,rec=nrec) work_gi4
            work_g1 = real(work_gi4,kind=dbl_kind)
            deallocate(work_gi4)
         elseif (atype == 'ida8') then
            allocate(work_gi8(nx_global,ny_global))
            read(nu,rec=nrec) work_gi8
            work_g1 = real(work_gi8,kind=dbl_kind)
            deallocate(work_gi8)
         elseif (atype == 'rda4') then
            allocate(work_gr(nx_global,ny_global))
            read(nu,rec=nrec) work_gr
            work_g1 = work_gr
            deallocate(work_gr)
         elseif (atype == 'rda8') then
            read(nu,rec=nrec) work_g1
         elseif (atype == 'ruf8') then
            if (present(ignore_eof)) then
               ignore_eof_use = ignore_eof
            else
               ignore_eof_use = .false.
            endif
            if (ignore_eof_use) then
             ! Read line from file, checking for end-of-file
               read(nu, iostat=ios) ((work_g1(i,j),i=1,nx_global), &
                                                   j=1,ny_global)
               if (present(hit_eof)) hit_eof = ios < 0
            else
               read(nu) ((work_g1(i,j),i=1,nx_global),j=1,ny_global)
            endif
         else
            write(nu_diag,*) ' ERROR: reading unknown atype ',atype
         endif
      endif                     ! my_task = master_task

      if (present(hit_eof)) then
         call broadcast_scalar(hit_eof,master_task)
         if (hit_eof) then
            deallocate(work_g1)
            return
         endif
      endif


    !-------------------------------------------------------------------
    ! optional diagnostics
    !-------------------------------------------------------------------
      if (my_task==master_task .and. diag) then
         amin = minval(work_g1)
         amax = maxval(work_g1, mask = work_g1 /= spval_dbl)
         write(nu_diag,*) ' read_global ',nu, nrec, amin, amax
      endif

    !-------------------------------------------------------------------
    ! Scatter data to individual processors.
    ! NOTE: Ghost cells are not updated unless field_loc is present.
    !-------------------------------------------------------------------

      if (present(field_loc)) then
         call scatter_global(work, work_g1, master_task, distrb_info, &
                             field_loc, field_type)
      else
         call scatter_global(work, work_g1, master_task, distrb_info, &
                             field_loc_noupdate, field_type_noupdate)
      endif

      deallocate(work_g1)

      end subroutine ice_read

!=======================================================================
!BOP
!
! !IROUTINE: ice_read_global - read an unformatted file
!
! !INTERFACE:
!
      subroutine ice_read_global (nu,  nrec,  work_g, atype, diag, &
                                  ignore_eof, hit_eof)
!
! !DESCRIPTION:
!
! Read an unformatted file \\
! Just like ice_read except that it returns a global array \\
! work_g is a real array, atype indicates the format of the data
!
! !REVISION HISTORY:
! Adapted by William Lipscomb, LANL, from ice_read
!
! !USES:
!
      use ice_work, only: work_gr, work_gi4, work_gi8
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
           nu            , & ! unit number
           nrec              ! record number (0 for sequential access)

      real (kind=dbl_kind), dimension(:,:), &
           intent(out) :: &
           work_g            ! output array (real, 8-byte)

      character (len=4) :: &
           atype             ! format for input array
                             ! (real/integer, 4-byte/8-byte)

      logical (kind=log_kind) :: &
           diag              ! if true, write diagnostic output

      logical (kind=log_kind), optional, intent(in)  :: ignore_eof
      logical (kind=log_kind), optional, intent(out) :: hit_eof
!
!EOP
!
      integer (kind=int_kind) :: i, j, ios

      real (kind=dbl_kind) :: &
         amin, amax         ! min and max values of input array

      logical (kind=log_kind) :: ignore_eof_use

      work_g(:,:) = c0

      if (my_task == master_task) then

    !-------------------------------------------------------------------
    ! Read global array according to format atype
    !-------------------------------------------------------------------
         if (present(hit_eof)) hit_eof = .false.

         if (atype == 'ida4') then
            allocate(work_gi4(nx_global,ny_global))
            read(nu,rec=nrec) work_gi4
            work_g = real(work_gi4,kind=dbl_kind)
            deallocate(work_gi4)
         elseif (atype == 'ida8') then
            allocate(work_gi8(nx_global,ny_global))
            read(nu,rec=nrec) work_gi8
            work_g = real(work_gi8,kind=dbl_kind)
            deallocate(work_gi8)
         elseif (atype == 'rda4') then
            allocate(work_gr(nx_global,ny_global))
            read(nu,rec=nrec) work_gr
            work_g = work_gr
            deallocate(work_gr)
         elseif (atype == 'rda8') then
            read(nu,rec=nrec) work_g
         elseif (atype == 'ruf8') then
            if (present(ignore_eof)) then
               ignore_eof_use = ignore_eof
            else
               ignore_eof_use = .false.
            endif
            if (ignore_eof_use) then
               ! Read line from file, checking for end-of-file
               read(nu, iostat=ios) ((work_g(i,j),i=1,nx_global), &
                                                  j=1,ny_global)
               if (present(hit_eof)) hit_eof = ios < 0
            else
               read(nu) ((work_g(i,j),i=1,nx_global),j=1,ny_global)
            endif
         else
            write(nu_diag,*) ' ERROR: reading unknown atype ',atype
         endif
      endif                     ! my_task = master_task

      if (present(hit_eof)) then
         call broadcast_scalar(hit_eof,master_task)
         if (hit_eof) return
      endif

    !-------------------------------------------------------------------
    ! optional diagnostics
    !-------------------------------------------------------------------
      if (my_task == master_task .and. diag) then
         amin = minval(work_g)
         amax = maxval(work_g, mask = work_g /= spval_dbl)
         write(nu_diag,*) ' read_global ',nu, nrec, amin, amax
      endif

      end subroutine ice_read_global

!=======================================================================
!BOP
!
! !IROUTINE: ice_write - writes an unformatted file
!
! !INTERFACE:
!
      subroutine ice_write(nu, nrec, work, atype, diag)
!
! !DESCRIPTION:
!
! Writes an unformatted file \\
! work is a real array, atype indicates the format of the data
!
! !REVISION HISTORY:
!
! author: Tony Craig, NCAR
!
! !USES:
!
      use ice_gather_scatter
      use ice_domain
      use ice_work, only: work_g1, work_gr, work_gi4, work_gi8
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
           nu            , & ! unit number
           nrec              ! record number (0 for sequential access)

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks), &
           intent(in) :: &
           work              ! input array (real, 8-byte)

      character (len=4) :: &
           atype             ! format for output array
                             ! (real/integer, 4-byte/8-byte)

      logical (kind=log_kind) :: &
           diag              ! if true, write diagnostic output
!
!EOP
!
      integer (kind=int_kind) :: i, j

      real (kind=dbl_kind) :: &
         amin, amax     ! min and max values of ouput array

    !-------------------------------------------------------------------
    ! Gather data from individual processors
    !-------------------------------------------------------------------

      if (my_task == master_task) then
         allocate(work_g1(nx_global,ny_global))
      else
         allocate(work_g1(1,1)) ! to save memory
      endif

      call gather_global(work_g1, work, master_task, distrb_info)

      if (my_task == master_task) then

    !-------------------------------------------------------------------
    ! Write global array according to format atype
    !-------------------------------------------------------------------
         if (atype == 'ida4') then
            allocate(work_gi4(nx_global,ny_global))
            work_gi4 = nint(work_g1)
            write(nu,rec=nrec) work_gi4
            deallocate(work_gi4)
         elseif (atype == 'ida8') then
            allocate(work_gi8(nx_global,ny_global))
            work_gi8 = nint(work_g1)
            write(nu,rec=nrec) work_gi8           
            deallocate(work_gi8)
         elseif (atype == 'rda4') then
            allocate(work_gr(nx_global,ny_global))
            work_gr = work_g1
            write(nu,rec=nrec) work_gr
            deallocate(work_gr)
         elseif (atype == 'rda8') then
            write(nu,rec=nrec) work_g1
         elseif (atype == 'ruf8') then
            write(nu) ((work_g1(i,j),i=1,nx_global),j=1,ny_global)
         else
            write(nu_diag,*) ' ERROR: writing unknown atype ',atype
         endif

    !-------------------------------------------------------------------
    ! diagnostics
    !-------------------------------------------------------------------
         if (diag) then
            amin = minval(work_g1)
            amax = maxval(work_g1, mask = work_g1 /= spval_dbl)
            write(nu_diag,*) ' write_global ', nu, nrec, amin, amax
         endif

      endif                     ! my_task = master_task

      deallocate(work_g1)

      end subroutine ice_write

!=======================================================================
!BOP
!
! !IROUTINE: ice_write_nc - writes a field to a netcdf file
!
! !INTERFACE:
!
      subroutine ice_write_nc(fid, nrec, varname, work, atype, diag)
!
! !DESCRIPTION:
!
! Writes a field to a netcdf file \\
! work is a real array, atype indicates the format of the data
!
! !REVISION HISTORY:
!
! author: David A Bailey, NCAR
!
! !USES:
!
      use ice_gather_scatter
      use ice_domain
      use ice_work, only: work_g1, work_gr, work_gi4, work_gi8
      use ice_exit
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
           fid          ,&   ! netcdf file id
           nrec              ! record number

      character (len=*), intent(in) :: varname

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks), &
           intent(in) :: &
           work              ! input array (real, 8-byte)

      character (len=4), intent(in) :: &
           atype             ! format for output array
                             ! (real/integer, 4-byte/8-byte)

      logical (kind=log_kind), intent(in) :: &
           diag              ! if true, write diagnostic output
!
!EOP
!
      integer (kind=int_kind) :: i, j, varid, numDims

      integer (kind=int_kind) :: &
         status        ! status variable from netCDF routine 

      real (kind=dbl_kind) :: &
         amin, amax     ! min and max values of ouput array

      integer (kind=int_kind), allocatable :: &
         start_arr(:), count_arr(:)

    !-------------------------------------------------------------------
    ! Gather data from individual processors
    !-------------------------------------------------------------------

      if (my_task == master_task) then
         allocate(work_g1(nx_global,ny_global))
      else
         allocate(work_g1(1,1)) ! to save memory
      endif

      call gather_global(work_g1, work, master_task, distrb_info)

      if (my_task == master_task) then

         status = nf90_inq_varid(fid, trim(varname), varid)

         if (status /= nf90_noerr) then
           call abort_ice ( &
               'ice_write_nc: Cannot find variable '//trim(varname) )
         endif

         status = nf90_inquire_variable(fid, varid, ndims = numDims)

         if (status /= nf90_noerr) then
           call abort_ice ( &
               'ice_write_nc: Cannot find dimensions for '//trim(varname) )
         endif

         allocate(start_arr(numDims))
         allocate(count_arr(numDims))

         if (numDims > 2) then
            start_arr(1) = 1
            start_arr(2) = 1
            start_arr(3) = nrec
            count_arr(1) = nx_global
            count_arr(2) = ny_global
            count_arr(3) = 1
         else
            start_arr(1) = 1
            start_arr(2) = 1
            count_arr(1) = nx_global
            count_arr(2) = ny_global
         endif

    !-------------------------------------------------------------------
    ! Write global array according to format atype
    !-------------------------------------------------------------------
         if (atype == 'ida4') then
            allocate(work_gi4(nx_global,ny_global))
            work_gi4 = nint(work_g1)
            status = nf90_put_var(fid,varid,work_gi4, &
                                  start=start_arr,  &
                                  count=count_arr)
            deallocate(work_gi4)
         elseif (atype == 'ida8') then
            allocate(work_gi8(nx_global,ny_global))
            work_gi8 = nint(work_g1)
            status = nf90_put_var(fid,varid,work_gi8, &
                                  start=start_arr,  &
                                  count=count_arr)
            deallocate(work_gi8)
         elseif (atype == 'rda4') then
            allocate(work_gr(nx_global,ny_global))
            work_gr = work_g1
            status = nf90_put_var(fid,varid,work_gr, &
                                  start=start_arr,  &
                                  count=count_arr)
            deallocate(work_gr)
         elseif (atype == 'rda8') then
            status = nf90_put_var(fid,varid,work_g1, &
                                  start=start_arr,  &
                                  count=count_arr)
         else
            write(nu_diag,*) ' ERROR: writing unknown atype ',atype
         endif

    !-------------------------------------------------------------------
    ! diagnostics
    !-------------------------------------------------------------------
         if (diag) then
            amin = minval(work_g1)
            amax = maxval(work_g1, mask = work_g1 /= spval_dbl)
            write(nu_diag,*) ' write_global ', fid, varid, nrec, amin, amax
         endif

         deallocate(start_arr)
         deallocate(count_arr)

      endif                     ! my_task = master_task

      deallocate(work_g1)

      end subroutine ice_write_nc
      
!=======================================================================
!
!BOP
!
! !IROUTINE: ice_open_nc - opens a netCDF file for reading
!
! !INTERFACE:
!
      subroutine ice_open_nc(filename, fid)
!
! !DESCRIPTION:
!
! Opens a netCDF file for reading
!
! !REVISION HISTORY:
!
! Adapted by Alison McLaren, Met Office from ice_open
!
! !USES:
 
      use ice_exit
!
! !INPUT/OUTPUT PARAMETERS:
!

      character (char_len_long), intent(in) :: & 
           filename      ! netCDF filename

      integer (kind=int_kind), intent(out) :: &
           fid           ! unit number
!
!EOP
!
#ifdef ncdf
      integer (kind=int_kind) :: &
        status        ! status variable from netCDF routine 


      if (my_task == master_task) then

          status = nf90_open(filename, NF90_NOWRITE, fid)
          if (status /= nf90_noerr) then
             call abort_ice ( & 
                   'ice_open_nc: Cannot open '//trim(filename) )
          endif

      endif                      ! my_task = master_task

#else
      fid = -999 ! to satisfy intent(out) attribute
#endif
      end subroutine ice_open_nc

!=======================================================================
!BOP
!
! !IROUTINE: ice_read_nc - read and scatter one field from a netCDF file
!
! !INTERFACE:
!
      subroutine ice_read_nc(fid,  nrec,  varname, work,  diag, &
                             field_loc, field_type)
!
! !DESCRIPTION:
!
! Read a netCDF file and scatter to processors\\
! If the optional variables field_loc and field_type are present \\
! the ghost cells are filled using values from the global array.\\
! This prevents them from being filled with zeroes in land cells \\
! (subroutine ice_HaloUpdate need not be called).
!
! !REVISION HISTORY:
!
! Adapted by Alison McLaren, Met Office from ice_read
!
! !USES:
!
      use ice_domain
      use ice_gather_scatter
#ifdef ORCA_GRID
      use ice_work, only: work_g1, work_g2
#else
      use ice_work, only: work_g1
#endif
      use ice_exit
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
           fid           , & ! file id
           nrec              ! record number 

      logical (kind=log_kind), intent(in) :: &
           diag              ! if true, write diagnostic output

      character (len=*), intent(in) :: & 
           varname           ! field name in netcdf file

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks), &
           intent(out) :: &
           work              ! output array (real, 8-byte)

      integer (kind=int_kind), optional, intent(in) :: &
           field_loc, &      ! location of field on staggered grid
           field_type        ! type of field (scalar, vector, angle)
!
!EOP
!
#ifdef ncdf
! netCDF file diagnostics:
      integer (kind=int_kind) :: & 
         varid,           & ! netcdf id for field
         status,          & ! status output from netcdf routines
         ndim, nvar,      & ! sizes of netcdf file
         id,              & ! dimension index
         dimlen             ! size of dimension

      real (kind=dbl_kind) :: &
         amin, amax         ! min and max values of input array

      character (char_len) :: &
         dimname            ! dimension name            
!
      if (my_task == master_task) then
         allocate(work_g1(nx_global,ny_global))
      else
         allocate(work_g1(1,1))   ! to save memory
      endif

#ifdef ORCA_GRID
      if (my_task == master_task) then
         allocate(work_g2(nx_global+2,ny_global+1))
      else
         allocate(work_g2(1,1))   ! to save memory
      endif
#endif

      if (my_task == master_task) then

        !-------------------------------------------------------------
        ! Find out ID of required variable
        !-------------------------------------------------------------

         status = nf90_inq_varid(fid, trim(varname), varid)
 
         if (status /= nf90_noerr) then
           call abort_ice ( & 
               'ice_read_nc: Cannot find variable '//trim(varname) )
         endif

       !--------------------------------------------------------------
       ! Read global array 
       !--------------------------------------------------------------

#ifndef ORCA_GRID
         status = nf90_get_var( fid, varid, work_g1, &
               start=(/1,1,nrec/), & 
               count=(/nx_global,ny_global,1/) )
#else
         status = nf90_get_var( fid, varid, work_g2, &
               start=(/1,1,nrec/), &
               count=(/nx_global+2,ny_global+1,1/) )
         work_g1=work_g2(2:nx_global+1,1:ny_global)
#endif

      endif                     ! my_task = master_task

    !-------------------------------------------------------------------
    ! optional diagnostics
    !-------------------------------------------------------------------

      if (my_task==master_task .and. diag) then

!        write(nu_diag,*) & 
!          'ice_read_nc, fid= ',fid, ', nrec = ',nrec, & 
!          ', varname = ',trim(varname)
         status = nf90_inquire(fid, nDimensions=ndim, nVariables=nvar)
!        write(nu_diag,*) 'ndim= ',ndim,', nvar= ',nvar
         do id=1,ndim
           status = nf90_inquire_dimension(fid,id,name=dimname,len=dimlen)
!          write(nu_diag,*) 'Dim name = ',trim(dimname),', size = ',dimlen
         enddo
         amin = minval(work_g1)
         amax = maxval(work_g1, mask = work_g1 /= spval_dbl)
         write(nu_diag,*) ' read_global ',fid, varid, nrec, amin, amax

      endif

    !-------------------------------------------------------------------
    ! Scatter data to individual processors.
    ! NOTE: Ghost cells are not updated unless field_loc is present.
    !-------------------------------------------------------------------

      if (present(field_loc)) then
         call scatter_global(work, work_g1, master_task, distrb_info, &
                             field_loc, field_type)
      else
         call scatter_global(work, work_g1, master_task, distrb_info, &
                             field_loc_noupdate, field_type_noupdate)
      endif

      deallocate(work_g1)
#ifdef ORCA_GRID
      deallocate(work_g2)
#endif

#else
      work = c0 ! to satisfy intent(out) attribute
#endif
      end subroutine ice_read_nc
!
!=======================================================================
!BOP
!
! !IROUTINE: ice_read_global_nc - read one field from a netcdf file
!
! !INTERFACE:
!
      subroutine ice_read_global_nc_dbl (fid,  nrec, varname, work_g, diag)
!
! !DESCRIPTION:
!
! Read a netcdf file \\
! Just like ice_read_nc except that it returns a global array \\
! work_g is a real array
!
! !REVISION HISTORY:
! Adapted by William Lipscomb, LANL, from ice_read
! Adapted by Ann Keen, Met Office, to read from a netcdf file 
!
! !USES:
! 
      use ice_exit
#ifdef ORCA_GRID
      use ice_work, only: work_g3
#endif
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
           fid           , & ! file id
           nrec              ! record number 

     character (len=*), intent(in) :: & 
           varname           ! field name in netcdf file        

      real (kind=dbl_kind), dimension(:,:), &
           intent(out) :: &
           work_g            ! output array (real, 8-byte)

      logical (kind=log_kind) :: &
           diag              ! if true, write diagnostic output
!
!EOP
!
#ifdef ncdf
! netCDF file diagnostics:
      integer (kind=int_kind) :: & 
         varid,           & ! netcdf id for field
         status,          & ! status output from netcdf routines
         ndim, nvar,      & ! sizes of netcdf file
         id,              & ! dimension index
         dimlen             ! size of dimension      

      real (kind=dbl_kind) :: &
         amin, amax         ! min and max values of input array

     character (char_len) :: &
         dimname            ! dimension name            

!
#ifdef ORCA_GRID
      if (my_task == master_task) then
          allocate(work_g3(nx_global+2,ny_global+1))
       else
          allocate(work_g3(1,1))   ! to save memory
       endif

      work_g3(:,:) = c0
#endif
      work_g(:,:) = c0

      if (my_task == master_task) then

        !-------------------------------------------------------------
        ! Find out ID of required variable
        !-------------------------------------------------------------

         status = nf90_inq_varid(fid, trim(varname), varid)

         if (status /= nf90_noerr) then
           call abort_ice ( & 
            'ice_read_global_nc: Cannot find variable '//trim(varname) )
         endif

       !--------------------------------------------------------------
       ! Read global array 
       !--------------------------------------------------------------
 
#ifndef ORCA_GRID
         status = nf90_get_var( fid, varid, work_g, &
               start=(/1,1,nrec/), & 
               count=(/nx_global,ny_global,1/) )
#else
         status = nf90_get_var( fid, varid, work_g3, &
               start=(/1,1,nrec/), &
               count=(/nx_global+2,ny_global+1,1/) )
         work_g=work_g3(2:nx_global+1,1:ny_global)
#endif

      endif                     ! my_task = master_task

    !-------------------------------------------------------------------
    ! optional diagnostics
    !-------------------------------------------------------------------

      if (my_task == master_task .and. diag) then

!         write(nu_diag,*) & 
!           'ice_read_global_nc, fid= ',fid, ', nrec = ',nrec, & 
!           ', varname = ',trim(varname)
          status = nf90_inquire(fid, nDimensions=ndim, nVariables=nvar)
!         write(nu_diag,*) 'ndim= ',ndim,', nvar= ',nvar
          do id=1,ndim
            status = nf90_inquire_dimension(fid,id,name=dimname,len=dimlen)
!           write(nu_diag,*) 'Dim name = ',trim(dimname),', size = ',dimlen
         enddo
         amin = minval(work_g)
         amax = maxval(work_g, mask = work_g /= spval_dbl)
!        write(nu_diag,*) 'min and max = ', amin, amax
!        write(nu_diag,*) ''
         write(nu_diag,*) ' read_global ',fid, varid, nrec, amin, amax

      endif

#ifdef ORCA_GRID
      deallocate(work_g3)
#endif

#else
      work_g = c0 ! to satisfy intent(out) attribute
#endif
      end subroutine ice_read_global_nc_dbl
!=======================================================================
!BOP
!
! !IROUTINE: ice_read_global_nc - read one field from a netcdf file
!
! !INTERFACE:
!
      subroutine ice_read_global_nc_r4 (fid,  nrec, varname, work_g, diag)
!
! !DESCRIPTION:
!
! Read a netcdf file \\
! Just like ice_read_nc except that it returns a global array \\
! work_g is a real array
!
! !REVISION HISTORY:
! Adapted by William Lipscomb, LANL, from ice_read
! Adapted by Ann Keen, Met Office, to read from a netcdf file 
!
! !USES:
! 
      use ice_exit
#ifdef ORCA_GRID
      use ice_work, only: work_g3
#endif
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
           fid           , & ! file id
           nrec              ! record number 

     character (len=*), intent(in) :: & 
           varname           ! field name in netcdf file        

      real (kind=real_kind), dimension(:,:), &
           intent(out) :: &
           work_g            ! output array (real, 8-byte)

      logical (kind=log_kind) :: &
           diag              ! if true, write diagnostic output
!
!EOP
!
#ifdef ncdf
! netCDF file diagnostics:
      integer (kind=int_kind) :: & 
         varid,           & ! netcdf id for field
         status,          & ! status output from netcdf routines
         ndim, nvar,      & ! sizes of netcdf file
         id,              & ! dimension index
         dimlen             ! size of dimension      

      real (kind=dbl_kind) :: &
         amin, amax         ! min and max values of input array

     character (char_len) :: &
         dimname            ! dimension name            

!
#ifdef ORCA_GRID
      if (my_task == master_task) then
          allocate(work_g3(nx_global+2,ny_global+1))
       else
          allocate(work_g3(1,1))   ! to save memory
       endif

      work_g3(:,:) = c0
#endif
      work_g(:,:) = c0

      if (my_task == master_task) then

        !-------------------------------------------------------------
        ! Find out ID of required variable
        !-------------------------------------------------------------

         status = nf90_inq_varid(fid, trim(varname), varid)

         if (status /= nf90_noerr) then
           call abort_ice ( & 
            'ice_read_global_nc: Cannot find variable '//trim(varname) )
         endif

       !--------------------------------------------------------------
       ! Read global array 
       !--------------------------------------------------------------
 
#ifndef ORCA_GRID
         status = nf90_get_var( fid, varid, work_g, &
               start=(/1,1,nrec/), & 
               count=(/nx_global,ny_global,1/) )
#else
         status = nf90_get_var( fid, varid, work_g3, &
               start=(/1,1,nrec/), &
               count=(/nx_global+2,ny_global+1,1/) )
         work_g=work_g3(2:nx_global+1,1:ny_global)
#endif

      endif                     ! my_task = master_task

    !-------------------------------------------------------------------
    ! optional diagnostics
    !-------------------------------------------------------------------

      if (my_task == master_task .and. diag) then

!         write(nu_diag,*) & 
!           'ice_read_global_nc, fid= ',fid, ', nrec = ',nrec, & 
!           ', varname = ',trim(varname)
          status = nf90_inquire(fid, nDimensions=ndim, nVariables=nvar)
!         write(nu_diag,*) 'ndim= ',ndim,', nvar= ',nvar
          do id=1,ndim
            status = nf90_inquire_dimension(fid,id,name=dimname,len=dimlen)
!           write(nu_diag,*) 'Dim name = ',trim(dimname),', size = ',dimlen
         enddo
         amin = minval(work_g)
         amax = maxval(work_g, mask = work_g /= spval_dbl)
!        write(nu_diag,*) 'min and max = ', amin, amax
!        write(nu_diag,*) ''
         write(nu_diag,*) ' read_global ',fid, varid, nrec, amin, amax

      endif

#ifdef ORCA_GRID
      deallocate(work_g3)
#endif

#else
      work_g = c0 ! to satisfy intent(out) attribute
#endif
      end subroutine ice_read_global_nc_r4

!=======================================================================
!BOP
!
! !IROUTINE: ice_close_nc - closes a netCDF file
!
! !INTERFACE:
!
      subroutine ice_close_nc(fid)
!
! !DESCRIPTION:
!
! Closes a netCDF file
!
! !REVISION HISTORY:
!
! author: Alison McLaren, Met Office
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
           fid           ! unit number
!
!EOP
!
#ifdef ncdf
      integer (kind=int_kind) :: &
        status        ! status variable from netCDF routine 

      if (my_task == master_task) then

         status = nf90_close(fid)

      endif

#endif
      end subroutine ice_close_nc

!=======================================================================

      end module ice_read_write

!=======================================================================
