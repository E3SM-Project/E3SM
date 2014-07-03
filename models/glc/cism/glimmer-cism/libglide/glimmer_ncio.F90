!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glimmer_ncio.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

#define NCO outfile%nc
#define NCI infile%nc

module glimmer_ncio
  !*FD module for common netCDF I/O
  !*FD written by Magnus Hagdorn, 2004

  use glimmer_ncdf

  implicit none

  integer,parameter,private :: msglen=512
  
contains
  !*****************************************************************************
  ! netCDF output
  !*****************************************************************************  
  subroutine openall_out(model,outfiles)
    !*FD open all netCDF files for output
    use glide_types
    use glimmer_ncdf
    implicit none
    type(glide_global_type) :: model
    type(glimmer_nc_output),pointer,optional :: outfiles
    
    ! local variables
    type(glimmer_nc_output), pointer :: oc

    if (present(outfiles)) then
       oc => outfiles
    else
       oc=>model%funits%out_first
    end if

    do while(associated(oc))
       if (oc%append) then
          call glimmer_nc_openappend(oc,model)
       else
          call glimmer_nc_createfile(oc,model)
       end if
       oc=>oc%next
    end do
  end subroutine openall_out

  subroutine closeall_out(model,outfiles)
    !*FD close all netCDF files for output
    use glide_types
    use glimmer_ncdf
    implicit none
    type(glide_global_type) :: model
    type(glimmer_nc_output),pointer,optional :: outfiles

    ! local variables
    type(glimmer_nc_output), pointer :: oc

    if (present(outfiles)) then
       oc => outfiles
    else
       oc=>model%funits%out_first
    end if

    do while(associated(oc))
       oc=>delete(oc)
    end do
    if (.not.present(outfiles)) model%funits%out_first=>NULL()
  end subroutine closeall_out

  subroutine glimmer_nc_openappend(outfile,model)
    !*FD open netCDF file for appending
    use parallel
    use glimmer_log
    use glide_types
    use glimmer_map_CFproj
    use glimmer_map_types
    use glimmer_filenames
    implicit none
    type(glimmer_nc_output), pointer :: outfile
    !*FD structure containg output netCDF descriptor
    type(glide_global_type) :: model
    !*FD the model instance

    ! local variables
    integer :: status,timedimid,ntime,timeid
    real(dp),dimension(1) :: last_time
    character(len=msglen) :: message

    ! open existing netCDF file
    status = parallel_open(process_path(NCO%filename),NF90_WRITE,NCO%id)
    call nc_errorhandle(__FILE__,__LINE__,status)
    call write_log_div
    write(message,*) 'Reopening file ',trim(process_path(NCO%filename)),' for output; '
    call write_log(trim(message))
    ! Find out when last time-slice was
    status = parallel_inq_dimid(NCO%id,'time',timedimid)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_inquire_dimension(NCO%id,timedimid,len=ntime)
    call nc_errorhandle(__FILE__,__LINE__,status)
    ! Set timecounter
    outfile%timecounter=ntime+1
    write(message,*) '  Starting output at ',outfile%next_write,' and write every ',outfile%freq,' years'
    call write_log(trim(message))
    
    ! Get time varid
    status = parallel_inq_varid(NCO%id,'time',NCO%timevar)
    call nc_errorhandle(__FILE__,__LINE__,status)

    ! Put dataset into define mode
    status = parallel_redef(NCO%id)
    call nc_errorhandle(__FILE__,__LINE__,status)

  end subroutine glimmer_nc_openappend

  subroutine glimmer_nc_createfile(outfile,model)
    !*FD create a new netCDF file
    use parallel
    use glimmer_log
    use glide_types
    use glimmer_map_CFproj
    use glimmer_map_types
    use glimmer_filenames
    implicit none
    type(glimmer_nc_output), pointer :: outfile
    !*FD structure containg output netCDF descriptor
    type(glide_global_type) :: model
    !*FD the model instance

    ! local variables
    integer status
    integer mapid
    character(len=msglen) message

    ! create new netCDF file
    status = parallel_create(process_path(NCO%filename),NF90_CLOBBER,NCO%id)
    call nc_errorhandle(__FILE__,__LINE__,status)
    call write_log_div
    write(message,*) 'Opening file ',trim(process_path(NCO%filename)),' for output; '
    call write_log(trim(message))
    write(message,*) '  Starting output at ',outfile%next_write,' and write every ',outfile%freq,' years'
    call write_log(trim(message))
    if (outfile%end_write < glimmer_nc_max_time) then
       write(message,*) '  Stop writing at ',outfile%end_write
       call write_log(trim(message))
    end if
    NCO%define_mode=.TRUE.

    ! writing meta data
    status = parallel_put_att(NCO%id, NF90_GLOBAL, 'Conventions', "CF-1.3")
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_put_att(NCO%id, NF90_GLOBAL,'title',trim(outfile%metadata%title))
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_put_att(NCO%id, NF90_GLOBAL,'institution',trim(outfile%metadata%institution))
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_put_att(NCO%id, NF90_GLOBAL,'source',trim(outfile%metadata%source))
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_put_att(NCO%id, NF90_GLOBAL,'history',trim(outfile%metadata%history))
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_put_att(NCO%id, NF90_GLOBAL,'references',trim(outfile%metadata%references))
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_put_att(NCO%id, NF90_GLOBAL,'comment',trim(outfile%metadata%comment))
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_put_att(NCO%id, NF90_GLOBAL,'configuration',trim(outfile%metadata%config))
    call nc_errorhandle(__FILE__,__LINE__,status)
  
    ! defining time dimension and variable
    status = parallel_def_dim(NCO%id,'time',NF90_UNLIMITED,NCO%timedim)
    call nc_errorhandle(__FILE__,__LINE__,status)
    !     time -- Model time
    call write_log('Creating variable time')
    !EIB! lanl version
    !status = nf90_def_var(NCO%id,'time',NF90_FLOAT,(/NCO%timedim/),NCO%timevar)
    !EIB! gc2 version
    status = parallel_def_var(NCO%id,'time',outfile%default_xtype,(/NCO%timedim/),NCO%timevar)
    !EIB! pick one and consistant
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_put_att(NCO%id, NCO%timevar, 'long_name', 'Model time')
    status = parallel_put_att(NCO%id, NCO%timevar, 'standard_name', 'time')
    status = parallel_put_att(NCO%id, NCO%timevar, 'units', 'year since 1-1-1 0:0:0')
    status = parallel_put_att(NCO%id, NCO%timevar, 'calendar', 'none')

    ! adding projection info
    if (glimmap_allocated(model%projection)) then
       status = parallel_def_var(NCO%id,glimmer_nc_mapvarname,NF90_CHAR,mapid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       call glimmap_CFPutProj(NCO%id,mapid,model%projection)
    end if

    ! setting the size of the level and staglevel dimension
    NCO%nlevel = model%general%upn
    NCO%nstaglevel = model%general%upn-1
    NCO%nstagwbndlevel = model%general%upn ! MJH this is the max index, not the size
  end subroutine glimmer_nc_createfile

  subroutine glimmer_nc_checkwrite(outfile,model,forcewrite,time)
    !*FD check if we should write to file
    use parallel
    use glimmer_log
    use glide_types
    use glimmer_filenames
    implicit none
    type(glimmer_nc_output), pointer :: outfile    
    type(glide_global_type) :: model
    logical forcewrite
    real(dp),optional :: time

    character(len=msglen) :: message
    integer status
    real(dp) :: sub_time

    real(dp), parameter :: eps = 1.d-11

    ! Check for optional time argument
    if (present(time)) then
       sub_time=time
    else
       sub_time=model%numerics%time
    end if

    ! check if we are still in define mode and if so leave it
    if (NCO%define_mode) then
       status = parallel_enddef(NCO%id)
       call nc_errorhandle(__FILE__,__LINE__,status)
       NCO%define_mode = .FALSE.
    end if

    if (sub_time > NCO%processsed_time) then
       if (NCO%just_processed) then
          ! finished writing during last time step, need to increase counter...
          
          outfile%timecounter = outfile%timecounter + 1
          status = parallel_sync(NCO%id)
          call nc_errorhandle(__FILE__,__LINE__,status)
          NCO%just_processed = .FALSE.
       end if
    end if

    !WHL - Allow for small roundoff error in computing the time
!!    if (sub_time >= outfile%next_write .or. (forcewrite .and. sub_time > outfile%next_write-outfile%freq)) then  ! prone to roundoff error
    if (sub_time + eps >= outfile%next_write .or. (forcewrite .and. sub_time > outfile%next_write-outfile%freq)) then
       if (sub_time <= outfile%end_write .and. .not.NCO%just_processed) then
          call write_log_div
          write(message,*) 'Writing to file ', trim(process_path(NCO%filename)), ' at time ', sub_time
          call write_log(trim(message))
          ! increase next_write
          outfile%next_write = outfile%next_write + outfile%freq
          NCO%processsed_time = sub_time
          ! write time
          status = parallel_put_var(NCO%id,NCO%timevar,sub_time,(/outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          NCO%just_processed = .TRUE.         
       end if
    end if

  end subroutine glimmer_nc_checkwrite

  !*****************************************************************************
  ! netCDF input
  !*****************************************************************************  
  subroutine openall_in(model)
    !*FD open all netCDF files for input
    use glide_types
    use glimmer_ncdf
    implicit none
    type(glide_global_type) :: model
    
    ! local variables
    type(glimmer_nc_input), pointer :: ic

    ic=>model%funits%in_first
    do while(associated(ic))
       call glimmer_nc_openfile(ic,model)
       ic=>ic%next
    end do
  end subroutine openall_in

  subroutine closeall_in(model)
    !*FD close all netCDF files for input
    use glide_types
    use glimmer_ncdf
    implicit none
    type(glide_global_type) :: model
    
    ! local variables
    type(glimmer_nc_input), pointer :: ic

    ic=>model%funits%in_first
    do while(associated(ic))
       ic=>delete(ic)
    end do
    model%funits%in_first=>NULL()
  end subroutine closeall_in

  subroutine glimmer_nc_openfile(infile,model)
    !*FD open an existing netCDF file
    use glide_types
    use glimmer_map_CFproj
    use glimmer_map_types
    use glimmer_log
    use glimmer_paramets, only: len0
    use glimmer_filenames
    use parallel
    implicit none
    type(glimmer_nc_input), pointer :: infile
    !*FD structure containg input netCDF descriptor
    type(glide_global_type) :: model
    !*FD the model instance

    ! local variables
    integer dimsize, dimid, varid
    real, dimension(2) :: delta
    integer status    
    character(len=msglen) message
    
    real,parameter :: small = 1.e-6

    ! open netCDF file
    status = parallel_open(process_path(NCI%filename),NF90_NOWRITE,NCI%id)
    if (status /= NF90_NOERR) then
       call write_log('Error opening file '//trim(process_path(NCI%filename))//': '//nf90_strerror(status),&
            type=GM_FATAL,file=__FILE__,line=__LINE__)
    end if
    call write_log_div
    call write_log('opening file '//trim(process_path(NCI%filename))//' for input')

    ! getting projection, if none defined already
    if (.not.glimmap_allocated(model%projection)) model%projection = glimmap_CFGetProj(NCI%id)

    ! getting time dimension
    status = parallel_inq_dimid(NCI%id, 'time', NCI%timedim)
    call nc_errorhandle(__FILE__,__LINE__,status)
    ! get id of time variable
    status = parallel_inq_varid(NCI%id,'time',NCI%timevar)
    call nc_errorhandle(__FILE__,__LINE__,status)
    
    ! getting length of time dimension and allocating memory for array containing times
    status = parallel_inquire_dimension(NCI%id,NCI%timedim,len=dimsize)
    call nc_errorhandle(__FILE__,__LINE__,status)
    allocate(infile%times(dimsize))
    infile%nt=dimsize
    status = parallel_get_var(NCI%id,NCI%timevar,infile%times)

    ! setting the size of the level and staglevel dimension
    NCI%nlevel = model%general%upn
    NCI%nstaglevel = model%general%upn-1
    NCI%nstagwbndlevel = model%general%upn !MJH This is the max index, not size

    ! checking if dimensions and grid spacing are the same as in the configuration file
    ! x1
    status = parallel_inq_dimid(NCI%id,'x1',dimid)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_inquire_dimension(NCI%id,dimid,len=dimsize)
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (dimsize /= global_ewn) then
       write(message,*) 'Dimension x1 of file '//trim(process_path(NCI%filename))// &
            ' does not match with config dimension: ', dimsize, global_ewn
       call write_log(message,type=GM_FATAL)
    end if
    status = parallel_inq_varid(NCI%id,'x1',varid)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_get_var(NCI%id,varid,delta)
    call nc_errorhandle(__FILE__,__LINE__,status)

    if (abs(delta(2)-delta(1) - model%numerics%dew*len0) > small) then
       write(message,*) 'deltax1 of file '//trim(process_path(NCI%filename))// &
            ' does not match with config deltax: ', delta(2)-delta(1),model%numerics%dew*len0
       call write_log(message,type=GM_FATAL)
    end if

    ! x0
    !status = nf90_inq_dimid(NCI%id,'x0',dimid)
    !call nc_errorhandle(__FILE__,__LINE__,status)
    !status = nf90_inquire_dimension(NCI%id,dimid,len=dimsize)
    !call nc_errorhandle(__FILE__,__LINE__,status)
    !if (dimsize /= model%general%ewn-1) then
    !   write(message,*) 'Dimension x0 of file ',trim(process_path(NCI%filename)),' does not match with config dimension: ', &
    !        dimsize, model%general%ewn-1
    !   call write_log(message,type=GM_FATAL)
    !end if
    !status = nf90_inq_varid(NCI%id,'x0',varid)
    !call nc_errorhandle(__FILE__,__LINE__,status)
    !status = nf90_get_var(NCI%id,varid,delta)
    !call nc_errorhandle(__FILE__,__LINE__,status)
    !if (abs(delta(2)-delta(1) - model%numerics%dew*len0) > small) then
    !   write(message,*) 'deltax0 of file '//trim(process_path(NCI%filename))//' does not match with config deltax: ', &
    !        delta(2)-delta(1),model%numerics%dew*len0
    !   call write_log(message,type=GM_FATAL)
    !end if

    ! y1
    status = parallel_inq_dimid(NCI%id,'y1',dimid)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_inquire_dimension(NCI%id,dimid,len=dimsize)
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (dimsize /= global_nsn) then
       write(message,*) 'Dimension y1 of file '//trim(process_path(NCI%filename))// &
            ' does not match with config dimension: ', dimsize, global_nsn
       call write_log(message,type=GM_FATAL)
    end if
    status = parallel_inq_varid(NCI%id,'y1',varid)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_get_var(NCI%id,varid,delta)
    call nc_errorhandle(__FILE__,__LINE__,status)

    if (abs(delta(2)-delta(1) - model%numerics%dns*len0) > small) then
       write(message,*) 'deltay1 of file '//trim(process_path(NCI%filename))// &
            ' does not match with config deltay: ', delta(2)-delta(1),model%numerics%dns*len0
       call write_log(message,type=GM_FATAL)
    end if
    
    ! y0
    !status = nf90_inq_dimid(NCI%id,'y0',dimid)
    !call nc_errorhandle(__FILE__,__LINE__,status)
    !status = nf90_inquire_dimension(NCI%id,dimid,len=dimsize)
    !call nc_errorhandle(__FILE__,__LINE__,status)
    !if (dimsize /= model%general%nsn-1) then
    !   write(message,*) 'Dimension y0 of file '//trim(process_path(NCI%filename))//' does not match with config dimension: ',&
    !        dimsize, model%general%nsn-1
    !   call write_log(message,type=GM_FATAL)
    !end if
    !status = nf90_inq_varid(NCI%id,'y0',varid)
    !call nc_errorhandle(__FILE__,__LINE__,status)
    !status = nf90_get_var(NCI%id,varid,delta)
    !call nc_errorhandle(__FILE__,__LINE__,status)
    !if (abs(delta(2)-delta(1) - model%numerics%dns*len0) > small) then
    !   write(message,*) 'deltay0 of file '//trim(process_path(NCI%filename))//' does not match with config deltay: ',&
    !        delta(2)-delta(1),model%numerics%dns*len0
    !   call write_log(message,type=GM_FATAL)
    !end if
  
  ! Check that the number of vertical layers is the same, though it's asking for trouble
  ! to check whether the spacing is the same (don't want to put that burden on setup,
  ! plus f.p. compare has been known to cause problems here)
  status = parallel_inq_dimid(NCI%id,'level',dimid)
  ! If we couldn't find the 'level' dimension fail with a warning.
  ! We don't want to throw an error, as input files are only required to have it if they
  ! include 3D data fields.
  if (status == NF90_NOERR) then
        status = parallel_inquire_dimension(NCI%id, dimid, len=dimsize)
        call nc_errorhandle(__FILE__, __LINE__, status)
        if (dimsize /= model%general%upn .and. dimsize  /=  1) then
            write(message,*) 'Dimension level of file '//trim(process_path(NCI%filename))//&
                ' does not match with config dimension: ', &
                dimsize, model%general%upn
            call write_log(message,type=GM_FATAL)
        end if
  else
        call write_log("Input file contained no level dimension.  This is not necessarily a problem.", type=GM_WARNING)
  end if
  
  end subroutine glimmer_nc_openfile

  subroutine glimmer_nc_checkread(infile,model,time)
    !*FD check if we should read from file
    use glimmer_log
    use glide_types
    use glimmer_filenames
    implicit none
    type(glimmer_nc_input), pointer :: infile  !*FD structure containg output netCDF descriptor
    type(glide_global_type) :: model    !*FD the model instance
    real(dp),optional :: time           !*FD Optional alternative time

    character(len=msglen) :: message
    real(dp) :: sub_time

    integer :: pos  ! to identify restart files

    real(dp) :: restart_time   ! time of restart (yr)

    if (present(time)) then
       sub_time = time
    else
       sub_time = model%numerics%time
    end if

    if (infile%current_time <= infile%nt) then
       if (.not.NCI%just_processed) then
          call write_log_div
          !EIB! added form gc2, needed?
          ! Reset model%numerics%tstart if reading a restart file
          write(message,*) 'Check for restart:', trim(infile%nc%filename)
          call write_log(message)
          pos = index(infile%nc%filename,'.r.')  ! use CESM naming convention for restart files
          if (pos /= 0) then   ! get the start time based on the current time slice
             restart_time = infile%times(infile%current_time)      ! years
             model%numerics%tstart = restart_time
             model%numerics%time = restart_time
             write(message,*) 'Restart: New tstart =', model%numerics%tstart
             call write_log(message)
          endif
          !EIB! end add
          write(message,*) 'Reading time slice ',infile%current_time,'(',infile%times(infile%current_time),') from file ', &
               trim(process_path(NCI%filename)), ' at time ', sub_time
          call write_log(message)
          NCI%just_processed = .TRUE.
          NCI%processsed_time = sub_time
       end if
    end if

    if (sub_time > NCI%processsed_time) then
       if (NCI%just_processed) then
          ! finished reading during last time step, need to increase counter...
          infile%current_time = infile%current_time + 1
          NCI%just_processed = .FALSE.
       end if
    end if

  end subroutine glimmer_nc_checkread

!------------------------------------------------------------------------------

    subroutine check_for_tempstag(whichdycore, nc)
      ! Check for the need to output tempstag and update the output variables if needed.
      !
      ! For the glam/glissade dycore, the vertical temperature grid has an extra level.
      ! In that case, the netCDF output file should include a variable
      ! called tempstag(0:nz) instead of temp(1:nz). This subroutine is added for
      ! convenience to allow the variable "temp" to be specified in the config
      ! file in all cases and have it converted to "tempstag" when appropriate.
      ! MJH

      use glimmer_log
      use glide_types

      implicit none
      integer, intent(in) :: whichdycore
      type(glimmer_nc_stat) :: nc

      ! Locals
      integer :: i

      ! Check if tempstag should be output

      ! TODO If both temp and tempstag are specified, should one be removed?
      ! TODO Modify this to work if multiple output files are specified?

      !print *, "Original varstring:", varstring

      if (whichdycore==DYCORE_GLAM .or. whichdycore==DYCORE_GLISSADE) then 
          ! We want temp to become tempstag
          i = index(nc%vars, " temp ")
          if (i > 0) then
            ! temp was specified - change it to tempstag
            ! If temp is listed more than once, this just changes the first instance
            nc%vars = nc%vars(1:i-1) // " tempstag " // nc%vars(i+6:len(nc%vars))
            call write_log('Temperature remapping option uses temperature on a staggered grid.' // &
              '  The netCDF output variable "temp" has been changed to "tempstag".' )
          endif
          ! Now check if flwa needs to be changed to flwastag
          i = index(nc%vars, " flwa ") ! Look for flwa
          if (i > 0) then
            ! flwa was specified - change to flwastag
            nc%vars = nc%vars(1:i-1) // " flwastag " // nc%vars(i+6:len(nc%vars))
            call write_log('Temperature remapping option uses flwa on a staggered grid.' // &
            '  The netCDF output variable "flwa" has been changed to "flwastag".' )
          endif
      else  ! glide dycore
          ! We want tempstag to become temp
          i = index(nc%vars, " tempstag ")
          if (i > 0) then
            !Change tempstag to temp
            nc%vars = nc%vars(1:i-1) // " temp " // nc%vars(i+10:len(nc%vars))
            call write_log('The netCDF output variable "tempstag" should only be used when remapping temperature.' // &
              '  The netCDF output variable "tempstag" has been changed to "temp".' )
          endif
          ! We want flwastag to become flwa
          i = index(nc%vars, " flwastag ")
          if (i > 0) then
            !Change flwastag to flwa
            nc%vars = nc%vars(1:i-1) // " flwa " // nc%vars(i+10:len(nc%vars))
            call write_log('The netCDF output variable "flwastag" should only be used when remapping temperature.' // &
              '  The netCDF output variable "flwastag" has been changed to "flwa".' )
          endif
      endif  ! whichdycore

      ! Copy any changes to vars_copy
      nc%vars_copy = nc%vars

    end subroutine check_for_tempstag

!------------------------------------------------------------------------------

end module glimmer_ncio

!------------------------------------------------------------------------------
