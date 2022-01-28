!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WARNING: this file was automatically generated on
! Fri, 08 Jan 2016 19:15:11 +0000
! from ncdf_template.F90.in
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! WJS (1-30-12): The following (turning optimization off) is needed as a workaround for an
! xlf compiler bug, at least in IBM XL Fortran for AIX, V12.1 on bluefire
#ifdef CPRIBM
@PROCESS OPT(0)
#endif

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   ncdf_template.F90.in - part of the Community Ice Sheet Model (CISM)  
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

#define NCO outfile%nc
#define NCI infile%nc


module glad_io
  ! template for creating subsystem specific I/O routines
  ! written by Magnus Hagdorn, 2004

  use glad_type

  implicit none

  private :: get_xtype, is_enabled, is_enabled_0dint, is_enabled_1dint, &
       is_enabled_2dint, is_enabled_0dreal, is_enabled_1dreal, is_enabled_2dreal, is_enabled_3dreal

  character(310), save :: restart_variable_list=''    ! list of variables needed for a restart
!TODO change 310 to a variable - see glimmer_ncdf.F90 in the definition for type glimmer_nc_stat for other instances of this value.

  interface is_enabled  ! MJH 10/21/13: Interface needed for determining if arrays have been enabled.  See notes below in glad_io_create.
    module procedure is_enabled_0dint
    module procedure is_enabled_1dint
    module procedure is_enabled_2dint
    module procedure is_enabled_0dreal
    module procedure is_enabled_1dreal
    module procedure is_enabled_2dreal
    module procedure is_enabled_3dreal
  end interface is_enabled

contains

  !*****************************************************************************
  ! netCDF output
  !*****************************************************************************
  subroutine glad_io_createall(model,data,outfiles)
    ! open all netCDF files for output
    use glad_type
    use glide_types
    use glimmer_ncdf
    use glimmer_ncio
    implicit none
    type(glide_global_type) :: model
    type(glad_instance) :: data ! MJH 10/21/13: Making 'data' mandatory.  See notes below in glad_io_create
    type(glimmer_nc_output),optional,pointer :: outfiles
    
    ! local variables
    type(glimmer_nc_output), pointer :: oc

    if (present(outfiles)) then
       oc => outfiles
    else
       oc=>model%funits%out_first
    end if

    do while(associated(oc))
       call glad_io_create(oc,model,data)
       oc=>oc%next
    end do
  end subroutine glad_io_createall

  subroutine glad_io_writeall(data,model,atend,outfiles,time)
    ! if necessary write to netCDF files
    use glad_type
    use glide_types
    use glimmer_ncdf
    use glimmer_ncio
    implicit none
    type(glad_instance) :: data
    type(glide_global_type) :: model
    logical, optional :: atend
    type(glimmer_nc_output),optional,pointer :: outfiles
    real(dp),optional :: time

    ! local variables
    type(glimmer_nc_output), pointer :: oc
    logical :: forcewrite=.false.

    if (present(outfiles)) then
       oc => outfiles
    else
       oc=>model%funits%out_first
    end if

    if (present(atend)) then
       forcewrite = atend
    end if

    do while(associated(oc))
#ifdef HAVE_AVG
       if (oc%do_averages) then
          call glad_avg_accumulate(oc,data,model)
       end if
#endif
       call glimmer_nc_checkwrite(oc,model,forcewrite,time)
       if (oc%nc%just_processed) then
          ! write standard variables
          call glad_io_write(oc,data)
#ifdef HAVE_AVG
          if (oc%do_averages) then
             call glad_avg_reset(oc,data)
          end if
#endif
       end if
       oc=>oc%next
    end do
  end subroutine glad_io_writeall
  
  subroutine glad_io_create(outfile,model,data)
    use parallel
    use glide_types
    use glad_type
    use glimmer_ncdf
    use glimmer_ncio
    use glimmer_map_types
    use glimmer_log
    use glimmer_paramets
    use glimmer_scales
    use glimmer_log
    implicit none
    type(glimmer_nc_output), pointer :: outfile
    type(glide_global_type) :: model
    type(glad_instance) :: data    ! MJH 10/21/13: Making 'data' mandatory.  See note below

    integer status,varid,pos

    ! MJH 10/21/13: Local variables needed for checking if a variable is enabled.
    real(dp) :: tavgf
    integer :: up

    integer :: time_dimid
    integer :: x1_dimid
    integer :: y1_dimid

    ! defining dimensions
    status = parallel_inq_dimid(NCO%id,'time',time_dimid)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_inq_dimid(NCO%id,'x1',x1_dimid)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_inq_dimid(NCO%id,'y1',y1_dimid)
    call nc_errorhandle(__FILE__,__LINE__,status)

    ! Expanding restart variables: if 'restart' or 'hot' is present, we remove that
    ! word from the variable list, and flip the restartfile flag.  
    ! In CISM 2.0, 'restart' is the preferred name to represent restart variables, 
    ! but 'hot' is supported for backward compatibility.  Thus, we check for both.
    NCO%vars = ' '//trim(adjustl(NCO%vars))//' '  ! Need to maintain a space at beginning and end of list
    ! expanding restart variables
    pos = index(NCO%vars,' restart ') 
    if (pos.ne.0) then
       NCO%vars = NCO%vars(:pos)//NCO%vars(pos+8:)
       NCO%restartfile = .true.
    end if
    pos = index(NCO%vars,' hot ') 
    if (pos.ne.0) then
       NCO%vars = NCO%vars(:pos)//NCO%vars(pos+4:)
       NCO%restartfile = .true.
    end if
    ! Now apply necessary changes if the file is a restart file.
    if (NCO%restartfile) then
       if ((len_trim(NCO%vars) + len_trim(restart_variable_list) + 2) >= len(NCO%vars) ) then
          call write_log('Adding restart variables has made the list of output variables too long for file ' // NCO%filename, &
               GM_FATAL)
       else
          ! Expand the restart variable list 
          ! Need to maintain a space at beginning and end of list
          NCO%vars = trim(NCO%vars) // ' ' // trim(restart_variable_list) // ' ' ! (a module variable)  
          ! Set the xtype to be double (required for an exact restart)
          outfile%default_xtype = NF90_DOUBLE   
       endif
    end if

    ! Convert temp and flwa to versions on stag grid, if needed
    ! Note: this check must occur after restart variables are expanded which happens in glimmer_nc_readparams
    call check_for_tempstag(model%options%whichdycore,NCO)

    ! checking if we need to handle time averages
    pos = index(NCO%vars,"_tavg")
    if (pos.ne.0) then
       outfile%do_averages = .True.
    end if    

    ! Now that the output variable list is finalized, make sure we aren't truncating what the user intends to be output.
    ! Note: this only checks that the text in the variable list does not extend to within one character of the end of the variable.
    ! It does not handle the case where the user exactly fills the allowable length with variables or has a too-long list with more than one space between variable names.
    if ((len_trim(NCO%vars) + 1 ) >= len(NCO%vars)) then 
       call write_log('The list of output variables is too long for file ' // NCO%filename, GM_FATAL)
    endif


    ! MJH, 10/21/13: In the auto-generated code below, the creation of each output variable is wrapped by a check if the data for that 
    !   variable has a size greater than 0.  This is because of recently added checks in glide_types.F90 that don't fully allocate
    !   some variables if certain model options are disabled.  This is to lower memory requirements while running the model.
    !   The reason they have to be allocated with size zero rather than left unallocated is because the data for
    !   some netCDF output variables is defined with math, which causes an error if the operands are unallocated.
    !   Note that if a variable is not created, then it will not be subsequently written to.
    !   Also note that this change requires that data be a mandatory argument to this subroutine.

    ! Some output variables will need tavgf.  The value does not matter, but it must exist.  
    ! Nonetheless, for completeness give it the proper value that it has in glad_io_write.
    tavgf = outfile%total_time
    if (tavgf.ne.0.d0) then
       tavgf = 1.d0/tavgf
    end if
    ! Similarly, some output variables use the variable up.  Give it value of 0 here.
    up = 0

    !     hflx_tavg -- heat flux to ice surface
    pos = index(NCO%vars,' hflx_tavg ')
    status = parallel_inq_varid(NCO%id,'hflx_tavg',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+9) = '         '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%hflx_tavg)) then
       call write_log('Creating variable hflx_tavg')
       status = parallel_def_var(NCO%id,'hflx_tavg',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'heat flux to ice surface')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'W m-2')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable hflx_tavg was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     lat -- latitude
    pos = index(NCO%vars,' lat ')
    status = parallel_inq_varid(NCO%id,'lat',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+3) = '   '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%lat)) then
       call write_log('Creating variable lat')
       status = parallel_def_var(NCO%id,'lat',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'latitude')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'latitude')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'degreeN')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable lat was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     lon -- longitude
    pos = index(NCO%vars,' lon ')
    status = parallel_inq_varid(NCO%id,'lon',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+3) = '   '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%lon)) then
       call write_log('Creating variable lon')
       status = parallel_def_var(NCO%id,'lon',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'longitude')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'longitude')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'degreeE')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable lon was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     rofi_tavg -- solid calving flux
    pos = index(NCO%vars,' rofi_tavg ')
    status = parallel_inq_varid(NCO%id,'rofi_tavg',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+9) = '         '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%rofi_tavg)) then
       call write_log('Creating variable rofi_tavg')
       status = parallel_def_var(NCO%id,'rofi_tavg',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'solid calving flux')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'kg m-2 s-1')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable rofi_tavg was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     rofl_tavg -- liquid runoff flux
    pos = index(NCO%vars,' rofl_tavg ')
    status = parallel_inq_varid(NCO%id,'rofl_tavg',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+9) = '         '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%rofl_tavg)) then
       call write_log('Creating variable rofl_tavg')
       status = parallel_def_var(NCO%id,'rofl_tavg',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'liquid runoff flux')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'kg m-2 s-1')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable rofl_tavg was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

  end subroutine glad_io_create

  subroutine glad_io_write(outfile,data)
    use parallel
    use glad_type
    use glimmer_ncdf
    use glimmer_paramets
    use glimmer_scales
    implicit none
    type(glimmer_nc_output), pointer :: outfile
    ! structure containg output netCDF descriptor
    type(glad_instance) :: data
    ! the model instance

    ! local variables
    real(dp) :: tavgf
    integer status, varid
    integer up
     
    tavgf = outfile%total_time
    if (tavgf.ne.0.d0) then
       tavgf = 1.d0/tavgf
    end if

    ! write variables
    status = parallel_inq_varid(NCO%id,'hflx_tavg',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%hflx_tavg, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'lat',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%lat, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'lon',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%lon, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'rofi_tavg',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%rofi_tavg, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'rofl_tavg',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%rofl_tavg, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

  end subroutine glad_io_write


  subroutine glad_add_to_restart_variable_list(vars_to_add)
    ! This subroutine adds variables to the list of variables needed for a restart.
    ! It is a public subroutine that allows other parts of the model to modify the list, 
    ! which is a module level variable.   MJH 1/17/2013

    use glimmer_log
    implicit none

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------
    character(len=*), intent (in) :: vars_to_add  ! list of variable(s) to be added to the list of restart variables 
    !character(*), intent (inout) :: restart_variable_list  ! list of variables needed to perform an exact restart - module variable

    !------------------------------------------------------------------------------------
    ! Internal variables
    !------------------------------------------------------------------------------------

    !------------------------------------------------------------------------------------

    ! Add the variables to the list so long as they don't make the list too long.
    if ( (len_trim(restart_variable_list) + 1 + len_trim(vars_to_add)) > len(restart_variable_list)) then
       call write_log('Adding restart variables has made the restart variable list too long.',GM_FATAL)
    else
       restart_variable_list = trim(adjustl(restart_variable_list)) // ' ' // trim(vars_to_add)
       !call write_log('Adding to glad restart variable list: ' // trim(vars_to_add) )
    endif

  end subroutine glad_add_to_restart_variable_list


  ! Functions for the interface 'is_enabled'.  These are needed by the auto-generated code in glad_io_create
  !   to determine if a variable is 'turned on', and should be written.

  function is_enabled_0dint(var)
    integer, intent(in) :: var
    logical :: is_enabled_0dint
    is_enabled_0dint = .true.  ! scalars are always enabled
    return
  end function is_enabled_0dint

  function is_enabled_1dint(var)
    integer, dimension(:), pointer, intent(in) :: var
    logical :: is_enabled_1dint
    if (associated(var)) then
      is_enabled_1dint = .true.
    else
      is_enabled_1dint = .false.
    endif
    return
  end function is_enabled_1dint

  function is_enabled_2dint(var)
    integer, dimension(:,:), pointer, intent(in) :: var
    logical :: is_enabled_2dint
    if (associated(var)) then
      is_enabled_2dint = .true.
    else
      is_enabled_2dint = .false.
    endif
    return
  end function is_enabled_2dint

  function is_enabled_0dreal(var)
    real(dp), intent(in) :: var
    logical :: is_enabled_0dreal
    is_enabled_0dreal = .true.  ! scalars are always enabled
    return
  end function is_enabled_0dreal

  function is_enabled_1dreal(var)
    real(dp), dimension(:), pointer, intent(in) :: var
    logical :: is_enabled_1dreal
    if (associated(var)) then
      is_enabled_1dreal = .true.
    else
      is_enabled_1dreal = .false.
    endif
    return
  end function is_enabled_1dreal

  function is_enabled_2dreal(var)
    real(dp), dimension(:,:), pointer, intent(in) :: var
    logical :: is_enabled_2dreal
    if (associated(var)) then
      is_enabled_2dreal = .true.
    else
      is_enabled_2dreal = .false.
    endif
    return
  end function is_enabled_2dreal

  function is_enabled_3dreal(var)
    real(dp), dimension(:,:,:), pointer, intent(in) :: var
    logical :: is_enabled_3dreal
    if (associated(var)) then
      is_enabled_3dreal = .true.
    else
      is_enabled_3dreal = .false.
    endif
    return
  end function is_enabled_3dreal


  !*****************************************************************************
  ! netCDF input
  !*****************************************************************************  
  subroutine glad_io_readall(data, model, filetype)
    ! read from netCDF file
    use glad_type
    use glide_types
    use glimmer_ncdf
    use glimmer_ncio
    implicit none
    type(glad_instance) :: data
    type(glide_global_type) :: model
    integer, intent(in), optional :: filetype  ! 0 for input, 1 for forcing; defaults to input

    ! local variables
    type(glimmer_nc_input), pointer :: ic
    integer :: filetype_local

    if (present(filetype)) then
      filetype_local = filetype
    else
      filetype_local = 0 ! default to input type
    end if

    if (filetype_local == 0) then
      ic=>model%funits%in_first
    else
      ic=>model%funits%frc_first
    endif
    do while(associated(ic))
       call glimmer_nc_checkread(ic,model)
       if (ic%nc%just_processed) then
          call glad_io_read(ic,data)
       end if
       ic=>ic%next
    end do
  end subroutine glad_io_readall


  subroutine glad_read_forcing(data, model)
    ! Read data from forcing files
    use glimmer_log
    use glide_types
    use glimmer_ncdf

    implicit none
    type(glad_instance) :: data
    type(glide_global_type), intent(inout) :: model

    ! Locals
    type(glimmer_nc_input), pointer :: ic
    integer :: t
    real(dp) :: eps ! a tolerance to use for stepwise constant forcing

    ! Make eps a fraction of the time step.
    eps = model%numerics%tinc * 1.0d-4

    ! read forcing files
    ic=>model%funits%frc_first
    do while(associated(ic))

       !print *, 'possible forcing times', ic%times

       ic%nc%just_processed = .true. ! until we find an acceptable time, set this to true which will prevent the file from being read.
       ! Find the current time in the file
       do t = ic%nt, 1, -1  ! look through the time array backwards
          if ( ic%times(t) <= model%numerics%time + eps) then
            ! use the largest time that is smaller or equal to the current time (stepwise forcing)

            ! Set the desired time to be read
            ic%current_time = t
            ic%nc%just_processed = .false.  ! set this to false so file will be read.
            !print *, 'time, forcing index, forcing time', model%numerics%time, ic%current_time, ic%times(ic%current_time)
            exit  ! once we find the time, exit the loop
          endif
          end do  ! if we get to end of loop without exiting, then this file will not be read at this time.

       ! move on to the next forcing file
       ic=>ic%next
    end do

    ! Now that we've updated metadata for each forcing file, actually perform the read.
    ! This call will only read forcing files where just_processed=.false.
    call glad_io_readall(data, model, filetype=1)

  end subroutine glad_read_forcing


!------------------------------------------------------------------------------


  subroutine glad_io_read(infile,data)
    ! read variables from a netCDF file
    use parallel
    use glimmer_log
    use glimmer_ncdf
    use glad_type
    use glimmer_paramets
    use glimmer_scales
    implicit none
    type(glimmer_nc_input), pointer :: infile
    ! structure containg output netCDF descriptor
    type(glad_instance) :: data
    ! the model instance

    ! local variables
    integer status,varid
    integer up
    real(dp) :: scaling_factor

    ! read variables
    status = parallel_inq_varid(NCI%id,'hflx_tavg',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%hflx_tavg)) then
       call write_log('  Loading hflx_tavg')
       status = distributed_get_var(NCI%id, varid, &
            data%hflx_tavg, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling hflx_tavg",GM_DIAGNOSTIC)
          data%hflx_tavg = &
               data%hflx_tavg*scaling_factor
       end if
    else
    call write_log('Variable hflx_tavg was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'lat',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%lat)) then
       call write_log('  Loading lat')
       status = distributed_get_var(NCI%id, varid, &
            data%lat, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling lat",GM_DIAGNOSTIC)
          data%lat = &
               data%lat*scaling_factor
       end if
    else
    call write_log('Variable lat was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'lon',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%lon)) then
       call write_log('  Loading lon')
       status = distributed_get_var(NCI%id, varid, &
            data%lon, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling lon",GM_DIAGNOSTIC)
          data%lon = &
               data%lon*scaling_factor
       end if
    else
    call write_log('Variable lon was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'rofi_tavg',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%rofi_tavg)) then
       call write_log('  Loading rofi_tavg')
       status = distributed_get_var(NCI%id, varid, &
            data%rofi_tavg, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling rofi_tavg",GM_DIAGNOSTIC)
          data%rofi_tavg = &
               data%rofi_tavg*scaling_factor
       end if
    else
    call write_log('Variable rofi_tavg was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'rofl_tavg',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%rofl_tavg)) then
       call write_log('  Loading rofl_tavg')
       status = distributed_get_var(NCI%id, varid, &
            data%rofl_tavg, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling rofl_tavg",GM_DIAGNOSTIC)
          data%rofl_tavg = &
               data%rofl_tavg*scaling_factor
       end if
    else
    call write_log('Variable rofl_tavg was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

  end subroutine glad_io_read

  subroutine glad_io_checkdim(infile,model,data)
    ! check if dimension sizes in file match dims of model
    use parallel
    use glimmer_log
    use glimmer_ncdf
    use glide_types
    use glad_type
    implicit none
    type(glimmer_nc_input), pointer :: infile
    ! structure containg output netCDF descriptor
    type(glide_global_type) :: model
    type(glad_instance), optional :: data

    integer status,dimid,dimsize
    character(len=150) message

    ! check dimensions
  end subroutine glad_io_checkdim

  !*****************************************************************************
  ! calculating time averages
  !*****************************************************************************  
#ifdef HAVE_AVG
  subroutine glad_avg_accumulate(outfile,data,model)
    use parallel
    use glide_types
    use glad_type
    use glimmer_ncdf
    implicit none
    type(glimmer_nc_output), pointer :: outfile
    ! structure containg output netCDF descriptor
    type(glide_global_type) :: model
    type(glad_instance) :: data

    ! local variables
    real(dp) :: factor
    integer status, varid

    ! increase total time
    outfile%total_time = outfile%total_time + model%numerics%tinc
    factor = model%numerics%tinc

  end subroutine glad_avg_accumulate

  subroutine glad_avg_reset(outfile,data)
    use parallel
    use glad_type
    use glimmer_ncdf
    implicit none
    type(glimmer_nc_output), pointer :: outfile
    ! structure containg output netCDF descriptor
    type(glad_instance) :: data

    ! local variables
    integer status, varid

    ! reset total time
    outfile%total_time = 0.d0

  end subroutine glad_avg_reset
#endif

  !*********************************************************************
  ! some private procedures
  !*********************************************************************

  !> apply default type to be used in netCDF file
  integer function get_xtype(outfile,xtype)
    use glimmer_ncdf
    implicit none
    type(glimmer_nc_output), pointer :: outfile !< derived type holding information about output file
    integer, intent(in) :: xtype                !< the external netCDF type

    get_xtype = xtype
    
    if (xtype.eq.NF90_REAL .and. outfile%default_xtype.eq.NF90_DOUBLE) then
       get_xtype = NF90_DOUBLE
    end if
    if (xtype.eq.NF90_DOUBLE .and. outfile%default_xtype.eq.NF90_REAL) then
       get_xtype = NF90_REAL
    end if
  end function get_xtype

  !*********************************************************************
  ! lots of accessor subroutines follow
  !*********************************************************************
  subroutine glad_get_lat(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glad_type
    implicit none
    type(glad_instance) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%lat
  end subroutine glad_get_lat

  subroutine glad_set_lat(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glad_type
    implicit none
    type(glad_instance) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%lat = inarray
  end subroutine glad_set_lat

  subroutine glad_get_lon(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glad_type
    implicit none
    type(glad_instance) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%lon
  end subroutine glad_get_lon

  subroutine glad_set_lon(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glad_type
    implicit none
    type(glad_instance) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%lon = inarray
  end subroutine glad_set_lon


end module glad_io
