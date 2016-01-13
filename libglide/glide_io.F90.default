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


module glide_io
  ! template for creating subsystem specific I/O routines
  ! written by Magnus Hagdorn, 2004

  use glide_types

  implicit none

  private :: get_xtype, is_enabled, is_enabled_0dint, is_enabled_1dint, &
       is_enabled_2dint, is_enabled_0dreal, is_enabled_1dreal, is_enabled_2dreal, is_enabled_3dreal

  character(310), save :: restart_variable_list=''    ! list of variables needed for a restart
!TODO change 310 to a variable - see glimmer_ncdf.F90 in the definition for type glimmer_nc_stat for other instances of this value.

  interface is_enabled  ! MJH 10/21/13: Interface needed for determining if arrays have been enabled.  See notes below in glide_io_create.
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
  subroutine glide_io_createall(model,data,outfiles)
    ! open all netCDF files for output
    use glide_types
    use glide_types
    use glimmer_ncdf
    use glimmer_ncio
    implicit none
    type(glide_global_type) :: model
    type(glide_global_type) :: data ! MJH 10/21/13: Making 'data' mandatory.  See notes below in glide_io_create
    type(glimmer_nc_output),optional,pointer :: outfiles
    
    ! local variables
    type(glimmer_nc_output), pointer :: oc

    if (present(outfiles)) then
       oc => outfiles
    else
       oc=>model%funits%out_first
    end if

    do while(associated(oc))
       call glide_io_create(oc,model,data)
       oc=>oc%next
    end do
  end subroutine glide_io_createall

  subroutine glide_io_writeall(data,model,atend,outfiles,time)
    ! if necessary write to netCDF files
    use glide_types
    use glide_types
    use glimmer_ncdf
    use glimmer_ncio
    implicit none
    type(glide_global_type) :: data
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
          call glide_avg_accumulate(oc,data,model)
       end if
#endif
       call glimmer_nc_checkwrite(oc,model,forcewrite,time)
       if (oc%nc%just_processed) then
          ! write standard variables
          call glide_io_write(oc,data)
#ifdef HAVE_AVG
          if (oc%do_averages) then
             call glide_avg_reset(oc,data)
          end if
#endif
       end if
       oc=>oc%next
    end do
  end subroutine glide_io_writeall
  
  subroutine glide_io_create(outfile,model,data)
    use parallel
    use glide_types
    use glide_types
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
    type(glide_global_type) :: data    ! MJH 10/21/13: Making 'data' mandatory.  See note below

    integer status,varid,pos

    ! MJH 10/21/13: Local variables needed for checking if a variable is enabled.
    real(dp) :: tavgf
    integer :: up

    integer :: level_dimid
    integer :: lithoz_dimid
    integer :: staglevel_dimid
    integer :: stagwbndlevel_dimid
    integer :: time_dimid
    integer :: x0_dimid
    integer :: x1_dimid
    integer :: y0_dimid
    integer :: y1_dimid

    ! defining dimensions
    if (.not.outfile%append) then
       status = parallel_def_dim(NCO%id,'level',model%general%upn,level_dimid)
    else
       status = parallel_inq_dimid(NCO%id,'level',level_dimid)
    endif
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (.not.outfile%append) then
       status = parallel_def_dim(NCO%id,'lithoz',model%lithot%nlayer,lithoz_dimid)
    else
       status = parallel_inq_dimid(NCO%id,'lithoz',lithoz_dimid)
    endif
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (.not.outfile%append) then
       status = parallel_def_dim(NCO%id,'staglevel',model%general%upn-1,staglevel_dimid)
    else
       status = parallel_inq_dimid(NCO%id,'staglevel',staglevel_dimid)
    endif
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (.not.outfile%append) then
       status = parallel_def_dim(NCO%id,'stagwbndlevel',model%general%upn+1,stagwbndlevel_dimid)
    else
       status = parallel_inq_dimid(NCO%id,'stagwbndlevel',stagwbndlevel_dimid)
    endif
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_inq_dimid(NCO%id,'time',time_dimid)
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (.not.outfile%append) then
       status = parallel_def_dim(NCO%id,'x0',global_ewn-1,x0_dimid)
    else
       status = parallel_inq_dimid(NCO%id,'x0',x0_dimid)
    endif
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (.not.outfile%append) then
       status = parallel_def_dim(NCO%id,'x1',global_ewn,x1_dimid)
    else
       status = parallel_inq_dimid(NCO%id,'x1',x1_dimid)
    endif
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (.not.outfile%append) then
       status = parallel_def_dim(NCO%id,'y0',global_nsn-1,y0_dimid)
    else
       status = parallel_inq_dimid(NCO%id,'y0',y0_dimid)
    endif
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (.not.outfile%append) then
       status = parallel_def_dim(NCO%id,'y1',global_nsn,y1_dimid)
    else
       status = parallel_inq_dimid(NCO%id,'y1',y1_dimid)
    endif
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
    ! Nonetheless, for completeness give it the proper value that it has in glide_io_write.
    tavgf = outfile%total_time
    if (tavgf.ne.0.d0) then
       tavgf = 1.d0/tavgf
    end if
    ! Similarly, some output variables use the variable up.  Give it value of 0 here.
    up = 0

    !     level -- sigma layers
    if (.not.outfile%append) then
       call write_log('Creating variable level')
       status = parallel_def_var(NCO%id,'level',get_xtype(outfile,NF90_FLOAT), &
            (/level_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'positive', &
            'down')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'sigma layers')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_sigma_coordinate')
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1')
     end if

    !     lithoz -- vertical coordinate of lithosphere layer
    if (.not.outfile%append) then
       call write_log('Creating variable lithoz')
       status = parallel_def_var(NCO%id,'lithoz',get_xtype(outfile,NF90_FLOAT), &
            (/lithoz_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'vertical coordinate of lithosphere layer')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter')
     end if

    !     staglevel -- stag sigma layers
    if (.not.outfile%append) then
       call write_log('Creating variable staglevel')
       status = parallel_def_var(NCO%id,'staglevel',get_xtype(outfile,NF90_FLOAT), &
            (/staglevel_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'positive', &
            'down')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'stag sigma layers')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_stag_sigma_coordinate')
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1')
     end if

    !     stagwbndlevel -- stag sigma layers with boundaries
    if (.not.outfile%append) then
       call write_log('Creating variable stagwbndlevel')
       status = parallel_def_var(NCO%id,'stagwbndlevel',get_xtype(outfile,NF90_FLOAT), &
            (/stagwbndlevel_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'positive', &
            'down')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'stag sigma layers with boundaries')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_stag_sigma_coordinate_with_bnd')
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1')
     end if

    !     x0 -- Cartesian x-coordinate, velocity grid
    if (.not.outfile%append) then
       call write_log('Creating variable x0')
       status = parallel_def_var(NCO%id,'x0',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'Cartesian x-coordinate, velocity grid')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter')
       status = parallel_put_att(NCO%id, varid, 'axis', &
            'X')
     end if

    !     x1 -- Cartesian x-coordinate
    if (.not.outfile%append) then
       call write_log('Creating variable x1')
       status = parallel_def_var(NCO%id,'x1',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'Cartesian x-coordinate')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter')
       status = parallel_put_att(NCO%id, varid, 'axis', &
            'X')
     end if

    !     y0 -- Cartesian y-coordinate, velocity grid
    if (.not.outfile%append) then
       call write_log('Creating variable y0')
       status = parallel_def_var(NCO%id,'y0',get_xtype(outfile,NF90_FLOAT), &
            (/y0_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'Cartesian y-coordinate, velocity grid')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter')
       status = parallel_put_att(NCO%id, varid, 'axis', &
            'Y')
     end if

    !     y1 -- Cartesian y-coordinate
    if (.not.outfile%append) then
       call write_log('Creating variable y1')
       status = parallel_def_var(NCO%id,'y1',get_xtype(outfile,NF90_FLOAT), &
            (/y1_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'Cartesian y-coordinate')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter')
       status = parallel_put_att(NCO%id, varid, 'axis', &
            'Y')
     end if

    !     C_space_factor -- spatial factor for basal shear stress
    pos = index(NCO%vars,' C_space_factor ')
    status = parallel_inq_varid(NCO%id,'C_space_factor',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+14) = '              '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%basal_physics%C_space_factor)) then
       call write_log('Creating variable C_space_factor')
       status = parallel_def_var(NCO%id,'C_space_factor',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'spatial factor for basal shear stress')
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable C_space_factor was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     acab -- accumulation, ablation rate
    pos = index(NCO%vars,' acab ')
    status = parallel_inq_varid(NCO%id,'acab',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%climate%acab)) then
       call write_log('Creating variable acab')
       status = parallel_def_var(NCO%id,'acab',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_acab))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'accumulation, ablation rate')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_surface_specific_mass_balance')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable acab was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     adv_cfl_dt -- advective CFL maximum time step
    pos = index(NCO%vars,' adv_cfl_dt ')
    status = parallel_inq_varid(NCO%id,'adv_cfl_dt',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+10) = '          '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%numerics%adv_cfl_dt)) then
       call write_log('Creating variable adv_cfl_dt')
       status = parallel_def_var(NCO%id,'adv_cfl_dt',get_xtype(outfile,NF90_FLOAT), &
            (/time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'advective CFL maximum time step')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'years')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable adv_cfl_dt was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     artm -- annual mean air temperature
    pos = index(NCO%vars,' artm ')
    status = parallel_inq_varid(NCO%id,'artm',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%climate%artm)) then
       call write_log('Creating variable artm')
       status = parallel_def_var(NCO%id,'artm',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'annual mean air temperature')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'surface_temperature')
       status = parallel_put_att(NCO%id, varid, 'cell_methods', &
            'time: mean')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'degree_Celsius')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable artm was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     beta -- higher-order bed stress coefficient
    pos = index(NCO%vars,' beta ')
    status = parallel_inq_varid(NCO%id,'beta',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%beta)) then
       call write_log('Creating variable beta')
       status = parallel_def_var(NCO%id,'beta',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_beta))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'higher-order bed stress coefficient')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa yr/m')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable beta was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     beta_internal -- weighted higher-order bed stress coefficient
    pos = index(NCO%vars,' beta_internal ')
    status = parallel_inq_varid(NCO%id,'beta_internal',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+13) = '             '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%beta_internal)) then
       call write_log('Creating variable beta_internal')
       status = parallel_def_var(NCO%id,'beta_internal',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_beta))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'weighted higher-order bed stress coefficient')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa yr/m')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable beta_internal was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     bfricflx -- basal friction heat flux
    pos = index(NCO%vars,' bfricflx ')
    status = parallel_inq_varid(NCO%id,'bfricflx',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+8) = '        '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%temper%bfricflx)) then
       call write_log('Creating variable bfricflx')
       status = parallel_def_var(NCO%id,'bfricflx',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(1.0))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal friction heat flux')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'watt/meter2')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable bfricflx was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     bheatflx -- upward basal heat flux
    pos = index(NCO%vars,' bheatflx ')
    status = parallel_inq_varid(NCO%id,'bheatflx',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+8) = '        '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%temper%bheatflx)) then
       call write_log('Creating variable bheatflx')
       status = parallel_def_var(NCO%id,'bheatflx',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_bflx))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'upward basal heat flux')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'watt/meter2')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable bheatflx was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     bmlt_float -- basal melt rate for floating ice
    pos = index(NCO%vars,' bmlt_float ')
    status = parallel_inq_varid(NCO%id,'bmlt_float',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+10) = '          '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%temper%bmlt_float)) then
       call write_log('Creating variable bmlt_float')
       status = parallel_def_var(NCO%id,'bmlt_float',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_acab))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal melt rate for floating ice')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_basal_melt_rate_floating')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable bmlt_float was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     bmlt_float_mask -- mask for basal melting of floating ice
    pos = index(NCO%vars,' bmlt_float_mask ')
    status = parallel_inq_varid(NCO%id,'bmlt_float_mask',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+15) = '               '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%temper%bmlt_float_mask)) then
       call write_log('Creating variable bmlt_float_mask')
       status = parallel_def_var(NCO%id,'bmlt_float_mask',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'mask for basal melting of floating ice')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_basal_melt_rate_floating_mask')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable bmlt_float_mask was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     bmlt_ground -- basal melt rate for grounded ice
    pos = index(NCO%vars,' bmlt_ground ')
    status = parallel_inq_varid(NCO%id,'bmlt_ground',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+11) = '           '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%temper%bmlt_ground)) then
       call write_log('Creating variable bmlt_ground')
       status = parallel_def_var(NCO%id,'bmlt_ground',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_acab))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal melt rate for grounded ice')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_basal_melt_rate_grounded')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable bmlt_ground was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     btemp -- basal ice temperature
    pos = index(NCO%vars,' btemp ')
    status = parallel_inq_varid(NCO%id,'btemp',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+5) = '     '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%temper%temp)) then
       call write_log('Creating variable btemp')
       status = parallel_def_var(NCO%id,'btemp',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal ice temperature')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_temperature')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'degree_Celsius')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable btemp was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     btractx -- basal traction (x-direction comp)
    pos = index(NCO%vars,' btractx ')
    status = parallel_inq_varid(NCO%id,'btractx',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+7) = '       '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%stress%btractx)) then
       call write_log('Creating variable btractx')
       status = parallel_def_var(NCO%id,'btractx',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_tau))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal traction (x-direction comp)')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable btractx was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     btractx_extend -- basal traction (x-direction comp)
    pos = index(NCO%vars,' btractx_extend ')
    status = parallel_inq_varid(NCO%id,'btractx_extend',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+14) = '              '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%stress%btractx_extend)) then
       call write_log('Creating variable btractx_extend')
       status = parallel_def_var(NCO%id,'btractx_extend',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_tau))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal traction (x-direction comp)')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable btractx_extend was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     btracty -- basal traction (y-direction comp)
    pos = index(NCO%vars,' btracty ')
    status = parallel_inq_varid(NCO%id,'btracty',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+7) = '       '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%stress%btracty)) then
       call write_log('Creating variable btracty')
       status = parallel_def_var(NCO%id,'btracty',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_tau))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal traction (y-direction comp)')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable btracty was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     btracty_extend -- basal traction (y-direction comp)
    pos = index(NCO%vars,' btracty_extend ')
    status = parallel_inq_varid(NCO%id,'btracty_extend',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+14) = '              '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%stress%btracty_extend)) then
       call write_log('Creating variable btracty_extend')
       status = parallel_def_var(NCO%id,'btracty_extend',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_tau))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal traction (y-direction comp)')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable btracty_extend was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     btrc -- basal slip coefficient
    pos = index(NCO%vars,' btrc ')
    status = parallel_inq_varid(NCO%id,'btrc',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%btrc)) then
       call write_log('Creating variable btrc')
       status = parallel_def_var(NCO%id,'btrc',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_btrc))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal slip coefficient')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/pascal/year')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable btrc was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     bwat -- basal water depth
    pos = index(NCO%vars,' bwat ')
    status = parallel_inq_varid(NCO%id,'bwat',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%temper%bwat)) then
       call write_log('Creating variable bwat')
       status = parallel_def_var(NCO%id,'bwat',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(thk0))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal water depth')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable bwat was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     bwatflx -- basal water flux
    pos = index(NCO%vars,' bwatflx ')
    status = parallel_inq_varid(NCO%id,'bwatflx',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+7) = '       '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%temper%bwatflx)) then
       call write_log('Creating variable bwatflx')
       status = parallel_def_var(NCO%id,'bwatflx',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(thk0))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal water flux')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter3/year')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable bwatflx was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     calving -- ice margin calving
    pos = index(NCO%vars,' calving ')
    status = parallel_inq_varid(NCO%id,'calving',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+7) = '       '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%calving%calving_thck)) then
       call write_log('Creating variable calving')
       status = parallel_def_var(NCO%id,'calving',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(thk0))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'ice margin calving')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable calving was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     damage -- ice damage
    pos = index(NCO%vars,' damage ')
    status = parallel_inq_varid(NCO%id,'damage',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+6) = '      '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%calving%damage)) then
       call write_log('Creating variable damage')
       status = parallel_def_var(NCO%id,'damage',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, staglevel_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'ice damage')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'unitless [0,1]')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable damage was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     damage_column -- vertically integrated ice damage
    pos = index(NCO%vars,' damage_column ')
    status = parallel_inq_varid(NCO%id,'damage_column',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+13) = '             '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%calving%damage_column)) then
       call write_log('Creating variable damage_column')
       status = parallel_def_var(NCO%id,'damage_column',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'vertically integrated ice damage')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'unitless [0,1]')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable damage_column was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     diff_cfl_dt -- diffusive CFL maximum time step
    pos = index(NCO%vars,' diff_cfl_dt ')
    status = parallel_inq_varid(NCO%id,'diff_cfl_dt',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+11) = '           '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%numerics%diff_cfl_dt)) then
       call write_log('Creating variable diff_cfl_dt')
       status = parallel_def_var(NCO%id,'diff_cfl_dt',get_xtype(outfile,NF90_FLOAT), &
            (/time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'diffusive CFL maximum time step')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'years')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable diff_cfl_dt was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     diffu -- apparent diffusivity
    pos = index(NCO%vars,' diffu ')
    status = parallel_inq_varid(NCO%id,'diffu',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+5) = '     '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%diffu)) then
       call write_log('Creating variable diffu')
       status = parallel_def_var(NCO%id,'diffu',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_diffu))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'apparent diffusivity')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter2/year')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable diffu was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     dissip -- dissipation rate (W m-3) divided by rhoi Ci
    pos = index(NCO%vars,' dissip ')
    status = parallel_inq_varid(NCO%id,'dissip',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+6) = '      '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%temper%dissip)) then
       call write_log('Creating variable dissip')
       status = parallel_def_var(NCO%id,'dissip',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, level_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scyr))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'dissipation rate (W m-3) divided by rhoi Ci')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'deg C/yr')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable dissip was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     dissipstag -- dissipation rate (W m-3) divided by rhoi Ci
    pos = index(NCO%vars,' dissipstag ')
    status = parallel_inq_varid(NCO%id,'dissipstag',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+10) = '          '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%temper%dissip)) then
       call write_log('Creating variable dissipstag')
       status = parallel_def_var(NCO%id,'dissipstag',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, staglevel_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scyr))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'dissipation rate (W m-3) divided by rhoi Ci')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'deg C/yr')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable dissipstag was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     dthckdtm -- tendency of ice thickness (NOTE: Glide only)
    pos = index(NCO%vars,' dthckdtm ')
    status = parallel_inq_varid(NCO%id,'dthckdtm',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+8) = '        '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geomderv%dthckdtm)) then
       call write_log('Creating variable dthckdtm')
       status = parallel_def_var(NCO%id,'dthckdtm',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_acab))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'tendency of ice thickness (NOTE: Glide only)')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable dthckdtm was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     dusrfdtm -- rate of upper ice surface elevation change (NOTE: Glide only)
    pos = index(NCO%vars,' dusrfdtm ')
    status = parallel_inq_varid(NCO%id,'dusrfdtm',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+8) = '        '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geomderv%dusrfdtm)) then
       call write_log('Creating variable dusrfdtm')
       status = parallel_def_var(NCO%id,'dusrfdtm',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_acab))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'rate of upper ice surface elevation change (NOTE: Glide only)')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable dusrfdtm was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     dynbcmask -- 2d array of higher-order model boundary condition mask values (NOTE: Glam ONLY)
    pos = index(NCO%vars,' dynbcmask ')
    status = parallel_inq_varid(NCO%id,'dynbcmask',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+9) = '         '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%dynbcmask)) then
       call write_log('Creating variable dynbcmask')
       status = parallel_def_var(NCO%id,'dynbcmask',get_xtype(outfile,NF90_INT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            '2d array of higher-order model boundary condition mask values (NOTE: Glam ONLY)')
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable dynbcmask was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     effecpress -- effective pressure
    pos = index(NCO%vars,' effecpress ')
    status = parallel_inq_varid(NCO%id,'effecpress',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+10) = '          '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%basal_physics%effecpress)) then
       call write_log('Creating variable effecpress')
       status = parallel_def_var(NCO%id,'effecpress',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'effective pressure')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable effecpress was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     efvs -- effective viscosity
    pos = index(NCO%vars,' efvs ')
    status = parallel_inq_varid(NCO%id,'efvs',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%stress%efvs)) then
       call write_log('Creating variable efvs')
       status = parallel_def_var(NCO%id,'efvs',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, staglevel_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_efvs))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'effective viscosity')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pascal * years')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable efvs was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     enthalpy -- specific enthalpy
    pos = index(NCO%vars,' enthalpy ')
    status = parallel_inq_varid(NCO%id,'enthalpy',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+8) = '        '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%temper%enthalpy)) then
       call write_log('Creating variable enthalpy')
       status = parallel_def_var(NCO%id,'enthalpy',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, stagwbndlevel_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'specific enthalpy')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'J/m^3')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable enthalpy was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     eus -- global average sea level
    pos = index(NCO%vars,' eus ')
    status = parallel_inq_varid(NCO%id,'eus',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+3) = '   '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%climate%eus)) then
       call write_log('Creating variable eus')
       status = parallel_def_var(NCO%id,'eus',get_xtype(outfile,NF90_FLOAT), &
            (/time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(thk0))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'global average sea level')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'global_average_sea_level_change')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable eus was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     f_flotation -- flotation function
    pos = index(NCO%vars,' f_flotation ')
    status = parallel_inq_varid(NCO%id,'f_flotation',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+11) = '           '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%f_flotation)) then
       call write_log('Creating variable f_flotation')
       status = parallel_def_var(NCO%id,'f_flotation',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'flotation function')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'unitless')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable f_flotation was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     f_ground -- grounded ice fraction
    pos = index(NCO%vars,' f_ground ')
    status = parallel_inq_varid(NCO%id,'f_ground',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+8) = '        '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%f_ground)) then
       call write_log('Creating variable f_ground')
       status = parallel_def_var(NCO%id,'f_ground',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'grounded ice fraction')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'grounded_fraction')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'unitless [0,1]')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable f_ground was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     flwa -- Pre-exponential flow law parameter
    pos = index(NCO%vars,' flwa ')
    status = parallel_inq_varid(NCO%id,'flwa',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%temper%flwa)) then
       call write_log('Creating variable flwa')
       status = parallel_def_var(NCO%id,'flwa',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, level_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_flwa))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'Pre-exponential flow law parameter')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'pascal**(-n) year**(-1)')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable flwa was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     flwastag -- Pre-exponential flow law parameter
    pos = index(NCO%vars,' flwastag ')
    status = parallel_inq_varid(NCO%id,'flwastag',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+8) = '        '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%temper%flwa)) then
       call write_log('Creating variable flwastag')
       status = parallel_def_var(NCO%id,'flwastag',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, staglevel_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_flwa))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'Pre-exponential flow law parameter')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'pascal**(-n) year**(-1)')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable flwastag was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     gravity -- gravitational acceleration
    pos = index(NCO%vars,' gravity ')
    status = parallel_inq_varid(NCO%id,'gravity',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+7) = '       '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(grav)) then
       call write_log('Creating variable gravity')
       status = parallel_def_var(NCO%id,'gravity',get_xtype(outfile,NF90_FLOAT), &
            (/time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',1.0)
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'gravitational acceleration')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'gravity')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/s/s')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable gravity was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     iarea -- area covered by ice
    pos = index(NCO%vars,' iarea ')
    status = parallel_inq_varid(NCO%id,'iarea',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+5) = '     '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%iarea)) then
       call write_log('Creating variable iarea')
       status = parallel_def_var(NCO%id,'iarea',get_xtype(outfile,NF90_FLOAT), &
            (/time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(len0*len0*1.e-6))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'area covered by ice')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'km2')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable iarea was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     iareaf -- area covered by floating ice
    pos = index(NCO%vars,' iareaf ')
    status = parallel_inq_varid(NCO%id,'iareaf',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+6) = '      '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%iareaf)) then
       call write_log('Creating variable iareaf')
       status = parallel_def_var(NCO%id,'iareaf',get_xtype(outfile,NF90_FLOAT), &
            (/time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(len0*len0*1.e-6))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'area covered by floating ice')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'km2')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable iareaf was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     iareag -- area covered by grounded ice
    pos = index(NCO%vars,' iareag ')
    status = parallel_inq_varid(NCO%id,'iareag',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+6) = '      '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%iareag)) then
       call write_log('Creating variable iareag')
       status = parallel_def_var(NCO%id,'iareag',get_xtype(outfile,NF90_FLOAT), &
            (/time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(len0*len0*1.e-6))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'area covered by grounded ice')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'km2')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable iareag was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     ice_age -- ice age
    pos = index(NCO%vars,' ice_age ')
    status = parallel_inq_varid(NCO%id,'ice_age',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+7) = '       '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%ice_age)) then
       call write_log('Creating variable ice_age')
       status = parallel_def_var(NCO%id,'ice_age',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, staglevel_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(tim0/scyr))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'ice age')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_age')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'year')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable ice_age was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     ice_mask -- real-valued mask denoting ice (1) or no ice (0)
    pos = index(NCO%vars,' ice_mask ')
    status = parallel_inq_varid(NCO%id,'ice_mask',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+8) = '        '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%ice_mask)) then
       call write_log('Creating variable ice_mask')
       status = parallel_def_var(NCO%id,'ice_mask',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(1.0))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'real-valued mask denoting ice (1) or no ice (0)')
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable ice_mask was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     ice_specific_heat -- ice specific heat
    pos = index(NCO%vars,' ice_specific_heat ')
    status = parallel_inq_varid(NCO%id,'ice_specific_heat',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+17) = '                 '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(shci)) then
       call write_log('Creating variable ice_specific_heat')
       status = parallel_def_var(NCO%id,'ice_specific_heat',get_xtype(outfile,NF90_FLOAT), &
            (/time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',1.0)
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'ice specific heat')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'ice_specific_heat')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'J/kg/K')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable ice_specific_heat was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     ice_thermal_conductivity -- ice thermal conductivity
    pos = index(NCO%vars,' ice_thermal_conductivity ')
    status = parallel_inq_varid(NCO%id,'ice_thermal_conductivity',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+24) = '                        '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(coni)) then
       call write_log('Creating variable ice_thermal_conductivity')
       status = parallel_def_var(NCO%id,'ice_thermal_conductivity',get_xtype(outfile,NF90_FLOAT), &
            (/time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',1.0)
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'ice thermal conductivity')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'ice_thermal_conductivity')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'J/(K kg)')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable ice_thermal_conductivity was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     ivol -- ice volume
    pos = index(NCO%vars,' ivol ')
    status = parallel_inq_varid(NCO%id,'ivol',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%ivol)) then
       call write_log('Creating variable ivol')
       status = parallel_def_var(NCO%id,'ivol',get_xtype(outfile,NF90_FLOAT), &
            (/time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(thk0*len0*len0*1.e-9))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'ice volume')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'km3')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable ivol was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     kinbcmask -- Mask of locations where uvel, vvel value should be held constant
    pos = index(NCO%vars,' kinbcmask ')
    status = parallel_inq_varid(NCO%id,'kinbcmask',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+9) = '         '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%kinbcmask)) then
       call write_log('Creating variable kinbcmask')
       status = parallel_def_var(NCO%id,'kinbcmask',get_xtype(outfile,NF90_INT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'Mask of locations where uvel, vvel value should be held constant')
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable kinbcmask was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     litho_temp -- lithosphere temperature
    pos = index(NCO%vars,' litho_temp ')
    status = parallel_inq_varid(NCO%id,'litho_temp',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+10) = '          '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%lithot%temp)) then
       call write_log('Creating variable litho_temp')
       status = parallel_def_var(NCO%id,'litho_temp',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, lithoz_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'lithosphere temperature')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'degree_Celsius')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable litho_temp was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     lsurf -- ice lower surface elevation
    pos = index(NCO%vars,' lsurf ')
    status = parallel_inq_varid(NCO%id,'lsurf',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+5) = '     '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%lsrf)) then
       call write_log('Creating variable lsurf')
       status = parallel_def_var(NCO%id,'lsurf',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(thk0))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'ice lower surface elevation')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable lsurf was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     relx -- relaxed bedrock topography
    pos = index(NCO%vars,' relx ')
    status = parallel_inq_varid(NCO%id,'relx',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%isostasy%relx)) then
       call write_log('Creating variable relx')
       status = parallel_def_var(NCO%id,'relx',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(thk0))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'relaxed bedrock topography')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable relx was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     resid_u -- u component of residual Ax - b (NOTE: Glam only)
    pos = index(NCO%vars,' resid_u ')
    status = parallel_inq_varid(NCO%id,'resid_u',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+7) = '       '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%resid_u)) then
       call write_log('Creating variable resid_u')
       status = parallel_def_var(NCO%id,'resid_u',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, level_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_resid))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'u component of residual Ax - b (NOTE: Glam only)')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa/m')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable resid_u was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     resid_v -- v component of residual Ax - b (NOTE: Glam only)
    pos = index(NCO%vars,' resid_v ')
    status = parallel_inq_varid(NCO%id,'resid_v',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+7) = '       '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%resid_v)) then
       call write_log('Creating variable resid_v')
       status = parallel_def_var(NCO%id,'resid_v',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, level_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_resid))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'v component of residual Ax - b (NOTE: Glam only)')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa/m')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable resid_v was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     rho_ice -- ice density
    pos = index(NCO%vars,' rho_ice ')
    status = parallel_inq_varid(NCO%id,'rho_ice',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+7) = '       '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(rhoi)) then
       call write_log('Creating variable rho_ice')
       status = parallel_def_var(NCO%id,'rho_ice',get_xtype(outfile,NF90_FLOAT), &
            (/time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',1.0)
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'ice density')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'rho_ice')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'kg/meter3')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable rho_ice was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     rho_seawater -- seawater density
    pos = index(NCO%vars,' rho_seawater ')
    status = parallel_inq_varid(NCO%id,'rho_seawater',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+12) = '            '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(rhoo)) then
       call write_log('Creating variable rho_seawater')
       status = parallel_def_var(NCO%id,'rho_seawater',get_xtype(outfile,NF90_FLOAT), &
            (/time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',1.0)
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'seawater density')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'rho_seawater')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'kg/meter3')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable rho_seawater was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     rhs_u -- u component of b in Ax = b
    pos = index(NCO%vars,' rhs_u ')
    status = parallel_inq_varid(NCO%id,'rhs_u',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+5) = '     '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%rhs_u)) then
       call write_log('Creating variable rhs_u')
       status = parallel_def_var(NCO%id,'rhs_u',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, level_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_resid))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'u component of b in Ax = b')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa/m')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable rhs_u was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     rhs_v -- v component of b in Ax = b
    pos = index(NCO%vars,' rhs_v ')
    status = parallel_inq_varid(NCO%id,'rhs_v',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+5) = '     '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%rhs_v)) then
       call write_log('Creating variable rhs_v')
       status = parallel_def_var(NCO%id,'rhs_v',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, level_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_resid))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'v component of b in Ax = b')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa/m')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable rhs_v was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     seconds_per_year -- seconds per year
    pos = index(NCO%vars,' seconds_per_year ')
    status = parallel_inq_varid(NCO%id,'seconds_per_year',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+16) = '                '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(scyr)) then
       call write_log('Creating variable seconds_per_year')
       status = parallel_def_var(NCO%id,'seconds_per_year',get_xtype(outfile,NF90_FLOAT), &
            (/time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',1.0)
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'seconds per year')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'seconds_per_year')
       status = parallel_put_att(NCO%id, varid, 'units', &
            's/yr')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable seconds_per_year was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     soft -- bed softness parameter
    pos = index(NCO%vars,' soft ')
    status = parallel_inq_varid(NCO%id,'soft',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%bed_softness)) then
       call write_log('Creating variable soft')
       status = parallel_def_var(NCO%id,'soft',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_btrc))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'bed softness parameter')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/pascal/year')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable soft was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     stagthk -- staggered ice thickness
    pos = index(NCO%vars,' stagthk ')
    status = parallel_inq_varid(NCO%id,'stagthk',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+7) = '       '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geomderv%stagthck)) then
       call write_log('Creating variable stagthk')
       status = parallel_def_var(NCO%id,'stagthk',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(thk0))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'staggered ice thickness')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'stag_land_ice_thickness')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable stagthk was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     surftemp -- annual mean surface temperature
    pos = index(NCO%vars,' surftemp ')
    status = parallel_inq_varid(NCO%id,'surftemp',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+8) = '        '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%climate%artm)) then
       call write_log('Creating variable surftemp')
       status = parallel_def_var(NCO%id,'surftemp',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'annual mean surface temperature')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'surface_temperature')
       status = parallel_put_att(NCO%id, varid, 'cell_methods', &
            'time: mean')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'degree_Celsius')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable surftemp was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     tau_eff -- effective stress
    pos = index(NCO%vars,' tau_eff ')
    status = parallel_inq_varid(NCO%id,'tau_eff',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+7) = '       '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%stress%tau%scalar)) then
       call write_log('Creating variable tau_eff')
       status = parallel_def_var(NCO%id,'tau_eff',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, staglevel_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_tau))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'effective stress')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable tau_eff was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     tau_xx -- x component of horiz. normal stress
    pos = index(NCO%vars,' tau_xx ')
    status = parallel_inq_varid(NCO%id,'tau_xx',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+6) = '      '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%stress%tau%xx)) then
       call write_log('Creating variable tau_xx')
       status = parallel_def_var(NCO%id,'tau_xx',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, staglevel_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_tau))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'x component of horiz. normal stress')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable tau_xx was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     tau_xy -- horiz. shear stress
    pos = index(NCO%vars,' tau_xy ')
    status = parallel_inq_varid(NCO%id,'tau_xy',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+6) = '      '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%stress%tau%xy)) then
       call write_log('Creating variable tau_xy')
       status = parallel_def_var(NCO%id,'tau_xy',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, staglevel_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_tau))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'horiz. shear stress')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable tau_xy was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     tau_xz -- X component vertical shear stress
    pos = index(NCO%vars,' tau_xz ')
    status = parallel_inq_varid(NCO%id,'tau_xz',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+6) = '      '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%stress%tau%xz)) then
       call write_log('Creating variable tau_xz')
       status = parallel_def_var(NCO%id,'tau_xz',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, staglevel_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_tau))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'X component vertical shear stress')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable tau_xz was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     tau_yy -- y component of horiz. normal stress
    pos = index(NCO%vars,' tau_yy ')
    status = parallel_inq_varid(NCO%id,'tau_yy',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+6) = '      '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%stress%tau%yy)) then
       call write_log('Creating variable tau_yy')
       status = parallel_def_var(NCO%id,'tau_yy',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, staglevel_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_tau))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'y component of horiz. normal stress')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable tau_yy was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     tau_yz -- Y component vertical shear stress
    pos = index(NCO%vars,' tau_yz ')
    status = parallel_inq_varid(NCO%id,'tau_yz',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+6) = '      '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%stress%tau%yz)) then
       call write_log('Creating variable tau_yz')
       status = parallel_def_var(NCO%id,'tau_yz',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, staglevel_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_tau))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'Y component vertical shear stress')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable tau_yz was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     taudx -- driving stress (x-direction comp)
    pos = index(NCO%vars,' taudx ')
    status = parallel_inq_varid(NCO%id,'taudx',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+5) = '     '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%stress%taudx)) then
       call write_log('Creating variable taudx')
       status = parallel_def_var(NCO%id,'taudx',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_tau))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'driving stress (x-direction comp)')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable taudx was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     taudy -- driving stress (y-direction comp)
    pos = index(NCO%vars,' taudy ')
    status = parallel_inq_varid(NCO%id,'taudy',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+5) = '     '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%stress%taudy)) then
       call write_log('Creating variable taudy')
       status = parallel_def_var(NCO%id,'taudy',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_tau))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'driving stress (y-direction comp)')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable taudy was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     tauf -- higher-order basal yield stress
    pos = index(NCO%vars,' tauf ')
    status = parallel_inq_varid(NCO%id,'tauf',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%basalproc%mintauf)) then
       call write_log('Creating variable tauf')
       status = parallel_def_var(NCO%id,'tauf',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_tau))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'higher-order basal yield stress')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable tauf was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     taux -- basal shear stress in x direction (NOTE: Glide only)
    pos = index(NCO%vars,' taux ')
    status = parallel_inq_varid(NCO%id,'taux',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%tau_x)) then
       call write_log('Creating variable taux')
       status = parallel_def_var(NCO%id,'taux',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(1e-3*thk0*thk0/len0))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal shear stress in x direction (NOTE: Glide only)')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'kilopascal')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable taux was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     tauy -- basal shear stress in y direction
    pos = index(NCO%vars,' tauy ')
    status = parallel_inq_varid(NCO%id,'tauy',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%tau_y)) then
       call write_log('Creating variable tauy')
       status = parallel_def_var(NCO%id,'tauy',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(1e-3*thk0*thk0/len0))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal shear stress in y direction')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'kilopascal')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable tauy was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     temp -- ice temperature
    pos = index(NCO%vars,' temp ')
    status = parallel_inq_varid(NCO%id,'temp',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%temper%temp)) then
       call write_log('Creating variable temp')
       status = parallel_def_var(NCO%id,'temp',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, level_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'ice temperature')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_temperature')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'degree_Celsius')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable temp was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     tempstag -- ice temperature on staggered vertical levels with boundaries
    pos = index(NCO%vars,' tempstag ')
    status = parallel_inq_varid(NCO%id,'tempstag',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+8) = '        '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%temper%temp)) then
       call write_log('Creating variable tempstag')
       status = parallel_def_var(NCO%id,'tempstag',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, stagwbndlevel_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'ice temperature on staggered vertical levels with boundaries')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_temperature_stag')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'degree_Celsius')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable tempstag was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     thk -- ice thickness
    pos = index(NCO%vars,' thk ')
    status = parallel_inq_varid(NCO%id,'thk',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+3) = '   '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%thck)) then
       call write_log('Creating variable thk')
       status = parallel_def_var(NCO%id,'thk',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(thk0))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'ice thickness')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_thickness')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable thk was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     thkmask -- mask
    pos = index(NCO%vars,' thkmask ')
    status = parallel_inq_varid(NCO%id,'thkmask',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+7) = '       '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%thkmask)) then
       call write_log('Creating variable thkmask')
       status = parallel_def_var(NCO%id,'thkmask',get_xtype(outfile,NF90_INT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'mask')
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable thkmask was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     topg -- bedrock topography
    pos = index(NCO%vars,' topg ')
    status = parallel_inq_varid(NCO%id,'topg',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%topg)) then
       call write_log('Creating variable topg')
       status = parallel_def_var(NCO%id,'topg',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(thk0))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'bedrock topography')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'bedrock_altitude')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable topg was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     ubas -- basal slip velocity in x direction
    pos = index(NCO%vars,' ubas ')
    status = parallel_inq_varid(NCO%id,'ubas',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%ubas)) then
       call write_log('Creating variable ubas')
       status = parallel_def_var(NCO%id,'ubas',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_uvel))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal slip velocity in x direction')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_basal_x_velocity')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable ubas was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     uflx -- flux in x direction (NOTE: Glide and Glam only)
    pos = index(NCO%vars,' uflx ')
    status = parallel_inq_varid(NCO%id,'uflx',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%uflx)) then
       call write_log('Creating variable uflx')
       status = parallel_def_var(NCO%id,'uflx',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_uflx))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'flux in x direction (NOTE: Glide and Glam only)')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter2/year')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable uflx was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     unstagbeta -- higher-order bed stress coefficient on the unstaggered grid (NOTE: this will overwrite beta if both are input)
    pos = index(NCO%vars,' unstagbeta ')
    status = parallel_inq_varid(NCO%id,'unstagbeta',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+10) = '          '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%unstagbeta)) then
       call write_log('Creating variable unstagbeta')
       status = parallel_def_var(NCO%id,'unstagbeta',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_beta))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'higher-order bed stress coefficient on the unstaggered grid (NOTE: this will overwrite beta if both are input)')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa yr/m')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable unstagbeta was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     usurf -- ice upper surface elevation
    pos = index(NCO%vars,' usurf ')
    status = parallel_inq_varid(NCO%id,'usurf',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+5) = '     '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%usrf)) then
       call write_log('Creating variable usurf')
       status = parallel_def_var(NCO%id,'usurf',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(thk0))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'ice upper surface elevation')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'surface_altitude')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable usurf was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     uvel -- ice velocity in x direction
    pos = index(NCO%vars,' uvel ')
    status = parallel_inq_varid(NCO%id,'uvel',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%uvel)) then
       call write_log('Creating variable uvel')
       status = parallel_def_var(NCO%id,'uvel',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, level_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_uvel))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'ice velocity in x direction')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_x_velocity')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable uvel was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     uvel_2d -- vertically averaged ice velocity in x direction
    pos = index(NCO%vars,' uvel_2d ')
    status = parallel_inq_varid(NCO%id,'uvel_2d',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+7) = '       '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%uvel_2d)) then
       call write_log('Creating variable uvel_2d')
       status = parallel_def_var(NCO%id,'uvel_2d',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_uvel))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'vertically averaged ice velocity in x direction')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_x_velocity_2d')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable uvel_2d was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     uvel_2d_extend -- vertically averaged ice velocity in x direction (extended grid)
    pos = index(NCO%vars,' uvel_2d_extend ')
    status = parallel_inq_varid(NCO%id,'uvel_2d_extend',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+14) = '              '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%uvel_2d_extend)) then
       call write_log('Creating variable uvel_2d_extend')
       status = parallel_def_var(NCO%id,'uvel_2d_extend',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_uvel))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'vertically averaged ice velocity in x direction (extended grid)')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_x_velocity_2d')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable uvel_2d_extend was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     uvel_extend -- ice velocity in x direction (extended grid)
    pos = index(NCO%vars,' uvel_extend ')
    status = parallel_inq_varid(NCO%id,'uvel_extend',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+11) = '           '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%uvel_extend)) then
       call write_log('Creating variable uvel_extend')
       status = parallel_def_var(NCO%id,'uvel_extend',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, level_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_uvel))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'ice velocity in x direction (extended grid)')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_x_velocity')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable uvel_extend was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     vbas -- basal slip velocity in y direction
    pos = index(NCO%vars,' vbas ')
    status = parallel_inq_varid(NCO%id,'vbas',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%vbas)) then
       call write_log('Creating variable vbas')
       status = parallel_def_var(NCO%id,'vbas',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_uvel))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal slip velocity in y direction')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_basal_y_velocity')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable vbas was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     velnorm -- Horizontal ice velocity magnitude
    pos = index(NCO%vars,' velnorm ')
    status = parallel_inq_varid(NCO%id,'velnorm',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+7) = '       '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%velnorm)) then
       call write_log('Creating variable velnorm')
       status = parallel_def_var(NCO%id,'velnorm',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, level_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_uvel))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'Horizontal ice velocity magnitude')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable velnorm was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     vflx -- flux in x direction (NOTE: Glide and Glam only)
    pos = index(NCO%vars,' vflx ')
    status = parallel_inq_varid(NCO%id,'vflx',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%vflx)) then
       call write_log('Creating variable vflx')
       status = parallel_def_var(NCO%id,'vflx',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_uflx))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'flux in x direction (NOTE: Glide and Glam only)')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter2/year')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable vflx was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     vvel -- ice velocity in y direction
    pos = index(NCO%vars,' vvel ')
    status = parallel_inq_varid(NCO%id,'vvel',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%vvel)) then
       call write_log('Creating variable vvel')
       status = parallel_def_var(NCO%id,'vvel',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, level_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_uvel))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'ice velocity in y direction')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_y_velocity')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable vvel was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     vvel_2d -- vertically averaged ice velocity in y direction
    pos = index(NCO%vars,' vvel_2d ')
    status = parallel_inq_varid(NCO%id,'vvel_2d',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+7) = '       '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%vvel_2d)) then
       call write_log('Creating variable vvel_2d')
       status = parallel_def_var(NCO%id,'vvel_2d',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_uvel))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'vertically averaged ice velocity in y direction')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_y_velocity_2d')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable vvel_2d was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     vvel_2d_extend -- vertically averaged ice velocity in y direction (extended grid)
    pos = index(NCO%vars,' vvel_2d_extend ')
    status = parallel_inq_varid(NCO%id,'vvel_2d_extend',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+14) = '              '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%vvel_2d_extend)) then
       call write_log('Creating variable vvel_2d_extend')
       status = parallel_def_var(NCO%id,'vvel_2d_extend',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_uvel))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'vertically averaged ice velocity in y direction (extended grid)')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_y_velocity_2d')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable vvel_2d_extend was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     vvel_extend -- ice velocity in y direction (extended grid)
    pos = index(NCO%vars,' vvel_extend ')
    status = parallel_inq_varid(NCO%id,'vvel_extend',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+11) = '           '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%vvel_extend)) then
       call write_log('Creating variable vvel_extend')
       status = parallel_def_var(NCO%id,'vvel_extend',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, level_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_uvel))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'ice velocity in y direction (extended grid)')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_y_velocity')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable vvel_extend was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     waterfrac -- internal water fraction
    pos = index(NCO%vars,' waterfrac ')
    status = parallel_inq_varid(NCO%id,'waterfrac',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+9) = '         '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%temper%waterfrac)) then
       call write_log('Creating variable waterfrac')
       status = parallel_def_var(NCO%id,'waterfrac',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, staglevel_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'internal water fraction')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'unitless [0,1]')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable waterfrac was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     wgrd -- Vertical grid velocity
    pos = index(NCO%vars,' wgrd ')
    status = parallel_inq_varid(NCO%id,'wgrd',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%wgrd)) then
       call write_log('Creating variable wgrd')
       status = parallel_def_var(NCO%id,'wgrd',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, level_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_wvel))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'Vertical grid velocity')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable wgrd was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     wvel -- vertical ice velocity
    pos = index(NCO%vars,' wvel ')
    status = parallel_inq_varid(NCO%id,'wvel',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%wvel)) then
       call write_log('Creating variable wvel')
       status = parallel_def_var(NCO%id,'wvel',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, level_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_wvel))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'vertical ice velocity')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable wvel was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     wvel_ho -- vertical ice velocity relative to ice sheet base from higher-order model (NOTE: Glam only)
    pos = index(NCO%vars,' wvel_ho ')
    status = parallel_inq_varid(NCO%id,'wvel_ho',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+7) = '       '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%wvel_ho)) then
       call write_log('Creating variable wvel_ho')
       status = parallel_def_var(NCO%id,'wvel_ho',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, level_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_wvel))
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'vertical ice velocity relative to ice sheet base from higher-order model (NOTE: Glam only)')
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
          status = parallel_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       end if
     else
     call write_log('Variable wvel_ho was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

  end subroutine glide_io_create

  subroutine glide_io_write(outfile,data)
    use parallel
    use glide_types
    use glimmer_ncdf
    use glimmer_paramets
    use glimmer_scales
    implicit none
    type(glimmer_nc_output), pointer :: outfile
    ! structure containg output netCDF descriptor
    type(glide_global_type) :: data
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
    status = parallel_inq_varid(NCO%id,'C_space_factor',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%basal_physics%C_space_factor, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'acab',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%climate%acab, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'adv_cfl_dt',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%numerics%adv_cfl_dt, (/outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'artm',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%climate%artm, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'beta',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%beta, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'beta_internal',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%beta_internal, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'bfricflx',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%temper%bfricflx, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'bheatflx',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%temper%bheatflx, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'bmlt_float',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%temper%bmlt_float, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'bmlt_float_mask',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%temper%bmlt_float_mask, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'bmlt_ground',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%temper%bmlt_ground, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'btemp',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%temper%temp(data%general%upn,1:data%general%ewn,1:data%general%nsn), (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'btractx',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%stress%btractx(:,:), (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'btractx_extend',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%stress%btractx_extend(:,:), (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'btracty',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%stress%btracty(:,:), (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'btracty_extend',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%stress%btracty_extend(:,:), (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'btrc',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%btrc, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'bwat',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%temper%bwat, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'bwatflx',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%temper%bwatflx, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'calving',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%calving%calving_thck, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'damage',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nstaglevel
          status = distributed_put_var(NCO%id, varid, &
               data%calving%damage(up,:,:), (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'damage_column',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%calving%damage_column, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'diff_cfl_dt',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%numerics%diff_cfl_dt, (/outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'diffu',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%diffu, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'dissip',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nlevel
          status = distributed_put_var(NCO%id, varid, &
               data%temper%dissip(up,:,:), (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'dissipstag',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nstaglevel
          status = distributed_put_var(NCO%id, varid, &
               data%temper%dissip(up,:,:), (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'dthckdtm',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geomderv%dthckdtm, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'dusrfdtm',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geomderv%dusrfdtm, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'dynbcmask',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%dynbcmask, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'effecpress',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%basal_physics%effecpress, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'efvs',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nstaglevel
          status = distributed_put_var(NCO%id, varid, &
               data%stress%efvs(up,:,:), (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'enthalpy',varid)
    if (status .eq. nf90_noerr) then
       do up=0,NCO%nstagwbndlevel
          status = distributed_put_var(NCO%id, varid, &
               data%temper%enthalpy(up,:,:), (/1,1,up+1,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'eus',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%climate%eus, (/outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'f_flotation',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%f_flotation, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'f_ground',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%f_ground, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'flwa',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nlevel
          status = distributed_put_var(NCO%id, varid, &
               data%temper%flwa(up,:,:), (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'flwastag',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nstaglevel
          status = distributed_put_var(NCO%id, varid, &
               data%temper%flwa(up,:,:), (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'gravity',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            grav, (/outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'iarea',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%iarea, (/outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'iareaf',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%iareaf, (/outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'iareag',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%iareag, (/outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'ice_age',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nstaglevel
          status = distributed_put_var(NCO%id, varid, &
               data%geometry%ice_age(up,:,:), (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'ice_mask',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%ice_mask, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'ice_specific_heat',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            shci, (/outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'ice_thermal_conductivity',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            coni, (/outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'ivol',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%ivol, (/outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'kinbcmask',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%kinbcmask(:,:), (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'litho_temp',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%lithot%temp, (/1,1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'lsurf',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%lsrf, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'relx',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%isostasy%relx, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'resid_u',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nlevel
          status = distributed_put_var(NCO%id, varid, &
               data%velocity%resid_u(up,:,:), (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'resid_v',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nlevel
          status = distributed_put_var(NCO%id, varid, &
               data%velocity%resid_v(up,:,:), (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'rho_ice',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            rhoi, (/outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'rho_seawater',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            rhoo, (/outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'rhs_u',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nlevel
          status = distributed_put_var(NCO%id, varid, &
               data%velocity%rhs_u(up,:,:), (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'rhs_v',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nlevel
          status = distributed_put_var(NCO%id, varid, &
               data%velocity%rhs_v(up,:,:), (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'seconds_per_year',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            scyr, (/outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'soft',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%bed_softness, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'stagthk',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geomderv%stagthck, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'surftemp',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%climate%artm, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'tau_eff',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nstaglevel
          status = distributed_put_var(NCO%id, varid, &
               data%stress%tau%scalar(up,:,:), (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'tau_xx',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nstaglevel
          status = distributed_put_var(NCO%id, varid, &
               data%stress%tau%xx(up,:,:), (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'tau_xy',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nstaglevel
          status = distributed_put_var(NCO%id, varid, &
               data%stress%tau%xy(up,:,:), (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'tau_xz',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nstaglevel
          status = distributed_put_var(NCO%id, varid, &
               data%stress%tau%xz(up,:,:), (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'tau_yy',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nstaglevel
          status = distributed_put_var(NCO%id, varid, &
               data%stress%tau%yy(up,:,:), (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'tau_yz',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nstaglevel
          status = distributed_put_var(NCO%id, varid, &
               data%stress%tau%yz(up,:,:), (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'taudx',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%stress%taudx(:,:), (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'taudy',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%stress%taudy(:,:), (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'tauf',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%basalproc%mintauf, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'taux',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%tau_x, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'tauy',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%tau_y, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'temp',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nlevel
          status = distributed_put_var(NCO%id, varid, &
               data%temper%temp(up,1:data%general%ewn,1:data%general%nsn), (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'tempstag',varid)
    if (status .eq. nf90_noerr) then
       do up=0,NCO%nstagwbndlevel
          status = distributed_put_var(NCO%id, varid, &
               data%temper%temp(up,1:data%general%ewn,1:data%general%nsn), (/1,1,up+1,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'thk',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%thck, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'thkmask',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%thkmask, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'topg',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%topg, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'ubas',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%ubas, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'uflx',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%uflx, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'unstagbeta',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%unstagbeta, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'usurf',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%usrf, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'uvel',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nlevel
          status = distributed_put_var(NCO%id, varid, &
               data%velocity%uvel(up,:,:), (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'uvel_2d',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%uvel_2d(:,:), (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'uvel_2d_extend',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%uvel_2d_extend(:,:), (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'uvel_extend',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nlevel
          status = distributed_put_var(NCO%id, varid, &
               data%velocity%uvel_extend(up,:,:), (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'vbas',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%vbas, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'velnorm',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nlevel
          status = distributed_put_var(NCO%id, varid, &
               data%velocity%velnorm(up,:,:), (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'vflx',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%vflx, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'vvel',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nlevel
          status = distributed_put_var(NCO%id, varid, &
               data%velocity%vvel(up,:,:), (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'vvel_2d',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%vvel_2d(:,:), (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'vvel_2d_extend',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%vvel_2d_extend(:,:), (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'vvel_extend',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nlevel
          status = distributed_put_var(NCO%id, varid, &
               data%velocity%vvel_extend(up,:,:), (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'waterfrac',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nstaglevel
          status = distributed_put_var(NCO%id, varid, &
               data%temper%waterfrac(up,:,:), (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'wgrd',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nlevel
          status = distributed_put_var(NCO%id, varid, &
               data%velocity%wgrd(up,:,:), (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'wvel',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nlevel
          status = distributed_put_var(NCO%id, varid, &
               data%velocity%wvel(up,:,:), (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'wvel_ho',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nlevel
          status = distributed_put_var(NCO%id, varid, &
               data%velocity%wvel_ho(up,:,:), (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

  end subroutine glide_io_write


  subroutine glide_add_to_restart_variable_list(vars_to_add)
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
       !call write_log('Adding to glide restart variable list: ' // trim(vars_to_add) )
    endif

  end subroutine glide_add_to_restart_variable_list


  ! Functions for the interface 'is_enabled'.  These are needed by the auto-generated code in glide_io_create
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
  subroutine glide_io_readall(data, model, filetype)
    ! read from netCDF file
    use glide_types
    use glide_types
    use glimmer_ncdf
    use glimmer_ncio
    implicit none
    type(glide_global_type) :: data
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
          call glide_io_read(ic,data)
       end if
       ic=>ic%next
    end do
  end subroutine glide_io_readall


  subroutine glide_read_forcing(data, model)
    ! Read data from forcing files
    use glimmer_log
    use glide_types
    use glimmer_ncdf

    implicit none
    type(glide_global_type) :: data
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
    call glide_io_readall(data, model, filetype=1)

  end subroutine glide_read_forcing


!------------------------------------------------------------------------------


  subroutine glide_io_read(infile,data)
    ! read variables from a netCDF file
    use parallel
    use glimmer_log
    use glimmer_ncdf
    use glide_types
    use glimmer_paramets
    use glimmer_scales
    implicit none
    type(glimmer_nc_input), pointer :: infile
    ! structure containg output netCDF descriptor
    type(glide_global_type) :: data
    ! the model instance

    ! local variables
    integer status,varid
    integer up
    real(dp) :: scaling_factor

    ! read variables
    status = parallel_inq_varid(NCI%id,'x1',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%general%x1)) then
       call write_log('  Loading x1')
       status = distributed_get_var(NCI%id, varid, &
            data%general%x1, (/1/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling x1",GM_DIAGNOSTIC)
          data%general%x1 = &
               data%general%x1*scaling_factor
       end if
    else
    call write_log('Variable x1 was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'y1',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%general%y1)) then
       call write_log('  Loading y1')
       status = distributed_get_var(NCI%id, varid, &
            data%general%y1, (/1/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling y1",GM_DIAGNOSTIC)
          data%general%y1 = &
               data%general%y1*scaling_factor
       end if
    else
    call write_log('Variable y1 was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'C_space_factor',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%basal_physics%C_space_factor)) then
       call write_log('  Loading C_space_factor')
       status = distributed_get_var(NCI%id, varid, &
            data%basal_physics%C_space_factor, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling C_space_factor",GM_DIAGNOSTIC)
          data%basal_physics%C_space_factor = &
               data%basal_physics%C_space_factor*scaling_factor
       end if
    else
    call write_log('Variable C_space_factor was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'acab',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%climate%acab)) then
       call write_log('  Loading acab')
       status = distributed_get_var(NCI%id, varid, &
            data%climate%acab, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_acab)
       else
          scaling_factor = scaling_factor/(scale_acab)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling acab",GM_DIAGNOSTIC)
          data%climate%acab = &
               data%climate%acab*scaling_factor
       end if
    else
    call write_log('Variable acab was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'artm',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%climate%artm)) then
       call write_log('  Loading artm')
       status = distributed_get_var(NCI%id, varid, &
            data%climate%artm, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling artm",GM_DIAGNOSTIC)
          data%climate%artm = &
               data%climate%artm*scaling_factor
       end if
    else
    call write_log('Variable artm was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'beta',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%velocity%beta)) then
       call write_log('  Loading beta')
       status = distributed_get_var(NCI%id, varid, &
            data%velocity%beta, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_beta)
       else
          scaling_factor = scaling_factor/(scale_beta)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling beta",GM_DIAGNOSTIC)
          data%velocity%beta = &
               data%velocity%beta*scaling_factor
       end if
    else
    call write_log('Variable beta was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'bfricflx',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%temper%bfricflx)) then
       call write_log('  Loading bfricflx')
       status = distributed_get_var(NCI%id, varid, &
            data%temper%bfricflx, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(1.0)
       else
          scaling_factor = scaling_factor/(1.0)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling bfricflx",GM_DIAGNOSTIC)
          data%temper%bfricflx = &
               data%temper%bfricflx*scaling_factor
       end if
    else
    call write_log('Variable bfricflx was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'bheatflx',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%temper%bheatflx)) then
       call write_log('  Loading bheatflx')
       status = distributed_get_var(NCI%id, varid, &
            data%temper%bheatflx, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_bflx)
       else
          scaling_factor = scaling_factor/(scale_bflx)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling bheatflx",GM_DIAGNOSTIC)
          data%temper%bheatflx = &
               data%temper%bheatflx*scaling_factor
       end if
    else
    call write_log('Variable bheatflx was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'bmlt_float',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%temper%bmlt_float)) then
       call write_log('  Loading bmlt_float')
       status = distributed_get_var(NCI%id, varid, &
            data%temper%bmlt_float, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_acab)
       else
          scaling_factor = scaling_factor/(scale_acab)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling bmlt_float",GM_DIAGNOSTIC)
          data%temper%bmlt_float = &
               data%temper%bmlt_float*scaling_factor
       end if
    else
    call write_log('Variable bmlt_float was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'bmlt_float_mask',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%temper%bmlt_float_mask)) then
       call write_log('  Loading bmlt_float_mask')
       status = distributed_get_var(NCI%id, varid, &
            data%temper%bmlt_float_mask, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling bmlt_float_mask",GM_DIAGNOSTIC)
          data%temper%bmlt_float_mask = &
               data%temper%bmlt_float_mask*scaling_factor
       end if
    else
    call write_log('Variable bmlt_float_mask was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'bmlt_ground',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%temper%bmlt_ground)) then
       call write_log('  Loading bmlt_ground')
       status = distributed_get_var(NCI%id, varid, &
            data%temper%bmlt_ground, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_acab)
       else
          scaling_factor = scaling_factor/(scale_acab)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling bmlt_ground",GM_DIAGNOSTIC)
          data%temper%bmlt_ground = &
               data%temper%bmlt_ground*scaling_factor
       end if
    else
    call write_log('Variable bmlt_ground was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'btractx',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%stress%btractx)) then
       call write_log('  Loading btractx')
       status = distributed_get_var(NCI%id, varid, &
            data%stress%btractx(:,:), (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_tau)
       else
          scaling_factor = scaling_factor/(scale_tau)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling btractx",GM_DIAGNOSTIC)
          data%stress%btractx(:,:) = &
               data%stress%btractx(:,:)*scaling_factor
       end if
    else
    call write_log('Variable btractx was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'btractx_extend',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%stress%btractx_extend)) then
       call write_log('  Loading btractx_extend')
       status = distributed_get_var(NCI%id, varid, &
            data%stress%btractx_extend(:,:), (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_tau)
       else
          scaling_factor = scaling_factor/(scale_tau)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling btractx_extend",GM_DIAGNOSTIC)
          data%stress%btractx_extend(:,:) = &
               data%stress%btractx_extend(:,:)*scaling_factor
       end if
    else
    call write_log('Variable btractx_extend was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'btracty',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%stress%btracty)) then
       call write_log('  Loading btracty')
       status = distributed_get_var(NCI%id, varid, &
            data%stress%btracty(:,:), (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_tau)
       else
          scaling_factor = scaling_factor/(scale_tau)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling btracty",GM_DIAGNOSTIC)
          data%stress%btracty(:,:) = &
               data%stress%btracty(:,:)*scaling_factor
       end if
    else
    call write_log('Variable btracty was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'btracty_extend',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%stress%btracty_extend)) then
       call write_log('  Loading btracty_extend')
       status = distributed_get_var(NCI%id, varid, &
            data%stress%btracty_extend(:,:), (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_tau)
       else
          scaling_factor = scaling_factor/(scale_tau)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling btracty_extend",GM_DIAGNOSTIC)
          data%stress%btracty_extend(:,:) = &
               data%stress%btracty_extend(:,:)*scaling_factor
       end if
    else
    call write_log('Variable btracty_extend was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'bwat',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%temper%bwat)) then
       call write_log('  Loading bwat')
       status = distributed_get_var(NCI%id, varid, &
            data%temper%bwat, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(thk0)
       else
          scaling_factor = scaling_factor/(thk0)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling bwat",GM_DIAGNOSTIC)
          data%temper%bwat = &
               data%temper%bwat*scaling_factor
       end if
    else
    call write_log('Variable bwat was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'damage',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%calving%damage)) then
       call write_log('  Loading damage')
       do up=1,NCI%nstaglevel
          status = distributed_get_var(NCI%id, varid, &
               data%calving%damage(up,:,:), (/1,1,up,infile%current_time/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
          if (status.ne.NF90_NOERR) then
             scaling_factor = 1.0d0
          end if
          if (abs(scaling_factor-1.0d0).gt.1.d-17) then
             call write_log("scaling damage",GM_DIAGNOSTIC)
             data%calving%damage(up,:,:) = &
                  data%calving%damage(up,:,:)*scaling_factor
          end if
       end do
    else
    call write_log('Variable damage was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'dissip',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%temper%dissip)) then
       call write_log('  Loading dissip')
       do up=1,NCI%nlevel
          status = distributed_get_var(NCI%id, varid, &
               data%temper%dissip(up,:,:), (/1,1,up,infile%current_time/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
          if (status.ne.NF90_NOERR) then
             scaling_factor = 1.0d0/(scyr)
          else
             scaling_factor = scaling_factor/(scyr)
          end if
          if (abs(scaling_factor-1.0d0).gt.1.d-17) then
             call write_log("scaling dissip",GM_DIAGNOSTIC)
             data%temper%dissip(up,:,:) = &
                  data%temper%dissip(up,:,:)*scaling_factor
          end if
       end do
    else
    call write_log('Variable dissip was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'dissipstag',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%temper%dissip)) then
       call write_log('  Loading dissipstag')
       do up=1,NCI%nstaglevel
          status = distributed_get_var(NCI%id, varid, &
               data%temper%dissip(up,:,:), (/1,1,up,infile%current_time/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
          if (status.ne.NF90_NOERR) then
             scaling_factor = 1.0d0/(scyr)
          else
             scaling_factor = scaling_factor/(scyr)
          end if
          if (abs(scaling_factor-1.0d0).gt.1.d-17) then
             call write_log("scaling dissipstag",GM_DIAGNOSTIC)
             data%temper%dissip(up,:,:) = &
                  data%temper%dissip(up,:,:)*scaling_factor
          end if
       end do
    else
    call write_log('Variable dissipstag was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'effecpress',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%basal_physics%effecpress)) then
       call write_log('  Loading effecpress')
       status = distributed_get_var(NCI%id, varid, &
            data%basal_physics%effecpress, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling effecpress",GM_DIAGNOSTIC)
          data%basal_physics%effecpress = &
               data%basal_physics%effecpress*scaling_factor
       end if
    else
    call write_log('Variable effecpress was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'efvs',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%stress%efvs)) then
       call write_log('  Loading efvs')
       do up=1,NCI%nstaglevel
          status = distributed_get_var(NCI%id, varid, &
               data%stress%efvs(up,:,:), (/1,1,up,infile%current_time/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
          if (status.ne.NF90_NOERR) then
             scaling_factor = 1.0d0/(scale_efvs)
          else
             scaling_factor = scaling_factor/(scale_efvs)
          end if
          if (abs(scaling_factor-1.0d0).gt.1.d-17) then
             call write_log("scaling efvs",GM_DIAGNOSTIC)
             data%stress%efvs(up,:,:) = &
                  data%stress%efvs(up,:,:)*scaling_factor
          end if
       end do
    else
    call write_log('Variable efvs was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'flwa',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%temper%flwa)) then
       call write_log('  Loading flwa')
       do up=1,NCI%nlevel
          status = distributed_get_var(NCI%id, varid, &
               data%temper%flwa(up,:,:), (/1,1,up,infile%current_time/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
          if (status.ne.NF90_NOERR) then
             scaling_factor = 1.0d0/(scale_flwa)
          else
             scaling_factor = scaling_factor/(scale_flwa)
          end if
          if (abs(scaling_factor-1.0d0).gt.1.d-17) then
             call write_log("scaling flwa",GM_DIAGNOSTIC)
             data%temper%flwa(up,:,:) = &
                  data%temper%flwa(up,:,:)*scaling_factor
          end if
       end do
    else
    call write_log('Variable flwa was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'flwastag',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%temper%flwa)) then
       call write_log('  Loading flwastag')
       do up=1,NCI%nstaglevel
          status = distributed_get_var(NCI%id, varid, &
               data%temper%flwa(up,:,:), (/1,1,up,infile%current_time/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
          if (status.ne.NF90_NOERR) then
             scaling_factor = 1.0d0/(scale_flwa)
          else
             scaling_factor = scaling_factor/(scale_flwa)
          end if
          if (abs(scaling_factor-1.0d0).gt.1.d-17) then
             call write_log("scaling flwastag",GM_DIAGNOSTIC)
             data%temper%flwa(up,:,:) = &
                  data%temper%flwa(up,:,:)*scaling_factor
          end if
       end do
    else
    call write_log('Variable flwastag was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'ice_age',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%geometry%ice_age)) then
       call write_log('  Loading ice_age')
       do up=1,NCI%nstaglevel
          status = distributed_get_var(NCI%id, varid, &
               data%geometry%ice_age(up,:,:), (/1,1,up,infile%current_time/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
          if (status.ne.NF90_NOERR) then
             scaling_factor = 1.0d0/(tim0/scyr)
          else
             scaling_factor = scaling_factor/(tim0/scyr)
          end if
          if (abs(scaling_factor-1.0d0).gt.1.d-17) then
             call write_log("scaling ice_age",GM_DIAGNOSTIC)
             data%geometry%ice_age(up,:,:) = &
                  data%geometry%ice_age(up,:,:)*scaling_factor
          end if
       end do
    else
    call write_log('Variable ice_age was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'kinbcmask',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%velocity%kinbcmask)) then
       call write_log('  Loading kinbcmask')
       status = distributed_get_var(NCI%id, varid, &
            data%velocity%kinbcmask(:,:), (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling kinbcmask",GM_DIAGNOSTIC)
          data%velocity%kinbcmask(:,:) = &
               data%velocity%kinbcmask(:,:)*scaling_factor
       end if
    else
    call write_log('Variable kinbcmask was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'litho_temp',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%lithot%temp)) then
       call write_log('  Loading litho_temp')
       status = distributed_get_var(NCI%id, varid, &
            data%lithot%temp, (/1,1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling litho_temp",GM_DIAGNOSTIC)
          data%lithot%temp = &
               data%lithot%temp*scaling_factor
       end if
    else
    call write_log('Variable litho_temp was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'relx',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%isostasy%relx)) then
       call write_log('  Loading relx')
       status = distributed_get_var(NCI%id, varid, &
            data%isostasy%relx, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(thk0)
       else
          scaling_factor = scaling_factor/(thk0)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling relx",GM_DIAGNOSTIC)
          data%isostasy%relx = &
               data%isostasy%relx*scaling_factor
       end if
    else
    call write_log('Variable relx was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'soft',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%velocity%bed_softness)) then
       call write_log('  Loading soft')
       status = distributed_get_var(NCI%id, varid, &
            data%velocity%bed_softness, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_btrc)
       else
          scaling_factor = scaling_factor/(scale_btrc)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling soft",GM_DIAGNOSTIC)
          data%velocity%bed_softness = &
               data%velocity%bed_softness*scaling_factor
       end if
    else
    call write_log('Variable soft was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'surftemp',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%climate%artm)) then
       call write_log('  Loading surftemp')
       status = distributed_get_var(NCI%id, varid, &
            data%climate%artm, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling surftemp",GM_DIAGNOSTIC)
          data%climate%artm = &
               data%climate%artm*scaling_factor
       end if
    else
    call write_log('Variable surftemp was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'taudx',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%stress%taudx)) then
       call write_log('  Loading taudx')
       status = distributed_get_var(NCI%id, varid, &
            data%stress%taudx(:,:), (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_tau)
       else
          scaling_factor = scaling_factor/(scale_tau)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling taudx",GM_DIAGNOSTIC)
          data%stress%taudx(:,:) = &
               data%stress%taudx(:,:)*scaling_factor
       end if
    else
    call write_log('Variable taudx was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'taudy',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%stress%taudy)) then
       call write_log('  Loading taudy')
       status = distributed_get_var(NCI%id, varid, &
            data%stress%taudy(:,:), (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_tau)
       else
          scaling_factor = scaling_factor/(scale_tau)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling taudy",GM_DIAGNOSTIC)
          data%stress%taudy(:,:) = &
               data%stress%taudy(:,:)*scaling_factor
       end if
    else
    call write_log('Variable taudy was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'tauf',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%basalproc%mintauf)) then
       call write_log('  Loading tauf')
       status = distributed_get_var(NCI%id, varid, &
            data%basalproc%mintauf, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_tau)
       else
          scaling_factor = scaling_factor/(scale_tau)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling tauf",GM_DIAGNOSTIC)
          data%basalproc%mintauf = &
               data%basalproc%mintauf*scaling_factor
       end if
    else
    call write_log('Variable tauf was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'temp',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%temper%temp)) then
       call write_log('  Loading temp')
       do up=1,NCI%nlevel
          status = distributed_get_var(NCI%id, varid, &
               data%temper%temp(up,1:data%general%ewn,1:data%general%nsn), (/1,1,up,infile%current_time/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
          if (status.ne.NF90_NOERR) then
             scaling_factor = 1.0d0
          end if
          if (abs(scaling_factor-1.0d0).gt.1.d-17) then
             call write_log("scaling temp",GM_DIAGNOSTIC)
             data%temper%temp(up,1:data%general%ewn,1:data%general%nsn) = &
                  data%temper%temp(up,1:data%general%ewn,1:data%general%nsn)*scaling_factor
          end if
       end do
    else
    call write_log('Variable temp was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'tempstag',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%temper%temp)) then
       call write_log('  Loading tempstag')
       do up=0,NCI%nstagwbndlevel
          status = distributed_get_var(NCI%id, varid, &
               data%temper%temp(up,1:data%general%ewn,1:data%general%nsn), (/1,1,up+1,infile%current_time/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
          if (status.ne.NF90_NOERR) then
             scaling_factor = 1.0d0
          end if
          if (abs(scaling_factor-1.0d0).gt.1.d-17) then
             call write_log("scaling tempstag",GM_DIAGNOSTIC)
             data%temper%temp(up,1:data%general%ewn,1:data%general%nsn) = &
                  data%temper%temp(up,1:data%general%ewn,1:data%general%nsn)*scaling_factor
          end if
       end do
    else
    call write_log('Variable tempstag was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'thk',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%geometry%thck)) then
       call write_log('  Loading thk')
       status = distributed_get_var(NCI%id, varid, &
            data%geometry%thck, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(thk0)
       else
          scaling_factor = scaling_factor/(thk0)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling thk",GM_DIAGNOSTIC)
          data%geometry%thck = &
               data%geometry%thck*scaling_factor
       end if
    else
    call write_log('Variable thk was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'thkmask',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%geometry%thkmask)) then
       call write_log('  Loading thkmask')
       status = distributed_get_var(NCI%id, varid, &
            data%geometry%thkmask, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling thkmask",GM_DIAGNOSTIC)
          data%geometry%thkmask = &
               data%geometry%thkmask*scaling_factor
       end if
    else
    call write_log('Variable thkmask was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'topg',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%geometry%topg)) then
       call write_log('  Loading topg')
       status = distributed_get_var(NCI%id, varid, &
            data%geometry%topg, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(thk0)
       else
          scaling_factor = scaling_factor/(thk0)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling topg",GM_DIAGNOSTIC)
          data%geometry%topg = &
               data%geometry%topg*scaling_factor
       end if
    else
    call write_log('Variable topg was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'ubas',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%velocity%ubas)) then
       call write_log('  Loading ubas')
       status = distributed_get_var(NCI%id, varid, &
            data%velocity%ubas, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_uvel)
       else
          scaling_factor = scaling_factor/(scale_uvel)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling ubas",GM_DIAGNOSTIC)
          data%velocity%ubas = &
               data%velocity%ubas*scaling_factor
       end if
    else
    call write_log('Variable ubas was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'unstagbeta',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%velocity%unstagbeta)) then
       call write_log('  Loading unstagbeta')
       status = distributed_get_var(NCI%id, varid, &
            data%velocity%unstagbeta, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_beta)
       else
          scaling_factor = scaling_factor/(scale_beta)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling unstagbeta",GM_DIAGNOSTIC)
          data%velocity%unstagbeta = &
               data%velocity%unstagbeta*scaling_factor
       end if
    else
    call write_log('Variable unstagbeta was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'usurf',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%geometry%usrf)) then
       call write_log('  Loading usurf')
       status = distributed_get_var(NCI%id, varid, &
            data%geometry%usrf, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(thk0)
       else
          scaling_factor = scaling_factor/(thk0)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling usurf",GM_DIAGNOSTIC)
          data%geometry%usrf = &
               data%geometry%usrf*scaling_factor
       end if
    else
    call write_log('Variable usurf was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'uvel',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%velocity%uvel)) then
       call write_log('  Loading uvel')
       do up=1,NCI%nlevel
          status = distributed_get_var(NCI%id, varid, &
               data%velocity%uvel(up,:,:), (/1,1,up,infile%current_time/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
          if (status.ne.NF90_NOERR) then
             scaling_factor = 1.0d0/(scale_uvel)
          else
             scaling_factor = scaling_factor/(scale_uvel)
          end if
          if (abs(scaling_factor-1.0d0).gt.1.d-17) then
             call write_log("scaling uvel",GM_DIAGNOSTIC)
             data%velocity%uvel(up,:,:) = &
                  data%velocity%uvel(up,:,:)*scaling_factor
          end if
       end do
    else
    call write_log('Variable uvel was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'uvel_2d',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%velocity%uvel_2d)) then
       call write_log('  Loading uvel_2d')
       status = distributed_get_var(NCI%id, varid, &
            data%velocity%uvel_2d(:,:), (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_uvel)
       else
          scaling_factor = scaling_factor/(scale_uvel)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling uvel_2d",GM_DIAGNOSTIC)
          data%velocity%uvel_2d(:,:) = &
               data%velocity%uvel_2d(:,:)*scaling_factor
       end if
    else
    call write_log('Variable uvel_2d was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'uvel_2d_extend',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%velocity%uvel_2d_extend)) then
       call write_log('  Loading uvel_2d_extend')
       status = distributed_get_var(NCI%id, varid, &
            data%velocity%uvel_2d_extend(:,:), (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_uvel)
       else
          scaling_factor = scaling_factor/(scale_uvel)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling uvel_2d_extend",GM_DIAGNOSTIC)
          data%velocity%uvel_2d_extend(:,:) = &
               data%velocity%uvel_2d_extend(:,:)*scaling_factor
       end if
    else
    call write_log('Variable uvel_2d_extend was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'uvel_extend',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%velocity%uvel_extend)) then
       call write_log('  Loading uvel_extend')
       do up=1,NCI%nlevel
          status = distributed_get_var(NCI%id, varid, &
               data%velocity%uvel_extend(up,:,:), (/1,1,up,infile%current_time/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
          if (status.ne.NF90_NOERR) then
             scaling_factor = 1.0d0/(scale_uvel)
          else
             scaling_factor = scaling_factor/(scale_uvel)
          end if
          if (abs(scaling_factor-1.0d0).gt.1.d-17) then
             call write_log("scaling uvel_extend",GM_DIAGNOSTIC)
             data%velocity%uvel_extend(up,:,:) = &
                  data%velocity%uvel_extend(up,:,:)*scaling_factor
          end if
       end do
    else
    call write_log('Variable uvel_extend was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'vbas',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%velocity%vbas)) then
       call write_log('  Loading vbas')
       status = distributed_get_var(NCI%id, varid, &
            data%velocity%vbas, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_uvel)
       else
          scaling_factor = scaling_factor/(scale_uvel)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling vbas",GM_DIAGNOSTIC)
          data%velocity%vbas = &
               data%velocity%vbas*scaling_factor
       end if
    else
    call write_log('Variable vbas was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'vvel',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%velocity%vvel)) then
       call write_log('  Loading vvel')
       do up=1,NCI%nlevel
          status = distributed_get_var(NCI%id, varid, &
               data%velocity%vvel(up,:,:), (/1,1,up,infile%current_time/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
          if (status.ne.NF90_NOERR) then
             scaling_factor = 1.0d0/(scale_uvel)
          else
             scaling_factor = scaling_factor/(scale_uvel)
          end if
          if (abs(scaling_factor-1.0d0).gt.1.d-17) then
             call write_log("scaling vvel",GM_DIAGNOSTIC)
             data%velocity%vvel(up,:,:) = &
                  data%velocity%vvel(up,:,:)*scaling_factor
          end if
       end do
    else
    call write_log('Variable vvel was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'vvel_2d',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%velocity%vvel_2d)) then
       call write_log('  Loading vvel_2d')
       status = distributed_get_var(NCI%id, varid, &
            data%velocity%vvel_2d(:,:), (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_uvel)
       else
          scaling_factor = scaling_factor/(scale_uvel)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling vvel_2d",GM_DIAGNOSTIC)
          data%velocity%vvel_2d(:,:) = &
               data%velocity%vvel_2d(:,:)*scaling_factor
       end if
    else
    call write_log('Variable vvel_2d was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'vvel_2d_extend',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%velocity%vvel_2d_extend)) then
       call write_log('  Loading vvel_2d_extend')
       status = distributed_get_var(NCI%id, varid, &
            data%velocity%vvel_2d_extend(:,:), (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_uvel)
       else
          scaling_factor = scaling_factor/(scale_uvel)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling vvel_2d_extend",GM_DIAGNOSTIC)
          data%velocity%vvel_2d_extend(:,:) = &
               data%velocity%vvel_2d_extend(:,:)*scaling_factor
       end if
    else
    call write_log('Variable vvel_2d_extend was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'vvel_extend',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%velocity%vvel_extend)) then
       call write_log('  Loading vvel_extend')
       do up=1,NCI%nlevel
          status = distributed_get_var(NCI%id, varid, &
               data%velocity%vvel_extend(up,:,:), (/1,1,up,infile%current_time/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
          if (status.ne.NF90_NOERR) then
             scaling_factor = 1.0d0/(scale_uvel)
          else
             scaling_factor = scaling_factor/(scale_uvel)
          end if
          if (abs(scaling_factor-1.0d0).gt.1.d-17) then
             call write_log("scaling vvel_extend",GM_DIAGNOSTIC)
             data%velocity%vvel_extend(up,:,:) = &
                  data%velocity%vvel_extend(up,:,:)*scaling_factor
          end if
       end do
    else
    call write_log('Variable vvel_extend was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'waterfrac',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%temper%waterfrac)) then
       call write_log('  Loading waterfrac')
       do up=1,NCI%nstaglevel
          status = distributed_get_var(NCI%id, varid, &
               data%temper%waterfrac(up,:,:), (/1,1,up,infile%current_time/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
          if (status.ne.NF90_NOERR) then
             scaling_factor = 1.0d0
          end if
          if (abs(scaling_factor-1.0d0).gt.1.d-17) then
             call write_log("scaling waterfrac",GM_DIAGNOSTIC)
             data%temper%waterfrac(up,:,:) = &
                  data%temper%waterfrac(up,:,:)*scaling_factor
          end if
       end do
    else
    call write_log('Variable waterfrac was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'wgrd',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%velocity%wgrd)) then
       call write_log('  Loading wgrd')
       do up=1,NCI%nlevel
          status = distributed_get_var(NCI%id, varid, &
               data%velocity%wgrd(up,:,:), (/1,1,up,infile%current_time/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
          if (status.ne.NF90_NOERR) then
             scaling_factor = 1.0d0/(scale_wvel)
          else
             scaling_factor = scaling_factor/(scale_wvel)
          end if
          if (abs(scaling_factor-1.0d0).gt.1.d-17) then
             call write_log("scaling wgrd",GM_DIAGNOSTIC)
             data%velocity%wgrd(up,:,:) = &
                  data%velocity%wgrd(up,:,:)*scaling_factor
          end if
       end do
    else
    call write_log('Variable wgrd was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'wvel',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%velocity%wvel)) then
       call write_log('  Loading wvel')
       do up=1,NCI%nlevel
          status = distributed_get_var(NCI%id, varid, &
               data%velocity%wvel(up,:,:), (/1,1,up,infile%current_time/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
          if (status.ne.NF90_NOERR) then
             scaling_factor = 1.0d0/(scale_wvel)
          else
             scaling_factor = scaling_factor/(scale_wvel)
          end if
          if (abs(scaling_factor-1.0d0).gt.1.d-17) then
             call write_log("scaling wvel",GM_DIAGNOSTIC)
             data%velocity%wvel(up,:,:) = &
                  data%velocity%wvel(up,:,:)*scaling_factor
          end if
       end do
    else
    call write_log('Variable wvel was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'wvel_ho',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%velocity%wvel_ho)) then
       call write_log('  Loading wvel_ho')
       do up=1,NCI%nlevel
          status = distributed_get_var(NCI%id, varid, &
               data%velocity%wvel_ho(up,:,:), (/1,1,up,infile%current_time/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
          if (status.ne.NF90_NOERR) then
             scaling_factor = 1.0d0/(scale_wvel)
          else
             scaling_factor = scaling_factor/(scale_wvel)
          end if
          if (abs(scaling_factor-1.0d0).gt.1.d-17) then
             call write_log("scaling wvel_ho",GM_DIAGNOSTIC)
             data%velocity%wvel_ho(up,:,:) = &
                  data%velocity%wvel_ho(up,:,:)*scaling_factor
          end if
       end do
    else
    call write_log('Variable wvel_ho was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

  end subroutine glide_io_read

  subroutine glide_io_checkdim(infile,model,data)
    ! check if dimension sizes in file match dims of model
    use parallel
    use glimmer_log
    use glimmer_ncdf
    use glide_types
    use glide_types
    implicit none
    type(glimmer_nc_input), pointer :: infile
    ! structure containg output netCDF descriptor
    type(glide_global_type) :: model
    type(glide_global_type), optional :: data

    integer status,dimid,dimsize
    character(len=150) message

    ! check dimensions
    status = parallel_inq_dimid(NCI%id,'level',dimid)
    if (dimid.gt.0) then
       status = parallel_inquire_dimension(NCI%id,dimid,len=dimsize)
       if (dimsize.ne.model%general%upn) then
          write(message,*) 'Error, reading file ',trim(NCI%filename),' size level does not match: ', &
               model%general%upn
          call write_log(message,GM_FATAL)
       end if
    end if
    status = parallel_inq_dimid(NCI%id,'lithoz',dimid)
    if (dimid.gt.0) then
       status = parallel_inquire_dimension(NCI%id,dimid,len=dimsize)
       if (dimsize.ne.model%lithot%nlayer) then
          write(message,*) 'Error, reading file ',trim(NCI%filename),' size lithoz does not match: ', &
               model%lithot%nlayer
          call write_log(message,GM_FATAL)
       end if
    end if
    status = parallel_inq_dimid(NCI%id,'staglevel',dimid)
    if (dimid.gt.0) then
       status = parallel_inquire_dimension(NCI%id,dimid,len=dimsize)
       if (dimsize.ne.model%general%upn-1) then
          write(message,*) 'Error, reading file ',trim(NCI%filename),' size staglevel does not match: ', &
               model%general%upn-1
          call write_log(message,GM_FATAL)
       end if
    end if
    status = parallel_inq_dimid(NCI%id,'stagwbndlevel',dimid)
    if (dimid.gt.0) then
       status = parallel_inquire_dimension(NCI%id,dimid,len=dimsize)
       if (dimsize.ne.model%general%upn+1) then
          write(message,*) 'Error, reading file ',trim(NCI%filename),' size stagwbndlevel does not match: ', &
               model%general%upn+1
          call write_log(message,GM_FATAL)
       end if
    end if
    status = parallel_inq_dimid(NCI%id,'x0',dimid)
    if (dimid.gt.0) then
       status = parallel_inquire_dimension(NCI%id,dimid,len=dimsize)
       if (dimsize.ne.global_ewn-1) then
          write(message,*) 'Error, reading file ',trim(NCI%filename),' size x0 does not match: ', &
               global_ewn-1
          call write_log(message,GM_FATAL)
       end if
    end if
    status = parallel_inq_dimid(NCI%id,'x1',dimid)
    if (dimid.gt.0) then
       status = parallel_inquire_dimension(NCI%id,dimid,len=dimsize)
       if (dimsize.ne.global_ewn) then
          write(message,*) 'Error, reading file ',trim(NCI%filename),' size x1 does not match: ', &
               global_ewn
          call write_log(message,GM_FATAL)
       end if
    end if
    status = parallel_inq_dimid(NCI%id,'y0',dimid)
    if (dimid.gt.0) then
       status = parallel_inquire_dimension(NCI%id,dimid,len=dimsize)
       if (dimsize.ne.global_nsn-1) then
          write(message,*) 'Error, reading file ',trim(NCI%filename),' size y0 does not match: ', &
               global_nsn-1
          call write_log(message,GM_FATAL)
       end if
    end if
    status = parallel_inq_dimid(NCI%id,'y1',dimid)
    if (dimid.gt.0) then
       status = parallel_inquire_dimension(NCI%id,dimid,len=dimsize)
       if (dimsize.ne.global_nsn) then
          write(message,*) 'Error, reading file ',trim(NCI%filename),' size y1 does not match: ', &
               global_nsn
          call write_log(message,GM_FATAL)
       end if
    end if
  end subroutine glide_io_checkdim

  !*****************************************************************************
  ! calculating time averages
  !*****************************************************************************  
#ifdef HAVE_AVG
  subroutine glide_avg_accumulate(outfile,data,model)
    use parallel
    use glide_types
    use glide_types
    use glimmer_ncdf
    implicit none
    type(glimmer_nc_output), pointer :: outfile
    ! structure containg output netCDF descriptor
    type(glide_global_type) :: model
    type(glide_global_type) :: data

    ! local variables
    real(dp) :: factor
    integer status, varid

    ! increase total time
    outfile%total_time = outfile%total_time + model%numerics%tinc
    factor = model%numerics%tinc

  end subroutine glide_avg_accumulate

  subroutine glide_avg_reset(outfile,data)
    use parallel
    use glide_types
    use glimmer_ncdf
    implicit none
    type(glimmer_nc_output), pointer :: outfile
    ! structure containg output netCDF descriptor
    type(glide_global_type) :: data

    ! local variables
    integer status, varid

    ! reset total time
    outfile%total_time = 0.d0

  end subroutine glide_avg_reset
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
  subroutine glide_get_C_space_factor(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%basal_physics%C_space_factor
  end subroutine glide_get_C_space_factor

  subroutine glide_set_C_space_factor(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%basal_physics%C_space_factor = inarray
  end subroutine glide_set_C_space_factor

  subroutine glide_get_acab(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_acab)*(data%climate%acab)
  end subroutine glide_get_acab

  subroutine glide_set_acab(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%climate%acab = inarray/(scale_acab)
  end subroutine glide_set_acab

  subroutine glide_get_adv_cfl_dt(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(out) :: outarray

    outarray = data%numerics%adv_cfl_dt
  end subroutine glide_get_adv_cfl_dt

  subroutine glide_set_adv_cfl_dt(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(in) :: inarray

    data%numerics%adv_cfl_dt = inarray
  end subroutine glide_set_adv_cfl_dt

  subroutine glide_get_artm(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%climate%artm
  end subroutine glide_get_artm

  subroutine glide_set_artm(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%climate%artm = inarray
  end subroutine glide_set_artm

  subroutine glide_get_beta(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_beta)*(data%velocity%beta)
  end subroutine glide_get_beta

  subroutine glide_set_beta(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%velocity%beta = inarray/(scale_beta)
  end subroutine glide_set_beta

  subroutine glide_get_beta_internal(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_beta)*(data%velocity%beta_internal)
  end subroutine glide_get_beta_internal

  subroutine glide_set_beta_internal(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%velocity%beta_internal = inarray/(scale_beta)
  end subroutine glide_set_beta_internal

  subroutine glide_get_bfricflx(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (1.0)*(data%temper%bfricflx)
  end subroutine glide_get_bfricflx

  subroutine glide_set_bfricflx(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%temper%bfricflx = inarray/(1.0)
  end subroutine glide_set_bfricflx

  subroutine glide_get_bheatflx(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_bflx)*(data%temper%bheatflx)
  end subroutine glide_get_bheatflx

  subroutine glide_set_bheatflx(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%temper%bheatflx = inarray/(scale_bflx)
  end subroutine glide_set_bheatflx

  subroutine glide_get_bmlt_float(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_acab)*(data%temper%bmlt_float)
  end subroutine glide_get_bmlt_float

  subroutine glide_set_bmlt_float(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%temper%bmlt_float = inarray/(scale_acab)
  end subroutine glide_set_bmlt_float

  subroutine glide_get_bmlt_float_mask(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%temper%bmlt_float_mask
  end subroutine glide_get_bmlt_float_mask

  subroutine glide_set_bmlt_float_mask(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%temper%bmlt_float_mask = inarray
  end subroutine glide_set_bmlt_float_mask

  subroutine glide_get_bmlt_ground(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_acab)*(data%temper%bmlt_ground)
  end subroutine glide_get_bmlt_ground

  subroutine glide_set_bmlt_ground(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%temper%bmlt_ground = inarray/(scale_acab)
  end subroutine glide_set_bmlt_ground

  subroutine glide_get_btemp(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%temper%temp(data%general%upn,1:data%general%ewn,1:data%general%nsn)
  end subroutine glide_get_btemp

  subroutine glide_get_btractx(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_tau)*(data%stress%btractx(:,:))
  end subroutine glide_get_btractx

  subroutine glide_set_btractx(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%stress%btractx(:,:) = inarray/(scale_tau)
  end subroutine glide_set_btractx

  subroutine glide_get_btractx_extend(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_tau)*(data%stress%btractx_extend(:,:))
  end subroutine glide_get_btractx_extend

  subroutine glide_set_btractx_extend(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%stress%btractx_extend(:,:) = inarray/(scale_tau)
  end subroutine glide_set_btractx_extend

  subroutine glide_get_btracty(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_tau)*(data%stress%btracty(:,:))
  end subroutine glide_get_btracty

  subroutine glide_set_btracty(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%stress%btracty(:,:) = inarray/(scale_tau)
  end subroutine glide_set_btracty

  subroutine glide_get_btracty_extend(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_tau)*(data%stress%btracty_extend(:,:))
  end subroutine glide_get_btracty_extend

  subroutine glide_set_btracty_extend(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%stress%btracty_extend(:,:) = inarray/(scale_tau)
  end subroutine glide_set_btracty_extend

  subroutine glide_get_btrc(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_btrc)*(data%velocity%btrc)
  end subroutine glide_get_btrc

  subroutine glide_set_btrc(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%velocity%btrc = inarray/(scale_btrc)
  end subroutine glide_set_btrc

  subroutine glide_get_bwat(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (thk0)*(data%temper%bwat)
  end subroutine glide_get_bwat

  subroutine glide_set_bwat(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%temper%bwat = inarray/(thk0)
  end subroutine glide_set_bwat

  subroutine glide_get_bwatflx(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (thk0)*(data%temper%bwatflx)
  end subroutine glide_get_bwatflx

  subroutine glide_set_bwatflx(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%temper%bwatflx = inarray/(thk0)
  end subroutine glide_set_bwatflx

  subroutine glide_get_calving(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (thk0)*(data%calving%calving_thck)
  end subroutine glide_get_calving

  subroutine glide_set_calving(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%calving%calving_thck = inarray/(thk0)
  end subroutine glide_set_calving

  subroutine glide_get_damage_column(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%calving%damage_column
  end subroutine glide_get_damage_column

  subroutine glide_set_damage_column(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%calving%damage_column = inarray
  end subroutine glide_set_damage_column

  subroutine glide_get_diff_cfl_dt(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(out) :: outarray

    outarray = data%numerics%diff_cfl_dt
  end subroutine glide_get_diff_cfl_dt

  subroutine glide_set_diff_cfl_dt(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(in) :: inarray

    data%numerics%diff_cfl_dt = inarray
  end subroutine glide_set_diff_cfl_dt

  subroutine glide_get_diffu(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_diffu)*(data%velocity%diffu)
  end subroutine glide_get_diffu

  subroutine glide_set_diffu(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%velocity%diffu = inarray/(scale_diffu)
  end subroutine glide_set_diffu

  subroutine glide_get_dthckdtm(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_acab)*(data%geomderv%dthckdtm)
  end subroutine glide_get_dthckdtm

  subroutine glide_set_dthckdtm(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%geomderv%dthckdtm = inarray/(scale_acab)
  end subroutine glide_set_dthckdtm

  subroutine glide_get_dusrfdtm(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_acab)*(data%geomderv%dusrfdtm)
  end subroutine glide_get_dusrfdtm

  subroutine glide_set_dusrfdtm(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%geomderv%dusrfdtm = inarray/(scale_acab)
  end subroutine glide_set_dusrfdtm

  subroutine glide_get_dynbcmask(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    integer, dimension(:,:), intent(out) :: outarray

    outarray = data%velocity%dynbcmask
  end subroutine glide_get_dynbcmask

  subroutine glide_set_dynbcmask(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    integer, dimension(:,:), intent(in) :: inarray

    data%velocity%dynbcmask = inarray
  end subroutine glide_set_dynbcmask

  subroutine glide_get_effecpress(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%basal_physics%effecpress
  end subroutine glide_get_effecpress

  subroutine glide_set_effecpress(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%basal_physics%effecpress = inarray
  end subroutine glide_set_effecpress

  subroutine glide_get_eus(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(out) :: outarray

    outarray = (thk0)*(data%climate%eus)
  end subroutine glide_get_eus

  subroutine glide_set_eus(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(in) :: inarray

    data%climate%eus = inarray/(thk0)
  end subroutine glide_set_eus

  subroutine glide_get_f_flotation(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%geometry%f_flotation
  end subroutine glide_get_f_flotation

  subroutine glide_set_f_flotation(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%geometry%f_flotation = inarray
  end subroutine glide_set_f_flotation

  subroutine glide_get_f_ground(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%geometry%f_ground
  end subroutine glide_get_f_ground

  subroutine glide_set_f_ground(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%geometry%f_ground = inarray
  end subroutine glide_set_f_ground

  subroutine glide_get_gravity(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(out) :: outarray

    outarray = grav
  end subroutine glide_get_gravity

  subroutine glide_set_gravity(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(in) :: inarray

!  no rescaling here
  end subroutine glide_set_gravity

  subroutine glide_get_iarea(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(out) :: outarray

    outarray = (len0*len0*1.e-6)*(data%geometry%iarea)
  end subroutine glide_get_iarea

  subroutine glide_set_iarea(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(in) :: inarray

    data%geometry%iarea = inarray/(len0*len0*1.e-6)
  end subroutine glide_set_iarea

  subroutine glide_get_iareaf(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(out) :: outarray

    outarray = (len0*len0*1.e-6)*(data%geometry%iareaf)
  end subroutine glide_get_iareaf

  subroutine glide_set_iareaf(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(in) :: inarray

    data%geometry%iareaf = inarray/(len0*len0*1.e-6)
  end subroutine glide_set_iareaf

  subroutine glide_get_iareag(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(out) :: outarray

    outarray = (len0*len0*1.e-6)*(data%geometry%iareag)
  end subroutine glide_get_iareag

  subroutine glide_set_iareag(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(in) :: inarray

    data%geometry%iareag = inarray/(len0*len0*1.e-6)
  end subroutine glide_set_iareag

  subroutine glide_get_ice_mask(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (1.0)*(data%geometry%ice_mask)
  end subroutine glide_get_ice_mask

  subroutine glide_set_ice_mask(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%geometry%ice_mask = inarray/(1.0)
  end subroutine glide_set_ice_mask

  subroutine glide_get_ice_specific_heat(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(out) :: outarray

    outarray = shci
  end subroutine glide_get_ice_specific_heat

  subroutine glide_set_ice_specific_heat(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(in) :: inarray

!  no rescaling here
  end subroutine glide_set_ice_specific_heat

  subroutine glide_get_ice_thermal_conductivity(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(out) :: outarray

    outarray = coni
  end subroutine glide_get_ice_thermal_conductivity

  subroutine glide_set_ice_thermal_conductivity(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(in) :: inarray

!  no rescaling here
  end subroutine glide_set_ice_thermal_conductivity

  subroutine glide_get_ivol(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(out) :: outarray

    outarray = (thk0*len0*len0*1.e-9)*(data%geometry%ivol)
  end subroutine glide_get_ivol

  subroutine glide_set_ivol(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(in) :: inarray

    data%geometry%ivol = inarray/(thk0*len0*len0*1.e-9)
  end subroutine glide_set_ivol

  subroutine glide_get_kinbcmask(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    integer, dimension(:,:), intent(out) :: outarray

    outarray = data%velocity%kinbcmask(:,:)
  end subroutine glide_get_kinbcmask

  subroutine glide_set_kinbcmask(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    integer, dimension(:,:), intent(in) :: inarray

    data%velocity%kinbcmask(:,:) = inarray
  end subroutine glide_set_kinbcmask

  subroutine glide_get_lsurf(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (thk0)*(data%geometry%lsrf)
  end subroutine glide_get_lsurf

  subroutine glide_set_lsurf(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%geometry%lsrf = inarray/(thk0)
  end subroutine glide_set_lsurf

  subroutine glide_get_relx(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (thk0)*(data%isostasy%relx)
  end subroutine glide_get_relx

  subroutine glide_set_relx(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%isostasy%relx = inarray/(thk0)
  end subroutine glide_set_relx

  subroutine glide_get_rho_ice(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(out) :: outarray

    outarray = rhoi
  end subroutine glide_get_rho_ice

  subroutine glide_set_rho_ice(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(in) :: inarray

!  no rescaling here
  end subroutine glide_set_rho_ice

  subroutine glide_get_rho_seawater(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(out) :: outarray

    outarray = rhoo
  end subroutine glide_get_rho_seawater

  subroutine glide_set_rho_seawater(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(in) :: inarray

!  no rescaling here
  end subroutine glide_set_rho_seawater

  subroutine glide_get_seconds_per_year(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(out) :: outarray

    outarray = scyr
  end subroutine glide_get_seconds_per_year

  subroutine glide_set_seconds_per_year(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(in) :: inarray

!  no rescaling here
  end subroutine glide_set_seconds_per_year

  subroutine glide_get_soft(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_btrc)*(data%velocity%bed_softness)
  end subroutine glide_get_soft

  subroutine glide_set_soft(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%velocity%bed_softness = inarray/(scale_btrc)
  end subroutine glide_set_soft

  subroutine glide_get_stagthk(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (thk0)*(data%geomderv%stagthck)
  end subroutine glide_get_stagthk

  subroutine glide_set_stagthk(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%geomderv%stagthck = inarray/(thk0)
  end subroutine glide_set_stagthk

  subroutine glide_get_surftemp(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%climate%artm
  end subroutine glide_get_surftemp

  subroutine glide_set_surftemp(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%climate%artm = inarray
  end subroutine glide_set_surftemp

  subroutine glide_get_taudx(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_tau)*(data%stress%taudx(:,:))
  end subroutine glide_get_taudx

  subroutine glide_set_taudx(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%stress%taudx(:,:) = inarray/(scale_tau)
  end subroutine glide_set_taudx

  subroutine glide_get_taudy(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_tau)*(data%stress%taudy(:,:))
  end subroutine glide_get_taudy

  subroutine glide_set_taudy(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%stress%taudy(:,:) = inarray/(scale_tau)
  end subroutine glide_set_taudy

  subroutine glide_get_tauf(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_tau)*(data%basalproc%mintauf)
  end subroutine glide_get_tauf

  subroutine glide_set_tauf(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%basalproc%mintauf = inarray/(scale_tau)
  end subroutine glide_set_tauf

  subroutine glide_get_taux(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (1e-3*thk0*thk0/len0)*(data%velocity%tau_x)
  end subroutine glide_get_taux

  subroutine glide_set_taux(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%velocity%tau_x = inarray/(1e-3*thk0*thk0/len0)
  end subroutine glide_set_taux

  subroutine glide_get_tauy(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (1e-3*thk0*thk0/len0)*(data%velocity%tau_y)
  end subroutine glide_get_tauy

  subroutine glide_set_tauy(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%velocity%tau_y = inarray/(1e-3*thk0*thk0/len0)
  end subroutine glide_set_tauy

  subroutine glide_get_thk(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (thk0)*(data%geometry%thck)
  end subroutine glide_get_thk

  subroutine glide_set_thk(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%geometry%thck = inarray/(thk0)
  end subroutine glide_set_thk

  subroutine glide_get_thkmask(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    integer, dimension(:,:), intent(out) :: outarray

    outarray = data%geometry%thkmask
  end subroutine glide_get_thkmask

  subroutine glide_set_thkmask(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    integer, dimension(:,:), intent(in) :: inarray

    data%geometry%thkmask = inarray
  end subroutine glide_set_thkmask

  subroutine glide_get_topg(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (thk0)*(data%geometry%topg)
  end subroutine glide_get_topg

  subroutine glide_set_topg(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%geometry%topg = inarray/(thk0)
  end subroutine glide_set_topg

  subroutine glide_get_ubas(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_uvel)*(data%velocity%ubas)
  end subroutine glide_get_ubas

  subroutine glide_set_ubas(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%velocity%ubas = inarray/(scale_uvel)
  end subroutine glide_set_ubas

  subroutine glide_get_uflx(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_uflx)*(data%velocity%uflx)
  end subroutine glide_get_uflx

  subroutine glide_set_uflx(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%velocity%uflx = inarray/(scale_uflx)
  end subroutine glide_set_uflx

  subroutine glide_get_unstagbeta(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_beta)*(data%velocity%unstagbeta)
  end subroutine glide_get_unstagbeta

  subroutine glide_set_unstagbeta(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%velocity%unstagbeta = inarray/(scale_beta)
  end subroutine glide_set_unstagbeta

  subroutine glide_get_usurf(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (thk0)*(data%geometry%usrf)
  end subroutine glide_get_usurf

  subroutine glide_set_usurf(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%geometry%usrf = inarray/(thk0)
  end subroutine glide_set_usurf

  subroutine glide_get_uvel_2d(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_uvel)*(data%velocity%uvel_2d(:,:))
  end subroutine glide_get_uvel_2d

  subroutine glide_set_uvel_2d(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%velocity%uvel_2d(:,:) = inarray/(scale_uvel)
  end subroutine glide_set_uvel_2d

  subroutine glide_get_uvel_2d_extend(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_uvel)*(data%velocity%uvel_2d_extend(:,:))
  end subroutine glide_get_uvel_2d_extend

  subroutine glide_set_uvel_2d_extend(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%velocity%uvel_2d_extend(:,:) = inarray/(scale_uvel)
  end subroutine glide_set_uvel_2d_extend

  subroutine glide_get_vbas(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_uvel)*(data%velocity%vbas)
  end subroutine glide_get_vbas

  subroutine glide_set_vbas(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%velocity%vbas = inarray/(scale_uvel)
  end subroutine glide_set_vbas

  subroutine glide_get_vflx(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_uflx)*(data%velocity%vflx)
  end subroutine glide_get_vflx

  subroutine glide_set_vflx(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%velocity%vflx = inarray/(scale_uflx)
  end subroutine glide_set_vflx

  subroutine glide_get_vvel_2d(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_uvel)*(data%velocity%vvel_2d(:,:))
  end subroutine glide_get_vvel_2d

  subroutine glide_set_vvel_2d(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%velocity%vvel_2d(:,:) = inarray/(scale_uvel)
  end subroutine glide_set_vvel_2d

  subroutine glide_get_vvel_2d_extend(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_uvel)*(data%velocity%vvel_2d_extend(:,:))
  end subroutine glide_get_vvel_2d_extend

  subroutine glide_set_vvel_2d_extend(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use glide_types
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%velocity%vvel_2d_extend(:,:) = inarray/(scale_uvel)
  end subroutine glide_set_vvel_2d_extend


end module glide_io
