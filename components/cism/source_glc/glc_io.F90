!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module glc_io

!BOP
! !MODULE: glc_io

! !DESCRIPTION:
!  Contains routines for specialized glc IO
!
! !REVISION HISTORY:
!
! !USES:

   use glc_time_management, only: iyear, imonth, iday, ihour, iminute, isecond, &
                                  runtype, cesm_date_stamp, elapsed_days, elapsed_days0
   use glc_communicate,     only: my_task, master_task
   use glimmer_ncdf,        only: add_output, delete_output, nc_errorhandle
   use glc_broadcast,       only: broadcast_scalar
   use glimmer_ncio,        only: glimmer_nc_checkwrite, &
                                  glimmer_nc_createfile
   use glimmer_global,      only: fname_length
   use glc_constants
   use glc_kinds_mod
   use esmf,            only: ESMF_Clock
   use seq_timemgr_mod,     only: seq_timemgr_EClockGetData
   use shr_sys_mod
   use shr_kind_mod,        only: CL=>SHR_KIND_CL, CX=>SHR_KIND_CX, &
                                  IN=>SHR_KIND_IN
   use shr_file_mod,        only: shr_file_getunit, shr_file_freeunit
   use netcdf

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: glc_io_read_restart_time,  &
             glc_io_write_history,      &
             glc_io_write_restart

!----------------------------------------------------------------------
!
!   module variables
!
!----------------------------------------------------------------------

  character(CL), public ::  &
      history_vars  ! Name of the CISM variables to be output in cesm
                    ! history files

!EOP
!BOC
!EOC
!***********************************************************************
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: glc_io_read_restart_time
! !INTERFACE:

   subroutine glc_io_read_restart_time(nhour_glint, filename)

    use glc_files, only : ptr_filename

    implicit none
    integer(IN),             intent(inout) :: nhour_glint
    character(fname_length), intent(inout) :: filename

    ! local variables
    character(fname_length) :: filename0
    integer(IN)             :: cesmYMD           ! cesm model date
    integer(IN)             :: cesmTOD           ! cesm model sec
    integer(IN)             :: glcYMD            ! glc model date
    integer(IN)             :: glcTOD            ! glc model sec
    integer(IN)             :: rst_elapsed_days  ! 
    integer(IN)             :: ptr_unit          ! unit for pointer file
    integer(IN)             :: rst_unit          ! unit for restart file
    integer(IN)             :: status            !

!-----------------------------------------------------------------------

    if (my_task == master_task) then
       
       ! get restart filename from rpointer file
       ptr_unit = shr_file_getUnit()
       open(ptr_unit,file=ptr_filename)
       read(ptr_unit,'(a)') filename0
       filename = trim(filename0)
       close(ptr_unit)
       write(stdout,*) &
            'glc_io_read_restart_time: using dumpfile for restart = ', filename
       call shr_sys_flush(stdout)
       call shr_file_freeunit(ptr_unit)

       ! read time from the restart file, since glimmer needs this to initialize
       rst_unit = shr_file_getUnit()
       status = nf90_open(filename,0,rst_unit)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_get_att(rst_unit, NF90_GLOBAL, 'cesmYMD', cesmYMD)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_get_att(rst_unit, NF90_GLOBAL, 'cesmTOD', cesmTOD)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_get_att(rst_unit, NF90_GLOBAL, 'glcYMD', glcYMD)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_get_att(rst_unit, NF90_GLOBAL, 'glcTOD', glcTOD)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_get_att(rst_unit, NF90_GLOBAL, 'elapsed_days', rst_elapsed_days)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_close(rst_unit)
       call nc_errorhandle(__FILE__,__LINE__,status)

    end if

    call broadcast_scalar (cesmYMD         , master_task)
    call broadcast_scalar (cesmTOD         , master_task)
    call broadcast_scalar (glcYMD          , master_task)
    call broadcast_scalar (glcTOD          , master_task)
    call broadcast_scalar (rst_elapsed_days, master_task)

    ! calculate nhour_glint for return
    nhour_glint = rst_elapsed_days * 24

  end subroutine glc_io_read_restart_time

!***********************************************************************
!BOP
! !IROUTINE: glc_io_write_history
! !INTERFACE:

   subroutine glc_io_write_history(instance, EClock)

    use glint_type
    use glide_io, only : glide_io_create, glide_io_write
    use glint_io, only : glint_io_create, glint_io_write
    use glide_nc_custom, only: glide_nc_filldvars
    implicit none
    type(glint_instance), intent(inout) :: instance
    type(ESMF_Clock),     intent(in)    :: EClock

    ! local variables
    type(glimmer_nc_output),  pointer :: oc => null()
    character(CL) :: filename
    integer(IN)   :: cesmYMD           ! cesm model date
    integer(IN)   :: cesmTOD           ! cesm model sec
    integer(IN)   :: cesmYR            ! cesm model year
    integer(IN)   :: cesmMON           ! cesm model month
    integer(IN)   :: cesmDAY           ! cesm model day
    integer(IN)   :: glcYMD            ! cism model date
    integer(IN)   :: glcTOD            ! cism model sec
    integer(IN)   :: rst_elapsed_days  ! 
    integer(IN)   :: ptr_unit          ! unit for pointer file
    integer(IN)   :: status            !

!-----------------------------------------------------------------------

    ! figure out history filename
    call seq_timemgr_EClockGetData(EClock, curr_ymd=cesmYMD, curr_tod=cesmTOD, &
                                   curr_yr=cesmYR, curr_mon=cesmMON, curr_day=cesmDAY)
    filename = glc_filename(cesmYR, cesmMON, cesmDAY, cesmTOD, 'history')

    if (my_task == master_task) then
       write(stdout,*) &
            'glc_io_write_history: calling dumpfile for history filename= ', filename
       call shr_sys_flush(stdout)
    endif

    allocate(oc)
    oc%freq          = 1
    oc%append        = .false.
    oc%default_xtype = NF90_DOUBLE
    oc%nc%filename   = ''
    oc%nc%filename   = trim(filename)
!jw    oc%nc%vars       = ' acab artm thk usurf uvel vvel uflx vflx temp '
    oc%nc%vars       = trim(history_vars)
    oc%nc%vars_copy  = oc%nc%vars
!jw TO DO: fill out the rest of the metadata
!jw    oc%metadata%title =
!jw    oc%metadata%institution =
!jw    oc%metadata%source =
!jw    oc%metadata%history =
!jw    oc%metadata%references =
!jw    oc%metadata%comment = 

    ! create the output unit
    call glimmer_nc_createfile(oc, instance%model)
    call glide_io_create(oc, instance%model)
    call glint_io_create(oc, instance%model)

    if (my_task == master_task) then
       ! write time to the file
       glcYMD = iyear*10000 + imonth*100 + iday
       glcTOD = ihour*3600 + iminute*60 + isecond
       status = nf90_put_att(oc%nc%id, NF90_GLOBAL, 'cesmYMD', cesmYMD)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(oc%nc%id, NF90_GLOBAL, 'cesmTOD', cesmTOD)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(oc%nc%id, NF90_GLOBAL, 'glcYMD', glcYMD)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(oc%nc%id, NF90_GLOBAL, 'glcTOD', glcTOD)
       call nc_errorhandle(__FILE__,__LINE__,status)
       rst_elapsed_days = elapsed_days - elapsed_days0
       status = nf90_put_att(oc%nc%id, NF90_GLOBAL, 'elapsed_days', rst_elapsed_days)
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if
    
    call glide_nc_filldvars(oc, instance%model)
    call glimmer_nc_checkwrite(oc, instance%model, forcewrite=.true., &
                               time=instance%glide_time)
    call glide_io_write(oc, instance%model)
    call glint_io_write(oc, instance)

    if (my_task == master_task) then
       status = nf90_close(oc%nc%id)
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    oc => null()
!jw TO DO: figure out why deallocate statement crashes the code
!jw    deallocate(oc)

  end subroutine glc_io_write_history

!***********************************************************************
!BOP
! !IROUTINE: glc_io_write_restart
! !INTERFACE:

   subroutine glc_io_write_restart(instance, EClock)

    use glc_files, only : ptr_filename
    use glint_type
    use glide_io, only : glide_io_create, glide_io_write
    use glint_io, only : glint_io_create, glint_io_write
    use glide_nc_custom, only: glide_nc_filldvars
    implicit none
    type(glint_instance), intent(inout) :: instance
    type(ESMF_Clock),     intent(in)    :: EClock

    ! local variables
    type(glimmer_nc_output),  pointer :: oc => null()
    character(CL) :: filename
    integer(IN)   :: cesmYMD           ! cesm model date
    integer(IN)   :: cesmTOD           ! cesm model sec
    integer(IN)   :: cesmYR            ! cesm model year
    integer(IN)   :: cesmMON           ! cesm model month
    integer(IN)   :: cesmDAY           ! cesm model day
    integer(IN)   :: glcYMD            ! cism model date
    integer(IN)   :: glcTOD            ! cism model sec
    integer(IN)   :: rst_elapsed_days  ! 
    integer(IN)   :: ptr_unit          ! unit for pointer file
    integer(IN)   :: status            !

!-----------------------------------------------------------------------

    ! figure out restart filename
    call seq_timemgr_EClockGetData(EClock, curr_ymd=cesmYMD, curr_tod=cesmTOD, &
                                   curr_yr=cesmYR, curr_mon=cesmMON, curr_day=cesmDAY)
    filename = glc_filename(cesmYR, cesmMON, cesmDAY, cesmTOD, 'restart')

    if (my_task == master_task) then
       write(stdout,*) &
            'glc_io_write_restart: calling dumpfile for restart filename= ', filename
       call shr_sys_flush(stdout)
    endif

    allocate(oc)
    oc%freq          = 1
    oc%append        = .false.
    oc%default_xtype = NF90_DOUBLE
    oc%nc%filename   = ''
    oc%nc%filename   = trim(filename)
    oc%nc%vars       = ' restart '
    oc%nc%vars_copy  = oc%nc%vars
!jw TO DO: fill out the rest of the metadata
!jw    oc%metadata%title =
!jw    oc%metadata%institution =
!jw    oc%metadata%source =
!jw    oc%metadata%history =
!jw    oc%metadata%references =
!jw    oc%metadata%comment = 

    ! create the output unit
    call glimmer_nc_createfile(oc, instance%model)
    call glide_io_create(oc, instance%model)
    call glint_io_create(oc, instance%model)

    if (my_task == master_task) then
       ! write time to the file
       glcYMD = iyear*10000 + imonth*100 + iday
       glcTOD = ihour*3600 + iminute*60 + isecond
       status = nf90_put_att(oc%nc%id, NF90_GLOBAL, 'cesmYMD', cesmYMD)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(oc%nc%id, NF90_GLOBAL, 'cesmTOD', cesmTOD)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(oc%nc%id, NF90_GLOBAL, 'glcYMD', glcYMD)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(oc%nc%id, NF90_GLOBAL, 'glcTOD', glcTOD)
       call nc_errorhandle(__FILE__,__LINE__,status)
       rst_elapsed_days = elapsed_days - elapsed_days0
       status = nf90_put_att(oc%nc%id, NF90_GLOBAL, 'elapsed_days', rst_elapsed_days)
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if
    
    call glide_nc_filldvars(oc, instance%model)
    call glimmer_nc_checkwrite(oc, instance%model, forcewrite=.true., &
                               time=instance%glide_time)
    call glide_io_write(oc, instance%model)
    call glint_io_write(oc, instance)

    if (my_task == master_task) then
       status = nf90_close(oc%nc%id)
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    oc => null()
!jw TO DO: figure out why deallocate statement crashes the code
!jw    deallocate(oc)

    ! write pointer to restart file
    if (my_task == master_task) then
       ptr_unit = shr_file_getUnit()
       open(ptr_unit,file=ptr_filename)
       write(ptr_unit,'(a)') filename
       close(ptr_unit)
       call shr_file_freeunit(ptr_unit)
    endif

  end subroutine glc_io_write_restart

!***********************************************************************
! BOP
!
! !ROUTINE: glc_filename
!
! !INTERFACE:
  character(CL) function glc_filename( yr_spec, mon_spec, day_spec, sec_spec, file_type )
!
! !DESCRIPTION: Create a filename from a filename specifier. Interpret filename specifier
! string with:
! %c for case
! %i for instance suffix
! %y for year
! %m for month
! %d for day
! %s for second
! %% for the "%" character
! If the filename specifier has spaces " ", they will be trimmed out
! of the resulting filename.
!
! !USES:
    use glc_time_management, only: runid
    use glc_ensemble       , only: get_inst_suffix
!
! !INPUT/OUTPUT PARAMETERS:
  integer,      intent(in)  :: yr_spec         ! Simulation year
  integer,      intent(in)  :: mon_spec        ! Simulation month
  integer,      intent(in)  :: day_spec        ! Simulation day
  integer,      intent(in)  :: sec_spec        ! Seconds into current simulation day
  character(7), intent(in)  :: file_type       ! file type, either history or restart
!
! EOP
!
  integer       :: i, n           ! Loop variables
  integer       :: year           ! Simulation year
  integer       :: month          ! Simulation month
  integer       :: day            ! Simulation day
  integer       :: ncsec          ! Seconds into current simulation day
  character(CX) :: string         ! Temporary character string
  character(CL) :: format         ! Format character string
  character(CL) :: filename_spec  ! cism filename specifier

  !---------------------------------------------------------------------------
  ! Determine what the file tpye is and set the filename specifier accordingly
  !---------------------------------------------------------------------------

  filename_spec = ' '
  if (file_type.eq.'history') then
     filename_spec = '%c.cism%i.h.%y-%m-%d-%s'
  else if (file_type.eq.'restart') then
     filename_spec = '%c.cism%i.r.%y-%m-%d-%s'
  else
     call shr_sys_abort ('glc_filename: file_type specifier is invalid')
  endif

  !-----------------------------------------------------------------
  ! Determine year, month, day and sec to put in filename
  !-----------------------------------------------------------------

 if ( len_trim(filename_spec) == 0 )then
     call shr_sys_abort ('glc_filename: filename specifier is empty')
  end if
  if ( index(trim(filename_spec)," ") /= 0 )then
     call shr_sys_abort ('glc_filename: filename specifier can not contain a space:'//trim(filename_spec))
  end if

  year  = yr_spec
  month = mon_spec
  day   = day_spec
  ncsec = sec_spec

  ! Go through each character in the filename specifier and interpret if special string

  i = 1
  glc_filename = ''
  string = ''
  do while ( i <= len_trim(filename_spec) )
     if ( filename_spec(i:i) == "%" )then
        i = i + 1
        select case( filename_spec(i:i) )
        case( 'c' )   ! runid
           string = trim(runid)
        case( 'i' )   ! instance suffix
           call get_inst_suffix(string)
        case( 'y' )   ! year
           if ( year > 99999   ) then
              format = '(i6.6)'
           else if ( year > 9999    ) then
              format = '(i5.5)'
           else
              format = '(i4.4)'
           end if
           write(string,format) year
        case( 'm' )   ! month
           write(string,'(i2.2)') month
        case( 'd' )   ! day
           write(string,'(i2.2)') day
        case( 's' )   ! second
           write(string,'(i5.5)') ncsec
        case( '%' )   ! percent character
           string = "%"
        case default
           call shr_sys_abort ('glc_filename: Invalid expansion character: '//filename_spec(i:i))
        end select
     else
       n = index( filename_spec(i:), "%" )
        if ( n == 0 ) n = len_trim( filename_spec(i:) ) + 1
        if ( n == 0 ) exit
        string = filename_spec(i:n+i-2)
        i = n + i - 2
     end if
     if ( len_trim(glc_filename) == 0 )then
        glc_filename = trim(string)
     else
        if ( (len_trim(glc_filename)+len_trim(string)) >= CL )then
           call shr_sys_abort ('glc_filename Resultant filename too long')
        end if
        glc_filename = trim(glc_filename) // trim(string)
     end if
     i = i + 1
  end do
  if ( len_trim(glc_filename) == 0 )then
     call shr_sys_abort ('glc_filename: Resulting filename is empty')
  end if

  ! add ".nc" to tail end
  glc_filename = trim(glc_filename) // '.nc'

end function glc_filename

end module glc_io
