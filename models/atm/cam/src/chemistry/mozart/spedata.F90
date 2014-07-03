module spedata
!----------------------------------------------------------------------- 
!
! BOP
!
! !MODULE: spedata
!
! !DESCRIPTION
! Handles reading and interpolating solar proton ionization data.
!
! !USES
  use shr_kind_mod, only: r8 => shr_kind_r8,r4 => shr_kind_r4, shr_kind_cl
  use shr_cal_mod,  only: shr_cal_gregorian
  use time_manager, only: get_curr_date, get_curr_calday, get_step_size, timemgr_is_caltype
  use spmd_utils,   only: masterproc
  use ppgrid,       only: pcols, pver, begchunk, endchunk
  use abortutils,   only: endrun
  use pio,          only: var_desc_t, file_desc_t, pio_get_var, pio_get_att, &
                          pio_setdebuglevel, pio_seterrorhandling, pio_bcast_error, &
                          pio_internal_error, pio_noerr, pio_inq_varid, pio_put_var, &
                          pio_put_att, pio_inq_dimid, pio_char, pio_def_dim, pio_def_var, &
                          pio_inq_dimlen, pio_closefile
  use cam_pio_utils,only: cam_pio_openfile
  use perf_mod,     only: t_startf, t_stopf
  use cam_logfile,  only: iulog

#if ( defined SPMD )
  use mpishorthand, only: mpicom, mpir8, mpiint, mpichar
#endif

  implicit none

  private  ! all unless made public
  save 

! !PUBLIC MEMBERS

  public spedata_init            ! subroutine to open files, allocate blocked arrays, etc
  public advance_spedata         ! subroutine to read more data and interpolate
  public get_ionpairs_profile    ! interface to get ionization profile
  public spedata_defaultopts
  public spedata_setopts

  !------------------------------------------------------------------
  ! Interface to access the meteorology fields.  Possible invocations
  ! are as follows:
  !   call get_met_fields( physics_state, us, vs , tend, dt )
  !   call get_met_fields( u, v )
  !   call get_met_fields( cam_in )
  !------------------------------------------------------------------
  Interface get_ionpairs_profile                        ! overload accessors
     Module Procedure get_ionpairs
  End Interface
  
  logical :: remove_spe_file = .false.  ! delete data file when finished with it

! !REVISION HISTORY:
!   25 Oct 2005  Francis Vitt     Creation
!
! EOP
!----------------------------------------------------------------------- 
! $Id: spedata.F90,v 1.1.2.1 2006/05/03 20:53:09 stacy Exp $
! $Author: stacy $
!----------------------------------------------------------------------- 

  type input_profile
     real(r8), dimension(:), pointer :: data
  endtype input_profile


  real(r8), allocatable :: ionpairs(:)  ! interpolated ionization profile
  real(r8), allocatable :: prod_pressures(:) 

  type(input_profile) :: ionpairs_i(2)

  integer :: prod_id ! var id of the data in the netCDF
  integer :: nprod_press

  integer :: dateid           ! var id of the date in the netCDF
  integer :: secid            ! var id of the sec data 
  real(r8) :: datatimem = -1.e36_r8     ! time of prv. values read in
  real(r8) :: datatimep = -1.e36_r8     ! time of nxt. values read in

  integer, parameter :: nm=1    ! array index for previous (minus) data
  integer, parameter :: np=2    ! array indes for next (plus) data

  real(r8) :: curr_mod_time ! model time - calendar day
  real(r8) :: next_mod_time ! model time - calendar day - next time step

  character(len=shr_kind_cl) :: curr_filename = ' '
  character(len=shr_kind_cl) :: next_filename = ' '
  character(len=shr_kind_cl) :: spedata_file = ' '
  type(file_desc_t) :: curr_fileid, next_fileid     ! the pio id of the NetCDF file
  type(var_desc_t), pointer :: currfnameid, nextfnameid  ! the id of the filename strings in the restart file
  real(r8), pointer, dimension(:)  :: curr_data_times, next_data_times
  character(len=shr_kind_cl) :: filenames_list = ''

  character(len=16) :: calendar

  logical :: spe_run     = .false.

contains

!--------------------------------------------------------------------------
! Get the default runtime options
!--------------------------------------------------------------------------
  subroutine spedata_defaultopts( spe_data_file_out, &
                                  spe_remove_file_out , &
                                  spe_filenames_list_out )

    implicit none

    character(len=shr_kind_cl), intent(out), optional :: spe_data_file_out
    character(len=shr_kind_cl), intent(out), optional :: spe_filenames_list_out
    logical, intent(out), optional :: spe_remove_file_out

    if ( present( spe_data_file_out ) ) then
       spe_data_file_out = spedata_file
    endif

    if ( present( spe_remove_file_out ) ) then
       spe_remove_file_out = remove_spe_file
    endif

    if ( present( spe_filenames_list_out ) ) then
       spe_filenames_list_out =  filenames_list
    endif

  end subroutine spedata_defaultopts

!--------------------------------------------------------------------------
! Set runtime options
!--------------------------------------------------------------------------
  subroutine spedata_setopts( spe_data_file_in, &
                              spe_remove_file_in, &
                              spe_filenames_list_in ) 

    implicit none

    character(len=shr_kind_cl), intent(in), optional :: spe_data_file_in
    character(len=shr_kind_cl), intent(in), optional :: spe_filenames_list_in
    logical, intent(in), optional :: spe_remove_file_in

    integer :: ierr

    if ( present( spe_data_file_in ) ) then
       spedata_file = spe_data_file_in 
    endif

    if ( present( spe_remove_file_in ) ) then
       remove_spe_file = spe_remove_file_in
    endif

    if ( present( spe_filenames_list_in ) ) then
       filenames_list = spe_filenames_list_in 
    endif

    if (len_trim(spedata_file) > 0) then 
      spe_run = .true.
    else
      spe_run = .false.
    endif
      
    if (masterproc) then
       write(iulog,*) 'This run includes NOx and HOx production from solor proton events: ',spe_run
       if ( spe_run ) then
         write(iulog,*) 'Time-variant solar proton ionization dataset (spedata_file) is: ', trim(spedata_file)
         write(iulog,*) 'Solar proton ionization data file will be removed (remove_spe_file): ', remove_spe_file
         write(iulog,*) 'Solar proton ionization data file names list file: ', trim(filenames_list) 
       endif
    endif

  end subroutine spedata_setopts

  subroutine spedata_init()
!--------------------------------------------------------------------------
! Opens file, allocates arrays
!--------------------------------------------------------------------------

  use cam_control_mod, only : nsrest
  use string_utils,    only : to_upper

    implicit none

!--------------------------------------------------------------------------
! 	... local variables
!--------------------------------------------------------------------------
    integer :: astat
    integer :: ierr

    if (.not. spe_run) return

!--------------------------------------------------------------------------
! 	... initial run ?
!--------------------------------------------------------------------------
!    if (nsrest == 0) then
       curr_filename = spedata_file
       next_filename = ''
!    endif

    call open_spe_datafile( curr_filename, curr_fileid, curr_data_times, read_pressures=.true. )

    if ( len_trim(next_filename) > 0 ) &
         call open_spe_datafile( next_filename, next_fileid, next_data_times )
    
    ierr = pio_inq_varid( curr_fileid, 'Prod', prod_id )
    

!--------------------------------------------------------------------------
! allocate space for data arrays ...
!--------------------------------------------------------------------------
    allocate( ionpairs(nprod_press),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'spedata_init: failed to allocate ionpairs; error = ',astat
       call endrun
    end if
    allocate( ionpairs_i(nm)%data(nprod_press),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'spedata_init: failed to allocate ionpairs_i; error = ',astat
       call endrun
    end if
    allocate( ionpairs_i(np)%data(nprod_press),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'spedata_init: failed to allocate ionpairs_i; error = ',astat
       call endrun
    end if

 end subroutine spedata_init


!-----------------------------------------------------------------------
! Reads more data if needed and interpolates data to current model time 
!-----------------------------------------------------------------------
 subroutine advance_spedata()

    implicit none

    if (.not. spe_run) return

    call t_startf('MET__advance')

    call get_model_time()

    if ( ( curr_mod_time > datatimep ) ) then
       call check_files()
    endif

    if ( curr_mod_time > datatimep ) then
       call read_next_spedata()
    end if

! need to inperpolate the data, regardless !
! each mpi tasks needs to interpolate
    call interpolate_spedata()

    call t_stopf('MET__advance')

  end subroutine advance_spedata

!------------------------------------------------------------------------
! get the ion pair production profile
!------------------------------------------------------------------------
  subroutine get_ionpairs( model_pressures, pairs )
    
    use interpolate_data, only : lininterp

    implicit none

    real(r8), intent(in)  :: model_pressures(:)
    real(r8), intent(out) :: pairs(:)

    integer :: npress, ub,lb

    pairs(:) = 0._r8

    if (.not. spe_run) return

    npress = size(model_pressures)

    ub = ubound( prod_pressures,1 )
    lb = lbound( prod_pressures,1 )

    ! interpolate to model levels
    call lininterp( ionpairs(ub:lb:-1), prod_pressures(ub:lb:-1)*1.e2_r8, nprod_press, pairs, model_pressures, npress )

  end subroutine get_ionpairs

!------------------------------------------------------------------------------
! internal methods :
!------------------------------------------------------------------------------

  subroutine get_model_time()

    implicit none

    integer yr, mon, day, ncsec  ! components of a date

    call t_startf('MET__get_model_time')

    call get_curr_date(yr, mon, day, ncsec)

    curr_mod_time = get_time_float( yr, mon, day, ncsec )
    next_mod_time = curr_mod_time + get_step_size()/86400._r8

    call t_stopf('MET__get_model_time')

  end subroutine get_model_time

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  subroutine check_files()

    use shr_sys_mod, only: shr_sys_system
    implicit none

!-----------------------------------------------------------------------
! 	... local variables
!-----------------------------------------------------------------------
    character(len=128) :: ctmp
    integer ::  istat


    if (next_mod_time > curr_data_times(size(curr_data_times))) then
       if ( .not. associated(next_data_times) ) then
          ! open next file...
          next_filename = incr_filename( curr_filename )
          call open_spe_datafile( next_filename, next_fileid, next_data_times )
       endif
    endif

    if ( associated(next_data_times) ) then
       if (curr_mod_time >= next_data_times(1)) then

          ! close current file ...
          call pio_closefile( curr_fileid )
          ! remove if requested
          if(masterproc) then
             if( remove_spe_file ) then
                write(iulog,*) 'check_files: removing file = ',trim(curr_filename) 
                ctmp = 'rm -f ' // trim(curr_filename) 
                write(iulog,*) 'check_files: fsystem issuing command - '
                write(iulog,*) trim(ctmp)
                call shr_sys_system( ctmp, istat )
             end if
          endif

          curr_filename = next_filename
          curr_fileid = next_fileid

          deallocate( curr_data_times )
          allocate( curr_data_times( size( next_data_times ) ) )
          curr_data_times(:) = next_data_times(:)

          next_filename = ''

          deallocate( next_data_times )
          nullify(  next_data_times )

       endif
    endif

  end subroutine check_files

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  function incr_filename( filename )

    !-----------------------------------------------------------------------
    ! 	... Increment or decrement a date string withing a filename
    !           the filename date section is assumed to be of the form
    !           yyyy-dd-mm
    !-----------------------------------------------------------------------

    use string_utils,  only : incstr

    implicit none


    character(len=*), intent(in) :: filename ! present dynamical dataset filename
    character(len=shr_kind_cl) :: incr_filename      ! next filename in the sequence

    ! set new next_filename ...

    !-----------------------------------------------------------------------
    !	... local variables
    !-----------------------------------------------------------------------
    integer :: pos, pos1, istat
    character(len=shr_kind_cl) :: fn_new, line
    character(len=6)   :: seconds
    character(len=5)   :: num

    integer :: ios

    if ( len_trim(filenames_list) .eq. 0) then

       ! increment the number in the filename to the next number
       pos = len_trim( filename )
       fn_new = filename(:pos)
       write(iulog,*) 'incr_flnm: old filename = ',trim(fn_new)
       if( fn_new(pos-2:) == '.nc' ) then
          pos = pos - 3
       end if
       istat = incstr( fn_new(:pos), 1 )
       if( istat /= 0 ) then
          write(iulog,*) 'incr_flnm: incstr returned ', istat
          write(iulog,*) '           while trying to decrement ',trim( fn_new )
          call endrun
       end if

    else

       ! open filenames_list
       write(iulog,*) 'incr_flnm: old filename = ',trim(filename)
       write(iulog,*) 'incr_flnm: open filenames_list : ',filenames_list 
       open( unit=9, file=filenames_list, iostat=ios, status="OLD")
       if (ios /= 0) then
          call endrun('not able to open filenames_list file: '//filenames_list)
       endif

       ! read file names
       read( unit=9, fmt='(A)', iostat=ios ) line 
       if (ios /= 0) then
          call endrun('not able to increment file name from filenames_list file: '//filenames_list)
       endif
       do while( trim(line) /= trim(filename) )
          read( unit=9, fmt='(A)', iostat=ios ) line 
          if (ios /= 0) then
             call endrun('not able to increment file name from filenames_list file: '//filenames_list)
          endif
       enddo

       read( unit=9, fmt='(A)', iostat=ios ) line 
       if (ios /= 0) then
          call endrun('not able to increment file name from filenames_list file: '//filenames_list)
       endif
       fn_new = trim(line)

       close(unit=9)

    endif

    incr_filename = trim(fn_new)
    write(iulog,*) 'incr_flnm: new filename = ',incr_filename

  end function incr_filename

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  subroutine find_times( itms, fids, datatm, datatp, time )

    implicit none

!------------------------------------------------------------------------------
!	... dummy arguments
!------------------------------------------------------------------------------
    integer, intent(out) :: itms(2) ! record numbers that bracket time
    type(file_desc_t), intent(out) :: fids(2) ! ids of files that contains these recs
    real(r8), intent(in) :: time    ! time of interest
    real(r8), intent(out):: datatm, datatp

!------------------------------------------------------------------------------
!	... local variables
!------------------------------------------------------------------------------
    integer np1        ! current forward time index of dataset
    integer n,i      ! 
    integer :: curr_tsize, next_tsize, all_tsize
    real(r8), allocatable, dimension(:):: all_data_times
    logical :: found_time

    curr_tsize = size(curr_data_times)
    next_tsize = 0
    if ( associated(next_data_times)) next_tsize = size(next_data_times)

    all_tsize = curr_tsize + next_tsize

    allocate( all_data_times( all_tsize ) )

    all_data_times(:curr_tsize) = curr_data_times(:)
    if (next_tsize > 0) all_data_times(curr_tsize+1:all_tsize) = next_data_times(:)

    ! find bracketing times 
    found_time = .false.
    do n = 1,all_tsize-1
       np1 = n + 1
       datatm = all_data_times(n)
       datatp = all_data_times(np1)
       if ( (time .ge. datatm) .and. (time .le. datatp) ) then
          found_time = .true.
          exit
       endif
    enddo
  
    if( found_time ) then
       deallocate( all_data_times )
  
       itms(1) = n
       itms(2) = np1
       fids(:) = curr_fileid
  
       do i=1,2
          if ( itms(i) > curr_tsize ) then 
             itms(i) = itms(i) - curr_tsize 
             fids(i) = next_fileid
          endif
       enddo
    else
       write(iulog,*)'FIND_TIMES: Failed to find dates bracketing desired time =', time
       write(iulog,*)' datatm = ',datatm
       write(iulog,*)' datatp = ',datatp
       write(iulog,*)' all_data_times = ',all_data_times
       call endrun
    end if

  end subroutine find_times

!------------------------------------------------------------------------
!------------------------------------------------------------------------
  subroutine read_next_spedata()
    implicit none

!------------------------------------------------------------------------
!	... local variables
!------------------------------------------------------------------------
    integer :: i
    integer :: recnos(2)
    type(file_desc_t) :: fids(2)
    integer :: cnt(2)            ! array of counts for each dimension
    integer :: strt(2)           ! array of starting indices
    integer :: ierr
    call t_startf('MET__read_next_spedata')

    call find_times( recnos, fids, datatimem, datatimep, curr_mod_time )

    cnt(1)  = nprod_press
    cnt(2)  = 1
    strt(1) = 1

    do i = 1,2
       strt(2) = recnos(i)
       ierr = pio_get_var( fids(i), prod_id, strt, cnt,  ionpairs_i(i)%data )
    enddo

    if (masterproc) write(iulog,*)'READ_NEXT_SPEDATA: Read soloar proton ionization data '

    call t_stopf('MET__read_next_spedata')

  end subroutine read_next_spedata


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  subroutine interpolate_spedata()

    implicit none

    real(r4) fact1, fact2
    real(r8) deltat 

    call t_startf('MET__interpolate_spedata')

    deltat = datatimep - datatimem

    fact1 = (datatimep - curr_mod_time)/deltat
    fact2 = 1._r8-fact1

    ionpairs(:) = fact1*ionpairs_i(nm)%data(:) + fact2*ionpairs_i(np)%data(:)

    call t_stopf('MET__interpolate_spedata')

  end subroutine interpolate_spedata

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine get_dimension( fid, dname, dsize )
    implicit none
    type(file_desc_t), intent(in) :: fid
    character(*), intent(in) :: dname
    integer, intent(out) :: dsize

    integer :: dimid, ierr

    ierr = pio_inq_dimid( fid, dname, dimid )
    ierr = pio_inq_dimlen( fid, dimid, dsize )


  end subroutine get_dimension

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine open_spe_datafile( fname, fileid, times, read_pressures )

    use ioFileMod, only : getfil

    implicit none

    character(*), intent(in) :: fname
    type(file_desc_t), intent(inout) :: fileid
    real(r8), pointer :: times(:)
    logical, optional, intent(in) :: read_pressures

    character(len=shr_kind_cl) :: filen   
    integer :: year, month, day, dsize, i, timesize
    integer :: dateid,secid,press_id
    integer, allocatable , dimension(:) :: dates, datesecs
    integer :: ierr
    !
    ! open file and get fileid
    !
    call getfil( fname, filen, 0 )
    call cam_pio_openfile( fileid, filen, 0 )
    if(masterproc)   write(iulog,*)'open_met_datafile: ',trim(filen)

    call get_dimension( fileid, 'time', timesize )

    if ( associated(times) ) deallocate(times)
    allocate( times(timesize) )

    allocate( dates(timesize) )
    allocate( datesecs(timesize) )

    ierr = pio_inq_varid( fileid, 'date',    dateid  )
    ierr = pio_inq_varid( fileid, 'datesec', secid  )

    ierr = pio_get_var( fileid, dateid, dates )
    ierr = pio_get_var( fileid, secid,  datesecs  )

    do i=1,timesize
       year = dates(i) / 10000
       month = mod(dates(i),10000)/100
       day = mod(dates(i),100)
       times(i) = get_time_float( year, month, day, datesecs(i) )
    enddo

    deallocate( dates )
    deallocate( datesecs )   

    if (present( read_pressures ) ) then
       call get_dimension( fileid, 'pressure', nprod_press )

       allocate( prod_pressures( nprod_press ) )

       ierr = pio_inq_varid( fileid, 'pressure',    press_id  )
       ierr = pio_get_var( fileid, press_id, (/ 1 /), (/ nprod_press /), prod_pressures )
    endif    

  end subroutine open_spe_datafile

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  function get_time_float( year, month, day, sec )


! returns float representation of time -- number of days
! since 1 jan 0001 00:00:00.000

    implicit none

    integer, intent(in) :: year, month, day
    integer, intent(in) :: sec
    real(r8) :: get_time_float

! ref date is 1 jan 0001

    integer  :: refyr, refmn, refdy
    real(r8) :: refsc, fltdy
    integer  :: doy(12)

!              jan feb mar apr may jun jul aug sep oct nov dec
!              31  28  31  30  31  30  31  31  31  31  30  31
    data doy /  1, 32, 60, 91,121,152,182,213,244,274,305,335 /

    refyr = 1
    refmn = 1
    refdy = 1
    refsc = 0._r8

    if ( timemgr_is_caltype(trim(shr_cal_gregorian))) then
       fltdy = greg2jday(year, month, day) - greg2jday(refyr,refmn,refdy)
    else ! assume no_leap (all years are 365 days)
       fltdy = (year - refyr)*365._r8 + &
               (doy(month)-doy(refmn)) + &
               (day-refdy) 
    endif

    get_time_float = fltdy + ((sec-refsc)/86400._r8)

  endfunction get_time_float

!-----------------------------------------------------------------------
! 	... Return Julian day number given Gregorian date.
!
! Algorithm from Hatcher,D.A., Simple Formulae for Julian Day Numbers
! and Calendar Dates, Q.Jl.R.astr.Soc. (1984) v25, pp 53-55.
!-----------------------------------------------------------------------
  function greg2jday( year, month, day )

    implicit none

    integer, intent(in) :: year, month, day
    integer :: greg2jday

    !-----------------------------------------------------------------------
    !	... Local variables
    !-----------------------------------------------------------------------
    integer :: ap, mp
    integer :: y, d, n, g

    !-----------------------------------------------------------------------
    !     	... Modify year and month numbers
    !-----------------------------------------------------------------------
    ap = year - (12 - month)/10
    mp = MOD( month-3,12 )
    if( mp < 0 ) then
       mp = mp + 12
    end if

    !-----------------------------------------------------------------------
    !     	... Julian day
    !-----------------------------------------------------------------------
    y = INT( 365.25_r8*( ap + 4712 ) )
    d = INT( 30.6_r8*mp + .5_r8 )
    n = y + d + day  + 59
    g = INT( .75_r8*INT( ap/100 + 49 ) ) - 38
    greg2jday = n - g

  end function greg2jday

end module spedata
