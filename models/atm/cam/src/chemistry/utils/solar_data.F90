!-----------------------------------------------------------------------
! Solar irradiance / photon flux data
!-----------------------------------------------------------------------
module solar_data
  use shr_kind_mod, only: r8 => shr_kind_r8
  use spmd_utils,   only: masterproc
  use abortutils,   only: endrun
  use pio
  use cam_pio_utils, only : cam_pio_openfile
  use cam_logfile,  only: iulog
  use infnan,       only: nan, assignment(=)

  implicit none

  save
  private
  public :: solar_data_readnl
  public :: solar_data_init
  public :: solar_data_advance

  public :: nbins, we  ! number of wavelength samples of spectrum, wavelength endpoints
  public :: sol_etf
  public :: sol_irrad
  public :: sol_tsi
  public :: do_spctrl_scaling
  public :: has_spectrum
  public :: has_ref_spectrum
  public :: ssi_ref
  public :: ref_tsi

  integer :: nbins
  integer :: ntimes
  real(r8), allocatable :: sol_etf(:)
  real(r8), allocatable :: irradi(:,:)
  real(r8)              :: itsi(2)
  real(r8), allocatable :: ssi_ref(:)  ! a reference spectrum constructed from 3 solar cycles of data

  real(r8), allocatable :: we(:)
  real(r8), allocatable :: data_times(:)

  integer :: last_index = 1
  type(file_desc_t) :: file_id
  integer :: ssi_vid
  integer :: tsi_vid
  integer :: ref_vid
  integer :: tsi_ref_vid

  logical :: initialized = .false.
  logical :: has_spectrum = .false.
  logical :: has_ref_spectrum = .false.
  logical :: has_tsi = .false.

  real(r8), allocatable :: sol_irrad(:)
  real(r8)              :: sol_tsi = -1.0_r8
  real(r8)              :: ref_tsi

  real(r8), allocatable :: irrad_fac(:)
  real(r8), allocatable :: etf_fac(:)
  real(r8), allocatable :: dellam(:)

  logical  :: fixed_scon
  logical  :: fixed_solar
  real(r8) :: offset_time

  logical            :: do_spctrl_scaling = .false.

! namelist vars
  character(len=256) :: solar_data_file = ''
  character(len=8)   :: solar_data_type = 'SERIAL' ! "FIXED" or "SERIAL"
  integer            :: solar_data_ymd = 0         ! YYYYMMDD for "FIXED" type
  integer            :: solar_data_tod = 0         ! seconds of day for "FIXED" type
  real(r8)           :: solar_const = -9999._r8         ! constant TSI (W/m2)
  logical            :: solar_htng_spctrl_scl = .false. ! do rad heating spectral scaling

contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine solar_data_readnl( nlfile )
    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit

    ! arguments
    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! local vars
    integer :: unitn, ierr

    namelist /solar_inparm/ &
         solar_data_file, solar_data_type, solar_data_ymd, solar_data_tod, solar_const, &
         solar_htng_spctrl_scl

    unitn = getunit()
    open( unitn, file=trim(nlfile), status='old' )
    call find_group_name(unitn, 'solar_inparm', status=ierr)
    if (ierr == 0) then
       read(unitn, solar_inparm, iostat=ierr)
       if (ierr /= 0) then
          call endrun('solar_data_readnl: ERROR reading namelist')
       end if
    end if
    close(unitn)
    call freeunit(unitn)

    if ( (len_trim(solar_data_file) > 0) .and. (solar_const>0._r8) ) then
       call endrun('solar_data_readnl: ERROR cannot specify both solar_data_file and solar_const')      
    endif

    if ( len_trim(solar_data_file) > 0 ) then
       fixed_scon = .false.
    else
       fixed_scon = .true.
    endif
    
    if ( solar_const>0._r8 ) then
       sol_tsi = solar_const
    endif
    ref_tsi = nan
    
    if (masterproc) then
       write(iulog,*) 'solar_data_readnl: solar_const (W/m2) = ', solar_const
       if ( .not.fixed_scon ) then
          write(iulog,*) 'solar_data_readnl: solar_data_file = ',trim(solar_data_file)
          write(iulog,*) 'solar_data_readnl: solar_data_type = ',trim(solar_data_type)
          write(iulog,*) 'solar_data_readnl: solar_data_ymd  = ',solar_data_ymd
          write(iulog,*) 'solar_data_readnl: solar_data_tod  = ',solar_data_tod
       endif
    endif

  end subroutine solar_data_readnl

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine solar_data_init

    use ioFileMod, only : getfil
    use physconst, only : c0, planck

    integer :: astat, dimid, vid
    character(len=256) :: filen   
    real(r8), allocatable :: lambda(:)
    integer,  allocatable :: dates(:)
    integer,  allocatable :: datesecs(:)

    integer :: i, wvl_vid
    real(r8), parameter :: c = c0     ! speed of light (m/s)
    real(r8), parameter :: h = planck ! Planck's constant (Joule sec)
    real(r8), parameter :: fac = 1._r8/(h*c)

    real(r8) :: model_time, time

    integer :: ierr

    has_spectrum = .false.
    if ( fixed_scon ) return

    fixed_solar = trim(solar_data_type) == 'FIXED'

    call getfil( solar_data_file, filen, 0 )
    call cam_pio_openfile( file_id, filen, 0 )
    if(masterproc)   write(iulog,*)'solar_data_init: data file = ',trim(filen)
    call pio_seterrorhandling(file_id, pio_bcast_error)
    ierr = pio_inq_varid( file_id, 'ssi', ssi_vid )
    has_spectrum = ierr==PIO_NOERR

    ierr = pio_inq_varid( file_id, 'tsi', tsi_vid )
    has_tsi      = ierr==PIO_NOERR .and. solar_const<0._r8

    ierr = pio_inq_varid( file_id, 'ssi_ref', ref_vid )
    has_ref_spectrum = ierr==PIO_NOERR
    call pio_seterrorhandling(file_id, pio_internal_error)

    ierr = pio_inq_dimid( file_id, 'time', dimid )
    ierr = pio_inq_dimlen( file_id, dimid, ntimes )

    if ( has_spectrum ) then
       call pio_seterrorhandling(file_id, pio_bcast_error)
       ierr = pio_inq_varid( file_id, 'wavelength', wvl_vid )
       call pio_seterrorhandling(file_id, pio_internal_error)
       
       if ( ierr==PIO_NOERR ) then
          ierr = pio_inq_dimid( file_id, 'wavelength', dimid )
       else ! for backwards compatibility
          ierr = pio_inq_varid( file_id, 'wvl', wvl_vid  )
          ierr = pio_inq_dimid( file_id, 'wvl', dimid )
       endif
       ierr = pio_inq_dimlen( file_id, dimid, nbins )
       if ( has_ref_spectrum ) then
          ierr = pio_inq_varid( file_id, 'tsi_ref', tsi_ref_vid )
       endif
    endif

    do_spctrl_scaling = has_spectrum .and. solar_htng_spctrl_scl

    allocate(lambda(nbins), stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'solar_data_init: failed to allocate lambda; error = ',astat
       call endrun('solar_data_init')
    end if
    allocate(dellam(nbins), stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'solar_data_init: failed to allocate dellam; error = ',astat
       call endrun('solar_data_init')
    end if
    allocate(data_times(ntimes), stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'solar_data_init: failed to allocate data_times; error = ',astat
       call endrun('solar_data_init')
    end if

    allocate(dates(ntimes), stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'solar_data_init: failed to allocate dates; error = ',astat
       call endrun('solar_data_init')
    end if

    allocate(datesecs(ntimes), stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'solar_data_init: failed to allocate datesecs; error = ',astat
       call endrun('solar_data_init')
    end if

    allocate(irrad_fac(nbins), stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'solar_data_init: failed to allocate irrad_fac; error = ',astat
       call endrun('solar_data_init')
    end if
    allocate(etf_fac(nbins), stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'solar_data_init: failed to allocate etf_fac; error = ',astat
       call endrun('solar_data_init')
    end if
    allocate(sol_irrad(nbins), stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'solar_data_init: failed to allocate sol_irrad; error = ',astat
       call endrun('solar_data_init')
    end if
    allocate(ssi_ref(nbins), stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'solar_data_init: failed to allocate ssi_ref; error = ',astat
       call endrun('solar_data_init')
    end if
    ssi_ref(:) = nan

    ierr = pio_inq_varid( file_id, 'date', vid  )
    ierr = pio_get_var( file_id, vid, dates )
    call pio_seterrorhandling(file_id, pio_bcast_error)
    ierr = pio_inq_varid( file_id, 'datesec', vid )
    call pio_seterrorhandling(file_id, pio_internal_error)
    if (ierr==PIO_NOERR) then
       ierr = pio_get_var( file_id, vid, datesecs )
    else
       datesecs(:) = 0
    endif
       
    if (has_spectrum) then
       ierr = pio_get_var( file_id, wvl_vid, lambda )
       ierr = pio_inq_varid( file_id, 'band_width', vid  )
       ierr = pio_get_var( file_id, vid, dellam )
    endif

    if(masterproc)   write(iulog,*)'solar_data_init: has_ref_spectrum',has_ref_spectrum
    if ( has_ref_spectrum ) then
       ierr = pio_inq_varid( file_id, 'ssi_ref', vid  )
       ierr = pio_get_var( file_id, vid, ssi_ref )
       ierr = pio_get_var( file_id, tsi_ref_vid, ref_tsi )
    endif

    if ( has_spectrum ) then
       allocate(sol_etf(nbins), stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'solar_data_init: failed to allocate sol_etf; error = ',astat
          call endrun('solar_data_init')
       end if
       allocate(irradi(nbins,2), stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'solar_data_init: failed to allocate irradi; error = ',astat
          call endrun('solar_data_init')
       end if

       allocate(we(nbins+1), stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'solar_data_init: failed to allocate we; error = ',astat
          call endrun('solar_data_init')
       end if

       we(:nbins)  = lambda(:nbins) - 0.5_r8*dellam(:nbins)
       we(nbins+1) = lambda(nbins)  + 0.5_r8*dellam(nbins)
       do i = 1,nbins
          irrad_fac(i) = 1.e-3_r8                ! mW/m2/nm --> W/m2/nm
          etf_fac(i)   = 1.e-16_r8*lambda(i)*fac ! mW/m2/nm --> photons/cm2/sec/nm
       enddo
       if(has_ref_spectrum) then
          ssi_ref = ssi_ref * 1.e-3_r8        ! mW/m2/nm --> W/m2/nm
       endif
    endif

    offset_time = 0._r8

    if ( solar_data_ymd > 0 ) then
      call get_model_time( model_time )
      call convert_date( solar_data_ymd, solar_data_tod, time )
      offset_time = time - model_time
    endif

    call convert_dates( dates, datesecs, data_times )

    data_times = data_times - offset_time

    deallocate(lambda)
    deallocate(dates)

    ! need to force data loading when the model starts at a time =/ 00:00:00.000
    ! -- may occur in restarts also
    call solar_data_advance()
    initialized = .true.

  end subroutine solar_data_init

!-----------------------------------------------------------------------
! Reads in the ETF data for the current date.  
!-----------------------------------------------------------------------
  subroutine solar_data_advance( )

    use physconst,   only : cday

    integer  :: year, month, day, sec
    integer  :: index, i, nt
    integer  :: offset(2), count(2)
    logical  :: do_adv, read_data
    real(r8) :: time, delt
    real(r8) :: data(nbins)
    integer  :: ierr

    if ( fixed_scon ) return
    if ( fixed_solar .and. initialized ) return

    index = -1
    call get_model_time( time, year=year, month=month, day=day, seconds=sec )
    read_data = time > data_times(last_index) .or. .not.initialized

    if ( read_data ) then

       find_ndx: do i = last_index, ntimes
          if ( data_times(i) > time ) then
             index = i-1
             exit find_ndx
          endif
       enddo find_ndx

       last_index = index+1

       nt = 2

       if ( index < 1 ) then
          write(iulog,102) year,month,day,sec
          call endrun('solar_data_advance: failed to read data from '//trim(solar_data_file))
       endif

       ! get the surrounding time slices
       offset = (/ 1, index /)
       count =  (/ nbins, nt /)

   
       if (has_spectrum) then
          ierr = pio_get_var( file_id, ssi_vid, offset, count, irradi )
       endif
       if (has_tsi .and. (.not.do_spctrl_scaling)) then
          ierr = pio_get_var( file_id, tsi_vid, (/index/), (/nt/), itsi )
          if ( any(itsi(:nt) < 0._r8) ) then
             call endrun( 'solar_data_advance: invalid or missing tsi data  ' )
          endif
       endif
    else
       index = last_index - 1
    endif

    delt = ( time - data_times(index) ) / ( data_times(index+1) - data_times(index) )

    ! this assures that FIXED data are b4b on restarts
    if ( fixed_solar ) then
       delt = dble(int(delt*cday+.5_r8))/dble(cday)
    endif

    if (has_spectrum) then
       data(:) = irradi(:,1) + delt*( irradi(:,2) - irradi(:,1) )
       do i = 1,nbins
          sol_irrad(i) = data(i)*irrad_fac(i) ! W/m2/nm
          sol_etf(i)   = data(i)*etf_fac(i)   ! photons/cm2/sec/nm 
       enddo
    endif
    if (has_tsi .and. (.not.do_spctrl_scaling)) then
       sol_tsi = itsi(1) + delt*( itsi(2) - itsi(1) )
       if ( masterproc ) then
          if (day == 1 .and. sec == 0) then
             write(iulog,101) year, month, day, sec, sol_tsi
          end if
       endif
    endif

 101 FORMAT('solar_data_advance: date, TSI : ',i4.4,'-',i2.2,'-',i2.2,'-',i5.5,',  ',f12.6)
 102 FORMAT('solar_data_advance: not able to find data for : ',i4.4,'-',i2.2,'-',i2.2,'-',i5.5)

  end subroutine solar_data_advance
  
  !---------------------------------------------------------------------------
  ! private methods
  !---------------------------------------------------------------------------
  subroutine convert_dates( dates, secs, times )

    use time_manager, only: set_time_float_from_date

    integer,  intent(in)  :: dates(:)
    integer,  intent(in)  :: secs(:)

    real(r8), intent(out) :: times(:)

    integer :: year, month, day, sec,n ,i

    n = size( dates ) 

    do i=1,n
       year = dates(i)/10000
       month = (dates(i)-year*10000)/100
       day = dates(i)-year*10000-month*100
       sec = secs(i)
       call set_time_float_from_date( times(i), year, month, day, sec )
    enddo

  end subroutine convert_dates

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  subroutine convert_date( date, sec, time )

    integer,  intent(in)  :: date
    integer,  intent(in)  :: sec
    real(r8), intent(out) :: time

    integer :: dates(1), secs(1)
    real(r8) :: times(1)
    dates(1) = date
    secs(1) = sec
    call convert_dates( dates, secs, times )
    time = times(1)
  end subroutine convert_date

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine get_model_time( time, year, month, day, seconds )

    use time_manager, only: get_curr_date

    real(r8), intent(out) :: time
    integer, optional, intent(out) :: year, month, day, seconds

    integer  :: yr, mn, dy, sc, date

    call get_curr_date(yr, mn, dy, sc)
    date = yr*10000 + mn*100 + dy
    call convert_date( date, sc, time )

    if (present(year))    year = yr
    if (present(month))   month = mn
    if (present(day))     day = dy
    if (present(seconds)) seconds = sc

  end subroutine get_model_time

end module solar_data
