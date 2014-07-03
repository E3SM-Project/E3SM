module mo_strato_sad
!---------------------------------------------------------------
! 	... prescribed strat aero surf area density module
!---------------------------------------------------------------

  use m_types,      only : time_ramp
  use shr_kind_mod, only : r8 => shr_kind_r8
  use ppgrid,       only : pver, pcols, pverp, begchunk, endchunk
  use time_manager, only : get_curr_date
  use time_utils,   only : flt_date, moz_findplb
  use abortutils,   only : endrun
  use ioFileMod,    only : getfil
  use spmd_utils,   only : masterproc
  use cam_logfile,  only : iulog

  implicit none

  private
  public  :: strato_sad_inti
  public  :: strato_sad_timestep_init
  public  :: strato_sad_set

  save

  integer, parameter    :: time_span = 1

  integer               :: ntimes
  integer               :: nlon
  integer               :: nlat
  integer               :: nlev
  integer               :: tim_ndx(2)
  integer,  allocatable :: dates(:)
  real(r8), allocatable :: times(:)
  real(r8), allocatable :: sad_lats(:)
  real(r8), allocatable :: sad_lons(:)
  real(r8), allocatable :: sad_levs(:)
  real(r8), allocatable :: sage_sad(:,:,:,:)
  character(len=256)     :: filename

  logical :: has_sulfate_rxts
  type(time_ramp)   :: strato_sad_timing

contains

  subroutine strato_sad_inti( sad_file, sad_timing )
!-----------------------------------------------------------------------
! 	... initialize the strato sad module
!-----------------------------------------------------------------------
    
    use mo_constants,  only : d2r
    use mo_chem_utls,  only : get_rxt_ndx
    use mo_strato_rates,only : has_strato_chem
    use string_utils,  only : to_upper
    use pio,           only : file_desc_t, pio_closefile, pio_inq_dimid, &
         pio_inq_dimlen, pio_inq_varid, pio_get_var, pio_nowrite
    use cam_pio_utils, only : cam_pio_openfile
    implicit none

!-----------------------------------------------------------------------
! 	... dummy arguments
!-----------------------------------------------------------------------
    character(len=*), intent(in) :: sad_file
    type(time_ramp),  intent(in) :: sad_timing

!-----------------------------------------------------------------------
! 	... local variables
!-----------------------------------------------------------------------
    real(r8), parameter :: hPa2Pa = 100._r8
    integer :: astat
    integer :: j, l, m, n                           ! Indices
    integer :: t1, t2
    type(file_desc_t) :: ncid
    integer :: dimid
    integer :: varid
    integer :: yr, wrk_date, wrk_sec
    real(r8)    :: seq
    real(r8)    :: wrk_time
    character(len=8)  :: time_type

    integer :: mon, day, tod, ncdate, ncsec
    integer :: ierr
    integer :: usrrxt_ndx(4)

    usrrxt_ndx(1) = get_rxt_ndx( 'usr_N2O5_aer' )
    usrrxt_ndx(2) = get_rxt_ndx( 'usr_NO3_aer' )
    usrrxt_ndx(3) = get_rxt_ndx( 'usr_NO2_aer' )
    usrrxt_ndx(4) = get_rxt_ndx( 'usr_HO2_aer' )

    has_sulfate_rxts = any( usrrxt_ndx > 0 )
    has_sulfate_rxts = has_sulfate_rxts .or. has_strato_chem

    if ( .not. has_sulfate_rxts) return

    call get_curr_date( yr, mon, day, tod )
    ncdate = yr*10000 + mon*100 + day
    ncsec = tod

!-----------------------------------------------------------------------
! 	... check timing
!-----------------------------------------------------------------------
    strato_sad_timing = sad_timing 
    strato_sad_timing%type = to_upper(strato_sad_timing%type)
    time_type = strato_sad_timing%type 


    if( time_type /= 'SERIAL' .and. time_type /= 'CYCLICAL' .and. time_type /= 'FIXED' ) then
       write(iulog,*) 'strato_sad_inti: time type ',trim(time_type),' is not SERIAL, CYCLICAL, or FIXED'
       call endrun
    end if

    wrk_sec  = ncsec
    if( (time_type == 'SERIAL') ) then
       wrk_date = ncdate
    else if( time_type == 'CYCLICAL' ) then
       wrk_date = strato_sad_timing%cycle_yr*10000 + mod( ncdate,10000 )
    else
       wrk_date = strato_sad_timing%fixed_ymd
       wrk_sec  = strato_sad_timing%fixed_tod
    end if
    wrk_time = flt_date( wrk_date, wrk_sec )
    if (masterproc) then
       write(iulog,*) 'strato_sad_inti: wrk_date,wrk_sec,wrk_time = ',wrk_date,wrk_sec,wrk_time

       !-----------------------------------------------------------------------
       ! 	... diagnostics
       !-----------------------------------------------------------------------
       write(iulog,*) ' '
       write(iulog,*) 'strato_sad_inti: diagnostics'
       write(iulog,*) ' '
       write(iulog,*) 'strato sad timing specs'
       write(iulog,*) 'type = ',strato_sad_timing%type
       if( time_type == 'CYCLICAL' ) then
          write(iulog,*) 'cycle year = ',strato_sad_timing%cycle_yr
       else
          write(iulog,*) 'fixed date = ',strato_sad_timing%fixed_ymd
          write(iulog,*) 'fixed time = ',strato_sad_timing%fixed_tod
       end if
    endif

!-----------------------------------------------------------------------
! 	... open netcdf file
!-----------------------------------------------------------------------
    call getfil ( sad_file, filename, 0 )
    call cam_pio_openfile (ncid, trim(filename), PIO_NOWRITE)
!-----------------------------------------------------------------------
! 	... get timing information, allocate arrays, and read in dates
!-----------------------------------------------------------------------
    ierr = pio_inq_dimid( ncid, 'time', dimid )
    ierr = pio_inq_dimlen( ncid, dimid, ntimes )
    allocate( dates(ntimes),stat=astat )
    if( astat/= 0 ) then
       write(iulog,*) 'strato_sad_inti: failed to allocate dates array; error = ',astat
       call endrun
    end if
    allocate( times(ntimes),stat=astat )
    if( astat/= 0 ) then
       write(iulog,*) 'strato_sad_inti: failed to allocate times array; error = ',astat
       call endrun
    end if
    ierr = pio_inq_varid( ncid, 'date', varid )
    ierr = pio_get_var( ncid, varid, dates )
    do n = 1,ntimes
       times(n) = flt_date( dates(n), 0 )
    end do
    if( time_type /= 'CYCLICAL' ) then
       if( wrk_time < times(1) .or. wrk_time > times(ntimes) ) then
          write(iulog,*) 'strato_sad_inti: time out of bounds for dataset = ',trim(filename)
          call endrun
       end if
       do n = 2,ntimes
          if( wrk_time <= times(n) ) then
             exit
          end if
       end do
       tim_ndx(1) = n - 1
    else
       yr = strato_sad_timing%cycle_yr
       do n = 1,ntimes
          if( yr == dates(n)/10000 ) then
             exit
          end if
       end do
       if( n >= ntimes ) then
          write(iulog,*) 'strato_sad_inti: time out of bounds for dataset = ',trim(filename)
          call endrun
       end if
       tim_ndx(1) = n
    end if
    select case( time_type )
    case( 'FIXED' )
       tim_ndx(2) = n
    case( 'CYCLICAL' )
       do n = tim_ndx(1),ntimes
          if( yr /= dates(n)/10000 ) then
             exit
          end if
       end do
       tim_ndx(2) = n - 1
       if( (tim_ndx(2) - tim_ndx(1)) < 2 ) then
          write(iulog,*) 'strato_sad_inti: cyclical sad require at least two time points'
          call endrun
       end if
    case( 'SERIAL' )
       tim_ndx(2) = min( ntimes,tim_ndx(1) + time_span )
    end select
    t1 = tim_ndx(1)
    t2 = tim_ndx(2)
    
    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) 'strato_sad time cnt = ',ntimes
       write(iulog,*) 'strato_sad times'
       write(iulog,'(10i10)') dates(:)
       write(iulog,'(1p,5g15.7)') times(:)
       write(iulog,*) 'strato_sad time indicies = ',tim_ndx(:)
       write(iulog,'(10i10)') dates(tim_ndx(1):tim_ndx(2))
       write(iulog,*) ' '
    endif

!-----------------------------------------------------------------------
!     	... inquire about latitudes
!-----------------------------------------------------------------------
    ierr = pio_inq_dimid( ncid, 'lat', dimid )
    ierr = pio_inq_dimlen( ncid, dimid, nlat )
!-----------------------------------------------------------------------
!     	... allocate space for latitude coordinates
!-----------------------------------------------------------------------
    allocate( sad_lats(nlat),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) ' strato_sad_inti : failed to allocate latitudes array; error = ',astat
       call endrun
    end if
!-----------------------------------------------------------------------
!     	... get latitudes
!-----------------------------------------------------------------------
    ierr = pio_inq_varid( ncid, 'lat', varid )
    ierr = pio_get_var( ncid, varid, sad_lats  )

!-----------------------------------------------------------------------
!     	... inquire about longitudes
!-----------------------------------------------------------------------
    ierr = pio_inq_dimid( ncid, 'lon', dimid )
    ierr = pio_inq_dimlen( ncid, dimid, nlon )
!-----------------------------------------------------------------------
!     	... allocate space for longitude coordinates
!-----------------------------------------------------------------------
    allocate( sad_lons(nlon),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) ' strato_sad_inti : failed to allocate longitudes array; error = ',astat
       call endrun
    end if
!-----------------------------------------------------------------------
!     	... get longitudes
!-----------------------------------------------------------------------
    ierr = pio_inq_varid( ncid, 'lon', varid )
    ierr = pio_get_var( ncid, varid, sad_lons )

!-----------------------------------------------------------------------
!     	... convert to radians and setup regridding
!-----------------------------------------------------------------------
    sad_lats(:) = d2r * sad_lats(:)
    sad_lons(:) = d2r * sad_lons(:)

!-----------------------------------------------------------------------
!     	... inquire about levels
!-----------------------------------------------------------------------
    ierr = pio_inq_dimid( ncid, 'lev', dimid )
    ierr = pio_inq_dimlen( ncid, dimid, nlev )
!-----------------------------------------------------------------------
!     	... allocate space for level coordinates
!-----------------------------------------------------------------------
    allocate( sad_levs(nlev),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) ' strato_sad_inti : failed to allocate levels array; error = ',astat
       call endrun
    end if
!-----------------------------------------------------------------------
!     	... get levels
!-----------------------------------------------------------------------
    ierr = pio_inq_varid( ncid, 'lev', varid )
    ierr = pio_get_var( ncid, varid, sad_levs )
    sad_levs(:) = hPa2Pa * sad_levs(:)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) 'strato_sad_inti: sad_levs'
       write(iulog,'(10g12.5)') sad_levs(:)
       write(iulog,*) ' '
    endif

!-----------------------------------------------------------------------
!     	... allocate module sad variable
!-----------------------------------------------------------------------
    allocate( sage_sad(pcols,nlev,begchunk:endchunk,t1:t2),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) ' strato_sad_inti : failed to allocate sage_sad array; error = ',astat
       call endrun
    end if
!-----------------------------------------------------------------------
! 	... read and interpolate sage sad data
!-----------------------------------------------------------------------
    call strato_sad_get( ncid )
!-----------------------------------------------------------------------
! 	... close netcdf file
!-----------------------------------------------------------------------
    call pio_closefile( ncid )

  end subroutine strato_sad_inti

  subroutine strato_sad_timestep_init( )
!-----------------------------------------------------------------------
!       ... check serial case for time span
!-----------------------------------------------------------------------
    use cam_pio_utils, only : cam_pio_openfile
    use pio, only : pio_closefile, file_desc_t, pio_nowrite
    implicit none

!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
    integer                     :: m
    integer                     :: t1, t2, tcnt
    integer                     :: astat
    type(file_desc_t)                     :: ncid
    real(r8)                    :: wrk_time

    integer ::  yr, mon, day, tod, ncdate, ncsec

    if ( .not. has_sulfate_rxts) return

    call get_curr_date( yr, mon, day, tod )
    ncdate = yr*10000 + mon*100 + day
    ncsec = tod

	
    if( strato_sad_timing%type == 'SERIAL' ) then
       wrk_time = flt_date( ncdate, ncsec )
       if( wrk_time > times(tim_ndx(2)) ) then
          tcnt = tim_ndx(2) - tim_ndx(1)
          tim_ndx(1) = tim_ndx(2)
          tim_ndx(2) = min( ntimes,tim_ndx(1) + time_span )
          t1 = tim_ndx(1)
          t2 = tim_ndx(2)
!!$          if( tcnt /= (t2 - t1) ) then
!-----------------------------------------------------------------------
! 	... allocate array
!-----------------------------------------------------------------------
             if( allocated( sage_sad ) ) then
                deallocate( sage_sad,stat=astat )
                if( astat/= 0 ) then
                   write(iulog,*) 'strato_sad_timestep_init: failed to deallocate strato_sad vmr; error = ',astat
                   call endrun
                end if
             end if
             allocate( sage_sad(pcols,nlev,begchunk:endchunk,t1:t2),stat=astat )
             if( astat/= 0 ) then
                write(iulog,*) 'strato_sad_timestep_init: failed to allocate sage_sad; error = ',astat
                call endrun
             end if
!!$          end if
!-----------------------------------------------------------------------
! 	... open netcdf file
!-----------------------------------------------------------------------
          call cam_pio_openfile (ncid, trim(filename), PIO_NOWRITE)
!-----------------------------------------------------------------------
! 	... read and interpolate sage sad data
!-----------------------------------------------------------------------
          call strato_sad_get( ncid )
!-----------------------------------------------------------------------
! 	... close netcdf file
!-----------------------------------------------------------------------
          call pio_closefile( ncid )
       end if
    end if

  end subroutine strato_sad_timestep_init

  subroutine strato_sad_get( ncid )
!-----------------------------------------------------------------------
!       ... read sad values
!-----------------------------------------------------------------------
    
    use interpolate_data, only : lininterp_init, lininterp, lininterp_finish, &
         interp_type
    use pio, only : file_desc_t, pio_inq_varid, pio_inq_varndims, pio_get_var
    use phys_grid, only : get_ncols_p, get_rlat_all_p, get_rlon_all_p
    use mo_constants, only : pi
    implicit none

!-----------------------------------------------------------------------
!       ... dummy arguments
!-----------------------------------------------------------------------
    type(file_desc_t), intent(in)           :: ncid

!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
    integer               :: i, c, j, k, n                            ! Indices
    integer               :: t1, t2, tcnt
    integer               :: ierr
    integer               :: vid
    real(r8), allocatable :: sage_sad_in(:,:,:,:)
    real(r8), allocatable :: tmp_3d(:,:,:)
    integer :: ndims, ncols
    real(r8) :: to_lats(pcols), to_lons(pcols)
    type(interp_type) :: lon_wgts, lat_wgts

    real(r8), parameter :: zero=0._r8, twopi=2._r8*pi
    
!-----------------------------------------------------------------------
!       ... read sage data
!-----------------------------------------------------------------------
    t1 = tim_ndx(1)
    t2 = tim_ndx(2)
    tcnt = t2 - t1 + 1

!-----------------------------------------------------------------------
!       ... allocate local wrk arrays
!-----------------------------------------------------------------------
    allocate( sage_sad_in(nlon,nlat,nlev,t1:t2), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'strato_sad_get: wrk allocation error = ',ierr
       call endrun
    end if

    ierr = pio_inq_varid(  ncid, 'sad_sage', vid )
    ierr = pio_inq_varndims (ncid, vid, ndims)

    if ( ndims == 4 ) then
       ierr = pio_get_var( ncid, vid, &
            (/ 1, 1, 1, t1/), &                     ! start
            (/ nlon, nlat, nlev, tcnt /), &  ! count
            sage_sad_in   )
    else
       allocate( tmp_3d(nlat,nlev,t1:t2), stat=ierr )
       if( ierr /= 0 ) then
          write(iulog,*) 'strato_sad_get: tmp_3d allocation error = ',ierr
          call endrun
       end if

       ierr = pio_get_var( ncid, vid, &
            (/ 1, 1, t1/), &                     ! start
            (/ nlat, nlev, tcnt /), &  ! count
            tmp_3d   )

       do i = 1,nlon
          sage_sad_in( i, :,:,:) = tmp_3d(:,:,:)
       enddo

       deallocate( tmp_3d, stat=ierr )
       if( ierr /= 0 ) then
          write(iulog,*) 'strato_sad_get: failed to deallocate tmp_3d, ierr = ',ierr
          call endrun
       end if

    endif


    do c=begchunk,endchunk
       ncols = get_ncols_p(c)
       call get_rlat_all_p(c, pcols, to_lats)
       call get_rlon_all_p(c, pcols, to_lons)
       call lininterp_init(sad_lons, nlon, to_lons, ncols, 2, lon_wgts, zero, twopi)
       call lininterp_init(sad_lats, nlat, to_lats, ncols, 1, lat_wgts)


       do n = t1,t2
          do k = 1,nlev
             call lininterp(sage_sad_in(:,:,k,n), nlon, nlat, sage_sad(:,k,c,n),&
                  ncols,lon_wgts,lat_wgts)
          end do
       end do
       call lininterp_finish(lat_wgts)
       call lininterp_finish(lon_wgts)
    end do
    deallocate( sage_sad_in, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'strato_sad_get: failed to deallocate sage_sad_in, ierr = ',ierr
       call endrun
    end if

  end subroutine strato_sad_get

  subroutine strato_sad_set( pmid, sad, ncol, lchnk )
!--------------------------------------------------------
!	... set the sad values
!--------------------------------------------------------
    
    implicit none

!--------------------------------------------------------
!	... dummy arguments
!--------------------------------------------------------
    integer, intent(in)     :: ncol, lchnk
    real(r8), intent(in)    :: pmid(pcols,pver) ! midpoint pressure (Pa)
    real(r8), intent(inout) :: sad(pcols, pver) ! stratospheric aerosol surface area density (1/cm)

!--------------------------------------------------------
!	... local variables
!--------------------------------------------------------
    integer  ::  i, k, n
    integer  ::  last, next
    integer  ::  wrk_date, wrk_sec
    integer  ::  tcnt
    integer  ::  astat
    real(r8) ::  dels
    real(r8) ::  wrk_time
    real(r8) ::  sad_int(ncol,pver,2)
    integer  ::  yr, mon, day, tod, ncdate, ncsec

    if ( .not. has_sulfate_rxts) return

    call get_curr_date( yr, mon, day, tod )
    ncdate = yr*10000 + mon*100 + day
    ncsec = tod

!--------------------------------------------------------
!	... setup the time interpolation
!--------------------------------------------------------
    wrk_sec  = ncsec
    select case( strato_sad_timing%type )
    case( 'SERIAL' )
       wrk_date = ncdate
    case( 'CYCLICAL' )
       wrk_date = strato_sad_timing%cycle_yr*10000 + mod( ncdate,10000 )
    case( 'FIXED' )
       wrk_date = strato_sad_timing%fixed_ymd
       wrk_sec  = strato_sad_timing%fixed_tod
    end select
    wrk_time = flt_date( wrk_date, wrk_sec )

!--------------------------------------------------------
!	... set time interpolation factor
!--------------------------------------------------------
    if( strato_sad_timing%type /= 'CYCLICAL' ) then
       do n = tim_ndx(1)+1,tim_ndx(2)
          if( wrk_time <= times(n) ) then
             last = n - 1
             next = n
             exit
          end if
       end do
       if( n > ntimes ) then
          write(iulog,*) 'strato_sad_set: interp time is out of bounds'
          call endrun
       end if
       dels = (wrk_time - times(last))/(times(next) - times(last))
    else
       tcnt = tim_ndx(2) - tim_ndx(1) + 1
       call moz_findplb( times(tim_ndx(1)), tcnt, wrk_time, n )
       if( n < tcnt ) then
          last = tim_ndx(1) + n - 1
          next = last + 1
          dels = (wrk_time - times(last))/(times(next) - times(last))
       else
          next = tim_ndx(1)
          last = tim_ndx(2)
          dels = wrk_time - times(last)
          if( dels < 0._r8 ) then
             dels = 365._r8 + dels
          end if
          dels = dels/(365._r8 + times(next) - times(last))
       end if
    end if
    
    dels = max( min( 1._r8,dels ),0._r8 )

#ifdef DEBUG
    write(iulog,*) ' '
    write(iulog,*) 'strato_sad_set: last,next,dels,ncdate,ncsec = ',last,next,dels,ncdate,ncsec
    write(iulog,*) 'strato_sad_set: dates(last),dates(next)     = ',dates(last),dates(next)
#endif

    call vinterp( pmid, sage_sad(:,:,lchnk,last), sad_int(:,:,1), ncol )
    call vinterp( pmid, sage_sad(:,:,lchnk,next), sad_int(:,:,2), ncol )

#ifdef DEBUG
    write(iulog,*) 'strato_sad_set: pmid(1,:)'
    write(iulog,'(1p,5g15.7)') pmid(1,:)
    write(iulog,*) 'strato_sad_set: sad_levs'
    write(iulog,'(1p,5g15.7)') sad_levs(:)
    write(iulog,*) 'strato_sad_set: sage_sad(last)'
    write(iulog,'(1p,5g15.7)') sage_sad(1,:,lat,ip,last)
    write(iulog,*) 'strato_sad_set: sage_sad(next)'
    write(iulog,'(1p,5g15.7)') sage_sad(1,:,lat,ip,next)
#endif

    do i=1,ncol
       sad(i,:) = sad_int(i,:,1) + dels * (sad_int(i,:,2) - sad_int(i,:,1))
    enddo

#ifdef DEBUG
    write(iulog,*) 'strato_sad_set: sad'
    write(iulog,'(1p,5g15.7)') sad(1,:)
    write(iulog,*) ' '
    !call endrun 
#endif

  end subroutine strato_sad_set

  subroutine vinterp( pmid, sad_src, sad_int, ncol )
!-----------------------------------------------------------------------
!   	... vertically interpolate input data
!-----------------------------------------------------------------------
    implicit none

!-----------------------------------------------------------------------
!   	... dummy arguments
!-----------------------------------------------------------------------
    integer,  intent(in)  :: ncol
    real(r8), intent(in)  :: pmid(:,:)
    real(r8), intent(in)  :: sad_src(:,:)
    real(r8), intent(out) :: sad_int(:,:)

!-----------------------------------------------------------------------
!   	... local variables
!-----------------------------------------------------------------------
    integer :: i
    integer :: k, kl, ku
    real(r8)    :: delp, pinterp

    level_loop : do k = 1,pver
       long_loop : do i = 1,ncol
          pinterp = pmid(i,k)
          if( pinterp <= sad_levs(1) ) then
             sad_int(i,k) = sad_src(i,1)
          else if( pinterp > sad_levs(nlev) ) then
             sad_int(i,k) = sad_src(i,nlev)
          else
             do ku = 2,nlev
                if( pinterp <= sad_levs(ku) ) then
                   kl = ku - 1
                   delp = log( pinterp/sad_levs(kl) ) &
                        / log( sad_levs(ku)/sad_levs(kl) )
                   sad_int(i,k) = sad_src(i,kl) + delp * (sad_src(i,ku) - sad_src(i,kl))
                   exit
                end if
             end do
          end if
       end do long_loop
    end do level_loop

  end subroutine vinterp


  





end module mo_strato_sad
