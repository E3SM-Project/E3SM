module mo_flbc
  !---------------------------------------------------------------
  ! 	... lower boundary module
  !---------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8
  use m_types,      only : time_ramp
  use spmd_utils,   only : masterproc,iam
  use cam_abortutils,   only : endrun
  use ioFileMod,    only : getfil
  use ppgrid,       only : pcols, begchunk, endchunk, pver
  use time_manager, only : get_curr_date, get_curr_calday
  use time_utils,   only : flt_date
  use cam_logfile,  only : iulog
  use constituents,  only : pcnst
  use constituents,  only : tracnam=>cnst_name

  implicit none

  type :: flbc
     integer            :: spc_ndx = -1
     real(r8), pointer  :: vmr(:,:,:)
     character(len=16)  :: species = ' '
     logical            :: has_mean
     real(r8), pointer  :: vmr_mean(:)
  end type flbc

  private
  public  :: flbc_inti, flbc_set, flbc_chk, has_flbc
  public  :: flbc_gmean_vmr

  save

  integer, parameter :: time_span = 1

  integer :: ntimes
  integer :: flbc_cnt
  integer :: gndx
  integer :: tim_ndx(2)
  integer :: jlim(2)
  integer, allocatable  :: dates(:)
  real(r8), allocatable     :: times(:)
  logical :: has_flbc(pcnst)
  character(len=256) :: filename, lpath, mspath

  type(time_ramp) :: flbc_timing
  integer ::  ncdate, ncsec

  integer, parameter :: nghg = 5
  integer, parameter :: max_nflbc = pcnst+nghg

  integer, parameter :: co2_ndx = 1
  integer, parameter :: ch4_ndx = 2
  integer, parameter :: n2o_ndx = 3
  integer, parameter :: f11_ndx = 4
  integer, parameter :: f12_ndx = 5
  character(len=5)  :: ghg_names(nghg) = (/ 'CO2  ','CH4  ','N2O  ','CFC11','CFC12' /)
  integer :: ghg_indices(nghg) = -1

  type(flbc) :: flbcs(max_nflbc)

  logical, parameter :: debug = .false.

contains

  subroutine flbc_inti( flbc_file, flbc_list, flbc_timing_in, co2vmr, ch4vmr, n2ovmr, f11vmr, f12vmr )
    !-----------------------------------------------------------------------
    ! 	... initialize the fixed lower bndy cond
    !-----------------------------------------------------------------------

    use mo_constants,  only : d2r, pi, rearth
    use string_utils,  only : to_upper
    use constituents,  only : cnst_get_ind
    use cam_pio_utils, only : cam_pio_openfile
    use pio,           only : pio_get_var,pio_inq_varid,pio_inq_dimid, pio_inq_dimlen
    use pio,           only : file_desc_t, pio_closefile, pio_nowrite

    implicit none

    !-----------------------------------------------------------------------
    ! 	... dummy arguments
    !-----------------------------------------------------------------------
    character(len=*), intent(in) :: flbc_file
    character(len=*), intent(in) :: flbc_list(:)
    type(time_ramp),  intent(in) :: flbc_timing_in
    real(r8),         intent(in) :: co2vmr, ch4vmr, n2ovmr, f11vmr, f12vmr

    !-----------------------------------------------------------------------
    ! 	... local variables
    !-----------------------------------------------------------------------
    integer :: astat
    integer :: j, l, m, n                     ! Indices
    integer :: t1, t2
    type(file_desc_t) :: ncid
    integer :: dimid
    integer :: varid
    integer :: yr, mon, day, wrk_date, wrk_sec
    real(r8)    :: seq
    real(r8)    :: wrk_time
    character(len=16)  :: species
    character(len=16)  :: spc_name
    character(len=8)   :: time_type
    integer :: ierr

    if ( len_trim( flbc_file ) == 0 ) return

    call get_curr_date( yr, mon, day, ncsec )
    ncdate = yr*10000 + mon*100 + day

    !-----------------------------------------------------------------------
    ! 	... check timing
    !-----------------------------------------------------------------------
    flbc_timing = flbc_timing_in
    time_type = to_upper(flbc_timing%type)
    flbc_timing%type = time_type
    if( time_type /= 'SERIAL' .and. time_type /= 'CYCLICAL' &
         .and. time_type /= 'FIXED' ) then
       write(iulog,*) 'flbc_inti: time type ',trim(time_type),' is not SERIAL,CYCLICAL, or FIXED'
       call endrun('flbc_inti: invalid time type ')
    end if

    if ( (flbc_timing%cycle_yr>0) .and. (time_type/='CYCLICAL') ) then
       call endrun('flbc_inti: cannot specify  flbc_cycle_yr if flbc_type is not CYCLICAL')
    endif
    if ( ((flbc_timing%fixed_ymd>0).or.(flbc_timing%fixed_tod>0)).and.(time_type/='FIXED') ) then
       call endrun('flbc_inti: cannot specify  flbc_fixed_ymd or flbc_fixed_tod if flbc_type is not FIXED')
    endif

    wrk_sec  = ncsec
    if( time_type == 'SERIAL' ) then
       wrk_date = ncdate 
    else if( time_type == 'CYCLICAL' ) then

    	! If this is a leap-day, we have to avoid asking for a non-leap-year
    	! on a cyclical dataset. When this happens, just use Feb 28 instead
    	if (( mon .eq. 2 ) .and. ( day.eq.29 )) then
	   ncdate = yr*10000 + mon*100 + (day-1)
           write(iulog,*)'WARNING: flbc_inti using Feb 28 instead of Feb 29 for cyclical dataset'
        endif 	
       wrk_date = flbc_timing%cycle_yr*10000 + mod(ncdate,10000)
    else
       wrk_date = flbc_timing%fixed_ymd
       wrk_sec  = flbc_timing%fixed_tod
    end if
    wrk_time = flt_date( wrk_date, wrk_sec )
    if (masterproc) write(iulog,*) 'flbc_inti: wrk_date,wrk_sec,wrk_time = ',wrk_date,wrk_sec,wrk_time

    !-----------------------------------------------------------------------
    ! 	... species with fixed lbc ?
    !-----------------------------------------------------------------------
    has_flbc(:) = .false.
    flbc_cnt = 0
    
    do m = 1,max_nflbc

       if ( len_trim(flbc_list(m))==0 ) exit

       flbc_cnt = flbc_cnt + 1

       call cnst_get_ind (flbc_list(m), n, abort=.false.)

       if (n > 0) then
          has_flbc(n) = .true.
          flbcs(flbc_cnt)%spc_ndx = n
       else ! must be one of the GHGs which is not prognosted
          if( .not. any( ghg_names(:) == flbc_list(m) ) ) then
             call endrun('flbc_inti: flbc_list member '// trim(flbc_list(m)) //' is not allowed')
          endif
          flbcs(flbc_cnt)%spc_ndx = -1
       endif

       flbcs(flbc_cnt)%species = trim( flbc_list(m) )

       where( ghg_names(:) == flbc_list(m) )
          ghg_indices = m
       endwhere

       if( trim(flbcs(flbc_cnt)%species) == 'CFC11' ) then
          flbcs(flbc_cnt)%species = 'CFCL3'
       elseif( trim(flbcs(flbc_cnt)%species) == 'CFC12' ) then
          flbcs(flbc_cnt)%species = 'CF2CL2'
       endif

    enddo

    ! check that user has not set vmr namelist values... 
    if ( ghg_indices(co2_ndx) > 0 .and. co2vmr>1.e-6_r8) then
       call endrun('flbc_inti: cannot specify both co2vmr and CO2 in flbc_file')
    endif
    if ( ghg_indices(ch4_ndx) > 0 .and. ch4vmr > 0._r8) then
       call endrun('flbc_inti: cannot specify both ch4vmr and CH4 in flbc_file')
    endif
    if ( ghg_indices(n2o_ndx) > 0 .and. n2ovmr > 0._r8) then
       call endrun('flbc_inti: cannot specify both n2ovmr and N2O in flbc_file')
    endif
    if ( ghg_indices(f11_ndx) > 0 .and. f11vmr > 0._r8) then
       call endrun('flbc_inti: cannot specify both f11vmr and CFC11 in flbc_file')
    endif
    if ( ghg_indices(f12_ndx) > 0 .and. f12vmr > 0._r8) then
       call endrun('flbc_inti: cannot specify both f12vmr and CFC12 in flbc_file')
    endif
    
    if( flbc_cnt == 0 ) then
       return
    end if

    if(masterproc) then
       write(iulog,*) ' '
       if( flbc_cnt > 0 ) then
          write(iulog,*) 'flbc_inti: Species with specified lower boundary values'
          do n = 1,flbc_cnt
             write(iulog,*) trim(flbcs(n)%species)
          enddo
       else
          write(iulog,*) 'There are no species with specified lower boundary values'
       end if
       write(iulog,*) ' '

       !-----------------------------------------------------------------------
       ! 	... diagnostics
       !-----------------------------------------------------------------------
       write(iulog,*) ' '
       write(iulog,*) 'flbc_inti: diagnostics'
       write(iulog,*) ' '
       write(iulog,*) 'lower bndy timing specs'
       write(iulog,*) 'type = ',flbc_timing%type
       if( time_type == 'CYCLICAL' ) then
          write(iulog,*) 'cycle year = ',flbc_timing%cycle_yr
       else
          write(iulog,*) 'fixed date = ',flbc_timing%fixed_ymd
          write(iulog,*) 'fixed time = ',flbc_timing%fixed_tod
       end if
       write(iulog,*) ' '
       write(iulog,*) 'there are ',flbc_cnt,' species with specified lower bndy values'
       write(iulog,*) ' '
    end if
    !-----------------------------------------------------------------------
    ! 	... get timing information, allocate arrays, and read in dates
    !-----------------------------------------------------------------------
    call getfil ( flbc_file, filename, 0)
    call cam_pio_openfile (ncid, trim(filename), PIO_NOWRITE)
    ierr = pio_inq_dimid( ncid, 'time', dimid )
    ierr = pio_inq_dimlen( ncid, dimid, ntimes )

    allocate( dates(ntimes),stat=astat )
    if( astat/= 0 ) then
       write(iulog,*) 'flbc_inti: failed to allocate dates array; error = ',astat
       call endrun
    end if
    allocate( times(ntimes),stat=astat )
    if( astat/= 0 ) then
       write(iulog,*) 'flbc_inti: failed to allocate times array; error = ',astat
       call endrun
    end if

    ierr = pio_inq_varid( ncid, 'date', varid )
    ierr = pio_get_var( ncid, varid, dates )

    do n = 1,ntimes
       times(n) = flt_date( dates(n), 0 )
    end do
    if( time_type /= 'CYCLICAL' ) then
       if( wrk_time < times(1) .or. wrk_time > times(ntimes) ) then
          write(iulog,*) 'flbc_inti: time out of bounds for dataset = ',trim(filename)
          call endrun
       end if
       do n = 2,ntimes
          if( wrk_time <= times(n) ) then
             exit
          end if
       end do
       tim_ndx(1) = n - 1
    else
       yr = flbc_timing%cycle_yr
       do n = 1,ntimes
          if( yr == dates(n)/10000 ) then
             exit
          end if
       end do
       if( n >= ntimes ) then
          write(iulog,*) 'flbc_inti: time out of bounds for dataset = ',trim(filename)
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
          write(iulog,*) 'flbc_inti: cyclical lb conds require at least two time points'
          call endrun
       end if
    case( 'SERIAL' )
       tim_ndx(2) = min( ntimes,tim_ndx(1) + time_span )
    end select
    t1 = tim_ndx(1)
    t2 = tim_ndx(2)

    if( masterproc .and. debug ) then
       write(iulog,*) ' '
       write(iulog,*) 'flbc time cnt = ',ntimes
       write(iulog,*) 'flbc times'
       write(iulog,'(10i10)') dates(:)
       write(iulog,'(1p,5g15.7)') times(:)
       write(iulog,*) 'flbc time indicies = ',tim_ndx(:)
       write(iulog,'(10i10)') dates(tim_ndx(1):tim_ndx(2))
       write(iulog,*) ' '
    endif

    do m = 1,flbc_cnt
       !-----------------------------------------------------------------------
       ! 	... allocate array
       !-----------------------------------------------------------------------
       allocate( flbcs(m)%vmr(pcols,begchunk:endchunk,t1:t2),stat=astat )
       if( astat/= 0 ) then
          write(iulog,*) 'flbc_inti: failed to allocate lbc vmr; error = ',astat
          call endrun
       end if
       flbcs(m)%has_mean = file_has_gmean(ncid,flbcs(m)%species)
       if ( flbcs(m)%has_mean) then
          allocate( flbcs(m)%vmr_mean(t1:t2),stat=astat )
          if( astat/= 0 ) then
             write(iulog,*) 'flbc_inti: failed to allocate lbc vmr_mean; error = ',astat
             call endrun
          end if
       endif
       !-----------------------------------------------------------------------
       ! 	... readin the flbc vmr
       !-----------------------------------------------------------------------
       call flbc_get( ncid, flbcs(m), .true., read_gmean=flbcs(m)%has_mean )
    end do

    !-----------------------------------------------------------------------
    ! 	... close the file
    !-----------------------------------------------------------------------
    call pio_closefile( ncid )

  end subroutine flbc_inti

  subroutine flbc_chk( )
    use cam_pio_utils, only : cam_pio_openfile
    use pio,           only : file_desc_t, pio_closefile, pio_nowrite
    !-----------------------------------------------------------------------
    !       ... check serial case for time span
    !-----------------------------------------------------------------------

    implicit none

    !-----------------------------------------------------------------------
    !       ... dummy arguments
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !       ... local variables
    !-----------------------------------------------------------------------
    integer                     :: m
    integer                     :: t1, t2, tcnt
    integer                     :: astat
    type(file_desc_t)           :: ncid
    real(r8)                        :: wrk_time
    integer ::  yr, mon, day

    call get_curr_date( yr, mon, day, ncsec )
    ncdate = yr*10000 + mon*100 + day

    if( flbc_cnt > 0 .and. flbc_timing%type == 'SERIAL' ) then
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
          do m = 1,flbc_cnt
             if( associated( flbcs(m)%vmr ) ) then
                deallocate( flbcs(m)%vmr,stat=astat )
                if( astat/= 0 ) then
                   write(iulog,*) 'flbc_chk: failed to deallocate flbc vmr; error = ',astat
                   call endrun
                end if
             end if
             allocate( flbcs(m)%vmr(pcols,begchunk:endchunk,t1:t2),stat=astat )
             if( astat/= 0 ) then
                write(iulog,*) 'flbc_chk: failed to allocate flbc vmr; error = ',astat
                call endrun
             end if
                
             if (flbcs(m)%has_mean) then
                if( associated( flbcs(m)%vmr_mean ) ) then
                   deallocate( flbcs(m)%vmr_mean,stat=astat )
                   if( astat/= 0 ) then
                      write(iulog,*) 'flbc_chk: failed to deallocate flbc vmr; error = ',astat
                      call endrun
                   end if
                end if
                allocate( flbcs(m)%vmr_mean(t1:t2),stat=astat )
                if( astat/= 0 ) then
                   write(iulog,*) 'flbc_chk: failed to allocate flbc vmr; error = ',astat
                   call endrun
                end if

             endif
          end do
!!$          end if

          call cam_pio_openfile (ncid, trim(filename), PIO_NOWRITE)
          !-----------------------------------------------------------------------
          ! 	... readin the lb concentrations
          !-----------------------------------------------------------------------
          do m = 1,flbc_cnt
             call flbc_get( ncid, flbcs(m), .true., read_gmean=flbcs(m)%has_mean )
          end do

          !-----------------------------------------------------------------------
          ! 	... close the file
          !-----------------------------------------------------------------------
          call pio_closefile( ncid )

       end if
    end if

  end subroutine flbc_chk
  
  ! checks for global mean in input file
  function file_has_gmean(ncid,species)
    use pio, only : file_desc_t, pio_inq_varid, pio_noerr, pio_seterrorhandling, &
         pio_bcast_error, pio_internal_error
    implicit none

    type(file_desc_t),      intent(inout) :: ncid
    character(*), intent(in) :: species
    logical :: file_has_gmean

    integer :: varid, ierr

    ! Allow pio to return the potential error and handle it locally
    call pio_seterrorhandling(ncid, PIO_BCAST_ERROR)
    ierr = pio_inq_varid( ncid, trim(species)//'_LBC_mean', varid)
    call pio_seterrorhandling(ncid, PIO_INTERNAL_ERROR)


    file_has_gmean = (ierr==PIO_NOERR)

  endfunction file_has_gmean

  subroutine flbc_get( ncid, lbcs, initial, read_gmean )
    !-----------------------------------------------------------------------
    !       ... read lower bndy values
    !-----------------------------------------------------------------------
    use mo_constants,  only : d2r, pi
    use phys_grid,     only: get_ncols_p, get_rlat_all_p, get_rlon_all_p
    use pio,           only: file_desc_t, pio_get_var, pio_inq_varndims, &
         pio_max_name, pio_inq_varid, pio_inq_dimlen, pio_inq_dimid
    use interpolate_data, only : interp_type, lininterp_init, lininterp_finish, lininterp

    implicit none

    !-----------------------------------------------------------------------
    !       ... dummy arguments
    !-----------------------------------------------------------------------
    type(file_desc_t), intent(in)           :: ncid
    logical, intent(in)           :: initial
    type(flbc), intent(inout) :: lbcs

    logical, intent(in), optional :: read_gmean

    !-----------------------------------------------------------------------
    !       ... local variables
    !-----------------------------------------------------------------------
    integer                     :: j, m               ! Indices
    integer                     :: t1, t2, tcnt
    integer                     :: ierr
    integer                     :: vid, nlat, nlon
    integer                     :: dimid_lat, dimid_lon
    integer                     :: plon, plat
    real(r8), allocatable           :: lat(:)
    real(r8), allocatable           :: lon(:)
    real(r8), allocatable           :: wrk(:,:,:), wrk_zonal(:,:)
    real(r8), allocatable           :: wrk2d(:,:)
    character(len=pio_max_name)  :: varname
    real(r8), allocatable       :: locl_vmr(:,:,:)
    integer :: ndims, t, c, ncols
    type(interp_type) :: lon_wgts, lat_wgts
    real(r8) :: to_lats(pcols), to_lons(pcols)
    real(r8), parameter :: twopi=2._r8*pi, zero=0._r8

    t1 = tim_ndx(1)
    t2 = tim_ndx(2)
    tcnt = t2 - t1 + 1
    allocate( locl_vmr(pcols,begchunk:endchunk,tcnt), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'srf_emis_get: locl_emis allocation error = ',ierr
       call endrun
    end if

    locl_vmr(:,:,:) = 0._r8

    initialization : if( initial ) then
       !-----------------------------------------------------------------------
       !       ... get grid dimensions from file
       !-----------------------------------------------------------------------
       !           latitudes
       !-----------------------------------------------------------------------
       ierr = pio_inq_dimid( ncid, 'lat', dimid_lat )
       ierr = pio_inq_dimlen( ncid, dimid_lat, nlat )
       allocate( lat(nlat),stat=ierr )
       if( ierr /= 0 ) then
          write(iulog,*) 'flbc_get: lat allocation error = ',ierr
          call endrun
       end if
       ierr = pio_inq_varid( ncid, 'lat', vid )
       ierr = pio_get_var( ncid, vid, lat )
       lat(:nlat) = lat(:nlat) * d2r
       
       !-----------------------------------------------------------------------
       !           longitudes
       !-----------------------------------------------------------------------
       ierr = pio_inq_dimid( ncid, 'lon', dimid_lon )
       ierr = pio_inq_dimlen( ncid, dimid_lon, nlon )
       allocate( lon(nlon),stat=ierr )
       if( ierr /= 0 ) then
          write(iulog,*) 'flbc_get: lon allocation error = ',ierr
          call endrun
       end if
       ierr = pio_inq_varid( ncid, 'lon', vid )
       ierr = pio_get_var( ncid, vid, lon )
       lon(:nlon) = lon(:nlon) * d2r
    end if initialization
        
    allocate( wrk(nlon,nlat,tcnt), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'flbc_get: wrk allocation error = ',ierr
       call endrun
    end if

    !-----------------------------------------------------------------------
    !       ... read data
    !-----------------------------------------------------------------------
    varname = trim(lbcs%species) // '_LBC'
    ierr = pio_inq_varid( ncid, trim(varname), vid )
    ierr = pio_inq_varndims (ncid, vid, ndims)
    
    if (ndims==2) then
       allocate( wrk_zonal(nlat,tcnt), stat=ierr )
       if( ierr /= 0 ) then
          write(iulog,*) 'flbc_get: wrk_zonal allocation error = ',ierr
          call endrun
       end if
    endif

    if (ndims==2) then
       ierr = pio_get_var( ncid, vid, (/ 1, t1/), &
            (/ nlat, tcnt /), wrk_zonal )
       do t = 1,tcnt
          do j = 1,nlat
             wrk(:nlon,j,t) = wrk_zonal(j,t)
          enddo
       enddo
    else
       ierr = pio_get_var( ncid, vid, (/ 1, 1, t1/), &
            (/ nlon, nlat, tcnt /), wrk )
    endif

    do c=begchunk,endchunk
       ncols = get_ncols_p(c)
       call get_rlat_all_p(c, pcols, to_lats)
       call get_rlon_all_p(c, pcols, to_lons)
       call lininterp_init(lon, nlon, to_lons, ncols, 2, lon_wgts, zero, twopi)
       call lininterp_init(lat, nlat, to_lats, ncols, 1, lat_wgts)
          
       do m = 1,tcnt
          call lininterp(wrk(:,:,m), nlon, nlat, locl_vmr(:,c,m), ncols, lon_wgts, lat_wgts) 
       end do
          

       call lininterp_finish(lon_wgts)
       call lininterp_finish(lat_wgts)


    end do

    deallocate(wrk, stat=ierr)
    if( ierr /= 0 ) then
       write(iulog,*) 'flbc_get: Failed to deallocate wrk, ierr = ',ierr
       call endrun
    end if

    if (ndims==2) then
       deallocate( wrk_zonal,stat=ierr )
       if( ierr /= 0 ) then
          write(iulog,*) 'flbc_get: Failed to deallocate wrk_zonal, ierr = ',ierr
          call endrun
       end if
    end if
    if (read_gmean) then
       varname = trim(lbcs%species) // '_LBC_mean'
       ierr = pio_inq_varid( ncid, trim(varname), vid )
       ierr = pio_get_var( ncid, vid, (/t1/), (/tcnt/), lbcs%vmr_mean(t1:t2) )
    endif


    do m = t1,t2
       lbcs%vmr(:,:,m) = locl_vmr(:,:,m-t1+1)
    enddo

    deallocate(locl_vmr, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'flbc_get: Failed to deallocate locl_vmr; ierr = ',ierr
       call endrun
    end if

  end subroutine flbc_get

  subroutine flbc_set( vmr, ncol, lchnk, map )
    !--------------------------------------------------------
    !	... set the lower bndy values
    !--------------------------------------------------------

    implicit none

    !--------------------------------------------------------
    !	... dummy arguments
    !--------------------------------------------------------
    integer,  intent(in)    ::   ncol
    integer,  intent(in)    ::   lchnk
    integer,  intent(in)    ::   map(:)
    real(r8), intent(inout) ::   vmr(:,:,:)    ! lower bndy concentrations( mol/mol )

    !--------------------------------------------------------
    !	... local variables
    !--------------------------------------------------------
    integer  :: m, n
    integer  :: last, next
    real(r8) :: dels

    if( flbc_cnt < 1 ) then
       return
    end if

    call get_dels( dels, last, next )

    do m = 1,flbc_cnt
       if ( flbcs(m)%spc_ndx > 0 ) then
          n = map( flbcs(m)%spc_ndx )
          vmr(:ncol,pver,n) = flbcs(m)%vmr(:ncol,lchnk,last) &
               + dels * (flbcs(m)%vmr(:ncol,lchnk,next) - flbcs(m)%vmr(:ncol,lchnk,last))
       endif
    end do

  end subroutine flbc_set

  subroutine get_dels( dels, last, next )

    implicit none

    real(r8), intent(out) :: dels
    integer,  intent(out) :: last
    integer,  intent(out) :: next

    !--------------------------------------------------------
    !	... local variables
    !--------------------------------------------------------
    integer  ::  wrk_date, wrk_sec
    integer  ::  tcnt, n
    real(r8)     ::  wrk_time

    !--------------------------------------------------------
    !	... setup the time interpolation
    !--------------------------------------------------------
    wrk_sec  = ncsec
    select case( flbc_timing%type )
    case( 'SERIAL' )
       wrk_date = ncdate
    case( 'CYCLICAL' )
       wrk_date = flbc_timing%cycle_yr*10000 + mod( ncdate,10000 )
    case( 'FIXED' )
       wrk_date = flbc_timing%fixed_ymd
       wrk_sec  = flbc_timing%fixed_tod
    end select

    wrk_time = flt_date( wrk_date, wrk_sec )

    !--------------------------------------------------------
    !	... set time interpolation factor
    !--------------------------------------------------------
    if( flbc_timing%type /= 'CYCLICAL' ) then
       do n = tim_ndx(1)+1,tim_ndx(2)
          if( wrk_time <= times(n) ) then
             last = n - 1
             next = n
             exit
          end if
       end do
       if( n > ntimes ) then
          write(iulog,*) 'flbc_set: interp time is out of bounds'
          call endrun
       end if
       dels = (wrk_time - times(last))/(times(next) - times(last))
       !        write(iulog,*) ' '
       !        write(iulog,*) 'flbc_set: last,next,dels,ncdate,ncsec = ',last,next,dels,ncdate,ncsec
    else
       tcnt = tim_ndx(2) - tim_ndx(1) + 1
       call findplb( times(tim_ndx(1)), tcnt, wrk_time, n )
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
       !        write(iulog,*) ' '
       !        write(iulog,*) 'flbc_set: last,next,dels,ncdate,ncsec = ',last,next,dels,ncdate,ncsec
    end if

    dels = max( min( 1._r8,dels ),0._r8 )

  end subroutine get_dels

  subroutine flbc_gmean_vmr(co2vmr,ch4vmr,n2ovmr,f11vmr,f12vmr)

     implicit none

     real(r8), intent(inout) :: co2vmr
     real(r8), intent(inout) :: ch4vmr
     real(r8), intent(inout) :: n2ovmr
     real(r8), intent(inout) :: f11vmr
     real(r8), intent(inout) :: f12vmr

     integer  :: last, next
     real(r8) :: dels

     if( flbc_cnt < 1 ) return

     call get_dels( dels, last, next )

     if (ghg_indices(co2_ndx)>0) &
          co2vmr = global_mean_vmr(flbcs(ghg_indices(co2_ndx)), dels, last, next )
     if (ghg_indices(ch4_ndx)>0) &
          ch4vmr = global_mean_vmr(flbcs(ghg_indices(ch4_ndx)), dels, last, next )
     if (ghg_indices(n2o_ndx)>0) &
          n2ovmr = global_mean_vmr(flbcs(ghg_indices(n2o_ndx)), dels, last, next )
     if (ghg_indices(f11_ndx)>0) &
          f11vmr = global_mean_vmr(flbcs(ghg_indices(f11_ndx)), dels, last, next )
     if (ghg_indices(f12_ndx)>0) &
          f12vmr = global_mean_vmr(flbcs(ghg_indices(f12_ndx)), dels, last, next )

  end subroutine flbc_gmean_vmr

  function global_mean_vmr( flbcs, dels, last, next  )
    !use phys_gmean, only: gmean!BSINGH - Commented out due to circular dependency, see below
    use phys_grid,  only: get_ncols_p

    implicit none

    type(flbc), intent(in) :: flbcs
    real(r8), intent(in) :: dels
    integer, intent(in) :: last
    integer, intent(in) :: next
    real(r8) :: global_mean_vmr
    real(r8) :: vmr_arr(pcols,begchunk:endchunk)

    integer  :: lchnk, ncol !, n

    if (flbcs%has_mean) then
       global_mean_vmr = flbcs%vmr_mean(last) &
            + dels * (flbcs%vmr_mean(next) - flbcs%vmr_mean(last))
    else 
       do lchnk = begchunk, endchunk
          ncol = get_ncols_p(lchnk)
          vmr_arr(:ncol,lchnk) = flbcs%vmr(:ncol,lchnk,last) &
               + dels * (flbcs%vmr(:ncol,lchnk,next) - flbcs%vmr(:ncol,lchnk,last))
       enddo
       call endrun ('BALLI- circular dependency[physics_types->cam_history->mo_flbc->phys_gmean->physics_type]')
       !BALLIcall gmean (vmr_arr, global_mean_vmr)
    endif

  endfunction global_mean_vmr

end module mo_flbc
