module read_volc_radiation_data

!To do list:
!1. Handle ALL allocate errors
!2. Add checks for bands so that model and read in input file bands match!
!3. Logic for cyclic time interpolation where we substracted "one_yr" should be generalized
!4. cleanup
!5. TESTING TESTING TESTING!!

!-------------------------------------------------------------------------------------------------------
! What these codes do:
! ===================
! These F90 codes read radiation data directly from the prescribed volcanic file. CMIP6 
! provided this file so that extinction (ext), single scattering albedos(ssa) and asymmetric factors (af)
! can be directly read from this prescribed file and fed into the model rather than 
! deriving these quantities from volcanic aerosol mixing ratios (VOLC_MMR). These quantites
! are used to update model variables only above the tropopause (see aer_rad_props.F90 
! and modal_aer_opt.F90, look for is_cmip6_volc logical variable) and only 50% contributions
! at the trpopause layer
!
!-------------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------------
!
! Some shorthand explained:
! sw  = short wave
! lw  = long wave
! ext = volcanic aerosols extinction
! ssa = volcanic aerosols single scattering albedo
! af  = volcanic aerosols asymmetric factors
! "input file" = prescribed volcanic input file from where 
!-------------------------------------------------------------------------------------------------------

  use pio,            only: pio_inq_dimid, pio_inq_varid, pio_nowrite, file_desc_t, &
       pio_inq_dimlen, pio_get_var
  
  use radconstants,   only: nswbands, nlwbands                                    !Number of bands in sw and lw
  use shr_kind_mod,   only: shr_kind_cl,r8 => shr_kind_r8
  use cam_abortutils, only: endrun
  use time_manager,   only: set_time_float_from_date
  use cam_logfile,    only: iulog
  use phys_grid,      only: pcols, pver, begchunk, endchunk, get_ncols_p, get_rlat_all_p

  implicit none

  private
  save

  public :: read_volc_radiation_data_init
  public :: advance_volc_radiation_data

  type(file_desc_t)  :: piofile                          !input file handle

  character(len=shr_kind_cl) :: curr_filename            !name of the input file(to use in error messages)
  logical :: iscyclic
  integer, parameter :: ntslc = 2                        !number of time slices to use for time interpolation

  integer :: mxnflds, mxnflds_sw, mxnflds_lw, nlats 
  integer :: nalts, nalts_int, ntimes                    !lengths of various fields of input file
  integer :: cnt_sw(4), cnt_lw(4)                        !dimension of data to be read from input file for one time slice
  integer :: cyc_ndx_beg, cyc_ndx_end, cyc_tsize         !variables for cyclical input file
  integer :: cyc_yr_in
  integer :: strt_t(4), strt_tp1(4)
  integer, allocatable :: pbuf_idx_sw(:), pbuf_idx_lw(:) !buffer index of all the radiation fields in pbuf

  real(r8) :: neg_huge, one_yr                           !largest possible negative number;days in a year
  real(r8) :: deltat, datatimem, datatimep               !two time indices for time interpolation and their difference

  real(r8), allocatable :: wrk_sw(:,:,:,:,:)               !work arrays for sw and lw
  real(r8), allocatable :: wrk_lw(:,:,:,:,:)               !PMC identified a bug here with a wrong order of dims, fixed now

  real(r8), allocatable :: lats(:), alts(:), alts_int(:) !input file dimension values
  real(r8), allocatable :: times(:) 
  
  
contains
  
  subroutine read_volc_radiation_data_init (specifier_sw, specifier_lw, filename, datapath, &
       data_type, cyc_yr)

    use mo_constants,    only: d2r        !degree to radians conversion
    use spmd_utils,      only: masterproc 
    use ioFileMod,       only: getfil
    use cam_pio_utils,   only: cam_pio_openfile
    use physics_buffer,  only: pbuf_get_index

    implicit none
    
    character(len=*), intent(in) :: specifier_sw(:), specifier_lw(:) !names of fields in sw and lw to be read from the input file
    character(len=*), intent(in) :: filename                         !name of the input file
    character(len=*), intent(in) :: datapath                         !path to input file
    character(len=*), intent(in) :: data_type                        !SERAIL or CYCLICAL
    integer,          intent(in) :: cyc_yr                           !cycle year, only used if data_type is CYCLICAL

    !Local variables
    character(len=shr_kind_cl) :: filen, filepath
    character(len = 512) :: msg
    logical :: need_first_ndx 

    integer :: ierr, ifld, errcode
    integer :: lat_dimid, lev_dimid, tim_dimid, alt_dimid, alt_int_dimid, sw_dimid, lw_dimid, old_dimid !dimension ids
    integer :: itimes
    integer :: nswbnds_prscb, nlwbnds_prscb !band numbers from the precribed input file
    integer :: year, month, day, idates

    integer, allocatable :: dates(:)

    real(r8) :: time1, time2                !times to compute "one_yr"


    !BALLI- add comment

    !Currently only handles cyclic or serial data types; exit if some other data_type detected
    iscyclic = .false.
    if (trim(data_type).ne. 'SERIAL' .and. trim(data_type).ne. 'CYCLICAL') then
       msg = 'read_volc_radiation_data.F90:Only SERIAL or CYCLICAL data types supported, current data type is:'//trim(data_type)
       call endrun (msg)
    endif
    if (trim(data_type).eq. 'CYCLICAL') then
       iscyclic = .true.
    endif
       
        
    curr_filename = trim(filename)
    cyc_yr_in     = cyc_yr
    neg_huge      = -1.0_r8* huge(1.0_r8)
    
    datatimem     = neg_huge
    datatimep     = neg_huge

    strt_t   = (/ -1, -1, -1, -1 /)
    strt_tp1 = (/ -1, -1, -1, -1 /)
    
    mxnflds_sw = size( specifier_sw )!size of sw fields in input file
    mxnflds_lw = size( specifier_lw )!size of lw fields in input file
    mxnflds = mxnflds_sw + mxnflds_lw

    if (mxnflds < 1) return          !return if no field needs to be read
    
    if (len_trim(datapath) == 0) then
       filepath = trim(filename)
    else
       filepath = trim(datapath) // '/' // trim(filename)
    end if

    ! open file and get fileid 
    call getfil( filepath , filen, 0 )
    call cam_pio_openfile( piofile, filen, PIO_NOWRITE)
    if(masterproc) write(iulog,*)'open_read_volc_radiation_datafile: ',trim(filen)


    !Polulate dates from the volc file, dates stamps are stored in variable "month" in the file
    ierr = pio_inq_dimid( piofile, 'month', old_dimid)
    ! Hack to work with weird netCDF and old gcc or NAG bug.
    tim_dimid = old_dimid
    ierr = pio_inq_dimlen( piofile, old_dimid, ntimes)
    allocate(dates(ntimes))
    allocate(times(ntimes))

    ierr = pio_inq_varid( piofile, 'date', old_dimid )
    ierr = pio_get_var( piofile, old_dimid, dates )
    
    call set_time_float_from_date( time2, 2, 1, 1, 0 )
    call set_time_float_from_date( time1, 1, 1, 1, 0 )
    one_yr = time2-time1

    if (iscyclic) then
       cyc_ndx_beg = -1
       cyc_ndx_end = -1
       if ( cyc_yr /= 0 ) then
          call set_time_float_from_date( time1, cyc_yr  , 1, 1, 0 )
          call set_time_float_from_date( time2, cyc_yr+1, 1, 1, 0 )
          one_yr = time2-time1
       endif
    endif
    
    !compute times array
    need_first_ndx = .true. !for cyclic year only
    
    do itimes = 1, ntimes
       idates = dates(itimes)
       year   = idates / 10000
       month  = mod(idates, 10000)/100
       day    = mod(idates, 100)
       call set_time_float_from_date( times(itimes), year, month, day, 0 )
       if (iscyclic) then
          if ( year == cyc_yr ) then
             if ( need_first_ndx ) then
                cyc_ndx_beg = itimes
                need_first_ndx = .false.
             endif
             cyc_ndx_end = itimes     
          endif
       endif
    enddo
    
    if(iscyclic) then
       if (cyc_ndx_beg < 0) then
          write(iulog,*) 'read_volc_radiation_data.F90: subr read_volc_radiation_data_init: cycle year not found : ' , cyc_yr
          call endrun('read_volc_radiation_data_init:: cycle year not found')
       endif

       cyc_tsize = cyc_ndx_end - cyc_ndx_beg + 1 
    endif


    !Polulate lats from the volc file
    ierr = pio_inq_dimid( piofile, 'latitude', old_dimid)
    lat_dimid = old_dimid
    ierr = pio_inq_dimlen( piofile, old_dimid, nlats)
    allocate(lats(nlats))

    ierr = pio_inq_varid( piofile, 'latitude', old_dimid )
    ierr = pio_get_var( piofile, old_dimid, lats )
    lats = lats * d2r !degree to radians

    !Polulate alts from the volc file
    ierr = pio_inq_dimid( piofile, 'altitude', old_dimid)
    alt_dimid = old_dimid
    ierr = pio_inq_dimlen( piofile, old_dimid, nalts)
    allocate(alts(nalts))

    ierr = pio_inq_varid( piofile, 'altitude', old_dimid )
    ierr = pio_get_var( piofile, old_dimid, alts )

    ierr = pio_inq_dimid( piofile, 'altitude_int', old_dimid)
    alt_int_dimid = old_dimid
    ierr = pio_inq_dimlen( piofile, old_dimid, nalts_int)
    allocate(alts_int(nalts_int))

    ierr = pio_inq_varid( piofile, 'altitude_int', old_dimid )
    ierr = pio_get_var( piofile, old_dimid, alts_int )

    allocate(wrk_sw(2,mxnflds_sw,nalts,nlats,nswbands))!BALLI- handle erros for allocate everywhere!!
    allocate(wrk_lw(2,mxnflds_lw,nalts,nlats,nlwbands))


    ierr = pio_inq_dimid( piofile, 'solar_bands', old_dimid)
    sw_dimid = old_dimid
    ierr = pio_inq_dimlen( piofile, old_dimid, nswbnds_prscb)

    ierr = pio_inq_dimid( piofile, 'terrestrial_bands', old_dimid)
    lw_dimid = old_dimid
    ierr = pio_inq_dimlen( piofile, old_dimid, nlwbnds_prscb)

    if(nswbands .ne. nswbnds_prscb .or. nlwbands .ne. nlwbnds_prscb) then
       write(iulog,*)'Bands in volc input file do not match the model bands'
       write(iulog,*)'Bands in volc input file:nswbands,nlwbands:',nswbnds_prscb,nlwbnds_prscb
       write(iulog,*)'Bands in model:nswbands,nlwbands:',nswbands,nlwbands
       call endrun ("read_volc_radiation_data.F90: shortwave and longwave bands mismatch ")
    endif


    cnt_sw  = (/ 1, nalts, nlats, nswbands /) 
    cnt_lw  = (/ 1, nalts, nlats, nlwbands /) 

    allocate(pbuf_idx_sw(mxnflds_sw),pbuf_idx_lw(mxnflds_lw))

    do ifld = 1, mxnflds_sw
       pbuf_idx_sw(ifld) = pbuf_get_index(trim(specifier_sw(ifld)),errcode)!BALLI handle error?
    enddo 
    do ifld = 1, mxnflds_lw
       pbuf_idx_lw(ifld) = pbuf_get_index(trim(specifier_lw(ifld)),errcode)!BALLI handle error?
    enddo

  end subroutine read_volc_radiation_data_init

  !------------------------------------------------------------------------------------------------


  subroutine advance_volc_radiation_data (specifier_sw, specifier_lw, state, pbuf2d )
    !------------------------------------------------------------------------------------------------
    !Logic:
    !Data is read from the prescribed input file and advanced in time. We do not need to read data each
    !time step but we do need to interpolate the data at each time step as "current model time" changes 
    !with each time step. "current model time" (variable name:curr_mdl_time) is needed for time 
    !interpolation, so we have to do time interpolation every time step. Vertical interpolation follows
    !time interpolation, so we have to do vertical interpolation each time step as well.
    
    use physics_buffer,   only: physics_buffer_desc
    use physics_types,    only: physics_state
    use time_manager,     only: get_curr_date
    use spmd_utils,       only: masterproc

    !Arguments    
    character(len=*),    intent(in)    :: specifier_sw(:),specifier_lw(:)
    type(physics_state), intent(in)    :: state(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    !Local variables
    logical :: read_data
    integer :: ifld, ierr, var_id
    integer :: errcode, yr, mon, day, ncsec, it, itp1, itp1_old, banddim

    real(r8) :: data_time, curr_mdl_time, fact1, fact2, deltat_cyc

  !------------------------------------------------------------------------------------------------

    !The pupose here is to find time indices to read data for two consecutive time indices to interpolate in time

    !Get current model date. 
    call get_curr_date(yr, mon, day, ncsec)
    !if cyclic year, set year to cyclic year
    if(iscyclic) yr = cyc_yr_in
 
    !get current model time in "days".Note: (year will be cyclic year if iscyclic is TRUE)
    !Output from the following call is "curr_mdl_time", which is in days, it can be fraction 
    !based on how far we are in a given day
    call set_time_float_from_date( curr_mdl_time, yr, mon, day, ncsec )

    !Generally data files has data for each month (it can also be more of less frequent)
    !"times" variable has all the times (day,month,year) for which we have data in the file.
    !Note that values of "times" array are also "number of days" starting from 0001/01/01

    !We only need to read data when we are at a "data time boundary" where we need next
    !time slice to interpolate

    !"datatimep" below is "plus" end of 2 time slices (i.e. 2nd time slice) to interpolate
    !Lets assign "data_time" to "datatimep" and check later if "curr_mdl_time" is greater
    ! that this slice or not. If it is, we would like to read next time sslice
    data_time = datatimep  
    !if cyclic, we might need to treat the yearly boundaries in a special manner
    if ( iscyclic) then
       ! wrap around                                                                                                                                              
       if ( (datatimep<datatimem) .and. (curr_mdl_time>datatimem) ) then
          !This "if" condition is here so that we DON'T read data when it isn't necessary at the
          !yearly boundary. To understand it, lets assume, we start from year 1 and we have data
          !at every 15th of the month in the file. At 12/15, we will read data from slices 12/15 and 01/15 of year 1 
          !and interpolate. Note that at 12/15, this "if" condition will not be executed as curr_mdl_time
          !is == datatimem (not greater than). At 12/16, we DO NOT want to read next time slice, so this
          !if condition will be executed to increment "data_time" by "one_yr" so that curr_mdl_time is
          !NOT greater than data_time (see next if condition). NOTE: from 12/16 to 12/31, this if 
          !condition will be executed because (datatimep<datatimem) is TRUE and (curr_mdl_time>datatimem)
          ! is also TRUE. Here datatimep and datatimem are days at 01/15 and 12/15 for year 1 respectively. 
          ! At 01/01, curr_mdl_time will start from 0, therefore, (curr_mdl_time>datatimem) condition
          ! will be false, so this if condition will not be executed. Note that next time slice will 
          !still not be read (as it should not be) at 01/01 as curr_mdl_time is 0 and data_time is 
          !datatimep, which is 01/15 (see next if condition). Therefore curr_mdl_time is NOT greater than
          !data_time through 01/01 to 01/14
          data_time = data_time + one_yr
       endif
    endif

    !if current model time is greater than or equal to second time index, we must move forward 
    !to increment time indices (first and second) and read the input file again
    read_data = .false.
    if (curr_mdl_time >= data_time ) then
       ! Save previous itp1 (will be -1 the first time step)
       itp1_old = strt_tp1(1)
       ! Move forward and get new values for datatimep,datatimem, it and itp1
       call find_times_to_interpolate(curr_mdl_time, datatimem, datatimep, it, itp1)
       deltat = datatimep - datatimem
       ! Two time indices for reading data for time interpolation
       strt_t   = (/ it,   1, 1, 1 /) !index for first time stamp
       strt_tp1 = (/ itp1, 1, 1, 1 /) !index for second time stamp

       !set read_data to true so that new data is read from the file
       read_data = .true.
    endif

    fact1 = (datatimep - curr_mdl_time)/deltat    
    if (iscyclic .and. deltat < 0.0_r8) then
       !if deltat is negative (it generally happens at yearly boundaries or at start of the simulation)
       ! recompute deltat (as deltat_cyc) and fact1:

       !Add a year, this will make deltat_cyc positive and it will be exact numder of days difference between
       !two time slices at yearly boundary (say between 12/15 to 01/15)
       deltat_cyc = deltat + one_yr
       if ( datatimep >= curr_mdl_time ) then
          ! if datatimep is greater than curr_mdl_time, it will be true between 01/01 to 01/15 (considering
          ! the example above), compute the fact1 using following similar formula except use deltat_cyc
          fact1 = (datatimep - curr_mdl_time)/deltat_cyc
       else
          ! if datatimep is less than curr_mdl_time, we are between 12/15 to 12/31, we have to add
          ! one_yr to datatimep so that fact1 is positive and interpolation is done assuming datatimep
          !is in next year (if it makes sense). Divide by deltat_cyc here as well
          fact1 = (datatimep+one_yr - curr_mdl_time)/deltat_cyc
       endif
    endif


    fact2 = 1._r8-fact1

    do ifld = 1,mxnflds_sw
       if (read_data) then
          !Get netcdf variable id for the field
          ierr = pio_inq_varid(piofile, trim(adjustl(specifier_sw(ifld))),var_id) !get id of the variable to read from iput file
          ! If current itp matches previous itp1, copy rather than reading from file
          if (it == itp1_old) then
             wrk_sw(1,ifld,:,:,:) = wrk_sw(2,ifld,:,:,:)
          else
             ierr = pio_get_var( piofile, var_id, strt_t, cnt_sw, wrk_sw(1,ifld,:,:,:) )!BALLI: handle error?
          endif
          ierr = pio_get_var( piofile, var_id, strt_tp1, cnt_sw, wrk_sw(2,ifld,:,:,:) )!BALLI: handle error?
       endif

       !we always have to do interpolation as current model time changes time factors-fact1 and fact2
       !interpolate in lats, time and vertical
       call interpolate_lats_time_vert(state, wrk_sw(:,ifld,:,:,:), nswbands, pbuf_idx_sw(ifld), fact1, fact2, pbuf2d )
    enddo

    do ifld = 1,mxnflds_lw
       !Note that we have to reverse (compared with how they are mentioned in the netcdf file) the array dim sizes
       if(read_data) then
          ierr = pio_inq_varid(piofile, trim(adjustl(specifier_lw(ifld))),var_id) !get id of the variable to read from iput file
          ! If current itp matches previous itp1, copy rather than reading from file
          if (it == itp1_old) then
             wrk_lw(1,ifld,:,:,:) = wrk_lw(2,ifld,:,:,:)
          else
             ierr = pio_get_var( piofile, var_id, strt_t, cnt_lw, wrk_lw(1,ifld,:,:,:) )!BALLI: handle error?
          endif
          ierr = pio_get_var( piofile, var_id, strt_tp1, cnt_lw, wrk_lw(2,ifld,:,:,:) )!BALLI: handle error?
       endif
       !we always have to do interpolation as current model time changes time factors-fact1 and fact2 
       !interpolate in lats, time and vertical
       call interpolate_lats_time_vert(state, wrk_lw(:,ifld,:,:,:), nlwbands, pbuf_idx_lw(ifld), fact1, fact2, pbuf2d )
    enddo

    ! Chris Golaz: for debugging
    !if (masterproc) write(iulog,'(a,i5,2i3,i6,3f9.2,l2,2f8.4)')'advance_volc_radiation_data:', &
    !  yr,mon,day,ncsec,curr_mdl_time,datatimem,datatimep,read_data,fact1,fact2
        
  end subroutine advance_volc_radiation_data

  !------------------------------------------------------------------------------------------------

  subroutine find_times_to_interpolate(curr_mdl_time, datatimem, datatimep, it, itp1)
    
    !Arguments
    real(r8), intent(in)  :: curr_mdl_time

    integer,  intent(out) :: it, itp1
    real(r8), intent(out) :: datatimem, datatimep
    
    !Local variables
    logical :: times_found
    integer :: n_out, np1

    datatimem = neg_huge
    datatimep = neg_huge

    times_found = .false.

    if(.not. iscyclic) then
       if ( all( times(:) > curr_mdl_time ) ) then
          write(iulog,*) 'FIND_TIMES: ALL data times are after ', curr_mdl_time
          write(iulog,*) 'FIND_TIMES: data times: ',times(:)
          write(iulog,*) 'FIND_TIMES: time: ',curr_mdl_time
          call endrun('find_times: all(times(:) > curr_mdl_time) '// trim(curr_filename) )
       endif
       find_times_loop : do it=1, ntimes-1
          itp1 = it + 1
          datatimem = times(it)
          datatimep = times(itp1)
          if ( (curr_mdl_time .ge. datatimem) .and. (curr_mdl_time .lt. datatimep) ) then
             times_found = .true.
             exit find_times_loop
          endif
       enddo find_times_loop
    else !if cyclic
       if ( cyc_tsize > 1 ) then
          call findplb(times(cyc_ndx_beg:cyc_ndx_end),cyc_tsize, curr_mdl_time, n_out )
          if (n_out == cyc_tsize) then
             np1 = 1
          else
             np1 = n_out+1
          endif
          it   = n_out + cyc_ndx_beg-1
          itp1 = np1   + cyc_ndx_beg-1
          datatimem = times(it)
          datatimep = times(itp1)
          times_found = .true.
       else
          call endrun('cyclic method for cyc_tsize<=1 not implemented yet')
       endif
    endif
       
    if(.not.times_found) call endrun ('read_volc_radiation_data.F90: sub find_times_to_interpolate: times not found!')

  end subroutine find_times_to_interpolate

  !------------------------------------------------------------------------------------------------

  subroutine interpolate_lats_time_vert(state, wrk, banddim, pbuf_idx, fact1, fact2, pbuf2d)

    use ppgrid,           only: pverp
    use physics_buffer,   only: pbuf_get_field, physics_buffer_desc
    use physics_types,    only: physics_state
    use interpolate_data, only: lininterp_init, lininterp, interp_type, lininterp_finish    
    use mo_util,          only: rebin

    !Arguments
    !--intent(in)
    type(physics_state), intent(in) :: state(begchunk:endchunk)
    integer,  intent(in) :: banddim, pbuf_idx
    real(r8), intent(in) :: fact1, fact2, wrk(ntslc,nalts,nlats,banddim)

    !--intent (out)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    !Local variables
    type(interp_type) :: lat_wgts
    integer  :: ichnk, ncols, k, iband, icol

    real(r8), parameter :: m2km  = 1.e-3_r8

    real(r8) :: wrk_1d(ntslc,pcols), to_lats(pcols)
    real(r8) :: datain(pcols,nalts,banddim)
    real(r8) :: model_z(pverp),data_out_tmp(pcols,pver,banddim)    

    real(r8), pointer :: data_out(:,:,:)

    do ichnk = begchunk, endchunk
       call pbuf_get_field(pbuf2d, ichnk, pbuf_idx, data_out)
       ncols = get_ncols_p(ichnk)
       call get_rlat_all_p(ichnk, pcols, to_lats)
       call lininterp_init(lats, nlats, to_lats, ncols, 1, lat_wgts)
       do iband = 1 , banddim 
          do k = 1, nalts
             !lats interpolation
             call lininterp(wrk(1,k,:,iband), nlats, wrk_1d(1,1:ncols), ncols, lat_wgts)
             call lininterp(wrk(2,k,:,iband), nlats, wrk_1d(2,1:ncols), ncols, lat_wgts)

             !time interpolation
             datain(1:ncols,k,iband) = fact1*wrk_1d(1,1:ncols) + fact2*wrk_1d(2,1:ncols)
          end do
       enddo       
       call lininterp_finish(lat_wgts)
       !vertical interpolation
       do iband = 1, banddim
          do icol = 1, ncols
             !convert model's vertical coordinate from m to km and flip it to match data
             model_z(1:pverp) = m2km * state(ichnk)%zi(icol,pverp:1:-1)
             call rebin( nalts, pver, alts_int, model_z, datain(icol,:,iband), data_out_tmp(icol,:,iband) )
             !flip in vertical to match model
             data_out(icol,:,iband) = data_out_tmp(icol,pver:1:-1,iband)
          enddo
       enddo
    end do

  end subroutine interpolate_lats_time_vert
  
end module read_volc_radiation_data
