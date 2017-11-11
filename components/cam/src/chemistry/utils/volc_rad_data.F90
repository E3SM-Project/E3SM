module volc_rad_data

  use pio,            only: pio_inq_dimid, pio_inq_varid, pio_nowrite, file_desc_t, &
       pio_inq_dimlen, pio_get_var
  
  use radconstants,   only: nswbands, nlwbands
  use shr_kind_mod,   only: shr_kind_cl,r8 => shr_kind_r8
  use cam_abortutils, only: endrun
  use time_manager,   only: set_time_float_from_date
  use cam_logfile,    only: iulog
  use phys_grid,      only: pcols, pver, begchunk, endchunk, get_ncols_p, get_rlat_all_p

  implicit none

  private
  save

  public :: volc_rad_data_init
  public :: advance_volc_rad_data

  type(file_desc_t)  :: piofile                       !input file handle

  character(len=shr_kind_cl) :: curr_filename         !name of the input file(to use in error messages)
  logical :: iscyclic
  integer, parameter :: ntslc = 2                  !number of time slices to use for time interpolation

  integer :: mxnflds, mxnflds_sw, mxnflds_lw, nlats, nalts, nalts_int, ntimes !lengths of various fields of input file
  integer :: cnt_sw(4), cnt_lw(4)                        !dimension of data to be read from input file for one time slice
  integer :: cyc_ndx_beg, cyc_ndx_end, cyc_tsize, cyc_yr_in
  integer, allocatable :: pbuf_idx_sw(:), pbuf_idx_lw(:) !buffer index of all the radiation fields in pbuf

  real(r8) :: neg_huge                                   !largest possible negative number

  real(r8), allocatable :: lats(:), alts(:), alts_int(:), times(:) !input file dimension values
  
contains
  
  subroutine volc_rad_data_init (specifier_sw, specifier_lw, filename, datapath, &
       data_type, cyc_yr)

    use mo_constants,    only: d2r
    use spmd_utils,      only: masterproc
    use ioFileMod,       only: getfil
    use cam_pio_utils,   only: cam_pio_openfile
    use physics_buffer,  only: pbuf_get_index

    implicit none
    
    character(len=*), intent(in) :: specifier_sw(:), specifier_lw(:)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: datapath
    character(len=*), intent(in) :: data_type
    integer,          intent(in) :: cyc_yr

    !Local variables
    character(len=shr_kind_cl) :: filen, filepath
    logical :: need_first_ndx

    integer :: ierr, ifld, errcode
    integer :: lat_dimid, lev_dimid, tim_dimid, alt_dimid, alt_int_dimid, sw_dimid, lw_dimid, old_dimid
    integer :: itimes
    integer :: nswbnds_prscb, nlwbnds_prscb, year, month, day, idates

    integer, allocatable :: dates(:)


    !BALLI- add comment

    !Currently only handles cyclic or serial data types
    iscyclic = .false.
    if (trim(data_type).ne. 'SERIAL' .and. trim(data_type).ne. 'CYCLICAL') then
       call endrun ('volc_rad_data.F90:Only SERIAL or CYCLICAL data types supported, current data type is:'//trim(data_type))
    endif
    if (trim(data_type).eq. 'CYCLICAL') then
       iscyclic = .true.
    endif
       
        
    curr_filename = trim(filename)
    cyc_yr_in     = cyc_yr
    neg_huge = -1.0_r8* huge(1.0)
    
    mxnflds_sw = size( specifier_sw )
    mxnflds_lw = size( specifier_lw )
    mxnflds = mxnflds_sw + mxnflds_lw

    if (mxnflds < 1) return
    
    if (len_trim(datapath) == 0) then
       filepath = trim(filename)
    else
       filepath = trim(datapath) // '/' // trim(filename)
    end if

    ! open file and get fileid 
    call getfil( filepath , filen, 0 )
    call cam_pio_openfile( piofile, filen, PIO_NOWRITE)
    if(masterproc) write(iulog,*)'open_volc_rad_datafile: ',trim(filen)


    !Polulate dates from the volc file
    ierr = pio_inq_dimid( piofile, 'month', old_dimid)
    ! Hack to work with weird netCDF and old gcc or NAG bug.
    tim_dimid = old_dimid
    ierr = pio_inq_dimlen( piofile, old_dimid, ntimes)
    allocate(dates(ntimes))
    allocate(times(ntimes))

    ierr = pio_inq_varid( piofile, 'month', old_dimid )
    ierr = pio_get_var( piofile, old_dimid, dates )
    
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
    
    if (cyc_ndx_beg < 0) then
       write(iulog,*) 'volc_rad_data.F90: subr volc_rad_data_init: cycle year not found : ' , cyc_yr
       call endrun('volc_rad_data_init:: cycle year not found')
    endif

    cyc_tsize = cyc_ndx_end - cyc_ndx_beg + 1 


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

    ierr = pio_inq_dimid( piofile, 'solar_bands', old_dimid)
    sw_dimid = old_dimid
    ierr = pio_inq_dimlen( piofile, old_dimid, nswbnds_prscb)

    ierr = pio_inq_dimid( piofile, 'terrestrial_bands', old_dimid)
    lw_dimid = old_dimid
    ierr = pio_inq_dimlen( piofile, old_dimid, nlwbnds_prscb)

    if(nswbands .ne. nswbnds_prscb .and. nlwbands .ne. nlwbnds_prscb) then
       write(iulog,*)'Bands in volc input file do not match the model bands'
       write(iulog,*)'Bands in volc input file:nswbands,nlwbands:',nswbnds_prscb,nlwbnds_prscb
       write(iulog,*)'Bands in model:nswbands,nlwbands:',nswbands,nlwbands
       call endrun ("volc_rad_data.F90: shortwave and longwave bands mismatch ")
    endif


    cnt_sw  = (/ 1, nalts, nlats, nswbands /) 
    cnt_lw  = (/ 1, nalts, nlats, nlwbands /) 

    allocate(pbuf_idx_sw(mxnflds_sw),pbuf_idx_lw(mxnflds_lw))

    do ifld = 1,mxnflds_sw
       pbuf_idx_sw(ifld) = pbuf_get_index(trim(specifier_sw(ifld)),errcode)!BALLI handle error?
    enddo 
    do ifld = 1,mxnflds_lw
       pbuf_idx_lw(ifld) = pbuf_get_index(trim(specifier_lw(ifld)),errcode)!BALLI handle error?
    enddo

  end subroutine volc_rad_data_init

  !------------------------------------------------------------------------------------------------


  subroutine advance_volc_rad_data (specifier_sw, specifier_lw, state, pbuf2d )
    
    use physics_buffer,   only: physics_buffer_desc
    use physics_types,    only: physics_state
    use time_manager,     only: get_curr_date

    !Arguments    
    character(len=*),    intent(in)    :: specifier_sw(:),specifier_lw(:)
    type(physics_state), intent(in)    :: state(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    !Local variables
    integer :: strt_t(4), strt_tp1(4), ifld, ierr, var_id
    integer :: errcode, yr, mon, day, ncsec, it, itp1, banddim

    real(r8) :: datatimem, datatimep, curr_mdl_time, deltat, fact1, fact2
    real(r8) :: wrk_sw(2,nswbands,nlats,nalts)
    real(r8) :: wrk_lw(2,nlwbands,nlats,nalts)

  !------------------------------------------------------------------------------------------------

    !find time indices to read data for two consecutive time indices to interpolate in time
    !get current model date
    call get_curr_date(yr, mon, day, ncsec)
    if(iscyclic) yr = cyc_yr_in
 
    call set_time_float_from_date( curr_mdl_time, yr, mon, day, ncsec )
    
    call find_times_to_interpolate(curr_mdl_time, datatimem, datatimep, it, itp1)

    else
       deltat = datatimep - datatimem
       fact1 = (datatimep - curr_mdl_time)/deltat
       fact2 = 1._r8-fact1
    end if

    !see if we need to read data (i.e. we are at a time where new data is needed)
    strt_t = (/ it, 1, 1, 1 /)
    strt_tp1 = (/ itp1, 1, 1, 1 /)


    do ifld = 1,mxnflds_sw
       !Get netcdf variable id for the field
       ierr = pio_inq_varid(piofile, trim(adjustl(specifier_sw(ifld))),var_id) !get id of the variable to read from iput file
       ierr = pio_get_var( piofile, var_id, strt_t  , cnt_sw, wrk_sw(1,:,:,:) )!BALLI: handle error?
       ierr = pio_get_var( piofile, var_id, strt_tp1, cnt_sw, wrk_sw(2,:,:,:) )!BALLI: handle error?
       !interpolate in lats, time and vertical
       call interpolate_lats_time_vert(state, wrk_sw, nswbands, pbuf_idx_sw(ifld), fact1, fact2, pbuf2d )
    enddo

    do ifld = 1,mxnflds_lw
       !Note that we have to reverse (compared with how they are mentioned in the netcdf file) the array dim sizes
       ierr = pio_inq_varid(piofile, trim(adjustl(specifier_lw(ifld))),var_id) !get id of the variable to read from iput file
       ierr = pio_get_var( piofile, var_id, strt_t  , cnt_lw, wrk_lw(1,:,:,:) )!BALLI: handle error?
       ierr = pio_get_var( piofile, var_id, strt_tp1, cnt_lw, wrk_lw(2,:,:,:) )!BALLI: handle error?
       !interpolate in lats, time and vertical
       call interpolate_lats_time_vert(state, wrk_lw, nlwbands, pbuf_idx_lw(ifld), fact1, fact2, pbuf2d )
    enddo
    
        
  end subroutine advance_volc_rad_data

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
             !PMC bugfix: in this case, target x is smaller than the 1st source x,
             !so it=1 and itp1=cyc_tsize. This was causing weird values for interpolation
             !solution: allow datatimep to be a negative number.
             it   = n_out + cyc_ndx_beg-1
             itp1 = cyc_ndx_beg          
             datatimem = times(it)
             datatimep = times(itp1) - 365.
             !KLUDGE: line above hardcodes assumption that times are in days!
             times_found = .true.
          else
             np1 = n_out+1
             it   = n_out + cyc_ndx_beg-1
             itp1 = n_out  + cyc_ndx_beg          
             datatimem = times(it)
             datatimep = times(itp1)
             times_found = .true.
          endif

       else
          call endrun('cyclic method for cyc_tsize<=1 not implemented yet')
       endif
    endif
       
    if(.not.times_found) call endrun ('volc_rad_data.F90: sub find_times_to_interpolate: times not found!')

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
    real(r8), intent(in) :: fact1, fact2, wrk(ntslc,banddim,nlats,nalts)

    !--intent (out)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    !Local variables
    type(interp_type) :: lat_wgts
    integer  :: ichnk, ncols, k, iband, icol

    real(r8), parameter :: m2km  = 1.e-3_r8

    real(r8) :: wrk_1d(ntslc,pcols), to_lats(pcols)
    real(r8) :: loc_arr(ntslc,banddim,pcols,nalts), datain(banddim,pcols,nalts)
    real(r8) :: model_z(pverp)

    real(r8), pointer :: data_out(:,:,:)    

    do ichnk = begchunk, endchunk
       call pbuf_get_field(pbuf2d, ichnk, pbuf_idx, data_out)
       ncols = get_ncols_p(ichnk)
       call get_rlat_all_p(ichnk, pcols, to_lats)
       call lininterp_init(lats, nlats, to_lats, ncols, 1, lat_wgts)
       do k = 1, nalts
          do iband = 1 , banddim 
             !lats interpolation
             call lininterp(wrk(1,iband,:,k), nlats, wrk_1d(1,1:ncols), ncols, lat_wgts)
             call lininterp(wrk(2,iband,:,k), nlats, wrk_1d(2,1:ncols), ncols, lat_wgts)
             loc_arr(1,iband,1:ncols,k) = wrk_1d(1,1:ncols)
             loc_arr(2,iband,1:ncols,k) = wrk_1d(2,1:ncols)
             !time interpolation                                                                                                                                                                                              
             datain(iband,1:ncols,k) = fact1*loc_arr(1,iband,1:ncols,k) + fact2*loc_arr(2,iband,1:ncols,k)
          end do
       enddo       
       call lininterp_finish(lat_wgts)
       !vertical interpolation                                                                                                                                                                                 
       do icol = 1, ncols
          model_z(1:pverp) = m2km * state(ichnk)%zi(icol,pverp:1:-1)
          do iband = 1, banddim
             call rebin( nalts, pver, alts_int, model_z, datain(iband,icol,:), data_out(iband,icol,:) )
          enddo
       enddo
    end do

  end subroutine interpolate_lats_time_vert
  
end module volc_rad_data
