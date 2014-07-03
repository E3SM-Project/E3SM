!----------------------------------------------------------------------- 
!
! BOP
!
! !MODULE: sst_data
!
! !DESCRIPTION:	Module to handle dealing with the Sea-Surface Temperature 
!	        datasets.  This module also figures out the location of
!	        sea-ice from these datasets where it is assumed that 
! 	        seawater at freezing or below is a flag for the existence of sea-ice.
!	        SST datasets that are created for use with the stand-alone CCM should
!	        take this into account and set grid-points where sea-ice fraction is
! 	        greater than 50% to -1.8C and ensure that other grid points where sea-ice
!	        is less than 50% have SST's greater than -1.8C.
!
! Public interfaces:
!
!	sstini -- Initialization and reading of dataset.
!	sstint -- Interpolate dataset SST to current time.
!
!----------------------------------------------------------------------- 

module sst_data
!
! USES:
!
  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_scam_mod, only: shr_scam_GetCloseLatLon

  use rgrid,           only: nlon, fullgrid
  use pmgrid,          only: plon, plat
  use ppgrid,          only: pcols, begchunk, endchunk
  use phys_grid,       only: get_ncols_p, get_rlat_all_p, get_rlon_all_p,&
                             scatter_field_to_chunk
  use abortutils,      only: endrun
  use cam_control_mod, only: aqua_planet
  use error_messages,  only: alloc_err
  use interpolate_data,only: get_timeinterp_factors

  use scamMod,         only: scmlon,scmlat,isrestart,single_column

  use ocn_spmd
  use ocn_time_manager,only: get_curr_date, get_curr_calday, &
                             is_perpetual, get_perp_date, &
                             get_step_size, is_first_step
  use physconst,       only: pi
  use cam_logfile,     only: iulog
  use pio

  implicit none
  private ! By default all data is private to this module
  save
!
! ! PUBLIC DATA:
!
  logical, public :: sstcyc  ! If false, do not cycle sst dataset(assume multiyear)
  public          :: sst
!
!
! ! PUBLIC MEMBER FUNCTIONS:
  !
  public sstini   ! Initialization
  public sstint   ! Time interpolation of SST data

!===============================================================================
!EOP
!===============================================================================
!----------------------------------------------------------------------- 
! PRIVATE: Everything else is private to this module
!----------------------------------------------------------------------- 
  real(r8), parameter :: daysperyear = 365.0_r8  ! Number of days in a year

  real(r8), allocatable, dimension(:,:,:) :: &
      sstbdy         ! SST values on boundary dataset (pcols,begchunk:endchunk,2)
  real(r8), allocatable, dimension(:,:) :: &
      sst            ! Interpolated model sst values (pcols,begchunk:endchunk)
  real(r8) :: cdaysstm   ! Calendar day for prv. month SST values read in
  real(r8) :: cdaysstp   ! Calendar day for nxt. month SST values read in
      
  integer :: nm,np   ! Array indices for prv., nxt month sst data
  integer :: lonsiz  ! size of longitude dimension on sst dataset
  integer :: levsiz  ! size of level dimension on sst dataset
  integer :: latsiz  ! size of latitude dimension on sst dataset
  integer :: timesiz ! size of time dimension on sst dataset
  integer :: np1     ! current forward time index of sst dataset
  integer, allocatable :: date_sst(:)! Date on sst dataset (YYYYMMDD)
  integer, allocatable :: sec_sst(:) ! seconds of date on sst dataset (0-86399)
  integer :: ret
  integer :: closelatidx,closelonidx
  real(r8):: srfdata
  real(r8):: closelat,closelon

  real(r8), parameter :: tsice = -1.7999_r8 ! Freezing point of sea ice degrees C
                                            ! Use this with global sst data
  character(len=*), parameter :: fieldname='SST_cpl' 
!===============================================================================
CONTAINS
!===============================================================================

!======================================================================
! PUBLIC ROUTINES: Following routines are publically accessable
!======================================================================

!----------------------------------------------------------------------- 
! 
! BOP
!
! !IROUTINE: sstini
!
! !DESCRIPTION:
!
! Initialize the procedure for specifying sea surface temperatures
! Do initial read of time-varying sst boundary dataset, reading two
! consecutive months on either side of the current model date.
!
! Method: 
! 
! Author: L.Bath
! 
!-----------------------------------------------------------------------
!
! !INTERFACE
!
subroutine sstini(ncid_sst)
  use ncdio_atm, only : infld
!
! EOP
!
!---------------------------Common blocks-------------------------------
  type(file_desc_t), intent(inout) :: ncid_sst
!---------------------------Local variables-----------------------------
  integer dtime                 ! timestep size [seconds]
  integer dateid                ! netcdf id for date variable
  integer secid                 ! netcdf id for seconds variable
  integer londimid              ! netcdf id for longitude variable
  integer latdimid              ! netcdf id for latitude variable
  integer latid                 ! netcdf id for latitude variable
  integer timeid                ! netcdf id for time variable
  integer nlonid                ! netcdf id for nlon variable (rgrid)
  integer cnt3(3)               ! array of counts for each dimension
  integer strt3(3)              ! array of starting indices
  integer n                     ! indices
  integer nlon_sst(plat)        ! number of lons per lat on bdy dataset
  integer j                     ! latitude index
  integer istat                 ! error return
  integer :: yr, mon, day       ! components of a date
  integer :: ncdate             ! current date in integer format [yyyymmdd]
  integer :: ncsec              ! current time of day [seconds]
  integer :: ret                ! return code
  real(r8) calday               ! calendar day (includes yr if no cycling)
  real(r8) caldayloc            ! calendar day (includes yr if no cycling)
  real(r8) xvar(plon,plat,2)    ! work space 
  integer :: ierr
  logical :: readvar

!-----------------------------------------------------------------------
  !
  ! Initialize time indices
  !
  nm = 1
  np = 2
  !
  ! Allocate space for data.
  !
  !
  if (.not.allocated(sst)) then 
     allocate( sst(pcols,begchunk:endchunk), stat=istat )
     call alloc_err( istat, 'sstini', 'sst', &
	  pcols*(endchunk-begchunk+1) )
  endif

  if(aqua_planet) then
     call prescribed_sst()
     return
  end if
  if (.not.allocated(sstbdy)) then 
     allocate( sstbdy(pcols,begchunk:endchunk,2), stat=istat )
     call alloc_err( istat, 'sstini', 'sstbdy', &
	  pcols*(endchunk-begchunk+1)*2 )
  endif

  !
  ! Use year information only if not cycling sst dataset
  !
  if (is_first_step()) then
     dtime = get_step_size()
     dtime = -dtime
     calday = get_curr_calday(offset=dtime)
  else
     calday = get_curr_calday()
  endif
  if ( is_perpetual() ) then
     call get_perp_date(yr, mon, day, ncsec)
  else
     if (is_first_step()) then
        call get_curr_date(yr, mon, day, ncsec,offset=dtime)
     else
        call get_curr_date(yr, mon, day, ncsec)
     endif
  end if
  ncdate = yr*10000 + mon*100 + day
  if (sstcyc) then
     calday = get_curr_calday()
     caldayloc = calday
  else
     caldayloc = calday + yr*daysperyear
  end if
  !
  ! Get and check dimension info
  !
  ierr = pio_inq_dimid( ncid_sst, 'lon', londimid   )
  ierr = pio_inq_dimid( ncid_sst, 'time', timeid  )
  ierr = pio_inq_dimid( ncid_sst, 'lat', latdimid   )
  ierr = pio_inq_dimlen( ncid_sst, londimid, lonsiz   )

  if (.not. single_column .and. lonsiz /= plon ) then
     write(iulog,*)'SSTINI: lonsiz=',lonsiz,' must = plon=',plon
     call endrun
  end if
  ierr = pio_inq_dimlen( ncid_sst, latdimid, latsiz   )
  if (.not. single_column .and. latsiz /= plat ) then
     write(iulog,*)'SSTINI: latsiz=',latsiz,' must = plat=',plat
     call endrun	
  endif
  ierr = pio_inq_dimlen( ncid_sst, timeid, timesiz   )

  if (.not. single_column) then
     !
     ! Check to ensure reduced or not grid of dataset matches that of model
     !
     if (fullgrid) then
        call pio_seterrorhandling(ncid_sst, pio_bcast_error)
        ierr = pio_inq_varid (ncid_sst, 'nlon', nlonid)
        call pio_seterrorhandling(ncid_sst, pio_internal_error)

        if (ierr == PIO_NOERR) then
           ierr = pio_get_var (ncid_sst, nlonid, nlon_sst)
           do j=1,plat
              if (nlon_sst(j) /= plon) then
                 call endrun ('SSTINI: model grid does not match dataset grid')
              end if
           end do
        end if
     else
        ierr = pio_inq_varid (ncid_sst, 'nlon', nlonid)
        ierr = pio_get_var(ncid_sst, nlonid, nlon_sst)
        do j=1,plat
           if (nlon_sst(j) /= nlon(j)) then
              call endrun ('SSTINI: model grid does not match dataset grid')
           end if
        end do
     end if
  endif
     
  ierr = pio_inq_varid( ncid_sst, 'date', dateid   )
  ierr = pio_inq_varid( ncid_sst, 'datesec', secid   )
  ierr = pio_inq_varid( ncid_sst, 'lat', latid   )
     !
     ! Retrieve entire date and sec variables.
     !

  allocate(date_sst(timesiz), sec_sst(timesiz))
  ierr = pio_get_var (ncid_sst,dateid,date_sst)
  ierr = pio_get_var (ncid_sst,secid,sec_sst)
  if (sstcyc) then
     if (timesiz<12) then 
        write(iulog,*)'SSTINI: ERROR' 
        write(iulog,*)'When cycling sst, sst data set must have 12' 
        write(iulog,*)'consecutive months of data starting with Jan'
        write(iulog,*)'Current dataset has only ',timesiz,' months'
        call endrun
     end if
     do n = 1,12
        if (mod(date_sst(n),10000)/100/=n) then
           write(iulog,*)'SSTINI: ERROR' 
           write(iulog,*)'When cycling sst, sst data set must have 12' 
           write(iulog,*)'consecutive months of data starting with Jan'
           write(iulog,*)'Month ',n,' of sst data set is out of order'
           call endrun
        end if
     end do
  end if


  if (single_column) then
     call shr_scam_GetCloseLatLon(ncid_sst%fh,scmlat,scmlon,closelat,closelon,closelatidx,closelonidx)
     strt3(1) = closelonidx
     strt3(2) = closelatidx
     strt3(3) = 1
     cnt3(1)  = 1
     cnt3(2)  = 1
     cnt3(3)  = 1
  else
     strt3(1) = 1
     strt3(2) = 1
     strt3(3) = 1
     cnt3(1)  = lonsiz
     cnt3(2)  = latsiz
     cnt3(3)  = 1
  endif
 
     !
     ! Special code for interpolation between December and January
     !
  if (sstcyc) then
     n = 12
     np1 = 1
     call bnddyi(date_sst(n  ), sec_sst(n  ), cdaysstm)
     call bnddyi(date_sst(np1), sec_sst(np1), cdaysstp)
     
     if (caldayloc<=cdaysstp .or. caldayloc>cdaysstm) then
        call infld(fieldname, ncid_sst, 'lon', 'lat', 1, pcols, begchunk,&
             endchunk, sstbdy(:,:,nm) ,readvar, grid_map='PHYS', timelevel=n)
        call infld(fieldname, ncid_sst, 'lon', 'lat', 1, pcols, begchunk,&
             endchunk, sstbdy(:,:,np) ,readvar, grid_map='PHYS', timelevel=np1)
        goto 10
     end if
  end if
  !
  ! Normal interpolation between consecutive time slices.
  !
  do n=1,timesiz-1
     np1 = n + 1
     call bnddyi(date_sst(n  ), sec_sst(n  ), cdaysstm)
     call bnddyi(date_sst(np1), sec_sst(np1), cdaysstp)
     if (.not.sstcyc) then
        yr = date_sst(n)/10000
        cdaysstm = cdaysstm + yr*daysperyear
        yr = date_sst(np1)/10000
        cdaysstp = cdaysstp + yr*daysperyear
     end if
     if (caldayloc>cdaysstm .and. caldayloc<=cdaysstp) then
        call infld(fieldname, ncid_sst, 'lon', 'lat', 1, pcols, begchunk,&
             endchunk, sstbdy(:,:,nm) ,readvar, grid_map='PHYS', timelevel=n)
        call infld(fieldname, ncid_sst, 'lon', 'lat', 1, pcols, begchunk,&
             endchunk, sstbdy(:,:,np) ,readvar, grid_map='PHYS', timelevel=np1)
        goto 10
     end if
  end do
  write(iulog,*)'SSTINI: Failed to find dates bracketing ncdate, ncsec=', ncdate, ncsec
  call endrun
10 continue
  write(iulog,*)'SSTINI: Read sst data for dates ',date_sst(n),sec_sst(n), &
       ' and ',date_sst(np1),sec_sst(np1)
  
  return
end subroutine sstini

!----------------------------------------------------------------------- 
! 
! BOP
!
! !IROUTINE: sstint
!
! !DESCRIPTION:
!
! if "aqua_planet", specify SST's analytically (Jerry Olson).
! Otherwise, time interpolate SST's to current time, reading in new monthly data if
! necessary.
! Method: 
! 
! Author: L.Bath
! 
!-----------------------------------------------------------------------
!
! !INTERFACE:
!
subroutine sstint(ncid_sst, prev_timestep)
  use ncdio_atm, only : infld
  !
  ! !INPUT PARAMETERS:
  !
  type(file_desc_t), intent(inout) :: ncid_sst
  logical, intent(in), optional :: prev_timestep ! If using previous timestep, set to true
  !
  ! EOP
  !
  !---------------------------Local variables-----------------------------
  integer dtime          ! timestep size [seconds]
  integer cnt3(3)        ! array of counts for each dimension
  integer strt3(3)       ! array of starting indices
  integer i,j,lchnk      ! indices
  integer ncol           ! number of columns in current chunk
  integer ntmp           ! temporary
  real(r8) fact1, fact2  ! time interpolation factors
  integer :: yr, mon, day! components of a date
  integer :: ncdate      ! current date in integer format [yyyymmdd]
  integer :: ncsec       ! current time of day [seconds]
  real(r8) :: calday     ! current calendar day
  real(r8) caldayloc     ! calendar day (includes yr if no cycling)
  real(r8) deltat        ! time (days) between interpolating sst data
  real(r8) xvar(plon,plat,2)    ! work space 
  integer :: ierr
  logical :: readvar
  logical :: previous              ! If using previous timestep, set to true
  !
  !-----------------------------------------------------------------------
  !
  if (aqua_planet) return

  !
  ! Use year information only if a multiyear dataset
  !
  if ( .not. present(prev_timestep) ) then
     previous = .false.
  else
     previous = prev_timestep
  end if
  if (previous .and. is_first_step()) then
     dtime = get_step_size()
     dtime = -dtime
     calday = get_curr_calday(offset=dtime)
  else
     calday = get_curr_calday()
  endif
  if ( is_perpetual() ) then
     call get_perp_date(yr, mon, day, ncsec)
  else
     if (previous .and. is_first_step()) then
        call get_curr_date(yr, mon, day, ncsec,offset=dtime)
     else
        call get_curr_date(yr, mon, day, ncsec)
     endif
  end if
  ncdate = yr*10000 + mon*100 + day
  if (sstcyc) then
     calday = get_curr_calday() 
     caldayloc = calday
  else
     caldayloc = calday + yr*daysperyear
  end if

  if (masterproc) then
     if (single_column) then
        call shr_scam_GetCloseLatLon(ncid_sst%fh,scmlat,scmlon,closelat,closelon,closelatidx,closelonidx)
        strt3(1) = closelonidx
        strt3(2) = closelatidx
        strt3(3) = 1
        cnt3(1)  = 1
        cnt3(2)  = 1
        cnt3(3)  = 1
     else
        strt3(1) = 1
        strt3(2) = 1
        strt3(3) = 1
        cnt3(1)  = lonsiz
        cnt3(2)  = latsiz
        cnt3(3)  = 1
     endif
  endif
  !
  ! If model time is past current forward sst timeslice, read in the next
  ! timeslice for time interpolation.  Messy logic is for sstcyc = .true. 
  ! interpolation between December and January (np1==1).  Note that 
  ! np1 is never 1 when sstcyc is .false.
  !
  if (caldayloc > cdaysstp .and. .not. (np1==1 .and. caldayloc>cdaysstm)) then
     if (sstcyc) then
        np1 = mod(np1,12) + 1
     else
        np1 = np1 + 1
     end if
     if (np1 > timesiz) then
        call endrun ('SSTINT: Attempt to read past end of SST dataset')
     end if
     cdaysstm = cdaysstp
     call bnddyi(date_sst(np1), sec_sst(np1), cdaysstp)

     if (.not.sstcyc) then
        yr = date_sst(np1)/10000
        cdaysstp = cdaysstp + yr*daysperyear
     end if

     if (.not. (np1 == 1 .or. caldayloc <= cdaysstp)) then
        if (masterproc) then
           write(iulog,*)'SSTINT: Input sst for date', date_sst(np1), ' sec ', sec_sst(np1), &
                ' does not exceed model date', ncdate, ' sec ', ncsec, ' Stopping.'
        end if
        call endrun ()
     end if

     ntmp = nm
     nm = np
     np = ntmp

     strt3(3) = np1

     call infld(fieldname, ncid_sst, 'lon', 'lat', 1, pcols, begchunk,&
          endchunk, sstbdy(:,:,np) ,readvar, grid_map='PHYS', timelevel=np1)
     if (masterproc) write(iulog,*)'SSTINT: Read sst for date (yyyymmdd) ',date_sst(np1), ' sec ',sec_sst(np1)
  end if
  !
  ! Determine time interpolation factors.
  !
  call get_timeinterp_factors (sstcyc, np1, cdaysstm, cdaysstp, caldayloc, &
       fact1, fact2, 'SSTINT:')

  do lchnk=begchunk,endchunk
     ncol = get_ncols_p(lchnk)
     do i=1,ncol
        sst(i,lchnk) = sstbdy(i,lchnk,nm)*fact1 + sstbdy(i,lchnk,np)*fact2
     end do
  end do

end subroutine sstint

subroutine prescribed_sst
  implicit none
  integer,parameter :: sst_option = 1

  real(r8), parameter :: pio180     = pi/180._r8
  !
  ! Parameters for zonally symmetric experiments
  !

  real(r8), parameter ::   t0_max     = 27._r8
  real(r8), parameter ::   t0_min     = 0._r8
  real(r8), parameter ::   maxlat     = 60._r8*pio180
  real(r8), parameter ::   shift      = 5._r8*pio180
  real(r8), parameter ::   shift9     = 10._r8*pio180
  real(r8), parameter ::   shift10    = 15._r8*pio180
  !
  ! Parameters for zonally asymmetric experiments
  !
  real(r8), parameter ::   t0_max6    = 1._r8
  real(r8), parameter ::   t0_max7    = 3._r8
  real(r8), parameter ::   latcen     = 0._r8*pio180
  real(r8), parameter ::   loncen     = 0._r8*pio180
  real(r8), parameter ::   latrad6    = 15._r8*pio180
  real(r8), parameter ::   latrad8    = 30._r8*pio180
  real(r8), parameter ::   lonrad     = 30._r8*pio180

  integer :: lchnk, i, ncols
  real(r8) :: tmp, tmp1, rlat(pcols), rlon(pcols)
  !
  ! Control
  !
  if(sst_option .lt. 1 .or. sst_option .gt. 10) then
     call endrun ('SSTINT:  sst_option must be between 1 and 10')
  endif
  if(sst_option == 1 .or. sst_option == 6 .or. &
       sst_option == 7 .or. sst_option == 8     ) then

     do lchnk=begchunk,endchunk
        call get_rlat_all_p(lchnk,pcols,rlat)
        ncols = get_ncols_p(lchnk)
        do i=1,ncols
           if(abs(rlat(i)) > maxlat) then
              sst(i,lchnk)=t0_min
           else
              tmp = sin(rlat(i)*pi*0.5_r8/maxlat)
              tmp = 1._r8 - tmp*tmp
              sst(i,lchnk)=tmp*(t0_max - t0_min) + t0_min
           end if
        end do
     end do
  end if
  !
  ! Flat
  !
  if(sst_option == 2) then
     do lchnk=begchunk,endchunk
        call get_rlat_all_p(lchnk,pcols,rlat)
        ncols = get_ncols_p(lchnk)
        do i=1,ncols
           if(abs(rlat(i)) > maxlat) then
              sst(i,lchnk)=t0_min
           else
              tmp = sin(rlat(i)*pi*0.5_r8/maxlat)
              tmp = 1._r8 - tmp*tmp*tmp*tmp
              sst(i,lchnk)=tmp*(t0_max - t0_min) + t0_min
           end if
        end do
     end do
  end if
  !
  ! Qobs
  !
  if(sst_option == 3) then
     do lchnk=begchunk,endchunk
        call get_rlat_all_p(lchnk,pcols,rlat)
        ncols = get_ncols_p(lchnk)
        do i=1,ncols
           if(abs(rlat(i)) > maxlat) then
              sst(i,lchnk)=t0_min
           else
              tmp = sin(rlat(i)*pi*0.5_r8/maxlat)
              tmp = (2._r8 - tmp*tmp*tmp*tmp - tmp*tmp)*0.5_r8
              sst(i,lchnk)=tmp*(t0_max - t0_min) + t0_min
           end if
        end do
     end do
  end if
  !
  ! Peaked
  !
  if(sst_option == 4) then
     do lchnk=begchunk,endchunk
        call get_rlat_all_p(lchnk,pcols,rlat)
        ncols = get_ncols_p(lchnk)
        do i=1,ncols
           if(abs(rlat(i)) > maxlat) then
              sst(i,lchnk)=t0_min
           else
              tmp = (maxlat - abs(rlat(i)))/maxlat
              tmp1 = 1._r8 - tmp
              sst(i,lchnk)= t0_max*tmp + t0_min*tmp1
           end if
        end do
     end do
  end if
  !
  ! Control-5N
  !
  if(sst_option == 5) then
     do lchnk=begchunk,endchunk
        call get_rlat_all_p(lchnk,pcols,rlat)
        ncols = get_ncols_p(lchnk)
        do i=1,ncols
           if(abs(rlat(i)) > maxlat) then
              sst(i,lchnk)=t0_min
           else if(rlat(i) > shift) then
              tmp = sin((rlat(i)-shift)*pi*0.5_r8/(maxlat-shift))
              tmp = 1._r8 - tmp*tmp
              sst(i,lchnk) = tmp*(t0_max - t0_min) + t0_min
           else
              tmp = sin((rlat(i)-shift)*pi*0.5_r8/(maxlat+shift))	
              tmp = 1._r8 - tmp*tmp
              sst(i,lchnk) = tmp*(t0_max - t0_min) + t0_min
           end if
        end do
     end do
  end if
  !
  ! 1KEQ
  !
  if(sst_option == 6) then
     do lchnk=begchunk,endchunk
        call get_rlat_all_p(lchnk,pcols,rlat)
        call get_rlon_all_p(lchnk,pcols,rlon)
        ncols = get_ncols_p(lchnk)
        do i=1,ncols
           if(abs(rlat(i)-latcen) <= latrad6) then
              tmp1 = cos((rlat(i)-latcen)*pi*0.5_r8/latrad6)
              tmp1 = tmp1*tmp1
              tmp = abs(rlon(i)-loncen)
              tmp = min(tmp , 2._r8*pi-tmp)
              if(tmp <= lonrad) then
                 tmp = cos(tmp*pi*0.5_r8/lonrad)
                 tmp = tmp*tmp
                 sst(i,lchnk) = sst(i,lchnk) + t0_max6*tmp*tmp1
              end if
           end if
        end do
     end do
  end if
  !
  ! 3KEQ
  !
  if(sst_option == 7) then
     do lchnk=begchunk,endchunk
        call get_rlat_all_p(lchnk,pcols,rlat)
        call get_rlon_all_p(lchnk,pcols,rlon)
        ncols = get_ncols_p(lchnk)
        do i=1,ncols
           if(abs(rlat(i)-latcen) <= latrad6) then
              tmp1 = cos((rlat(i)-latcen)*pi*0.5_r8/latrad6)
              tmp1 = tmp1*tmp1
              tmp = abs(rlon(i)-loncen)
              tmp = min(tmp , 2._r8*pi-tmp)
              if(tmp <= lonrad) then
                 tmp = cos(tmp*pi*0.5_r8/lonrad)
                 tmp = tmp*tmp
                 sst(i,lchnk) = sst(i,lchnk) + t0_max7*tmp*tmp1
              end if
           end if
        end do
     end do
  end if
  !
  ! 3KW1
  !
  if(sst_option == 8) then
     do lchnk=begchunk,endchunk
        call get_rlat_all_p(lchnk,pcols,rlat)
        call get_rlon_all_p(lchnk,pcols,rlon)
        ncols = get_ncols_p(lchnk)
        do i=1,ncols
           if(abs(rlat(i)-latcen) <= latrad8) then
              tmp1 = cos((rlat(i)-latcen)*pi*0.5_r8/latrad8)
              tmp1 = tmp1*tmp1
              tmp = cos(rlon(i)-loncen)
              sst(i,lchnk) = sst(i,lchnk) + t0_max7*tmp*tmp1
           end if
        end do
     end do
  end if
  !
  ! Control-10N
  !
  if(sst_option == 9) then
     do lchnk=begchunk,endchunk
        call get_rlat_all_p(lchnk,pcols,rlat)
        ncols = get_ncols_p(lchnk)
        do i=1,ncols
           if(abs(rlat(i)) > maxlat) then
              sst(i,lchnk) = t0_min
           elseif(rlat(i) > shift9) then
              tmp = sin((rlat(i)-shift9)*pi*0.5_r8/(maxlat-shift9))
              tmp = 1._r8 - tmp*tmp
              sst(i,lchnk) = tmp*(t0_max - t0_min) + t0_min
           else
              tmp = sin((rlat(i)-shift9)*pi*0.5_r8/(maxlat+shift9))
              tmp = 1._r8 - tmp*tmp
              sst(i,lchnk) = tmp*(t0_max - t0_min) + t0_min
           end if
        end do
     end do
  end if
  !
  ! Control-15N
  !
  if(sst_option == 10) then
     do lchnk=begchunk,endchunk
        call get_rlat_all_p(lchnk,pcols,rlat)
        ncols = get_ncols_p(lchnk)
        do i=1,ncols
           if(abs(rlat(i)) > maxlat) then
              sst(i,lchnk) = t0_min
           elseif(rlat(i) > shift10) then
              tmp = sin((rlat(i)-shift10)*pi*0.5_r8/(maxlat-shift10))
              tmp = 1._r8 - tmp*tmp
              sst(i,lchnk) = tmp*(t0_max - t0_min) + t0_min
           else
              tmp = sin((rlat(i)-shift10)*pi*0.5_r8/(maxlat+shift10))
              tmp = 1._r8 - tmp*tmp
              sst(i,lchnk) = tmp*(t0_max - t0_min) + t0_min
           end if
        end do
     end do
  endif
end subroutine prescribed_sst

end module sst_data

