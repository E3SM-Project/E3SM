module getinterpnetcdfdata

! Description:
!   Routines for extracting a column from a netcdf file
!
! Author: 
!   
! Modules Used:
!
  use abortutils,    only: endrun
  use pmgrid,        only: plev
  use scamMod,       only: scm_crm_mode
  use cam_logfile,   only: iulog
  implicit none
  private
!
! Public Methods:
!
  public getinterpncdata

contains

subroutine getinterpncdata( NCID, camlat, camlon, TimeIdx, &
   varName, have_surfdat, surfdat, fill_ends, &
   press, npress, ps, outData, STATUS )

!     getinterpncdata: extracts the entire level dimension for a 
!     particular lat,lon,time from a netCDF file
!     and interpolates it onto the input pressure levels, placing
!     result in outData, and the error status inx STATUS

   use shr_kind_mod, only: r8 => shr_kind_r8, i8 => shr_kind_i8
   use shr_scam_mod, only: shr_scam_GetCloseLatLon
   use netcdf
   implicit none
!-----------------------------------------------------------------------


!     ---------- inputs ------------

   integer, intent(in)  :: NCID          ! NetCDF ID
   integer, intent(in)  :: TimeIdx       ! time index
   real(r8), intent(in) :: camlat,camlon ! target lat and lon to be extracted  
   logical, intent(in)  :: have_surfdat  ! is surfdat provided
   logical, intent(in)  :: fill_ends ! extrapolate the end values
   integer, intent(in)  :: npress        ! number of dataset pressure levels
   real(r8), intent(in) :: press(npress) ! dataset pressure levels
   real(r8), intent(in) :: ps ! dataset pressure levels

!     ---------- outputs ----------

   real(r8), intent(inout)  :: outData(:)    ! interpolated output data
   integer, intent(out)   :: STATUS        ! return status of netcdf calls

!     -------  locals ---------

   real(r8)  surfdat       ! surface value to be added before interpolation
   integer nlev          ! number of levels in dataset
   integer     latIdx        ! latitude index
   integer     lonIdx        ! longitude index
   real(r8),allocatable :: tmp(:)
   real(r8)        closelat,closelon
   real(r8)        dx, dy, m             ! slope for interpolation of surfdat
   integer     varID
   integer     var_ndims
   integer     dims_set
   integer     i
   integer     var_dimIDs( NF90_MAX_VAR_DIMS )
   integer     start( NF90_MAX_VAR_DIMS ) 
   integer     count( NF90_MAX_VAR_DIMS )

   character   varName*(*)
   character   dim_name*( 256 )
   real(r8)        missing_val
   logical     usable_var

!     -------  code ---------

   call shr_scam_GetCloseLatLon(ncid,camlat,camlon,closelat,closelon,latidx,lonidx)

!
! Check mode: double or single precision
!

!
! Get var ID.  Check to make sure it's there.
!
   STATUS = NF90_INQ_VARID( NCID, varName, varID )

   if ( STATUS .NE. NF90_NOERR )  return

!
! Check the var variable's information with what we are expecting
! it to be.
!

   STATUS = NF90_INQUIRE_VARIABLE( NCID, varID, ndims=var_ndims )
   if ( var_ndims .GT. 4 ) then
      write(iulog,* ) 'ERROR - extractdata.F: The input var',varName, &
         'has', var_ndims, 'dimensions'
      STATUS = -1
   endif

!
!     surface variables
!
   if ( var_ndims .EQ. 0 ) then
      STATUS = NF90_GET_VAR( NCID, varID, outData )
      return
   endif

   STATUS = NF90_INQUIRE_VARIABLE( NCID, varID, dimids=var_dimIDs )
   if ( STATUS .NE. NF90_NOERR ) then
      write(iulog,* ) 'ERROR - extractdata.F:Cant get dimension IDs for', varName
      return
   endif
!     
!     Initialize the start and count arrays 
!     
   dims_set = 0
   nlev = 1
   do i =  var_ndims, 1, -1

      usable_var = .false.
      STATUS = NF90_INQUIRE_DIMENSION( NCID, var_dimIDs( i ), dim_name )

      if ( dim_name .EQ. 'lat' ) then
         start( i ) =  latIdx
         count( i ) = 1           ! Extract a single value 
         dims_set = dims_set + 1
         usable_var = .true.
      endif

      if ( dim_name .EQ. 'lon' ) then
         start( i ) = lonIdx
         count( i ) = 1           ! Extract a single value
         dims_set = dims_set + 1
         usable_var = .true.
      endif

      if ( dim_name .EQ. 'lev' ) then
         STATUS = NF90_INQUIRE_DIMENSION( NCID, var_dimIDs( i ), len=nlev )
         start( i ) = 1
         count( i ) = nlev       ! Extract all levels
         dims_set = dims_set + 1
         usable_var = .true.
      endif

      if ( dim_name .EQ. 'ilev' ) then
         STATUS = NF90_INQUIRE_DIMENSION( NCID, var_dimIDs( i ), len=nlev )
         start( i ) = 1
         count( i ) = nlev        ! Extract all levels
         dims_set = dims_set + 1
         usable_var = .true.
      endif

      if ( dim_name .EQ. 'time' .OR. dim_name .EQ. 'tsec' ) then 
         start( i ) = TimeIdx
         count( i ) = 1           ! Extract a single value 
         dims_set = dims_set + 1   
         usable_var = .true.
      endif

      if ( usable_var .EQV. .false. ) then
         write(iulog,* )'ERROR - extractdata.F: The input var ',varName, &
            ' has an unusable dimension ', dim_name
         STATUS = 1
      endif
   end do

   if ( dims_set .NE. var_ndims ) then
      write(iulog,* )'ERROR - extractdata.F: Could not find all the', &
         ' dimensions for input var ', varName
      write(iulog,* )'Found ',dims_set, ' of ',var_ndims
      STATUS = 1
   endif

   allocate(tmp(nlev+1))

   STATUS = NF90_GET_VAR( NCID, varID, tmp, start, count )

   if ( STATUS .NE. NF90_NOERR ) then
      write(iulog,* )'ERROR - extractdata.F: Could not get data for input var ', varName
      return
   endif

   if ( nlev .eq. 1 ) then
      outdata(1) = tmp(1)
      return                 ! no need to do interpolation 
   endif
!   if ( use_camiop .and. nlev.eq.plev) then
   if ( nlev.eq.plev .or. nlev.eq.plev+1) then
      outData(:nlev)= tmp(:nlev)! no need to do interpolation 
   else
!
!     add the surface data if available, else
!     fill in missing surface data by extrapolation
!
      if(.not.scm_crm_mode) then
         if ( have_surfdat ) then
            tmp(npress) = surfdat
         else
            dy = press(npress-1) - press(npress-2)
            dx = tmp(npress-1) - tmp(npress-2)
            if ( dx .ne. 0.0_r8 ) then
               m = dy/dx
               tmp(npress) = ((press(npress) - press(npress-1)) / m ) + tmp(npress-1)
            else
               tmp(npress) = tmp(npress-1)
            endif
            surfdat = tmp(npress)
         endif
      endif

#if DEBUG > 1
!
!     check data for missing values
!

      STATUS = NF90_GET_ATT( NCID, varID, 'missing_value', missing_val )
      if ( STATUS .NE. NF90_NOERR ) then
         missing_val = -9999999.0_r8
      endif
!
! reset status to zero
!     
      STATUS = 0
!
      do i=1, npress
         if ( tmp(i) .eq. missing_val ) then
            write(iulog,*) 'ERROR - missing value found in ', varname
            write(iulog,*) 'time,lat,lon,lev = ' ,timeidx, latidx, lonidx, i
            stop
         endif
      enddo
#endif
!
      call interplevs( tmp(:npress), press, npress, ps, fill_ends,outdata )

   endif

   deallocate(tmp)
   return
 end subroutine getinterpncdata

subroutine interplevs( inputdata,   dplevs,   nlev, &
                       ps, fill_ends,   outdata)

   use shr_kind_mod, only: r8 => shr_kind_r8, i8 => shr_kind_i8
   use hycoef, only: hyam, hybm
   use interpolate_data, only: lininterp
   implicit none

!
!     WARNING: ps, siga and sigb must be initialized before calling this routine
!

!------------------------------Commons----------------------------------

!-----------------------------------------------------------------------


!     ------- inputs -----------
   integer, intent(in) :: nlev                 ! num press levels in dataset

   real(r8), intent(in) :: ps     ! surface pressure
   real(r8), intent(in) :: inputdata(nlev)     ! data from netcdf dataset
   real(r8), intent(in) :: dplevs(nlev)         ! input data pressure levels 

   logical, intent(in) :: fill_ends            ! fill in missing end values(used for
                                            ! global model datasets)


! ------- outputs ----------
   real(r8), intent(inout) :: outdata(:)      ! interpolated column data

! ------- locals -----------

   real(r8) mplevs( PLEV )
   real(r8) interpdata( PLEV )


   integer dstart_lev, dend_lev 
   integer mstart_lev, mend_lev
   integer data_nlevs, model_nlevs, i
   integer STATUS

!
!     Initialize  model_pressure_levels.  ps should be set in the calling
!     routine to the value in the dataset
!
   do i = 1, plev
      mplevs( i ) = 1000.0_r8 * hyam( i ) + ps * hybm( i ) / 100.0_r8
   end do
!     
!     the following algorithm assumes that pressures are increasing in the
!     arrays
!     
!     
!     Find the data pressure levels that are just outside the range
!     of the model pressure levels, and that contain valid values
!     
   dstart_lev = 1
   do i= 1, nlev
      if ( dplevs(i) .LE. mplevs(1) ) dstart_lev  = i
   end do

   dend_lev = nlev
   do i= nlev, 1, -1
      if ( dplevs(i) .GE. mplevs(plev) ) then
         dend_lev  = i
      endif
   end do
!         
!     Find the model pressure levels that are just inside the range
!     of the data pressure levels
!
   mstart_lev = 1
   do i=plev, 1, -1
      if ( mplevs( i ) .GE. dplevs( dstart_lev ) )  mstart_lev = i
   end do

   mend_lev = plev
   do i=1,plev
      if ( mplevs( i ) .LE. dplevs( dend_lev ) ) mend_lev = i
   end do

   data_nlevs = dend_lev - dstart_lev +1
   model_nlevs = mend_lev - mstart_lev +1

   call lininterp (inputdata(dstart_lev:dend_lev),dplevs(dstart_lev:dend_lev),data_nlevs, &
		   interpdata,mplevs(mstart_lev:mend_lev),model_nlevs)
!
!     interpolate data onto the model pressure levels
!
!!$   call lininterp (inputdata,dplevs,nlev, &
!!$		   outdata(:plev),mplevs,plev)
   do i=1 , model_nlevs
      outdata( i+mstart_lev-1 ) = interpdata( i )
   end do
!
!     fill in the missing end values 
!           (usually  done if this is global model dataset)
!
   if ( fill_ends ) then 
      do i=1, mstart_lev
         outdata(i) = inputdata(1)
      end do
      do i= mend_lev, plev
         outdata(i) = inputdata(nlev)
      end do
   end if

   return
end subroutine interplevs
end module getinterpnetcdfdata

