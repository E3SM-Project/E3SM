subroutine driver (ncprec, rgridnl,isncol,nxi,nyi,nxo,nyo,nz,ntime,nxof,nyof)
!------------------------------------------------------------------
! Purpose:
!   Horizontally interpolate file containing aerosol masses to CAM grid.
!   Sum the values in the vertical such that output arrays contain the
!   column sum from the bottom of the atmosphere to that level.
!   Convert monthly averages to mid month values.
!   Input file assumed to be MATCH ncep runs, averaged by month.
!        and backsolved to provide Mid-month values)
!     
!  Method:
!    read data from file
!    interpolate data onto CAM horizontal grid
!  Modified: by Jim Edwards  3/2006
!    CAM horizontal grid may either be a rectangular lat/lon grid, or a 1d list of columns
!    in the case of a 1d list the output arrays are of size f(ncol,nz,1) and lonout(ncol) 
!    and latout(ncol) contain a lon/lat pair for each point of the output grid.
!    to achieve this a few variables in this file are overloaded.  
!    Specifically in the lat/lon case we have nxof=nxo and nyof=nyo and in the ncols case:
!    nxo=nyo=ncols
!    nxof=ncols
!    nyof=1
!------------------------------------------------------------------
   use shr_kind_mod, only : r8=>shr_kind_r8
   use globals
   use preserve_mean, only: monthly_to_midmonth
   use interpolate_data, only: bilin, lininterp_init, lininterp, lininterp_finish, interp_type
   use netcdf
   use error_messages, only : handle_ncerr  

   implicit none


!
! Arguments
!
   integer, intent(in) :: ncprec     ! specify 32-bit or 64-bit precision for output variables
   character(len=*), intent(in) :: rgridnl ! reduced grid namelist (if applicable)
   logical, intent(in) :: isncol ! full grid list of lat/lon pairs instead of rectangular lat by lon grid
   integer, intent(in) :: nxi ! first dimension of input array (longitude)
   integer, intent(in) :: nz  ! second dimension of input and output array (level)
   integer, intent(in) :: nyi ! third dimension of input array (latitude)
   integer, intent(in) :: nxo ! size of output longitude or ncol array
   integer, intent(in) :: nyo ! size of output latitude or ncol array
   integer, intent(in) :: ntime ! number of time levels to interpolate
   integer, intent(in) :: nxof ! first dimension of output array (longitude or ncol)
   integer, intent(in) :: nyof ! third dimension of output array (latitude or 1)
!
! Local workspace
!
   character(len=8), parameter :: aerosol_name(naer) =  &
        (/"MSUL    "&
         ,"MSSLT   "&
         ,"MDUST1  "&
         ,"MDUST2  "&
         ,"MDUST3  "&
         ,"MDUST4  "&
         ,"MOCPHO  "&
         ,"MBCPHO  "&
         ,"MOCPHI  "&
         ,"MBCPHI  "/)
   
   integer, parameter :: idxsslt = 2   ! index of SSLT
   
   integer :: lonidi = -1              ! longitude id input file
   integer :: latidi = -1              ! latitude id input file
   integer :: dateidi = -1             ! date id input file
   integer :: datesecidi = -1          ! datesec id input file
   integer :: mhybiid = -1             ! MATCH hybi id input file
   integer :: mpsid = -1               ! MATCH PS id input file
   integer :: species_id(naer)         ! aerosol ids input file

   integer :: lonido = -1              ! longitude id output file
   integer :: latido = -1              ! latitude id output file
   integer :: dateido = -1             ! date id output file
   integer :: datesecido = -1          ! datesec id output file
   integer :: hybiido = -1             ! MATCH hybi id output file
   integer :: psido = -1               ! MATCH PS id output file
   integer :: varido(naer)             ! aerosol ids output file
   integer :: timeido                  ! time id output filep

   integer :: nlonid                   ! nlon id output file (if present)
   integer :: rlonid                   ! rlon id output file (if present)

   real(r8) :: tempin(nxi,nz,nyi)      ! temp variable for bilinear interpolation
   real(r8) :: lonin(nxi)              ! longitude on input grid (rectangular)
   real(r8) :: latin(nyi)              ! latitude on input grid
   real(r8) :: m_hybi(nz+1)            ! MATCH hybi
   real(r8) :: m_ps(nxi,nyi,ntime)     ! surface pressure from MATCH

   real(r8) :: tempout(nxof,nz,nyof)     ! temp variable for bilinear interpolation
   real(r8) :: rlonout(nxof,nyof)        ! longitude on output grid (2-d for reduced grid)
   real(r8) :: lonout(nxo)             ! longitude on output grid (if rectangular)
   real(r8) :: latout(nyo)             ! latitude on input grid
   real(r8) :: M_ps_cam(nxof,nyof,ntime) ! surface pressure from MATCH on cam grid

   integer :: istat                    ! status return
   integer :: i, j, k                  ! x,y,z indices
   integer :: m                        ! constituent index
   integer :: mo                       ! month index
   integer :: dimids(nf90_max_var_dims)  ! variable shape on output netcdf file
   integer :: date(ntime)              ! date (yyyymmdd)
   integer :: datesec(ntime)           ! seconds of date
   integer :: nlonout(nyo)             ! number of longitudes per latitude
   integer :: start(4)                 ! starting position in netcdf file
   integer :: kount(4)                 ! number of values to take from netcdf file
   integer :: dimcnt
   type(interp_type) :: lon_wgts, lat_wgts

!
! Reduced grid namelist.  Fortran insists namelist arrs not be dimensioned
! dynamically 
!
   integer, parameter :: maxlat = 10000
   integer :: nlon(maxlat) = -1
   namelist /reduced/ nlon
   real(r8), parameter :: fillvalue = 1.d36
!
! a temporary place to store mmr's from files (generated from MATCH runs)
!
   real(r8) :: fspecies(nxi,nyi,nz)          ! aerosol mmr's from MATCH file
   real(r8) :: aerosol(nxof,nyof,nz,ntime) ! aerosol mmr's from MATCH file on CAM grid
!!!   real(r8), allocatable :: aerosol(:,:,:,:) ! aerosol mmr's from MATCH file on CAM grid
                                             ! allocate instead of dimension directly to prevent stack overflow

!!!   allocate (aerosol(nxof,nyof,nz,ntime))
!
! Get required info from input file.
!
   call handle_ncerr( nf90_inq_varid (ncidi, 'lon', lonidi),&
     'driver.f90:127')
   call handle_ncerr( nf90_inq_varid (ncidi, 'lat', latidi),&
     'driver.f90:129')
   call handle_ncerr( nf90_inq_varid (ncidi, 'hybi', mhybiid),&
     'driver.f90:131')
   call handle_ncerr( nf90_inq_varid (ncidi, 'PS', mpsid),&
     'driver.f90:133')
   call handle_ncerr( nf90_inq_varid (ncidi, 'date', dateidi),&
     'driver.f90:135')
   call handle_ncerr( nf90_inq_varid (ncidi, 'datesec', datesecidi),&
     'driver.f90:137')

   call handle_ncerr( nf90_get_var (ncidi, lonidi, lonin),&
     'driver.f90:140')
   call handle_ncerr( nf90_get_var (ncidi, latidi, latin),&
     'driver.f90:142')
   call handle_ncerr( nf90_get_var (ncidi, mhybiid, m_hybi),&
     'driver.f90:144')
   call handle_ncerr( nf90_get_var (ncidi, mpsid, m_ps),&
     'driver.f90:146')
   call handle_ncerr( nf90_get_var (ncidi, dateidi, date),&
     'driver.f90:148')
   call handle_ncerr( nf90_get_var (ncidi, datesecidi, datesec),&
     'driver.f90:150')
!
! If time variable is not found on output file, create it.
! Then create date and datesec variables.
!
   if (nf90_inq_varid (ncido, 'time', timeido) /= nf90_noerr) then
      call handle_ncerr( nf90_def_var (ncido, 'time', nf90_double,  timedimido, timeido),&
        'driver.f90:157')
   end if
   call handle_ncerr( nf90_def_var (ncido, 'date', nf90_int, timedimido, dateido),&
     'driver.f90:160')
   call handle_ncerr( nf90_def_var (ncido, 'datesec', nf90_int, timedimido, datesecido),&
     'driver.f90:162')
!
! Define hybi and PS on output grid
!
   call handle_ncerr( nf90_def_var (ncido, 'hybi', nf90_double, ilevdimido, hybiido),&
     'driver.f90:167')
   if(isncol) then
      dimids(1:4) = (/ncoldimido, timedimido, -1, -1/)
      dimcnt=2
   else
      dimids(1:4) = (/londimido, latdimido, timedimido, -1/)
      dimcnt=3
   end if
   call handle_ncerr( nf90_def_var (ncido, 'PS', ncprec, dimids(1:dimcnt), psido),&
     'driver.f90:176')
   call handle_ncerr( nf90_put_att (ncido, psido, '_FillValue',  fillvalue),&
     'driver.f90:178')
!
! Read input variable names and define output names accordingly
! Append '_V' to indicate field has been vertically summed from sfc to each level
!
   dimcnt=dimcnt+1
   dimids(dimcnt-1:dimcnt+1) = (/levdimido, timedimido, -1/)
   do m = 1, naer
      call handle_ncerr( nf90_inq_varid (ncidi, trim (aerosol_name(m)), species_id(m)),&
        'driver.f90:187')
      call handle_ncerr( nf90_def_var (ncido, trim (aerosol_name(m))//'_V', ncprec, dimids(1:dimcnt), varido(m)),&
        'driver.f90:189')
      call handle_ncerr( nf90_put_att (ncido, varido(m), '_FillValue',  fillvalue),&
        'driver.f90:191')
   end do
!
! Define global attribute "cam-ready" which will be checked for by CAM to prevent
! the use of datasets not run through the interpaerosols procedure.
!
   call handle_ncerr( nf90_put_att (ncido, nf90_global, 'cam-ready', 'yes'),&
     'driver.f90:198')
!
! End define mode on output file.  Copy required data from input file to output file
!
   call handle_ncerr( nf90_enddef (ncido),&
     'driver.f90:203')
!
! Retrieve output grid definition
!
   call handle_ncerr( nf90_inq_varid (ncido, 'lat', latido),&
     'driver.f90:208')
   call handle_ncerr( nf90_get_var (ncido, latido, latout),&
     'driver.f90:210')
!
! Define nlon and rlon if output file is on reduced grid
!
   if (rgridnl /= ' ') then
      call handle_ncerr( nf90_redef (ncido),&
        'driver.f90:216')
      call handle_ncerr( nf90_def_var (ncido, 'nlon', NF90_INT, latdimido, nlonid),&
        'driver.f90:218')
      dimids(:3) = (/londimido, latdimido, -1/);
      dimcnt=2
      call handle_ncerr( nf90_def_var (ncido, 'rlon', NF90_DOUBLE,  dimids(1:dimcnt), rlonid),&
        'driver.f90:221')
      call handle_ncerr( nf90_put_att (ncido, rlonid, '_FillValue',  fillvalue),&
        'driver.f90:223')
      call handle_ncerr( nf90_enddef (ncido),&
        'driver.f90:225')
      open (unit=7, form='formatted', status='old', file=rgridnl, iostat=istat)
      if (istat /= 0) then
         write(6,*)'Namelist file ', rgridnl, ' cannot be opened for reading'
         stop 999
      end if
      if (nyo > maxlat) then
         write(6,*)'Parameter maxlat must be at least nyo=', nyo
         write(6,*)'Should just need to change parameter setting in regrid.f90'
         stop 999
      end if
      read(7,reduced)
      nlonout(:nyo) = nlon(:nyo)
      do j=1,nyo
         if (nlonout(j) < 1 .or. nlonout(j) > nxo .or. &
             nlonout(j) /= nlonout(nyo-j+1)) then
            write(6,*)'Bad nlonout value=', nlonout(j)
            stop 999
         end if
      end do
      do j=1,nyo
         rlonout(:,j) = fillvalue
         do i=1,nlonout(j)
            rlonout(i,j) = (i-1)*360./nlonout(j)
         end do
      end do
   else                                                 ! full grid: define nlon and rlon
      call handle_ncerr( nf90_inq_varid (ncido, 'lon', lonido),&
        'driver.f90:253')
      call handle_ncerr( nf90_get_var (ncido, lonido, lonout),&
        'driver.f90:255')
      do j=1,nyof
         rlonout(:,j) = lonout(:)
         nlonout(j) = nxo
      end do
   end if
!
! interpolate match's surface pressure and get mid-month values
!
   if(isncol) then
      ! 2 cyclic boundaries, 1 means set to boundary value.
      call lininterp_init(lonin, nxi, lonout, nxo, 2, lon_wgts)
      call lininterp_init(latin, nyi, latout, nyo, 1, lat_wgts)
      do mo=1,ntime
         call lininterp(M_ps(:,:,mo),nxi,nyi,M_ps_cam(:,1,mo),nxo,lon_wgts,lat_wgts)
 !        print *, 'mean ps month: ',mo,' input: ',sum(M_ps(:,:,mo))/real(nxi*nyi,kind=r8), &
 !             ' output: ', sum(M_ps_cam(:,1,mo))/real(nxo,kind=r8)
 !        call simplemean(nxi,nyi,nyi,M_ps(:,:,mo),latin)
 !        call simplemean(nxo,1,nxo,M_PS_cam(:,:,mo),latout)
 !        print *,minval(M_ps(:,:,mo)), maxval(M_ps(:,:,mo)), &
 !             minval(M_ps_cam(:,1,mo)), maxval(M_ps_cam(:,1,mo))
      end do
   else
      do mo=1,ntime
         call bilin (M_ps(1,1,mo), lonin, latin, nxi, nxi, &
              1, 1, nyi, M_ps_cam(1,1,mo), rlonout, &
              latout, nxo, nlonout, 1, nyo)
         print *, 'mean ps month: ',mo,' input: ',sum(M_ps(:,:,mo))/real(nxi*nyi,kind=r8), &
              ' output: ', sum(M_ps_cam(:,:,mo))/real(nxo*nyo,kind=r8)
      end do
   end if
   call monthly_to_midmonth (M_ps_cam, 1, nlonout,nxof,nyof,ntime)
!
! Retrieve Aerosol Masses (kg/m^2 in each layer)
!
   do m=1,naer
      do mo=1,ntime
         start(:) = (/1,1,1,mo/)
         kount(:) = (/nxi,nyi,nz,1/)
         call handle_ncerr( nf90_get_var (ncidi, species_id(m), fspecies, start,kount),&
           'driver.f90:297')
!
! Accumulate mass below each interface level
!  note that lowest level (nz+1) is assumed to be zero 
!  but there isn't even storage for this value in the array
!
         do k=nz-1,1,-1
            fspecies(:,:,k) = fspecies(:,:,k) + fspecies(:,:,k+1)
         end do
!
! Transpose coords of retrieved aerosols to enable using CAM's bi-linear
! interpolation code.  Interpolate onto CAM horizontal grid
!
         if(isncol) then
            do k=1,nz
               call lininterp(fspecies(:,:,k),nxi,nyi,tempout(:,k,1),nxo,lon_wgts,lat_wgts)
            end do
         else
            do k=1,nz
               tempin(:,k,:) = fspecies(:,:,k)
            end do
            call bilin (tempin, lonin, latin, nxi, nxi, &
                 nz, nz, nyi, tempout, rlonout, &
                 latout, nxo, nlonout, nz, nyo)
         end if
         do k=1,nz
            do j=1,nyof
               do i=1,nlonout(j)
                  aerosol(i,j,k,mo) = tempout(i,k,j)
               end do
            end do
         end do
!
! Sea Salt over land is minuscule.  After interpolation, differencing
! different levels of the cumulative mass can lead to underflow errors.
! To solve this problem, set sea salt to 0 for any column where
! the total mass is less than a threashold which looks like roundoff
! (and is unmeasurable?)
!
         if (m == idxSSLT) then
            do j=1,nyof
               do i=1,nlonout(j)
                  if (aerosol(i,j,1,mo) < 1.e-24) then
                     aerosol(i,j,:,mo) = 0.
                  end if
               end do
            end do
         end if
      end do ! mo
!
! convert from monthly average to mid-month values
!
      call monthly_to_midmonth (aerosol, nz, nlonout,nxof,nyof,ntime)

      do mo=1,ntime
         do j=1,nyof
!
! make sure total column mass total is not negative
!
            do i=1,nlonout(j)
               aerosol(i,j,1,mo) = max (aerosol(i,j,1,mo), 0._r8)
            end do
!
! make function non-increasing and positive
!
            do k=2,nz
               do i=1,nlonout(j)
                  aerosol(i,j,k,mo) = min (aerosol(i,j,k,mo), aerosol(i,j,k-1,mo))
                  aerosol(i,j,k,mo) = max (aerosol(i,j,k,mo), 0._r8)
               end do
            end do
            do k=1,nz
               do i=nlonout(j)+1,nxo
                  aerosol(i,j,k,mo) = fillvalue
               end do
            end do
         end do
!
! Write interpolated data to output file
!
         if(isncol) then
            start(:) = (/1,1,mo,-1/)
            kount(:) = (/nxo,nz,1,-1/)
         else
            start(:) = (/1,1,1,mo/)
            kount(:) = (/nxo,nyo,nz,1/)
         end if
         call handle_ncerr( nf90_put_var (ncido, varido(m), aerosol(:,:,:,mo),start,kount),&
           'driver.f90:385')
      end do     ! loop over months (mo)
   end do        ! loop over constituents (m)

   call handle_ncerr( nf90_put_var (ncido, hybiido, m_hybi),&
     'driver.f90:390')

   if (rgridnl /= ' ') then
      call handle_ncerr( nf90_put_var (ncido, nlonid, nlon),&
        'driver.f90:394')
      call handle_ncerr( nf90_put_var (ncido, rlonid, rlonout),&
        'driver.f90:396')
   end if

   do mo=1,ntime
      do j=1,nyof
         do i=nlonout(j)+1,nxo
            M_ps_cam(i,j,mo) = fillvalue
         end do
      end do
      if(isncol) then
         start(:) = (/1,mo,-1,-1/)
         kount(:) = (/nxo,1,-1,-1/)
      else
         start(:) = (/1,1,mo,-1/)
         kount(:) = (/nxo,nyo,1,-1/)
      end if
      call handle_ncerr( nf90_put_var (ncido, psido, M_ps_cam(:,:,mo),start,kount),&
        'driver.f90:413')

      start(:) = (/mo,-1,-1,-1/)
      kount(:) = (/1,-1,-1,-1/)
      call handle_ncerr( nf90_put_var (ncido, dateido, date(mo:mo),start,kount),&
        'driver.f90:418')
      call handle_ncerr( nf90_put_var (ncido, datesecido, datesec(mo:mo),start,kount),&
        'driver.f90:420')
   end do

!!!   deallocate (aerosol)
   return
end subroutine driver

subroutine simplemean(n1,n2,n3,f,lat)
  use shr_kind_mod, only : r8=>shr_kind_r8
  integer, intent(in) :: n1, n2, n3
  real(r8), intent(in) :: f(n1,n2)
  real(r8), intent(in) :: lat(n3)
  real(r8) :: mean, wgt
  integer :: i, j


  if(n2==1) then
     mean = sum(f(:,1)*cos(lat))/sum(cos(lat))
  else
     mean = 0.
     wgt = 0.
     do j=1,n2
        do i=1,n1
           wgt = wgt+cos(lat(j))
           mean = mean+f(i,j)*cos(lat(j))
        end do
     end do
     mean = mean/wgt
  end if
  print *, 'weighted mean = ',mean

end subroutine simplemean
