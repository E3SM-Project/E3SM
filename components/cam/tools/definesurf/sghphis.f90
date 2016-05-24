subroutine sghphis (plon, plat, numlons, mlatcnts, mloncnts, & 
                    topofile, verbose, sgh, sgh30, have_sgh30, phis, fland )

!-----------------------------------------------------------------------
!
! Read high resolution topo dataset and calculate values of phis and sgh
! for the model resolution this model
!
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8

   implicit none
   include 'netcdf.inc'
!
!-----------------------------------------------------------------------
!
! parameters
!
   integer , parameter :: ntopolon = 2160
   integer , parameter :: ntopolat = 1080
   integer , parameter :: n2x2lon  = 180
   integer , parameter :: n2x2lat  = 90
   integer , parameter :: n3x3lon  = 120
   integer , parameter :: n3x3lat  = 60
   real(r8), parameter :: r8_360   = 360. ! For argument compatibility to mod
!
! arguments
!
   integer , intent(in) :: plon                ! maximum number of model longitudes
   integer , intent(in) :: plat                ! number of model latitudes
   integer , intent(in) :: numlons(plat)       ! number of model longitudes per latitude
   real(r8), intent(in) :: mlatcnts(plat)      ! model cell center latitudes 
   real(r8), intent(in) :: mloncnts(plon,plat) ! model cell ceneter longitudes
   logical , intent(in) :: verbose             ! true => verbose output 
   character(len=*), intent(in) :: topofile    ! high resolution topo file    
   real(r8), intent(out):: phis(plon,plat)     ! model geopotention height 
   real(r8), intent(out):: sgh(plon,plat)      ! model standard dev of geopotential height above 10min
   real(r8), intent(out):: sgh30(plon,plat)    ! model standard dev of geopotential height from 30s to 10m
   logical , intent(out):: have_sgh30          ! true => variance is on topofile, sgh30 will be output
   real(r8), intent(out):: fland(plon,plat)    ! model fractional land
!
! Local workspace : note that anything with plon or plat in its dimension is dynamic
!
   real(r8) wt                          ! weight for area averaging
   real(r8) dx,dy                       ! increments for definition of intermed grid   

! high resolution topo grid 

   integer lonid_topo, latid_topo       ! input topo file vars
   integer htopoid,ftopoid,ret,varianceid ! input topo file vars
   real(r8) tloncnts(ntopolon)          ! topo cell center lon boundaries
   real(r8) tlatcnts(ntopolat)          ! topo cell center lat boundaries
   real(r8) tlons(ntopolon+1,ntopolat)  ! topo cell W lon boundaries
   real(r8) tlats(ntopolat+1)           ! topo cell N lat boundaries
   real(r8) ftopo(ntopolon,ntopolat)    ! Land fraction array 
   real(r8) htopo(ntopolon,ntopolat)    ! Topographic heights
   real(r8) variance(ntopolon,ntopolat)    ! Variance of elev at 30sec

! intermediate grid

   real(r8) lons3x3(n3x3lon+1,n3x3lat)  ! list of topo cell W lon boundaries
   real(r8) lats3x3(n3x3lat+1)          ! list of topo cell N lat boundaries
   integer num3x3lons(n3x3lat)          ! number if longitudes per latitude
   real(r8) mnhgt3x3(n3x3lon,n3x3lat)   ! intermediate topo height 
   real(r8) varhgt3x3(n3x3lon,n3x3lat)  ! intermediate topovariance

! model grid

   real(r8) mlons(plon+1,plat)          ! model cell W lon boundaries
   real(r8) mlats(plat+1)               ! model cell N lat boundaries
   real(r8) mnhgt(plon,plat)            ! model topographic height
   real(r8) varhgt(plon,plat)           ! model topographic variance
   real(r8) summn, sumvar              ! use only for pole point calculations

! other vars

   real(r8) xmax                       ! temporary variable 
   real(r8), parameter :: eps = 1.e-6  ! eps criterion for pole point
   integer imax, jmax                  ! indices
   integer i,j,ii,ji,io,jo,n           ! indices     
   integer ncid_topo                   ! topographic netcdf id
   integer ioe
   integer mxovr ! max number of fine grid points used in area calculation of model grid point
!
! Space needed in 3 dimensions to store the initial data. This space is
! required because the input data file does not have a predetermined 
! ordering of the latitude records. A  specific order is imposed in the 
! transforms so that the results will be reproducible.
!
! Dynamic
!
   integer , allocatable :: iovr(:,:,:) ! lon index of overlap input cell
   integer , allocatable :: jovr(:,:,:) ! lat index of overlap input cell
   real(r8), allocatable :: wovr(:,:,:) ! weight of overlap input cell
!
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------------
! Read in navy topo cell locations and determine cell edges (Uniform grid)
!----------------------------------------------------------------------------
!
   ret = nf_open (topofile, nf_nowrite, ncid_topo)
   if (ret == nf_noerr) then
      if (verbose) write(6,*)'Successfully opened netcdf topofile ',trim(topofile)
      ret = nf_inq_varid (ncid_topo, 'variance', varianceid)
      if (ret == NF_NOERR) then
         if (verbose) write(6,*)'Found a new style topofile.'
         call wrap_get_var8 (ncid_topo, varianceid, variance  )
         call wrap_inq_varid (ncid_topo, 'landfract', ftopoid   )
         have_sgh30 = .true.
      else
         if (verbose) write(6,*)'Found an old style topofile.'
         call wrap_inq_varid (ncid_topo, 'ftopo', ftopoid   )
         have_sgh30 = .false.
      end if
      call wrap_get_var8 (ncid_topo, ftopoid, ftopo)
      call wrap_inq_varid (ncid_topo, 'htopo', htopoid   )
      call wrap_get_var8 (ncid_topo, htopoid, htopo)
   else
      write(6,*)'cannot open topo file successfully'
      call endrun
   endif

   call wrap_inq_varid (ncid_topo, 'lon', lonid_topo)
   call wrap_inq_varid (ncid_topo, 'lat', latid_topo)
   
   call wrap_get_var8 (ncid_topo, latid_topo, tlatcnts)
   call wrap_get_var8 (ncid_topo, lonid_topo, tloncnts)
   ret = nf_close (ncid_topo)
 
   tloncnts(:) = mod(tloncnts(:)+r8_360,r8_360) 

   tlats(:) = 1.e36
   tlats(1) = -90.                                    ! south pole
   do j = 2, ntopolat                               
      tlats(j) = (tlatcnts(j-1) + tlatcnts(j)) / 2.   ! southern edges
   end do
   tlats(ntopolat+1) = 90.                            ! north pole

   tlons(:,:) = 1.e36
   do j = 1,ntopolat
      dx = 360./ntopolon
      tlons(1,j) = tloncnts(1) - dx/2.
      do i = 2, ntopolon
         tlons(i,j) = tloncnts(i) - dx/2.
      end do
      tlons(ntopolon+1,j) = tloncnts(ntopolon) + dx/2.
   end do
!
!----------------------------------------------------------------------------
! Determine model cell edges
!----------------------------------------------------------------------------
!
   mlats(:) = 1.e36
   mlats(1) = -90.                                    ! south pole
   do j = 2,plat                               
      mlats(j) = (mlatcnts(j-1) + mlatcnts(j)) / 2.   ! southern edges
   end do
   mlats(plat+1) = 90.                                ! north pole

   do j = 1,plat
      dx = 360./(numlons(j))
      do i = 1,plon+1
         mlons(i,j) = -dx/2. + (i-1)*dx
      end do
   end do

!
!----------------------------------------------------------------------------
! Calculate fractional land
!----------------------------------------------------------------------------
!
   call binf2c(tloncnts ,tlatcnts  ,ntopolon  ,ntopolat ,ftopo, &
               mlons    ,mlats     ,plon      ,plat     ,fland)
!
!----------------------------------------------------------------------------
! Calculate standard deviation of elevation from 30sec to 10min
!----------------------------------------------------------------------------

   if (have_sgh30) then
      call binf2c(tloncnts ,tlatcnts  ,ntopolon  ,ntopolat ,variance, &
           mlons    ,mlats     ,plon      ,plat     ,sgh30)
   else
      sgh30 = -1
   endif
!-------------------------------------------------------------------------
! Calculate determine mean and variance of topographic height, plon >=128 
!-------------------------------------------------------------------------
!
   if (plon >= 128) then
      call binf2c(tloncnts ,tlatcnts  ,ntopolon  ,ntopolat ,htopo, &
                  mlons    ,mlats     ,plon      ,plat     ,mnhgt)

      call varf2c(tloncnts ,tlatcnts  ,ntopolon  ,ntopolat ,htopo  , &
                  mlons   ,mlats     ,plon      ,plat     ,mnhgt  , &
                  varhgt  )
   end if

!-------------------------------------------------------------------------
! Calculate determine mean and variance of topographic height, plon < 128 
!-------------------------------------------------------------------------

   if (plon < 128) then
!
! bin to uniform 3x3 deg grid then area avg to output grid
! get 3x3 cell boundaries for binning routine
!
      dy = 180./n3x3lat
      do j = 1, n3x3lat+1
         lats3x3(j) = -90.0 + (j-1)*dy
      end do

      num3x3lons(:) = n3x3lon
      do j = 1,n3x3lat
         dx = 360./(num3x3lons(j))
         do i = 1, num3x3lons(j)+1
            lons3x3(i,j) =  0. + (i-1)*dx
         end do
      end do
!
! bin mean height to intermed grid
!
      call binf2c (tloncnts, tlatcnts, ntopolon, ntopolat, htopo,  &
                   lons3x3 , lats3x3 , n3x3lon , n3x3lat , mnhgt3x3)
!
! get variation of topography mean height over the intermed grid
!
      call varf2c (tloncnts, tlatcnts, ntopolon, ntopolat, htopo   , &
                   lons3x3 , lats3x3 , n3x3lon , n3x3lat , mnhgt3x3, &
                   varhgt3x3 )
!
! get maximum number of 3x3 cells which will to be used in area average
! for each model cell
!
      call max_ovr (n3x3lon, n3x3lat, num3x3lons, plon  , plat, numlons, &
                    lons3x3, lats3x3, mlons     , mlats , mxovr   )
!
! do area average from intermediate regular grid to gauss grid
! get memory for pointer based arrays
!
      allocate(iovr(plon,plat,mxovr))
      allocate(jovr(plon,plat,mxovr))
      allocate(wovr(plon,plat,mxovr))

      call map_i (n3x3lon, n3x3lat, num3x3lons, lons3x3, lats3x3, &
                  plon   , plat   , numlons   , mlons  , mlats  , &
                  mxovr  , iovr   , jovr      , wovr   )

      do jo = 1, plat
         do io = 1, numlons(jo)
            mnhgt(io,jo)  = 0.
            varhgt(io,jo) = 0.
            do n = 1, mxovr             ! overlap cell index
               ii = iovr(io,jo,n)       ! lon index (input grid) of overlap cell
               ji = jovr(io,jo,n)       ! lat index (input grid) of overlap cell
               wt = wovr(io,jo,n)       ! overlap weight
               mnhgt(io,jo)  = mnhgt(io,jo)  + mnhgt3x3(ii,ji)  * wt
               varhgt(io,jo) = varhgt(io,jo) + varhgt3x3(ii,ji) * wt
            end do
         end do
      end do
    
! If model grid contains pole points, then overwrite above values of phis and sgh at the
! poles with average of values of nearest 2x2 band - this is a fair approximation and
! is done so that above mapping routines do not have to be rewritten to correctly evaulte
! the area average of the pole points

      if (mlatcnts(1)-eps < -90.0 .and. mlatcnts(plat)+eps > 90.0) then
         write(6,*)' determining sgh and phis at poles'
         summn = 0
         sumvar = 0
         do io = 1,numlons(2)
            summn  = summn + mnhgt(io,2)
            sumvar = sumvar + varhgt(io,2)
         end do
         do io = 1,numlons(1)
            mnhgt(io,1)  = summn/numlons(2)
            varhgt(io,1) = sumvar/numlons(2)
         end do
         summn = 0
         sumvar = 0
         do io = 1,numlons(plat-1)
            summn  = summn + mnhgt(io,plat-1)
            sumvar = sumvar + varhgt(io,plat-1)
         end do
         do io = 1,numlons(plat)
            mnhgt(io,plat)  = summn/numlons(plat-1)
            varhgt(io,plat) = sumvar/numlons(plat-1)
         end do
      endif
     
      deallocate(iovr)
      deallocate(jovr)
      deallocate(wovr)

   end if

! 1-2-1 smoothing for variation height

   call sm121(varhgt,plon,plat,numlons)
   call sm121(varhgt,plon,plat,numlons)
   if (have_sgh30) then
      call sm121(sgh30,plon,plat,numlons)
      call sm121(sgh30,plon,plat,numlons)
   end if
!
! get standard deviation for smoothed height field
!
! determine geopotential height field.  The multiplication by 9.80616 
! causes phis to be only accurate to 32-bit roundoff on some machines
!
   xmax = -1.d99
   do jo=1,plat
      do io=1,numlons(jo)
         if (varhgt(io,jo) < 0.5) then
            sgh(io,jo) = 0.
         else
            sgh(io,jo) = sqrt(varhgt(io,jo))
         end if
         if (have_sgh30) then
            if (sgh30(io,jo) < 0.5) then
               sgh30(io,jo) = 0.
            else
               sgh30(io,jo) = sqrt(sgh30(io,jo))
            end if
         end if
         if (sgh(io,jo) > xmax) then
            xmax = sgh(io,jo)
            imax = io
            jmax = jo
         end if
         phis(io,jo) = mnhgt(io,jo) * 9.80616
      end do
   end do
   
   if (verbose) write(6,*)'Max SGH =',xmax,' at i,j=', imax, jmax
  
   return
end subroutine sghphis
