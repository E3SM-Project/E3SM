module horizontal_interpolate

!   
! Modules Used:
!
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use shr_const_mod,  only: SHR_CONST_PI
  use abortutils,     only: endrun
  use scamMod,        only: single_column
  use cam_logfile,    only: iulog
  implicit none
  private
  save

  real(r8) :: gw1(1000), gw2(1000)

  public :: xy_interp_init, xy_interp

contains
  subroutine xy_interp_init(im1,jm1,lon0,lat0,im2,jm2,weight_x,weight_y)
!------------------------------------------------------------------------------------------------------------
! This program computes weighting functions to map a variable of (im1,jm1) resolution to (im2,jm2) resolution
! weight_x(im2,im1) is the weighting function for zonal interpolation
! weight_y(jm2,jm1) is the weighting function for meridional interpolation
! 
! Author: Chih-Chieh (Jack) Chen  -- May 2010
!
!------------------------------------------------------------------------------------------------------------
  implicit none
  integer,  intent(in)  :: im1, jm1, im2, jm2
  real(r8), intent(in)  :: lon0(im1), lat0(jm1)
  real(r8), intent(out) :: weight_x(im2,im1), weight_y(jm2,jm1)

  real(r8) :: lon1(im1), lat1(jm1)
  real(r8) :: lon2(im2), lat2(jm2)
  real(r8) :: slon1(im1+1), slon2(im2+1), slat1(jm1+1), slat2(jm2+1)  
  real(r8) :: x1_west, x1_east, x2_west, x2_east
  real(r8) :: y1_south, y1_north, y2_south, y2_north
  integer  :: i1, j1, i2, j2

  weight_x(:,:) = 0.0_r8
  weight_y(:,:) = 0.0_r8

! lon0 & lat0 are longitude & latitude on the source mesh in radians
! convert lon1, lat1 from radians to degrees
  lon1(:) = lon0(:)/SHR_CONST_PI*180.0_r8
  lat1(:) = lat0(:)/SHR_CONST_PI*180.0_r8

! set up lon2, lat2 (target mesh), in CAM convention
  do i2=1,im2
   lon2(i2) = (float(i2)-1.0_r8)*360.0_r8/float(im2)
  enddo
  do j2=1,jm2
   lat2(j2) = -90.0_r8+(float(j2)-1.0_r8)*180.0_r8/(float(jm2)-1.0_r8)
  enddo


! set up staggered longitudes (cell edges in x)
  do i1=2,im1
	slon1(i1) = (lon1(i1-1)+lon1(i1))/2.0_r8
  enddo
  slon1(1) = lon1(1)-(lon1(2)-lon1(1))/2.0_r8
  slon1(im1+1) = lon1(im1)+(lon1(im1)-lon1(im1-1))/2.0_r8

  do i2=2,im2
	slon2(i2) = (lon2(i2-1)+lon2(i2))/2.0_r8
  enddo
  slon2(1) = lon2(1)-(lon2(2)-lon2(1))/2.0_r8
  slon2(im2+1) = lon2(im2)+(lon2(im2)-lon2(im2-1))/2.0_r8

! set up staggered lattiudes (cell edges in y)
  slat1(1)=-90.0_r8
  do j1=2,jm1
	slat1(j1) = (lat1(j1-1)+lat1(j1))/2.0_r8
  enddo
  slat1(jm1+1)=90.0_r8

  slat2(1)=-90.0_r8
  do j2=2,jm2
	slat2(j2)=(lat2(j2-1)+lat2(j2))/2.0_r8
  enddo	
  slat2(jm2+1)=90.0_r8

! compute Guassian weight for two meshes (discrete form of cos(lat).)
  do j1=1,jm1
   	gw1(j1) = sin(slat1(j1+1)/180.0_r8*SHR_CONST_PI)-sin(slat1(j1)/180.0_r8*SHR_CONST_PI)
  enddo

  do j2=1,jm2
	gw2(j2) = sin(slat2(j2+1)/180.0_r8*SHR_CONST_PI)-sin(slat2(j2)/180.0_r8*SHR_CONST_PI)
  enddo


! add 360 to slon1 and slon2 
  slon1(:) = slon1(:)+360.0_r8
  slon2(:) = slon2(:)+360.0_r8

 do i2=1,im2

! target grid east-west boundaries
  x2_west=slon2(i2)
  x2_east=slon2(i2+1)

  do i1=1,im1

! source grid east-west boundaries
   x1_west=slon1(i1)
   x1_east=slon1(i1+1)

! check if there is any overlap between the source grid and the target grid
! if no overlap, then weighting is zero
! there are three scenarios overlaps can take place 
   if( (x1_west.ge.x2_west).and.(x1_east.le.x2_east) ) then
! case 1: 
!                x1_west             x1_east
!                  |-------------------|
!            |---------------------------------|
!          x2_west                           x2_east
     weight_x(i2,i1) =  (x1_east-x1_west)/(x2_east-x2_west)
   elseif ( (x1_west.ge.x2_west).and.(x1_west.lt.x2_east) ) then
! case 2: 
!                x1_west                          x1_east
!                  |--------------------------------|
!            |---------------------------------|
!          x2_west                           x2_east
     weight_x(i2,i1) = (x2_east-x1_west)/(x2_east-x2_west)
   elseif ( (x1_east>x2_west).and.(x1_east.le.x2_east) ) then
! case 3: 
!       x1_west                          x1_east
!         |--------------------------------|
!                |---------------------------------|
!              x2_west                           x2_east
     weight_x(i2,i1) = (x1_east-x2_west)/(x2_east-x2_west)
   endif

   enddo	
  enddo


! consider end points
      if(slon1(im1+1).gt.slon2(im2+1)) then
! case 1:
!           slon1(im1)                slon1(im1+1) <--- end point
!              |-------------------------|
!           |----------------|......................|
!        slon2(im2)         slon2(im2+1)        slon2(2)  (note: slon2(im2+1) = slon2(1))
     	weight_x(1,im1)= weight_x(1,im1)+(slon1(im1+1)-slon2(im2+1))/(slon2(2)-slon2(1))
      endif	

      if(slon1(im1+1).lt.slon2(im2+1)) then
! case 1:
!           slon1(im1)                slon1(im1+1)                  slon1(2)    (note: slon1(im1+1) = slon1(1))
!              |-------------------------|.............................|
!                   |-------------------------------|
!               slon2(im2)                        slon2(im2+1) <--- end point
        weight_x(im2,1) = weight_x(im2,1)+(slon2(1)-slon1(1))/(slon2(2)-slon2(1)) 
      endif



      do j2=1,jm2
! target grid north-south boundaries
         y2_south=slat2(j2)
         y2_north=slat2(j2+1)

         do j1=1,jm1

! source grid north-south boundaries
            y1_south=slat1(j1)
            y1_north=slat1(j1+1)
 
! check if there is any overlap between the source grid and the target grid
! if no overlap, then weighting is zero
! there are three scenarios overlaps can take place 
! note: there is Guassian weight to consider in the meridional direction!

            if( (y1_south.ge.y2_south).and.(y1_north.le.y2_north) ) then
! case 1: 
!                y1_south             y1_north
!                  |-------------------|
!            |---------------------------------|
!          y2_south                           y2_north
                 weight_y(j2,j1) =  gw1(j1)/gw2(j2)
            elseif ( (y1_south.ge.y2_south).and.(y1_south.lt.y2_north) ) then
! case 2: 
!                y1_south                          y1_north
!                  |--------------------------------|
!            |---------------------------------|
!          y2_south                           y2_north
                 weight_y(j2,j1) = (y2_north-y1_south)/(y1_north-y1_south)*gw1(j1)/gw2(j2)
            elseif ( (y1_north.gt.y2_south).and.(y1_north.le.y2_north) ) then
! case 3: 
!       y1_south                          y1_north
!         |--------------------------------|
!                |---------------------------------|
!              y2_south                           y2_north
                 weight_y(j2,j1) = (y1_north-y2_south)/(y1_north-y1_south)*gw1(j1)/gw2(j2)
            endif

          enddo
        enddo

  end subroutine xy_interp_init

  subroutine xy_interp(im1,jm1,km1,im2,jm2,pcols,ncols,weight_x,weight_y,var_src,var_trg,lons,lats,count_x,count_y,index_x,index_y)
!-------------------------------------------------------------------------------------------------------------
! This program interpolates var_src(im1,jm1,km1) to var_trg(im2,jm2,km1) based on weighting functions weight_x & weight_y.
!-------------------------------------------------------------------------------------------------------------
  implicit none
  integer,  intent(in)  :: im1   ! source number of longitudes
  integer,  intent(in)  :: jm1   ! source number of latitudes
  integer,  intent(in)  :: km1   ! source/target number of levels
  integer,  intent(in)  :: im2   ! target number of longitudes
  integer,  intent(in)  :: jm2   ! target number of latitudes
  integer,  intent(in)  :: pcols
  integer,  intent(in)  :: ncols
  real(r8), intent(in)  :: weight_x(im2,im1), weight_y(jm2,jm1)
  real(r8), intent(in)  :: var_src(im1,jm1,km1)
  integer,  intent(in)  :: lons(pcols), lats(pcols)
  integer,  intent(in)  :: count_x(im2), count_y(jm2)
  integer,  intent(in)  :: index_x(im2,im1), index_y(jm2,jm1)
  real(r8), intent(out) :: var_trg(pcols,km1)
  integer  :: n, i1, j1, k1, i2, j2, ii, jj
  real(r8) :: sum_x

  var_trg(:,:) = 0.0_r8
 

  do k1=1,km1
     do n=1,ncols
! interpolate in x
        do jj=1,count_y(lats(n))
           sum_x = 0.0_r8
           do ii=1,count_x(lons(n))
              sum_x = sum_x + var_src(index_x(lons(n),ii),index_y(lats(n),jj),k1)* &
                             weight_x(lons(n),index_x(lons(n),ii))
           enddo
           var_trg(n,k1) = var_trg(n,k1)+sum_x*weight_y(lats(n),index_y(lats(n),jj))
        enddo
     enddo
  enddo

  end subroutine xy_interp


end module horizontal_interpolate
