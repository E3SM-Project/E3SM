!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     This module reads in and initializes two grids for remapping.
!     NOTE: grid1 must be the master grid -- the grid that determines
!           which cells participate (e.g. land mask) and the fractional
!           area of grid2 cells that participate in the remapping.
!
!-----------------------------------------------------------------------
!
!     CVS:$Id: grids.f,v 1.6 2001/08/21 21:06:41 pwjones Exp $
!
!     Copyright (c) 1997, 1998 the Regents of the University of 
!       California.
!
!     This software and ancillary information (herein called software) 
!     called SCRIP is made available under the terms described here.  
!     The software has been approved for release with associated 
!     LA-CC Number 98-45.
!
!     Unless otherwise indicated, this software has been authored
!     by an employee or employees of the University of California,
!     operator of the Los Alamos National Laboratory under Contract
!     No. W-7405-ENG-36 with the U.S. Department of Energy.  The U.S.
!     Government has rights to use, reproduce, and distribute this
!     software.  The public may copy and use this software without
!     charge, provided that this Notice and any statement of authorship
!     are reproduced on all copies.  Neither the Government nor the
!     University makes any warranty, express or implied, or assumes
!     any liability or responsibility for the use of this software.
!
!     If software is modified to produce derivative works, such modified
!     software should be clearly marked, so as not to confuse it with 
!     the version available from Los Alamos National Laboratory.
!
!     this grids.f has been modified from the version available from 
!     Los Alamos National Laboratory.
!     modifications are marked with "NRL"
!     list of modifications:
!     - netcdf operations for reading from netcdf file removed
!           (major change)
!     - print statements added
!     - allocations removed (moved to scrip_wrapper subroutine)
!     - imask removed (this was an intermediate step from ncdf to 
!           the logicals grid1_mask grid2_mask; creation of these 
!           logicals is now handled by scrip_wrapper)
!
!***********************************************************************

      module scrip_grids

!-----------------------------------------------------------------------

      use SCRIP_KindsMod    ! defines data types
      use SCRIP_ErrorMod    ! error tracking
      use SCRIP_IOUnitsMod  ! manages I/O units
      use SCRIP_constants    ! common constants
!      use SCRIP_NetcdfMod   ! netCDF stuff
!      use netcdf            ! netCDF library module

      implicit none

!-----------------------------------------------------------------------
!
!     variables that describe each grid
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4), save ::
     &             grid1_size, grid2_size, ! total points on each grid
     &             grid1_rank, grid2_rank, ! rank of each grid
     &             grid1_corners, grid2_corners ! number of corners
                                                ! for each grid cell

      integer (SCRIP_i4), dimension(:), allocatable, save ::
     &             grid1_dims, grid2_dims  ! size of each grid dimension

      character(SCRIP_charLength), save :: 
     &             grid1_name, grid2_name  ! name for each grid

      character (SCRIP_charLength), save :: 
     &             grid1_units, ! units for grid coords (degs/radians)
     &             grid2_units  ! units for grid coords

      real (SCRIP_r8), parameter ::
     &      deg2rad = pi/180.   ! conversion for deg to rads

!-----------------------------------------------------------------------
!
!     grid coordinates and masks
!
!-----------------------------------------------------------------------

      logical (SCRIP_logical), dimension(:), allocatable, target,save ::
     &             grid1_mask,        ! flag which cells participate
     &             grid2_mask,        ! flag which cells participate
     &     special_polar_cell1,       ! cell with only 1 corner at pole
     &     special_polar_cell2        !

      real (SCRIP_r8), dimension(:), allocatable, target, save ::
     &             grid1_center_lat,  ! lat/lon coordinates for
     &             grid1_center_lon,  ! each grid center in radians
     &             grid2_center_lat, 
     &             grid2_center_lon,
     &             grid1_area,        ! tot area of each grid1 cell
     &             grid2_area,        ! tot area of each grid2 cell
     &             grid1_area_in,     ! area of grid1 cell from file
     &             grid2_area_in,     ! area of grid2 cell from file
     &             grid1_frac,        ! fractional area of grid cells
     &             grid2_frac,        ! participating in remapping
     &             grid1_centroid_lat,! Centroid of grid1 cell
     &             grid1_centroid_lon,!
     &             grid2_centroid_lat,! Centroid of grid2 cell
     &             grid2_centroid_lon !


      real (SCRIP_r8), dimension(:,:), allocatable, target, save ::
     &             grid1_corner_lat,  ! lat/lon coordinates for
     &             grid1_corner_lon,  ! each grid corner in radians
     &             grid2_corner_lat, 
     &             grid2_corner_lon

      logical (SCRIP_logical), save ::
     &             luse_grid_centers ! use centers for bounding boxes
     &,            luse_grid1_area   ! use area from grid file
     &,            luse_grid2_area   ! use area from grid file

      real (SCRIP_r8), dimension(:,:), allocatable, target, save ::
     &             grid1_bound_box,  ! lat/lon bounding box for use
     &             grid2_bound_box   ! in restricting grid searches

      integer (SCRIP_i4), save ::    ! Cells overlapping the poles 
                                     ! (may be 0)
     &     grid1_npole_cell,        
     &     grid1_spole_cell,         
     &     grid2_npole_cell,
     &     grid2_spole_cell
      

!-----------------------------------------------------------------------
!
!     bins for restricting searches
!
!-----------------------------------------------------------------------

      character (SCRIP_charLength), save ::
     &        restrict_type  ! type of bins to use

      integer (SCRIP_i4), save ::
     &        num_srch_bins  ! num of bins for restricted srch

      integer (SCRIP_i4), dimension(:,:), allocatable, save ::
     &        bin_addr1, ! min,max adds for grid1 cells in this lat bin
     &        bin_addr2  ! min,max adds for grid2 cells in this lat bin

      real(SCRIP_r8), dimension(:,:), allocatable, save ::
     &        bin_lats   ! min,max latitude for each search bin
     &,       bin_lons   ! min,max longitude for each search bin

!-----------------------------------------------------------------------
!
!     Parameters for dealing with intersections in the polar region
!
!-----------------------------------------------------------------------

      real (SCRIP_r8), save :: 
     &     north_thresh,  ! threshold for coord transf.
     &     south_thresh   ! threshold for coord transf.


      !*** Number of subsegments used to represents edges near
      !*** the polar regions - choose an odd number to avoid obvious
      !*** degeneracies in intersection

      integer (SCRIP_i4), save ::
     &     npseg                    

!***********************************************************************

      contains

!***********************************************************************

      subroutine grid_init(errorCode,l_master,l_test)

!-----------------------------------------------------------------------
!
!     this routine reads grid info from grid files and makes any
!     necessary changes (e.g. for 0,2pi longitude range)
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     input variables
!
!-----------------------------------------------------------------------

      logical(SCRIP_Logical), intent(in) :: l_master  ! Am I the master
                                                  ! processor (do I/O)?
      logical(SCRIP_Logical), intent(in) :: l_test ! Whether to include
                                                   ! test output

!-----------------------------------------------------------------------
!
!     output variables
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4), intent(out) ::
     &   errorCode                 ! returned error code

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4) :: 
     &  n,      ! loop counter
     &  nele,   ! element loop counter
     &  i,j,
     &  ip1,jp1,
     &  n_add, e_add, ne_add,
     &  nx, ny, ncorners_at_pole

      integer (SCRIP_i4) ::
     &     zero_crossing, pi_crossing,
     &     grid1_add, grid2_add,
     &     corner, next_corn
     
      real (SCRIP_r8) ::
     &     beglon, beglat, endlon, endlat

      logical (SCRIP_logical) ::
     &     found

!NRL      integer (SCRIP_i4), dimension(:), allocatable :: 
!NRL     &                            imask ! integer mask read from file

      real (SCRIP_r8) :: 
     &  dlat,dlon           ! lat/lon intervals for search bins

      real (SCRIP_r8), dimension(4) ::
     &  tmp_lats, tmp_lons  ! temps for computing bounding boxes

      character (9), parameter ::
     &   rtnName = 'grid_init'

      if(l_master.and.l_test)write(SCRIP_stdout,*)'subroutine grid_init'

!NRL #####################################################
!NRL .....Begin part of code formerly handled by netcdf read
!NRL #####################################################

!-----------------------------------------------------------------------
!
!     allocate grid coordinates/masks and read data
!
!-----------------------------------------------------------------------

      allocate( 
!NRL &          grid1_mask      (grid1_size),
!NRL &          grid2_mask      (grid2_size),
     &          special_polar_cell1(grid1_size),
     &          special_polar_cell2(grid2_size),
!NRL &          grid1_center_lat(grid1_size), 
!NRL &          grid1_center_lon(grid1_size),
!NRL &          grid2_center_lat(grid2_size), 
!NRL &          grid2_center_lon(grid2_size),
     &          grid1_area      (grid1_size),
     &          grid1_area_in   (grid1_size),
     &          grid2_area      (grid2_size),
     &          grid2_area_in   (grid2_size),
     &          grid1_frac      (grid1_size),
     &          grid2_frac      (grid2_size),
!NRL &          grid1_corner_lat(grid1_corners, grid1_size),
!NRL &          grid1_corner_lon(grid1_corners, grid1_size),
!NRL &          grid2_corner_lat(grid2_corners, grid2_size),
!NRL &          grid2_corner_lon(grid2_corners, grid2_size),
     &          grid1_bound_box (4            , grid1_size),
     &          grid2_bound_box (4            , grid2_size),
     &          grid1_centroid_lat(grid1_size),
     &          grid1_centroid_lon(grid1_size),
     &          grid2_centroid_lat(grid2_size),
     &          grid2_centroid_lon(grid2_size))

!NRL      allocate(imask(grid1_size))

      grid1_area = zero
      grid1_frac = zero
      grid1_centroid_lat = zero
      grid1_centroid_lon = zero

!-----------------------------------------------------------------------
!
!     initialize logical mask and convert lat/lon units if required
!
!-----------------------------------------------------------------------

!NRL      where (imask == 1)
!NRL        grid1_mask = .true.
!NRL      elsewhere
!NRL        grid1_mask = .false.
!NRL      endwhere
!NRL      deallocate(imask)

      select case (grid1_units(1:7))
      case ('degrees')

        grid1_center_lat = grid1_center_lat*deg2rad
        grid1_center_lon = grid1_center_lon*deg2rad

      case ('radians')

        !*** no conversion necessary

      case default

        print *,'unknown units supplied for grid1 center lat/lon: '
        print *,'proceeding assuming radians'

      end select

      select case (grid1_units(1:7))
      case ('degrees')

        grid1_corner_lat = grid1_corner_lat*deg2rad
        grid1_corner_lon = grid1_corner_lon*deg2rad

      case ('radians')

        !*** no conversion necessary

      case default

        print *,'unknown units supplied for grid1 corner lat/lon: '
        print *,'proceeding assuming radians'

      end select

!-----------------------------------------------------------------------
!
!     read data for grid 2
!
!-----------------------------------------------------------------------

!NRL      allocate(imask(grid2_size))

      grid2_area = zero
      grid2_frac = zero
      grid2_centroid_lat = zero
      grid2_centroid_lon = zero

!-----------------------------------------------------------------------
!
!     initialize logical mask and convert lat/lon units if required
!
!-----------------------------------------------------------------------

!NRL      where (imask == 1)
!NRL        grid2_mask = .true.
!NRL      elsewhere
!NRL        grid2_mask = .false.
!NRL      endwhere
!NRL      deallocate(imask)

      select case (grid2_units(1:7))
      case ('degrees')

        grid2_center_lat = grid2_center_lat*deg2rad
        grid2_center_lon = grid2_center_lon*deg2rad

      case ('radians')

        !*** no conversion necessary

      case default

        print *,'unknown units supplied for grid2 center lat/lon: '
        print *,'proceeding assuming radians'

      end select

      select case (grid2_units(1:7))
      case ('degrees')

        grid2_corner_lat = grid2_corner_lat*deg2rad
        grid2_corner_lon = grid2_corner_lon*deg2rad

      case ('radians')

        !*** no conversion necessary

      case default

        print *,'no units supplied for grid2 corner lat/lon: '
        print *,'proceeding assuming radians'

      end select

!NRL #####################################################
!NRL .....End part of code formerly handled by netcdf read
!NRL #####################################################


!-----------------------------------------------------------------------
!
!     convert longitudes to 0,2pi interval
!
!-----------------------------------------------------------------------

      where (grid1_center_lon .gt. pi2)  grid1_center_lon =
     &                                   grid1_center_lon - pi2
      where (grid1_center_lon .lt. zero) grid1_center_lon =
     &                                   grid1_center_lon + pi2
      where (grid2_center_lon .gt. pi2)  grid2_center_lon =
     &                                   grid2_center_lon - pi2
      where (grid2_center_lon .lt. zero) grid2_center_lon =
     &                                   grid2_center_lon + pi2
      where (grid1_corner_lon .gt. pi2)  grid1_corner_lon =
     &                                   grid1_corner_lon - pi2
      where (grid1_corner_lon .lt. zero) grid1_corner_lon =
     &                                   grid1_corner_lon + pi2
      where (grid2_corner_lon .gt. pi2)  grid2_corner_lon =
     &                                   grid2_corner_lon - pi2
      where (grid2_corner_lon .lt. zero) grid2_corner_lon =
     &                                   grid2_corner_lon + pi2

!-----------------------------------------------------------------------
!
!     make sure input latitude range is within the machine values
!     for +/- pi/2 
!
!-----------------------------------------------------------------------

      where (grid1_center_lat >  pih) grid1_center_lat =  pih
      where (grid1_corner_lat >  pih) grid1_corner_lat =  pih
      where (grid1_center_lat < -pih) grid1_center_lat = -pih
      where (grid1_corner_lat < -pih) grid1_corner_lat = -pih

      where (grid2_center_lat >  pih) grid2_center_lat =  pih
      where (grid2_corner_lat >  pih) grid2_corner_lat =  pih
      where (grid2_center_lat < -pih) grid2_center_lat = -pih
      where (grid2_corner_lat < -pih) grid2_corner_lat = -pih


!-----------------------------------------------------------------------
!
!     also, different grids consider the pole to be a slightly different
!     values (1.570796326789 vs 1.5707963267977). Find the closest
!     approach to the pole and if it is within a tolerance, move such
!     points that are practically at the pole to the pole to avoid
!     problems
!
!-----------------------------------------------------------------------
      
      where (abs(grid1_corner_lat-pih) < 1.e-03) grid1_corner_lat =  pih
      where (abs(grid1_corner_lat+pih) < 1.e-03) grid1_corner_lat = -pih
      where (abs(grid2_corner_lat-pih) < 1.e-03) grid2_corner_lat =  pih
      where (abs(grid2_corner_lat+pih) < 1.e-03) grid2_corner_lat = -pih


!-----------------------------------------------------------------------
!
!     compute bounding boxes for restricting future grid searches
!
!-----------------------------------------------------------------------

      if (.not. luse_grid_centers) then
        grid1_bound_box(1,:) = minval(grid1_corner_lat, DIM=1)
        grid1_bound_box(2,:) = maxval(grid1_corner_lat, DIM=1)
        grid1_bound_box(3,:) = minval(grid1_corner_lon, DIM=1)
        grid1_bound_box(4,:) = maxval(grid1_corner_lon, DIM=1)

        grid2_bound_box(1,:) = minval(grid2_corner_lat, DIM=1)
        grid2_bound_box(2,:) = maxval(grid2_corner_lat, DIM=1)
        grid2_bound_box(3,:) = minval(grid2_corner_lon, DIM=1)
        grid2_bound_box(4,:) = maxval(grid2_corner_lon, DIM=1)

      else

        nx = grid1_dims(1)
        ny = grid1_dims(2)

        do n=1,grid1_size

          !*** find N,S and NE points to this grid point

          j = (n - 1)/nx +1
          i = n - (j-1)*nx

          if (i < nx) then
            ip1 = i + 1
          else
            !*** assume cyclic
            ip1 = 1
            !*** but if it is not, correct
            e_add = (j - 1)*nx + ip1
            if (abs(grid1_center_lat(e_add) - 
     &              grid1_center_lat(n   )) > pih) then
              ip1 = i
            endif
          endif

          if (j < ny) then
            jp1 = j+1
          else
            !*** assume cyclic
            jp1 = 1
            !*** but if it is not, correct
            n_add = (jp1 - 1)*nx + i
            if (abs(grid1_center_lat(n_add) - 
     &              grid1_center_lat(n   )) > pih) then
              jp1 = j
            endif
          endif

          n_add = (jp1 - 1)*nx + i
          e_add = (j - 1)*nx + ip1
          ne_add = (jp1 - 1)*nx + ip1

          !*** find N,S and NE lat/lon coords and check bounding box

          tmp_lats(1) = grid1_center_lat(n)
          tmp_lats(2) = grid1_center_lat(e_add)
          tmp_lats(3) = grid1_center_lat(ne_add)
          tmp_lats(4) = grid1_center_lat(n_add)

          tmp_lons(1) = grid1_center_lon(n)
          tmp_lons(2) = grid1_center_lon(e_add)
          tmp_lons(3) = grid1_center_lon(ne_add)
          tmp_lons(4) = grid1_center_lon(n_add)

          grid1_bound_box(1,n) = minval(tmp_lats)
          grid1_bound_box(2,n) = maxval(tmp_lats)
          grid1_bound_box(3,n) = minval(tmp_lons)
          grid1_bound_box(4,n) = maxval(tmp_lons)
        end do

        nx = grid2_dims(1)
        ny = grid2_dims(2)

        do n=1,grid2_size

          !*** find N,S and NE points to this grid point

          j = (n - 1)/nx +1
          i = n - (j-1)*nx

          if (i < nx) then
            ip1 = i + 1
          else
            !*** assume cyclic
            ip1 = 1
            !*** but if it is not, correct
            e_add = (j - 1)*nx + ip1
            if (abs(grid2_center_lat(e_add) - 
     &              grid2_center_lat(n   )) > pih) then
              ip1 = i
            endif
          endif

          if (j < ny) then
            jp1 = j+1
          else
            !*** assume cyclic
            jp1 = 1
            !*** but if it is not, correct
            n_add = (jp1 - 1)*nx + i
            if (abs(grid2_center_lat(n_add) - 
     &              grid2_center_lat(n   )) > pih) then
              jp1 = j
            endif
          endif

          n_add = (jp1 - 1)*nx + i
          e_add = (j - 1)*nx + ip1
          ne_add = (jp1 - 1)*nx + ip1

          !*** find N,S and NE lat/lon coords and check bounding box

          tmp_lats(1) = grid2_center_lat(n)
          tmp_lats(2) = grid2_center_lat(e_add)
          tmp_lats(3) = grid2_center_lat(ne_add)
          tmp_lats(4) = grid2_center_lat(n_add)

          tmp_lons(1) = grid2_center_lon(n)
          tmp_lons(2) = grid2_center_lon(e_add)
          tmp_lons(3) = grid2_center_lon(ne_add)
          tmp_lons(4) = grid2_center_lon(n_add)

          grid2_bound_box(1,n) = minval(tmp_lats)
          grid2_bound_box(2,n) = maxval(tmp_lats)
          grid2_bound_box(3,n) = minval(tmp_lons)
          grid2_bound_box(4,n) = maxval(tmp_lons)
        end do

      endif

      where (abs(grid1_bound_box(4,:) - grid1_bound_box(3,:)) > pi)
        grid1_bound_box(3,:) = zero
        grid1_bound_box(4,:) = pi2
      end where

      where (abs(grid2_bound_box(4,:) - grid2_bound_box(3,:)) > pi)
        grid2_bound_box(3,:) = zero
        grid2_bound_box(4,:) = pi2
      end where

      where (grid1_center_lat > grid1_bound_box(2,:))
     &  grid1_bound_box(2,:) = pih

      where (grid1_center_lat < grid1_bound_box(1,:))
     &  grid1_bound_box(1,:) = -pih

      where (grid2_center_lat > grid2_bound_box(2,:))
     &  grid2_bound_box(2,:) = pih

      where (grid2_center_lat < grid2_bound_box(1,:))
     &  grid2_bound_box(1,:) = -pih


      !***
      !*** Check for cells that overlap poles and explicitly
      !*** store their addresses
      !***

      grid1_npole_cell = 0
      grid1_spole_cell = 0
         
      do grid1_add = 1, grid1_size

         found = .false.
         do corner = 1, grid1_corners
            endlat = grid1_corner_lat(corner,grid1_add)
            if (abs(abs(endlat)-pih) .lt. 1e-5) then
               found = .true.   ! cell has polar pnt; so pole is 
                                ! not in the interior of the cell
               exit
            endif
         enddo

         if (found) cycle


         beglon = grid1_corner_lon(1,grid1_add)
         zero_crossing = 0
         pi_crossing = 0

         do corner = 1, grid1_corners
            next_corn = mod(corner,grid1_corners) + 1
            endlon = grid1_corner_lon(next_corn,grid1_add)

            if (abs(beglon-endlon) .gt. pi) then
               zero_crossing = 1
            else
               if ((beglon .lt. pi .and. endlon .ge. pi) .or.
     &              (endlon .lt. pi .and. beglon .ge. pi)) then
                  pi_crossing = 1
               endif
            endif

            beglon = endlon
         enddo
         
         if (zero_crossing .eq. 1 .and. pi_crossing .eq. 1) then

            !***
            !*** We have a polar cell
            !***

            if (grid1_center_lat(grid1_add) .gt. 0) then
               grid1_npole_cell = grid1_add
            else if (grid1_center_lat(grid1_add) .lt. 0) then
               grid1_spole_cell = grid1_add
            endif

            if (grid1_npole_cell .ne. 0 .and.
     &           grid1_spole_cell .ne. 0) then
               exit
            endif

         endif

      enddo



      grid2_npole_cell = 0
      grid2_spole_cell = 0

      do grid2_add = 1, grid2_size

         found = .false.
         do corner = 1, grid2_corners
            endlat = grid2_corner_lat(corner,grid2_add)
            if (abs(abs(endlat)-pih) .lt. 1e-5) then
               found = .true.   ! cell has polar pnt; so pole is 
                                ! not in the interior of the cell
               exit
            endif
         enddo

         if (found) cycle

         beglon = grid2_corner_lon(1,grid2_add)
         zero_crossing = 0
         pi_crossing = 0

         do corner = 1, grid2_corners
            next_corn = mod(corner,grid2_corners) + 1
            endlon = grid2_corner_lon(next_corn,grid2_add)
            
            if (abs(beglon-endlon) > pi) then
               zero_crossing = 1
            else
               if ((beglon .lt. pi .and. endlon .ge. pi) .or.
     &              (endlon .lt. pi .and. beglon .ge. pi)) then
                  pi_crossing = 1
               endif
            endif

            beglon = endlon
         enddo
         
         if (zero_crossing .eq. 1 .and. pi_crossing .eq. 1) then

            !***
            !*** We have a polar cell
            !***

            if (grid2_center_lat(grid2_add) .gt. 0) then
               grid2_npole_cell = grid2_add
            else if (grid2_center_lat(grid2_add) .lt. 0) then
               grid2_spole_cell = grid2_add
            endif

            if (grid2_npole_cell .ne. 0 .and.
     &           grid2_spole_cell .ne. 0) then
               exit
            endif

         endif

      enddo


      special_polar_cell1 = .false.
      do grid1_add = 1, grid1_size

         ncorners_at_pole = 0
         do i = 1, grid1_corners
            beglat = grid1_corner_lat(i,grid1_add)
            if (abs(abs(beglat)-pih) .le. 1.e-5) 
     &           ncorners_at_pole = ncorners_at_pole + 1
         enddo

         if (ncorners_at_pole .eq. 1) 
     &        special_polar_cell1(grid1_add) = .true.
            
      enddo

      special_polar_cell2 = .false.
      do grid2_add = 1, grid2_size

         ncorners_at_pole = 0
         do i = 1, grid2_corners
            beglat = grid2_corner_lat(i,grid2_add)
            if (abs(abs(beglat)-pih) .le. 1.e-5) 
     &              ncorners_at_pole = ncorners_at_pole + 1
         enddo

         if (ncorners_at_pole .eq. 1) 
     &        special_polar_cell2(grid2_add) = .true.
            
      enddo

      if(l_master)print *, ' '
      if(l_master)print *, 'Grid 1 size', grid1_size
      if(l_master)print *, 'Grid 2 size', grid2_size


!     if(l_master)print *, 'grid1_npole_cell',grid1_npole_cell
      if (grid1_npole_cell .gt. 0) then
         do i = 1, grid1_corners
            print *, grid1_corner_lat(i,grid1_npole_cell),
     &           grid1_corner_lon(i,grid1_npole_cell)
         enddo
      endif
!     if(l_master)print *, 'grid1_spole_cell',grid1_spole_cell
      if (grid1_spole_cell .gt. 0) then
         do i = 1, grid1_corners
            print *, grid1_corner_lat(i,grid1_spole_cell), 
     &           grid1_corner_lon(i,grid1_spole_cell)
         enddo
      endif
!     if(l_master)print *, 'grid2_npole_cell',grid2_npole_cell
      if (grid2_npole_cell .gt. 0) then
         do i = 1, grid2_corners
            print *, grid2_corner_lat(i,grid2_npole_cell), 
     &           grid2_corner_lon(i,grid2_npole_cell)
         enddo
      endif
!     if(l_master)print *, 'grid2_spole_cell',grid2_spole_cell
      if (grid2_spole_cell .gt. 0) then
         do i = 1, grid2_corners
            print *, grid2_corner_lat(i,grid2_spole_cell), 
     &           grid2_corner_lon(i,grid2_spole_cell)
         enddo
      endif
      if(l_master)print *


!-----------------------------------------------------------------------
!
!     set up and assign address ranges to search bins in order to 
!     further restrict later searches
!
!-----------------------------------------------------------------------

      select case (restrict_type)

      case ('latitude')
        if(l_master.and.l_test)write(SCRIP_stdout,*)
     &        'Using latitude bins to restrict search.'

        allocate(bin_addr1(2,num_srch_bins))
        allocate(bin_addr2(2,num_srch_bins))
        allocate(bin_lats (2,num_srch_bins))
        allocate(bin_lons (2,num_srch_bins))

        dlat = pi/num_srch_bins

        do n=1,num_srch_bins
          bin_lats(1,n) = (n-1)*dlat - pih
          bin_lats(2,n) =     n*dlat - pih
          bin_lons(1,n) = zero
          bin_lons(2,n) = pi2
          bin_addr1(1,n) = grid1_size + 1
          bin_addr1(2,n) = 0
          bin_addr2(1,n) = grid2_size + 1
          bin_addr2(2,n) = 0
        end do

        do nele=1,grid1_size
          do n=1,num_srch_bins
            if (grid1_bound_box(1,nele) <= bin_lats(2,n) .and.
     &          grid1_bound_box(2,nele) >= bin_lats(1,n)) then
              bin_addr1(1,n) = min(nele,bin_addr1(1,n))
              bin_addr1(2,n) = max(nele,bin_addr1(2,n))
            endif
          end do
        end do

        do nele=1,grid2_size
          do n=1,num_srch_bins
            if (grid2_bound_box(1,nele) <= bin_lats(2,n) .and.
     &          grid2_bound_box(2,nele) >= bin_lats(1,n)) then
              bin_addr2(1,n) = min(nele,bin_addr2(1,n))
              bin_addr2(2,n) = max(nele,bin_addr2(2,n))
            endif
          end do
        end do

      case ('latlon')
        if(l_master.and.l_test)write(SCRIP_stdout,*)
     &         'Using lat/lon boxes to restrict search.'

        dlat = pi /num_srch_bins
        dlon = pi2/num_srch_bins

        allocate(bin_addr1(2,num_srch_bins*num_srch_bins))
        allocate(bin_addr2(2,num_srch_bins*num_srch_bins))
        allocate(bin_lats (2,num_srch_bins*num_srch_bins))
        allocate(bin_lons (2,num_srch_bins*num_srch_bins))

        n = 0
        do j=1,num_srch_bins
        do i=1,num_srch_bins
          n = n + 1

          bin_lats(1,n) = (j-1)*dlat - pih
          bin_lats(2,n) =     j*dlat - pih
          bin_lons(1,n) = (i-1)*dlon
          bin_lons(2,n) =     i*dlon
          bin_addr1(1,n) = grid1_size + 1
          bin_addr1(2,n) = 0
          bin_addr2(1,n) = grid2_size + 1
          bin_addr2(2,n) = 0
        end do
        end do

        num_srch_bins = num_srch_bins**2

        do nele=1,grid1_size
          do n=1,num_srch_bins
            if (grid1_bound_box(1,nele) <= bin_lats(2,n) .and.
     &          grid1_bound_box(2,nele) >= bin_lats(1,n) .and.
     &          grid1_bound_box(3,nele) <= bin_lons(2,n) .and.
     &          grid1_bound_box(4,nele) >= bin_lons(1,n)) then
              bin_addr1(1,n) = min(nele,bin_addr1(1,n))
              bin_addr1(2,n) = max(nele,bin_addr1(2,n))
            endif
          end do
        end do

        do nele=1,grid2_size
          do n=1,num_srch_bins
            if (grid2_bound_box(1,nele) <= bin_lats(2,n) .and.
     &          grid2_bound_box(2,nele) >= bin_lats(1,n) .and.
     &          grid2_bound_box(3,nele) <= bin_lons(2,n) .and.
     &          grid2_bound_box(4,nele) >= bin_lons(1,n)) then
              bin_addr2(1,n) = min(nele,bin_addr2(1,n))
              bin_addr2(2,n) = max(nele,bin_addr2(2,n))
            endif
          end do
        end do

      case default
        stop 'unknown search restriction method'
      end select

!-----------------------------------------------------------------------
!
!     if area not read in, compute an area
!
!-----------------------------------------------------------------------

      if (.not. luse_grid1_area) then
         call SCRIP_GridComputeArea(grid1_area_in, grid1_corner_lat,
     &                              grid1_corner_lon, errorCode)

         if (SCRIP_ErrorCheck(errorCode, rtnName, 
     &                        'error computing grid1 area')) return
      endif

      if (.not. luse_grid2_area) then
         call SCRIP_GridComputeArea(grid2_area_in, grid2_corner_lat,
     &                              grid2_corner_lon, errorCode)

         if (SCRIP_ErrorCheck(errorCode, rtnName, 
     &                        'error computing grid2 area')) return
      endif

!-----------------------------------------------------------------------

      end subroutine grid_init

!***********************************************************************
!BOP
! !IROUTINE: SCRIP_GridComputeArea -- computes grid cell areas
! !INTERFACE:

      subroutine SCRIP_GridComputeArea(area, cornerLat, cornerLon,
     &                                 errorCode)

! !DESCRIPTION:
!  This routine computes a grid cell area based on corner lat/lon
!  coordinates.  It is provided in the case that a user supplied
!  area is not available.
!
! !REVISION HISTORY:
!  same as module

! !OUTPUT PARAMETERS:

      real (SCRIP_r8), dimension(:), intent(out) ::
     &   area              ! computed area for each grid cell

      integer (SCRIP_i4), intent(out) ::
     &   errorCode         ! returned error code

! !INPUT PARAMETERS:

      real (SCRIP_r8), dimension(:,:), intent(in) ::
     &   cornerLat,        ! latitude  of each cell corner
     &   cornerLon         ! longitude of each cell corner

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4) ::
     &   numCells,           ! number of grid cells
     &   numCorners,         ! number of corners in each cell
     &   nCell,              ! loop index for grid cells
     &   nCorner,            ! loop index for corners in each cell
     &   nextCorner          ! next corner around cell perimeter

      real (SCRIP_r8) ::
     &   dphi                ! delta(longitude) for this segment
     
!-----------------------------------------------------------------------
!
!  determine size of grid and initialize
!
!-----------------------------------------------------------------------

      errorCode = SCRIP_Success

      numCells   = size(CornerLat, dim=2)
      numCorners = size(CornerLat, dim=1)

!-----------------------------------------------------------------------
!
!  compute area for each cell by integrating around cell edge
!
!-----------------------------------------------------------------------

      do nCell=1,numCells

         Area(nCell) = 0.0_SCRIP_r8

         do nCorner=1,numCorners
            nextCorner = mod(nCorner,numCorners) + 1

            !*** trapezoid rule - delta(Lon) is -0.5*dx
            dphi = CornerLon(   nCorner,nCell) - 
     &             CornerLon(nextCorner,nCell)
            if (dphi > pi) then
               dphi = dphi - pi2
            else if (dphi < -pi) then
               dphi = dphi + pi2
            endif
            dphi = 0.5_SCRIP_r8*dphi

            Area(nCell) = Area(nCell) + 
     &                    dphi*(sin(CornerLat(   nCorner,nCell)) +
     &                          sin(CornerLat(nextCorner,nCell)))
         end do

      end do

!-----------------------------------------------------------------------
!EOC

      end subroutine SCRIP_GridComputeArea

!***********************************************************************

      end module scrip_grids

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

