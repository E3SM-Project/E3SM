! Copyright (c) 2016,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module mpas_init_atm_gwd

   use mpas_framework
   use mpas_timekeeping

   public :: compute_gwd_fields

   private

   real (kind=RKIND), parameter :: Re = 6371229.0_RKIND        ! Earth radius in MPAS-Atmosphere
   real (kind=RKIND), parameter :: Pi = 2.0_RKIND * asin(1.0_RKIND)
   real (kind=RKIND), parameter :: rad2deg = 180.0_RKIND / Pi

   integer, parameter :: topo_x = 43200            ! x-dimension of global 30-arc-second topography array
   integer, parameter :: topo_y = 21600            ! y-dimension of global 30-arc-second topography array
   real (kind=RKIND), parameter :: pts_per_degree = real(topo_x,RKIND) / 360.0_RKIND
   real (kind=RKIND), parameter :: start_lat = -90.0_RKIND
   real (kind=RKIND), parameter :: start_lon = -180.0_RKIND

   ! Nominal delta-x (in meters) for sub-grid topography cells
   real (kind=RKIND), parameter ::  sg_delta = 2.0 * Pi * Re / (360.0_RKIND * real(pts_per_degree,RKIND))

   real (kind=R4KIND), dimension(:,:), pointer :: topo ! Global 30-arc-second topography
   real (kind=RKIND), dimension(:,:), pointer :: box   ! Subset of topography covering a grid cell
   real (kind=RKIND), dimension(:,:), pointer :: dxm   ! Size (meters) in zonal direction of a grid cell
   real (kind=RKIND) :: box_mean                       ! Mean value of topography in box
   integer :: nx, ny                                   ! Dimensions of box covering grid cell


   contains


   !***********************************************************************
   !
   !  function compute_gwd_fields
   !
   !> \brief   Main routine for computing GWDO fields on an MPAS mesh
   !> \author  Michael Duda
   !> \date    29 May 2015
   !> \details 
   !>  This is the main routine for computing GWDO statistics on an MPAS mesh.
   !>  Currently computed fields are:
   !>  var2d
   !>  con
   !>  ol{1,2,3,4}
   !>  oa{1,2,3,4}
   !
   !-----------------------------------------------------------------------
   function compute_gwd_fields(domain) result(iErr)
   
      use mpas_derived_types
      use mpas_kind_types
      use mpas_timer
      use mpas_stream_manager
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      integer :: iErr

      type (mpas_pool_type), pointer :: mesh, state
      integer :: iCell, i
      real (kind=RKIND) :: dc
      integer, pointer :: nCells
      integer, pointer :: nEdges
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnCell
      logical :: onUnitSphere
      real (kind=RKIND), pointer :: sphere_radius
      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell, dcEdge
      real (kind=RKIND), dimension(:), pointer :: var2d, con, oa1, oa2, oa3, oa4, ol1, ol2, ol3, ol4
      real (kind=RKIND), dimension(:), pointer :: elvmax, htheta, hgamma, hsigma
      character(len=StrKIND), pointer :: config_geog_data_path


      allocate(topo(topo_x,topo_y))

      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', mesh)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', state)

      call mpas_pool_get_config(mesh, 'sphere_radius', sphere_radius)
      call mpas_pool_get_config(domain % configs, 'config_geog_data_path', config_geog_data_path)

      !
      ! Retrieve pointers to arrays holding the latitudes and longitudes of
      ! cells, and arrays that will hold the computed GWDO statistics
      !
      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call mpas_pool_get_array(mesh, 'latCell', latCell)
      call mpas_pool_get_array(mesh, 'lonCell', lonCell)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(mesh, 'dcEdge', dcEdge)
      call mpas_pool_get_array(mesh, 'var2d', var2d)
      call mpas_pool_get_array(mesh, 'con', con)
      call mpas_pool_get_array(mesh, 'ol1', ol1)
      call mpas_pool_get_array(mesh, 'ol2', ol2)
      call mpas_pool_get_array(mesh, 'ol3', ol3)
      call mpas_pool_get_array(mesh, 'ol4', ol4)
      call mpas_pool_get_array(mesh, 'oa1', oa1)
      call mpas_pool_get_array(mesh, 'oa2', oa2)
      call mpas_pool_get_array(mesh, 'oa3', oa3)
      call mpas_pool_get_array(mesh, 'oa4', oa4)
!      call mpas_pool_get_array(mesh, 'elvmax', elvmax)
!      call mpas_pool_get_array(mesh, 'theta', htheta)
!      call mpas_pool_get_array(mesh, 'gamma', hgamma)
!      call mpas_pool_get_array(mesh, 'sigma', hsigma)

      iErr = read_global_30s_topo(config_geog_data_path)
      if (iErr /= 0) then
         write(stderrUnit,*) 'Error reading global 30-arc-sec topography for GWD statistics'
         return
      end if

      !
      ! It is possible that this code is called before the mesh fields have been scaled
      ! up to "Earth-sized". Because we need "Earth" distances to cut out bounding
      ! boxes from topography, we try here to detect whether we are on an unscaled
      ! unit sphere or not: if the maximum dcEdge value is less than 1.0, assume this
      ! is the case.
      !
      if (maxval(dcEdge(1:nEdges)) < 1.0_RKIND) then
         write(stderrUnit,*) 'Computing GWD statistics on a unit sphere'
         onUnitSphere = .true.
      else
         onUnitSphere = .false.
      end if

      !
      ! Main loop to compute each of the GWDO fields for every horizontal
      ! grid cell in the mesh.
      !
      do iCell=1,nCells

         !
         ! First, get an estimate of the mean diameter (in meters) of the grid 
         ! cell by averaging the distances to each of the neighboring cells
         !
         dc = 0.0
         do i=1,nEdgesOnCell(iCell)
            dc = dc + dcEdge(edgesOnCell(i,iCell))
         end do
         dc = dc / real(nEdgesOnCell(iCell),RKIND)
         if (onUnitSphere) then
            dc = dc * sphere_radius
         end if

         !
         ! Cut out a rectangular piece of the global 30-arc-second topography
         ! data that is centered at the lat/lon of the current cell being
         ! processed and that is just large enough to cover the cell. The
         ! rectangular array of topography data is stored in the module
         ! variable 'box', and the dimensions of this array are given by the
         ! module variables 'nx' and 'ny'. The get_box() routine also
         ! computes the mean elevation in the array and stores that value in
         ! the module variable 'box_mean'.
         !
         call get_box(latCell(iCell)*rad2deg, lonCell(iCell)*rad2deg, dc)

         !
         ! With a box of 30-arc-second data for the current grid cell, call
         ! subroutines to compute each sub-grid orography statistic
         !
         var2d(iCell) = get_var()
         con(iCell) = get_con()
         oa1(iCell) = get_oa1()
         oa2(iCell) = get_oa2()
         oa3(iCell) = get_oa3()
         oa4(iCell) = get_oa4()
         ol1(iCell) = get_ol1()
         ol2(iCell) = get_ol2()
         ol3(iCell) = get_ol3()
         ol4(iCell) = get_ol4()
!         elvmax(iCell) = get_elvmax()
!         htheta(iCell) = get_htheta()
!         hgamma(iCell) = get_hgamma()
!         hsigma(iCell) = get_hsigma()
      end do

      deallocate(topo)

      iErr = 0

   end function compute_gwd_fields


   !***********************************************************************
   !
   !  function read_global_30s_topo
   !
   !> \brief   Reads global 30-arc-second topography into 'topo' module variable
   !> \author  Michael Duda
   !> \date    31 October 2016
   !> \details 
   !>  This subroutine reads the global 30-arc-second topography from the subdirectory 
   !>  'topo_30s' of the path provided as an argument.
   !
   !-----------------------------------------------------------------------
   function read_global_30s_topo(path) result(iErr)

      implicit none

      character(len=*), intent(in) :: path

      integer :: iErr

      integer, parameter :: tile_x = 1200       ! x-dimension of each tile of global 30-arc-second topography
      integer, parameter :: tile_y = 1200       ! y-dimension of each tile of global 30-arc-second topography
      integer, parameter :: tile_bdr = 3        ! number of layers of border/halo points surrounding each tile

      integer :: istatus
      integer :: ix, iy
      integer :: isigned, endian, wordsize, nx, ny, nz
      real (kind=R4KIND) :: scalefactor
      real (kind=R4KIND), dimension(:,:,:), allocatable :: tile
      character(len=StrKIND) :: filename

      allocate(tile(tile_x+2*tile_bdr,tile_y+2*tile_bdr,1))

      isigned  = 1
      endian   = 0
      wordsize = 2
      scalefactor = 1.0
      nx = tile_x + 2*tile_bdr
      ny = tile_y + 2*tile_bdr
      nz = 1

      do iy=1,topo_y,tile_y
      do ix=1,topo_x,tile_x
         write(filename,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(path)//'/topo_30s/', ix, '-', (ix+tile_x-1), '.', &
                                                                                      iy, '-', (iy+tile_y-1)
         call read_geogrid(filename, len_trim(filename), tile, nx, ny, nz, isigned, endian, &
                           scalefactor, wordsize, istatus)
         if (istatus /= 0) then
            write(stderrUnit,*) 'Error reading topography tile '//trim(filename)
            iErr = 1
            return
         end if

         topo(ix:(ix+tile_x-1),iy:(iy+tile_y-1)) = tile((tile_bdr+1):(tile_x+tile_bdr),(tile_bdr+1):(tile_y+tile_bdr),1)

      end do
      end do

      deallocate(tile)

      iErr = 0

   end function read_global_30s_topo


   !***********************************************************************
   !
   !  subroutine get_box
   !
   !> \brief   Cuts out a rectangular box of data centered at a given (lat,lon)
   !> \author  Michael Duda
   !> \date    29 May 2015
   !> \details 
   !>  This subroutine extracts a rectangular sub-array of the 30-arc-second
   !>  global topography dataset, stored in the module variable 'topo'; the
   !>  sub-array will be centered at the (lat,lon) specified in the argument
   !>  list, and will have a width and height large enough to span 'dx' meters.
   !>  The extracted sub-array is stored in the module variable 'box', and the
   !>  dimensions of this sub-array are stored in the module variables 'nx' and
   !>  'ny'.
   !>  Since the mean value of the terrain in a grid cell is needed by many of
   !>  the GWDO statistics computations, this mean value is also computed by
   !>  this subroutine and stored in the module variable 'box_mean'.
   !
   !-----------------------------------------------------------------------
   subroutine get_box(lat, lon, dx)

      implicit none

      real (kind=RKIND), intent(in) :: lat, lon, dx

      integer :: i, j, ii, jj, ic, jc
      real (kind=RKIND) :: sg_lat

      !
      ! Get number of points to extract in the zonal direction
      !
      if (cos(lat/rad2deg) > (2.0 * pts_per_degree * dx * 180.0) / (real(topo_x,RKIND) * Pi * Re)) then
         nx = ceiling((180.0 * dx * pts_per_degree) / (Pi * Re * cos(lat/rad2deg)))
      else
         nx = topo_x / 2
      end if

      !
      ! Get number of points to extract in the meridional direction
      !
      ny = ceiling((180.0 * dx * pts_per_degree) / (Pi * Re))

      !
      ! Find coordinates in global topography array of the box center
      !
      ic = nint((lon - start_lon) * pts_per_degree) + 1
      jc = nint((lat - start_lat) * pts_per_degree) + 1

      if (ic <= 0) ic = ic + topo_x
      if (ic > topo_x) ic = ic - topo_x


      if (associated(box)) deallocate(box)
      allocate(box(nx,ny))

      if (associated(dxm)) deallocate(dxm)
      allocate(dxm(nx,ny))

      !
      ! Extract sub-array (box) from global array; must properly account for 
      ! the periodicity in the longitude coordinate, as well as the poles
      !
      box_mean = 0.0
      do j=1,ny
      do i=1,nx

         ii = i - nx/2 + ic
         jj = j - ny/2 + jc

         if (jj <= 0) then
            jj = -jj + 1
            ii = ii + topo_y
         end if
         if (jj > topo_y) then
            jj = topo_y - (jj - topo_y - 1)
            ii = ii + topo_y
         end if
         do while (ii <= 0)
            ii = ii + topo_x
         end do
         do while (ii > topo_x)
            ii = ii - topo_x
         end do
 
         box(i,j) = topo(ii,jj)
         sg_lat = (start_lat + (real(jj-1,RKIND) + 0.5) / pts_per_degree) / rad2deg  ! Add 0.5 for cell center
         dxm(i,j) = sg_delta * cos(sg_lat)
         box_mean = box_mean + box(i,j)

      end do
      end do

      !
      ! Compute mean topography in the extracted box
      !
      box_mean = box_mean / real(nx*ny, RKIND)

   end subroutine get_box


   !***********************************************************************
   !
   !  function get_var
   !
   !> \brief   Computes standard deviation of sub-grid-scale terrain
   !> \author  Michael Duda
   !> \date    29 May 2015
   !> \details 
   !
   !-----------------------------------------------------------------------
   real (kind=RKIND) function get_var()

      implicit none

      integer :: i, j
      real (kind=RKIND) :: s2

      s2 = 0.0

      do j=1,ny
         do i=1,nx
            s2 = s2 + (box(i,j) - box_mean)**2
         end do
      end do

      get_var = sqrt(s2 / real(nx*ny,RKIND))

   end function get_var


   !***********************************************************************
   !
   !  function get_con
   !
   !> \brief   Computes orographic convexityof sub-grid-scale terrain
   !> \author  Michael Duda
   !> \date    29 May 2015
   !> \details 
   !
   !-----------------------------------------------------------------------
   real (kind=RKIND) function get_con()

      implicit none

      integer :: i, j
      real (kind=RKIND) :: s2, s4, var

      s2 = 0.0
      s4 = 0.0

      do j=1,ny
         do i=1,nx
            s2 = s2 + (box(i,j) - box_mean)**2
            s4 = s4 + (box(i,j) - box_mean)**4
         end do
      end do

      var = s2 / real(nx*ny,RKIND)
 
      if (abs(var) < 1.0e-5) then
         get_con = 0.0
      else
         get_con = s4 / (var**2 * real(nx*ny,RKIND))
      end if

   end function get_con


   !***********************************************************************
   !
   !  function get_oa1
   !
   !> \brief   Computes orographic asymmetry in the West direction
   !> \author  Michael Duda
   !> \date    29 May 2015
   !> \details 
   !>  This function computes the sub-grid orographic asymmetry following 
   !>  the comment from N. Wood in the footnote of Kim and Doyle (QRJMS, 2005).
   !
   !-----------------------------------------------------------------------
   real (kind=RKIND) function get_oa1()

      implicit none

      integer :: i, j
      integer :: nu, nd

      nu = 0
      nd = 0
      do j=1,ny
         do i=1,nx/2
            if (box(i,j) > box_mean) nu = nu + 1
         end do
         do i=nx/2+1,nx
            if (box(i,j) > box_mean) nd = nd + 1
         end do
      end do

      if (nu + nd > 0) then
         get_oa1 = real((nu - nd),RKIND) / real((nu + nd),RKIND)
      else
         get_oa1 = 0.0
      end if

   end function get_oa1


   !***********************************************************************
   !
   !  function get_oa2
   !
   !> \brief   Computes orographic asymmetry in the South direction
   !> \author  Michael Duda
   !> \date    29 May 2015
   !> \details 
   !>  This function computes the sub-grid orographic asymmetry following 
   !>  the comment from N. Wood in the footnote of Kim and Doyle (QRJMS, 2005).
   !
   !-----------------------------------------------------------------------
   real (kind=RKIND) function get_oa2()

      implicit none

      integer :: i, j
      integer :: nu, nd

      nu = 0
      nd = 0
      do j=1,ny/2
         do i=1,nx
            if (box(i,j) > box_mean) nu = nu + 1
         end do
      end do
      do j=ny/2+1,ny
         do i=1,nx
            if (box(i,j) > box_mean) nd = nd + 1
         end do
      end do

      if (nu + nd > 0) then
         get_oa2 = real((nu - nd),RKIND) / real((nu + nd),RKIND)
      else
         get_oa2 = 0.0
      end if

   end function get_oa2


   !***********************************************************************
   !
   !  function get_oa3
   !
   !> \brief   Computes orographic asymmetry in the South-West direction
   !> \author  Michael Duda
   !> \date    29 May 2015
   !> \details 
   !>  This function computes the sub-grid orographic asymmetry following 
   !>  the comment from N. Wood in the footnote of Kim and Doyle (QRJMS, 2005).
   !
   !-----------------------------------------------------------------------
   real (kind=RKIND) function get_oa3()

      implicit none

      integer :: i, j
      integer :: nu, nd
      real (kind=RKIND) :: ratio

      nu = 0
      nd = 0
      ratio = real(ny,RKIND)/real(nx,RKIND)
      do j=1,ny
         do i=1,nx
            if (nint(real(i,RKIND) * ratio) < (ny - j)) then
               if (box(i,j) > box_mean) nu = nu + 1
            else
               if (box(i,j) > box_mean) nd = nd + 1
            end if
         end do
      end do

      if (nu + nd > 0) then
         get_oa3 = real((nu - nd),RKIND) / real((nu + nd),RKIND)
      else
         get_oa3 = 0.0
      end if

   end function get_oa3


   !***********************************************************************
   !
   !  function get_oa4
   !
   !> \brief   Computes orographic asymmetry in the North-West direction
   !> \author  Michael Duda
   !> \date    29 May 2015
   !> \details 
   !>  This function computes the sub-grid orographic asymmetry following 
   !>  the comment from N. Wood in the footnote of Kim and Doyle (QRJMS, 2005).
   !
   !-----------------------------------------------------------------------
   real (kind=RKIND) function get_oa4()

      implicit none

      integer :: i, j
      integer :: nu, nd
      real (kind=RKIND) :: ratio

      nu = 0
      nd = 0
      ratio = real(ny,RKIND)/real(nx,RKIND)
      do j=1,ny
         do i=1,nx
            if (nint(real(i,RKIND) * ratio) < j) then
               if (box(i,j) > box_mean) nu = nu + 1
            else
               if (box(i,j) > box_mean) nd = nd + 1
            end if
         end do
      end do

      if (nu + nd > 0) then
         get_oa4 = real((nu - nd),RKIND) / real((nu + nd),RKIND)
      else
         get_oa4 = 0.0
      end if

   end function get_oa4


   !***********************************************************************
   !
   !  function get_ol1
   !
   !> \brief   Computes orographic effective length for Westerly flow
   !> \author  Michael Duda
   !> \date    29 May 2015
   !> \details 
   !
   !-----------------------------------------------------------------------
   real (kind=RKIND) function get_ol1()

      implicit none

      integer :: i, j
      integer :: nw
      integer :: nt

      nw = 0
      nt = 0

      do j=ny/4,3*ny/4
         do i=1,nx
            if (box(i,j) > box_mean) nw = nw + 1
            nt = nt + 1
         end do
      end do

      get_ol1 = real(nw,RKIND) / real(nt,RKIND)

   end function get_ol1


   !***********************************************************************
   !
   !  function get_ol2
   !
   !> \brief   Computes orographic effective length for Southerly flow
   !> \author  Michael Duda
   !> \date    29 May 2015
   !> \details 
   !
   !-----------------------------------------------------------------------
   real (kind=RKIND) function get_ol2()

      implicit none

      integer :: i, j
      integer :: nw
      integer :: nt

      nw = 0
      nt = 0

      do j=1,ny
         do i=nx/4,3*nx/4
            if (box(i,j) > box_mean) nw = nw + 1
            nt = nt + 1
         end do
      end do

      get_ol2 = real(nw,RKIND) / real(nt,RKIND)

   end function get_ol2


   !***********************************************************************
   !
   !  function get_ol3
   !
   !> \brief   Computes orographic effective length for South-Westerly flow
   !> \author  Michael Duda
   !> \date    29 May 2015
   !> \details 
   !
   !-----------------------------------------------------------------------
   real (kind=RKIND) function get_ol3()

      implicit none

      integer :: i, j
      integer :: nw
      integer :: nt

      nw = 0
      nt = 0

      do j=1,ny/2
         do i=1,nx/2
            if (box(i,j) > box_mean) nw = nw + 1
            nt = nt + 1
         end do
      end do
      do j=ny/2+1,ny
         do i=nx/2+1,nx
            if (box(i,j) > box_mean) nw = nw + 1
            nt = nt + 1
         end do
      end do

      get_ol3 = real(nw,RKIND) / real(nt,RKIND)

   end function get_ol3


   !***********************************************************************
   !
   !  function get_ol4
   !
   !> \brief   Computes orographic effective length for North-Westerly flow
   !> \author  Michael Duda
   !> \date    29 May 2015
   !> \details 
   !
   !-----------------------------------------------------------------------
   real (kind=RKIND) function get_ol4()

      implicit none

      integer :: i, j
      integer :: nw
      integer :: nt

      nw = 0
      nt = 0

      do j=ny/2+1,ny
         do i=1,nx/2
            if (box(i,j) > box_mean) nw = nw + 1
            nt = nt + 1
         end do
      end do
      do j=1,ny/2
         do i=nx/2+1,nx
            if (box(i,j) > box_mean) nw = nw + 1
            nt = nt + 1
         end do
      end do

      get_ol4 = real(nw,RKIND) / real(nt,RKIND)

   end function get_ol4


   !***********************************************************************
   !
   !  function get_elvmax
   !
   !> \brief   Computes maximum subgrid orography height
   !> \author  Michael Duda
   !> \date    20 December 2015
   !> \details 
   !
   !-----------------------------------------------------------------------
   real (kind=RKIND) function get_elvmax()

      implicit none

      integer :: i, j

      get_elvmax = box(1,1)

      do j=1,ny
         do i=1,nx
            if (box(i,j) > get_elvmax) then
               get_elvmax = box(i,j)
            end if
         end do
      end do

   end function get_elvmax


   !***********************************************************************
   !
   !  function get_htheta
   !
   !> \brief   Computes angle of principle axis of the gradient correlation tensor
   !> \author  Michael Duda
   !> \date    20 December 2015
   !> \details Computation following Lott and Miller (QJRMS 1997)
   !
   !-----------------------------------------------------------------------
   real (kind=RKIND) function get_htheta()

      implicit none

      integer :: i, j
      real (kind=RKIND) :: dx, dy
      real (kind=RKIND) :: xfp, yfp
      real (kind=RKIND) :: hx2, hy2, hxy
      real (kind=RKIND) :: hk, hl

      hx2 = 0.0
      hy2 = 0.0
      hxy = 0.0

      do j=2,ny-1
      do i=2,nx-1
         dx = dxm(i,j)
         dy = sg_delta
         xfp = (box(i+1,j) - box(i-1,j)) / (2.0 * dx)
         yfp = (box(i,j+1) - box(i,j-1)) / (2.0 * dy)
         hx2 = hx2 + xfp * xfp
         hy2 = hy2 + yfp * yfp
         hxy = hxy + xfp * yfp
      end do
      end do

      hx2 = hx2 / real((nx-2)*(ny-2),RKIND)
      hy2 = hy2 / real((nx-2)*(ny-2),RKIND)
      hxy = hxy / real((nx-2)*(ny-2),RKIND)

      hk = 0.5 * (hx2 + hy2)
      hl = 0.5 * (hx2 - hy2)

      get_htheta = 0.5 * atan2(hxy, hl)

   end function get_htheta


   !***********************************************************************
   !
   !  function get_hgamma
   !
   !> \brief   Computes anisotropy of subgrid orography
   !> \author  Michael Duda
   !> \date    20 December 2015
   !> \details Computation following Lott and Miller (QJRMS 1997)
   !
   !-----------------------------------------------------------------------
   real (kind=RKIND) function get_hgamma()

      implicit none

      integer :: i, j
      real (kind=RKIND) :: dx, dy
      real (kind=RKIND) :: xfp, yfp
      real (kind=RKIND) :: hx2, hy2, hxy
      real (kind=RKIND) :: hk, hl, hlp

      hx2 = 0.0
      hy2 = 0.0
      hxy = 0.0

      do j=2,ny-1
      do i=2,nx-1
         dx = dxm(i,j)
         dy = sg_delta
         xfp = (box(i+1,j) - box(i-1,j)) / (2.0 * dx)
         yfp = (box(i,j+1) - box(i,j-1)) / (2.0 * dy)
         hx2 = hx2 + xfp * xfp
         hy2 = hy2 + yfp * yfp
         hxy = hxy + xfp * yfp
      end do
      end do

      hx2 = hx2 / real((nx-2)*(ny-2),RKIND)
      hy2 = hy2 / real((nx-2)*(ny-2),RKIND)
      hxy = hxy / real((nx-2)*(ny-2),RKIND)

      hk = 0.5 * (hx2 + hy2)
      hl = 0.5 * (hx2 - hy2)
      hlp = sqrt(hl*hl + hxy*hxy)

      if ((hk + hlp) > 0.0 .and. (hk - hlp) >= 0.0) then
         get_hgamma = sqrt((hk - hlp) / (hk + hlp))
      else
         get_hgamma = 0.0
      end if

   end function get_hgamma


   !***********************************************************************
   !
   !  function get_hsigma
   !
   !> \brief   Computes mean slope of subgrid orography
   !> \author  Michael Duda
   !> \date    20 December 2015
   !> \details Computation following Lott and Miller (QJRMS 1997)
   !
   !-----------------------------------------------------------------------
   real (kind=RKIND) function get_hsigma()

      implicit none

      integer :: i, j
      real (kind=RKIND) :: dx, dy
      real (kind=RKIND) :: xfp, yfp
      real (kind=RKIND) :: hx2, hy2, hxy
      real (kind=RKIND) :: hk, hl, hlp

      hx2 = 0.0
      hy2 = 0.0
      hxy = 0.0

      do j=2,ny-1
      do i=2,nx-1
         dx = dxm(i,j)
         dy = sg_delta
         xfp = (box(i+1,j) - box(i-1,j)) / (2.0 * dx)
         yfp = (box(i,j+1) - box(i,j-1)) / (2.0 * dy)
         hx2 = hx2 + xfp * xfp
         hy2 = hy2 + yfp * yfp
         hxy = hxy + xfp * yfp
      end do
      end do

      hx2 = hx2 / real((nx-2)*(ny-2),RKIND)
      hy2 = hy2 / real((nx-2)*(ny-2),RKIND)
      hxy = hxy / real((nx-2)*(ny-2),RKIND)

      hk = 0.5 * (hx2 + hy2)
      hl = 0.5 * (hx2 - hy2)
      hlp = sqrt(hl*hl + hxy*hxy)

      get_hsigma = sqrt(hk + hlp)

   end function get_hsigma

end module mpas_init_atm_gwd
