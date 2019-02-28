!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     this module contains necessary routines for computing addresses
!     and weights for a conservative interpolation  between any two 
!     grids on a sphere.  the weights are computed by performing line 
!     integrals around all overlap regions of the two grids.  see 
!     Dukowicz and Kodis, SIAM J. Sci. Stat. Comput. 8, 305 (1987) and
!     Jones, P.W. Monthly Weather Review (submitted).
!
!-----------------------------------------------------------------------
!
!     CVS:$Id: remap_conserv.f,v 1.10 2001/08/21 21:05:13 pwjones Exp $
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
!     This code has been modified from the version available from 
!     Los Alamos National Laboratory, for the purpose of running it
!     within WW3. Primary modifications:
!     - renamed many variables to be unique across the code
!     - "save" variables moved from subroutine to module so that
!        we can "clear" them later.
!     - print statements added.
!     - phi_or_theta = 2 instead of phi_or_theta = 1 (important!)
!
!***********************************************************************
!  Modifications introduced by M. Dutour (MD) for 
!  running with WAVEWATCH III  ... see below
!
!     
!     BE CAREFUL ABOUT EXPLICIT INITIALIZATION OF VARIABLES IN
!     MULTI-THREADED VERSION OF THE CODE - INLINE INITIALIZATION OF
!     A VARIABLE IN FORTRAN 90/95 MAKES THE VARIABLE IMPLICITLY STATIC.
!     OPENMP FORCES _ALL_ FORTRAN IMPLEMENTATIONS TO MAKE THE VARIABLE
!     STATIC (OR OF THE TYPE SAVE) IF IT IS INITIALIZED IN THE 
!     DECLARATION LINE
!
!

      module scrip_remap_conservative

!-----------------------------------------------------------------------

      use SCRIP_KindsMod ! defines common data types
      use SCRIP_constants    ! defines common constants
      use scrip_timers       ! module for timing
      use scrip_grids        ! module containing grid information
      use scrip_remap_vars   ! module containing remap information
      use omp_lib

      implicit none

      integer (SCRIP_i4) :: nthreads=1  ! Number of parallel threads

!............variables that needed to be moved from "local level" to
!............ "module level" in order that we can clear them later.
!............These are all local variables that had the "save" attribute
!............in the standard version of SCRIP

      integer (SCRIP_i4), save :: 
     &     avoid_pole_count = 0  ! count attempts to avoid pole

      real (SCRIP_r8), save :: 
     &     avoid_pole_offset = tiny  ! endpoint offset to avoid pole

      integer (SCRIP_i4), dimension(:,:), allocatable, save ::
     &        link_add1,  ! min,max link add to restrict search
     &        link_add2   ! min,max link add to restrict search

      logical (SCRIP_logical), save :: 
     &        first_call_store_link_cnsrv = .true.

      logical (SCRIP_logical), save ::
     &     first_call_locate_segstart= .true.

      integer (SCRIP_i4), save :: 
     &     last_cell_locate_segstart=0,     ! save the search parameters
     &     last_cell_grid_num_locate_segstart=0,   ! if unchanged, reuse
                                                   ! search lists
     &     last_srch_grid_num_locate_segstart=0

      integer (SCRIP_i4), save ::
     &     num_srch_cells_locate_segstart=0,
     &     srch_corners_locate_segstart      ! number of corners for 
                                             ! each cell

      integer (SCRIP_i4), dimension(:), allocatable, save :: 
     &        srch_add_locate_segstart       ! global address of cells
                                             ! in srch arrays

      real (SCRIP_r8), dimension(:,:), allocatable, save ::
     &     srch_corner_lat_locate_segstart,  ! lat of each corner of
                                             ! srch cells
     &     srch_corner_lon_locate_segstart   ! lon of each corner of
                                             ! srch cells

      real(SCRIP_r8), dimension(:), allocatable, save ::
     &     srch_center_lat_locate_segstart,! lat of center of srch cells
     &     srch_center_lon_locate_segstart ! lon of center of srch cells
     
      logical (SCRIP_logical), save ::
     &     first_call_locate_point= .true.

      integer (SCRIP_i4), save :: 
     &     last_cell_locate_point=0,       ! save the search parameters
     &     last_cell_grid_num_locate_point=0,     ! if unchanged, reuse
                                                  ! search lists
     &     last_srch_grid_num_locate_point=0

      integer (SCRIP_i4), save ::
     &     num_srch_cell_locate_points=0,
     &     srch_corners_locate_point  ! number of corners for each cell

      integer (SCRIP_i4), dimension(:), allocatable, save :: 
     &        srch_add_locate_point       ! global address of cells in
                                          ! srch arrays

      real (SCRIP_r8), dimension(:,:), allocatable, save ::
     &     srch_corner_lat_locate_point,  ! lat of each corner of srch
                                          ! cells
     &     srch_corner_lon_locate_point   ! lon of each corner of srch
                                          ! cells

      real (SCRIP_r8), dimension(:), allocatable, save ::
     &     srch_center_lat_locate_point,  ! lat of center of srch cells
     &     srch_center_lon_locate_point   ! lon of center of srch cells

      integer (SCRIP_i4), save ::
     &     num_srch_cells_loc_get_srch_cells,  ! Number of srch cells
                                               ! found
     &     srch_corners_loc_get_srch_cells     ! Number of corners for
                                               ! search cells

      integer (SCRIP_i4), dimension(:), allocatable, save ::
     &     srch_add_loc_get_srch_cells         ! Global addresses of
                                               ! search cells

      real (SCRIP_r8), dimension(:,:), allocatable, save ::
     &     srch_corner_lat_loc_get_srch_cells, 
     &     srch_corner_lon_loc_get_srch_cells
      
      real (SCRIP_r8), dimension(:), allocatable, save ::
     &     srch_center_lat_loc_get_srch_cells,
     &     srch_center_lon_loc_get_srch_cells

      integer (SCRIP_i4), save ::
     &     last_cell_add_get_srch_cells, 
     &     last_cell_grid_num_get_srch_cells, 
     &     last_srch_grid_num_get_srch_cells

      logical (SCRIP_logical), save ::
     &     first_call_get_srch_cells=.true.

      logical (SCRIP_logical), save ::
     &     first_call_find_adj_cell=.true.

      logical (SCRIP_logical), private :: is_master
           ! module's equivalent of "l_master"

      integer (SCRIP_i4), save :: 
     &     last_cell_find_adj_cell,
     &     last_cell_grid_num_find_adj_cell,
     &     num_srch_cells_find_adj_cell, 
     &     srch_corners_find_adj_cell

      integer (SCRIP_i4), dimension(:), allocatable, save :: 
     &     srch_add_find_adj_cell
      real (SCRIP_r8), dimension(:,:), allocatable, save ::
     &     srch_corner_lat_find_adj_cell, srch_corner_lon_find_adj_cell

      real (SCRIP_r8), dimension(:), allocatable, save ::
     &     srch_center_lat_find_adj_cell, srch_center_lon_find_adj_cell

C$OMP THREADPRIVATE(last_cell_grid_num_get_srch_cells,
C$OMP& last_srch_grid_num_get_srch_cells,
C$OMP& first_call_get_srch_cells,
C$OMP& last_cell_add_get_srch_cells,
C$OMP& num_srch_cells_loc_get_srch_cells,
C$OMP& srch_corners_loc_get_srch_cells,
C$OMP& srch_add_loc_get_srch_cells,
C$OMP& srch_corner_lat_loc_get_srch_cells,
C$OMP& srch_corner_lon_loc_get_srch_cells,
C$OMP& srch_center_lat_loc_get_srch_cells,
C$OMP& srch_center_lon_loc_get_srch_cells) 

C$OMP THREADPRIVATE(first_call_locate_segstart,
C$OMP& last_cell_locate_segstart,
C$OMP& last_cell_grid_num_locate_segstart,
C$OMP& last_srch_grid_num_locate_segstart,
C$OMP& num_srch_cells_locate_segstart,
C$OMP& srch_corners_locate_segstart,
C$OMP& srch_add_locate_segstart,
C$OMP& srch_corner_lat_locate_segstart,
C$OMP& srch_corner_lon_locate_segstart,
C$OMP& srch_center_lat_locate_segstart,
C$OMP& srch_center_lon_locate_segstart)

C$OMP THREADPRIVATE(first_call_locate_point,
C$OMP& last_cell_locate_point,
C$OMP& last_cell_grid_num_locate_point,
C$OMP& last_srch_grid_num_locate_point,
C$OMP& num_srch_cell_locate_points,
C$OMP& srch_add_locate_point,srch_corner_lat_locate_point,
C$OMP& srch_corner_lon_locate_point,
C$OMP& srch_center_lat_locate_point,
C$OMP& srch_center_lon_locate_point)

C$OMP THREADPRIVATE(first_call_find_adj_cell,
C$OMP& last_cell_find_adj_cell,
C$OMP& last_cell_grid_num_find_adj_cell,
C$OMP& num_srch_cells_find_adj_cell,
C$OMP& srch_corners_find_adj_cell,
C$OMP& srch_add_find_adj_cell,
C$OMP& srch_corner_lat_find_adj_cell,
C$OMP& srch_corner_lon_find_adj_cell,
C$OMP& srch_center_lat_find_adj_cell,
C$OMP& srch_center_lon_find_adj_cell)

!***********************************************************************

      contains

!***********************************************************************

      subroutine remap_conserv(l_master, l_test)

!-----------------------------------------------------------------------
!
!     this routine traces the perimeters of every grid cell on each
!     grid checking for intersections with the other grid and computing
!     line integrals for each subsegment.
!
!-----------------------------------------------------------------------

      logical(SCRIP_Logical), intent(in) :: l_master   ! Am I the master
                                                   ! processor (do I/O)?
      logical(SCRIP_Logical), intent(in) :: l_test     ! Whether to 
                                                    !include test output

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4), parameter :: 
     &     phi_or_theta = 2      ! integrate w.r.t. phi (1) or theta (2)
      

      integer (SCRIP_i4) :: 
     &     i, inext,            !
     &     n, nwgt, 
     &     grid1_add,           ! Current linear address for grid1 cell
     &     grid2_add,           ! Current linear address for grid2 cell
     &     grid_num,            ! Index (1,2) of grid that we are
                                ! processing
     &     opp_grid_num,        ! Index of opposite grid (2,1)
     &     maxrd_cell,          ! cell with the max. relative difference
                                ! in area
     &     progint              ! Intervals at which progress is to be
                                ! printed
     &     ,icount              ! for counting

      real (SCRIP_r8) ::
     &     norm_factor          ! factor for normalizing wts

      real (SCRIP_r8), dimension(6) :: 
     &     weights             ! Weights array

      real (SCRIP_r8) ::
     &     beglat, beglon,
     &     endlat, endlon,
     &     ave_reldiff,         ! Average rel. diff. in areas
     &     max_reldiff,         ! Maximum rel. diff in areas
     &     maxrd_area,          ! Computed area for cell with max rel
                                ! diff
     &     maxrd_true           ! True area for cell with max rel diff

      real (SCRIP_r8), dimension(:), allocatable ::
     &     reldiff,             ! Relative difference in computed
                                ! and true area
     &     ref_area             ! Area of cell as computed by direct 
                                ! integration around its boundaries

!      call OMP_SET_DYNAMIC(.FALSE.)

!-----------------------------------------------------------------------
!
!     integrate around each cell on grid1
!
!-----------------------------------------------------------------------

      is_master=l_master ! set module variable using subroutine input
                         ! argument variable. 
                         ! Use the former subsequently.

      if(is_master)print *,'grid1 sweep'

!NRL  Progress is slow when the other grid (grid 2) is large, so we use
!NRL    that. Really, it would be a better to do this with a timer...
      if (grid2_size >     500000) then
         progint =         1000
      elseif (grid2_size > 250000) then
         progint =         2000 
      elseif (grid2_size > 100000) then
         progint =         5000 
      else
         progint =         10000
      endif

      grid_num = 1
      opp_grid_num = 2

      call timer_start(1)

C$OMP PARALLEL DEFAULT(SHARED) PRIVATE(grid1_add) NUM_THREADS(nthreads)

C$OMP DO SCHEDULE(DYNAMIC) 

      do grid1_add = 1,grid1_size

         if (mod(grid1_add,progint) .eq. 0 .and. is_master) then
            print *, grid1_add,' of ',grid1_size,' cells processed ...'
         endif

         call cell_integrate(grid1_add, grid_num, phi_or_theta)

      end do                    ! do grid1_add=...

C$OMP END DO

C$OMP END PARALLEL

!-----------------------------------------------------------------------
!
!     integrate around each cell on grid2
!
!-----------------------------------------------------------------------

      if(is_master)print *,'grid2 sweep '

!NRL  Progress is slow when the other grid (grid 1) is large, so we use
!NRL    that.
      if (grid1_size >     500000) then
         progint =         1000
      elseif (grid1_size > 250000) then
         progint =         2000 
      elseif (grid1_size > 100000) then
         progint =         5000 
      else
         progint =         10000
      endif

      grid_num = 2
      opp_grid_num = 1

      call timer_start(2)

C$OMP PARALLEL DEFAULT(SHARED) PRIVATE(grid2_add) NUM_THREADS(nthreads)

C$OMP DO SCHEDULE(DYNAMIC)

      do grid2_add = 1,grid2_size

         if (mod(grid2_add,progint) .eq. 0 .and. is_master) then
            print *, grid2_add,' of ',grid2_size,' cells processed ...'
         endif

         call cell_integrate(grid2_add, grid_num, phi_or_theta)

      end do                    ! do grid2_add=...

C$OMP END DO

C$OMP END PARALLEL

      call timer_stop(2)

!-----------------------------------------------------------------------
!
!     correct for situations where N/S pole not explicitly included in
!     grid (i.e. as a grid corner point). if pole is missing from only
!     one grid, need to correct only the area and centroid of that 
!     grid.  if missing from both, do complete weight calculation.
!     This is necessary only when integrating w.r.t. phi (longitude)
!
!-----------------------------------------------------------------------

      if (phi_or_theta .eq. 1) then

         !*** North Pole
         weights(1) =  pi2
         weights(2) =  pi*pi
         weights(3) =  zero
         weights(4) =  pi2
         weights(5) =  pi*pi
         weights(6) =  zero
         
         if (grid1_npole_cell /=0) then
            grid1_area(grid1_npole_cell) = grid1_area(grid1_npole_cell) 
     &           + weights(1)
            grid1_centroid_lat(grid1_npole_cell) = 
     &           grid1_centroid_lat(grid1_npole_cell) + weights(2)
            grid1_centroid_lon(grid1_npole_cell) =
     &           grid1_centroid_lon(grid1_npole_cell) + weights(3)
         endif
         
         if (grid2_npole_cell /=0) then
            grid2_area(grid2_npole_cell) = grid2_area(grid2_npole_cell) 
     &           + weights(num_wts+1)
            grid2_centroid_lat(grid2_npole_cell) = 
     &           grid2_centroid_lat(grid2_npole_cell) + 
     &           weights(num_wts+2)
            grid2_centroid_lon(grid2_npole_cell) =
     &           grid2_centroid_lon(grid2_npole_cell) + 
     &           weights(num_wts+3)
         endif
         
         if (grid1_npole_cell /= 0 .and. grid2_npole_cell /=0) then
            call store_link_cnsrv(grid1_npole_cell, 
     &           grid2_npole_cell, weights)
            
            grid1_frac(grid1_npole_cell) = grid1_frac(grid1_npole_cell) 
     &           + weights(1)
            grid2_frac(grid2_npole_cell) = grid2_frac(grid2_npole_cell) 
     &           + weights(num_wts+1)
         endif

         
         !*** South Pole
         weights(1) =  pi2
         weights(2) = -pi*pi
         weights(3) =  zero
         weights(4) =  pi2
         weights(5) = -pi*pi
         weights(6) =  zero
         
         if (grid1_spole_cell /=0) then
            grid1_area(grid1_spole_cell) = grid1_area(grid1_spole_cell) 
     &           + weights(1)
            grid1_centroid_lat(grid1_spole_cell) = 
     &           grid1_centroid_lat(grid1_spole_cell) + weights(2)
            grid1_centroid_lon(grid1_spole_cell) =
     &           grid1_centroid_lon(grid1_spole_cell) + weights(3)
         endif
         
         if (grid2_spole_cell /=0) then
            grid2_area(grid2_spole_cell) = grid2_area(grid2_spole_cell) 
     &           + weights(num_wts+1)
            grid2_centroid_lat(grid2_spole_cell) = 
     &           grid2_centroid_lat(grid2_spole_cell) + 
     &           weights(num_wts+2)
            grid2_centroid_lon(grid2_spole_cell) =
     &           grid2_centroid_lon(grid2_spole_cell) + 
     &           weights(num_wts+3)
         endif

         if (grid1_spole_cell /= 0 .and. grid2_spole_cell /=0) then
            call store_link_cnsrv(grid1_spole_cell, 
     &           grid2_spole_cell, weights)
            
            grid1_frac(grid1_spole_cell) = grid1_frac(grid1_spole_cell) 
     &           + weights(1)
            grid2_frac(grid2_spole_cell) = grid2_frac(grid2_spole_cell) 
     &           + weights(num_wts+1)
         endif
      endif


      
      if(is_master)print *, 'Grid sweeps completed'
      

!-----------------------------------------------------------------------
!
!     finish centroid computation
!
!-----------------------------------------------------------------------

      call timer_start(3)

C$OMP PARALLEL
C$OMP WORKSHARE
      where (grid1_area /= zero)
        grid1_centroid_lat = grid1_centroid_lat/grid1_area
        grid1_centroid_lon = grid1_centroid_lon/grid1_area
      end where
C$OMP END WORKSHARE

C$OMP WORKSHARE
      where (grid2_area /= zero)
        grid2_centroid_lat = grid2_centroid_lat/grid2_area
        grid2_centroid_lon = grid2_centroid_lon/grid2_area
      end where
C$OMP END WORKSHARE
C$OMP END PARALLEL


!-----------------------------------------------------------------------
!
!     include centroids in weights and normalize using destination
!     area if requested
!
!-----------------------------------------------------------------------

C$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(nthreads)
C$OMP& PRIVATE(n,grid1_add,grid2_add,nwgt,weights,norm_factor)

C$OMP DO SCHEDULE(DYNAMIC)

      do n=1,num_links_map1
        grid1_add = grid1_add_map1(n)
        grid2_add = grid2_add_map1(n)
        do nwgt=1,num_wts
          weights(        nwgt) = wts_map1(nwgt,n)
          if (num_maps > 1) then
            weights(num_wts+nwgt) = wts_map2(nwgt,n)
          endif
        end do

        select case(norm_opt)
        case (norm_opt_dstarea)
          if (grid2_area(grid2_add) /= zero) then
            if (luse_grid2_area) then
              norm_factor = one/grid2_area_in(grid2_add)
            else
              norm_factor = one/grid2_area(grid2_add)
            endif
          else
            norm_factor = zero
          endif
        case (norm_opt_frcarea)
          if (grid2_frac(grid2_add) /= zero) then
            if (luse_grid2_area) then
              norm_factor = grid2_area(grid2_add)/
     &                     (grid2_frac(grid2_add)*
     &                      grid2_area_in(grid2_add))
            else
              norm_factor = one/grid2_frac(grid2_add)
            endif
          else
            norm_factor = zero
          endif
        case (norm_opt_none)
          norm_factor = one
        end select

        wts_map1(1,n) =  weights(1)*norm_factor
        wts_map1(2,n) = (weights(2) - weights(1)*
     &                              grid1_centroid_lat(grid1_add))*
     &                              norm_factor
        wts_map1(3,n) = (weights(3) - weights(1)*
     &                              grid1_centroid_lon(grid1_add))*
     &                              norm_factor

        if (num_maps > 1) then
          select case(norm_opt)
          case (norm_opt_dstarea)
            if (grid1_area(grid1_add) /= zero) then
              if (luse_grid1_area) then
                norm_factor = one/grid1_area_in(grid1_add)
              else
                norm_factor = one/grid1_area(grid1_add)
              endif
            else
              norm_factor = zero
            endif
          case (norm_opt_frcarea)
            if (grid1_frac(grid1_add) /= zero) then
              if (luse_grid1_area) then
                norm_factor = grid1_area(grid1_add)/
     &                       (grid1_frac(grid1_add)*
     &                        grid1_area_in(grid1_add))
              else
                norm_factor = one/grid1_frac(grid1_add)
              endif
            else
              norm_factor = zero
            endif
          case (norm_opt_none)
            norm_factor = one
          end select

          wts_map2(1,n) =  weights(num_wts+1)*norm_factor
          wts_map2(2,n) = (weights(num_wts+2) - weights(num_wts+1)*
     &                                grid2_centroid_lat(grid2_add))*
     &                                norm_factor
          wts_map2(3,n) = (weights(num_wts+3) - weights(num_wts+1)*
     &                                grid2_centroid_lon(grid2_add))*
     &                                norm_factor
        endif

      end do

C$OMP END DO

C$OMP END PARALLEL

      if(is_master)print *, 'Total number of links = ',num_links_map1

C$OMP PARALLEL
C$OMP WORKSHARE
      where (grid1_area /= zero) grid1_frac = grid1_frac/grid1_area
C$OMP END WORKSHARE
C$OMP WORKSHARE
      where (grid2_area /= zero) grid2_frac = grid2_frac/grid2_area
C$OMP END WORKSHARE
C$OMP END PARALLEL

      call timer_stop(3)

!-----------------------------------------------------------------------
!
!     perform some error checking on final weights
!
!-----------------------------------------------------------------------

      allocate(ref_area(grid1_size))
      allocate(reldiff(grid1_size))

C$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(nthreads)
C$OMP&  PRIVATE(n, i, inext, beglat, beglon, endlat, endlon, weights)
C$OMP DO SCHEDULE(DYNAMIC)

      do n=1,grid1_size
        if (grid1_area(n) < -.01 .and. is_master) then
          print *,'Grid 1 area error: ',n,grid1_area(n)
        endif
        if ((grid1_centroid_lat(n) < -pih-.01 .or.
     &      grid1_centroid_lat(n) >  pih+.01) .and. is_master) then
          print *,'Grid 1 centroid lat error: ',n,grid1_centroid_lat(n)
        endif

        ref_area(n) = 0.0
        do i = 1, grid1_corners
           inext = 1 + mod(i,grid1_corners)

           beglat = grid1_corner_lat(i,n)
           beglon = grid1_corner_lon(i,n)
           endlat = grid1_corner_lat(inext,n)
           endlon = grid1_corner_lon(inext,n)

           if ((phi_or_theta .eq. 1 .and. beglon .eq. endlon) .or.
     &          (phi_or_theta .eq. 2 .and. beglat .eq. endlat)) cycle

           call line_integral(phi_or_theta, weights, num_wts, beglon, 
     &          endlon, beglat, endlat, grid1_center_lat(n), 
     &          grid1_center_lon(n), grid1_center_lat(n),
     &          grid1_center_lon(n))

           ref_area(n) = ref_area(n) + weights(1)
        enddo
      enddo
C$OMP END DO
C$OMP END PARALLEL


!     Correct for polar cells

      if (phi_or_theta .eq. 1) then

         !*** North Pole
         weights(1) =  pi2

         if (grid1_npole_cell /=0) then
            ref_area(grid1_npole_cell) = ref_area(grid1_npole_cell) 
     &           + weights(1)
         endif
         
         !*** South Pole
         weights(1) =  pi2
         
         if (grid1_spole_cell /=0) then
            ref_area(grid1_spole_cell) = ref_area(grid1_spole_cell) 
     &           + weights(1)
         endif
         
      endif


      ave_reldiff = 0.0
      max_reldiff = -1.0

      do n = 1, grid1_size
         if(ref_area(n).gt.0.0)then ! added May 21 2013
            reldiff(n) = abs(ref_area(n)-grid1_area(n))/abs(ref_area(n))
         endif
         ave_reldiff = ave_reldiff + reldiff(n)
         if (reldiff(n) > max_reldiff) then
            max_reldiff = reldiff(n)
            maxrd_cell = n
            maxrd_area = grid1_area(n)
            maxrd_true = ref_area(n)
         endif
      end do
      
      ave_reldiff = ave_reldiff/grid1_size

      if(is_master.and.l_test)then
         print *
         print *
         print *,'Grid 1: Ave. rel. diff. in areas: ',
     &        ave_reldiff
         print *,'        rel. diff. = abs(area-refarea)/refarea'
         print *
         print *,'Grid 1: Max. rel. diff. in areas: ',
     &        max_reldiff
         print *, 'Max rel. diff. is in cell ',maxrd_cell
         print *, 'Computed Area: ', maxrd_area
         print *, 'Reference Area: ',maxrd_true 
         print *
      endif

      deallocate(ref_area, reldiff)



      allocate(ref_area(grid2_size))
      allocate(reldiff(grid2_size))

C$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(nthreads) 
C$OMP&   PRIVATE(n, i, inext, beglat, beglon, endlat, endlon, weights)
C$OMP DO SCHEDULE(DYNAMIC)

      do n=1,grid2_size
        if (grid2_area(n) < -.01 .and. is_master) then
          print *,'Grid 2 area error: ',n,grid2_area(n)
        endif
        if ((grid2_centroid_lat(n) < -pih-.01 .or.
     &      grid2_centroid_lat(n) >  pih+.01) .and. is_master) then
          print *,'Grid 2 centroid lat error: ',n,grid2_centroid_lat(n)
        endif

        ref_area(n) = 0.0
        do i = 1, grid2_corners
           inext = 1 + mod(i,grid2_corners)

           beglat = grid2_corner_lat(i,n)
           beglon = grid2_corner_lon(i,n)
           endlat = grid2_corner_lat(inext,n)
           endlon = grid2_corner_lon(inext,n)

           if ((phi_or_theta .eq. 1 .and. beglon .eq. endlon) .or.
     &          (phi_or_theta .eq. 2 .and. beglat .eq. endlat)) cycle

           call line_integral(phi_or_theta, weights, num_wts, beglon, 
     &          endlon, beglat, endlat, grid2_center_lat(n), 
     &          grid2_center_lon(n), grid2_center_lat(n),
     &          grid2_center_lon(n))

           ref_area(n) = ref_area(n) + weights(1)
        enddo
      enddo
C$OMP END DO
C$OMP END PARALLEL


!     Correct for polar cells

      if (phi_or_theta .eq. 1) then

         !*** North Pole
         weights(1) =  pi2
         
         if (grid2_npole_cell /=0) then
            ref_area(grid2_npole_cell) = ref_area(grid2_npole_cell) 
     &           + weights(1)
         endif
         
         !*** South Pole
         weights(1) =  pi2
         
         if (grid2_spole_cell /=0) then
            ref_area(grid2_spole_cell) = ref_area(grid2_spole_cell) 
     &           + weights(1)
         endif
         
      endif


      ave_reldiff = 0.0
      max_reldiff = -1.0

      do n = 1, grid2_size
         reldiff(n) = abs(ref_area(n)-grid2_area(n))/abs(ref_area(n))
         ave_reldiff = ave_reldiff + reldiff(n)
         if (reldiff(n) > max_reldiff) then
            max_reldiff = reldiff(n)
            maxrd_cell = n
            maxrd_area = grid2_area(n)
            maxrd_true = ref_area(n)
         endif
      end do

      ave_reldiff = ave_reldiff/grid2_size

      if(is_master.and.l_test)then
         print *
         print *,'Grid 2: Ave. rel. diff. in areas: ',
     &        ave_reldiff
         print *,'        rel. diff. = abs(area-refarea)/refarea'
         print *
         print *,'Grid 2: Max. rel. diff. in areas: ',
     &        max_reldiff
         print *, 'Max rel. diff. is in cell ',maxrd_cell
         print *, 'Computed Area: ', maxrd_area
         print *, 'Reference Area: ',maxrd_true 
         print *
      endif

      deallocate(ref_area,reldiff)

      if(is_master.and.l_test)then
        print *, 'Computed area = Area of cell computed by adding areas'
        print *, '                of intersection with other cells'
        print *, 'Reference area = Area of cell by direct integration'
        print *
      endif

      !***
      !*** In the following code, gridN_centroid_lat is being used to 
      !*** store running tallies of the cell areas - so it is a 
      !*** misnomer used to avoid allocation of a new variable
      !***

      grid1_centroid_lat = zero
      grid2_centroid_lat = zero
      icount=0

C$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(nthreads)
C$OMP&   PRIVATE(n,grid1_add,grid2_add,nwgt,weights)
C$OMP DO SCHEDULE(DYNAMIC)

      do n=1,num_links_map1
        grid1_add = grid1_add_map1(n)
        grid2_add = grid2_add_map1(n)

        do nwgt=1,num_wts
          weights(        nwgt) = wts_map1(nwgt,n)
          if (num_maps > 1) then
            weights(num_wts+nwgt) = wts_map2(nwgt,n)
          endif
        end do

! count warnings about weights that will be excluded
        if (grid2_frac(grid2_add).gt.frac_lowest .and. 
     &       grid2_frac(grid2_add).lt.frac_highest .and. is_master) then
           if ( (wts_map1(1,n) < wt_lowest) )then
              icount=icount+1
! print statements that were here have been moved to another routine...
           endif
           if (norm_opt /= norm_opt_none .and. wts_map1(1,n) > 
     &          wt_highest)then
              icount=icount+1
! print statements that were here have been moved to another routine...
           endif
        endif
C$OMP   CRITICAL
        grid2_centroid_lat(grid2_add) = 
     &  grid2_centroid_lat(grid2_add) + wts_map1(1,n)
C$OMP   END CRITICAL

        if (num_maps > 1) then
          if (wts_map2(1,n) < -.01 .and. is_master) then
            print *,'Map 2 weight < 0 ',grid1_add,grid2_add,
     &                                  wts_map2(1,n)
          endif
          if (norm_opt /= norm_opt_none .and. wts_map2(1,n) > 1.01
     &         .and. is_master) then
            print *,'Map 2 weight > 1 ',grid1_add,grid2_add,
     &                                  wts_map2(1,n)
          endif
C$OMP     CRITICAL
          grid1_centroid_lat(grid1_add) = 
     &    grid1_centroid_lat(grid1_add) + wts_map2(1,n)
C$OMP     END CRITICAL
        endif
      end do

C$OMP END DO
C$OMP END PARALLEL

      if(icount.gt.0.and.is_master)then
         print *,'We had problems in ',icount,' points.'
      endif
! stop condition was here...has been moved to another routine...

      !***
      !*** If grid1 has masks, links between some cells of grid1 and
      !*** grid2 do not exist even though they overlap. In such a case,
      !*** the following code will generate errors even though nothing
      !*** is wrong (grid1_centroid_lat or grid2_centroid_lat are never
      !*** updated in the above loop)
      !*** 

      do n=1,grid2_size
        select case(norm_opt)
        case (norm_opt_dstarea)
          norm_factor = grid2_frac(n)
        case (norm_opt_frcarea)
          norm_factor = one
        case (norm_opt_none)
          if (luse_grid2_area) then
            norm_factor = grid2_area_in(n)
          else
            norm_factor = grid2_area(n)
          endif
        end select
!       if (abs(grid2_centroid_lat(n)-norm_factor) > .01 
!    &     .and. is_master) then
!         print *,'Warning: sum of wts for map1 ',n,
!    &            grid2_centroid_lat(n),norm_factor
!       endif
!       write(501,*)n,grid2_centroid_lat(n)
      end do


      if (num_maps > 1) then
        do n=1,grid1_size
          select case(norm_opt)
          case (norm_opt_dstarea)
            norm_factor = grid1_frac(n)
          case (norm_opt_frcarea)
            norm_factor = one
          case (norm_opt_none)
            if (luse_grid1_area) then
              norm_factor = grid1_area_in(n)
            else
              norm_factor = grid1_area(n)
            endif
          end select
          if (abs(grid1_centroid_lat(n)-norm_factor) > .01
     &      .and. is_master) then
            print *,'Error: sum of wts for map2 ',n,
     &              grid1_centroid_lat(n),norm_factor
          endif
        end do
      endif
!-----------------------------------------------------------------------

      call timer_stop(4)

      if(is_master)print *, 'Finished Conservative Remapping'

      if(l_test)then
         call timer_print(1)
         call timer_print(2)
         call timer_print(3)
         call timer_print(4)
      endif

      end subroutine remap_conserv

!***********************************************************************



!***********************************************************************

      subroutine cellblock_integrate(ibegin, iend, grid_num, 
     &     phi_or_theta) 

      integer (SCRIP_i4) :: ibegin, iend, grid_num, phi_or_theta

      integer (SCRIP_i4) :: cell_add

      
      do cell_add = ibegin, iend

         call cell_integrate(cell_add, grid_num, phi_or_theta)

      enddo



      end subroutine cellblock_integrate

!***********************************************************************



!***********************************************************************

      subroutine cell_integrate(cell_add, grid_num, phi_or_theta)

!-----------------------------------------------------------------------
!
!     Integrate around cell while finding intersecting with opposite
!     grid cells and finding segments of cell boundary lying in cells
!     of opposite grid
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     Input variables
!
!-----------------------------------------------------------------------
      
      integer (SCRIP_i4) ::
     &     cell_add,            ! cell to be processed
     &     grid_num,            ! grid that the cell belongs to
     &     phi_or_theta         ! Integration var :
                                ! phi (lon) or theta (lat)


!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4), parameter ::
     &     max_subseg = 100     ! max number of subsegments per segment
                                ! to prevent infinite loop


      integer (SCRIP_i4) :: 
     &     i, inext,            !
     &     j, jnext,            ! generic counters
     &     ic, k, ns,           !
     &     n, next_n,           !
     &     nwgt, it,            !
     &     oppcell_add,         ! Cell from opposite grid we are 
                                ! intersecting
     &     opp_grid_num,        ! Index of opposite grid (2,1)
     &     min_add,             ! addresses for restricting search of
     &     max_add,             ! destination grid
     &     corner,              ! corner of cell that segment starts 
                                ! from
     &     next_corn,           ! corner of cell that segment ends on
     &     nseg,                ! number of segments to use to represent
                                ! edges near the pole                   
     &     num_subseg,          ! number of subsegments
     &     bedgeid1,            !
     &     bedgeid2,            ! ID of edge that a point is on
     &     bedgeid3,            !
     &     intedge,             ! ID of intersected edge
     &     last_add,            ! Address of last cell we were in
     &     next_add,            ! Address of next cell we will go into
     &     adj_add              ! Address of cell adjacent to current 
                                ! one

      logical (SCRIP_logical) :: 
     &     lcoinc,              ! Are segments coincident?
     &     lrevers,             ! Are we integrating segment in reverse?
     &     lboundary1, 
     &     lboundary2,          ! Is point is on cell boundary?
     &     lboundary3,          
     &     last_lboundary,      ! Is last point is on cell bdry?
     &     loutside,            ! Is point outside the grid?
     &     lthresh,             ! Has segment crossed threshold?
     &     srch_success,        ! Was search for segment start 
                                ! successful?
     &     intrsct_success,     ! Was intersection of segment with 
                                ! opposite grid successful?
     &     inpoly,              ! Is point is in polygon
     &     last_endpt_inpoly,   ! Was end point of last segment in cell
     &     last_endpt_onedge,   ! Was end point of last segment on edge
                                ! of cell
     &     lstuck,              ! Is the walk stuck inside a cell
     &     seg_outside,         ! Is segment completely outside the grid
     &     bndedge,             ! Is segment on the boundary of the grid
     &     search,              ! Do we have to search to locate point 
                                ! in grid
     &     inpolar,             ! Are we in the polar region?
     &     special_cell,        ! Is this a special cell 
                                ! (only 1 vtx at pole)
     &     L_exit_do            ! Do we need to escape from infinite 
                                ! loop? (NRL)

      real (SCRIP_r8) ::
     &     intrsct_lat,         ! lat of next intersection point
     &     intrsct_lon,         ! lon of next intersection point
     &     beglat, beglon,      ! start point of current sub seg
     &     endlat, endlon,      ! endpoint of current seg  
                                ! (chg to endseg?)
     &     endlat1, endlon1,    ! endpoint of current subseg
     &     norm_factor          ! factor for normalizing wts

      real (SCRIP_r8), dimension(2) :: 
     &     begseg               ! begin lat/lon for full segment

      real (SCRIP_r8), dimension(6) :: 
     &     weights,             ! local wgt array
     &     rev_weights          ! Weights for grid1 and grid2 flipped

      real (SCRIP_r8) ::
     &     vec1_lat, vec1_lon,  ! vectors, products
     &     vec2_lat, vec2_lon,  ! used in grid search
     &     vec1_len, dp,
     &     midlat, midlon,      ! Midpoint of segment
     &     tmplat, tmplon,
     &     srchpt_lat,          ! Search point (offset from seg. start)
     &     srchpt_lon,          
     &     offset, delta,       ! Offset and offset increase for search
     &     sinang2,             ! Square of sine of angle b/w two  
                                ! segments
     &     dist2,               ! Square of distance b/w two points
     &     fullseg_len2,        ! Square of full segment length
     &     partseg_len2,        ! Square of length of segment integrated
                                ! so far
     &     fullseg_dlat,        ! Lat diff of full segment endpoints
     &     fullseg_dlon,        ! Lon diff of full segment endpoints
     &     prevlon,
     &     nextlon,
     &     pole_lat,
     &     cell_center_lat,
     &     cell_center_lon,
     &     oppcell_center_lat,
     &     oppcell_center_lon

      real (SCRIP_r8), dimension(:), allocatable ::
     &     cell_corner_lat,        ! Local copies of cell coordinates
     &     cell_corner_lon,        ! May be augmented for computational
                                   ! reasons
     &     oppcell_corner_lat,
     &     oppcell_corner_lon

      integer (SCRIP_i4) ::
     &     ncorners,        ! Number of corners in local copy of cell
     &     ncorners_opp,    ! Number of corners in local copy of oppcell
     &     nalloc,          ! Allocation for the cell_corner_* array
     &     nalloc_opp       ! Allocation for the oppcell_corner_* array

      real (SCRIP_r8) ::
     &     tmpwt1, tmpwt2

      integer (SCRIP_i4) ::
     &     ncorners_at_pole,
     &     previdx,
     &     nextidx


      if (grid_num .eq. 1) then

         !***
         !*** Set up a local copy of the cell with room to add
         !*** degenerate edges
         !***

         ncorners = grid1_corners
         nalloc = ncorners+2
         allocate (cell_corner_lat(nalloc),
     &        cell_corner_lon(nalloc))

         do corner = 1, ncorners
            cell_corner_lat(corner) = grid1_corner_lat(corner,cell_add)
            cell_corner_lon(corner) = grid1_corner_lon(corner,cell_add)
         enddo

         cell_center_lat = grid1_center_lat(cell_add)
         cell_center_lon = grid1_center_lon(cell_add)

         special_cell = special_polar_cell1(cell_add)


         !***
         !*** Also, allocate storage for the cell from the opposite grid
         !***

         opp_grid_num = 2
         ncorners_opp = grid2_corners
         nalloc_opp = ncorners_opp+2
         allocate (oppcell_corner_lat(nalloc_opp),
     &        oppcell_corner_lon(nalloc_opp))

      else

         !***
         !*** Set up the cell info with room to add degenerate edges
         !***

         ncorners = grid2_corners
         nalloc = ncorners + 2
         allocate (cell_corner_lat(nalloc),
     &        cell_corner_lon(nalloc))

         do corner = 1, ncorners
            cell_corner_lat(corner) = grid2_corner_lat(corner,cell_add)
            cell_corner_lon(corner) = grid2_corner_lon(corner,cell_add)
         enddo

         cell_center_lat = grid2_center_lat(cell_add)
         cell_center_lon = grid2_center_lon(cell_add)

         special_cell = special_polar_cell2(cell_add)

         !***
         !*** Also, allocate storage for the cell from the opposite grid
         !***

         opp_grid_num = 1
         ncorners_opp = grid1_corners
         nalloc_opp = ncorners_opp + 2
         allocate (oppcell_corner_lat(nalloc_opp),
     &        oppcell_corner_lon(nalloc_opp))

      endif

      if (special_cell) then

         !***
         !*** Special cell with only one corner at the pole Such cells
         !*** can have an artificially extreme distortion of the edges
         !*** when mapped to the Lambert projection because of the span
         !*** of longitudes on the edges So we will augment such cells
         !*** with degenerate edges at the pole so that the edges coming
         !*** off the pole will actually have the same longitude values
         !*** at both ends
         !***
         !***           lon_p           lon_p+    lon_p   lon_p-
         !***           pi/2            pi/2      pi/2     pi/2
         !***            *                 *--------*------*
         !***           / \                |               |
         !***          /   \               |               |
         !***         /     \              |               |
         !***        *       *             *               *
         !***    lon_p+      lon_p-      lon_p+           lon_p-
         !***    lat_p+      lat_p-      lat_p+           lat_p-
         !***

         call modify_polar_cell(ncorners,nalloc,cell_corner_lat,
     &        cell_corner_lon)

      endif

      !***
      !*** Cell info set up - Now process the cell
      !***

      do corner = 1, ncorners
         next_corn = mod(corner,ncorners) + 1

         !***
         !*** define endpoints of the current segment
         !***

         beglat = cell_corner_lat(corner)
         beglon = cell_corner_lon(corner)
         endlat = cell_corner_lat(next_corn)
         endlon = cell_corner_lon(next_corn)
         lrevers = .false.

         !***
         !*** if this is a constant-longitude segment, skip the rest 
         !*** since the line integral contribution will be zero.
         !***

         if ((phi_or_theta == 1 .and. endlon == beglon) .or.
     &        (phi_or_theta == 2 .and. endlat == beglat)) cycle

         !***
         !*** to ensure exact path taken during both
         !*** sweeps, always integrate segments in the same 
         !*** direction (SW to NE).
         !***

         if ((endlat < beglat) .or.
     &        (endlat == beglat .and. endlon < beglon)) then
            tmplat = beglat
            beglat = endlat
            endlat = tmplat
            tmplon = beglon
            beglon = endlon
            endlon = tmplon
            lrevers = .not. lrevers
         endif

         !*** But if one of the segment ends is in the polar region,
         !*** we want to start from that (makes some logic easier)

         if ((beglat < north_thresh .and. endlat > north_thresh) .or.
     &       (beglat > south_thresh .and. endlat < south_thresh)) 
     &        then
            tmplat = beglat
            beglat = endlat
            endlat = tmplat
            tmplon = beglon
            beglon = endlon
            endlon = tmplon
            lrevers = .not. lrevers
         endif

         begseg(1) = beglat
         begseg(2) = beglon

         fullseg_dlat = endlat-beglat
         fullseg_dlon = endlon-beglon
         if (fullseg_dlon >  pi) fullseg_dlon = fullseg_dlon - pi2
         if (fullseg_dlon < -pi) fullseg_dlon = fullseg_dlon + pi2
         fullseg_len2 = fullseg_dlat*fullseg_dlat + 
     &        fullseg_dlon*fullseg_dlon
         
         partseg_len2 = 0.0

         !***
         !*** Is this an edge on the boundary of the grid or 
         !*** on the boundary of the active cells
         !*** 

! Commented out by MD 
!         call find_adj_cell(cell_add, corner, grid_num, adj_add)
!         if (grid_num .eq. 1) then
!            if (adj_add .eq. 0 .or. .not. grid1_mask(adj_add)) then
!               bndedge = .true.
!            else
!               bndedge = .false.
!            endif
!         else
!            if (adj_add .eq. 0 .or. .not. grid2_mask(adj_add)) then
!               bndedge = .true.
!            else
!              bndedge = .false.
!            endif
!         endif

         call find_adj_cell(cell_add, corner, grid_num, adj_add)
         bndedge = .false.
         if (grid_num .eq. 1) then
            if (adj_add .eq. 0) then
               bndedge = .true.
            else
               if (.not. grid1_mask(adj_add)) then
                 bndedge = .true.
               endif
            endif
         else
            if (adj_add .eq. 0) then
               bndedge = .true.
            else
               if (.not. grid2_mask(adj_add)) then
                 bndedge = .true.
               endif
            endif
         endif

         !***
         !*** integrate along this segment, detecting intersections 
         !*** and computing the line integral for each sub-segment
         !***

         if (beglat .gt. north_thresh .or. beglat .lt. south_thresh)
     &        then
            nseg = npseg     ! Use multiple subsegments near the pole
            inpolar = .true.
         else
            nseg = 1
            inpolar = .false.
         endif


         last_add = 0
         last_lboundary = .false.
         last_endpt_inpoly = .false.
         last_endpt_onedge = .false.
         next_add = 0
         search = .true.
         ns = 1

! outer "do while"

         do while (beglat /= endlat .or. beglon /= endlon)
            
            L_exit_do=.false.         !NRL

            if ((ns .eq. nseg) .or. (inpolar .eqv. .false.)) then
               !
               ! Last subseg or out of the polar region
               ! Go directly to end of segment
               !
               endlat1 = endlat
               endlon1 = endlon
            else
               endlat1 = begseg(1) + ns*(fullseg_dlat)/nseg
               endlon1 = begseg(2) + ns*(fullseg_dlon)/nseg
            endif
            
            num_subseg = 0

! inner "do while"

            do while (beglat /= endlat1 .or. beglon /= endlon1)
               
               !*** 
               !*** If we integrated to the end or just past it (due to 
               !*** numerical errors), we are done with this segment
               !***

!NRL see notes below re: infinite "do while" loop
               L_exit_do=.false.                  !NRL
               if (partseg_len2 .ge. fullseg_len2) then
                  write(*,*)'partseg_len2 .ge. fullseg_len2'
                  write(*,*)'beglat,beglon = ',beglat,beglon
                  write(*,*)'endlat,endlon = ',endlat,endlon
                  write(*,*)'endlat1,endlon1 = ',endlat1,endlon1
                  write(*,*)'exiting inner do while loop'
                  L_exit_do=.true.                !NRL
                  exit
               end if

               !******************************************************
               !*** Try to find which cell of the opposite grid this
               !*** segment is starting in and where it is exiting this
               !*** cell
               !******************************************************
         
               vec1_lat = endlat1-beglat
               vec1_lon = endlon1-beglon
               if (vec1_lon > pi) vec1_lon = vec1_lon - pi2
               if (vec1_lon < -pi) vec1_lon = vec1_lon + pi2
               vec1_len = sqrt(vec1_lat*vec1_lat+vec1_lon*vec1_lon)
               vec1_lat = vec1_lat/vec1_len
               vec1_lon = vec1_lon/vec1_len
               
               offset = 100.0*tiny
               oppcell_add = 0
               delta = 10*tiny
               intrsct_success = .false.
               loutside = .false.
               lstuck = .false.
               lboundary1 = .false.
               lboundary2 = .false.
               lboundary3 = .false.
               
               do while (.not. intrsct_success) 

                  !*************************************************
                  !*** Find out which cell the segment starts in
                  !*************************************************

                  srch_success = .false.
                  if (search) then

                     !***
                     !*** Offset the start point in ever increasing 
                     !*** amounts until we are able to reliably locate 
                     !*** the point in a cell of grid2. Inability to locate 
                     !*** the point causes the offset amount to increase

                     it = 0
                     do while (.not. srch_success) 
                        
                        srchpt_lat = beglat + offset*vec1_lat
                        srchpt_lon = beglon + offset*vec1_lon
                        
                        call locate_point(srchpt_lat, srchpt_lon,
     &                       cell_add, grid_num, opp_grid_num, 
     &                       oppcell_add, lboundary1, bedgeid1)
                        
                        if (oppcell_add .eq. 0) then
                           loutside = .true.
! lcoinc added by MD
                           lcoinc = .false.
                           exit ! exit the search loop
                        else
                           if (oppcell_add .ne. last_add .or. lthresh) 
     &                          then
                              srch_success = .true.
                           else
                              offset = offset + delta
                              if (offset .ge. vec1_len) then
                                 exit
                              endif                               
                              if (it .gt. 3) then
                                 delta = 2.0*delta
                                 it = 0
                              endif
                           endif                        
                        endif
                        
                        it = it + 1                         
                     enddo   ! do while (.not. srch_success)
                     
                  else
                     if (last_endpt_inpoly) then

                        !*** We know the grid cell the end of the last 
                        !*** segment (which is the beginning of this 
                        !*** segment)

                        oppcell_add = last_add
                        lboundary1 = last_lboundary

                     else if (next_add .ne. 0) then

                        !*** We know the edge of the grid2 cell that the
                        !*** last segment intersected, so we move into
                        !*** the adjacent cell

                        oppcell_add = next_add
                        lboundary1 = .true.
                        
                     endif
                     
                     srch_success = .true.
                     
                  endif
                  
                  !*****************************************************
                  !*** Find where the segment exits this cell, if at all
                  !*****************************************************

                  if (srch_success) then 

                     !***
                     !*** First setup local copy of oppcell with room for
                     !*** adding degenerate edges
                     !***

                     if (grid_num .eq. 1) then
                        ncorners_opp = grid2_corners
                        do i = 1, ncorners_opp
                           oppcell_corner_lat(i) = 
     &                          grid2_corner_lat(i,oppcell_add)
                           oppcell_corner_lon(i) = 
     &                          grid2_corner_lon(i,oppcell_add)
                        enddo
                        oppcell_center_lat = 
     &                       grid2_center_lat(oppcell_add)
                        oppcell_center_lon = 
     &                       grid2_center_lon(oppcell_add)

                        special_cell = special_polar_cell2(oppcell_add)
                     else
                        ncorners_opp = grid1_corners
                        do i = 1, ncorners_opp
                           oppcell_corner_lat(i) = 
     &                          grid1_corner_lat(i,oppcell_add)
                           oppcell_corner_lon(i) =
     &                          grid1_corner_lon(i,oppcell_add)
                        enddo
                        oppcell_center_lat = 
     &                       grid1_center_lat(oppcell_add)
                        oppcell_center_lon = 
     &                       grid1_center_lon(oppcell_add)

                        special_cell = special_polar_cell1(oppcell_add)
                     endif

                     if (special_cell) then
                        call modify_polar_cell(ncorners_opp, nalloc_opp,
     &                       oppcell_corner_lat, oppcell_corner_lon)
                     endif
                     
                     !***
                     !*** First see if the segment end is in the same cell
                     !***

                     call ptincell(endlat1,endlon1, oppcell_add, 
     &                    ncorners_opp,
     &                    oppcell_corner_lat,oppcell_corner_lon,
     &                    oppcell_center_lat,oppcell_center_lon,
     &                    opp_grid_num,inpoly,
     &                    lboundary2,bedgeid2) 
                     
                     if (inpoly) then
                        intrsct_lat = endlat1
                        intrsct_lon = endlon1
                        intrsct_success = .true.                  
                        search = .false.
                        next_add = 0
                        last_add = oppcell_add        ! for next subseg
                        last_lboundary = lboundary2
                        last_endpt_inpoly = .true.
                        
                        if (lboundary1 .and. lboundary2) then
                  
                           !*** This is a edge on the boundary of the 
                           !*** active mesh and both of its endpoints 
                           !*** are on the boundary of the containing 
                           !*** cell. Check if the the segment is also 
                           !*** on the boundary

                           midlat = (beglat+endlat1)/2.0
                           if (abs(beglon-endlon1) .ge. pi) then
                              midlon = (beglon+endlon1)/2.0 - pi
                           else
                              midlon = (beglon+endlon1)/2.0
                           endif 
                           
                           call ptincell(midlat,midlon, oppcell_add, 
     &                          ncorners_opp,
     &                          oppcell_corner_lat, oppcell_corner_lon,
     &                          oppcell_center_lat, oppcell_center_lon,
     &                          opp_grid_num, inpoly, lboundary3,
     &                          bedgeid3) 
                           
                           if (inpoly .and. lboundary3) then
                              lcoinc = .true.
                              intedge = bedgeid3
                           endif
                           
                        else
                           lcoinc = .false.
                        endif
                        
                     else
                      
                        !***
                        !*** Do an intersection to find out where the 
                        !*** segment exits the cell
                        !***

                        call intersection(cell_add,grid_num,
     &                       beglat, beglon, endlat1, endlon1, 
     &                       begseg,
     &                       bedgeid1,
     &                       oppcell_add, ncorners_opp,
     &                       oppcell_corner_lat, oppcell_corner_lon,
     &                       opp_grid_num,
     &                       intrsct_lat, intrsct_lon, intedge,
     &                       sinang2, lcoinc, lthresh)
                        
                        if (intedge /= 0) then
                           intrsct_success = .true.
                           last_add = oppcell_add     ! for next subseg
                           last_endpt_onedge = .true.
                           last_endpt_inpoly = .false.
                           last_lboundary = .true.

                           if (.not. lthresh) then
                              call find_adj_cell(oppcell_add,intedge,
     &                             opp_grid_num,next_add)
                              if (next_add .ne. 0) then
                                 search = .false.
                              else
                                 search = .true.
                              endif
                           else
                              search = .true.
                           endif
                        endif
                        
                     endif
                     
                     if (.not. intrsct_success) then
                        
                        !*** Offset point and try again

                        search = .true.
                        delta = 2.0*delta
                        offset = offset + delta
                        if (offset .gt. vec1_len) then

                           ! Punt - exit the intersection loop

                           intrsct_lat = endlat1
                           intrsct_lon = endlon1
                           last_add = 0
                           last_lboundary = .false.
                           exit              

                        endif
                     endif

!NRL                 if (lcoinc .and. .not. bndedge) then

                     if (lcoinc .and. .not. bndedge  !NRL
     &                   .and. intedge /= 0) then    !NRL

                        !***
                        !*** Segment is coincident with edge of other grid
                        !*** which means it could belong to one of 2 cells
                        !*** Choose the cell such that edge that is 
                        !*** coincident with the segment is in the same
                        !*** dir as the segment

                        i = intedge
                        inext = mod(i,ncorners_opp)+1
                        vec2_lat = oppcell_corner_lat(inext) -
     &                          oppcell_corner_lat(i)
                        vec2_lon = oppcell_corner_lon(inext) - 
     &                          oppcell_corner_lon(i)

                        if (vec2_lon >  pi) vec2_lon = vec2_lon - pi2
                        if (vec2_lon < -pi) vec2_lon = vec2_lon + pi2
                        
                        dp = vec1_lat*vec2_lat + vec1_lon*vec2_lon
                        
                        if ((.not. lrevers .and. dp .lt. 0) .or.
     &                       (lrevers .and. dp .gt. 0)) then
                           
                           !*** Integrals from this segment must be
                           !*** assigned to the adjacent cell of
                           !*** opcell_add but only if such an adjacent
                           !*** cell exists

                           call find_adj_cell(oppcell_add, intedge, 
     &                          opp_grid_num, adj_add)

                           if (adj_add .gt. 0) then
                              oppcell_add = adj_add

                              if (grid_num .eq. 1) then
                                 ncorners_opp = grid2_corners
                                 do i = 1, ncorners_opp
                                    oppcell_corner_lat(i) = 
     &                                   grid2_corner_lat(i,oppcell_add)
                                    oppcell_corner_lon(i) = 
     &                                   grid2_corner_lon(i,oppcell_add)
                                 enddo
                                 oppcell_center_lat = 
     &                                grid2_center_lat(oppcell_add)
                                 oppcell_center_lon = 
     &                                grid2_center_lon(oppcell_add)
                                 
                                 special_cell = 
     &                                special_polar_cell2(oppcell_add)
                              else
                                 ncorners_opp = grid1_corners
                                 do i = 1, ncorners_opp
                                    oppcell_corner_lat(i) = 
     &                                   grid1_corner_lat(i,oppcell_add)
                                    oppcell_corner_lon(i) =
     &                                   grid1_corner_lon(i,oppcell_add)
                                 enddo
                                 oppcell_center_lat = 
     &                                grid1_center_lat(oppcell_add)
                                 oppcell_center_lon = 
     &                                grid1_center_lon(oppcell_add)

                                 special_cell = 
     &                                special_polar_cell1(oppcell_add)
                              endif
                              
                              if (special_cell) then
                                 call modify_polar_cell(ncorners_opp, 
     &                                nalloc_opp, oppcell_corner_lat, 
     &                                oppcell_corner_lon)
                              endif

                           endif

                        endif
                        
                     endif
                     
                  else

                     !***
                     !*** Could not locate a viable cell for the segment
                     !*** start
                     !***

                     if (oppcell_add .eq. 0) then                     
                        loutside = .true.
! lcoinc added by MD
                        lcoinc = .false.
                      
                        !***
                        !*** Take baby steps to see if any part of the
                        !*** segment is inside a cell of the other grid
                        !***

                        seg_outside = .false.
                        delta = vec1_len/100.00
                        offset = delta
                        do while (.not. srch_success)
                         
                           srchpt_lat = beglat + offset*vec1_lat
                           srchpt_lon = beglon + offset*vec1_lon
                           
                           call locate_point(srchpt_lat, srchpt_lon, 
     &                          cell_add, grid_num, opp_grid_num, 
     &                          oppcell_add, lboundary1, bedgeid1)
                           
                           if (oppcell_add /= 0) then
                              srch_success = .true.

                              !***
                              !*** Found a point of the segment in the
                              !*** cell. Do a bisection method to find
                              !*** the starting point of the segment 
                              !*** in the cell
                              !*** 

                              call converge_to_bdry(oppcell_add, 
     &                             opp_grid_num, ncorners_opp, 
     &                             oppcell_corner_lat, 
     &                             oppcell_corner_lon,
     &                             oppcell_center_lat,
     &                             oppcell_center_lon,
     &                             srchpt_lat, srchpt_lon,
     &                             beglat, beglon,
     &                             intrsct_lat, intrsct_lon, 
     &                             bedgeid1)

                              search = .false.
                              last_endpt_onedge = .true.
                              next_add = oppcell_add
                              last_lboundary = .true.
                                                         
                              oppcell_add = 0  

                           else

                              offset = offset + delta
                              
                              if (offset .ge. vec1_len) then
!                                 print *, 
!     &                                'Segment fully outside grid2'
!                                 print *, 'Segment of grid1_add',
!     &                                grid1_add
!                                 print *, beglat,beglon
!                                 print *, endlat1,endlon1

                                 seg_outside = .true.

                                 intrsct_lat = endlat1
                                 intrsct_lon = endlon1 
                                 
                                 search = .true.
                                 last_add = 0
                                 last_lboundary = .false.

                                 exit ! leave search loop
                              endif
                           endif
                         
                        enddo

                        ! int. loop
                        if (srch_success .or. seg_outside) exit 
                        
                     else
                        
                        if(is_master)then
                           print *, 'Unable to move out of last cell'
                           print *, 'Segment of edge ',corner,
     &                          ' of grid cell ',cell_add
                           print *, 'Stuck in opposite grid cell ',
     &                          oppcell_add
                           dist2 = 
     &                      (endlat1-begseg(1))*(endlat1-begseg(1)) +
     &                      (endlon1-begseg(2))*(endlon1-begseg(2))
                           print *, 'Fraction of segment left ',
     &                          vec1_len/sqrt(dist2)
                        endif
                        lstuck = .true.
                        
                        !***
                        !*** Punt - just assign the rest of the segment 
                        !*** to the current cell it is stuck in by
                        !*** tagging the segment endpoint as the 
                        !*** intersection point
                        !*** 

                        intrsct_lat = endlat1
                        intrsct_lon = endlon1

                        search = .true.
                        last_add = 0
                        last_lboundary = .false.
                        
                     endif

                     exit    ! exit the intersection loop
                     
                  endif      ! if (srch_success) then ... else ....

               end do        ! do while (.not. intrsct_success)

               !********************************************************
               !*** Compute the line integrals for this subsegment
               !********************************************************
            
               if (oppcell_add /= 0) then
                  call line_integral(phi_or_theta, weights, num_wts,
     &                 beglon, intrsct_lon, beglat, intrsct_lat,
     &                 cell_center_lat, cell_center_lon,
     &                 oppcell_center_lat, oppcell_center_lon)
               else
                  call line_integral(phi_or_theta, weights, num_wts,
     &                 beglon, intrsct_lon, beglat, intrsct_lat,
     &                 cell_center_lat, cell_center_lon,
     &                 cell_center_lat, cell_center_lon)
               endif
             
               !***
               !*** if integrating in reverse order, change
               !*** sign of weights
               !***

               if (lrevers) then
                  weights = -weights
               endif

               !***
               !*** store the appropriate addresses and weights. 
               !*** also add contributions to cell areas and centroids.
               !***

               if (grid_num .eq. 1) then

                  if (oppcell_add /= 0) then
                     if (grid1_mask(cell_add)) then
                        call store_link_cnsrv(cell_add, oppcell_add, 
     &                       weights)

C$OMP CRITICAL(block1)
!
!     Could have another thread that found an intersection between that
!     cell address and oppcell_add in which case it will try to write
!     into this address - we have to block that until we are finished
!
                        grid1_frac(cell_add) = 
     &                       grid1_frac(cell_add) + weights(1)

                        grid2_frac(oppcell_add) = 
     &                     grid2_frac(oppcell_add) + weights(num_wts+1)
C$OMP END CRITICAL(block1)
                     endif
                     
                  endif

C$OMP CRITICAL(block2)                  
                  grid1_area(cell_add) = grid1_area(cell_add) + 
     &                 weights(1)
                  grid1_centroid_lat(cell_add) = 
     &                 grid1_centroid_lat(cell_add) + weights(2)
                  grid1_centroid_lon(cell_add) = 
     &                 grid1_centroid_lon(cell_add) + weights(3)
C$OMP END CRITICAL(block2)

               else

                  !*** swap weights because in store_link_cnsrv
                  !*** we are always sending in grid1 weights first
                  !*** and then grid2 weights

                  do i = 1, num_wts                     
                     rev_weights(num_wts+i) = weights(i)
                     rev_weights(i) = weights(num_wts+i)
                  enddo

                  if (.not. lcoinc .and. oppcell_add /= 0) then
                     if (grid1_mask(oppcell_add)) then
                        call store_link_cnsrv(oppcell_add, cell_add, 
     &                       rev_weights)

C$OMP CRITICAL(block3)
!
!     Could have another thread that found an intersection between that
!     cell address and oppcell_add in which case it will try to write
!     into this address - we have to block that until we are finished
!
                        grid2_frac(cell_add) = 
     &                       grid2_frac(cell_add) + weights(1)

                        grid1_frac(oppcell_add) = 
     &                      grid1_frac(oppcell_add) + weights(num_wts+1)
C$OMP END CRITICAL(block3)

                     endif
                     
                  endif
                  
C$OMP CRITICAL(block4)
                  grid2_area(cell_add) = grid2_area(cell_add) + 
     &                 weights(1)
                  grid2_centroid_lat(cell_add) = 
     &                 grid2_centroid_lat(cell_add) + weights(2)
                  grid2_centroid_lon(cell_add) = 
     &                 grid2_centroid_lon(cell_add) + weights(3)
C$OMP END CRITICAL(block4)
               endif

               !***
               !*** reset beglat and beglon for next subsegment.
               !***

               beglat = intrsct_lat
               beglon = intrsct_lon

               !***
               !*** How far have we come from the start of the segment
               !***
             
               vec2_lat = intrsct_lat-begseg(1)
               vec2_lon = intrsct_lon-begseg(2)
               if (vec2_lon >  pi) vec2_lon = vec2_lon - pi2
               if (vec2_lon < -pi) vec2_lon = vec2_lon + pi2
               
               partseg_len2 = vec2_lat*vec2_lat + vec2_lon*vec2_lon

               !***
               !*** prevent infinite loops if integration gets stuck
               !*** near cell or threshold boundary
               !***

               num_subseg = num_subseg + 1
               if (num_subseg > max_subseg) then
                  print *, 
     &               'integration stalled: num_subseg exceeded limit'
                  print *, 'Cell ',cell_add
                  print *, 'Edge ',corner
                  print *, 'Grid ',1
                  dist2 = (endlat1-begseg(1))*(endlat1-begseg(1)) +
     &                 (endlon1-begseg(2))*(endlon1-begseg(2))
                  print *, 'Fraction of segment left ',
     &                 vec1_len/sqrt(dist2)
!                 exit       ! Give up and exit
                  stop       ! Give up and stop
               endif

! inner "do while"
            end do           ! do while (beglat /= endlat1 ...

!NRL We add an exit to outer do similar to exit of inner do:
!NRL   This was an apparent bug: exit statement would escape 
!NRL   inner do but then computation could not get out of 
!NRL   outer do since beglat, beglon controlling outer do
!NRL   never changed b/c it never gets to the part of the
!NRL   code that changes beglat, beglon, b/c it keeps 
!NRL   exiting inner do.

!NRL This should happen very rarely, so we have a print
!NRL   statement to notify user.

            if (L_exit_do)then                           ! NRL
               write(*,*)'partseg_len2,fullseg_len2 = ', ! NRL
     &                    partseg_len2,fullseg_len2      ! NRL
               write(*,*)'exiting outer do while loop'   ! NRL
               exit                                      ! NRL
            endif                                        ! NRL

            ns = ns + 1
            if ((beglat > 0 .and. beglat < north_thresh) .or. 
     &           (beglat < 0 .and. beglat > south_thresh)) 
     &           then
               inpolar = .false.
            endif

! outer "do while"
         end do              ! do while (beglat /= endlat ....

         call line_integral(phi_or_theta, weights, num_wts,
     &        begseg(2), endlon, begseg(1), endlat,
     &        cell_center_lat, 
     &        cell_center_lon,
     &        cell_center_lat, 
     &        cell_center_lon)

         !***
         !*** end of segment
         !***

      end do                    ! do corner=....

      end subroutine cell_integrate
!***********************************************************************


!***********************************************************************

      subroutine modify_polar_cell(ncorners, nalloc, cell_corner_lat, 
     &     cell_corner_lon)

      !*** Input variables

      integer (SCRIP_i4), intent(in) ::
     &     nalloc

      !*** In/Out Variables 

      integer (SCRIP_i4), intent(inout) ::
     &     ncorners
      real (SCRIP_r8), dimension(:), intent(inout) ::
     &     cell_corner_lat(nalloc),
     &     cell_corner_lon(nalloc)

      !*** Local variables

      integer (SCRIP_i4) ::
     &     npcorners,      ! Number of polar corners
     &     pcorner,        ! Index of the polar corner 
                           ! (if only 1 is found)
     &     corner,         ! Corner iterator variable
     &     previdx,        ! Index of previous corner to polar corner
     &     nextidx         ! Index of next corner to polar corner

      real (SCRIP_r8) ::
     &     pole_lat,       ! Latitude considered to be pole
     &     prevlon,        ! Latitude of previous corner to polar corner
     &     nextlon         ! Latitude of next corner to polar corner

      
      !***
      !*** Modify special cell with only one corner at the pole. Such
      !*** cells can have an artificially extreme distortion of the
      !*** edges when mapped to the Lambert projection because of the
      !*** span of longitudes on the edges So we will augment such
      !*** cells with degenerate edges at the pole so that the edges
      !*** coming off the pole will actually have the same longitude
      !*** values at both ends
      !***
      !***           lon_p           lon_p+    lon_p   lon_p-
      !***           pi/2            pi/2      pi/2     pi/2
      !***            *                 *--------*------*
      !***           / \                |               |
      !***          /   \               |               |
      !***         /     \              |               |
      !***        *       *             *               *
      !***    lon_p+      lon_p-      lon_p+           lon_p-
      !***    lat_p+      lat_p-      lat_p+           lat_p-
      !***

      
      !***
      !*** MAJOR ASSUMPTION HERE IS THAT CELL_CORNER_LAT AND 
      !*** CELL_CORNER_LON HAVE ROOM TO GROW
      !***


      pcorner = 0
      npcorners = 0
      do corner = 1, ncorners
         if (abs(abs(cell_corner_lat(corner))-pih) .le. 1.0e-05) then
            pcorner = corner
            pole_lat = cell_corner_lat(corner)
            npcorners = npcorners + 1
         endif
      enddo


      if (npcorners .ne. 1) return !*** Not the kind of cell we want
      
      previdx = mod((pcorner-1)-1+ncorners,ncorners) + 1
      prevlon = cell_corner_lon(previdx)
      
      nextidx = mod(pcorner,ncorners) + 1
      nextlon = cell_corner_lon(nextidx)
      
      !*** Move entries from pcorner+1 on back by one

      do corner = ncorners, pcorner+1, -1
         cell_corner_lat(corner+1) = cell_corner_lat(corner)
         cell_corner_lon(corner+1) = cell_corner_lon(corner)
      enddo

      !*** Add a corner after pcorner

      cell_corner_lat(pcorner+1) = pole_lat
      cell_corner_lon(pcorner+1) = nextlon
      
      ncorners = ncorners+1

      !*** Move entries from pcorner on back by one

      do corner = ncorners, pcorner, -1
         cell_corner_lat(corner+1) = cell_corner_lat(corner)
         cell_corner_lon(corner+1) = cell_corner_lon(corner)
      enddo

      !*** Add a corner before pcorner

      cell_corner_lat(pcorner) = pole_lat
      cell_corner_lon(pcorner) = prevlon

      ncorners = ncorners+1
      
      end subroutine modify_polar_cell


!***********************************************************************

      subroutine intersection(seg_cell_id, seg_grid_id,
     &     beglat, beglon, endlat, endlon, begseg, begedge,
     &     cell_id, ncorners, cell_corner_lat,
     &     cell_corner_lon, cell_grid_id, intrsct_lat, intrsct_lon,
     &     intedge, sinang2, lcoinc, lthresh)

!-----------------------------------------------------------------------
!
!     this routine finds the intersection of a line segment given by 
!     beglon, endlon, etc. with a cell from another grid
!     A coincidence flag is returned if the segment is entirely 
!     coincident with an edge of the opposite.  
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     intent(in): 
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4), intent(in) ::
     &     seg_cell_id    ! ID of cell that intersecting segment is from

      integer (SCRIP_i4), intent(in) ::
     &     seg_grid_id    ! ID of grid that intersecting segment is from

      real (SCRIP_r8), intent(in) :: 
     &     beglat, beglon,     ! beginning lat/lon endpoints for segment
     &     endlat, endlon      ! ending    lat/lon endpoints for segment

      real (SCRIP_r8), dimension(2), intent(inout) :: 
     &     begseg              ! begin lat/lon of full segment

      integer (SCRIP_i4), intent(in) ::
     &     begedge             ! edge that beginning point is on (can be 0)

      integer (SCRIP_i4), intent(in) ::
     &     cell_id             ! cell to intersect with
      
      integer (SCRIP_i4), intent(in) ::
     &     ncorners            ! number of corners of cell

      real (SCRIP_r8), dimension(ncorners), intent(in) ::
     &     cell_corner_lat,    ! coordinates of cell corners
     &     cell_corner_lon

      integer (SCRIP_i4), intent(in) ::
     &     cell_grid_id        ! which grid is the cell from?


!-----------------------------------------------------------------------
!
!     intent(out): 
!
!-----------------------------------------------------------------------

      real (SCRIP_r8), intent(out) ::
     &     intrsct_lat, 
     &     intrsct_lon          ! lat/lon coords of intersection

      real (SCRIP_r8), intent(out) ::
     &     sinang2              ! square of sine of angle between
                                ! intersecting lines

      integer (SCRIP_i4), intent(out) ::
     &     intedge              ! edge that is intersected

      logical (SCRIP_logical), intent(out) ::
     &     lcoinc               ! True if segment is coincident with 
                                ! a cell edge

      logical (SCRIP_logical), intent(out) ::
     &     lthresh              ! True if segment crosses threshold

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4) :: 
     &     n, next_n

      logical (SCRIP_logical) ::
     &     found, first

      real (SCRIP_r8) ::
     &     lon1, lon2,         ! local longitude variables for segment
     &     lat1, lat2,         ! local latitude  variables for segment
     &     grdlon1, grdlon2,   ! local longitude variables for grid cell
     &     grdlat1, grdlat2,   ! local latitude  variables for grid cell
     &     vec1_lat, vec1_lon, 
     &     vec2_lat, vec2_lon, !
     &     vec3_lat, vec3_lon, ! vectors and vector products used
     &     cross_product,      ! during grid search
     &     dot_product,        !
     &     lensqr1, lensqr2,   !
     &     lensqr3,            !
     &     s1, s2, determ,      
     &     mat1, mat2,         ! variables used for linear solve to
     &     mat3, mat4,         ! find intersection
     &     rhs1, rhs2,         !
     &     denom,               
     &     begsegloc(2),       ! local copy of full segment start
     &     dist2,              ! distance from start pt to intersection 
                               ! pt
     &     maxdist2,           ! max dist from start pt to any 
                               ! intersection pt
     &     max_intrsct_lat,    ! latitude of farthest intersection point
     &     max_intrsct_lon,    ! longitude of farthest intersection
                               ! point
     &     minlat, maxlat,     ! min and max latitudes of segment
     &     minlon, maxlon,     ! min and max longitudes of segment
     &     tmplat, tmplon


!-----------------------------------------------------------------------
!
!     initialize defaults, flags, etc.
!
!-----------------------------------------------------------------------

      lcoinc = .false.
      lthresh = .false.
      intedge = 0
      first = .true.

      lat1 = beglat
      lon1 = beglon
      lat2 = endlat
      lon2 = endlon

      ! No edge is allowed to span more than pi radians
      ! Accordingly transform one or the other end point

      if ((lon2-lon1) > pi) then
        lon2 = lon2 - pi2
      else if ((lon2-lon1) < -pi) then
        lon1 = lon1 - pi2
      endif
      s1 = zero

!-----------------------------------------------------------------------
!
!     loop over sides of the cell to find intersection with side
!     must check all sides for coincidences or intersections
!
!-----------------------------------------------------------------------

      if (beglat > north_thresh .or. beglat < south_thresh) then

         !*** Special intersection routine for cells near the pole
         !*** Intersection is done in a transformed space using 
         !*** multi-segmented representation of the cell

        call pole_intersection(cell_id,ncorners,
     &               cell_corner_lat,cell_corner_lon,cell_grid_id,
     &               beglat, beglon, endlat, 
     &               endlon, begseg, begedge,
     &               intedge,intrsct_lat,intrsct_lon,
     &               sinang2,lcoinc,lthresh)

        return

      endif


      maxdist2 = -9999999.0 

      begsegloc(1) = begseg(1)
      begsegloc(2) = begseg(2)

      lthresh = .false.
      intrsct_loop: do n=1,ncorners
         next_n = mod(n,ncorners) + 1
         
         grdlat1 = cell_corner_lat(n)
         grdlon1 = cell_corner_lon(n)
         grdlat2 = cell_corner_lat(next_n)
         grdlon2 = cell_corner_lon(next_n)

         lensqr2 = (grdlat1-grdlat2)*(grdlat1-grdlat2) +
     &        (grdlon1-grdlon2)*(grdlon1-grdlon2)

         if (lensqr2 .le. tiny*tiny) cycle       ! degenerate edge

         ! No edge can span more than pi radians

         if (grdlon2-grdlon1 > pi) then
            grdlon2 = grdlon2 - pi2
         else if (grdlon2-grdlon1 < -pi) then
            grdlon1 = grdlon1 - pi2
         endif

         ! Also the two intersecting segments together 
         ! cannot span more than 2*pi radians

         minlon = min(lon1,lon2)
         maxlon = max(grdlon1,grdlon2)
         if (maxlon-minlon > pi2) then
            grdlon1 = grdlon1 - pi2
            grdlon2 = grdlon2 - pi2
         else
            minlon = min(grdlon1,grdlon2)
            maxlon = max(lon1,lon2)
            if (maxlon-minlon > pi2) then
               grdlon1 = grdlon1 + pi2
               grdlon2 = grdlon2 + pi2
            endif
         endif     
         

        !***
        !*** set up linear system to solve for intersection
        !***

        mat1 = lat2 - lat1
        mat2 = grdlat1 - grdlat2
        mat3 = lon2 - lon1
        mat4 = grdlon1 - grdlon2
        rhs1 = grdlat1 - lat1
        rhs2 = grdlon1 - lon1

        determ = mat1*mat4 - mat2*mat3

        !***
        !*** if the determinant is zero, the segments are either 
        !***   parallel or coincident.  coincidences were detected 
        !***   above so do nothing.

        if (abs(determ) > tiny*tiny) then

          !*** if the determinant is non-zero, solve for the linear 
          !***   parameters s for the intersection point on each line 
          !***   segment.
          !*** if 0<s1,s2<1 then the segment intersects with this side.
          !***   return the point of intersection (adding a small
          !***   number so the intersection is off the grid line).
          !***

          s1 = (rhs1*mat4 - mat2*rhs2)/determ
          s2 = (mat1*rhs2 - rhs1*mat3)/determ

          if (s2 >= zero .and. s2 <= one .and.
     &        s1 >  zero .and. s1 <= one) then

             !***
             !*** recompute intersection based on full segment
             !*** so intersections are consistent for both sweeps
             !***

             if (lon2-begsegloc(2) > pi) then
                lon2 = lon2 - pi2
             else if (lon2-begsegloc(2) < -pi) then
                begsegloc(2) = begsegloc(2) - pi2
             endif


            ! Also the two intersecting segments together 
            ! cannot span more than 2*pi radians

             minlon = min(begsegloc(2),lon2)
             maxlon = max(grdlon1,grdlon2)
             if (maxlon-minlon > pi2) then
                grdlon1 = grdlon1 - pi2
                grdlon2 = grdlon2 - pi2
             else
                minlon = min(grdlon1,grdlon2)
                maxlon = max(begsegloc(2),lon2)
                if (maxlon-minlon > pi2) then
                   grdlon1 = grdlon1 + pi2
                   grdlon2 = grdlon2 + pi2
                endif
             endif


             mat1 = lat2 - begsegloc(1)
             mat3 = lon2 - begsegloc(2)
             rhs1 = grdlat1 - begsegloc(1)
             rhs2 = grdlon1 - begsegloc(2)

             determ = mat1*mat4 - mat2*mat3

             !***
             !*** sometimes due to roundoff, the previous 
             !*** determinant is non-zero, but the lines
             !*** are actually coincident.  if this is the
             !*** case, skip the rest.
             !***

             if (determ /= zero) then
                s1 = (rhs1*mat4 - mat2*rhs2)/determ
                s2 = (mat1*rhs2 - rhs1*mat3)/determ

                intrsct_lat = begsegloc(1) + mat1*s1
                intrsct_lon = begsegloc(2) + mat3*s1

                if (intrsct_lon < 0.0) then
                   intrsct_lon = intrsct_lon + pi2
                else if (intrsct_lon > pi2) then
                   intrsct_lon = intrsct_lon - pi2
                endif

                !***
                !*** Make sure the intersection point is not within
                !*** tolerance of the starting point
                !***

                if (first) then
                   max_intrsct_lat = intrsct_lat
                   max_intrsct_lon = intrsct_lon

                   vec1_lat = intrsct_lat-beglat
                   vec1_lon = intrsct_lon-beglon
                   if (vec1_lon > pi) then
                      vec1_lon = vec1_lon - pi2
                   else if (vec1_lon < -pi) then
                      vec1_lon = vec1_lon + pi2
                   endif

                   maxdist2 = vec1_lat*vec1_lat + vec1_lon*vec1_lon
                   dist2 = maxdist2
                   
                   denom = (mat1*mat1+mat2*mat2)*(mat3*mat3+mat4*mat4)
                   sinang2 = determ*determ/denom
                   intedge = n
                   first = .false.
                else                   
                   vec1_lat = intrsct_lat-beglat
                   vec1_lon = intrsct_lon-beglon
                   if (vec1_lon > pi) then
                      vec1_lon = vec1_lon - pi2
                   else if (vec1_lon < -pi) then
                      vec1_lon = vec1_lon + pi2
                   endif

                   dist2 = vec1_lat*vec1_lat + vec1_lon*vec1_lon

                   if (dist2 > maxdist2) then
                      if (begedge .eq. 0 .or. begedge .ne. n) then
                         max_intrsct_lat = intrsct_lat
                         max_intrsct_lon = intrsct_lon
                         maxdist2 = dist2
                         
                         denom = 
     &                       (mat1*mat1+mat2*mat2)*(mat3*mat3+mat4*mat4)
                         sinang2 = determ*determ/denom
                         intedge = n
                      endif
                   endif
                endif

             else
                print *, 'DEBUG: zero determ'
                stop
             endif
             
          endif

       else

          !***
          !*** Coincident lines or parallel lines
          !*** 

          cross_product = mat2*rhs2 - mat4*rhs1

          !***
          !*** If area of triangle formed by endlat,endlon and
          !*** the gridline is negligible then the lines are coincident
          !***


          if (abs(cross_product) < tiny) then

             dot_product = mat1*(-mat2) + mat3*(-mat4)

             lensqr1 = mat1*mat1 + mat3*mat3 ! length sqrd of input 
                                             ! segment

             if (dot_product < zero) then

                !***
                !*** Segments oriented in the same direction
                !***
             

                tmplat = grdlat2
                tmplon = grdlon2
                grdlat2 = grdlat1
                grdlon2 = grdlon1
                grdlat1 = tmplat
                grdlon1 = tmplon

             endif


             vec2_lat = grdlat1 - lat1
             vec2_lon = grdlon1 - lon1
             if (vec2_lon >  pi) vec2_lon = vec2_lon - pi2
             if (vec2_lon < -pi) vec2_lon = vec2_lon + pi2

             lensqr2 = vec2_lat*vec2_lat + vec2_lon*vec2_lon

             if (vec2_lat*mat1 + vec2_lon*mat3 < 0) then
                lensqr2 = -lensqr2
             endif

             vec3_lat = grdlat2 - lat1
             vec3_lon = grdlon2 - lon1
             if (vec3_lon >  pi) vec3_lon = vec3_lon - pi2
             if (vec3_lon < -pi) vec3_lon = vec3_lon + pi2             

             lensqr3 = (vec3_lat*vec3_lat+vec3_lon*vec3_lon)

             if (vec3_lat*mat1 + vec3_lon*mat3 < 0) then
                lensqr3 = -lensqr3
             endif
             
             found = .false.

             if (lensqr2 > 0) then
                if (lensqr2 <= lensqr1) then
                   intrsct_lat = grdlat1
                   intrsct_lon = grdlon1
                   found = .true.
                endif
             else
                if (lensqr3 > 0) then
                   if (lensqr3 > lensqr1) then
                      intrsct_lat = lat2
                      intrsct_lon = lon2
                      found = .true.
                   else
                      intrsct_lat = grdlat2
                      intrsct_lon = grdlon2
                      found = .true.
                   endif
                endif
             endif

             if (found) then

                dist2 = (intrsct_lat-beglat)*(intrsct_lat-beglat)+
     &               (intrsct_lon-beglon)*(intrsct_lon-beglon)

               !*** Coincidence intersection always wins

                max_intrsct_lat = intrsct_lat
                max_intrsct_lon = intrsct_lon
                maxdist2 = dist2
                sinang2 = 0
                intedge = n
                lcoinc = .true.

                exit intrsct_loop
             endif
                
          endif

       endif

       !*** restore lon1 and lon2 in case it got modified
         
       lon1 = beglon
       lon2 = endlon
       begsegloc(2) = begseg(2)
       if ((lon2-lon1) > pi) then
          lon2 = lon2 - pi2
       else if ((lon2-lon1) < -pi) then
          lon1 = lon1 - pi2
       endif

      end do intrsct_loop

      if (intedge .eq. 0) then
         return
      else
         if (maxdist2 < 1e6*tiny*tiny) then
            intedge = 0
            return
         else
            intrsct_lat = max_intrsct_lat
            intrsct_lon = max_intrsct_lon
         endif
      endif

!-----------------------------------------------------------------------
!
!     if the segment crosses a pole threshold, reset the intersection
!     to be the threshold latitude.  only check if this was not a
!     threshold segment since sometimes coordinate transform can end
!     up on other side of threshold again.
!
!-----------------------------------------------------------------------

      if (lthresh) then
        if (intrsct_lat < north_thresh .or. intrsct_lat > south_thresh)
     &      lthresh = .false.
      else if (lat1 > zero .and. intrsct_lat > north_thresh) then
!        intrsct_lat = north_thresh + tiny
        intrsct_lat = north_thresh
        mat1 = lat2 - begsegloc(1)
        mat3 = lon2 - begsegloc(2)
        s1 = (intrsct_lat - begsegloc(1))/mat1
        intrsct_lon     = begsegloc(2) + s1*mat3
        lthresh = .true.
      else if (lat1 < zero .and. intrsct_lat < south_thresh) then
!        intrsct_lat = south_thresh - tiny
        intrsct_lat = south_thresh
        mat1 = lat2 - begsegloc(1)
        mat3 = lon2 - begsegloc(2)
        s1 = (intrsct_lat - begsegloc(1))/mat1
        intrsct_lon     = begsegloc(2) + s1*mat3
        lthresh = .true.
      endif

      if (intrsct_lon < 0.0) then
         intrsct_lon = intrsct_lon + pi2
      else if (intrsct_lon > pi2) then
         intrsct_lon = intrsct_lon - pi2
      endif


!-----------------------------------------------------------------------

      end subroutine intersection

!***********************************************************************


      subroutine pole_intersection(location,ncorners,
     &                 cell_corners_lat,cell_corners_lon,cell_grid_id,
     &                 beglat, beglon, endlat, endlon, begseg, begedge,
     &                 intedge, intrsct_lat, intrsct_lon,
     &                 sinang2, lcoinc, lthresh)

!-----------------------------------------------------------------------
!
!     Special intersection routine for line segment in cell close to 
!     poles
!     A coordinate transformation (using a Lambert azimuthal
!     equivalent projection) is performed to perform the intersection
!     Also, since a straight line in lat-lon space is a curve in this
!     transformed space, we represent each edge of the cell as having 
!     'npseg' segments whose endpoints are mapped using the Lambert 
!     projection
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     intent(in): 
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4), intent(in) ::
     &     location             ! cell to intersect segment with

      integer (SCRIP_i4), intent(in) ::
     &     ncorners             ! Number of cell corners

      real    (SCRIP_r8), dimension(ncorners), intent(in) ::
     &     cell_corners_lat,    ! Cell corner coordinates
     &     cell_corners_lon

      integer (SCRIP_i4), intent(in) ::
     &     cell_grid_id         ! which grid is the cell from?

      real (SCRIP_r8), intent(in) :: 
     &     beglat, beglon,      ! beginning lat/lon coords for segment
     &     endlat, endlon       ! ending    lat/lon coords for segment

      real (SCRIP_r8), dimension(2), intent(inout) :: 
     &     begseg               ! begin lat/lon of full segment

      integer (SCRIP_i4), intent(in) ::
     &     begedge              ! edge on which segment start is on 
                                ! (can be 0)

!-----------------------------------------------------------------------
!
!     intent(out): 
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4), intent(out) ::
     &     intedge              ! Edge that segment intersects

      real (SCRIP_r8), intent(out) ::
     &     intrsct_lat,         ! lat/lon coords of intersection
     &     intrsct_lon          

      real (SCRIP_r8), intent(out) ::
     &     sinang2              ! square of sine of angle between
                                ! intersecting line segments

      logical (SCRIP_logical), intent(out) ::
     &     lcoinc               ! True if segment is coincident with 
                                ! a cell edge

      logical (SCRIP_logical), intent(inout) ::
     &     lthresh              ! True if segment crosses threshold 


!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4) :: 
     &     n, n1, next_n, prev_n,
     &     it, i, j, 
     &     ncorners2, 
     &     intedge1

      logical (SCRIP_logical) :: 
     &     first, 
     &     found

      real (SCRIP_r8) :: 
     &     pi4, rns,             ! north/south conversion
     &     x1, x2,               ! local x variables for segment
     &     y1, y2,               ! local y variables for segment
     &     grdx1, grdx2,         ! local x variables for grid cell
     &     grdy1, grdy2,         ! local y variables for grid cell
     &     grdlat1, grdlat2,     ! latitude vars for grid cell
     &     grdlon1, grdlon2,     ! longitude vars for grid cell
     &     vec1_y, vec1_x,       !
     &     vec2_y, vec2_x,       ! vectors and cross products used
     &     vec3_y, vec3_x,       ! 
     &     vec1_lat, vec1_lon,   !
     &     vec2_lat, vec2_lon,   !
     &     vec3_lon,             !
     &     cross_product,        !
     &     dot_product,          !
     &     s1, s2, determ,       ! variables used for linear solve to
     &     mat1, mat2,           !
     &     mat3, mat4,           ! find intersection
     &     rhs1, rhs2,           !
     &     denom,                !
     &     intrsct_x, intrsct_y, ! intersection coordinates in 
                                 ! transformed space
     &     max_intrsct_lat,      ! intersection point at max distance
     &     max_intrsct_lon,      ! from the start point
     &     dist2,                ! dist of intersection point from start
                                 ! point
     &     maxdist2,             ! max dist of intersection point from
                                 ! start pnt
     &     lensqr1, lensqr2,     ! various segment lengths
     &     lensqr3,             
     &     tmpx, tmpy, 
     &     tmplat, tmplon,
     &     ldummy

      !***
      !*** variables necessary if segment manages to hit pole
      !***

      real (SCRIP_r8), dimension(npseg*ncorners) ::
     &     cell_corners_lat_loc,! Lat/Lon coordinates of multi-segmented
     &     cell_corners_lon_loc ! version of cell 

      

!-----------------------------------------------------------------------
!
!     initialize defaults, flags, etc.
!
!-----------------------------------------------------------------------

      max_intrsct_lat = pi     ! intersection point at max distance
      max_intrsct_lon = 4*pi   ! from the start point

      intedge = 0
      first = .true.
      maxdist2 = -999999.00

      s1 = zero

!-----------------------------------------------------------------------
!
!     convert coordinates
!
!-----------------------------------------------------------------------

      if (beglat > zero) then
        pi4 = quart*pi
        rns = one
      else
        pi4 = -quart*pi
        rns = -one
      endif

      x1 = rns*two*sin(pi4 - half*beglat)*cos(beglon)
      y1 =     two*sin(pi4 - half*beglat)*sin(beglon)
      x2 = rns*two*sin(pi4 - half*endlat)*cos(endlon)
      y2 =     two*sin(pi4 - half*endlat)*sin(endlon)

      intrsct_x = x2
      intrsct_y = y2


!-----------------------------------------------------------------------
!
!     now that a cell is found, search for the next intersection.
!     loop over sides of the cell to find intersection with side
!     must check all sides for coincidences or intersections
!
!-----------------------------------------------------------------------


      if (abs(x1) .le. tiny .and. abs(y1) .le. tiny .and. 
     &     abs(x2) .le. tiny .and. abs(y2) .le. tiny) then

         !***
         !*** The segment is a polar segment which is degenerate
         !*** in the transformed Lambert space. Find out which
         !*** cell edge it is coincident with and find the 
         !*** point where the segment exits this cell (if at all)
         !*** NOTE 1: THIS MUST BE DONE IN LAT-LON SPACE
         !*** NOTE 2: CODE RELEVANT ONLY FOR INTEGRATION W.R.T. phi
         !***

         intrsct_loop1: do n = 1, ncorners
           next_n = mod(n,ncorners) + 1
           
           grdlat1  = cell_corners_lat(n)
           grdlon1  = cell_corners_lon(n)
           grdlat2  = cell_corners_lat(next_n)
           grdlon2  = cell_corners_lon(next_n)
           grdx1 = rns*two*sin(pi4 - half*grdlat1)*cos(grdlon1)
           grdy1 =     two*sin(pi4 - half*grdlat1)*sin(grdlon1)
           grdx2 = rns*two*sin(pi4 - half*grdlat2)*cos(grdlon2)
           grdy2 =     two*sin(pi4 - half*grdlat2)*sin(grdlon2)
          
           if (abs(grdx1) .le. tiny .and. abs(grdy1) .le. tiny .and.
     &          abs(grdx2) .le. tiny .and. abs(grdy2) .le. tiny) then
              
              !***
              !*** Found polar segment in cell
              !***

              vec1_lon = endlon-beglon
              if (vec1_lon .gt.  pi) vec1_lon = vec1_lon - pi2
              if (vec1_lon .lt. -pi) vec1_lon = vec1_lon + pi2

              vec2_lon = grdlon2-grdlon1
              if (vec2_lon .gt.  pi) vec2_lon = vec2_lon - pi2
              if (vec2_lon .lt. -pi) vec2_lon = vec2_lon + pi2

              if (vec1_lon*vec2_lon .lt. 0) then

                 !*** switch coordinates to simplify logic below

                 tmplat = grdlat2
                 tmplon = grdlon2
                 grdlat2 = grdlat1
                 grdlon2 = grdlon1
                 grdlat1 = tmplat
                 grdlon1 = tmplon                 
              endif

              vec2_lon = grdlon1 - beglon
              if (vec2_lon .gt. pi) vec2_lon = vec2_lon - pi2
              if (vec2_lon .lt. -pi) vec2_lon = vec2_lon + pi2

              vec3_lon = grdlon2 - beglon
              if (vec3_lon .gt. pi) vec3_lon = vec3_lon - pi2
              if (vec3_lon .lt. -pi) vec3_lon = vec3_lon + pi2

              found = .false.

              if (vec2_lon*vec1_lon > 0) then
                 if (abs(vec3_lon) < abs(vec1_lon)) then
                    intrsct_lon = grdlon2
                    found = .true.
                 else if (abs(vec2_lon) < abs(vec1_lon)) then
                    intrsct_lon = grdlon1   ! Shouldn't be here
                    found = .true.
                 endif
              else
                 if (vec3_lon*vec1_lon > 0) then
                    if (abs(vec3_lon) < abs(vec1_lon)) then
                       intrsct_lon = grdlon2
                       found = .true.
                    endif
                 endif

              endif

              if (found) then
                 intrsct_lat = endlat
                 lcoinc = .true.
                 sinang2 = 0.0
                 intedge = n
                 exit intrsct_loop1
              endif
 
           endif
         
        end do intrsct_loop1

        return
      endif

      


      !****
      !**** General intersection
      !****



      !***
      !*** Construct multi-segmented version of the cell
      !***

      i = 0
      do n = ncorners, 1, -1
         i = i+1
         n1 = mod(n,ncorners)+1         
         cell_corners_lat_loc(i) = cell_corners_lat(n1)
         cell_corners_lon_loc(i) = cell_corners_lon(n1)

         prev_n = n1-1
         if (prev_n .eq. 0) prev_n = ncorners ! how do we do (j-1+n)%n
                                              ! in F90 ?

         vec1_lat = cell_corners_lat(prev_n)-cell_corners_lat(n1)
         vec1_lon = cell_corners_lon(prev_n)-cell_corners_lon(n1)
         if (vec1_lon > pi) then
            vec1_lon = vec1_lon - pi2
         else if (vec1_lon < -pi) then
            vec1_lon = vec1_lon + pi2
         endif

         do j = 1, npseg-1
            i = i+1         
            cell_corners_lat_loc(i) = 
     &           cell_corners_lat(n1) + j*vec1_lat/npseg
            cell_corners_lon_loc(i) =
     &           cell_corners_lon(n1) + j*vec1_lon/npseg
         enddo
      enddo

      ncorners2 = npseg*ncorners



      !***
      !*** Now intersect segment with multi-segmented version of cell
      !***


      intrsct_loop2: do n= 1, ncorners2

        next_n = mod(n,ncorners2) + 1

        grdlat1  = cell_corners_lat_loc(n)
        grdlon1  = cell_corners_lon_loc(n)
        grdlat2  = cell_corners_lat_loc(next_n)
        grdlon2  = cell_corners_lon_loc(next_n)
        grdx1 = rns*two*sin(pi4 - half*grdlat1)*cos(grdlon1)
        grdy1 =     two*sin(pi4 - half*grdlat1)*sin(grdlon1)
        grdx2 = rns*two*sin(pi4 - half*grdlat2)*cos(grdlon2)
        grdy2 =     two*sin(pi4 - half*grdlat2)*sin(grdlon2)

        if ((grdx1-grdx2)*(grdx1-grdx2)+(grdy1-grdy2)*(grdy1-grdy2) .le.
     &       tiny*tiny) cycle


        !***
        !*** set up linear system to solve for intersection
        !***

        mat1 = x2 - x1
        mat2 = grdx1 - grdx2
        mat3 = y2 - y1
        mat4 = grdy1 - grdy2
        rhs1 = grdx1 - x1
        rhs2 = grdy1 - y1

        determ = mat1*mat4 - mat2*mat3

        !***
        !*** if the determinant is zero, the segments are either 
        !***   parallel or coincident or one segment has zero length.  

        !*** if the determinant is non-zero, solve for the linear 
        !***   parameters s for the intersection point on each line 
        !***   segment.
        !*** if 0<s1,s2<1 then the segment intersects with this side.
        !***   return the point of intersection (adding a small
        !***   number so the intersection is off the grid line).
        !***

        if (abs(determ) > 1.e-30) then

           s1 = (rhs1*mat4 - mat2*rhs2)/determ
           s2 = (mat1*rhs2 - rhs1*mat3)/determ
           
           if (s2 >= zero .and. s2 <= one .and.
     &          s1 >  tiny .and. s1 <= one) then

              intrsct_x = x1 + s1*mat1
              intrsct_y = y1 + s1*mat3
              
              !***
              !*** convert back to lat/lon coordinates
              !***

              if (abs(intrsct_x) .gt. tiny .or.
     &             abs(intrsct_y) .gt. tiny) then

                 intrsct_lon = rns*atan2(intrsct_y,intrsct_x)

              else

                 !*** Degenerate case - we don't have a good way of
                 !*** finding out what the longitude corresponding
                 !*** to a (0,0) intersection is. So we take the 
                 !*** the intersection as one of the two endpoints of
                 !*** the grid segment
                 
                 if (abs(abs(grdlat1)-pih) .lt. 1e-5 .and.
     &                abs(abs(grdlat2)-pih) .lt. 1e-5) then

                 !*** Both endpoints of the grid segment are at the pole
                 !*** but at different longitudes

                    vec1_lat = grdlat1-beglat
                    vec1_lon = grdlon1-beglon
                    if (vec1_lon > pi) then
                       vec1_lon = vec1_lon - pi2
                    else if (vec1_lon < -pi) then
                       vec1_lon = vec1_lon + pi2
                    endif
                    dist2 = vec1_lat*vec1_lat + vec1_lon*vec1_lon

                    vec2_lat = grdlat2-beglat
                    vec2_lon = grdlon2-beglon
                    if (vec2_lon > pi) then
                       vec2_lon = vec2_lon - pi2
                    else if (vec2_lon < -pi) then
                       vec2_lon = vec2_lon + pi2
                    endif
                    
                    !*** pick the endpoint of the grid segment that is
                    !*** farthest from the beg point of the segment

                    if ((vec1_lat*vec1_lat + vec1_lon*vec1_lon) .ge.
     &                   (vec2_lat*vec2_lat + vec2_lon*vec2_lon)) then
                       intrsct_lon = grdlon1
                    else
                       intrsct_lon = grdlon2
                    endif

                 else if (abs(abs(grdlat1)-pih) .lt. 1e-5) then
                    intrsct_lon = grdlon1
                 else if (abs(abs(grdlat2)-pih) .lt. 1e-5) then
                    intrsct_lon = grdlon2
                 endif

                 !*** Make sure this longitude is not outside the
                 !*** beglon,endlon range

                 vec1_lon = endlon-intrsct_lon
                 if (vec1_lon > pi) then
                    vec1_lon = vec1_lon - pi2
                 else if (vec1_lon < -pi) then
                    vec1_lon = vec1_lon + pi2
                 endif

                 vec2_lon = beglon-intrsct_lon
                 if (vec2_lon > pi) then
                    vec2_lon = vec2_lon - pi2
                 else if (vec2_lon < -pi) then
                    vec2_lon = vec2_lon + pi2
                 endif

                 !*** if vec1_lon and vec2_lon are of the same sign
                 !*** then intrsct_lon is outside the beglon,endlon
                 !*** range

                 if (vec1_lon*vec2_lon > 0) cycle

              endif

              if (intrsct_lon < zero) 
     &             intrsct_lon = intrsct_lon + pi2
           
              if (abs(intrsct_x) > 1.d-10) then
                 intrsct_lat = (pi4 - 
     &                asin(rns*half*intrsct_x/cos(intrsct_lon)))*two
                 ldummy = two*(pi4 -
     &           asin(sqrt(intrsct_x*intrsct_x+intrsct_y*intrsct_y)/2.))
              else if (abs(intrsct_y) > 1.d-10) then
                 intrsct_lat = (pi4 - 
     &                asin(half*intrsct_y/sin(intrsct_lon)))*two
                 ldummy = two*(pi4 -
     &           asin(sqrt(intrsct_x*intrsct_x+intrsct_y*intrsct_y)/2.))
              else
                 intrsct_lat = two*pi4
              endif


              !***
              !*** If there are multiple intersection points, accept the
              !*** one that is not on the edge we started from but is
              !*** closest to the start point - need this for 
              !*** intersection to work for non-convex edges
              !***
            
              if (first) then

                 intedge1 = (n-1)/npseg + 1
                 intedge1 = ncorners - intedge1 + 1  ! dir of edges was
                                                     ! reversed
                 if (intedge1 .ne. begedge) then

                    max_intrsct_lat = intrsct_lat
                    max_intrsct_lon = intrsct_lon
                    
                    vec1_lat = intrsct_lat-beglat
                    vec1_lon = intrsct_lon-beglon
                    if (vec1_lon > pi) then
                       vec1_lon = vec1_lon - pi2
                    else if (vec1_lon < -pi) then
                       vec1_lon = vec1_lon + pi2
                    endif
                    maxdist2 = vec1_lat*vec1_lat + vec1_lon*vec1_lon
                    dist2 = maxdist2
                    
                    denom = (mat1*mat1+mat2*mat2)*(mat3*mat3+mat4*mat4)
                    sinang2 = determ*determ/denom
                    intedge = intedge1
                 
                    first = .false.
                 endif

              else
                 vec1_lat = intrsct_lat-beglat
                 vec1_lon = intrsct_lon-beglon
                 if (vec1_lon > pi) then
                    vec1_lon = vec1_lon - pi2
                 else if (vec1_lon < -pi) then
                    vec1_lon = vec1_lon + pi2
                 endif
                 dist2 = vec1_lat*vec1_lat + vec1_lon*vec1_lon
                 
                 !*** if the first intersection was on the same edge
                 !*** as the starting edge or 
                 !*** the current intersection point is not on the 
                 !*** starting edge and the distance to the beginning 
                 !*** point is less than that of the previous 
                 !*** intersection accept this intersection

                 intedge1 = (n-1)/npseg + 1
                 intedge1 = ncorners - intedge1 + 1 ! dir of edges was
                                                    ! reversed 
                 if (dist2 > maxdist2) then
                    if (begedge == 0 .or. intedge1 .ne. begedge) then
                       max_intrsct_lat = intrsct_lat
                       max_intrsct_lon = intrsct_lon
                       maxdist2 = dist2
                       
                       denom = 
     &                      (mat1*mat1+mat2*mat2)*(mat3*mat3+mat4*mat4)
                       sinang2 = determ*determ/denom
                       intedge = intedge1
                    endif
                 endif
              endif
           endif   

        else

          !***
          !*** Coincident lines or parallel lines
          !*** 

           cross_product = mat2*rhs2 - mat4*rhs1
         
           if (abs(cross_product) < tiny) then
              
              dot_product = mat1*(-mat2) + mat3*(-mat4)
              
              !***
              !*** If area of triangle formed by x2,y2 and the gridline
              !*** is negligible then the lines are coincident
              !***

              lensqr1 = mat1*mat1 + mat3*mat3 ! length sqrd of input
                                              ! segment
              
              if (dot_product < zero) then
                tmpx = grdx2
                tmpy = grdy2
                tmplat = grdlat2
                tmplon = grdlon2
                grdx2 = grdx1
                grdy2 = grdy1
                grdlat2 = grdlat1
                grdlon2 = grdlon1
                grdx1 = tmpx
                grdy1 = tmpy
                grdlat1 = tmplat
                grdlon1 = tmplon
             endif

             
             vec2_x = grdx1 - x1
             vec2_y = grdy1 - y1
             lensqr2 = vec2_x*vec2_x + vec2_y*vec2_y
             if (vec2_x*mat1+vec2_y*mat3 < 0) then
                lensqr2 = -lensqr2
             endif
             

             vec3_x = grdx2 - x1
             vec3_y = grdy2 - y1             
             lensqr3 = (vec3_x*vec3_x+vec3_y*vec3_y)
             if (vec3_x*mat1+vec3_y*mat3 < 0) then
                lensqr3 = -lensqr3
             endif

             found = .false.

             if (lensqr2 > 0) then
                if (lensqr2 <= lensqr1) then
                   intrsct_x = grdx1
                   intrsct_y = grdy1
                   intrsct_lat = grdlat1
                   intrsct_lon = grdlon1
                   found = .true.
                endif
             else
                if (lensqr3 > 0) then
                   if (lensqr3 > lensqr1) then
                      intrsct_x = x2
                      intrsct_y = y2
                      intrsct_lat = endlat
                      intrsct_lon = endlon
                      found = .true.
                   else
                      intrsct_x = grdx2
                      intrsct_y = grdy2
                      intrsct_lat = grdlat2
                      intrsct_lon = grdlon2
                      found = .true.
                   endif
                endif
             endif
             
             if (found) then
                dist2 = (intrsct_lat-beglat)*(intrsct_lat-beglat)+
     &               (intrsct_lon-beglon)*(intrsct_lon-beglon)
                
                if (dist2 > tiny*tiny) then
                   
                   !*** Coincidence intersection always wins
                   
                   max_intrsct_lat = intrsct_lat
                   max_intrsct_lon = intrsct_lon
                   maxdist2 = dist2
                   sinang2 = 0
                   intedge = (n-1)/npseg + 1
                   intedge = ncorners - intedge + 1
                   lcoinc = .true.
                   
                   exit intrsct_loop2
                endif
             endif


           endif            ! if (abs(cross_product) < tiny)
              
        endif               ! if (abs(determ) > 1.e-30) .. else .. endif

      end do intrsct_loop2

      if (maxdist2 < 1e6*tiny*tiny)  then
         intedge = 0
         return
      else
         intrsct_lat = max_intrsct_lat
         intrsct_lon = max_intrsct_lon
      endif

!-----------------------------------------------------------------------
!
!     if segment manages to cross over pole, shift the beginning 
!     endpoint in order to avoid hitting pole directly
!     (it is ok for endpoint to be pole point)
!
!-----------------------------------------------------------------------

      if (abs(intrsct_x) < 1.e-10 .and. abs(intrsct_y) < 1.e-10 .and.
     &    (x2 /= zero .and. y2 /=0)) then
        if (avoid_pole_count > 2) then
           avoid_pole_count = 0
           avoid_pole_offset = 10.*avoid_pole_offset
        endif

        cross_product = x1*(y2-y1) - y1*(x2-x1)
        intrsct_lat = beglat
        if (cross_product*intrsct_lat > zero) then
          intrsct_lon = beglon    + avoid_pole_offset
        else
          intrsct_lon = beglon    - avoid_pole_offset
        endif

        avoid_pole_count = avoid_pole_count + 1
      else
        avoid_pole_count = 0
        avoid_pole_offset = tiny
      endif

!-----------------------------------------------------------------------
!
!     if the segment crosses a pole threshold, reset the intersection
!     to be the threshold latitude and do not reuse x,y intersect
!     on next entry.  only check if did not cross threshold last
!     time - sometimes the coordinate transformation can place a
!     segment on the other side of the threshold again
!
!-----------------------------------------------------------------------

      if (lthresh) then
        if (intrsct_lat > north_thresh .or. intrsct_lat < south_thresh)
     &    lthresh = .false.
      else if (beglat > zero .and. intrsct_lat < north_thresh) then
        mat4 = endlat - begseg(1)
        mat3 = endlon - begseg(2)
        if (mat3 >  pi) mat3 = mat3 - pi2
        if (mat3 < -pi) mat3 = mat3 + pi2
!        intrsct_lat = north_thresh - tiny
        intrsct_lat = north_thresh
        s1 = (north_thresh - begseg(1))/mat4
        intrsct_lon = begseg(2) + s1*mat3
        lthresh = .true.
      else if (beglat < zero .and. intrsct_lat > south_thresh) then
        mat4 = endlat - begseg(1)
        mat3 = endlon - begseg(2)
        if (mat3 >  pi) mat3 = mat3 - pi2
        if (mat3 < -pi) mat3 = mat3 + pi2
!        intrsct_lat = south_thresh + tiny
        intrsct_lat = south_thresh
        s1 = (south_thresh - begseg(1))/mat4
        intrsct_lon = begseg(2) + s1*mat3
        lthresh = .true.
      endif


!-----------------------------------------------------------------------

      end subroutine pole_intersection

!***********************************************************************



      subroutine line_integral(phi_or_theta, weights, num_wts, 
     &                       in_phi1, in_phi2, theta1, theta2,
     &                       grid1_lat, grid1_lon, grid2_lat, grid2_lon)

!-----------------------------------------------------------------------
!
!     this routine computes the line integral of the flux function 
!     that results in the interpolation weights.  the line is defined
!     by the input lat/lon of the endpoints.
!
!-----------------------------------------------------------------------
      
!-----------------------------------------------------------------------
!
!     intent(in):
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4), intent(in) :: 
     &     phi_or_theta         ! Integration variable (lat or lon)

      integer (SCRIP_i4), intent(in) ::
     &     num_wts              ! number of weights to compute

      real (SCRIP_r8), intent(in) :: 
     &     in_phi1, in_phi2,     ! longitude endpoints for the segment
     &     theta1, theta2,       ! latitude  endpoints for the segment
     &     grid1_lat, grid1_lon, ! reference coordinates for each
     &     grid2_lat, grid2_lon  ! grid (to ensure correct 0,2pi interv.

!-----------------------------------------------------------------------
!
!     intent(out):
!
!-----------------------------------------------------------------------

      real (SCRIP_r8), dimension(2*num_wts), intent(out) ::
     &     weights   ! line integral contribution to weights


!     write(*,*)'subroutine line_integral'
      if (phi_or_theta .eq. 1) then
         call line_integral_phi(weights, num_wts, in_phi1, in_phi2,
     &        theta1, theta2, grid1_lat, grid1_lon, 
     &        grid2_lat, grid2_lon)
      else
         call line_integral_theta(weights, num_wts,in_phi1,in_phi2,
     &        theta1, theta2, grid1_lat, grid1_lon, 
     &        grid2_lat, grid2_lon)
      endif


      return

!-----------------------------------------------------------------------

      end subroutine line_integral

!***********************************************************************



      subroutine line_integral_phi(weights, num_wts, 
     &                       in_phi1, in_phi2, theta1, theta2,
     &                       grid1_lat, grid1_lon, grid2_lat, grid2_lon)

!-----------------------------------------------------------------------
!
!     this routine computes the line integral of the flux function 
!     that results in the interpolation weights.  the line is defined
!     by the input lat/lon of the endpoints. Integration is w.r.t. lon
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     intent(in):
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4), intent(in) ::
     &        num_wts  ! number of weights to compute

      real (SCRIP_r8), intent(in) :: 
     &     in_phi1, in_phi2,     ! longitude endpoints for the segment
     &     theta1, theta2,       ! latitude  endpoints for the segment
     &     grid1_lat, grid1_lon, ! reference coordinates for each
     &     grid2_lat, grid2_lon  ! grid (to ensure correct 0,2pi interv.

!-----------------------------------------------------------------------
!
!     intent(out):
!
!-----------------------------------------------------------------------

      real (SCRIP_r8), dimension(2*num_wts), intent(out) ::
     &     weights   ! line integral contribution to weights

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      real (SCRIP_r8) :: dphi, sinth1, sinth2, costh1, costh2, fac,
     &                        phi1, phi2
      real (SCRIP_r8) :: f1, f2, fint

!-----------------------------------------------------------------------
!
!     weights for the general case based on a trapezoidal approx to
!     the integrals.
!
!-----------------------------------------------------------------------


!     write(*,*)'subroutine line_integral_phi'

      sinth1 = SIN(theta1)
      sinth2 = SIN(theta2)
      costh1 = COS(theta1)
      costh2 = COS(theta2)

      dphi = in_phi1 - in_phi2
      if (dphi >  pi) then
        dphi = dphi - pi2
      else if (dphi < -pi) then
        dphi = dphi + pi2
      endif
      dphi = half*dphi

!-----------------------------------------------------------------------
!
!     the first weight is the area overlap integral. the second and
!     fourth are second-order latitude gradient weights.
!
!-----------------------------------------------------------------------

      weights(        1) = dphi*(sinth1 + sinth2)
      write(401,*)weights(1),' % A'
      weights(num_wts+1) = dphi*(sinth1 + sinth2)
      weights(        2) = dphi*(costh1 + costh2 + (theta1*sinth1 +
     &                                              theta2*sinth2))
      weights(num_wts+2) = dphi*(costh1 + costh2 + (theta1*sinth1 +
     &                                              theta2*sinth2))

!-----------------------------------------------------------------------
!
!     the third and fifth weights are for the second-order phi gradient
!     component.  must be careful of longitude range.
!
!-----------------------------------------------------------------------

      f1 = half*(costh1*sinth1 + theta1)
      f2 = half*(costh2*sinth2 + theta2)

      phi1 = in_phi1 - grid1_lon
      if (phi1 >  pi) then
        phi1 = phi1 - pi2
      else if (phi1 < -pi) then
        phi1 = phi1 + pi2
      endif

      phi2 = in_phi2 - grid1_lon
      if (phi2 >  pi) then
        phi2 = phi2 - pi2
      else if (phi2 < -pi) then
        phi2 = phi2 + pi2
      endif

      if ((phi2-phi1) <  pi .and. (phi2-phi1) > -pi) then
        weights(3) = dphi*(phi1*f1 + phi2*f2)
      else
        if (phi1 > zero) then
          fac = pi
        else
          fac = -pi
        endif
        fint = f1 + (f2-f1)*(fac-phi1)/abs(dphi)
        weights(3) = half*phi1*(phi1-fac)*f1 -
     &               half*phi2*(phi2+fac)*f2 +
     &               half*fac*(phi1+phi2)*fint
      endif

      phi1 = in_phi1 - grid2_lon
      if (phi1 >  pi) then
        phi1 = phi1 - pi2
      else if (phi1 < -pi) then
        phi1 = phi1 + pi2
      endif

      phi2 = in_phi2 - grid2_lon
      if (phi2 >  pi) then
        phi2 = phi2 - pi2
      else if (phi2 < -pi) then
        phi2 = phi2 + pi2
      endif

      if ((phi2-phi1) <  pi .and. (phi2-phi1) > -pi) then
        weights(num_wts+3) = dphi*(phi1*f1 + phi2*f2)
      else
        if (phi1 > zero) then
          fac = pi
        else
          fac = -pi
        endif
        fint = f1 + (f2-f1)*(fac-phi1)/abs(dphi)
        weights(num_wts+3) = half*phi1*(phi1-fac)*f1 -
     &                       half*phi2*(phi2+fac)*f2 +
     &                       half*fac*(phi1+phi2)*fint
      endif

!-----------------------------------------------------------------------

      end subroutine line_integral_phi

!***********************************************************************



!***********************************************************************

      subroutine line_integral_theta(weights, num_wts, 
     &                       in_phi1, in_phi2, theta1, theta2,
     &                       grid1_lat, grid1_lon, grid2_lat, grid2_lon)

!-----------------------------------------------------------------------
!
!    this routine computes the line integral of the flux function 
!    that results in the interpolation weights.  the line is defined
!    by the input lat/lon of the endpoints. Integration is w.r.t. lat
!
!    Needed to use Simpson rule for this integration to get lower errors
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     intent(in):
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4), intent(in) ::
     &        num_wts  ! number of weights to compute

      real (SCRIP_r8), intent(in) :: 
     &     in_phi1, in_phi2,     ! longitude endpoints for the segment
     &     theta1, theta2,       ! latitude  endpoints for the segment
     &     grid1_lat, grid1_lon, ! reference coordinates for each
     &     grid2_lat, grid2_lon  ! grid (to ensure correct 0,2pi interv.

!-----------------------------------------------------------------------
!
!     intent(out):
!
!-----------------------------------------------------------------------

      real (SCRIP_r8), dimension(2*num_wts), intent(out) ::
     &     weights   ! line integral contribution to weights

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      real (SCRIP_r8) :: dtheta, dtheta2, costh1, costh2, costhpi,
     &                   phi1, phi2, theta_pi, f1, f2, fpi, 
     &                   fm, costhm, part1, part2

!-----------------------------------------------------------------------
!
!     weights for the general case based on a trapezoidal approx to
!     the integrals.
!
!-----------------------------------------------------------------------

      costh1 = COS(theta1)
      costh2 = COS(theta2)
      costhm = COS(half*(theta1+theta2))

      dtheta = theta2 - theta1
      dtheta2 = half*dtheta

!      write(*,*)' subroutine line_integral_theta'


!-----------------------------------------------------------------------
!
!     Need to account for double value of longitude in calculations of
!     all the weights.  First we transform all the phis to be relative
!     to the grid center This takes care of a good number of cases where
!     the the phis span the periodic boundary in the longitudinal
!     direction. If we still have a line that spans the periodic
!     boundary then we have to integrate along the line in two parts -
!     from point 1 to the periodic boundary and from the periodic
!     boundary to the second point
!
!     Example: Consider a line which has points at phi1 = -100 and phi2
!     = 100 degrees and say the grid center is at phi_c = 0
!     degrees. Then phi1-phi_c > -180 and phi2-phi_c < 180. But
!     phi2-phi1 > 180.
!
!     *********************************************!!!!!!!!!!!
!     If we are doing the second step anyway, why are we normalizing the
!     coordinates with respect to the grid centers?
!
!     We need it particularly in this integration because phi figures
!     explicitly in the expressions - so if a cell straddles the 0,2pi
!     boundary, we integrate some edges with phi values close to zero
!     and others with phi values close to 2pi leading to errors
!     *********************************************!!!!!!!!!!!
!
!-----------------------------------------------------------------------

      phi1 = in_phi1 - grid1_lon
      if (phi1 >  pi) then
        phi1 = phi1 - pi2
      else if (phi1 < -pi) then
        phi1 = phi1 + pi2
      endif

      phi2 = in_phi2 - grid1_lon
      if (phi2 >  pi) then
        phi2 = phi2 - pi2
      else if (phi2 < -pi) then
        phi2 = phi2 + pi2
      endif

      f1 = phi1*costh1
      f2 = phi2*costh2

      if ((phi2-phi1) <  pi .and. (phi2-phi1) > -pi) then

         fm = half*(phi1+phi2)*costhm

        weights(1) = dtheta*(f1 + 4*fm + f2)/6.0
!        write(401,*)weights(1),' % A'

        weights(2) = dtheta2*(theta1*f1 + theta2*f2)

        weights(3) = half*dtheta2*(f1*f1 + f2*f2)

      else
         if (phi1 > zero) then  ! Means phi2-phi1 < -pi

!           theta at phi = pi
            theta_pi = theta1 + (pi - phi1)*dtheta/(phi2 + pi2 - phi1)
!            print *, ''
!            print *, 'phi1',phi1,'    phi2',phi2
!            print *, 'theta1',theta1,'    theta2',theta2
!            print *, 'theta_pi',theta_pi
            
            costhpi = COS(theta_pi)
            fpi = pi*costhpi

            fm = half*(phi1+pi)*cos(half*(theta1+theta_pi))            
            part1 = (theta_pi - theta1)*(f1 + 4*fm + fpi)/6.0

            fm = half*(phi2-pi)*cos(half*(theta1+theta_pi))
            part2 = 0.5*(theta2 - theta_pi)*(-fpi + 4*fm + f2)/6.0

            weights(1) = part1 + part2 
!            write(401,*)weights(1),' % B'
            
            part1 = 0.5*(theta_pi - theta1)*(theta1*f1 + theta_pi*fpi)
            part2 = 0.5*(theta2 - theta_pi)*(-theta_pi*fpi + theta2*f2)
            weights(2) = part1 + part2  
            
            
         else                   ! Means phi2-phi1 > pi
       
!           theta at phi = -pi
            theta_pi = theta1 + (-pi - phi1)*dtheta/(phi2 - pi2 - phi1)
!            print *, ''
!            print *, 'phi1',phi1,'    phi2',phi2
!            print *, 'theta1',theta1,'    theta2',theta2
!            print *, 'theta_pi',theta_pi
            
            costhpi = COS(theta_pi)
            fpi = pi*costhpi
            
            fm = half*(phi1-pi)*cos(half*(theta1+theta_pi))
            part1 = 0.5*(theta_pi - theta1)*(f1 + 4*fm - fpi)/6.0

            fm = half*(pi+phi2)*cos(half*(theta2+theta_pi))
            part2 = 0.5*(theta2 - theta_pi)*(fpi + 4*fm + f2)/6.0
            weights(1) = part1 + part2
!            write(401,*)weights(1),' % C'
            
            part1 = 0.5*(theta_pi - theta1)*(theta1*f1 - theta_pi*fpi)
            part2 = 0.5*(theta2 - theta_pi)*(theta_pi*fpi + theta2*f2)
            weights(2) = part1 + part2
            
            
         endif
         
         part1 = 0.25*(theta_pi - theta1)*(f1*f1 + fpi*fpi)
         part2 = 0.25*(theta2 - theta_pi)*(fpi*fpi + f2*f2)
         weights(3) = part1 + part2
         
      endif
      
      
      phi1 = in_phi1 - grid2_lon
      if (phi1 >  pi) then
         phi1 = phi1 - pi2
      else if (phi1 < -pi) then
         phi1 = phi1 + pi2
      endif
      
      phi2 = in_phi2 - grid2_lon
      if (phi2 >  pi) then
         phi2 = phi2 - pi2
      else if (phi2 < -pi) then
         phi2 = phi2 + pi2
      endif
      
      
      f1 = phi1*costh1
      f2 = phi2*costh2
      
      if ((phi2-phi1) <  pi .and. (phi2-phi1) > -pi) then

         fm = half*(phi1+phi2)*costhm

         weights(num_wts+1) = dtheta2*(f1 + f2)
         
         weights(num_wts+2) = dtheta2*(theta1*f1 + theta2*f2)
         
         weights(num_wts+3) = half*dtheta2*(f1*f1 + f2*f2)
         
      else
         if (phi1 > zero) then
            
            theta_pi = theta1 + (pi - phi1)*dtheta/(phi2 + pi2 - phi1)
!            print *, ''
!            print *, 'phi1',phi1,'    phi2',phi2
!            print *, 'theta1',theta1,'    theta2',theta2
!            print *, 'theta_pi',theta_pi
            
            costhpi = COS(theta_pi)
            fpi = pi*costhpi

            fm = half*(phi1+pi)*cos(half*(theta1+theta_pi))
            part1 = (theta_pi - theta1)*(f1 + 4*fm + fpi)/6.0

            fm = half*(-pi+phi2)*cos(half*(theta2+theta_pi))
            part2 = (theta2 - theta_pi)*(-fpi + 4*fm + f2)/6.0
            weights(num_wts+1) = part1 + part2
            
            part1 = 0.5*(theta_pi - theta1)*(theta1*f1 + theta_pi*fpi)
            part2 = 0.5*(theta2 - theta_pi)*(-theta_pi*fpi + theta2*f2)
            weights(num_wts+2) = part1 + part2
            
            
         else
            
            theta_pi = theta1 + (-pi - phi1)*dtheta/(phi2 - pi2 - phi1)
!            print *, ''
!            print *, 'phi1',phi1,'    phi2',phi2
!            print *, 'theta1',theta1,'    theta2',theta2
!            print *, 'theta_pi',theta_pi
            
            costhpi = COS(theta_pi)
            fpi = pi*costhpi
            
            fm = half*(phi1-pi)*cos(half*(theta1+theta_pi))
            part1 = (theta_pi - theta1)*(f1 +4*fm - fpi)/6.0

            fm = half*(phi2+pi)*cos(half*(theta2+theta_pi))
            part2 = 0.5*(theta2 - theta_pi)*(fpi + 4*fm + f2)/6.0
            weights(num_wts+1) = part1 + part2
            
            part1 = 0.5*(theta_pi - theta1)*(theta1*f1 - theta_pi*fpi)
            part2 = 0.5*(theta2 - theta_pi)*(theta_pi*fpi + theta2*f2)
            weights(num_wts+2) = part1 + part2
            
         endif
         
         part1 = 0.25*(theta_pi - theta1)*(f1*f1 + fpi*fpi)
         part2 = 0.25*(theta2 - theta_pi)*(fpi*fpi + f2*f2)
         weights(num_wts+3) = part1 + part2
         
      endif

!-----------------------------------------------------------------------

      end subroutine line_integral_theta

!***********************************************************************



      subroutine store_link_cnsrv(add1, add2, weights)

!-----------------------------------------------------------------------
!
!     this routine stores the address and weight for this link in
!     the appropriate address and weight arrays and resizes those
!     arrays if necessary.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     input variables
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4), intent(in) ::
     &        add1,  ! address on grid1
     &        add2   ! address on grid2

      real (SCRIP_r8), dimension(:), intent(in) ::
     &        weights ! array of remapping weights for this link

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4) :: nlink, min_link, max_link ! link index

      logical (SCRIP_logical) :: found


!-----------------------------------------------------------------------
!
!     if all weights are zero, do not bother storing the link
!
!-----------------------------------------------------------------------

      if (all(weights == zero)) return

!-----------------------------------------------------------------------
!
!     restrict the range of links to search for existing links
!
!-----------------------------------------------------------------------

C$OMP CRITICAL(block5)
!     first_call should be within critical block or else multiple 
!     threads will see it as true the first time around

      if (first_call_store_link_cnsrv) then
        allocate(link_add1(2,grid1_size), link_add2(2,grid2_size))
        link_add1 = 0
        link_add2 = 0
        first_call_store_link_cnsrv = .false.
        min_link = 1
        max_link = 0
      else
        min_link = min(link_add1(1,add1),link_add2(1,add2))
        max_link = max(link_add1(2,add1),link_add2(2,add2))
        if (min_link == 0) then
          min_link = 1
          max_link = 0
        endif
      endif
C$OMP END CRITICAL(block5)

!-----------------------------------------------------------------------
!
!     if the link already exists, add the weight to the current weight
!     arrays
!
!-----------------------------------------------------------------------

      found = .false.

      do nlink=min_link,max_link
        if (add1 == grid1_add_map1(nlink)) then
        if (add2 == grid2_add_map1(nlink)) then

C$OMP CRITICAL(block3a)
          wts_map1(:,nlink) = wts_map1(:,nlink) + weights(1:num_wts)
          if (num_maps == 2) then
            wts_map2(:,nlink) = wts_map2(:,nlink) + 
     &                                  weights(num_wts+1:2*num_wts)
          endif
C$OMP END CRITICAL(block3a)
          found = .true.
          exit

        endif
        endif
      end do


      if (found) return

!-----------------------------------------------------------------------
!
!     if the link does not yet exist, increment number of links and 
!     check to see if remap arrays need to be increased to accomodate 
!     the new link.  then store the link.
!
!-----------------------------------------------------------------------

C$OMP CRITICAL(block6)

      num_links_map1  = num_links_map1 + 1
      if (num_links_map1 > max_links_map1) 
     &   call resize_remap_vars(1,resize_increment)

      grid1_add_map1(num_links_map1) = add1
      grid2_add_map1(num_links_map1) = add2
      wts_map1    (:,num_links_map1) = weights(1:num_wts)

      if (num_maps > 1) then
        num_links_map2  = num_links_map2 + 1
        if (num_links_map2 > max_links_map2) 
     &     call resize_remap_vars(2,resize_increment)

        grid1_add_map2(num_links_map2) = add1
        grid2_add_map2(num_links_map2) = add2
        wts_map2    (:,num_links_map2) = weights(num_wts+1:2*num_wts)
      endif

      if (link_add1(1,add1) == 0) link_add1(1,add1) = num_links_map1
      if (link_add2(1,add2) == 0) link_add2(1,add2) = num_links_map1
      link_add1(2,add1) = num_links_map1
      link_add2(2,add2) = num_links_map1

C$OMP END CRITICAL(block6)

!-----------------------------------------------------------------------

      end subroutine store_link_cnsrv

!***********************************************************************




      subroutine locate_segstart(cell_grid_num, cell,
     &     beglat, beglon, endlat, endlon, offset,
     &     srch_grid_num, cont_cell, lboundary, edgeid)

!-----------------------------------------------------------------------
!
!     Find the cell containing the given point
!
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!
!     intent(in): 
!
!-----------------------------------------------------------------------

      real (SCRIP_r8), intent(in) :: 
     &     beglat, beglon,      ! beginning and end points of segment
     &     endlat, endlon       ! on which the point to be located lies

      real (SCRIP_r8), intent(in) ::
     &     offset               ! Offset to calculate the search point

      integer (SCRIP_i4), intent(in) ::
     &     cell,                ! Cell from which point originates
                                ! Point will be on boundary of orig_cell
     &     cell_grid_num        ! Index of grid to which cell belongs

      integer (SCRIP_i4), intent(in) ::
     &     srch_grid_num        ! num indicating if we are locating a 
                                ! grid1 point in a cell of grid2 (num=2)
                                ! or a grid2 point in a cell of grid1 
                                ! (num=1)

!-----------------------------------------------------------------------
!
!     intent(out): 
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4), intent(out) ::
     &     cont_cell            ! grid cell containing this point

      logical (SCRIP_logical), intent(out) ::
     &     lboundary            ! flag points that lie on the boundary 
                                ! of the cell

      integer (SCRIP_i4), intent(out) ::
     &     edgeid               ! if point is on boundary, which local
                                ! edge is it on? (0 otherwise)

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4) :: i, j, k, n, ic
      integer (SCRIP_i4) :: whichpole, srch_cell_add,
     &                      grid1_add, grid2_add, min_add, max_add

      real (SCRIP_r8), dimension(:), allocatable ::
     &     cell_corner_x, cell_corner_y

      logical (SCRIP_logical) :: inpoly, latlon      

      real (SCRIP_r8) ::
     &     vec1_x, vec1_y, vec1_lenx, vec1_lat, vec1_lon, vec1_len,
     &     begx, begy, endx, endy, ptx, pty, rns, pi4, ptlat, ptlon,
     &     lat, lon, cell_center_x, cell_center_y

!-----------------------------------------------------------------------
!
!     initialize defaults, flags, etc.
!
!-----------------------------------------------------------------------

      lboundary = .false.
      edgeid = 0
      cont_cell = 0
      

      if (cell /= last_cell_locate_segstart .or. 
     &     cell_grid_num /= last_cell_grid_num_locate_segstart 
     &     .or. srch_grid_num /= last_srch_grid_num_locate_segstart)
     &      then

         last_cell_locate_segstart = cell
         last_cell_grid_num_locate_segstart = cell_grid_num
         last_srch_grid_num_locate_segstart = srch_grid_num

         if (first_call_locate_segstart) then
            first_call_locate_segstart = .false.
            last_cell_locate_segstart = 0
            last_cell_grid_num_locate_segstart = 0
            last_srch_grid_num_locate_segstart = 0
            num_srch_cells_locate_segstart = 0
         else
            if (num_srch_cells_locate_segstart .gt. 0) then
               deallocate(srch_add_locate_segstart,
     &         srch_corner_lat_locate_segstart,
     &         srch_corner_lon_locate_segstart,
     &         srch_center_lat_locate_segstart, 
     &         srch_center_lon_locate_segstart)
            endif
         endif

         call get_srch_cells(cell, cell_grid_num, srch_grid_num, 
     &        num_srch_cells_locate_segstart, srch_add_locate_segstart, 
     &        srch_corners_locate_segstart, 
     &        srch_corner_lat_locate_segstart, 
     &        srch_corner_lon_locate_segstart,
     &        srch_center_lat_locate_segstart, 
     &        srch_center_lon_locate_segstart)

      endif

      if (num_srch_cells_locate_segstart == 0) return


      do ic=1,num_srch_cells_locate_segstart
         
         srch_cell_add = srch_add_locate_segstart(ic)
         


         !**** CAN WE ACCOMPLISH THE FOLLOWING THROUGH A SUBROUTINE 
         !**** CALLED SEGSTART_INCELL ?? 


    !*** IF POINT IS IN POLAR REGION, CHECK IN A TRANSFORMED SPACE
    !*** HOWEVER, POINTS THAT ARE PRACTICALLY AT THE POLE CANNOT
    !*** BE CORRECTLY LOCATED THIS WAY BECAUSE THE POLE IS A SINGULARITY
    !*** AND CONTAINMENT IN ANY CELL INCIDENT ON THE POLE WILL GIVE US A 
    !*** POSITIVE ANSWER. FOR THESE POINTS REVERT TO THE LATLON SPACE
    !***



         vec1_lat = endlat-beglat
         vec1_lon = endlon-beglon
         if (vec1_lon > pi) then
            vec1_lon = vec1_lon - pi2
         else if (vec1_lon < -pi) then
            vec1_lon = vec1_lon + pi2
         endif
         vec1_len = sqrt(vec1_lat*vec1_lat+vec1_lon*vec1_lon)
         vec1_lat = vec1_lat/vec1_len
         vec1_lon = vec1_lon/vec1_len
         
         ptlat = beglat + offset*vec1_lat
         ptlon = beglon + offset*vec1_lon
            

         if ((ptlat .gt. north_thresh .and. abs(ptlat-pih) .ge. 0.001) 
     &        .or. 
     &        (ptlat .lt. south_thresh .and. abs(ptlat+pih) .ge. 0.001))
     &        then

            if (ptlat > zero) then
               pi4 = quart*pi
               rns = one
            else
               pi4 = -quart*pi
               rns = -one
            endif

        
            begx = rns*two*sin(pi4 - half*beglat)*cos(beglon)
            begy =     two*sin(pi4 - half*beglat)*sin(beglon)
            endx = rns*two*sin(pi4 - half*endlat)*cos(endlon)
            endy =     two*sin(pi4 - half*endlat)*sin(endlon)
        
            vec1_x = endx-begx
            vec1_y = endy-begy

            vec1_lenx = sqrt(vec1_x*vec1_x + vec1_y*vec1_y)
            vec1_x = vec1_x/vec1_lenx
            vec1_y = vec1_y/vec1_lenx


            !*** Must calculate ptx and pty as an offset on straight 
            !*** line in polar space rather than calculating it on a
            !*** straight line in latlon space an offset point in latlon
            !***  space will be off the straight line in polar space

            ptx = begx + offset*vec1_x
            pty = begy + offset*vec1_y

            latlon = .false.

         ! Since we want greater fidelity for locating the points 
         ! we send in the mid-points of the polygon edges too
         ! BUT THAT MAKES THE POLYGON NON-CONVEX SOMETIMES AND 
         ! THE CROSS-PRODUCT CHECK FAILS. SO USE CODE TO CHECK GENERAL
         ! POLYGONS


            allocate(cell_corner_x(npseg*srch_corners_locate_segstart),
     &           cell_corner_y(npseg*srch_corners_locate_segstart))


            k = 0
            do i = srch_corners_locate_segstart, 1, -1
               k = k+1
               lat = srch_corner_lat_locate_segstart(i,ic)
               lon = srch_corner_lon_locate_segstart(i,ic)
               cell_corner_x(k) = rns*two*sin(pi4-half*lat)*cos(lon)
               cell_corner_y(k) =     two*sin(pi4-half*lat)*sin(lon)
           
               j = i-1
               if (j .eq. 0) j = srch_corners_locate_segstart ! how do 
                                              ! we do (j-1+n)%n in F90?
           
               vec1_lat = srch_corner_lat_locate_segstart(j,ic)
     &              -srch_corner_lat_locate_segstart(i,ic)
               vec1_lon = srch_corner_lon_locate_segstart(j,ic)
     &              -srch_corner_lon_locate_segstart(i,ic)
               if (vec1_lon > pi) then
                  vec1_lon = vec1_lon - pi2
               else if (vec1_lon < -pi) then
                  vec1_lon = vec1_lon + pi2
               endif

               do j = 1, npseg-1
                  k = k+1
                  lat = srch_corner_lat_locate_segstart(i,ic)
     &               + j*vec1_lat/npseg
                  lon = srch_corner_lon_locate_segstart(i,ic)
     &               + j*vec1_lon/npseg
                  cell_corner_x(k) = rns*two*sin(pi4-half*lat)*cos(lon)
                  cell_corner_y(k) =     two*sin(pi4-half*lat)*sin(lon)
               enddo
            enddo


            call ptinpolygen2(ptx, pty, k, cell_corner_x,
     &           cell_corner_y, latlon, inpoly, lboundary, edgeid)

            if (lboundary) then
               edgeid = (edgeid-1)/npseg + 1 ! convert from index in 
                                       ! multi-segmented to regular cell
            endif

            deallocate(cell_corner_x, cell_corner_y)

         else

            latlon = .true.
            
            whichpole = 0
            if (srch_grid_num .eq. 1 .and. 
     &           srch_cell_add .eq. grid1_spole_cell) then
               
               whichpole = -1   ! S pole
               call ptinpolarpoly(ptlat, ptlon, 
     &              srch_corners_locate_segstart,
     &              srch_corner_lat_locate_segstart(:,ic), 
     &              srch_corner_lon_locate_segstart(:,ic),
     &              latlon, whichpole, inpoly, lboundary, edgeid)
               
            else if (srch_grid_num .eq. 1 .and.
     &              srch_cell_add .eq. grid1_npole_cell) then
               
               whichpole =  1   ! N pole
               call ptinpolarpoly(ptlat, ptlon, 
     &              srch_corners_locate_segstart,
     &              srch_corner_lat_locate_segstart(:,ic), 
     &              srch_corner_lon_locate_segstart(:,ic),
     &              latlon, whichpole, inpoly, lboundary, edgeid)
               
            else if (srch_grid_num .eq. 2 .and.
     &              srch_cell_add .eq. grid2_spole_cell) then
               
               whichpole = -1   ! S pole
               call ptinpolarpoly(ptlat, ptlon, 
     &              srch_corners_locate_segstart,
     &              srch_corner_lat_locate_segstart(:,ic), 
     &              srch_corner_lon_locate_segstart(:,ic),
     &              latlon, whichpole, inpoly, lboundary, edgeid)
               
            else if (srch_grid_num .eq. 2 .and.
     &              srch_cell_add .eq. grid2_npole_cell) then
               
               whichpole =  1   ! N pole
               call ptinpolarpoly(ptlat, ptlon, 
     &              srch_corners_locate_segstart,
     &              srch_corner_lat_locate_segstart(:,ic), 
     &              srch_corner_lon_locate_segstart(:,ic),
     &              latlon, whichpole, inpoly, lboundary, edgeid)
               
            else
            
      !***
      !***  General cell
      !***
               
               call ptinpoly(ptlat, ptlon, srch_corners_locate_segstart,
     &              srch_corner_lat_locate_segstart(:,ic), 
     &              srch_corner_lon_locate_segstart(:,ic),
     &              latlon, inpoly, lboundary, edgeid)
               
            endif
            
         endif

         if (inpoly) then
            cont_cell = srch_cell_add
            exit
         endif         
            
      end do

      return

!----------------------------------------------------------------------

      end subroutine locate_segstart

!**********************************************************************




!**********************************************************************

      subroutine locate_point(ptlat, ptlon, cell, cell_grid_num,
     &     srch_grid_num, cont_cell, lboundary, edgeid)

!-----------------------------------------------------------------------
!
!     Find the cell containing the given point
!
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!
!     intent(in): 
!
!-----------------------------------------------------------------------

      real (SCRIP_r8), intent(in) :: 
     &     ptlat, ptlon         ! Point to locate

      integer (SCRIP_i4), intent(in) ::
     &     cell,                ! Cell from which point originates
                                ! Point will be on boundary of orig_cell
     &     cell_grid_num        ! Index of grid to which cell belongs

      integer (SCRIP_i4), intent(in) ::
     &     srch_grid_num        ! num indicating if we are locating a 
                                ! grid1 point in a cell of grid2 (num=2)
                                ! or a grid2 point in a cell of grid1
                                ! (num=1)

!-----------------------------------------------------------------------
!
!     intent(out): 
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4), intent(out) ::
     &     cont_cell            ! grid cell containing this point

      logical (SCRIP_logical), intent(out) ::
     &     lboundary            ! flag points that lie on the boundary 
                                ! of the cell

      integer (SCRIP_i4), intent(out) ::
     &     edgeid               ! if point is on boundary, which local
                                ! edge is it on? (0 otherwise)

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4) :: i, j, n, ic
      integer (SCRIP_i4) :: whichpole, srch_cell_add,
     &                      grid1_add, grid2_add, min_add, max_add,
     &                      previdx, nextidx, pcorner, corner, 
     &                      ncorners, nalloc

       real (SCRIP_r8), dimension(:), allocatable ::
     &     cell_corner_lat,
     &     cell_corner_lon

      real (SCRIP_r8) ::
     &     prevlon,
     &     nextlon,
     &     polelat,
     &     cell_center_lat,
     &     cell_center_lon


      logical (SCRIP_logical) :: inpoly, latlon      
      logical (SCRIP_logical) :: test

!-----------------------------------------------------------------------
!
!     initialize defaults, flags, etc.
!
!-----------------------------------------------------------------------

      lboundary = .false.
      edgeid = 0
      cont_cell = 0

      if (cell /= last_cell_locate_point .or. cell_grid_num /=
     &      last_cell_grid_num_locate_point 
     &     .or. srch_grid_num /= last_srch_grid_num_locate_point) then

         last_cell_locate_point = cell
         last_cell_grid_num_locate_point = cell_grid_num
         last_srch_grid_num_locate_point = srch_grid_num

         if (first_call_locate_point) then
            first_call_locate_point = .false.
            last_cell_locate_point = 0
            last_cell_grid_num_locate_point = 0
            last_srch_grid_num_locate_point = 0
            num_srch_cell_locate_points = 0
         else
            if (num_srch_cell_locate_points .gt. 0) then
               deallocate(srch_add_locate_point,
     &         srch_corner_lat_locate_point,
     &         srch_corner_lon_locate_point)
            endif
         endif

         call get_srch_cells(cell, cell_grid_num, srch_grid_num, 
     &        num_srch_cell_locate_points, srch_add_locate_point,
     &        srch_corners_locate_point, 
     &        srch_corner_lat_locate_point,srch_corner_lon_locate_point,
     &        srch_center_lat_locate_point,srch_center_lon_locate_point)

      endif

      if (num_srch_cell_locate_points == 0) return


      ncorners = srch_corners_locate_point
      nalloc = ncorners+2
      allocate(cell_corner_lat(nalloc),
     &     cell_corner_lon(nalloc))


      do ic=1,num_srch_cell_locate_points
         
         srch_cell_add = srch_add_locate_point(ic)

         do i = 1, ncorners
            cell_corner_lat(i) = srch_corner_lat_locate_point(i,ic)
            cell_corner_lon(i) = srch_corner_lon_locate_point(i,ic)
         enddo
         
         cell_center_lat = srch_center_lat_locate_point(ic)
         cell_center_lon = srch_center_lon_locate_point(ic)

!         if ((srch_grid_num .eq. 1 .and. 
!     &        (special_polar_cell1(srch_cell_add))) .or.
!     &        (srch_grid_num .eq. 2 .and.
!     &        (special_polar_cell2(srch_cell_add)))) then
!
! Modified by MD
         test=.false.
         if (srch_grid_num .eq. 1) then
           if (special_polar_cell1(srch_cell_add)) then
             test=.true.
           endif
         else
           if (special_polar_cell2(srch_cell_add)) then
             test=.true.
           endif
         endif
         if (test) then
            call modify_polar_cell(ncorners, nalloc, cell_corner_lat,
     &           cell_corner_lon)

         endif

         call ptincell(ptlat, ptlon, srch_cell_add, ncorners,
     &        cell_corner_lat, cell_corner_lon,
     &        cell_center_lat, cell_center_lon,
     &        srch_grid_num, inpoly, lboundary, edgeid)

         
         if (inpoly) then
            cont_cell = srch_cell_add
            exit
         endif         

         ncorners = srch_corners_locate_point ! reset it for other srch 
                                              !cells
      end do

!----------------------------------------------------------------------

      end subroutine locate_point

!**********************************************************************



!**********************************************************************

      subroutine ptincell(ptlat, ptlon, cell_add, ncorners, 
     &     cell_corner_lat, cell_corner_lon,
     &     cell_center_lat, cell_center_lon,
     &     cell_grid_id, inpoly, lboundary, edgeid)

!----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!
!     intent(in): 
!
!-----------------------------------------------------------------------

      real (SCRIP_r8), intent(in) :: 
     &     ptlat, ptlon         ! Point to locate

      integer (SCRIP_i4), intent(in) ::
     &     cell_add             ! ID of cell

      integer (SCRIP_i4), intent(in) ::
     &     ncorners 

      real (SCRIP_r8), dimension(ncorners), intent(in) ::
     &     cell_corner_lat, cell_corner_lon

      real (SCRIP_r8), intent(in) ::
     &     cell_center_lat, 
     &     cell_center_lon

      integer (SCRIP_i4), intent(in) ::
     &     cell_grid_id     ! num indicating if we are locating a grid1
                            ! point in a cell of grid2 (num = 2) or 
                            ! a grid2 point in a cell of grid1 (num = 1)


!-----------------------------------------------------------------------
!
!     intent(out): 
!
!-----------------------------------------------------------------------

      logical (SCRIP_logical), intent(out) ::
     &     inpoly               ! is point in polygon?

      logical (SCRIP_logical), intent(out) ::
     &     lboundary            ! flag points that lie on the boundary 
                                ! of the cell

      integer (SCRIP_i4), intent(out) ::
     &     edgeid               ! if point is on boundary, which local
                                ! edge is it on? (0 otherwise)
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4) :: i, j, k, ic
      integer (SCRIP_i4) :: whichpole

      real (SCRIP_r8) :: rns, pi4, ptx, pty, lat, lon,
     &     cell_center_x, cell_center_y, vec1_lat, vec1_lon

      logical (SCRIP_logical) ::
     &     latlon

      real (kind=SCRIP_r8), dimension(npseg*ncorners) ::
     &     cell_corner_x,  ! x of each corner of cell
     &     cell_corner_y   ! y of each corner of cell

!----------------------------------------------------------------------

      edgeid = 0


    !*** IF POINTS ARE ABOVE THE THRESHOLD, CHECK THEM IN A TRANSFORMED
    !*** SPACE
    !*** HOWEVER, POINTS THAT ARE PRACTICALLY AT THE POLE CANNOT
    !*** BE CORRECTLY LOCATED THIS WAY BECAUSE THE POLE IS A SINGULARITY
    !*** AND CONTAINMENT IN ANY CELL INCIDENT ON THE POLE WILL GIVE US A 
    !*** POSITIVE ANSWER. FOR THESE POINTS REVERT TO THE LATLON SPACE
    !***

      if ((ptlat .gt. north_thresh .and. abs(ptlat-pih) .ge. 0.001) .or.
     &    (ptlat .lt. south_thresh .and. abs(ptlat+pih) .ge. 0.001)) 
     &     then

         if (ptlat > zero) then
            pi4 = quart*pi
            rns = one
         else
            pi4 = -quart*pi
            rns = -one
         endif

         ptx = rns*two*sin(pi4 - half*ptlat)*cos(ptlon)
         pty =     two*sin(pi4 - half*ptlat)*sin(ptlon)
         
         latlon = .false.

         ! Since we want greater fidelity for locating the points 
         ! we send in the mid-points of the polygon edges too
         ! BUT THAT MAKES THE POLYGON NON-CONVEX SOMETIMES AND 
         ! THE CROSS-PRODUCT CHECK FAILS. SO USE CODE TO CHECK GENERAL
         ! POLYGONS


         k = 0
         do i = ncorners, 1, -1
            k = k+1
            lat = cell_corner_lat(i)
            lon = cell_corner_lon(i)
            cell_corner_x(k) = rns*two*sin(pi4-half*lat)*cos(lon)
            cell_corner_y(k) =     two*sin(pi4-half*lat)*sin(lon)
            
            j = i-1
            if (j .eq. 0) j = ncorners  ! how do we do (j-1+n)%n in F90?

            vec1_lat = cell_corner_lat(j)-cell_corner_lat(i)
            vec1_lon = cell_corner_lon(j)-cell_corner_lon(i)
            if (vec1_lon > pi) then
               vec1_lon = vec1_lon - pi2
            else if (vec1_lon < -pi) then
               vec1_lon = vec1_lon + pi2
            endif

            do j = 1, npseg-1
               k = k+1
               lat = cell_corner_lat(i) + j*vec1_lat/npseg
               lon = cell_corner_lon(i) + j*vec1_lon/npseg
               cell_corner_x(k) = rns*two*sin(pi4-half*lat)*cos(lon)
               cell_corner_y(k) =     two*sin(pi4-half*lat)*sin(lon)
            enddo
         enddo

         !*** cell is so non-convex that no feasible center exists
         !*** we have to fall back on a different algorithm

         call ptinpolygen2(ptx, pty, k, cell_corner_x,
     &        cell_corner_y, latlon, inpoly, lboundary, edgeid)

         if (lboundary) then
            edgeid = (edgeid-1)/npseg + 1 ! convert from index in 
                                          ! multi-segmented cell to 
                                          ! regular cell
         endif
      else

         latlon = .true.
         
         whichpole = 0
         if (cell_grid_id .eq. 1 .and. 
     &        cell_add .eq. grid1_spole_cell) then
            
            whichpole = -1      ! S pole
            call ptinpolarpoly(ptlat, ptlon, ncorners,
     &           cell_corner_lat, cell_corner_lon,
     &           latlon, whichpole, inpoly, lboundary, edgeid)
            
         else if (cell_grid_id .eq. 1 .and.
     &           cell_add .eq. grid1_npole_cell) then
            
            whichpole =  1      ! N pole
            call ptinpolarpoly(ptlat, ptlon, ncorners,
     &           cell_corner_lat, cell_corner_lon,
     &           latlon, whichpole, inpoly, lboundary, edgeid)
            
         else if (cell_grid_id .eq. 2 .and.
     &           cell_add .eq. grid2_spole_cell) then
            
            whichpole = -1      ! S pole
            call ptinpolarpoly(ptlat, ptlon, ncorners,
     &           cell_corner_lat, cell_corner_lon,
     &           latlon, whichpole, inpoly, lboundary, edgeid)
            
         else if (cell_grid_id .eq. 2 .and.
     &           cell_add .eq. grid2_npole_cell) then
            
            whichpole =  1      ! N pole
            call ptinpolarpoly(ptlat, ptlon, ncorners,
     &           cell_corner_lat, cell_corner_lon,
     &           latlon, whichpole, inpoly, lboundary, edgeid)
            
         else
            
      !***
      !***  General cell
      !***
               
            call ptinpoly(ptlat, ptlon, ncorners, 
     &           cell_corner_lat, cell_corner_lon,
     &           latlon, inpoly, lboundary, edgeid)

         endif

      endif

      return

!----------------------------------------------------------------------

      end subroutine ptincell

!**********************************************************************

!**********************************************************************

      subroutine ptinpoly(ptx, pty, ncorners, cell_corner_x,
     &     cell_corner_y, latlon, inpoly, lboundary, edgeid)

!----------------------------------------------------------------------
!
!     Check if point is in (convex) polygonal cell 
!
!----------------------------------------------------------------------

      implicit none

!----------------------------------------------------------------------
!
!     Input arguments
!
!----------------------------------------------------------------------

      real (SCRIP_r8), intent(in) ::
     &     ptx, pty             ! Point to check

      integer (SCRIP_i4), intent(in) ::
     &     ncorners             ! Number of polygon corners

      real (SCRIP_r8), dimension(ncorners), intent(in) ::
     &     cell_corner_x,       ! Coordinates of cell corners
     &     cell_corner_y        ! Could be x-y or lat-lon or ...

      logical (SCRIP_logical), intent(in) ::
     &     latlon               ! Are coordinates in latlon space?

!----------------------------------------------------------------------
!
!     Output arguments
!
!----------------------------------------------------------------------

      logical (SCRIP_logical), intent(out) ::
     &     inpoly             ! Is point in the polygon?

      logical (SCRIP_logical), intent(out) ::
     &     lboundary          ! Is point on the boundary of the polygon?

      integer (SCRIP_i4), intent(out) ::
     &     edgeid             ! if point is on boundary, which local
                              ! edge is it on? (0 otherwise)

!----------------------------------------------------------------------
!
!     Local variables
!
!----------------------------------------------------------------------

      integer (SCRIP_i4) :: n, next_n

      real (SCRIP_r8) :: x1, y1, x2, y2, vec1_x, vec1_y, vec2_x, vec2_y,
     &     cross_product, minlon, maxlon, ptx_loc, pty_loc

      real (SCRIP_r8), dimension(ncorners) ::
     &     cell_corner_lat_loc, cell_corner_lon_loc


      !***********************************************************
      !*** We should just remove the latlon argument since that is
      !*** the only coordinate system we are using it for
      !***********************************************************      
      

      !***
      !*** here we take the cross product of the vector making 
      !*** up each cell side with the vector formed by the vertex
      !*** and search point.  if all the cross products are 
      !*** positive, the point is contained in the cell.
      !***

      inpoly = .false.
      lboundary = .false.
      edgeid = 0
        
      if (.not. latlon) then
         
         do n = 1, ncorners
            next_n = MOD(n,ncorners) + 1
            
            x1 = cell_corner_x(n)
            y1 = cell_corner_y(n)
            x2 = cell_corner_x(next_n)
            y2 = cell_corner_y(next_n)
            
            vec1_x = x2 - x1
            vec1_y = y2 - y1
            vec2_x = ptx - x1
            vec2_y = pty - y1
                        
            cross_product = vec1_y*vec2_x - vec2_y*vec1_x
            
            !***
            !***   if the cross product for a side is zero, the point 
            !***   lies exactly on the side or the side is degenerate
            !***   (zero length).  if degenerate, set the cross 
            !***   product to a positive number.  
            !***

            if (abs(cross_product) < tiny) then
               if (vec1_x*vec1_x + vec1_y*vec1_y .le. tiny*tiny) then
                  cross_product = one
               else
                  lboundary = .true.
                  edgeid = n
               endif
            else

            !***
            !*** if cross product is less than zero, this cell
            !*** doesn't work
            !***
            !*** Should we say "if (cp < zero .and. abs(cp) > tiny)" ?

               if (cross_product < zero) then
                  inpoly = .false.
                  lboundary = .false.
                  return
               endif
            endif
            
         end do

      else

         !*** Checking in latlon space
         !*** If the grid cell coordinates spans more than pi radians
         !*** transform the coordinates so that they don't

         cell_corner_lat_loc = cell_corner_x
         cell_corner_lon_loc = cell_corner_y
            
         minlon =  9999.0
         maxlon = -9999.0
         do n = 1, ncorners
            if (cell_corner_lon_loc(n) < minlon) then
               minlon = cell_corner_lon_loc(n)
            endif
            if (cell_corner_lon_loc(n) > maxlon) then
               maxlon = cell_corner_lon_loc(n)
            endif
         enddo

         if (maxlon-minlon > pi) then

            do n = 1, ncorners
               if (cell_corner_lon_loc(n)-minlon > pi) then
                  cell_corner_lon_loc(n) = cell_corner_lon_loc(n)-pi2
               endif
            enddo

         endif

         ptx_loc = ptx
         pty_loc = pty
         if (pty_loc - minlon > pi) then
            pty_loc = pty_loc - pi2
         else if (pty_loc - minlon < -pi) then
            pty_loc = pty_loc + pi2
         endif


         do n = 1, ncorners
            next_n = MOD(n,ncorners) + 1
            
            x1 = cell_corner_lat_loc(n)
            y1 = cell_corner_lon_loc(n)
            x2 = cell_corner_lat_loc(next_n)
            y2 = cell_corner_lon_loc(next_n)
            
            vec1_x = x2 - x1
            vec1_y = y2 - y1
            vec2_x = ptx_loc - x1
            vec2_y = pty_loc - y1
                        
            cross_product = vec1_y*vec2_x - vec2_y*vec1_x
            
            !***
            !***   if the cross product for a side is zero, the point 
            !***   lies exactly on the side or the side is degenerate
            !***   (zero length).  if degenerate, set the cross 
            !***   product to a positive number.  
            !***

            if (abs(cross_product) < tiny) then
               if (vec1_x*vec1_x + vec1_y*vec1_y .le. tiny*tiny) then
                  cross_product = one
               else
                  lboundary = .true.
                  edgeid = n
               endif
            else

            !***
            !*** if cross product is less than zero, this cell
            !*** doesn't work
            !***
            !*** Should we say "if (cp < zero .and. abs(cp) > tiny)" ?

               if (cross_product < zero) then
                  inpoly = .false.
                  lboundary = .false.
                  return
               endif
            endif
            
         end do

      endif
      !***
      !*** if cross products all positive, we found the location
      !***

      inpoly = .true.
      return

!----------------------------------------------------------------------

      end subroutine ptinpoly

!**********************************************************************



      subroutine ptinpolarpoly(ptx, pty, ncorners, cell_corner_x,
     &     cell_corner_y, latlon, whichpole, inpoly, lboundary, edgeid)

!----------------------------------------------------------------------
!
!     Check if point is in polygonal cell overlapping the pole
!     Cannot check the containment as is in latlon space - We have 
!     to check by connecting each edge of the polygon to the pole
!     and check containment in the resulting quadrilateral in latlon 
!     space
!     The cell can be non-convex as long as the pole is 'visible' to
!     all the edges of the polygon, i.e., we can connect the pole to 
!     each edge of the polygon and form a triangle with positive area
!
!----------------------------------------------------------------------

      implicit none

!----------------------------------------------------------------------
!
!     Input arguments
!
!----------------------------------------------------------------------

      real (SCRIP_r8), intent(in) ::
     &     ptx, pty             ! Point to check

      integer (SCRIP_i4), intent(in) ::
     &     ncorners             ! Number of polygon corners

      real (SCRIP_r8), dimension(ncorners), intent(in) ::
     &     cell_corner_x,       ! Coordinates of cell corners
     &     cell_corner_y        ! Could be x-y or lat-lon or ...

      logical (SCRIP_logical), intent(in) ::
     &     latlon               ! Are coordinates in latlon space?

      integer (SCRIP_i4), intent(in) ::
     &     whichpole            ! South or North pole      

!----------------------------------------------------------------------
!
!     Output arguments
!
!----------------------------------------------------------------------

      logical (SCRIP_logical), intent(out) ::
     &     inpoly             ! Is point in the polygon?

      logical (SCRIP_logical), intent(out) ::
     &     lboundary          ! Is point on the boundary of the polygon?

      integer (SCRIP_i4), intent(out) ::
     &     edgeid             ! if point is on boundary, which local
                              ! edge is it on? (0 otherwise)

!----------------------------------------------------------------------
!
!     Local variables
!
!----------------------------------------------------------------------

      integer (SCRIP_i4) :: n, next_n, ledgeid

      real (SCRIP_r8), dimension(4) ::
     &     pquad_corner_x,       ! Coordinates of polar quad
     &     pquad_corner_y        

      real (SCRIP_r8) :: x1, y1, x2, y2, vec1_x, vec1_y, vec2_x, vec2_y,
     &     cross_product, pole_lat

      pole_lat = whichpole*pih

      !***
      !*** This is a polygon that overlaps the pole
      !*** A normal point in polygon check could fail
      !*** So, with each edge of the polygon form a quadrilateral 
      !*** in latlon space using the polar latitude and the longitude
      !*** values of the endpoints of the edge. Then check containment
      !*** of the point in this quadrilateral
      !***
      
      inpoly = .false.
      lboundary = .false.
      
      do n = 1, ncorners
         next_n = MOD(n,ncorners) + 1

         pquad_corner_x(1) = cell_corner_x(n)
         pquad_corner_y(1) = cell_corner_y(n)
         pquad_corner_x(2) = cell_corner_x(next_n)
         pquad_corner_y(2) = cell_corner_y(next_n)
         pquad_corner_x(3) = pole_lat
         pquad_corner_y(3) = cell_corner_y(next_n)
         pquad_corner_x(4) = pole_lat
         pquad_corner_y(4) = cell_corner_y(n)

         
         call ptinpoly(ptx,pty,4,pquad_corner_x,pquad_corner_y, 
     &        latlon,inpoly,lboundary, ledgeid)

         if (inpoly) then

            if (lboundary) then

               !***
               !*** Check to see if the lboundary flag is being
               !*** triggered by the outer edge of the polygon or 
               !*** by one of the artificial internal edges            
               !***

               vec1_x = pquad_corner_x(2) - pquad_corner_x(1)
               vec1_y = pquad_corner_y(2) - pquad_corner_y(1)
               vec2_x = ptx - pquad_corner_x(1)
               vec2_y = pty - pquad_corner_y(1)


               if (latlon) then

                  !***
                  !*** check for 0,2pi crossings
                  !***

                  if (vec1_y >  pi) vec1_y = vec1_y - pi2
                  if (vec1_y < -pi) vec1_y = vec1_y + pi2
                  if (vec2_y >  pi) vec2_y = vec2_y - pi2
                  if (vec2_y < -pi) vec2_y = vec2_y + pi2
               
               endif

               cross_product = vec1_y*vec2_x - vec2_y*vec1_x

               !***
               !***   if the cross product for a side is zero, the point
               !***   lies exactly on the side or the side is degenerate
               !***   (zero length).  if degenerate, set the cross 
               !***   product to a positive number.  
               !***

               if (abs(cross_product) < tiny) then
                  if (vec1_x .eq. zero .and. vec1_y .eq. zero) then
                     cross_product = one
                     lboundary = .false.
                  else
                     edgeid = n
                     lboundary = .true.
                  endif
               else
                  lboundary = .false.
               endif
            endif               ! if (lboundary)

            return              ! pt in polygon

         endif                  ! if (inpoly)

      end do

      return                    ! pt outside polygon

!----------------------------------------------------------------------

      end subroutine ptinpolarpoly

!**********************************************************************



      subroutine ptinpolygen(ptx, pty, ncorners, cell_corner_x,
     &     cell_corner_y, cell_center_x, cell_center_y,
     &     latlon, inpoly, lboundary, edgeid)

!----------------------------------------------------------------------
!
!     Check if point is in general (convex or mildly non-convex) 
!     polygonal cell by connecting each edge of the polygon to a
!     a central point (average of vertices) and check containment in 
!     the resulting triangle
!
!     The cell can be non-convex as long as the 'center' is 'visible' to
!     all the edges of the polygon, i.e., we can connect the 'center' to
!     each edge of the polygon and form a triangle with positive area
!
!----------------------------------------------------------------------

      implicit none

!----------------------------------------------------------------------
!
!     Input arguments
!
!----------------------------------------------------------------------

      real (SCRIP_r8), intent(in) ::
     &     ptx, pty             ! Point to check

      integer (SCRIP_i4), intent(in) ::
     &     ncorners             ! Number of polygon corners

      real (SCRIP_r8), dimension(ncorners), intent(in) ::
     &     cell_corner_x,       ! Coordinates of cell corners
     &     cell_corner_y        ! Could be x-y or lat-lon or ...

      real (SCRIP_r8), intent(in) ::
     &     cell_center_x,
     &     cell_center_y

      logical (SCRIP_logical), intent(in) ::
     &     latlon               ! Are coordinates in latlon space?

!----------------------------------------------------------------------
!
!     Output arguments
!
!----------------------------------------------------------------------

      logical (SCRIP_logical), intent(out) ::
     &     inpoly             ! Is point in the polygon?

      logical (SCRIP_logical), intent(out) ::
     &     lboundary          ! Is point on the boundary of the polygon?

      integer (SCRIP_i4), intent(out) ::
     &     edgeid             ! if point is on boundary, which local
                              ! edge is it on? (0 otherwise)

!----------------------------------------------------------------------
!
!     Local variables
!
!----------------------------------------------------------------------

      integer (SCRIP_i4) :: n, next_n, ledgeid

      real (SCRIP_r8), dimension(3) ::
     &     tri_corner_x,       ! Coordinates of triangle
     &     tri_corner_y        

      real (SCRIP_r8) :: x1, y1, x2, y2, vec1_x, vec1_y, vec2_x, vec2_y,
     &     cross_product


      !***
      !*** So, with each edge of the polygon form a triangle
      !*** by connecting a 'central' point to the endpoints of 
      !*** the edge. Then check containment of the point in this tri
      !***
      
      inpoly = .false.
      lboundary = .false.
      
      do n = 1, ncorners
         next_n = MOD(n,ncorners) + 1

         tri_corner_x(1) = cell_corner_x(n)
         tri_corner_y(1) = cell_corner_y(n)
         tri_corner_x(2) = cell_corner_x(next_n)
         tri_corner_y(2) = cell_corner_y(next_n)
         tri_corner_x(3) = cell_center_x
         tri_corner_y(3) = cell_center_y

         vec1_x = tri_corner_x(2) - tri_corner_x(1)
         vec1_y = tri_corner_y(2) - tri_corner_y(1)

         !*** Skip triangles arising from degenerate edges

         if (vec1_x*vec1_x+vec1_y*vec1_y .le. tiny*tiny) cycle
         
         call ptinpoly(ptx,pty,3,tri_corner_x,tri_corner_y, 
     &        latlon,inpoly,lboundary, ledgeid)

         if (inpoly) then

            if (lboundary) then

               !***
               !*** Check to see if the lboundary flag is being
               !*** triggered by the outer edge of the polygon or 
               !*** by one of the artificial internal edges            
               !***

               vec2_x = ptx - tri_corner_x(1)
               vec2_y = pty - tri_corner_y(1)


               if (latlon) then

                  !***
                  !*** check for 0,2pi crossings
                  !***

                  if (vec1_y >  pi) vec1_y = vec1_y - pi2
                  if (vec1_y < -pi) vec1_y = vec1_y + pi2
                  if (vec2_y >  pi) vec2_y = vec2_y - pi2
                  if (vec2_y < -pi) vec2_y = vec2_y + pi2
               
               endif

               cross_product = vec1_y*vec2_x - vec2_y*vec1_x

               !***
               !***   if the cross product for a side is zero, the point 
               !***   lies exactly on the side or the side is degenerate
               !***   (zero length).  if degenerate, set the cross 
               !***   product to a positive number.  
               !***

               if (abs(cross_product) < tiny) then
                  if (vec1_x*vec1_x+vec1_y*vec1_y .le. tiny*tiny) then
                     cross_product = one
                     lboundary = .false.
                  else
                     edgeid = n
                     lboundary = .true.
                  endif
               else
                  lboundary = .false.
               endif
            endif               ! if (lboundary)

            return              ! pt in polygon

         endif                  ! if (inpoly)

      end do

      return                    ! pt outside polygon

!----------------------------------------------------------------------

      end subroutine ptinpolygen

!**********************************************************************



      subroutine ptinpolygen2(ptx, pty, ncorners, cell_corner_x,
     &     cell_corner_y, latlon, inpoly, lboundary, edgeid)

!----------------------------------------------------------------------
!
!     Check if point is in general (convex or mildly non-convex) 
!     polygonal cell by connecting each edge of the polygon to a
!     a central point (average of vertices) and check containment in 
!     the resulting triangle
!
!     The cell can be non-convex as long as the 'center' is 'visible' to
!     all the edges of the polygon, i.e., we can connect the 'center' to 
!     each edge of the polygon and form a triangle with positive area
!
!----------------------------------------------------------------------

      implicit none

!----------------------------------------------------------------------
!
!     Input arguments
!
!----------------------------------------------------------------------

      real (SCRIP_r8), intent(in) ::
     &     ptx, pty             ! Point to check

      integer (SCRIP_i4), intent(in) ::
     &     ncorners             ! Number of polygon corners

      real (SCRIP_r8), dimension(ncorners), intent(in) ::
     &     cell_corner_x,       ! Coordinates of cell corners
     &     cell_corner_y        ! Could be x-y or lat-lon or ...

      logical (SCRIP_logical), intent(in) ::
     &     latlon               ! Are coordinates in latlon space?

!----------------------------------------------------------------------
!
!     Output arguments
!
!----------------------------------------------------------------------

      logical (SCRIP_logical), intent(out) ::
     &     inpoly             ! Is point in the polygon?

      logical (SCRIP_logical), intent(out) ::
     &     lboundary          ! Is point on the boundary of the polygon?

      integer (SCRIP_i4), intent(out) ::
     &     edgeid             ! if point is on boundary, which local
                              ! edge is it on? (0 otherwise)

!----------------------------------------------------------------------
!
!     Local variables
!
!----------------------------------------------------------------------

      integer (SCRIP_i4) :: c, n, next_n

      real (SCRIP_r8) :: x1, y1, x2, y2, vec1_x, vec1_y, vec2_x, vec2_y,
     &     vec3_x, vec3_y, vec1_len, vec2_len, vec3_len,
     &     cross_product, dot_product


      !***
      !*** So, with each edge of the polygon form a triangle
      !*** by connecting a 'central' point to the endpoints of 
      !*** the edge. Then check containment of the point in this tri
      !***
      
      inpoly = .false.
      lboundary = .false.
      
      c = 0
      do n = 1, ncorners
         next_n = MOD(n,ncorners) + 1

         x1 = cell_corner_x(n)
         y1 = cell_corner_y(n)
         x2 = cell_corner_x(next_n)
         y2 = cell_corner_y(next_n)

         if (((y1 > pty .and. y2 <= pty) .or.
     &        (y2 > pty .and. y1 <= pty)) .and.
     &        (ptx <= (x1 + (pty-y1)*(x2-x1)/(y2-y1)))) then

            c = 1 - c
            
         endif
      enddo

      if (c .eq. 1) inpoly = .true.
         

      !*** Check if the point is on the boundary of the polygon

      do n = 1, ncorners

         next_n = MOD(n,ncorners) + 1

         x1 = cell_corner_x(n)
         y1 = cell_corner_y(n)
         x2 = cell_corner_x(next_n)
         y2 = cell_corner_y(next_n)

         vec1_x = x2 - x1
         vec1_y = y2 - y1
         vec1_len = sqrt(vec1_x*vec1_x + vec1_y*vec1_y)
         vec1_x = vec1_x/vec1_len
         vec1_y = vec1_y/vec1_len

         vec2_x = ptx - x1
         vec2_y = pty - y1
         vec2_len = sqrt(vec2_x*vec2_x + vec2_y*vec2_y)

         cross_product = vec1_x*vec2_y - vec2_x*vec1_y
         if (abs(cross_product) > tiny .and. vec2_len > tiny) then
            cross_product = cross_product/vec2_len
         endif

         if (abs(cross_product) < 1e5*tiny .and. 
     &        abs(cross_product) > 10*tiny) then

            !*** Sometimes when the point is too close to a vertex
            !*** then the cross product computation has errors due
            !*** to subtraction of two small numbers - So check w.r.t.
            !*** other vertex of the segment as well

            vec3_x = ptx - x2
            vec3_y = pty - y2
            vec3_len = sqrt(vec3_x*vec3_x + vec3_y*vec3_y)

            cross_product = -vec1_x*vec3_y + vec1_y*vec3_x
            if (abs(cross_product) > tiny .and. vec3_len > tiny) then
               !***
               !*** Normalize only if we won't be dividing two small 
               !*** numbers
               cross_product = cross_product/vec3_len
            endif
         endif

         if (abs(cross_product) < 10*tiny) then

            if (vec1_x*vec1_x+vec1_y*vec1_y .le. tiny*tiny) then
               cross_product = one
            else
               dot_product = vec1_x*vec2_x + vec1_y*vec2_y

               if (dot_product >= 0 .and. dot_product <= vec1_len) then
                  inpoly = .true.
                  lboundary = .true.
                  edgeid = n
                  exit
               endif
            endif
         endif

      enddo
         
      return   

!----------------------------------------------------------------------

      end subroutine ptinpolygen2

!**********************************************************************





!----------------------------------------------------------------------
!
!----------------------------------------------------------------------

      subroutine get_srch_cells(cell_add, cell_grid_num, srch_grid_num, 
     &     num_srch_cells, srch_add, srch_corners, 
     &     srch_corner_lat, srch_corner_lon,
     &     srch_center_lat, srch_center_lon)

!----------------------------------------------------------------------
!
!     Input arguments
!
!----------------------------------------------------------------------

      integer (SCRIP_i4), intent(in) ::
     &     cell_add,         ! cell in whose nbrhood we must find other 
     &     cell_grid_num,    ! cells grid number from which 'cell_add' 
     &     srch_grid_num     ! is grid number in which we must find
                             ! search cells 

!----------------------------------------------------------------------
!
!     Output arguments
!
!----------------------------------------------------------------------

      integer (SCRIP_i4), intent(out) ::
     &     num_srch_cells,
     &     srch_corners         ! Number of corners for search cells

      integer (SCRIP_i4), dimension(:), allocatable, intent(out) ::
     &     srch_add             ! Global addresses of search cells

      real (SCRIP_r8), dimension(:,:), allocatable, intent(out) ::
     &     srch_corner_lat, srch_corner_lon

      real (SCRIP_r8), dimension(:), allocatable, intent(out) ::
     &     srch_center_lat, srch_center_lon


!-----------------------------------------------------------------------
!
!     Local arguments
!
!-----------------------------------------------------------------------

      logical (SCRIP_logical), dimension(:), allocatable ::
     &     srch_mask

      integer (SCRIP_i4) :: grid1_add, grid2_add, max_add, min_add,
     &     n

!-----------------------------------------------------------------------

      num_srch_cells = 0
         
      !***
      !*** restrict searches first using search bins
      !***

      if (last_cell_add_get_srch_cells /= cell_add .or. 
     &     last_cell_grid_num_get_srch_cells /= cell_grid_num .or.
     &     last_srch_grid_num_get_srch_cells /= srch_grid_num) then

         if (first_call_get_srch_cells) then
            first_call_get_srch_cells = .false.
            num_srch_cells_loc_get_srch_cells = 0
            srch_corners_loc_get_srch_cells = 0
            last_cell_add_get_srch_cells = 0
            last_cell_grid_num_get_srch_cells = 0
            last_srch_grid_num_get_srch_cells = 0
         else
            if (num_srch_cells_loc_get_srch_cells .gt. 0) then
               deallocate(srch_add_loc_get_srch_cells,
     &              srch_corner_lat_loc_get_srch_cells,
     &              srch_corner_lon_loc_get_srch_cells,
     &              srch_center_lat_loc_get_srch_cells,
     &              srch_center_lon_loc_get_srch_cells)
            endif

         endif


         last_cell_add_get_srch_cells = cell_add
         last_cell_grid_num_get_srch_cells = cell_grid_num
         last_srch_grid_num_get_srch_cells = srch_grid_num


         if (cell_grid_num == 1) then

            if (srch_grid_num == 1) then

               !*** Grid 1 neighbors of grid 1 cell

               allocate(srch_mask(grid1_size))

               min_add = grid1_size
               max_add = 1
               do n=1,num_srch_bins
                  if (cell_add >= bin_addr1(1,n) .and.
     &                 cell_add <= bin_addr1(2,n)) then
                     min_add = min(min_add, bin_addr1(1,n))
                     max_add = max(max_add, bin_addr1(2,n))
                  endif
               end do

               !***
               !*** further restrict searches using bounding boxes
               !***

               num_srch_cells_loc_get_srch_cells = 0
               do grid1_add = min_add,max_add
                  srch_mask(grid1_add) = 
     &                 (grid1_bound_box(1,grid1_add) <= 
     &                 grid1_bound_box(2,cell_add)) .and.
     &                 (grid1_bound_box(2,grid1_add) >= 
     &                 grid1_bound_box(1,cell_add)) .and.
     &                 (grid1_bound_box(3,grid1_add) <= 
     &                 grid1_bound_box(4,cell_add)) .and.
     &                 (grid1_bound_box(4,grid1_add) >= 
     &                 grid1_bound_box(3,cell_add))
                  
                  if (srch_mask(grid1_add)) 
     &                 num_srch_cells_loc_get_srch_cells =
     &                 num_srch_cells_loc_get_srch_cells+1
               end do

               if (num_srch_cells_loc_get_srch_cells /= 0) then

                  !***
                  !*** create search arrays
                  !***

                  allocate(srch_add_loc_get_srch_cells
     &            (num_srch_cells_loc_get_srch_cells),
     &            srch_corner_lat_loc_get_srch_cells
     &            (grid1_corners,num_srch_cells_loc_get_srch_cells),
     &            srch_corner_lon_loc_get_srch_cells
     &            (grid1_corners,num_srch_cells_loc_get_srch_cells),
     &            srch_center_lat_loc_get_srch_cells
     &            (num_srch_cells_loc_get_srch_cells),
     &            srch_center_lon_loc_get_srch_cells
     &            (num_srch_cells_loc_get_srch_cells))

                  n = 0
                  do grid1_add = min_add,max_add
                     if (srch_mask(grid1_add)) then
                        n = n+1
                        srch_add_loc_get_srch_cells(n) = grid1_add
                        srch_corner_lat_loc_get_srch_cells(:,n) = 
     &                       grid1_corner_lat(:,grid1_add)
                        srch_corner_lon_loc_get_srch_cells(:,n) = 
     &                       grid1_corner_lon(:,grid1_add)
                        srch_center_lat_loc_get_srch_cells(n) = 
     &                       grid1_center_lat(grid1_add)
                        srch_center_lon_loc_get_srch_cells(n) = 
     &                       grid1_center_lon(grid1_add)
                     endif
                  end do

                  srch_corners_loc_get_srch_cells = grid1_corners
               endif

               deallocate(srch_mask)

            else
            
               !*** Grid 2 neighbors of grid 1 cell

               allocate(srch_mask(grid2_size))

               min_add = grid2_size
               max_add = 1
               do n=1,num_srch_bins
                  if (cell_add >= bin_addr1(1,n) .and.
     &                 cell_add <= bin_addr1(2,n)) then
                     min_add = min(min_add, bin_addr2(1,n))
                     max_add = max(max_add, bin_addr2(2,n))
                  endif
               end do

               !***
               !*** further restrict searches using bounding boxes
               !***

               num_srch_cells_loc_get_srch_cells = 0
               do grid2_add = min_add,max_add
                  srch_mask(grid2_add) = 
     &                 (grid2_bound_box(1,grid2_add) <= 
     &                 grid1_bound_box(2,cell_add)) .and.
     &                 (grid2_bound_box(2,grid2_add) >= 
     &                 grid1_bound_box(1,cell_add)) .and.
     &                 (grid2_bound_box(3,grid2_add) <= 
     &                 grid1_bound_box(4,cell_add)) .and.
     &                 (grid2_bound_box(4,grid2_add) >= 
     &                 grid1_bound_box(3,cell_add))


                  
                  if (srch_mask(grid2_add)) 
     &                 num_srch_cells_loc_get_srch_cells = 
     &                 num_srch_cells_loc_get_srch_cells+1
               end do


               if (num_srch_cells_loc_get_srch_cells /= 0) then

                  !***
                  !*** create search arrays
                  !***

                  allocate(srch_add_loc_get_srch_cells(
     &            num_srch_cells_loc_get_srch_cells),
     &            srch_corner_lat_loc_get_srch_cells(
     &            grid2_corners,num_srch_cells_loc_get_srch_cells),
     &            srch_corner_lon_loc_get_srch_cells(
     &            grid2_corners,num_srch_cells_loc_get_srch_cells),
     &            srch_center_lat_loc_get_srch_cells(
     &            num_srch_cells_loc_get_srch_cells),
     &            srch_center_lon_loc_get_srch_cells(
     &            num_srch_cells_loc_get_srch_cells))

                  n = 0
                  do grid2_add = min_add,max_add
                     if (srch_mask(grid2_add)) then
                        n = n+1
                        srch_add_loc_get_srch_cells(n) = grid2_add
                        srch_corner_lat_loc_get_srch_cells(:,n) = 
     &                       grid2_corner_lat(:,grid2_add)
                        srch_corner_lon_loc_get_srch_cells(:,n) = 
     &                       grid2_corner_lon(:,grid2_add)
                        srch_center_lat_loc_get_srch_cells(n) = 
     &                       grid2_center_lat(grid2_add)
                        srch_center_lon_loc_get_srch_cells(n) = 
     &                       grid2_center_lon(grid2_add)
                     endif
                  end do
                  
                  srch_corners_loc_get_srch_cells = grid2_corners
               endif

               deallocate(srch_mask)
            endif

         else

            if (srch_grid_num == 1) then

               !*** Grid 1 neighbors of grid 2 cell

               allocate(srch_mask(grid1_size))

               min_add = grid1_size
               max_add = 1
               do n=1,num_srch_bins
                  if (cell_add >= bin_addr2(1,n) .and.
     &                 cell_add <= bin_addr2(2,n)) then
                     min_add = min(min_add, bin_addr1(1,n))
                     max_add = max(max_add, bin_addr1(2,n))
                  endif
               end do

              !***
              !*** further restrict searches using bounding boxes
              !***

               num_srch_cells_loc_get_srch_cells = 0
               do grid1_add = min_add,max_add
                  srch_mask(grid1_add) = 
     &                 (grid1_bound_box(1,grid1_add) <= 
     &                 grid2_bound_box(2,cell_add)) .and.
     &                 (grid1_bound_box(2,grid1_add) >= 
     &                 grid2_bound_box(1,cell_add)) .and.
     &                 (grid1_bound_box(3,grid1_add) <= 
     &                 grid2_bound_box(4,cell_add)) .and.
     &                 (grid1_bound_box(4,grid1_add) >= 
     &                 grid2_bound_box(3,cell_add))

                  if (srch_mask(grid1_add)) 
     &                 num_srch_cells_loc_get_srch_cells = 
     &                 num_srch_cells_loc_get_srch_cells+1
               end do


               if (num_srch_cells_loc_get_srch_cells /= 0) then

                  !***
                  !*** create search arrays
                  !***

                  allocate(srch_add_loc_get_srch_cells(
     &            num_srch_cells_loc_get_srch_cells),
     &            srch_corner_lat_loc_get_srch_cells(
     &            grid1_corners,num_srch_cells_loc_get_srch_cells),
     &            srch_corner_lon_loc_get_srch_cells(
     &            grid1_corners,num_srch_cells_loc_get_srch_cells),
     &            srch_center_lat_loc_get_srch_cells(
     &            num_srch_cells_loc_get_srch_cells),
     &            srch_center_lon_loc_get_srch_cells(
     &            num_srch_cells_loc_get_srch_cells))

                  n = 0
                  do grid1_add = min_add,max_add
                     if (srch_mask(grid1_add)) then
                        n = n+1
                        srch_add_loc_get_srch_cells(n) = grid1_add
                        srch_corner_lat_loc_get_srch_cells(:,n) = 
     &                       grid1_corner_lat(:,grid1_add)
                        srch_corner_lon_loc_get_srch_cells(:,n) = 
     &                       grid1_corner_lon(:,grid1_add)
                        srch_center_lat_loc_get_srch_cells(n) = 
     &                       grid1_center_lat(grid1_add)
                        srch_center_lon_loc_get_srch_cells(n) = 
     &                       grid1_center_lon(grid1_add)
                     endif
                  end do

                  srch_corners_loc_get_srch_cells = grid1_corners
               endif

               deallocate(srch_mask)

            else

               !*** Grid 2 neighbors of grid 2 cell

               allocate(srch_mask(grid2_size))

               min_add = grid2_size
               max_add = 1
               do n=1,num_srch_bins
                  if (cell_add >= bin_addr2(1,n) .and.
     &                 cell_add <= bin_addr2(2,n)) then
                     min_add = min(min_add, bin_addr2(1,n))
                     max_add = max(max_add, bin_addr2(2,n))
                  endif
               end do

               !***
               !*** further restrict searches using bounding boxes
               !***

               num_srch_cells_loc_get_srch_cells = 0
               do grid2_add = min_add,max_add
                  srch_mask(grid2_add) = 
     &                 (grid2_bound_box(1,grid2_add) <= 
     &                 grid2_bound_box(2,cell_add)) .and.
     &                 (grid2_bound_box(2,grid2_add) >= 
     &                 grid2_bound_box(1,cell_add)) .and.
     &                 (grid2_bound_box(3,grid2_add) <= 
     &                 grid2_bound_box(4,cell_add)) .and.
     &                 (grid2_bound_box(4,grid2_add) >= 
     &                 grid2_bound_box(3,cell_add))

                  if (srch_mask(grid2_add)) 
     &                 num_srch_cells_loc_get_srch_cells = 
     &                 num_srch_cells_loc_get_srch_cells+1
               end do
               

               if (num_srch_cells_loc_get_srch_cells /= 0) then

                  !***
                  !*** create search arrays
                  !***

                  allocate(srch_add_loc_get_srch_cells(
     &            num_srch_cells_loc_get_srch_cells),
     &            srch_corner_lat_loc_get_srch_cells(
     &            grid2_corners,num_srch_cells_loc_get_srch_cells),
     &            srch_corner_lon_loc_get_srch_cells(
     &            grid2_corners,num_srch_cells_loc_get_srch_cells),
     &            srch_center_lat_loc_get_srch_cells(
     &            num_srch_cells_loc_get_srch_cells),
     &            srch_center_lon_loc_get_srch_cells(
     &            num_srch_cells_loc_get_srch_cells))

                  n = 0
                  do grid2_add = min_add,max_add
                     if (srch_mask(grid2_add)) then
                        n = n+1
                        srch_add_loc_get_srch_cells(n) = grid2_add
                        srch_corner_lat_loc_get_srch_cells(:,n) = 
     &                       grid2_corner_lat(:,grid2_add)
                        srch_corner_lon_loc_get_srch_cells(:,n) = 
     &                       grid2_corner_lon(:,grid2_add)
                        srch_center_lat_loc_get_srch_cells(n) = 
     &                       grid2_center_lat(grid2_add)
                        srch_center_lon_loc_get_srch_cells(n) = 
     &                       grid2_center_lon(grid2_add)
                     endif
                  end do
                  
                  srch_corners_loc_get_srch_cells = grid2_corners
               endif

               deallocate(srch_mask)

            endif

         endif

      endif


      num_srch_cells = num_srch_cells_loc_get_srch_cells

      if (num_srch_cells .eq. 0) then
         return
      endif

      srch_corners = srch_corners_loc_get_srch_cells
      allocate(srch_add(num_srch_cells),
     &     srch_corner_lat(srch_corners,num_srch_cells),
     &     srch_corner_lon(srch_corners,num_srch_cells),
     &     srch_center_lat(num_srch_cells),
     &     srch_center_lon(num_srch_cells))
      srch_add = srch_add_loc_get_srch_cells
      srch_corner_lat = srch_corner_lat_loc_get_srch_cells
      srch_corner_lon = srch_corner_lon_loc_get_srch_cells
      srch_center_lat = srch_center_lat_loc_get_srch_cells
      srch_center_lon = srch_center_lon_loc_get_srch_cells

      end subroutine get_srch_cells


!**********************************************************************


!----------------------------------------------------------------------
!
!     Find cell adjacent to edge (edge_id) of given cell (cell_add)
!
!----------------------------------------------------------------------

      subroutine find_adj_cell(cell_add, edge_id, cell_grid_num, 
     &     adj_add) 

!----------------------------------------------------------------------
!
!     Input variables
!
!----------------------------------------------------------------------

      integer (SCRIP_i4), intent(in) ::
     &     cell_add,              ! cell whose edge we are checking
     &     edge_id,               ! index of edge that we are check
     &     cell_grid_num          ! grid to which cell belongs

!----------------------------------------------------------------------
!
!     Output variables
!
!----------------------------------------------------------------------

      integer (SCRIP_i4), intent(out) :: adj_add

!----------------------------------------------------------------------
!
!     Local variables
!
!----------------------------------------------------------------------

      integer (SCRIP_i4) :: i, inx, n, global_add
      logical (SCRIP_logical) :: found
      real (SCRIP_r8) :: lat1, lon1, lat2, lon2

      adj_add = 0

      if (cell_grid_num .eq. 1) then

         i = edge_id
         inx = 1 + mod(edge_id,grid1_corners)

         lat1 = grid1_corner_lat(i,cell_add)
         lon1 = grid1_corner_lon(i,cell_add)
         lat2 = grid1_corner_lat(inx,cell_add)
         lon2 = grid1_corner_lon(inx,cell_add)

         !***
         !*** Often the cell with the next or previous index is 
         !*** the adjacent cell. Check that first
         !***

         if (cell_add .lt. grid1_size) then

            global_add = cell_add + 1

            do i = 1, grid1_corners
               inx = mod(i,grid1_corners)+1
               if (abs(grid1_corner_lat(inx,global_add)-lat1) .le. tiny 
     &              .and.
     &              abs(grid1_corner_lat(i,global_add)-lat2) .le. tiny 
     &              .and.
     &              abs(grid1_corner_lon(inx,global_add)-lon1) .le. tiny
     &              .and.
     &              abs(grid1_corner_lon(i,global_add)-lon2) .le. tiny) 
     &              then
                  
                  adj_add = global_add
                  return
               endif
            enddo

         endif
         
         if (cell_add .gt. 1) then
            
            global_add = cell_add - 1

            do i = 1, grid1_corners
               inx = mod(i,grid1_corners)+1
               if (abs(grid1_corner_lat(inx,global_add)-lat1) .le. tiny 
     &              .and.
     &              abs(grid1_corner_lat(i,global_add)-lat2) .le. tiny 
     &              .and.
     &              abs(grid1_corner_lon(inx,global_add)-lon1) .le. tiny
     &              .and.
     &              abs(grid1_corner_lon(i,global_add)-lon2) .le. tiny) 
     &              then
                  
                  adj_add = global_add
                  return
               endif
            enddo

         endif
         
         

      else

         i = edge_id
         inx = 1 + mod(edge_id,grid2_corners)

         lat1 = grid2_corner_lat(i,cell_add)
         lon1 = grid2_corner_lon(i,cell_add)
         lat2 = grid2_corner_lat(inx,cell_add)
         lon2 = grid2_corner_lon(inx,cell_add)


         !***
         !*** Often the cell with the next or previous index is 
         !*** the adjacent cell. Check that first
         !***

         if (cell_add .lt. grid2_size) then

            global_add = cell_add + 1

            do i = 1, grid2_corners
               inx = mod(i,grid2_corners)+1
               if (abs(grid2_corner_lat(inx,global_add)-lat1) .le. tiny 
     &              .and.
     &              abs(grid2_corner_lat(i,global_add)-lat2) .le. tiny 
     &              .and.
     &              abs(grid2_corner_lon(inx,global_add)-lon1) .le. tiny
     &              .and.
     &              abs(grid2_corner_lon(i,global_add)-lon2) .le. tiny) 
     &              then
                  
                  adj_add = global_add
                  return
               endif
            enddo

         endif
         
         if (cell_add .gt. 1) then
            
            global_add = cell_add - 1

            do i = 1, grid2_corners
               inx = mod(i,grid2_corners)+1
               if (abs(grid2_corner_lat(inx,global_add)-lat1) .le. tiny 
     &              .and.
     &              abs(grid2_corner_lat(i,global_add)-lat2) .le. tiny 
     &              .and.
     &              abs(grid2_corner_lon(inx,global_add)-lon1) .le. tiny
     &              .and.
     &              abs(grid2_corner_lon(i,global_add)-lon2) .le. tiny) 
     &              then
                  
                  adj_add = global_add
                  return
               endif
            enddo

         endif
         

      endif


      if (cell_add /= last_cell_find_adj_cell .or. 
     &     cell_grid_num /= last_cell_grid_num_find_adj_cell) then

         last_cell_find_adj_cell = cell_add
         last_cell_grid_num_find_adj_cell = cell_grid_num

         if (first_call_find_adj_cell) then
            first_call_find_adj_cell = .false.
            last_cell_find_adj_cell = 0
            last_cell_grid_num_find_adj_cell = 0
         else
            if (num_srch_cells_find_adj_cell .gt. 0) then
               deallocate(srch_add_find_adj_cell, 
     &         srch_corner_lat_find_adj_cell, 
     &         srch_corner_lon_find_adj_cell,
     &         srch_center_lat_find_adj_cell, 
     &         srch_center_lon_find_adj_cell)
            endif
         endif

         call get_srch_cells(cell_add, cell_grid_num, cell_grid_num,
     &        num_srch_cells_find_adj_cell, srch_add_find_adj_cell, 
     &        srch_corners_find_adj_cell, srch_corner_lat_find_adj_cell,
     &        srch_corner_lon_find_adj_cell, 
     &        srch_center_lat_find_adj_cell, 
     &        srch_center_lon_find_adj_cell)

      endif


      found = .false.
      do n = 1, num_srch_cells_find_adj_cell
         
         global_add = srch_add_find_adj_cell(n)

         do i = 1, srch_corners_find_adj_cell
            inx = mod(i,srch_corners_find_adj_cell)+1
            if (abs(srch_corner_lat_find_adj_cell(inx,n)-lat1) .le. tiny
     &           .and.
     &           abs(srch_corner_lat_find_adj_cell(i,n)-lat2) .le. tiny 
     &           .and.
     &           abs(srch_corner_lon_find_adj_cell(inx,n)-lon1) .le.tiny
     &           .and.
     &           abs(srch_corner_lon_find_adj_cell(i,n)-lon2) .le. tiny)
     &           then

               adj_add = global_add
               found = .true.

               exit
            endif
         enddo

         if (found) exit
         
      enddo

      return
      end subroutine find_adj_cell


!----------------------------------------------------------------------
!
!     Given points inside and outside a cell, converge to the boundary
!
!----------------------------------------------------------------------


      subroutine converge_to_bdry(cell_add, cell_grid_num,
     &     ncorners, cell_corner_lat,
     &     cell_corner_lon, cell_center_lat, cell_center_lon, 
     &     inpt_lat, inpt_lon, outpt_lat, outpt_lon,
     &     bpt_lat, bpt_lon, bedgeid)

!----------------------------------------------------------------------
!
!     Input variables
!
!----------------------------------------------------------------------

      integer (SCRIP_i4), intent(in) ::
     &     cell_add,            ! Cell in which we are operating
     &     cell_grid_num,       ! Grid to which cell belongs
     &     ncorners             ! Number of corners in cell

      real (SCRIP_r8), dimension(ncorners), intent(in) ::
     &     cell_corner_lat,     ! Latitude values of cell corners
     &     cell_corner_lon      ! Longitude values of cell corners

      real (SCRIP_r8), intent(in) ::
     &     cell_center_lat,     ! Latitude of cell center
     &     cell_center_lon,     ! Longitude of cell center,
     &     inpt_lat,            ! Latitude of inside point
     &     inpt_lon,            ! Longitude of inside point
     &     outpt_lat,           ! Latitude of outside point
     &     outpt_lon            ! Longitude of outside point


!----------------------------------------------------------------------
!
!     Output variables
!
!----------------------------------------------------------------------

      real (SCRIP_r8), intent(out) ::
     &     bpt_lat,             ! Latitude of boundary point
     &     bpt_lon              ! Longitude of boundary point

      integer (SCRIP_i4), intent(out) ::
     &     bedgeid              ! ID of edge that point converged to

!----------------------------------------------------------------------
!
!     Local variables
!
!----------------------------------------------------------------------

      logical (SCRIP_logical) ::
     &     converged, 
     &     lboundary,
     &     inpoly

      integer (SCRIP_i4) ::
     &     it

      real (SCRIP_r8) ::
     &     lat1, lon1,
     &     lat2, lon2,
     &     midlat, midlon

      bedgeid = 0

      lat1 = inpt_lat
      lon1 = inpt_lon
      lat2 = outpt_lat
      lon2 = outpt_lon      


      converged = .false.
      it = 0
      do while (.not. converged) 
         
         midlat = (lat1+lat2)/2.0
         if (abs(lon1-lon2) < pi) then            
            midlon = (lon1+lon2)/2.0
         else
            midlon = (lon1+lon2)/2.0 - pi2
         endif
                               
                               
         call ptincell(midlat, midlon, 
     &        cell_add, ncorners, 
     &        cell_corner_lat, cell_corner_lon,
     &        cell_center_lat, cell_center_lon,
     &        cell_grid_num,
     &        inpoly, lboundary, bedgeid) 
                               
         if (inpoly) then
            lat1 = midlat
            lon1 = midlon
         else
            lat2 = midlat
            lon2 = midlon
         endif
                               
         if (abs(lat1-lat2) < tiny .and.
     &        abs(lon1-lon2) < tiny .and. lboundary) then
            converged = .true.
         endif                 
                               
         if (it > 100) then
            exit
         endif
                               
         it = it + 1
      enddo                     ! do while (not converged)


      bpt_lat = midlat
      bpt_lon = midlon

      end subroutine converge_to_bdry




      end module scrip_remap_conservative

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
