!/ ------------------------------------------------------------------- /
      module scrip_interface

!  1. Original author :
!
!     Erick Rogers, NRL
!
!  2. Last update :
!
!     See revisions.
!
!  3. Revisions :
!
!     29-Apr-2011 : Origination                         ( version 4.01 )
!
!  4. Copyright :
!
!  5. Purpose :
!
!     Routines to provide interface between WMGHGH and SCRIP
!
!  6. Variables and types :
!
!  7. Subroutines and functions :
!
!  8. Subroutines and functions used :
!
!  9. Remarks :
!     On parallelization:
!         When WW3 is working on WMGHGH, each MPI task performs identical computations
!         and therefore creates its own copy of the an identical solution. Thus there
!         is no MPI parallelization of WMGHGH. Presumably, this was designed with the
!         expectation that WMGHGH computations are fairly quick. Unfortunately, SCRIP 
!         can be slow for large grids, making WMGHGH slow. It is possible to run SCRIP 
!         with OpenMP (with pgf90, it simply requires that we compile with the -mp flag). 
!         Unfortunately, this doesn't change the situation whereby each MPI task calls 
!         SCRIP and creates its own copy of an identical solution. Thus, if we want 
!         to use OpenMP to speed up SCRIP, it means that we have to leave some cores idle 
!         during the actual WW3 compuation. Example: using pgf90, we execute ww3_multi 
!         with mpirun -n 4. Thus, four MPI tasks are created. And say we compiled SCRIP 
!         using -mp and have nthreads=2 specified in scrip_interface.f. When multi gets 
!         to SCRIP, each MPI task will create 2 OpenMP threads, which means that we'd 
!         better have 8 threads available to do this, and also note that 4 threads will 
!         be idle after SCRIP is finished. This isn't an ideal solution. A better 
!         prior to the SCRIP call, is instructed to only compute using one MPI thread. 
!         Then, SCRIP does its thing using OpenMPI and returns to WMGHGH. Then, WMGHGH
!         shares all the information with all the nodes and continues.
!      On connection to WW3: 
!         This is handled by "scrip_wrapper" routines. These routines formerly resided 
!         in this file (scrip_interface.f*) but are now in wmscrpmd.ftn
!
! 10. Switches :
!
! 11. Source code :

      use SCRIP_KindsMod    ! defines data types

!.....notes: since the calling subroutine (wmghgh) does not know a priori the size 
!............of the weights arrays, we manage the communication via this module.
      implicit none

      type weight_data
         integer (SCRIP_i4)              :: n    ! number of weights for dst cell
                                                 ! n is equivalent to NR1 and NLOC in original WMGHGH
                                                 ! NR1 is the counter of |MAPSTA|=1 (indicates sea point)
         integer (SCRIP_i4)              :: NR0  ! counter of MAPSTA=0 (indicates excluded point)
         integer (SCRIP_i4)              :: NR2  ! counter of |MAPSTA|=2 (indicates boundary point)
         integer (SCRIP_i4)              :: NRL  ! counter of MAPSTA=0 (indicates excluded point) and MAPST2=0 (indicates land)
         real    (SCRIP_r8), allocatable :: w(:) ! weights, sized by n, formerly wxwy(:,:)
         integer (SCRIP_i4), allocatable :: k(:) ! source grid cells, sized by n, formerly ksrc(:,:)
      end type weight_data

      type(weight_data), allocatable :: wgtdata(:)

      contains
!/ ------------------------------------------------------------------- /



!#######################################################################
      subroutine scrip_clear
!#######################################################################

!  1. Original author :
!
!     Erick Rogers, NRL
!
!  2. Last update :
!
!     See revisions.
!
!  3. Revisions :
!
!     5-May-2011 : Origination                         ( version 4.01 )
!
!  4. Copyright :
!
!  5. Purpose : 
!
!     "Clear" all variables declared at module level of SCRIP routines
!     (clear "common block" equivalent)
!
!  6. Method :
!
!     rules: 
!        - if not an array (scalar), set to zero or other start value
!        - if dimensioned array, set to zero
!        - if allocatable array, deallocate
!        - private variables: ignore, since we would need to 
!             make them public in order to clear them, which may do more
!             harm than good.
!
!  7. Parameters, Variables and types :
!  
!  8. Called by : 
!
!     Subroutine SCRIP_interface
!
!  9. Subroutines and functions used :
!
!     None
!
! 10. Error messages: 
!
! 11. Remarks :
!
!     We "clear" all variables with "save" attribute, 
!     both "module variables" and "subroutine variables"
!     including all variables that are initialized with a value
!     in the type declaration, e.g. "real :: x=5.0"
!
! 12. Structure :
!
! 13. Switches :
!
! 14. Source code :

      use SCRIP_KindsMod
      use scrip_timers
      use scrip_remap_vars
      use scrip_remap_conservative
      use scrip_iounitsmod
      use scrip_grids

      implicit none

      call timers_init ! takes care of all variables in timers.f

!.....scrip_remap_vars.f :
      max_links_map1=0
      num_links_map1=0
      max_links_map2=0
      num_links_map2=0
      num_maps=0 
      num_wts=0   
      map_type=0  
      norm_opt=0
      resize_increment=0
      if(allocated(grid1_add_map1))deallocate(grid1_add_map1)
      if(allocated(grid2_add_map1))deallocate(grid2_add_map1)
      if(allocated(grid1_add_map2))deallocate(grid1_add_map2)
      if(allocated(grid2_add_map2))deallocate(grid2_add_map2)
      if(allocated(wts_map1))deallocate(wts_map1)
      if(allocated(wts_map2))deallocate(wts_map2)

!.....remap_conserv.f :
!.....scalars: 
      first_call_store_link_cnsrv = .true.
      first_call_locate_segstart= .true.
      first_call_locate_point= .true.
      first_call_get_srch_cells=.true.
      first_call_find_adj_cell=.true.
      avoid_pole_count = 0
      avoid_pole_offset = tiny
      last_cell_locate_segstart=0
      last_cell_grid_num_locate_segstart=0
      last_srch_grid_num_locate_segstart=0
      num_srch_cells_locate_segstart=0
      last_cell_locate_point=0
      last_cell_grid_num_locate_point=0
      last_srch_grid_num_locate_point=0
      num_srch_cell_locate_points=0
      srch_corners_locate_point=0
      srch_corners_find_adj_cell=0
      srch_corners_locate_segstart=0
      srch_corners_loc_get_srch_cells=0
      num_srch_cells_loc_get_srch_cells=0
      num_srch_cells_find_adj_cell=0
      last_cell_add_get_srch_cells=0
      last_cell_grid_num_get_srch_cells=0
      last_srch_grid_num_get_srch_cells=0
      last_cell_find_adj_cell=0
      last_cell_grid_num_find_adj_cell=0
!.....arrays :
      if(allocated(link_add1))deallocate(link_add1)
      if(allocated(link_add2))deallocate(link_add2)

      if(allocated(srch_add_loc_get_srch_cells))deallocate(srch_add_loc_get_srch_cells)
      if(allocated(srch_corner_lat_loc_get_srch_cells))deallocate(srch_corner_lat_loc_get_srch_cells)
      if(allocated(srch_corner_lon_loc_get_srch_cells))deallocate(srch_corner_lon_loc_get_srch_cells)
      if(allocated(srch_center_lat_loc_get_srch_cells))deallocate(srch_center_lat_loc_get_srch_cells)
      if(allocated(srch_center_lon_loc_get_srch_cells))deallocate(srch_center_lon_loc_get_srch_cells)

      if(allocated(srch_add_find_adj_cell))deallocate(srch_add_find_adj_cell)
      if(allocated(srch_corner_lat_find_adj_cell))deallocate(srch_corner_lat_find_adj_cell)
      if(allocated(srch_corner_lon_find_adj_cell))deallocate(srch_corner_lon_find_adj_cell)
      if(allocated(srch_center_lat_find_adj_cell))deallocate(srch_center_lat_find_adj_cell)
      if(allocated(srch_center_lon_find_adj_cell))deallocate(srch_center_lon_find_adj_cell)

      if(allocated(srch_add_locate_segstart))deallocate(srch_add_locate_segstart)
      if(allocated(srch_corner_lat_locate_segstart))deallocate(srch_corner_lat_locate_segstart)
      if(allocated(srch_corner_lon_locate_segstart))deallocate(srch_corner_lon_locate_segstart)
      if(allocated(srch_center_lat_locate_segstart))deallocate(srch_center_lat_locate_segstart)
      if(allocated(srch_center_lon_locate_segstart))deallocate(srch_center_lon_locate_segstart)

      if(allocated(srch_add_locate_point))deallocate(srch_add_locate_point)
      if(allocated(srch_corner_lat_locate_point))deallocate(srch_corner_lat_locate_point)
      if(allocated(srch_corner_lon_locate_point))deallocate(srch_corner_lon_locate_point)
      if(allocated(srch_center_lat_locate_point))deallocate(srch_center_lat_locate_point)
      if(allocated(srch_center_lon_locate_point))deallocate(srch_center_lon_locate_point)

!.....scrip_grids.f :
      grid1_size=0
      grid2_size=0
      grid1_rank=0
      grid2_rank=0
      grid1_corners=0
      grid2_corners=0
      grid1_name=''
      grid2_name=''
      grid1_units=''
      grid2_units=''
      luse_grid_centers=.false.
      luse_grid1_area=.false.
      luse_grid2_area=.false.
      restrict_type=''
      num_srch_bins=0
      if(allocated(bin_addr1))deallocate(bin_addr1)
      if(allocated(bin_addr2))deallocate(bin_addr2)
      if(allocated(bin_lats))deallocate(bin_lats)
      if(allocated(bin_lons))deallocate(bin_lons)
      if(allocated(grid1_dims))deallocate(grid1_dims)
      if(allocated(grid2_dims))deallocate(grid2_dims)
      if(allocated(grid1_mask))deallocate(grid1_mask)
      if(allocated(grid2_mask))deallocate(grid2_mask)
      if(allocated(grid1_center_lat))deallocate(grid1_center_lat)
      if(allocated(grid1_center_lon))deallocate(grid1_center_lon)
      if(allocated(grid2_center_lat))deallocate(grid2_center_lat)
      if(allocated(grid2_center_lon))deallocate(grid2_center_lon)
      if(allocated(grid1_area))deallocate(grid1_area)
      if(allocated(grid2_area))deallocate(grid2_area)
      if(allocated(grid1_area_in))deallocate(grid1_area_in)
      if(allocated(grid2_area_in))deallocate(grid2_area_in)
      if(allocated(grid1_frac))deallocate(grid1_frac)
      if(allocated(grid2_frac))deallocate(grid2_frac)
      if(allocated(grid1_corner_lat))deallocate(grid1_corner_lat)
      if(allocated(grid1_corner_lon))deallocate(grid1_corner_lon)
      if(allocated(grid2_corner_lat))deallocate(grid2_corner_lat)
      if(allocated(grid2_corner_lon))deallocate(grid2_corner_lon)
      if(allocated(grid1_bound_box))deallocate(grid1_bound_box)
      if(allocated(grid2_bound_box))deallocate(grid2_bound_box)
      if(allocated(special_polar_cell1))deallocate(special_polar_cell1)
      if(allocated(special_polar_cell2))deallocate(special_polar_cell2)
      if(allocated(grid1_centroid_lat))deallocate(grid1_centroid_lat)
      if(allocated(grid1_centroid_lon))deallocate(grid1_centroid_lon)
      if(allocated(grid2_centroid_lat))deallocate(grid2_centroid_lat)
      if(allocated(grid2_centroid_lon))deallocate(grid2_centroid_lon)

!#######################################################################
      end subroutine scrip_clear
!#######################################################################

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     subroutine scrip:
!     This routine is the driver for computing the addresses and weights 
!     for interpolating between two grids on a sphere.
!
!-----------------------------------------------------------------------
!
!     CVS:$Id: scrip.f,v 1.6 2001/08/21 21:06:44 pwjones Exp $
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
!     within WW3. Here is a list of modifications:
!     - changed from standalone program to a subroutine, to be called
!          by scrip_wrapper
!     - all modules are now prepended with "scrip_" to make them easier
!          to distinguish from WW3 code.
!     - code changed to free format style (e.g. continuation characters)
!     - print statements added
!     - initial values for settings are changed
!     - initial value for "map_method" added (is set to "conservative")
!     - read of settings from "scrip_in" is removed, and file "scrip_in"
!          is thus no longer used
!     - in context of scrip_wrapper, these initial values for the settings
!          are never changed. Thus, these settings are basically hardwired
!          in here.
!     - lines associated with remap_distance_weight, remap_bilinear, 
!          remap_bicubic are removed. Motivation: we do not need these
!          routines, so we opt to exclude them in our compile
!          (fewer .f files, fewer compiles).
!
!***********************************************************************

      subroutine scrip(src_num, dst_num, l_master, l_read, l_test)

!-----------------------------------------------------------------------

      use SCRIP_KindsMod                  ! module defining data types
      use scrip_constants                  ! module for common constants
      use scrip_iounitsmod                    ! I/O unit manager
      use scrip_timers                     ! CPU timers
      use scrip_grids                      ! module with grid information
      use scrip_remap_vars                 ! common remapping variables
      use scrip_remap_conservative         ! routines for conservative remap

!/SCRIPNC      use scrip_remap_write                ! routines for remap output
!/SCRIPNC      use scrip_remap_read                 ! routines for remap input

      use scrip_errormod

      implicit none

!-----------------------------------------------------------------------
!
!     input variables formerly part of namelist
!
!-----------------------------------------------------------------------

      character (SCRIP_charLength) :: &
                 interp_file1,& ! filename for output remap data (map1)
                 interp_file2,& ! filename for output remap data (map2)
                 map1_name,   & ! name for mapping from grid1 to grid2
                 map2_name,   & ! name for mapping from grid2 to grid1
                 map_method,  & ! choice for mapping method
                 normalize_opt,&! option for normalizing weights
                 output_opt    ! option for output conventions

      integer (SCRIP_i4) :: &
                 nmap          ! number of mappings to compute (1 or 2)

!-----------------------------------------------------------------------
!
!     input variables not part of namelist
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4), intent(in)      :: dst_num   ! number of destination grid GDST
      integer (SCRIP_i4), intent(in)      :: src_num   ! number of source grid GSRC
      logical (SCRIP_Logical), intent(in) :: l_master  ! Am I the master processor (do I/O)?
      logical( SCRIP_Logical), intent(in) :: l_read    ! Do I read the remap file?
      logical(SCRIP_Logical), intent(in) :: l_test     ! Whether to include test output
                                                       ! in subroutines

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4) :: n, &   ! dummy counter
                                 iunit  ! unit number for namelist file

      integer (SCRIP_i4) :: &
           errorCode      ! error flag

      character (12), parameter :: &
           rtnName = 'SCRIP_driver'

!/SCRIPNC      character (LEN=3) :: cdst  ! 3 character number of destination map
!/SCRIPNC      character (LEN=3) :: csrc  ! 3 character number of source map
!/T38      CHARACTER (LEN=10) :: CDATE_TIME(3)
!/T38      INTEGER            :: DATE_TIME(8)
!/T38      INTEGER            :: ELAPSED_TIME, BEG_TIME, END_TIME

!-----------------------------------------------------------------------
!
!     initialize timers and errorcode
!
!-----------------------------------------------------------------------

!/T38 if(l_master)write(SCRIP_stdout,*)'subroutine scrip'

      call timers_init
      do n=1,max_timers
        call timer_clear(n)
      end do

      errorCode = SCRIP_Success

!-----------------------------------------------------------------------
!
!     set variables that were previously read in as a namelist
!
!-----------------------------------------------------------------------

      num_maps = 1
!/SCRIPNC!     Note: Only master does I/O, but all processors need to know about
!/SCRIPNC!     file existence
!/SCRIPNC      interp_file1 = "rmp_src_to_dst_conserv_XXX_XXX.nc"
!/SCRIPNC      interp_file2 = 'not_used.nc'
!/SCRIPNC      map1_name = 'source to destination Conservative Mapping'
!/SCRIPNC      map2_name = 'map not used'
!/SCRIPNC      write(cdst, "(i3.3)") dst_num
!/SCRIPNC      write(csrc, "(i3.3)") src_num
!/SCRIPNC      interp_file1(24:26) = csrc
!/SCRIPNC      interp_file1(28:30) = cdst

      map_method = 'conservative'
      normalize_opt = 'fracarea'
      output_opt = 'scrip'
      restrict_type = 'latitude'
      num_srch_bins = 90 
      luse_grid1_area = .false.
      luse_grid2_area = .false.
      npseg=11 ! or num_polar_segs
      north_thresh=1.5_SCRIP_r8 ! or npole_threshold
      south_thresh=-1.5_SCRIP_r8 ! or spole_threshold
      nthreads=2 ! or num_threads
     
      select case(map_method)
      case ('conservative')
         map_type = map_type_conserv
         luse_grid_centers = .false.
      case ('bilinear')
         map_type = map_type_bilinear
         luse_grid_centers = .true.
      case ('bicubic')
         map_type = map_type_bicubic
         luse_grid_centers = .true.
      case ('distwgt')
         map_type = map_type_distwgt
         luse_grid_centers = .true.
      case ('particle')
         map_type = map_type_particle
         luse_grid_centers = .false.
      case default
         call SCRIP_ErrorSet(errorCode, rtnName, 'unknown mapping method')
         call SCRIP_driverExit(errorCode, 'unknown mapping method')
      end select

      select case(normalize_opt(1:4))
      case ('none')
         norm_opt = norm_opt_none
      case ('frac')
         norm_opt = norm_opt_frcarea
      case ('dest')
         norm_opt = norm_opt_dstarea
      case default
         call SCRIP_ErrorSet(errorCode, rtnName, 'unknown normalization option')
         call SCRIP_driverExit(errorCode, 'unknown normalization option')
      end select

!-----------------------------------------------------------------------
!
!     initialize grid information for both grids
!
!-----------------------------------------------------------------------

!/T38 if(l_master)write(SCRIP_stdout,*)'calling grid_init'

      call grid_init( errorCode,l_master,l_test)

!/T38 if(l_master)write(SCRIP_stdout, *) 'Computing remappings between: ',grid1_name
!/T38 if(l_master)write(SCRIP_stdout, *) '                         and  ',grid2_name

!-----------------------------------------------------------------------
!
!     initialize some remapping variables.
!
!-----------------------------------------------------------------------

      call init_remap_vars

!-----------------------------------------------------------------------
!
!     call appropriate interpolation setup routine based on type of
!     remapping requested. or read in remapping data.
!/SCRIPNC!or, read in remapping data.
!
!-----------------------------------------------------------------------

!/T38      call date_and_time (CDATE_TIME(1), CDATE_TIME(2), CDATE_TIME(3), DATE_TIME)
!/T38      beg_time = ((date_time(5)*60 + date_time(6))*60 +date_time(7))*1000 + date_time(8)
      
!/SCRIPNC      if (l_read) then
!/SCRIPNC         if(l_master)write(SCRIP_stdout, *) 'Reading remapping data from ', interp_file1
!/SCRIPNC         call read_remap_ww3(map1_name, interp_file1, errorCode)
!/T38         call date_and_time (CDATE_TIME(1), CDATE_TIME(2), CDATE_TIME(3), DATE_TIME)
!/T38         end_time = ((date_time(5)*60 + date_time(6))*60 +date_time(7))*1000 + date_time(8)
!/T38         elapsed_time = end_time - beg_time
!/T38         write(0,*) "SCRIP: READING ", elapsed_time, " MSEC"
!/SCRIPNC      else
         select case(map_type)
         case(map_type_conserv)
!/T38       if(l_master)write(SCRIP_stdout,*)'calling remap_conserv'
            call remap_conserv(l_master,l_test)
!/T38       if(l_master)write(SCRIP_stdout,*)'back from remap_conserv'
!/T38            call date_and_time (CDATE_TIME(1), CDATE_TIME(2), CDATE_TIME(3), DATE_TIME)
!/T38            end_time = ((date_time(5)*60 + date_time(6))*60 +date_time(7))*1000 + date_time(8)
!/T38            elapsed_time = end_time - beg_time
!/T38            write(0,*) "SCRIP: CALCULATING ", elapsed_time, " MSEC"
         case default
            call SCRIP_ErrorSet(errorCode, rtnName, 'Invalid Map Type')
            call SCRIP_driverExit(errorCode, 'Invalid Map Type')
         end select

!-----------------------------------------------------------------------
!
!     reduce size of remapping arrays
!
!-----------------------------------------------------------------------
      
!/T38         call date_and_time (CDATE_TIME(1), CDATE_TIME(2), CDATE_TIME(3), DATE_TIME)
!/T38         beg_time = ((date_time(5)*60 + date_time(6))*60 +date_time(7))*1000 + date_time(8)

         if (num_links_map1 /= max_links_map1) then
            call resize_remap_vars(1, num_links_map1-max_links_map1)
         endif
         if ((num_maps > 1) .and. (num_links_map2 /= max_links_map2)) then
            call resize_remap_vars(2, num_links_map2-max_links_map2)
         endif

         call sort_add_v2(grid2_add_map1, grid1_add_map1, wts_map1)

!-----------------------------------------------------------------------
!
!/SCRIPNC!     write remapping info to a file.
!
!-----------------------------------------------------------------------

!/SCRIPNC         if (l_master) then
!/SCRIPNC            write(SCRIP_stdout, *) 'Writing remapping data to ', interp_file1
!/SCRIPNC         endif

!/T38        call date_and_time (CDATE_TIME(1), CDATE_TIME(2), CDATE_TIME(3), DATE_TIME)
!/T38        end_time = ((date_time(5)*60 + date_time(6))*60 +date_time(7))*1000 + date_time(8)
!/T38        elapsed_time = end_time - beg_time
!/T38        write(0,*) "SCRIP: RESIZING ", elapsed_time, " MSEC"
!/T38        call date_and_time (CDATE_TIME(1), CDATE_TIME(2), CDATE_TIME(3), DATE_TIME)
!/T38        beg_time = ((date_time(5)*60 + date_time(6))*60 +date_time(7))*1000 + date_time(8)

! Use  write_remap if you want the extra variables in the .nc files for diagnostics
! Use  write_remap_ww3 if you don't want any extra variables in the .nc files

!/SCRIPNC         if(l_test)then
!/SCRIPNC            call write_remap(map1_name, map2_name, interp_file1, interp_file2, &
!/SCRIPNC            output_opt, l_master, errorCode)
!/SCRIPNC         else
!/SCRIPNC            call write_remap_ww3(map1_name, interp_file1, output_opt, &
!/SCRIPNC            l_master, errorCode)
!/SCRIPNC         endif

!/SCRIPNC      end if

!-----------------------------------------------------------------------

    end subroutine scrip

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine SCRIP_driverExit(errorCode,errormsg)

! !DESCRIPTION:
!  This routine exits the SCRIP driver program. It first calls the 
!  SCRIP error print function to print any errors encountered and then
!  exits the message environment before stopping.
!
! !USES:

   use SCRIP_KindsMod

! !INPUT PARAMETERS:

   integer (SCRIP_i4), intent(in) :: &
      errorCode        ! error flag to detect any errors encountered

   CHARACTER*(*), INTENT(IN)  :: errormsg
!-----------------------------------------------------------------------
!
!  call SCRIP error print function to output any logged errors that
!  were encountered during execution.  Then stop.
!
!-----------------------------------------------------------------------

   write(*,*)'error encountered : ',errorcode
   write(*,*)errormsg

   stop

!-----------------------------------------------------------------------

   end subroutine SCRIP_driverExit

!#######################################################################
   subroutine sort_add_v2(add1, add2, weights)
!#######################################################################

!-----------------------------------------------------------------------
!
!     this routine sorts address and weight arrays based on the
!     destination address with the source address as a secondary
!     sorting criterion.  the method is a standard heap sort.
!
!     sort_add_v2 is identical to subroutine sort_add, but is moved into
!     scrip_interface.ftn and converted to f90 format (line continuations)
!
!-----------------------------------------------------------------------
     
     implicit none

!-----------------------------------------------------------------------
!
!     Input and Output arrays
!
!-----------------------------------------------------------------------
     
     integer (SCRIP_i4), intent(inout), dimension(:) :: &
          add1,  &    ! destination address array (num_links)
          add2        ! source      address array

     real (SCRIP_r8), intent(inout), dimension(:,:) :: &
          weights     ! remapping weights (num_wts, num_links)

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------
     
     integer (SCRIP_i4) ::   &
          num_links,         & ! num of links for this mapping
          num_wts,           & ! num of weights for this mapping
          add1_tmp, add2_tmp,& ! temp for addresses during swap
          lvl, final_lvl,    & ! level indexes for heap sort levels
          chk_lvl1, chk_lvl2, max_lvl

     real (SCRIP_r8), dimension(SIZE(weights,DIM=1)) :: &
          wgttmp              ! temp for holding wts during swap

!-----------------------------------------------------------------------
!
!     determine total number of links to sort and number of weights
!
!-----------------------------------------------------------------------
     
     num_links = SIZE(add1)
     num_wts   = SIZE(weights, DIM=1)

!-----------------------------------------------------------------------
!
!     start at the lowest level (N/2) of the tree and sift lower 
!     values to the bottom of the tree, promoting the larger numbers
!
!-----------------------------------------------------------------------
     
     do lvl=num_links/2,1,-1

        final_lvl = lvl
        add1_tmp = add1(lvl)
        add2_tmp = add2(lvl)
        wgttmp(:) = weights(:,lvl)

        !***
        !*** loop until proper level is found for this link, or reach
        !*** bottom
        !***

        sift_loop1: do

           !***
           !*** find the largest of the two daughters
           !***

           chk_lvl1 = 2*final_lvl
           chk_lvl2 = 2*final_lvl+1
           if (chk_lvl1 .EQ. num_links) chk_lvl2 = chk_lvl1

           if ((add1(chk_lvl1) >  add1(chk_lvl2)) .OR.    &
                ((add1(chk_lvl1) == add1(chk_lvl2)) .AND. &
                (add2(chk_lvl1) >  add2(chk_lvl2)))) then
              max_lvl = chk_lvl1
           else 
              max_lvl = chk_lvl2
           endif

           !***
           !*** if the parent is greater than both daughters,
           !*** the correct level has been found
           !***

           if ((add1_tmp .GT. add1(max_lvl)) .OR.    &
                ((add1_tmp .EQ. add1(max_lvl)) .AND. &
                (add2_tmp .GT. add2(max_lvl)))) then
              add1(final_lvl) = add1_tmp
              add2(final_lvl) = add2_tmp
              weights(:,final_lvl) = wgttmp(:)
              exit sift_loop1

              !***
              !*** otherwise, promote the largest daughter and push
              !*** down one level in the tree.  if haven't reached
              !*** the end of the tree, repeat the process.  otherwise
              !*** store last values and exit the loop
              !***

           else 
              add1(final_lvl) = add1(max_lvl)
              add2(final_lvl) = add2(max_lvl)
              weights(:,final_lvl) = weights(:,max_lvl)

              final_lvl = max_lvl
              if (2*final_lvl > num_links) then
                 add1(final_lvl) = add1_tmp
                 add2(final_lvl) = add2_tmp
                 weights(:,final_lvl) = wgttmp(:)
                 exit sift_loop1
              endif
           endif
        end do sift_loop1
     end do

!-----------------------------------------------------------------------
!
!     now that the heap has been sorted, strip off the top (largest)
!     value and promote the values below
!
!-----------------------------------------------------------------------
     
     do lvl=num_links,3,-1

        !***
        !*** move the top value and insert it into the correct place
        !***

        add1_tmp = add1(lvl)
        add1(lvl) = add1(1)

        add2_tmp = add2(lvl)
        add2(lvl) = add2(1)

        wgttmp(:) = weights(:,lvl)
        weights(:,lvl) = weights(:,1)

        !***
        !*** as above this loop sifts the tmp values down until proper 
        !*** level is reached
        !***

        final_lvl = 1

        sift_loop2: do

           !***
           !*** find the largest of the two daughters
           !***

           chk_lvl1 = 2*final_lvl
           chk_lvl2 = 2*final_lvl+1
           if (chk_lvl2 >= lvl) chk_lvl2 = chk_lvl1

           if ((add1(chk_lvl1) >  add1(chk_lvl2)) .OR.    &
                ((add1(chk_lvl1) == add1(chk_lvl2)) .AND. &
                (add2(chk_lvl1) >  add2(chk_lvl2)))) then
              max_lvl = chk_lvl1
           else 
              max_lvl = chk_lvl2
           endif

           !***
           !*** if the parent is greater than both daughters,
           !*** the correct level has been found
           !***

           if ((add1_tmp >  add1(max_lvl)) .OR.    &
                ((add1_tmp == add1(max_lvl)) .AND. &
                (add2_tmp >  add2(max_lvl)))) then
              add1(final_lvl) = add1_tmp
              add2(final_lvl) = add2_tmp
              weights(:,final_lvl) = wgttmp(:)
              exit sift_loop2

              !***
              !*** otherwise, promote the largest daughter and push
              !*** down one level in the tree.  if haven't reached
              !*** the end of the tree, repeat the process.  otherwise
              !*** store last values and exit the loop
              !***

           else 
              add1(final_lvl) = add1(max_lvl)
              add2(final_lvl) = add2(max_lvl)
              weights(:,final_lvl) = weights(:,max_lvl)

              final_lvl = max_lvl
              if (2*final_lvl >= lvl) then
                 add1(final_lvl) = add1_tmp
                 add2(final_lvl) = add2_tmp
                 weights(:,final_lvl) = wgttmp(:)
                 exit sift_loop2
              endif
           endif
        end do sift_loop2
     end do

     !***
     !*** swap the last two entries
     !***


     add1_tmp = add1(2)
     add1(2)  = add1(1)
     add1(1)  = add1_tmp

     add2_tmp = add2(2)
     add2(2)  = add2(1)
     add2(1)  = add2_tmp

     wgttmp (:)   = weights(:,2)
     weights(:,2) = weights(:,1)
     weights(:,1) = wgttmp (:)

!#######################################################################
  end subroutine sort_add_v2
!#######################################################################

!/
!/ End of module SCRIP_INTERFACE -------------------------------------------- /
!/
      END MODULE SCRIP_INTERFACE


