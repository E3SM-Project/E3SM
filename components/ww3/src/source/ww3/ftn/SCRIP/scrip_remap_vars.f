!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     this module contains necessary variables for remapping between
!     two grids.  also routines for resizing and initializing these
!     variables.
!
!-----------------------------------------------------------------------
!
!     CVS:$Id: remap_vars.f,v 1.5 2000/04/19 21:56:26 pwjones Exp $
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
!     within WW3.
!
!***********************************************************************

      module scrip_remap_vars

      use SCRIP_KindsMod ! defines common data types
      use SCRIP_constants
      use scrip_grids

      implicit none

!-----------------------------------------------------------------------
!
!     module variables
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4), parameter ::
     &      norm_opt_none    = 1
     &,     norm_opt_dstarea = 2
     &,     norm_opt_frcarea = 3

      integer (SCRIP_i4), parameter ::
     &      map_type_conserv  = 1
     &,     map_type_bilinear = 2
     &,     map_type_bicubic  = 3
     &,     map_type_distwgt  = 4
     &,     map_type_particle = 5

      integer (SCRIP_i4), save :: 
     &      max_links_map1  ! current size of link arrays
     &,     num_links_map1  ! actual number of links for remapping
     &,     max_links_map2  ! current size of link arrays
     &,     num_links_map2  ! actual number of links for remapping
     &,     num_maps        ! num of remappings for this grid pair
     &,     num_wts         ! num of weights used in remapping
     &,     map_type        ! identifier for remapping method
     &,     norm_opt        ! option for normalization (conserv only)
     &,     resize_increment ! default amount to increase array size

      integer (SCRIP_i4), dimension(:), allocatable, save ::
     &      grid1_add_map1, ! grid1 address for each link in mapping 1
     &      grid2_add_map1, ! grid2 address for each link in mapping 1
     &      grid1_add_map2, ! grid1 address for each link in mapping 2
     &      grid2_add_map2  ! grid2 address for each link in mapping 2

      real (SCRIP_r8), dimension(:,:), allocatable, save ::
     &      wts_map1, ! map weights for each link (num_wts,max_links)
     &      wts_map2  ! map weights for each link (num_wts,max_links)

      real (kind = SCRIP_r8) :: frac_lowest,frac_highest
      real (kind = SCRIP_r8) :: wt_lowest,wt_highest

!***********************************************************************

      contains

!***********************************************************************

      subroutine init_remap_vars

!-----------------------------------------------------------------------
!
!     this routine initializes some variables and provides an initial
!     allocation of arrays (fairly large so frequent resizing 
!     unnecessary).
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     determine the number of weights
!
!-----------------------------------------------------------------------

      select case (map_type)
      case(map_type_conserv)
        num_wts = 3
      case(map_type_bilinear)
        num_wts = 1
      case(map_type_bicubic)
        num_wts = 4
      case(map_type_distwgt)
        num_wts = 1
      case(map_type_particle)
        num_wts = 1
      end select

!-----------------------------------------------------------------------
!
!     initialize num_links and set max_links to four times the largest 
!     of the destination grid sizes initially (can be changed later).
!     set a default resize increment to increase the size of link
!     arrays if the number of links exceeds the initial size
!   
!-----------------------------------------------------------------------

      num_links_map1 = 0
      max_links_map1 = 4*grid2_size
      if (num_maps > 1) then
        num_links_map2 = 0
        max_links_map1 = max(4*grid1_size,4*grid2_size)
        max_links_map2 = max_links_map1
      endif

      resize_increment = 0.1*max(grid1_size,grid2_size)

!-----------------------------------------------------------------------
!
!     allocate address and weight arrays for mapping 1
!   
!-----------------------------------------------------------------------

      allocate (grid1_add_map1(max_links_map1),
     &          grid2_add_map1(max_links_map1),
     &          wts_map1(num_wts, max_links_map1))

!-----------------------------------------------------------------------
!
!     allocate address and weight arrays for mapping 2 if necessary 
!   
!-----------------------------------------------------------------------

      if (num_maps > 1) then
        allocate (grid1_add_map2(max_links_map2),
     &            grid2_add_map2(max_links_map2),
     &            wts_map2(num_wts, max_links_map2))
      endif

!-----------------------------------------------------------------------

      end subroutine init_remap_vars

!***********************************************************************

      subroutine resize_remap_vars(nmap, increment)

!-----------------------------------------------------------------------
!
!     this routine resizes remapping arrays by increasing(decreasing)
!     the max_links by increment
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     input variables
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4), intent(in) ::
     &     nmap,      ! identifies which mapping array to resize
     &     increment  ! the number of links to add(subtract) to arrays

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4) ::
     &   ierr,     ! error flag
     &   mxlinks   ! size of link arrays

      integer (SCRIP_i4), dimension(:), allocatable ::
     &   add1_tmp, ! temp array for resizing address arrays
     &   add2_tmp  ! temp array for resizing address arrays

      real (SCRIP_r8), dimension(:,:), allocatable ::
     &   wts_tmp   ! temp array for resizing weight arrays

!-----------------------------------------------------------------------
!
!     resize map 1 arrays if required.
!
!-----------------------------------------------------------------------

      select case (nmap)
      case(1)

        !***
        !*** allocate temporaries to hold original values
        !***

        mxlinks = size(grid1_add_map1)
        allocate (add1_tmp(mxlinks), add2_tmp(mxlinks), 
     &            wts_tmp(num_wts,mxlinks))

        add1_tmp = grid1_add_map1
        add2_tmp = grid2_add_map1
        wts_tmp  = wts_map1
        
        !***
        !*** deallocate originals and increment max_links then
        !*** reallocate arrays at new size
        !***

        deallocate (grid1_add_map1, grid2_add_map1, wts_map1)
        max_links_map1 = mxlinks + increment
        allocate (grid1_add_map1(max_links_map1),
     &            grid2_add_map1(max_links_map1),
     &            wts_map1(num_wts,max_links_map1))

        !***
        !*** restore original values from temp arrays and
        !*** deallocate temps
        !***

        mxlinks = min(mxlinks, max_links_map1)
        grid1_add_map1(1:mxlinks) = add1_tmp (1:mxlinks)
        grid2_add_map1(1:mxlinks) = add2_tmp (1:mxlinks)
        wts_map1    (:,1:mxlinks) = wts_tmp(:,1:mxlinks)
        deallocate(add1_tmp, add2_tmp, wts_tmp)

!-----------------------------------------------------------------------
!
!     resize map 2 arrays if required.
!
!-----------------------------------------------------------------------

      case(2)

        !***
        !*** allocate temporaries to hold original values
        !***

        mxlinks = size(grid1_add_map2)
        allocate (add1_tmp(mxlinks), add2_tmp(mxlinks), 
     &            wts_tmp(num_wts,mxlinks),stat=ierr)
        if (ierr .ne. 0) then
          print *,'error allocating temps in resize: ',ierr
          stop
        endif

        add1_tmp = grid1_add_map2
        add2_tmp = grid2_add_map2
        wts_tmp  = wts_map2
        
        !***
        !*** deallocate originals and increment max_links then
        !*** reallocate arrays at new size
        !***

        deallocate (grid1_add_map2, grid2_add_map2, wts_map2)
        max_links_map2 = mxlinks + increment
        allocate (grid1_add_map2(max_links_map2),
     &            grid2_add_map2(max_links_map2),
     &            wts_map2(num_wts,max_links_map2),stat=ierr)
        if (ierr .ne. 0) then
          print *,'error allocating new arrays in resize: ',ierr
          stop
        endif


        !***
        !*** restore original values from temp arrays and
        !*** deallocate temps
        !***

        mxlinks = min(mxlinks, max_links_map2)
        grid1_add_map2(1:mxlinks) = add1_tmp (1:mxlinks)
        grid2_add_map2(1:mxlinks) = add2_tmp (1:mxlinks)
        wts_map2    (:,1:mxlinks) = wts_tmp(:,1:mxlinks)
        deallocate(add1_tmp, add2_tmp, wts_tmp)

      end select

!-----------------------------------------------------------------------

      end subroutine resize_remap_vars

!***********************************************************************

      end module scrip_remap_vars

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
