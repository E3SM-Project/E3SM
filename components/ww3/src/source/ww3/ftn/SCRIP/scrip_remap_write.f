!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     This module contains routines for writing the remapping data to 
!     a file.  Before writing the data for each mapping, the links are 
!     sorted by destination grid address.
!
!-----------------------------------------------------------------------
!
!     CVS:$Id: remap_write.f,v 1.7 2001/08/21 21:06:42 pwjones Exp $
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

      module scrip_remap_write

!-----------------------------------------------------------------------

      use SCRIP_KindsMod  ! defines common data types
      use SCRIP_ErrorMod  ! SCRIP  error handler
      use SCRIP_NetcdfMod ! netCDF error handler
      use netcdf          ! netCDF library

      use SCRIP_constants     ! defines common scalar constants
      use scrip_grids         ! module containing grid information
      use scrip_remap_vars    ! module containing remap information

      implicit none

!-----------------------------------------------------------------------
!
!     module variables
!
!-----------------------------------------------------------------------

      character(SCRIP_charLength), private :: 
     &   map_method       ! character string for map_type
     &,  normalize_opt    ! character string for normalization option
     &,  history          ! character string for history information
     &,  convention       ! character string for output convention

      character(8), private :: 
     &   cdate            ! character date string

      integer (SCRIP_i4), dimension(:), allocatable, private ::
     &   src_mask_int     ! integer masks to determine
     &,  dst_mask_int     ! cells that participate in map

!-----------------------------------------------------------------------
!
!     various netCDF identifiers used by output routines
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4), private ::
     &   ncstat               ! error flag for netCDF calls 
     &,  nc_file_id           ! id for netCDF file
     &,  nc_srcgrdsize_id     ! id for source grid size
     &,  nc_dstgrdsize_id     ! id for destination grid size
     &,  nc_srcgrdcorn_id     ! id for number of source grid corners
     &,  nc_dstgrdcorn_id     ! id for number of dest grid corners
     &,  nc_srcgrdrank_id     ! id for source grid rank
     &,  nc_dstgrdrank_id     ! id for dest grid rank
     &,  nc_numlinks_id       ! id for number of links in mapping
     &,  nc_numwgts_id        ! id for number of weights for mapping
     &,  nc_srcgrddims_id     ! id for source grid dimensions
     &,  nc_dstgrddims_id     ! id for dest grid dimensions
     &,  nc_srcgrdcntrlat_id  ! id for source grid center latitude
     &,  nc_dstgrdcntrlat_id  ! id for dest grid center latitude
     &,  nc_srcgrdcntrlon_id  ! id for source grid center longitude
     &,  nc_dstgrdcntrlon_id  ! id for dest grid center longitude
     &,  nc_srcgrdimask_id    ! id for source grid mask
     &,  nc_dstgrdimask_id    ! id for dest grid mask
     &,  nc_srcgrdcrnrlat_id  ! id for latitude of source grid corners
     &,  nc_srcgrdcrnrlon_id  ! id for longitude of source grid corners
     &,  nc_dstgrdcrnrlat_id  ! id for latitude of dest grid corners
     &,  nc_dstgrdcrnrlon_id  ! id for longitude of dest grid corners
     &,  nc_srcgrdarea_id     ! id for area of source grid cells
     &,  nc_dstgrdarea_id     ! id for area of dest grid cells
     &,  nc_srcgrdfrac_id     ! id for area fraction on source grid
     &,  nc_dstgrdfrac_id     ! id for area fraction on dest grid
     &,  nc_srcadd_id         ! id for map source address
     &,  nc_dstadd_id         ! id for map destination address
     &,  nc_rmpmatrix_id      ! id for remapping matrix

      integer (SCRIP_i4), dimension(2), private ::
     &   nc_dims2_id  ! netCDF ids for 2d array dims

!***********************************************************************

      contains

!***********************************************************************

      subroutine write_remap(map1_name, map2_name, interp_file1, 
     &                    interp_file2, output_opt, l_master, errorCode)

!-----------------------------------------------------------------------
!
!     calls correct output routine based on output format choice
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     input variables
!
!-----------------------------------------------------------------------

      character(SCRIP_charLength), intent(in) ::
     &            map1_name,    ! name for mapping grid1 to grid2
     &            map2_name,    ! name for mapping grid2 to grid1
     &            interp_file1, ! filename for map1 remap data
     &            interp_file2, ! filename for map2 remap data
     &            output_opt    ! option for output conventions

      logical, intent(in) ::
     &            l_master      ! Am I the master processor?

!-----------------------------------------------------------------------
!
!     output variables
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4), intent(out) ::
     &   errorCode              ! returned error code

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      character (11), parameter :: rtnName = 'write_remap'

!-----------------------------------------------------------------------
!
!     define some common variables to be used in all routines
!
!-----------------------------------------------------------------------

      errorCode = SCRIP_Success

      select case(norm_opt)
      case (norm_opt_none)
        normalize_opt = 'none'
      case (norm_opt_frcarea)
        normalize_opt = 'fracarea'
      case (norm_opt_dstarea)
        normalize_opt = 'destarea'
      end select

      select case(map_type)
      case(map_type_conserv)
        map_method = 'Conservative remapping'
      case(map_type_bilinear)
        map_method = 'Bilinear remapping'
      case(map_type_distwgt)
        map_method = 'Distance weighted avg of nearest neighbors'
      case(map_type_bicubic)
        map_method = 'Bicubic remapping'
      case(map_type_particle)
        map_method = 'Particle remapping'
      case default
         call SCRIP_ErrorSet(errorCode, rtnName, 'Invalid Map Type')
         return
      end select

      call date_and_time(date=cdate)
      write (history,1000) cdate(5:6),cdate(7:8),cdate(1:4)
 1000 format('Created: ',a2,'-',a2,'-',a4)

!-----------------------------------------------------------------------
!
!     sort address and weight arrays
!
!-----------------------------------------------------------------------

! New Apr4 2013: sort_add is instead called from scrip_interface.ftn
!     prior to entering this routine...no need to call it again here.

!     call sort_add(grid2_add_map1, grid1_add_map1, wts_map1)
      if (num_maps > 1) then
         call sort_add(grid1_add_map2, grid2_add_map2, wts_map2)
      endif

!-----------------------------------------------------------------------
!
!     call appropriate output routine
!
!-----------------------------------------------------------------------

      select case(output_opt)
      case ('scrip')
         if (l_master)  
     &   call write_remap_scrip(map1_name, interp_file1, 1, errorCode)
         if (SCRIP_ErrorCheck(errorCode, rtnName,
     &       'error in write_remap_scrip')) return
      case ('ncar-csm')
         call write_remap_csm  (map1_name, interp_file1, 1, errorCode)
         if (SCRIP_ErrorCheck(errorCode, rtnName,
     &       'error in write_remap_csm')) return
      case default
         call SCRIP_ErrorSet(errorCode, rtnName, 
     &                       'unknown output file convention')
         return
      end select

!-----------------------------------------------------------------------
!
!     call appropriate output routine for second mapping if required
!
!-----------------------------------------------------------------------

      if (num_maps > 1) then
        select case(output_opt)
        case ('scrip')
          if (l_master)  
     &    call write_remap_scrip(map2_name, interp_file2, 2, errorCode)
          if (SCRIP_ErrorCheck(errorCode, rtnName,
     &        'error in write_remap_scrip')) return
        case ('ncar-csm')
          call write_remap_csm  (map2_name, interp_file2, 2, errorCode)
          if (SCRIP_ErrorCheck(errorCode, rtnName,
     &        'error in write_remap_csm')) return
        case default
          call SCRIP_ErrorSet(errorCode, rtnName,
     &                        'unknown output file convention')
          return
        end select
      endif

!-----------------------------------------------------------------------

      end subroutine write_remap

!***********************************************************************

      subroutine write_remap_scrip(map_name, interp_file, direction,
     &                             errorCode)

!-----------------------------------------------------------------------
!
!     writes remap data to a netCDF file using SCRIP conventions
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     input variables
!
!-----------------------------------------------------------------------

      character(SCRIP_charLength), intent(in) ::
     &            map_name     ! name for mapping 
     &,           interp_file  ! filename for remap data

      integer (SCRIP_i4), intent(in) ::
     &  direction              ! direction of map (1=grid1 to grid2
                               !                   2=grid2 to grid1)

!-----------------------------------------------------------------------
!
!     output variables
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4), intent(out) ::
     &   errorCode              ! returned error code

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      character(SCRIP_charLength) ::
     &  grid1_ctmp        ! character temp for grid1 names
     &, grid2_ctmp        ! character temp for grid2 names

      integer (SCRIP_i4) ::
     &  itmp1             ! integer temp
     &, itmp2             ! integer temp
     &, itmp3             ! integer temp
     &, itmp4             ! integer temp

      character (17), parameter :: rtnName = 'write_remap_scrip'

!-----------------------------------------------------------------------
!
!     create netCDF file for mapping and define some global attributes
!
!-----------------------------------------------------------------------

      errorCode = SCRIP_Success

      ncstat = nf90_create (interp_file, NF90_CLOBBER, nc_file_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &                           'error creating remap file')) return

      !***
      !*** map name
      !***
      ncstat = nf90_put_att(nc_file_id, NF90_GLOBAL, 'title', map_name)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &                           'error writing remap name')) return

      !***
      !*** normalization option
      !***
      ncstat = nf90_put_att(nc_file_id, NF90_GLOBAL, 'normalization',
     &                      normalize_opt)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &                     'error writing normalize option')) return

      !***
      !*** map method
      !***
      ncstat = nf90_put_att(nc_file_id, NF90_GLOBAL, 'map_method',
     &                      map_method)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &                         'error writing remap method')) return

      !***
      !*** history
      !***
      ncstat = nf90_put_att(nc_file_id, NF90_GLOBAL, 'history',
     &                      history)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &                           'error writing history')) return

      !***
      !*** file convention
      !***
      convention = 'SCRIP'
      ncstat = nf90_put_att(nc_file_id, NF90_GLOBAL, 'conventions',
     &                                                convention)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &                  'error writing output convention')) return

      !***
      !*** source and destination grid names
      !***

      if (direction == 1) then
        grid1_ctmp = 'source_grid'
        grid2_ctmp = 'dest_grid'
      else
        grid1_ctmp = 'dest_grid'
        grid2_ctmp = 'source_grid'
      endif

      ncstat = nf90_put_att(nc_file_id, NF90_GLOBAL, trim(grid1_ctmp),
     &                      grid1_name)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &                      'error writing source grid name')) return

      ncstat = nf90_put_att(nc_file_id, NF90_GLOBAL, trim(grid2_ctmp),
     &                      grid2_name)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &                   'error writing destination grid name')) return

!-----------------------------------------------------------------------
!
!     prepare netCDF dimension info
!
!-----------------------------------------------------------------------

      !***
      !*** define grid size dimensions
      !***

      if (direction == 1) then
        itmp1 = grid1_size
        itmp2 = grid2_size
      else
        itmp1 = grid2_size
        itmp2 = grid1_size
      endif

      ncstat = nf90_def_dim(nc_file_id, 'src_grid_size', itmp1, 
     &                      nc_srcgrdsize_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &                 'error defining source grid size')) return

      ncstat = nf90_def_dim(nc_file_id, 'dst_grid_size', itmp2, 
     &                      nc_dstgrdsize_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &            'error defining destination grid size')) return

      !***
      !*** define grid corner dimension
      !***

      if (direction == 1) then
        itmp1 = grid1_corners
        itmp2 = grid2_corners
      else
        itmp1 = grid2_corners
        itmp2 = grid1_corners
      endif

      ncstat = nf90_def_dim(nc_file_id, 'src_grid_corners', 
     &                      itmp1, nc_srcgrdcorn_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining num corners on source grid')) return

      ncstat = nf90_def_dim(nc_file_id, 'dst_grid_corners', 
     &                      itmp2, nc_dstgrdcorn_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &     'error defining num corners on destination grid')) return

      !***
      !*** define grid rank dimension
      !***

      if (direction == 1) then
        itmp1 = grid1_rank
        itmp2 = grid2_rank
      else
        itmp1 = grid2_rank
        itmp2 = grid1_rank
      endif

      ncstat = nf90_def_dim(nc_file_id, 'src_grid_rank', 
     &                      itmp1, nc_srcgrdrank_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &                    'error defining source grid rank')) return

      ncstat = nf90_def_dim(nc_file_id, 'dst_grid_rank', 
     &                      itmp2, nc_dstgrdrank_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &               'error defining destination grid rank')) return

      !***
      !*** define map size dimensions
      !***

      if (direction == 1) then
        itmp1 = num_links_map1
      else
        itmp1 = num_links_map2
      endif

      ncstat = nf90_def_dim(nc_file_id, 'num_links', 
     &                      itmp1, nc_numlinks_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &                           'error defining remap size')) return

      ncstat = nf90_def_dim(nc_file_id, 'num_wgts', 
     &                      num_wts, nc_numwgts_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &                    'error defining number of weights')) return

      !***
      !*** define grid dimensions
      !***

      ncstat = nf90_def_var(nc_file_id, 'src_grid_dims', NF90_INT,
     &                      nc_srcgrdrank_id, nc_srcgrddims_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining source grid dims')) return

      ncstat = nf90_def_var(nc_file_id, 'dst_grid_dims', NF90_INT,
     &                      nc_dstgrdrank_id, nc_dstgrddims_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining destination grid dims')) return

!-----------------------------------------------------------------------
!
!     define all arrays for netCDF descriptors
!
!-----------------------------------------------------------------------

      !***
      !*** define grid center latitude array
      !***

      ncstat = nf90_def_var(nc_file_id, 'src_grid_center_lat', 
     &                      NF90_DOUBLE, nc_srcgrdsize_id, 
     &                      nc_srcgrdcntrlat_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining source grid center lat')) return

      ncstat = nf90_def_var(nc_file_id, 'dst_grid_center_lat', 
     &                      NF90_DOUBLE, nc_dstgrdsize_id, 
     &                      nc_dstgrdcntrlat_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining destination grid center lat')) return

      !***
      !*** define grid center longitude array
      !***

      ncstat = nf90_def_var(nc_file_id, 'src_grid_center_lon', 
     &                      NF90_DOUBLE, nc_srcgrdsize_id, 
     &                      nc_srcgrdcntrlon_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining source grid center lon')) return

      ncstat = nf90_def_var(nc_file_id, 'dst_grid_center_lon', 
     &                      NF90_DOUBLE, nc_dstgrdsize_id, 
     &                      nc_dstgrdcntrlon_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining destination grid center lon')) return

      !***
      !*** define grid corner lat/lon arrays
      !***

      nc_dims2_id(1) = nc_srcgrdcorn_id
      nc_dims2_id(2) = nc_srcgrdsize_id

      ncstat = nf90_def_var(nc_file_id, 'src_grid_corner_lat', 
     &                      NF90_DOUBLE, nc_dims2_id, 
     &                      nc_srcgrdcrnrlat_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining source grid corner lats')) return

      ncstat = nf90_def_var(nc_file_id, 'src_grid_corner_lon', 
     &                      NF90_DOUBLE, nc_dims2_id, 
     &                      nc_srcgrdcrnrlon_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining source grid corner lons')) return

      nc_dims2_id(1) = nc_dstgrdcorn_id
      nc_dims2_id(2) = nc_dstgrdsize_id

      ncstat = nf90_def_var(nc_file_id, 'dst_grid_corner_lat', 
     &                      NF90_DOUBLE, nc_dims2_id, 
     &                      nc_dstgrdcrnrlat_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining destination grid corner lats')) return

      ncstat = nf90_def_var(nc_file_id, 'dst_grid_corner_lon', 
     &                      NF90_DOUBLE, nc_dims2_id, 
     &                      nc_dstgrdcrnrlon_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining destination grid corner lons')) return

      !***
      !*** define units for all coordinate arrays
      !***

      if (direction == 1) then
        grid1_ctmp = grid1_units
        grid2_ctmp = grid2_units
      else
        grid1_ctmp = grid2_units
        grid2_ctmp = grid1_units
      endif

      ncstat = nf90_put_att(nc_file_id, nc_srcgrdcntrlat_id, 
     &                      'units', grid1_ctmp)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing source grid units')) return

      ncstat = nf90_put_att(nc_file_id, nc_dstgrdcntrlat_id, 
     &                      'units', grid2_ctmp)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing destination grid units')) return

      ncstat = nf90_put_att(nc_file_id, nc_srcgrdcntrlon_id, 
     &                      'units', grid1_ctmp)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing source grid units')) return

      ncstat = nf90_put_att(nc_file_id, nc_dstgrdcntrlon_id, 
     &                      'units', grid2_ctmp)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing destination grid units')) return

      ncstat = nf90_put_att(nc_file_id, nc_srcgrdcrnrlat_id, 
     &                      'units', grid1_ctmp)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing source grid units')) return

      ncstat = nf90_put_att(nc_file_id, nc_srcgrdcrnrlon_id, 
     &                      'units', grid1_ctmp)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing source grid units')) return

      ncstat = nf90_put_att(nc_file_id, nc_dstgrdcrnrlat_id, 
     &                      'units', grid2_ctmp)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing destination grid units')) return

      ncstat = nf90_put_att(nc_file_id, nc_dstgrdcrnrlon_id, 
     &                      'units', grid2_ctmp)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing destination grid units')) return

      !***
      !*** define grid mask
      !***

      ncstat = nf90_def_var(nc_file_id, 'src_grid_imask', NF90_INT,
     &                      nc_srcgrdsize_id, nc_srcgrdimask_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining source grid mask')) return

      ncstat = nf90_put_att(nc_file_id, nc_srcgrdimask_id, 
     &                      'units', 'unitless')
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing source grid mask units')) return

      ncstat = nf90_def_var(nc_file_id, 'dst_grid_imask', NF90_INT,
     &                      nc_dstgrdsize_id, nc_dstgrdimask_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining destination grid mask')) return

      ncstat = nf90_put_att(nc_file_id, nc_dstgrdimask_id, 
     &                      'units', 'unitless')
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing destination grid mask units')) return

      !***
      !*** define grid area arrays
      !***

      ncstat = nf90_def_var(nc_file_id, 'src_grid_area', 
     &                      NF90_DOUBLE, nc_srcgrdsize_id, 
     &                      nc_srcgrdarea_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining source grid area')) return

      ncstat = nf90_put_att(nc_file_id, nc_srcgrdarea_id, 
     &                      'units', 'square radians')
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing source area units')) return

      ncstat = nf90_def_var(nc_file_id, 'dst_grid_area', 
     &                      NF90_DOUBLE, nc_dstgrdsize_id, 
     &                      nc_dstgrdarea_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining destination grid area')) return

      ncstat = nf90_put_att(nc_file_id, nc_dstgrdarea_id, 
     &                      'units', 'square radians')
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing destination area units')) return

      !***
      !*** define grid fraction arrays
      !***

      ncstat = nf90_def_var(nc_file_id, 'src_grid_frac', 
     &                      NF90_DOUBLE, nc_srcgrdsize_id, 
     &                      nc_srcgrdfrac_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining source grid fraction')) return

      ncstat = nf90_put_att(nc_file_id, nc_srcgrdfrac_id, 
     &                      'units', 'unitless')
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing source fraction units')) return

      ncstat = nf90_def_var(nc_file_id, 'dst_grid_frac', 
     &                      NF90_DOUBLE, nc_dstgrdsize_id, 
     &                      nc_dstgrdfrac_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining destination fraction')) return

      ncstat = nf90_put_att(nc_file_id, nc_dstgrdfrac_id, 
     &                      'units', 'unitless')
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing destination frac units')) return

      !***
      !*** define mapping arrays
      !***

      ncstat = nf90_def_var(nc_file_id, 'src_address', 
     &                      NF90_INT, nc_numlinks_id, 
     &                      nc_srcadd_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining source addresses')) return

      ncstat = nf90_def_var(nc_file_id, 'dst_address', 
     &                      NF90_INT, nc_numlinks_id, 
     &                      nc_dstadd_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining destination addresses')) return

      nc_dims2_id(1) = nc_numwgts_id
      nc_dims2_id(2) = nc_numlinks_id

      ncstat = nf90_def_var(nc_file_id, 'remap_matrix', 
     &                      NF90_DOUBLE, nc_dims2_id, 
     &                      nc_rmpmatrix_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining remapping weights')) return

      !***
      !*** end definition stage
      !***

      ncstat = nf90_enddef(nc_file_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &                         'error ending definition phase')) return

!-----------------------------------------------------------------------
!
!     compute integer masks
!
!-----------------------------------------------------------------------

      if (direction == 1) then
        allocate (src_mask_int(grid1_size),
     &            dst_mask_int(grid2_size))

        where (grid2_mask)
          dst_mask_int = 1
        elsewhere
          dst_mask_int = 0
        endwhere

        where (grid1_mask)
          src_mask_int = 1
        elsewhere
          src_mask_int = 0
        endwhere
      else
        allocate (src_mask_int(grid2_size),
     &            dst_mask_int(grid1_size))

        where (grid1_mask)
          dst_mask_int = 1
        elsewhere
          dst_mask_int = 0
        endwhere

        where (grid2_mask)
          src_mask_int = 1
        elsewhere
          src_mask_int = 0
        endwhere
      endif

!-----------------------------------------------------------------------
!
!     change units of lat/lon coordinates if input units different
!     from radians
!
!-----------------------------------------------------------------------

      if (grid1_units(1:7) == 'degrees' .and. direction == 1) then
        grid1_center_lat = grid1_center_lat/deg2rad
        grid1_center_lon = grid1_center_lon/deg2rad
        grid1_corner_lat = grid1_corner_lat/deg2rad
        grid1_corner_lon = grid1_corner_lon/deg2rad
      endif

      if (grid2_units(1:7) == 'degrees' .and. direction == 1) then
        grid2_center_lat = grid2_center_lat/deg2rad
        grid2_center_lon = grid2_center_lon/deg2rad
        grid2_corner_lat = grid2_corner_lat/deg2rad
        grid2_corner_lon = grid2_corner_lon/deg2rad
      endif

!-----------------------------------------------------------------------
!
!     write mapping data
!
!-----------------------------------------------------------------------

      if (direction == 1) then
        itmp1 = nc_srcgrddims_id
        itmp2 = nc_dstgrddims_id
      else
        itmp2 = nc_srcgrddims_id
        itmp1 = nc_dstgrddims_id
      endif

      ncstat = nf90_put_var(nc_file_id, itmp1, grid1_dims)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing source grid dims')) return

      ncstat = nf90_put_var(nc_file_id, itmp2, grid2_dims)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing destination grid dims')) return

      ncstat = nf90_put_var(nc_file_id, nc_srcgrdimask_id, src_mask_int)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing source grid mask')) return

      ncstat = nf90_put_var(nc_file_id, nc_dstgrdimask_id, dst_mask_int)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing destination grid mask')) return

      deallocate(src_mask_int, dst_mask_int)

      if (direction == 1) then
        itmp1 = nc_srcgrdcntrlat_id
        itmp2 = nc_srcgrdcntrlon_id
        itmp3 = nc_srcgrdcrnrlat_id
        itmp4 = nc_srcgrdcrnrlon_id
      else
        itmp1 = nc_dstgrdcntrlat_id
        itmp2 = nc_dstgrdcntrlon_id
        itmp3 = nc_dstgrdcrnrlat_id
        itmp4 = nc_dstgrdcrnrlon_id
      endif

      ncstat = nf90_put_var(nc_file_id, itmp1, grid1_center_lat)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing source grid center lats')) return

      ncstat = nf90_put_var(nc_file_id, itmp2, grid1_center_lon)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing source grid center lons')) return

      ncstat = nf90_put_var(nc_file_id, itmp3, grid1_corner_lat)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing source grid corner lats')) return

      ncstat = nf90_put_var(nc_file_id, itmp4, grid1_corner_lon)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing source grid corner lons')) return

      if (direction == 1) then
        itmp1 = nc_dstgrdcntrlat_id
        itmp2 = nc_dstgrdcntrlon_id
        itmp3 = nc_dstgrdcrnrlat_id
        itmp4 = nc_dstgrdcrnrlon_id
      else
        itmp1 = nc_srcgrdcntrlat_id
        itmp2 = nc_srcgrdcntrlon_id
        itmp3 = nc_srcgrdcrnrlat_id
        itmp4 = nc_srcgrdcrnrlon_id
      endif

      ncstat = nf90_put_var(nc_file_id, itmp1, grid2_center_lat)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing destination grid center lats')) return

      ncstat = nf90_put_var(nc_file_id, itmp2, grid2_center_lon)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing destination grid center lons')) return

      ncstat = nf90_put_var(nc_file_id, itmp3, grid2_corner_lat)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing destination grid corner lats')) return

      ncstat = nf90_put_var(nc_file_id, itmp4, grid2_corner_lon)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing destination grid corner lons')) return

      if (direction == 1) then
        itmp1 = nc_srcgrdarea_id
        itmp2 = nc_srcgrdfrac_id
        itmp3 = nc_dstgrdarea_id
        itmp4 = nc_dstgrdfrac_id
      else
        itmp1 = nc_dstgrdarea_id
        itmp2 = nc_dstgrdfrac_id
        itmp3 = nc_srcgrdarea_id
        itmp4 = nc_srcgrdfrac_id
      endif

      if (luse_grid1_area) then
        ncstat = nf90_put_var(nc_file_id, itmp1, grid1_area_in)
      else
        ncstat = nf90_put_var(nc_file_id, itmp1, grid1_area)
      endif
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing grid1 area')) return

      ncstat = nf90_put_var(nc_file_id, itmp2, grid1_frac)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing grid1 frac')) return

      if (luse_grid2_area) then
        ncstat = nf90_put_var(nc_file_id, itmp3, grid2_area_in)
      else
        ncstat = nf90_put_var(nc_file_id, itmp3, grid2_area)
      endif
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing grid2 area')) return

      ncstat = nf90_put_var(nc_file_id, itmp4, grid2_frac)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing grid2 frac')) return

      if (direction == 1) then
        ncstat = nf90_put_var(nc_file_id, nc_srcadd_id, 
     &                        grid1_add_map1)
        if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &           'error writing source addresses')) return

        ncstat = nf90_put_var(nc_file_id, nc_dstadd_id, 
     &                        grid2_add_map1)
        if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &           'error writing destination addresses')) return

        ncstat = nf90_put_var(nc_file_id, nc_rmpmatrix_id, wts_map1)
        if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &           'error writing weights')) return
      else
        ncstat = nf90_put_var(nc_file_id, nc_srcadd_id, 
     &                        grid2_add_map2)
        if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &           'error writing source addresses')) return

        ncstat = nf90_put_var(nc_file_id, nc_dstadd_id, 
     &                          grid1_add_map2)
        if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &           'error writing destination addresses')) return

        ncstat = nf90_put_var(nc_file_id, nc_rmpmatrix_id, wts_map2)
        if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &           'error writing weights')) return
      endif

      ncstat = nf90_close(nc_file_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &                           'error closing file')) return

!-----------------------------------------------------------------------

      end subroutine write_remap_scrip

!***********************************************************************

      subroutine write_remap_ww3(map1_name, interp_file1, 
     &                           output_opt, l_master, errorCode)

!-----------------------------------------------------------------------
!
!     calls correct output routine for WW3 based on output format choice
!
!     only output variables needed for WW3 remapping
!
!-----------------------------------------------------------------------
      use netcdf

!-----------------------------------------------------------------------
!
!     input variables
!
!-----------------------------------------------------------------------

      character(SCRIP_charLength), intent(in) ::
     &            map1_name,    ! name for mapping grid1 to grid2
     &            interp_file1, ! filename for map1 remap data
     &            output_opt    ! option for output conventions
      logical, intent(in) ::
     &            l_master      ! Am I the master processor?

!-----------------------------------------------------------------------
!
!     output variables
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4), intent(out) ::
     &   errorCode              ! returned error code

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      character (15), parameter :: rtnName = 'write_remap_ww3'
      character(SCRIP_charLength) ::
     &            map1_name_pass,    ! name for mapping grid1 to grid2
     &            interp_file1_pass  ! filename for map1 remap data

!-----------------------------------------------------------------------
!
!     define some common variables to be used in all routines
!
!-----------------------------------------------------------------------

      errorCode = SCRIP_Success
      map1_name_pass = map1_name
      interp_file1_pass = interp_file1

      select case(norm_opt)
      case (norm_opt_none)
        normalize_opt = 'none'
      case (norm_opt_frcarea)
        normalize_opt = 'fracarea'
      case (norm_opt_dstarea)
        normalize_opt = 'destarea'
      end select

      select case(map_type)
      case(map_type_conserv)
        map_method = 'Conservative remapping'
      case(map_type_bilinear)
        map_method = 'Bilinear remapping'
      case(map_type_distwgt)
        map_method = 'Distance weighted avg of nearest neighbors'
      case(map_type_bicubic)
        map_method = 'Bicubic remapping'
      case(map_type_particle)
        map_method = 'Particle remapping'
      case default
         call SCRIP_ErrorSet(errorCode, rtnName, 'Invalid Map Type')
         return
      end select

      call date_and_time(date=cdate)
      write (history,1000) cdate(5:6),cdate(7:8),cdate(1:4)
 1000 format('Created: ',a2,'-',a2,'-',a4)

!-----------------------------------------------------------------------
!
!     sort address and weight arrays
!
!-----------------------------------------------------------------------

      call sort_add(grid2_add_map1, grid1_add_map1, wts_map1)
      if (num_maps > 1) then
        call sort_add(grid1_add_map2, grid2_add_map2, wts_map2)
      endif

!-----------------------------------------------------------------------
!
!     call appropriate output routine
!
!-----------------------------------------------------------------------

      select case(output_opt)
      case ('scrip')
         if (l_master) then
            call write_remap_scrip_ww3(map1_name_pass, 
     &           interp_file1_pass, errorCode)
         endif
         if (SCRIP_ErrorCheck(errorCode, rtnName,
     &       'error in write_remap_scrip')) return
      case default
         call SCRIP_ErrorSet(errorCode, rtnName, 
     &                       'unknown output file convention')
         return
      end select

!-----------------------------------------------------------------------

      end subroutine write_remap_ww3

!***********************************************************************

      subroutine write_remap_scrip_ww3(map_name, interp_file, errorCode)

!-----------------------------------------------------------------------
!
!     writes remap data to a netCDF file using SCRIP conventions
!     only writes variables needed for WW3 remap
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     input variables
!
!-----------------------------------------------------------------------

      character(SCRIP_charLength), intent(in) ::
     &            map_name     ! name for mapping 
     &,           interp_file  ! filename for remap data

!-----------------------------------------------------------------------
!
!     output variables
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4), intent(out) ::
     &   errorCode              ! returned error code

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      character(SCRIP_charLength) ::
     &  grid1_ctmp        ! character temp for grid1 names
     &, grid2_ctmp        ! character temp for grid2 names
     &, map_name_ctmp     ! character temp for name for mapping 
     &, interp_file_ctmp  ! character temp filename for remap data

      integer (SCRIP_i4) ::
     &  nc_file_id        ! netCDF file id
     &, itmp1             ! integer temp
     &, itmp2             ! integer temp
     &, itmp3             ! integer temp
     &, itmp4             ! integer temp

      character (21), parameter :: rtnName = 'write_remap_scrip_ww3'

!-----------------------------------------------------------------------
!
!     create netCDF file for mapping and define some global attributes
!
!-----------------------------------------------------------------------

      errorCode = SCRIP_Success
      map_name_ctmp = map_name
      interp_file_ctmp = trim(interp_file)

      itmp1 = len(interp_file_ctmp)
      itmp2 = len_trim(interp_file_ctmp)
      ncstat = nf90_create(interp_file_ctmp, NF90_CLOBBER, nc_file_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &                           'error creating remap file')) return

      !***
      !*** map name
      !***
      ncstat = nf90_put_att(nc_file_id, NF90_GLOBAL, 'title', map_name)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &                           'error writing remap name')) return

      !***
      !*** normalization option
      !***
      ncstat = nf90_put_att(nc_file_id, NF90_GLOBAL, 'normalization',
     &                      normalize_opt)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &                     'error writing normalize option')) return

      !***
      !*** map method
      !***
      ncstat = nf90_put_att(nc_file_id, NF90_GLOBAL, 'map_method',
     &                      map_method)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &                         'error writing remap method')) return

      !***
      !*** history
      !***
      ncstat = nf90_put_att(nc_file_id, NF90_GLOBAL, 'history',
     &                      history)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &                           'error writing history')) return

      !***
      !*** file convention
      !***
      convention = 'SCRIP'
      ncstat = nf90_put_att(nc_file_id, NF90_GLOBAL, 'conventions',
     &                                                convention)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &                  'error writing output convention')) return

      !***
      !*** source and destination grid names
      !***

      grid1_ctmp = 'source_grid'
      grid2_ctmp = 'dest_grid'

      ncstat = nf90_put_att(nc_file_id, NF90_GLOBAL, trim(grid1_ctmp),
     &                      grid1_name)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &                      'error writing source grid name')) return

      ncstat = nf90_put_att(nc_file_id, NF90_GLOBAL, trim(grid2_ctmp),
     &                      grid2_name)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &                   'error writing destination grid name')) return

!-----------------------------------------------------------------------
!
!     prepare netCDF dimension info
!
!-----------------------------------------------------------------------

      !***
      !*** define grid size dimensions
      !***

      itmp2 = grid2_size

      ncstat = nf90_def_dim(nc_file_id, 'dst_grid_size', itmp2, 
     &                      nc_dstgrdsize_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &            'error defining destination grid size')) return

      !***
      !*** define map size dimensions
      !***

      itmp1 = num_links_map1

      ncstat = nf90_def_dim(nc_file_id, 'num_links', 
     &                      itmp1, nc_numlinks_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &                           'error defining remap size')) return

      ncstat = nf90_def_dim(nc_file_id, 'num_wgts', 
     &                      num_wts, nc_numwgts_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &                    'error defining number of weights')) return

!-----------------------------------------------------------------------
!
!     define all arrays for netCDF descriptors
!
!-----------------------------------------------------------------------


      !***
      !*** define grid fraction arrays
      !***


      ncstat = nf90_def_var(nc_file_id, 'dst_grid_frac', 
     &                      NF90_DOUBLE, nc_dstgrdsize_id, 
     &                      nc_dstgrdfrac_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining destination fraction')) return

      ncstat = nf90_put_att(nc_file_id, nc_dstgrdfrac_id, 
     &                      'units', 'unitless')
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing destination frac units')) return

      !***
      !*** define mapping arrays
      !***

      ncstat = nf90_def_var(nc_file_id, 'src_address', 
     &                      NF90_INT, nc_numlinks_id, 
     &                      nc_srcadd_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining source addresses')) return

      ncstat = nf90_def_var(nc_file_id, 'dst_address', 
     &                      NF90_INT, nc_numlinks_id, 
     &                      nc_dstadd_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining destination addresses')) return

      nc_dims2_id(1) = nc_numwgts_id
      nc_dims2_id(2) = nc_numlinks_id

      ncstat = nf90_def_var(nc_file_id, 'remap_matrix', 
     &                      NF90_DOUBLE, nc_dims2_id, 
     &                      nc_rmpmatrix_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining remapping weights')) return

      !***
      !*** end definition stage
      !***

      ncstat = nf90_enddef(nc_file_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &                         'error ending definition phase')) return


      !***
      !*** write variable arrays
      !***

      itmp4 = nc_dstgrdfrac_id

      ncstat = nf90_put_var(nc_file_id, itmp4, grid2_frac)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing grid2 frac')) return

      ncstat = nf90_put_var(nc_file_id, nc_srcadd_id, 
     &                        grid1_add_map1)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &           'error writing source addresses')) return

      ncstat = nf90_put_var(nc_file_id, nc_dstadd_id, 
     &                        grid2_add_map1)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &           'error writing destination addresses')) return

      ncstat = nf90_put_var(nc_file_id, nc_rmpmatrix_id, wts_map1)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &           'error writing weights')) return


      ncstat = nf90_close(nc_file_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &                           'error closing file')) return


!-----------------------------------------------------------------------

      end subroutine write_remap_scrip_ww3

!***********************************************************************

      subroutine write_remap_csm(map_name, interp_file, direction,
     &                           errorCode)

!-----------------------------------------------------------------------
!
!     writes remap data to a netCDF file using NCAR-CSM conventions
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     input variables
!
!-----------------------------------------------------------------------

      character(SCRIP_charLength), intent(in) ::
     &            map_name     ! name for mapping 
     &,           interp_file  ! filename for remap data

      integer (SCRIP_i4), intent(in) ::
     &  direction              ! direction of map (1=grid1 to grid2
                               !                   2=grid2 to grid1)

!-----------------------------------------------------------------------
!
!     output variables
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4), intent(out) ::
     &   errorCode              ! returned error code

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      character(SCRIP_charLength) ::
     &  grid1_ctmp        ! character temp for grid1 names
     &, grid2_ctmp        ! character temp for grid2 names

      integer (SCRIP_i4) ::
     &  itmp1             ! integer temp
     &, itmp2             ! integer temp
     &, itmp3             ! integer temp
     &, itmp4             ! integer temp
     &, nc_numwgts1_id    ! extra netCDF id for additional weights
     &, nc_src_isize_id   ! extra netCDF id for ni_a
     &, nc_src_jsize_id   ! extra netCDF id for nj_a
     &, nc_dst_isize_id   ! extra netCDF id for ni_b
     &, nc_dst_jsize_id   ! extra netCDF id for nj_b
     &, nc_rmpmatrix2_id  ! extra netCDF id for high-order remap matrix

      real (SCRIP_r8), dimension(:),allocatable ::
     &  wts1              ! CSM wants single array for 1st-order wts

      real (SCRIP_r8), dimension(:,:),allocatable ::
     &  wts2              ! write remaining weights in different array

      character (15), parameter :: rtnName = 'write_remap_csm'

!-----------------------------------------------------------------------
!
!     create netCDF file for mapping and define some global attributes
!
!-----------------------------------------------------------------------

      errorCode = SCRIP_Success

      ncstat = nf90_create(interp_file, NF90_CLOBBER, nc_file_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &                           'error opening file')) return

      !***
      !*** map name
      !***
      ncstat = nf90_put_att(nc_file_id, NF90_GLOBAL, 'title', map_name)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &                           'error writing remap name')) return

      !***
      !*** normalization option
      !***
      ncstat = nf90_put_att(nc_file_id, NF90_GLOBAL, 'normalization',
     &                      normalize_opt)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing normalization option')) return

      !***
      !*** map method
      !***
      ncstat = nf90_put_att(nc_file_id, NF90_GLOBAL, 'map_method',
     &                      map_method)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing remap method')) return

      !***
      !*** history
      !***
      ncstat = nf90_put_att(nc_file_id, NF90_GLOBAL, 'history',
     &                      history)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &                           'error writing history')) return

      !***
      !*** file convention
      !***
      convention = 'NCAR-CSM'
      ncstat = nf90_put_att(nc_file_id, NF90_GLOBAL, 'conventions',
     &                      convention)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing output convention')) return

      !***
      !*** source and destination grid names
      !***

      if (direction == 1) then
        grid1_ctmp = 'domain_a'
        grid2_ctmp = 'domain_b'
      else
        grid1_ctmp = 'domain_b'
        grid2_ctmp = 'domain_a'
      endif

      ncstat = nf90_put_att(nc_file_id, NF90_GLOBAL, trim(grid1_ctmp),
     &                      grid1_name)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &                           'error writing grid1 name')) return

      ncstat = nf90_put_att(nc_file_id, NF90_GLOBAL, trim(grid2_ctmp),
     &                      grid2_name)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &                           'error writing grid2 name')) return

!-----------------------------------------------------------------------
!
!     prepare netCDF dimension info
!
!-----------------------------------------------------------------------

      !***
      !*** define grid size dimensions
      !***

      if (direction == 1) then
        itmp1 = grid1_size
        itmp2 = grid2_size
      else
        itmp1 = grid2_size
        itmp2 = grid1_size
      endif

      ncstat = nf90_def_dim(nc_file_id, 'n_a', itmp1, nc_srcgrdsize_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining source grid size')) return

      ncstat = nf90_def_dim(nc_file_id, 'n_b', itmp2, nc_dstgrdsize_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining destination grid size')) return

      !***
      !*** define grid corner dimension
      !***

      if (direction == 1) then
        itmp1 = grid1_corners
        itmp2 = grid2_corners
      else
        itmp1 = grid2_corners
        itmp2 = grid1_corners
      endif

      ncstat = nf90_def_dim(nc_file_id, 'nv_a', itmp1, nc_srcgrdcorn_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &       'error defining number of corners on source grid')) return

      ncstat = nf90_def_dim(nc_file_id, 'nv_b', itmp2, nc_dstgrdcorn_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &  'error defining number of corners on destination grid')) return

      !***
      !*** define grid rank dimension
      !***

      if (direction == 1) then
        itmp1 = grid1_rank
        itmp2 = grid2_rank
      else
        itmp1 = grid2_rank
        itmp2 = grid1_rank
      endif

      ncstat = nf90_def_dim(nc_file_id, 'src_grid_rank', 
     &                      itmp1, nc_srcgrdrank_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining source grid rank')) return

      ncstat = nf90_def_dim(nc_file_id, 'dst_grid_rank', 
     &                      itmp2, nc_dstgrdrank_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining destination grid rank')) return

      !***
      !*** define first two dims as if 2-d cartesian domain
      !***

      if (direction == 1) then
        itmp1 = grid1_dims(1)
        if (grid1_rank > 1) then
          itmp2 = grid1_dims(2)
        else
          itmp2 = 0
        endif
        itmp3 = grid2_dims(1)
        if (grid2_rank > 1) then
          itmp4 = grid2_dims(2)
        else
          itmp4 = 0
        endif
      else
        itmp1 = grid2_dims(1)
        if (grid2_rank > 1) then
          itmp2 = grid2_dims(2)
        else
          itmp2 = 0
        endif
        itmp3 = grid1_dims(1)
        if (grid1_rank > 1) then
          itmp4 = grid1_dims(2)
        else
          itmp4 = 0
        endif
      endif

      ncstat = nf90_def_dim(nc_file_id, 'ni_a', itmp1, nc_src_isize_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining source isize')) return

      ncstat = nf90_def_dim(nc_file_id, 'nj_a', itmp2, nc_src_jsize_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining source jsize')) return

      ncstat = nf90_def_dim(nc_file_id, 'ni_b', itmp3, nc_dst_isize_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining destination isize')) return

      ncstat = nf90_def_dim(nc_file_id, 'nj_b', itmp4, nc_dst_jsize_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining destination jsize')) return

      !***
      !*** define map size dimensions
      !***

      if (direction == 1) then
        itmp1 = num_links_map1
      else
        itmp1 = num_links_map2
      endif

      ncstat = nf90_def_dim(nc_file_id, 'n_s', itmp1, nc_numlinks_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining remap size')) return

      ncstat = nf90_def_dim(nc_file_id, 'num_wgts', 
     &                     num_wts, nc_numwgts_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining number of weights')) return

      if (num_wts > 1) then
        ncstat = nf90_def_dim(nc_file_id, 'num_wgts1', 
     &                       num_wts-1, nc_numwgts1_id)
        if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &           'error defining number of weights1')) return
      endif

      !***
      !*** define grid dimensions
      !***

      ncstat = nf90_def_var(nc_file_id, 'src_grid_dims', NF90_INT,
     &                      nc_srcgrdrank_id, nc_srcgrddims_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining source grid dims')) return

      ncstat = nf90_def_var(nc_file_id, 'dst_grid_dims', NF90_INT,
     &                      nc_dstgrdrank_id, nc_dstgrddims_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining destination grid dims')) return

!-----------------------------------------------------------------------
!
!     define all arrays for netCDF descriptors
!
!-----------------------------------------------------------------------

      !***
      !*** define grid center latitude array
      !***

      ncstat = nf90_def_var(nc_file_id, 'yc_a',
     &                      NF90_DOUBLE, nc_srcgrdsize_id, 
     &                      nc_srcgrdcntrlat_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining source grid center lats')) return

      ncstat = nf90_def_var(nc_file_id, 'yc_b', 
     &                      NF90_DOUBLE, nc_dstgrdsize_id, 
     &                      nc_dstgrdcntrlat_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining destination grid center lats')) return

      !***
      !*** define grid center longitude array
      !***

      ncstat = nf90_def_var(nc_file_id, 'xc_a', 
     &                      NF90_DOUBLE, nc_srcgrdsize_id, 
     &                      nc_srcgrdcntrlon_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining source grid center lons')) return

      ncstat = nf90_def_var(nc_file_id, 'xc_b', 
     &                      NF90_DOUBLE, nc_dstgrdsize_id, 
     &                      nc_dstgrdcntrlon_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining destination grid center lons')) return

      !***
      !*** define grid corner lat/lon arrays
      !***

      nc_dims2_id(1) = nc_srcgrdcorn_id
      nc_dims2_id(2) = nc_srcgrdsize_id

      ncstat = nf90_def_var(nc_file_id, 'yv_a', 
     &                      NF90_DOUBLE, nc_dims2_id, 
     &                      nc_srcgrdcrnrlat_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining source grid corner lats')) return

      ncstat = nf90_def_var(nc_file_id, 'xv_a', 
     &                      NF90_DOUBLE, nc_dims2_id, 
     &                      nc_srcgrdcrnrlon_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining source grid corner lons')) return

      nc_dims2_id(1) = nc_dstgrdcorn_id
      nc_dims2_id(2) = nc_dstgrdsize_id

      ncstat = nf90_def_var(nc_file_id, 'yv_b', 
     &                      NF90_DOUBLE, nc_dims2_id, 
     &                      nc_dstgrdcrnrlat_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining destination grid corner lats')) return

      ncstat = nf90_def_var(nc_file_id, 'xv_b', 
     &                      NF90_DOUBLE, nc_dims2_id, 
     &                      nc_dstgrdcrnrlon_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining destination grid corner lons')) return

      !***
      !*** CSM wants all in degrees
      !***

      grid1_units = 'degrees'
      grid2_units = 'degrees'

      if (direction == 1) then
        grid1_ctmp = grid1_units
        grid2_ctmp = grid2_units
      else
        grid1_ctmp = grid2_units
        grid2_ctmp = grid1_units
      endif

      ncstat = nf90_put_att(nc_file_id, nc_srcgrdcntrlat_id, 
     &                      'units', grid1_ctmp)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing grid units')) return

      ncstat = nf90_put_att(nc_file_id, nc_dstgrdcntrlat_id, 
     &                      'units', grid2_ctmp)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing grid units')) return

      ncstat = nf90_put_att(nc_file_id, nc_srcgrdcntrlon_id, 
     &                      'units', grid1_ctmp)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing grid units')) return

      ncstat = nf90_put_att(nc_file_id, nc_dstgrdcntrlon_id, 
     &                      'units', grid2_ctmp)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing grid units')) return

      ncstat = nf90_put_att(nc_file_id, nc_srcgrdcrnrlat_id, 
     &                      'units', grid1_ctmp)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing grid units')) return

      ncstat = nf90_put_att(nc_file_id, nc_srcgrdcrnrlon_id, 
     &                      'units', grid1_ctmp)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing grid units')) return

      ncstat = nf90_put_att(nc_file_id, nc_dstgrdcrnrlat_id, 
     &                      'units', grid2_ctmp)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing grid units')) return

      ncstat = nf90_put_att(nc_file_id, nc_dstgrdcrnrlon_id, 
     &                      'units', grid2_ctmp)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing grid units')) return

      !***
      !*** define grid mask
      !***

      ncstat = nf90_def_var(nc_file_id, 'mask_a', NF90_INT,
     &                      nc_srcgrdsize_id, nc_srcgrdimask_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining source grid mask')) return

      ncstat = nf90_put_att(nc_file_id, nc_srcgrdimask_id, 
     &                      'units', 'unitless')
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing source mask units')) return

      ncstat = nf90_def_var(nc_file_id, 'mask_b', NF90_INT,
     &                      nc_dstgrdsize_id, nc_dstgrdimask_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining destination grid mask')) return

      ncstat = nf90_put_att(nc_file_id, nc_dstgrdimask_id, 
     &                      'units', 'unitless')
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing destination mask units')) return

      !***
      !*** define grid area arrays
      !***

      ncstat = nf90_def_var(nc_file_id, 'area_a', 
     &                      NF90_DOUBLE, nc_srcgrdsize_id, 
     &                      nc_srcgrdarea_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining source grid area')) return

      ncstat = nf90_put_att(nc_file_id, nc_srcgrdarea_id, 
     &                      'units', 'square radians')
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining source area units')) return

      ncstat = nf90_def_var(nc_file_id, 'area_b', 
     &                      NF90_DOUBLE, nc_dstgrdsize_id, 
     &                      nc_dstgrdarea_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining destination grid area')) return

      ncstat = nf90_put_att(nc_file_id, nc_dstgrdarea_id, 
     &                      'units', 'square radians')
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining destination area units')) return

      !***
      !*** define grid fraction arrays
      !***

      ncstat = nf90_def_var(nc_file_id, 'frac_a', 
     &                      NF90_DOUBLE, nc_srcgrdsize_id, 
     &                      nc_srcgrdfrac_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining source grid frac')) return

      ncstat = nf90_put_att(nc_file_id, nc_srcgrdfrac_id, 
     &                      'units', 'unitless')
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining source frac units')) return

      ncstat = nf90_def_var(nc_file_id, 'frac_b', 
     &                      NF90_DOUBLE, nc_dstgrdsize_id, 
     &                      nc_dstgrdfrac_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining destination grid frac')) return

      ncstat = nf90_put_att(nc_file_id, nc_dstgrdfrac_id, 
     &                      'units', 'unitless')
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining destination frac units')) return

      !***
      !*** define mapping arrays
      !***

      ncstat = nf90_def_var(nc_file_id, 'col', 
     &                      NF90_INT, nc_numlinks_id, 
     &                      nc_srcadd_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining source addresses')) return

      ncstat = nf90_def_var(nc_file_id, 'row', 
     &                      NF90_INT, nc_numlinks_id, 
     &                      nc_dstadd_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error defining destination addresses')) return

      ncstat = nf90_def_var(nc_file_id, 'S', 
     &                      NF90_DOUBLE, nc_numlinks_id, 
     &                      nc_rmpmatrix_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &                      'error defining weights')) return

      if (num_wts > 1) then
        nc_dims2_id(1) = nc_numwgts1_id
        nc_dims2_id(2) = nc_numlinks_id

        ncstat = nf90_def_var(nc_file_id, 'S2', 
     &                     NF90_DOUBLE, nc_dims2_id, 
     &                     nc_rmpmatrix2_id)
        if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &                     'error defining weights2')) return
      endif

      !***
      !*** end definition stage
      !***

      ncstat = nf90_enddef(nc_file_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &                       'error ending definition phase')) return

!-----------------------------------------------------------------------
!
!     compute integer masks
!
!-----------------------------------------------------------------------

      if (direction == 1) then
        allocate (src_mask_int(grid1_size),
     &            dst_mask_int(grid2_size))

        where (grid2_mask)
          dst_mask_int = 1
        elsewhere
          dst_mask_int = 0
        endwhere

        where (grid1_mask)
          src_mask_int = 1
        elsewhere
          src_mask_int = 0
        endwhere
      else
        allocate (src_mask_int(grid2_size),
     &            dst_mask_int(grid1_size))

        where (grid1_mask)
          dst_mask_int = 1
        elsewhere
          dst_mask_int = 0
        endwhere

        where (grid2_mask)
          src_mask_int = 1
        elsewhere
          src_mask_int = 0
        endwhere
      endif

!-----------------------------------------------------------------------
!
!     change units of lat/lon coordinates if input units different
!     from radians. if this is the second mapping, the conversion has
!     alread been done.
!
!-----------------------------------------------------------------------

      if (grid1_units(1:7) == 'degrees' .and. direction == 1) then
        grid1_center_lat = grid1_center_lat/deg2rad
        grid1_center_lon = grid1_center_lon/deg2rad
        grid1_corner_lat = grid1_corner_lat/deg2rad
        grid1_corner_lon = grid1_corner_lon/deg2rad
      endif

      if (grid2_units(1:7) == 'degrees' .and. direction == 1) then
        grid2_center_lat = grid2_center_lat/deg2rad
        grid2_center_lon = grid2_center_lon/deg2rad
        grid2_corner_lat = grid2_corner_lat/deg2rad
        grid2_corner_lon = grid2_corner_lon/deg2rad
      endif

!-----------------------------------------------------------------------
!
!     write mapping data
!
!-----------------------------------------------------------------------

      if (direction == 1) then
        itmp1 = nc_srcgrddims_id
        itmp2 = nc_dstgrddims_id
      else
        itmp2 = nc_srcgrddims_id
        itmp1 = nc_dstgrddims_id
      endif

      ncstat = nf90_put_var(nc_file_id, itmp1, grid1_dims)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing grid1 dims')) return

      ncstat = nf90_put_var(nc_file_id, itmp2, grid2_dims)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing grid2 dims')) return

      ncstat = nf90_put_var(nc_file_id, nc_srcgrdimask_id, 
     &                      src_mask_int)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing source grid mask')) return

      ncstat = nf90_put_var(nc_file_id, nc_dstgrdimask_id,
     &                      dst_mask_int)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing destination grid mask')) return

      deallocate(src_mask_int, dst_mask_int)

      if (direction == 1) then
        itmp1 = nc_srcgrdcntrlat_id
        itmp2 = nc_srcgrdcntrlon_id
        itmp3 = nc_srcgrdcrnrlat_id
        itmp4 = nc_srcgrdcrnrlon_id
      else
        itmp1 = nc_dstgrdcntrlat_id
        itmp2 = nc_dstgrdcntrlon_id
        itmp3 = nc_dstgrdcrnrlat_id
        itmp4 = nc_dstgrdcrnrlon_id
      endif

      ncstat = nf90_put_var(nc_file_id, itmp1, grid1_center_lat)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing grid1 center lats')) return

      ncstat = nf90_put_var(nc_file_id, itmp2, grid1_center_lon)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing grid1 center lons')) return

      ncstat = nf90_put_var(nc_file_id, itmp3, grid1_corner_lat)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing grid1 corner lats')) return

      ncstat = nf90_put_var(nc_file_id, itmp4, grid1_corner_lon)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing grid1 corner lons')) return

      if (direction == 1) then
        itmp1 = nc_dstgrdcntrlat_id
        itmp2 = nc_dstgrdcntrlon_id
        itmp3 = nc_dstgrdcrnrlat_id
        itmp4 = nc_dstgrdcrnrlon_id
      else
        itmp1 = nc_srcgrdcntrlat_id
        itmp2 = nc_srcgrdcntrlon_id
        itmp3 = nc_srcgrdcrnrlat_id
        itmp4 = nc_srcgrdcrnrlon_id
      endif

      ncstat = nf90_put_var(nc_file_id, itmp1, grid2_center_lat)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing grid2 center lats')) return

      ncstat = nf90_put_var(nc_file_id, itmp2, grid2_center_lon)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing grid2 center lons')) return

      ncstat = nf90_put_var(nc_file_id, itmp3, grid2_corner_lat)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing grid2 corner lats')) return

      ncstat = nf90_put_var(nc_file_id, itmp4, grid2_corner_lon)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing grid2 corner lons')) return

      if (direction == 1) then
        itmp1 = nc_srcgrdarea_id
        itmp2 = nc_srcgrdfrac_id
        itmp3 = nc_dstgrdarea_id
        itmp4 = nc_dstgrdfrac_id
      else
        itmp1 = nc_dstgrdarea_id
        itmp2 = nc_dstgrdfrac_id
        itmp3 = nc_srcgrdarea_id
        itmp4 = nc_srcgrdfrac_id
      endif

      if (luse_grid1_area) then
        ncstat = nf90_put_var(nc_file_id, itmp1, grid1_area_in)
      else
        ncstat = nf90_put_var(nc_file_id, itmp1, grid1_area)
      endif
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing grid1 area')) return

      ncstat = nf90_put_var(nc_file_id, itmp2, grid1_frac)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing grid1 frac')) return

      if (luse_grid2_area) then
        ncstat = nf90_put_var(nc_file_id, itmp3, grid2_area)
      else
        ncstat = nf90_put_var(nc_file_id, itmp3, grid2_area)
      endif
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing grid2 area')) return

      ncstat = nf90_put_var(nc_file_id, itmp4, grid2_frac)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error writing grid2 frac')) return

      if (direction == 1) then
        ncstat = nf90_put_var(nc_file_id, nc_srcadd_id, 
     &                        grid1_add_map1)
        if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &           'error writing source addresses')) return

        ncstat = nf90_put_var(nc_file_id, nc_dstadd_id, 
     &                        grid2_add_map1)
        if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &           'error writing destination addresses')) return

        if (num_wts == 1) then
          ncstat = nf90_put_var(nc_file_id, nc_rmpmatrix_id, 
     &                               wts_map1(1,:))
          if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &             'error writing weights')) return
        else
          allocate(wts1(num_links_map1),wts2(num_wts-1,num_links_map1))

          wts1 = wts_map1(1,:)
          wts2 = wts_map1(2:,:)

          ncstat = nf90_put_var(nc_file_id, nc_rmpmatrix_id, wts1)
          if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &             'error writing weights1')) return
          ncstat = nf90_put_var(nc_file_id, nc_rmpmatrix2_id, wts2)
          if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &             'error writing weights2')) return
          deallocate(wts1,wts2)
        endif
      else
        ncstat = nf90_put_var(nc_file_id, nc_srcadd_id, 
     &                        grid2_add_map2)
        if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &           'error writing source addresses')) return

        ncstat = nf90_put_var(nc_file_id, nc_dstadd_id, 
     &                        grid1_add_map2)
        if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &           'error writing destination addresses')) return

        if (num_wts == 1) then
          ncstat = nf90_put_var(nc_file_id, nc_rmpmatrix_id, 
     &                          wts_map2(1,:))
          if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &             'error writing weights')) return
        else
          allocate(wts1(num_links_map2),wts2(num_wts-1,num_links_map2))

          wts1 = wts_map2(1,:)
          wts2 = wts_map2(2:,:)

          ncstat = nf90_put_var(nc_file_id, nc_rmpmatrix_id, wts1)
          if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &             'error writing weights1')) return
          ncstat = nf90_put_var(nc_file_id, nc_rmpmatrix2_id, wts2)
          if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &             'error writing weights2')) return
          deallocate(wts1,wts2)
        endif
      endif

      ncstat = nf90_close(nc_file_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName,
     &         'error closing file')) return

!-----------------------------------------------------------------------

      end subroutine write_remap_csm

!***********************************************************************

      subroutine sort_add(add1, add2, weights)

!-----------------------------------------------------------------------
!
!     this routine sorts address and weight arrays based on the
!     destination address with the source address as a secondary
!     sorting criterion.  the method is a standard heap sort.
!
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!
!     Input and Output arrays
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4), intent(inout), dimension(:) ::
     &        add1,       ! destination address array (num_links)
     &        add2        ! source      address array

      real (SCRIP_r8), intent(inout), dimension(:,:) ::
     &        weights     ! remapping weights (num_wts, num_links)

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4) ::
     &          num_links,          ! num of links for this mapping
     &          num_wts,            ! num of weights for this mapping
     &          add1_tmp, add2_tmp, ! temp for addresses during swap
     &          lvl, final_lvl,     ! level indexes for heap sort levels
     &          chk_lvl1, chk_lvl2, max_lvl

      real (SCRIP_r8), dimension(SIZE(weights,DIM=1)) ::
     &          wgttmp              ! temp for holding wts during swap

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

          if ((add1(chk_lvl1) >  add1(chk_lvl2)) .OR.
     &       ((add1(chk_lvl1) == add1(chk_lvl2)) .AND.
     &        (add2(chk_lvl1) >  add2(chk_lvl2)))) then
            max_lvl = chk_lvl1
          else 
            max_lvl = chk_lvl2
          endif

          !***
          !*** if the parent is greater than both daughters,
          !*** the correct level has been found
          !***

          if ((add1_tmp .GT. add1(max_lvl)) .OR.
     &       ((add1_tmp .EQ. add1(max_lvl)) .AND.
     &        (add2_tmp .GT. add2(max_lvl)))) then
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

          if ((add1(chk_lvl1) >  add1(chk_lvl2)) .OR.
     &       ((add1(chk_lvl1) == add1(chk_lvl2)) .AND.
     &        (add2(chk_lvl1) >  add2(chk_lvl2)))) then
            max_lvl = chk_lvl1
          else 
            max_lvl = chk_lvl2
          endif

          !***
          !*** if the parent is greater than both daughters,
          !*** the correct level has been found
          !***

          if ((add1_tmp >  add1(max_lvl)) .OR.
     &       ((add1_tmp == add1(max_lvl)) .AND.
     &        (add2_tmp >  add2(max_lvl)))) then
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

!-----------------------------------------------------------------------

      end subroutine sort_add

!***********************************************************************

      end module scrip_remap_write

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
