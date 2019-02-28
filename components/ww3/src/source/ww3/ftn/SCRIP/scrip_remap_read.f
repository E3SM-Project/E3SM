!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     This routine reads remapping information from files written
!     by remap_setup.  If remapping in both directions are required,
!     two input files must be specified.
!
!-----------------------------------------------------------------------
!
!     CVS:$Id: remap_read.f,v 1.6 2000/04/19 21:56:26 pwjones Exp $
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

      module scrip_remap_read 

!-----------------------------------------------------------------------
!
!     contains routines for reading a remap file
!
!-----------------------------------------------------------------------

      use SCRIP_KindsMod   ! defines common data types
      use SCRIP_ErrorMod   ! SCRIP error handler
      use SCRIP_NetcdfMod  ! netCDF error handler
      use netcdf           ! netCDF library

      use SCRIP_constants  ! defines common scalar constants
      use scrip_grids      ! module containing grid information
      use scrip_remap_vars ! module containing remap information

      implicit none

!-----------------------------------------------------------------------
!
!     module variables
!
!-----------------------------------------------------------------------


      character(SCRIP_charLength), private :: 
     &   map_method       ! character string for map_type
     &,  normalize_opt    ! character string for normalization option
     &,  convention       ! character string for output convention

!-----------------------------------------------------------------------
!
!     various netCDF ids for files variables
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4), private ::
     &   ncstat                ! error flag for netCDF calls
     &,  nc_file_id            ! id for netCDF file
     &,  nc_srcgrdsize_id      ! id for source grid size
     &,  nc_dstgrdsize_id      ! id for destination grid size
     &,  nc_srcgrdcorn_id      ! id for number of source grid corners
     &,  nc_dstgrdcorn_id      ! id for number of dest grid corners
     &,  nc_srcgrdrank_id      ! id for source grid rank
     &,  nc_dstgrdrank_id      ! id for dest grid rank
     &,  nc_numlinks_id        ! id for number of links in mapping
     &,  nc_numwgts_id         ! id for number of weights for mapping
     &,  nc_srcgrddims_id      ! id for source grid dimensions
     &,  nc_dstgrddims_id      ! id for dest grid dimensions,
     &,  nc_srcgrdcntrlat_id   ! id for source grid center latitude
     &,  nc_dstgrdcntrlat_id   ! id for dest grid center latitude
     &,  nc_srcgrdcntrlon_id   ! id for source grid center longitude
     &,  nc_dstgrdcntrlon_id   ! id for dest grid center longitude
     &,  nc_srcgrdimask_id     ! id for source grid mask
     &,  nc_dstgrdimask_id     ! id for dest grid mask
     &,  nc_srcgrdcrnrlat_id   ! id for latitude of source grid corners
     &,  nc_srcgrdcrnrlon_id   ! id for longitude of source grid corners
     &,  nc_dstgrdcrnrlat_id   ! id for latitude of dest grid corners
     &,  nc_dstgrdcrnrlon_id   ! id for longitude of dest grid corners
     &,  nc_srcgrdarea_id      ! id for area of source grid cells
     &,  nc_dstgrdarea_id      ! id for area of dest grid cells
     &,  nc_srcgrdfrac_id      ! id for area fraction on source grid
     &,  nc_dstgrdfrac_id      ! id for area fraction on dest grid
     &,  nc_srcgrdadd_id       ! id for map source address
     &,  nc_dstgrdadd_id       ! id for map dest address
     &,  nc_rmpmatrix_id       ! id for remapping matrix

!***********************************************************************

      contains

!***********************************************************************

      subroutine read_remap_ww3(map_name, interp_file, errorCode)

!-----------------------------------------------------------------------
!
!     this driver routine reads some global attributes and then
!     calls a specific read routine based on file conventions
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     input variables
!
!-----------------------------------------------------------------------

      character(SCRIP_CharLength), intent(in) ::
     &            interp_file  ! filename for remap data

!-----------------------------------------------------------------------
!
!     output variables
!
!-----------------------------------------------------------------------

      character(SCRIP_CharLength), intent(out) ::
     &            map_name     ! name for mapping

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (SCRIP_i4) ::  
     &  errorCode            ! error code for SCRIP routine

      character (14), parameter :: 
     &  rtnName = 'read_remap_ww3'

!-----------------------------------------------------------------------
!
!     open file and read some global information
!
!-----------------------------------------------------------------------

      ncstat = nf90_open(interp_file, NF90_NOCLOBBER, nc_file_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, 'error
     &   opening interpolation file')) return

      !***
      !*** map name
      !***
      map_name = ' '
      ncstat = nf90_get_att(nc_file_id, NF90_GLOBAL, 'title',
     &                         map_name)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, 'error
     &   reading map name')) return


      !***
      !*** normalization option
      !***
      normalize_opt = ' '
      ncstat = nf90_get_att(nc_file_id, NF90_GLOBAL, 'normalization',
     &                         normalize_opt)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, 'error
     &   reading interpolation option')) return

      select case(normalize_opt)
      case ('none')
        norm_opt = norm_opt_none
      case ('fracarea')
        norm_opt = norm_opt_frcarea
      case ('destarea')
        norm_opt = norm_opt_dstarea
      case default
        print *,'normalize_opt = ',normalize_opt
        stop 'Invalid normalization option'
      end select

      !***
      !*** map method
      !***
      map_method = ' '
      ncstat = nf90_get_att (nc_file_id, NF90_GLOBAL, 'map_method',
     &                          map_method)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, 'error
     &   reading map method')) return

      select case(map_method)
      case('Conservative remapping')
        map_type = map_type_conserv
      case('Bilinear remapping')
        map_type = map_type_bilinear
      case('Distance weighted avg of nearest neighbors')
        map_type = map_type_distwgt
      case('Bicubic remapping')
        map_type = map_type_bicubic
      case default
        print *,'map_type = ',map_method
        stop 'Invalid Map Type'
      end select

      !***
      !*** file convention
      !***
      convention = ' '
      ncstat = nf90_get_att (nc_file_id, NF90_GLOBAL, 'conventions',
     &                          convention)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, 'error
     &   reading file convention')) return

!-----------------------------------------------------------------------
!
!     call appropriate read routine based on output convention
!
!-----------------------------------------------------------------------

      select case(convention)
      case ('SCRIP')
        call read_remap_scrip_ww3
      case ('NCAR-CSM')
        print *,'convention = NCAR-CSM not supported for WW3'
        stop 'unsupported file convention'
!       call read_remap_csm
      case default
        print *,'convention = ',convention
        stop 'unknown output file convention'
      end select

!-----------------------------------------------------------------------

      end subroutine read_remap_ww3

!***********************************************************************

      subroutine read_remap_scrip_ww3

!-----------------------------------------------------------------------
!
!     the routine reads a netCDF file to extract remapping info
!     in SCRIP format
!
!     Only read variables needed by WW3
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      character (SCRIP_charLength) ::
     &  grid1_name           ! grid name for source grid
     &, grid2_name           ! grid name for dest   grid

      integer (SCRIP_i4) ::  
     &  n                    ! dummy index
     &, errorCode            ! error code for SCRIP routine

      integer (SCRIP_i4), dimension(:), allocatable ::
     &  grid1_mask_int,      ! integer masks to determine
     &  grid2_mask_int       ! cells that participate in map

      character (20), parameter :: 
     &  rtnName = 'read_remap_scrip_ww3'

!-----------------------------------------------------------------------
!
!     read some additional global attributes
!
!-----------------------------------------------------------------------

      !***
      !*** source and destination grid names
      !***

      grid1_name = ' '
      grid2_name = ' '
      errorCode = SCRIP_Success
      ncstat = nf90_get_att (nc_file_id, NF90_GLOBAL, 'source_grid',
     &                          grid1_name)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, 'error 
     &   reading source grid name')) return

      ncstat = nf90_get_att (nc_file_id, NF90_GLOBAL, 'dest_grid',
     &                          grid2_name)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, 'error
     &   reading destination grid name')) return

! Let's not include routine write statements w/out check on improc, l_master, etc.
!     print *,' '
!     print *,'Remapping between:',trim(grid1_name)
!     print *,'and ',trim(grid2_name)
!     print *,' '

!-----------------------------------------------------------------------
!
!     read dimension information
!
!-----------------------------------------------------------------------

      ncstat = nf90_inq_dimid(nc_file_id, 'dst_grid_size', 
     &                      nc_dstgrdsize_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, 'error
     &   reading destination grid id')) return
!     ncstat = nf90_inq_dimlen(nc_file_id, nc_dstgrdsize_id, grid2_size)
      ncstat = nf90_inquire_dimension(nc_file_id, nc_dstgrdsize_id, 
     &                                len = grid2_size)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, 'error
     &   reading destination grid size')) return

      ncstat = nf90_inq_dimid(nc_file_id, 'num_links', 
     &                      nc_numlinks_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, 'error
     &   reading number of links id')) return
!     ncstat = nf90_inq_dimlen(nc_file_id, nc_numlinks_id, 
!    &                       num_links_map1)
      ncstat = nf90_inquire_dimension(nc_file_id, nc_numlinks_id, 
     &                                len = num_links_map1)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, 'error
     &   reading number of links')) return
      max_links_map1 = num_links_map1

      ncstat = nf90_inq_dimid(nc_file_id, 'num_wgts', 
     &                      nc_numwgts_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, 'error
     &   reading number of weights id')) return
!     ncstat = nf90_inq_dimlen(nc_file_id, nc_numwgts_id, num_wts)
      ncstat = nf90_inquire_dimension(nc_file_id, nc_numwgts_id, 
     &                                len = num_wts)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, 'error
     &   reading number of weights')) return

!-----------------------------------------------------------------------
!
!     allocate arrays
!
!-----------------------------------------------------------------------

! grid2_frac is allocated in grid_init
      if(allocated(grid2_frac))deallocate(grid2_frac)
      allocate( grid2_frac      (grid2_size) )

! grid1_add_map1, grid2_add_map1, wts_map1 are allocated in init_remap_vars
      if(allocated(grid1_add_map1))deallocate(grid1_add_map1)
      if(allocated(grid2_add_map1))deallocate(grid2_add_map1)
      if(allocated(wts_map1))deallocate(wts_map1)
      allocate( grid1_add_map1(num_links_map1),
     &          grid2_add_map1(num_links_map1),
     &          wts_map1(num_wts,num_links_map1) )

!-----------------------------------------------------------------------
!
!     get variable ids
!
!-----------------------------------------------------------------------

      ncstat = nf90_inq_varid(nc_file_id, 'dst_grid_frac', 
     &                                   nc_dstgrdfrac_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, 'error
     &   reading destination grid fraction id')) return

      ncstat = nf90_inq_varid(nc_file_id, 'src_address', 
     &                                   nc_srcgrdadd_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, 'error
     &   reading source grid address id')) return

      ncstat = nf90_inq_varid(nc_file_id, 'dst_address', 
     &                                   nc_dstgrdadd_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, 'error
     &   reading destination grid address id')) return

      ncstat = nf90_inq_varid(nc_file_id, 'remap_matrix', 
     &                                   nc_rmpmatrix_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, 'error
     &   reading remap matrix id')) return

!-----------------------------------------------------------------------
!
!     read all variables
!
!-----------------------------------------------------------------------

      ncstat = nf90_get_var(nc_file_id, nc_dstgrdfrac_id, 
     &                                       grid2_frac)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, 'error
     &   reading destination grid fraction')) return

      ncstat = nf90_get_var(nc_file_id, nc_srcgrdadd_id, 
     &                        grid1_add_map1)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, 'error
     &   reading source grid address')) return

      ncstat = nf90_get_var(nc_file_id, nc_dstgrdadd_id, 
     &                        grid2_add_map1)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, 'error
     &   reading destination grid address')) return

      ncstat = nf90_get_var(nc_file_id, nc_rmpmatrix_id, 
     &                                       wts_map1)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, 'error
     &   reading remap weights')) return


      ncstat = nf90_close(nc_file_id)
      if (SCRIP_NetcdfErrorCheck(ncstat, errorCode, rtnName, 'error
     &   closing netCDF file')) return

!-----------------------------------------------------------------------

      end subroutine read_remap_scrip_ww3

!***********************************************************************

      end module scrip_remap_read

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
