! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_init_atm_static
!==================================================================================================
 use atm_advection
 use mpas_dmpar
 use mpas_pool_routines
 use mpas_derived_types, only : MPAS_LOG_WARN, MPAS_LOG_ERR, MPAS_LOG_CRIT
 use mpas_log, only : mpas_log_write
 use init_atm_hinterp
 use init_atm_llxy
 use mpas_c_interfacing, only : mpas_f_to_c_string

 use mpas_geometry_utils, only : mpas_in_cell
 use mpas_atmphys_utilities

 use mpas_kd_tree, only : mpas_kd_type, mpas_kd_construct, mpas_kd_free

 use iso_c_binding, only : c_char, c_int, c_float, c_loc, c_ptr

 implicit none
 private
 public:: init_atm_static,           &
          init_atm_check_read_error, &
          nearest_cell,              &
          sphere_distance

 interface
    subroutine read_geogrid(fname, rarray, nx, ny, nz, isigned, endian, &
                            scalefactor, wordsize, status) bind(C)
       use iso_c_binding, only : c_char, c_int, c_float, c_ptr
       character (c_char), dimension(*), intent(in) :: fname
       type (c_ptr), value :: rarray
       integer (c_int), intent(in), value :: nx
       integer (c_int), intent(in), value :: ny
       integer (c_int), intent(in), value :: nz
       integer (c_int), intent(in), value :: isigned
       integer (c_int), intent(in), value :: endian
       real (c_float), intent(in), value :: scalefactor
       integer (c_int), intent(in), value :: wordsize
       integer (c_int), intent(inout) :: status
    end subroutine read_geogrid
 end interface

 contains

!==================================================================================================
 subroutine init_atm_static(mesh, dims, configs)
!==================================================================================================

!inout arguments:
 type (mpas_pool_type), intent(inout) :: mesh
 type (mpas_pool_type), intent(in) :: dims
 type (mpas_pool_type), intent(in) :: configs

!constants
 integer, parameter :: nBdyLayers = 7   ! The number of relaxation layers plus the number of specified layers
                                        ! This value is used in determining whether extra checks are needed
                                        ! in the remapping of terrain, land use, and soil category pixels

!local variables:
 type(proj_info):: proj

 character(len=StrKIND) :: fname
 character(kind=c_char), dimension(StrKIND+1) :: c_fname
 character(len=StrKIND), pointer :: config_geog_data_path
 character(len=StrKIND), pointer :: config_landuse_data
 character(len=StrKIND), pointer :: config_topo_data
 character(len=StrKIND), pointer :: config_vegfrac_data
 character(len=StrKIND), pointer :: config_albedo_data
 character(len=StrKIND), pointer :: config_maxsnowalbedo_data
 character(len=StrKIND+1) :: geog_data_path      ! same as config_geog_data_path, but guaranteed to have a trailing slash
 character(len=StrKIND+1) :: geog_sub_path       ! subdirectory names in config_geog_data_path, with trailing slash

 integer:: ismax_lu

 integer(c_int):: nx,ny,nz
 integer(c_int):: endian,isigned,istatus,wordsize
 integer:: i,j,k
 integer :: ii, jj
 integer:: iCell,iEdge,iVtx,iPoint,iTileStart,iTileEnd,jTileStart,jTileEnd
 integer,dimension(5) :: interp_list
 integer,dimension(:),allocatable  :: nhs
 integer,dimension(:,:),allocatable:: ncat
      
 real(kind=c_float):: scalefactor
 real(kind=c_float),dimension(:,:,:),pointer,contiguous :: rarray
 type(c_ptr) :: rarray_ptr

 real(kind=RKIND):: start_lat
 real(kind=RKIND):: start_lon

 integer, pointer :: supersample_fac

 real(kind=RKIND):: lat,lon,x,y,z
 real(kind=RKIND):: lat_pt,lon_pt
 real(kind=RKIND),dimension(:,:),allocatable  :: soiltemp_1deg
 real(kind=RKIND),dimension(:,:),allocatable  :: maxsnowalb
 real(kind=RKIND),dimension(:,:,:),allocatable:: vegfra

 integer, pointer :: isice_lu, iswater_lu
 integer, pointer :: nCells, nEdges, nVertices, maxEdges
 logical, pointer :: on_a_sphere
 real (kind=RKIND), pointer :: sphere_radius
 
 integer, dimension(:), pointer :: nEdgesOnCell
 integer, dimension(:,:), pointer :: cellsOnCell
 integer, dimension(:,:), pointer :: verticesOnCell
 real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell
 real (kind=RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex
 real (kind=RKIND), dimension(:), pointer :: xEdge, yEdge, zEdge
 real (kind=RKIND), dimension(:), pointer :: dvEdge, dcEdge
 real (kind=RKIND), dimension(:), pointer :: areaCell, areaTriangle
 real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex
 real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
 real (kind=RKIND), dimension(:), pointer :: latVertex, lonVertex
 real (kind=RKIND), dimension(:), pointer :: latEdge, lonEdge
 real (kind=RKIND), dimension(:), pointer :: fEdge, fVertex

 real (kind=RKIND), dimension(:), pointer :: ter
 real (kind=RKIND), dimension(:), pointer :: soiltemp
 real (kind=RKIND), dimension(:), pointer :: snoalb
 real (kind=RKIND), dimension(:), pointer :: shdmin, shdmax
 real (kind=RKIND), dimension(:,:), pointer :: greenfrac
 real (kind=RKIND), dimension(:,:), pointer :: albedo12m
 real (kind=RKIND) :: msgval, fillval
 integer, dimension(:), pointer :: lu_index
 integer, dimension(:), pointer :: soilcat_top
 integer, dimension(:), pointer :: landmask
 integer, dimension(:), pointer :: bdyMaskCell
 character(len=StrKIND), pointer :: mminlu

 real (kind=RKIND) :: xPixel, yPixel, zPixel

 type (mpas_kd_type), dimension(:), pointer :: kd_points
 type (mpas_kd_type), pointer :: tree

!--------------------------------------------------------------------------------------------------


 call mpas_log_write('')
 call mpas_log_write('--- enter subroutine init_atm_static:')

 call mpas_pool_get_config(configs, 'config_geog_data_path', config_geog_data_path)
 call mpas_pool_get_config(configs, 'config_landuse_data', config_landuse_data)
 call mpas_pool_get_config(configs, 'config_topo_data', config_topo_data)
 call mpas_pool_get_config(configs, 'config_vegfrac_data', config_vegfrac_data)
 call mpas_pool_get_config(configs, 'config_albedo_data', config_albedo_data)
 call mpas_pool_get_config(configs, 'config_maxsnowalbedo_data', config_maxsnowalbedo_data)
 call mpas_pool_get_config(configs, 'config_supersample_factor', supersample_fac)

 write(geog_data_path, '(a)') config_geog_data_path
 i = len_trim(geog_data_path)
 if (geog_data_path(i:i) /= '/') then
    geog_data_path(i+1:i+1) = '/'
 end if

!
! Scale all distances and areas from a unit sphere to one with radius sphere_radius
!


 call mpas_pool_get_array(mesh, 'xCell', xCell)
 call mpas_pool_get_array(mesh, 'yCell', yCell)
 call mpas_pool_get_array(mesh, 'zCell', zCell)
 call mpas_pool_get_array(mesh, 'xVertex', xVertex)
 call mpas_pool_get_array(mesh, 'yVertex', yVertex)
 call mpas_pool_get_array(mesh, 'zVertex', zVertex)
 call mpas_pool_get_array(mesh, 'xEdge', xEdge)
 call mpas_pool_get_array(mesh, 'yEdge', yEdge)
 call mpas_pool_get_array(mesh, 'zEdge', zEdge)
 call mpas_pool_get_array(mesh, 'dcEdge', dcEdge)
 call mpas_pool_get_array(mesh, 'dvEdge', dvEdge)
 call mpas_pool_get_array(mesh, 'areaCell', areaCell)
 call mpas_pool_get_array(mesh, 'areaTriangle', areaTriangle)
 call mpas_pool_get_array(mesh, 'kiteAreasOnVertex', kiteAreasOnVertex)
 call mpas_pool_get_array(mesh, 'latCell', latCell)
 call mpas_pool_get_array(mesh, 'lonCell', lonCell)
 call mpas_pool_get_array(mesh, 'latEdge', latEdge)
 call mpas_pool_get_array(mesh, 'lonEdge', lonEdge)
 call mpas_pool_get_array(mesh, 'latVertex', latVertex)
 call mpas_pool_get_array(mesh, 'lonVertex', lonVertex)
 call mpas_pool_get_array(mesh, 'fEdge', fEdge)
 call mpas_pool_get_array(mesh, 'fVertex', fVertex)
 call mpas_pool_get_array(mesh, 'bdyMaskCell', bdyMaskCell)
 
 call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
 call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
 call mpas_pool_get_array(mesh, 'verticesOnCell', verticesOnCell)

 call mpas_pool_get_array(mesh, 'ter', ter)
 call mpas_pool_get_array(mesh, 'lu_index', lu_index)
 call mpas_pool_get_array(mesh, 'mminlu', mminlu)
 call mpas_pool_get_array(mesh, 'isice_lu', isice_lu)
 call mpas_pool_get_array(mesh, 'iswater_lu', iswater_lu)
 call mpas_pool_get_array(mesh, 'soilcat_top', soilcat_top)
 call mpas_pool_get_array(mesh, 'landmask', landmask)
 call mpas_pool_get_array(mesh, 'soiltemp', soiltemp)
 call mpas_pool_get_array(mesh, 'snoalb', snoalb)
 call mpas_pool_get_array(mesh, 'greenfrac', greenfrac)
 call mpas_pool_get_array(mesh, 'albedo12m', albedo12m)
 call mpas_pool_get_array(mesh, 'shdmin', shdmin)
 call mpas_pool_get_array(mesh, 'shdmax', shdmax)

 call mpas_pool_get_config(mesh, 'on_a_sphere', on_a_sphere)
 call mpas_pool_get_config(mesh, 'sphere_radius', sphere_radius)

 call mpas_pool_get_dimension(dims, 'nCells', nCells)
 call mpas_pool_get_dimension(dims, 'nEdges', nEdges)
 call mpas_pool_get_dimension(dims, 'nVertices', nVertices)
 call mpas_pool_get_dimension(dims, 'maxEdges', maxEdges)

 xCell = xCell * sphere_radius
 yCell = yCell * sphere_radius
 zCell = zCell * sphere_radius
 xVertex = xVertex * sphere_radius
 yVertex = yVertex * sphere_radius
 zVertex = zVertex * sphere_radius
 xEdge = xEdge * sphere_radius
 yEdge = yEdge * sphere_radius
 zEdge = zEdge * sphere_radius
 dvEdge = dvEdge * sphere_radius
 dcEdge = dcEdge * sphere_radius
 areaCell = areaCell * sphere_radius**2.0
 areaTriangle = areaTriangle * sphere_radius**2.0
 kiteAreasOnVertex = kiteAreasOnVertex * sphere_radius**2.0

!
! Initialize the KD-Tree
!
 allocate(kd_points(nCells))
 do i = 1, nCells
     allocate(kd_points(i) % point(3))
     kd_points(i) % point = (/xCell(i), yCell(i), zCell(i)/)
     kd_points(i) % id = i ! Cell ID
 enddo
 tree => null()
 tree => mpas_kd_construct(kd_points, 3)
 if (.not. associated(tree)) then
     call mpas_log_write('Error creating the KD-Tree for static interpolation', messageType=MPAS_LOG_CRIT)
 endif


!
! Initialize Coriolis parameter field on edges and vertices
!
 do iEdge=1,nEdges
    fEdge(iEdge)  = 2.0 * omega * sin(latEdge(iEdge))
 end do
 do iVtx=1,nVertices
    fVertex(iVtx) = 2.0 * omega * sin(latVertex(iVtx))
 end do


!
! Compute weights used in advection and deformation calculation
!
 call atm_initialize_advection_rk(mesh, nCells, nEdges, maxEdges, on_a_sphere, sphere_radius) 
 call atm_initialize_deformation_weights(mesh, nCells, on_a_sphere, sphere_radius) 


!
! Set land use and soil category parameters for water and ice
!
 surface_input_select0: select case(trim(config_landuse_data))
    case('USGS')
       isice_lu = 24
       iswater_lu = 16
       ismax_lu = 24
       write(mminlu,'(a)') 'USGS'
    case('MODIFIED_IGBP_MODIS_NOAH')
       isice_lu = 15
       iswater_lu = 17
       ismax_lu = 20
       write(mminlu,'(a)') 'MODIFIED_IGBP_MODIS_NOAH'
    case default
         call mpas_log_write('*****************************************************************', messageType=MPAS_LOG_ERR)
         call mpas_log_write('Invalid land use dataset '''//trim(config_landuse_data) &
                                       //''' selected for config_landuse_data',                   messageType=MPAS_LOG_ERR)
         call mpas_log_write('   Possible options are: ''USGS'', ''MODIFIED_IGBP_MODIS_NOAH''',   messageType=MPAS_LOG_ERR)
         call mpas_log_write('*****************************************************************', messageType=MPAS_LOG_ERR)
         call mpas_log_write('Please correct the namelist.', messageType=MPAS_LOG_CRIT)
 end select surface_input_select0


!
! Interpolate HGT
!
 select case(trim(config_topo_data))
    case('GTOPO30')
       call mpas_log_write('Using GTOPO30 terrain dataset')
       geog_sub_path = 'topo_30s/'
    case('GMTED2010')
       call mpas_log_write('Using GMTED2010 terrain dataset')
       geog_sub_path = 'topo_gmted2010_30s/'
    case('default')
       call mpas_log_write('*****************************************************************', messageType=MPAS_LOG_ERR)
       call mpas_log_write('Invalid topography dataset '''//trim(config_topo_data) &
                                     //''' selected for config_topo_data',                      messageType=MPAS_LOG_ERR)
       call mpas_log_write('   Possible options are: ''GTOPO30'', ''GMTED2010''',               messageType=MPAS_LOG_ERR)
       call mpas_log_write('*****************************************************************', messageType=MPAS_LOG_ERR)
       call mpas_log_write('Please correct the namelist.', messageType=MPAS_LOG_CRIT)
 end select

 ter(:) = 0.0_RKIND
 call mpas_log_write('--- end interpolate TER')


!
! Interpolate LU_INDEX
!
 surface_input_select1: select case(trim(config_landuse_data))
    case('USGS')
       call mpas_log_write('Using 24-class USGS 30-arc-second land cover dataset')
       geog_sub_path = 'landuse_30s/'
    case('MODIFIED_IGBP_MODIS_NOAH')
       call mpas_log_write('Using 20-class MODIS 30-arc-second land cover dataset')
       geog_sub_path = 'modis_landuse_20class_30s/'
    case default
       call mpas_log_write('*****************************************************************', messageType=MPAS_LOG_ERR)
       call mpas_log_write('Invalid land use dataset '''//trim(config_landuse_data) &
                                     //''' selected for config_landuse_data',                   messageType=MPAS_LOG_ERR)
       call mpas_log_write('   Possible options are: ''USGS'', ''MODIFIED_IGBP_MODIS_NOAH''',   messageType=MPAS_LOG_ERR)
       call mpas_log_write('*****************************************************************', messageType=MPAS_LOG_ERR)
       call mpas_log_write('Please correct the namelist.', messageType=MPAS_LOG_CRIT)
 end select surface_input_select1
 nx = 1200
 ny = 1200
 nz = 1
 isigned  = 1
 endian   = 0
 wordsize = 1
 scalefactor = 1.0
 allocate(rarray(nx,ny,nz))
 allocate(ncat(ismax_lu,nCells))
 ncat(:,:) = 0
 lu_index(:) = 0.0

 rarray_ptr = c_loc(rarray)

 do jTileStart = 1,20401,ny
    jTileEnd = jTileStart + ny - 1

    do iTileStart = 1,42001,nx
       iTileEnd = iTileStart + nx - 1
       write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
             trim(geog_sub_path),iTileStart,'-',iTileEnd,'.',jTileStart,'-',jTileEnd
       call mpas_log_write(trim(fname))
       call mpas_f_to_c_string(fname, c_fname)

       call read_geogrid(c_fname,rarray_ptr,nx,ny,nz,isigned,endian, &
                         scalefactor,wordsize,istatus)
       call init_atm_check_read_error(istatus, fname)

       iPoint = 1
       do j=1,ny
       do i=1,nx
!
! The MODIS dataset appears to have zeros at the South Pole, possibly other places, too
!
if (rarray(i,j,1) == 0) cycle

          lat_pt = -89.99583  + (jTileStart + j - 2) * 0.0083333333
          lon_pt = -179.99583 + (iTileStart + i - 2) * 0.0083333333
          lat_pt = lat_pt * PI / 180.0
          lon_pt = lon_pt * PI / 180.0

          iPoint = nearest_cell(lat_pt,lon_pt,iPoint,nCells,maxEdges, &
                                nEdgesOnCell,cellsOnCell, &
                                latCell,lonCell)

          !
          ! For all but the outermost boundary cells, we can safely assume that the nearest
          ! model grid cell contains the pixel (else, a different cell would be nearest)
          !
          if (bdyMaskCell(iPoint) < nBdyLayers) then
             ncat(int(rarray(i,j,1)),iPoint) = ncat(int(rarray(i,j,1)),iPoint) + 1

          ! For outermost boundary cells, additional work is needed to verify that the pixel
          ! actually lies within the nearest cell
          else
             zPixel = sphere_radius * sin(lat_pt)                  ! Model cell coordinates assume a "full" sphere radius
             xPixel = sphere_radius * cos(lon_pt) * cos(lat_pt)    ! at this point, so we need to ues the same radius
             yPixel = sphere_radius * sin(lon_pt) * cos(lat_pt)    ! for source pixel coordinates

             if (mpas_in_cell(xPixel, yPixel, zPixel, xCell(iPoint), yCell(iPoint), zCell(iPoint), &
                         nEdgesOnCell(iPoint), verticesOnCell(:,iPoint), xVertex, yVertex, zVertex)) then

                ncat(int(rarray(i,j,1)),iPoint) = ncat(int(rarray(i,j,1)),iPoint) + 1

             end if
          end if
       end do
       end do

    end do
 end do

 do iCell = 1,nCells
    lu_index(iCell) = 1
    do i = 2,ismax_lu
       if(ncat(i,iCell) > ncat(lu_index(iCell),iCell)) then
          lu_index(iCell) = i
       end if
    end do
 end do
 deallocate(rarray)
 deallocate(ncat)
 call mpas_log_write('--- end interpolate LU_INDEX')


!
! Interpolate SOILCAT_TOP
!
 nx = 1200
 ny = 1200
 nz = 1
 isigned     = 1
 endian      = 0
 wordsize    = 1
 scalefactor = 1.0
 allocate(rarray(nx,ny,nz))
 allocate(ncat(16,nCells))
 ncat(:,:) = 0
 soilcat_top(:) = 0.0

 rarray_ptr = c_loc(rarray)

 do jTileStart = 1,20401,ny
    jTileEnd = jTileStart + ny - 1

    do iTileStart = 1,42001,nx
       iTileEnd = iTileStart + nx - 1
       write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
             'soiltype_top_30s/',iTileStart,'-',iTileEnd,'.',jTileStart,'-',jTileEnd
       call mpas_log_write(trim(fname))
       call mpas_f_to_c_string(fname, c_fname)

       call read_geogrid(c_fname,rarray_ptr,nx,ny,nz,isigned,endian, &
                         scalefactor,wordsize,istatus)
       call init_atm_check_read_error(istatus, fname)

       iPoint = 1
       do j=1,ny
       do i=1,nx
          lat_pt = -89.99583  + (jTileStart + j - 2) * 0.0083333333
          lon_pt = -179.99583 + (iTileStart + i - 2) * 0.0083333333
          lat_pt = lat_pt * PI / 180.0
          lon_pt = lon_pt * PI / 180.0

          iPoint = nearest_cell(lat_pt,lon_pt,iPoint,nCells,maxEdges, &
                                nEdgesOnCell,cellsOnCell, &
                                latCell,lonCell)

          !
          ! For all but the outermost boundary cells, we can safely assume that the nearest
          ! model grid cell contains the pixel (else, a different cell would be nearest)
          !
          if (bdyMaskCell(iPoint) < nBdyLayers) then
             ncat(int(rarray(i,j,1)),iPoint) = ncat(int(rarray(i,j,1)),iPoint) + 1

          ! For outermost boundary cells, additional work is needed to verify that the pixel
          ! actually lies within the nearest cell
          else
             zPixel = sphere_radius * sin(lat_pt)                  ! Model cell coordinates assume a "full" sphere radius
             xPixel = sphere_radius * cos(lon_pt) * cos(lat_pt)    ! at this point, so we need to ues the same radius
             yPixel = sphere_radius * sin(lon_pt) * cos(lat_pt)    ! for source pixel coordinates

             if (mpas_in_cell(xPixel, yPixel, zPixel, xCell(iPoint), yCell(iPoint), zCell(iPoint), &
                         nEdgesOnCell(iPoint), verticesOnCell(:,iPoint), xVertex, yVertex, zVertex)) then

                ncat(int(rarray(i,j,1)),iPoint) = ncat(int(rarray(i,j,1)),iPoint) + 1

             end if
          end if
       end do
       end do

    end do
 end do

 do iCell = 1,nCells
    soilcat_top(iCell) = 1
    do i = 2,16
       if(ncat(i,iCell) > ncat(soilcat_top(iCell),iCell)) then
          soilcat_top(iCell) = i
       end if
    end do
 end do
 deallocate(rarray)
 deallocate(ncat)
 call mpas_log_write('--- end interpolate SOILCAT_TOP')


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! KLUDGE TO FIX SOIL TYPE OVER ANTARCTICA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 where (lu_index == isice_lu) soilcat_top = 16

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CORRECT INCONSISTENT SOIL AND LAND USE DATA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 do iCell = 1,nCells
    if (lu_index(iCell) == iswater_lu .or. &
        soilcat_top(iCell) == 14) then
        if (lu_index(iCell) /= iswater_lu) then
            call mpas_log_write('Turning lu_index into water at $i', intArgs=(/iCell/))
            lu_index(iCell) = iswater_lu
        end if
        if (soilcat_top(iCell) /= 14) then
            call mpas_log_write('Turning soilcat_top into water at $i', intArgs=(/iCell/))
            soilcat_top(iCell) = 14
        end if
    end if
 end do


!
! Derive LANDMASK
!
 landmask(:) = 0
 do iCell=1, nCells
    if (lu_index(iCell) /= iswater_lu) landmask(iCell) = 1
 end do
 call mpas_log_write('--- end interpolate LANDMASK')


!
! Interpolate SOILTEMP:
!
 nx = 186
 ny = 186
 nz = 1
 isigned  = 0
 endian   = 0
 wordsize = 2
 scalefactor = 0.01
 allocate(rarray(nx,ny,nz))
 allocate(soiltemp_1deg(-2:363,-2:183))
 soiltemp(:) = 0.0

 rarray_ptr = c_loc(rarray)

 call map_set(PROJ_LATLON, proj,  &
              latinc = 1.0_RKIND, &
              loninc = 1.0_RKIND, &
              knowni = 1.0_RKIND, &
              knownj = 1.0_RKIND, &
              lat1 = -89.5_RKIND, &
              lon1 = -179.5_RKIND)

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
       'soiltemp_1deg/',1,'-',180,'.',1,'-',180
 call mpas_log_write(trim(fname))
 call mpas_f_to_c_string(fname, c_fname)

 call read_geogrid(c_fname,rarray_ptr,nx,ny,nz,isigned, endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus, fname)
 soiltemp_1deg(-2:180,-2:183) = rarray(1:183,1:186,1)

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
            'soiltemp_1deg/',181,'-',360,'.',1,'-',180
 call mpas_log_write(trim(fname))
 call mpas_f_to_c_string(fname, c_fname)

 call read_geogrid(c_fname,rarray_ptr,nx,ny,nz,isigned,endian, &
                        scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 soiltemp_1deg(181:363,-2:183) = rarray(4:186,1:186,1)

 interp_list(1) = FOUR_POINT
 interp_list(2) = W_AVERAGE4
 interp_list(3) = W_AVERAGE16
 interp_list(4) = SEARCH
 interp_list(5) = 0

 do iCell = 1,nCells
  
    if(landmask(iCell) == 1) then
       lat = latCell(iCell) * DEG_PER_RAD
       lon = lonCell(iCell) * DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       if(x < 0.5) then
          lon = lon + 360.0
          call latlon_to_ij(proj, lat, lon, x, y)
       else if (x >= 360.5) then
          lon = lon - 360.0
          call latlon_to_ij(proj, lat, lon, x, y)
       end if
       if (y < 1.0) y = 1.0
       if (y > 179.0) y = 179.0
       soiltemp(iCell) = interp_sequence(x,y,1,soiltemp_1deg,-2,363,-2,183, &
                                           1,1,0.0_RKIND,interp_list,1)
    else
       soiltemp(iCell) = 0.0
    end if

 end do
 deallocate(rarray)
 deallocate(soiltemp_1deg)
 call mpas_log_write('--- end interpolate SOILTEMP')


!
! Interpolate SNOALB
!
 if (trim(config_maxsnowalbedo_data) == 'MODIS') then

    call mpas_log_write('Using MODIS 0.05-deg data for maximum snow albedo')
    if (supersample_fac > 1) then
       call mpas_log_write('   Dataset will be supersampled by a factor of $i', intArgs=(/supersample_fac/))
    end if

    nx = 1206
    ny = 1206
    nz = 1
    isigned  = 1
    endian   = 0
    wordsize = 2
    scalefactor = 0.01
    msgval = real(-999.0,kind=R4KIND)*real(0.01,kind=R4KIND)
    fillval = 0.0
    allocate(rarray(nx,ny,nz))
    allocate(nhs(nCells))
    nhs(:) = 0
    snoalb(:) = 0.0

    rarray_ptr = c_loc(rarray)

    start_lat = 90.0 - 0.05 * 0.5 / supersample_fac
    start_lon = -180.0 + 0.05 * 0.5 / supersample_fac
    geog_sub_path = 'maxsnowalb_modis/'

    do jTileStart = 1,02401,ny-6
       jTileEnd = jTileStart + ny - 1 - 6

       do iTileStart=1,06001,nx-6
          iTileEnd = iTileStart + nx - 1 - 6
          write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)//trim(geog_sub_path), &
                       iTileStart,'-',iTileEnd,'.',jTileStart,'-',jTileEnd
          call mpas_log_write(trim(fname))
          call mpas_f_to_c_string(fname, c_fname)

          call read_geogrid(c_fname,rarray_ptr,nx,ny,nz,isigned,endian, &
                            scalefactor,wordsize,istatus)
          call init_atm_check_read_error(istatus, fname)

          iPoint = 1
          do j=supersample_fac * 3 + 1, supersample_fac * (ny-3)
          do i=supersample_fac * 3 + 1, supersample_fac * (nx-3)
             ii = (i - 1) / supersample_fac + 1
             jj = (j - 1) / supersample_fac + 1

             lat_pt = start_lat - (supersample_fac*(jTileStart-1) + j - (supersample_fac*3+1)) * 0.05 / supersample_fac
             lon_pt = start_lon + (supersample_fac*(iTileStart-1) + i - (supersample_fac*3+1)) * 0.05 / supersample_fac
             lat_pt = lat_pt * PI / 180.0
             lon_pt = lon_pt * PI / 180.0

             iPoint = nearest_cell(lat_pt,lon_pt,iPoint,nCells,maxEdges, &
                                   nEdgesOnCell,cellsOnCell, &
                                   latCell,lonCell)
             if (rarray(ii,jj,1) /= msgval) then

                !
                ! This field only matters for land cells, and for all but the outermost boundary cells,
                ! we can safely assume that the nearest model grid cell contains the pixel (else, a different
                ! cell would be nearest)
                !
                if (landmask(iPoint) == 1 .and. bdyMaskCell(iPoint) < nBdyLayers) then
                   snoalb(iPoint) = snoalb(iPoint) + rarray(ii,jj,1)
                   nhs(iPoint) = nhs(iPoint) + 1

                ! For outermost land cells, additional work is needed to verify that the pixel
                ! actually lies within the nearest cell
                else if (landmask(iPoint) == 1) then
                   zPixel = sphere_radius * sin(lat_pt)                  ! Model cell coordinates assume a "full" sphere radius
                   xPixel = sphere_radius * cos(lon_pt) * cos(lat_pt)    ! at this point, so we need to ues the same radius
                   yPixel = sphere_radius * sin(lon_pt) * cos(lat_pt)    ! for source pixel coordinates

                   if (mpas_in_cell(xPixel, yPixel, zPixel, xCell(iPoint), yCell(iPoint), zCell(iPoint), &
                               nEdgesOnCell(iPoint), verticesOnCell(:,iPoint), xVertex, yVertex, zVertex)) then
                      snoalb(iPoint) = snoalb(iPoint) + rarray(ii,jj,1)
                      nhs(iPoint) = nhs(iPoint) + 1
                   end if
                end if
             end if
          end do
          end do

       end do
    end do

    do iCell = 1,nCells
       !
       ! Mismatches in land mask can lead to MPAS land points with no maximum snow albedo.
       ! Ideally, we would perform a search for nearby valid albedos, but for now using
       ! the fill value will at least allow the model to run. In general, the number of cells
       ! to be treated in this way tends to be a very small fraction of the total number of cells.
       !
       if (nhs(iCell) == 0) then
          snoalb(iCell) = fillval
       else
          snoalb(iCell) = snoalb(iCell) / real(nhs(iCell))
       end if
       snoalb(iCell) = 0.01_RKIND * snoalb(iCell)        ! Convert from percent to fraction
    end do
    deallocate(rarray)
    deallocate(nhs)

 else if (trim(config_maxsnowalbedo_data) == 'NCEP') then

    call mpas_log_write('Using NCEP 1.0-deg data for maximum snow albedo')

    nx = 186
    ny = 186
    nz = 1
    isigned     = 0
    endian      = 0
    wordsize    = 1
    scalefactor = 1.0
    allocate(rarray(nx,ny,nz))
    allocate(maxsnowalb(-2:363,-2:183))
    snoalb(:) = 0.0

    rarray_ptr = c_loc(rarray)

    call map_set(PROJ_LATLON, proj,  &
                 latinc = 1.0_RKIND, &
                 loninc = 1.0_RKIND, &
                 knowni = 1.0_RKIND, &
                 knownj = 1.0_RKIND, &
                 lat1 = -89.5_RKIND, &
                 lon1 = -179.5_RKIND)

    write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
          'maxsnowalb/',1,'-',180,'.',1,'-',180
    call mpas_log_write(trim(fname))
    call mpas_f_to_c_string(fname, c_fname)

    call read_geogrid(c_fname,rarray_ptr,nx,ny,nz,isigned,endian, &
                      scalefactor,wordsize,istatus)
    call init_atm_check_read_error(istatus,fname)
    maxsnowalb(-2:180,-2:183) = rarray(1:183,1:186,1)

    write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
          'maxsnowalb/',181,'-',360,'.',1,'-',180
    call mpas_log_write(trim(fname))
    call mpas_f_to_c_string(fname, c_fname)

    call read_geogrid(c_fname,rarray_ptr,nx,ny,nz,isigned,endian, &
                      scalefactor,wordsize,istatus)
    call init_atm_check_read_error(istatus, fname)
    maxsnowalb(181:363,-2:183) = rarray(4:186,1:186,1)

    interp_list(1) = FOUR_POINT
    interp_list(2) = W_AVERAGE4
    interp_list(3) = W_AVERAGE16
    interp_list(4) = SEARCH
    interp_list(5) = 0

    do iCell = 1,nCells

       if(landmask(iCell) == 1) then
          lat = latCell(iCell) * DEG_PER_RAD
          lon = lonCell(iCell) * DEG_PER_RAD
          call latlon_to_ij(proj, lat, lon, x, y)
          if(x < 0.5) then
             lon = lon + 360.0
             call latlon_to_ij(proj, lat, lon, x, y)
          else if (x >= 360.5) then
             lon = lon - 360.0
             call latlon_to_ij(proj, lat, lon, x, y)
          end if
          if (y < 1.0) y = 1.0
          if (y > 179.0) y = 179.0
          snoalb(iCell) = interp_sequence(x,y,1,maxsnowalb,-2,363,-2,183, &
                                            1,1,0.0_RKIND,interp_list,1)
       else
          snoalb(iCell) = 0.0
       end if

    end do
    snoalb(:) = snoalb(:) / 100.0
    deallocate(rarray)
    deallocate(maxsnowalb)

 else

    call mpas_log_write('*****************************************************************', messageType=MPAS_LOG_ERR)
    call mpas_log_write('Invalid maximum snow albedo dataset '''//trim(config_maxsnowalbedo_data) &
                                  //''' selected for config_maxsnowalbedo_data',             messageType=MPAS_LOG_ERR)
    call mpas_log_write('   Possible options are: ''MODIS'', ''NCEP''',                      messageType=MPAS_LOG_ERR)
    call mpas_log_write('*****************************************************************', messageType=MPAS_LOG_ERR)
    call mpas_log_write('Please correct the namelist.', messageType=MPAS_LOG_CRIT)

 end if

 call mpas_log_write('--- end interpolate SNOALB')


!
! Interpolate GREENFRAC
!
 if (trim(config_vegfrac_data) == 'MODIS') then

    call mpas_log_write('Using MODIS FPAR 30-arc-second data for climatological monthly vegetation fraction')

    nx = 1200
    ny = 1200
    nz = 12
    isigned  = 0
    endian   = 0
    wordsize = 1
    scalefactor = 1.0
    msgval = 200.0
    fillval = 0.0
    allocate(rarray(nx,ny,nz))
    allocate(nhs(nCells))
    nhs(:) = 0
    greenfrac(:,:) = 0.0

    rarray_ptr = c_loc(rarray)

    start_lat = -89.99583
    start_lon = -179.99583
    geog_sub_path = 'greenfrac_fpar_modis/'

    do jTileStart = 1,20401,ny
       jTileEnd = jTileStart + ny - 1

       do iTileStart=1,42001,nx
          iTileEnd = iTileStart + nx - 1
          write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)//trim(geog_sub_path), &
                       iTileStart,'-',iTileEnd,'.',jTileStart,'-',jTileEnd
          call mpas_log_write(trim(fname))
          call mpas_f_to_c_string(fname, c_fname)

          call read_geogrid(c_fname,rarray_ptr,nx,ny,nz,isigned,endian, &
                            scalefactor,wordsize,istatus)
          call init_atm_check_read_error(istatus, fname)

          iPoint = 1
          do j=1,ny
          do i=1,nx
             lat_pt = start_lat + (jTileStart + j - 2) * 0.0083333333
             lon_pt = start_lon + (iTileStart + i - 2) * 0.0083333333
             lat_pt = lat_pt * PI / 180.0
             lon_pt = lon_pt * PI / 180.0

             iPoint = nearest_cell(lat_pt,lon_pt,iPoint,nCells,maxEdges, &
                                   nEdgesOnCell,cellsOnCell, &
                                   latCell,lonCell)

             !
             ! This field only matters for land cells, and for all but the outermost boundary cells,
             ! we can safely assume that the nearest model grid cell contains the pixel (else, a different
             ! cell would be nearest)
             !
             if (landmask(iPoint) == 1 .and. bdyMaskCell(iPoint) < nBdyLayers) then
                do k=1,nz
                   if (rarray(i,j,k) == msgval) then
                      rarray(i,j,k) = fillval
                   end if
                   greenfrac(k,iPoint) = greenfrac(k,iPoint) + rarray(i,j,k)
                end do
                nhs(iPoint) = nhs(iPoint) + 1

             ! For outermost land cells, additional work is needed to verify that the pixel
             ! actually lies within the nearest cell
             else if (landmask(iPoint) == 1) then
                zPixel = sphere_radius * sin(lat_pt)                  ! Model cell coordinates assume a "full" sphere radius
                xPixel = sphere_radius * cos(lon_pt) * cos(lat_pt)    ! at this point, so we need to ues the same radius
                yPixel = sphere_radius * sin(lon_pt) * cos(lat_pt)    ! for source pixel coordinates

                if (mpas_in_cell(xPixel, yPixel, zPixel, xCell(iPoint), yCell(iPoint), zCell(iPoint), &
                            nEdgesOnCell(iPoint), verticesOnCell(:,iPoint), xVertex, yVertex, zVertex)) then
                   do k=1,nz
                      if (rarray(i,j,k) == msgval) then
                         rarray(i,j,k) = fillval
                      end if
                      greenfrac(k,iPoint) = greenfrac(k,iPoint) + rarray(i,j,k)
                   end do
                   nhs(iPoint) = nhs(iPoint) + 1
                end if
             end if
          end do
          end do

       end do
    end do

    do iCell = 1,nCells
       ! For land points that have no overlap with valid data, and for water points,
       ! just use the fill value...
       if (nhs(iCell) == 0) then
          greenfrac(:,iCell) = fillval
       else
          greenfrac(:,iCell) = greenfrac(:,iCell) / real(nhs(iCell))
       end if
       shdmin(iCell) = minval(greenfrac(:,iCell))
       shdmax(iCell) = maxval(greenfrac(:,iCell))
    end do
    deallocate(rarray)
    deallocate(nhs)

 else if (trim(config_vegfrac_data) == 'NCEP') then

    call mpas_log_write('Using NCEP 0.144-deg data for climatological monthly vegetation fraction')

    nx = 1256
    ny = 1256
    nz = 12
    isigned     = 0
    endian      = 0
    wordsize    = 1
    scalefactor = 1.0
    allocate(rarray(nx,ny,nz))
    allocate(vegfra(-2:2503,-2:1253,12))
    greenfrac(:,:) = 0.0

    rarray_ptr = c_loc(rarray)

    call map_set(PROJ_LATLON, proj,    &
                 latinc = 0.144_RKIND, &
                 loninc = 0.144_RKIND, &
                 knowni = 1.0_RKIND,   &
                 knownj = 1.0_RKIND,   &
                 lat1 = -89.928_RKIND, &
                 lon1 = -179.928_RKIND)

    write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
          'greenfrac/',1,'-',1250,'.',1,'-',1250
    call mpas_log_write(trim(fname))
    call mpas_f_to_c_string(fname, c_fname)

    call read_geogrid(c_fname,rarray_ptr,nx,ny,nz,isigned,endian, &
                      scalefactor,wordsize,istatus)
    call init_atm_check_read_error(istatus,fname)
    vegfra(-2:1250,-2:1253,1:12) = rarray(1:1253,1:1256,1:12)

    write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
          'greenfrac/',1251,'-',2500,'.',1,'-',1250
    call mpas_log_write(trim(fname))
    call mpas_f_to_c_string(fname, c_fname)

    call read_geogrid(c_fname,rarray_ptr,nx,ny,nz,isigned,endian, &
                      scalefactor,wordsize,istatus)
    call init_atm_check_read_error(istatus,fname)
    vegfra(1251:2503,-2:1253,1:12) = rarray(4:1256,1:1256,1:12)

    do iCell = 1,nCells

       if (landmask(iCell) == 1) then
          lat = latCell(iCell) * DEG_PER_RAD
          lon = lonCell(iCell) * DEG_PER_RAD
          call latlon_to_ij(proj, lat, lon, x, y)
          if(x < 0.5) then
             lon = lon + 360.0
             call latlon_to_ij(proj, lat, lon, x, y)
          else if(x >= 2500.5) then
             lon = lon - 360.0
             call latlon_to_ij(proj, lat, lon, x, y)
          end if
          if (y < 1.0) y = 1.0
          if (y > 1249.0) y = 1249.0
          do k = 1,12
             greenfrac(k,iCell) = interp_sequence(x,y,k,vegfra,-2,2503,-2,1253, &
                                                    1,12,-1.e30_RKIND,interp_list,1)
          end do
       else
          greenfrac(:,iCell) = 0.0
       end if
       shdmin(iCell) = minval(greenfrac(:,iCell))
       shdmax(iCell) = maxval(greenfrac(:,iCell))

    end do
    deallocate(rarray)
    deallocate(vegfra)

 else

    call mpas_log_write('*****************************************************************', messageType=MPAS_LOG_ERR)
    call mpas_log_write('Invalid monthly vegetation fraction dataset '''//trim(config_vegfrac_data) &
                                  //''' selected for config_vegfrac_data',                   messageType=MPAS_LOG_ERR)
    call mpas_log_write('   Possible options are: ''MODIS'', ''NCEP''',                      messageType=MPAS_LOG_ERR)
    call mpas_log_write('*****************************************************************', messageType=MPAS_LOG_ERR)
    call mpas_log_write('Please correct the namelist.', messageType=MPAS_LOG_CRIT)

 end if

 call mpas_log_write('--- end interpolate GREENFRAC')


!
! Interpolate ALBEDO12M
!
 if (trim(config_albedo_data) == 'MODIS') then

    call mpas_log_write('Using MODIS 0.05-deg data for climatological monthly albedo')
    if (supersample_fac > 1) then
       call mpas_log_write('   Dataset will be supersampled by a factor of $i', intArgs=(/supersample_fac/))
    end if

    nx = 1206
    ny = 1206
    nz = 12
    isigned  = 1
    endian   = 0
    wordsize = 2
    scalefactor = 0.01
    msgval = real(-999.0,kind=R4KIND)*real(0.01,kind=R4KIND)
    fillval = 8.0
    allocate(rarray(nx,ny,nz))
    allocate(nhs(nCells))
    nhs(:) = 0
    albedo12m(:,:) = 0.0

    rarray_ptr = c_loc(rarray)

    start_lat = 90.0 - 0.05 * 0.5 / supersample_fac
    start_lon = -180.0 + 0.05 * 0.5 / supersample_fac
    geog_sub_path = 'albedo_modis/'

    do jTileStart = 1,02401,ny-6
       jTileEnd = jTileStart + ny - 1 - 6

       do iTileStart=1,06001,nx-6
          iTileEnd = iTileStart + nx - 1 - 6
          write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)//trim(geog_sub_path), &
                       iTileStart,'-',iTileEnd,'.',jTileStart,'-',jTileEnd
          call mpas_log_write(trim(fname))
          call mpas_f_to_c_string(fname, c_fname)

          call read_geogrid(c_fname,rarray_ptr,nx,ny,nz,isigned,endian, &
                            scalefactor,wordsize,istatus)
          call init_atm_check_read_error(istatus, fname)

          iPoint = 1
          do j=supersample_fac * 3 + 1, supersample_fac * (ny-3)
          do i=supersample_fac * 3 + 1, supersample_fac * (nx-3)
             ii = (i - 1) / supersample_fac + 1
             jj = (j - 1) / supersample_fac + 1

             lat_pt = start_lat - (supersample_fac*(jTileStart-1) + j - (supersample_fac*3+1)) * 0.05 / supersample_fac
             lon_pt = start_lon + (supersample_fac*(iTileStart-1) + i - (supersample_fac*3+1)) * 0.05 / supersample_fac
             lat_pt = lat_pt * PI / 180.0
             lon_pt = lon_pt * PI / 180.0

             iPoint = nearest_cell(lat_pt,lon_pt,iPoint,nCells,maxEdges, &
                                   nEdgesOnCell,cellsOnCell, &
                                   latCell,lonCell)

             !
             ! This field only matters for land cells, and for all but the outermost boundary cells,
             ! we can safely assume that the nearest model grid cell contains the pixel (else, a different
             ! cell would be nearest)
             !
             if (landmask(iPoint) == 1 .and. bdyMaskCell(iPoint) < nBdyLayers) then
                do k=1,nz
                   if (rarray(ii,jj,k) == msgval) then
                      rarray(ii,jj,k) = fillval
                   end if
                   albedo12m(k,iPoint) = albedo12m(k,iPoint) + rarray(ii,jj,k)
                end do
                nhs(iPoint) = nhs(iPoint) + 1

             ! For outermost land cells, additional work is needed to verify that the pixel
             ! actually lies within the nearest cell
             else if (landmask(iPoint) == 1) then
                zPixel = sphere_radius * sin(lat_pt)                  ! Model cell coordinates assume a "full" sphere radius
                xPixel = sphere_radius * cos(lon_pt) * cos(lat_pt)    ! at this point, so we need to ues the same radius
                yPixel = sphere_radius * sin(lon_pt) * cos(lat_pt)    ! for source pixel coordinates

                if (mpas_in_cell(xPixel, yPixel, zPixel, xCell(iPoint), yCell(iPoint), zCell(iPoint), &
                            nEdgesOnCell(iPoint), verticesOnCell(:,iPoint), xVertex, yVertex, zVertex)) then
                   do k=1,nz
                      if (rarray(ii,jj,k) == msgval) then
                         rarray(ii,jj,k) = fillval
                      end if
                      albedo12m(k,iPoint) = albedo12m(k,iPoint) + rarray(ii,jj,k)
                   end do
                   nhs(iPoint) = nhs(iPoint) + 1
                end if
             end if
          end do
          end do

       end do
    end do

    do iCell = 1,nCells
       ! For land points that have no overlap with valid data, and for water points,
       ! just use the fill value...
       if (nhs(iCell) == 0) then
          albedo12m(:,iCell) = fillval
       else
          albedo12m(:,iCell) = albedo12m(:,iCell) / real(nhs(iCell))
       end if
       if (lu_index(iCell) == isice_lu) then
          albedo12m(:,iCell) = 70.0
       end if
    end do
    deallocate(rarray)
    deallocate(nhs)

 else if (trim(config_albedo_data) == 'NCEP') then

    call mpas_log_write('Using NCEP 0.144-deg data for climatological monthly albedo')

    nx = 1256
    ny = 1256
    nz = 12
    isigned     = 0
    endian      = 0
    wordsize    = 1
    scalefactor = 1.0
    allocate(rarray(nx,ny,nz))
    allocate(vegfra(-2:2503,-2:1253,12))
    albedo12m(:,:) = 0.0

    rarray_ptr = c_loc(rarray)

    call map_set(PROJ_LATLON, proj,    &
                 latinc = 0.144_RKIND, &
                 loninc = 0.144_RKIND, &
                 knowni = 1.0_RKIND,   &
                 knownj = 1.0_RKIND,   &
                 lat1 = -89.928_RKIND, &
                 lon1 = -179.928_RKIND)

    write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
          'albedo_ncep/',1,'-',1250,'.',1,'-',1250
    call mpas_log_write(trim(fname))
    call mpas_f_to_c_string(fname, c_fname)

    call read_geogrid(c_fname,rarray_ptr,nx,ny,nz,isigned,endian, &
                      scalefactor, wordsize, istatus)
    call init_atm_check_read_error(istatus,fname)
    vegfra(-2:1250,-2:1253,1:12) = rarray(1:1253,1:1256,1:12)

    write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
          'albedo_ncep/',1251,'-',2500,'.',1,'-',1250
    call mpas_log_write(trim(fname))
    call mpas_f_to_c_string(fname, c_fname)

    call read_geogrid(c_fname,rarray_ptr,nx,ny,nz,isigned,endian, &
                      scalefactor,wordsize,istatus)
    call init_atm_check_read_error(istatus,fname)
    vegfra(1251:2503,-2:1253,1:12) = rarray(4:1256,1:1256,1:12)

    do iCell = 1,nCells

       if (landmask(iCell) == 1) then
          lat = latCell(iCell) * DEG_PER_RAD
          lon = lonCell(iCell) * DEG_PER_RAD
          call latlon_to_ij(proj, lat, lon, x, y)
          if(x < 0.5) then
             lon = lon + 360.0
             call latlon_to_ij(proj, lat, lon, x, y)
          else if(x >= 2500.5) then
             lon = lon - 360.0
             call latlon_to_ij(proj, lat, lon, x, y)
          end if
          if (y < 1.0) y = 1.0
          if (y > 1249.0) y = 1249.0
          do k = 1,12
             albedo12m(k,iCell) = interp_sequence(x,y,k,vegfra,-2,2503,-2,1253, &
                                                    1,12,0.0_RKIND,interp_list,1)
          end do
       else
          albedo12m(:,iCell) = 8.0
       end if
    end do
    deallocate(rarray)
    deallocate(vegfra)

 else

    call mpas_log_write('*****************************************************************', messageType=MPAS_LOG_ERR)
    call mpas_log_write('Invalid monthly albedo dataset '''//trim(config_albedo_data) &
                                  //''' selected for config_albedo_data',                    messageType=MPAS_LOG_ERR)
    call mpas_log_write('   Possible options are: ''MODIS'', ''NCEP''',                      messageType=MPAS_LOG_ERR)
    call mpas_log_write('*****************************************************************', messageType=MPAS_LOG_ERR)
    call mpas_log_write('Please correct the namelist.', messageType=MPAS_LOG_CRIT)

 end if

 call mpas_log_write('--- end interpolate ALBEDO12M')

!
! Deallocate and free the KD Tree
!
 call mpas_kd_free(tree)
 deallocate(kd_points)


 end subroutine init_atm_static

!==================================================================================================
 subroutine init_atm_check_read_error(istatus, fname)
!==================================================================================================
 implicit none

 integer, intent(in) :: istatus
 character (len=*), intent(in) :: fname

 if (istatus /= 0) then
     call mpas_log_write('Could not read file '//trim(fname), messageType=MPAS_LOG_CRIT)
 end if

 end subroutine init_atm_check_read_error

!==================================================================================================
 integer function nearest_cell(target_lat, target_lon, start_cell, nCells, maxEdges, &
                               nEdgesOnCell, cellsOnCell, latCell, lonCell)
!==================================================================================================
 implicit none

 real (kind=RKIND), intent(in) :: target_lat, target_lon
 integer, intent(in) :: start_cell
 integer, intent(in) :: nCells, maxEdges
 integer, dimension(nCells), intent(in) :: nEdgesOnCell
 integer, dimension(maxEdges,nCells), intent(in) :: cellsOnCell
 real (kind=RKIND), dimension(nCells), intent(in) :: latCell, lonCell

 integer :: i
 integer :: iCell
 integer :: current_cell
 real (kind=RKIND) :: current_distance, d
 real (kind=RKIND) :: nearest_distance

 nearest_cell = start_cell
 current_cell = -1

 do while (nearest_cell /= current_cell)
    current_cell = nearest_cell
    current_distance = sphere_distance(latCell(current_cell), lonCell(current_cell), target_lat, &
                                       target_lon, 1.0_RKIND)
    nearest_cell = current_cell
    nearest_distance = current_distance
    do i = 1, nEdgesOnCell(current_cell)
       iCell = cellsOnCell(i,current_cell)
       if (iCell <= nCells) then
          d = sphere_distance(latCell(iCell), lonCell(iCell), target_lat, target_lon, 1.0_RKIND)
          if (d < nearest_distance) then
             nearest_cell = iCell
             nearest_distance = d
          end if
       end if
    end do
 end do

 end function nearest_cell

!==================================================================================================
 real (kind=RKIND) function sphere_distance(lat1, lon1, lat2, lon2, radius)

!Compute the great-circle distance between (lat1, lon1) and (lat2, lon2) on a
!sphere with given radius.
!==================================================================================================
 implicit none

 real (kind=RKIND), intent(in) :: lat1, lon1, lat2, lon2, radius
 real (kind=RKIND) :: arg1

 arg1 = sqrt( sin(0.5*(lat2-lat1))**2 +  &
              cos(lat1)*cos(lat2)*sin(0.5*(lon2-lon1))**2 )
 sphere_distance = 2.*radius*asin(arg1)

 end function sphere_distance


!==================================================================================================
 end module mpas_init_atm_static
!==================================================================================================
