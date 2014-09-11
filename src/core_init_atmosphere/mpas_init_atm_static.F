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
! use mpas_configure
 use mpas_dmpar
 use init_atm_hinterp
 use init_atm_llxy

 use mpas_atmphys_utilities

 implicit none
 private
 public:: init_atm_static,           &
          init_atm_static_orogwd,    &
          init_atm_check_read_error, &
          nearest_cell,              &
          sphere_distance

 contains

!==================================================================================================
 subroutine init_atm_static(mesh, dims, configs)
!==================================================================================================

!inout arguments:
 type (mpas_pool_type), intent(inout) :: mesh
 type (mpas_pool_type), intent(in) :: dims
 type (mpas_pool_type), intent(in) :: configs

!local variables:
 type(proj_info):: proj

 character(len=StrKIND) :: fname
 character(len=StrKIND), pointer :: config_geog_data_path
 character(len=StrKIND+1) :: geog_data_path      ! same as config_geog_data_path, but guaranteed to have a trailing slash

 integer:: nx,ny,nz
 integer:: endian,isigned,istatus,wordsize
 integer:: i,j,k
 integer:: iCell,iEdge,iVtx,iPoint,iTileStart,iTileEnd,jTileStart,jTileEnd
 integer,dimension(5) :: interp_list
 integer,dimension(:),allocatable  :: nhs
 integer,dimension(:,:),allocatable:: ncat
      
 real(kind=4):: scalefactor
 real(kind=4),dimension(:,:,:),allocatable:: rarray

 real(kind=RKIND):: lat,lon,x,y
 real(kind=RKIND):: lat_pt,lon_pt
 real(kind=RKIND),dimension(:,:),allocatable  :: soiltemp_1deg
 real(kind=RKIND),dimension(:,:),allocatable  :: maxsnowalb
 real(kind=RKIND),dimension(:,:,:),allocatable:: vegfra

 integer, pointer :: nCells, nEdges, nVertices, maxEdges
 logical, pointer :: on_a_sphere
 real (kind=RKIND), pointer :: sphere_radius
 
 integer, dimension(:), pointer :: nEdgesOnCell
 integer, dimension(:,:), pointer :: cellsOnCell
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
 integer, dimension(:), pointer :: lu_index
 integer, dimension(:), pointer :: soilcat_top
 integer, dimension(:), pointer :: soilcat_bot
 integer, dimension(:), pointer :: landmask

!--------------------------------------------------------------------------------------------------


 write(0,*)
 write(0,*) '--- enter subroutine init_atm_static:'

 call mpas_pool_get_config(configs, 'config_geog_data_path', config_geog_data_path)

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
 
 call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
 call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)

 call mpas_pool_get_array(mesh, 'ter', ter)
 call mpas_pool_get_array(mesh, 'lu_index', lu_index)
 call mpas_pool_get_array(mesh, 'soilcat_top', soilcat_top)
 call mpas_pool_get_array(mesh, 'soilcat_bot', soilcat_bot)
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
! Interpolate HGT
!
!nx = 126
!ny = 126
 nx = 1206
 ny = 1206
 nz = 1
 isigned  = 1
 endian   = 0
 wordsize = 2
 scalefactor = 1.0
 allocate(rarray(nx,ny,nz))
 allocate(nhs(nCells))
 nhs(:) = 0
 ter(:) = 0.0

 do jTileStart = 1,20401,ny-6
    jTileEnd = jTileStart + ny - 1 - 6

    do iTileStart=1,42001,nx-6
       iTileEnd = iTileStart + nx - 1 - 6
       write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
             'topo_30s/',iTileStart,'-',iTileEnd,'.',jTileStart,'-',jTileEnd
       write(0,*) trim(fname)

       call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                         scalefactor,wordsize,istatus)
       call init_atm_check_read_error(istatus, fname)

       iPoint = 1
       do j=4,ny-3
       do i=4,nx-3
          lat_pt = -89.99583  + (jTileStart + j - 5) * 0.0083333333
          lon_pt = -179.99583 + (iTileStart + i - 5) * 0.0083333333
          lat_pt = lat_pt * PI / 180.0
          lon_pt = lon_pt * PI / 180.0

          iPoint = nearest_cell(lat_pt,lon_pt,iPoint,nCells,maxEdges, &
                                nEdgesOnCell,cellsOnCell, &
                                latCell,lonCell)
          ter(iPoint) = ter(iPoint) + rarray(i,j,1)
          nhs(iPoint) = nhs(iPoint) + 1
       end do
       end do

    end do
 end do

 do iCell = 1,nCells
    ter(iCell) = ter(iCell) / real(nhs(iCell))
 end do
 deallocate(rarray)
 deallocate(nhs)
 write(0,*) '--- end interpolate TER'


!
! Interpolate LU_INDEX
!
 nx = 1200
 ny = 1200
 nz = 1
 isigned  = 1
 endian   = 0
 wordsize = 1
 scalefactor = 1.0
 allocate(rarray(nx,ny,nz))
 allocate(ncat(24,nCells))
 ncat(:,:) = 0
 lu_index(:) = 0.0

 do jTileStart = 1,20401,ny
    jTileEnd = jTileStart + ny - 1

    do iTileStart = 1,42001,nx
       iTileEnd = iTileStart + nx - 1
       write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
             'landuse_30s/',iTileStart,'-',iTileEnd,'.',jTileStart,'-',jTileEnd
       write(0,*) trim(fname)

       call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
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
          ncat(int(rarray(i,j,1)),iPoint) = ncat(int(rarray(i,j,1)),iPoint) + 1
       end do
       end do

    end do
 end do

 do iCell = 1,nCells
    lu_index(iCell) = 1
    do i = 2,24
       if(ncat(i,iCell) > ncat(lu_index(iCell),iCell)) then
          lu_index(iCell) = i
       end if
    end do
 end do
 deallocate(rarray)
 deallocate(ncat)
 write(0,*) '--- end interpolate LU_INDEX'


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

 do jTileStart = 1,20401,ny
    jTileEnd = jTileStart + ny - 1

    do iTileStart = 1,42001,nx
       iTileEnd = iTileStart + nx - 1
       write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
             'soiltype_top_30s/',iTileStart,'-',iTileEnd,'.',jTileStart,'-',jTileEnd
       write(0,*) trim(fname)

       call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
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
          ncat(int(rarray(i,j,1)),iPoint) = ncat(int(rarray(i,j,1)),iPoint) + 1
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
 write(0,*) '--- end interpolate SOILCAT_TOP'


!
! Interpolate SOILCAT_BOT
!
 nx = 1200
 ny = 1200
 nz = 1
 isigned  = 1
 endian   = 0
 wordsize = 1
 scalefactor = 1.0
 allocate(rarray(nx,ny,nz))
 allocate(ncat(16,nCells))
 ncat(:,:) = 0
 soilcat_bot(:) = 0.0

 do jTileStart = 1,20401,ny
    jTileEnd = jTileStart + ny - 1

    do iTileStart = 1,42001,nx
       iTileEnd = iTileStart + nx - 1
       write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
             'soiltype_bot_30s/',iTileStart,'-',iTileEnd,'.',jTileStart,'-',jTileEnd
       write(0,*) trim(fname)

       call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
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
          ncat(int(rarray(i,j,1)),iPoint) = ncat(int(rarray(i,j,1)),iPoint) + 1
       end do
       end do

    end do
 end do

 do iCell =1,nCells
    soilcat_bot(iCell) = 1
    do i = 2,16
       if(ncat(i,iCell) > ncat(soilcat_bot(iCell),iCell)) then
          soilcat_bot(iCell) = i
       end if
    end do
 end do
 deallocate(rarray)
 deallocate(ncat)
 write(0,*) '--- end interpolate SOILCAT_BOT'


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! KLUDGE TO FIX SOIL TYPE OVER ANTARCTICA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 where (lu_index == 24) soilcat_top = 16
 where (lu_index == 24) soilcat_bot = 16

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CORRECT INCONSISTENT SOIL AND LAND USE DATA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 do iCell = 1,nCells
    if (lu_index(iCell) == 16 .or. &
        soilcat_top(iCell) == 14 .or. &
        soilcat_bot(iCell) == 14) then
        if (lu_index(iCell) /= 16) then
            write(0,*) 'Turning lu_index into water at ', iCell
            lu_index(iCell) = 16
        end if
        if (soilcat_top(iCell) /= 14) then
            write(0,*) 'Turning soilcat_top into water at ', iCell
            soilcat_top(iCell) = 14
        end if
        if (soilcat_bot(iCell) /= 14) then
            write(0,*) 'Turning soilcat_bot into water at ', iCell
            soilcat_bot(iCell) = 14
        end if
    end if
 end do


!
! Derive LANDMASK
!
 landmask(:) = 0
 do iCell=1, nCells
    if (lu_index(iCell) /= 16) landmask(iCell) = 1
 end do
 write(0,*) '--- end interpolate LANDMASK'


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

 call map_set(PROJ_LATLON, proj,  &
              latinc = 1.0_RKIND, &
              loninc = 1.0_RKIND, &
              knowni = 1.0_RKIND, &
              knownj = 1.0_RKIND, &
              lat1 = -89.5_RKIND, &
              lon1 = -179.5_RKIND)

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
       'soiltemp_1deg/',1,'-',180,'.',1,'-',180
 write(0,*) trim(fname)

 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned, endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus, fname)
 soiltemp_1deg(-2:180,-2:183) = rarray(1:183,1:186,1)

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
            'soiltemp_1deg/',181,'-',360,'.',1,'-',180
 write(0,*) trim(fname)

 call read_geogrid(fname, len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
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
 write(0,*) '--- end interpolate SOILTEMP'


!
! Interpolate SNOALB
!
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

 call map_set(PROJ_LATLON, proj,  &
              latinc = 1.0_RKIND, &
              loninc = 1.0_RKIND, &
              knowni = 1.0_RKIND, &
              knownj = 1.0_RKIND, &
              lat1 = -89.5_RKIND, &
              lon1 = -179.5_RKIND)

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
       'maxsnowalb/',1,'-',180,'.',1,'-',180
 write(0,*) trim(fname)

 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, & 
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 maxsnowalb(-2:180,-2:183) = rarray(1:183,1:186,1)

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
       'maxsnowalb/',181,'-',360,'.',1,'-',180
 write(0,*) trim(fname)

 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
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
 write(0,*) '--- end interpolate SNOALB'


!
! Interpolate GREENFRAC
!
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

 call map_set(PROJ_LATLON, proj,    &
              latinc = 0.144_RKIND, &
              loninc = 0.144_RKIND, &
              knowni = 1.0_RKIND,   &
              knownj = 1.0_RKIND,   &
              lat1 = -89.928_RKIND, &
              lon1 = -179.928_RKIND)

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
       'greenfrac/',1,'-',1250,'.',1,'-',1250
 write(0,*) trim(fname)

 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, & 
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 vegfra(-2:1250,-2:1253,1:12) = rarray(1:1253,1:1256,1:12)

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
       'greenfrac/',1251,'-',2500,'.',1,'-',1250
 write(0,*) trim(fname)

 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
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
 write(0,*) '--- end interpolate GREENFRAC'


!
! Interpolate ALBEDO12M
!
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

 call map_set(PROJ_LATLON, proj,    &
              latinc = 0.144_RKIND, &
              loninc = 0.144_RKIND, &
              knowni = 1.0_RKIND,   &
              knownj = 1.0_RKIND,   &
              lat1 = -89.928_RKIND, &
              lon1 = -179.928_RKIND)

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
       'albedo_ncep/',1,'-',1250,'.',1,'-',1250
 write(0,*) trim(fname)

 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor, wordsize, istatus)
 call init_atm_check_read_error(istatus,fname)
 vegfra(-2:1250,-2:1253,1:12) = rarray(1:1253,1:1256,1:12)

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
       'albedo_ncep/',1251,'-',2500,'.',1,'-',1250
 write(0,*) trim(fname)

 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, & 
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
 write(0,*) '--- end interpolate ALBEDO12M'


 end subroutine init_atm_static

!==================================================================================================
 subroutine init_atm_static_orogwd(mesh, dims, configs)
!==================================================================================================

!inout arguments:
 type (mpas_pool_type), intent(inout) :: mesh
 type (mpas_pool_type), intent(in) :: dims
 type (mpas_pool_type), intent(in) :: configs

!local variables:
 type(proj_info):: proj

 character(len=StrKIND) :: mess
 character(len=StrKIND) :: fname
 character(len=StrKIND) :: dir_gwdo
 character(len=StrKIND), pointer :: config_geog_data_path
 character(len=StrKIND+1) :: geog_data_path      ! same as config_geog_data_path, but guaranteed to have a trailing slash

 integer, pointer :: nCells, nEdges, maxEdges

 integer, dimension(:), pointer :: nEdgesOnCell
 integer, dimension(:,:), pointer :: cellsOnCell
 integer, dimension(:), pointer :: landmask

 integer:: nx,ny,nz
 integer:: endian,isigned,istatus,wordsize
 integer:: i,j
 integer:: iCell,iPoint,iTileStart,iTileEnd,jTileStart,jTileEnd
 integer,dimension(5) :: interp_list
 integer,dimension(:),allocatable:: nhs

 real(kind=4):: scalefactor
 real(kind=4),dimension(:,:,:),allocatable:: rarray

 real(kind=RKIND):: lat,lon,x,y
 real(kind=RKIND):: lat_pt,lon_pt
 real(kind=RKIND):: dx,dy,known_lat,known_lon,known_x,known_y
 real(kind=RKIND):: minMeshD,maxMeshD
 real(kind=RKIND):: mindcEdge,maxdcEdge
 real(kind=RKIND),dimension(:,:),allocatable:: xarray

 real(kind=RKIND), dimension(:), pointer :: latCell, lonCell
 real(kind=RKIND), dimension(:), pointer :: meshDensity
 real(kind=RKIND), dimension(:), pointer :: dcEdge
 real(kind=RKIND), dimension(:), pointer :: varsso
 real(kind=RKIND), dimension(:), pointer :: con, oa1, oa2, oa3, oa4, ol1, ol2, ol3, ol4, var2d


 call mpas_pool_get_dimension(dims, 'nCells', nCells)
 call mpas_pool_get_dimension(dims, 'nEdges', nEdges)
 call mpas_pool_get_dimension(dims, 'maxEdges', maxEdges)

 call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
 call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
 call mpas_pool_get_array(mesh, 'landmask', landmask)
 call mpas_pool_get_array(mesh, 'latCell', latCell)
 call mpas_pool_get_array(mesh, 'lonCell', lonCell)
 call mpas_pool_get_array(mesh, 'varsso', varsso)
 call mpas_pool_get_array(mesh, 'meshDensity', meshDensity)
 call mpas_pool_get_array(mesh, 'dcEdge', dcEdge)
 call mpas_pool_get_array(mesh, 'con', con)
 call mpas_pool_get_array(mesh, 'oa1', oa1)
 call mpas_pool_get_array(mesh, 'oa2', oa2)
 call mpas_pool_get_array(mesh, 'oa3', oa3)
 call mpas_pool_get_array(mesh, 'oa4', oa4)
 call mpas_pool_get_array(mesh, 'ol1', ol1)
 call mpas_pool_get_array(mesh, 'ol2', ol2)
 call mpas_pool_get_array(mesh, 'ol3', ol3)
 call mpas_pool_get_array(mesh, 'ol4', ol4)
 call mpas_pool_get_array(mesh, 'var2d', var2d)


 write(0,*)
 write(0,*) '--- enter subroutine init_atm_static_orogwd:'

 call mpas_pool_get_config(configs, 'config_geog_data_path', config_geog_data_path)

 write(geog_data_path, '(a)') config_geog_data_path
 i = len_trim(geog_data_path)
 if (geog_data_path(i:i) /= '/') then
    geog_data_path(i+1:i+1) = '/'
 end if

!
! Interpolate VARSSO:
 varsso(:) = 0.0_RKIND
 nx = 600
 ny = 600
 nz = 1
 isigned     = 0
 endian      = 0
 wordsize    = 4
 scalefactor = 1.0

 dx = 0.00833333
 dy = 0.00833333
 known_x = 1.0
 known_y = 1.0
 known_lat = -59.99583
 known_lon = -179.99583
  
 allocate(rarray(nx,ny,nz))
 allocate(nhs(nCells))
 nhs(:) = 0
 rarray(:,:,:) = 0._RKIND
 do jTileStart = 1,13801,ny
    jTileEnd = jTileStart + ny - 1

    do iTileStart = 1,42601,nx
       iTileEnd = iTileStart + nx -1
       write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)//'varsso/', &
             iTileStart,'-',iTileEnd,'.',jTileStart,'-',jTileEnd
       write(0,*) trim(fname)

       call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                         scalefactor,wordsize,istatus)
       call init_atm_check_read_error(istatus,fname)

       iPoint = 1
       do j = 1,ny
       do i = 1,nx
          lat_pt = known_lat + (jTileStart + j - 2) * dy
          lon_pt = known_lon + (iTileStart + i - 2) * dx
          lat_pt = lat_pt * PI / 180.0
          lon_pt = lon_pt * PI / 180.0

          iPoint = nearest_cell(lat_pt,lon_pt,iPoint,nCells,maxEdges, &
                                nEdgesOnCell,cellsOnCell, &
                                latCell,lonCell)
          varsso(iPoint) = varsso(iPoint) + rarray(i,j,1)
          nhs(iPoint) = nhs(iPoint) + 1
       enddo
       enddo

    enddo
 enddo

 do iCell = 1,nCells
    if(nhs(iCell) .gt. 0) &
       varsso(iCell) = varsso(iCell) / real(nhs(iCell))
 enddo
 deallocate(rarray)
 deallocate(nhs)
 write(0,*) '--- end interpolate VARSSO'

!... statistic fields needed for the parameterization of gravity wavwe drag over orography. The
!input directory depends on the mesh resolution, and the mesh must be a uniform mesh.
 minMeshD  = minval(meshDensity(1:nCells))
 maxMeshD  = maxval(meshDensity(1:nCells))
 mindcEdge = minval(dcEdge(1:nEdges))
 maxdcEdge = maxval(dcEdge(1:nEdges))

 write(0,*)
 write(0,*) 'BEGIN INTERPOLATION OF STATISTICAL FIELDS FOR GRAVITY WAVE DRAG OVER OROGRAPHY'
 write(0,*) 'min MeshD  =', minMeshD
 write(0,*) 'max MeshD  =', maxMeshD
 write(0,*) 'min dcEdge =', mindcEdge
 write(0,*) 'max dcEdge =', maxdcEdge

 dir_gwdo = '   '
 if(minMeshD == 1.0_RKIND .and. maxMeshD == 1.0_RKIND) then
    !... uniform 10242 mesh:
    if(mindcEdge .ge. 200000._RKIND .and. maxdcEdge .lt. 260000._RKIND) then
       dir_gwdo = 'orogwd_2deg'
    elseif(mindcEdge .ge. 90000._RKIND .and. maxdcEdge .lt. 150000_RKIND) then
       dir_gwdo = 'orogwd_1deg'
    elseif(mindcEdge .ge. 40000._RKIND .and. maxdcEdge .lt. 70000._RKIND) then
       dir_gwdo = 'orogwd_30m'
    else
       write(0,*)
!      write(mess,*) 'GWDO: Interpolation not available. The initialization will abort'
!      call physics_error_fatal(mess)
       write(mess,*) 'GWDO: Interpolation not available. Set config_gwdo_scheme = .false.'
       return
    endif
 else
    write(0,*)
!   write(mess,*) 'GWDO: The input mesh must be a uniform mesh. The initialization will abort'
!   call physics_error_fatal(mess)
    write(mess,*) 'GWDO: The input mesh must be a uniform mesh. Set config_gwdo_scheme = .false.'
    return
 endif
 write(0,*) 'dir_gwdo   =    ', trim(dir_gwdo)
 write(0,*)

!
! Interpolate CON:
!
 con(:) = 0.0_RKIND

 con_select: select case(dir_gwdo)
    case("orogwd_2deg")
       nx = 180
       ny =  90
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.025
       dx = 2.0
       dy = 2.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.0
       known_lon =   1.0
    case("orogwd_1deg")
       nx = 360
       ny = 180
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.025
       dx = 1.0
       dy = 1.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.5
       known_lon =   0.5
    case("orogwd_30m")
       nx = 720
       ny = 360
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.025
       dx = 0.5
       dy = 0.5
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.75
       known_lon = 0.25
    case("orogwd_10m")
       nx = 2160
       ny = 1080
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.025
       dx = 0.16666667
       dy = 0.16666667
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.916667
       known_lon = 0.0833333
    case default
 end select con_select

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') &
       trim(geog_data_path)//trim(dir_gwdo)//'/con/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 xarray(1:nx,1:ny) = rarray(1:nx,1:ny,1)

 call map_set(PROJ_LATLON, proj,  &
              latinc = dy,        &
              loninc = dx,        &
              knowni = known_x,   &
              knownj = known_y,   &
              lat1   = known_lat, &
              lon1   = known_lon)

 interp_list(1) = AVERAGE4
 interp_list(2) = AVERAGE4
 interp_list(3) = AVERAGE4
 interp_list(4) = AVERAGE4
 interp_list(5) = 0

 do iCell = 1,nCells
    if(landmask(iCell) == 1) then
       lat = latCell(iCell) * DEG_PER_RAD
       lon = lonCell(iCell) * DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       con(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                                    0.0_RKIND,interp_list,1)
    endif
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate CON'

!
! Interpolate OA1:
!
 oa1(:) = 0.0_RKIND

 oa1_select: select case(dir_gwdo)
    case("orogwd_2deg")
       nx = 180
       ny =  90
       nz = 1
       isigned     = 1
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 2.0
       dy = 2.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.0
       known_lon =   1.0
    case("orogwd_1deg")
       nx = 360
       ny = 180
       nz = 1
       isigned     = 1
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 1.0
       dy = 1.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.5
       known_lon =   0.5
    case("orogwd_30m")
       nx = 720
       ny = 360
       nz = 1
       isigned     = 1
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 0.5
       dy = 0.5
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.75
       known_lon = 0.25
    case("orogwd_10m")
       nx = 2160
       ny = 1080
       nz = 1
       isigned     = 1
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 0.16666667
       dy = 0.16666667
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.916667
       known_lon = 0.0833333
    case default
 end select oa1_select

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') &
       trim(geog_data_path)//trim(dir_gwdo)//'/oa1/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 xarray(1:nx,1:ny) = rarray(1:nx,1:ny,1)

 call map_set(PROJ_LATLON, proj,  &
              latinc = dy,        &
              loninc = dx,        &
              knowni = known_x,   &
              knownj = known_y,   &
              lat1   = known_lat, &
              lon1   = known_lon)

 interp_list(1) = AVERAGE4
 interp_list(2) = AVERAGE4
 interp_list(3) = AVERAGE4
 interp_list(4) = AVERAGE4
 interp_list(5) = 0

 do iCell = 1,nCells
    if(landmask(iCell) == 1) then
       lat = latCell(iCell) * DEG_PER_RAD
       lon = lonCell(iCell) * DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       oa1(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                                    0.0_RKIND,interp_list,1)
    endif
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate OA1'

!
! Interpolate OA2:
 oa2(:) = 0.0_RKIND

 oa2_select: select case(dir_gwdo)
    case("orogwd_2deg")
       nx = 180
       ny =  90
       nz = 1
       isigned     = 1
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 2.0
       dy = 2.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.0
       known_lon =   1.0
    case("orogwd_1deg")
       nx = 360
       ny = 180
       nz = 1
       isigned     = 1
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 1.0
       dy = 1.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.5
       known_lon =   0.5
    case("orogwd_30m")
       nx = 720
       ny = 360
       nz = 1
       isigned     = 1
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 0.5
       dy = 0.5
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.75
       known_lon = 0.25
    case("orogwd_10m")
       nx = 2160
       ny = 1080
       nz = 1
       isigned     = 1
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 0.16666667
       dy = 0.16666667
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.916667
       known_lon = 0.0833333
    case default
 end select oa2_select

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') &
       trim(geog_data_path)//trim(dir_gwdo)//'/oa2/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 xarray(1:nx,1:ny) = rarray(1:nx,1:ny,1)

 call map_set(PROJ_LATLON, proj,  &
              latinc = dy,        &
              loninc = dx,        &
              knowni = known_x,   &
              knownj = known_y,   &
              lat1   = known_lat, &
              lon1   = known_lon)

 interp_list(1) = AVERAGE4
 interp_list(2) = AVERAGE4
 interp_list(3) = AVERAGE4
 interp_list(4) = AVERAGE4
 interp_list(5) = 0

 do iCell = 1,nCells
    if(landmask(iCell) == 1) then
       lat = latCell(iCell) * DEG_PER_RAD
       lon = lonCell(iCell) * DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       oa2(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                                    0.0_RKIND,interp_list,1)
     endif
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate OA2'

!
! Interpolate OA3:
!
 oa3(:) = 0.0_RKIND

 oa3_select: select case(dir_gwdo)
    case("orogwd_2deg")
       nx = 180
       ny =  90
       nz = 1
       isigned     = 1
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 2.0
       dy = 2.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.0
       known_lon =   1.0
    case("orogwd_1deg")
       nx = 360
       ny = 180
       nz = 1
       isigned     = 1
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 1.0
       dy = 1.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.5
       known_lon =   0.5
    case("orogwd_30m")
       nx = 720
       ny = 360
       nz = 1
       isigned     = 1
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 0.5
       dy = 0.5
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.75
       known_lon = 0.25
    case("orogwd_10m")
       nx = 2160
       ny = 1080
       nz = 1
       isigned     = 1
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 0.16666667
       dy = 0.16666667
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.916667
       known_lon = 0.0833333
    case default
 end select oa3_select

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') &
       trim(geog_data_path)//trim(dir_gwdo)//'/oa3/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 xarray(1:nx,1:ny) = rarray(1:nx,1:ny,1)

 call map_set(PROJ_LATLON, proj,  &
              latinc = dy,        &
              loninc = dx,        &
              knowni = known_x,   &
              knownj = known_y,   &
              lat1   = known_lat, &
              lon1   = known_lon)

 interp_list(1) = AVERAGE4
 interp_list(2) = AVERAGE4
 interp_list(3) = AVERAGE4
 interp_list(4) = AVERAGE4
 interp_list(5) = 0

 do iCell = 1,nCells
    if(landmask(iCell) == 1) then
       lat = latCell(iCell) * DEG_PER_RAD
       lon = lonCell(iCell) * DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       oa3(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                                    0.0_RKIND,interp_list,1)
    endif
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate OA3'

!
! Interpolate OA4:
!
 oa4(:) = 0.0_RKIND

 oa4_select: select case(dir_gwdo)
    case("orogwd_2deg")
       nx = 180
       ny =  90
       nz = 1
       isigned     = 1
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 2.0
       dy = 2.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.0
       known_lon =   1.0
    case("orogwd_1deg")
       nx = 360
       ny = 180
       nz = 1
       isigned     = 1
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 1.0
       dy = 1.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.5
       known_lon =   0.5
    case("orogwd_30m")
       nx = 720
       ny = 360
       nz = 1
       isigned     = 1
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 0.5
       dy = 0.5
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.75
       known_lon = 0.25
    case("orogwd_10m")
       nx = 2160
       ny = 1080
       nz = 1
       isigned     = 1
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 0.16666667
       dy = 0.16666667
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.916667
       known_lon = 0.0833333
    case default
 end select oa4_select

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') &
       trim(geog_data_path)//trim(dir_gwdo)//'/oa4/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 xarray(1:nx,1:ny) = rarray(1:nx,1:ny,1)

 call map_set(PROJ_LATLON, proj,  &
              latinc = dy,        &
              loninc = dx,        &
              knowni = known_x,   &
              knownj = known_y,   &
              lat1   = known_lat, &
              lon1   = known_lon)

 interp_list(1) = AVERAGE4
 interp_list(2) = AVERAGE4
 interp_list(3) = AVERAGE4
 interp_list(4) = AVERAGE4
 interp_list(5) = 0

 do iCell = 1,nCells
    if(landmask(iCell) == 1) then
       lat = latCell(iCell) * DEG_PER_RAD
       lon = lonCell(iCell) * DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       oa4(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                                    0.0_RKIND,interp_list,1)
    endif
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate OA4'

!
! Interpolate OL1:
!
 ol1(:) = 0.0_RKIND

 ol1_select: select case(dir_gwdo)
    case("orogwd_2deg")
       nx = 180
       ny =  90
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 2.0
       dy = 2.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.0
       known_lon =   1.0
    case("orogwd_1deg")
       nx = 360
       ny = 180
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 1.0
       dy = 1.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.5
       known_lon =   0.5
    case("orogwd_30m")
       nx = 720
       ny = 360
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 0.5
       dy = 0.5
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.75
       known_lon = 0.25
    case("orogwd_10m")
       nx = 2160
       ny = 1080
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 0.16666667
       dy = 0.16666667
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.916667
       known_lon = 0.0833333
    case default
 end select ol1_select

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') &
       trim(geog_data_path)//trim(dir_gwdo)//'/ol1/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 xarray(1:nx,1:ny) = rarray(1:nx,1:ny,1)

 call map_set(PROJ_LATLON, proj,  &
              latinc = dy,        &
              loninc = dx,        &
              knowni = known_x,   &
              knownj = known_y,   &
              lat1   = known_lat, &
              lon1   = known_lon)

 interp_list(1) = AVERAGE4
 interp_list(2) = AVERAGE4
 interp_list(3) = AVERAGE4
 interp_list(4) = AVERAGE4
 interp_list(5) = 0

 do iCell = 1,nCells
    if(landmask(iCell) == 1) then
       lat = latCell(iCell) * DEG_PER_RAD
       lon = lonCell(iCell) * DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       ol1(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                                    0.0_RKIND,interp_list,1)
    endif
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate OL1'

!
! Interpolate OL2:
!
 ol2(:) = 0.0_RKIND

 ol2_select: select case(dir_gwdo)
    case("orogwd_2deg")
       nx = 180
       ny =  90
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 2.0
       dy = 2.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.0
       known_lon =   1.0
    case("orogwd_1deg")
       nx = 360
       ny = 180
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 1.0
       dy = 1.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.5
       known_lon =   0.5
    case("orogwd_30m")
       nx = 720
       ny = 360
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 0.5
       dy = 0.5
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.75
       known_lon = 0.25
    case("orogwd_10m")
       nx = 2160
       ny = 1080
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 0.16666667
       dy = 0.16666667
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.916667
       known_lon = 0.0833333
    case default
 end select ol2_select

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') &
       trim(geog_data_path)//trim(dir_gwdo)//'/ol2/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 xarray(1:nx,1:ny) = rarray(1:nx,1:ny,1)

 call map_set(PROJ_LATLON, proj,  &
              latinc = dy,        &
              loninc = dx,        &
              knowni = known_x,   &
              knownj = known_y,   &
              lat1   = known_lat, &
              lon1   = known_lon)

 interp_list(1) = AVERAGE4
 interp_list(2) = AVERAGE4
 interp_list(3) = AVERAGE4
 interp_list(4) = AVERAGE4
 interp_list(5) = 0

 do iCell = 1,nCells
    if(landmask(iCell) == 1) then
       lat = latCell(iCell) * DEG_PER_RAD
       lon = lonCell(iCell) * DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       ol2(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                                    0.0_RKIND,interp_list,1)
    endif
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate OL2'

!
! Interpolate OL3:
!
 ol3(:) = 0.0_RKIND

 ol3_select: select case(dir_gwdo)
    case("orogwd_2deg")
       nx = 180
       ny =  90
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 2.0
       dy = 2.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.0
       known_lon =   1.0
    case("orogwd_1deg")
       nx = 360
       ny = 180
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 1.0
       dy = 1.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.5
       known_lon =   0.5
    case("orogwd_30m")
       nx = 720
       ny = 360
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 0.5
       dy = 0.5
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.75
       known_lon = 0.25
    case("orogwd_10m")
       nx = 2160
       ny = 1080
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 0.16666667
       dy = 0.16666667
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.916667
       known_lon = 0.0833333
    case default
 end select ol3_select

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') &
       trim(geog_data_path)//trim(dir_gwdo)//'/ol3/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 xarray(1:nx,1:ny) = rarray(1:nx,1:ny,1)

 call map_set(PROJ_LATLON, proj,  &
              latinc = dy,        &
              loninc = dx,        &
              knowni = known_x,   &
              knownj = known_y,   &
              lat1   = known_lat, &
              lon1   = known_lon)

 interp_list(1) = AVERAGE4
 interp_list(2) = AVERAGE4
 interp_list(3) = AVERAGE4
 interp_list(4) = AVERAGE4
 interp_list(5) = 0

 do iCell = 1,nCells
    if(landmask(iCell) == 1) then
       lat = latCell(iCell) * DEG_PER_RAD
       lon = lonCell(iCell) * DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       ol3(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                                    0.0_RKIND,interp_list,1)
    endif
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate OL3'

!
! Interpolate OL4:
!
 ol4(:) = 0.0_RKIND

 ol4_select: select case(dir_gwdo)
    case("orogwd_2deg")
       nx = 180
       ny =  90
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 2.0
       dy = 2.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.0
       known_lon =   1.0
    case("orogwd_1deg")
       nx = 360
       ny = 180
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 1.0
       dy = 1.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.5
       known_lon =   0.5
    case("orogwd_30m")
       nx = 720
       ny = 360
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 0.5
       dy = 0.5
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.75
       known_lon = 0.25
    case("orogwd_10m")
       nx = 2160
       ny = 1080
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.0001
       dx = 0.16666667
       dy = 0.16666667
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.916667
       known_lon = 0.0833333
    case default
 end select ol4_select

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') &
       trim(geog_data_path)//trim(dir_gwdo)//'/ol4/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 xarray(1:nx,1:ny) = rarray(1:nx,1:ny,1)

 call map_set(PROJ_LATLON, proj,  &
              latinc = dy,        &
              loninc = dx,        &
              knowni = known_x,   &
              knownj = known_y,   &
              lat1   = known_lat, &
              lon1   = known_lon)

 interp_list(1) = AVERAGE4
 interp_list(2) = AVERAGE4
 interp_list(3) = AVERAGE4
 interp_list(4) = AVERAGE4
 interp_list(5) = 0

 do iCell = 1,nCells
    if(landmask(iCell) == 1) then
       lat = latCell(iCell) * DEG_PER_RAD
       lon = lonCell(iCell) * DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       ol4(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                                    0.0_RKIND,interp_list,1)
    endif
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate OL4'

!
! Interpolate VAR2D:
!
 var2d(:) = 0.0_RKIND

 var2d_select: select case(dir_gwdo)
    case("orogwd_2deg")
       nx = 180
       ny = 90
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 4
       scalefactor = 0.02
       dx = 2.0
       dy = 2.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.0
       known_lon =   1.0
    case("orogwd_1deg")
       nx = 360
       ny = 180
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 4
       scalefactor = 0.02
       dx = 1.0
       dy = 1.0
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.5
       known_lon =   0.5
    case("orogwd_30m")
       nx = 720
       ny = 360
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 4
       scalefactor = 0.02
       dx = 0.5
       dy = 0.5
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.75
       known_lon = 0.25
    case("orogwd_10m")
       nx = 2160
       ny = 1080
       nz = 1
       isigned     = 0
       endian      = 0
       wordsize    = 2
       scalefactor = 0.02
       dx = 0.16666667
       dy = 0.16666667
       known_x = 1.0
       known_y = 1.0
       known_lat = -89.916667
       known_lon = 0.0833333
    case default
 end select var2d_select

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') &
       trim(geog_data_path)//trim(dir_gwdo)//'/var/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)


 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 xarray(1:nx,1:ny) = rarray(1:nx,1:ny,1)

 call map_set(PROJ_LATLON, proj,  &
              latinc = dy,        &
              loninc = dx,        &
              knowni = known_x,   &
              knownj = known_y,   &
              lat1   = known_lat, &
              lon1   = known_lon)

 interp_list(1) = AVERAGE4
 interp_list(2) = AVERAGE4
 interp_list(3) = AVERAGE4
 interp_list(4) = AVERAGE4
 interp_list(5) = 0

 do iCell = 1,nCells
    if(landmask(iCell) == 1) then
       lat = latCell(iCell) * DEG_PER_RAD
       lon = lonCell(iCell) * DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       var2d(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                                      0.0_RKIND,interp_list,1)
    endif
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate VAR2D'

 end subroutine init_atm_static_orogwd

!==================================================================================================
 subroutine init_atm_check_read_error(istatus, fname)
!==================================================================================================
 implicit none

 integer, intent(in) :: istatus
 character (len=*), intent(in) :: fname

 if (istatus /= 0) then
     call mpas_dmpar_global_abort('ERROR: Could not read file '//trim(fname))
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
