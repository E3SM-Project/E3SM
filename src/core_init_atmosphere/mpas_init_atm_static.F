! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
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
 use mpas_configure
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
 subroutine init_atm_static(mesh)
!==================================================================================================

!inout arguments:
 type(mesh_type),intent(inout):: mesh

!local variables:
 type(proj_info):: proj
 type(dm_info),pointer :: dminfo

 character(len=StrKIND):: fname

 integer:: nx,ny,nz
 integer:: endian,isigned,istatus,wordsize
 integer:: i,j,k
 integer:: iCell,iEdge,iVtx,iPoint,iTileStart,iTileEnd,jTileStart,jTileEnd
 integer,dimension(5) :: interp_list
 integer,dimension(:),allocatable  :: nhs
 integer,dimension(:,:),allocatable:: ncat
      
 real(kind=4):: scalefactor
 real(kind=4),dimension(:,:,:),allocatable:: rarray

 real(kind=RKIND):: r_earth
 real(kind=RKIND):: lat,lon,x,y
 real(kind=RKIND):: lat_pt,lon_pt
 real(kind=RKIND),dimension(:,:),allocatable  :: soiltemp_1deg
 real(kind=RKIND),dimension(:,:),allocatable  :: maxsnowalb
 real(kind=RKIND),dimension(:,:,:),allocatable:: vegfra


!--------------------------------------------------------------------------------------------------

 dminfo => mesh % block % domain % dminfo

 write(0,*)
 write(0,*) '--- enter subroutine init_atm_static:'

!
! Scale all distances and areas from a unit sphere to one with radius sphere_radius
!

 r_earth = mesh % sphere_radius

 mesh % xCell % array = mesh % xCell % array * r_earth
 mesh % yCell % array = mesh % yCell % array * r_earth
 mesh % zCell % array = mesh % zCell % array * r_earth
 mesh % xVertex % array = mesh % xVertex % array * r_earth
 mesh % yVertex % array = mesh % yVertex % array * r_earth
 mesh % zVertex % array = mesh % zVertex % array * r_earth
 mesh % xEdge % array = mesh % xEdge % array * r_earth
 mesh % yEdge % array = mesh % yEdge % array * r_earth
 mesh % zEdge % array = mesh % zEdge % array * r_earth
 mesh % dvEdge % array = mesh % dvEdge % array * r_earth
 mesh % dcEdge % array = mesh % dcEdge % array * r_earth
 mesh % areaCell % array = mesh % areaCell % array * r_earth**2.0
 mesh % areaTriangle % array = mesh % areaTriangle % array * r_earth**2.0
 mesh % kiteAreasOnVertex % array = mesh % kiteAreasOnVertex % array * r_earth**2.0


 !
 ! Initialize Coriolis parameter field on edges and vertices
 !
 do iEdge=1,mesh % nEdges
    mesh % fEdge % array(iEdge)  = 2.0 * omega * sin(mesh % latEdge % array(iEdge))
 end do
 do iVtx=1,mesh % nVertices
    mesh % fVertex % array(iVtx) = 2.0 * omega * sin(mesh % latVertex % array(iVtx))
 end do


 !
 ! Compute weights used in advection and deformation calculation
 !
 call atm_initialize_advection_rk(mesh) 
 call atm_initialize_deformation_weights(mesh) 


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
 allocate(nhs(mesh%nCells))
 nhs(:) = 0
 mesh%ter%array(:) = 0.0

 do jTileStart = 1,20401,ny-6
    jTileEnd = jTileStart + ny - 1 - 6

    do iTileStart=1,42001,nx-6
       iTileEnd = iTileStart + nx - 1 - 6
       write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)// &
             'topo_30s/',iTileStart,'-',iTileEnd,'.',jTileStart,'-',jTileEnd
       write(0,*) trim(fname)

       call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                         scalefactor,wordsize,istatus)
       call init_atm_check_read_error(istatus, fname, dminfo)

       iPoint = 1
       do j=4,ny-3
       do i=4,nx-3
          lat_pt = -89.99583  + (jTileStart + j - 5) * 0.0083333333
          lon_pt = -179.99583 + (iTileStart + i - 5) * 0.0083333333
          lat_pt = lat_pt * PI / 180.0
          lon_pt = lon_pt * PI / 180.0

          iPoint = nearest_cell(lat_pt,lon_pt,iPoint,mesh%nCells,mesh%maxEdges, &
                                mesh%nEdgesOnCell%array,mesh%cellsOnCell%array, &
                                mesh%latCell%array,mesh%lonCell%array)
          mesh%ter%array(iPoint) = mesh%ter%array(iPoint) + rarray(i,j,1)
          nhs(iPoint) = nhs(iPoint) + 1
       end do
       end do

    end do
 end do

 do iCell = 1,mesh%nCells
    mesh%ter%array(iCell) = mesh%ter%array(iCell) / real(nhs(iCell))
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
 allocate(ncat(24,mesh%nCells))
 ncat(:,:) = 0
 mesh%lu_index%array(:) = 0.0

 do jTileStart = 1,20401,ny
    jTileEnd = jTileStart + ny - 1

    do iTileStart = 1,42001,nx
       iTileEnd = iTileStart + nx - 1
       write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)// &
             '/landuse_30s/',iTileStart,'-',iTileEnd,'.',jTileStart,'-',jTileEnd
       write(0,*) trim(fname)

       call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                         scalefactor,wordsize,istatus)
       call init_atm_check_read_error(istatus, fname, dminfo)

       iPoint = 1
       do j=1,ny
       do i=1,nx
          lat_pt = -89.99583  + (jTileStart + j - 2) * 0.0083333333
          lon_pt = -179.99583 + (iTileStart + i - 2) * 0.0083333333
          lat_pt = lat_pt * PI / 180.0
          lon_pt = lon_pt * PI / 180.0

          iPoint = nearest_cell(lat_pt,lon_pt,iPoint,mesh%nCells,mesh%maxEdges, &
                                mesh%nEdgesOnCell%array,mesh%cellsOnCell%array, &
                                mesh%latCell%array,mesh%lonCell%array)
          ncat(int(rarray(i,j,1)),iPoint) = ncat(int(rarray(i,j,1)),iPoint) + 1
       end do
       end do

    end do
 end do

 do iCell = 1,mesh%nCells
    mesh%lu_index%array(iCell) = 1
    do i = 2,24
       if(ncat(i,iCell) > ncat(mesh%lu_index%array(iCell),iCell)) then
          mesh%lu_index%array(iCell) = i
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
 allocate(ncat(16,mesh%nCells))
 ncat(:,:) = 0
 mesh%soilcat_top%array(:) = 0.0

 do jTileStart = 1,20401,ny
    jTileEnd = jTileStart + ny - 1

    do iTileStart = 1,42001,nx
       iTileEnd = iTileStart + nx - 1
       write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)// &
             '/soiltype_top_30s/',iTileStart,'-',iTileEnd,'.',jTileStart,'-',jTileEnd
       write(0,*) trim(fname)

       call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                         scalefactor,wordsize,istatus)
       call init_atm_check_read_error(istatus, fname, dminfo)

       iPoint = 1
       do j=1,ny
       do i=1,nx
          lat_pt = -89.99583  + (jTileStart + j - 2) * 0.0083333333
          lon_pt = -179.99583 + (iTileStart + i - 2) * 0.0083333333
          lat_pt = lat_pt * PI / 180.0
          lon_pt = lon_pt * PI / 180.0

          iPoint = nearest_cell(lat_pt,lon_pt,iPoint,mesh%nCells,mesh%maxEdges, &
                                mesh%nEdgesOnCell%array,mesh%cellsOnCell%array, &
                                mesh%latCell%array,mesh%lonCell%array)
          ncat(int(rarray(i,j,1)),iPoint) = ncat(int(rarray(i,j,1)),iPoint) + 1
       end do
       end do

    end do
 end do

 do iCell = 1,mesh%nCells
    mesh%soilcat_top%array(iCell) = 1
    do i = 2,16
       if(ncat(i,iCell) > ncat(mesh%soilcat_top%array(iCell),iCell)) then
          mesh%soilcat_top%array(iCell) = i
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
 allocate(ncat(16,mesh%nCells))
 ncat(:,:) = 0
 mesh%soilcat_bot%array(:) = 0.0

 do jTileStart = 1,20401,ny
    jTileEnd = jTileStart + ny - 1

    do iTileStart = 1,42001,nx
       iTileEnd = iTileStart + nx - 1
       write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)// &
             '/soiltype_bot_30s/',iTileStart,'-',iTileEnd,'.',jTileStart,'-',jTileEnd
       write(0,*) trim(fname)

       call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                         scalefactor,wordsize,istatus)
       call init_atm_check_read_error(istatus, fname, dminfo)

       iPoint = 1
       do j=1,ny
       do i=1,nx
          lat_pt = -89.99583  + (jTileStart + j - 2) * 0.0083333333
          lon_pt = -179.99583 + (iTileStart + i - 2) * 0.0083333333
          lat_pt = lat_pt * PI / 180.0
          lon_pt = lon_pt * PI / 180.0

          iPoint = nearest_cell(lat_pt,lon_pt,iPoint,mesh%nCells,mesh%maxEdges, &
                                mesh%nEdgesOnCell%array,mesh%cellsOnCell%array, &
                                mesh%latCell%array,mesh%lonCell%array)
          ncat(int(rarray(i,j,1)),iPoint) = ncat(int(rarray(i,j,1)),iPoint) + 1
       end do
       end do

    end do
 end do

 do iCell =1, mesh%nCells
    mesh%soilcat_bot%array(iCell) = 1
    do i = 2,16
       if(ncat(i,iCell) > ncat(mesh%soilcat_bot%array(iCell),iCell)) then
          mesh%soilcat_bot%array(iCell) = i
       end if
    end do
 end do
 deallocate(rarray)
 deallocate(ncat)
 write(0,*) '--- end interpolate SOILCAT_BOT'


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! KLUDGE TO FIX SOIL TYPE OVER ANTARCTICA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 where (mesh%lu_index%array == 24) mesh%soilcat_top%array = 16
 where (mesh%lu_index%array == 24) mesh%soilcat_bot%array = 16

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CORRECT INCONSISTENT SOIL AND LAND USE DATA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 do iCell = 1,mesh%nCells
    if (mesh%lu_index%array(iCell) == 16 .or. &
        mesh%soilcat_top%array(iCell) == 14 .or. &
        mesh%soilcat_bot%array(iCell) == 14) then
        if (mesh%lu_index%array(iCell) /= 16) then
            write(0,*) 'Turning lu_index into water at ', iCell
            mesh%lu_index%array(iCell) = 16
        end if
        if (mesh%soilcat_top%array(iCell) /= 14) then
            write(0,*) 'Turning soilcat_top into water at ', iCell
            mesh%soilcat_top%array(iCell) = 14
        end if
        if (mesh%soilcat_bot%array(iCell) /= 14) then
            write(0,*) 'Turning soilcat_bot into water at ', iCell
            mesh%soilcat_bot%array(iCell) = 14
        end if
    end if
 end do


!
! Derive LANDMASK
!
 mesh%landmask%array(:) = 0
 do iCell=1, mesh%nCells
    if (mesh%lu_index%array(iCell) /= 16) mesh%landmask%array(iCell) = 1
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
 mesh%soiltemp%array(:) = 0.0

 call map_set(PROJ_LATLON, proj,  &
              latinc = 1.0_RKIND, &
              loninc = 1.0_RKIND, &
              knowni = 1.0_RKIND, &
              knownj = 1.0_RKIND, &
              lat1 = -89.5_RKIND, &
              lon1 = -179.5_RKIND)

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)// &
       'soiltemp_1deg/',1,'-',180,'.',1,'-',180
 write(0,*) trim(fname)

 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned, endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus, fname, dminfo)
 soiltemp_1deg(-2:180,-2:183) = rarray(1:183,1:186,1)

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)// &
            'soiltemp_1deg/',181,'-',360,'.',1,'-',180
 write(0,*) trim(fname)

 call read_geogrid(fname, len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                        scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname,dminfo)
 soiltemp_1deg(181:363,-2:183) = rarray(4:186,1:186,1)

 interp_list(1) = FOUR_POINT
 interp_list(2) = W_AVERAGE4
 interp_list(3) = W_AVERAGE16
 interp_list(4) = SEARCH
 interp_list(5) = 0

 do iCell = 1,mesh%nCells
  
    if(mesh%landmask%array(iCell) == 1) then
       lat = mesh % latCell % array(iCell) * DEG_PER_RAD
       lon = mesh % lonCell % array(iCell) * DEG_PER_RAD
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
       mesh%soiltemp%array(iCell) = interp_sequence(x,y,1,soiltemp_1deg,-2,363,-2,183, &
                                           1,1,0.0_RKIND,interp_list,1)
    else
       mesh%soiltemp%array(iCell) = 0.0
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
 mesh%snoalb%array(:) = 0.0

 call map_set(PROJ_LATLON, proj,  &
              latinc = 1.0_RKIND, &
              loninc = 1.0_RKIND, &
              knowni = 1.0_RKIND, &
              knownj = 1.0_RKIND, &
              lat1 = -89.5_RKIND, &
              lon1 = -179.5_RKIND)

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)// &
       'maxsnowalb/',1,'-',180,'.',1,'-',180
 write(0,*) trim(fname)

 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, & 
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname,dminfo)
 maxsnowalb(-2:180,-2:183) = rarray(1:183,1:186,1)

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)// &
       'maxsnowalb/',181,'-',360,'.',1,'-',180
 write(0,*) trim(fname)

 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus, fname, dminfo)
 maxsnowalb(181:363,-2:183) = rarray(4:186,1:186,1)

 interp_list(1) = FOUR_POINT
 interp_list(2) = W_AVERAGE4
 interp_list(3) = W_AVERAGE16
 interp_list(4) = SEARCH
 interp_list(5) = 0

 do iCell = 1,mesh%nCells
  
    if(mesh%landmask%array(iCell) == 1) then
       lat = mesh % latCell % array(iCell) * DEG_PER_RAD
       lon = mesh % lonCell % array(iCell) * DEG_PER_RAD
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
       mesh%snoalb%array(iCell) = interp_sequence(x,y,1,maxsnowalb,-2,363,-2,183, &
                                         1,1,0.0_RKIND,interp_list,1)
    else
       mesh%snoalb%array(iCell) = 0.0
    end if

 end do
 mesh%snoalb%array(:) = mesh%snoalb%array(:) / 100.0
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
 mesh%greenfrac%array(:,:) = 0.0

 call map_set(PROJ_LATLON, proj,    &
              latinc = 0.144_RKIND, &
              loninc = 0.144_RKIND, &
              knowni = 1.0_RKIND,   &
              knownj = 1.0_RKIND,   &
              lat1 = -89.928_RKIND, &
              lon1 = -179.928_RKIND)

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)// &
       'greenfrac/',1,'-',1250,'.',1,'-',1250
 write(0,*) trim(fname)

 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, & 
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname,dminfo)
 vegfra(-2:1250,-2:1253,1:12) = rarray(1:1253,1:1256,1:12)

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)// &
       'greenfrac/',1251,'-',2500,'.',1,'-',1250
 write(0,*) trim(fname)

 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname,dminfo)
 vegfra(1251:2503,-2:1253,1:12) = rarray(4:1256,1:1256,1:12)

 do iCell = 1,mesh%nCells

    if (mesh%landmask%array(iCell) == 1) then
       lat = mesh % latCell % array(iCell) * DEG_PER_RAD
       lon = mesh % lonCell % array(iCell) * DEG_PER_RAD
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
          mesh%greenfrac%array(k,iCell) = interp_sequence(x,y,k,vegfra,-2,2503,-2,1253, &
                                                 1,12,-1.e30_RKIND,interp_list,1)
       end do
    else
       mesh%greenfrac%array(:,iCell) = 0.0
    end if
    mesh%shdmin%array(iCell) = minval(mesh%greenfrac%array(:,iCell))
    mesh%shdmax%array(iCell) = maxval(mesh%greenfrac%array(:,iCell))
      
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
 mesh%albedo12m%array(:,:) = 0.0

 call map_set(PROJ_LATLON, proj,    &
              latinc = 0.144_RKIND, &
              loninc = 0.144_RKIND, &
              knowni = 1.0_RKIND,   &
              knownj = 1.0_RKIND,   &
              lat1 = -89.928_RKIND, &
              lon1 = -179.928_RKIND)

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)// &
       'albedo_ncep/',1,'-',1250,'.',1,'-',1250
 write(0,*) trim(fname)

 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor, wordsize, istatus)
 call init_atm_check_read_error(istatus,fname, dminfo)
 vegfra(-2:1250,-2:1253,1:12) = rarray(1:1253,1:1256,1:12)

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)// &
       'albedo_ncep/',1251,'-',2500,'.',1,'-',1250
 write(0,*) trim(fname)

 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, & 
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname,dminfo)
 vegfra(1251:2503,-2:1253,1:12) = rarray(4:1256,1:1256,1:12)

 do iCell = 1,mesh%nCells

    if (mesh%landmask%array(iCell) == 1) then
       lat = mesh % latCell % array(iCell) * DEG_PER_RAD
       lon = mesh % lonCell % array(iCell) * DEG_PER_RAD
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
          mesh%albedo12m%array(k,iCell) = interp_sequence(x,y,k,vegfra,-2,2503,-2,1253, &
                                                 1,12,0.0_RKIND,interp_list,1)
       end do
    else
       mesh%albedo12m%array(:,iCell) = 8.0
    end if
 end do
 deallocate(rarray)
 deallocate(vegfra)
 write(0,*) '--- end interpolate ALBEDO12M'


 end subroutine init_atm_static

!==================================================================================================
 subroutine init_atm_static_orogwd(mesh)
!==================================================================================================

!inout arguments:
 type(mesh_type),intent(inout):: mesh

!local variables:
 type(proj_info):: proj
 type(dm_info),pointer :: dminfo

 character(len=StrKIND):: mess
 character(len=StrKIND):: fname
 character(len=StrKIND):: dir_gwdo

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


!--------------------------------------------------------------------------------------------------

 dminfo => mesh % block % domain % dminfo

 write(0,*)
 write(0,*) '--- enter subroutine init_atm_static_orogwd:'

!
! Interpolate VARSSO:
 mesh%varsso%array(:) = 0.0_RKIND
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
 allocate(nhs(mesh%nCells))
 nhs(:) = 0
 rarray(:,:,:) = 0._RKIND
 do jTileStart = 1,13801,ny
    jTileEnd = jTileStart + ny - 1

    do iTileStart = 1,42601,nx
       iTileEnd = iTileStart + nx -1
       write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)//'varsso/', &
             iTileStart,'-',iTileEnd,'.',jTileStart,'-',jTileEnd
       write(0,*) trim(fname)

       call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                         scalefactor,wordsize,istatus)
       call init_atm_check_read_error(istatus,fname,dminfo)

       iPoint = 1
       do j = 1,ny
       do i = 1,nx
          lat_pt = known_lat + (jTileStart + j - 2) * dy
          lon_pt = known_lon + (iTileStart + i - 2) * dx
          lat_pt = lat_pt * PI / 180.0
          lon_pt = lon_pt * PI / 180.0

          iPoint = nearest_cell(lat_pt,lon_pt,iPoint,mesh%nCells,mesh%maxEdges, &
                                mesh%nEdgesOnCell%array,mesh%cellsOnCell%array, &
                                mesh%latCell%array,mesh%lonCell%array)
          mesh%varsso%array(iPoint) = mesh%varsso%array(iPoint) + rarray(i,j,1)
          nhs(iPoint) = nhs(iPoint) + 1
       enddo
       enddo

    enddo
 enddo

 do iCell = 1,mesh%nCells
    if(nhs(iCell) .gt. 0) &
       mesh%varsso%array(iCell) = mesh%varsso%array(iCell) / real(nhs(iCell))
 enddo
 deallocate(rarray)
 deallocate(nhs)
 write(0,*) '--- end interpolate VARSSO'

!... statistic fields needed for the parameterization of gravity wavwe drag over orography. The
!input directory depends on the mesh resolution, and the mesh must be a uniform mesh.
 minMeshD  = minval(mesh%meshDensity%array(1:mesh%nCells))
 maxMeshD  = maxval(mesh%meshDensity%array(1:mesh%nCells))
 mindcEdge = minval(mesh%dcEdge%array(1:mesh%nEdges))
 maxdcEdge = maxval(mesh%dcEdge%array(1:mesh%nEdges))

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
 mesh%con%array(:) = 0.0_RKIND

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
       trim(config_geog_data_path)//trim(dir_gwdo)//'/con/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname,dminfo)
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

 do iCell = 1,mesh%nCells
    if(mesh % landmask % array(iCell) == 1) then
       lat = mesh % latCell % array(iCell) * DEG_PER_RAD
       lon = mesh % lonCell % array(iCell) * DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       mesh % con % array(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                                          0.0_RKIND,interp_list,1)
    endif
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate CON'

!
! Interpolate OA1:
!
 mesh%oa1%array(:) = 0.0_RKIND

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
       trim(config_geog_data_path)//trim(dir_gwdo)//'/oa1/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname,dminfo)
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

 do iCell = 1,mesh%nCells
    if(mesh % landmask % array(iCell) == 1) then
       lat = mesh % latCell % array(iCell) * DEG_PER_RAD
       lon = mesh % lonCell % array(iCell) * DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       mesh % oa1 % array(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                                          0.0_RKIND,interp_list,1)
    endif
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate OA1'

!
! Interpolate OA2:
 mesh%oa2%array(:) = 0.0_RKIND

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
       trim(config_geog_data_path)//trim(dir_gwdo)//'/oa2/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname,dminfo)
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

 do iCell = 1,mesh%nCells
    if(mesh % landmask % array(iCell) == 1) then
       lat = mesh % latCell % array(iCell) * DEG_PER_RAD
       lon = mesh % lonCell % array(iCell) * DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       mesh % oa2 % array(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                                          0.0_RKIND,interp_list,1)
     endif
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate OA2'

!
! Interpolate OA3:
!
 mesh%oa3%array(:) = 0.0_RKIND

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
       trim(config_geog_data_path)//trim(dir_gwdo)//'/oa3/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname,dminfo)
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

 do iCell = 1,mesh%nCells
    if(mesh % landmask % array(iCell) == 1) then
       lat = mesh % latCell % array(iCell) * DEG_PER_RAD
       lon = mesh % lonCell % array(iCell) * DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       mesh % oa3 % array(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                                          0.0_RKIND,interp_list,1)
    endif
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate OA3'

!
! Interpolate OA4:
!
 mesh%oa4%array(:) = 0.0_RKIND

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
       trim(config_geog_data_path)//trim(dir_gwdo)//'/oa4/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname,dminfo)
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

 do iCell = 1,mesh%nCells
    if(mesh % landmask % array(iCell) == 1) then
       lat = mesh % latCell % array(iCell) * DEG_PER_RAD
       lon = mesh % lonCell % array(iCell) * DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       mesh % oa4 % array(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                                          0.0_RKIND,interp_list,1)
    endif
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate OA4'

!
! Interpolate OL1:
!
 mesh%ol1%array(:) = 0.0_RKIND

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
       trim(config_geog_data_path)//trim(dir_gwdo)//'/ol1/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname,dminfo)
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

 do iCell = 1,mesh%nCells
    if(mesh % landmask % array(iCell) == 1) then
       lat = mesh % latCell % array(iCell) * DEG_PER_RAD
       lon = mesh % lonCell % array(iCell) * DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       mesh % ol1 % array(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                                          0.0_RKIND,interp_list,1)
    endif
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate OL1'

!
! Interpolate OL2:
!
 mesh%ol2%array(:) = 0.0_RKIND

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
       trim(config_geog_data_path)//trim(dir_gwdo)//'/ol2/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname,dminfo)
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

 do iCell = 1,mesh%nCells
    if(mesh % landmask % array(iCell) == 1) then
       lat = mesh % latCell % array(iCell) * DEG_PER_RAD
       lon = mesh % lonCell % array(iCell) * DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       mesh % ol2 % array(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                                          0.0_RKIND,interp_list,1)
    endif
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate OL2'

!
! Interpolate OL3:
!
 mesh%ol3%array(:) = 0.0_RKIND

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
       trim(config_geog_data_path)//trim(dir_gwdo)//'/ol3/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname,dminfo)
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

 do iCell = 1,mesh%nCells
    if(mesh % landmask % array(iCell) == 1) then
       lat = mesh % latCell % array(iCell) * DEG_PER_RAD
       lon = mesh % lonCell % array(iCell) * DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       mesh % ol3 % array(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                                          0.0_RKIND,interp_list,1)
    endif
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate OL3'

!
! Interpolate OL4:
!
 mesh%ol4%array(:) = 0.0_RKIND

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
       trim(config_geog_data_path)//trim(dir_gwdo)//'/ol4/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname,dminfo)
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

 do iCell = 1,mesh%nCells
    if(mesh % landmask % array(iCell) == 1) then
       lat = mesh % latCell % array(iCell) * DEG_PER_RAD
       lon = mesh % lonCell % array(iCell) * DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       mesh % ol4 % array(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                                          0.0_RKIND,interp_list,1)
    endif
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate OL4'

!
! Interpolate VAR2D:
!
 mesh%var2d%array(:) = 0.0_RKIND

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
       trim(config_geog_data_path)//trim(dir_gwdo)//'/var/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)


 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname,dminfo)
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

 do iCell = 1,mesh%nCells
    if(mesh % landmask % array(iCell) == 1) then
       lat = mesh % latCell % array(iCell) * DEG_PER_RAD
       lon = mesh % lonCell % array(iCell) * DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       mesh % var2d % array(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                                            0.0_RKIND,interp_list,1)
    endif
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate VAR2D'

 end subroutine init_atm_static_orogwd

!==================================================================================================
 subroutine init_atm_check_read_error(istatus, fname, dminfo)
!==================================================================================================
 implicit none

 integer, intent(in) :: istatus
 character (len=*), intent(in) :: fname
 type (dm_info), intent(in) :: dminfo

 if (istatus /= 0) then
     write(0,*) 'ERROR: Could not read file '//trim(fname)
     call mpas_dmpar_abort(dminfo)
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
