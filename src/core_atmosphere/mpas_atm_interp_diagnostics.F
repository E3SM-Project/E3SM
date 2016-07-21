! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_atm_interp_diagnostics
 use mpas_dmpar
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines
 use mpas_constants

 implicit none
 private
 public:: interp_diagnostics

 contains

!==================================================================================================
 subroutine interp_diagnostics(mesh, state, time_lev, diag)
!==================================================================================================

!input arguments:
 type (mpas_pool_type), intent(in)  :: mesh
 type (mpas_pool_type), intent(in) :: state
 integer, intent(in) :: time_lev              ! which time level to use from state

!inout arguments:
 type (mpas_pool_type), intent(inout) :: diag

!local variables:
 integer :: iCell,iVert,iVertD,k,kk
 integer, pointer :: nCells, nCellsSolve, nVertLevels, nVertices, vertexDegree, nIsoLevelsT, nIsoLevelsZ
 integer :: nVertLevelsP1
 integer, pointer :: index_qv, num_scalars
 integer, dimension(:,:), pointer :: cellsOnVertex

 type (field2DReal), pointer:: pressure_p_field

 real (kind=RKIND), dimension(:), pointer :: areaTriangle
 real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex
 
 real (kind=RKIND), dimension(:,:), pointer :: exner, height
 real (kind=RKIND), dimension(:,:), pointer :: pressure_b, pressure_p 
 real (kind=RKIND), dimension(:,:), pointer :: relhum, theta_m, vorticity
 real (kind=RKIND), dimension(:,:), pointer :: umeridional, uzonal, vvel
 real (kind=RKIND), dimension(:,:,:), pointer :: scalars

 real (kind=RKIND), dimension(:), pointer :: t_iso_levels
 real (kind=RKIND), dimension(:), pointer :: z_iso_levels
 real (kind=RKIND), dimension(:,:), pointer :: t_isobaric
 real (kind=RKIND), dimension(:,:), pointer :: z_isobaric
 real (kind=RKIND), dimension(:), pointer :: meanT_500_300

 real (kind=RKIND), dimension(:), pointer :: temperature_200hPa
 real (kind=RKIND), dimension(:), pointer :: temperature_250hPa
 real (kind=RKIND), dimension(:), pointer :: temperature_500hPa
 real (kind=RKIND), dimension(:), pointer :: temperature_700hPa
 real (kind=RKIND), dimension(:), pointer :: temperature_850hPa
 real (kind=RKIND), dimension(:), pointer :: temperature_925hPa

 real (kind=RKIND), dimension(:), pointer :: relhum_200hPa
 real (kind=RKIND), dimension(:), pointer :: relhum_250hPa
 real (kind=RKIND), dimension(:), pointer :: relhum_500hPa
 real (kind=RKIND), dimension(:), pointer :: relhum_700hPa
 real (kind=RKIND), dimension(:), pointer :: relhum_850hPa
 real (kind=RKIND), dimension(:), pointer :: relhum_925hPa

 real (kind=RKIND), dimension(:), pointer :: dewpoint_200hPa
 real (kind=RKIND), dimension(:), pointer :: dewpoint_250hPa
 real (kind=RKIND), dimension(:), pointer :: dewpoint_500hPa
 real (kind=RKIND), dimension(:), pointer :: dewpoint_700hPa
 real (kind=RKIND), dimension(:), pointer :: dewpoint_850hPa
 real (kind=RKIND), dimension(:), pointer :: dewpoint_925hPa

 real (kind=RKIND), dimension(:), pointer :: uzonal_200hPa
 real (kind=RKIND), dimension(:), pointer :: uzonal_250hPa
 real (kind=RKIND), dimension(:), pointer :: uzonal_500hPa
 real (kind=RKIND), dimension(:), pointer :: uzonal_700hPa
 real (kind=RKIND), dimension(:), pointer :: uzonal_850hPa
 real (kind=RKIND), dimension(:), pointer :: uzonal_925hPa

 real (kind=RKIND), dimension(:), pointer :: umeridional_200hPa
 real (kind=RKIND), dimension(:), pointer :: umeridional_250hPa
 real (kind=RKIND), dimension(:), pointer :: umeridional_500hPa
 real (kind=RKIND), dimension(:), pointer :: umeridional_700hPa
 real (kind=RKIND), dimension(:), pointer :: umeridional_850hPa
 real (kind=RKIND), dimension(:), pointer :: umeridional_925hPa

 real (kind=RKIND), dimension(:), pointer :: height_200hPa
 real (kind=RKIND), dimension(:), pointer :: height_250hPa
 real (kind=RKIND), dimension(:), pointer :: height_500hPa
 real (kind=RKIND), dimension(:), pointer :: height_700hPa
 real (kind=RKIND), dimension(:), pointer :: height_850hPa
 real (kind=RKIND), dimension(:), pointer :: height_925hPa

 real (kind=RKIND), dimension(:), pointer :: w_200hPa
 real (kind=RKIND), dimension(:), pointer :: w_250hPa
 real (kind=RKIND), dimension(:), pointer :: w_500hPa
 real (kind=RKIND), dimension(:), pointer :: w_700hPa
 real (kind=RKIND), dimension(:), pointer :: w_850hPa
 real (kind=RKIND), dimension(:), pointer :: w_925hPa

 real (kind=RKIND), dimension(:), pointer :: vorticity_200hPa
 real (kind=RKIND), dimension(:), pointer :: vorticity_250hPa
 real (kind=RKIND), dimension(:), pointer :: vorticity_500hPa
 real (kind=RKIND), dimension(:), pointer :: vorticity_700hPa
 real (kind=RKIND), dimension(:), pointer :: vorticity_850hPa
 real (kind=RKIND), dimension(:), pointer :: vorticity_925hPa

!--------------------
! fields added specifically as convective diagnostics
! WCS 201503

 real (kind=RKIND), dimension(:), pointer :: cape
 real (kind=RKIND), dimension(:), pointer :: cin
 real (kind=RKIND), dimension(:), pointer :: lcl
 real (kind=RKIND), dimension(:), pointer :: lfc
 real (kind=RKIND), dimension(:), pointer :: srh_0_1km
 real (kind=RKIND), dimension(:), pointer :: srh_0_3km
 real (kind=RKIND), dimension(:), pointer :: uzonal_surface
 real (kind=RKIND), dimension(:), pointer :: uzonal_1km
 real (kind=RKIND), dimension(:), pointer :: uzonal_6km
 real (kind=RKIND), dimension(:), pointer :: umeridional_surface
 real (kind=RKIND), dimension(:), pointer :: umeridional_1km
 real (kind=RKIND), dimension(:), pointer :: umeridional_6km
 real (kind=RKIND), dimension(:), pointer :: temperature_surface
 real (kind=RKIND), dimension(:), pointer :: dewpoint_surface


 real (kind=RKIND), parameter :: dev_motion = 7.5, z_bunker_bot = 0., z_bunker_top = 6000.
 real (kind=RKIND)            :: u_storm, v_storm, u_srh_bot, v_srh_bot, u_srh_top, v_srh_top
 real (kind=RKIND)            :: u_mean, v_mean, u_shear, v_shear, shear_magnitude
 real (kind=RKIND)            :: b_term, cape_out, cin_out
 real (kind=RKIND)            :: tempC, logrelhum
 real (kind=RKIND), dimension(:), allocatable :: dudz, dvdz, zp
 real (kind=RKIND), dimension(:), allocatable :: zrel, srh
 real (kind=RKIND), dimension(:), allocatable :: p_in, t_in, td_in

!--------------------

 real (kind=RKIND), dimension(:), pointer :: mslp

 real (kind=RKIND), dimension(:,:), allocatable :: pressure, pressureCp1, pressure2, pressure_v, temperature
 real (kind=RKIND), dimension(:,:), allocatable :: dewpoint

!local interpolated fields:
 integer :: nIntP
 real (kind=RKIND) :: w1,w2,z0,z1,z2
 real (kind=RKIND), dimension(:,:), allocatable :: field_in,press_in
 real (kind=RKIND), dimension(:,:), allocatable :: field_interp,press_interp
 
!--------------------------------------------------------------------------------------------------

! write(0,*)
! write(0,*) '--- enter subroutine interp_diagnostics:'

 call mpas_pool_get_dimension(mesh, 'nCells', nCells)
 call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
 call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
 call mpas_pool_get_dimension(mesh, 'nVertices', nVertices)
 call mpas_pool_get_dimension(mesh, 'vertexDegree', vertexDegree)
 call mpas_pool_get_dimension(mesh, 'nIsoLevelsT', nIsoLevelsT)
 call mpas_pool_get_dimension(mesh, 'nIsoLevelsZ', nIsoLevelsZ)
 call mpas_pool_get_dimension(state, 'index_qv', index_qv)
 call mpas_pool_get_dimension(state, 'num_scalars', num_scalars)

 nVertLevelsP1 = nVertLevels + 1

 call mpas_pool_get_array(mesh, 'cellsOnVertex', cellsOnVertex)
 call mpas_pool_get_array(mesh, 'areaTriangle', areaTriangle)
 call mpas_pool_get_array(mesh, 'kiteAreasOnVertex', kiteAreasOnVertex)

 call mpas_pool_get_array(mesh, 'zgrid', height)
 call mpas_pool_get_array(state, 'w', vvel, time_lev)
 call mpas_pool_get_array(state, 'theta_m', theta_m, time_lev)
 call mpas_pool_get_array(state, 'scalars', scalars, time_lev)

 call mpas_pool_get_field(diag, 'pressure_p', pressure_p_field)
 call mpas_dmpar_exch_halo_field(pressure_p_field)

 call mpas_pool_get_array(diag, 'exner', exner)
 call mpas_pool_get_array(diag, 'pressure_base', pressure_b)
 call mpas_pool_get_array(diag, 'pressure_p', pressure_p)
 call mpas_pool_get_array(diag, 'vorticity', vorticity)
 call mpas_pool_get_array(diag, 'uReconstructMeridional', umeridional)
 call mpas_pool_get_array(diag, 'uReconstructZonal', uzonal)
 call mpas_pool_get_array(diag, 'relhum', relhum)

 call mpas_pool_get_array(mesh, 't_iso_levels', t_iso_levels)
 call mpas_pool_get_array(mesh, 'z_iso_levels', z_iso_levels)
 call mpas_pool_get_array(diag, 't_isobaric', t_isobaric)
 call mpas_pool_get_array(diag, 'z_isobaric', z_isobaric)
 call mpas_pool_get_array(diag, 'meanT_500_300', meanT_500_300)

 call mpas_pool_get_array(diag, 'temperature_200hPa', temperature_200hPa)
 call mpas_pool_get_array(diag, 'temperature_250hPa', temperature_250hPa)
 call mpas_pool_get_array(diag, 'temperature_500hPa', temperature_500hPa)
 call mpas_pool_get_array(diag, 'temperature_700hPa', temperature_700hPa)
 call mpas_pool_get_array(diag, 'temperature_850hPa', temperature_850hPa)
 call mpas_pool_get_array(diag, 'temperature_925hPa', temperature_925hPa)

 call mpas_pool_get_array(diag, 'relhum_200hPa', relhum_200hPa)
 call mpas_pool_get_array(diag, 'relhum_250hPa', relhum_250hPa)
 call mpas_pool_get_array(diag, 'relhum_500hPa', relhum_500hPa)
 call mpas_pool_get_array(diag, 'relhum_700hPa', relhum_700hPa)
 call mpas_pool_get_array(diag, 'relhum_850hPa', relhum_850hPa)
 call mpas_pool_get_array(diag, 'relhum_925hPa', relhum_925hPa)

 call mpas_pool_get_array(diag, 'dewpoint_200hPa', dewpoint_200hPa)
 call mpas_pool_get_array(diag, 'dewpoint_250hPa', dewpoint_250hPa)
 call mpas_pool_get_array(diag, 'dewpoint_500hPa', dewpoint_500hPa)
 call mpas_pool_get_array(diag, 'dewpoint_700hPa', dewpoint_700hPa)
 call mpas_pool_get_array(diag, 'dewpoint_850hPa', dewpoint_850hPa)
 call mpas_pool_get_array(diag, 'dewpoint_925hPa', dewpoint_925hPa)

 call mpas_pool_get_array(diag, 'uzonal_200hPa', uzonal_200hPa)
 call mpas_pool_get_array(diag, 'uzonal_250hPa', uzonal_250hPa)
 call mpas_pool_get_array(diag, 'uzonal_500hPa', uzonal_500hPa)
 call mpas_pool_get_array(diag, 'uzonal_700hPa', uzonal_700hPa)
 call mpas_pool_get_array(diag, 'uzonal_850hPa', uzonal_850hPa)
 call mpas_pool_get_array(diag, 'uzonal_925hPa', uzonal_925hPa)

 call mpas_pool_get_array(diag, 'umeridional_200hPa', umeridional_200hPa)
 call mpas_pool_get_array(diag, 'umeridional_250hPa', umeridional_250hPa)
 call mpas_pool_get_array(diag, 'umeridional_500hPa', umeridional_500hPa)
 call mpas_pool_get_array(diag, 'umeridional_700hPa', umeridional_700hPa)
 call mpas_pool_get_array(diag, 'umeridional_850hPa', umeridional_850hPa)
 call mpas_pool_get_array(diag, 'umeridional_925hPa', umeridional_925hPa)

 call mpas_pool_get_array(diag, 'height_200hPa', height_200hPa)
 call mpas_pool_get_array(diag, 'height_250hPa', height_250hPa)
 call mpas_pool_get_array(diag, 'height_500hPa', height_500hPa)
 call mpas_pool_get_array(diag, 'height_700hPa', height_700hPa)
 call mpas_pool_get_array(diag, 'height_850hPa', height_850hPa)
 call mpas_pool_get_array(diag, 'height_925hPa', height_925hPa)

 call mpas_pool_get_array(diag, 'w_200hPa', w_200hPa)
 call mpas_pool_get_array(diag, 'w_250hPa', w_250hPa)
 call mpas_pool_get_array(diag, 'w_500hPa', w_500hPa)
 call mpas_pool_get_array(diag, 'w_700hPa', w_700hPa)
 call mpas_pool_get_array(diag, 'w_850hPa', w_850hPa)
 call mpas_pool_get_array(diag, 'w_925hPa', w_925hPa)

 call mpas_pool_get_array(diag, 'vorticity_200hPa', vorticity_200hPa)
 call mpas_pool_get_array(diag, 'vorticity_250hPa', vorticity_250hPa)
 call mpas_pool_get_array(diag, 'vorticity_500hPa', vorticity_500hPa)
 call mpas_pool_get_array(diag, 'vorticity_700hPa', vorticity_700hPa)
 call mpas_pool_get_array(diag, 'vorticity_850hPa', vorticity_850hPa)
 call mpas_pool_get_array(diag, 'vorticity_925hPa', vorticity_925hPa)

 call mpas_pool_get_array(diag, 'mslp', mslp)

 if(.not.allocated(pressure)    ) allocate(pressure(nVertLevels,nCells)      )
 if(.not.allocated(pressureCp1) ) allocate(pressureCp1(nVertLevels,nCells+1) )
 if(.not.allocated(pressure2)   ) allocate(pressure2(nVertLevelsP1,nCells)   )
 if(.not.allocated(pressure_v)  ) allocate(pressure_v(nVertLevels,nVertices) )
 if(.not.allocated(temperature) ) allocate(temperature(nVertLevels,nCells)   )
 if(.not.allocated(dewpoint) ) allocate(dewpoint(nVertLevels,nCells)   )

 t_iso_levels(1) = 30000.0
 t_iso_levels(2) = 35000.0
 t_iso_levels(3) = 40000.0
 t_iso_levels(4) = 45000.0
 t_iso_levels(5) = 50000.0

 z_iso_levels(1)  = 30000.0
 z_iso_levels(2)  = 35000.0
 z_iso_levels(3)  = 40000.0
 z_iso_levels(4)  = 45000.0
 z_iso_levels(5)  = 50000.0
 z_iso_levels(6)  = 55000.0
 z_iso_levels(7)  = 60000.0
 z_iso_levels(8)  = 65000.0
 z_iso_levels(9)  = 70000.0
 z_iso_levels(10) = 75000.0
 z_iso_levels(11) = 80000.0
 z_iso_levels(12) = 85000.0
 z_iso_levels(13) = 90000.0

!calculation of total pressure at cell centers (at mass points):
 do iCell = 1, nCells
 do k = 1, nVertLevels
    pressure(k,iCell)    = (pressure_p(k,iCell) + pressure_b(k,iCell)) / 100._RKIND
    pressureCp1(k,iCell) = pressure(k,iCell)
 enddo
 enddo
 do iCell = nCells+1, nCells+1
 do k = 1, nVertLevels
    pressureCp1(k,iCell)   = (pressure_p(k,iCell) + pressure_b(k,iCell)) / 100._RKIND
 enddo
 enddo

!calculation of total pressure at cell centers (at vertical velocity points):
 k = nVertLevelsP1
 do iCell = 1, nCells
    z0 = height(k,iCell)
    z1 = 0.5*(height(k,iCell)+height(k-1,iCell)) 
    z2 = 0.5*(height(k-1,iCell)+height(k-2,iCell))
    w1 = (z0-z2)/(z1-z2)
    w2 = 1.-w1
    !use log of pressure to avoid occurrences of negative top-of-the-model pressure.
    pressure2(k,iCell) = exp(w1*log(pressure(k-1,iCell))+w2*log(pressure(k-2,iCell)))
 enddo
 do k = 2, nVertLevels
 do iCell = 1, nCells
    w1 = (height(k,iCell)-height(k-1,iCell)) / (height(k+1,iCell)-height(k-1,iCell))
    w2 = (height(k+1,iCell)-height(k,iCell)) / (height(k+1,iCell)-height(k-1,iCell))
    pressure2(k,iCell) = w1*pressure(k,iCell) + w2*pressure(k-1,iCell)
 enddo
 enddo
 k = 1
 do iCell = 1, nCells
    z0 = height(k,iCell)
    z1 = 0.5*(height(k,iCell)+height(k+1,iCell)) 
    z2 = 0.5*(height(k+1,iCell)+height(k+2,iCell))
    w1 = (z0-z2)/(z1-z2)
    w2 = 1.-w1
    pressure2(k,iCell) = w1*pressure(k,iCell)+w2*pressure(k+1,iCell)
 enddo

!calculation of total pressure at cell vertices (at mass points):
 do iVert = 1, nVertices
    pressure_v(:,iVert) = 0._RKIND

    do k = 1, nVertLevels
    do iVertD = 1, vertexDegree
       pressure_v(k,iVert) = pressure_v(k,iVert) &
               + kiteAreasOnVertex(iVertD,iVert)*pressureCp1(k,cellsOnVertex(iVertD,iVert))
    enddo
    pressure_v(k,iVert) = pressure_v(k,iVert) / areaTriangle(iVert)
    enddo
 enddo

!calculation of temperature at cell centers:
 do iCell = 1,nCells
 do k = 1,nVertLevels
     temperature(k,iCell) = (theta_m(k,iCell)/(1._RKIND+rvord*scalars(index_qv,k,iCell)))*exner(k,iCell) 
     tempC = temperature(k,iCell) - 273.15
     relhum(k,iCell) = max(1.e-08,min(1.,relhum(k,iCell)))
     logrelhum = log(relhum(k,iCell))
     dewpoint(k,iCell) = 243.04*(logrelhum+(17.625*tempC/(243.04+tempC))) & 
                         /(17.625-logrelhum-((17.625*tempC)/(243.04+tempC))) 
 enddo
 enddo

!interpolation to fixed pressure levels for fields located at cells centers and at mass points:
 nIntP = 6
 if(.not.allocated(field_interp)) allocate(field_interp(nCells,nIntP) )
 if(.not.allocated(press_interp)) allocate(press_interp(nCells,nIntP) )
 do iCell = 1, nCells
    press_interp(iCell,1) = 200.0_RKIND
    press_interp(iCell,2) = 250.0_RKIND
    press_interp(iCell,3) = 500.0_RKIND
    press_interp(iCell,4) = 700.0_RKIND
    press_interp(iCell,5) = 850.0_RKIND
    press_interp(iCell,6) = 925.0_RKIND
 enddo

 if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevels))
 do iCell = 1, nCells
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    press_in(iCell,kk) = pressure(k,iCell)
 enddo
 enddo

 if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))
!... temperature:
 do iCell = 1, nCells
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    field_in(iCell,kk) = temperature(k,iCell)
 enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 temperature_200hPa(1:nCells) = field_interp(1:nCells,1)
 temperature_250hPa(1:nCells) = field_interp(1:nCells,2)
 temperature_500hPa(1:nCells) = field_interp(1:nCells,3)
 temperature_700hPa(1:nCells) = field_interp(1:nCells,4)
 temperature_850hPa(1:nCells) = field_interp(1:nCells,5)
 temperature_925hPa(1:nCells) = field_interp(1:nCells,6)
! write(0,*) '--- end interpolate temperature:'


!... relative humidity:
 do iCell = 1, nCells
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    field_in(iCell,kk) = relhum(k,iCell)
 enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 relhum_200hPa(1:nCells) = field_interp(1:nCells,1)
 relhum_250hPa(1:nCells) = field_interp(1:nCells,2)
 relhum_500hPa(1:nCells) = field_interp(1:nCells,3)
 relhum_700hPa(1:nCells) = field_interp(1:nCells,4)
 relhum_850hPa(1:nCells) = field_interp(1:nCells,5)
 relhum_925hPa(1:nCells) = field_interp(1:nCells,6)
! write(0,*) '--- end interpolate relative humidity:'

!... dewpoint
 do iCell = 1, nCells
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    field_in(iCell,kk) = dewpoint(k,iCell)
 enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 dewpoint_200hPa(1:nCells) = field_interp(1:nCells,1)
 dewpoint_250hPa(1:nCells) = field_interp(1:nCells,2)
 dewpoint_500hPa(1:nCells) = field_interp(1:nCells,3)
 dewpoint_700hPa(1:nCells) = field_interp(1:nCells,4)
 dewpoint_850hPa(1:nCells) = field_interp(1:nCells,5)
 dewpoint_925hPa(1:nCells) = field_interp(1:nCells,6)
! write(0,*) '--- end interpolate relative humidity:'

!... u zonal wind:
 do iCell = 1, nCells
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    field_in(iCell,kk) = uzonal(k,iCell)
 enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 uzonal_200hPa(1:nCells) = field_interp(1:nCells,1)
 uzonal_250hPa(1:nCells) = field_interp(1:nCells,2)
 uzonal_500hPa(1:nCells) = field_interp(1:nCells,3)
 uzonal_700hPa(1:nCells) = field_interp(1:nCells,4)
 uzonal_850hPa(1:nCells) = field_interp(1:nCells,5)
 uzonal_925hPa(1:nCells) = field_interp(1:nCells,6)
! write(0,*) '--- end interpolate zonal wind:'

!... u meridional wind:
 do iCell = 1, nCells
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    field_in(iCell,kk) = umeridional(k,iCell)
 enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 umeridional_200hPa(1:nCells) = field_interp(1:nCells,1)
 umeridional_250hPa(1:nCells) = field_interp(1:nCells,2)
 umeridional_500hPa(1:nCells) = field_interp(1:nCells,3)
 umeridional_700hPa(1:nCells) = field_interp(1:nCells,4)
 umeridional_850hPa(1:nCells) = field_interp(1:nCells,5)
 umeridional_925hPa(1:nCells) = field_interp(1:nCells,6)
! write(0,*) '--- end interpolate meridional wind:'

 if(allocated(field_in)) deallocate(field_in)
 if(allocated(press_in)) deallocate(press_in)

!interpolation to fixed pressure levels for fields located at cells centers and at vertical
!velocity points:
 if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevelsP1))
 do iCell = 1, nCells
 do k = 1, nVertLevelsP1
    kk = nVertLevelsP1+1-k
    press_in(iCell,kk) = pressure2(k,iCell)
 enddo
 enddo

 if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevelsP1))
 !... height:
 do iCell = 1, nCells
 do k = 1, nVertLevelsP1
    kk = nVertLevelsP1+1-k
    field_in(iCell,kk) = height(k,iCell)
 enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevelsP1,nIntP,press_in,field_in,press_interp,field_interp)
 height_200hPa(1:nCells) = field_interp(1:nCells,1)
 height_250hPa(1:nCells) = field_interp(1:nCells,2)
 height_500hPa(1:nCells) = field_interp(1:nCells,3)
 height_700hPa(1:nCells) = field_interp(1:nCells,4)
 height_850hPa(1:nCells) = field_interp(1:nCells,5)
 height_925hPa(1:nCells) = field_interp(1:nCells,6)
! write(0,*) '--- end interpolate height:'
 
!... vertical velocity
 do iCell = 1, nCells
 do k = 1, nVertLevelsP1
    kk = nVertLevelsP1+1-k
    field_in(iCell,kk) = vvel(k,iCell)
 enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevelsP1,nIntP,press_in,field_in,press_interp,field_interp)
 w_200hPa(1:nCells) = field_interp(1:nCells,1)
 w_250hPa(1:nCells) = field_interp(1:nCells,2)
 w_500hPa(1:nCells) = field_interp(1:nCells,3)
 w_700hPa(1:nCells) = field_interp(1:nCells,4)
 w_850hPa(1:nCells) = field_interp(1:nCells,5)
 w_925hPa(1:nCells) = field_interp(1:nCells,6)
! write(0,*) '--- end interpolate vertical velocity:'
 
 if(allocated(field_interp)) deallocate(field_interp)
 if(allocated(press_interp)) deallocate(press_interp)

!interpolation to fixed pressure levels for fields located at cell vertices and at mass points:
 nIntP = 6
 if(.not.allocated(field_interp)) allocate(field_interp(nVertices,nIntP) )
 if(.not.allocated(press_interp)) allocate(press_interp(nVertices,nIntP) )
 do iVert = 1, nVertices
    press_interp(iVert,1) = 200.0_RKIND
    press_interp(iVert,2) = 250.0_RKIND
    press_interp(iVert,3) = 500.0_RKIND
    press_interp(iVert,4) = 700.0_RKIND
    press_interp(iVert,5) = 850.0_RKIND
    press_interp(iVert,6) = 925.0_RKIND
 enddo

 if(allocated(field_in)) deallocate(field_in)
 if(allocated(press_in)) deallocate(press_in)

 if(.not.allocated(press_in)) allocate(press_in(nVertices,nVertLevels))
 do iVert = 1, nVertices
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    press_in(iVert,kk) = pressure_v(k,iVert)
 enddo
 enddo

 if(.not.allocated(field_in)) allocate(field_in(nVertices,nVertLevels))
!... relative vorticity:
 do iVert = 1, nVertices
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    field_in(iVert,kk) = vorticity(k,iVert)
 enddo
 enddo
 call interp_tofixed_pressure(nVertices,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 vorticity_200hPa(1:nVertices) = field_interp(1:nVertices,1)
 vorticity_250hPa(1:nVertices) = field_interp(1:nVertices,2)
 vorticity_500hPa(1:nVertices) = field_interp(1:nVertices,3)
 vorticity_700hPa(1:nVertices) = field_interp(1:nVertices,4)
 vorticity_850hPa(1:nVertices) = field_interp(1:nVertices,5)
 vorticity_925hPa(1:nVertices) = field_interp(1:nVertices,6)
! write(0,*) '--- end interpolate relative vorticity:'

!... compute SLP (requires temp, height, pressure, qvapor)
 call compute_slp(nCells, nVertLevels, num_scalars, temperature, height, pressure, index_qv, scalars, mslp)
 mslp(:) = mslp(:) * 100.0   ! Convert from hPa to Pa
!... alternative way
!do iCell = 1, nCells
!   mslp(iCell) = diag % surface_pressure % array(iCell) + 11.38*height(1,iCell)
!   mslp(iCell) = mslp(iCell)/100.
!enddo

!--------------------------------
! convective diagnostics, WCS 201503

 call mpas_pool_get_array(diag, 'cape', cape)
 call mpas_pool_get_array(diag, 'cin',  cin)
 call mpas_pool_get_array(diag, 'lcl',  lcl)
 call mpas_pool_get_array(diag, 'lfc',  lfc)
 call mpas_pool_get_array(diag, 'srh_0_1km', srh_0_1km)
 call mpas_pool_get_array(diag, 'srh_0_3km', srh_0_3km)
 call mpas_pool_get_array(diag, 'uzonal_surface', uzonal_surface)
 call mpas_pool_get_array(diag, 'uzonal_1km', uzonal_1km)
 call mpas_pool_get_array(diag, 'uzonal_6km', uzonal_6km)
 call mpas_pool_get_array(diag, 'umeridional_surface', umeridional_surface)
 call mpas_pool_get_array(diag, 'umeridional_1km', umeridional_1km)
 call mpas_pool_get_array(diag, 'umeridional_6km', umeridional_6km)
 call mpas_pool_get_array(diag, 'temperature_surface', temperature_surface)
 call mpas_pool_get_array(diag, 'dewpoint_surface', dewpoint_surface)

 allocate(dudz(nVertLevels))
 allocate(dvdz(nVertLevels))
 allocate(zp(nVertLevels))
 allocate(zrel(nVertLevels+1))
 allocate(srh(nVertLevels+1))

 allocate(p_in(nVertLevels))
 allocate(t_in(nVertLevels))
 allocate(td_in(nVertLevels))

! first the shear values.  We will use lowest model level velocity for surface velocity

 do iCell=1,nCellsSolve

   zp(1:nVertLevels) = 0.5*(height(1:nVertLevels,iCell)+height(2:nVertlevels+1,iCell)) - height(1,iCell)
   zrel(1:nVertLevels+1) = height(1:nVertLevels+1,iCell) - height(1,iCell)

   uzonal_surface(iCell) = uzonal(1,iCell)
   umeridional_surface(iCell) = umeridional(1,iCell)
   temperature_surface(iCell) = temperature(1,iCell)
   dewpoint_surface(iCell) = dewpoint(1,iCell)
   uzonal_1km(iCell)      = column_height_value(uzonal(1:nVertLevels,iCell),      zp, 1000., nVertLevels)
   umeridional_1km(iCell) = column_height_value(umeridional(1:nVertLevels,iCell), zp, 1000., nVertLevels)
   uzonal_6km(iCell)      = column_height_value(uzonal(1:nVertLevels,iCell),      zp, 6000., nVertLevels)
   umeridional_6km(iCell) = column_height_value(umeridional(1:nVertLevels,iCell), zp, 6000., nVertLevels)

!  storm-relative helicity
!  first, calculate storm motion, using Bunkers formula for right-moving storms

    u_srh_bot = uzonal(1,iCell)
    v_srh_bot = umeridional(1,iCell)
    if(z_bunker_bot .gt. zp(1)) then
      u_srh_bot = column_height_value( uzonal(1:nVertLevels,iCell), zp, z_bunker_bot, nVertLevels)
      v_srh_bot = column_height_value( umeridional(1:nVertLevels,iCell), zp, z_bunker_bot, nVertLevels)
    end if
    u_srh_top = column_height_value( uzonal(1:nVertLevels,iCell), zp, z_bunker_top, nVertLevels)
    v_srh_top = column_height_value( umeridional(1:nVertLevels,iCell), zp, z_bunker_top, nVertLevels)
    u_shear = u_srh_top - u_srh_bot
    v_shear = v_srh_top - v_srh_bot
    u_mean = integral_zstaggered(uzonal(1:nVertLevels,iCell),zrel,z_bunker_bot,z_bunker_top,nVertLevels,nVertLevelsP1)/(z_bunker_top-z_bunker_bot)
    v_mean = integral_zstaggered(umeridional(1:nVertLevels,iCell),zrel,z_bunker_bot,z_bunker_top,nVertLevels,nVertLevelsP1)/(z_bunker_top-z_bunker_bot)
    shear_magnitude = max(0.0001,sqrt(u_shear**2 + v_shear**2))
    u_storm = u_mean + dev_motion * v_shear/shear_magnitude
    v_storm = v_mean - dev_motion * u_shear/shear_magnitude

    !  calculate horizontal vorticity

    do k=2, nVertLevels-1
      dudz(k) = (uzonal(k,iCell)     -uzonal(k-1,iCell)     )/(0.5*(height(k+1,iCell)-height(k-1,iCell)))      
      dvdz(k) = (umeridional(k,iCell)-umeridional(k-1,iCell))/(0.5*(height(k+1,iCell)-height(k-1,iCell)))
    enddo
    dudz(1) = dudz(2)
    dvdz(1) = dvdz(2)
    dudz(nVertLevels) = dudz(nVertLevels-1)
    dvdz(nVertLevels) = dvdz(nVertLevels-1)

    do k=2,nVertLevels
      srh(k) = - (0.5*(uzonal(k,iCell)      + uzonal(k-1,iCell)      )-u_storm)*dvdz(k) &
               + (0.5*(umeridional(k,iCell) + umeridional(k-1,iCell) )-v_storm)*dudz(k)
    end do
    srh(1) = - (uzonal(1,iCell)     - u_storm)*dvdz(1) &
             + (umeridional(1,iCell)- v_storm)*dudz(1)
    srh(nVertLevelsP1) = srh(nVertLevels)

    do k=1, nVertLevels+1
      srh(k) = max(0.,srh(k))  ! discounting negative SRH
    end do

    srh_0_1km(iCell) = integral_zpoint(srh, zrel, 0., 1000., nVertLevelsP1)
    srh_0_3km(iCell) = integral_zpoint(srh, zrel, 0., 3000., nVertLevelsP1)

 end do

!  calculate cape and cin

 do iCell=1, nCellsSolve

   p_in(1:nVertLevels) = pressure(1:nVertLevels,iCell)
   t_in(1:nVertLevels) = temperature(1:nVertLevels,iCell) - 273.15
   td_in(1:nVertLevels) = dewpoint(1:nVertLevels,iCell)

!   do k=1,nVertLevels
!     relhum(k,iCell) = max(1.e-08,min(1.,relhum(k,iCell)))
!     td_in(k) = 243.04*(log(relhum(k,iCell))+((17.625*t_in(k))/(243.04+t_in(k)))) & 
!                      /(17.625-log(relhum(k,iCell))-((17.625*t_in(k))/(243.04+t_in(k)))) 
!   end do

   call getcape( nVertLevels, p_in, t_in, td_in, cape_out, cin_out )

   cape(iCell) = cape_out
   cin(iCell) = cin_out

 end do

 if(allocated(dudz)  ) deallocate(dudz)
 if(allocated(dvdz)  ) deallocate(dvdz)
 if(allocated(zp  )  ) deallocate(zp)
 if(allocated(zrel)  ) deallocate(zrel)
 if(allocated(srh )  ) deallocate(srh)
 if(allocated(p_in)  ) deallocate(p_in)
 if(allocated(t_in)  ) deallocate(t_in)
 if(allocated(td_in) ) deallocate(td_in)

 if(allocated(field_in    )) deallocate(field_in)
 if(allocated(press_in    )) deallocate(press_in)
 if(allocated(field_interp)) deallocate(field_interp)
 if(allocated(press_interp)) deallocate(press_interp)
 if(allocated(pressureCp1) ) deallocate(pressureCp1 )
 if(allocated(pressure_v)  ) deallocate(pressure_v  )


   !!!!!!!!!!! Additional temperature levels for vortex tracking !!!!!!!!!!!

   allocate(field_in(nCells, nVertLevels))
   allocate(press_in(nCells, nVertLevels))
   allocate(field_interp(nCells, nIsoLevelsT))
   allocate(press_interp(nCells, nIsoLevelsT))

   do k=1,nIsoLevelsT
      press_interp(:,k) = t_iso_levels(k)
   end do

   ! Additional temperature levels for vortex tracking
   do iCell=1,nCells
   do k=1,nVertLevels
      kk = nVertLevels+1-k
      field_in(iCell,kk) = temperature(k,iCell)
   end do
   end do

   do iCell=1,nCells
   do k=1,nVertLevels
      kk = nVertLevels+1-k
      press_in(iCell,kk) = pressure(k,iCell) * 100.0
   end do
   end do

   call interp_tofixed_pressure(nCells, nVertLevels, nIsoLevelsT, press_in, field_in, press_interp, field_interp)

   do k=1,nIsoLevelsT
      t_isobaric(k,1:nCells) = field_interp(1:nCells,k)
   end do


   !!!!!!!!!!! Calculate mean temperature in 500 hPa - 300 hPa layer !!!!!!!!!!!

   call compute_layer_mean(meanT_500_300, 50000.0, 30000.0, field_in, press_in)


   deallocate(field_in)
   deallocate(field_interp)
   deallocate(press_in)
   deallocate(press_interp)


   !!!!!!!!!!! Additional height levels for vortex tracking !!!!!!!!!!!

   allocate(field_in(nCells, nVertLevelsP1))
   allocate(press_in(nCells, nVertLevelsP1))
   allocate(field_interp(nCells, nIsoLevelsZ))
   allocate(press_interp(nCells, nIsoLevelsZ))

   do k=1,nIsoLevelsZ
      press_interp(:,k) = z_iso_levels(k)
   end do

   do iCell=1,nCells
   do k=1,nVertLevelsP1
      kk = nVertLevelsP1+1-k
      field_in(iCell,kk) = height(k,iCell)
   end do
   end do

   do iCell=1,nCells
   do k=1,nVertLevelsP1
      kk = nVertLevelsP1+1-k
      press_in(iCell,kk) = pressure2(k,iCell) * 100.0
   end do
   end do

   call interp_tofixed_pressure(nCells, nVertLevelsP1, nIsoLevelsZ, press_in, field_in, press_interp, field_interp)

   do k=1,nIsoLevelsZ
      z_isobaric(k,1:nCells) = field_interp(1:nCells,k)
   end do

   deallocate(field_in)
   deallocate(field_interp)
   deallocate(press_in)
   deallocate(press_interp)

 if(allocated(temperature) ) deallocate(temperature )
 if(allocated(pressure2)   ) deallocate(pressure2   )
 if(allocated(pressure)    ) deallocate(pressure    )
 if(allocated(dewpoint)    ) deallocate(dewpoint )

!formats:
! 201 format(i5,4(1x,e15.8))

 end subroutine interp_diagnostics

!==================================================================================================
 subroutine interp_tofixed_pressure(ncol,nlev_in,nlev_out,pres_in,field_in,pres_out,field_out)
!==================================================================================================

!input arguments:
 integer,intent(in):: ncol,nlev_in,nlev_out

 real(kind=RKIND),intent(in),dimension(ncol,nlev_in) :: pres_in,field_in
 real(kind=RKIND),intent(in),dimension(ncol,nlev_out):: pres_out

!output arguments:
 real(kind=RKIND),intent(out),dimension(ncol,nlev_out):: field_out

!local variables:
! integer:: i1,i2,icol,k,kk
 integer:: icol,k,kk
 integer:: kkstart,kount
 integer,dimension(ncol):: kupper

 real(kind=RKIND):: dpl,dpu

!--------------------------------------------------------------------------------------------------

!formats:
! 201 format(i5,8(1x,e15.8))

!write(0,*)
!write(0,*) '--- enter subroutine interp_tofixed_pressure:'
!write(0,*) '... ncol     = ',ncol
!write(0,*) '... nlev_in  = ',nlev_in
!write(0,*) '... nlev_out = ',nlev_out
!i1=1 ; i2=ncol
!do k = 1, nlev_in
!   write(0,201) k,pres_in(i1,k),field_in(i1,k),pres_in(i2,k),field_in(i2,k)
!enddo
!write(0,*)

 do icol = 1, ncol
    kupper(icol) = 1
 enddo

 do k = 1, nlev_out 

    kkstart = nlev_in
    do icol = 1, ncol
       kkstart = min0(kkstart,kupper(icol))
    enddo
    kount = 0

    do kk = kkstart, nlev_in-1
       do icol = 1, ncol
          if(pres_out(icol,k).gt.pres_in(icol,kk).and.pres_out(icol,k).le.pres_in(icol,kk+1)) then
             kupper(icol) = kk
             kount = kount + 1
!            write(0,201) kupper(icol),pres_out(icol,k),pres_in(icol,kk),pres_in(icol,kk+1)
          endif
       enddo

       if(kount.eq.ncol) then
          do icol = 1, ncol
             dpu = pres_out(icol,k) - pres_in(icol,kupper(icol))
             dpl = pres_in(icol,kupper(icol)+1) - pres_out(icol,k)
             field_out(icol,k) = (field_in(icol,kupper(icol))*dpl &
                               + field_in(icol,kupper(icol)+1)*dpu)/(dpl + dpu)
          end do
          goto 35
        end if
    enddo

    do icol = 1, ncol
       if(pres_out(icol,k) .lt. pres_in(icol,1)) then
          field_out(icol,k) = field_in(icol,1)*pres_out(icol,k)/pres_in(icol,1)
       elseif(pres_out(icol,k) .gt. pres_in(icol,nlev_in)) then
          field_out(icol,k) = field_in(icol,nlev_in)
       else
          dpu = pres_out(icol,k) - pres_in(icol,kupper(icol))
          dpl = pres_in(icol,kupper(icol)+1) - pres_out(icol,k)
          field_out(icol,k) = (field_in(icol,kupper(icol))*dpl &
                            + field_in(icol,kupper(icol)+1)*dpu)/(dpl + dpu)
       endif
    enddo

 35 continue
!   write(0,201) kupper(i1),pres_out(i1,k),pres_in(i1,kupper(i1)),pres_in(i1,kupper(i1)+1), &
!                field_out(i1,k),field_in(i1,kupper(i1)),field_in(i1,kupper(i1)+1)
!   write(0,201) kupper(i2),pres_out(i2,k),pres_in(i2,kupper(i2)),pres_in(i2,kupper(i2)+1), &
!                field_out(i2,k),field_in(i2,kupper(i2)),field_in(i2,kupper(i2)+1)

 enddo

 end subroutine interp_tofixed_pressure

 subroutine compute_slp(ncol,nlev_in,nscalars,t,height,p,index_qv,scalars,slp)

    implicit none

   !input arguments:
    integer, intent(in) :: ncol, nlev_in, nscalars
   
   !p: in mb
   !t: in K
   !scalars: in kg/kg
   !height: in m
    real(kind=RKIND), intent(in), dimension(nlev_in,ncol) :: p,t
    real(kind=RKIND), intent(in), dimension(nlev_in+1,ncol) :: height
    integer, intent(in) :: index_qv
    real(kind=RKIND), intent(in), dimension(nscalars,nlev_in,ncol) :: scalars
   
   !output arguments:
    real(kind=RKIND), intent(out), dimension(ncol) :: slp
   
   !local variables:
    integer :: icol, k, kcount
    integer :: klo, khi
   
    real(kind=RKIND) :: gamma, rr, grav
    parameter (rr=287.0, grav=9.80616, gamma=0.0065)
   
    real(kind=RKIND) :: tc, pconst
    parameter (tc=273.16+17.5, pconst=100.)
   
    logical mm5_test
    parameter (mm5_test=.true.)
   
    integer, dimension(:), allocatable :: level
    real(kind=RKIND), dimension(:), allocatable :: t_surf, t_msl
    real(kind=RKIND) :: plo , phi , tlo, thi , zlo , zhi
    real(kind=RKIND) :: p_at_pconst , t_at_pconst , z_at_pconst, z_half_lowest
   
    logical :: l1, l2, l3, found
   
   ! Find least zeta level that is PCONST Pa above the surface.  We later use this
   ! level to extrapolate a surface pressure and temperature, which is supposed
   ! to reduce the effect of the diurnal heating cycle in the pressure field.
   
    if (.not.allocated(level))  allocate(level(ncol))
    if (.not.allocated(t_surf)) allocate(t_surf(ncol))
    if (.not.allocated(t_msl))  allocate(t_msl(ncol))
   
    do icol = 1 , ncol
       level(icol) = -1
   
       k = 1
       found = .false.
       do while ( (.not. found) .and. (k.le.nlev_in))
             if ( p(k,icol) .lt. p(1,icol)-pconst ) then
                level(icol) = k
                found = .true.
             end if
             k = k+1
       end do
   
       if ( level(icol) .eq. -1 ) then
          write(0,*) 'Troubles finding level ', pconst,' above ground.'
          write(0,*) 'Problems first occur at (',icol,')'
          write(0,*) 'Surface pressure = ',p(1,icol),' hPa.'
          write(0,*) '*** MSLP field will not be computed'
          slp(:) = 0.0
          return
       end if
   
    end do
   
   ! Get temperature PCONST hPa above surface.  Use this to extrapolate
   ! the temperature at the surface and down to sea level.
   
    do icol = 1 , ncol
   
       klo = max ( level(icol) - 1 , 1      )
       khi = min ( klo + 1        , nlev_in - 1 )
   
       if ( klo .eq. khi ) then
          write(0,*) 'Trapping levels are weird.'
          write(0,*) 'klo = ',klo,', khi = ',khi, ': and they should not be equal.'
          call mpas_dmpar_global_abort('ERROR: Error_trapping_levels')
       end if
   
       plo = p(klo,icol)
       phi = p(khi,icol)
       tlo = t(klo,icol) * (1. + 0.608 * scalars(index_qv,klo,icol))
       thi = t(khi,icol) * (1. + 0.608 * scalars(index_qv,khi,icol))
       zlo = 0.5*(height(klo,icol)+height(klo+1,icol))
       zhi = 0.5*(height(khi,icol)+height(khi+1,icol))
   
       p_at_pconst = p(1,icol) - pconst
       t_at_pconst = thi-(thi-tlo)*log(p_at_pconst/phi)*log(plo/phi)
       z_at_pconst = zhi-(zhi-zlo)*log(p_at_pconst/phi)*log(plo/phi)
   
       t_surf(icol) = t_at_pconst*(p(1,icol)/p_at_pconst)**(gamma*rr/grav)
       t_msl(icol) = t_at_pconst+gamma*z_at_pconst
   !   if (icol.eq.500) then
   !      write(0,*) plo,phi,tlo,thi,zlo,zhi,p_at_pconst,t_at_pconst,z_at_pconst
   !      write(0,*) t_surf(icol),t_msl(icol),level(icol),klo,khi
   !      write(0,*) height(klo,icol),height(khi,icol),height(khi+1,icol)
   !   endif
   
    end do
   
   ! If we follow a traditional computation, there is a correction to the sea level
   ! temperature if both the surface and sea level temnperatures are *too* hot.
   
    if ( mm5_test ) then
       kcount = 0
       do icol = 1 , ncol
             l1 = t_msl(icol) .lt. tc
             l2 = t_surf(icol) .le. tc
             l3 = .not. l1
             if ( l2 .and. l3 ) then
                t_msl(icol) = tc
             else
                t_msl(icol) = tc - 0.005*(t_surf(icol)-tc)**2
                kcount = kcount+1
             end if
       end do
   !   write(0,*) 'These number of points had t_msl adjusted ', kcount
    end if
   
    do icol = 1 , ncol
       z_half_lowest=0.5*(height(1,icol)+height(2,icol))
       slp(icol) = p(1,icol) * exp((2.*grav*z_half_lowest)/ &
                                 (rr*(t_msl(icol)+t_surf(icol))))
    end do
   
    if (allocated(level))  deallocate(level)
    if (allocated(t_surf)) deallocate(t_surf)
    if (allocated(t_msl))  deallocate(t_msl)

 end subroutine compute_slp

!***********************************************************************
!
!  routine compute_layer_mean
!
!> \brief   Computes the mean of a field in the specified layer.
!> \author  Michael Duda
!> \date    3 July 2014
!> \details
!>  Given a 3d pressure field, press_in(nCells,nVertLevels), with pressure 
!>  increasing with vertical index, and a 3d field, 
!>  field_in(nCells,nVertLevels) with levels in the same order, this routine
!>  will compute the mean of the field for each column between pressures
!>  p1 and p2.
!
!----------------------------------------------------------------------- 
 subroutine compute_layer_mean(layerMean, p1, p2, field_in, press_in)

    implicit none

    real(kind=RKIND), dimension(:), intent(out) :: layerMean
    real(kind=RKIND), intent(in) :: p1, p2
    real(kind=RKIND), dimension(:,:), intent(in) :: field_in
    real(kind=RKIND), dimension(:,:), intent(in) :: press_in

    integer :: nCells, nVertLevels
    integer :: iCell, k
    integer :: k_bot, k_top
    real(kind=RKIND) :: p_bot, p_top
    real(kind=RKIND) :: wtop_p, wtop_m
    real(kind=RKIND) :: wbot_p, wbot_m
    real(kind=RKIND) :: wtotal, w
    real(kind=RKIND) :: temp


    !
    ! Get dimensions of input arrays
    !
    nCells = size(field_in, 1) 
    nVertLevels = size(field_in, 2) 


    !
    ! Check that pressure is increasing with index
    !
    if (press_in(1,1) > press_in(1,nVertLevels)) then
        write(stderrUnit,*) 'Error in compute_layer_mean: pressure should increase with index'
        layerMean(:) = 0.0
        return
    end if

    
    !
    ! Set the pressure at the top and bottom of the layer
    !
    if (p1 < p2) then
       p_top = p1
       p_bot = p2
    else
       p_top = p2
       p_bot = p1
    end if


    !
    ! For each column, compute the mean value of the field between p_bot and
    ! p_top, with the field weighted by delta-p in each layer
    !
    do iCell=1,nCells
       k_bot = -1
       k_top = -1

       ! Search for trapping levels: k_top is the index just above (or equal to)
       ! p_top, and k_bot is the index just above (or equal to) p_bot.
       do k=1,nVertLevels-1
          if (press_in(iCell,k) <= p_top .and. press_in(iCell,k+1) > p_top) then
             k_top = k
             wtop_p = (p_top - press_in(iCell,k)) / (press_in(iCell,k+1) - press_in(iCell,k))
             wtop_m = (press_in(iCell,k+1) - p_top) / (press_in(iCell,k+1) - press_in(iCell,k))
          end if
          if (press_in(iCell,k) <= p_bot .and. press_in(iCell,k+1) > p_bot) then
             k_bot = k
             wbot_m = (p_bot - press_in(iCell,k)) / (press_in(iCell,k+1) - press_in(iCell,k))
             wbot_p = (press_in(iCell,k+1) - p_bot) / (press_in(iCell,k+1) - press_in(iCell,k))
          end if
       end do

       if (k_top == -1 .or. k_bot == -1) then      ! Layer intersects top or bottom boundary

          layerMean(iCell) = 0.0  

       else if (k_top == k_bot) then               ! Layer lies entirely within a single model layer

          layerMean(iCell) = wtop_m * field_in(iCell,k_top) + wtop_p * field_in(iCell,k_top+1)
          layerMean(iCell) = layerMean(iCell) + wbot_m * field_in(iCell,k_bot) + wbot_p * field_in(iCell,k_bot+1)
          layerMean(iCell) = 0.5 * layerMean(iCell)

       else

          ! First layer: from p_top down to press_in(iCell,k_top+1)
          wtotal = press_in(iCell,k_top+1) - p_top  
          temp = wtop_m * field_in(iCell,k_top) + wtop_p * field_in(iCell,k_top+1)
          layerMean(iCell) = wtotal * 0.5 * (field_in(iCell,k_top+1) + temp)

          ! Middle layers
          do k=k_top+1,k_bot-1
             w = press_in(iCell,k+1) - press_in(iCell,k)
             wtotal = wtotal + w
             layerMean(iCell) = layerMean(iCell) + w * 0.5 * (field_in(iCell,k) + field_in(iCell,k+1))
          end do

          ! Last layer: from press_in(iCell,k_bot) down to p_bot
          w = p_bot - press_in(iCell,k_bot)
          wtotal = wtotal + w
          temp = wbot_m * field_in(iCell,k_bot) + wbot_p * field_in(iCell,k_bot+1)
          layerMean(iCell) = layerMean(iCell) + w * 0.5 * (field_in(iCell,k_bot) + temp)

          layerMean(iCell) = layerMean(iCell) / wtotal
       end if

    end do

 end subroutine compute_layer_mean

!-----------------------

      real (kind=RKIND) function column_height_value( column_values, z, z_interp, n )
      implicit none
      integer n
      real (kind=RKIND) :: column_values(n), z(n), z_interp, wz, wzp1
      integer :: kz, k
!  we assume height increases monotonically with n
      kz = 1
      do k=1,n
        if(z(k) <= z_interp) kz = k
      end do
      kz = min(kz,n-1)
      
      wz = (z(kz+1)-z_interp)/(z(kz+1)-z(kz))
      wzp1 = 1. - wz
      column_height_value = wz*column_values(kz) + wzp1*column_values(kz+1)

      end function column_height_value

!---------------------------

      real (kind=RKIND) function integral_zstaggered( column_values, z, zbot, ztop, n, np1 )
      implicit none
      integer n, np1
      real (kind=RKIND) :: column_values(n), z(np1), zbot, ztop
      real (kind=RKIND) :: zb, zt

      integer :: k

!  integral from z_bot to z_top, assume cell-average values (first-order integration)
!  z increases monotonically

      integral_zstaggered = 0.
      do k=1,n
        zb = max(z(k), zbot)
        zt = min(z(k+1), ztop)
        integral_zstaggered = integral_zstaggered + column_values(k)*max(0.,(zt-zb))
      end do
      end function integral_zstaggered

!---------------------------------

      real (kind=RKIND) function integral_zpoint( column_values, z, zbot, ztop, n )
      implicit none
      integer n
      real :: column_values(n), z(n), zbot, ztop
      real :: zb, zt, dz, zr_midpoint, midpoint_value

      integer :: k

!  integral from z_bot to z_top, assume point values (second-order integration)
!  z increases monotonically

      integral_zpoint = 0.
      do k=1,n-1
        zb = max(z(k), zbot)
        zt = min(z(k+1), ztop)
        dz = max(0.,zt-zb)
        zr_midpoint = (0.5*(zt+zb) - z(k))/(z(k+1)-z(k))
        midpoint_value = column_values(k) + (column_values(k+1)-column_values(k))*zr_midpoint
        integral_zpoint = integral_zpoint + dz*midpoint_value
      end do
      end function integral_zpoint

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

    subroutine getcape( nk , p_in , t_in , td_in, cape , cin )
    implicit none

    integer, intent(in) :: nk
    real (kind=RKIND), dimension(nk), intent(in) :: p_in,t_in,td_in
    real (kind=RKIND), intent(out) :: cape,cin

!-----------------------------------------------------------------------
!
!  getcape - a fortran90 subroutine to calculate Convective Available
!            Potential Energy (CAPE) from a sounding.
!
!  Version 1.02                           Last modified:  10 October 2008
!
!  Author:  George H. Bryan
!           Mesoscale and Microscale Meteorology Division
!           National Center for Atmospheric Research
!           Boulder, Colorado, USA
!           gbryan@ucar.edu
!
!  Disclaimer:  This code is made available WITHOUT WARRANTY.
!
!  References:  Bolton (1980, MWR, p. 1046) (constants and definitions)
!               Bryan and Fritsch (2004, MWR, p. 2421) (ice processes)
!
!-----------------------------------------------------------------------
!
!  Input:     nk - number of levels in the sounding (integer)
!
!           p_in - one-dimensional array of pressure (mb) (real)
!
!           t_in - one-dimensional array of temperature (C) (real)
!
!          td_in - one-dimensional array of dewpoint temperature (C) (real)
!
!  Output:  cape - Convective Available Potential Energy (J/kg) (real)
!
!            cin - Convective Inhibition (J/kg) (real)
!
!-----------------------------------------------------------------------
!  User options:

    real (kind=RKIND), parameter :: pinc = 100.0   ! Pressure increment (Pa)
                                      ! (smaller number yields more accurate
                                      !  results,larger number makes code 
                                      !  go faster)

    integer, parameter :: source = 2    ! Source parcel:
                                        ! 1 = surface
                                        ! 2 = most unstable (max theta-e)
                                        ! 3 = mixed-layer (specify ml_depth)

    real (kind=RKIND), parameter :: ml_depth =  200.0  ! depth (m) of mixed layer 
                                          ! for source=3

    integer, parameter :: adiabat = 1   ! Formulation of moist adiabat:
                                        ! 1 = pseudoadiabatic, liquid only
                                        ! 2 = reversible, liquid only
                                        ! 3 = pseudoadiabatic, with ice
                                        ! 4 = reversible, with ice

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
!            No need to modify anything below here:
!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

    logical :: doit,ice,cloud,not_converged
    integer :: k,kmax,n,nloop,i,orec
    real (kind=RKIND), dimension(nk) :: p,t,td,pi,q,th,thv,z,pt,pb,pc,pn,ptv

    real (kind=RKIND) :: the,maxthe,parea,narea,lfc
    real (kind=RKIND) :: th1,p1,t1,qv1,ql1,qi1,b1,pi1,thv1,qt,dp,dz,ps,frac
    real (kind=RKIND) :: th2,p2,t2,qv2,ql2,qi2,b2,pi2,thv2
    real (kind=RKIND) :: thlast,fliq,fice,tbar,qvbar,qlbar,qibar,lhv,lhs,lhf,rm,cpm
    real*8 :: avgth,avgqv
!    real (kind=RKIND) :: getqvs,getqvi,getthe

!-----------------------------------------------------------------------

    real (kind=RKIND), parameter :: g     = 9.81
    real (kind=RKIND), parameter :: p00   = 100000.0
    real (kind=RKIND), parameter :: cp    = 1005.7
    real (kind=RKIND), parameter :: rd    = 287.04
    real (kind=RKIND), parameter :: rv    = 461.5
    real (kind=RKIND), parameter :: xlv   = 2501000.0
    real (kind=RKIND), parameter :: xls   = 2836017.0
    real (kind=RKIND), parameter :: t0    = 273.15
    real (kind=RKIND), parameter :: cpv   = 1875.0
    real (kind=RKIND), parameter :: cpl   = 4190.0
    real (kind=RKIND), parameter :: cpi   = 2118.636
    real (kind=RKIND), parameter :: lv1   = xlv+(cpl-cpv)*t0
    real (kind=RKIND), parameter :: lv2   = cpl-cpv
    real (kind=RKIND), parameter :: ls1   = xls+(cpi-cpv)*t0
    real (kind=RKIND), parameter :: ls2   = cpi-cpv

    real (kind=RKIND), parameter :: rp00  = 1.0/p00
    real (kind=RKIND), parameter :: eps   = rd/rv
    real (kind=RKIND), parameter :: reps  = rv/rd
    real (kind=RKIND), parameter :: rddcp = rd/cp
    real (kind=RKIND), parameter :: cpdrd = cp/rd
    real (kind=RKIND), parameter :: cpdg  = cp/g

    real (kind=RKIND), parameter :: converge = 0.0002

    integer, parameter :: debug_level =   0

!-----------------------------------------------------------------------

!---- convert p,t,td to mks units; get pi,q,th,thv ----!

    do k=1,nk
        p(k) = 100.0*p_in(k)
        t(k) = 273.15+t_in(k)
       td(k) = 273.15+td_in(k)
       pi(k) = (p(k)*rp00)**rddcp
        q(k) = getqvs(p(k),td(k))
       th(k) = t(k)/pi(k)
      thv(k) = th(k)*(1.0+reps*q(k))/(1.0+q(k))
    enddo

!---- get height using the hydrostatic equation ----!

    z(1) = 0.0
    do k=2,nk
      dz = -cpdg*0.5*(thv(k)+thv(k-1))*(pi(k)-pi(k-1))
      z(k) = z(k-1) + dz
    enddo

!---- find source parcel ----!

  IF(source.eq.1)THEN
    ! use surface parcel
    kmax = 1

  ELSEIF(source.eq.2)THEN
    ! use most unstable parcel (max theta-e)

    IF(p(1).lt.50000.0)THEN
      ! first report is above 500 mb ... just use the first level reported
      kmax = 1
      maxthe = getthe(p(1),t(1),td(1),q(1))
    ELSE
      ! find max thetae below 500 mb
      maxthe = 0.0
      do k=1,nk
        if(p(k).ge.50000.0)then
          the = getthe(p(k),t(k),td(k),q(k))
          if( the.gt.maxthe )then
            maxthe = the
            kmax = k
          endif
        endif
      enddo
    ENDIF
    if(debug_level.ge.100) print *,'  kmax,maxthe = ',kmax,maxthe

  ELSEIF(source.eq.3)THEN
    ! use mixed layer

    IF( (z(2)-z(1)).gt.ml_depth )THEN
      ! the second level is above the mixed-layer depth:  just use the
      ! lowest level

      avgth = th(1)
      avgqv = q(1)
      kmax = 1

    ELSEIF( z(nk).lt.ml_depth )THEN
      ! the top-most level is within the mixed layer:  just use the
      ! upper-most level

      avgth = th(nk)
      avgqv = q(nk)
      kmax = nk

    ELSE
      ! calculate the mixed-layer properties:

      avgth = 0.0
      avgqv = 0.0
      k = 2
      if(debug_level.ge.100) print *,'  ml_depth = ',ml_depth
      if(debug_level.ge.100) print *,'  k,z,th,q:'
      if(debug_level.ge.100) print *,1,z(1),th(1),q(1)

      do while( (z(k).le.ml_depth) .and. (k.le.nk) )

        if(debug_level.ge.100) print *,k,z(k),th(k),q(k)

        avgth = avgth + 0.5*(z(k)-z(k-1))*(th(k)+th(k-1))
        avgqv = avgqv + 0.5*(z(k)-z(k-1))*(q(k)+q(k-1))

        k = k + 1

      enddo

      th2 = th(k-1)+(th(k)-th(k-1))*(ml_depth-z(k-1))/(z(k)-z(k-1))
      qv2 =  q(k-1)+( q(k)- q(k-1))*(ml_depth-z(k-1))/(z(k)-z(k-1))

      if(debug_level.ge.100) print *,999,ml_depth,th2,qv2

      avgth = avgth + 0.5*(ml_depth-z(k-1))*(th2+th(k-1))
      avgqv = avgqv + 0.5*(ml_depth-z(k-1))*(qv2+q(k-1))

      if(debug_level.ge.100) print *,k,z(k),th(k),q(k)

      avgth = avgth/ml_depth
      avgqv = avgqv/ml_depth

      kmax = 1

    ENDIF

    if(debug_level.ge.100) print *,avgth,avgqv

  ELSE

!    print *
!    print *,'  Unknown value for source'
!    print *
!    print *,'  source = ',source
!    print *
!    stop
    write(0,*) 'getcape: unknown value for source'
    return

  ENDIF

!---- define parcel properties at initial location ----!
    narea = 0.0

  if( (source.eq.1).or.(source.eq.2) )then
    k    = kmax
    th2  = th(kmax)
    pi2  = pi(kmax)
    p2   = p(kmax)
    t2   = t(kmax)
    thv2 = thv(kmax)
    qv2  = q(kmax)
    b2   = 0.0
  elseif( source.eq.3 )then
    k    = kmax
    th2  = avgth
    qv2  = avgqv
    thv2 = th2*(1.0+reps*qv2)/(1.0+qv2)
    pi2  = pi(kmax)
    p2   = p(kmax)
    t2   = th2*pi2
    b2   = g*( thv2-thv(kmax) )/thv(kmax)
  endif

    ql2 = 0.0
    qi2 = 0.0
    qt  = qv2

    cape = 0.0
    cin  = 0.0
    lfc  = 0.0

    doit = .true.
    cloud = .false.
    if(adiabat.eq.1.or.adiabat.eq.2)then
      ice = .false.
    else
      ice = .true.
    endif

      the = getthe(p2,t2,t2,qv2)
      if(debug_level.ge.100) print *,'  the = ',the

!---- begin ascent of parcel ----!

      if(debug_level.ge.100)then
        print *,'  Start loop:'
        print *,'  p2,th2,qv2 = ',p2,th2,qv2
      endif

    do while( doit .and. (k.lt.nk) )

        k = k+1
       b1 =  b2

       dp = p(k-1)-p(k)

      if( dp.lt.pinc )then
        nloop = 1
      else
        nloop = 1 + int( dp/pinc )
        dp = dp/float(nloop)
      endif

      do n=1,nloop

         p1 =  p2
         t1 =  t2
        pi1 = pi2
        th1 = th2
        qv1 = qv2
        ql1 = ql2
        qi1 = qi2
        thv1 = thv2

        p2 = p2 - dp
        pi2 = (p2*rp00)**rddcp

        thlast = th1
        i = 0
        not_converged = .true.

        do while( not_converged )
          i = i + 1
          t2 = thlast*pi2
          if(ice)then
            fliq = max(min((t2-233.15)/(273.15-233.15),1.0),0.0)
            fice = 1.0-fliq
          else
            fliq = 1.0
            fice = 0.0
          endif
          qv2 = min( qt , fliq*getqvs(p2,t2) + fice*getqvi(p2,t2) )
          qi2 = max( fice*(qt-qv2) , 0.0 )
          ql2 = max( qt-qv2-qi2 , 0.0 )

          tbar  = 0.5*(t1+t2)
          qvbar = 0.5*(qv1+qv2)
          qlbar = 0.5*(ql1+ql2)
          qibar = 0.5*(qi1+qi2)

          lhv = lv1-lv2*tbar
          lhs = ls1-ls2*tbar
          lhf = lhs-lhv

          rm=rd+rv*qvbar
          cpm=cp+cpv*qvbar+cpl*qlbar+cpi*qibar
          th2=th1*exp(  lhv*(ql2-ql1)/(cpm*tbar)     &
                       +lhs*(qi2-qi1)/(cpm*tbar)     &
                       +(rm/cpm-rd/cp)*alog(p2/p1) )

          if(i .gt. 90 .and. debug_level .gt. 0) print *,i,th2,thlast,th2-thlast
          if(i .gt. 100)then
!            print *
!            print *,'  Error:  lack of convergence'
!            print *
!            print *,'  ... stopping iteration '
!            print *
!            stop 1001
            if (debug_level .gt. 0) then
               write(0,*) 'getcape: lack of convergence'
            end if
            return
          endif
          if( abs(th2-thlast).gt.converge )then
            thlast=thlast+0.3*(th2-thlast)
          else
            not_converged = .false.
          endif
        enddo

        ! Latest pressure increment is complete.  Calculate some
        ! important stuff:

        if( ql2.ge.1.0e-10 ) cloud = .true.

        IF(adiabat.eq.1.or.adiabat.eq.3)THEN
          ! pseudoadiabat
          qt  = qv2
          ql2 = 0.0
          qi2 = 0.0
        ELSEIF(adiabat.le.0.or.adiabat.ge.5)THEN
!          print *
!          print *,'  Undefined adiabat'
!          print *
!          stop 10000
          write(0,*) 'getcape: Undefined adiabat'
          return
        ENDIF

      enddo

      thv2 = th2*(1.0+reps*qv2)/(1.0+qv2+ql2+qi2)
        b2 = g*( thv2-thv(k) )/thv(k)
        dz = -cpdg*0.5*(thv(k)+thv(k-1))*(pi(k)-pi(k-1))

      the = getthe(p2,t2,t2,qv2)

      ! Get contributions to CAPE and CIN:

      if( (b2.ge.0.0) .and. (b1.lt.0.0) )then
        ! first trip into positive area
        ps = p(k-1)+(p(k)-p(k-1))*(0.0-b1)/(b2-b1)
        frac = b2/(b2-b1)
        parea =  0.5*b2*dz*frac
        narea = narea-0.5*b1*dz*(1.0-frac)
        if(debug_level.ge.200)then
          print *,'      b1,b2 = ',b1,b2
          print *,'      p1,ps,p2 = ',p(k-1),ps,p(k)
          print *,'      frac = ',frac
          print *,'      parea = ',parea
          print *,'      narea = ',narea
        endif
        cin  = cin  + narea
        narea = 0.0
      elseif( (b2.lt.0.0) .and. (b1.gt.0.0) )then
        ! first trip into neg area
        ps = p(k-1)+(p(k)-p(k-1))*(0.0-b1)/(b2-b1)
        frac = b1/(b1-b2)
        parea =  0.5*b1*dz*frac
        narea = -0.5*b2*dz*(1.0-frac)
        if(debug_level.ge.200)then
          print *,'      b1,b2 = ',b1,b2
          print *,'      p1,ps,p2 = ',p(k-1),ps,p(k)
          print *,'      frac = ',frac
          print *,'      parea = ',parea
          print *,'      narea = ',narea
        endif
      elseif( b2.lt.0.0 )then
        ! still collecting negative buoyancy
        parea =  0.0
        narea = narea-0.5*dz*(b1+b2)
      else
        ! still collecting positive buoyancy
        parea =  0.5*dz*(b1+b2)
        narea =  0.0
      endif

      cape = cape + max(0.0,parea)

      if(debug_level.ge.200)then
        write(6,102) p2,b1,b2,cape,cin,cloud
102     format(5(f13.4),2x,l1)
      endif

      if( (p(k).le.10000.0).and.(b2.lt.0.0) )then
        ! stop if b < 0 and p < 100 mb
        doit = .false.
      endif

    enddo

!---- All done ----!

    return
    end subroutine getcape

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

    real (kind=RKIND) function getqvs(p,t)
    implicit none

    real (kind=RKIND) :: p,t,es

    real (kind=RKIND), parameter :: eps = 287.04/461.5

    es = 611.2*exp(17.67*(t-273.15)/(t-29.65))
    getqvs = eps*es/(p-es)

    return
    end function getqvs

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

    real (kind=RKIND) function getqvi(p,t)
    implicit none

    real (kind=RKIND) :: p,t,es

    real (kind=RKIND), parameter :: eps = 287.04/461.5

    es = 611.2*exp(21.8745584*(t-273.15)/(t-7.66))
    getqvi = eps*es/(p-es)

    return
    end function getqvi

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

    real (kind=RKIND) function getthe(p,t,td,q)
    implicit none

    real (kind=RKIND) :: p,t,td,q
    real (kind=RKIND) :: tlcl

    if( (td-t).ge.-0.1 )then
      tlcl = t
    else
      tlcl = 56.0 + ( (td-56.0)**(-1) + 0.00125*alog(t/td) )**(-1)
    endif

    getthe=t*( (100000.0/p)**(0.2854*(1.0-0.28*q)) )   &
            *exp( ((3376.0/tlcl)-2.54)*q*(1.0+0.81*q) )

    return
    end function getthe

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

!==================================================================================================
 end module mpas_atm_interp_diagnostics
!==================================================================================================
