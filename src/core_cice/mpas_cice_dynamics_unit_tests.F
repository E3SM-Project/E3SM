module cice_dynamics_unit_tests

  use mpas_grid_types

  implicit none

  private
  public :: &
       strain_rate_operator_unit_test, &
       stress_divergence_operator_unit_test

contains

  !--------------------------------------------------------------------------
  ! Spherical strain rate unit test 
  !--------------------------------------------------------------------------

  subroutine strain_rate_operator_unit_test(block, unit_test_subtype)

    use cice_dyn_hex_weak, only: strain_tensor_hex_weak

    use cice_diagnostics, only: plot_pop_grid_cells, rms_difference

    use cice_mesh, only: latlon_from_xyz, grid_rotation_forward, grid_rotation_forward_nochange

    use cice_constants, only: deg_to_rad, rad_to_deg
    
    type(block_type), intent(inout) :: block

    character(len=*), intent(in) :: &
         unit_test_subtype

    type (MPAS_pool_type), pointer :: mesh
    type (MPAS_pool_type), pointer :: normal
    type (MPAS_pool_type), pointer :: hexdyn
    type (MPAS_pool_type), pointer :: hexweak
    
    real(kind=RKIND), dimension(:), allocatable :: &
         strain11_test, &
         strain22_test, &
         strain12_test, &
         strain11_diff, &
         strain22_diff, &
         strain12_diff, &
         latVertexRotated0, &
         lonVertexRotated0, &
         latCellRotated0, &
         lonCellRotated0

    integer :: &
         iCell, &
         iCellOnCell, &
         iVertex

    real(kind=RKIND) :: &
         up, vp, ut, vt, &
         xp, yp, zp, &
         rms_strain11, &
         rms_strain22, &
         rms_strain12

    integer, pointer :: &
         nCells, &
         nVertices

    integer, dimension(:), pointer :: &
         solveStress

    real(kind=RKIND), pointer :: &
         sphere_radius

    real(kind=RKIND), dimension(:), pointer :: &
         xVertex, &
         yVertex, &
         zVertex, &
         xCell, &
         yCell, &
         zCell, &
         latCell, &
         lonCell, &
         uVelocity, &
         vVelocity, &
         latCellRotated, &
         strain11, &
         strain22, &
         strain12
  
    real(kind=RKIND), dimension(:,:,:), pointer :: &
         normalVectorPolygon

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "normal", normal)
    call MPAS_pool_get_subpool(block % structs, "hexdyn", hexdyn)
    call MPAS_pool_get_subpool(block % structs, "hexweak", hexweak)

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)
    call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)

    call MPAS_pool_get_array(mesh, "sphere_radius", sphere_radius)
    call MPAS_pool_get_array(mesh, "xVertex", xVertex)
    call MPAS_pool_get_array(mesh, "yVertex", yVertex)
    call MPAS_pool_get_array(mesh, "zVertex", zVertex)
    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)
    call MPAS_pool_get_array(mesh, "zCell", zCell)
    call MPAS_pool_get_array(mesh, "latCell", latCell)
    call MPAS_pool_get_array(mesh, "lonCell", lonCell)

    call MPAS_pool_get_array(hexdyn, "solveStress", solveStress)
    call MPAS_pool_get_array(hexdyn, "uVelocity", uVelocity)
    call MPAS_pool_get_array(hexdyn, "vVelocity", vVelocity)

    call MPAS_pool_get_array(normal, "latCellRotated", latCellRotated)
    call MPAS_pool_get_array(normal, "normalVectorPolygon", normalVectorPolygon)

    call MPAS_pool_get_array(hexweak, "strain11", strain11)
    call MPAS_pool_get_array(hexweak, "strain22", strain22)
    call MPAS_pool_get_array(hexweak, "strain12", strain12)

    open(54,file="strain_rate_test.txt")

    solveStress(:) = 1

    allocate(strain11_test(nCells))
    allocate(strain22_test(nCells))
    allocate(strain12_test(nCells))
    allocate(strain11_diff(nCells))
    allocate(strain22_diff(nCells))
    allocate(strain12_diff(nCells))

    allocate(lonVertexRotated0(nVertices))
    allocate(latVertexRotated0(nVertices))
    allocate(lonCellRotated0(nCells))
    allocate(latCellRotated0(nCells))

    do iVertex = 1, nVertices
    
       call grid_rotation_forward(&
            xp,               yp,               zp, &
            xVertex(iVertex), yVertex(iVertex), zVertex(iVertex))

       call latlon_from_xyz(latVertexRotated0(iVertex), lonVertexRotated0(iVertex), &
            xp, yp, zp, sphere_radius)

    enddo ! iVertex

    do iCell = 1, nCells

       call grid_rotation_forward(&
            xp,           yp,           zp, &
            xCell(iCell), yCell(iCell), zCell(iCell))

       call latlon_from_xyz(latCellRotated0(iCell), lonCellRotated0(iCell), &
            xp, yp, zp, sphere_radius)

    enddo ! iCell
    
    ! set test velocities and expected outputs
    call spherical_test_strain(mesh, &
         uVelocity,              &
         vVelocity,              &
         lonVertexRotated0,      &
         latVertexRotated0,      &
         strain11_test,          &
         strain22_test,          &
         strain12_test,          &
         lonCellRotated0,        &
         latCellRotated0,        &
         trim(unit_test_subtype))

    write(54,*) 
    write(54,*) "Analytic solution:"
    write(54,*) "strain11 (min/max):", minval(strain11_test), maxval(strain11_test)
    write(54,*) "strain22 (min/max):", minval(strain22_test), maxval(strain22_test)
    write(54,*) "strain12 (min/max):", minval(strain12_test), maxval(strain12_test)

    ! calculate strain
    call strain_tensor_hex_weak(mesh, &
         strain11,            &
         strain22,            &
         strain12,            &
         uVelocity,           &
         vVelocity,           &
         normalVectorPolygon, &
         latCellRotated,      &
         solveStress)

    write(54,*) 
    write(54,*) "Numerical solution:"
    write(54,*) "strain11 (min/max):", minval(strain11(1:nCells)), &
                                       maxval(strain11(1:nCells))
    write(54,*) "strain22 (min/max):", minval(strain22(1:nCells)), &
                                       maxval(strain22(1:nCells))
    write(54,*) "strain12 (min/max):", minval(strain12(1:nCells)), &
                                       maxval(strain12(1:nCells))

    ! calculate difference
    strain11_diff = strain11(1:nCells) - strain11_test
    strain22_diff = strain22(1:nCells) - strain22_test
    strain12_diff = strain12(1:nCells) - strain12_test

    write(54,*) 
    write(54,*) "Difference:"
    write(54,*) "strain11 (min/max):", minval(strain11_diff), maxval(strain11_diff)
    write(54,*) "strain22 (min/max):", minval(strain22_diff), maxval(strain22_diff)
    write(54,*) "strain12 (min/max):", minval(strain12_diff), maxval(strain12_diff)

    open(55,file="strains.txt")
    do iCell = 1, nCells

       write(55,*) iCell, solveStress(iCell), &
            strain11(iCell), strain11_test(iCell), &
            strain22(iCell), strain22_test(iCell), &
            strain12(iCell), strain12_test(iCell), &
            strain12(iCell) - strain12_test(iCell)

    enddo ! iCell
    close(55)

    ! set cells to solve
    do iCell = 1, nCells
       solveStress(iCell) = 1
       ! cut off low latitude cells
       if (latCell(iCell) < 30.0 * deg_to_rad .and. latCell(iCell) > -30.0 * deg_to_rad) then
          solveStress(iCell) = 0
       endif
    enddo ! iCell

    open(55,file='strain11.txt')
    open(56,file='strain22.txt')
    open(57,file='strain12.txt')
    do iCell = 1, nCells
       write(55,*) iCell, lonCell(iCell) * rad_to_deg, latCell(iCell) * rad_to_deg, &
            strain11(iCell), strain11_test(iCell), strain11_diff(iCell), &
            strain11_diff(iCell) * real(solveStress(iCell),RKIND)
       write(56,*) iCell, lonCell(iCell) * rad_to_deg, latCell(iCell) * rad_to_deg, &
            strain22(iCell), strain22_test(iCell), strain22_diff(iCell), &
            strain22_diff(iCell) * real(solveStress(iCell),RKIND)
       write(57,*) iCell, lonCell(iCell) * rad_to_deg, latCell(iCell) * rad_to_deg, &
            strain12(iCell), strain12_test(iCell), strain12_diff(iCell), &
            strain12_diff(iCell) * real(solveStress(iCell),RKIND)
    enddo ! iCell
    close(55)
    close(56)
    close(57)

    ! rms difference
    call rms_difference( &
         strain11(1:nCells),    &
         strain11_test,         &
         solveStress(1:nCells), &
         nCells,                &
         rms_strain11)

    call rms_difference( &
         strain22(1:nCells),    &
         strain22_test,         &
         solveStress(1:nCells), &
         nCells,                &
         rms_strain22)

    call rms_difference( &
         strain12(1:nCells),    &
         strain12_test,         &
         solveStress(1:nCells), &
         nCells,                &
         rms_strain12)

    write(54,*)
    write(54,*) "RMS:"
    write(54,*) "strain11:", rms_strain11
    write(54,*) "strain22:", rms_strain22
    write(54,*) "strain12:", rms_strain12

    deallocate(strain11_test)
    deallocate(strain22_test)
    deallocate(strain12_test)
    deallocate(strain11_diff)
    deallocate(strain22_diff)
    deallocate(strain12_diff)

    deallocate(lonVertexRotated0)
    deallocate(latVertexRotated0)
    deallocate(lonCellRotated0)
    deallocate(latCellRotated0)

    close(54)

  end subroutine strain_rate_operator_unit_test

  !-------------------------------------------------------------

  subroutine spherical_test_strain(mesh, &
       uVelocity,         &
       vVelocity,         &
       longitudeVelocity, &
       latitudeVelocity,  &
       strain11,          &
       strain22,          &
       strain12,          &
       longitudeStrain,   &
       latitudeStrain,    &
       test_type)

    type(MPAS_pool_type), pointer, intent(in) :: mesh

    real(kind=RKIND), dimension(:), intent(out) :: &
         uVelocity, &
         vVelocity, &
         strain11,  &
         strain22,  &
         strain12

    real(kind=RKIND), dimension(:), intent(in) :: &
         longitudeVelocity, &
         latitudeVelocity,  &
         longitudeStrain,   &
         latitudeStrain

    character(len=*), intent(in) :: &
         test_type

    real(kind=RKIND) :: &
         uVelocity_test, &
         vVelocity_test, &
         du_dlon, &
         du_dlat, &
         dv_dlon, &
         dv_dlat, &
         radius

    integer :: &
         nPointsVelocity, &
         nPointsStrain, &
         iPoint

    real(kind=RKIND), pointer :: &
         sphere_radius

    call MPAS_pool_get_array(mesh, "sphere_radius", sphere_radius)

    nPointsVelocity = size(uVelocity)
    nPointsStrain   = size(strain11)

    radius = sphere_radius

    ! set velocity points
    do iPoint = 1, nPointsVelocity

       call spherical_test_strain_velocities( &
            uVelocity(iPoint),         &
            vVelocity(iPoint),         &
            du_dlon,                   &
            du_dlat,                   &
            dv_dlon,                   &
            dv_dlat,                   &
            longitudeVelocity(iPoint), &
            latitudeVelocity(iPoint),  &
            test_type)

    enddo ! iPoint

    ! set strain points
    do iPoint = 1, nPointsStrain

       call spherical_test_strain_velocities( &
            uVelocity_test,          &
            vVelocity_test,          &
            du_dlon,                 &
            du_dlat,                 &
            dv_dlon,                 &
            dv_dlat,                 &
            longitudeStrain(iPoint), &
            latitudeStrain(iPoint),  &
            test_type)

       strain11(iPoint) = strain11_component(uVelocity_test, vVelocity_test, &
                                             du_dlon, du_dlat, dv_dlon, dv_dlat, &
                                             radius, longitudeStrain(iPoint), latitudeStrain(iPoint))
       strain22(iPoint) = strain22_component(uVelocity_test, vVelocity_test, &
                                             du_dlon, du_dlat, dv_dlon, dv_dlat, &
                                             radius, longitudeStrain(iPoint), latitudeStrain(iPoint))
       strain12(iPoint) = strain12_component(uVelocity_test, vVelocity_test, &
                                             du_dlon, du_dlat, dv_dlon, dv_dlat, &
                                             radius, longitudeStrain(iPoint), latitudeStrain(iPoint))

    enddo ! iPoint

  end subroutine spherical_test_strain

  !-------------------------------------------------------------
  
  subroutine spherical_test_strain_velocities( &
       u,       &
       v,       &
       du_dlon, &
       du_dlat, &
       dv_dlon, &
       dv_dlat, &
       lon,     &
       lat,     &
       test_type)

    real(kind=RKIND), intent(out) :: &
         u,       &
         v,       &
         du_dlon, &
         du_dlat, &
         dv_dlon, &
         dv_dlat
    
    real(kind=RKIND), intent(in) :: &
         lon, &
         lat
    
    character(len=*), intent(in) :: &
         test_type

    if (trim(test_type) == "zero") then

       u = 0.0_RKIND
       v = 0.0_RKIND

       du_dlon = 0.0_RKIND
       du_dlat = 0.0_RKIND
       dv_dlon = 0.0_RKIND
       dv_dlat = 0.0_RKIND

    else if (trim(test_type) == "zonal") then

       u = 1.0_RKIND
       v = 0.0_RKIND

       du_dlon = 0.0_RKIND
       du_dlat = 0.0_RKIND
       dv_dlon = 0.0_RKIND
       dv_dlat = 0.0_RKIND

    else if (trim(test_type) == "meridonal") then

       u = 0.0_RKIND
       v = 1.0_RKIND

       du_dlon = 0.0_RKIND
       du_dlat = 0.0_RKIND
       dv_dlon = 0.0_RKIND
       dv_dlat = 0.0_RKIND

    else if (trim(test_type) == "solid_body") then

       u = cos(lat)
       v = 0.0_RKIND

       du_dlon = 0.0_RKIND
       du_dlat = -sin(lat)
       dv_dlon = 0.0_RKIND
       dv_dlat = 0.0_RKIND

    else if (trim(test_type) == "sinusoidal1") then

       u = cos(lon) * (1.0_RKIND + cos(2.0_RKIND * lat))
       v = 0.0_RKIND

       du_dlon = -sin(lon) * (1.0_RKIND + cos(2.0_RKIND * lat))
       du_dlat = -2.0_RKIND * cos(lon) * sin(2.0_RKIND * lat)
       dv_dlon = 0.0_RKIND
       dv_dlat = 0.0_RKIND

    else if (trim(test_type) == "sinusoidal2") then

       u = 0.0_RKIND
       v = cos(lon) * (1.0_RKIND + cos(2.0_RKIND * lat))

       du_dlon = 0.0_RKIND
       du_dlat = 0.0_RKIND
       dv_dlon = -sin(lon) * (1.0_RKIND + cos(2.0_RKIND * lat))
       dv_dlat = -2.0_RKIND * cos(lon) * sin(2.0_RKIND * lat)

    endif

  end subroutine spherical_test_strain_velocities

  !-------------------------------------------------------------

  function strain11_component(u, v, du_dlon, du_dlat, dv_dlon, dv_dlat, r, lon, lat) result(e11)

    real(kind=RKIND), intent(in) :: &
         u,       &
         v,       &
         du_dlon, &
         du_dlat, &
         dv_dlon, &
         dv_dlat, &                                  
         r,       &
         lon,     &
         lat

    real(kind=RKIND) :: e11

    e11 = (1.0_RKIND / (r * cos(lat))) * (du_dlon - v * sin(lat)) 

  end function strain11_component

  !-------------------------------------------------------------

  function strain22_component(u, v, du_dlon, du_dlat, dv_dlon, dv_dlat, r, lon, lat) result(e22)

    real(kind=RKIND), intent(in) :: &
         u,       &
         v,       &
         du_dlon, &
         du_dlat, &
         dv_dlon, &
         dv_dlat, &                                  
         r,       &
         lon,     &
         lat

    real(kind=RKIND) :: e22

    e22 = dv_dlat / r

  end function strain22_component

  !-------------------------------------------------------------

  function strain12_component(u, v, du_dlon, du_dlat, dv_dlon, dv_dlat, r, lon, lat) result(e12)

    real(kind=RKIND), intent(in) :: &
         u,       &
         v,       &
         du_dlon, &
         du_dlat, &
         dv_dlon, &
         dv_dlat, &                                  
         r,       &
         lon,     &
         lat

    real(kind=RKIND) :: e12

    e12 = (0.5_RKIND / r) * (du_dlat + u * tan(lat) + dv_dlon / cos(lat))

  end function strain12_component

  !--------------------------------------------------------------------------
  ! Spherical stress divergence unit test
  !--------------------------------------------------------------------------

  subroutine stress_divergence_operator_unit_test(block, unit_test_subtype)

    use cice_dyn_hex_weak, only: stress_divergence_hex_weak

    use cice_diagnostics, only: rms_difference

    use cice_mesh, only: grid_rotation_forward, latlon_from_xyz

    use cice_constants, only: deg_to_rad, rad_to_deg

    type(block_type), intent(inout) :: block

    character(len=*), intent(in) :: &
         unit_test_subtype

    type (MPAS_pool_type), pointer :: mesh
    type (MPAS_pool_type), pointer :: normal
    type (MPAS_pool_type), pointer :: hexdyn
    type (MPAS_pool_type), pointer :: hexweak

    real(kind=RKIND), dimension(:), allocatable :: &
         stressDivergenceU_test, &
         stressDivergenceV_test, &
         stressDivergenceU_diff, &
         stressDivergenceV_diff, &
         latVertexRotated0, &
         lonVertexRotated0, &
         latCellRotated0, &
         lonCellRotated0
    
    real(kind=RKIND) :: &
         xp, yp, zp, &
         rms_divergenceU, &
         rms_divergenceV

    integer :: &
         iCell, &
         iVertex

    integer, pointer :: &
         nCells, &
         nVertices

    integer, dimension(:), pointer :: &
         solveVelocity

    real(kind=RKIND), pointer :: &
         sphere_radius

    real(kind=RKIND), dimension(:), pointer :: &
         xVertex, &
         yVertex, &
         zVertex, &
         xCell, &
         yCell, &
         zCell, &
         latVertex, &
         lonVertex, &
         latVertexRotated, &
         stress11, &
         stress22, &
         stress12, &
         stressDivergenceU, &
         stressDivergenceV

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         normalVectorTriangle

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "normal", normal)
    call MPAS_pool_get_subpool(block % structs, "hexdyn", hexdyn)
    call MPAS_pool_get_subpool(block % structs, "hexweak", hexweak)

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)
    call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)

    call MPAS_pool_get_array(mesh, "sphere_radius", sphere_radius)
    call MPAS_pool_get_array(mesh, "xVertex", xVertex)
    call MPAS_pool_get_array(mesh, "yVertex", yVertex)
    call MPAS_pool_get_array(mesh, "zVertex", zVertex)
    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)
    call MPAS_pool_get_array(mesh, "zCell", zCell)
    call MPAS_pool_get_array(mesh, "latVertex", latVertex)
    call MPAS_pool_get_array(mesh, "lonVertex", lonVertex)

    call MPAS_pool_get_array(hexdyn, "solveVelocity", solveVelocity)
    call MPAS_pool_get_array(hexdyn, "stressDivergenceU", stressDivergenceU)
    call MPAS_pool_get_array(hexdyn, "stressDivergenceV", stressDivergenceV)

    call MPAS_pool_get_array(normal, "normalVectorTriangle", normalVectorTriangle)
    call MPAS_pool_get_array(normal, "latVertexRotated", latVertexRotated)

    call MPAS_pool_get_array(hexweak, "stress11", stress11)
    call MPAS_pool_get_array(hexweak, "stress22", stress22)
    call MPAS_pool_get_array(hexweak, "stress12", stress12)

    open(54,file="stress_divergence_test.txt")

    solveVelocity(:) = 1

    allocate(stressDivergenceU_test(nVertices))
    allocate(stressDivergenceV_test(nVertices))
    allocate(stressDivergenceU_diff(nVertices))
    allocate(stressDivergenceV_diff(nVertices))

    allocate(lonVertexRotated0(nVertices))
    allocate(latVertexRotated0(nVertices))
    allocate(lonCellRotated0(nCells))
    allocate(latCellRotated0(nCells))

    do iVertex = 1, nVertices
    
       call grid_rotation_forward(&
            xp,               yp,               zp, &
            xVertex(iVertex), yVertex(iVertex), zVertex(iVertex))

       call latlon_from_xyz(latVertexRotated0(iVertex), lonVertexRotated0(iVertex), &
            xp, yp, zp, sphere_radius)

    enddo ! iVertex

    do iCell = 1, nCells

       call grid_rotation_forward(&
            xp,           yp,           zp, &
            xCell(iCell), yCell(iCell), zCell(iCell))

       call latlon_from_xyz(latCellRotated0(iCell), lonCellRotated0(iCell), &
            xp, yp, zp, sphere_radius)

    enddo ! iCell

    ! set test velocities and expected outputs
    call spherical_test_divergence_stress(mesh, &
         stress11,               &
         stress22,               &
         stress12,               &
         lonCellRotated0,        &
         latCellRotated0,        &
         stressDivergenceU_test, &
         stressDivergenceV_test, &
         lonVertexRotated0,      &
         latVertexRotated0,      &
         trim(unit_test_subtype))

    write(54,*)
    write(54,*) "Analytic solution:"
    write(54,*) "stressDivergenceU:", minval(stressDivergenceU_test), maxval(stressDivergenceU_test)
    write(54,*) "stressDivergenceV:", minval(stressDivergenceV_test), maxval(stressDivergenceV_test)

    ! calculate divergence of stress
    call stress_divergence_hex_weak(mesh, &
         stressDivergenceU,    &
         stressDivergenceV,    &
         stress11,             &
         stress22,             &
         stress12,             &
         normalVectorTriangle, &
         latVertexRotated,     &
         solveVelocity)

    write(54,*)
    write(54,*) "Numerical solution:"
    write(54,*) "stressDivergenceU:", minval(stressDivergenceU), &
                                      maxval(stressDivergenceU)
    write(54,*) "stressDivergenceV:", minval(stressDivergenceV), &
                                      maxval(stressDivergenceV)

    ! calculate difference
    stressDivergenceU_diff = stressDivergenceU(1:nVertices) - &
                             stressDivergenceU_test
    stressDivergenceV_diff = stressDivergenceV(1:nVertices) - &
                             stressDivergenceV_test

    write(54,*)
    write(54,*) "Difference:"
    write(54,*) "stressDivergenceU:", minval(stressDivergenceU_diff), maxval(stressDivergenceU_diff)
    write(54,*) "stressDivergenceV:", minval(stressDivergenceV_diff), maxval(stressDivergenceV_diff)

    open(55, file="divergences")
    do iVertex = 1, nVertices

       write(55,*) iVertex, solveVelocity(iVertex), &
            stressDivergenceU(iVertex), stressDivergenceU_test(iVertex), &
            stressDivergenceV(iVertex), stressDivergenceV_test(iVertex), &
            stressDivergenceU(iVertex) - stressDivergenceU_test(iVertex), &
            stressDivergenceV(iVertex) - stressDivergenceV_test(iVertex)

    enddo ! iVertex
    close(55)

    ! set vertices to solve
    do iVertex = 1, nVertices
       solveVelocity(iVertex) = 1
       ! cut off low latitude cells
       if (latVertex(iVertex) < 30.0 * deg_to_rad .and. latVertex(iVertex) > -30.0 * deg_to_rad) then
          solveVelocity(iVertex) = 0
       endif
    enddo ! iVertex

    open(55,file='divergenceu.txt')
    open(56,file='divergencev.txt')
    do iVertex = 1, nVertices
       write(55,*) iVertex, lonVertex(iVertex) * rad_to_deg, latVertex(iVertex) * rad_to_deg, &
            stressDivergenceU(iVertex), stressDivergenceU_test(iVertex), stressDivergenceU_diff(iVertex), &
            stressDivergenceU_diff(iVertex) * real(solveVelocity(iVertex),RKIND)
       write(56,*) iVertex, lonVertex(iVertex) * rad_to_deg, latVertex(iVertex) * rad_to_deg, &
            stressDivergenceV(iVertex), stressDivergenceV_test(iVertex), stressDivergenceV_diff(iVertex), &
            stressDivergenceV_diff(iVertex) * real(solveVelocity(iVertex),RKIND)
    enddo ! iVertex
    close(55)
    close(56)

    ! rms difference
    call rms_difference( &
         stressDivergenceU(1:nVertices), &
         stressDivergenceU_test,         &
         solveVelocity(1:nVertices),     &
         nVertices,                      &
         rms_divergenceU)

    call rms_difference( &
         stressDivergenceV(1:nVertices), &
         stressDivergenceV_test,         &
         solveVelocity(1:nVertices),     &
         nVertices,                      &
         rms_divergenceV)

    write(54,*) 
    write(54,*) "RMS:"
    write(54,*) "rms_divergenceU:", rms_divergenceU
    write(54,*) "rms_divergenceV:", rms_divergenceV

    ! cleanup
    deallocate(stressDivergenceU_test)
    deallocate(stressDivergenceV_test)
    deallocate(stressDivergenceU_diff)
    deallocate(stressDivergenceV_diff)

    deallocate(lonVertexRotated0)
    deallocate(latVertexRotated0)
    deallocate(lonCellRotated0)
    deallocate(latCellRotated0)

    close(54)

  end subroutine stress_divergence_operator_unit_test

  !--------------------------------------------------------------------------

  subroutine spherical_test_divergence_stress(mesh, &
       stress11,            &
       stress22,            &
       stress12,            &
       longitudeStress,     &
       latitudeStress,      &
       stressDivergenceU,   &
       stressDivergenceV,   &
       longitudeDivergence, &
       latitudeDivergence,  &
       test_type)

    type(MPAS_pool_type), pointer, intent(in) :: mesh

    real(kind=RKIND), dimension(:), intent(out) :: &
         stressDivergenceU, &
         stressDivergenceV, &
         stress11,          &
         stress22,          &
         stress12
    
    real(kind=RKIND), dimension(:), intent(in) :: &
         longitudeStress,     &
         latitudeStress,      &
         longitudeDivergence, &
         latitudeDivergence

    character(len=*), intent(in) :: &
         test_type

    real(kind=RKIND) :: &
         stress11_test,  &
         stress22_test,  &
         stress12_test,  &
         dstress11_dlon, &
         dstress11_dlat, &
         dstress22_dlon, &
         dstress22_dlat, &
         dstress12_dlon, &
         dstress12_dlat, &
         radius

    integer :: &
         nPointsStress, &
         nPointsDivergence, &
         iPoint

    real(kind=RKIND), pointer :: &
         sphere_radius
    
    call MPAS_pool_get_array(mesh, "sphere_radius", sphere_radius)

    radius = sphere_radius

    nPointsStress     = size(stress11)
    nPointsDivergence = size(stressDivergenceU)

    ! set stress points
    do iPoint = 1, nPointsStress

       call spherical_test_divergence_stress_stresses( &
         stress11(iPoint),        &
         stress22(iPoint),        &
         stress12(iPoint),        &
         dstress11_dlon,          &
         dstress11_dlat,          &
         dstress22_dlon,          &
         dstress22_dlat,          &
         dstress12_dlon,          &
         dstress12_dlat,          &
         longitudeStress(iPoint), &
         latitudeStress(iPoint),  &
         test_type)

    enddo ! iPoint

    ! set divergence poinys
    do iPoint = 1, nPointsDivergence

       call spherical_test_divergence_stress_stresses( &
         stress11_test,               &
         stress22_test,               &
         stress12_test,               &
         dstress11_dlon,              &
         dstress11_dlat,              &
         dstress22_dlon,              &
         dstress22_dlat,              &
         dstress12_dlon,              &
         dstress12_dlat,              &
         longitudeDivergence(iPoint), &
         latitudeDivergence(iPoint),  &
         test_type)

       stressDivergenceU(iPoint) = divergenceStressU(stress11_test, stress22_test, stress12_test, &
                                                     dstress11_dlon, dstress11_dlat, dstress22_dlon, dstress22_dlat, dstress12_dlon, dstress12_dlat, &
                                                     radius, longitudeDivergence(iPoint), latitudeDivergence(iPoint))
       stressDivergenceV(iPoint) = divergenceStressV(stress11_test, stress22_test, stress12_test, &
                                                     dstress11_dlon, dstress11_dlat, dstress22_dlon, dstress22_dlat, dstress12_dlon, dstress12_dlat, &
                                                     radius, longitudeDivergence(iPoint), latitudeDivergence(iPoint))

       write(55,*) iPoint, longitudeDivergence(iPoint), latitudeDivergence(iPoint), stressDivergenceU(iPoint), stressDivergenceV(iPoint), radius

    enddo ! iPoint

  end subroutine spherical_test_divergence_stress

  !-------------------------------------------------------------

  subroutine spherical_test_divergence_stress_stresses( &
       stress11,       &
       stress22,       &
       stress12,       &
       dstress11_dlon, &
       dstress11_dlat, &
       dstress22_dlon, &
       dstress22_dlat, &
       dstress12_dlon, &
       dstress12_dlat, &
       lon,            &
       lat,            &
       test_type)
    
    real(kind=RKIND), intent(out) :: &
         stress11,       &
         stress22,       &
         stress12,       &
         dstress11_dlon, &
         dstress11_dlat, &
         dstress22_dlon, &
         dstress22_dlat, &
         dstress12_dlon, &
         dstress12_dlat

    real(kind=RKIND), intent(in) :: &
         lon, &
         lat

    character(len=*), intent(in) :: &
         test_type

    if (trim(test_type) == "zero") then
       
       stress11 = 0.0_RKIND
       stress22 = 0.0_RKIND
       stress12 = 0.0_RKIND
       
       dstress11_dlon = 0.0_RKIND
       dstress11_dlat = 0.0_RKIND
       dstress22_dlon = 0.0_RKIND
       dstress22_dlat = 0.0_RKIND
       dstress12_dlon = 0.0_RKIND
       dstress12_dlat = 0.0_RKIND

    else if (trim(test_type) == "const11") then

       stress11 = 1.0_RKIND
       stress22 = 0.0_RKIND
       stress12 = 0.0_RKIND
       
       dstress11_dlon = 0.0_RKIND
       dstress11_dlat = 0.0_RKIND
       dstress22_dlon = 0.0_RKIND
       dstress22_dlat = 0.0_RKIND
       dstress12_dlon = 0.0_RKIND
       dstress12_dlat = 0.0_RKIND

    else if (trim(test_type) == "const22") then

       stress11 = 0.0_RKIND
       stress22 = 1.0_RKIND
       stress12 = 0.0_RKIND
       
       dstress11_dlon = 0.0_RKIND
       dstress11_dlat = 0.0_RKIND
       dstress22_dlon = 0.0_RKIND
       dstress22_dlat = 0.0_RKIND
       dstress12_dlon = 0.0_RKIND
       dstress12_dlat = 0.0_RKIND

    else if (trim(test_type) == "const12") then

       stress11 = 0.0_RKIND
       stress22 = 0.0_RKIND
       stress12 = 1.0_RKIND
       
       dstress11_dlon = 0.0_RKIND
       dstress11_dlat = 0.0_RKIND
       dstress22_dlon = 0.0_RKIND
       dstress22_dlat = 0.0_RKIND
       dstress12_dlon = 0.0_RKIND
       dstress12_dlat = 0.0_RKIND

    else if (trim(test_type) == "test1") then

       stress11 = cos(lon) * (1.0_RKIND + cos(2.0_RKIND * lat))
       stress22 = 0.0_RKIND
       stress12 = 0.0_RKIND
       
       dstress11_dlon = -sin(lon) * (1.0_RKIND + cos(2.0_RKIND * lat))
       dstress11_dlat = -cos(lon) * 2.0_RKIND * sin(2.0_RKIND * lat)
       dstress22_dlon = 0.0_RKIND
       dstress22_dlat = 0.0_RKIND
       dstress12_dlon = 0.0_RKIND
       dstress12_dlat = 0.0_RKIND

    else if (trim(test_type) == "test2") then

       stress11 = 0.0_RKIND
       stress22 = cos(lon) * (1.0_RKIND + cos(2.0_RKIND * lat))
       stress12 = 0.0_RKIND
       
       dstress11_dlon = 0.0_RKIND
       dstress11_dlat = 0.0_RKIND
       dstress22_dlon = -sin(lon) * (1.0_RKIND + cos(2.0_RKIND * lat))
       dstress22_dlat = -cos(lon) * 2.0_RKIND * sin(2.0_RKIND * lat)
       dstress12_dlon = 0.0_RKIND
       dstress12_dlat = 0.0_RKIND

    else if (trim(test_type) == "test3") then

       stress11 = 0.0_RKIND
       stress22 = 0.0_RKIND
       stress12 = cos(lon) * (1.0_RKIND + cos(2.0_RKIND * lat))
       
       dstress11_dlon = 0.0_RKIND
       dstress11_dlat = 0.0_RKIND
       dstress22_dlon = 0.0_RKIND
       dstress22_dlat = 0.0_RKIND
       dstress12_dlon = -sin(lon) * (1.0_RKIND + cos(2.0_RKIND * lat))
       dstress12_dlat = -cos(lon) * 2.0_RKIND * sin(2.0_RKIND * lat)

    else if (trim(test_type) == "test4") then

       stress11 = cos(lon) * (1.0_RKIND + cos(2.0_RKIND * lat))
       stress22 = cos(lon) * (1.0_RKIND + cos(2.0_RKIND * lat))
       stress12 = cos(lon) * (1.0_RKIND + cos(2.0_RKIND * lat))
       
       dstress11_dlon = -sin(lon) * (1.0_RKIND + cos(2.0_RKIND * lat))
       dstress11_dlat = -cos(lon) * 2.0_RKIND * sin(2.0_RKIND * lat)
       dstress22_dlon = -sin(lon) * (1.0_RKIND + cos(2.0_RKIND * lat))
       dstress22_dlat = -cos(lon) * 2.0_RKIND * sin(2.0_RKIND * lat)
       dstress12_dlon = -sin(lon) * (1.0_RKIND + cos(2.0_RKIND * lat))
       dstress12_dlat = -cos(lon) * 2.0_RKIND * sin(2.0_RKIND * lat)

    else
       
       write(*,*) "Unknown test case"
       stop

    endif

  end subroutine spherical_test_divergence_stress_stresses

  !-------------------------------------------------------------

  function divergenceStressU(stress11, stress22, stress12, &
                             dstress11_dlon, dstress11_dlat, dstress22_dlon, dstress22_dlat, dstress12_dlon, dstress12_dlat, &
                             r, lon, lat) result(divu)

    real(kind=RKIND), intent(in) :: &
         stress11,       &
         stress22,       &
         stress12,       &
         dstress11_dlon, &
         dstress11_dlat, &
         dstress22_dlon, &
         dstress22_dlat, &
         dstress12_dlon, &
         dstress12_dlat, &
         r,              &
         lon,            &
         lat

    real(kind=RKIND) :: divu

    divu = (1.0_RKIND / (r * cos(lat))) * dstress11_dlon + &
           (1.0_RKIND / r)              * dstress12_dlat - &
           (2.0_RKIND / r) * tan(lat)   * stress12

  end function divergenceStressU

  !-------------------------------------------------------------

  function divergenceStressV(stress11, stress22, stress12, &
                             dstress11_dlon, dstress11_dlat, dstress22_dlon, dstress22_dlat, dstress12_dlon, dstress12_dlat, &
                             r, lon, lat) result(divv)

    real(kind=RKIND), intent(in) :: &
         stress11,       &
         stress22,       &
         stress12,       &
         dstress11_dlon, &
         dstress11_dlat, &
         dstress22_dlon, &
         dstress22_dlat, &
         dstress12_dlon, &
         dstress12_dlat, &
         r,              &
         lon,            &
         lat

    real(kind=RKIND) :: divv

    divv = (1.0_RKIND / (r * cos(lat))) * dstress12_dlon + &
           (1.0_RKIND / r)              * dstress22_dlat + &
           (1.0_RKIND / r) * tan(lat)   * stress11       - &
           (1.0_RKIND / r) * tan(lat)   * stress22

  end function divergenceStressV

  !--------------------------------------------------------------------------

end module cice_dynamics_unit_tests



