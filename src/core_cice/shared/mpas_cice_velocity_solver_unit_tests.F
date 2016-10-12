!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_velocity_solver_unit_tests
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

module cice_velocity_solver_unit_tests

  use mpas_derived_types
  use mpas_pool_routines
  use mpas_io_units

  implicit none

  private
  save

  public :: &
       cice_strain_rate_operator_unit_test, &
       cice_stress_divergence_operator_unit_test, &
       cice_constitutive_relationship_unit_test

contains

!-----------------------------------------------------------------------
! Spherical strain rate unit test
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_strain_rate_operator_unit_test
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_strain_rate_operator_unit_test(&
       domain, &
       unit_test_subtype)!{{{

    type(domain_type), intent(inout) :: &
         domain !< Input/Output:

    character(len=*), intent(in) :: &
         unit_test_subtype !< Input:

    integer, parameter :: ntests = 6
    character(len=200), dimension(ntests), parameter :: &
         test_names = (/"zero       ", &
                        "zonal      ", &
                        "meridonal  ", &
                        "solid_body ", &
                        "sinusoidal1", &
                        "sinusoidal2"/)

    integer :: &
         itest

    if (trim(unit_test_subtype) == "all") then

       do itest = 1, ntests

          call cice_strain_rate_operator_unit_test_individual(&
               domain, &
               trim(test_names(itest)))

       enddo ! itest

    else

       call cice_strain_rate_operator_unit_test_individual(&
            domain, &
            unit_test_subtype)

    endif

  end subroutine cice_strain_rate_operator_unit_test

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_strain_rate_operator_unit_test
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_strain_rate_operator_unit_test_individual(&
       domain, &
       unit_test_subtype)!{{{

    use cice_velocity_solver_weak, only: &
         cice_strain_tensor_weak

    use cice_velocity_solver_variational, only: &
         cice_strain_tensor_variational

    use cice_debug, only: &
         cice_rms_difference, &
         cice_plot_spherical_latlon

    use cice_mesh, only: &
         cice_latlon_from_xyz, &
         cice_grid_rotation_forward

    use cice_constants, only: &
         ciceRadiansToDegrees

    use cice_constants, only: &
         pii

    type(domain_type) :: &
         domain !< Input/Output:

    character(len=*), intent(in) :: &
         unit_test_subtype !< Input:

    type(block_type), pointer :: &
         block

    type (MPAS_pool_type), pointer :: &
         mesh, &
         boundary, &
         velocity_solver, &
         velocity_weak, &
         velocity_variational

    real(kind=RKIND), dimension(:), allocatable :: &
         strain11, &
         strain22, &
         strain12, &
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
         solveStress, &
         interiorCell

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
         strain11_weak, &
         strain22_weak, &
         strain12_weak, &
         tanLatVertexRotatedOverRadius

    real(kind=RKIND), dimension(:,:), pointer :: &
         strain11_variational, &
         strain22_variational, &
         strain12_variational

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         basisGradientU, &
         basisGradientV

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         normalVectorPolygon

    logical, pointer :: &
         config_rotate_cartesian_grid

    character(len=strKIND), pointer :: &
         config_stress_divergence_scheme

    real(kind=RKIND), parameter :: &
         polarAngleLimit = 0.95_RKIND * 0.5_RKIND * pii

    call MPAS_pool_get_config(domain % configs, "config_rotate_cartesian_grid", config_rotate_cartesian_grid)
    call MPAS_pool_get_config(domain % configs, "config_stress_divergence_scheme", config_stress_divergence_scheme)

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
       call MPAS_pool_get_subpool(block % structs, "boundary", boundary)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocity_solver)
       call MPAS_pool_get_subpool(block % structs, "velocity_weak", velocity_weak)
       call MPAS_pool_get_subpool(block % structs, "velocity_variational", velocity_variational)

       call MPAS_pool_get_dimension(mesh, "nCells", nCells)
       call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)
       call MPAS_pool_get_config(mesh, "sphere_radius", sphere_radius)
       sphere_radius = 6371220.0_RKIND

       call MPAS_pool_get_array(mesh, "xVertex", xVertex)
       call MPAS_pool_get_array(mesh, "yVertex", yVertex)
       call MPAS_pool_get_array(mesh, "zVertex", zVertex)
       call MPAS_pool_get_array(mesh, "xCell", xCell)
       call MPAS_pool_get_array(mesh, "yCell", yCell)
       call MPAS_pool_get_array(mesh, "zCell", zCell)
       call MPAS_pool_get_array(mesh, "latCell", latCell)
       call MPAS_pool_get_array(mesh, "lonCell", lonCell)

       call MPAS_pool_get_array(velocity_solver, "solveStress", solveStress)
       call MPAS_pool_get_array(velocity_solver, "uVelocity", uVelocity)
       call MPAS_pool_get_array(velocity_solver, "vVelocity", vVelocity)

       call MPAS_pool_get_array(boundary, "interiorCell", interiorCell)

       if (trim(config_stress_divergence_scheme) == "weak") then
          call MPAS_pool_get_array(velocity_weak, "latCellRotated", latCellRotated)
          call MPAS_pool_get_array(velocity_weak, "normalVectorPolygon", normalVectorPolygon)
          call MPAS_pool_get_array(velocity_weak, "strain11", strain11_weak)
          call MPAS_pool_get_array(velocity_weak, "strain22", strain22_weak)
          call MPAS_pool_get_array(velocity_weak, "strain12", strain12_weak)
       else if (trim(config_stress_divergence_scheme) == "variational") then
          call MPAS_pool_get_array(velocity_variational, "strain11", strain11_variational)
          call MPAS_pool_get_array(velocity_variational, "strain22", strain22_variational)
          call MPAS_pool_get_array(velocity_variational, "strain12", strain12_variational)
          call MPAS_pool_get_array(velocity_variational, "tanLatVertexRotatedOverRadius", tanLatVertexRotatedOverRadius)
          call MPAS_pool_get_array(velocity_variational, "basisGradientU", basisGradientU)
          call MPAS_pool_get_array(velocity_variational, "basisGradientV", basisGradientV)
       endif

       open(54,file="strain_rate_test_"//trim(unit_test_subtype)//".txt",position="append")

       write(54,*) "-----------------------------------------------------------"
       write(54,*)
       write(54,*) "Test type: ", trim(unit_test_subtype)
       write(54,*)

       allocate(strain11(nCells))
       allocate(strain22(nCells))
       allocate(strain12(nCells))
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

          call cice_grid_rotation_forward(&
               xp,               yp,               zp, &
               xVertex(iVertex), yVertex(iVertex), zVertex(iVertex), &
               config_rotate_cartesian_grid)

          call cice_latlon_from_xyz(&
               latVertexRotated0(iVertex), lonVertexRotated0(iVertex), &
               xp, yp, zp, sphere_radius)

       enddo ! iVertex

       do iCell = 1, nCells

          call cice_grid_rotation_forward(&
               xp,           yp,           zp, &
               xCell(iCell), yCell(iCell), zCell(iCell), &
               config_rotate_cartesian_grid)

          call cice_latlon_from_xyz(&
               latCellRotated0(iCell), lonCellRotated0(iCell), &
               xp, yp, zp, sphere_radius)

       enddo ! iCell

       ! set solveStrain away from rotated pole
       do iCell = 1, nCells

          solveStress(iCell) = interiorCell(iCell)

          if (abs(latCellRotated0(iCell)) > 0.8_RKIND * pii * 0.5_RKIND) then
             solveStress(iCell) = 0
          endif

       enddo ! iCell

       ! set test velocities and expected outputs
       call spherical_test_strain(&
            mesh, &
            uVelocity,              &
            vVelocity,              &
            lonVertexRotated0,      &
            latVertexRotated0,      &
            strain11_test,          &
            strain22_test,          &
            strain12_test,          &
            lonCellRotated0,        &
            latCellRotated0,        &
            nVertices,              &
            nCells,                 &
            trim(unit_test_subtype))

       write(54,*)
       write(54,*) "Analytic solution:"
       write(54,*) "strain11 (min/max):", minval(strain11_test,mask=(solveStress(1:nCells)==1)), &
            maxval(strain11_test,mask=(solveStress(1:nCells)==1))
       write(54,*) "strain22 (min/max):", minval(strain22_test,mask=(solveStress(1:nCells)==1)), &
            maxval(strain22_test,mask=(solveStress(1:nCells)==1))
       write(54,*) "strain12 (min/max):", minval(strain12_test,mask=(solveStress(1:nCells)==1)), &
            maxval(strain12_test,mask=(solveStress(1:nCells)==1))

       ! calculate strain
       if (trim(config_stress_divergence_scheme) == "weak") then

          call cice_strain_tensor_weak(&
               mesh, &
               strain11_weak,       &
               strain22_weak,       &
               strain12_weak,       &
               uVelocity,           &
               vVelocity,           &
               normalVectorPolygon, &
               latCellRotated,      &
               solveStress)

          strain11(1:nCells) = strain11_weak(1:nCells)
          strain22(1:nCells) = strain22_weak(1:nCells)
          strain12(1:nCells) = strain12_weak(1:nCells)

       else if (trim(config_stress_divergence_scheme) == "variational") then

          call cice_strain_tensor_variational(&
               mesh, &
               strain11_variational, &
               strain22_variational, &
               strain12_variational, &
               uVelocity, &
               vVelocity, &
               basisGradientU, &
               basisGradientV, &
               tanLatVertexRotatedOverRadius, &
               solveStress)

          call average_variational_stress(&
               mesh, &
               strain11_variational, &
               strain22_variational, &
               strain12_variational, &
               strain11, &
               strain22, &
               strain12)

       endif

       write(54,*)
       write(54,*) "Numerical solution:"
       write(54,*) "strain11 (min/max):", minval(strain11(1:nCells),mask=(solveStress(1:nCells)==1)), &
            maxval(strain11(1:nCells),mask=(solveStress(1:nCells)==1))
       write(54,*) "strain22 (min/max):", minval(strain22(1:nCells),mask=(solveStress(1:nCells)==1)), &
            maxval(strain22(1:nCells),mask=(solveStress(1:nCells)==1))
       write(54,*) "strain12 (min/max):", minval(strain12(1:nCells),mask=(solveStress(1:nCells)==1)), &
            maxval(strain12(1:nCells),mask=(solveStress(1:nCells)==1))

       ! calculate difference
       strain11_diff = strain11(1:nCells) - strain11_test
       strain22_diff = strain22(1:nCells) - strain22_test
       strain12_diff = strain12(1:nCells) - strain12_test

       write(54,*)
       write(54,*) "Difference:"
       write(54,*) "strain11 (min/max):", minval(strain11_diff,mask=(solveStress(1:nCells)==1)), &
            maxval(strain11_diff,mask=(solveStress(1:nCells)==1))
       write(54,*) "strain22 (min/max):", minval(strain22_diff,mask=(solveStress(1:nCells)==1)), &
            maxval(strain22_diff,mask=(solveStress(1:nCells)==1))
       write(54,*) "strain12 (min/max):", minval(strain12_diff,mask=(solveStress(1:nCells)==1)), &
            maxval(strain12_diff,mask=(solveStress(1:nCells)==1))

       open(55,file="strains_"//trim(unit_test_subtype)//".txt")
       do iCell = 1, nCells

          write(55,*) iCell, solveStress(iCell), &
               strain11(iCell), strain11_test(iCell), &
               strain22(iCell), strain22_test(iCell), &
               strain12(iCell), strain12_test(iCell), &
               strain12(iCell) - strain12_test(iCell)

       enddo ! iCell
       close(55)

       open(55,file="strain11_"//trim(unit_test_subtype)//".txt")
       open(56,file="strain22_"//trim(unit_test_subtype)//".txt")
       open(57,file="strain12_"//trim(unit_test_subtype)//".txt")
       do iCell = 1, nCells
          write(55,*) iCell, lonCell(iCell) * ciceRadiansToDegrees, latCell(iCell) * ciceRadiansToDegrees, &
               strain11(iCell), strain11_test(iCell), strain11_diff(iCell), &
               strain11_diff(iCell) * real(solveStress(iCell),RKIND)
          write(56,*) iCell, lonCell(iCell) * ciceRadiansToDegrees, latCell(iCell) * ciceRadiansToDegrees, &
               strain22(iCell), strain22_test(iCell), strain22_diff(iCell), &
               strain22_diff(iCell) * real(solveStress(iCell),RKIND)
          write(57,*) iCell, lonCell(iCell) * ciceRadiansToDegrees, latCell(iCell) * ciceRadiansToDegrees, &
               strain12(iCell), strain12_test(iCell), strain12_diff(iCell), &
               strain12_diff(iCell) * real(solveStress(iCell),RKIND)
       enddo ! iCell
       close(55)
       close(56)
       close(57)

       ! rms difference
       call cice_rms_difference( &
            strain11(1:nCells),    &
            strain11_test,         &
            solveStress(1:nCells), &
            nCells,                &
            rms_strain11)

       call cice_rms_difference( &
            strain22(1:nCells),    &
            strain22_test,         &
            solveStress(1:nCells), &
            nCells,                &
            rms_strain22)

       call cice_rms_difference( &
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

       do iCell = 1, nCells

          if (solveStress(iCell) == 0) then
             strain11_test(iCell) = 0.0_RKIND
             strain22_test(iCell) = 0.0_RKIND
             strain12_test(iCell) = 0.0_RKIND
             strain11_diff(iCell) = 0.0_RKIND
             strain22_diff(iCell) = 0.0_RKIND
             strain12_diff(iCell) = 0.0_RKIND
          endif

       enddo ! iCell

       ! plot spatial distribution
       call cice_plot_spherical_latlon(mesh, real(interiorCell,RKIND), "interiorCell_"//trim(unit_test_subtype)//".txt")
       call cice_plot_spherical_latlon(mesh, real(solveStress,RKIND), "solveStress_"//trim(unit_test_subtype)//".txt")
       call cice_plot_spherical_latlon(mesh, latCellRotated0, "latCellRotated0_"//trim(unit_test_subtype)//".txt")

       call cice_plot_spherical_latlon(mesh, strain11,      "strain11_plot_"//trim(unit_test_subtype)//"_model.txt")
       call cice_plot_spherical_latlon(mesh, strain11_test, "strain11_plot_"//trim(unit_test_subtype)//"_analytical.txt")
       call cice_plot_spherical_latlon(mesh, strain11_diff, "strain11_plot_"//trim(unit_test_subtype)//"_diff.txt")

       call cice_plot_spherical_latlon(mesh, strain22,      "strain22_plot_"//trim(unit_test_subtype)//"_model.txt")
       call cice_plot_spherical_latlon(mesh, strain22_test, "strain22_plot_"//trim(unit_test_subtype)//"_analytical.txt")
       call cice_plot_spherical_latlon(mesh, strain22_diff, "strain22_plot_"//trim(unit_test_subtype)//"_diff.txt")

       call cice_plot_spherical_latlon(mesh, strain12,      "strain12_plot_"//trim(unit_test_subtype)//"_model.txt")
       call cice_plot_spherical_latlon(mesh, strain12_test, "strain12_plot_"//trim(unit_test_subtype)//"_analytical.txt")
       call cice_plot_spherical_latlon(mesh, strain12_diff, "strain12_plot_"//trim(unit_test_subtype)//"_diff.txt")

       deallocate(strain11)
       deallocate(strain22)
       deallocate(strain12)
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

       write(54,*)

       close(54)

       block => block % next
    enddo

  end subroutine cice_strain_rate_operator_unit_test_individual!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  average_variational_stress
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date December 12th 2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine average_variational_stress(&
       mesh, &
       strain11_variational, &
       strain22_variational, &
       strain12_variational, &
       strain11, &
       strain22, &
       strain12)

    type(mpas_pool_type), intent(in) :: &
         mesh

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         strain11_variational, &
         strain22_variational, &
         strain12_variational

    real(kind=RKIND), dimension(:), intent(out) :: &
         strain11, &
         strain22, &
         strain12

    integer, pointer :: &
         nCells

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer :: &
         iCell, &
         iVertexOnCell

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)
    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)

    do iCell = 1, nCells

       strain11(iCell) = 0.0_RKIND
       strain22(iCell) = 0.0_RKIND
       strain12(iCell) = 0.0_RKIND

       do iVertexOnCell = 1, nEdgesOnCell(iCell)

          strain11(iCell) = strain11(iCell) + strain11_variational(iVertexOnCell,iCell)
          strain22(iCell) = strain22(iCell) + strain22_variational(iVertexOnCell,iCell)
          strain12(iCell) = strain12(iCell) + strain12_variational(iVertexOnCell,iCell)

       enddo ! iVertexOnCell

       strain11(iCell) = strain11(iCell) / real(nEdgesOnCell(iCell), RKIND)
       strain22(iCell) = strain22(iCell) / real(nEdgesOnCell(iCell), RKIND)
       strain12(iCell) = strain12(iCell) / real(nEdgesOnCell(iCell), RKIND)

    enddo ! iCell

  end subroutine average_variational_stress

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  spherical_test_strain
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

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
       nPointsVelocity,   &
       nPointsStrain,     &
       test_type)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    real(kind=RKIND), dimension(:), intent(out) :: &
         uVelocity, & !< Output:
         vVelocity, & !< Output:
         strain11,  & !< Output:
         strain22,  & !< Output:
         strain12     !< Output:

    real(kind=RKIND), dimension(:), intent(in) :: &
         longitudeVelocity, & !< Input:
         latitudeVelocity,  & !< Input:
         longitudeStrain,   & !< Input:
         latitudeStrain       !< Input:

    integer, intent(in) :: &
         nPointsVelocity, & !< Input:
         nPointsStrain      !< Input:

    character(len=*), intent(in) :: &
         test_type !< Input:

    real(kind=RKIND) :: &
         uVelocity_test, &
         vVelocity_test, &
         du_dlon, &
         du_dlat, &
         dv_dlon, &
         dv_dlat

    integer :: &
         iPoint

    real(kind=RKIND), pointer :: &
         sphere_radius

    call MPAS_pool_get_config(mesh, "sphere_radius", sphere_radius)

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
                                             sphere_radius, longitudeStrain(iPoint), latitudeStrain(iPoint))
       strain22(iPoint) = strain22_component(uVelocity_test, vVelocity_test, &
                                             du_dlon, du_dlat, dv_dlon, dv_dlat, &
                                             sphere_radius, longitudeStrain(iPoint), latitudeStrain(iPoint))
       strain12(iPoint) = strain12_component(uVelocity_test, vVelocity_test, &
                                             du_dlon, du_dlat, dv_dlon, dv_dlat, &
                                             sphere_radius, longitudeStrain(iPoint), latitudeStrain(iPoint))

    enddo ! iPoint

  end subroutine spherical_test_strain!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  spherical_test_strain_velocities
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine spherical_test_strain_velocities( &
       u,       &
       v,       &
       du_dlon, &
       du_dlat, &
       dv_dlon, &
       dv_dlat, &
       lon,     &
       lat,     &
       test_type)!{{{

    real(kind=RKIND), intent(out) :: &
         u,       & !< Output:
         v,       & !< Output:
         du_dlon, & !< Output:
         du_dlat, & !< Output:
         dv_dlon, & !< Output:
         dv_dlat    !< Output:

    real(kind=RKIND), intent(in) :: &
         lon, & !< Input:
         lat    !< Input:

    character(len=*), intent(in) :: &
         test_type !< Input:

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

    else

       write(stderrUnit,*) "Unknown test case: ", trim(test_type)
       stop

    endif

  end subroutine spherical_test_strain_velocities!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  strain11_component
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  function strain11_component(&
       u, &
       v, &
       du_dlon, &
       du_dlat, &
       dv_dlon, &
       dv_dlat, &
       r, &
       lon, &
       lat) &
       result(e11)!{{{

    real(kind=RKIND), intent(in) :: &
         u,       & !< Input: 
         v,       & !< Input: 
         du_dlon, & !< Input: 
         du_dlat, & !< Input: 
         dv_dlon, & !< Input: 
         dv_dlat, & !< Input:                                
         r,       & !< Input: 
         lon,     & !< Input: 
         lat        !< Input: 

    real(kind=RKIND) :: e11

    e11 = (1.0_RKIND / (r * cos(lat))) * (du_dlon - v * sin(lat)) 

  end function strain11_component!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  strain22_component
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  function strain22_component(&
       u, &
       v, &
       du_dlon, &
       du_dlat, &
       dv_dlon, &
       dv_dlat, &
       r, &
       lon, &
       lat) &
       result(e22)!{{{

    real(kind=RKIND), intent(in) :: &
         u,       & !< Input:
         v,       & !< Input:
         du_dlon, & !< Input:
         du_dlat, & !< Input:
         dv_dlon, & !< Input:
         dv_dlat, & !< Input:
         r,       & !< Input:
         lon,     & !< Input:
         lat        !< Input:

    real(kind=RKIND) :: e22

    e22 = dv_dlat / r

  end function strain22_component!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  strain12_component
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  function strain12_component(&
       u, &
       v, &
       du_dlon, &
       du_dlat, &
       dv_dlon, &
       dv_dlat, &
       r, &
       lon, &
       lat) &
       result(e12)!{{{

    real(kind=RKIND), intent(in) :: &
         u,       & !< Input:
         v,       & !< Input:
         du_dlon, & !< Input:
         du_dlat, & !< Input:
         dv_dlon, & !< Input:
         dv_dlat, & !< Input:
         r,       & !< Input:
         lon,     & !< Input:
         lat        !< Input:

    real(kind=RKIND) :: e12

    e12 = (0.5_RKIND / r) * (du_dlat + u * tan(lat) + dv_dlon / cos(lat))

  end function strain12_component!}}}

!-----------------------------------------------------------------------
! Spherical stress divergence unit test
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_stress_divergence_operator_unit_test
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_stress_divergence_operator_unit_test(&
       domain, &
       unit_test_subtype)!{{{

    type(domain_type) :: &
         domain !< Input/Output:

    character(len=*), intent(in) :: &
         unit_test_subtype !< Input:

    integer, parameter :: ntests = 8
    character(len=200), dimension(ntests), parameter :: &
         test_names = (/"zero   ", &
                        "const11", &
                        "const22", &
                        "const12", &
                        "test1  ", &
                        "test2  ", &
                        "test3  ", &
                        "test4  "/)

    integer :: &
         itest

    if (trim(unit_test_subtype) == "all") then

       do itest = 1, ntests

          call cice_stress_divergence_operator_unit_test_individual(&
               domain, &
               trim(test_names(itest)))

       enddo ! itest

    else

       call cice_stress_divergence_operator_unit_test_individual(&
            domain, &
            unit_test_subtype)

    endif

  end subroutine cice_stress_divergence_operator_unit_test

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_stress_divergence_operator_unit_test_individual
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_stress_divergence_operator_unit_test_individual(&
       domain, &
       unit_test_subtype)!{{{

    use cice_velocity_solver_weak, only: &
         cice_stress_divergence_weak

    use cice_velocity_solver_variational, only: &
         cice_stress_divergence_variational

    use cice_debug, only: &
         cice_rms_difference, &
         cice_plot_spherical_latlon

    use cice_mesh, only: &
         cice_grid_rotation_forward, &
         cice_latlon_from_xyz

    use cice_constants, only: &
         ciceRadiansToDegrees

    use cice_constants, only: &
         pii

    type(domain_type), intent(inout) :: &
         domain !< Input/Output:

    character(len=*), intent(in) :: &
         unit_test_subtype !< Input:

    type(block_type), pointer :: &
         block

    type (MPAS_pool_type), pointer :: &
         mesh, &
         boundary, &
         velocity_solver, &
         velocity_weak, &
         velocity_variational

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
         solveVelocity, &
         interiorVertex

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
         stress11_weak, &
         stress22_weak, &
         stress12_weak, &
         stressDivergenceU, &
         stressDivergenceV

    integer, dimension(:,:), pointer :: &
         cellVerticesAtVertex

    real(kind=RKIND), dimension(:), pointer :: &
         tanLatVertexRotatedOverRadius

    real(kind=RKIND), dimension(:,:), pointer :: &
         stress11_variational, &
         stress22_variational, &
         stress12_variational

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         basisGradientU, &
         basisGradientV, &
         basisIntegralsU, &
         basisIntegralsV, &
         basisIntegralsMetric

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         normalVectorTriangle

    logical, pointer :: &
         config_rotate_cartesian_grid

    character(len=strKIND), pointer :: &
         config_stress_divergence_scheme

    call MPAS_pool_get_config(domain % configs, "config_rotate_cartesian_grid", config_rotate_cartesian_grid)
    call MPAS_pool_get_config(domain % configs, "config_stress_divergence_scheme", config_stress_divergence_scheme)

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
       call MPAS_pool_get_subpool(block % structs, "boundary", boundary)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocity_solver)
       call MPAS_pool_get_subpool(block % structs, "velocity_weak", velocity_weak)
       call MPAS_pool_get_subpool(block % structs, "velocity_variational", velocity_variational)

       call MPAS_pool_get_dimension(mesh, "nCells", nCells)
       call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)
       call MPAS_pool_get_config(mesh, "sphere_radius", sphere_radius)

       call MPAS_pool_get_array(mesh, "xVertex", xVertex)
       call MPAS_pool_get_array(mesh, "yVertex", yVertex)
       call MPAS_pool_get_array(mesh, "zVertex", zVertex)
       call MPAS_pool_get_array(mesh, "xCell", xCell)
       call MPAS_pool_get_array(mesh, "yCell", yCell)
       call MPAS_pool_get_array(mesh, "zCell", zCell)
       call MPAS_pool_get_array(mesh, "latVertex", latVertex)
       call MPAS_pool_get_array(mesh, "lonVertex", lonVertex)

       call MPAS_pool_get_array(boundary, "interiorVertex", interiorVertex)

       call MPAS_pool_get_array(velocity_solver, "solveVelocity", solveVelocity)
       call MPAS_pool_get_array(velocity_solver, "stressDivergenceU", stressDivergenceU)
       call MPAS_pool_get_array(velocity_solver, "stressDivergenceV", stressDivergenceV)

       if (trim(config_stress_divergence_scheme) == "weak") then
          call MPAS_pool_get_array(velocity_weak, "normalVectorTriangle", normalVectorTriangle)
          call MPAS_pool_get_array(velocity_weak, "latVertexRotated", latVertexRotated)
          call MPAS_pool_get_array(velocity_weak, "stress11", stress11_weak)
          call MPAS_pool_get_array(velocity_weak, "stress22", stress22_weak)
          call MPAS_pool_get_array(velocity_weak, "stress12", stress12_weak)
       else if (trim(config_stress_divergence_scheme) == "variational") then
          call MPAS_pool_get_array(velocity_variational, "stress11", stress11_variational)
          call MPAS_pool_get_array(velocity_variational, "stress22", stress22_variational)
          call MPAS_pool_get_array(velocity_variational, "stress12", stress12_variational)
          call MPAS_pool_get_array(velocity_variational, "basisIntegralsU", basisIntegralsU)
          call MPAS_pool_get_array(velocity_variational, "basisIntegralsV", basisIntegralsV)
          call MPAS_pool_get_array(velocity_variational, "basisIntegralsMetric", basisIntegralsMetric)
          call MPAS_pool_get_array(velocity_variational, "cellVerticesAtVertex", cellVerticesAtVertex)
          call MPAS_pool_get_array(velocity_variational, "tanLatVertexRotatedOverRadius", tanLatVertexRotatedOverRadius)
       endif

       open(54,file="stress_divergence_test_"//trim(unit_test_subtype)//".txt")

       write(54,*) "-----------------------------------------------------------"
       write(54,*)
       write(54,*) "Test type: ", trim(unit_test_subtype)
       write(54,*)

       allocate(stressDivergenceU_test(nVertices))
       allocate(stressDivergenceV_test(nVertices))
       allocate(stressDivergenceU_diff(nVertices))
       allocate(stressDivergenceV_diff(nVertices))

       allocate(lonVertexRotated0(nVertices))
       allocate(latVertexRotated0(nVertices))
       allocate(lonCellRotated0(nCells))
       allocate(latCellRotated0(nCells))

       do iVertex = 1, nVertices

          call cice_grid_rotation_forward(&
               xp,               yp,               zp, &
               xVertex(iVertex), yVertex(iVertex), zVertex(iVertex), &
               config_rotate_cartesian_grid)

          call cice_latlon_from_xyz(&
               latVertexRotated0(iVertex), lonVertexRotated0(iVertex), &
               xp, yp, zp, sphere_radius)

       enddo ! iVertex

       do iCell = 1, nCells

          call cice_grid_rotation_forward(&
               xp,           yp,           zp, &
               xCell(iCell), yCell(iCell), zCell(iCell), &
               config_rotate_cartesian_grid)

          call cice_latlon_from_xyz(&
               latCellRotated0(iCell), lonCellRotated0(iCell), &
               xp, yp, zp, sphere_radius)

       enddo ! iCell

       ! set solveStrain away from rotated pole
       do iVertex = 1, nVertices

          solveVelocity(iVertex) = interiorVertex(iVertex)

          if (abs(latVertexRotated0(iVertex)) > 0.8_RKIND * pii * 0.5_RKIND) then
             solveVelocity(iVertex) = 0
          endif

       enddo ! iCell

       ! set test velocities and expected outputs
       if (trim(config_stress_divergence_scheme) == "weak") then
          call spherical_test_divergence_stress_weak(mesh, &
               stress11_weak,          &
               stress22_weak,          &
               stress12_weak,          &
               lonCellRotated0,        &
               latCellRotated0,        &
               stressDivergenceU_test, &
               stressDivergenceV_test, &
               lonVertexRotated0,      &
               latVertexRotated0,      &
               trim(unit_test_subtype))
       else if (trim(config_stress_divergence_scheme) == "variational") then
          call spherical_test_divergence_stress_variational(mesh, &
               stress11_variational,          &
               stress22_variational,          &
               stress12_variational,          &
               lonVertexRotated0,        &
               latVertexRotated0,        &
               stressDivergenceU_test, &
               stressDivergenceV_test, &
               lonVertexRotated0,      &
               latVertexRotated0,      &
               trim(unit_test_subtype))
       endif

       write(54,*)
       write(54,*) "Analytic solution:"
       write(54,*) "stressDivergenceU:", minval(stressDivergenceU_test,mask=(solveVelocity(1:nVertices)==1)), &
            maxval(stressDivergenceU_test,mask=(solveVelocity(1:nVertices)==1))
       write(54,*) "stressDivergenceV:", minval(stressDivergenceV_test,mask=(solveVelocity(1:nVertices)==1)), &
            maxval(stressDivergenceV_test,mask=(solveVelocity(1:nVertices)==1))

       ! calculate divergence of stress
       if (trim(config_stress_divergence_scheme) == "weak") then

          call cice_stress_divergence_weak(&
               mesh, &
               stressDivergenceU,    &
               stressDivergenceV,    &
               stress11_weak,        &
               stress22_weak,        &
               stress12_weak,        &
               normalVectorTriangle, &
               latVertexRotated,     &
               solveVelocity)

       else if (trim(config_stress_divergence_scheme) == "variational") then

          call cice_stress_divergence_variational(&
               mesh, &
               stressDivergenceU, &
               stressDivergenceV, &
               stress11_variational, &
               stress22_variational, &
               stress12_variational, &
               basisIntegralsU, &
               basisIntegralsV, &
               basisIntegralsMetric, &
               tanLatVertexRotatedOverRadius, &
               cellVerticesAtVertex, &
               solveVelocity)

       endif

       write(54,*)
       write(54,*) "Numerical solution:"
       write(54,*) "stressDivergenceU:", minval(stressDivergenceU(1:nVertices),mask=(solveVelocity(1:nVertices)==1)), &
            maxval(stressDivergenceU(1:nVertices),mask=(solveVelocity(1:nVertices)==1))
       write(54,*) "stressDivergenceV:", minval(stressDivergenceV(1:nVertices),mask=(solveVelocity(1:nVertices)==1)), &
            maxval(stressDivergenceV(1:nVertices),mask=(solveVelocity(1:nVertices)==1))

       ! calculate difference
       stressDivergenceU_diff = stressDivergenceU(1:nVertices) - &
            stressDivergenceU_test
       stressDivergenceV_diff = stressDivergenceV(1:nVertices) - &
            stressDivergenceV_test

       write(54,*)
       write(54,*) "Difference:"
       write(54,*) "stressDivergenceU:", minval(stressDivergenceU_diff,mask=(solveVelocity(1:nVertices)==1)), &
            maxval(stressDivergenceU_diff,mask=(solveVelocity(1:nVertices)==1))
       write(54,*) "stressDivergenceV:", minval(stressDivergenceV_diff,mask=(solveVelocity(1:nVertices)==1)), &
            maxval(stressDivergenceV_diff,mask=(solveVelocity(1:nVertices)==1))

       open(55, file="divergences_"//trim(unit_test_subtype)//".txt")
       do iVertex = 1, nVertices

          write(55,*) iVertex, solveVelocity(iVertex), &
               stressDivergenceU(iVertex), stressDivergenceU_test(iVertex), &
               stressDivergenceV(iVertex), stressDivergenceV_test(iVertex), &
               stressDivergenceU(iVertex) - stressDivergenceU_test(iVertex), &
               stressDivergenceV(iVertex) - stressDivergenceV_test(iVertex)

       enddo ! iVertex
       close(55)

       open(55,file="divergenceu_"//trim(unit_test_subtype)//".txt")
       open(56,file="divergencev_"//trim(unit_test_subtype)//".txt")
       do iVertex = 1, nVertices
          write(55,*) iVertex, lonVertex(iVertex) * ciceRadiansToDegrees, latVertex(iVertex) * ciceRadiansToDegrees, &
               stressDivergenceU(iVertex), stressDivergenceU_test(iVertex), stressDivergenceU_diff(iVertex), &
               stressDivergenceU_diff(iVertex) * real(solveVelocity(iVertex),RKIND)
          write(56,*) iVertex, lonVertex(iVertex) * ciceRadiansToDegrees, latVertex(iVertex) * ciceRadiansToDegrees, &
               stressDivergenceV(iVertex), stressDivergenceV_test(iVertex), stressDivergenceV_diff(iVertex), &
               stressDivergenceV_diff(iVertex) * real(solveVelocity(iVertex),RKIND)
       enddo ! iVertex
       close(55)
       close(56)

       ! rms difference
       call cice_rms_difference( &
            stressDivergenceU(1:nVertices), &
            stressDivergenceU_test,         &
            solveVelocity(1:nVertices),     &
            nVertices,                      &
            rms_divergenceU)

       call cice_rms_difference( &
            stressDivergenceV(1:nVertices), &
            stressDivergenceV_test,         &
            solveVelocity(1:nVertices),     &
            nVertices,                      &
            rms_divergenceV)

       write(54,*)
       write(54,*) "RMS:"
       write(54,*) "rms_divergenceU:", rms_divergenceU
       write(54,*) "rms_divergenceV:", rms_divergenceV

       do iVertex = 1, nVertices

          if (solveVelocity(iVertex) == 0) then
             stressDivergenceU_test(iVertex) = 0.0_RKIND
             stressDivergenceU_test(iVertex) = 0.0_RKIND
             stressDivergenceU_diff(iVertex) = 0.0_RKIND
             stressDivergenceU_diff(iVertex) = 0.0_RKIND
          endif

       enddo ! iCell

       ! plot spatial distribution
       call cice_plot_spherical_latlon(mesh, real(interiorVertex,RKIND), "interiorVertex_"//trim(unit_test_subtype)//".txt")
       call cice_plot_spherical_latlon(mesh, real(solveVelocity,RKIND), "solveVelocity_"//trim(unit_test_subtype)//".txt")
       call cice_plot_spherical_latlon(mesh, latVertexRotated0, "latVertexRotated0_"//trim(unit_test_subtype)//".txt")

       call cice_plot_spherical_latlon(mesh, stressDivergenceU,      "stressDivergenceU_plot_"//trim(unit_test_subtype)//"_model.txt")
       call cice_plot_spherical_latlon(mesh, stressDivergenceU_test, "stressDivergenceU_plot_"//trim(unit_test_subtype)//"_analytical.txt")
       call cice_plot_spherical_latlon(mesh, stressDivergenceU_diff, "stressDivergenceU_plot_"//trim(unit_test_subtype)//"_diff.txt")

       call cice_plot_spherical_latlon(mesh, stressDivergenceV,      "stressDivergenceV_plot_"//trim(unit_test_subtype)//"_model.txt")
       call cice_plot_spherical_latlon(mesh, stressDivergenceV_test, "stressDivergenceV_plot_"//trim(unit_test_subtype)//"_analytical.txt")
       call cice_plot_spherical_latlon(mesh, stressDivergenceV_diff, "stressDivergenceV_plot_"//trim(unit_test_subtype)//"_diff.txt")

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

       block => block % next
    enddo

  end subroutine cice_stress_divergence_operator_unit_test_individual!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  spherical_test_divergence_stress_weak
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine spherical_test_divergence_stress_weak(mesh, &
       stress11,            &
       stress22,            &
       stress12,            &
       longitudeStress,     &
       latitudeStress,      &
       stressDivergenceU,   &
       stressDivergenceV,   &
       longitudeDivergence, &
       latitudeDivergence,  &
       test_type)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    real(kind=RKIND), dimension(:), intent(out) :: &
         stressDivergenceU, & !< Output:
         stressDivergenceV, & !< Output:
         stress11,          & !< Output:
         stress22,          & !< Output:
         stress12             !< Output:

    real(kind=RKIND), dimension(:), intent(in) :: &
         longitudeStress,     & !< Input:
         latitudeStress,      & !< Input:
         longitudeDivergence, & !< Input:
         latitudeDivergence     !< Input:

    character(len=*), intent(in) :: &
         test_type !< Input:

    real(kind=RKIND) :: &
         stress11_test,  &
         stress22_test,  &
         stress12_test,  &
         dstress11_dlon, &
         dstress11_dlat, &
         dstress22_dlon, &
         dstress22_dlat, &
         dstress12_dlon, &
         dstress12_dlat

    integer :: &
         nPointsStress, &
         nPointsDivergence, &
         iPoint

    real(kind=RKIND), pointer :: &
         sphere_radius

    call MPAS_pool_get_config(mesh, "sphere_radius", sphere_radius)

    nPointsStress     = size(stress11)-1
    nPointsDivergence = size(stressDivergenceU)-1

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
                                                     sphere_radius, longitudeDivergence(iPoint), latitudeDivergence(iPoint))
       stressDivergenceV(iPoint) = divergenceStressV(stress11_test, stress22_test, stress12_test, &
                                                     dstress11_dlon, dstress11_dlat, dstress22_dlon, dstress22_dlat, dstress12_dlon, dstress12_dlat, &
                                                     sphere_radius, longitudeDivergence(iPoint), latitudeDivergence(iPoint))

    enddo ! iPoint

  end subroutine spherical_test_divergence_stress_weak!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  spherical_test_divergence_stress_variational
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine spherical_test_divergence_stress_variational(mesh, &
       stress11,            &
       stress22,            &
       stress12,            &
       longitudeStress,     &
       latitudeStress,      &
       stressDivergenceU,   &
       stressDivergenceV,   &
       longitudeDivergence, &
       latitudeDivergence,  &
       test_type)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    real(kind=RKIND), dimension(:), intent(out) :: &
         stressDivergenceU, & !< Output:
         stressDivergenceV    !< Output:

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         stress11,          & !< Output:
         stress22,          & !< Output:
         stress12             !< Output:

    real(kind=RKIND), dimension(:), intent(in) :: &
         longitudeStress,     & !< Input:
         latitudeStress,      & !< Input:
         longitudeDivergence, & !< Input:
         latitudeDivergence     !< Input:

    character(len=*), intent(in) :: &
         test_type !< Input:

    real(kind=RKIND) :: &
         stress11_test,  &
         stress22_test,  &
         stress12_test,  &
         dstress11_dlon, &
         dstress11_dlat, &
         dstress22_dlon, &
         dstress22_dlat, &
         dstress12_dlon, &
         dstress12_dlat

    integer, pointer :: &
         nCells, &
         nVertices

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         verticesOnCell

    integer :: &
         iVertex, &
         iCell, &
         iVertexOnCell

    real(kind=RKIND), pointer :: &
         sphere_radius

    call MPAS_pool_get_config(mesh, "sphere_radius", sphere_radius)

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)
    call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)

    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "verticesOnCell", verticesOnCell)

    ! set stress points
    do iCell = 1, nCells

       do iVertexOnCell = 1, nEdgesOnCell(iCell)

          iVertex = verticesOnCell(iVertexOnCell,iCell)

          call spherical_test_divergence_stress_stresses( &
               stress11(iVertexOnCell,iCell),        &
               stress22(iVertexOnCell,iCell),        &
               stress12(iVertexOnCell,iCell),        &
               dstress11_dlon,          &
               dstress11_dlat,          &
               dstress22_dlon,          &
               dstress22_dlat,          &
               dstress12_dlon,          &
               dstress12_dlat,          &
               longitudeStress(iVertex), &
               latitudeStress(iVertex),  &
               test_type)

       enddo ! iVertexOnCell

    enddo ! iCell

    ! set divergence points
    do iVertex = 1, nVertices

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
         longitudeDivergence(iVertex), &
         latitudeDivergence(iVertex), &
         test_type)

       stressDivergenceU(iVertex) = divergenceStressU(stress11_test, stress22_test, stress12_test, &
                                                      dstress11_dlon, dstress11_dlat, dstress22_dlon, dstress22_dlat, dstress12_dlon, dstress12_dlat, &
                                                      sphere_radius, longitudeDivergence(iVertex), latitudeDivergence(iVertex))
       stressDivergenceV(iVertex) = divergenceStressV(stress11_test, stress22_test, stress12_test, &
                                                      dstress11_dlon, dstress11_dlat, dstress22_dlon, dstress22_dlat, dstress12_dlon, dstress12_dlat, &
                                                      sphere_radius, longitudeDivergence(iVertex), latitudeDivergence(iVertex))

    enddo ! iVertex

  end subroutine spherical_test_divergence_stress_variational!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  spherical_test_divergence_stress_stresses
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

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
       test_type)!{{{

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
         lon, & !< Input:
         lat    !< Input:

    character(len=*), intent(in) :: &
         test_type !< Input:

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

       write(stderrUnit,*) "Unknown test case: ", trim(test_type)
       stop

    endif

  end subroutine spherical_test_divergence_stress_stresses!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  divergenceStressU
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  function divergenceStressU(&
       stress11, &
       stress22, &
       stress12, &
       dstress11_dlon, &
       dstress11_dlat, &
       dstress22_dlon, &
       dstress22_dlat, &
       dstress12_dlon, &
       dstress12_dlat, &
       r, &
       lon, &
       lat) &
       result(divu)!{{{

    real(kind=RKIND), intent(in) :: &
         stress11,       & !< Input:
         stress22,       & !< Input:
         stress12,       & !< Input:
         dstress11_dlon, & !< Input:
         dstress11_dlat, & !< Input:
         dstress22_dlon, & !< Input:
         dstress22_dlat, & !< Input:
         dstress12_dlon, & !< Input:
         dstress12_dlat, & !< Input:
         r,              & !< Input:
         lon,            & !< Input:
         lat               !< Input:

    real(kind=RKIND) :: divu

    divu = (1.0_RKIND / (r * cos(lat))) * dstress11_dlon + &
           (1.0_RKIND / r)              * dstress12_dlat - &
           (2.0_RKIND / r) * tan(lat)   * stress12

  end function divergenceStressU!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  divergenceStressV
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  function divergenceStressV(&
       stress11, &
       stress22, &
       stress12, &
       dstress11_dlon, &
       dstress11_dlat, &
       dstress22_dlon, &
       dstress22_dlat, &
       dstress12_dlon, &
       dstress12_dlat, &
       r, &
       lon, &
       lat) &
       result(divv)!{{{

    real(kind=RKIND), intent(in) :: &
         stress11,       & !< Input:
         stress22,       & !< Input:
         stress12,       & !< Input:
         dstress11_dlon, & !< Input:
         dstress11_dlat, & !< Input:
         dstress22_dlon, & !< Input:
         dstress22_dlat, & !< Input:
         dstress12_dlon, & !< Input:
         dstress12_dlat, & !< Input:
         r,              & !< Input:
         lon,            & !< Input:
         lat               !< Input:

    real(kind=RKIND) :: divv

    divv = (1.0_RKIND / (r * cos(lat))) * dstress12_dlon + &
           (1.0_RKIND / r)              * dstress22_dlat + &
           (1.0_RKIND / r) * tan(lat)   * stress11       - &
           (1.0_RKIND / r) * tan(lat)   * stress22

  end function divergenceStressV!}}}

!-----------------------------------------------------------------------
! EVP Constitutive relationship unit test
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_constitutive_relationship_unit_test
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_constitutive_relationship_unit_test(&
       domain)!{{{

    use cice_velocity_solver_constitutive_relation, only: &
         cice_evp_constitutive_relation, &
         cice_evp_constitutive_relation_revised, &
         cice_init_evp

    type(domain_type), intent(inout) :: &
         domain !< Input/Output:

    type(block_type), pointer :: &
         block

    integer, parameter :: &
         nTests = 2

    real(kind=RKIND), dimension(nTests) :: &
         icePressure, &
         stress11, &
         stress22, &
         stress12, &
         strain11, &
         strain22, &
         strain12

    real(kind=RKIND) :: &
         replacementPressure, &
         areaCell

    logical, pointer :: &
         config_revised_evp

    type(MPAS_pool_type), pointer :: &
         velocitySolver

    real(kind=RKIND), pointer :: &
         elasticTimeStep

    integer :: &
         iTest

    call MPAS_pool_get_config(domain % configs, "config_revised_evp", config_revised_evp)

    call cice_init_evp(domain)

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolver)
       call MPAS_pool_get_array(velocitySolver, "elasticTimeStep", elasticTimeStep)

       call constitutive_relationship_set_strains(&
            stress11, &
            stress22, &
            stress12, &
            icePressure, &
            strain11, &
            strain22, &
            strain12)

       areaCell = 1.0_RKIND

       do iTest = 1, nTests

          if (.not. config_revised_evp) then

             call cice_evp_constitutive_relation(&
                  stress11(iTest), &
                  stress22(iTest), &
                  stress12(iTest), &
                  strain11(iTest), &
                  strain22(iTest), &
                  strain12(iTest), &
                  icePressure(iTest), &
                  replacementPressure,   &
                  areaCell, &
                  elasticTimeStep)

          else

             call cice_evp_constitutive_relation_revised(&
                  stress11(iTest), &
                  stress22(iTest), &
                  stress12(iTest), &
                  strain11(iTest), &
                  strain22(iTest), &
                  strain12(iTest), &
                  icePressure(iTest), &
                  replacementPressure,   &
                  areaCell)

          endif

       enddo ! iTest

       call constitutive_relationship_writeout(&
            icePressure, &
            strain11, &
            strain22, &
            strain12, &
            stress11, &
            stress22, &
            stress12)

       block => block % next
    enddo

  end subroutine cice_constitutive_relationship_unit_test!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  constitutive_relationship_set_strains
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine constitutive_relationship_set_strains(&
       stress11, &
       stress22, &
       stress12, &
       icePressure, &
       strain11, &
       strain22, &
       strain12)!{{{

    real(kind=RKIND), dimension(:), intent(out) :: &
         stress11, &    !< Output:
         stress22, &    !< Output:
         stress12, &    !< Output:
         icePressure, & !< Output:
         strain11, &    !< Output:
         strain22, &    !< Output:
         strain12       !< Output:

    integer :: &
         iTest

    ! test 1
    iTest = 1
    stress11(iTest) = 0.0_RKIND
    stress22(iTest) = 0.0_RKIND
    stress12(iTest) = 0.0_RKIND
    icePressure(iTest) = 0.0_RKIND
    strain11(iTest) = 0.0_RKIND
    strain22(iTest) = 0.0_RKIND
    strain12(iTest) = 0.0_RKIND

    ! test 2
    iTest = 2
    stress11(iTest) = 0.0_RKIND
    stress22(iTest) = 0.0_RKIND
    stress12(iTest) = 0.0_RKIND
    icePressure(iTest) = 1000000000.0_RKIND
    strain11(iTest) = 1.0_RKIND
    strain22(iTest) = 1.0_RKIND
    strain12(iTest) = 1.0_RKIND

  end subroutine constitutive_relationship_set_strains!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  constitutive_relationship_writeout
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine constitutive_relationship_writeout(&
       icePressure, &
       strain11, &
       strain22, &
       strain12, &
       stress11, &
       stress22, &
       stress12)!{{{

    real(kind=RKIND), dimension(:), intent(in) :: &
         icePressure, & !< Input:
         strain11, &    !< Input:
         strain22, &    !< Input:
         strain12, &    !< Input:
         stress11, &    !< Input:
         stress22, &    !< Input:
         stress12       !< Input:

    integer :: &
         iTest

    write(*,*) "Constitutive tests"
    write(*,*)

    do iTest = 1, size(icePressure)

       write(*,*) "Constitutive test ", iTest

       write(*,*) "strength: ", icePressure(iTest)
       write(*,*) "strain11: ", strain11(iTest)
       write(*,*) "strain22: ", strain22(iTest)
       write(*,*) "strain12: ", strain12(iTest)
       write(*,*) "stress11: ", stress11(iTest)
       write(*,*) "stress22: ", stress22(iTest)
       write(*,*) "stress12: ", stress12(iTest)

    enddo ! iTest

  end subroutine constitutive_relationship_writeout!}}}

  !--------------------------------------------------------------------------

end module cice_velocity_solver_unit_tests
