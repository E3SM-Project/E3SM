module cice_dynamics_tri

  use mpas_grid_types

  implicit none

  private
  public :: init_dynamics_tri, &
            run_dynamics_tri

contains

  !-------------------------------------------------------------
  ! Initialization
  !-------------------------------------------------------------
  
  subroutine init_dynamics_tri(mesh, &
                               boundary, normal, &
                               dt, divergenceFormulation)

    use cice_dynamics_shared, only: init_evp
    use cice_dyn_tri_weak, only: init_dynamics_tri_weak

    type(mesh_type), intent(inout) :: mesh

    type(boundary_type), pointer :: boundary
    type(normal_type),   pointer :: normal

    real(kind=RKIND), intent(in) :: dt

    character(len=*), intent(in) :: divergenceFormulation

    call init_evp(dt)
    
    call init_dynamics_tri_weak(mesh, boundary, normal)
    
  end subroutine init_dynamics_tri

  !-------------------------------------------------------------
  ! Time step
  !-------------------------------------------------------------

  subroutine run_dynamics_tri(mesh, &
                              icestate, boundary,        &
                              normal,   atmos_forcing,   &
                              tridyn,   triweak,         &
                              trifor,                    &
                              dt, divergenceFormulation, &
                              evpDamping)

    use cice_dynamics_shared, only: ice_strength, &
                                    interpolate_cell_to_vertex, &
                                    air_stress, &
                                    ocean_stress, &
                                    surface_tilt, &
                                    principal_stresses, &
                                    dtDynamics, &
                                    dtElastic, &
                                    elapsedTime

    use cice_testing, only: gnuplot_triangle, gnuplot_cell, writeout_minmax

    type(mesh_type), intent(in) :: mesh

    type(icestate_type),      pointer :: icestate
    type(boundary_type),      pointer :: boundary
    type(normal_type),        pointer :: normal
    type(atmos_forcing_type), pointer :: atmos_forcing    
    type(tridyn_type),        pointer :: tridyn
    type(triweak_type),       pointer :: triweak
    type(trifor_type),        pointer :: trifor

    real(kind=RKIND), intent(in) :: &
         dt

    character(len=*), intent(in) :: &
         divergenceFormulation

    logical, intent(in) :: &
         evpDamping

    call ice_strength(mesh, &     
         tridyn % icePressure % array,    &
         icestate % iceAreaCell % array,  &
         icestate % iceVolumeCell % array, &
         boundary % interiorCell % array)

    call interpolate_cell_to_vertex(mesh, &
         triweak % icePressureVertex % array, &
         tridyn % icePressure % array)

    call air_stress(&
         tridyn % airStressU % array,          &
         tridyn % airStressV % array,          &
         atmos_forcing % uAirVelocity % array, & 
         atmos_forcing % vAirVelocity % array, &
         icestate % iceAreaCell % array)

    call ocean_stress(& 
         tridyn % oceanStressU % array,   &
         tridyn % oceanStressV % array,   &
         trifor % uOceanVelocity % array, & 
         trifor % vOceanVelocity % array, &
         mesh % fCell % array,            &
         boundary % interiorCell % array)
    
    call surface_tilt(&
         tridyn % surfaceTiltForceU % array, &
         tridyn % surfaceTiltForceV % array, &
         trifor % uOceanVelocity % array,    & 
         trifor % vOceanVelocity % array,    &
         icestate % totalMassCell % array,   & 
         mesh % fVertex % array,             &
         boundary % interiorCell % array)

    call subcycle_dynamics_tri(mesh, &
                               icestate,   boundary,  &
                               normal,                &
                               tridyn,     triweak,   &
                               trifor,                &
                               dtDynamics, dtElastic, &
                               divergenceFormulation, &
                               evpDamping)

    elapsedTime = elapsedTime + dt

  end subroutine run_dynamics_tri

  !-------------------------------------------------------------

  subroutine subcycle_dynamics_tri(mesh, &
                                   icestate,   boundary,  &
                                   normal,                &
                                   tridyn,     triweak,   &
                                   trifor,                &
                                   dtDynamics, dtElastic, &
                                   divergenceFormulation, &
                                   evpDamping)

    use cice_dynamics_shared, only: nElasticSubcycle

    use cice_testing, only: gnuplot_triangle, gnuplot_cell, writeout_minmax, gnuplot_vertexvector

    type(mesh_type), intent(in) :: mesh

    type(icestate_type), pointer :: icestate
    type(boundary_type), pointer :: boundary
    type(normal_type),   pointer :: normal
    type(tridyn_type),   pointer :: tridyn
    type(triweak_type),  pointer :: triweak
    type(trifor_type),   pointer :: trifor

    real(kind=RKIND), intent(in) :: &
         dtDynamics, &
         dtElastic

    character(len=*), intent(in) :: &
         divergenceFormulation

    logical, intent(in) :: &
         evpDamping

    integer :: &
         iElasticSubcycle

    do iElasticSubcycle = 1, nElasticSubcycle

       write(*,*) "subcycle: ", iElasticSubcycle
       call single_subcycle_dynamics_tri(mesh, &
                                         icestate,   boundary,  &
                                         normal,                &
                                         tridyn,     triweak,   &
                                         trifor,                &
                                         dtDynamics, dtElastic, &
                                         divergenceFormulation, &
                                         evpDamping)

    enddo

  end subroutine subcycle_dynamics_tri

  !-------------------------------------------------------------

  subroutine single_subcycle_dynamics(mesh, &
                                      icestate,   boundary,  &
                                      normal,                &
                                      tridyn,     triweak,   &
                                      trifor,                &
                                      dtDynamics, dtElastic, &
                                      divergenceFormulation, &
                                      evpDamping)

    use cice_dynamics_shared, only: interpolate_cell_to_vertex, &
                                    ocean_stress_coefficient, &
                                    solve_velocity

    use cice_dyn_tri_weak, only: strain_tensor_tri_weak, &
                                 stress_tensor_tri_weak, &
                                 stress_divergence_tri_weak

    use cice_testing, only: gnuplot_triangle, gnuplot_cell, writeout_minmax, gnuplot_vertexvector, &
                            divergence_stress_test_velocity_set

    type(mesh_type), intent(in) :: mesh

    type(icestate_type), pointer :: icestate
    type(boundary_type), pointer :: boundary
    type(normal_type),   pointer :: normal
    type(tridyn_type),   pointer :: tridyn
    type(triweak_type),  pointer :: triweak
    type(trifor_type),   pointer :: trifor

    real(kind=RKIND), intent(in) :: &
         dtDynamics, &
         dtElastic

    character(len=*), intent(in) :: &
         divergenceFormulation
    
    logical, intent(in) :: &
         evpDamping

    call strain_tensor_tri_weak(mesh, &
         triweak % strain11 % array,            &
         triweak % strain22 % array,            &
         triweak % strain12 % array,            &
         tridyn % uVelocity % array,            &
         tridyn % vVelocity % array,            &
         normal % normalVectorTriangle % array, &
         boundary % interiorVertex2 % array)
    
    call stress_tensor_tri_weak(mesh, &
         triweak % stress11 % array,          &
         triweak % stress22 % array,          &
         triweak % stress12 % array,          &
         triweak % strain11 % array,          &
         triweak % strain22 % array,          &
         triweak % strain12 % array,          &
         triweak % icePressureVertex % array, &
         dtElastic, evpDamping,               &
         boundary % interiorVertex % array)
    
    call stress_divergence_tri_weak(mesh, &
         tridyn % stressDivergenceU % array,   &
         tridyn % stressDivergenceV % array,   &
         triweak % stress11 % array,           &
         triweak % stress22 % array,           &
         triweak % stress12 % array,           &
         normal % normalVectorPolygon % array)
    
    call interpolate_cell_to_vertex(mesh, &
         icestate % totalMassVertex % array, &
         icestate % totalMassCell % array)

    call ocean_stress_coefficient(& 
         tridyn % oceanStressCoeff % array, &
         trifor % uOceanVelocity % array,   & 
         trifor % vOceanVelocity % array,   &
         tridyn % uVelocity % array,        &   
         tridyn % vVelocity % array,        &
         icestate % iceAreaCell % array,    &
         boundary % interiorCell % array)

    call solve_velocity(mesh % nVertices, &
         boundary % interiorCell % array,    &
         tridyn % uVelocity % array,         &   
         tridyn % vVelocity % array,         &
         icestate % totalMassCell % array,   & 
         mesh % fCell % array,               &
         tridyn % stressDivergenceU % array, & 
         tridyn % stressDivergenceV % array, &
         tridyn % airStressU % array,        & 
         tridyn % airStressV % array,        &
         tridyn % surfaceTiltForceU % array, &
         tridyn % surfaceTiltForceV % array, &
         tridyn % oceanStressU % array,      &
         tridyn % oceanStressV % array,      &
         tridyn % oceanStressCoeff % array,  &
         dtElastic)

  end subroutine single_subcycle_dynamics

  !-------------------------------------------------------------

end module cice_dynamics_tri

