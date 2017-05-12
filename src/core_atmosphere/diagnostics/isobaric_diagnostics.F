! Copyright (c) 2016,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module isobaric_diagnostics

    use mpas_dmpar
    use mpas_kind_types
    use mpas_derived_types
    use mpas_pool_routines
    use mpas_constants

    type (MPAS_pool_type), pointer :: mesh
    type (MPAS_pool_type), pointer :: state
    type (MPAS_pool_type), pointer :: diag

    type (MPAS_clock_type), pointer :: clock

    public :: isobaric_diagnostics_setup, &
              isobaric_diagnostics_compute

    private

    logical :: need_mslp, &
               need_relhum_200, need_relhum_250, need_relhum_500, need_relhum_700, need_relhum_850, need_relhum_925, &
               need_dewpoint_200, need_dewpoint_250, need_dewpoint_500, need_dewpoint_700, need_dewpoint_850, need_dewpoint_925, &
               need_temp_200, need_temp_250, need_temp_500, need_temp_700, need_temp_850, need_temp_925, &
               need_height_200, need_height_250, need_height_500, need_height_700, need_height_850, need_height_925, &
               need_uzonal_200, need_uzonal_250, need_uzonal_500, need_uzonal_700, need_uzonal_850, need_uzonal_925, &
               need_umeridional_200, need_umeridional_250, need_umeridional_500, need_umeridional_700, need_umeridional_850, need_umeridional_925, &
               need_w_200, need_w_250, need_w_500, need_w_700, need_w_850, need_w_925, &
               need_vorticity_200, need_vorticity_250, need_vorticity_500, need_vorticity_700, need_vorticity_850, need_vorticity_925, &
               need_t_isobaric, need_z_isobaric, need_meanT_500_300
    logical :: need_temp, need_relhum, need_dewpoint, need_w, need_uzonal, need_umeridional, need_vorticity, need_height


    contains


    !-----------------------------------------------------------------------
    !  routine isobaric_diagnostics_setup
    !
    !> \brief Set up the isobaric diagnostics module
    !> \author Michael Duda
    !> \date   21 October 2016
    !> \details
    !>  This routine sets up the isobaric diagnostics module, principally by
    !>  saving pointers to pools that are used in the computation of diagnostics.
    !
    !-----------------------------------------------------------------------
    subroutine isobaric_diagnostics_setup(all_pools, simulation_clock)

        use mpas_derived_types, only : MPAS_pool_type, MPAS_clock_type
        use mpas_pool_routines, only : mpas_pool_get_subpool

        implicit none

        type (MPAS_pool_type), pointer :: all_pools
        type (MPAS_clock_type), pointer :: simulation_clock

        clock => simulation_clock

        call mpas_pool_get_subpool(all_pools, 'mesh', mesh)
        call mpas_pool_get_subpool(all_pools, 'state', state)
        call mpas_pool_get_subpool(all_pools, 'diag', diag)
   
    end subroutine isobaric_diagnostics_setup


    !-----------------------------------------------------------------------
    !  routine isobaric_diagnostics_compute
    !
    !> \brief Compute isobaric diagnostic before model output is written
    !> \author Michael Duda
    !> \date   21 October 2016
    !> \details
    !>  Compute isobaric diagnostic before model output is written. Code called
    !>  from here was previously in mpas_atm_interp_diagnostics.F.
    !
    !-----------------------------------------------------------------------
    subroutine isobaric_diagnostics_compute()

        use mpas_atm_diagnostics_utils, only : MPAS_field_will_be_written

        implicit none

        logical :: need_any_diags

        need_any_diags = .false.

        need_temp = .false.
        need_dewpoint = .false.
        need_relhum = .false.
        need_w = .false.
        need_uzonal = .false.
        need_umeridional = .false.
        need_vorticity = .false.
        need_height = .false.

        need_mslp = MPAS_field_will_be_written('mslp')
        need_any_diags = need_any_diags .or. need_mslp
        need_relhum_200 = MPAS_field_will_be_written('relhum_200hPa')
        need_relhum = need_relhum .or. need_relhum_200
        need_any_diags = need_any_diags .or. need_relhum_200
        need_relhum_250 = MPAS_field_will_be_written('relhum_250hPa')
        need_relhum = need_relhum .or. need_relhum_250
        need_any_diags = need_any_diags .or. need_relhum_250
        need_relhum_500 = MPAS_field_will_be_written('relhum_500hPa')
        need_relhum = need_relhum .or. need_relhum_500
        need_any_diags = need_any_diags .or. need_relhum_500
        need_relhum_700 = MPAS_field_will_be_written('relhum_700hPa')
        need_relhum = need_relhum .or. need_relhum_700
        need_any_diags = need_any_diags .or. need_relhum_700
        need_relhum_850 = MPAS_field_will_be_written('relhum_850hPa')
        need_relhum = need_relhum .or. need_relhum_850
        need_any_diags = need_any_diags .or. need_relhum_850
        need_relhum_925 = MPAS_field_will_be_written('relhum_925hPa')
        need_relhum = need_relhum .or. need_relhum_925
        need_any_diags = need_any_diags .or. need_relhum_925
        need_dewpoint_200 = MPAS_field_will_be_written('dewpoint_200hPa')
        need_dewpoint = need_dewpoint .or. need_dewpoint_200
        need_any_diags = need_any_diags .or. need_dewpoint_200
        need_dewpoint_250 = MPAS_field_will_be_written('dewpoint_250hPa')
        need_dewpoint = need_dewpoint .or. need_dewpoint_250
        need_any_diags = need_any_diags .or. need_dewpoint_250
        need_dewpoint_500 = MPAS_field_will_be_written('dewpoint_500hPa')
        need_dewpoint = need_dewpoint .or. need_dewpoint_500
        need_any_diags = need_any_diags .or. need_dewpoint_500
        need_dewpoint_700 = MPAS_field_will_be_written('dewpoint_700hPa')
        need_dewpoint = need_dewpoint .or. need_dewpoint_700
        need_any_diags = need_any_diags .or. need_dewpoint_700
        need_dewpoint_850 = MPAS_field_will_be_written('dewpoint_850hPa')
        need_dewpoint = need_dewpoint .or. need_dewpoint_850
        need_any_diags = need_any_diags .or. need_dewpoint_850
        need_dewpoint_925 = MPAS_field_will_be_written('dewpoint_925hPa')
        need_dewpoint = need_dewpoint .or. need_dewpoint_925
        need_any_diags = need_any_diags .or. need_dewpoint_925
        need_temp_200 = MPAS_field_will_be_written('temperature_200hPa')
        need_temp = need_temp .or. need_temp_200
        need_any_diags = need_any_diags .or. need_temp_200
        need_temp_250 = MPAS_field_will_be_written('temperature_250hPa')
        need_temp = need_temp .or. need_temp_250
        need_any_diags = need_any_diags .or. need_temp_250
        need_temp_500 = MPAS_field_will_be_written('temperature_500hPa')
        need_temp = need_temp .or. need_temp_500
        need_any_diags = need_any_diags .or. need_temp_500
        need_temp_700 = MPAS_field_will_be_written('temperature_700hPa')
        need_temp = need_temp .or. need_temp_700
        need_any_diags = need_any_diags .or. need_temp_700
        need_temp_850 = MPAS_field_will_be_written('temperature_850hPa')
        need_temp = need_temp .or. need_temp_850
        need_any_diags = need_any_diags .or. need_temp_850
        need_temp_925 = MPAS_field_will_be_written('temperature_925hPa')
        need_temp = need_temp .or. need_temp_925
        need_any_diags = need_any_diags .or. need_temp_925
        need_height_200 = MPAS_field_will_be_written('height_200hPa')
        need_height = need_height .or. need_height_200
        need_any_diags = need_any_diags .or. need_height_200
        need_height_250 = MPAS_field_will_be_written('height_250hPa')
        need_height = need_height .or. need_height_250
        need_any_diags = need_any_diags .or. need_height_250
        need_height_500 = MPAS_field_will_be_written('height_500hPa')
        need_height = need_height .or. need_height_500
        need_any_diags = need_any_diags .or. need_height_500
        need_height_700 = MPAS_field_will_be_written('height_700hPa')
        need_height = need_height .or. need_height_700
        need_any_diags = need_any_diags .or. need_height_700
        need_height_850 = MPAS_field_will_be_written('height_850hPa')
        need_height = need_height .or. need_height_850
        need_any_diags = need_any_diags .or. need_height_850
        need_height_925 = MPAS_field_will_be_written('height_925hPa')
        need_height = need_height .or. need_height_925
        need_any_diags = need_any_diags .or. need_height_925
        need_uzonal_200 = MPAS_field_will_be_written('uzonal_200hPa')
        need_uzonal = need_uzonal .or. need_uzonal_200
        need_any_diags = need_any_diags .or. need_uzonal_200
        need_uzonal_250 = MPAS_field_will_be_written('uzonal_250hPa')
        need_uzonal = need_uzonal .or. need_uzonal_250
        need_any_diags = need_any_diags .or. need_uzonal_250
        need_uzonal_500 = MPAS_field_will_be_written('uzonal_500hPa')
        need_uzonal = need_uzonal .or. need_uzonal_500
        need_any_diags = need_any_diags .or. need_uzonal_500
        need_uzonal_700 = MPAS_field_will_be_written('uzonal_700hPa')
        need_uzonal = need_uzonal .or. need_uzonal_700
        need_any_diags = need_any_diags .or. need_uzonal_700
        need_uzonal_850 = MPAS_field_will_be_written('uzonal_850hPa')
        need_uzonal = need_uzonal .or. need_uzonal_850
        need_any_diags = need_any_diags .or. need_uzonal_850
        need_uzonal_925 = MPAS_field_will_be_written('uzonal_925hPa')
        need_uzonal = need_uzonal .or. need_uzonal_925
        need_any_diags = need_any_diags .or. need_uzonal_925
        need_umeridional_200 = MPAS_field_will_be_written('umeridional_200hPa')
        need_umeridional = need_umeridional .or. need_umeridional_200
        need_any_diags = need_any_diags .or. need_umeridional_200
        need_umeridional_250 = MPAS_field_will_be_written('umeridional_250hPa')
        need_umeridional = need_umeridional .or. need_umeridional_250
        need_any_diags = need_any_diags .or. need_umeridional_250
        need_umeridional_500 = MPAS_field_will_be_written('umeridional_500hPa')
        need_umeridional = need_umeridional .or. need_umeridional_500
        need_any_diags = need_any_diags .or. need_umeridional_500
        need_umeridional_700 = MPAS_field_will_be_written('umeridional_700hPa')
        need_umeridional = need_umeridional .or. need_umeridional_700
        need_any_diags = need_any_diags .or. need_umeridional_700
        need_umeridional_850 = MPAS_field_will_be_written('umeridional_850hPa')
        need_umeridional = need_umeridional .or. need_umeridional_850
        need_any_diags = need_any_diags .or. need_umeridional_850
        need_umeridional_925 = MPAS_field_will_be_written('umeridional_925hPa')
        need_umeridional = need_umeridional .or. need_umeridional_925
        need_any_diags = need_any_diags .or. need_umeridional_925
        need_w_200 = MPAS_field_will_be_written('w_200hPa')
        need_w = need_w .or. need_w_200
        need_any_diags = need_any_diags .or. need_w_200
        need_w_250 = MPAS_field_will_be_written('w_250hPa')
        need_w = need_w .or. need_w_250
        need_any_diags = need_any_diags .or. need_w_250
        need_w_500 = MPAS_field_will_be_written('w_500hPa')
        need_w = need_w .or. need_w_500
        need_any_diags = need_any_diags .or. need_w_500
        need_w_700 = MPAS_field_will_be_written('w_700hPa')
        need_w = need_w .or. need_w_700
        need_any_diags = need_any_diags .or. need_w_700
        need_w_850 = MPAS_field_will_be_written('w_850hPa')
        need_w = need_w .or. need_w_850
        need_any_diags = need_any_diags .or. need_w_850
        need_w_925 = MPAS_field_will_be_written('w_925hPa')
        need_w = need_w .or. need_w_925
        need_any_diags = need_any_diags .or. need_w_925
        need_vorticity_200 = MPAS_field_will_be_written('vorticity_200hPa')
        need_vorticity = need_vorticity .or. need_vorticity_200
        need_any_diags = need_any_diags .or. need_vorticity_200
        need_vorticity_250 = MPAS_field_will_be_written('vorticity_250hPa')
        need_vorticity = need_vorticity .or. need_vorticity_250
        need_any_diags = need_any_diags .or. need_vorticity_250
        need_vorticity_500 = MPAS_field_will_be_written('vorticity_500hPa')
        need_vorticity = need_vorticity .or. need_vorticity_500
        need_any_diags = need_any_diags .or. need_vorticity_500
        need_vorticity_700 = MPAS_field_will_be_written('vorticity_700hPa')
        need_vorticity = need_vorticity .or. need_vorticity_700
        need_any_diags = need_any_diags .or. need_vorticity_700
        need_vorticity_850 = MPAS_field_will_be_written('vorticity_850hPa')
        need_vorticity = need_vorticity .or. need_vorticity_850
        need_any_diags = need_any_diags .or. need_vorticity_850
        need_vorticity_925 = MPAS_field_will_be_written('vorticity_925hPa')
        need_vorticity = need_vorticity .or. need_vorticity_925
        need_any_diags = need_any_diags .or. need_vorticity_925
        need_t_isobaric = MPAS_field_will_be_written('t_isobaric')
        need_any_diags = need_any_diags .or. need_t_isobaric
        need_z_isobaric = MPAS_field_will_be_written('z_isobaric')
        need_any_diags = need_any_diags .or. need_z_isobaric
        need_meanT_500_300 = MPAS_field_will_be_written('meanT_500_300')
        need_any_diags = need_any_diags .or. need_meanT_500_300

        if (need_any_diags) then
            call interp_diagnostics(mesh, state, 1, diag)
        end if
   
    end subroutine isobaric_diagnostics_compute


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
       
        real (kind=RKIND) :: evp
       
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
       
        call mpas_pool_get_array(diag, 't_iso_levels', t_iso_levels)
        call mpas_pool_get_array(diag, 'z_iso_levels', z_iso_levels)
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
       
        if (need_t_isobaric) then
            t_iso_levels(1) = 30000.0
            t_iso_levels(2) = 35000.0
            t_iso_levels(3) = 40000.0
            t_iso_levels(4) = 45000.0
            t_iso_levels(5) = 50000.0
        end if
       
        if (need_z_isobaric) then
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
       end if
       
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
       
        if (NEED_TEMP .or. NEED_RELHUM .or. NEED_DEWPOINT .or. need_mslp) then
           !calculation of temperature at cell centers:
            do iCell = 1,nCells
            do k = 1,nVertLevels
                temperature(k,iCell) = (theta_m(k,iCell)/(1._RKIND+rvord*scalars(index_qv,k,iCell)))*exner(k,iCell) 

                ! Vapor pressure (NB: pressure here is already in hPa)
                evp = pressure(k,iCell) * scalars(index_qv,k,iCell) / (scalars(index_qv,k,iCell) + 0.622_RKIND)
                evp = max(evp, 1.0e-8_RKIND)

                ! Dewpoint temperature following Bolton (1980)
                dewpoint(k,iCell) = (243.5_RKIND * log(evp/6.112_RKIND)) / (17.67_RKIND - log(evp/6.112_RKIND))
                dewpoint(k,iCell) = dewpoint(k,iCell) + 273.15
            enddo
            enddo
        end if
       
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

        if (NEED_TEMP) then
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
        end if
       
       
        if (NEED_RELHUM) then
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
        end if
       
        if (NEED_DEWPOINT) then
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
        end if
       
        if (NEED_UZONAL) then
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
        end if
       
        if (NEED_UMERIDIONAL) then
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
        end if
       
        if(allocated(field_in)) deallocate(field_in)
        if(allocated(press_in)) deallocate(press_in)
       
        if (NEED_W .or. NEED_HEIGHT) then
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
       
            if(allocated(field_in)) deallocate(field_in)
            if(allocated(press_in)) deallocate(press_in)
           ! write(0,*) '--- end interpolate vertical velocity:'
        end if

        if(allocated(field_interp)) deallocate(field_interp)
        if(allocated(press_interp)) deallocate(press_interp)
       
        if (NEED_VORTICITY) then
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

            if(allocated(field_interp)) deallocate(field_interp)
            if(allocated(press_interp)) deallocate(press_interp)
       
            if(allocated(field_in    )) deallocate(field_in)
            if(allocated(press_in    )) deallocate(press_in)
        end if

        if(allocated(pressureCp1) ) deallocate(pressureCp1 )
        if(allocated(pressure_v)  ) deallocate(pressure_v  )
       
        if (need_mslp) then
            !... compute SLP (requires temp, height, pressure, qvapor)
             call compute_slp(nCells, nVertLevels, num_scalars, temperature, height, pressure, index_qv, scalars, mslp)
             mslp(:) = mslp(:) * 100.0   ! Convert from hPa to Pa
            !... alternative way
            !do iCell = 1, nCells
            !   mslp(iCell) = diag % surface_pressure % array(iCell) + 11.38*height(1,iCell)
            !   mslp(iCell) = mslp(iCell)/100.
            !enddo
        end if
    
    
        !!!!!!!!!!! Additional temperature levels for vortex tracking !!!!!!!!!!!
        if (need_t_isobaric .or. need_meanT_500_300) then
     
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
     
            if (need_t_isobaric) then
                call interp_tofixed_pressure(nCells, nVertLevels, nIsoLevelsT, press_in, field_in, press_interp, field_interp)
         
                do k=1,nIsoLevelsT
                   t_isobaric(k,1:nCells) = field_interp(1:nCells,k)
                end do
            end if
     
     
            !!!!!!!!!!! Calculate mean temperature in 500 hPa - 300 hPa layer !!!!!!!!!!!
     
            if (need_meanT_500_300) then
                call compute_layer_mean(meanT_500_300, 50000.0, 30000.0, field_in, press_in)
            end if
     
     
            deallocate(field_in)
            deallocate(field_interp)
            deallocate(press_in)
            deallocate(press_interp)
        end if
     
     
        !!!!!!!!!!! Additional height levels for vortex tracking !!!!!!!!!!!
        if (need_z_isobaric) then
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
        end if
    
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
             write(0,*) 'icol = ',icol
             write(0,*) 'klo = ',klo,', khi = ',khi, ': and they should not be equal.'
             write(0,*) '*** MSLP field will not be computed'
             slp(:) = 0.0
             return
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
   
end module isobaric_diagnostics
