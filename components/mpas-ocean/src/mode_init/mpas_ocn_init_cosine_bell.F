! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the
! LICENSE file
! distributed with this code, or at
! http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_cosine_bell
!
!> \brief MPAS ocean initialize case -- Cosine Bell
!> \author Luke Van Roekel
!> \date   09/01/2020
!> \details
!>  This module contains the routines for initializing the
!>  the cosine bell test case
!>  Reference: Section 2a in Skamarock, W.C. and A. Gassmann, 2011: Conservative
!>  Transport Schemes for Spherical Geodesic Grids: High-Order Flux Operators for
!>  ODE-Based Time Integration. Mon. Wea. Rev., 139, 2962–2975,
!>  https://doi.org/10.1175/MWR-D-10-05056.1
!
!-----------------------------------------------------------------------

module ocn_init_cosine_bell

   use mpas_kind_types
   use mpas_io_units
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_stream_manager
   use mpas_dmpar

   use ocn_init_cell_markers
   use ocn_init_vertical_grids

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_init_setup_cosine_bell, &
             ocn_init_validate_cosine_bell

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_setup_cosine_bell
!
!> \brief   Set-up for test case
!> \author  Luke Van Roekel
!> \date    09/01/2020
!> \details
!>  Reference: Section 2a in Skamarock, W.C. and A. Gassmann, 2011: Conservative
!>  Transport Schemes for Spherical Geodesic Grids: High-Order Flux Operators for
!>  ODE-Based Time Integration. Mon. Wea. Rev., 139, 2962–2975,
!>  https://doi.org/10.1175/MWR-D-10-05056.1
!
!-----------------------------------------------------------------------

   subroutine ocn_init_setup_cosine_bell(domain, iErr)!{{{

   !--------------------------------------------------------------------

      type (domain_type), intent(inout) :: domain
      integer, intent(out) :: iErr
      real (kind=RKIND) :: temperature, salinity

      type (block_type), pointer :: block_ptr

      type (mpas_pool_type), pointer :: meshPool, verticalMeshPool,statePool
      type (mpas_pool_type), pointer :: forcingPool

      type (mpas_pool_type), pointer :: tracersPool

      integer, pointer :: nVertLevels, nVertLevelsP1, nCellsSolve,nEdgesSolve, nEdge, nVerticesSolve
      integer, pointer :: index_temperature, index_salinity,index_tracer1

      integer, dimension(:), pointer :: maxLevelCell
      real (kind=RKIND), dimension(:), pointer :: refBottomDepth,refZMid
      real (kind=RKIND), dimension(:), pointer :: lonCell, latEdge, latCell, latVertex
      real (kind=RKIND), dimension(:), pointer :: bottomDepth, angleEdge
      real (kind=RKIND), dimension(:), pointer :: fCell, fEdge, fVertex,xCell, yCell
      real (kind=RKIND), dimension(:, :), pointer :: layerThickness,restingThickness,normalVelocity
      real (kind=RKIND), dimension(:, :, :), pointer :: activeTracers,debugTracers

      real (kind=RKIND), dimension(:), pointer :: interfaceLocations

      integer :: iCell, iEdge, iVertex, k, kML

      real (kind=RKIND) :: temp, velocity

      character (len=StrKIND), pointer :: config_init_configuration

      real (kind=RKIND), pointer ::config_cosine_bell_temperature, &
                                   config_cosine_bell_salinity,    &
                                   config_cosine_bell_radius,  &
                                   config_cosine_bell_lat_center, &
                                   config_cosine_bell_lon_center, &
                                   config_cosine_bell_psi0,  &
                                   config_cosine_bell_vel_pd, &
                                   sphere_radius
      ! assume no error
      iErr = 0

      ! get and test if this is the configuration specified
      call mpas_pool_get_config(domain % configs,'config_init_configuration', config_init_configuration)
      if(config_init_configuration .ne. trim('cosine_bell')) return

      ! load the remaining configuration parameters
      call mpas_pool_get_config(domain % configs,'config_cosine_bell_temperature', &
                                config_cosine_bell_temperature)
      call mpas_pool_get_config(domain % configs,'config_cosine_bell_salinity', &
                    config_cosine_bell_salinity)
      call mpas_pool_get_config(domain % configs,'config_cosine_bell_radius', &
                                config_cosine_bell_radius)
      call mpas_pool_get_config(domain % configs,'config_cosine_bell_psi0', &
                                config_cosine_bell_psi0)
      call mpas_pool_get_config(domain % configs,'config_cosine_bell_vel_pd', &
                                config_cosine_bell_vel_pd)
      call mpas_pool_get_config(domain % configs,'config_cosine_bell_lat_center', &
                                config_cosine_bell_lat_center)
      call mpas_pool_get_config(domain % configs,'config_cosine_bell_lon_center', &
                                config_cosine_bell_lon_center)
      !load data that required to initialize the ocean simulation
      block_ptr => domain % blocklist
      do while(associated(block_ptr))
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh',meshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh',verticalMeshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'state',statePool)
        call mpas_pool_get_subpool(block_ptr % structs, 'forcing',forcingPool)

        call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

        call mpas_pool_get_dimension(meshPool, 'nVertLevels',nVertLevels)
        call mpas_pool_get_dimension(meshPool, 'nCellsSolve',nCellsSolve)
        call mpas_pool_get_dimension(meshPool, 'nEdge', nEdge)
        call mpas_pool_get_dimension(meshPool, 'nEdgesSolve',nEdgesSolve)
        call mpas_pool_get_dimension(meshPool, 'nVerticesSolve',nVerticesSolve)
        call mpas_pool_get_config(meshPool, 'sphere_radius', sphere_radius)

        call mpas_pool_get_dimension(tracersPool, 'index_temperature',index_temperature)
        call mpas_pool_get_dimension(tracersPool, 'index_salinity',index_salinity)
        call mpas_pool_get_dimension(tracersPool, 'index_tracer1',index_tracer1)

        call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
        call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
        call mpas_pool_get_array(meshPool, 'angleEdge', angleEdge)
        call mpas_pool_get_array(meshPool, 'latCell', latCell)
        call mpas_pool_get_array(meshPool, 'latEdge', latEdge)
        call mpas_pool_get_array(meshPool, 'latVertex', latVertex)
        call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
        call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
        call mpas_pool_get_array(meshPool, 'fCell', fCell)
        call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
        call mpas_pool_get_array(meshPool, 'fVertex', fVertex)
        call mpas_pool_get_array(meshPool, 'xCell', xCell)
        call mpas_pool_get_array(meshPool, 'yCell', yCell)


        call mpas_pool_get_array(verticalMeshPool, 'refZMid', refZMid)
        call mpas_pool_get_array(verticalMeshPool, 'restingThickness',restingThickness)

        call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, 1)
        call mpas_pool_get_array(tracersPool, 'activeTracers',activeTracers, 1)
        call mpas_pool_get_array(tracersPool, 'debugTracers',debugTracers, 1)
        call mpas_pool_get_array(statePool, 'layerThickness',layerThickness, 1)

        ! Set refBottomDepth and refBottomDepthTopOfCell
        do k = 1, nVertLevels
           refBottomDepth(k) = 100.0_RKIND
           refZMid(k) = 50.0_RKIND
        end do

        do iCell = 1, nCellsSolve
           if(associated(activeTracers) ) then

           ! Loop from surface through surface layer depth

           activeTracers(index_temperature, :, iCell) = config_cosine_bell_temperature
           activeTracers(index_salinity, :, iCell) = config_cosine_bell_salinity
          endif

        if(associated(debugTracers)) then
          temp = sphere_radius*acos(sin(config_cosine_bell_lat_center)*sin(latCell(iCell)) + &
                  cos(config_cosine_bell_lat_center)*cos(latCell(iCell))*cos(lonCell(iCell) - &
                  config_cosine_bell_lon_center))

          if( temp < config_cosine_bell_radius ) then
            debugTracers(:,:,iCell) = config_cosine_bell_psi0/2.0_RKIND * ( 1.0_RKIND + &
              cos(3.1415926_RKIND*temp/config_cosine_bell_radius))
          else
            debugTracers(:,:,iCell) = 0.0_RKIND
          endif
        endif
        ! Set layerThickness
        do k = 1, nVertLevels
           layerThickness(k, iCell) = 100.0_RKIND
           restingThickness(k, iCell) = layerThickness(k, iCell)
        end do

        ! Set Coriolis parameter
        fCell(iCell) = 0.0_RKIND*omega*sin(latCell(iCell))

        ! Set bottomDepth
        bottomDepth(iCell) = 100.0_RKIND

        ! Set maxLevelCell
        maxLevelCell(iCell) = nVertLevels

     end do  ! do iCell

     do iEdge = 1, nEdgesSolve
        fEdge(iEdge) = 0.0_RKIND*omega*sin(latEdge(iEdge))
     end do

     do iVertex=1, nVerticesSolve
        fVertex(iVertex) = 0.0_RKIND*omega*sin(latVertex(iVertex))
     end do

     !set normal velocity
     do iEdge = 1,nEdgesSolve
       velocity = 2.0*3.14159265*sphere_radius*cos(latEdge(iEdge)) / &
         (86400.0_RKIND*config_cosine_bell_vel_pd)
       normalVelocity(:,iEdge) = velocity*cos(angleEdge(iEdge))
     enddo

     block_ptr => block_ptr % next
   end do

  !--------------------------------------------------------------------

   end subroutine ocn_init_setup_cosine_bell!}}}

!***********************************************************************
!
!  routine ocn_init_validate_cosine_bell
!
!> \brief   Validation for test case
!> \author  Luke Van Roekel
!> \date    09/01/2020
!> \details
!
!-----------------------------------------------------------------------

   subroutine ocn_init_validate_cosine_bell(configPool, packagePool,iocontext, iErr)!{{{

   !--------------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: configPool
      type (mpas_pool_type), intent(inout) :: packagePool
      type (mpas_io_context_type), intent(inout) :: iocontext

      integer, intent(out) :: iErr

      character (len=StrKIND), pointer :: config_init_configuration
      integer, pointer :: config_vert_levels

      iErr = 0

      call mpas_pool_get_config(configPool, 'config_init_configuration',config_init_configuration)

      if(config_init_configuration .ne. trim('cosine_bell')) return

      call mpas_pool_get_config(configPool, 'config_vert_levels',config_vert_levels)

      config_vert_levels = 3

   !--------------------------------------------------------------------

   end subroutine ocn_init_validate_cosine_bell!}}}

!***********************************************************************

end module ocn_init_cosine_bell

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
