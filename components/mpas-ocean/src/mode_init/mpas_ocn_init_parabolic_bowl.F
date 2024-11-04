! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_parabolic_bowl
!
!> \brief Initialize the parabolic bowl test case
!> \author S. Brus, D. Wirasaet  
!> \date   May-June 2022
!> \details
!>  This module contains the routines for initializing the
!>  parabolic bowl initial condition (Thacker's problem)
!>
!
!-----------------------------------------------------------------------

module ocn_init_parabolic_bowl

   use mpas_kind_types
   use mpas_io_units
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_dmpar
   use mpas_vector_operations

   use ocn_constants
   use ocn_config
   use ocn_init_vertical_grids
   use ocn_init_cell_markers
   use ocn_subgrid
   use ocn_init_subgrid

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

   public :: ocn_init_setup_parabolic_bowl, &
             ocn_init_validate_parabolic_bowl

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   real (kind=RKIND):: LL, CC
   real (kind=RKIND):: oneMC2, sqrtOneMC2, oneMC
   real (kind=RKIND),  parameter :: eps=1.0e-10

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_setup_parabolic_bowl
!
!> \brief   Setup for this initial condition
!> \author  D. Wirasaet
!> \date    May-June 2022
!> \details
!>  This routine sets up the initial conditions for this case.
!>   To be run in sigma vertical coordinates and single-layer
!
!-----------------------------------------------------------------------

  subroutine ocn_init_setup_parabolic_bowl(domain, iErr)!{{{
    
    implicit none
    !--------------------------------------------------------------------

    type (domain_type), intent(inout) :: domain
    integer, intent(out) :: iErr

    type (block_type), pointer :: block_ptr
    type (mpas_pool_type), pointer :: meshPool
    type (mpas_pool_type), pointer :: statePool
    type (mpas_pool_type), pointer :: tracersPool
    type (mpas_pool_type), pointer :: verticalMeshPool

    ! local variables
    integer :: iCell, iEdge, iVertex, k, idx
    real (kind=RKIND) :: yMin, yMax, xMin, xMax, dcEdgeMin, dcEdgeMinGlobal
    real (kind=RKIND) :: yMinGlobal, yMaxGlobal, yMidGlobal, xMinGlobal, xMaxGlobal
    real (kind=RKIND) :: localVar1, localVar2
    real (kind=RKIND), dimension(:), pointer :: interfaceLocations

    ! Define dimension pointers
    integer, pointer :: nCellsSolve, nEdgesSolve, nVerticesSolve, nVertLevels, nVertLevelsP1
    integer, pointer :: index_temperature, index_salinity
    integer, pointer :: maxEdges

    ! Define variable pointers
    logical, pointer :: on_a_sphere
    integer, dimension(:), pointer :: minLevelCell, maxLevelCell
    integer, dimension(:), pointer :: nEdgesOnCell
    integer, dimension(:,:), pointer :: verticesOnCell, verticesOnEdge
    integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnVertex
    real (kind=RKIND), dimension(:), pointer :: xCell, yCell, refBottomDepth, refZMid, &
         vertCoordMovementWeights, bottomDepth, fCell, fEdge, fVertex, dcEdge
    real (kind=RKIND), dimension(:,:), pointer:: zMid

    real (kind=RKIND), dimension(:), pointer:: xEdge, yEdge, xVertex, yVertex
    real (kind=RKIND) :: minBottomDepth, maxBottomDepth, globalMaxBottomDepth, globalMinBottomDepth
    real (kind=RKIND), dimension(:,:), pointer :: layerThickness, restingThickness
    real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers
    
    real (kind=RKIND), dimension(:), pointer :: ssh
    real (kind=RKIND), dimension(:), pointer :: areaCell
    real (kind=RKIND), dimension(:,:), pointer :: edgeNormalVectors
    real (kind=RKIND), dimension(:,:), pointer :: normalVelocity

  
    real (kind=RKIND):: HH, uu, vv
    real (kind=RKIND):: RR, num, den
    real (kind=RKIND):: xshift = 0.0, yshift = 0.0 
    real (kind=RKIND) :: layerThicknessEdgeAverage
    real (kind=RKIND), dimension(:,:), allocatable :: rSubgridPoints, sSubgridPoints
    real (kind=RKIND), dimension(:), allocatable :: subgridBathymetryValues, subgridAreas
    real (kind=RKIND), dimension(:), allocatable :: subgridSshValues
    real (kind=RKIND), dimension(:), allocatable :: subgridUValues, subgridVValues
    real (kind=RKIND), dimension(:), allocatable :: uVelocityAverage, vVelocityAverage
    integer :: nSubgridCell, nSubgridEdge, nSubgridVertex
    integer :: nSubgridTriPerSlice
    integer :: v1, v2
    integer :: c1, c2
    real (kind=RKIND) :: x(3), y(3)
    integer :: slice, nSlice
    integer :: i,j
    real (kind=RKIND) :: deltaZ

    integer:: nsubgridCellEdge, iEdgeSegment
    real (kind=RKIND), dimension(:,:), allocatable :: cellEdgeBathymetryValues
    real (kind=RKIND), dimension(:), allocatable:: dsEdge
    integer:: jj
    real (kind=RKIND), dimension(:), allocatable:: xSubgridCell, ySubgridCell

    iErr = 0

    if(config_init_configuration .ne. trim('parabolic_bowl')) return

    ! Determine vertical grid for configuration
    call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
    call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
    call mpas_pool_get_dimension(meshPool, 'nVertLevelsP1', nVertLevelsP1)
    call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)

    nVertLevels  = config_parabolic_bowl_vert_levels ;
    nVertLevelsP1 = nVertLevels + 1

    allocate(interfaceLocations(nVertLevelsP1))
    call ocn_generate_vertical_grid( config_vertical_grid, interfaceLocations, ocnConfigs ) ;

    ! Initalize min/max values to large positive and negative values
    yMin = 1.0E10_RKIND
    yMax = -1.0E10_RKIND
    xMin = 1.0E10_RKIND
    xMax = -1.0E10_RKIND
    dcEdgeMin = 1.0E10_RKIND

    ! Determine local min and max values.
    block_ptr => domain % blocklist
    do while(associated(block_ptr))
       call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)

       call mpas_pool_get_dimension( meshPool, 'nCellsSolve', nCellsSolve )
       call mpas_pool_get_dimension( meshPool, 'nEdgesSolve', nEdgesSolve )

       call mpas_pool_get_array(meshPool, 'xCell', xCell)
       call mpas_pool_get_array(meshPool, 'yCell', yCell)
       call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)

       yMin = min( yMin, minval(yCell(1:nCellsSolve)))
       yMax = max( yMax, maxval(yCell(1:nCellsSolve)))
       xMin = min( xMin, minval(xCell(1:nCellsSolve)))
       xMax = max( xMax, maxval(xCell(1:nCellsSolve)))
       dcEdgeMin = min( dcEdgeMin, minval(dcEdge(1:nEdgesSolve)))

       block_ptr => block_ptr % next
    end do

    ! Determine global min and max values.
    call mpas_dmpar_min_real(domain % dminfo, yMin, yMinGlobal)
    call mpas_dmpar_max_real(domain % dminfo, yMax, yMaxGlobal)
    call mpas_dmpar_min_real(domain % dminfo, xMin, xMinGlobal)
    call mpas_dmpar_max_real(domain % dminfo, xMax, xMaxGlobal)
    call mpas_dmpar_min_real(domain % dminfo, dcEdgeMin, dcEdgeMinGlobal)

    ! Calculate constants
    xshift = 0.5_RKIND*(xMin + xMax) ;
    yshift = 0.5_RKIND*(YMin + yMax) ;  
    
    den = ((config_parabolic_bowl_omega**2.0) - (config_parabolic_bowl_Coriolis_parameter**2.0)) ; 
    num = (8.0_RKIND*config_parabolic_bowl_gravity*config_parabolic_bowl_b0) ;
    LL  = sqrt(num/den) ; 

    den = (config_parabolic_bowl_b0 + config_parabolic_bowl_eta0)**2.0 + config_parabolic_bowl_b0**2.0 ;
    num = (config_parabolic_bowl_b0 + config_parabolic_bowl_eta0)**2.0 - config_parabolic_bowl_b0**2.0 ; 

    CC = num/den ;

    oneMC2 = 1.0_RKIND - CC**2.0 ;
    sqrtOneMC2 = sqrt(oneMC2) ; 
    oneMC = 1.0_RKIND - CC ;


    block_ptr => domain % blocklist
    call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
    call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
    call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)
    call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

    call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels) ;
    call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve) ;
    call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve) ; 
    call mpas_pool_get_dimension(meshPool, 'nVerticesSolve', nVerticesSolve) ; 
    call mpas_pool_get_dimension(meshPool, 'maxEdges', maxEdges)
    
    call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_temperature)
    call mpas_pool_get_dimension(tracersPool, 'index_salinity', index_salinity)

    call mpas_pool_get_array(meshPool, 'xCell', xCell)
    call mpas_pool_get_array(meshPool, 'yCell', yCell)
    call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
    call mpas_pool_get_array(meshPool, 'vertCoordMovementWeights', vertCoordMovementWeights)
    call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
    call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)
    call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
    call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
    call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

    call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
    call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
    call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
    call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)

    call mpas_pool_get_array(meshPool, 'fCell', fCell)
    call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
    call mpas_pool_get_array(meshPool, 'fVertex', fVertex)
    
    call mpas_pool_get_array(meshPool, 'xEdge', xEdge ) 
    call mpas_pool_get_array(meshPool, 'yEdge', yEdge ) 
    call mpas_pool_get_array(meshPool, 'xVertex', xVertex ) 
    call mpas_pool_get_array(meshPool, 'yVertex', yVertex ) 

    call mpas_pool_get_array(statePool, 'zMid', zMid, 1) ;

    call mpas_pool_get_array(statePool, 'ssh', ssh, 1)
    call mpas_pool_get_array(meshPool, 'edgeNormalVectors', edgeNormalVectors ) ;
    call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity ) ;

    call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
    call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)

    call mpas_pool_get_array(verticalMeshPool, 'refZMid', refZMid)
    call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)


    ! if config_parabolic_bowl_adjust_domain_center == .true.,
    ! Adjust center of the mesh so that its center is located at (0,0)
    if ( config_parabolic_bowl_adjust_domain_center ) then
       xCell = xCell - xshift ; 
       yCell = yCell - yshift ;

       xEdge = xEdge - xshift ;
       yEdge = yEdge - yshift ;

       xVertex = xVertex - xshift ;
       yVertex = yVertex - yshift ; 
    end if


    ! Initlialze vector
    call mpas_initialize_vectors(meshPool) ; 


    minLevelCell(:) = 1
    do iCell = 1, nCellsSolve
      ! Set up vertical grid
      maxLevelCell(iCell) = nVertLevels ; ! sigma coordinates 
    end do


    do iCell = 1, nCellsSolve
      
       ! Set temperature
       activeTracers(index_temperature, :, iCell) = 10.0_RKIND

       ! Set salinity
       activeTracers(index_salinity, :, iCell) = 30.0_RKIND

      ! Set Coriolis parameters, if other than zero
      fCell(iCell) = config_parabolic_bowl_coriolis_parameter ; 
    end do

    do iEdge = 1, nEdgesSolve
      fEdge(iEdge) = config_parabolic_bowl_coriolis_parameter ;  
    end do

    do iVertex = 1, nVerticesSolve
      fVertex(iVertex) =  config_parabolic_bowl_coriolis_parameter ;
    end do

    allocate(uVelocityAverage(nEdgesSolve))
    allocate(vVelocityAverage(nEdgesSolve))
    
    if (config_use_subgrid_wetting_drying) then

       call ocn_subgrid_init(domain,iErr)
       call ocn_init_subgrid_calculations(domain, &
                                          ocn_init_parabolic_bowl_bathymetry, &
                                          ocn_init_parabolic_bowl_velocity, &
                                          ocn_init_parabolic_bowl_ssh, &
                                          config_parabolic_bowl_subgrid_refinement_level, &
                                          config_parabolic_bowl_subgrid_edge_bathymetry_max_pixel, &
                                          config_parabolic_bowl_subgrid_use_thin_layer, &
                                          uVelocityAverage, &
                                          vVelocityAverage, &
                                          iErr)
    end if

    ! Find max bottom depth
    maxBottomDepth = maxval( bottomDepth ) ; 
    minBottomDepth = minval( bottomDepth ) ; 
    call mpas_dmpar_max_real( domain % dminfo, maxBottomDepth, globalMaxBottomDepth ) ;
    call mpas_dmpar_min_real( domain % dminfo, minBottomDepth, globalMinBottomDepth ) ;

    ! Set refBottomDepth and refZMid
    do k = 1, nVertLevels
      refBottomDepth(k) = globalMaxBottomDepth*interfaceLocations(k+1) ; 
      refZMid(k) = -0.5_RKIND*( interfaceLocations(k+1) + interfaceLocations(k))*globalMaxBottomDepth ; 
    end do

    ! Set vertCoordMovementWeights
    vertCoordMovementWeights(:) = 1.0_RKIND

    ! Set velocity
    do iEdge = 1, nEdgesSolve

      if (config_use_subgrid_wetting_drying) then
         do k = 1, nVertLevels 
            normalVelocity(k,iEdge) = uVelocityAverage(iEdge)*edgeNormalVectors(1,iEdge) &
                                    + vVelocityAverage(iEdge)*edgeNormalVectors(2,iEdge) ;
         end do
      else
         call ocn_init_parabolic_bowl_velocity(xEdge(iEdge), yEdge(iEdge), uu, vv)
         do k = 1, nVertLevels 
            normalVelocity(k,iEdge) = uu*edgeNormalVectors(1,iEdge) + vv*edgeNormalVectors(2,iEdge) ;
         end do
      end if
    end do

    ! Set layer thickness and ssh
    if (config_use_wetting_drying) then

      do iCell = 1, nCellsSolve
        ! Set up vertical grid
        maxLevelCell(iCell) = nVertLevels ; ! sigma coordinates 
       end do

       do iCell = 1, nCellsSolve
         !
         ! make sure depth is thick enough via ssh = TOTAL_DEPTH - bottomDepth
         ! add a thin layer of nlayer*config_drying_min_cellhight 
         ! 
 
         if (config_use_subgrid_wetting_drying) then

           call ocn_subgrid_ssh_lookup(layerThickness(1,iCell),&
                                       subgridWetVolumeCellTable(:,iCell),&
                                       subgridSshCellTableRange(:,iCell),&
                                       bottomDepth(iCell),&
                                       subgridCellBathymetryMin(iCell),&
                                       ssh(iCell))
           !call ocn_subgrid_layer_thickness_lookup(ssh(iCell), & 
           !                                subgridWetVolumeCellTable(:,iCell), &
           !                                subgridSshCellTableRange(:,iCell),&
           !                                bottomDepth(iCell),&
           !                                LayerThickness(1,iCell))

         else
           call ocn_init_parabolic_bowl_bathymetry(xCell(iCell),yCell(iCell),bottomDepth(iCell))
           call ocn_init_parabolic_bowl_ssh(xCell(iCell),yCell(iCell),bottomDepth(iCell),ssh(iCell))
           ssh(iCell) = - bottomDepth(iCell) + &
             max(ssh(iCell) + bottomDepth(iCell), &
             maxLevelCell(iCell)*(config_drying_min_cell_height + eps))

           do k = 1, maxLevelCell(iCell)
             layerThickness(k,iCell) = max(config_drying_min_cell_height + eps, &
               (ssh(iCell) + bottomDepth(iCell))/real(maxLevelCell(iCell),RKIND))
            
             if (layerThickness(k,iCell) < config_drying_min_cell_height) then
               call mpas_log_write('layerThickness($i,$i)=$r', MPAS_LOG_CRIT, &
                 intArgs=(/k,iCell/), &
                 realArgs=(/layerThickness(k,iCell)/))
             end if
           end do
         endif
       
       
         do k = 1, maxLevelCell(iCell)
           restingThickness(k,iCell) = bottomDepth(iCell)/maxLevelCell(iCell)
         end do
       end do
       
    end if 

    if (config_use_subgrid_wetting_drying) then
    do iCell = 1,nCellsSolve
      call ocn_subgrid_layer_thickness_lookup(ssh(iCell), & 
                                              subgridWetVolumeCellTable(:,iCell), &
                                              subgridSshCellTableRange(:,iCell),&
                                              bottomDepth(iCell),&
                                              subgridLayerThicknessDebug(iCell))
    enddo
    endif

    deallocate(interfaceLocations)
    !--------------------------------------------------------------------

    print*, "****** End parabolic bowl init *****" ;

    return ;  
  end subroutine ocn_init_setup_parabolic_bowl!}}}


!***********************************************************************
!
!  routine ocn_init_parabolic_bowl_bathymetry
!
!> \brief   Compute bathymetry 
!> \author  Steven Brus, D. Wirasaet
!> \date    November 2022
!> \details Return the value of the bathymetry at a given x,y point
!>  
!
!-----------------------------------------------------------------------

   subroutine ocn_init_parabolic_bowl_bathymetry(x, y, depth)!{{{
           
      implicit none 

      real (kind=RKIND), intent(in) :: x, y
      real (kind=RKIND), intent(out) :: depth
      real (kind=RKIND) :: RR

   !--------------------------------------------------------------------

      RR = sqrt(x**2 + y**2) 
      depth  = config_parabolic_bowl_b0*(1.0_RKIND - (RR**2)/(LL**2))

   !--------------------------------------------------------------------

   end subroutine ocn_init_parabolic_bowl_bathymetry!}}}

!***********************************************************************
!
!  routine ocn_init_parabolic_bowl_ssh
!
!> \brief   Compute initial ssh field
!> \author  Steven Brus, D. Wirasaet
!> \date    November 2022
!> \details Use exact solution to compute ssh field for initial conditions
!>  
!
!-----------------------------------------------------------------------

   subroutine ocn_init_parabolic_bowl_ssh(x, y, bottomDepth, ssh)!{{{
           
      implicit none 

      real (kind=RKIND), intent(in) :: x, y
      real (kind=RKIND), intent(in) :: bottomDepth
      real (kind=RKIND), intent(out) :: ssh
      real (kind=RKIND) :: RR

   !--------------------------------------------------------------------

      RR = sqrt(x**2 + y**2) ;
 
      ssh = (sqrtOneMC2/oneMC) - 1.0_RKIND - ((RR**2)/(LL**2))*( (oneMC2/(oneMC**2)) - 1.0_RKIND ) ;
      ssh = config_parabolic_bowl_b0*ssh ; 
      !ssh = - bottomDepth + max(ssh + bottomDepth, config_drying_min_cell_height + eps)

   !--------------------------------------------------------------------

   end subroutine ocn_init_parabolic_bowl_ssh!}}}

!***********************************************************************
!
!  routine ocn_init_parabolic_bowl_velocity
!
!> \brief   Compute initial velocity field
!> \author  Steven Brus, D. Wirasaet
!> \date    November 2022
!> \details Use exact solution to comupte velocity field for initial conditions
!>  
!
!-----------------------------------------------------------------------

   subroutine ocn_init_parabolic_bowl_velocity(x, y, u, v)!{{{
           
      implicit none 

      real (kind=RKIND), intent(in) :: x, y
      real (kind=RKIND), intent(out) :: u, v 
      real (kind=RKIND) :: RR, HH

   !--------------------------------------------------------------------

   RR = sqrt(x**2 + y**2) ;
   
   HH = (sqrtOneMC2/oneMC) - ( (RR**2)/(LL**2))*( (oneMC2/(oneMC**2))  ) ;
   HH = config_parabolic_bowl_b0*HH  ;
    
   u =  - config_parabolic_bowl_Coriolis_parameter*y*( sqrtOneMC2 + CC - 1.0_RKIND  ) ;
   v =    config_parabolic_bowl_Coriolis_parameter*x*( sqrtOneMc2 + CC - 1.0_RKIND  ) ; 

   u = 0.5_RKIND*u/oneMC ; 
   v = 0.5_RKIND*v/oneMC ; 
   if ( HH < eps ) then
       u = 0.0_RKIND ;
       v = 0.0_RKIND ;
   end if

   !--------------------------------------------------------------------

   end subroutine ocn_init_parabolic_bowl_velocity!}}}

!***********************************************************************
!
!  routine ocn_init_validate_parabolic_bowl
!
!> \brief   Validation for this initial condition
!> \author  Steven Brus
!> \date   November 2022
!> \details
!>  This routine validates the configuration options for this case.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_validate_parabolic_bowl(configPool, packagePool, iocontext, iErr)!{{{
           
      implicit none 

   !--------------------------------------------------------------------
      type (mpas_pool_type), intent(inout) :: configPool, packagePool
      type (mpas_io_context_type), intent(inout) :: iocontext

      integer, intent(out) :: iErr

      character (len=StrKIND), pointer :: config_init_configuration
      integer, pointer :: config_vert_levels, config_parabolic_bowl_vert_levels
      integer, pointer :: config_subgrid_table_levels, config_parabolic_bowl_subgrid_table_levels
   
      iErr = 0

      call mpas_pool_get_config(configPool, 'config_init_configuration', config_init_configuration)

      if(config_init_configuration .ne. trim('parabolic_bowl')) return

      call mpas_pool_get_config(configPool, 'config_vert_levels', config_vert_levels)
      call mpas_pool_get_config(configPool, 'config_parabolic_bowl_vert_levels', config_parabolic_bowl_vert_levels)

      if(config_vert_levels <= 0 .and. config_parabolic_bowl_vert_levels > 0) then
         config_vert_levels = config_parabolic_bowl_vert_levels
      else if (config_vert_levels <= 0) then
         call mpas_log_write( 'Validation failed for parabolic_bowl. Not given a usable value for vertical levels.', MPAS_LOG_CRIT)
         iErr = 1
      end if

      call mpas_pool_get_config(configPool, 'config_subgrid_table_levels', config_subgrid_table_levels)
      call mpas_pool_get_config(configPool, 'config_parabolic_bowl_subgrid_table_levels', config_parabolic_bowl_subgrid_table_levels)

      if (config_subgrid_table_levels <= 0 .and. config_parabolic_bowl_subgrid_table_levels > 0) then
         config_subgrid_table_levels = config_parabolic_bowl_subgrid_table_levels  
      else if (config_subgrid_table_levels <= 0) then
         call mpas_log_write( 'Validation failed for parabolic_bowl. Not given a usable value for subgrid table levels.', MPAS_LOG_CRIT)
         iErr = 1
      end if

      print*, "in ocn_init_validate_parabolic_bowl()" ; 
   !--------------------------------------------------------------------

   end subroutine ocn_init_validate_parabolic_bowl!}}}


!***********************************************************************

end module ocn_init_parabolic_bowl

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
