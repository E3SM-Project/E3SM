! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_subgrid
!
!> \brief   MPAS ocean subgrid wetting and drying
!> \authors Steven Brus, Damrongsak Wirasaet
!> \date    September 2022
!> \details
!>  This module contains routines for subgrid wetting and drying
!>  These routines primarily determine the relationship between the
!>  wet volume per unit area (which can be thought of a as a layer thickness
!>  for a completely wet cell) of a control volume and a the ssh
!>  defined over the wet fraction of that control volume. When a
!>  control volume is wet, this relationship is linear, i.e. H = \eta + b.
!>  However, for partially wet cells, the nonlinear relationship between
!>  these values can be pre-computed based on subgrid scale high-resolution data.
!>  The routines in this module are responsible for providing interpolated
!>  values for the forward and inverse relationship to be used in
!>  computing the tendency terms from pre-computed look-up tables.
!>
!>  See other details in Kennedy et al. (2019):
!>  https://doi.org/10.1016/j.ocemod.2019.101491
!
!-----------------------------------------------------------------------

module ocn_subgrid

   use mpas_kind_types
   use mpas_constants
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timekeeping
   use mpas_timer
   use ocn_constants
   use ocn_config
   use ocn_mesh
   use ocn_diagnostics_variables

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   real(kind=RKIND), public, dimension(:,:), pointer :: &
      subgridWetVolumeCellTable, &
      subgridWetVolumeEdgeTable, &
      subgridWetVolumeVertexTable, &
      subgridWetFractionCellTable, &
      subgridWetFractionEdgeTable, &
      subgridWetFractionVertexTable, &
      subgridSshCellTableRange, &
      subgridSshEdgeTableRange, &
      subgridSshVertexTableRange
   real(kind=RKIND), public, dimension(:), pointer :: &
      subgridEdgeBathymetryMean, &
      subgridVertexBathymetryMean, &
      subgridCellBathymetryMin, &
      subgridEdgeBathymetryMin, &
      subgridVertexBathymetryMin, &
      subgridLayerThicknessDebug

   integer, public, pointer :: nSubgridTableLevels

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_subgrid_layerThickEdgeFlux_center, &
             ocn_subgrid_layerThickEdgeFlux_upwind, &
             ocn_subgrid_vorticity

   public :: ocn_subgrid_layer_thickness_lookup, &
             ocn_subgrid_wet_fraction_lookup, &
             ocn_subgrid_ssh_lookup, &
             ocn_subgrid_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------



contains

!***********************************************************************
!
!  routine ocn_subgrid_layerThickEdgeFlux_center
!
!> \brief   Copmutes centered layer thickness edge flux using subgrid info
!> \author  Steven Brus, Damrongsak Wirasaet
!> \date    September 2022
!> \details
!>
!>
!
!-----------------------------------------------------------------------

   subroutine ocn_subgrid_layerThickEdgeFlux_center(ssh, layerThickEdgeFlux)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(in) :: ssh

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: layerThickEdgeFlux

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: iEdge
      integer :: cell1, cell2
      real (kind=RKIND) :: sshMean
      real(kind=RKIND):: eps = 1.0e-10_RKIND


      do iEdge = 1, nEdgesAll
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         sshMean = 0.5_RKIND * (ssh(cell1) + ssh(cell2))
         call ocn_subgrid_layer_thickness_lookup(sshMean, &
                                                 subgridWetVolumeEdgeTable(:,iEdge), &
                                                 subgridSshEdgeTableRange(:,iEdge), &
                                                 subgridEdgeBathymetryMean(iEdge), &
                                                 layerThickEdgeFlux(1,iEdge) )

         if ( layerThickEdgeFlux(1,iEdge) < eps ) then
             layerThickEdgeFlux(1,iEdge) = layerThickEdgeMean(1,iEdge) ;
         end if
      end do

   end subroutine ocn_subgrid_layerThickEdgeFlux_center!}}}

!***********************************************************************
!
!  routine ocn_subgrid_layerThickEdgeFlux_upwind
!
!> \brief   Copmutes upwind layer thickness edge flux using subgrid info
!> \author  Steven Brus, Damrongsak Wirasaet
!> \date    September 2022
!> \details
!>
!>
!
!-----------------------------------------------------------------------

   subroutine ocn_subgrid_layerThickEdgeFlux_upwind(ssh, normalVelocity, layerThickness, layerThickEdgeFlux)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(in) :: ssh
      real (kind=RKIND), dimension(:,:), intent(in) :: normalVelocity
      real (kind=RKIND), dimension(:,:), intent(in) :: layerThickness

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: layerThickEdgeFlux

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: iEdge
      integer :: cell1, cell2
      real (kind=RKIND) :: sshMean
      real(kind=RKIND):: eps = 1.0e-10_RKIND

      do iEdge = 1, nEdgesAll
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         if (normalVelocity(1,iEdge) > 0.0_RKIND) then
            call ocn_subgrid_layer_thickness_lookup(ssh(cell1), &
                                                    subgridWetVolumeEdgeTable(:,iEdge), &
                                                    subgridSshEdgeTableRange(:,iEdge), &
                                                    subgridEdgeBathymetryMean(iEdge), &
                                                    layerThickEdgeFlux(1,iEdge))
         elseif (normalVelocity(1,iEdge) < 0.0_RKIND) then
            call ocn_subgrid_layer_thickness_lookup(ssh(cell2), &
                                                    subgridWetVolumeEdgeTable(:,iEdge), &
                                                    subgridSshEdgeTableRange(:,iEdge), &
                                                    subgridEdgeBathymetryMean(iEdge), &
                                                    layerThickEdgeFlux(1,iEdge))
         else
            sshMean = 0.5_RKIND*(ssh(cell1) + ssh(cell2))
            call ocn_subgrid_layer_thickness_lookup(sshMean, &
                                                    subgridWetVolumeEdgeTable(:,iEdge), &
                                                    subgridSshEdgeTableRange(:,iEdge), &
                                                    subgridEdgeBathymetryMean(iEdge), &
                                                    layerThickEdgeFlux(1,iEdge))
         end if

         if ( layerThickEdgeFlux(1,iEdge) < eps ) then
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)

            layerThickEdgeFlux(1,iEdge) = 0.5_RKIND * &
                                     ( layerThickness(1,cell1) + &
                                       layerThickness(1,cell2) )
         end if
      end do

   end subroutine ocn_subgrid_layerThickEdgeFlux_upwind!}}}

!***********************************************************************
!
!  routine ocn_subgrid_vorticity
!
!> \brief   Copmutes vorticity using subgrid info
!> \author  Steven Brus, Damrongsak Wirasaet
!> \date    September 2022
!> \details
!>
!>
!
!-----------------------------------------------------------------------

   subroutine ocn_subgrid_vorticity(ssh, normalizedRelativeVorticityVertex, normalizedPlanetaryVorticityVertex)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(in) :: ssh

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: normalizedRelativeVorticityVertex
      real (kind=RKIND), dimension(:,:), intent(out) :: normalizedPlanetaryVorticityVertex

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: iVertex, vertex
      integer :: i, k
      real (kind=RKIND) :: sshVertex
      real (kind=RKIND) :: invAreaTri1
      real (kind=RKIND) :: layerThicknessVertex

      do iVertex = 1, nVerticesHalo(2)
         invAreaTri1 = 1.0_RKIND / areaTriangle(iVertex)
         do k = 1, maxLevelVertexBot(iVertex)
            sshVertex = 0.0_RKIND
            do i = 1, vertexDegree
               sshVertex = sshVertex + ssh(cellsOnVertex(i,iVertex)) &
                                    * kiteAreasOnVertex(i,iVertex)
            end do
            sshVertex = sshVertex * invAreaTri1

            call ocn_subgrid_layer_thickness_lookup(sshVertex, &
                                                    subgridWetVolumeVertexTable(:,iVertex), &
                                                    subgridSshVertexTableRange(:,iVertex), &
                                                    subgridVertexBathymetryMean(iVertex), &
                                                    layerThicknessVertex)
            if (layerThicknessVertex == 0) cycle

            normalizedRelativeVorticityVertex(k,iVertex) = relativeVorticity(k,iVertex) / layerThicknessVertex
            normalizedPlanetaryVorticityVertex(k,iVertex) = fVertex(iVertex) / layerThicknessVertex
         end do
      end do

   end subroutine ocn_subgrid_vorticity!}}}

!***********************************************************************
!
!  routine ocn_subgrid_layer_thickness_lookup
!
!> \brief   Forward subgrid lookup (ssh -> layerThicknes)
!> \author  Steven Brus, Damrongsak Wirasaet
!> \date    September 2022
!> \details Returns the wet volume per unit area (layerThick) for a
!>          given ssh (defined over the wet-fraction based on subgrid
!>          lookup table information
!
!-----------------------------------------------------------------------

   subroutine ocn_subgrid_layer_thickness_lookup(zeta, &
                                                 subgridTable, &
                                                 subgridTableRange, &
                                                 bathymetry, &
                                                 layerThick)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), intent(in) :: zeta
      real (kind=RKIND), dimension(:), intent(in) :: subgridTable
      real (kind=RKIND), dimension(:), intent(in) :: subgridTableRange
      real (kind=RKIND), intent(in) :: bathymetry

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), intent(out) :: layerThick

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: lev

      real (kind=RKIND) :: tableMin
      real (kind=RKIND) :: tableMax
      real (kind=RKIND) :: deltaZ
      real (kind=RKIND) :: zeta0, zeta1

      tableMin = subgridTableRange(1)
      tableMax = subgridTableRange(2)
      deltaZ = (tableMax - tableMin)/real(nSubgridTableLevels-1,RKIND)

      if (zeta >= tableMax) then
        layerThick = zeta + bathymetry
      else if (zeta <= tableMin) then
        layerThick = 0.0_RKIND
      else
        do lev = 1, nSubgridTableLevels-1
          zeta0 = (real(lev,RKIND)-1.0_RKIND)*deltaZ + tableMin
          zeta1 = zeta0 + deltaZ

          if ((zeta <= zeta1) .and. (zeta >= zeta0)) then
            layerThick = ((zeta-zeta0)*subgridTable(lev+1) - (zeta-zeta1)*subgridTable(lev))/deltaZ
            return
          end if

        end do
      end if

   end subroutine ocn_subgrid_layer_thickness_lookup!}}}

!***********************************************************************
!
!  routine ocn_subgrid_wet_fraction_lookup
!
!> \brief   Wet fraction lookup
!> \author  Steven Brus, Damrongsak Wirasaet
!> \date    September 2022
!> \details Returns the wet fraction for a given ssh (defined over the
!>          wet fraction) value based on subgrid lookup table information
!
!-----------------------------------------------------------------------

   subroutine ocn_subgrid_wet_fraction_lookup(zeta, &
                                              subgridTable, &
                                              subgridTableRange, &
                                              wetFraction)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), intent(in) :: zeta
      real (kind=RKIND), dimension(:), intent(in) :: subgridTable
      real (kind=RKIND), dimension(:), intent(in) :: subgridTableRange

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), intent(out) :: wetFraction

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: lev

      real (kind=RKIND) :: tableMin
      real (kind=RKIND) :: tableMax
      real (kind=RKIND) :: deltaZ
      real (kind=RKIND) :: zeta0, zeta1

      tableMin = subgridTableRange(1)
      tableMax = subgridTableRange(2)
      deltaZ = (tableMax - tableMin)/real(nSubgridTableLevels-1,RKIND)

      if (zeta >= tableMax) then
        wetFraction = 1.0_RKIND
      else if (zeta <= tableMin) then
        wetFraction = 0.0_RKIND
      else
        do lev = 1, nSubgridTableLevels-1
          zeta0 = (real(lev,RKIND)-1.0_RKIND)*deltaZ + tableMin
          zeta1 = zeta0 + deltaZ

          if ((zeta <= zeta1) .and. (zeta >= zeta0)) then
            wetFraction = ((zeta-zeta0)*subgridTable(lev+1) - (zeta-zeta1)*subgridTable(lev))/deltaZ
            return
          end if

        end do
      end if

   end subroutine ocn_subgrid_wet_fraction_lookup!}}}

!***********************************************************************
!
!  routine ocn_subgrid_ssh_lookup
!> \brief   Inverse subgrid lookup (layerThickness -> ssh)
!> \author  Steven Brus, Damrongsak Wirasaet
!> \date    September 2022
!> \details Returns the ssh (defined over the wet fraction) for a given
!>          wet volume per unit area value (layerThick) based on subgrid
!>          lookup table information
!
!-----------------------------------------------------------------------

   subroutine ocn_subgrid_ssh_lookup(layerThick, &
                                     subgridTable, &
                                     subgridTableRange, &
                                     bathymetryMean, &
                                     bathymetryMin, &
                                     zeta)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), intent(in) :: layerThick
      real (kind=RKIND), dimension(:), intent(in) :: subgridTable
      real (kind=RKIND), dimension(:), intent(in) :: subgridTableRange
      real (kind=RKIND), intent(in) :: bathymetryMean
      real (kind=RKIND), intent(in) :: bathymetryMin

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), intent(inout) :: zeta

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: lev

      real (kind=RKIND) :: deltaZ
      real (kind=RKIND) :: tableMin, tableMax
      real (kind=RKIND) :: layerThickMin, layerThickMax
      real (kind=RKIND) :: layerThick0, layerThick1
      real (kind=RKIND) :: zeta0, zeta1
      real (kind=RKIND) :: phi0, phi1

      tableMin = subgridTableRange(1)
      tableMax = subgridTableRange(2)
      deltaZ = (tableMax - tableMin)/real(nSubgridTableLevels-1,RKIND)

      layerThickMin = subgridTable(1)
      layerThickMax = subgridTable(nSubgridTableLevels)

      if (layerThick >= layerThickMax) then
        zeta = layerThick - bathymetryMean

      else if (layerThick <= layerThickMin) then
        zeta = - bathymetryMin ! prevent_drying likely fails to ensure positive water columbn
      else

       do lev = 1, nSubgridTableLevels-1
          zeta0 = (real(lev,RKIND)-1.0_RKIND)*deltaZ + tableMin
          zeta1 = zeta0 + deltaZ

          layerThick0 = subgridTable(lev)
          layerThick1 = subgridTable(lev+1)

          if ((layerThick <= layerThick1) .and. (layerThick >= layerThick0)) then
            phi0 = (layerThick-layerThick1)/(layerThick0-layerThick1)
            phi1 = (layerThick-layerThick0)/(layerThick1-layerThick0)
            zeta = phi0*zeta0 + phi1*zeta1
            return
          end if

        end do
      end if

   end subroutine ocn_subgrid_ssh_lookup!}}}

!***********************************************************************
!
!  routine ocn_subgrid_init
!
!> \brief   Initializes subgrid wetting and drying module.
!> \author  Steven Brus, Damrongsak Wirasaet
!> \date    September 2022
!> \details
!>  This routine initializes the subgrid wetting and drying module
!
!-----------------------------------------------------------------------

   subroutine ocn_subgrid_init(domain,err)!{{{

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      type(block_type), pointer :: block
      type(mpas_pool_type), pointer :: meshPool

      err = 0
      if (.not. config_use_subgrid_wetting_drying) then
        return
      end if

      if ((config_ocean_run_mode == 'forward') .and. (nVertLevels .ne. 1)) then
        call mpas_log_write('config_config_use_subgrid_wetting_drying = .true. requires single layer' , MPAS_LOG_CRIT)
      end if

      block => domain%blocklist
      call mpas_pool_get_subpool(block%structs, 'mesh', meshPool)

      call mpas_pool_get_array(meshPool, 'subgridWetVolumeCellTable', &
                               subgridWetVolumeCellTable)
      call mpas_pool_get_array(meshPool, 'subgridWetVolumeEdgeTable', &
                               subgridWetVolumeEdgeTable)
      call mpas_pool_get_array(meshPool, 'subgridWetVolumeVertexTable', &
                               subgridWetVolumeVertexTable)
      call mpas_pool_get_array(meshPool, 'subgridWetFractionCellTable', &
                               subgridWetFractionCellTable)
      call mpas_pool_get_array(meshPool, 'subgridWetFractionEdgeTable', &
                               subgridWetFractionEdgeTable)
      call mpas_pool_get_array(meshPool, 'subgridWetFractionVertexTable', &
                               subgridWetFractionVertexTable)
      call mpas_pool_get_array(meshPool, 'subgridSshCellTableRange', &
                               subgridSshCellTableRange)
      call mpas_pool_get_array(meshPool, 'subgridSshEdgeTableRange', &
                               subgridSshEdgeTableRange)
      call mpas_pool_get_array(meshPool, 'subgridSshVertexTableRange', &
                               subgridSshVertexTableRange)
      call mpas_pool_get_array(meshPool, 'subgridEdgeBathymetryMean', &
                               subgridEdgeBathymetryMean)
      call mpas_pool_get_array(meshPool, 'subgridVertexBathymetryMean', &
                               subgridVertexBathymetryMean)
      call mpas_pool_get_array(meshPool, 'subgridCellBathymetryMin', &
                               subgridCellBathymetryMin)
      call mpas_pool_get_array(meshPool, 'subgridEdgeBathymetryMin', &
                               subgridEdgeBathymetryMin)
      call mpas_pool_get_array(meshPool, 'subgridVertexBathymetryMin', &
                               subgridVertexBathymetryMin)
      call mpas_pool_get_array(meshPool, 'subgridLayerThicknessDebug', &
                               subgridLayerThicknessDebug)

      call mpas_pool_get_dimension(meshPool, 'nSubgridTableLevels', &
                                     nSubgridTableLevels)

   end subroutine ocn_subgrid_init!}}}

end module ocn_subgrid

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
