! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module ocn_gm

   use mpas_grid_types
   use mpas_timer
   use ocn_constants
   
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

   public :: ocn_gm_compute_uBolus

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------
   character (len=StrKIND), pointer :: config_vert_coord_movement
   real (kind=RKIND), pointer :: config_h_kappa, config_h_kappa_q

contains

   subroutine ocn_gm_compute_uBolus(statePool, diagnosticsPool, meshPool)!{{{
      implicit none
      type (mpas_pool_type), intent(in)  :: statePool !< Input/Output: State information
      type (mpas_pool_type), intent(inout)  :: diagnosticsPool !< Input/Output: Diagnostics information
      type (mpas_pool_type), intent(in)     :: meshPool !< Input: Mesh information

      real(kind=RKIND), dimension(:,:), pointer :: uBolusGM, hEddyFlux, layerThicknessEdge

      integer, dimension(:), pointer   :: maxLevelEdgeTop
      integer                          :: k, iEdge
      integer, pointer                 :: nEdges

      call mpas_pool_get_array(diagnosticsPool, 'uBolusGM', uBolusGM)
      call mpas_pool_get_array(diagnosticsPool, 'hEddyFlux', hEddyFlux)
      call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)

      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)

      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)

      call ocn_gm_compute_hEddyFlux(statePool, diagnosticsPool, meshPool)

      if (config_vert_coord_movement .EQ. 'impermeable_interfaces') then

         do iEdge = 1, nEdges
            do k = 1, maxLevelEdgeTop(iEdge)
               uBolusGM(k,iEdge) = hEddyFlux(k,iEdge) / layerThicknessEdge(k,iEdge)
            end do
         end do

      else

         ! Nothing for now for all other mesh types (zlevel, zstar, ztilde)
         uBolusGM(:,:) = 0.0

      end if

   end subroutine ocn_gm_compute_uBolus!}}}

   subroutine ocn_gm_compute_hEddyFlux(statePool, diagnosticsPool, meshPool, timeLevelIn)!{{{
      implicit none
      type (mpas_pool_type), intent(in) :: statePool !< Input/Output: State information
      type (mpas_pool_type), intent(inout)    :: diagnosticsPool !< Input: Diagnostics information
      type (mpas_pool_type), intent(in)    :: meshPool !< Input: Mesh information
      integer, intent(in), optional :: timeLevelIn

      real(kind=RKIND), dimension(:,:), pointer  :: hEddyFlux, layerThickness
      real(kind=RKIND), dimension(:), pointer    :: dcEdge
      integer, dimension(:,:), pointer           :: cellsOnEdge
      integer, dimension(:), pointer             :: maxLevelEdgeTop
      integer                                    :: k, cell1, cell2, iEdge
      integer, pointer                           :: nEdges
      integer :: timeLevel

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      call mpas_pool_get_array(diagnosticsPool, 'hEddyFlux', hEddyFlux)

      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)

      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)

      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)

      hEddyFlux(:,:) = 0.0

      if (config_vert_coord_movement .EQ. 'impermeable_interfaces') then
            do iEdge = 1, nEdges
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               do k = 1, maxLevelEdgeTop(iEdge)
                  hEddyFlux(k,iEdge) = -config_h_kappa * (layerThickness(k,cell2) - layerThickness(k,cell1)) / dcEdge(iEdge)
               end do
            end do
      else

         !Nothing for now for all other mesh types (zlevel, zstar, ztilde)

      end if
                  
   end subroutine ocn_gm_compute_hEddyFlux!}}}

   subroutine ocn_get_h_kappa(diagnosticsPool, meshPool)!{{{

      type (mpas_pool_type), intent(inout) :: diagnosticsPool !< Input/Output: Diagnostics information
      type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information

      real(kind=RKIND), dimension(:,:), pointer    :: hKappa


      call mpas_pool_get_array(diagnosticsPool, 'hKappa', hKappa)

      hKappa(:,:) = config_h_kappa

   end subroutine ocn_get_h_kappa!}}}

   subroutine ocn_get_h_kappa_q(diagnosticsPool, meshPool)!{{{

      type (mpas_pool_type), intent(inout) :: diagnosticsPool !< Input/Output: Diagnostics information
      type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information

      real(kind=RKIND), dimension(:,:), pointer    :: hKappaQ


      call mpas_pool_get_array(diagnosticsPool, 'hKappaQ', hKappaQ)

      hKappaQ(:,:) = config_h_kappa_q


   end subroutine ocn_get_h_kappa_q!}}}

   subroutine ocn_gm_init(err)!{{{
      integer, intent(out) :: err

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_vert_coord_movement', config_vert_coord_movement)
      call mpas_pool_get_config(ocnConfigs, 'config_h_kappa', config_h_kappa)
      call mpas_pool_get_config(ocnConfigs, 'config_h_kappa_q', config_h_kappa_q)

   end subroutine ocn_gm_init!}}}

end module ocn_gm
