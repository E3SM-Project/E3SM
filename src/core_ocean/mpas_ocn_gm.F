! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module ocn_gm

   use mpas_grid_types
   use mpas_configure
   use mpas_timer
   
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

contains

   subroutine ocn_gm_compute_uBolus(state, diagnostics, mesh)!{{{
      implicit none
      type(state_type), intent(inout)        :: state !< Input/Output: State information
      type(diagnostics_type), intent(inout)  :: diagnostics !< Input/Output: Diagnostics information
      type(mesh_type), intent(in)            :: mesh !< Input: Mesh information

      real(kind=RKIND), dimension(:,:), pointer :: uBolusGM, hEddyFlux, layerThicknessEdge

      integer, dimension(:), pointer   :: maxLevelEdgeTop
      integer                          :: k, iEdge, nEdges

      uBolusGM           => diagnostics % uBolusGM % array
      hEddyFlux          => diagnostics % hEddyFlux % array
      layerThicknessEdge => diagnostics % layerThicknessEdge % array

      maxLevelEdgeTop => mesh % maxLevelEdgeTop % array

      nEdges = mesh % nEdges

      call ocn_gm_compute_hEddyFlux(state, diagnostics, mesh)

      if (config_vert_coord_movement .EQ. 'impermeable_interfaces') then

         do iEdge = 1, nEdges
            do k = 1, maxLevelEdgeTop(iEdge)
               uBolusGM(k,iEdge) = hEddyFlux(k,iEdge)/layerThicknessEdge(k,iEdge)
            end do
         end do

      else

         ! Nothing for now for all other mesh types (zlevel, zstar, ztilde)
         uBolusGM(:,:) = 0.0

      end if

   end subroutine ocn_gm_compute_uBolus!}}}

   subroutine ocn_gm_compute_hEddyFlux(state, diagnostics, mesh)!{{{
      implicit none
      type(state_type), intent(inout)     :: state !< Input/Output: State information
      type(diagnostics_type), intent(in)  :: diagnostics !< Input: Diagnostics information
      type(mesh_type), intent(in)         :: mesh !< Input: Mesh information

      real(kind=RKIND), dimension(:,:), pointer  :: hEddyFlux, layerThickness
      real(kind=RKIND), dimension(:), pointer    :: dcEdge
      integer, dimension(:,:), pointer           :: cellsOnEdge
      integer, dimension(:), pointer             :: maxLevelEdgeTop
      integer                                    :: k, cell1, cell2, iEdge, nEdges

      hEddyFlux      => diagnostics % hEddyFlux % array
      layerThickness => state % layerThickness % array

      dcEdge         => mesh % dcEdge % array
      cellsOnEdge    => mesh % cellsOnEdge % array
      maxLevelEdgeTop => mesh % maxLevelEdgeTop % array

      nEdges = mesh % nEdges

      hEddyFlux(:,:) = 0.0

      if (config_vert_coord_movement .EQ. 'impermeable_interfaces') then
            do iEdge = 1,nEdges
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               do k=1,maxLevelEdgeTop(iEdge)
                  hEddyFlux(k,iEdge) = -config_h_kappa * (layerThickness(k,cell2) - layerThickness(k,cell1)) / dcEdge(iEdge)
               end do
            end do
      else

         !Nothing for now for all other mesh types (zlevel, zstar, ztilde)

      end if
                  
   end subroutine ocn_gm_compute_hEddyFlux!}}}

   subroutine ocn_get_h_kappa(diagnostics, mesh)!{{{

      type (diagnostics_type), intent(inout) :: diagnostics !< Input/Output: Diagnostics information
      type (mesh_type), intent(in) :: mesh !< Input: Mesh information

      real(kind=RKIND), dimension(:,:), pointer    :: hKappa


      hKappa  => diagnostics % hKappa % array

      hKappa(:,:) = config_h_kappa


   end subroutine ocn_get_h_kappa!}}}

   subroutine ocn_get_h_kappa_q(diagnostics, mesh)!{{{

      type (diagnostics_type), intent(inout) :: diagnostics !< Input/Output: Diagnostics information
      type (mesh_type), intent(in) :: mesh !< Input: Mesh information

      real(kind=RKIND), dimension(:,:), pointer    :: hKappaQ


      hKappaQ  => diagnostics % hKappaQ % array

      hKappaQ(:,:) = config_h_kappa_q


   end subroutine ocn_get_h_kappa_q!}}}

end module ocn_gm
