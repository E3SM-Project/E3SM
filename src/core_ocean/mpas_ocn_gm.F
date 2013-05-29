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

   subroutine ocn_gm_compute_uBolus(s, grid)!{{{
      implicit none
      type(state_type), intent(inout)        :: s
      type(mesh_type), intent(in)            :: grid

      real(kind=RKIND), dimension(:,:), pointer :: uBolusGM, hEddyFlux, layerThicknessEdge

      integer, dimension(:), pointer   :: maxLevelEdgeTop
      integer                          :: k, iEdge, nEdges

      uBolusGM         => s % uBolusGM % array
      layerThicknessEdge => s % layerThicknessEdge % array
      hEddyFlux      => s % hEddyFlux % array

      maxLevelEdgeTop => grid % maxLevelEdgeTop % array

      nEdges = grid % nEdges

      call ocn_gm_compute_hEddyFlux(s, grid)

      if (config_vert_coord_movement .EQ. 'isopycnal') then

         do iEdge = 1, nEdges
            do k = 1, maxLevelEdgeTop(iEdge)
               uBolusGM(k,iEdge) = hEddyFlux(k,iEdge)/layerThicknessEdge(k,iEdge)
            end do
         end do

      else

         ! Nothing for now for all other grid types (zlevel, zstar, ztilde)
         uBolusGM(:,:) = 0.0

      end if

   end subroutine ocn_gm_compute_uBolus!}}}

   subroutine ocn_gm_compute_hEddyFlux(s, grid)!{{{
      implicit none
      type(state_type), intent(inout)     :: s
      type(mesh_type), intent(in)         :: grid

      real(kind=RKIND), dimension(:,:), pointer  :: hEddyFlux, layerThickness
      real(kind=RKIND), dimension(:), pointer    :: dcEdge
      integer, dimension(:,:), pointer           :: cellsOnEdge
      integer, dimension(:), pointer             :: maxLevelEdgeTop
      integer                                    :: k, cell1, cell2, iEdge, nEdges

      hEddyFlux      => s % hEddyFlux % array
      layerThickness => s % layerThickness % array

      dcEdge         => grid % dcEdge % array
      cellsOnEdge    => grid % cellsOnEdge % array
      maxLevelEdgeTop => grid % maxLevelEdgeTop % array

      nEdges = grid % nEdges

      hEddyFlux(:,:) = 0.0

      if (config_vert_coord_movement .EQ. 'isopycnal') then
            do iEdge = 1,nEdges
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               do k=1,maxLevelEdgeTop(iEdge)
                  hEddyFlux(k,iEdge) = -config_h_kappa * (layerThickness(k,cell2) - layerThickness(k,cell1)) / dcEdge(iEdge)
               end do
            end do
      else

         !Nothing for now for all other grid types (zlevel, zstar, ztilde)

      end if
                  
   end subroutine ocn_gm_compute_hEddyFlux!}}}

   subroutine ocn_get_h_kappa(s, grid)!{{{

      type (state_type), intent(inout) :: s
      type (mesh_type), intent(in) :: grid

      real(kind=RKIND), dimension(:,:), pointer    :: h_kappa


      h_kappa  => s % h_kappa % array

      h_kappa(:,:) = config_h_kappa


   end subroutine ocn_get_h_kappa!}}}

   subroutine ocn_get_h_kappa_q(s, grid)!{{{

      type (state_type), intent(inout) :: s
      type (mesh_type), intent(in) :: grid

      real(kind=RKIND), dimension(:,:), pointer    :: h_kappa_q


      h_kappa_q  => s % h_kappa_q % array

      h_kappa_q(:,:) = config_h_kappa_q


   end subroutine ocn_get_h_kappa_q!}}}

end module ocn_gm
