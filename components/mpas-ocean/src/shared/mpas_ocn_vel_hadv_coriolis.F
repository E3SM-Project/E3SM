! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_hadv_coriolis
!
!> \brief MPAS ocean horizontal momentum mixing driver
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011
!> \details
!>  This module contains the routine for computing
!>  tendencies from the horizontal advection and coriolis force.
!>
!
!-----------------------------------------------------------------------

module ocn_vel_hadv_coriolis

   use mpas_timer
   use ocn_constants
   use ocn_mesh
   use ocn_config

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

   public :: ocn_vel_hadv_coriolis_tend, &
             ocn_vel_hadv_coriolis_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: &
      hadvCoriolisDisabled ! on/off switch for hadv/Coriolis tend
                           ! use disabled since default is on

   logical :: &
      usePlanetVorticity   ! multiplicative mask for including
                           ! planetary vorticity term

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vel_hadv_coriolis_tend
!
!> \brief   Computes tendency term for horizontal advection and coriolis force
!> \author  Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date    September 2011
!> \details This routine computes the horizontal advection and
!>          coriolis advection tendencies for momentum based on
!>          current state.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_hadv_coriolis_tend(normRelVortEdge, &
                                      normPlanetVortEdge, &
                                      layerThickEdgeFlux, normalVelocity, &
                                      kineticEnergyCell, tend, err)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normRelVortEdge,    &!< [in] relative vorticity/thickness at edge
         normPlanetVortEdge, &!< [in] planetary vorticity/thickness at edge
         layerThickEdgeFlux,     &!< [in] layer thickness on edge
         normalVelocity,     &!< [in] Horizontal velocity
         kineticEnergyCell    !< [in] Kinetic Energy

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend                 !< [inout] accumulated velocity tendency

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< [out] error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer ::       &
         iEdge, j, k,  &! loop indices for edges, vertical
         cell1, cell2, &! neighbor cell indices
         eoe,          &! edge on edge index
         kmin, kmax     ! min and max indices for vertical

      real (kind=RKIND) :: &
         avgVorticity, &! vorticity averaged across edge
         invLength,    &! temp variable for 1/dcedge
         edgeWeight     ! weight for summing over edges

      real (kind=RKIND), dimension(:,:), allocatable :: &
         qArr,         &! temporary vorticity array
         tmpVorticity   ! vorticity on edge

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      !*** Set error code and exit early if disabled
      !*** Otherwise, start timer

      err = 0
      if (hadvCoriolisDisabled) return
      call mpas_timer_start("coriolis")

      !*** allocate and transfer temporary arrays

      allocate(tmpVorticity(nVertLevels,nEdgesAll), &
                       qArr(nVertLevels,nEdgesAll))
      !$acc enter data create(tmpVorticity, qArr)

#ifndef MPAS_OPENACC
      !$omp parallel
#endif
      if (usePlanetVorticity) then
#ifdef MPAS_OPENACC
         !$acc parallel loop &
         !$acc    present(tmpVorticity, normRelVortEdge, &
         !$acc            normPlanetVortEdge) &
         !$acc    private(k)
#else
         !$omp do schedule(runtime) &
         !$omp    private(k)
#endif
         do iEdge = 1, nEdgesAll
         do k=1,nVertLevels
            tmpVorticity(k,iEdge) = normRelVortEdge(k,iEdge) + &
                                 normPlanetVortEdge(k,iEdge)
         end do
         end do
#ifndef MPAS_OPENACC
         !$omp end do
#endif
      else
#ifdef MPAS_OPENACC
         !$acc parallel loop &
         !$acc    present(tmpVorticity, normRelVortEdge) &
         !$acc    private(k)
#else
         !$omp do schedule(runtime) &
         !$omp    private(k)
#endif
         do iEdge = 1, nEdgesAll
         do k=1,nVertLevels
            tmpVorticity(k,iEdge) = normRelVortEdge(k,iEdge)
         end do
         end do
#ifndef MPAS_OPENACC
         !$omp end do
#endif
      endif
#ifndef MPAS_OPENACC
         !$omp end parallel
#endif

#ifdef MPAS_OPENACC
      !$acc parallel loop &
      !$acc    present(nEdgesOnEdge, edgesOnEdge, weightsOnEdge, &
      !$acc            minLevelEdgeBot, maxLevelEdgeTop, qArr, tmpVorticity, &
      !$acc            normalVelocity, layerThickEdgeFlux) &
      !$acc    private(k, j, eoe, edgeWeight, avgVorticity, kmin, kmax)
#else
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp    private(k, j, eoe, edgeWeight, avgVorticity, kmin, kmax)
#endif
      do iEdge = 1, nEdgesOwned
         kmin = minLevelEdgeBot(iEdge)
         kmax = maxLevelEdgeTop(iEdge)

         do k = 1, nVertLevels
            qArr(k,iEdge) = 0.0_RKIND
         end do

         do j = 1, nEdgesOnEdge(iEdge)
            eoe = edgesOnEdge(j, iEdge)
            edgeWeight = weightsOnEdge(j, iEdge)

            do k = kmin, kmax
               avgVorticity = 0.5_RKIND* &
                              (tmpVorticity(k,iEdge) + &
                               tmpVorticity(k,eoe))
               qArr(k,iEdge) = qArr(k,iEdge) + &
                               edgeWeight*normalVelocity(k,eoe)* &
                               avgVorticity*layerThickEdgeFlux(k,eoe)
            end do

         end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
#endif

#ifdef MPAS_OPENACC
      !$acc parallel loop &
      !$acc    present(cellsOnEdge, minLevelEdgeBot, maxLevelEdgeTop, edgeMask, &
      !$acc            dcEdge, tend, qArr, kineticEnergyCell) &
      !$acc    private(cell1, cell2, invLength, k, kmin, kmax)
#else
      !$omp do schedule(runtime) &
      !$omp    private(cell1, cell2, invLength, k, kmin, kmax)
#endif
      do iEdge = 1, nEdgesOwned
         kmin = minLevelEdgeBot(iEdge)
         kmax = maxLevelEdgeTop(iEdge)

         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         invLength = 1.0_RKIND / dcEdge(iEdge)

         do k = kmin, kmax
            tend(k,iEdge) = tend(k,iEdge) + &
                        edgeMask(k,iEdge)* (qArr(k,iEdge) - &
                           (kineticEnergyCell(k,cell2) &
                          - kineticEnergyCell(k,cell1))*invLength)
         end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      !$acc exit data delete(tmpVorticity, qArr)
      deallocate(qArr,tmpVorticity)

      call mpas_timer_stop("coriolis")

   !--------------------------------------------------------------------

   end subroutine ocn_vel_hadv_coriolis_tend!}}}

!***********************************************************************
!
!  routine ocn_vel_hadv_coriolis_init
!
!> \brief   Initializes ocean momentum horizontal advection and
!>          coriolis tendencies
!> \author  Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date    September 2011
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_hadv_coriolis_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      ! Output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< [out] error flag

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      !*** Initialize returned error code and default flags

      err = 0

      hadvCoriolisDisabled = .false.
      usePlanetVorticity   = .false.

      !*** Reset module variables based in input configuration

      if ( config_disable_vel_coriolis ) hadvCoriolisDisabled = .true.

      select case (trim(config_time_integrator))
      case ('RK4','rk4')
         ! For RK4, coriolis tendency term includes f: (eta+f)/h.
         usePlanetVorticity = .true.

      case ('split_explicit')
         ! For split explicit, Coriolis tendency uses eta/h because
         ! the Coriolis term is added separately to the momentum
         ! tendencies.
         usePlanetVorticity = .false.

      case ('unsplit_explicit')
         ! For unsplit explicit, Coriolis tendency uses eta/h because
         ! the Coriolis term is added separately to the momentum
         ! tendencies.
         usePlanetVorticity = .false.

      case ('split_implicit')
         ! For split implicit, Coriolis tendency uses eta/h because
         ! the Coriolis term is added separately to the momentum
         ! tendencies.
         usePlanetVorticity = .false.

      end select

   !--------------------------------------------------------------------

   end subroutine ocn_vel_hadv_coriolis_init!}}}

!***********************************************************************

end module ocn_vel_hadv_coriolis

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
