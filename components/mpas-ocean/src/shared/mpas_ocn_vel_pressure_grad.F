! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_pressure_grad
!
!> \brief MPAS ocean pressure gradient module
!> \author Mark Petersen
!> \date   September 2011
!> \details
!>  This module contains the routine for computing
!>  tendencies from the horizontal pressure gradient.
!>
!
!-----------------------------------------------------------------------

module ocn_vel_pressure_grad

   use mpas_timer
   use mpas_constants
   use mpas_log

   use ocn_constants
   use ocn_config
   use ocn_mesh

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   integer, public :: &
      pGradType   ! id for pressure gradient method selected

   integer, parameter, public :: &! ids for supported methods
      pGradTypeNone        = 0, &! none selected
      pGradTypeSSHgrad     = 1, &! ssh gradient
      pGradTypePZmid       = 2, &! pressure and zMid
      pGradTypeMontPot     = 3, &! Montgomery potential
      pGradTypeMontPotDens = 4, &! Montgomery potential and density
      pGradTypeJacobDens   = 5, &! Jacobian from density
      pGradTypeJacobTS     = 6, &! Jacobian from T,S
      pGradTypeConstForced = 7   ! constant forced

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_vel_pressure_grad_tend, &
             ocn_vel_pressure_grad_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: &
      pGradOff,         &! flag for turning pressure gradient on/off
      timeIntegratorLTS  ! flag for Local Time Stepping

   real (kind=RKIND) :: &! precomputed constants for efficiency
      density0Inv,      &! 1/density0
      gdensity0Inv,     &! g/density0
      inv12,            &! 1/12
      pGradLvlWgt        ! weighting for levels in Jacobian formulations

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vel_pressure_grad_tend
!
!> \brief   Computes tendency term for horizontal pressure gradient
!> \author  Mark Petersen
!> \date    February 2014
!> \details
!>  This routine computes the pressure gradient tendency for momentum
!>  based on current state.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_pressure_grad_tend(ssh, pressure, surfacePressure, &
                                  montgomeryPotential, zMid, &
                                  density, potentialDensity, &
                                  indxT, indxS, tracers, &
                                  thermExpCoeff, salineContractCoeff, &
                                  tend, err)!{{{

      implicit none

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      integer, intent(in) :: &
         indxT,              &!< [in] tracer array index for temperature
         indxS                !< [in] tracer array index for salt

      real (kind=RKIND), dimension(:), intent(in) :: &
         ssh,                &!< [in] sea surface height
         surfacePressure      !< [in] surface pressure

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         pressure,            &!< [in] Pressure field
         montgomeryPotential, &!< [in] Mongomery potential
         zMid,                &!< [in] z-coordinate at layer mid-depth
         density,             &!< [in] density
         potentialDensity      !< [in] potentialDensity

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         thermExpCoeff,       &!< [in] in situ thermal expansion coeff
         salineContractCoeff   !< [in] in situ saline contraction coeff

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         tracers               !< [in] array of active tracers

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend          !< [inout] accumulated velocity tendency

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< [out] error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer ::           &
         iEdge, k,  &! loop indices for edge, cell, vertical loops
         cell1, cell2,     &! neighbor cell indices across edge
         kMin, kMax         ! shallowest and deepest active layer

      real (kind=RKIND) :: &
         invdcEdge,        &! temporary for 1/dcEdge
         pGrad,            &! vertically integrated pgrad
         alpha, beta,      &! T,S expansion factors for Jacobian
         Area,             &! Area for Jacobian calculation
         zStar, zC, zGamma,&! z coordinate combination for Jacobian
         rhoL, rhoR,       &! density neighbors for Jacobian
         TL, TR,           &! temperature neighbors for Jacobian
         SL, SR             ! salinity neighbors for Jacobian

      real (kind=RKIND), dimension(:,:), allocatable :: &
         JacobianDxDs,     &! Jacobian in density * Dx * DS
         JacobianTz,       &! Jacobian associated with temperature
         JacobianSz         ! Jacobian associated with salinity

      real (kind=RKIND) :: &
         meridSshGrad,     &! config_meridional_ssh_grad
         zonalSshGrad       ! config_zonal_ssh_grad

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      !*** set error code and exit if turned off
      !*** start timer if turned on

      err = 0
      if (pGradOff) return
      call mpas_timer_start("pressure grad")

      meridSshGrad = config_meridional_ssh_grad
      zonalSshGrad = config_zonal_ssh_grad

      !*** Compute pGrad based on method selected
      select case (pGradType)

      case (pGradTypeSSHgrad)

         ! pressure for sea surface height = - g grad ssh

         if (timeIntegratorLTS) then
#ifdef MPAS_OPENACC
            !$acc parallel loop &
            !$acc    present(cellsOnEdge, minLevelEdgeBot, maxLevelEdgeTop, dcEdge, &
            !$acc            tend, edgeMask, ssh,surfacePressure) &
            !$acc    private(cell1, cell2, invdcEdge, k, kMin, kMax)
#else
            !$omp parallel
            !$omp do schedule(runtime) &
            !$omp    private(cell1, cell2, invdcEdge, k, kMin, kMax)
#endif
            do iEdge=1,nEdgesOwned
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               invdcEdge = 1.0_RKIND / dcEdge(iEdge)
               kMin = minLevelEdgeBot(iEdge)
               kMax = maxLevelEdgeTop(iEdge)

               do k=kMin,kMax
                  tend(k,iEdge) = tend(k,iEdge) - &
                                  edgeMask(k,iEdge)*invdcEdge* &
                                  (density0Inv*(  surfacePressure(cell2) - surfacePressure(cell1)) )
               end do
            end do
#ifndef MPAS_OPENACC
            !$omp end do
            !$omp end parallel
#endif
         else ! not LTS

#ifdef MPAS_OPENACC
            !$acc parallel loop &
            !$acc    present(cellsOnEdge, minLevelEdgeBot, maxLevelEdgeTop, dcEdge, &
            !$acc            tend, edgeMask, ssh,surfacePressure) &
            !$acc    private(cell1, cell2, invdcEdge, k, kMin, kMax)
#else
            !$omp parallel
            !$omp do schedule(runtime) &
            !$omp    private(cell1, cell2, invdcEdge, k, kMin, kMax)
#endif
            do iEdge=1,nEdgesOwned
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               invdcEdge = 1.0_RKIND / dcEdge(iEdge)
               kMin = minLevelEdgeBot(iEdge)
               kMax = maxLevelEdgeTop(iEdge)

               do k=kMin,kMax
                  tend(k,iEdge) = tend(k,iEdge) - &
                                  edgeMask(k,iEdge)*invdcEdge* &
                                  (  gravity*(ssh(cell2) - ssh(cell1)) &
                                   + density0Inv*(  surfacePressure(cell2) &
                                                  - surfacePressure(cell1)) )
               end do
            end do
#ifndef MPAS_OPENACC
            !$omp end do
            !$omp end parallel
#endif
         endif ! LTS

      case (pGradTypePZmid)

         ! pressure for generalized coordinates
         ! -1/density_0 (grad p_k + density g grad z_k^{mid})

#ifdef MPAS_OPENACC
         !$acc parallel loop &
         !$acc    present(cellsOnEdge, minLevelEdgeBot, maxLevelEdgeTop, dcEdge, zMid, &
         !$acc            tend, edgeMask, pressure, density) &
         !$acc    private(cell1, cell2, invdcEdge, k, kMin, kMax)
#else
         !$omp parallel
         !$omp do schedule(runtime) &
         !$omp    private(cell1, cell2, invdcEdge, k, kMin, kMax)
#endif
         do iEdge=1,nEdgesOwned
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            invdcEdge = 1.0_RKIND / dcEdge(iEdge)
            kMin = minLevelEdgeBot(iEdge)
            kMax = maxLevelEdgeTop(iEdge)

            do k=kMin,kMax
               tend(k,iEdge) = tend(k,iEdge) + &
                               edgeMask(k,iEdge)*invdcEdge*( &
                               - density0Inv*(pressure(k,cell2) - &
                                              pressure(k,cell1)) &
                               - gdensity0Inv*0.5_RKIND* &
                                 (density(k,cell1)+density(k,cell2))* &
                                 (zMid(k,cell2)-zMid(k,cell1)))
            end do
         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif

      case (pGradTypeMontPot)

         ! For pure isopycnal coordinates, this is just grad(M),
         ! the gradient of Montgomery Potential

#ifdef MPAS_OPENACC
         !$acc parallel loop &
         !$acc    present(cellsOnEdge, minLevelEdgeBot, maxLevelEdgeTop, dcEdge, &
         !$acc            tend, edgeMask, montgomeryPotential) &
         !$acc    private(cell1, cell2, invdcEdge, k, kMin, kMax)
#else
         !$omp parallel
         !$omp do schedule(runtime) &
         !$omp    private(cell1, cell2, invdcEdge, k, kMin, kMax)
#endif
         do iEdge=1,nEdgesOwned
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            invdcEdge = 1.0_RKIND / dcEdge(iEdge)
            kMin = minLevelEdgeBot(iEdge)
            kMax = maxLevelEdgeTop(iEdge)

            do k=kMin,kMax
               tend(k,iEdge) = tend(k,iEdge) + &
                               edgeMask(k,iEdge)*invdcEdge* ( &
                               - (montgomeryPotential(k,cell2) - &
                                  montgomeryPotential(k,cell1)))
            end do
         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif

      case (pGradTypeMontPotDens)

         ! This formulation has not been extensively tested and is not
         ! supported at this time.

         ! This is -grad(M)+p grad(1/rho)
         ! Where rho is the potential density.
         ! See Bleck (2002) equation 1, and last equation in Appendix A.

#ifdef MPAS_OPENACC
         !$acc parallel loop &
         !$acc    present(cellsOnEdge, dcEdge, minLevelEdgeBot, maxLevelEdgeTop, &
         !$acc            edgeMask, tend, montgomeryPotential, &
         !$acc            pressure, potentialDensity) &
         !$acc    private(cell1, cell2, invdcEdge, k, kMin, kMax)
#else
         !$omp parallel
         !$omp do schedule(runtime) &
         !$omp    private(cell1, cell2, invdcEdge, k, kMin, kMax)
#endif
         do iEdge=1,nEdgesOwned
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            invdcEdge = 1.0_RKIND / dcEdge(iEdge)
            kMin = minLevelEdgeBot(iEdge)
            kMax = maxLevelEdgeTop(iEdge)

            do k=kMin,kMax
               tend(k,iEdge) = tend(k,iEdge) + &
                               edgeMask(k,iEdge)*invdcEdge*( &
                               - (montgomeryPotential(k,cell2) - &
                                  montgomeryPotential(k,cell1)) &
                               +  0.5_RKIND*(pressure(k,cell1) + &
                                             pressure(k,cell2))* &
                              ( 1.0_RKIND/potentialDensity(k,cell2) &
                              - 1.0_RKIND/potentialDensity(k,cell1)))
            end do
         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif

      case (pGradTypeJacobDens)

         allocate(JacobianDxDs(nVertLevels,nEdgesAll))
         !$acc enter data create (JacobianDxDs)

#ifdef MPAS_OPENACC
         !$acc parallel loop &
         !$acc    present(cellsOnEdge, minLevelEdgeBot, maxLevelEdgeTop, dcEdge, &
         !$acc            edgeMask, pressure, density, zMid, tend, &
         !$acc            JacobianDxDs) &
         !$acc    private(cell1, cell2, invdcEdge, k, kMin, kMax, pGrad, &
         !$acc            Area, zStar, zC, zGamma, rhoL, rhoR)
#else
         !$omp parallel
         !$omp do schedule(runtime) &
         !$omp    private(cell1, cell2, invdcEdge, k, kMin, kMax, pGrad, &
         !$omp            Area, zStar, zC, zGamma, rhoL, rhoR)
#endif
         do iEdge=1,nEdgesOwned
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            invdcEdge = 1.0_RKIND / dcEdge(iEdge)
            kMin = minLevelEdgeBot(iEdge)
            kMax = maxLevelEdgeTop(iEdge)

            ! Compute the density-Jacobian in common_level form.
            ! See Shchepetkin and McWilliams (2003) Ocean Modeling,
            !   sections 2-4

            JacobianDxDs(kMin,iEdge) = 0.0_RKIND

            do k=kMin+1,kMax

               ! eqn 2.7 in Shchepetkin and McWilliams (2003)
               ! Note delta x was removed.  It must be an error in the
               ! paper, ! as it makes the units incorrect.
               Area = 0.5_RKIND*(zMid(k-1,cell1) - zMid(k,cell1) + &
                                 zMid(k-1,cell2) - zMid(k,cell2))

               ! eqn 2.8
               zStar = ( zMid(k-1,cell2)*zMid(k-1,cell1) - &
                         zMid(k  ,cell2)*zMid(k  ,cell1) )/ &
                       ( zMid(k-1,cell2)-zMid(k  ,cell2) + &
                         zMid(k-1,cell1)-zMid(k,cell1))

               ! eqn 3.2
               zC = 0.25_RKIND*(zMid(k,cell1) + zMid(k-1,cell1) + &
                                zMid(k,cell2) + zMid(k-1,cell2))

               ! eqn 4.1
               zGamma = (1.0_RKIND - pGradLvlWgt)*zStar + pGradLvlWgt*zC

               rhoL = (density(k  ,cell1)*(zMid(k-1,cell1)-zGamma) + &
                       density(k-1,cell1)*(zGamma-zMid(k  ,cell1)))/ &
                      (zMid(k-1,cell1) - zMid(k,cell1))
               rhoR = (density(k  ,cell2)*(zMid(k-1,cell2)-zGamma) +  &
                       density(k-1,cell2)*(zGamma-zMid(k  ,cell2)))/ &
                      (zMid(k-1,cell2) - zMid(k,cell2))

               ! eqn 2.6 in Shchepetkin and McWilliams (2003)
               JacobianDxDs(k,iEdge) = Area * (rhoL - rhoR)
            end do

            ! In the top layer, use pressure for generalized coordinates
            ! pGrad = -1/density_0 (grad p_k + density g grad z_k^{mid})
            k = kMin
            pGrad = edgeMask(k,iEdge)*invdcEdge*( &
                    - density0Inv*(pressure(k,cell2) - &
                                   pressure(k,cell1)) &
                    - gdensity0Inv*0.5_RKIND* &
                      (density(k,cell1)+density(k,cell2))* &
                      (zMid(k,cell2)- zMid(k,cell1) ) )

            tend(k,iEdge) = tend(k,iEdge) + pGrad

            do k=kMin+1,kMax

               ! note JacobianDxDs includes negative sign, so
               ! pGrad is - g/rho_0 dP/dx

               pGrad = pGrad + gdensity0Inv*JacobianDxDs(k,iEdge)*invdcEdge

               tend(k,iEdge) = tend(k,iEdge) + pGrad

            end do
         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif

         !$acc exit data delete(JacobianDxDs)
         deallocate(JacobianDxDs)

     case (pGradTypeJacobTS)

         allocate(JacobianDxDs(nVertLevels,nEdgesAll), &
                  JacobianTz(nVertLevels,nEdgesAll), &
                  JacobianSz(nVertLevels,nEdgesAll))
         !$acc enter data create(JacobianDxDs, &
         !$acc                   JacobianTz, JacobianSz)

#ifdef MPAS_OPENACC
         !$acc parallel loop &
         !$acc    present(cellsOnEdge, minLevelEdgeBot, maxLevelEdgeTop, dcEdge, &
         !$acc            zMid, tracers, edgeMask, pressure, &
         !$acc            density, tend, thermExpCoeff, &
         !$acc            salineContractCoeff, &
         !$acc            JacobianTz, JacobianSz, JacobianDxDs) &
         !$acc    private(cell1, cell2, invdcEdge, k, kMin, kMax, alpha, beta,&
         !$acc         TL, TR, SL, SR, Area, zStar, zC, zGamma, pGrad)
#else
         !$omp parallel
         !$omp do schedule(runtime) &
         !$omp private(cell1, cell2, invdcEdge, k, kMin, kMax, alpha, beta, &
         !$omp         TL, TR, SL, SR, Area, zStar, zC, zGamma, pGrad)
#endif
         do iEdge=1,nEdgesOwned
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            invdcEdge = 1.0_RKIND / dcEdge(iEdge)
            kMin = minLevelEdgeBot(iEdge)
            kMax = maxLevelEdgeTop(iEdge)

            ! compute J(T,z) and J(S,z)
            ! in Shchepetkin and McWilliams (2003) (7.16)

            JacobianTz(kMin,iEdge) = 0.0_RKIND
            JacobianSz(kMin,iEdge) = 0.0_RKIND

            do k=kMin+1,kMax

               ! eqn 2.7 in Shchepetkin and McWilliams (2003)
               ! Note delta x was removed.  It must be an error in the
               ! paper, ! as it makes the units incorrect.
               Area = 0.5_RKIND*(zMid(k-1,cell1) - zMid(k,cell1) + &
                                 zMid(k-1,cell2) - zMid(k,cell2))

               ! eqn 2.8
               zStar = ( zMid(k-1,cell2)*zMid(k-1,cell1) - &
                         zMid(k  ,cell2)*zMid(k  ,cell1) )/ &
                       ( zMid(k-1,cell2)-zMid(k  ,cell2) + &
                         zMid(k-1,cell1)-zMid(k,cell1))

               ! eqn 3.2
               zC = 0.25_RKIND*(zMid(k,cell1) + zMid(k-1,cell1) + &
                                zMid(k,cell2) + zMid(k-1,cell2))

               ! eqn 4.1
               zGamma = (1.0_RKIND - pGradLvlWgt)*zStar + pGradLvlWgt*zC


               TL = (tracers(indxT,k  ,cell1)*(zMid(k-1,cell1)-zGamma) + &
                     tracers(indxT,k-1,cell1)*(zGamma-zMid(k  ,cell1)))/ &
                    (zMid(k-1,cell1) - zMid(k,cell1))
               TR = (tracers(indxT,k  ,cell2)*(zMid(k-1,cell2)-zGamma) +  &
                     tracers(indxT,k-1,cell2)*(zGamma-zMid(k  ,cell2)))/ &
                    (zMid(k-1,cell2) - zMid(k,cell2))

               SL = (tracers(indxS,k  ,cell1)*(zMid(k-1,cell1)-zGamma) + &
                     tracers(indxS,k-1,cell1)*(zGamma-zMid(k  ,cell1)))/ &
                    (zMid(k-1,cell1) - zMid(k,cell1))
               SR = (tracers(indxS,k  ,cell2)*(zMid(k-1,cell2)-zGamma) +  &
                     tracers(indxS,k-1,cell2)*(zGamma-zMid(k  ,cell2)))/ &
                    (zMid(k-1,cell2) - zMid(k,cell2))


               ! eqn 2.6 in Shchepetkin and McWilliams (2003)
               JacobianTz(k,iEdge) = Area*(TL - TR)
               JacobianSz(k,iEdge) = Area*(SL - SR)
            end do

            ! In top layer, use pressure for generalized coordinates
            ! pGrad = -1/density_0 (grad p_k + density g grad z_k^{mid})

            k = kMin
            pGrad = edgeMask(k,iEdge)*invdcEdge*( &
                    - density0Inv*(pressure(k,cell2) - &
                                   pressure(k,cell1)) &
                    - gdensity0Inv*0.5_RKIND* &
                      (density(k,cell1)+density(k,cell2))* &
                      (zMid(k,cell2)- zMid(k,cell1) ) )

            tend(k,iEdge) = tend(k,iEdge) + pGrad

            do k=kMin+1,kMax

               ! Average alpha and beta over four data points of the Jacobian cell.
               ! Note that thermExpCoeff and salineContractCoeff include a 1/density factor,
               ! so must multiply by density here.
               alpha = 0.25_RKIND*( &
                       density(k  ,cell1)*thermExpCoeff (k  ,cell1) &
                     + density(k-1,cell1)*thermExpCoeff (k-1,cell1) &
                     + density(k  ,cell2)*thermExpCoeff (k  ,cell2) &
                     + density(k-1,cell2)*thermExpCoeff (k-1,cell2) )
               beta  = 0.25_RKIND*( &
                       density(k  ,cell1)*salineContractCoeff(k  ,cell1) &
                     + density(k-1,cell1)*salineContractCoeff(k-1,cell1) &
                     + density(k  ,cell2)*salineContractCoeff(k  ,cell2) &
                     + density(k-1,cell2)*salineContractCoeff(k-1,cell2) )

               ! Shchepetkin and McWilliams (2003) (7.16)
               JacobianDxDs(k,iEdge) = -alpha*JacobianTz(k,iEdge) + &
                                         beta*JacobianSz(k,iEdge)

               ! note JacobianDxDs includes negative sign, so
               ! pGrad is - g/rho_0 dP/dx

               pGrad = pGrad + gdensity0Inv * JacobianDxDs(k,iEdge) * invdcEdge

               tend(k,iEdge) = tend(k,iEdge) + pGrad

            end do
         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif

         !$acc exit data delete (JacobianDxDs, JacobianTz, JacobianSz)
         deallocate(JacobianDxDs, JacobianTz, JacobianSz)

      case (pGradTypeConstForced)

#ifdef MPAS_OPENACC
         !$acc parallel loop &
         !$acc    present(angleEdge,maxLevelEdgeTop, tend, edgeMask) &
         !$acc    private(pGrad, k)
#else
         !$omp parallel
         !$omp do schedule(runtime) &
         !$omp    private(pGrad, k)
#endif
         do iEdge=1,nEdgesOwned
            ! compute the pressure gradient to be applied at the edge
            pGrad = zonalSshGrad * cos(angleEdge(iEdge)) + meridSshGrad * sin(angleEdge(iEdge))

            do k=minLevelEdgeBot(iEdge),maxLevelEdgeTop(iEdge)
               tend(k,iEdge) = tend(k,iEdge) - gravity * edgeMask(k,iEdge) * pGrad
            end do
         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif

      case default

         !*** error already trapped during init

      end select ! pGradType

      call mpas_timer_stop("pressure grad")

   !--------------------------------------------------------------------

   end subroutine ocn_vel_pressure_grad_tend!}}}

!***********************************************************************
!
!  routine ocn_vel_pressure_grad_init
!
!> \brief   Initializes ocean momentum horizontal pressure gradient
!> \author  Mark Petersen
!> \date    September 2011
!> \details
!>  This routine initializes parameters required for the computation of
!>  the horizontal pressure gradient.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_pressure_grad_init(err)!{{{

      !-----------------------------------------------------------------
      ! Output Variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< [out] error flag

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      !*** Initialize error code and default values

      err = 0
      pGradOff    = .true.
      pGradType   = pGradTypeNone
      pGradLvlWgt = 0.0_RKIND

      density0Inv = 1.0_RKIND / rho_sw
      gdensity0Inv = gravity / rho_sw
      inv12 = 1.0_RKIND / 12.0_RKIND

      !*** Now reset variables based on input configuration

      if (config_disable_vel_pgrad) then
         pGradOff = .true.
         return
      else
         pGradOff = .false.
         pGradLvlWgt = config_common_level_weight
      endif

      select case (trim(config_pressure_gradient_type))

      case ('ssh_gradient')
         pGradType = pGradTypeSSHgrad
         call mpas_log_write(' Pressure type is: ssh_gradient')

      case ('pressure_and_zmid')
         pGradType = pGradTypePZmid
         call mpas_log_write(' Pressure type is: pressure_and_zmid')

      case ('MontgomeryPotential')
         pGradType = pGradTypeMontPot
         call mpas_log_write(' Pressure type is: MontgomeryPotential')

      case ('MontgomeryPotential_and_density')
         pGradType = pGradTypeMontPotDens
         call mpas_log_write( &
            ' Pressure type is: MontgomeryPotential_and_density')

      case ('Jacobian_from_density')
         pGradType = pGradTypeJacobDens
         call mpas_log_write(' Pressure type is: Jacobian_from_density')

      case ('Jacobian_from_TS')
         pGradType = pGradTypeJacobTS
         call mpas_log_write(' Pressure type is: Jacobian_from_TS')

      case ('constant_forced')
         pGradType = pGradTypeConstForced
         call mpas_log_write(' Pressure type is: constant_forced')

      case default
         call mpas_log_write( &
            ' Incorrect choice of config_pressure_gradient_type.', &
            MPAS_LOG_CRIT)
         err = 1

      end select

      if ( (config_time_integrator == 'LTS') &
           .or. (config_time_integrator == 'FB_LTS') ) then
         timeIntegratorLTS = .true.
      else
         timeIntegratorLTS = .false.
      endif
   !--------------------------------------------------------------------

   end subroutine ocn_vel_pressure_grad_init!}}}

!***********************************************************************

end module ocn_vel_pressure_grad

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
