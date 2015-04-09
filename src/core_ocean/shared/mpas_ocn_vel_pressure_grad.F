! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
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
!>  tendencie from the horizontal pressure gradient.
!>
!
!-----------------------------------------------------------------------

module ocn_vel_pressure_grad

   use mpas_grid_types
   use mpas_constants
   use mpas_io_units

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

   public :: ocn_vel_pressure_grad_tend, &
             ocn_vel_pressure_grad_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   character (len=StrKIND), pointer :: config_pressure_gradient_type
   real (kind=RKIND), pointer :: config_common_level_weight
   logical :: pgradOn
   real (kind=RKIND) :: density0Inv, gdensity0Inv, inv12


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

   subroutine ocn_vel_pressure_grad_tend(meshPool, pressure, montgomeryPotential, zMid, density, potentialDensity, &
      indexT, indexS, tracers, tend, err, inSituThermalExpansionCoeff,inSituSalineContractionCoeff)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      integer, intent(in) :: indexT, indexS

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         pressure, & !< Input: Pressure field
         montgomeryPotential, & !< Input: Mongomery potential
         zMid, &     !< Input: z-coordinate at mid-depth of layer
         density, &      !< Input: density
         potentialDensity         !< Input: potentialDensity

      real (kind=RKIND), dimension(:,:), intent(in), optional :: &
         inSituThermalExpansionCoeff, &  
         inSituSalineContractionCoeff

      real (kind=RKIND), dimension(:,:,:), intent(in) :: tracers

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend          !< Input/Output: velocity tendency

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iEdge, k, cell1, cell2, iCell, kMax
      integer, pointer :: nVertLevels, nCells, nEdgesSolve
      integer, dimension(:), pointer :: maxLevelEdgeTop, maxLevelCell
      integer, dimension(:,:), pointer :: cellsOnEdge, edgeMask

      real (kind=RKIND), dimension(:), pointer :: dcEdge
      real (kind=RKIND), dimension(:), allocatable :: JacobianDxDs,JacobianTz,JacobianSz,T1,T2,S1,S2
      real (kind=RKIND), dimension(:,:), allocatable :: FXTop, work1, work2
      real (kind=RKIND) :: invdcEdge, pGrad, sumAJTop, AJTop, FC, FCPrev, alpha, beta

      err = 0

      if (.not. pgradOn) return

      call mpas_pool_get_dimension(meshPool, 'nVertLevels',nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'nCells',nCells)
      call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)

      if (config_pressure_gradient_type.eq.'pressure_and_zmid') then

         ! pressure for generalized coordinates
         ! -1/density_0 (grad p_k + density g grad z_k^{mid})

         do iEdge=1,nEdgesSolve
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            invdcEdge = 1.0 / dcEdge(iEdge)

            do k=1,maxLevelEdgeTop(iEdge)
               tend(k,iEdge) = tend(k,iEdge) + edgeMask(k,iEdge) * invdcEdge * ( &
                 - density0Inv * ( pressure(k,cell2) - pressure(k,cell1) ) &
                 - gdensity0Inv * 0.5*(density(k,cell1)+density(k,cell2)) * ( zMid(k,cell2) - zMid(k,cell1) ) )
            end do
         end do

      elseif (config_pressure_gradient_type.eq.'MontgomeryPotential') then

         ! For pure isopycnal coordinates, this is just grad(M), 
         ! the gradient of Montgomery Potential

         do iEdge=1,nEdgesSolve
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            invdcEdge = 1.0 / dcEdge(iEdge)

            do k=1,maxLevelEdgeTop(iEdge)
               tend(k,iEdge) = tend(k,iEdge) + edgeMask(k,iEdge) * invdcEdge * ( &
                  - ( montgomeryPotential(k,cell2) - montgomeryPotential(k,cell1) ) )
            end do
         end do

      elseif (config_pressure_gradient_type.eq.'MontgomeryPotential_and_density') then

         ! This formulation has not been extensively tested and is not supported at this time.

         ! This is -grad(M)+p grad(1/rho)
         ! Where rho is the potential density.
         ! See Bleck (2002) equation 1, and last equation in Appendix A.

         do iEdge=1,nEdgesSolve
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            invdcEdge = 1.0 / dcEdge(iEdge)

            do k=1,maxLevelEdgeTop(iEdge)
               tend(k,iEdge) = tend(k,iEdge) + edgeMask(k,iEdge) * invdcEdge * ( &
                  - ( montgomeryPotential(k,cell2) - montgomeryPotential(k,cell1) ) &
                  +  0.5*(pressure(k,cell1)+pressure(k,cell2)) * ( 1.0/potentialDensity(k,cell2) - 1.0/potentialDensity(k,cell1) ) )
            end do
         end do

      elseif (config_pressure_gradient_type.eq.'Jacobian_from_density') then

         allocate(JacobianDxDs(nVertLevels))

         do iEdge=1,nEdgesSolve
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            invdcEdge = 1.0 / dcEdge(iEdge)

            call pGrad_Jacobian_common_level(density(:,cell1),density(:,cell2),zMid(:,cell1),zMid(:,cell2), &
               maxLevelEdgeTop(iEdge), config_common_level_weight, JacobianDxDs)

            ! In layer 1, use pressure for generalized coordinates
            ! pGrad = -1/density_0 (grad p_k + density g grad z_k^{mid})
            k = 1
            pGrad = edgeMask(k,iEdge) * invdcEdge * ( &
                 - density0Inv * ( pressure(k,cell2) - pressure(k,cell1) ) &
                 - gdensity0Inv * 0.5*(density(k,cell1)+density(k,cell2)) * ( zMid(k,cell2) - zMid(k,cell1) ) )

            tend(k,iEdge) = tend(k,iEdge) + pGrad

            do k=2,maxLevelEdgeTop(iEdge)

               ! note JacobianDxDs includes negative sign, so
               ! pGrad is - g/rho_0 dP/dx

               pGrad = pGrad + gdensity0Inv * JacobianDxDs(k) * invdcEdge

               tend(k,iEdge) = tend(k,iEdge) + pGrad

            end do
         end do

         deallocate(JacobianDxDs)

     elseif (config_pressure_gradient_type.eq.'Jacobian_from_TS') then

         allocate(JacobianDxDs(nVertLevels),JacobianTz(nVertLevels),JacobianSz(nVertLevels), T1(nVertLevels), T2(nVertLevels), S1(nVertLevels), S2(nVertLevels))

         do iEdge=1,nEdgesSolve
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            invdcEdge = 1.0 / dcEdge(iEdge)
            kMax = maxLevelEdgeTop(iEdge)

            ! copy T and S to local column arrays
            T1(1:kMax) = tracers(indexT,1:kMax,cell1)
            T2(1:kMax) = tracers(indexT,1:kMax,cell2)
            S1(1:kMax) = tracers(indexS,1:kMax,cell1)
            S2(1:kMax) = tracers(indexS,1:kMax,cell2)

            ! compute J(T,z) and J(S,z) in Shchepetkin and McWilliams (2003) (7.16)
            call pGrad_Jacobian_common_level(T1, T2 ,zMid(:,cell1),zMid(:,cell2),kMax,config_common_level_weight, JacobianTz)
            call pGrad_Jacobian_common_level(S1, S2 ,zMid(:,cell1),zMid(:,cell2),kMax,config_common_level_weight, JacobianSz)

            ! In layer 1, use pressure for generalized coordinates
            ! pGrad = -1/density_0 (grad p_k + density g grad z_k^{mid})
            k = 1
            pGrad = edgeMask(k,iEdge) * invdcEdge * ( &
                 - density0Inv * ( pressure(k,cell2) - pressure(k,cell1) ) &
                 - gdensity0Inv * 0.5*(density(k,cell1)+density(k,cell2)) * ( zMid(k,cell2) - zMid(k,cell1) ) )

            tend(k,iEdge) = tend(k,iEdge) + pGrad

            do k=2,kMax

               ! Average alpha and beta over four data points of the Jacobian cell.
               ! Note that inSituThermalExpansionCoeff and inSituSalineContractionCoeff include a 1/density factor,
               ! so must multiply by density here.
               alpha = 0.25*(  density(k,cell1)*inSituThermalExpansionCoeff (k,cell1) + density(k-1,cell1)*inSituThermalExpansionCoeff (k-1,cell1) &
                             + density(k,cell2)*inSituThermalExpansionCoeff (k,cell2) + density(k-1,cell2)*inSituThermalExpansionCoeff (k-1,cell2) )
               beta  = 0.25*(  density(k,cell1)*inSituSalineContractionCoeff(k,cell1) + density(k-1,cell1)*inSituSalineContractionCoeff(k-1,cell1) &
                             + density(k,cell2)*inSituSalineContractionCoeff(k,cell2) + density(k-1,cell2)*inSituSalineContractionCoeff(k-1,cell2) )

               ! Shchepetkin and McWilliams (2003) (7.16)
               JacobianDxDs(k) = -alpha*JacobianTz(k) + beta*JacobianSz(k)

               ! note JacobianDxDs includes negative sign, so
               ! pGrad is - g/rho_0 dP/dx

               pGrad = pGrad + gdensity0Inv * JacobianDxDs(k) * invdcEdge

               tend(k,iEdge) = tend(k,iEdge) + pGrad

            end do
         end do

         deallocate(JacobianDxDs,JacobianTz,JacobianSz, T1, T2, S1, S2)

      else

         write (stderrUnit,'(a,a)') ' Pressure type is: ',trim(config_pressure_gradient_type)
         write (stderrUnit,*) ' Incorrect choice of config_pressure_gradient_type.'
         err = 1

      endif

   !--------------------------------------------------------------------

   end subroutine ocn_vel_pressure_grad_tend!}}}

!***********************************************************************
!
!  routine pGrad_Jacobian_common_level
!
!> \brief   Computes density-Jacobian
!> \author  Mark Petersen
!> \date    February 2014
!> \details 
!>  This routine computes the density-Jacobian in common_level form.
!>  See Shchepetkin and McWilliams (2003) Ocean Modeling, sections 2-4
!
!-----------------------------------------------------------------------

   subroutine pGrad_Jacobian_common_level(rho1,rho2,z1,z2,kMax,gamma,JacobianDxDs)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(in) :: &
         rho1, & ! density of column 1 
         rho2, & ! density of column 2
         z1,   & ! z-coordinate at middle of cell, column 1 
         z2      ! z-coordinate at middle of cell, column 2 

      real (kind=RKIND), intent(in) :: &
         gamma   ! weight between zStar (original Jacobian) and z_C (weighted Jacobian)

      integer, intent(in) :: &
         kMax

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(out) :: &
         JacobianDxDs  ! - Delta x Delta s J(rho,z)

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: k
      real (kind=RKIND) :: Area, zStar, rhoL, rhoR, zC, zGamma

      JacobianDxDs = 0.0

      do k=2,kMax

         ! eqn 2.7 in Shchepetkin and McWilliams (2003)
         ! Note delta x was removed.  It must be an error in the paper,
         ! as it makes the units incorrect.
         Area = 0.5*(z1(k-1) - z1(k) + z2(k-1) - z2(k) ) 

         ! eqn 2.8
         zStar = ( z2(k-1)*z1(k-1) - z2(k)*z1(k) )/(z2(k-1)-z2(k) + z1(k-1)-z1(k))

         ! eqn 3.2
         zC = 0.25*( z1(k) + z1(k-1) + z2(k) + z2(k-1) )

         ! eqn 4.1
         zGamma = (1.0 - gamma)*zStar + gamma*zC

         rhoL = (rho1(k)*(z1(k-1)-zGamma) + rho1(k-1)*(zGamma-z1(k)))/(z1(k-1) - z1(k))
         rhoR = (rho2(k)*(z2(k-1)-zGamma) + rho2(k-1)*(zGamma-z2(k)))/(z2(k-1) - z2(k))

         ! eqn 2.6 in Shchepetkin and McWilliams (2003)
         JacobianDxDs(k) = Area * (rhoL - rhoR)
      end do

   end subroutine pGrad_Jacobian_common_level

!***********************************************************************
!
!  routine pGrad_Jacobian_POM_SCRUM
!
!> \brief   Computes density-Jacobian
!> \author  Mark Petersen
!> \date    February 2014
!> \details 
!>  This routine computes the density-Jacobian in POM/SCRUM form.
!>  See Shchepetkin and McWilliams (2003) Ocean Modeling, section 2.
!
!-----------------------------------------------------------------------

   subroutine pGrad_Jacobian_POM_SCRUM(rho1,rho2,z1,z2,kMax,JacobianDxDs)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(in) :: &
         rho1, & ! density of column 1 
         rho2, & ! density of column 2
         z1,   & ! z-coordinate at middle of cell, column 1 
         z2      ! z-coordinate at middle of cell, column 2 

      integer, intent(in) :: &
         kMax  ! maximum level

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(out) :: &
         JacobianDxDs  ! - Delta x Delta s J(rho,z)

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: k

      JacobianDxDs = 0.0

      do k=2,kMax

         ! eqn 2.3 in Shchepetkin and McWilliams (2003)
         JacobianDxDs(k) = 0.25*(&
              (rho1(k) + rho1(k-1) - rho2(k) - rho2(k-1) )*(z1(k-1) - z1(k) + z2(k-1) - z2(k) ) & 
            - (rho1(k-1) - rho1(k) + rho2(k-1) - rho2(k) )*(z1(k) + z1(k-1) - z2(k) - z2(k-1) ) )
      end do

   end subroutine pGrad_Jacobian_POM_SCRUM

!***********************************************************************
!
!  routine pGrad_Jacobian_diagonal
!
!> \brief   Computes density-Jacobian
!> \author  Mark Petersen
!> \date    February 2014
!> \details 
!>  This routine computes the density-Jacobian in diagonal form.
!>  See Shchepetkin and McWilliams (2003) Ocean Modeling, section 2.
!
!-----------------------------------------------------------------------

   subroutine pGrad_Jacobian_diagonal(rho1,rho2,z1,z2,kMax,JacobianDxDs)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(in) :: &
         rho1, & ! density of column 1 
         rho2, & ! density of column 2
         z1,   & ! z-coordinate at middle of cell, column 1 
         z2      ! z-coordinate at middle of cell, column 2 

      integer, intent(in) :: &
         kMax

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(out) :: &
         JacobianDxDs  ! - Delta x Delta s J(rho,z)

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: k

      JacobianDxDs = 0.0

      do k=2,kMax

         ! eqn 2.5 in Shchepetkin and McWilliams (2003)
         JacobianDxDs(k) = 0.5*( &
              (rho1(k-1) - rho2(k))*(z2(k-1) - z1(k) ) & 
            + (rho1(k) - rho2(k-1))*(z1(k-1) - z2(k)) )
      end do

   end subroutine pGrad_Jacobian_diagonal

!***********************************************************************
!
!  routine pGrad_Jacobian_pseudo_flux
!
!> \brief   Computes density-Jacobian
!> \author  Mark Petersen
!> \date    February 2014
!> \details 
!>  This routine computes the density-Jacobian in pseudo_flux form.
!>  See Shchepetkin and McWilliams (2003) Ocean Modeling, section 2.
!
!-----------------------------------------------------------------------

   subroutine pGrad_Jacobian_pseudo_flux(rho1,rho2,z1,z2,kMax,JacobianDxDs)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(in) :: &
         rho1, & ! density of column 1 
         rho2, & ! density of column 2
         z1,   & ! z-coordinate at middle of cell, column 1 
         z2      ! z-coordinate at middle of cell, column 2 

      integer, intent(in) :: &
         kMax

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(out) :: &
         JacobianDxDs  ! - Delta x Delta s J(rho,z)

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: k
      real (kind=RKIND) :: FLeft, FTop, FRight, FBottom 

      JacobianDxDs = 0.0

      do k=2,kMax

         FLeft   = 0.5*( rho1(k) + rho1(k-1) ) * (z1(k-1) - z1(k))
         FTop    = 0.5*( rho1(k-1) + rho2(k-1) ) * (z2(k-1) - z1(k-1))
         FRight  = 0.5*( rho2(k) + rho2(k-1) ) * (z2(k-1) - z2(k))
         FBottom = 0.5*( rho1(k) + rho2(k) ) * (z2(k) - z1(k))

         ! eqn 2.11 in Shchepetkin and McWilliams (2003)
         JacobianDxDs(k) = FLeft + FTop - FRight - FBottom 
      end do

   end subroutine pGrad_Jacobian_pseudo_flux

!***********************************************************************
!
!  routine ocn_vel_pressure_grad_init
!
!> \brief   Initializes ocean momentum horizontal pressure gradient
!> \author  Mark Petersen
!> \date    September 2011
!> \details 
!>  This routine initializes parameters required for the computation of the
!>  horizontal pressure gradient.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_pressure_grad_init(err)!{{{

   !--------------------------------------------------------------------


      !-----------------------------------------------------------------
      !
      ! Output Variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag


      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------
      real (kind=RKIND), pointer :: config_density0
      logical, pointer :: config_disable_vel_pgrad

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_pressure_gradient_type', config_pressure_gradient_type)
      call mpas_pool_get_config(ocnConfigs, 'config_common_level_weight', config_common_level_weight)
      call mpas_pool_get_config(ocnConfigs, 'config_density0', config_density0)
      call mpas_pool_get_config(ocnConfigs, 'config_disable_vel_pgrad', config_disable_vel_pgrad)

      pgradOn = .true.

      density0Inv = 1.0/config_density0
      gdensity0Inv = gravity/config_density0
      inv12 = 1.0/12.0

      if (config_disable_vel_pgrad) pgradOn = .false.

      write (stdoutUnit,'(a,a)') ' Pressure type is: ',trim(config_pressure_gradient_type)

   !--------------------------------------------------------------------

   end subroutine ocn_vel_pressure_grad_init!}}}

!***********************************************************************

end module ocn_vel_pressure_grad

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
