! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_equation_of_state_linear
!
!> \brief MPAS ocean linear equation of state
!> \author Mark Petersen, Todd Ringler, modified Phil Jones
!> \date   September 2011, modified May 2018
!> \details
!>  This module contains routines for computing density of sea water
!>  using a linear approximation.
!
!-----------------------------------------------------------------------

module ocn_equation_of_state_linear

   use mpas_kind_types
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

   public :: ocn_equation_of_state_linear_density, &
             ocn_equation_of_state_linear_init

   interface ocn_equation_of_state_linear_density
      module procedure ocn_equation_of_state_linear_density_only
      module procedure ocn_equation_of_state_linear_density_exp
   end interface

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   real (kind=RKIND) ::       &
      ocnEqStateLinearRhoRef, &! reference density
      ocnEqStateLinearAlpha,  &! scalar temperature expansion coeff
      ocnEqStateLinearBeta,   &! scalar salinity  contraction coeff
      ocnEqStateLinearTref,   &! reference density
      ocnEqStateLinearSref     ! reference density

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_equation_of_state_linear_density_only
!
!> \brief   Computes density using a linear equation of state
!> \author  Mark Petersen, Todd Ringler, modified by Phil Jones
!> \date    September 2011, modified May 2018
!> \details
!>  This routine uses a linear equation of state to update the density.
!>  The density is a linear perturbation from a reference density
!>  at reference T,S values using fixed (linear) expansion coefficients.
!>
!>  While somewhat unnecessary, we make the interface and capability
!>  of linear eos to be identical to nonlinear eos
!>
!>  Density can be computed in-situ using k_displaced=0 and
!>      displacement_type = 'relative'.
!>
!>  Potential density (referenced to top layer) can be computed
!>      using k_displaced=1 and displacement_type = 'absolute'.
!>
!>  The density of SST/SSS after adiabatic displacement to each layer
!>      can be computed using displacement_type = 'surfaceDisplaced'.
!>
!>  When using displacement_type = 'surfaceDisplaced', k_displaced is
!>      ignored and tracersSurfaceLayerValue must be present.
!
!-----------------------------------------------------------------------

   subroutine ocn_equation_of_state_linear_density_only(          &
              nVertLevels, nCells, kDisplaced, displacementType,  &
              indexT, indexS, tracers, density, err,              &
              tracersSurfaceLayerValue)
      !{{{
      !-----------------------------------------------------------------
      !
      ! Input variables
      !
      !-----------------------------------------------------------------

      integer, intent(in) ::     &
         nCells,                 &! num of horizontal cells
         nVertLevels,            &! num of vertical levels
         kDisplaced,             &! num of levels to displace parcel
         indexT, indexS           ! index into tracer array for T,S

      character(len=*), intent(in) :: &
         displacementType

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         tracers                  ! array containing T,S for calculation

      real (kind=RKIND), dimension(:,:), intent(in), optional :: &
         tracersSurfaceLayerValue ! optional sfc tracer values for sfc
                                  ! displacement option

      !-----------------------------------------------------------------
      !
      ! Output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: &
         err                   ! error flag

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         density               ! computed density

      !-----------------------------------------------------------------
      !
      ! Local variables
      !
      !-----------------------------------------------------------------

      integer :: &
         iCell, k    ! cell and level loop indices

      !-----------------------------------------------------------------

      !*** initialize error flag

      err = 0

      !*** check for displacement choice
      !*** for linear EOS, only choice that impacts the
      !*** density is surface displacement

      ! ignore this test for performance reasons
      ! test of request to address out of bounds
      !if (displacement_type_local == 'absolute' .and.   &
      !   (k_displaced_local <= 0 .or. k_displaced_local > nVertLevels) ) then
      !   call mpas_log_write('Abort: In equation_of_state_linear' // &
      !       ' k_displaced must be between 1 and nVertLevels for ' // &
      !       'displacement_type = absolute', MPAS_LOG_CRIT)
      !endif

      ! copy some intent(in) into local work space

      !*** if surfaceDisplaced, then compute density at all levels
      !*** based on input surface values

      if (trim(displacementType) == 'surfaceDisplaced') then

#ifdef MPAS_OPENACC
         !$acc parallel loop gang vector collapse(2) &
         !$acc& present(density, tracersSurfaceLayerValue)
#else
         !$omp parallel
         !$omp do schedule(runtime) private(k)
#endif
         do iCell = 1, nCells
         do k = 1, nVertLevels

            density(k,iCell) =  ocnEqStateLinearRhoRef &
                             - ocnEqStateLinearAlpha * &
                               (tracersSurfaceLayerValue(indexT,iCell) - &
                                ocnEqStateLinearTref)  &
                             + ocnEqStateLinearBeta  * &
                               (tracersSurfaceLayerValue(indexS,iCell) - &
                                ocnEqStateLinearSref)
         end do
         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif

      else  ! all other displacement types

#ifdef MPAS_OPENACC
         !$acc parallel loop gang vector collapse(2) &
         !$acc& present(density, tracers)
#else
         !$omp parallel
         !$omp do schedule(runtime) private(k)
#endif
         do iCell = 1, nCells
         do k = 1, nVertLevels

            density(k,iCell) = ocnEqStateLinearRhoRef     &
                             - ocnEqStateLinearAlpha *    &
                               (tracers(indexT,k,iCell) - &
                                ocnEqStateLinearTref)     &
                             + ocnEqStateLinearBeta  *    &
                               (tracers(indexS,k,iCell) - &
                                ocnEqStateLinearSref)
         end do
         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif

      endif

      !-----------------------------------------------------------------

   end subroutine ocn_equation_of_state_linear_density_only!}}}

!***********************************************************************
!
!  routine ocn_equation_of_state_linear_density_exp
!
!> \brief   Computes density using a linear equation of state
!> \author  Mark Petersen, Todd Ringler, modified by Phil Jones
!> \date    September 2011, modified May 2018
!> \details
!>  This routine uses a linear equation of state to update the density.
!>  This instance also returns the expansion coefficients.
!>
!-----------------------------------------------------------------------

   subroutine ocn_equation_of_state_linear_density_exp(           &
              nVertLevels, nCells, kDisplaced, displacementType,  &
              indexT, indexS, tracers, density, err,              &
              thermalExpansionCoeff, salineContractionCoeff,      &
              tracersSurfaceLayerValue)
      !{{{
      !-----------------------------------------------------------------
      !
      ! Input variables
      !
      !-----------------------------------------------------------------

      integer, intent(in) ::     &
         nCells,                 &! num of horizontal cells
         nVertLevels,            &! num of vertical levels
         kDisplaced,             &! num of levels to displace parcel
         indexT, indexS           ! index into tracer array for T,S

      character(len=*), intent(in) :: &
         displacementType

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         tracers                  ! array containing T,S for calculation

      real (kind=RKIND), dimension(:,:), intent(in), optional :: &
         tracersSurfaceLayerValue ! optional sfc tracer values for sfc
                                  ! displacement option

      !-----------------------------------------------------------------
      !
      ! Output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: &
         err                   ! error flag

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         density               ! computed density

      ! Optional expansion (contraction) coefficients
      ! alpha defined as $-1/\rho d\rho/dT$ (note negative sign)
      ! beta  defined as $1/\rho d\rho/dS$
      real (kind=RKIND), dimension(:,:), intent(out) :: &
         thermalExpansionCoeff,  &! Thermal expansion coeff (alpha)
         salineContractionCoeff   ! Saline contraction coeff (beta)

      !-----------------------------------------------------------------
      !
      ! Local variables
      !
      !-----------------------------------------------------------------

      integer :: &
         iCell, k    ! cell and level loop indices

      !-----------------------------------------------------------------

      !*** initialize error flag

      err = 0

      !*** check for displacement choice
      !*** for linear EOS, only choice that impacts the
      !*** density is surface displacement

      ! ignore this test for performance reasons
      ! test of request to address out of bounds
      !if (displacement_type_local == 'absolute' .and.   &
      !   (k_displaced_local <= 0 .or. k_displaced_local > nVertLevels) ) then
      !   call mpas_log_write('Abort: In equation_of_state_linear' // &
      !       ' k_displaced must be between 1 and nVertLevels for ' // &
      !       'displacement_type = absolute', MPAS_LOG_CRIT)
      !endif

      ! copy some intent(in) into local work space

      !*** if surfaceDisplaced, then compute density at all levels
      !*** based on input surface values
      !*** return fixed values of expansion coefficients

      if (trim(displacementType) == 'surfaceDisplaced') then

#ifdef MPAS_OPENACC
         !$acc parallel loop gang vector collapse(2) &
         !$acc& present(density, &
         !$acc&         tracersSurfaceLayerValue, &
         !$acc&         salineContractionCoeff, &
         !$acc&         thermalExpansionCoeff)
#else
         !$omp parallel
         !$omp do schedule(runtime) private(k)
#endif
         do iCell = 1, nCells
         do k = 1, nVertLevels

            density(k,iCell) =  ocnEqStateLinearRhoRef &
                             - ocnEqStateLinearAlpha * &
                               (tracersSurfaceLayerValue(indexT,iCell) - &
                                ocnEqStateLinearTref)  &
                             + ocnEqStateLinearBeta  * &
                               (tracersSurfaceLayerValue(indexS,iCell) - &
                                ocnEqStateLinearSref)

            thermalExpansionCoeff(k,iCell) = ocnEqStateLinearAlpha/ &
                                             density(k,iCell)

            salineContractionCoeff(k,iCell) = ocnEqStateLinearBeta/ &
                                              density(k,iCell)
         end do
         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif

      else  ! all other displacement types

#ifdef MPAS_OPENACC
         !$acc parallel loop gang vector collapse(2) &
         !$acc& present(density, &
         !$acc&         tracers, &
         !$acc&         salineContractionCoeff, &
         !$acc&         thermalExpansionCoeff)
#else
         !$omp parallel
         !$omp do schedule(runtime) private(k)
#endif
         do iCell = 1, nCells
         do k = 1, nVertLevels

            density(k,iCell) = ocnEqStateLinearRhoRef     &
                             - ocnEqStateLinearAlpha *    &
                               (tracers(indexT,k,iCell) - &
                                ocnEqStateLinearTref)     &
                             + ocnEqStateLinearBeta  *    &
                               (tracers(indexS,k,iCell) - &
                                ocnEqStateLinearSref)

            thermalExpansionCoeff(k,iCell) = ocnEqStateLinearAlpha/ &
                                             density(k,iCell)

            salineContractionCoeff(k,iCell) = ocnEqStateLinearBeta/ &
                                              density(k,iCell)
         end do
         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif

      endif

      !-----------------------------------------------------------------

   end subroutine ocn_equation_of_state_linear_density_exp!}}}

!***********************************************************************
!
!  routine ocn_equation_of_state_linear_init
!
!> \brief   Initializes linear ocean equation of state
!> \author  Mark Petersen, Todd Ringler, modified by Phil Jones
!> \date    September 2011, modified May 2018
!> \details
!>  This routine initializes a variety of quantities for computing
!>  density from temperature and salinity using a linear approximation.
!
!-----------------------------------------------------------------------

   subroutine ocn_equation_of_state_linear_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: &
         err                   ! error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------

      !*** Initialize error flag

      err = 0

      !***
      !*** Extract input namelist options into module variables
      !*** Needed both to rename and because scalar pointers can
      !*** cause problems for some compilers
      !***

      ocnEqStateLinearRhoRef = config_eos_linear_densityref
      ocnEqStateLinearAlpha  = config_eos_linear_alpha
      ocnEqStateLinearBeta   = config_eos_linear_beta
      ocnEqStateLinearTref   = config_eos_linear_Tref
      ocnEqStateLinearSref   = config_eos_linear_Sref

   !--------------------------------------------------------------------

   end subroutine ocn_equation_of_state_linear_init!}}}

!***********************************************************************

end module ocn_equation_of_state_linear

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
