! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_equation_of_state_wright
!
!> \brief MPAS ocean Wright (1997) equation of state
!> \author Xylar Asay-Davis
!> \date   June 30, 2019
!> \details
!>  This module contains the routines for computing density from
!>  temperature, salinity and depth using an equation of state by Wright (1997)
!>  that approximates the nonlinear equation of state using a set of
!>  rational polynomials.
!>  See Wright (1997), doi: 10.1175/1520-0426(1997)014<0735:AEOSFU>2.0.CO;2
!
!-----------------------------------------------------------------------

module ocn_equation_of_state_wright

   use mpas_kind_types
   use mpas_pool_routines
   use mpas_constants
   use ocn_constants
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

   public :: ocn_equation_of_state_wright_density, &
             ocn_equation_of_state_wright_init

   !*** generic interface for case of density only or density and
   !*** expansion coeffs
   interface ocn_equation_of_state_wright_density
      module procedure ocn_equation_of_state_wright_density_only
      module procedure ocn_equation_of_state_wright_density_exp
   end interface

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   !*** temporary array to hold modified T, S, p for EOS calculation
   !*** declared here to make it thread-shared - could be replaced
   !***   by a local subroutine temporary in a different thread model

   real (kind=RKIND), dimension(:,:), allocatable :: &
      tracerTemp, tracerSalt, boussinesqPres
#ifdef MPAS_OPENACC
   !$acc declare device_resident(tracerTemp, tracerSalt, boussinesqPres)
#endif

   !*** valid range of T,S for Wright (1997) EOS
   !***   allowing T to become more negative because of sub-ice-shelf cavities

   real (kind=RKIND), parameter :: &
      ocnEqStateTmin = -3.0_RKIND, &! valid pot. temp. range
      ocnEqStateTmax = 30.0_RKIND, &
      ocnEqStateSmin = 28.0_RKIND, &! valid salinity, in psu
      ocnEqStateSmax = 38.0_RKIND

   !***  Wright (1997) constants from Table 1, last column

   real (kind=RKIND), parameter ::      &
      a0 = 7.057924e-4_RKIND,           &
      a1 = 3.480336e-7_RKIND,           &
      a2 = -1.112733e-7_RKIND

   real (kind=RKIND), parameter ::      &
      b0 = 5.790749e8_RKIND,            &
      b1 = 3.516535e6_RKIND,            &
      b2 = -4.002714e4_RKIND,           &
      b3 = 2.084372e2_RKIND,            &
      b4 = 5.944068e5_RKIND,            &
      b5 = -9.643486e3_RKIND

   real (kind=RKIND), parameter ::      &
      c0 = 1.704853e5_RKIND,            &
      c1 = 7.904722e2_RKIND,            &
      c2 = -7.984422_RKIND,             &
      c3 = 5.140652e-2_RKIND,           &
      c4 = -2.302158e2_RKIND,           &
      c5 = -3.079464_RKIND

      real (kind=RKIND) eos_wright_ref_pressure !  = config_eos_wright_ref_pressure

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_equation_of_state_wright_density_only
!
!> \brief   Computes equation of state
!> \author  Xylar Asay-Davis
!> \date    June 30, 2019
!> \details
!>  This routine computes the density from model temperature,
!>  salinity, and depth using potential-temperature formulation from
!>  Wright (1997), doi: 10.1175/1520-0426(1997)014<0735:AEOSFU>2.0.CO;2
!>  and a simple linearization of pressure p = -rho0*g*z with the Boussinesq
!>  reference density
!>
!>  Density can be computed in-situ using kDisplaced=0 and
!>      displacementType = 'relative'.
!>
!>  Potential density (referenced to zero pressure) can be computed
!>      using displacementType = 'absolute'. kDisplaced is ignored.
!>
!>  The density of SST/SSS after adiabatic displacement to each layer
!>      can be computed using displacementType = 'surfaceDisplaced'.
!>      kDisplaced is ignored and tracersSurfaceLayerValue must be present.
!
!-----------------------------------------------------------------------

   subroutine ocn_equation_of_state_wright_density_only(nVertLevels,  &
                               nCells, kDisplaced, displacementType,  &
                               indexT, indexS, tracers, zMid,         &
                               maxLevelCell, density, err, &
                               tracersSurfaceLayerValue)

   !{{{
   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------

      integer, intent(in) :: &
         nCells,             &! number of horizontal cells
         nVertLevels,        &! max number of verical levels
         kDisplaced,         &! target layer for displacement
         indexT,             &! temperature index in tracer array
         indexS               ! salinity    index in tracer array

      character(len=*), intent(in) :: &
         displacementType     ! choice of displacement

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         tracers              ! array of tracers including T,S

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         zMid

      integer, dimension(:), intent(in) :: &
         maxLevelCell

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         tracersSurfaceLayerValue

      !-----------------------------------------------------------------
      ! Output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err  ! error flag

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         density         ! computed density

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------


      integer :: &
         iCell, k ! cell and vertical level loop indices

      real (kind=RKIND) :: &
         T, S, T2, T3,     &! adjusted T, S, T squared, T cubed,
         p,                &! Boussinesq pressure at reference level,
         alpha0, lambda0,  &! functions of theta and S defined by Wright (1997)
         p0

      !-----------------------------------------------------------------

      !*** initialize error flag

      err = 0

      allocate(tracerTemp(nVertLevels, nCells))
      allocate(tracerSalt(nVertLevels, nCells))
      allocate(boussinesqPres(nVertLevels, nCells))

      if (displacementType == 'surfaceDisplaced') then
         call compute_surface_displaced_T_S(nVertLevels, nCells, indexT,      &
                                            indexS, tracersSurfaceLayerValue, &
                                            tracerTemp, tracerSalt)
      else
         call compute_bounded_T_S(nVertLevels, nCells, indexT, indexS, tracers, &
                                  tracerTemp, tracerSalt)
      end if

      call compute_boussinesq_p(nVertLevels, nCells, kDisplaced,      &
                                displacementType, zMid, maxLevelCell, &
                                boussinesqPres)

#ifdef MPAS_OPENACC
      !$acc parallel loop gang vector collapse(2) &
      !$acc&   present(density)
#else
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp private(k, S, T, p, T2, T3, p0, lambda0, alpha0)
#endif

      do iCell=1,nCells
         do k=1,nVertLevels

            S = tracerSalt(k, iCell)
            T = tracerTemp(k, iCell)
            p = boussinesqPres(k, iCell)

            T2 = T*T
            T3 = T*T2

            p0 = b0 + b1*T + b2*T2 + b3*T3 + b4*S + b5*T*S
            lambda0 = c0 + c1*T + c2*T2 + c3*T3 + c4*S + c5*T*S
            alpha0 = a0 + a1*T + a2*S

            density(k, iCell) = (p + p0)/(lambda0 + alpha0*(p + p0))

         end do
      end do

#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      deallocate(tracerTemp)
      deallocate(tracerSalt)
      deallocate(boussinesqPres)

   !--------------------------------------------------------------------

   end subroutine ocn_equation_of_state_wright_density_only!}}}

!***********************************************************************
!
!  routine ocn_equation_of_state_wright_density_exp
!
!> \brief   Computes equation of state with expansion/contraction coeffs
!> \author  Xylar Asay-Davis
!> \date    July 1, 2019
!> \details
!>  This routine computes the density from model temperature,
!>  salinity, and depth using potential-temperature formulation from
!>  Wright (1997), doi: 10.1175/1520-0426(1997)014<0735:AEOSFU>2.0.CO;2
!>  and a simple linearization of pressure p = -rho0*g*z with the Boussinesq
!>  reference density
!>
!>  Density can be computed in-situ using kDisplaced=0 and
!>      displacementType = 'relative'.
!>
!>  Potential density (referenced to zero pressure) can be computed
!>      using displacementType = 'absolute'. kDisplaced is ignored.
!>
!>  The density of SST/SSS after adiabatic displacement to each layer
!>      can be computed using displacementType = 'surfaceDisplaced'.
!>      kDisplaced is ignored and tracersSurfaceLayerValue must be present.
!
!-----------------------------------------------------------------------

   subroutine ocn_equation_of_state_wright_density_exp(nVertLevels,         &
                               nCells, kDisplaced, displacementType,        &
                               indexT, indexS, tracers, zMid, maxLevelCell, &
                               density, err, thermalExpansionCoeff,         &
                               salineContractionCoeff,                      &
                               tracersSurfaceLayerValue)
   !{{{
   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------

      integer, intent(in) :: &
         nCells,             &! number of horizontal cells
         nVertLevels,        &! max number of verical levels
         kDisplaced,         &! target layer for displacement
         indexT,             &! temperature index in tracer array
         indexS               ! salinity    index in tracer array

      character(len=*), intent(in) :: &
         displacementType     ! choice of displacement

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         tracers              ! array of tracers including T,S

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         zMid

      integer, dimension(:), intent(in) :: &
         maxLevelCell

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         tracersSurfaceLayerValue

      !-----------------------------------------------------------------
      ! Output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err  ! error flag

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         density         ! computed density

      ! Thermal expansion coeff, $-1/\rho d\rho/dT$ (note negative sign)
      ! Saline contraction coeff, $1/\rho d\rho/dS$
      real (kind=RKIND), dimension(:,:), intent(out) :: &
         thermalExpansionCoeff,  &! Thermal expansion  coeff (alpha)
         salineContractionCoeff   ! Saline contraction coeff (beta)

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------

      integer :: &
         iCell, k ! cell and vertical level loop indices

      real (kind=RKIND) :: &
         T, S, T2, T3,     &! adjusted T, S, T squared, T cubed,
         p,                &! Boussinesq pressure at reference level,
         alpha0, lambda0,  &! functions of theta and S defined by Wright (1997)
         p0,               &
         denom, denom2,    &! denominator of rational polynomial for rho
         dalpha0dT,        &! d(alpha0)/d(temperature)
         dalpha0dS,        &! d(alpha0)/d(salinity)
         dlambda0dT,       &! d(lambda0)/d(temperature)
         dlambda0dS,       &! d(lambda0)/d(salinity)
         dp0dT,            &! d(p0)/d(temperature)
         dp0dS,            &! d(p0)/d(salinity)
         drhodT,           &! derivative of density with respect to temperature
         drhodS             ! derivative of density with respect to salinity

      !-----------------------------------------------------------------

      !*** initialize error flag

      err = 0


      allocate(tracerTemp(nVertLevels, nCells))
      allocate(tracerSalt(nVertLevels, nCells))
      allocate(boussinesqPres(nVertLevels, nCells))

      if (displacementType == 'surfaceDisplaced') then
         call compute_surface_displaced_T_S(nVertLevels, nCells, indexT,      &
                                            indexS, tracersSurfaceLayerValue, &
                                            tracerTemp, tracerSalt)
      else
         call compute_bounded_T_S(nVertLevels, nCells, indexT, indexS, tracers, &
                                  tracerTemp, tracerSalt)
      end if
      call compute_boussinesq_p(nVertLevels, nCells, kDisplaced,      &
                                displacementType, zMid, maxLevelCell, &
                                boussinesqPres)


#ifdef MPAS_OPENACC
      !$acc parallel loop gang vector collapse(2)     &
      !$acc&   present(density, &
      !$acc&           thermalExpansionCoeff,         &
      !$acc&           salineContractionCoeff)
#else
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp private(k, S, T, p, T2, T3, p0, dp0dT, dp0dS, lambda0, dlambda0dT, dlambda0dS, &
      !$omp         alpha0, dalpha0dT, dalpha0dS, denom, denom2, drhodT, drhodS)
#endif


      do iCell=1,nCells
         do k=1,nVertLevels

            S = tracerSalt(k, iCell)
            T = tracerTemp(k, iCell)
            p = boussinesqPres(k, iCell)

            T2 = T*T
            T3 = T*T2

            p0 = b0 + b1*T + b2*T2 + b3*T3 + b4*S + b5*T*S
            dp0dT = b1 + 2.0_RKIND*b2*T + 3.0_RKIND*b3*T2 + b5*S
            dp0dS = b4 + b5*T

            lambda0 = c0 + c1*T + c2*T2 + c3*T3 + c4*S + c5*T*S
            dlambda0dT = c1 + 2.0_RKIND*c2*T + 3.0_RKIND*c3*T2 + c5*S
            dlambda0dS = c4 + c5*T

            alpha0 = a0 + a1*T + a2*S
            dalpha0dT = a1
            dalpha0dS = a2

            denom = 1.0_RKIND/(lambda0 + alpha0*(p + p0))
            denom2 = denom*denom

            drhodT = (lambda0*dp0dT - (p + p0)*(dlambda0dT + (p + p0)*dalpha0dT))*denom2
            drhodS = (lambda0*dp0dS - (p + p0)*(dlambda0dS + (p + p0)*dalpha0dS))*denom2

            density(k, iCell) = (p + p0)*denom
            thermalExpansionCoeff(k, iCell) = -drhodT/density(k, iCell)
            salineContractionCoeff(k, iCell) = drhodS/density(k, iCell)

         end do
      end do

#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      deallocate(tracerTemp)
      deallocate(tracerSalt)
      deallocate(boussinesqPres)

   !--------------------------------------------------------------------

   end subroutine ocn_equation_of_state_wright_density_exp!}}}

!***********************************************************************
!
!  routine ocn_equation_of_state_wright_init
!
!> \brief   Initializes Wright equation of state
!> \author  Xylar Asay-Davis
!> \date    June 30, 2019
!> \details
!>  This routine initializes a variety of quantities related to
!>  the equation of state from Wright (1997).
!
!-----------------------------------------------------------------------

   subroutine ocn_equation_of_state_wright_init(domain, err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------

      type (domain_type), intent(in) :: &
         domain        ! domain containing all state, mesh info

      !-----------------------------------------------------------------
      ! Output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err  ! error flag

      !-----------------------------------------------------------------

      !*** initialize error flag

      err = 0

      eos_wright_ref_pressure  = config_eos_wright_ref_pressure

   !--------------------------------------------------------------------

   end subroutine ocn_equation_of_state_wright_init!}}}

!***********************************************************************
!
!  routine compute_bounded_T_S
!
!> \brief   Compute T, S within valid bounds of Wright EOS
!> \author  Xylar Asay-Davis
!> \date    July 8, 2019
!
!-----------------------------------------------------------------------

   subroutine compute_bounded_T_S(nVertLevels, nCells, indexT, indexS, &
                                  tracers, tracerTemp, tracerSalt)!{{{
   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------

      integer, intent(in) :: &
         nVertLevels,        &! max number of verical levels
         nCells,             &! number of horizontal cells
         indexT,             &! temperature index in tracer array
         indexS               ! salinity    index in tracer array

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         tracers              ! array of tracers including T,S

      !-----------------------------------------------------------------
      ! Output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         tracerTemp, tracerSalt

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------


      integer :: &
         iCell, k

      real (kind=RKIND) :: &
         T, S               ! adjusted T, S

#ifdef MPAS_OPENACC
      !$acc parallel loop gang vector collapse(2) &
      !$acc    present(tracers, tracerTemp, tracerSalt)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(k, T, S)
#endif
      do iCell=1,nCells
         do k = 1, nVertLevels
            T = min(tracers(indexT, k, iCell), ocnEqStateTmax)
            S = min(tracers(indexS, k, iCell), ocnEqStateSmax)
            tracerTemp(k, iCell) = max(T, ocnEqStateTmin)
            tracerSalt(k, iCell) = max(S, ocnEqStateSmin)
         end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

   !--------------------------------------------------------------------

   end subroutine compute_bounded_T_S!}}}


!***********************************************************************
!
!  routine compute_surface_displaced_T_S
!
!> \brief   Compute surface-displaced T, S within valid bounds of Wright EOS
!> \author  Xylar Asay-Davis
!> \date    Feb 4, 2021
!
!-----------------------------------------------------------------------

   subroutine compute_surface_displaced_T_S(nVertLevels, nCells, indexT,      &
                                            indexS, tracersSurfaceLayerValue, &
                                            tracerTemp, tracerSalt)!{{{
   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------

      integer, intent(in) :: &
         nVertLevels,        &! max number of verical levels
         nCells,             &! number of horizontal cells
         indexT,             &! temperature index in tracer array
         indexS               ! salinity    index in tracer array

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         tracersSurfaceLayerValue

      !-----------------------------------------------------------------
      ! Output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         tracerTemp, tracerSalt

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------


      integer :: &
         iCell, k

      real (kind=RKIND) :: &
         T, S               ! adjusted T, S

#ifdef MPAS_OPENACC
      !$acc parallel loop gang vector collapse(2), &
      !$acc    present(tracersSurfaceLayerValue, &
      !$acc&           tracerTemp, tracerSalt)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(k, T, S)
#endif
      do iCell=1,nCells
         do k = 1, nVertLevels
            T = min(tracersSurfaceLayerValue(indexT, iCell), ocnEqStateTmax)
            S = min(tracersSurfaceLayerValue(indexS, iCell), ocnEqStateSmax)
            tracerTemp(k, iCell) = max(T, ocnEqStateTmin)
            tracerSalt(k, iCell) = max(S, ocnEqStateSmin)
         end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

   !--------------------------------------------------------------------

   end subroutine compute_surface_displaced_T_S!}}}

!***********************************************************************
!
!  routine compute_boussinesq_p
!
!> \brief   Compute Boussinesq pressure, referenced to requested layer
!> \author  Xylar Asay-Davis
!> \date    July 8, 2019
!
!-----------------------------------------------------------------------

   subroutine compute_boussinesq_p(nVertLevels, nCells, kDisplaced,      &
                                   displacementType, zMid, maxLevelCell, &
                                   boussinesqPres)
   !{{{
   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------

      integer, intent(in) :: &
         nCells,             &! number of horizontal cells
         nVertLevels,        &! max number of verical levels
         kDisplaced           ! target layer for displacement

      character(len=*), intent(in) :: &
         displacementType     ! choice of displacement

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         zMid

      integer, dimension(:), intent(in) :: &
         maxLevelCell

      !-----------------------------------------------------------------
      ! Output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         boussinesqPres

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------

      integer ::          &
         iCell, k,        &! cell and vertical level loop indices
         kPressure         ! index to determine ref level for pressure

      !-----------------------------------------------------------------
      !  if displacementType = 'relative', density is calculated
      !     with pressure referenced to level k + kDisplaced
      !     If kDisplaced=0, in-situ density is returned (no displacement)
      !     If kDisplaced/=0, a displaced density is returned
      !
      !  if displacementType = 'absolute', potential density is calculated
      !     referenced to zero pressure

      if ((displacementType == 'relative') .or. &
          (displacementType == 'surfaceDisplaced')) then

#ifdef MPAS_OPENACC
         !$acc parallel loop gang vector collapse(2) &
         !$acc    present(boussinesqPres, maxLevelCell, zMid)
#else
         !$omp parallel
         !$omp do schedule(runtime) private(k, kPressure)
#endif
         do iCell=1,nCells
            do k = 1, nVertLevels
               kPressure = min(k + kDisplaced, maxLevelCell(iCell))
               kPressure = max(kPressure, 1)
               boussinesqPres(k, iCell) = -rho_sw*gravity*zMid(kPressure, iCell)
            end do
         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif

      else ! displacementType == 'absolute'

#ifdef MPAS_OPENACC
         !$acc parallel loop collapse(2) present(boussinesqPres)
#else
         !$omp parallel
         !$omp do schedule(runtime) private(k)
#endif
         do iCell=1,nCells
            do k = 1, nVertLevels
               boussinesqPres(k, iCell) = eos_wright_ref_pressure
            end do
         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif

      endif

   !--------------------------------------------------------------------

   end subroutine compute_boussinesq_p!}}}

!***********************************************************************

end module ocn_equation_of_state_wright

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
