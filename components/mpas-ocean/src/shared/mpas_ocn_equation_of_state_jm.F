! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_equation_of_state_jm
!
!> \brief MPAS ocean Jackett and McDougall equation of state
!> \author Mark Petersen, modified by Phil Jones
!> \date   September 2011, modified May 2018
!> \details
!>  This module contains the routines for computing density from
!>  temperature and salinity using an equation of state by Jackett
!>  and McDougall that approximates the nonlinear equation of state
!>  using a set of polynomials.
!>  See Jackett and McDougall, JTECH, Vol.12, pp 381-389, April, 1995.
!
!-----------------------------------------------------------------------

module ocn_equation_of_state_jm

   use mpas_kind_types
   use mpas_pool_routines

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

   public :: ocn_equation_of_state_jm_density, &
             ocn_equation_of_state_jm_init

   !*** generic interface for case of density only or density and
   !*** expansion coeffs
   interface ocn_equation_of_state_jm_density
      module procedure ocn_equation_of_state_jm_density_only
      module procedure ocn_equation_of_state_jm_density_exp
   end interface

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   !*** depth-based reference pressure for pressure term in eq state

   real (kind=RKIND), dimension(:), allocatable :: &
      ocnEqStatePRef   ! reference pressure at each layer

   !*** valid range of T,S for Jackett and McDougall

   real (kind=RKIND), parameter :: &
      ocnEqStateTmin = -2.0_RKIND, &! valid pot. temp. range
      ocnEqStateTmax = 40.0_RKIND, &
      ocnEqStateSmin =  0.0_RKIND, &! valid salinity, in psu
      ocnEqStateSmax = 42.0_RKIND

   !***  UNESCO EOS constants and JMcD bulk modulus constants

   !*** for density of fresh water (standard UNESCO)

   real (kind=RKIND), parameter ::         &
      unt0 =   999.842594_RKIND,           &
      unt1 =  6.793952e-2_RKIND,           &
      unt2 = -9.095290e-3_RKIND,           &
      unt3 =  1.001685e-4_RKIND,           &
      unt4 = -1.120083e-6_RKIND,           &
      unt5 =  6.536332e-9_RKIND

   !*** for dependence of surface density on salinity (UNESCO)

   real (kind=RKIND), parameter ::         &
      uns1t0 =  0.824493_RKIND ,           &
      uns1t1 = -4.0899e-3_RKIND,           &
      uns1t2 =  7.6438e-5_RKIND,           &
      uns1t3 = -8.2467e-7_RKIND,           &
      uns1t4 =  5.3875e-9_RKIND,           &
      unsqt0 = -5.72466e-3_RKIND,          &
      unsqt1 =  1.0227e-4_RKIND,           &
      unsqt2 = -1.6546e-6_RKIND,           &
      uns2t0 =  4.8314e-4_RKIND

   !*** from Table A1 of Jackett and McDougall

   real (kind=RKIND), parameter ::         &
      bup0s0t0 =  1.965933e+4_RKIND,       &
      bup0s0t1 =  1.444304e+2_RKIND,       &
      bup0s0t2 = -1.706103_RKIND   ,       &
      bup0s0t3 =  9.648704e-3_RKIND,       &
      bup0s0t4 = -4.190253e-5_RKIND

   real (kind=RKIND), parameter ::         &
      bup0s1t0 =  5.284855e+1_RKIND,       &
      bup0s1t1 = -3.101089e-1_RKIND,       &
      bup0s1t2 =  6.283263e-3_RKIND,       &
      bup0s1t3 = -5.084188e-5_RKIND

   real (kind=RKIND), parameter ::         &
      bup0sqt0 =  3.886640e-1_RKIND,       &
      bup0sqt1 =  9.085835e-3_RKIND,       &
      bup0sqt2 = -4.619924e-4_RKIND

   real (kind=RKIND), parameter ::         &
      bup1s0t0 =  3.186519_RKIND   ,       &
      bup1s0t1 =  2.212276e-2_RKIND,       &
      bup1s0t2 = -2.984642e-4_RKIND,       &
      bup1s0t3 =  1.956415e-6_RKIND

   real (kind=RKIND), parameter ::         &
      bup1s1t0 =  6.704388e-3_RKIND,       &
      bup1s1t1 = -1.847318e-4_RKIND,       &
      bup1s1t2 =  2.059331e-7_RKIND,       &
      bup1sqt0 =  1.480266e-4_RKIND

   real (kind=RKIND), parameter ::         &
      bup2s0t0 =  2.102898e-4_RKIND,       &
      bup2s0t1 = -1.202016e-5_RKIND,       &
      bup2s0t2 =  1.394680e-7_RKIND,       &
      bup2s1t0 = -2.040237e-6_RKIND,       &
      bup2s1t1 =  6.128773e-8_RKIND,       &
      bup2s1t2 =  6.207323e-10_RKIND

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_equation_of_state_jm_density_only
!
!> \brief   Calls JM equation of state
!> \author  Mark Petersen and Todd Ringler, modified Phil Jones
!> \date    September 2011, updated August 2013, May 2018
!> \details
!>  This routine computes the density from model temperature and
!>  salinity using a potential-temperature-based bulk modulus from
!>  Jackett and McDougall (JTECH, Vol.12, pp 381-389, April, 1995)
!>  to approximate the UNESCO equation of state.
!>
!>  Density can be computed in-situ using kDisplaced=0 and
!>      displacementType = 'relative'.
!>
!>  Potential density (referenced to top layer) can be computed
!>      using kDisplaced=1 and displacementType = 'absolute'.
!>
!>  The density of SST/SSS after adiabatic displacement to each layer
!>      can be computed using displacementType = 'surfaceDisplaced'.
!>
!>  When using displacementType = 'surfaceDisplaced', kDisplaced is
!>      ignored and tracersSurfaceLayerValue must be present.
!
!-----------------------------------------------------------------------

   subroutine ocn_equation_of_state_jm_density_only(nVertLevels,      &
                               nCells, kDisplaced, displacementType,  &
                               indexT, indexS, tracers, density, err, &
                               tracersSurfaceLayerValue)
   !{{{
   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------

      integer, intent(in) :: &
         nCells,             &! number of horizontal cells
         nVertLevels,        &! max number of vertical levels
         kDisplaced,         &! target layer for displacement
         indexT,             &! temperature index in tracer array
         indexS               ! salinity    index in tracer array

      character(len=*), intent(in) :: &
         displacementType     ! choice of displacement

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         tracers              ! array of tracers including T,S

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

      integer ::          &
         iCell, k,        &! cell and vertical level loop indices
         kTmp, kRef,      &! indices to determine ref level for pressure
         kDisplacedLocal   ! locally modified displacement target level

      character(len=60) :: &
         displacementTypeLocal ! locally modified displacement type

      real (kind=RKIND), dimension(:), allocatable :: &
         p, p2 ! temporary pressure scalars

      real (kind=RKIND), dimension(:,:), allocatable :: &
         tracerTemp,          &! modified temperature
         tracerSalt            ! modified salinity
#ifdef MPAS_OPENACC
      !$acc declare device_resident(tracerTemp, tracerSalt)
#endif

      real (kind=RKIND) :: &
         work1, work2,     &! temporary scalars for calc
         work3, work4,     &!
         tq, sq,           &! adjusted T,S
         t2, sqr,          &! temperature squared and square root salt
         bulkMod,          &! Bulk modulus
         rhosfc             ! density at the surface

      !-----------------------------------------------------------------

      !*** initialize error flag

      err = 0

      !*** allocate and compute pressure array

      allocate(p(nVertLevels),p2(nVertLevels))
#ifdef MPAS_OPENACC
      !$acc enter data create(p, p2)
#endif

      ! Determine pressure to use in density calculation
      !  If kDisplaced=0, in-situ density is returned (no displacement)
      !  If kDisplaced/=0, potential density is returned

      !  if displacementType = 'relative', potential density is calculated
      !     referenced to level k + kDisplaced
      !  if displacementType = 'absolute', potential density is calculated
      !     referenced to level kDisplaced for all k
      !  NOTE: kDisplaced = 0 or > nVertLevels is incompatible with 'absolute'
      !     so abort if necessary

      if (displacementType == 'surfaceDisplaced') then
         displacementTypeLocal = 'relative'
         kDisplacedLocal = 0
      else
         displacementTypeLocal = trim(displacementType)
         kDisplacedLocal = kDisplaced
      endif

      ! Eliminate this error check for performance reasons
      !if (displacementTypeLocal == 'absolute' .and.   &
      !   (kDisplacedLocal <= 0 .or. kDisplacedLocal > nVertLevels) ) then
      !
      !   call mpas_log_write('Abort: In equation_of_state_jm' // &
      !       ' kDisplaced must be between 1 and nVertLevels for ' // &
      !       'displacementType = absolute', MPAS_LOG_CRIT)
      !endif

      if (kDisplacedLocal == 0) then
         ! use pressure at in situ level
#ifdef MPAS_OPENACC
         !$acc parallel loop gang vector &
         !$acc    present(p, p2, ocnEqStatePRef)
#endif
         do k=1,nVertLevels
            p (k) = ocnEqStatePRef(k)
            p2(k) = ocnEqStatePRef(k)*ocnEqStatePRef(k)
         enddo
      else ! kDisplacedLocal /= 0
#ifdef MPAS_OPENACC
         !$acc parallel loop gang vector &
         !$acc    present(p, p2, ocnEqStatePRef)
#endif
         do k=1,nVertLevels
            ! determine level to which parcel is displaced
            !  and thus which pressure (depth) to use in EOS
            if (displacementTypeLocal == 'relative') then
               kTmp = min(k + kDisplacedLocal, nVertLevels)
            else
               kTmp = min(kDisplacedLocal, nVertLevels)
            endif
            kRef  = max(kTmp, 1)  ! make sure index bounded to 1
            p (k) = ocnEqStatePRef(kRef)
            p2(k) = ocnEqStatePRef(kRef)*ocnEqStatePRef(kRef)
         enddo
      endif

      !*** compute modified T,S to account for displacement and
      !*** valid range

      ! For OpenACC, these are device resident only
      allocate(tracerTemp(nVertLevels,nCells), &
               tracerSalt(nVertLevels,nCells))

      if (displacementType == 'surfaceDisplaced') then

#ifdef MPAS_OPENACC
         !$acc parallel loop gang vector collapse(2) &
         !$acc&  present(tracerssurfacelayervalue)
#else
         !$omp parallel
         !$omp do schedule(runtime) private(k, tq, sq)
#endif
         do iCell=1,nCells
         do k=1,nVertLevels
            tq = min(tracersSurfaceLayerValue(indexT,iCell), &
                     ocnEqStateTmax)
            sq = min(tracersSurfaceLayerValue(indexS,iCell), &
                     ocnEqStateSmax)
            tracerTemp(k,iCell) = max(tq,ocnEqStateTmin)
            tracerSalt(k,iCell) = max(sq,ocnEqStateSmin)
         enddo
         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif

      else

#ifdef MPAS_OPENACC
         !$acc parallel loop gang vector collapse(2) &
         !$acc&  present(tracers)
#else
         !$omp parallel
         !$omp do schedule(runtime) private(k, tq, sq)
#endif
         do iCell=1,nCells
         do k = 1, nVertLevels
            tq = min(tracers(indexT,k,iCell), ocnEqStateTmax)
            sq = min(tracers(indexS,k,iCell), ocnEqStateSmax)
            tracerTemp(k,iCell) = max(tq,ocnEqStateTmin)
            tracerSalt(k,iCell) = max(sq,ocnEqStateSmin)
         end do
         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif

      endif

#ifdef MPAS_OPENACC
      !$acc parallel loop gang vector collapse(2) &
      !$acc&   present(density, p, p2)
#else
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp private(k, sq, tq, sqr, t2, work1, work2, rhosfc, work3, work4, bulkMod)
#endif
      do iCell=1,nCells
      do k=1,nVertLevels

         sq  = tracerSalt(k,iCell)
         tq  = tracerTemp(k,iCell)

         sqr = sqrt(sq)
         t2  = tq*tq

         !***
         !*** first calculate surface (p=0) values from UNESCO eqns.
         !***

         work1 =      uns1t0 + uns1t1*tq + &
                     (uns1t2 + uns1t3*tq + uns1t4*t2)*t2
         work2 = sqr*(unsqt0 + unsqt1*tq + unsqt2*t2)

         rhosfc = unt1*tq + (unt2 + unt3*tq + (unt4 + unt5*tq)*t2)*t2 &
                          + (uns2t0*sq + work1 + work2)*sq

         !***
         !*** now calculate bulk modulus at pressure p from
         !*** Jackett and McDougall formula
         !***

         work3 = bup0s1t0 + bup0s1t1*tq +     &
                (bup0s1t2 + bup0s1t3*tq)*t2 + &
          p(k) *(bup1s1t0 + bup1s1t1*tq + bup1s1t2*t2) + &
          p2(k)*(bup2s1t0 + bup2s1t1*tq + bup2s1t2*t2)

         work4 = sqr*(bup0sqt0 + bup0sqt1*tq + bup0sqt2*t2 + &
                      bup1sqt0*p(k))

         bulkMod = bup0s0t0 + bup0s0t1*tq +                    &
                  (bup0s0t2 + bup0s0t3*tq + bup0s0t4*t2)*t2 + &
            p(k) *(bup1s0t0 + bup1s0t1*tq +                &
                  (bup1s0t2 + bup1s0t3*tq)*t2) +           &
            p2(k)*(bup2s0t0 + bup2s0t1*tq + bup2s0t2*t2) + &
                   sq*(work3 + work4)


         density(k,iCell) = (unt0 + rhosfc)*bulkMod/ &
                            (bulkMod - p(k))

      end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

#ifdef MPAS_OPENACC
      !$acc exit data delete(p, p2)
#endif
      deallocate(p,p2)
      deallocate(tracerTemp, tracerSalt)

   !--------------------------------------------------------------------

   end subroutine ocn_equation_of_state_jm_density_only!}}}

!***********************************************************************
!
!  routine ocn_equation_of_state_jm_density
!
!> \brief   Calls JM equation of state
!> \author  Mark Petersen and Todd Ringler, modified Phil Jones
!> \date    September 2011, updated August 2013, May 2018
!> \details
!>  This routine computes the density from model temperature and
!>  salinity using a potential-temperature-based bulk modulus from
!>  Jackett and McDougall (JTECH, Vol.12, pp 381-389, April, 1995)
!>  to approximate the UNESCO equation of state.
!>
!>  Density can be computed in-situ using kDisplaced=0 and
!>      displacementType = 'relative'.
!>
!>  Potential density (referenced to top layer) can be computed
!>      using kDisplaced=1 and displacementType = 'absolute'.
!>
!>  The density of SST/SSS after adiabatic displacement to each layer
!>      can be computed using displacementType = 'surfaceDisplaced'.
!>
!>  When using displacementType = 'surfaceDisplaced', kDisplaced is
!>      ignored and tracersSurfaceLayerValue must be present.
!
!-----------------------------------------------------------------------

   subroutine ocn_equation_of_state_jm_density_exp(nVertLevels,       &
                               nCells, kDisplaced, displacementType,  &
                               indexT, indexS, tracers, density, err, &
                               thermalExpansionCoeff,                 &
                               salineContractionCoeff,                &
                               tracersSurfaceLayerValue)
   !{{{
   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------

      integer, intent(in) :: &
         nCells,             &! number of horizontal cells
         nVertLevels,        &! max number of vertical levels
         kDisplaced,         &! target layer for displacement
         indexT,             &! temperature index in tracer array
         indexS               ! salinity    index in tracer array

      character(len=*), intent(in) :: &
         displacementType     ! choice of displacement

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         tracers              ! array of tracers including T,S

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

      integer ::          &
         iCell, k,        &! cell and vertical level loop indices
         kTmp, kRef,      &! indices to determine ref level for pressure
         kDisplacedLocal   ! locally modified displacement target level

      character(len=60) :: &
         displacementTypeLocal ! locally modified displacement type

      real (kind=RKIND) :: &
         tq, sq,           &! adjusted T,S
         t2, sqr,          &! T squared and square root salt
         rhosfc,           &! density at the surface
         bulkMod,          &! Bulk modulus
         denomk,           &! temp for avoiding division
         work1, work2,     &! temporary work space
         work3, work4,     &!
         drdt0,            &! d(density)/d(temperature), for surface
         drds0,            &! d(density)/d(salinity   ), for surface
         dkdt,             &! d(bulk modulus)/d(pot. temp.)
         dkds,             &! d(bulk modulus)/d(salinity  )
         drhodt,           &! derivative of density with respect to temperature
         drhods             ! derivative of density with respect to salinity

      real (kind=RKIND), dimension(:), allocatable :: &
         p, p2 ! temporary pressure scalars

      real (kind=RKIND), dimension(:,:), allocatable :: &
         tracerTemp,          &! modified temperature
         tracerSalt            ! modified salinity
#ifdef MPAS_OPENACC
      !$acc declare device_resident(tracerTemp, tracerSalt)
#endif

      !-----------------------------------------------------------------

      !*** initialize error flag

      err = 0

      !*** allocate and compute pressure array

      allocate(p(nVertLevels),p2(nVertLevels))
#ifdef MPAS_OPENACC
      !$acc enter data create(p, p2)
#endif

      ! Determine pressure to use in density calculation
      !  If kDisplaced=0, in-situ density is returned (no displacement)
      !  If kDisplaced/=0, potential density is returned

      !  if displacementType = 'relative', potential density is calculated
      !     referenced to level k + kDisplaced
      !  if displacementType = 'absolute', potential density is calculated
      !     referenced to level kDisplaced for all k
      !  NOTE: kDisplaced = 0 or > nVertLevels is incompatible with 'absolute'
      !     so abort if necessary

      if (displacementType == 'surfaceDisplaced') then
         displacementTypeLocal = 'relative'
         kDisplacedLocal = 0
      else
         displacementTypeLocal = trim(displacementType)
         kDisplacedLocal = kDisplaced
      endif

      ! Eliminate this error check for performance reasons
      !if (displacementTypeLocal == 'absolute' .and.   &
      !   (kDisplacedLocal <= 0 .or. kDisplacedLocal > nVertLevels) ) then
      !
      !   call mpas_log_write('Abort: In equation_of_state_jm' // &
      !       ' kDisplaced must be between 1 and nVertLevels for ' // &
      !       'displacementType = absolute', MPAS_LOG_CRIT)
      !endif

      if (kDisplacedLocal == 0) then
         ! use pressure at in situ level
#ifdef MPAS_OPENACC
         !$acc parallel loop gang vector &
         !$acc    present(p, p2, ocnEqStatePRef)
#endif
         do k=1,nVertLevels
            p (k) = ocnEqStatePRef(k)
            p2(k) = ocnEqStatePRef(k)*ocnEqStatePRef(k)
         enddo
      else ! kDisplacedLocal /= 0
#ifdef MPAS_OPENACC
         !$acc parallel loop gang vector &
         !$acc    present(p, p2, ocnEqStatePRef)
#endif
         do k=1,nVertLevels
            ! determine level to which parcel is displaced
            !  and thus which pressure (depth) to use in EOS
            if (displacementTypeLocal == 'relative') then
               kTmp = min(k + kDisplacedLocal, nVertLevels)
            else
               kTmp = min(kDisplacedLocal, nVertLevels)
            endif
            kRef  = max(kTmp, 1)  ! make sure index bounded to 1
            p (k) = ocnEqStatePRef(kRef)
            p2(k) = ocnEqStatePRef(kRef)*ocnEqStatePRef(kRef)
         enddo
      endif

      !*** compute modified T,S to account for displacement and
      !*** valid range

      ! For OpenACC, these are device resident only
      allocate(tracerTemp(nVertLevels,nCells), &
               tracerSalt(nVertLevels,nCells))

      if (displacementType == 'surfaceDisplaced') then

#ifdef MPAS_OPENACC
         !$acc parallel loop gang vector collapse(2) &
         !$acc&  present(tracerssurfacelayervalue)
#else
         !$omp parallel
         !$omp do schedule(runtime) private(k, tq, sq)
#endif
         do iCell=1,nCells
         do k=1,nVertLevels
            tq = min(tracersSurfaceLayerValue(indexT,iCell), &
                     ocnEqStateTmax)
            sq = min(tracersSurfaceLayerValue(indexS,iCell), &
                     ocnEqStateSmax)
            tracerTemp(k,iCell) = max(tq,ocnEqStateTmin)
            tracerSalt(k,iCell) = max(sq,ocnEqStateSmin)
         enddo
         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif

      else

#ifdef MPAS_OPENACC
         !$acc parallel loop gang vector collapse(2) &
         !$acc&  present(tracers)
#else
         !$omp parallel
         !$omp do schedule(runtime) private(k, tq, sq)
#endif
         do iCell=1,nCells
         do k = 1, nVertLevels
            tq = min(tracers(indexT,k,iCell), ocnEqStateTmax)
            sq = min(tracers(indexS,k,iCell), ocnEqStateSmax)
            tracerTemp(k,iCell) = max(tq,ocnEqStateTmin)
            tracerSalt(k,iCell) = max(sq,ocnEqStateSmin)
         end do
         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif

      endif

#ifdef MPAS_OPENACC
      !$acc parallel loop gang vector collapse(2)     &
      !$acc&   present(density, p, p2,                &
      !$acc&           thermalExpansionCoeff,         &
      !$acc&           salineContractionCoeff)
#else
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp private(k, sq, tq, sqr, t2, work1, work2, rhosfc, work3, work4, bulkMod, &
      !$omp         denomk, drdt0, dkdt, drhodt, drds0, dkds, drhods)
#endif
      do iCell=1,nCells
      do k=1,nVertLevels

         sq  = tracerSalt(k,iCell)
         tq  = tracerTemp(k,iCell)

         sqr = sqrt(sq)
         t2  = tq*tq

         !***
         !*** first calculate surface (p=0) values from UNESCO eqns.
         !***

         work1 = uns1t0 + uns1t1*tq + &
                (uns1t2 + uns1t3*tq + uns1t4*t2)*t2
         work2 = sqr*(unsqt0 + unsqt1*tq + unsqt2*t2)

         rhosfc = unt1*tq + (unt2 + unt3*tq + (unt4 + unt5*tq)*t2)*t2 &
                          + (uns2t0*sq + work1 + work2)*sq

         !***
         !*** now calculate bulk modulus at pressure p from
         !*** Jackett and McDougall formula
         !***

         work3 = bup0s1t0 + bup0s1t1*tq +                    &
                (bup0s1t2 + bup0s1t3*tq)*t2 +                &
          p(k) *(bup1s1t0 + bup1s1t1*tq + bup1s1t2*t2) + &
          p2(k)*(bup2s1t0 + bup2s1t1*tq + bup2s1t2*t2)

         work4 = sqr*(bup0sqt0 + bup0sqt1*tq + bup0sqt2*t2 + &
                                 bup1sqt0*p(k))

         bulkMod  = bup0s0t0 + bup0s0t1*tq +                    &
                   (bup0s0t2 + bup0s0t3*tq + bup0s0t4*t2)*t2 + &
             p(k) *(bup1s0t0 + bup1s0t1*tq +                &
                   (bup1s0t2 + bup1s0t3*tq)*t2) +           &
             p2(k)*(bup2s0t0 + bup2s0t1*tq + bup2s0t2*t2) + &
                                        sq*(work3 + work4)

         !***
         !*** compute density
         !***

         denomk = 1.0_RKIND/(bulkMod - p(k))

         density(k, iCell) = (unt0 + rhosfc)*bulkMod*denomk

         !***
         !*** compute temperature expansion coeff
         !***  by differentiating above formulae
         !***

         drdt0 =             unt1 + 2.0_RKIND*unt2*tq +        &
                  (3.0_RKIND*unt3 + 4.0_RKIND*unt4*tq +        &
                                    5.0_RKIND*unt5*t2)*t2 +    &
                          (uns1t1 + 2.0_RKIND*uns1t2*tq +      &
                (3.0_RKIND*uns1t3 + 4.0_RKIND*uns1t4*tq)*t2 +  &
                          (unsqt1 + 2.0_RKIND*unsqt2*tq)*sqr )*sq

         dkdt  =            bup0s0t1 + 2.0_RKIND*bup0s0t2*tq +      &
                 (3.0_RKIND*bup0s0t3 + 4.0_RKIND*bup0s0t4*tq)*t2 +  &
                     p(k) *(bup1s0t1 + 2.0_RKIND*bup1s0t2*tq +      &
                                       3.0_RKIND*bup1s0t3*t2) +     &
                     p2(k)*(bup2s0t1 + 2.0_RKIND*bup2s0t2*tq) +     &
                        sq*(bup0s1t1 + 2.0_RKIND*bup0s1t2*tq +      &
                                       3.0_RKIND*bup0s1t3*t2 +      &
                    p(k)  *(bup1s1t1 + 2.0_RKIND*bup1s1t2*tq) +     &
                    p2(k) *(bup2s1t1 + 2.0_RKIND*bup2s1t2*tq) +     &
                       sqr*(bup0sqt1 + 2.0_RKIND*bup0sqt2*tq))

         drhodt = (denomk*(drdt0*bulkMod -                    &
                  p(k)*(unt0+rhosfc)*dkdt*denomk))

         thermalExpansionCoeff(k,iCell) = -drhodt/density(k,iCell)

         !***
         !*** compute salinity contraction coeff
         !***  by differentiating above formulae
         !***

         drds0  = 2.0_RKIND*uns2t0*sq + work1 + 1.5_RKIND*work2
         dkds   = work3 + 1.5_RKIND*work4

         drhods = denomk*(drds0*bulkMod -                    &
                   p(k)*(unt0+rhosfc)*dkds*denomk)

         salineContractionCoeff(k,iCell) = drhods/density(k,iCell)

      end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

#ifdef MPAS_OPENACC
      !$acc exit data delete(p, p2)
#endif
      deallocate(p,p2)
      deallocate(tracerTemp, tracerSalt)

   !--------------------------------------------------------------------

   end subroutine ocn_equation_of_state_jm_density_exp!}}}

!***********************************************************************
!
!  routine ocn_equation_of_state_jm_init
!
!> \brief   Initializes ocean momentum horizontal mixing quantities
!> \author  Mark Petersen
!> \date    September 2011
!> \details
!>  This routine initializes a variety of quantities related to
!>  horizontal velocity mixing in the ocean. Since a variety of
!>  parameterizations are available, this routine primarily calls the
!>  individual init routines for each parameterization.
!
!-----------------------------------------------------------------------

   subroutine ocn_equation_of_state_jm_init(domain, err)!{{{

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
      ! Local variables
      !-----------------------------------------------------------------

      integer ::            &
         k,                 &! depth loop index
         nVertLevels         ! max number of vertical levels

      real (RKIND) :: depth  ! depth of layer midpoint

      integer, pointer ::   &
         nVertLevelsPtr      ! pointer for extracting max vert levels

      real (kind=RKIND), dimension(:), pointer :: &
         refBottomDepth      ! depth at bottom of each layer for
                             !  a reference depth profile

      type (block_type), pointer :: block ! pool info for each block

      type (mpas_pool_type), pointer :: meshPool

      !-----------------------------------------------------------------

      !*** initialize error flag

      err = 0

      !*** extract vertical mesh info for computing ref pressure
      !*** since vertical mesh info same for all blocks, only need first

      block => domain % blocklist
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

      call mpas_pool_get_array(meshPool,     'refBottomDepth', &
                                              refBottomDepth)

      call mpas_pool_get_dimension(meshPool, 'nVertLevels', &
                                              nVertLevelsPtr)
      nVertLevels = nVertLevelsPtr

      !*** compute reference pressure in bars from depth in meters
      !***  using a mean density derived from depth-dependent global
      !***  average temperatures and salinities from Levitus 1994, and
      !***  integrating using hydrostatic balance.
      !*** Note this calculation uses refBottomDepth, so pressure on
      !***  top level does not include SSH contribution.
      !*** Also note that this is only valid for case in which the
      !***  vertical mesh does not vary significantly in time. Future
      !***  Lagrangian meshes will likely need to use full
      !***  3-d dynamic depth/press

      allocate(ocnEqStatePRef(nVertLevels))

      depth = 0.5_RKIND*refBottomDepth(1)
      ocnEqStatePRef(1) = &
           0.059808_RKIND  *(exp(-0.025_RKIND*depth) - 1.0_RKIND) &
         + 0.100766_RKIND  *depth &
         + 2.28405e-7_RKIND*depth**2
      do k = 2,nVertLevels
         depth = 0.5_RKIND*(refBottomDepth(k) + refBottomDepth(k-1))
         ocnEqStatePRef(k) = &
              0.059808_RKIND  *(exp(-0.025_RKIND*depth) - 1.0_RKIND) &
            + 0.100766_RKIND  *depth                                 &
            + 2.28405e-7_RKIND*depth**2
      enddo
#ifdef MPAS_OPENACC
      !$acc enter data copyin(ocnEqStatePRef)
#endif

   !--------------------------------------------------------------------

   end subroutine ocn_equation_of_state_jm_init!}}}

!***********************************************************************

end module ocn_equation_of_state_jm

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
