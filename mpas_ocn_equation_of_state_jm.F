! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
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
!> \brief MPAS ocean equation of state driver
!> \author Mark Petersen
!> \date   September 2011
!> \details
!>  This module contains the main driver routine for calling
!>  the equation of state.
!
!-----------------------------------------------------------------------

module ocn_equation_of_state_jm

   use mpas_kind_types
   use mpas_grid_types
   use mpas_configure

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

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_equation_of_state_jm_density
!
!> \brief   Calls JM equation of state
!> \author  Mark Petersen
!> \date    September 2011
!> \details 
!>  This routine uses a JM equation of state to update the density
!
!-----------------------------------------------------------------------

   subroutine ocn_equation_of_state_jm_density(grid, k_displaced, displacement_type, indexT, indexS, tracers, density, err, tracersSurfaceValue)!{{{
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !  This module contains routines necessary for computing the density
   !  from model temperature and salinity using an equation of state.
   !
   !  The UNESCO equation of state computed using the
   !  potential-temperature-based bulk modulus from Jackett and
   !  McDougall, JTECH, Vol.12, pp 381-389, April, 1995.
   !
   ! Input: grid - grid metadata
   !        s - state: tracers
   !        k_displaced 

   !  If k_displaced<=0, density is returned with no displaced
   !  If k_displaced>0,the density returned is that for a parcel
   !  adiabatically displaced from its original level to level 
   !  k_displaced.

   !
   ! Output: s - state: computed density
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (mesh_type), intent(in) :: grid
      integer, intent(in) :: k_displaced, indexT, indexS
      character(len=*), intent(in) :: displacement_type
      real (kind=RKIND), dimension(:,:,:), intent(in) :: tracers
      real (kind=RKIND), dimension(:,:), intent(out) :: density
      integer, intent(out) :: err
      real (kind=RKIND), dimension(:,:), intent(in), optional :: tracersSurfaceValue

      type (dm_info) :: dminfo
      integer :: iEdge, iCell, iVertex, k

      integer :: nCells, nEdges, nVertices, nVertLevels
      integer :: k_displaced_local
      character(len=60) :: displacement_type_local

      real (kind=RKIND), dimension(:), pointer :: &
        refBottomDepth, pRefEOS

      integer, dimension(:), pointer :: maxLevelCell

      real (kind=RKIND) :: &
         TQ,SQ,             &! adjusted T,S
         BULK_MOD,          &! Bulk modulus
         density_S,             &! density at the surface
         DRDT0,             &! d(density)/d(temperature), for surface
         DRDS0,             &! d(density)/d(salinity   ), for surface
         DKDT,              &! d(bulk modulus)/d(pot. temp.)
         DKDS,              &! d(bulk modulus)/d(salinity  )
         SQR,DENOMK,        &! work arrays
         WORK1, WORK2, WORK3, WORK4, T2, depth

      real (kind=RKIND) :: & 
         tmin, tmax,        &! valid temperature range for level k
         smin, smax          ! valid salinity    range for level k

      real (kind=RKIND), dimension(:), allocatable :: &
         p, p2 ! temporary pressure scalars

      real (kind=RKIND), dimension(:,:,:), allocatable :: &
         tracerTS

!-----------------------------------------------------------------------
!
!  UNESCO EOS constants and JMcD bulk modulus constants
!
!-----------------------------------------------------------------------

      !*** for density of fresh water (standard UNESCO)

      real (kind=RKIND), parameter ::              &
         unt0 =   999.842594,           &
         unt1 =  6.793952e-2,           &
         unt2 = -9.095290e-3,           &
         unt3 =  1.001685e-4,           &
         unt4 = -1.120083e-6,           &
         unt5 =  6.536332e-9
   
      !*** for dependence of surface density on salinity (UNESCO)

      real (kind=RKIND), parameter ::              &
         uns1t0 =  0.824493 ,           &
         uns1t1 = -4.0899e-3,           &
         uns1t2 =  7.6438e-5,           &
         uns1t3 = -8.2467e-7,           &
         uns1t4 =  5.3875e-9,           &
         unsqt0 = -5.72466e-3,          &
         unsqt1 =  1.0227e-4,           &
         unsqt2 = -1.6546e-6,           &
         uns2t0 =  4.8314e-4
   
      !*** from Table A1 of Jackett and McDougall
   
      real (kind=RKIND), parameter ::              &
         bup0s0t0 =  1.965933e+4,       &
         bup0s0t1 =  1.444304e+2,       &
         bup0s0t2 = -1.706103   ,       &
         bup0s0t3 =  9.648704e-3,       &
         bup0s0t4 = -4.190253e-5
   
      real (kind=RKIND), parameter ::              &
         bup0s1t0 =  5.284855e+1,       &
         bup0s1t1 = -3.101089e-1,       &
         bup0s1t2 =  6.283263e-3,       &
         bup0s1t3 = -5.084188e-5
   
      real (kind=RKIND), parameter ::              &
         bup0sqt0 =  3.886640e-1,       &
         bup0sqt1 =  9.085835e-3,       &
         bup0sqt2 = -4.619924e-4
   
      real (kind=RKIND), parameter ::              &
         bup1s0t0 =  3.186519   ,       &
         bup1s0t1 =  2.212276e-2,       &
         bup1s0t2 = -2.984642e-4,       &
         bup1s0t3 =  1.956415e-6 
   
      real (kind=RKIND), parameter ::              &
         bup1s1t0 =  6.704388e-3,       &
         bup1s1t1 = -1.847318e-4,       &
         bup1s1t2 =  2.059331e-7,       &
         bup1sqt0 =  1.480266e-4 
   
      real (kind=RKIND), parameter ::              &
         bup2s0t0 =  2.102898e-4,       &
         bup2s0t1 = -1.202016e-5,       &
         bup2s0t2 =  1.394680e-7,       &
         bup2s1t0 = -2.040237e-6,       &
         bup2s1t1 =  6.128773e-8,       &
         bup2s1t2 =  6.207323e-10
   
      integer :: k_test, k_ref
   
      err = 0
   
      nCells      = grid % nCells
      maxLevelCell      => grid % maxLevelCell % array
      nVertLevels = grid % nVertLevels
      refBottomDepth => grid % refBottomDepth % array

!  allocate local T,S tracer field
      allocate(tracerTS(2,nVertLevels,nCells))

!  fill tracerTS
      if (displacement_type == 'surfaceDisplaced') then
        if(present(tracersSurfaceValue)) then
          do k=1,nVertLevels
            tracerTS(1,k,:) = tracersSurfaceValue(indexT,:)
            tracerTS(2,k,:) = tracersSurfaceValue(indexS,:)
          enddo
          displacement_type_local = 'relative'
          k_displaced_local = 0
        else
          write(0,*) 'Abort: tracersSurfaceValue must be present'
          call mpas_dmpar_abort(dminfo)
        endif
      else
        tracerTS(1,:,:) = tracers(indexT,:,:)
        tracerTS(2,:,:) = tracers(indexS,:,:)
        displacement_type_local = trim(displacement_type)
        k_displaced_local = k_displaced
      endif

!  Jackett and McDougall
      tmin = -2.0  ! valid pot. temp. range
      tmax = 40.0 
      smin =  0.0  ! valid salinity, in psu   
      smax = 42.0 

!  This function computes pressure in bars from depth in meters
!  using a mean density derived from depth-dependent global 
!  average temperatures and salinities from Levitus 1994, and 
!  integrating using hydrostatic balance.

      allocate(pRefEOS(nVertLevels),p(nVertLevels),p2(nVertLevels))

      ! This could be put in the init routine.
      ! Note I am using refBottomDepth, so pressure on top level does
      ! not include SSH contribution.  I am not sure if that matters, but
      ! POP does it the same way.
      depth = 0.5*refBottomDepth(1)
      pRefEOS(1) = 0.059808*(exp(-0.025*depth) - 1.0) &
          + 0.100766*depth + 2.28405e-7*depth**2
      do k = 2,nVertLevels
         depth = 0.5*(refBottomDepth(k)+refBottomDepth(k-1))
         pRefEOS(k) = 0.059808*(exp(-0.025*depth) - 1.0) &
             + 0.100766*depth + 2.28405e-7*depth**2
      enddo

      !  If k_displaced=0, in-situ density is returned (no displacement)
      !  If k_displaced/=0, potential density is returned

      !  if displacement_type = 'relative', potential density is calculated
      !     referenced to level k + k_displaced
      !  if displacement_type = 'absolute', potential density is calculated
      !     referenced to level k_displaced for all k
      !  NOTE: k_displaced = 0 or > nVertLevels is incompatible with 'absolute'
      !     so abort if necessary

      if (displacement_type_local == 'absolute' .and.   &
         (k_displaced_local <= 0 .or. k_displaced_local > nVertLevels) ) then

         write(0,*) 'Abort: In equation_of_state_jm', &
             ' k_displaced must be between 1 and nVertLevels for ', &
             'displacement_type = absolute'
         call mpas_dmpar_abort(dminfo)
      endif

      if (k_displaced_local == 0) then
         do k=1,nVertLevels
            p(k)   = pRefEOS(k)
            p2(k)  = p(k)*p(k)
         enddo
      else ! k_displaced_local /= 0
         do k=1,nVertLevels
            if (displacement_type_local == 'relative') then
               k_test = min(k + k_displaced_local, nVertLevels)
               k_ref  = max(k_test, 1)
            else
               k_test = min(k_displaced_local, nVertLevels)
               k_ref  = max(k_test, 1)
            endif
            p(k)   = pRefEOS(k_ref)
            p2(k)  = p(k)*p(k)
         enddo
      endif

      do iCell=1,nCells
         do k=1,maxLevelCell(iCell)
            SQ  = max(min(tracerTS(2,k,iCell),smax),smin)
            TQ  = max(min(tracerTS(1,k,iCell),tmax),tmin)
   
            SQR = sqrt(SQ)
            T2  = TQ*TQ

            !***
            !*** first calculate surface (p=0) values from UNESCO eqns.
            !***

            WORK1 = uns1t0 + uns1t1*TQ + & 
                   (uns1t2 + uns1t3*TQ + uns1t4*T2)*T2
            WORK2 = SQR*(unsqt0 + unsqt1*TQ + unsqt2*T2)

            density_S = unt1*TQ + (unt2 + unt3*TQ + (unt4 + unt5*TQ)*T2)*T2 &
                            + (uns2t0*SQ + WORK1 + WORK2)*SQ

            !***
            !*** now calculate bulk modulus at pressure p from 
            !*** Jackett and McDougall formula
            !***

            WORK3 = bup0s1t0 + bup0s1t1*TQ +                    &
                    (bup0s1t2 + bup0s1t3*TQ)*T2 +                &
                    p(k) *(bup1s1t0 + bup1s1t1*TQ + bup1s1t2*T2) + &
                    p2(k)*(bup2s1t0 + bup2s1t1*TQ + bup2s1t2*T2)
            WORK4 = SQR*(bup0sqt0 + bup0sqt1*TQ + bup0sqt2*T2 + &
                         bup1sqt0*p(k))
   
            BULK_MOD  = bup0s0t0 + bup0s0t1*TQ +                    &
                        (bup0s0t2 + bup0s0t3*TQ + bup0s0t4*T2)*T2 + &
                        p(k) *(bup1s0t0 + bup1s0t1*TQ +                &
                        (bup1s0t2 + bup1s0t3*TQ)*T2) +           &
                        p2(k)*(bup2s0t0 + bup2s0t1*TQ + bup2s0t2*T2) + &
                        SQ*(WORK3 + WORK4)
  
            DENOMK = 1.0/(BULK_MOD - p(k))
   
            density(k,iCell) = (unt0 + density_S)*BULK_MOD*DENOMK

         end do
      end do

      deallocate(pRefEOS,p,p2)
      deallocate(tracerTS)

   end subroutine ocn_equation_of_state_jm_density!}}}

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

   subroutine ocn_equation_of_state_jm_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err

      err = 0

   !--------------------------------------------------------------------

   end subroutine ocn_equation_of_state_jm_init!}}}

!***********************************************************************

end module ocn_equation_of_state_jm

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
