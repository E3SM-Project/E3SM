! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_tidal_potential
!
!> \brief MPAS ocean tidal potential forcing module
!> \author Steven Brus 
!> \date   September 2019
!> \details
!>  This module contains routines for the tidal potential forcing.
!>  Design document located in :
!>    MPAS-Model/docs/ocean/design_docs
!
!-----------------------------------------------------------------------

module ocn_vel_tidal_potential

   use mpas_kind_types
   use mpas_constants
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timekeeping
   use mpas_timer
   use ocn_constants
   use ocn_config
   use ocn_mesh
   use ocn_diagnostics_variables

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

   public :: ocn_vel_tidal_potential_tend, &
             ocn_compute_tidal_potential_forcing, &
             ocn_vel_tidal_potential_init, &
             tidal_constituent_factors

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: &
      tidalPotentialOff, &! on/off switch for tidal potential
      timeIntegratorLTS   ! flag for Local Time Stepping

   !*** physical constants
   real (kind=RKIND) ::  &
      betaSelfAttrLoad,  &! self-attraction and loading beta
      tidalPotRamp,      &! scale for ramping tidal potential
      rho0gInv            ! 1 / rho0 / gravity

   !*** eventually these will be private module arrays and 
   !*** not in the forcing pool so retain pointers here as placeholders
   real (kind=RKIND), dimension(:), pointer :: &
      tidalPotEta,                    &!
      tidalConstituentAmplitude,      &!
      tidalConstituentFrequency,      &!
      tidalConstituentLoveNumbers,    &!
      tidalConstituentNodalAmplitude, &!
      tidalConstituentAstronomical,   &!
      tidalConstituentNodalPhase       !

   real (kind=RKIND), dimension(:,:), pointer :: &
      latitudeFunction 

   integer, dimension(:), pointer :: &
      tidalConstituentType

   integer :: &
      nTidalConstituents    ! number of tidal constituents

   type :: char_array
     character(:), allocatable :: constituent
   end type
   type(char_array), dimension(37) :: constituentList
   public :: char_array

   ! Flag for pgf_sal
   real(kind=RKIND) :: pgf_sal_on

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vel_tidal_potential_tend
!
!> \brief   Computes tendency term for tidal potential
!> \author  Steven Brus
!> \date    April 2020
!> \details
!>  This routine computes the tidal potential tendency for momentum
!>  based on current state.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_tidal_potential_tend(ssh, surfacePressure, tend, err)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(in) :: &
         ssh, surfacePressure             !< [in] Sea surface height

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend            !< [inout] accumulated velocity tendency

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< [out] error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer ::       &
         iEdge, k,     &! loop indices for edge, vertical loops
         cell1, cell2, &! neighbor cell indices across edge
         kMax           ! deepest active layer on edge

      real (kind=RKIND) :: &
         invdcEdge,    &! 1/dcEdge
         potentialGrad  ! tmp potential gradient term

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      !*** Initialize error code and return early if not turned on
      err = 0
      if (tidalPotentialOff) return

      if (timeIntegratorLTS) then
#ifdef MPAS_OPENACC
         !$acc parallel loop &
         !$acc    present(cellsOnEdge, maxLevelEdgeTop, dcEdge, edgeMask, &
         !$acc            tidalPotEta, ssh, surfacePressure, tend, pgf_sal) &
         !$acc    private(cell1, cell2, invdcEdge, potentialGrad, k, kMax)
#else
         !$omp parallel do schedule(runtime) &
         !$omp    private(cell1, cell2, invdcEdge, potentialGrad, k, kMax)
#endif
         do iEdge=1,nEdgesOwned
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            invdcEdge = 1.0_RKIND / dcEdge(iEdge)
            kMax = maxLevelEdgeTop(iEdge)

            potentialGrad = - gravity*invdcEdge* &
               ( tidalPotEta(cell2) - tidalPotEta(cell1) &
                + (1.0_RKIND - pgf_sal_on) * betaSelfAttrLoad * ( rho0gInv*(surfacePressure(cell2) - surfacePressure(cell1))) &
                + pgf_sal_on * (pgf_sal(cell2) - pgf_sal(cell1)) )

            do k=1,kMax
               tend(k,iEdge) = tend(k,iEdge) - &
                               edgeMask(k,iEdge)*potentialGrad
            end do
         end do
#ifndef MPAS_OPENACC
         !$omp end parallel do
#endif

      else ! not LTS

#ifdef MPAS_OPENACC
         !$acc parallel loop &
         !$acc    present(cellsOnEdge, maxLevelEdgeTop, dcEdge, edgeMask, &
         !$acc            tidalPotEta, ssh, surfacePressure, tend, pgf_sal) &
         !$acc    private(cell1, cell2, invdcEdge, potentialGrad, k, kMax)
#else
         !$omp parallel do schedule(runtime) &
         !$omp    private(cell1, cell2, invdcEdge, potentialGrad, k, kMax)
#endif
         do iEdge=1,nEdgesOwned
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            invdcEdge = 1.0_RKIND / dcEdge(iEdge)
            kMax = maxLevelEdgeTop(iEdge)

            potentialGrad = - gravity*invdcEdge* &
               (  tidalPotEta(cell2) - tidalPotEta(cell1) &
               + (1.0_RKIND - pgf_sal_on) * betaSelfAttrLoad *( (ssh(cell2) - ssh(cell1)) & 
                                                        + rho0gInv*(surfacePressure(cell2) - surfacePressure(cell1))) &
               + pgf_sal_on * (pgf_sal(cell2) - pgf_sal(cell1)) )

            do k=1,kMax
               tend(k,iEdge) = tend(k,iEdge) - &
                               edgeMask(k,iEdge)*potentialGrad
            end do
         end do
#ifndef MPAS_OPENACC
         !$omp end parallel do
#endif

   endif

   end subroutine ocn_vel_tidal_potential_tend!}}}

!***********************************************************************
!
!  routine ocn_compute_tidal_potential_forcing
!
!> \brief   Computes equilibrium tidal potential 
!> \author  Steven Brus
!> \date    September 2019
!> \details
!>  This routine computes the equilibrium tidal potential due to all 
!>  tidal constituents 
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_tidal_potential_forcing(err)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: Error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: &
         iCell, jCon,  &! loop indices for cell, constituent loops
         conType        ! id for constituent type

      real (kind=RKIND) :: &
         lon,&
         tArg,&
         ramp,&
         t, &
         nCycles,&
         period

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      !*** Set error flag and return early if not turned on
      err = 0
      if (tidalPotentialOff) return

      ramp = tanh((2.0_RKIND*daysSinceStartOfSim)/tidalPotRamp)
      t = daysSinceStartOfSim*86400.0_RKIND

      do iCell = 1, nCellsAll
         tidalPotEta(iCell) = 0.0_RKIND
      end do

      !*** Compute eta by summing all constituent contributions
      do jCon = 1, nTidalConstituents
         period = 2.0_RKIND*pi/tidalConstituentFrequency(jCon)
         nCycles = real(int(t/period),RKIND)
         targ = tidalConstituentFrequency(jCon)*(t - nCycles*period) + &
                tidalConstituentNodalPhase(jCon) + &
                tidalConstituentAstronomical(jCon)
         conType = tidalConstituentType(jCon)
         do iCell = 1, nCellsAll
           lon = lonCell(iCell)
           tidalPotEta(iCell) = tidalPotEta(iCell) + &
                         ramp * tidalConstituentAmplitude(jCon) &
                              * tidalConstituentNodalAmplitude(jCon) &
                              * tidalConstituentLoveNumbers(jCon) &
                              * latitudeFunction(iCell,conType+1) &
                              * cos(tArg  + real(conType,RKIND)*lon)
         end do
        
      end do

   end subroutine ocn_compute_tidal_potential_forcing!}}}


!***********************************************************************
!
!  routine ocn_vel_tidal_potential_init
!
!> \brief   Initializes ocean tidal protential forcing module.
!> \author  Steven Brus
!> \date    September 2019
!> \details
!>  This routine initializes the ocean tidal potential forcing module 
!>  and variables.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_tidal_potential_init(domain,err)!{{{

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< [out] error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      type (MPAS_Time_Type) :: &
         refTime            ! internal reference time

      integer :: &
         iCell, iCon        ! loop indices for cell, constituent loops

      real (kind=RKIND) :: &
         lat                ! latitude temporary

      ! pointers for pool variables
      type (block_type), pointer :: block_ptr 
      type (mpas_pool_type), pointer :: forcingPool
      type (mpas_pool_type), pointer :: meshPool
      integer, pointer :: nCells
      real (kind=RKIND), dimension(:), pointer :: latCell

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      !*** Set error flags and default values
      err = 0
      tidalPotentialOff = .true.
      betaSelfAttrLoad  = 0.0_RKIND
      tidalPotRamp      = 1.0e20_RKIND
      rho0gInv = 1.0_RKIND / rho_sw / gravity
      !*** If tidal potential turned on, set all relevant variables

      if (config_use_tidal_potential_forcing) then

         tidalPotentialOff = .false.
         if (config_use_self_attraction_loading) then
            pgf_sal_on = 1.0_RKIND
         else
            pgf_sal_on = 0.0_RKIND
         endif

         betaSelfAttrLoad = config_self_attraction_and_loading_beta
         tidalPotRamp     = config_tidal_potential_ramp

         !*** retrieve variables from pool
         !*** eventually replace with allocated module private vars

         block_ptr => domain % blocklist

         call mpas_pool_get_subpool(block_ptr%structs, 'forcing', &
                                                      forcingPool)

         call mpas_pool_get_subpool(block_ptr%structs, 'mesh', &
                                                      meshPool)

         call mpas_pool_get_dimension(block_ptr % dimensions, 'nCells', nCells)

         call mpas_pool_get_array(meshPool, 'latCell', latCell)

         call mpas_pool_get_array(forcingPool, 'tidalPotentialEta', &
                                                tidalPotEta)

         call mpas_pool_get_array(forcingPool, &
                            'tidalPotentialConstituentAmplitude', & 
                             tidalConstituentAmplitude)
         call mpas_pool_get_array(forcingPool, &
                            'tidalPotentialConstituentFrequency', &
                            tidalConstituentFrequency)
         call mpas_pool_get_array(forcingPool, &
                            'tidalPotentialConstituentLoveNumbers', &
                            tidalConstituentLoveNumbers)
         call mpas_pool_get_array(forcingPool, &
                            'tidalPotentialConstituentNodalAmplitude', &
                            tidalConstituentNodalAmplitude)
         call mpas_pool_get_array(forcingPool, &
                            'tidalPotentialConstituentAstronomical', &
                            tidalConstituentAstronomical)
         call mpas_pool_get_array(forcingPool, &
                            'tidalPotentialConstituentNodalPhase', &
                            tidalConstituentNodalPhase)
         call mpas_pool_get_array(forcingPool, &
                            'tidalPotentialConstituentType', &
                            tidalConstituentType)
         call mpas_pool_get_array(forcingPool, &
                            'tidalPotentialLatitudeFunction', &
                            latitudeFunction)

         call mpas_set_time(refTime, &
               dateTimeString=config_tidal_potential_reference_time)

         do iCell = 1,nCells
            tidalPotEta(iCell) = 0.0_RKIND
         end do

         nTidalConstituents = 0
         if (config_use_tidal_potential_forcing_M2) then
            nTidalConstituents = nTidalConstituents + 1
            constituentList(nTidalConstituents)%constituent = 'M2'
         end if 

         if (config_use_tidal_potential_forcing_S2) then
            nTidalConstituents = nTidalConstituents + 1
            constituentList(nTidalConstituents)%constituent = 'S2'
         end if 
 
         if (config_use_tidal_potential_forcing_N2) then
            nTidalConstituents = nTidalConstituents + 1
            constituentList(nTidalConstituents)%constituent = 'N2'
         end if 

         if (config_use_tidal_potential_forcing_K2) then
            nTidalConstituents = nTidalConstituents + 1
            constituentList(nTidalConstituents)%constituent = 'K2'
         end if 

         if (config_use_tidal_potential_forcing_K1) then
            nTidalConstituents = nTidalConstituents + 1
            constituentList(nTidalConstituents)%constituent = 'K1'
         end if 

         if (config_use_tidal_potential_forcing_O1) then
            nTidalConstituents = nTidalConstituents + 1
            constituentList(nTidalConstituents)%constituent = 'O1'
         end if 

         if (config_use_tidal_potential_forcing_Q1) then
            nTidalConstituents = nTidalConstituents + 1
            constituentList(nTidalConstituents)%constituent = 'Q1'
         end if 

         if (config_use_tidal_potential_forcing_P1) then
            nTidalConstituents = nTidalConstituents + 1
            constituentList(nTidalConstituents)%constituent = 'P1'
         end if 

         call tidal_constituent_factors(constituentList, &
                                       nTidalConstituents, refTime, &
                                       tidalConstituentFrequency, &
                                       tidalConstituentAmplitude, &
                                       tidalConstituentLoveNumbers, &
                                       tidalConstituentNodalAmplitude, &
                                       tidalConstituentAstronomical, &
                                       tidalConstituentNodalPhase, &
                                       tidalConstituentType, &
                                       err)

         do iCell = 1,nCells
            lat = latCell(iCell)
            latitudeFunction(iCell,1) = 3.0_RKIND*sin(lat)**2 -1.0_RKIND 
            latitudeFunction(iCell,2) = sin(2.0_RKIND*lat)
            latitudeFunction(iCell,3) = cos(lat)**2
         end do
        
         do iCon = 1,nTidalConstituents
            call mpas_log_write( &
                 'Constituent '//constituentList(iCon)%constituent)
            call mpas_log_write('  Frequency $r', &
                 realArgs=(/ tidalConstituentFrequency(iCon) /))
            call mpas_log_write('  Amplitude $r', &
                 realArgs=(/ tidalConstituentAmplitude(iCon) /))
            call mpas_log_write('  LoveNumbers $r', &
                 realArgs=(/ tidalConstituentLoveNumbers(iCon) /))
            call mpas_log_write('  NodalAmplitude $r', &
                 realArgs=(/ tidalConstituentNodalAmplitude(iCon) /))
            call mpas_log_write('  Astronomical argument $r', &
                 realArgs=(/ tidalConstituentAstronomical(iCon) /))
            call mpas_log_write('  NodalPhase $r', &
                 realArgs=(/ tidalConstituentNodalPhase(iCon) /))
            call mpas_log_write('  Type $i', &
                 intArgs=(/ tidalConstituentType(iCon) /))
            call mpas_log_write(' ')
         end do

         if ( (config_time_integrator == 'LTS') & 
              .or. (config_time_integrator == 'FB_LTS') ) then
            timeIntegratorLTS = .true.
         else
            timeIntegratorLTS = .false.
         endif

      end if ! tidal potential on

      !-----------------------------------------------------------------

   end subroutine ocn_vel_tidal_potential_init!}}}


!***********************************************************************
!
!  routine tidal_constituent_factors
!
!> \brief   Initialized tidal constituent factors 
!> \author  Steven Brus
!> \date    September 2019
!> \details
!>  This routine initializes the amplitude, frequency, love numbers, 
!>  astronomical argument and nodal factors for each tidal constituent
!>  Nodal factor equations are from:
!>         "Manual of Harmonic Analysis and Prediction of Tides"
!>          by Paul Schureman, Special Publication #98, US Coast
!>          and Geodetic Survey, Department of Commerce (1958).
!
!-----------------------------------------------------------------------

   subroutine tidal_constituent_factors(constituentList,nTidalConstituents,refTime, &
                                     tidalConstituentFrequency, &
                                     tidalConstituentAmplitude, &
                                     tidalConstituentLoveNumbers, &
                                     tidalConstituentNodalAmplitude, &
                                     tidalConstituentAstronomical, &
                                     tidalConstituentNodalPhase, &
                                     tidalConstituentType, &
                                     iErr) !{{{

      type(char_array), dimension(:), intent(in) :: constituentList
      integer, intent(in) :: nTidalConstituents
      type (MPAS_Time_Type), intent(in) :: refTime
      real (kind=RKIND), dimension(:), intent(out) :: tidalConstituentAmplitude
      real (kind=RKIND), dimension(:), intent(out) :: tidalConstituentFrequency
      real (kind=RKIND), dimension(:), intent(out) :: tidalConstituentLoveNumbers
      real (kind=RKIND), dimension(:), intent(out) :: tidalConstituentNodalAmplitude
      real (kind=RKIND), dimension(:), intent(out) :: tidalConstituentAstronomical
      real (kind=RKIND), dimension(:), intent(out) :: tidalConstituentNodalPhase
      integer, dimension(:), intent(out) :: tidalConstituentType
      integer, intent(out) :: iErr

      integer :: j
      integer :: refDay
      integer :: refYear
      integer :: refHour
      !integer :: inty,time_mjd
      !real (kind=RKIND) :: T,s,h,p,N
      real (kind=RKIND) :: deg2rad
      real (kind=RKIND) :: s,p,h,p1,n,I,nu,xi,nup,nup2,pc,T

      iErr = 0

      call mpas_get_time(refTime, YYYY=refYear,DoY=refDay,H=refHour)
      refDay = refDay + 1

      ! Calculate locations of sun and moon
      call orbit(refYear,refDay,refHour, &
                 s,p,h,p1,N,I,nu,xi,nup,nup2,pc)


      !! From tide design doc
      !inty = floor((real(refYear,RKIND)-1857.0_RKIND)/4.0_RKIND) - 1
      !time_mjd = 365*(refYear-1858)+inty-(31+28+31+30+31+30+31+31+30+31+17) + 1
      !T =real(time_mjd,RKIND)-51544.4993_RKIND;
      !s = 218.3164_RKIND + 13.17639648_RKIND*T;
      !h = 280.4661_RKIND + 0.98564736_RKIND*T;
      !p = 83.3535_RKIND + 0.11140353_RKIND*T;
      !N = 125.0445_RKIND - 0.05295377_RKIND*T;
      !N = N*pi/180.0_RKIND
      
      !! M2
          !tidalConstituentAstronomical(j) = 2.0_RKIND*h - 2.0_RKIND*s 
          !tidalConstituentNodalPhase(j) = -2.1_RKIND*sin(N)
          !tidalConstituentNodalAmplitude(j) = 1.0_RKIND - 0.037_RKIND*cos(N) 

      !! S2
          !tidalConstituentAstronomical(j) = 0.0_RKIND
          !tidalConstituentNodalPhase(j) = 0.0_RKIND
          !tidalConstituentNodalAmplitude(j) = 1.0_RKIND

      !! N2
          !tidalConstituentAstronomical(j) = 2.0_RKIND*h - 3.0_RKIND*s + p
          !tidalConstituentNodalPhase(j) = -2.1_RKIND*sin(N)
          !tidalConstituentNodalAmplitude(j) = 1.0_RKIND - 0.037_RKIND*cos(N)

      !! K1
          !tidalConstituentAstronomical(j) = h + 90.0_RKIND
          !tidalConstituentNodalPhase(j) = -8.9_RKIND*sin(N)
          !tidalConstituentNodalAmplitude(j) = 1.006_RKIND + 0.115_RKIND*cos(N)

      !! O1
          !tidalConstituentAstronomical(j) = h - 2.0_RKIND*s - 90.0_RKIND 
          !tidalConstituentNodalPhase(j) = 10.8_RKIND*sin(N)
          !tidalConstituentNodalAmplitude(j) = 1.009_RKIND + 0.187_RKIND*cos(N)


      deg2rad = pi/180.0_RKIND
      T = adjust_angle(180.0_RKIND + real(refHour,RKIND)*(360.0_RKIND/24.0_RKIND))

      do j = 1,nTidalConstituents
        if (constituentList(j)%constituent == 'M2') then
          tidalConstituentAmplitude(j) = 0.242334_RKIND
          tidalConstituentFrequency(j) = 1.405189e-4_RKIND
          tidalConstituentLoveNumbers(j) = 0.693_RKIND
          tidalConstituentAstronomical(j) = 0.0_RKIND 
          tidalConstituentNodalAmplitude(j) = cos(0.5_RKIND*I*deg2rad)**4/0.91544_RKIND
          tidalConstituentNodalPhase(j) = 2.0_RKIND*(T-s+h) + 2.0_RKIND*(xi-nu)
          tidalConstituentType(j) = 2

        else if (constituentList(j)%constituent == 'S2') then
          tidalConstituentAmplitude(j) = 0.112743_RKIND
          tidalConstituentFrequency(j) = 1.454441e-4_RKIND
          tidalConstituentLoveNumbers(j) = 0.693_RKIND
          tidalConstituentAstronomical(j) = 0.0_RKIND 
          tidalConstituentNodalAmplitude(j) = 1.0_RKIND
          tidalConstituentNodalPhase(j) = 2.0_RKIND*T
          tidalConstituentType(j) = 2

        else if (constituentList(j)%constituent == 'N2') then
          tidalConstituentAmplitude(j) = 0.046397_RKIND
          tidalConstituentFrequency(j) = 1.378797e-4_RKIND
          tidalConstituentLoveNumbers(j) = 0.693_RKIND
          tidalConstituentAstronomical(j) = 0.0_RKIND 
          tidalConstituentNodalAmplitude(j) = cos(0.5_RKIND*I*deg2rad)**4/0.91544_RKIND
          tidalConstituentNodalPhase(j) = 2.0_RKIND*(T+h) - 3.0_RKIND*s+p + 2.0_RKIND*(xi-nu)
          tidalConstituentType(j) = 2

        else if (constituentList(j)%constituent == 'K2') then
          tidalConstituentAmplitude(j) = 0.030684_RKIND 
          tidalConstituentFrequency(j) = 1.458423e-4_RKIND
          tidalConstituentLoveNumbers(j) = 0.693_RKIND 
          tidalConstituentAstronomical(j) = 0.0_RKIND 
          tidalConstituentNodalAmplitude(j) = 0.001_RKIND+sqrt(19.0444_RKIND*sin(I*deg2rad)**4 + &
                                                               2.7702_RKIND*sin(I*deg2rad)**2*cos(2.0_RKIND*nu*deg2rad) + &
                                                               0.0981_RKIND) 
          tidalConstituentNodalPhase(j) = 2.0_RKIND*(T+h) - 2.0_RKIND*nup2 
          tidalConstituentType(j) = 2

        else if (constituentList(j)%constituent == 'K1') then
          tidalConstituentAmplitude(j) = 0.141565_RKIND
          tidalConstituentFrequency(j) = 0.7292117e-4_RKIND
          tidalConstituentLoveNumbers(j) = 0.736_RKIND
          tidalConstituentAstronomical(j) = 0.0_RKIND 
          tidalConstituentNodalAmplitude(j) = sqrt(0.8965_RKIND*sin(2.0_RKIND*I*deg2rad)**2 + &
                                                   0.6001_RKIND*sin(2.0_RKIND*I*deg2rad) * &
                                                   cos(nu*deg2rad) + 0.1006_RKIND)
          tidalConstituentNodalPhase(j) = T + h - 90.0_RKIND - nup
          tidalConstituentType(j) = 1

        else if (constituentList(j)%constituent == 'O1') then
          tidalConstituentAmplitude(j) = 0.100661_RKIND
          tidalConstituentFrequency(j) = 0.6759774e-4_RKIND
          tidalConstituentLoveNumbers(j) = 0.695_RKIND
          tidalConstituentAstronomical(j) = 0.0_RKIND 
          tidalConstituentNodalAmplitude(j) = sin(I*deg2rad)*cos(0.5_RKIND*I*deg2rad)**2/0.37988_RKIND
          tidalConstituentNodalPhase(j) = T - 2.0_RKIND*s + h + 90.0_RKIND + 2.0_RKIND*xi - nu
          tidalConstituentType(j) = 1

        else if (constituentList(j)%constituent == 'Q1') then
          tidalConstituentAmplitude(j) = 0.019273_RKIND 
          tidalConstituentFrequency(j) = 0.6495854e-4_RKIND
          tidalConstituentLoveNumbers(j) = 0.695_RKIND
          tidalConstituentAstronomical(j) = 0.0_RKIND 
          tidalConstituentNodalAmplitude(j) = sin(I*deg2rad)*cos(0.5_RKIND*I*deg2rad)**2/0.37988_RKIND
          tidalConstituentNodalPhase(j) = T - 3.0_RKIND*s + h + p + 90.0_RKIND + 2.0_RKIND*xi - nu 
          tidalConstituentType(j) = 1

        else if (constituentList(j)%constituent == 'P1') then
          tidalConstituentAmplitude(j) = 0.046848_RKIND
          tidalConstituentFrequency(j) = 0.7252295e-4_RKIND
          tidalConstituentLoveNumbers(j) = 0.706_RKIND
          tidalConstituentAstronomical(j) = 0.0_RKIND 
          tidalConstituentNodalAmplitude(j) = 1.0_RKIND 
          tidalConstituentNodalPhase(j) = T - h + 90.0_RKIND 
          tidalConstituentType(j) = 1

        else
          call mpas_log_write( 'Validation failed for tidal potential.' &
            // 'Not given a usable value for constituent name.',MPAS_LOG_CRIT)
          iErr = 1
        end if

        tidalConstituentNodalPhase(j) = adjust_angle(tidalConstituentNodalPhase(j))
        tidalConstituentNodalPhase(j) = tidalConstituentNodalPhase(j)*deg2rad
        tidalConstituentAstronomical(j) = tidalConstituentAstronomical(j)*deg2rad
      end do


  end subroutine tidal_constituent_factors!}}}

!***********************************************************************
!
!  routine orbit
!
!> \brief   Computes orbital factors 
!> \author  Steven Brus
!> \date    September 2019
!> \details
!>  Computes factors related to the orbits of the sun and moon for 
!>  calculating the nodal factors in the tidal potential.
!>  Adapted from tide_fac.f which can be found here:
!>  http://adcirc.org/home/related-software/adcirc-utility-programs/
!
!-----------------------------------------------------------------------
  subroutine orbit(year,julianDay,hour, &
                   s,p,h,p1,N,I,nu,xi,nup,nup2,pc) !{{{
      implicit none
      integer, intent(in) :: year
      integer, intent(in) :: julianDay
      integer, intent(in) :: hour
      real (kind = RKIND), intent(out) :: s
      real (kind = RKIND), intent(out) :: p
      real (kind = RKIND), intent(out) :: h
      real (kind = RKIND), intent(out) :: p1
      real (kind = RKIND), intent(out) :: N
      real (kind = RKIND), intent(out) :: I
      real (kind = RKIND), intent(out) :: nu
      real (kind = RKIND), intent(out) :: xi
      real (kind = RKIND), intent(out) :: nup
      real (kind = RKIND), intent(out) :: nup2
      real (kind = RKIND), intent(out) :: pc 

      integer :: x
      real (kind=RKIND) :: yr,day,hr
      real (kind=RKIND) :: deg2rad,rad2deg
      real (kind=RKIND) :: NRad,pRad,IRad,nuRad,xiRad,nupRad,nup2Rad

      deg2rad = pi/180.0_RKIND
      rad2deg = 180.0_RKIND/pi

      x = int((real(year,RKIND)-1901.0_RKIND)/4.0_RKIND)
      yr = real(year,RKIND) - 1900.0_RKIND
      day = real(julianDay,RKIND) + real(x,RKIND) - 1.0_RKIND
      hr = real(hour,RKIND)

      ! N is the longitude of the moon's node (capital N, Table 1, Schureman)
      N = 259.1560564_RKIND - 19.328185764_RKIND*yr - .0529539336_RKIND*day - .0022064139_RKIND*hr
      N = adjust_angle(N)
      NRad = N*deg2rad

      ! p is the lunar perigee (small p, Table 1)
      p = 334.3837214_RKIND + 40.66246584_RKIND*yr + .111404016_RKIND*day + .004641834_RKIND*hr
      p = adjust_angle(p)
      pRad = p*deg2rad

      IRad = ACOS(.9136949_RKIND - .0356926_RKIND*COS(NRad))
      I = adjust_angle(IRad*rad2deg)

      nuRad = ASIN(.0897056*SIN(NRad)/SIN(IRad))
      nu = nuRad*rad2deg

      xiRad = NRad-2.0_RKIND*ATAN(.64412_RKIND*TAN(NRad/2.0_RKIND))-nuRad
      xi = xiRad*rad2deg

      pc = adjust_angle(p-xi)

      ! h is the mean longitude of the sun (small h, Table 1)
      h = 280.1895014_RKIND - 0.238724988_RKIND*yr + 0.9856473288_RKIND*day + 0.0410686387_RKIND*hr
      h = adjust_angle(h)

      ! p1 is the solar perigee (small p1, Table 1)
      p1 = 281.2208569_RKIND + 0.01717836_RKIND*yr + 0.000047064_RKIND*day + 0.000001961_RKIND*hr
      p1 = adjust_angle(p1)

      ! s is the mean longitude of the moon (small s, Table 1)
      s = 277.0256206_RKIND + 129.38482032_RKIND*yr + 13.176396768_RKIND*day + 0.549016532_RKIND*hr
      s = adjust_angle(s)

      nupRad = ATAN(SIN(nuRad)/(COS(nuRad) + 0.334766_RKIND/SIN(2.0_RKIND*IRad)))
      nup = nupRad*rad2deg

      nup2Rad = ATAN(SIN(2.0_RKIND*nuRad)/(COS(2.0_RKIND*nuRad) + 0.0726184_RKIND/SIN(IRad)**2))/2.0_RKIND
      nup2 = nup2Rad*rad2deg

  end subroutine orbit!}}}

!***********************************************************************
!
!  function adjust_angle 
!
!> \brief   Puts angles into the 0-360 degree range
!> \author  Steven Brus
!> \date    September 2019
!> \details 
!
!-----------------------------------------------------------------------
  function adjust_angle(arg) result(angle) !{{{

      implicit none

      real (kind=RKIND), intent(in) :: arg
      real (kind=RKIND) :: angle

      integer :: cycles

      cycles = int(arg/360.0_RKIND)
      angle =  arg - real(cycles,RKIND)*360.0_RKIND
      if (angle <  0.0_RKIND) then
        angle = angle + 360.0_RKIND
      end if
  END


!***********************************************************************

end module ocn_vel_tidal_potential!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
