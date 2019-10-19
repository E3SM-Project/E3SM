! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tidal_potential_forcing
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

module ocn_tidal_potential_forcing

   use mpas_kind_types
   use mpas_constants
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timekeeping
   use mpas_timer
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

   public :: ocn_compute_tidal_potential_forcing, &
             ocn_tidal_potential_forcing_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: tidalPotentialOn
   type :: char_array
     character(:), allocatable :: constituent
   end type
   type(char_array), dimension(37) :: constituentList


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tidal_potential_forcing_layer_thickness
!
!> \brief   Computes equilibrium tidal potential 
!> \author  Steven Brus
!> \date    September 2019
!> \details
!>  This routine computes the equilibrium tidal potential due to all 
!>  tidal constituents 
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_tidal_potential_forcing(meshPool, forcingPool, diagnosticsPool, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: meshPool !< Input: mesh information
      type (mpas_pool_type), intent(in) :: diagnosticsPool 

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(inout) :: forcingPool 

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: Error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer, pointer :: nCells
      real (kind=RKIND), pointer :: tidalPotentialRamp
      real (kind=RKIND), dimension(:), pointer :: lonCell
      integer, dimension(:), pointer ::  maxLevelCell

      integer, pointer :: nTidalConstituents
      real (kind=RKIND), dimension(:), pointer :: amplitude, frequency, loveNumbers
      real (kind=RKIND), dimension(:), pointer :: nodalFactorAmplitude, nodalFactorPhase, astronomicalArgument
      real (kind=RKIND), dimension(:,:), pointer :: latitudeFunction
      integer, dimension(:), pointer :: constituentType
      real (kind=RKIND), dimension(:), pointer :: eta
      real (kind=RKIND), dimension(:,:), pointer :: zMid, tidalPotentialZMid
      real (kind=RKIND), pointer :: daysSinceStartOfSim
      
      integer :: iCell
      integer :: jCon, conType
      real (kind=RKIND) :: lon,tArg,ramp,t
      real (kind=RKIND) :: nCycles,period

      err = 0

      if ( .not. tidalPotentialOn ) return

      call mpas_pool_get_config(ocnConfigs, 'config_tidal_potential_ramp', tidalPotentialRamp)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
      call mpas_pool_get_array(forcingPool, 'nTidalPotentialConstituents', nTidalConstituents)
      call mpas_pool_get_array(forcingPool, 'tidalPotentialConstituentAmplitude', amplitude)
      call mpas_pool_get_array(forcingPool, 'tidalPotentialConstituentFrequency', frequency)
      call mpas_pool_get_array(forcingPool, 'tidalPotentialConstituentLoveNumbers', loveNumbers)
      call mpas_pool_get_array(forcingPool, 'tidalPotentialConstituentNodalAmplitude', nodalFactorAmplitude)
      call mpas_pool_get_array(forcingPool, 'tidalPotentialConstituentAstronomical', astronomicalArgument)
      call mpas_pool_get_array(forcingPool, 'tidalPotentialConstituentNodalPhase', nodalFactorPhase)
      call mpas_pool_get_array(forcingPool, 'tidalPotentialConstituentType', constituentType)
      call mpas_pool_get_array(forcingPool, 'tidalPotentialLatitudeFunction', latitudeFunction)
      call mpas_pool_get_array(forcingPool, 'tidalPotentialEta', eta)
      call mpas_pool_get_array(diagnosticsPool, "daysSinceStartOfSim", daysSinceStartOfSim)

      ramp = tanh((2.0_RKIND*daysSinceStartOfSim)/tidalPotentialRamp)
      t = daysSinceStartOfSim*86400.0_RKIND

      do iCell = 1, nCells
        eta(iCell) = 0.0_RKIND
      end do

      do jCon = 1, nTidalConstituents
        period = 2.0_RKIND*pii/frequency(jCon)
        nCycles = real(int(t/period),RKIND)
        targ = frequency(jCon)*(t - nCycles*period) + nodalFactorPhase(jCon) + astronomicalArgument(jCon)
        conType = constituentType(jCon)
        do iCell = 1, nCells
          lon = lonCell(iCell)
          eta(iCell) = eta(iCell) + ramp &
                                  * amplitude(jCon) &
                                  * nodalFactorAmplitude(jCon) &
                                  * loveNumbers(jCon) &
                                  * latitudeFunction(iCell,conType+1) &
                                  * cos(tArg  + real(conType,RKIND)*lon)
        end do
        
      end do

   end subroutine ocn_compute_tidal_potential_forcing!}}}


!***********************************************************************
!
!  routine ocn_tidal_potential_forcing_init
!
!> \brief   Initializes ocean tidal protential forcing module.
!> \author  Steven Brus
!> \date    September 2019
!> \details
!>  This routine initializes the ocean tidal potential forcing module 
!>  and variables.
!
!-----------------------------------------------------------------------

   subroutine ocn_tidal_potential_forcing_init(domain,err)!{{{

      type (domain_type), intent(inout) :: domain
      integer, intent(out) :: err !< Output: error flag

      logical, pointer :: config_use_tidal_potential_forcing
      logical, pointer :: config_use_tidal_potential_forcing_M2
      logical, pointer :: config_use_tidal_potential_forcing_S2
      logical, pointer :: config_use_tidal_potential_forcing_N2
      logical, pointer :: config_use_tidal_potential_forcing_K2
      logical, pointer :: config_use_tidal_potential_forcing_K1
      logical, pointer :: config_use_tidal_potential_forcing_O1
      logical, pointer :: config_use_tidal_potential_forcing_Q1
      logical, pointer :: config_use_tidal_potential_forcing_P1
      character (len=StrKIND), pointer :: config_tidal_potential_reference_time 

      type (block_type), pointer :: block_ptr 
      type (mpas_pool_type), pointer :: forcingPool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: diagnosticsPool 
      integer, pointer :: nTidalConstituents
      integer, pointer :: nCells 
      real (kind=RKIND), dimension(:), pointer :: tidalConstituentAmplitude
      real (kind=RKIND), dimension(:), pointer :: tidalConstituentFreqency
      real (kind=RKIND), dimension(:), pointer :: tidalConstituentLoveNumbers
      real (kind=RKIND), dimension(:), pointer :: tidalConstituentNodalAmplitude
      real (kind=RKIND), dimension(:), pointer :: tidalConstituentAstronomical
      real (kind=RKIND), dimension(:), pointer :: tidalConstituentNodalPhase
      real (kind=RKIND), dimension(:), pointer :: latCell 
      real (kind=RKIND), dimension(:), pointer :: eta
      real (kind=RKIND), dimension(:,:), pointer :: latitudeFunction 
      integer, dimension(:), pointer :: tidalConstituentType
      type (MPAS_Time_Type) :: refTime
      integer :: iCell,iCon
      real (kind=RKIND) :: lat


      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_use_tidal_potential_forcing', config_use_tidal_potential_forcing)
      call mpas_pool_get_config(ocnConfigs, 'config_use_tidal_potential_forcing_M2', config_use_tidal_potential_forcing_M2)
      call mpas_pool_get_config(ocnConfigs, 'config_use_tidal_potential_forcing_S2', config_use_tidal_potential_forcing_S2)
      call mpas_pool_get_config(ocnConfigs, 'config_use_tidal_potential_forcing_N2', config_use_tidal_potential_forcing_N2)
      call mpas_pool_get_config(ocnConfigs, 'config_use_tidal_potential_forcing_K2', config_use_tidal_potential_forcing_K2)
      call mpas_pool_get_config(ocnConfigs, 'config_use_tidal_potential_forcing_K1', config_use_tidal_potential_forcing_K1)
      call mpas_pool_get_config(ocnConfigs, 'config_use_tidal_potential_forcing_O1', config_use_tidal_potential_forcing_O1)
      call mpas_pool_get_config(ocnConfigs, 'config_use_tidal_potential_forcing_Q1', config_use_tidal_potential_forcing_Q1)
      call mpas_pool_get_config(ocnConfigs, 'config_use_tidal_potential_forcing_P1', config_use_tidal_potential_forcing_P1)
      call mpas_pool_get_config(ocnConfigs, 'config_tidal_potential_reference_time', config_tidal_potential_reference_time)

      block_ptr => domain % blocklist
      do while(associated(block_ptr))
        call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
        call mpas_pool_get_array(forcingPool, 'tidalPotentialEta', eta)

        do iCell = 1,nCells
          eta(iCell) = 0.0_RKIND
        end do

        block_ptr => block_ptr % next
      end do      

      tidalPotentialOn = .false.
      if (config_use_tidal_potential_forcing) then
        tidalPotentialOn = .true.
      else
        return
      end if


      block_ptr => domain % blocklist
      do while(associated(block_ptr))

        call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)
        call mpas_pool_get_array(forcingPool, 'nTidalPotentialConstituents', nTidalConstituents)
        call mpas_pool_get_array(forcingPool, 'tidalPotentialConstituentAmplitude', tidalConstituentAmplitude)
        call mpas_pool_get_array(forcingPool, 'tidalPotentialConstituentFrequency', tidalConstituentFreqency)
        call mpas_pool_get_array(forcingPool, 'tidalPotentialConstituentLoveNumbers', tidalConstituentLoveNumbers)
        call mpas_pool_get_array(forcingPool, 'tidalPotentialConstituentNodalAmplitude', tidalConstituentNodalAmplitude)
        call mpas_pool_get_array(forcingPool, 'tidalPotentialConstituentAstronomical', tidalConstituentAstronomical)
        call mpas_pool_get_array(forcingPool, 'tidalPotentialConstituentNodalPhase', tidalConstituentNodalPhase)
        call mpas_pool_get_array(forcingPool, 'tidalPotentialConstituentType', tidalConstituentType)
        call mpas_pool_get_array(forcingPool, 'tidalPotentialLatitudeFunction', latitudeFunction)
        call mpas_pool_get_array(meshPool, 'latCell', latCell)
        call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

        call mpas_set_time(refTime, dateTimeString=config_tidal_potential_reference_time)

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

        call tidal_constituent_factors(constituentList,nTidalConstituents,refTime, &
                                       tidalConstituentFreqency, &
                                       tidalConstituentAmplitude, &
                                       tidalConstituentLoveNumbers, &
                                       tidalConstituentNodalAmplitude, &
                                       tidalConstituentAstronomical, &
                                       tidalConstituentNodalPhase, &
                                       tidalConstituentType, &
                                       err)

        do iCell = 1,nCells
          lat = latCell(iCell)
          latitudeFunction(iCell,1) = 3.0_RKIND*sin(lat)**2 - 1.0_RKIND 
          latitudeFunction(iCell,2) = sin(2.0_RKIND*lat)
          latitudeFunction(iCell,3) = cos(lat)**2
        end do
        
        do iCon = 1,nTidalConstituents  
          call mpas_log_write('Constituent '//constituentList(iCon)%constituent)
          call mpas_log_write('  Frequency $r', realArgs=(/ tidalConstituentFreqency(iCon) /))
          call mpas_log_write('  Amplitude $r', realArgs=(/ tidalConstituentAmplitude(iCon) /))
          call mpas_log_write('  LoveNumbers $r', realArgs=(/ tidalConstituentLoveNumbers(iCon) /))
          call mpas_log_write('  NodalAmplitude $r', realArgs=(/ tidalConstituentNodalAmplitude(iCon) /))
          call mpas_log_write('  Astronomical argument $r', realArgs=(/ tidalConstituentAstronomical(iCon) /))
          call mpas_log_write('  NodalPhase $r', realArgs=(/ tidalConstituentNodalPhase(iCon) /))
          call mpas_log_write('  Type $i', intArgs=(/ tidalConstituentType(iCon) /))
          call mpas_log_write(' ')
        end do

        block_ptr => block_ptr % next
      end do      


   end subroutine ocn_tidal_potential_forcing_init!}}}


!***********************************************************************
!
!  routine tidal_constituent_factors
!
!> \brief   Initialized tidal constituent factors 
!> \author  Steven Brus
!> \date    September 2019
!> \details
!>  This routine initializes the ampiltude, frequency, love numbers, 
!>  astronomical argument and nodal factors for each tidal constituent
!>  Nodal factor equations are from:
!>         "Manual of Harmonic Analysis and Prediction of Tides"
!>          by Paul Schureman, Special Publication #98, US Coast
!>          and Geodetic Survey, Department of Commerce (1958).
!
!-----------------------------------------------------------------------

   subroutine tidal_constituent_factors(constituentList,nTidalConstituents,refTime, &
                                     tidalConstituentFreqency, &
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
      real (kind=RKIND), dimension(:), intent(out) :: tidalConstituentFreqency
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
      !N = N*pii/180.0_RKIND
      
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


      deg2rad = pii/180.0_RKIND
      T = adjust_angle(180.0_RKIND + real(refHour,RKIND)*(360.0_RKIND/24.0_RKIND))

      do j = 1,nTidalConstituents
        if (constituentList(j)%constituent == 'M2') then
          tidalConstituentAmplitude(j) = 0.242334_RKIND
          tidalConstituentFreqency(j) = 1.405189e-4_RKIND
          tidalConstituentLoveNumbers(j) = 0.693_RKIND
          tidalConstituentAstronomical(j) = 0.0_RKIND 
          tidalConstituentNodalAmplitude(j) = cos(0.5_RKIND*I*deg2rad)**4/0.91544_RKIND
          tidalConstituentNodalPhase(j) = 2.0_RKIND*(T-s+h) + 2.0_RKIND*(xi-nu)
          tidalConstituentType(j) = 2

        else if (constituentList(j)%constituent == 'S2') then
          tidalConstituentAmplitude(j) = 0.112743_RKIND
          tidalConstituentFreqency(j) = 1.454441e-4_RKIND
          tidalConstituentLoveNumbers(j) = 0.693_RKIND
          tidalConstituentAstronomical(j) = 0.0_RKIND 
          tidalConstituentNodalAmplitude(j) = 1.0_RKIND
          tidalConstituentNodalPhase(j) = 2.0_RKIND*T
          tidalConstituentType(j) = 2

        else if (constituentList(j)%constituent == 'N2') then
          tidalConstituentAmplitude(j) = 0.046397_RKIND
          tidalConstituentFreqency(j) = 1.378797e-4_RKIND
          tidalConstituentLoveNumbers(j) = 0.693_RKIND
          tidalConstituentAstronomical(j) = 0.0_RKIND 
          tidalConstituentNodalAmplitude(j) = cos(0.5_RKIND*I*deg2rad)**4/0.91544_RKIND
          tidalConstituentNodalPhase(j) = 2.0_RKIND*(T+h) - 3.0_RKIND*s+p + 2.0_RKIND*(xi-nu)
          tidalConstituentType(j) = 2

        else if (constituentList(j)%constituent == 'K2') then
          tidalConstituentAmplitude(j) = 0.030684_RKIND 
          tidalConstituentFreqency(j) = 1.458423e-4_RKIND
          tidalConstituentLoveNumbers(j) = 0.693_RKIND 
          tidalConstituentAstronomical(j) = 0.0_RKIND 
          tidalConstituentNodalAmplitude(j) = 0.001_RKIND+sqrt(19.0444_RKIND*sin(I*deg2rad)**4 + &
                                                               2.7702_RKIND*sin(I*deg2rad)**2*cos(2.0_RKIND*nu*deg2rad) + &
                                                               0.0981_RKIND) 
          tidalConstituentNodalPhase(j) = 2.0_RKIND*(T+h) - 2.0_RKIND*nup2 
          tidalConstituentType(j) = 2

        else if (constituentList(j)%constituent == 'K1') then
          tidalConstituentAmplitude(j) = 0.141565_RKIND
          tidalConstituentFreqency(j) = 0.7292117e-4_RKIND
          tidalConstituentLoveNumbers(j) = 0.736_RKIND
          tidalConstituentAstronomical(j) = 0.0_RKIND 
          tidalConstituentNodalAmplitude(j) = sqrt(0.8965_RKIND*sin(2.0_RKIND*I*deg2rad)**2 + &
                                                   0.6001_RKIND*sin(2.0_RKIND*I*deg2rad) * &
                                                   cos(nu*deg2rad) + 0.1006_RKIND)
          tidalConstituentNodalPhase(j) = T + h - 90.0_RKIND - nup
          tidalConstituentType(j) = 1

        else if (constituentList(j)%constituent == 'O1') then
          tidalConstituentAmplitude(j) = 0.100661_RKIND
          tidalConstituentFreqency(j) = 0.6759774e-4_RKIND
          tidalConstituentLoveNumbers(j) = 0.695_RKIND
          tidalConstituentAstronomical(j) = 0.0_RKIND 
          tidalConstituentNodalAmplitude(j) = sin(I*deg2rad)*cos(0.5_RKIND*I*deg2rad)**2/0.37988_RKIND
          tidalConstituentNodalPhase(j) = T - 2.0_RKIND*s + h + 90.0_RKIND + 2.0_RKIND*xi - nu
          tidalConstituentType(j) = 1

        else if (constituentList(j)%constituent == 'Q1') then
          tidalConstituentAmplitude(j) = 0.019273_RKIND 
          tidalConstituentFreqency(j) = 0.6495854e-4_RKIND
          tidalConstituentLoveNumbers(j) = 0.695_RKIND
          tidalConstituentAstronomical(j) = 0.0_RKIND 
          tidalConstituentNodalAmplitude(j) = sin(I*deg2rad)*cos(0.5_RKIND*I*deg2rad)**2/0.37988_RKIND
          tidalConstituentNodalPhase(j) = T - 3.0_RKIND*s + h + p + 90.0_RKIND + 2.0_RKIND*xi - nu 
          tidalConstituentType(j) = 1

        else if (constituentList(j)%constituent == 'P1') then
          tidalConstituentAmplitude(j) = 0.046848_RKIND
          tidalConstituentFreqency(j) = 0.7252295e-4_RKIND
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

      deg2rad = pii/180.0_RKIND
      rad2deg = 180.0_RKIND/pii

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

end module ocn_tidal_potential_forcing!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
