! copyright (c) 2013,  los alamos national security, llc (lans)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_CFC
!
!> \brief MPAS ocean CFC
!> \author Mathew Maltrud
!> \date   12/22/2016
!> \details
!>  This module contains routines for computing evolution of CFC11 and CFC12
!
!-----------------------------------------------------------------------

module ocn_tracer_CFC

   use mpas_timer
   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_forcing
   use ocn_framework_forcing
   use mpas_timekeeping
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

   public :: ocn_tracer_CFC_compute, &
             ocn_tracer_CFC_surface_flux_compute,  &
             ocn_tracer_CFC_init,  &
             ocn_get_CFCData,  &
             ocn_CFC_forcing_write_restart

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!  solubility coefficients
   real (kind=RKIND), dimension(11:12) :: a1sol, a2sol, a3sol, a4sol, b1sol, b2sol, b3sol

!  schmidt number coefficients
   real (kind=RKIND), dimension(11:12) :: a1schmidt, a2schmidt, a3schmidt, a4schmidt

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tracer_CFC_compute
!
!> \brief   computes a tracer tendency due to CFC
!> \author  Mathew Maltrud
!> \date    12/22/2016
!> \details
!>  This routine computes a tracer tendency for CFC11 and/or CFC12
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_CFC_compute!{{{

! no CFC sources/sinks

   end subroutine ocn_tracer_CFC_compute!}}}

!***********************************************************************
!
!  routine ocn_tracer_CFC_surface_flux_compute
!
!> \brief   computes a tracer tendency due to CFC surface fluxes
!> \author  Mathew Maltrud
!> \date    12/22/2016
!> \details
!>  This routine computes a tracer tendency due to CFC surface fluxes
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_CFC_surface_flux_compute(activeTracers, CFCTracers, forcingPool, nTracers,   &
      nCellsSolve, indexTemperature, indexSalinity, CFCSurfaceFlux, indexCFC11, indexCFC12, err)!{{{


 use ocn_config

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      ! two dimensional arrays
      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         CFCSurfaceFlux

      ! three dimensional arrays
      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         CFCTracers
      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         activeTracers

      ! scalars
      integer, intent(in) :: nTracers, nCellsSolve, indexTemperature, indexSalinity, indexCFC11, indexCFC12

      type (mpas_pool_type), intent(inout) :: forcingPool

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

      integer, intent(out) :: err !< Output: Error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), pointer :: CFC11FluxDiagnostics,  &
                                        CFC12FluxDiagnostics,  &
                                        CFCAuxiliary

      integer :: numColumns, column, iCell, iLevelSurface

      real (kind=RKIND), dimension(:), pointer :: &
         pCFC11,   &
         pCFC12,   &
         atmosphericPressure,   &
         iceFraction,          &
         windSpeedSquared10m,  &
         CFC11_flux_ifrac,  &
         CFC11_flux_xkw,    &
         CFC11_flux_atm_press,  &
         CFC11_flux_pv,     &
         CFC11_flux_schmidt,&
         CFC11_flux_sat,    &
         CFC11_flux_surf,   &
         CFC11_flux_ws,     &
         CFC12_flux_ifrac,  &
         CFC12_flux_xkw,    &
         CFC12_flux_atm_press,  &
         CFC12_flux_pv,     &
         CFC12_flux_schmidt,&
         CFC12_flux_sat,    &
         CFC12_flux_surf,   &
         CFC12_flux_ws

      real (kind=RKIND) :: &
!        xkw_coeff = 8.6e-7_RKIND,  &       ! xkw_coeff = 0.0031 m/hr s^2/m^2, older OCMIP value
         xkw_coeff = 9.6361e-7_RKIND,  &    ! xkw_coeff = 0.00337 m/hr s^2/m^2
         PascalsToAtmospheres = 1.0_RKIND/101.325e+3_RKIND,  &
         mSquared_to_cmSquared = 1.0e+4_RKIND

      real (kind=RKIND) :: &
         xkw_ice,   &
         sst,   &
         sss,   &
         saturation_CFC11_1atm,   &
         saturation_CFC12_1atm,   &
         CFC_flux_ifrac_nonNegative,   &
         CFC_flux_ifrac_clamped

      err = 0

      call mpas_timer_start("CFC surface flux")

      call mpas_pool_get_array(forcingPool, 'atmosphericPressure', atmosphericPressure)
      call mpas_pool_get_array(forcingPool, 'iceFraction', iceFraction)
      call mpas_pool_get_subpool(forcingPool, 'CFCAuxiliary', CFCAuxiliary)
      call mpas_pool_get_array(CFCAuxiliary, 'windSpeedSquared10mCFC', windSpeedSquared10m)

      if (config_use_CFC11) then
        call mpas_pool_get_array(CFCAuxiliary, 'pCFC11', pCFC11)
        call mpas_pool_get_subpool(forcingPool, 'CFC11FluxDiagnostics', CFC11FluxDiagnostics)
        call mpas_pool_get_array(CFC11FluxDiagnostics, 'CFC11_flux_xkw', CFC11_flux_xkw)
        call mpas_pool_get_array(CFC11FluxDiagnostics, 'CFC11_flux_atm_press', CFC11_flux_atm_press)
        call mpas_pool_get_array(CFC11FluxDiagnostics, 'CFC11_flux_pv', CFC11_flux_pv)
        call mpas_pool_get_array(CFC11FluxDiagnostics, 'CFC11_flux_schmidt', CFC11_flux_schmidt)
        call mpas_pool_get_array(CFC11FluxDiagnostics, 'CFC11_flux_sat', CFC11_flux_sat)
        call mpas_pool_get_array(CFC11FluxDiagnostics, 'CFC11_flux_surf', CFC11_flux_surf)
        call mpas_pool_get_array(CFC11FluxDiagnostics, 'CFC11_flux_ws', CFC11_flux_ws)
        call mpas_pool_get_array(CFC11FluxDiagnostics, 'CFC11_flux_ifrac', CFC11_flux_ifrac)
      endif

      if (config_use_CFC12) then
        call mpas_pool_get_array(CFCAuxiliary, 'pCFC12', pCFC12)
        call mpas_pool_get_subpool(forcingPool, 'CFC12FluxDiagnostics', CFC12FluxDiagnostics)
        call mpas_pool_get_array(CFC12FluxDiagnostics, 'CFC12_flux_xkw', CFC12_flux_xkw)
        call mpas_pool_get_array(CFC12FluxDiagnostics, 'CFC12_flux_atm_press', CFC12_flux_atm_press)
        call mpas_pool_get_array(CFC12FluxDiagnostics, 'CFC12_flux_pv', CFC12_flux_pv)
        call mpas_pool_get_array(CFC12FluxDiagnostics, 'CFC12_flux_schmidt', CFC12_flux_schmidt)
        call mpas_pool_get_array(CFC12FluxDiagnostics, 'CFC12_flux_sat', CFC12_flux_sat)
        call mpas_pool_get_array(CFC12FluxDiagnostics, 'CFC12_flux_surf', CFC12_flux_surf)
        call mpas_pool_get_array(CFC12FluxDiagnostics, 'CFC12_flux_ws', CFC12_flux_ws)
        call mpas_pool_get_array(CFC12FluxDiagnostics, 'CFC12_flux_ifrac', CFC12_flux_ifrac)
      endif

      numColumns = 1
      column = 1
      iLevelSurface = 1

      !DWJ 08/05/2016: This loop needs OpenMP added to it.
      do iCell=1,nCellsSolve

       CFC_flux_ifrac_nonNegative = max(iceFraction(iCell), 0.0_RKIND)
       CFC_flux_ifrac_clamped = min(CFC_flux_ifrac_nonNegative, 1.0_RKIND)

       xkw_ice = (1.0_RKIND - CFC_flux_ifrac_clamped) * xkw_coeff * windSpeedSquared10m(iCell)

       sst = activeTracers(indexTemperature,iLevelSurface,iCell)
       sss = activeTracers(indexSalinity,iLevelSurface,iCell)

       if (config_use_CFC11) then

          CFC11_flux_schmidt(iCell) = SCHMIDT_CFC(sst, 11)

          saturation_CFC11_1atm = pCFC11(iCell) * SOLUBILITY_CFC(sst, sss, 11)

          CFC11_flux_pv(iCell) = xkw_ice * sqrt(660.0_RKIND / CFC11_flux_schmidt(iCell))
          CFC11_flux_atm_press(iCell) = atmosphericPressure(iCell)*PascalsToAtmospheres
          CFC11_flux_sat(iCell) = CFC11_flux_atm_press(iCell) * saturation_CFC11_1atm
          CFC11_flux_surf(iCell) = max(0.0_RKIND, CFCTracers(indexCFC11,iLevelSurface,iCell))
          CFCSurfaceFlux(indexCFC11,iCell) = CFC11_flux_pv(iCell) *   &
             (CFC11_flux_sat(iCell) - CFC11_flux_surf(iCell))
          CFC11_flux_xkw(iCell) = xkw_ice
          CFC11_flux_ifrac(iCell) = CFC_flux_ifrac_clamped
          CFC11_flux_ws(iCell) = sqrt(windSpeedSquared10m(iCell))

       endif  !  config_use_CFC11

       if (config_use_CFC12) then

          CFC12_flux_schmidt(iCell) = SCHMIDT_CFC(sst, 12)

          saturation_CFC12_1atm = pCFC12(iCell) * SOLUBILITY_CFC(sst, sss, 12)

          CFC12_flux_pv(iCell) = xkw_ice * sqrt(660.0_RKIND / CFC12_flux_schmidt(iCell))
          CFC12_flux_atm_press(iCell) = atmosphericPressure(iCell)*PascalsToAtmospheres
          CFC12_flux_sat(iCell) = CFC12_flux_atm_press(iCell) * saturation_CFC12_1atm
          CFC12_flux_surf(iCell) = max(0.0_RKIND, CFCTracers(indexCFC12,iLevelSurface,iCell))
          CFCSurfaceFlux(indexCFC12,iCell) = CFC12_flux_pv(iCell) *   &
             (CFC12_flux_sat(iCell) - CFC12_flux_surf(iCell))
          CFC12_flux_xkw(iCell) = xkw_ice
          CFC12_flux_ifrac(iCell) = CFC_flux_ifrac_clamped
          CFC12_flux_ws(iCell) = sqrt(windSpeedSquared10m(iCell))

       endif  !  config_use_CFC12

      enddo  !  iCell

      call mpas_timer_stop("CFC surface flux")

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_CFC_surface_flux_compute!}}}

!***********************************************************************
!
!  routine ocn_tracer_CFC_init
!
!> \brief   Initializes ocean CFC11 and CFC12
!> \author  Mathew Maltrud
!> \date    12/28/2016
!> \details
!>  This routine initializes fields required for CFC tracers
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_CFC_init(domain,err)!{{{

!NOTE:  called from mpas_ocn_forward_mode.F

      type (domain_type), intent(inout) :: domain !< Input/Output: domain information

      integer, intent(out) :: err !< Output: error flag

      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: tracersPool
      type (mpas_pool_type), pointer :: forcingPool
      type (mpas_pool_type), pointer :: CFCAuxiliary
      type (mpas_pool_type), pointer :: CFCAnnualForcing

      ! three dimensional pointers
      real (kind=RKIND), dimension(:,:,:), pointer :: &
        CFCTracers

! input flux components in CFCAuxiliary
      real (kind=RKIND), dimension(:), pointer :: &
         pCFC11,    &
         pCFC12

! input flux components in CFCAnnualForcing
      real (kind=RKIND), dimension(:), pointer :: &
         atmCFC11,    &
         atmCFC12

      character(len=strKIND) :: &
         forcingIntervalAnnual,  &
         forcingReferenceTimeAnnual

      type (MPAS_Time_Type) :: currTime
      integer :: ierr
      character(len=StrKIND)  :: timeStamp

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      err = 0

      if (.not. config_use_CFCTracers) return

      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)
      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
      call mpas_pool_get_array(tracersPool, 'CFCTracers', CFCTracers, 1)

!-----------------------------------------------------------------------
!     schmidt number coefficients for CFC 11
!-----------------------------------------------------------------------
      a1schmidt(11) = 3501.8_RKIND
      a2schmidt(11) = -210.31_RKIND
      a3schmidt(11) =    6.1851_RKIND
      a4schmidt(11) =   -0.07513_RKIND

!-----------------------------------------------------------------------
!     schmidt number coefficients for CFC 12
!-----------------------------------------------------------------------
      a1schmidt(12) = 3845.4_RKIND
      a2schmidt(12) = -228.95_RKIND
      a3schmidt(12) =    6.1908_RKIND
      a4schmidt(12) =   -0.067430_RKIND

!-----------------------------------------------------------------------
!     solubility coefficients for CFC 11
!-----------------------------------------------------------------------
      a1sol ( 11) = -229.9261_RKIND
      a2sol ( 11) =  319.6552_RKIND
      a3sol ( 11) =  119.4471_RKIND
      a4sol ( 11) =   -1.39165_RKIND
      b1sol ( 11) =   -0.142382_RKIND
      b2sol ( 11) =    0.091459_RKIND
      b3sol ( 11) =   -0.0157274_RKIND

!-----------------------------------------------------------------------
!     solubility coefficients for CFC 12
!-----------------------------------------------------------------------
      a1sol ( 12) = -218.0971_RKIND
      a2sol ( 12) =  298.9702_RKIND
      a3sol ( 12) =  113.8049_RKIND
      a4sol ( 12) =   -1.39165_RKIND
      b1sol ( 12) =   -0.143566_RKIND
      b2sol ( 12) =    0.091015_RKIND
      b3sol ( 12) =   -0.0153924_RKIND

!-----------------------------------------------------------------------
! initialize annual forcing to be read from file
!-----------------------------------------------------------------------

      forcingIntervalAnnual = "0001-00-00_00:00:00"
      forcingReferenceTimeAnnual = "1936-07-01_00:00:00"

      currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)
      call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)

      call MPAS_forcing_init_group( forcingGroupHead,  &
                "CFCAnnualAtmosphericForcing", &
                domain, &
                trim(timeStamp), &
                '1936-07-01_00:00:00', &
                '0080-00-00_00:00:00', &
                config_do_restart)

      call mpas_pool_get_subpool(domain % blocklist % structs, 'forcing', forcingPool)
      call mpas_pool_get_subpool(forcingPool, 'CFCAuxiliary', CFCAuxiliary)
      if (config_use_CFC11) call mpas_pool_get_array(CFCAuxiliary, 'pCFC11', pCFC11)
      if (config_use_CFC12) call mpas_pool_get_array(CFCAuxiliary, 'pCFC12', pCFC12)

      call mpas_pool_get_subpool(domain % blocklist % structs, 'CFCAnnualForcing', CFCAnnualForcing)
      if (config_use_CFC11) call mpas_pool_get_array(CFCAnnualForcing, 'atmCFC11', atmCFC11)
      if (config_use_CFC12) call mpas_pool_get_array(CFCAnnualForcing, 'atmCFC12', atmCFC12)

      if (config_use_CFC11) call MPAS_forcing_init_field( domain % streamManager, &
                                      forcingGroupHead, &
                                      'CFCAnnualAtmosphericForcing', &
                                      'atmCFC11', &
                                      'CFC_annual_mole_fraction', &
                                      'CFCAnnualForcing',  &
                                      'atmCFC11',  &
                                      'linear',  &
                                      forcingReferenceTimeAnnual,  &
                                      forcingIntervalAnnual)

      if (config_use_CFC12) call MPAS_forcing_init_field( domain % streamManager, &
                                      forcingGroupHead, &
                                      'CFCAnnualAtmosphericForcing', &
                                      'atmCFC12', &
                                      'CFC_annual_mole_fraction', &
                                      'CFCAnnualForcing',  &
                                      'atmCFC12',  &
                                      'linear',  &
                                      forcingReferenceTimeAnnual,  &
                                      forcingIntervalAnnual)

      call MPAS_forcing_init_field_data( forcingGroupHead, &
        'CFCAnnualAtmosphericForcing',  &
         domain % streamManager,  &
         config_do_restart, &
         .false.)

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_CFC_init!}}}

!***********************************************************************

!***********************************************************************
!
!  function SCHMIDT_CFC(Temperature,CFCmolecule)
!
!> \brief   computes the Schmidt Number for CFC11 and CFC 12
!> \author  Mathew Maltrud
!> \date    12/28/2016
!> \details
!>  CFC 11 and 12 schmidt number as a function of temperature.
!>  ref: Zheng et al (1998), JGR, vol 103,No C1
!>  Temperature: temperature (degree Celcius)
!>  CFCmolecule: = 11 for CFC-11,  12 for CFC-12
!>  J-C Dutay - LSCE
!-----------------------------------------------------------------------

   function SCHMIDT_CFC(Temperature,CFCmolecule)

      implicit none

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), intent(in) :: Temperature
      integer, intent(in) :: CFCmolecule

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND) :: SCHMIDT_CFC

      SCHMIDT_CFC = a1schmidt(CFCmolecule) + a2schmidt(CFCmolecule) * Temperature + &
                    a3schmidt(CFCmolecule) *Temperature*Temperature +  &
                    a4schmidt(CFCmolecule) *Temperature*Temperature*Temperature

   end function SCHMIDT_CFC

!***********************************************************************
!
!  function SOLUBILITY_CFC(Temperature,Salinity,CFCmolecule)
!
!> \brief   computes CFC 11 and 12 Solubilities in seawater
!> \author  Mathew Maltrud
!> \date    12/28/2016
!> \details
!>    computes CFC 11 and 12 Solubilities in seawater
!>    ref: Warner & Weiss (1985) , Deep Sea Research, vol32
!>    Temperature:       temperature (degree Celcius)
!>    Salinity   :       salinity    (o/oo)
!>    CFCmolecule:       11 = CFC-11, 12 = CFC-12
!>    SOLUBILITY_CFC:  in mol/m3/pptv
!>              1 pptv = 1 part per trillion = 10^-12 atm = 1 picoatm
!>    J-C Dutay - LSCE
!
!-----------------------------------------------------------------------

   function SOLUBILITY_CFC(Temperature,Salinity,CFCmolecule)

      implicit none

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), intent(in) :: Temperature, Salinity
      integer, intent(in) :: CFCmolecule

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND) :: p01SSTKelvin
      real (kind=RKIND) :: SOLUBILITY_CFC

      p01SSTKelvin = (Temperature + T0_Kelvin)* 0.01_RKIND

      SOLUBILITY_CFC  &
          = exp ( a1sol ( CFCmolecule)   &
          +       a2sol ( CFCmolecule)/ p01SSTKelvin   &
          +       a3sol ( CFCmolecule)* log ( p01SSTKelvin )   &
          +       a4sol ( CFCmolecule)* p01SSTKelvin * p01SSTKelvin   &
          +       Salinity*(   &
                    ( b3sol( CFCmolecule)*p01SSTKelvin + b2sol( CFCmolecule) )*p01SSTKelvin + b1sol(CFCmolecule) )  &
                   )

!     conversion from mol/(l * atm) to mol/(m^3 * atm)
!     ------------------------------------------------
      SOLUBILITY_CFC = SOLUBILITY_CFC * 1000.0_RKIND

!     conversion from mol/(m^3 * atm) to mol/(m3 * pptv)
!     --------------------------------------------------
      SOLUBILITY_CFC = SOLUBILITY_CFC * 1.0e-12_RKIND

   end function SOLUBILITY_CFC

!***********************************************************************

!***********************************************************************
!
!  routine get_CFCData
!
!> \brief   retrieve data needed to compute CFC surface fluxes
!> \author  Mathew Maltrud
!> \date    08/31/21
!> \details
!>  This routine calls mpas_forcing routines to acquire needed CFC forcing data and interpolates
!>    between time levels.  directly copied from ocn_get_ecosysData.
!
!-----------------------------------------------------------------------

    subroutine ocn_get_CFCData( streamManager, &
        domain, &
        simulationClock, &
        firstTimeStep) !{{{

        type (MPAS_streamManager_type), intent(inout) :: streamManager

        type (domain_type) :: domain
        type (MPAS_timeInterval_type) :: timeStepCFC
        type (MPAS_clock_type) :: simulationClock

        logical, intent(in) :: firstTimeStep
        real(kind=RKIND) :: dt

        type (mpas_pool_type), pointer :: forcingPool
        type (mpas_pool_type), pointer :: meshPool
        type (mpas_pool_type), pointer :: CFCAuxiliary
        type (mpas_pool_type), pointer :: CFCAnnualForcing

        real (kind=RKIND), dimension(:), pointer :: &
         pCFC11,         &
         pCFC12

! input flux components in CFCAnnualForcing
        real (kind=RKIND), dimension(:), pointer :: &
         atmCFC11,         &
         atmCFC12

        integer, pointer :: nCells
        integer :: iCell

        call mpas_set_timeInterval(timeStepCFC,timeString=config_dt)
        call mpas_get_timeInterval(timeStepCFC,dt=dt)

        if (firstTimeStep .and. config_do_restart) then
          call MPAS_forcing_get_forcing(forcingGroupHead, &
             'CFCAnnualAtmosphericForcing', streamManager, 0.0_RKIND)
        else
          call MPAS_forcing_get_forcing(forcingGroupHead, &
             'CFCAnnualAtmosphericForcing', streamManager, dt)
        endif

        call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(domain % blocklist % structs, 'forcing', forcingPool)
        call mpas_pool_get_subpool(forcingPool, 'CFCAuxiliary', CFCAuxiliary)
        call mpas_pool_get_subpool(domain % blocklist % structs, 'CFCAnnualForcing', CFCAnnualForcing)

        call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

        if (config_use_CFC11) then
          call mpas_pool_get_array(CFCAuxiliary, 'pCFC11', pCFC11)
          call mpas_pool_get_array(CFCAnnualForcing, 'atmCFC11', atmCFC11)
          do iCell = 1, nCells
             pCFC11(iCell) = atmCFC11(iCell)
          enddo
        endif

        if (config_use_CFC12) then
          call mpas_pool_get_array(CFCAuxiliary, 'pCFC12', pCFC12)
          call mpas_pool_get_array(CFCAnnualForcing, 'atmCFC12', atmCFC12)
          do iCell = 1, nCells
             pCFC12(iCell) = atmCFC12(iCell)
          enddo
        endif

    end subroutine ocn_get_CFCData!}}}

!***********************************************************************
!
!  routine ocn_CFC_forcing_write_restart
!
!> \brief   writes restart timestamp for CFC data to be read in on future restart
!> \author  Mathew Maltrud
!> \date    03/07/2016

!
!-----------------------------------------------------------------------

   subroutine ocn_CFC_forcing_write_restart(domain)!{{{

      type(domain_type) :: domain

      call MPAS_forcing_write_restart_times(forcingGroupHead)

    end subroutine ocn_CFC_forcing_write_restart!}}}

end module ocn_tracer_CFC

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
