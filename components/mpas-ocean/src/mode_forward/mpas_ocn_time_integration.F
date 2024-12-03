! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integration
!
!> \brief MPAS ocean time integration driver
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011
!> \details
!>  This module contains the main driver routine for calling
!>  the time integration scheme
!
!-----------------------------------------------------------------------

module ocn_time_integration

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_timekeeping
   use mpas_timer
   use mpas_log

   use ocn_constants
   use ocn_config
   use ocn_diagnostics_variables
   use ocn_time_integration_rk4
   use ocn_time_integration_split
   use ocn_time_integration_si
   use ocn_time_integration_lts
   use ocn_time_integration_fblts
   use ocn_time_integration_split_ab2


   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

    ! Enum for selecting different time integrators
    integer, public :: timeIntegratorChoice

    integer, public, parameter :: &
       timeIntUnknown          = 0, &! unknown or undefined
       timeIntSplitExplicit    = 1, &! split-explicit
       timeIntUnsplitExplicit  = 2, &! unsplit-explicit
       timeIntSemiImplicit     = 3, &! Semi-implicit
       timeIntRK4              = 4, &! 4th-order Runge-Kutta
       timeIntLTS              = 5, &! local time-stepping
       timeIntFBLTS            = 6, &! forward-backward lts
       timeIntSplitExplicitAB2 = 7   ! split-explicit AB2 baroclinic

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_timestep, &
             ocn_timestep_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------


!***********************************************************************

   contains

!***********************************************************************
!
!  ocn_timestep
!
!> \brief MPAS ocean time integration driver
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011
!> \details
!>  This routine handles a single timestep for the ocean and advances
!>  the ocean by one time step dt. It determines the time integrator
!>  that will be used for the run, and calls the appropriate one.
!
!-----------------------------------------------------------------------

   subroutine ocn_timestep(domain, dt, timeStamp)!{{{

      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), intent(in) :: &
         dt              !< [in] time step in seconds to advance model

      character(len=*), intent(in) :: &
         timeStamp       !< [in] current time step in simulation

      !-----------------------------------------------------------------
      ! Input/Output variables
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: &
         domain !< [inout] model state to advance forward by dt seconds

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------

      type (MPAS_Time_type) :: &
         xtime_timeType,               &! current time in time type
         simulationStartTime_timeType, &! start time in time type
         referenceTime_timeType         ! reference time for Time coordinate

      ! End preamble
      !-------------
      ! Begin code

      ! Advance the model base on choice of time integration scheme
      select case (timeIntegratorChoice)
      case (timeIntSplitExplicit)
         call ocn_time_integrator_split(domain, dt)
      case (timeIntUnsplitExplicit)
         call ocn_time_integrator_split(domain, dt)
      case (timeIntSemiImplicit)
         call ocn_time_integrator_si(domain, dt)
      case (timeIntRK4)
         call ocn_time_integrator_rk4(domain, dt)
      case (timeIntLTS)
         call ocn_time_integrator_lts(domain, dt)
      case (timeIntFBLTS)
         call ocn_time_integrator_fblts(domain, dt)
      case (timeIntSplitExplicitAB2)
         call ocn_time_integrator_split_ab2(domain, dt)
      end select

      ! compute time since start of simulation, in days
      ! time variables stored in diagnostic variables
      xtime = timeStamp  ! current time (input)

      call mpas_set_time(xtime_timeType, dateTimeString=xtime)
      call mpas_set_time(simulationStartTime_timeType, &
                         dateTimeString=simulationStartTime)
      call mpas_get_timeInterval( &
           xtime_timeType - simulationStartTime_timeType, &
           dt=daysSinceStartOfSim)

      daysSinceStartOfSim = daysSinceStartOfSim*days_per_second

      call mpas_set_time(referenceTime_timeType, &
                         dateTimeString=config_output_reference_time)
      call mpas_get_timeInterval( &
           xtime_timeType - referenceTime_timeType, &
           dt=Time)

      Time = Time*days_per_second
      ! bounds on instantaneous time are just the time itself
      Time_bnds(1) = Time
      Time_bnds(2) = Time

   !--------------------------------------------------------------------

   end subroutine ocn_timestep!}}}

!***********************************************************************
!
!  ocn_timestep_init
!
!> \brief MPAS ocean time integration initialization driver
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011
!> \details
!>  This routine initializes variables need for the time integration
!>  scheme. It calls the relevant init routine base on choice of scheme.
!
!-----------------------------------------------------------------------

   subroutine ocn_timestep_init(domain, dt, err)!{{{

      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), intent(in) :: &
         dt              !< [in] time step (sec)

      !-----------------------------------------------------------------
      ! Input/Output variables
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: &
         domain          !< [inout] model state

      !-----------------------------------------------------------------
      ! Output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err ! returned error flag

      ! End preamble
      !-------------
      ! Begin code

      err = 0 ! initialize return code to success (0)

      ! Set enum for choice of time integrator to avoid character
      ! comparisons later and call relevant init routine

      select case (trim(config_time_integrator))
      case ('split_explicit')
         timeIntegratorChoice = timeIntSplitExplicit
         call ocn_time_integration_split_init(domain)

      case ('unsplit_explicit')
         timeIntegratorChoice = timeIntUnsplitExplicit
         call ocn_time_integration_split_init(domain)

      case ('split_implicit')
         timeIntegratorChoice = timeIntSemiImplicit
         call ocn_time_integration_si_init(domain, dt)

      case ('RK4')
         timeIntegratorChoice = timeIntRK4
         call ocn_time_integration_rk4_init(domain)

      case ('LTS')
         timeIntegratorChoice = timeIntLTS
         call ocn_time_integration_lts_init(domain)
      
      case ('FB_LTS')
         timeIntegratorChoice = timeIntFBLTS
         call ocn_time_integration_fblts_init(domain)

      case ('split_explicit_ab2')
         timeIntegratorChoice = timeIntSplitExplicitAB2
         call ocn_time_integration_split_ab2_init(domain)

      case default

         err = 1
         timeIntegratorChoice = timeIntUnknown
         call mpas_log_write( &
            'Incorrect choice for config_time_integrator:' // &
             trim(config_time_integrator) // &
             '   choices are: RK4, split_explicit, ' // &
                 'unsplit_explicit, split_implicit, LTS, ' // &
                 'split_explicit_ab2', MPAS_LOG_ERR)

      end select

   !--------------------------------------------------------------------

   end subroutine ocn_timestep_init!}}}

!***********************************************************************

end module ocn_time_integration

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
