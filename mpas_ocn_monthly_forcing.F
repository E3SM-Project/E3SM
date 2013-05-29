! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_monthly_forcing
!
!> \brief MPAS ocean monthly forcing
!> \author Doug Jacobsen
!> \date   04/25/12
!> \details
!>  This module contains routines for building the forcing arrays,
!>  if monthly forcing is used.
!
!-----------------------------------------------------------------------

module ocn_monthly_forcing

   use mpas_grid_types
   use mpas_configure
   use mpas_timekeeping

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

   public :: ocn_build_forcing_arrays, &
             ocn_monthly_forcing_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: monthlyForcingOn !< Flag to turn on/off resotring

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_build_forcing_arrays
!
!> \brief   Determines the forcing array used for the monthly forcing.
!> \author  Doug Jacobsen
!> \date    04/25/12
!> \details 
!>  This routine computes the forcing arrays used later in MPAS.
!
!-----------------------------------------------------------------------

   subroutine ocn_build_forcing_arrays(timeStamp, grid, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type(MPAS_Time_type), intent(in) :: timeStamp

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mesh_type), intent(inout) :: &
         grid          !< Input: grid information


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


      real (kind=RKIND), dimension(:,:), pointer :: temperatureRestoreMonthly 
      real (kind=RKIND), dimension(:,:), pointer :: salinityRestoreMonthly 
      real (kind=RKIND), dimension(:,:), pointer :: windStressMonthly 
      real (kind=RKIND), dimension(:), pointer :: temperatureRestore 
      real (kind=RKIND), dimension(:), pointer :: salinityRestore 
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocityForcing
      integer :: iCell, iEdge, nCells, nEdges, nMonths, k
      integer :: iMonth, iMonthP1, iDayInMonth, ierr
      real (kind=RKIND) :: data, dataP1, weight, weightP1

      err = 0

      if(.not.monthlyForcingOn) return

      nCells = grid % nCells
      nEdges = grid % nEdges
      nMonths = grid % nMonths

      temperatureRestore => grid % temperatureRestore % array
      salinityRestore => grid % salinityRestore % array
      normalVelocityForcing => grid % normalVelocityForcing % array

      temperatureRestoreMonthly => grid % temperatureRestoreMonthly % array
      salinityRestoreMonthly => grid % salinityRestoreMonthly % array
      windStressMonthly => grid % windStressMonthly % array

      call mpas_get_time(timeStamp, MM = iMonth, DD = iDayInMonth, ierr = ierr)

      err = ierr

      iMonthP1 = mod(iMonth, nMonths) + 1

      weight = 1.0 - (iDayInMonth-1) / 30.0
      weightP1 = 1.0 - weight

      do iCell=1,nCells
        ! Interpolate between iMonth and iMonthP1 records, using iDayInMonth
        data = temperatureRestoreMonthly(iMonth,iCell)
        dataP1 = temperatureRestoreMonthly(iMonthP1,iCell)
        temperatureRestore(iCell) = data * weight + dataP1 * weightP1
        data = salinityRestoreMonthly(iMonth,iCell)
        dataP1 = salinityRestoreMonthly(iMonthP1,iCell)
        salinityRestore(iCell) = data * weight + dataP1 * weightP1
      end do

      do iEdge=1,nEdges
        ! Interpolate between iMonth and iMonthP1 records, using iDayInMonth
        data = windStressMonthly(iMonth,iEdge)
        dataP1 = windStressMonthly(iMonthP1,iEdge)
        normalVelocityForcing(1,iEdge) = data * weight + dataP1 * weightP1
      end do

   !--------------------------------------------------------------------

   end subroutine ocn_build_forcing_arrays!}}}

!***********************************************************************
!
!  routine ocn_monthly_forcing_init
!
!> \brief   Initializes monthly forcing module
!> \author  Doug Jacobsen
!> \date    04/25/12
!> \details 
!>  This routine initializes the monthly forcing module.
!
!-----------------------------------------------------------------------

   subroutine ocn_monthly_forcing_init(err)!{{{

      integer, intent(out) :: err !< Output: error flag

      err = 0

      monthlyForcingOn = .false.

      if(config_use_monthly_forcing) then
        monthlyForcingOn = .true.

        write (0,'(a)') " Monthly forcing is on.  Make sure monthly forcing variables include iro in Registry, and are in your initial condition or restart file."
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_monthly_forcing_init!}}}

!***********************************************************************

end module ocn_monthly_forcing

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
