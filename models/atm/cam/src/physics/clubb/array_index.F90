!-----------------------------------------------------------------------
! $Id: array_index.F90 5216 2011-06-06 18:58:41Z dschanen@uwm.edu $
!-----------------------------------------------------------------------
module array_index

! Description:
!   Contains indices to variables in larger arrays.
!   Note that the 'ii' is necessary because 'i' is used in
!   statistics to track locations in the zt/zm/sfc derived types.

! References:
!   None
!-----------------------------------------------------------------------
  implicit none

  ! Variables
  ! Microphysics mixing ratios
  integer, public :: &
    iirrainm, iirsnowm, iiricem, iirgraupelm ! [kg/kg]
!$omp threadprivate(iirrainm, iirsnowm, iiricem, iirgraupelm)

  ! Microphysics number concentration
  integer, public :: &
    iiNrm, iiNsnowm, iiNim, iiNgraupelm, iiNcm ! [#/kg]
!$omp threadprivate(iiNrm, iiNsnowm, iiNim, iiNgraupelm, iiNcm)

  ! Scalar quantities
  integer, public :: & 
    iisclr_rt, iisclr_thl, iisclr_CO2, & ! [kg/kg]/[K]/[1e6 mol/mol]
    iiedsclr_rt, iiedsclr_thl, iiedsclr_CO2 ! "    "
!$omp threadprivate(iisclr_rt, iisclr_thl, iisclr_CO2, &
!$omp   iiedsclr_rt, iiedsclr_thl, iiedsclr_CO2)

  private ! Default Scope

end module array_index
!-----------------------------------------------------------------------
