!---------------------------------------------------------------------------
! $Id: array_index.F90 7118 2014-07-25 00:12:15Z raut@uwm.edu $
!===============================================================================
module array_index

  ! Description:
  ! Contains indices to variables in larger arrays.
  ! Note that the 'ii' is necessary because 'i' is used in
  ! statistics to track locations in the zt/zm/sfc derived types.

  ! References:
  !   None
  !-------------------------------------------------------------------------

  use clubb_precision, only: &
    core_rknd      ! Precision

  implicit none

  ! Variables
  ! Microphysics mixing ratios
  integer, public :: &
    iirrm,    & ! Hydrometeor array index for rain water mixing ratio, rr
    iirsm,    & ! Hydrometeor array index for snow mixing ratio, rs
    iirim,     & ! Hydrometeor array index for ice mixing ratio, ri
    iirgm    ! Hydrometeor array index for graupel mixing ratio, rg
!$omp threadprivate(iirrm, iirsm, iirim, iirgm)

  ! Microphysics concentrations
  integer, public :: &
    iiNrm,       & ! Hydrometeor array index for rain drop concentration, Nr
    iiNsm,    & ! Hydrometeor array index for snow concentration, Ns
    iiNim,       & ! Hydrometeor array index for ice concentration, Ni
    iiNgm    ! Hydrometeor array index for graupel concentration, Ng
!$omp threadprivate(iiNrm, iiNsm, iiNim, iiNgm)

  ! Scalar quantities
  integer, public :: & 
    iisclr_rt, iisclr_thl, iisclr_CO2, & ! [kg/kg]/[K]/[1e6 mol/mol]
    iiedsclr_rt, iiedsclr_thl, iiedsclr_CO2 ! "    "
!$omp threadprivate(iisclr_rt, iisclr_thl, iisclr_CO2, &
!$omp   iiedsclr_rt, iiedsclr_thl, iiedsclr_CO2)

  ! Logical fields
  logical, dimension(:), allocatable, public :: &
    l_frozen_hm, & ! if true, then the hydrometeor is frozen; otherwise liquid
    l_mix_rat_hm   ! if true, then the quantity is a hydrometeor mixing ratio
!$omp threadprivate(l_frozen_hm, l_mix_rat_hm)

  character(len=10), dimension(:), allocatable, public :: & 
    hydromet_list

!$omp threadprivate( hydromet_list )

  real( kind = core_rknd ), dimension(:), allocatable, public :: &
    hydromet_tol    ! Tolerance values for all hydrometeors    [units vary]

!$omp threadprivate( hydromet_tol )   

  private ! Default Scope

!===============================================================================

end module array_index
