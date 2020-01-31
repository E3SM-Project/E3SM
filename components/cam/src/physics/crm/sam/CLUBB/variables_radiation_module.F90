!---------------------------------------------------------------
! $Id: variables_radiation_module.F90 5982 2012-11-21 19:20:12Z raut@uwm.edu $
module variables_radiation_module

!   This module contains definitions of all radiation arrays
!   used in the single column model, as well as subroutines to
!   allocate, deallocate, and initialize them.
!---------------------------------------------------------------

  use clubb_precision, only: &
    core_rknd ! Variable(s)

  implicit none


  public :: &
    setup_radiation_variables, &
    cleanup_radiation_variables

  private ! Set Default Scoping

  integer, private, parameter :: dp = selected_real_kind( p=12 )

  real( kind = core_rknd ), public, dimension(:), allocatable :: &
    radht_LW, & ! LW heating rate   [K/s]
    radht_SW, & ! SW heating rate   [K/s]
    Frad_SW,  & ! SW radiative flux [W/m^2]
    Frad_LW     ! LW radiative flux [W/m^2]

!$omp threadprivate(radht_LW, radht_SW, Frad_SW, Frad_LW)

  real(kind = dp), public, dimension(:,:), allocatable :: &
    T_in_K,   & ! Temperature        [K]
    rcil,     & ! Ice mixing ratio   [kg/kg]
    o3l         ! Ozone mixing ratio [kg/kg]

!$omp threadprivate(T_in_K, rcil, o3l)

  real(kind = dp), public, dimension(:,:), allocatable :: &
    rsnowm_2d,& ! Two-dimensional copies of the input parameters
    rcm_in_cloud_2d, &
    cloud_frac_2d, &
    ice_supersat_frac_2d

!$omp threadprivate(rsnowm_2d, rcm_in_cloud_2d, cloud_frac_2d)

  real(kind = dp), public, dimension(:,:), allocatable :: &
    radht_SW_2d, & ! SW Radiative heating rate  [W/m^2]
    radht_LW_2d    ! LW Radiative heating rate  [W/m^2]

!$omp threadprivate(radht_SW_2d, radht_LW_2d)

  real(kind = dp), public, dimension(:,:), allocatable :: &
    Frad_uLW, & ! LW upwelling flux         [W/m^2]
    Frad_dLW, & ! LW downwelling flux       [W/m^2]
    Frad_uSW, & ! SW upwelling flux         [W/m^2]
    Frad_dSW    ! SW downwelling flux       [W/m^2]

!$omp threadprivate(Frad_uLW, Frad_dLW, Frad_uSW, Frad_dSW)

  real(kind = dp), public, dimension(:,:), allocatable :: &
     fdswcl, &  !Downward clear-sky SW flux                 (W/m^-2).
     fuswcl, &  !Upward clear-sky SW flux                   (W/m^-2).
     fdlwcl, &  !Downward clear-sky LW flux                 (W/m^-2).
     fulwcl     !Upward clear-sky LW flux                   (W/m^-2).

!$omp threadprivate(fdswcl, fuswcl, fdlwcl, fulwcl)

  ! Constant parameters
  integer, private, parameter :: &
    nlen = 1, &   ! Length of the total domain
    slen = 1      ! Length of the sub domain

  contains

  !---------------------------------------------------------------------
  subroutine setup_radiation_variables( nzmax, lin_int_buffer, &
                                        extend_atmos_range_size )
  ! Description:
  !   Allocates and initializes prognostic scalar and array variables
  !   for the CLUBB model code.
  !---------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nzmax, & ! Number of grid levels [-]
      lin_int_buffer,& ! Number of interpolated levels between the computational
                       ! grid and the extended atmosphere [-]
      extend_atmos_range_size ! The number of levels in the extended atmosphere [-]

    ! Local Variables

    integer :: rad_zt_dim, rad_zm_dim ! Dimensions of the radiation grid

    !----------------------------BEGIN CODE-------------------------------

    rad_zt_dim = (nzmax-1)+lin_int_buffer+extend_atmos_range_size
    rad_zm_dim = (nzmax-1)+lin_int_buffer+extend_atmos_range_size+1


    ! --- Allocation ---

    allocate( radht_SW(1:nzmax) )
    allocate( radht_LW(1:nzmax) )
    allocate( Frad_SW(1:nzmax) )
    allocate( Frad_LW(1:nzmax) )

    allocate( T_in_K(nlen, rad_zt_dim ) )
    allocate( rcil(nlen, rad_zt_dim ) )
    allocate( o3l(nlen, rad_zt_dim ) )

    allocate( rsnowm_2d(nlen, rad_zt_dim ) )
    allocate( rcm_in_cloud_2d(nlen, rad_zt_dim ) )
    allocate( cloud_frac_2d(nlen, rad_zt_dim ) )
    allocate( ice_supersat_frac_2d(nlen, rad_zt_dim ) )

    allocate( radht_SW_2d(nlen, rad_zt_dim ) )
    allocate( radht_LW_2d(nlen, rad_zt_dim ) )

    allocate( Frad_uLW(nlen, rad_zm_dim ) )
    allocate( Frad_dLW(nlen, rad_zm_dim ) )
    allocate( Frad_uSW(nlen, rad_zm_dim ) )
    allocate( Frad_dSW(nlen, rad_zm_dim ) )

    allocate( fdswcl(slen, rad_zm_dim ) )
    allocate( fuswcl(slen, rad_zm_dim ) )
    allocate( fdlwcl(slen, rad_zm_dim ) )
    allocate( fulwcl(slen, rad_zm_dim ) )


    ! --- Initialization ---

    radht_SW = 0.0_core_rknd
    radht_LW = 0.0_core_rknd
    Frad_SW = 0.0_core_rknd
    Frad_LW = 0.0_core_rknd
    T_in_K = 0.0_dp
    rcil = 0.0_dp
    o3l = 0.0_dp
    rsnowm_2d = 0.0_dp
    rcm_in_cloud_2d = 0.0_dp
    cloud_frac_2d = 0.0_dp
    ice_supersat_frac_2d = 0.0_dp
    radht_SW_2d = 0.0_dp
    radht_LW_2d = 0.0_dp
    Frad_uLW = 0.0_dp
    Frad_dLW = 0.0_dp
    Frad_uSW = 0.0_dp
    Frad_dSW = 0.0_dp
    fdswcl = 0.0_dp
    fuswcl = 0.0_dp
    fdlwcl = 0.0_dp
    fulwcl = 0.0_dp

  end subroutine setup_radiation_variables

  !---------------------------------------------------------------------
  subroutine cleanup_radiation_variables( )
  
  ! Description:
  !   Subroutine to deallocate variables defined in module global
  !---------------------------------------------------------------------
 
    implicit none

    ! --- Deallocate ---

    deallocate( radht_SW )
    deallocate( radht_LW )
    deallocate( Frad_SW )
    deallocate( Frad_LW )

    deallocate( T_in_K )
    deallocate( rcil )
    deallocate( o3l )

    deallocate( rsnowm_2d )
    deallocate( rcm_in_cloud_2d )
    deallocate( cloud_frac_2d )
    deallocate( ice_supersat_frac_2d )

    deallocate( radht_SW_2d )
    deallocate( radht_LW_2d )

    deallocate( Frad_uLW )
    deallocate( Frad_dLW )
    deallocate( Frad_uSW )
    deallocate( Frad_dSW )

    deallocate( fdswcl )
    deallocate( fuswcl )
    deallocate( fdlwcl )
    deallocate( fulwcl )

  end subroutine cleanup_radiation_variables
    
    
end module variables_radiation_module
