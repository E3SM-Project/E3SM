MODULE rrtm_params
  !-----------------------------------------------------------------------
  ! Define parameters for the RRTMG radiation code.
  !-----------------------------------------------------------------------

  USE parkind,   only: kind_rm
  USE rrtm_grid, only: nx, ny

  REAL (KIND=kind_rm), PARAMETER :: &
       cp = 1004.5_kind_rm, &      ! specific heat of dry air (J/kg/deg)
       ggr = 9.81_kind_rm,  &      ! gravitational constant (m/s2)
       secday = 86400.0_kind_rm    ! seconds in one day

  REAL (KIND=kind_rm) :: coszrs        ! cosine solar zenith angle

  REAL (KIND=kind_rm), allocatable, DIMENSION(:,:) :: &
       latitude, &   ! latitude (degrees)
       longitude     ! longitude (degrees)

  public :: allocate_rrtm_params

CONTAINS
  subroutine allocate_rrtm_params()
    allocate(latitude(nx, ny))
    allocate(longitude(nx, ny))
  end subroutine allocate_rrtm_params

END MODULE rrtm_params
