module eatmMod

  ! !USES:

  use shr_kind_mod, only: &
    IN=>SHR_KIND_IN, &
    R4=>SHR_KIND_R4, &
    R8=>SHR_KIND_R8, &
    CS=>SHR_KIND_CS, &
    CL=>SHR_KIND_CL

  ! !PUBLIC TYPES:

  implicit none
  private ! except

  !--------------------------------------------------------------------------
  ! Public module data
  !--------------------------------------------------------------------------
  integer, public           :: gsize, lsize, lsize_x, lsize_y
  character(CL), public     :: restart_file
  character(CL), public     :: case_name      ! case name
  character(len=16), public :: inst_name
  character(len=16), public :: inst_suffix = ""    ! char string associated with instance (ie. "_0001" or "")

  ! imported arrays first
  real(kind=R8), dimension(:,:), allocatable, public :: shf          ! sensible heat flux
  real(kind=R8), dimension(:,:), allocatable, public :: cflx         ! constituent flux (emissions)
  real(kind=R8), dimension(:,:), allocatable, public :: lhf          ! latent heat flux
  real(kind=R8), dimension(:,:), allocatable, public :: wsx          ! surface u-stress (N)
  real(kind=R8), dimension(:,:), allocatable, public :: wsy          ! surface v-stress (N)
  real(kind=R8), dimension(:,:), allocatable, public :: lwup         ! longwave up radiative flux
  real(kind=R8), dimension(:,:), allocatable, public :: asdir        ! albedo: shortwave, direct
  real(kind=R8), dimension(:,:), allocatable, public :: aldir        ! albedo: longwave, direct
  real(kind=R8), dimension(:,:), allocatable, public :: asdif        ! albedo: shortwave, diffuse
  real(kind=R8), dimension(:,:), allocatable, public :: aldif        ! albedo: longwave, diffuse
  real(kind=R8), dimension(:,:), allocatable, public :: ts           ! merged surface temp
  real(kind=R8), dimension(:,:), allocatable, public :: sst          ! sea surface temp
  real(kind=R8), dimension(:,:), allocatable, public :: snowhland    ! snow depth (liquid water equivalent) over land
  real(kind=R8), dimension(:,:), allocatable, public :: snowhice     ! snow depth over ice
  real(kind=R8), dimension(:,:), allocatable, public :: tref         ! ref height surface air temp
  real(kind=R8), dimension(:,:), allocatable, public :: qref         ! ref height specific humidity
  real(kind=R8), dimension(:,:), allocatable, public :: u10          ! 10m wind speed
  real(kind=R8), dimension(:,:), allocatable, public :: u10withgusts ! 10m wind speed with gustiness
  real(kind=R8), dimension(:,:), allocatable, public :: icefrac      ! sea-ice areal fraction
  real(kind=R8), dimension(:,:), allocatable, public :: ocnfrac      ! ocean areal fraction
  real(kind=R8), dimension(:,:), allocatable, public :: lndfrac      ! land area fraction

  ! exported arrays
  real(kind=R8), dimension(:,:), allocatable, public :: zbot         ! bot level height above surface
  real(kind=R8), dimension(:,:), allocatable, public :: ubot         ! bot level u wind
  real(kind=R8), dimension(:,:), allocatable, public :: vbot         ! bot level v wind
  real(kind=R8), dimension(:,:), allocatable, public :: tbot         ! bot level temperature
  real(kind=R8), dimension(:,:), allocatable, public :: ptem         ! bot level potential temperature
  real(kind=R8), dimension(:,:), allocatable, public :: shum         ! bot level specific humidity
  real(kind=R8), dimension(:,:), allocatable, public :: dens         ! bot level density
  real(kind=R8), dimension(:,:), allocatable, public :: pbot         ! bot level pressure
  real(kind=R8), dimension(:,:), allocatable, public :: pslv         ! sea level atm pressure
  real(kind=R8), dimension(:,:), allocatable, public :: lwdn         ! Down longwave flux at surface
  real(kind=R8), dimension(:,:), allocatable, public :: rainc        ! liquid "convective" precip
  real(kind=R8), dimension(:,:), allocatable, public :: rainl        ! liquid "large scale" precip
  real(kind=R8), dimension(:,:), allocatable, public :: snowc        ! frozen "convective" precip
  real(kind=R8), dimension(:,:), allocatable, public :: snowl        ! frozen "large scale" precip
  real(kind=R8), dimension(:,:), allocatable, public :: swndr        ! direct near-infrared incident solar radiation
  real(kind=R8), dimension(:,:), allocatable, public :: swvdr        ! direct visible incident solar radiation
  real(kind=R8), dimension(:,:), allocatable, public :: swndf        ! diffuse near-infrared incident solar radiation
  real(kind=R8), dimension(:,:), allocatable, public :: swvdf        ! diffuse visible incident solar radiation
  real(kind=R8), dimension(:,:), allocatable, public :: swnet        ! net shortwave radiation

  !
  real(kind=R4), dimension(:, :, :, :), allocatable, target, public :: net_inputs
  real(kind=R4), dimension(:, :, :, :), allocatable, target, public :: net_outputs

  character(CS), public :: myModelName = 'atm'   ! user defined model name

  character(len=*), parameter, public :: rpfile = 'rpointer.atm'

  type :: t_eatm_interpolator(kind)
    integer, kind :: kind
    real(kind=kind), dimension(:, :, :), allocatable :: t_im1
    real(kind=kind), dimension(:, :, :), allocatable :: t_ip1
  end type t_eatm_interpolator

  type(t_eatm_interpolator(kind=R4)), public :: eatm_intrp

  type, public :: t_normalization_struct
    real(kind=R4), dimension(:), allocatable :: means
    real(kind=R4), dimension(:), allocatable :: stds
  end type t_normalization_struct

  type, extends(t_normalization_struct) :: t_normalizer
    contains
      procedure :: normalize
  end type t_normalizer

  type, extends(t_normalization_struct) :: t_denormalizer
    contains
      procedure :: denormalize
  end type t_denormalizer

  type(t_normalizer), public :: normalizer
  type(t_denormalizer), public :: denormalizer

contains

  subroutine normalize(self, inputs)
    class(t_normalizer) :: self
    real(kind=R4), intent(inout) :: inputs(:, :, :, :)

    integer :: i, j, k
    integer :: nx, ny, nc

    nc = SIZE(inputs, dim=2)
    nx = SIZE(inputs, dim=3)
    ny = SIZE(inputs, dim=4)

    do k = 1, nc
      do j = 1, ny
        do i = 1, nx
          inputs(1, k, i, j) = (inputs(1, k, i, j) + self%means(k)) / self%stds(k)
        enddo
      enddo
    enddo

  end subroutine normalize

  subroutine denormalize(self, outputs)
    class(t_denormalizer) :: self
    real(kind=R4), intent(inout) :: outputs(:, :, :, :)

    integer :: i, j, k
    integer :: nx, ny, nc

    nc = SIZE(outputs, dim=2)
    nx = SIZE(outputs, dim=3)
    ny = SIZE(outputs, dim=4)

    do k = 1, nc
      do j = 1, ny
        do i = 1, nx
          outputs(1, k, i, j) = outputs(1, k, i, j) * self%stds(k) + self%means(k)
        enddo
      enddo
    enddo

  end subroutine denormalize

end module eatmMod
