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

  type t_eatm_interpolator(kind)
    integer, kind :: kind
    real(kind=kind), dimension(:, :, :), allocatable :: t_im1
    real(kind=kind), dimension(:, :, :), allocatable :: t_ip1
  end type t_eatm_interpolator

  type(t_eatm_interpolator(kind=R4)), public :: eatm_intrp

end module eatmMod
