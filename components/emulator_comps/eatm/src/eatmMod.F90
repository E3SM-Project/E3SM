module eatmMod

  ! !USES:

  use shr_kind_mod   , only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL

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

     !JW TODO: load up all arrays into a big 3D container? lots of 2D arrays?
     !         for now, using 2D arrays with EAM naming convention
     !         and EAM sign conventions
  ! imported arrays first
  real(kind=R8), dimension(:,:), allocatable, public :: shf
  real(kind=R8), dimension(:,:), allocatable, public :: cflx
  real(kind=R8), dimension(:,:), allocatable, public :: lhf
  real(kind=R8), dimension(:,:), allocatable, public :: wsx
  real(kind=R8), dimension(:,:), allocatable, public :: wsy
  real(kind=R8), dimension(:,:), allocatable, public :: lwup
  real(kind=R8), dimension(:,:), allocatable, public :: asdir
  real(kind=R8), dimension(:,:), allocatable, public :: aldir
  real(kind=R8), dimension(:,:), allocatable, public :: asdif
  real(kind=R8), dimension(:,:), allocatable, public :: aldif
  real(kind=R8), dimension(:,:), allocatable, public :: ts
  real(kind=R8), dimension(:,:), allocatable, public :: sst
  real(kind=R8), dimension(:,:), allocatable, public :: snowhland
  real(kind=R8), dimension(:,:), allocatable, public :: snowhice
  real(kind=R8), dimension(:,:), allocatable, public :: tref
  real(kind=R8), dimension(:,:), allocatable, public :: qref
  real(kind=R8), dimension(:,:), allocatable, public :: u10
  real(kind=R8), dimension(:,:), allocatable, public :: u10withgusts
  real(kind=R8), dimension(:,:), allocatable, public :: icefrac
  real(kind=R8), dimension(:,:), allocatable, public :: ocnfrac
  real(kind=R8), dimension(:,:), allocatable, public :: lndfrac

  ! exported arrays
  real(kind=R8), dimension(:,:), allocatable, public :: zbot
  real(kind=R8), dimension(:,:), allocatable, public :: ubot
  real(kind=R8), dimension(:,:), allocatable, public :: vbot
  real(kind=R8), dimension(:,:), allocatable, public :: tbot
  real(kind=R8), dimension(:,:), allocatable, public :: thbot
  real(kind=R8), dimension(:,:), allocatable, public :: qbot
  real(kind=R8), dimension(:,:), allocatable, public :: rho
  real(kind=R8), dimension(:,:), allocatable, public :: pbot
  real(kind=R8), dimension(:,:), allocatable, public :: psl
  real(kind=R8), dimension(:,:), allocatable, public :: flwds
  real(kind=R8), dimension(:,:), allocatable, public :: rainc
  real(kind=R8), dimension(:,:), allocatable, public :: rainl
  real(kind=R8), dimension(:,:), allocatable, public :: snowc
  real(kind=R8), dimension(:,:), allocatable, public :: snowl
  real(kind=R8), dimension(:,:), allocatable, public :: soll
  real(kind=R8), dimension(:,:), allocatable, public :: sols
  real(kind=R8), dimension(:,:), allocatable, public :: solld
  real(kind=R8), dimension(:,:), allocatable, public :: solsd
  real(kind=R8), dimension(:,:), allocatable, public :: netsw

  character(CS), public :: myModelName = 'atm'   ! user defined model name

  character(len=*), parameter, public :: rpfile = 'rpointer.atm'

end module eatmMod
