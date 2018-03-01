module BgcCentCnpDecompType

!
! DESCRIPTIONS
! contains data structure for ECA-century decomposition.

! USES

  use bshr_kind_mod  , only : r8 => shr_kind_r8
implicit none

  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  type, public :: DecompCent_type
  real(r8) :: o_scalar       ! fraction by which decomposition is limited by anoxia
  real(r8) :: w_scalar       ! fraction by which decomposition is limited by h2osoi_liqure availability
  real(r8) :: t_scalar       ! fraction by which decomposition is limited by temperature
  real(r8) :: depth_scalar   ! depth dependent factor for heteorotrophic respiration

  !parameters
  real(r8) :: Q10
  real(r8) :: froz_q10
  real(r8) :: decomp_depth_efolding
  contains
    procedure, public  :: Init
    procedure, public  :: set_decompk_scalar
    procedure, private :: InitCold
    procedure, private :: InitAllocate
    procedure, private :: initPar
  end type DecompCent_type

 contains

 !------------------------------------------------------------------------
 subroutine Init(this, biogeo_con)

  use BiogeoConType, only : BiogeoCon_type
  implicit none
  class(DecompCent_type), intent(inout) :: this
  type(BiogeoCon_type),intent(in) :: biogeo_con

  call this%InitAllocate ()

  call this%InitCold ()

  call this%InitPar(biogeo_con)

 end subroutine Init
 !------------------------------------------------------------------------
 subroutine InitAllocate(this)
   !
   ! !DESCRIPTION:
   ! Initialize module data structure
   !
   ! !USES:
   use bshr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
   use betr_varcon    , only : spval  => bspval
   !
   ! !ARGUMENTS:
   class(DecompCent_type), intent(inout) :: this

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine initCold(this)
    !
    ! !USES:

    implicit none
    !
    ! !ARGUMENTS:
    class(DecompCent_type), intent(inout) :: this

    this%o_scalar = 1._r8
    this%w_scalar  = 1._r8
    this%t_scalar    = 1._r8
    this%depth_scalar  = 1._r8

  end subroutine initCold

  !-----------------------------------------------------------------------
  subroutine InitPar(this, biogeo_con)

  use BiogeoConType, only : BiogeoCon_type
  implicit none
  class(DecompCent_type) , intent(inout) :: this
  type(BiogeoCon_type)   , intent(in) :: biogeo_con

  ! set "Q10" parameter
  this%Q10 = biogeo_con%Q10
  this%froz_q10 = biogeo_con%froz_q10
  this%decomp_depth_efolding = biogeo_con%decomp_depth_efolding
  end subroutine InitPar
  !-----------------------------------------------------------------------
  subroutine set_decompk_scalar(this, o2b, centuryeca_forc)

  use BgcCentCnpForcType , only : centuryeca_forc_type
  use bshr_const_mod     , only : SHR_CONST_TKFRZ
  implicit none
  ! !ARGUMENTS:
  class(DecompCent_type)     , intent(inout) :: this
  real(r8)                   , intent(in) :: o2b
  type(centuryeca_forc_type) , intent(in) :: centuryeca_forc

  ! !LOCAL VARIABLES:
  real(r8), parameter :: normalization_tref = 15._r8 ! reference temperature for normalizaion (degrees C)
  real(r8)            :: minpsi
  real(r8)            :: maxpsi
  real(r8)            :: normalization_factor
  real(r8)            :: catanf_30
  real(r8)            :: t1
  real(r8)            :: o2w
  real(r8)            :: psi

  associate(                                            &
    temp          => centuryeca_forc%temp    ,          &
    depz          => centuryeca_forc%depz    ,          &
    o2_w2b        => centuryeca_forc%o2_w2b  ,          &
    sucsat        => centuryeca_forc%sucsat  ,          & ! Input:  [real(r8) (:,:)] minimum soil suction [mm]
    soilpsi       => centuryeca_forc%soilpsi ,          & ! Input:  [real(r8) (:,:)] soilwater pontential in each soil layer [MPa]
    Q10           => this%Q10                ,          &
    froz_q10      => this%froz_q10           ,          &
    decomp_depth_efolding => this%decomp_depth_efolding &
  )

  catanf_30 = catanf(30._r8)

  !temperature scalar
  this%t_scalar     = 1._r8
  !use Charlie's Q10 based temperature scalar
  if (temp >= SHR_CONST_TKFRZ) then
    this%t_scalar= (Q10**((temp-(SHR_CONST_TKFRZ+25._r8))/10._r8))
  else
    this%t_scalar= (Q10**(-25._r8/10._r8))*(froz_q10**((temp-SHR_CONST_TKFRZ)/10._r8))
  endif

  ! scale all decomposition rates by a constant to compensate for offset between original CENTURY temp func and Q10
  normalization_factor = (catanf(normalization_tref)/catanf_30) / (q10**((normalization_tref-25._r8)/10._r8))
  this%t_scalar = this%t_scalar * normalization_factor

  !h2osoi_liqure scalar, also follows what Charlie has done
  minpsi = -10.0_r8
  this%w_scalar     = 1._r8
  maxpsi = sucsat * (-9.8e-6_r8)   !kg -> MPa
  psi = min(soilpsi,maxpsi)

  ! decomp only if soilpsi is higher than minpsi, some modification is needed for the following
  ! double check the paper by Wilson and Griffin, 1975
  if (psi > minpsi) then
    this%w_scalar = (log(minpsi/psi)/log(minpsi/maxpsi))
  else
    this%w_scalar = 0._r8
  end if

  !oxygen scalar, this is different from what CLM4.5bgc does, I use a M-M formulation to indicate O2 stress
  !and the O2 budget is done on the fly
  o2w = o2b / o2_w2b
  this%o_scalar = o2w/(o2w+0.02_r8)   !the value 0.22 mol O3/m3 is from Arah and Kirk, 2000

  !depth scalar, according to Koven et al. (2013), BG, the depth scalar is needed to resolve the radiocarbon profile
  this%depth_scalar = exp(-depz/decomp_depth_efolding)

  end associate
  end subroutine set_decompk_scalar

  !-----------------------------------------------------------------------
  function catanf(t1)result(ans)
  !DESCRIPTION:
  ! CENTURY T response function
  use bshr_const_mod, only : SHR_CONST_PI
  implicit none
  real(r8), intent(in) :: t1

  real(r8) :: ans
  ans = 11.75_r8 +(29.7_r8 / SHR_CONST_PI) * atan( SHR_CONST_PI * 0.031_r8  * ( t1 - 15.4_r8 ))

  end function catanf
end module BgcCentCnpDecompType
