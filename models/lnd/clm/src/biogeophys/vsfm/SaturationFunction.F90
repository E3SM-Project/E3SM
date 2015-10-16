#ifdef USE_PETSC_LIB


module SaturationFunction

  ! !USES:
  use clm_varctl         , only : iulog
  use abortutils         , only : endrun
  use shr_log_mod        , only : errMsg => shr_log_errMsg
  !
  ! !PUBLIC TYPES:
  implicit none
  private

#include "finclude/petscsys.h"


  ! Identify saturation function.
  ! These must have unique values (actual values are arbitrary).
  PetscInt, parameter, public :: SAT_FUNC_VAN_GENUCHTEN         = 1301
  PetscInt, parameter, public :: SAT_FUNC_BROOKS_COREY          = 1302
  PetscInt, parameter, public :: SAT_FUNC_SMOOTHED_BROOKS_COREY = 1303
  PetscInt, parameter, public :: RELPERM_FUNC_MUALEM            = 1304
  PetscInt, parameter, public :: RELPERM_FUNC_WEIBULL           = 1305


  type, public :: saturation_params_type
     PetscInt  :: sat_func_type                   ! Identify saturation function using, e.g., `SAT_FUNC_VAN_GENUCHTEN`.
     PetscInt  :: relperm_func_type               ! Idenitfy relative permeability function
     PetscReal :: sat_res                         ! Residual saturation [1].
     PetscReal :: alpha                           ! Empirical constant, \f$ \alpha = 1/p_c^0 \f$ [Pa^{-1}].
     PetscReal :: vg_m, vg_n                      ! Specific to Van Genuchten function [1].
     PetscReal :: bc_lambda                       ! Specific to Brooks-Corey functions [1].
     PetscReal :: sbc_pu, sbc_ps, sbc_b2, sbc_b3  ! Specific to smoothed Brooks-Corey function.
     PetscReal :: w_c, w_d                        ! Specific to relative permeability using Weibull parameterization
   contains
     procedure, public :: Copy    => SatFunc_Copy
     procedure, public :: Init    => SatFunc_Init
  end type saturation_params_type


  public ::                          &
       ! Parameterize a saturation function.
       SatFunc_Set_VG              , &
       SatFunc_Set_BC              , &
       SatFunc_Set_SBC             , &
       SatFunc_Set_SBC_bz2         , &
       SatFunc_Set_SBC_bz3         , &
       SatFunc_Set_Weibull_RelPerm , &
       
       ! Gateway routines (i.e. , generic functions).
       SatFunc_PressToSat          , &
       SatFunc_SatToPress          , &
       SatFunc_PressToRelPerm      , &
       
       ! Find saturation for a specific function.
       SatFunc_PcToSat_VG          , &
       SatFunc_PcToSat_BC          , &
       SatFunc_PcToSat_SBC         , &
       
       ! Find relative permeability for a specific function.
       SatFunc_PcToRelPerm_VG      , &
       SatFunc_PcToRelPerm_BC      , &
       SatFunc_PcToRelPerm_SBC     , &
       
       ! Find capillary pressure for a specific function.
       SatFunc_SatToPc_VG          , &
       SatFunc_SatToPc_BC          , &
       SatFunc_SatToPc_SBC


contains


  !------------------------------------------------------------------------
  subroutine SatFunc_Init(this)
    !
    ! !DESCRIPTION:
    ! Initializes the parameters for saturation function
    !
    implicit none
    !
    ! !ARGUMENTS
    class(saturation_params_type) :: this

    this%sat_func_type     = 0
    this%relperm_func_type = 0

    this%sat_res           = 0.d0
    this%alpha             = 0.d0
    this%vg_m              = 0.d0
    this%vg_n              = 0.d0
    this%bc_lambda         = 0.d0
    this%sbc_pu            = 0.d0
    this%sbc_ps            = 0.d0
    this%sbc_b2            = 0.d0
    this%sbc_b3            = 0.d0
    this%w_c               = 0.d0
    this%w_d               = 0.d0

  end subroutine SatFunc_Init


  !------------------------------------------------------------------------
  subroutine SatFunc_Set_VG(satParams, sat_res, alpha, vg_m)
    !
    ! !DESCRIPTION:
    !
    !
    implicit none
    !
    ! !ARGUMENTS
    type(saturation_params_type) , intent(out) :: satParams
    PetscReal                    , intent(in)  :: sat_res
    PetscReal                    , intent(in)  :: alpha
    PetscReal                    , intent(in)  :: vg_m

    ! Check arguments.
    if( sat_res<0.d0 .or. sat_res>0.5d0  &
         .or. alpha<=0.d0 .or. alpha>2.d0  &
         .or. vg_m<=0.d0 .or. vg_m>=1.d0 ) then
       write(iulog,*) 'SatFunc_Set_VG: bad param'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    call satParams%Init()

    satParams%sat_func_type     = SAT_FUNC_VAN_GENUCHTEN
    satParams%sat_res           = sat_res
    satParams%alpha             = alpha
    satParams%vg_m              = vg_m

    ! Apply Mualem model.
    satParams%relperm_func_type = RELPERM_FUNC_MUALEM
    satParams%vg_n              = 1.d0 / (1.d0 - vg_m)

  end subroutine SatFunc_Set_VG


  !------------------------------------------------------------------------
  subroutine SatFunc_Set_BC(satParams, sat_res, alpha, lambda)
    !
    ! !DESCRIPTION:
    !
    !
    implicit none
    !
    ! !ARGUMENTS
    type(saturation_params_type) , intent(out) :: satParams
    PetscReal                    , intent(in)  :: sat_res
    PetscReal                    , intent(in)  :: alpha
    PetscReal                    , intent(in)  :: lambda

    ! Check arguments.
    if( sat_res<0.d0 .or. sat_res>0.5d0  &
         .or. alpha<=0.d0 .or. alpha>2.d0  &
         .or. lambda<=0.d0 .or. lambda>=2.d0 ) then
       write(iulog,*) 'SatFunc_Set_BC: bad param'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    call satParams%Init()

    satParams%sat_func_type     = SAT_FUNC_BROOKS_COREY
    satParams%relperm_func_type = RELPERM_FUNC_MUALEM
    satParams%sat_res           = sat_res
    satParams%alpha             = alpha
    satParams%bc_lambda         = lambda

  end subroutine SatFunc_Set_BC


  !------------------------------------------------------------------------
  subroutine SatFunc_Set_SBC(satParams, sat_res, alpha, lambda, ps, pu, &
       smoothingRegimeSuperSat)
    !
    ! !DESCRIPTION:
    !
    !
    implicit none
    !
    ! !ARGUMENTS
    type(saturation_params_type) , intent(out) :: satParams
    PetscReal                    , intent(in)  :: sat_res
    PetscReal                    , intent(in)  :: alpha
    PetscReal                    , intent(in)  :: lambda
    PetscReal                    , intent(in)  :: ps
    PetscReal                    , intent(in)  :: pu
    logical                      , intent(out) :: smoothingRegimeSuperSat
    !
    ! !LOCAL VARIABLES:
    PetscReal :: bcAtPu, lambdaDeltaPuOnPu, oneOnDeltaPu
    PetscReal :: deltaPcCrit

    ! Check arguments.
    if( sat_res<0.d0 .or. sat_res>0.5d0  &
         .or. alpha<=0.d0 .or. alpha>2.d0  &
         .or. lambda<=0.d0 .or. lambda>=2.d0  &
         .or. ps<=-1.d0/alpha .or. ps>0.d0  &
         .or. pu>=-1.d0/alpha ) then
       write(iulog,*) 'SatFunc_Set_SBC: bad param'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    call satParams%Init()

    satParams%sat_func_type     = SAT_FUNC_SMOOTHED_BROOKS_COREY
    satParams%relperm_func_type = RELPERM_FUNC_MUALEM
    satParams%sat_res           = sat_res
    satParams%alpha             = alpha
    satParams%bc_lambda         = lambda
    satParams%sbc_pu            = pu
    satParams%sbc_ps            = ps

    ! Find helper constants.
    bcAtPu            = (-alpha*pu)**(-lambda)
    lambdaDeltaPuOnPu = lambda * (1.d0 - ps/pu)
    oneOnDeltaPu      = 1.d0 / (pu - ps)

    ! Store coefficients for cubic function.
    satParams%sbc_b2 = -(3.d0 - bcAtPu*(3.d0+lambdaDeltaPuOnPu)) * oneOnDeltaPu * oneOnDeltaPu
    satParams%sbc_b3 =  (2.d0 - bcAtPu*(2.d0+lambdaDeltaPuOnPu)) * oneOnDeltaPu * oneOnDeltaPu * oneOnDeltaPu

    ! Check for pathological behavior in the cubic smoothing regime.
    smoothingRegimeSuperSat = .false.
    if( satParams%sbc_b2 > 0.d0 ) then
       deltaPcCrit = -2 * satParams%sbc_b2 / (3*satParams%sbc_b3)
       if( deltaPcCrit<0.d0 .and. deltaPcCrit>(pu-ps) ) then
          ! Here, saturation curve has a local maximum inside the
          ! cubic smoothing regime.
          smoothingRegimeSuperSat = .true.
       endif
    endif

  end subroutine SatFunc_Set_SBC

  !------------------------------------------------------------------------
  subroutine SatFunc_Set_SBC_bz2(satParams, sat_res, alpha, lambda, ps)
    !
    ! !DESCRIPTION:
    !
    !
    implicit none
    !
    ! !ARGUMENTS
    type(saturation_params_type) , intent(out) :: satParams
    PetscReal                    , intent(in)  :: sat_res
    PetscReal                    , intent(in)  :: alpha
    PetscReal                    , intent(in)  :: lambda
    PetscReal                    , intent(in)  :: ps
    !
    ! !LOCAL VARIABLES:
    PetscReal                                  :: pu
    PetscReal                                  :: bcAtPu
    PetscReal                                  :: lambdaDeltaPuOnPu
    PetscReal                                  :: oneOnDeltaPu

    ! Check arguments.
    if( sat_res<0.d0 .or. sat_res>0.5d0  &
         .or. alpha<=0.d0 .or. alpha>2.d0  &
         .or. lambda<=0.d0 .or. lambda>=2.d0  &
         .or. ps<=-1.d0/alpha .or. ps>0.d0 ) then
       write(iulog,*) 'SatFunc_Set_SBC_bz2: bad param'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    call satParams%Init()

    satParams%sat_func_type     = SAT_FUNC_SMOOTHED_BROOKS_COREY
    satParams%relperm_func_type = RELPERM_FUNC_MUALEM
    satParams%sat_res           = sat_res
    satParams%alpha             = alpha
    satParams%bc_lambda         = lambda
    satParams%sbc_ps            = ps

    ! Choose `pu` that forces `sbc_b2 = 0`.
    pu               = findGu_SBC_zeroCoeff(lambda, 3, -alpha*ps) / (-alpha)
    satParams%sbc_pu = pu

    ! Find helper constants.
    bcAtPu            = (-alpha*pu)**(-lambda)
    lambdaDeltaPuOnPu = lambda * (1.d0 - ps/pu)
    oneOnDeltaPu      = 1.d0 / (pu - ps)

    ! Store coefficients for cubic function.
    satParams%sbc_b2 = 0.d0
    satParams%sbc_b3 = (2.d0 - bcAtPu*(2.d0+lambdaDeltaPuOnPu)) * oneOnDeltaPu * oneOnDeltaPu * oneOnDeltaPu
    if( satParams%sbc_b3 <= 0.d0 ) then
       write(iulog,*) 'SatFunc_Set_SBC_bz2: b3 <= 0'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

  end subroutine SatFunc_Set_SBC_bz2


  !------------------------------------------------------------------------
  subroutine SatFunc_Set_SBC_bz3(satParams, sat_res, alpha, lambda, ps)
    !
    ! !DESCRIPTION:
    !
    !
    implicit none
    !
    ! !ARGUMENTS
    type(saturation_params_type) , intent(out) :: satParams
    PetscReal                    , intent(in)  :: sat_res
    PetscReal                    , intent(in)  :: alpha
    PetscReal                    , intent(in)  :: lambda
    PetscReal                    , intent(in)  :: ps
    !
    ! !LOCAL VARIABLES:
    PetscReal                                  :: pu
    PetscReal                                  :: bcAtPu
    PetscReal                                  :: lambdaDeltaPuOnPu
    PetscReal                                  :: oneOnDeltaPu

    ! Check arguments.
    if( sat_res<0.d0 .or. sat_res>0.5d0  &
         .or. alpha<=0.d0 .or. alpha>2.d0  &
         .or. lambda<=0.d0 .or. lambda>=2.d0  &
         .or. ps<=-1.d0/alpha .or. ps>0.d0 ) then
       write(iulog,*) 'SatFunc_Set_SBC_bz3: bad param'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    satParams%sat_func_type     = SAT_FUNC_SMOOTHED_BROOKS_COREY
    satParams%relperm_func_type = RELPERM_FUNC_MUALEM
    satParams%sat_res           = sat_res
    satParams%alpha             = alpha
    satParams%bc_lambda         = lambda
    satParams%sbc_ps            = ps

    ! Choose `pu` that forces `sbc_b3 = 0`.
    pu               = findGu_SBC_zeroCoeff(lambda, 2, -alpha*ps) / (-alpha)
    satParams%sbc_pu = pu

    ! Find helper constants.
    bcAtPu            = (-alpha*pu)**(-lambda)
    lambdaDeltaPuOnPu = lambda * (1.d0 - ps/pu)
    oneOnDeltaPu      = 1.d0 / (pu - ps)

    ! Store coefficients for cubic function.
    satParams%sbc_b2 = -(3.d0 - bcAtPu*(3.d0+lambdaDeltaPuOnPu)) * oneOnDeltaPu * oneOnDeltaPu
    if( satParams%sbc_b2 >= 0.d0 ) then
       write(iulog,*) 'SatFunc_Set_SBC_bz3: b2 >= 0'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif
    satParams%sbc_b3 = 0.d0

  end subroutine SatFunc_Set_SBC_bz3


  !------------------------------------------------------------------------
  ! Find `pu` that forces a coefficient of the smoothing cubic polynomial to zero.
  !
  ! Work in terms of multipliers of `pc0`:
  !
  ! + Argument `gs` satisfies `ps = gs*pc0`.
  ! + Return `gu` such that `pu = gu*pc0`.
  !
  ! Argument `AA`:
  !
  ! + To set `b2 = 0`, let `A = 3`.
  ! + To set `b3 = 0`, let `A = 2`.
  !
  PetscReal function findGu_SBC_zeroCoeff(lambda, AA, gs)
    !
    ! !DESCRIPTION:
    !
    !
    implicit none
    !
    ! !ARGUMENTS
    PetscReal , intent(in) :: lambda
    PetscReal , intent(in) :: gs
    integer   , intent(in) :: AA
    !
    ! !LOCAL VARIABLES:
    PetscReal              :: guLeft, gu, guRight  ! Bracketed search.
    PetscReal              :: deltaGu, resid, dr_dGu  ! Newton-Raphson search.
    PetscReal              :: guInv, guToMinusLam, gsOnGu  ! Helper variables.
    PetscReal, parameter   :: relTol = 1.d-12

    ! Check arguments.
    !   Note this is more for documentation than anything else-- this
    ! fcn should only get used internally, by trusted callers.
    if( lambda<=0.d0 .or. lambda>=2.d0  &
         .or. (AA/=2 .and. AA/=3)  &
         .or. gs>=1.d0 .or. gs<0.d0 ) then
       write(iulog,*) 'findGu_SBC_zeroCoeff: bad param'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    ! Approach:
    ! + Bracketed Newton-Raphson search.
    ! + Note expect `1 < gu <= gu{gs=0}`.
    ! + Note if this was a critical inner loop, could try solving for
    !     `gui == 1/gu`, rather than for `gu`, in order to avoid division.
    !     Could also try using Ridder's method, since the residual here
    !     has a strong exponential component.

    ! Initialize.
    gu = (AA / (AA + lambda))**(-1.d0/lambda)  ! Solution if `gs = 0`.

    ! Search for root, using bracketed Newton-Raphson.
    !   Not necessary if `gs = 0`.
    if( gs > 0.d0 ) then
       guLeft = 1.d0
       guRight = gu
       do
          ! Here, assume:
          ! + Have an bracket on the root, between `guLeft` and `guRight`.
          ! + The derivative `dr/d{gu} > 0`.
          ! + The residual `r{guLeft} < 0`, and `r{guRight} > 0`.
          ! + Have a current guess `gu` at the root.  However, that guess
          !     might not lie in the bracket (and does not at first iteration).

          ! Reset `gu` using bisection if necessary.
          if( gu<=guLeft .or. gu>=guRight ) then
             gu = guLeft + 0.5d0*(guRight - guLeft)
          endif

          ! Find residual.
          guInv = 1.d0 / gu
          guToMinusLam = gu**(-lambda)  ! Could also do `guInv**lambda`; not sure if any numerical consequences.
          gsOnGu = gs * guInv
          resid = AA - guToMinusLam*(AA + lambda - lambda*gsOnGu)

          ! Update bracket.
          if( resid < 0.d0 ) then
             guLeft = gu
          else
             guRight = gu
          endif

          ! Find next guess using Newton-Raphson's method.
          dr_dGu = (1.d0 + lambda) * (1.d0 - gsOnGu) + (AA - 1)
          dr_dGu = lambda * guToMinusLam * guInv * dr_dGu
          deltaGu = resid / dr_dGu
          ! write(unit=*, fmt='("findGu_SBC_zeroCoeff, NR step: ", 6(a,g15.6))')  &
          !     'guLeft', guLeft, 'gu', gu, 'guRight', guRight,  &
          !     'deltaGu', deltaGu, 'resid', resid, 'dr_dGu', dr_dGu
          gu = gu - deltaGu

          ! Test for convergence.
          !   Note this test implicitly also tests `resid == 0`.
          if( abs(deltaGu) < relTol*abs(gu) ) then
             exit
          endif
       enddo
    endif

    ! Finish up.
    !   Note assuming the last Newton-Raphson step landed in the bracket,
    ! and had a smaller residual than either of the bracket points.  This
    ! seems a safe enough assumption, compared to cost of tracking residuals.
    findGu_SBC_zeroCoeff = gu

  end function findGu_SBC_zeroCoeff


  !------------------------------------------------------------------------
  subroutine SatFunc_Set_Weibull_RelPerm(satParams, d, c)
    !
    ! !DESCRIPTION:
    ! Sets parameters for Weibull model for relative permeability
    !
    implicit none
    !
    ! !ARGUMENTS
    type(saturation_params_type) , intent(inout) :: satParams
    PetscReal                    , intent(in)    :: d
    PetscReal                    , intent(in)    :: c
    !
    ! !LOCAL VARIABLES:

    satParams%relperm_func_type  = RELPERM_FUNC_WEIBULL
    satParams%w_d                = d
    satParams%w_c                = c

  end subroutine SatFunc_Set_Weibull_RelPerm


  !------------------------------------------------------------------------
  subroutine SatFunc_PressToSat(satParams, press, sat, dsat_dP)
    !
    ! !DESCRIPTION:
    !
    !
    use MultiPhysicsProbConstants, only : PRESSURE_REF
    !
    implicit none
    !
    ! !ARGUMENTS
    type(saturation_params_type) , intent(in)  :: satParams
    PetscReal                    , intent(in)  :: press
    PetscReal                    , intent(out) :: sat
    PetscReal                    , intent(out) :: dsat_dP
    !
    ! !LOCAL VARIABLES:
    PetscReal                                  :: pc  ! Capillary pressure

    pc = press - PRESSURE_REF

    select case(satParams%sat_func_type)
    case (SAT_FUNC_VAN_GENUCHTEN)
       call SatFunc_PcToSat_VG(satParams, pc, sat, dsat_dP)
    case (SAT_FUNC_BROOKS_COREY)
       call SatFunc_PcToSat_BC(satParams, pc, sat, dsat_dP)
    case (SAT_FUNC_SMOOTHED_BROOKS_COREY)
       call SatFunc_PcToSat_SBC(satParams, pc, sat, dsat_dP)
    case default
       write(iulog,*) 'SatFunc_PressToSat: Unknown type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine SatFunc_PressToSat


  !------------------------------------------------------------------------
  subroutine SatFunc_PressToRelPerm(satParams, press, frac_liq, kr, dkr_dP)
    !
    ! !DESCRIPTION:
    !
    !
    use MultiPhysicsProbConstants, only : PRESSURE_REF
    !
    implicit none
    !
    ! !ARGUMENTS
    type(saturation_params_type) , intent(in)  :: satParams
    PetscReal                    , intent(in)  :: press
    PetscReal                    , intent(in)  :: frac_liq
    PetscReal                    , intent(out) :: kr
    PetscReal                    , intent(out) :: dkr_dP
    !
    ! !LOCAL VARIABLES:
    PetscReal                                  :: pc  ! Capillary pressure

    pc = press - PRESSURE_REF

    select case(satParams%relperm_func_type)
    case (RELPERM_FUNC_MUALEM)
       select case(satParams%sat_func_type)
       case (SAT_FUNC_VAN_GENUCHTEN)
          call SatFunc_PcToRelPerm_VG(satParams, pc, kr, dkr_dP)
       case (SAT_FUNC_BROOKS_COREY)
          call SatFunc_PcToRelPerm_BC(satParams, pc, frac_liq, kr, dkr_dP)
       case (SAT_FUNC_SMOOTHED_BROOKS_COREY)
          call SatFunc_PcToRelPerm_SBC(satParams, pc, kr, dkr_dP)
       case default
          write(iulog,*) 'SatFunc_PressToRelPerm: Unknown type'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end select

    case(RELPERM_FUNC_WEIBULL)
       call SatFunc_PcToRelPerm_Weibull(satParams, pc, kr, dkr_dP)
    case default
       write(iulog,*)'SatFunc_PressToRelPerm: Unknown type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine SatFunc_PressToRelPerm


  !------------------------------------------------------------------------
  subroutine SatFunc_PcToRelPerm_Weibull(satParams, pc, kr, dkr_dP)
    !
    ! !DESCRIPTION:
    ! Computes relative permeability and its derivative w.r.t. pressure
    ! for a Weibull relative permeability model
    !
    implicit none
    !
    ! !ARGUMENTS
    type(saturation_params_type) , intent(in)  :: satParams
    PetscReal                    , intent(in)  :: pc
    PetscReal                    , intent(out) :: kr
    PetscReal                    , intent(out) :: dkr_dP
    !
    ! !LOCAL VARIABLES:
    PetscReal                                  :: AA

    if (pc >= 0.d0) then
       kr      = 1.d0
       dkr_dP  = 0.d0
    else
       AA     = (-pc/satParams%w_d)**satParams%w_c
       kr     = exp(-AA)
       dkr_dP = -satParams%w_c/pc*AA*kr
    endif

  end subroutine SatFunc_PcToRelPerm_Weibull


  !------------------------------------------------------------------------
  subroutine SatFunc_SatToPress(satParams, sat, press)
    !
    ! !DESCRIPTION:
    !
    !
    use MultiPhysicsProbConstants, only : PRESSURE_REF
    !
    implicit none
    !
    ! !ARGUMENTS
    type(saturation_params_type) , intent(in)  :: satParams
    PetscReal                    , intent(in)  :: sat
    PetscReal                    , intent(out) :: press
    !
    ! !LOCAL VARIABLES:
    PetscReal                                  :: pc  ! Capillary pressure.

    select case(satParams%sat_func_type)
    case (SAT_FUNC_VAN_GENUCHTEN)
       call SatFunc_SatToPc_VG(satParams, sat, pc)
    case (SAT_FUNC_BROOKS_COREY)
       call SatFunc_SatToPc_BC(satParams, sat, pc)
    case (SAT_FUNC_SMOOTHED_BROOKS_COREY)
       call SatFunc_SatToPc_SBC(satParams, sat, pc)
    case default
       write(iulog,*) 'SatFunc_SatToPress: Unrecognized type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    press = pc + PRESSURE_REF

  end subroutine SatFunc_SatToPress


  !------------------------------------------------------------------------
  subroutine SatFunc_PcToSat_VG(satParams, pc, sat, dsat_dP)
    !
    ! !DESCRIPTION:
    !
    !
    implicit none
    !
    ! !ARGUMENTS
    type(saturation_params_type) , intent(in)  :: satParams
    PetscReal                    , intent(in)  :: pc
    PetscReal                    , intent(out) :: sat
    PetscReal                    , intent(out) :: dsat_dP
    !
    ! !LOCAL VARIABLES:
    PetscReal                                  :: sat_res
    PetscReal                                  :: alpha
    PetscReal                                  :: mm
    PetscReal                                  :: nn
    PetscReal                                  :: pc_alpha_n
    PetscReal                                  :: one_plus_pc_alpha_n
    PetscReal                                  :: Se
    PetscReal                                  :: dSe_dpc
    PetscReal                                  :: AA

    sat_res = satParams%sat_res
    alpha   = satParams%alpha
    mm      = satParams%vg_m
    nn      = satParams%vg_n

    if( pc < 0.d0 ) then
       pc_alpha_n          = (-alpha*pc)**nn
       one_plus_pc_alpha_n = 1.d0 + pc_alpha_n

       Se                  = one_plus_pc_alpha_n**(-mm)
       sat                 = sat_res + (1.d0 - sat_res)*Se

       ! Find common subexpression
       ! AA = 1 - s_e^{1/m} = (-alpha*p_c)^{n} / (1 + (-alpha*pc)^n)
       AA = pc_alpha_n / one_plus_pc_alpha_n

       dSe_dpc = -mm*nn*Se* AA / pc
       dsat_dp = (1.d0 - sat_res)*dSe_dpc
    else
       ! Here, `pc >= 0`.
       sat     = 1.d0
       dsat_dP = 0.d0
    endif

  end subroutine SatFunc_PcToSat_VG


  !------------------------------------------------------------------------
  subroutine SatFunc_PcToRelPerm_VG(satParams, pc, kr, dkr_dP)
    !
    ! !DESCRIPTION:
    !
    !
    implicit none
    !
    ! !ARGUMENTS
    type(saturation_params_type) , intent(in)  :: satParams
    PetscReal                    , intent(in)  :: pc
    PetscReal                    , intent(out) :: kr
    PetscReal                    , intent(out) :: dkr_dP
    !
    ! !LOCAL VARIABLES:
    PetscReal                                  :: sat
    PetscReal                                  :: sat_res
    PetscReal                                  :: alpha
    PetscReal                                  :: mm
    PetscReal                                  :: nn
    PetscReal                                  :: pc_alpha_n
    PetscReal                                  :: one_plus_pc_alpha_n
    PetscReal                                  :: Se
    PetscReal                                  :: dSe_dpc
    PetscReal                                  :: AA
    PetscReal                                  :: BB
    PetscReal                                  :: dkr_dSe

    sat_res = satParams%sat_res
    alpha   = satParams%alpha
    mm      = satParams%vg_m
    nn      = satParams%vg_n

    if( pc < 0.d0 ) then
       pc_alpha_n          = (-alpha*pc)**nn
       one_plus_pc_alpha_n = 1.d0 + pc_alpha_n

       Se                  = one_plus_pc_alpha_n**(-mm)

       ! Find common subexpression
       ! AA = 1 - s_e^{1/m} = (-alpha*p_c)^{n} / (1 + (-alpha*pc)^n)
       AA = pc_alpha_n / one_plus_pc_alpha_n

       dSe_dpc = -mm*nn*Se* AA / pc

       ! Find common subexpression
       ! BB = 1 - (1 - s_e^{1/m})^{m}
       BB = 1.d0 - AA**mm

       kr      = sqrt(Se) * BB*BB
       dkr_dSe = 0.5d0*kr/Se +  &
            2.d0 * Se**(1d0/mm - 0.5d0) * AA**(mm - 1.d0) * BB
       dkr_dp  = dkr_dSe * dSe_dpc
    else
       ! Here, `pc > = 0`.
       kr     = 1.d0
       dkr_dP = 0.d0
    endif

  end subroutine SatFunc_PcToRelPerm_VG


  !------------------------------------------------------------------------
  subroutine SatFunc_SatToPc_VG(satParams, sat, pc)
    !
    ! !DESCRIPTION:
    !
    !
    implicit none
    !
    ! !ARGUMENTS
    type(saturation_params_type) , intent(in)  :: satParams
    PetscReal                    , intent(in)  :: sat
    PetscReal                    , intent(out) :: pc
    !
    ! !LOCAL VARIABLES:
    PetscReal                                  :: sat_res
    PetscReal                                  :: alpha
    PetscReal                                  :: mm
    PetscReal                                  :: nn
    PetscReal                                  :: Se

    sat_res = satParams%sat_res
    alpha   = satParams%alpha
    mm      = satParams%vg_m
    nn      = satParams%vg_n

    if( sat < 1.d0 ) then
       Se = (sat - sat_res)/(1.d0 - sat_res)
       if( Se < 0.d0 ) then
          ! Here, `sat < sat_res`.
          Se = 0.d0
       endif
       pc = -(Se**(-1.d0/mm) - 1.d0)**(1.d0/nn) / alpha
    else
       pc = 0.d0
    endif

  end subroutine SatFunc_SatToPc_VG


  !------------------------------------------------------------------------
  subroutine SatFunc_PcToSat_BC(satParams, pc, sat, dsat_dP)
    !
    ! !DESCRIPTION:
    !
    !
    implicit none
    !
    ! !ARGUMENTS
    type(saturation_params_type) , intent(in)  :: satParams
    PetscReal                    , intent(in)  :: pc
    PetscReal                    , intent(out) :: sat
    PetscReal                    , intent(out) :: dsat_dP
    !
    ! !LOCAL VARIABLES:
    PetscReal                                  :: sat_res
    PetscReal                                  :: alpha
    PetscReal                                  :: lambda
    PetscReal                                  :: pc_alpha
    PetscReal                                  :: Se
    PetscReal                                  :: dSe_dpc

    sat_res = satParams%sat_res
    alpha   = satParams%alpha
    lambda  = satParams%bc_lambda

    pc_alpha = -alpha*pc
    if( pc_alpha > 1.d0 ) then
       Se      = pc_alpha**(-lambda)
       sat     = sat_res + (1.d0 - sat_res)*Se

       dSe_dpc = -lambda*Se/pc
       dsat_dp = (1.d0 - sat_res)*dSe_dpc
    else
       ! Here, `pc >= pc0`.
       sat     = 1.d0
       dsat_dP = 0.d0
    endif

  end subroutine SatFunc_PcToSat_BC


  !------------------------------------------------------------------------
  subroutine SatFunc_PcToRelPerm_BC(satParams, pc, frac_liq, kr, dkr_dP)
    !
    ! !DESCRIPTION:
    !
    !
    implicit none
    !
    ! !ARGUMENTS
    type(saturation_params_type) , intent(in)  :: satParams
    PetscReal                    , intent(in)  :: pc
    PetscReal                    , intent(in)  :: frac_liq
    PetscReal                    , intent(out) :: kr
    PetscReal                    , intent(out) :: dkr_dP
    !
    ! !LOCAL VARIABLES:
    PetscReal                                  :: sat
    PetscReal                                  :: sat_res
    PetscReal                                  :: alpha
    PetscReal                                  :: lambda
    PetscReal                                  :: pc_alpha
    PetscReal                                  :: Se
    PetscReal                                  :: dSe_dpc
    PetscReal                                  :: dkr_dSe

    sat_res = satParams%sat_res
    alpha   = satParams%alpha
    lambda  = satParams%bc_lambda

    pc_alpha = -alpha*pc
    if( pc_alpha > 1.d0 ) then
       Se      = pc_alpha**(-lambda)
       sat     = sat_res + (1.d0 - sat_res)*Se

       dSe_dpc = -lambda*Se/pc

       kr      = Se ** (2.5d0 + 2.d0/lambda)

       dkr_dSe = (2.5d0 + 2.d0/lambda)*kr/Se
       dkr_dP  = dkr_dSe*dSe_dpc
    else
       ! Here, `pc >= pc0`.
       kr     = 1.d0
       dkr_dP = 0.d0
    endif

    kr     = frac_liq*kr
    dkr_dP = frac_liq*dkr_dP

  end subroutine SatFunc_PcToRelPerm_BC


  !------------------------------------------------------------------------
  subroutine SatFunc_SatToPc_BC(satParams, sat, pc)
    !
    ! !DESCRIPTION:
    !
    !
    implicit none
    !
    ! !ARGUMENTS
    type(saturation_params_type) , intent(in)  :: satParams
    PetscReal                    , intent(in)  :: sat
    PetscReal                    , intent(out) :: pc
    !
    ! !LOCAL VARIABLES:
    PetscReal                                  :: sat_res
    PetscReal                                  :: alpha
    PetscReal                                  :: lambda
    PetscReal                                  :: Se

    sat_res = satParams%sat_res
    alpha   = satParams%alpha
    lambda  = satParams%bc_lambda

    if( sat < 1.d0 ) then
       Se = (sat - sat_res)/(1.d0 - sat_res)
       pc = -Se**(-1.d0/lambda)/alpha
    else
       pc = 0.d0
    endif

  end subroutine SatFunc_SatToPc_BC


  !------------------------------------------------------------------------
  subroutine SatFunc_PcToSat_SBC(satParams, pc, sat, dsat_dP)
    !
    ! !DESCRIPTION:
    !
    !
    implicit none
    !
    ! !ARGUMENTS
    type(saturation_params_type) , intent(in)  :: satParams
    PetscReal                    , intent(in)  :: pc
    PetscReal                    , intent(out) :: sat
    PetscReal                    , intent(out) :: dsat_dP
    !
    ! !LOCAL VARIABLES:
    PetscReal                                  :: sat_res
    PetscReal                                  :: alpha
    PetscReal                                  :: lambda
    PetscReal                                  :: Se
    PetscReal                                  :: deltaPc
    PetscReal                                  :: dSe_dpc

    sat_res = satParams%sat_res
    alpha   = satParams%alpha
    lambda  = satParams%bc_lambda

    if( pc <= satParams%sbc_pu ) then
       ! Unsaturated full Brooks-Corey regime.
       ! Here, `pc <= pu < 0`.
       Se      = (-alpha*pc)**(-lambda)
       sat     = sat_res + (1.d0 - sat_res)*Se

       dSe_dpc = -lambda*Se/pc
       dsat_dp = (1.d0 - sat_res)*dSe_dpc
    elseif( pc < satParams%sbc_ps ) then
       ! Cubic smoothing regime.
       ! Here, `pu < pc < ps <= 0`.
       deltaPc = pc - satParams%sbc_ps
       Se      = 1.d0 + deltaPc*deltaPc*(satParams%sbc_b2 + deltaPc*satParams%sbc_b3)
       sat     = sat_res + (1.d0 - sat_res)*Se

       dSe_dpc = deltaPc*(2*satParams%sbc_b2 + 3*deltaPc*satParams%sbc_b3)
       dsat_dp = (1.d0 - sat_res)*dSe_dpc
    else
       ! Saturated regime.
       ! Here, `pc >= ps`.
       sat       = 1.d0
       dsat_dP   = 0.d0
    endif

  end subroutine SatFunc_PcToSat_SBC


  !------------------------------------------------------------------------
  subroutine SatFunc_PcToRelPerm_SBC(satParams, pc, kr, dkr_dP)
    !
    ! !DESCRIPTION:
    !
    !
    implicit none
    !
    ! !ARGUMENTS
    type(saturation_params_type) , intent(in)  :: satParams
    PetscReal                    , intent(in)  :: pc
    PetscReal                    , intent(out) :: kr
    PetscReal                    , intent(out) :: dkr_dP
    !
    ! !LOCAL VARIABLES:
    PetscReal                                  :: sat_res
    PetscReal                                  :: alpha
    PetscReal                                  :: lambda
    PetscReal                                  :: Se
    PetscReal                                  :: deltaPc
    PetscReal                                  :: dSe_dpc
    PetscReal                                  :: dkr_dSe

    sat_res = satParams%sat_res
    alpha   = satParams%alpha
    lambda  = satParams%bc_lambda

    if( pc <= satParams%sbc_pu ) then
       ! Unsaturated full Brooks-Corey regime.
       ! Here, `pc <= pu < 0`.
       Se      = (-alpha*pc)**(-lambda)

       dSe_dpc = -lambda*Se/pc

       kr      = Se ** (2.5d0 + 2.d0/lambda)

       dkr_dSe = (2.5d0 + 2.d0/lambda)*kr/Se
       dkr_dp  = dkr_dSe*dSe_dpc
    elseif( pc < satParams%sbc_ps ) then
       ! Cubic smoothing regime.
       ! Here, `pu < pc < ps <= 0`.
       deltaPc = pc - satParams%sbc_ps
       Se      = 1.d0 + deltaPc*deltaPc*(satParams%sbc_b2 + deltaPc*satParams%sbc_b3)

       dSe_dpc = deltaPc*(2*satParams%sbc_b2 + 3*deltaPc*satParams%sbc_b3)

       ! For `kr`, use same expressions as for full Brooks-Corey regime.
       ! Know this gives different `kr` than if integrate Mualem's formula;
       ! however, provided `satParams%sbc_ps` is sufficiently close to
       ! `pc0  = -1/alpha`, errors are fairly small.
       kr      = Se ** (2.5d0 + 2.d0/lambda)

       dkr_dSe = (2.5d0 + 2.d0/lambda)*kr/Se
       dkr_dp  = dkr_dSe*dSe_dpc
    else
       ! Saturated regime.
       ! Here, `pc >= ps`.
       kr        = 1.d0
       dkr_dP    = 0.d0
    endif

  end subroutine SatFunc_PcToRelPerm_SBC


  !------------------------------------------------------------------------
  subroutine SatFunc_SatToPc_SBC(satParams, sat, pc)
    !
    ! !DESCRIPTION:
    !
    !
    implicit none
    !
    ! !ARGUMENTS
    type(saturation_params_type) , intent(in)  :: satParams
    PetscReal                    , intent(in)  :: sat
    PetscReal                    , intent(out) :: pc
    !
    ! !LOCAL VARIABLES:
    PetscReal                                  :: sat_res
    PetscReal                                  :: alpha
    PetscReal                                  :: lambda
    PetscReal                                  :: Se
    PetscReal                                  :: xL
    PetscReal                                  :: xc
    PetscReal                                  :: xR
    PetscReal                                  :: resid
    PetscReal                                  :: dx
    PetscReal, parameter                       :: relTol = 1.d-9

    sat_res = satParams%sat_res
    alpha   = satParams%alpha
    lambda  = satParams%bc_lambda

    if( sat < 1.d0 ) then
       ! Find the `pc` that satisfies the unmodified Brooks-Corey function.
       Se = (sat - sat_res)/(1.d0 - sat_res)
       pc = -(Se**(-1.d0/lambda)) / alpha
       if( pc > satParams%sbc_pu ) then
          ! Here, solution is in the cubic smoothing regime.
          if( satParams%sbc_b2 == 0.d0 ) then
             ! Note know `b3 > 0`.
             pc = satParams%sbc_ps - ((1.d0 - Se) / satParams%sbc_b3)**(1.d0/3.d0)
          elseif( satParams%sbc_b3 == 0.d0 ) then
             ! Note know `b2 < 0`.
             pc = satParams%sbc_ps - sqrt((Se - 1.d0) / satParams%sbc_b2)
          else
             ! Here, want to solve general cubic
             ! `1 + b2*x^2 + b3*x^3 = Se`
             ! where `x = pc - pu`.
             ! Write as residual function
             ! `r = x^2 * (b2 + b3*x) + (1 - Se)`.
             ! Perform a Newton-Raphson search on `x`.
             ! Have
             ! `dr/dx = x*(2*b2 + 3*b3*x)`
             ! And Newton-Raphson sets
             ! `x[i+1] = x[i] - r[i]/(dr/dx[i])`.
             ! Note that r{0} = 1 - Se > 0.
             ! Therefore maintain the right bracket as having a positive
             ! residual, and the left bracket as having a negative residual.
             ! Note that it is possible, due to numerical effects with `pc`
             ! very close to `pu`, to get an `xL` with a positive residual.
             ! However, in this case also have `xc` very close to `xL`, and
             ! the Newton-Raphson search will converge after a single step.
             ! Therefore do not insert a special test to catch the case here.
             xL = satParams%sbc_pu - satParams%sbc_ps
             xR = 0.d0
             xc = pc - satParams%sbc_ps
             ! write(unit=*, fmt='("SatFunc_SatToPc_SBC: NR search:", 6(a,g15.6))')  &
             !     ' pu', satParams%sbc_pu, ' ps', satParams%sbc_ps,  &
             !     ' xL', xL, ' xR', xR,  &
             !     ' r{xL}', xL*xL*(satParams%sbc_b2 + satParams%sbc_b3*xL) + 1.d0 - Se,  &
             !     ' r{xR}', xR*xR*(satParams%sbc_b2 + satParams%sbc_b3*xR) + 1.d0 - Se
             do
                ! Here, assume:
                ! + Have a bracket on the root, between `xL` and `xR`.
                ! + The residual `r{xL} < 0` and `r{xR} > 0`.
                ! + Have a current guess `xc` at the root.  However, that guess
                !     might not lie in the bracket.

                ! Reset `xc` using bisection if necessary.
                if( xc<=xL .or. xc>=xR ) then
                   ! write(unit=*, fmt='("Bisecting")')
                   xc = xL + 0.5d0*(xR - xL)
                endif

                ! Find NR step.
                dx = satParams%sbc_b3 * xc
                resid = xc*xc*(satParams%sbc_b2 + dx) + 1.d0 - Se
                dx = resid / (xc*(2.d0*satParams%sbc_b2 + 3.d0*dx))

                ! Update bracket.
                if( resid > 0.d0 ) then
                   xR = xc
                else
                   xL = xc
                endif

                ! Take the Newton-Raphson step.
                xc = xc - dx
                ! write(unit=*, fmt='(6(a,g15.6))')  &
                !     ' xL', xL, ' xc', xc, ' xR', xR,  &
                !     ' r{xL}', xL*xL*(satParams%sbc_b2 + satParams%sbc_b3*xL) + 1.d0 - Se,  &
                !     ' r{xc}', xc*xc*(satParams%sbc_b2 + satParams%sbc_b3*xc) + 1.d0 - Se,  &
                !     ' r{xR}', xR*xR*(satParams%sbc_b2 + satParams%sbc_b3*xR) + 1.d0 - Se

                ! Test for convergence.
                !   Note this test implicitly also tests `resid == 0`.
                if( abs(dx) < -relTol*satParams%sbc_pu ) then
                   exit
                endif
             enddo

             ! Here, have `xc = pc - ps`.
             pc = xc + satParams%sbc_ps
          endif
       endif
    else
       pc = 0.d0
    endif

  end subroutine SatFunc_SatToPc_SBC


  !------------------------------------------------------------------------
  !
  !
  subroutine SatFunc_Copy(this, satParams)
    !
    ! !DESCRIPTION:
    !
    !
    implicit none
    !
    ! !ARGUMENTS
    class(saturation_params_type)             :: this
    class(saturation_params_type), intent(in) :: satParams

    ! Local variables.
    logical :: smoothingRegimeSuperSat

    select case(satParams%sat_func_type)
    case (SAT_FUNC_VAN_GENUCHTEN)

       call SatFunc_Set_VG(this,  &
            satParams%sat_res,    &
            satParams%alpha,      &
            satParams%vg_m)

    case (SAT_FUNC_BROOKS_COREY)

       call SatFunc_Set_BC(this,  &
            satParams%sat_res,    &
            satParams%alpha,      &
            satParams%bc_lambda)

    case (SAT_FUNC_SMOOTHED_BROOKS_COREY)

       call SatFunc_Set_SBC(this, &
            satParams%sat_res,    &
            satParams%alpha,      &
            satParams%bc_lambda,  &
            satParams%sbc_ps,     &
            satParams%sbc_pu,     &
            smoothingRegimeSuperSat)

    case default
       write(iulog,*) 'SaturationPressToSatSaturationFunctionCopy: Unknown type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine SatFunc_Copy


end module SaturationFunction
#endif
