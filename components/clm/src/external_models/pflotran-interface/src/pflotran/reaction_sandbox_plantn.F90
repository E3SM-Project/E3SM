module Reaction_Sandbox_PlantN_class

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Reaction_Sandbox_Base_class
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  use PFLOTRAN_Constants_module
  use Utility_module, only : HFunctionSmooth
  
  implicit none
  
  private
  
  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_plantn_type
    PetscInt  :: ispec_nh4
    PetscInt  :: ispec_no3
    PetscInt  :: ispec_plantn
    PetscInt  :: ispec_plantndemand
    PetscInt  :: ispec_plantnh4uptake
    PetscInt  :: ispec_plantno3uptake
    PetscReal :: rate_plantndemand
    PetscReal :: half_saturation_nh4
    PetscReal :: half_saturation_no3
    PetscReal :: inhibition_nh4_no3
    PetscReal :: x0eps_nh4
    PetscReal :: x0eps_no3

  contains
    procedure, public :: ReadInput => PlantNRead
    procedure, public :: Setup => PlantNSetup
    procedure, public :: Evaluate => PlantNReact
    procedure, public :: Destroy => PlantNDestroy
  end type reaction_sandbox_plantn_type

  public :: PlantNCreate

contains

! ************************************************************************** !
!
! PlantNCreate: Allocates plantn reaction object.
! author: Guoping Tang
! date: 09/09/2013
!
! ************************************************************************** !
function PlantNCreate()

  implicit none
  
  class(reaction_sandbox_plantn_type), pointer :: PlantNCreate

  allocate(PlantNCreate)
  PlantNCreate%ispec_nh4 = 0
  PlantNCreate%ispec_no3 = 0
  PlantNCreate%ispec_plantn = 0
  PlantNCreate%ispec_plantndemand = 0
  PlantNCreate%ispec_plantnh4uptake = 0
  PlantNCreate%ispec_plantno3uptake = 0
  PlantNCreate%rate_plantndemand = 0.d0
  PlantNCreate%half_saturation_nh4 = 1.d-15
  PlantNCreate%half_saturation_no3 = 1.d-15
  PlantNCreate%inhibition_nh4_no3  = 1.d0
  PlantNCreate%x0eps_nh4  = 1.d-20
  PlantNCreate%x0eps_no3  = 1.d-20
  nullify(PlantNCreate%next)
      
end function PlantNCreate

! ************************************************************************** !
!
! PlantNRead: Reads input deck for plantn reaction parameters
! author: Guoping Tang
! date: 09/09/2013
!
! ************************************************************************** !
subroutine PlantNRead(this,input,option)

  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module, only : UnitsConvertToInternal
  
  implicit none
  
  class(reaction_sandbox_plantn_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  PetscInt :: i
  character(len=MAXWORDLENGTH) :: word
  
  do 
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,REACTION_SANDBOX,PLANTN')
    call StringToUpper(word)   

    select case(trim(word))
      case('RATE_PLANTNDEMAND')
          call InputReadDouble(input,option,this%rate_plantndemand)
          call InputErrorMsg(input,option,'rate_plantndemand', &
                  'CHEMISTRY,REACTION_SANDBOX,PLANTN,REACTION')
      case('AMMONIUM_HALF_SATURATION')
          call InputReadDouble(input,option,this%half_saturation_nh4)
          call InputErrorMsg(input,option,'half saturation for ammonium', &
                     'CHEMISTRY,REACTION_SANDBOX,PLANTN,REACTION')
      case('NITRATE_HALF_SATURATION')
          call InputReadDouble(input,option,this%half_saturation_no3)
          call InputErrorMsg(input,option,'half saturation for nitrate', &
                     'CHEMISTRY,REACTION_SANDBOX,PLANTN,REACTION')
      case('AMMONIUM_INHIBITION_NITRATE')
          call InputReadDouble(input,option,this%inhibition_nh4_no3)
          call InputErrorMsg(input,option,'ammonium inhibition on nitrate', &
                     'CHEMISTRY,REACTION_SANDBOX,PLANTN,REACTION')
          if (this%inhibition_nh4_no3<0.d-20) then
            option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,PLANTN,' // &
              'AMMONIUM_INHIBITION_NITRATE cannot be too small to close to 0'
            call printErrMsg(option)
          endif
      case('X0EPS_NH4')
          call InputReadDouble(input,option,this%x0eps_nh4)
          call InputErrorMsg(input,option,'x0eps_nh4', &
                  'CHEMISTRY,REACTION_SANDBOX,PLANTN,REACTION')
      case('X0EPS_NO3')
          call InputReadDouble(input,option,this%x0eps_no3)
          call InputErrorMsg(input,option,'x0eps_no3', &
                  'CHEMISTRY,REACTION_SANDBOX,PLANTN,REACTION')
      case default
          option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,PLANTN,' // &
                'REACTION keyword: ' // trim(word) // ' not recognized.'
          call printErrMsg(option)
    end select
  enddo
  
end subroutine PlantNRead

! ************************************************************************** !
!
! PlantNSetup: Sets up the plantn reaction with parameters either
!                read from the input deck or hardwired.
! author: Guoping Tang
! date: 09/09/2013
!
! ************************************************************************** !
subroutine PlantNSetup(this,reaction,option)

  use Reaction_Aux_module
  use Option_module
  use Reaction_Immobile_Aux_module

  implicit none
  
  class(reaction_sandbox_plantn_type) :: this
  type(reaction_type) :: reaction
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word

!------------------------------------------------------------------------------------
  word = 'NH4+'
  this%ispec_nh4 = GetPrimarySpeciesIDFromName(word, reaction, PETSC_FALSE, option)

  word = 'NO3-'
  this%ispec_no3 = GetPrimarySpeciesIDFromName(word,reaction, PETSC_FALSE,option)

  if(this%ispec_nh4 < 0 .and. this%ispec_no3 < 0) then
     option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,PLANTN: ' // &
       ' at least one of NH4+ and NO3- must be specified as primary species in the input file.'
     call printErrMsg(option)
  endif

  word = 'PlantN'
  this%ispec_plantn = GetImmobileSpeciesIDFromName(word, reaction%immobile, &
                 PETSC_FALSE,option)

  if(this%ispec_plantn < 0) then
     option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,PLANTN: ' // &
       ' PlantN is not specified as immobile species in the input file.'
     call printErrMsg(option)
  endif

  word = 'Plantndemand'
  this%ispec_plantndemand = GetImmobileSpeciesIDFromName(word, reaction%immobile, &
                 PETSC_FALSE,option)

  word = 'Plantnh4uptake'
  this%ispec_plantnh4uptake = GetImmobileSpeciesIDFromName(word, reaction%immobile, &
                 PETSC_FALSE,option)
  word = 'Plantno3uptake'
  this%ispec_plantno3uptake = GetImmobileSpeciesIDFromName(word, reaction%immobile, &
                 PETSC_FALSE,option)

#ifdef CLM_PFLOTRAN
  if(this%ispec_plantnh4uptake < 0 .and. this%ispec_plantno3uptake < 0) then
     option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,PLANTN: ' // &
       'At least one of "Plantnh4uptake" or "Plantno3uptake " ' // &
       'must be specified as immobile species in the ' // &
       'input file, which required when coupled with CLM.'
     call printErrMsg(option)
  endif
#endif

#ifdef CLM_PFLOTRAN
  if(this%ispec_plantnh4uptake < 0 .and. this%ispec_plantno3uptake < 0) then
     option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,PLANTN: ' // &
       'At least one of "Plantnh4uptake" or "Plantno3uptake " ' // &
       'must be specified as immobile species in the ' // &
       'input file, which required when coupled with CLM.'
     call printErrMsg(option)
  endif
#endif

end subroutine PlantNSetup

! ************************************************************************** !
!
! PlantNReact: Evaluates reaction storing residual and/or Jacobian
! author: Guoping Tang
! date: 09/09/2013
!
! Rewritten by Fengming Yuan @Aug-14-2014. The orginal was messed-up with 'patches',
! which caused a lot of issues.
!
! ************************************************************************** !
subroutine PlantNReact(this,Residual,Jacobian,compute_derivative, &
                         rt_auxvar,global_auxvar,material_auxvar,reaction, &
                         option)

  use Option_module
  use Reaction_Aux_module
  use Reaction_Immobile_Aux_module
  use Material_Aux_class, only : material_auxvar_type
  use CLM_RspFuncs_module

#ifdef CLM_PFLOTRAN
#include "petsc/finclude/petscvec.h"
  use petscvec
  use clm_pflotran_interface_data
#endif
  
  implicit none

  class(reaction_sandbox_plantn_type) :: this
  type(option_type) :: option
  type(reaction_type) :: reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscBool :: compute_derivative

  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  PetscReal :: volume, porosity, saturation, tc
  PetscReal :: theta, L_water
  PetscInt :: ghosted_id
  PetscErrorCode :: ierr

  character(len=MAXWORDLENGTH) :: word

  PetscInt, parameter :: iphase = 1
  PetscInt :: ires_nh4, ires_no3, ires_plantn
  PetscInt :: ires_plantndemand, ires_plantnh4uptake, ires_plantno3uptake
  PetscInt :: ires

  PetscReal :: c_nh4         ! concentration (mole/m3)
  PetscReal :: fnh4          ! nh4 / (half_saturation + nh4): rate dependence on substrate
  PetscReal :: dfnh4_dnh4    ! d(fnh4)/d(nh4)

  PetscReal :: c_no3         ! concentration (mole/m3)
  PetscReal :: fno3          ! no3 / (half_saturation + no3): rate dependence on substrate
  PetscReal :: dfno3_dno3    ! d(fno3)/d(no3)

  PetscReal :: nratecap, dtmin  ! max. n uptake rate within allowable min. timestep
  PetscReal :: fnratecap        ! max. nratecap as function of c_nh4/c_no3 extracting rate vs. potential
  PetscReal :: dfnratecap_dnh4, dfnratecap_dno3

  ! nh4 inhibition on no3 uptake, or plant N uptake preference btw nh4 and no3
  ! (Currently it's similar function as microbial N immobilization)
  ! crate_nh4 = fnh4*fnh4_inhibit_no3, while crate_no3 = 1.-fnh4*fnh4_inhibition_no3
  ! by the following eq., if 'inhibition=1', uptake will be equal for both (if same conc.);
  ! higher coef, strong NH4 inhibition on NO3 (i.e., more NH4 uptake over NO3)
  PetscReal :: fnh4_inhibit_no3 ! inhibition_coef/(inhibition_coef + no3/nh4):
  PetscReal :: dfnh4_inhibit_no3_dnh4 ! d(fnh4_inhibit_no3)/dnh4
  PetscReal :: dfnh4_inhibit_no3_dno3 ! d(fnh4_inhibit_no3)/dno3

  PetscReal :: temp_real, feps0, dfeps0_dx

  PetscReal :: nrate_nh4
  PetscReal :: nrate_no3
  PetscReal :: dnrate_nh4_dnh4       !d(nrate_nh4)/d(nh4)
  PetscReal :: dnrate_nh4_dno3       !d(nrate_nh4)/d(no3)
  PetscReal :: dnrate_no3_dnh4       !d(nrate_no3)/d(nh4)
  PetscReal :: dnrate_no3_dno3       !d(nrate_no3)/d(no3)

#ifdef CLM_PFLOTRAN
  PetscScalar, pointer :: rate_plantndemand_pf_loc(:)   !
#endif 

  !--------------------------------------------------------------------------------------------
  volume = material_auxvar%volume
  porosity = material_auxvar%porosity
  saturation = global_auxvar%sat(1)
  if(saturation < 0.01d0) return

  theta = saturation * porosity
  L_water = theta * 1.0d3      ! Litres of H2O /m3 bulk soil

  tc = global_auxvar%temp
  if(tc < -0.1d0) return

  if(this%ispec_plantn>0) ires_plantn = this%ispec_plantn + reaction%offset_immobile
  if(this%ispec_nh4>0) ires_nh4 = this%ispec_nh4 + reaction%offset_aqueous
  if(this%ispec_no3>0) ires_no3 = this%ispec_no3 + reaction%offset_aqueous

  c_nh4 = 0.d0
  c_no3 = 0.d0
  fnh4       = 1.d0
  dfnh4_dnh4 = 0.d0
  fno3       = 1.d0
  dfno3_dno3 = 0.d0
  fnh4_inhibit_no3 = 1.d0
  dfnh4_inhibit_no3_dnh4 = 0.d0
  dfnh4_inhibit_no3_dno3 = 0.d0

  ! NOTE: the codes below are individually written for either NH4 or NO3 or both,
  !    so that, these two N species can be used separately for plant N uptake, or none, or both

  !--------------------------------------------------------------------------------------------
  !
  ! nh4 inhibition on no3 uptake, if any ('this%inhibition_nh4_no3')
  ! this is for quantifying plant N uptake preference btw NH4 and NO3
  ! (solely based on ratio of nh4/no3)
  if (this%ispec_nh4 >0 .and. this%ispec_no3 > 0) then
    c_nh4     = rt_auxvar%total(this%ispec_nh4, iphase) * L_water  ! mol/L (M) --> mole/m3
    c_no3     = rt_auxvar%total(this%ispec_no3, iphase) * L_water

    ! (DON'T change the 'rate' and 'derivatives' after this)
    if((c_nh4>this%x0eps_nh4 .and. c_no3>this%x0eps_no3) &
       .and. this%inhibition_nh4_no3>0.d0) then
       ! assuming that: f = c_nh4/(c_nh4+c_no3), or, f= (c_nh4/c_no3)/(c_nh4/c_no3+1)
       ! and adding 'preference kp - ratio of nh4/no3'
       ! f = kp*(c_nh4/c_no3)/(kp*(c_nh4/c_no3)+1), i.e.
       ! f = (c_nh4/c_no3)/(c_nh4/c_no3+1/kp)  - Monod type with half-saturation of 1/kp
      temp_real = c_nh4/c_no3
      fnh4_inhibit_no3 = funcMonod(temp_real, 1.0d0/this%inhibition_nh4_no3, PETSC_FALSE)

      ! the following appears troublesome (TODO - checking later on)
      ! symptoms: if both NH4 and NO3 available, the sum is NOT matching with single species of exactly same of total conc.
      !dfnh4_inhibit_no3_dnh4 = funcMonod(temp_real, 1.0d0/this%inhibition_nh4_no3, PETSC_TRUE)  ! over 'dtemp_real'
      !dfnh4_inhibit_no3_dnh4 = dfnh4_inhibit_no3_dnh4*(1.d0/c_no3)                              ! df_dtemp_real * dtemp_real_dnh4

      !dfnh4_inhibit_no3_dno3 = funcMonod(temp_real, 1.0d0/this%inhibition_nh4_no3, PETSC_TRUE)  ! over 'dtemp_real'
      !dfnh4_inhibit_no3_dno3 = dfnh4_inhibit_no3_dno3*(c_nh4/c_no3/c_no3)                       ! df_dtemp_real * dtemp_real_dno3

    else
      if (c_nh4>this%x0eps_nh4 .and. c_no3<=this%x0eps_no3) then
        fnh4_inhibit_no3 = 1.0d0
      elseif (c_nh4<=this%x0eps_nh4 .and. c_no3>this%x0eps_no3) then
        fnh4_inhibit_no3 = 0.0d0
      else
        return
        !fnh4_inhibit_no3 = 0.50d0
      endif
      dfnh4_inhibit_no3_dnh4 = 0.d0
      dfnh4_inhibit_no3_dno3 = 0.d0
    endif
    !
  endif

  !--------------------------------------------------------------------------------------------
  if (this%ispec_nh4 > 0) then
    c_nh4    = rt_auxvar%total(this%ispec_nh4, iphase) * L_water  ! mol/L (M) --> mole/m3

    ! half-saturation regulated rate (by input, it may be representing competetiveness)
    fnh4       = funcMonod(c_nh4, this%half_saturation_nh4, PETSC_FALSE)
    dfnh4_dnh4 = funcMonod(c_nh4, this%half_saturation_nh4, PETSC_TRUE)

    ! using the following for trailer smoothing
    ! note: 'x0eps' is different from 'half_saturation_nh4' above.
    !       'x0eps' is for mathematic reason in the code;
    !       'half_saturation_nh4' is for using 'monod-type' function to quantify plant
    !    NH4 uptake dependence on resources (NH4). So, physiologically it may be
    !    used as a method to quantify multiple consummer competition over resources.
    if(this%x0eps_nh4>0.d0) then
      ! GP's cut-off approach (sort of Heaviside function)
      call HfunctionSmooth(c_nh4, this%x0eps_nh4*10.d0, this%x0eps_nh4, feps0, dfeps0_dx)
    else
      feps0 = 1.0d0
      dfeps0_dx = 0.d0
    endif
    dfnh4_dnh4 = dfnh4_dnh4*feps0 + fnh4 * dfeps0_dx
    fnh4 = fnh4 * feps0

  endif

  !--------------------------------------------------------------------------------------------
  if (this%ispec_no3 > 0) then
    c_no3     = rt_auxvar%total(this%ispec_no3, iphase) * L_water       ! mol/L (M) --> mole/m3
    ! half-saturation regulated rate (by input, it may be representing competetiveness)
    fno3       = funcMonod(c_no3, this%half_saturation_no3, PETSC_FALSE)
    dfno3_dno3 = funcMonod(c_no3, this%half_saturation_no3, PETSC_TRUE)

    ! using the following for trailer smoothing
    ! note: 'x0eps' is different from 'half_saturation_no3' above.
    !       'x0eps' is for mathematic reason in the code;
    !       'half_saturation_no3' is for using 'monod-type' function to quantify plant
    !    NO3 uptake dependence on resources (NO3). So, physiologically it may be
    !    used as a method to quantify multiple consummer competition over resources.
    if(this%x0eps_no3>0.d0) then
      ! GP's cut-off approach (sort of Heaviside function)
      call HfunctionSmooth(c_no3, this%x0eps_no3*10.d0, this%x0eps_no3, feps0, dfeps0_dx)
    else
      feps0 = 1.0d0
      dfeps0_dx = 0.d0
    endif
    dfno3_dno3 = dfno3_dno3*feps0 + fno3 * dfeps0_dx
    fno3 = fno3 * feps0

  endif

  !--------------------------------------------------------------------------------------------
  ! plant N demand rates
#ifdef CLM_PFLOTRAN
  ghosted_id = option%iflag

  call VecGetArrayReadF90(clm_pf_idata%rate_plantndemand_pfs, &
       rate_plantndemand_pf_loc, ierr)

  this%rate_plantndemand = max(0.d0, rate_plantndemand_pf_loc(ghosted_id) * volume)          ! moles/m3/s * m3

  call VecRestoreArrayReadF90(clm_pf_idata%rate_plantndemand_pfs, &
       rate_plantndemand_pf_loc, ierr)

  if(this%rate_plantndemand<=0.d0) return

#else
!for testing
  this%rate_plantndemand = 1.d-2*volume

#endif

  if (this%ispec_plantndemand > 0) then  ! for tracking
    ires_plantndemand = this%ispec_plantndemand + reaction%offset_immobile
    Residual(ires_plantndemand) = Residual(ires_plantndemand) - this%rate_plantndemand
  endif

  ! constraining 'N demand rate' if too high compared to available within the allowable min. time-step
  ! It can be achieved by cutting time-step, but it may be taking a very small timestep finally
  ! - implying tiny timestep in model, which potentially crashes model
  if (this%rate_plantndemand > 0.d0) then

    ! in the following, '2.1' multiplier is chosen because that is slightly larger(5% to avoid numerical issue) than '2.0',
    ! which should be the previous-timestep before reaching the 'option%dt_min'
    ! (the time-cut used in PF is like dt=0.5*dt, when cutting)
    !dtmin = 2.1d0*option%dt_min
    dtmin = max(option%tran_dt, 2.1d0*option%dt_min)   ! this 'dtmin' may be accelerating the timing, but may not be appropriate to mulitple consummers

    if (this%ispec_nh4 > 0) then
       nratecap = this%rate_plantndemand * dtmin
       if (this%ispec_no3 > 0) nratecap = this%rate_plantndemand*fnh4_inhibit_no3*dtmin
       if (nratecap > c_nh4*volume) then
         ! assuming that monod type reduction used and 'nratecap' must be reduced to 'c_nh4', i.e.
         ! c_nh4/dt = nratecap/dt * (c_nh4/(c_nh4+c0))
         ! then, c0 = nratecap - c_nh4 (this is the 'half-saturation' term used above)
         fnratecap = funcMonod(c_nh4*volume, nratecap-c_nh4*volume, PETSC_FALSE)
         dfnratecap_dnh4 = funcMonod(c_nh4*volume, nratecap-c_nh4*volume, PETSC_TRUE)
       else
         fnratecap       = 1.d0
         dfnratecap_dnh4 = 0.d0
       endif
       ! modifying the 'fnh4' and 'dfnh4_dnh4' calculated above
       ! so that NO need to modify the major codes below
       dfnh4_dnh4 = dfnh4_dnh4*fnratecap + fnh4 * dfnratecap_dnh4   ! do the derivative first
       fnh4 = fnh4 * fnratecap
    endif
    !
    if (this%ispec_no3 > 0) then
       nratecap = this%rate_plantndemand * dtmin
       if (this%ispec_nh4 > 0) nratecap = this%rate_plantndemand*(1.-fnh4_inhibit_no3)*dtmin
       if (nratecap > c_no3*volume) then
         fnratecap = funcMonod(c_no3*volume, nratecap-c_no3*volume, PETSC_FALSE)
         dfnratecap_dno3 = funcMonod(c_no3*volume, nratecap-c_no3*volume, PETSC_TRUE)
       else
         fnratecap       = 1.d0
         dfnratecap_dno3 = 0.d0
       endif
       ! modifying the 'fno3' and 'dfno3_dno3' calculated above
       ! so that NO need to modify the major codes below
       dfno3_dno3 = dfno3_dno3*fnratecap + fno3 * dfnratecap_dno3
       fno3 = fno3 * fnratecap
    endif
    !
  endif

  !--------------------------------------------------------------------------------------------
  ! residuals and derivatives

  if(this%ispec_nh4 > 0) then

    ! rates
    nrate_nh4 = this%rate_plantndemand * fnh4
    if(this%ispec_no3 > 0) then
    ! splitting (fractioning) potential uptake rate by the 'fnh4_inhibition_no3' for NH4 uptake
      nrate_nh4 = this%rate_plantndemand * fnh4 * fnh4_inhibit_no3
    endif

    ! residuals
    Residual(ires_nh4) = Residual(ires_nh4) + nrate_nh4
    Residual(ires_plantn) = Residual(ires_plantn) - nrate_nh4

    if (this%ispec_plantnh4uptake>0) then   ! for tracking
      ires_plantnh4uptake = this%ispec_plantnh4uptake + reaction%offset_immobile
      Residual(ires_plantnh4uptake) = Residual(ires_plantnh4uptake) - nrate_nh4
    endif

    ! jacobians
    if(compute_derivative) then

      dnrate_nh4_dnh4 = this%rate_plantndemand * dfnh4_dnh4
      if(this%ispec_no3 > 0) then
        temp_real = fnh4 * dfnh4_inhibit_no3_dnh4 + &        ! d(fnh4*fnh4_inhibit_no3)/dnh4
                    fnh4_inhibit_no3 * dfnh4_dnh4
        dnrate_nh4_dnh4 = this%rate_plantndemand * temp_real

        temp_real = fnh4 * dfnh4_inhibit_no3_dno3            ! d(fnh4*fnh4_inhibit_no3)/dno3
        dnrate_nh4_dno3 = this%rate_plantndemand * temp_real
      endif

      Jacobian(ires_nh4,ires_nh4) = Jacobian(ires_nh4,ires_nh4) + &
        dnrate_nh4_dnh4 * &
        rt_auxvar%aqueous%dtotal(this%ispec_nh4,this%ispec_nh4,iphase)

      Jacobian(ires_plantn,ires_nh4) = Jacobian(ires_plantn,ires_nh4) - &
        dnrate_nh4_dnh4

#ifndef nojacobian_track_vars
      if (this%ispec_plantnh4uptake>0) then   ! for tracking
        Jacobian(ires_plantnh4uptake,ires_nh4) = Jacobian(ires_plantnh4uptake,ires_nh4) - &
          dnrate_nh4_dnh4
      endif
#endif

      !if(this%ispec_no3 > 0) then
      !  Jacobian(ires_nh4,ires_no3)=Jacobian(ires_nh4,ires_no3) - &       ! may need a checking of the sign (+/-) here
      !    dnrate_nh4_dno3 * &
      !    rt_auxvar%aqueous%dtotal(this%ispec_nh4,this%ispec_no3,iphase)  ! this actually is 0
      !endif

    endif ! if(compute_derivative)

  endif !if(this%ispec_nh4 > 0)

  !
  if(this%ispec_no3 > 0) then

    ! rates
    nrate_no3 = this%rate_plantndemand * fno3
    if(this%ispec_nh4 > 0) then
    ! splitting (fractioning) potential uptake rate by the rest of nrate_nh4,
    ! which adjusted by 'fnh4_inhibition_no3'
    ! i.e., 1.0-fnh4_inhibit_no3
      nrate_no3 = this%rate_plantndemand * fno3 * (1.0d0-fnh4_inhibit_no3)
    endif

    ! residuals
    Residual(ires_no3) = Residual(ires_no3) + nrate_no3
    Residual(ires_plantn) = Residual(ires_plantn) - nrate_no3

    if (this%ispec_plantno3uptake>0) then   ! for tracking
      ires_plantno3uptake = this%ispec_plantno3uptake + reaction%offset_immobile
      Residual(ires_plantno3uptake) = Residual(ires_plantno3uptake) - nrate_no3
    endif

    ! jacobians
    if (compute_derivative) then

      dnrate_no3_dno3 = this%rate_plantndemand * dfno3_dno3
      if(this%ispec_nh4 > 0) then
        temp_real = dfno3_dno3 * (1.d0-fnh4_inhibit_no3) + &
                    fno3 * (-1.0d0*dfnh4_inhibit_no3_dno3)
        dnrate_no3_dno3 = this%rate_plantndemand * temp_real

        temp_real = fno3 * (-1.0d0 * dfnh4_inhibit_no3_dnh4)                 ! 'dfno3_dnh4=0'
        dnrate_no3_dnh4 = this%rate_plantndemand * temp_real
      endif

      Jacobian(ires_no3,ires_no3) = Jacobian(ires_no3,ires_no3) + &
        dnrate_no3_dno3 * &
        rt_auxvar%aqueous%dtotal(this%ispec_no3,this%ispec_no3,iphase)

      Jacobian(ires_plantn,ires_no3) = Jacobian(ires_plantn,ires_no3) - &
        dnrate_no3_dno3

#ifndef nojacobian_track_vars
      if (this%ispec_plantno3uptake>0) then   ! for tracking
        Jacobian(ires_plantno3uptake,ires_no3) = Jacobian(ires_plantno3uptake,ires_no3) - &
          dnrate_no3_dno3
      endif
#endif

      !if(this%ispec_nh4 > 0) then
      !  Jacobian(ires_no3,ires_nh4)=Jacobian(ires_no3,ires_nh4) - &      ! may need a checking of sign (+/-) here
      !    dnrate_no3_dnh4 * &
      !    rt_auxvar%aqueous%dtotal(this%ispec_no3,this%ispec_nh4,iphase)  ! this actually is 0
      !endif

    endif


  endif

#ifdef DEBUG
  if (option%print_file_flag) then

    if(option%tran_dt<=option%dt_min) then
      if (this%rate_plantndemand*fnh4*fnh4_inhibit_no3*option%dt_min>=c_nh4 .or.    &
          this%rate_plantndemand*fno3*(1.d0-fnh4_inhibit_no3)*option%dt_min>=c_no3) then
        write(option%fid_out, *) '----------------------------------------------'
        write(option%fid_out, *) 'Reaction Sandbox: PLANT N UPTAKE'
        write(option%fid_out, *) 'dt=',option%tran_dt, ' dt_min=',option%dt_min
        write(option%fid_out, *) 'ghosted_id=',ghosted_id, &
          ' c_nh4=',c_nh4, ' c_no3=',c_no3, ' fnh4=',fnh4,' fno3=', fno3, &
          'uprate_nh4=',this%rate_plantndemand*fnh4*fnh4_inhibit_no3*option%dt, &
          'uprate_no3=',this%rate_plantndemand*fno3*(1.d0-fnh4_inhibit_no3)*option%dt
       endif
    endif

    do ires=1, reaction%ncomp
      temp_real = Residual(ires)

      if (abs(temp_real) > huge(temp_real)) then
        write(option%fid_out, *) 'infinity of Residual matrix checking at ires=', ires
        write(option%fid_out, *) 'Reaction Sandbox: PLANT N UPTAKE'
        option%io_buffer = ' checking infinity of Residuals matrix @ PlantNReact '
        call printErrMsg(option)
      endif
    enddo

    if (temp_real /= temp_real) then
      write(option%fid_out, *) 'NaN of Residual matrix checking at ires=', ires
      write(option%fid_out, *) 'Reaction Sandbox: PLANT N UPTAKE'
      option%io_buffer = ' checking NaN of Residuals matrix  @ PlantNReact '
      call printErrMsg(option)
    endif

  endif
#endif

end subroutine PlantNReact

! ************************************************************************** !
!
! PlantNDestroy: Destroys allocatable or pointer objects created in this
!                  module
! author: Guoping Tang
! date: 09/09/2013
!
! ************************************************************************** !
subroutine PlantNDestroy(this)

  implicit none
  
  class(reaction_sandbox_plantn_type) :: this

end subroutine PlantNDestroy

end module Reaction_Sandbox_PlantN_class
