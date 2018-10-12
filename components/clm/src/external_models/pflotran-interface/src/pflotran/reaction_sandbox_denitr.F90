module Reaction_Sandbox_Denitr_class

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Reaction_Sandbox_Base_class
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  use PFLOTRAN_Constants_module
  use Utility_module, only : HFunctionSmooth
  
  implicit none
  
  private
  
  PetscInt, parameter :: TEMPERATURE_RESPONSE_FUNCTION_CLMCN = 1
  PetscInt, parameter :: TEMPERATURE_RESPONSE_FUNCTION_Q10   = 2

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_denitr_type
    PetscInt :: ispec_no3
    PetscInt :: ispec_n2
    PetscInt :: ispec_n2o
    PetscInt :: ispec_ngasdeni

    PetscReal :: half_saturation
    PetscInt  :: temperature_response_function
    PetscReal :: Q10
    PetscReal :: k_deni_max                     ! denitrification rate (1/second)
    PetscReal :: x0eps

  contains
    procedure, public :: ReadInput => DenitrRead
    procedure, public :: Setup => DenitrSetup
    procedure, public :: Evaluate => DenitrReact
    procedure, public :: Destroy => DenitrDestroy
  end type reaction_sandbox_denitr_type

  public :: DenitrCreate

contains

! ************************************************************************** !
!
! DenitrCreate: Allocates denitrification reaction object.
!
! ************************************************************************** !
function DenitrCreate()

  implicit none
  
  class(reaction_sandbox_denitr_type), pointer :: DenitrCreate

! Add code to allocate object and initialized all variables to zero and
!    nullify all pointers. E.g.
  allocate(DenitrCreate)
  DenitrCreate%ispec_no3 = 0
  DenitrCreate%ispec_n2o = 0
  DenitrCreate%ispec_n2 = 0
  DenitrCreate%ispec_ngasdeni = 0

  DenitrCreate%half_saturation = 1.0d-15
  DenitrCreate%temperature_response_function = TEMPERATURE_RESPONSE_FUNCTION_CLMCN
  DenitrCreate%Q10 = 1.5d0
  DenitrCreate%k_deni_max = 2.5d-6  ! max. denitrification rate (1/sec)
  DenitrCreate%x0eps = 1.0d-20

  nullify(DenitrCreate%next)
      
end function DenitrCreate

! ************************************************************************** !
!
! DenitrRead: Reads input deck for denitrification reaction parameters (if any)
!
! ************************************************************************** !
subroutine DenitrRead(this,input,option)

  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module, only : UnitsConvertToInternal
  
  implicit none
  
  class(reaction_sandbox_denitr_type) :: this
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
                       'CHEMISTRY,REACTION_SANDBOX,DENITRIFICATION')
    call StringToUpper(word)   

    select case(trim(word))
      case('TEMPERATURE_RESPONSE_FUNCTION')
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit

          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword', &
            'CHEMISTRY,REACTION_SANDBOX,DENITRIFICATION,TEMPERATURE RESPONSE FUNCTION')
          call StringToUpper(word)   

          select case(trim(word))
            case('CLMCN')
              this%temperature_response_function = TEMPERATURE_RESPONSE_FUNCTION_CLMCN
            case('Q10')
              this%temperature_response_function = TEMPERATURE_RESPONSE_FUNCTION_Q10    
              call InputReadDouble(input,option,this%Q10)  
              call InputErrorMsg(input,option,'Q10', &
                'CHEMISTRY,REACTION_SANDBOX_DENITRIFICATION,TEMPERATURE RESPONSE FUNCTION')
            case default
              option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,DENITRIFICATION,TEMPERATURE RESPONSE FUNCTION keyword: ' // &
                trim(word) // ' not recognized - Valid keyword: "CLMCN","Q10" '
              call printErrMsg(option)
          end select
        enddo 

      case('DENITRIFICATION_RATE_COEF')
        call InputReadDouble(input,option,this%k_deni_max)
        call InputErrorMsg(input,option,'k_deni_max', &
                 'CHEMISTRY,REACTION_SANDBOX,DENITRIFICATION,REACTION')
      case('NITRATE_HALF_SATURATION')
        call InputReadDouble(input,option,this%half_saturation)
        call InputErrorMsg(input,option,'nitrate half-saturation', &
                 'CHEMISTRY,REACTION_SANDBOX,DENITRIFICATION,REACTION')
      case('X0EPS')
        call InputReadDouble(input,option,this%x0eps)
        call InputErrorMsg(input,option,'x0eps', &
                  'CHEMISTRY,REACTION_SANDBOX,NITRIFICATION,REACTION')
      case default
        option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,DENITRIFICATION,' // &
          'REACTION keyword: ' // trim(word) // ' not recognized.'
        call printErrMsg(option)
    end select
  enddo
  
end subroutine DenitrRead

! ************************************************************************** !
!
! DenitrSetup: Sets up the denitrification reaction either with parameters either
!                read from the input deck or hardwired.
!
! ************************************************************************** !
subroutine DenitrSetup(this,reaction,option)

  use Reaction_Aux_module, only : reaction_type, GetPrimarySpeciesIDFromName
  use Option_module
  use Reaction_Immobile_Aux_module, only : GetImmobileSpeciesIDFromName

  implicit none
  
  class(reaction_sandbox_denitr_type) :: this
  type(reaction_type) :: reaction
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word
 
  word = 'NO3-'
  this%ispec_no3 = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)
  if(this%ispec_no3 < 0) then
     option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,DENITRIFICATION: ' // &
                        ' NO3- is not specified in the input file.'
     call printErrMsg(option)
  endif

  word = 'N2O(aq)'
  this%ispec_n2o = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)
  if(this%ispec_n2o < 0) then
     option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,DENITRIFICATION: ' // &
                        ' N2O(aq) is not specified in the input file.'
     call printErrMsg(option)
  endif

  word = 'N2(aq)'
  this%ispec_n2 = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)

  if(this%ispec_n2 < 0) then
     option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,DENITRIFICATION: ' // &
                        ' N2(aq) is not specified in the input file.'
     call printErrMsg(option)
  endif

  word = 'NGASdeni'
  this%ispec_ngasdeni = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)
#ifdef CLM_PFLOTRAN
  if(this%ispec_ngasdeni < 0) then
     option%io_buffer = 'CHEMISTRY,REACTION_SANDBOX,DENITRIFICATION: ' // &
       ' NGASdeni is not specified as immobile species in the input file. ' // &
       ' It is required when coupled with CLM.'
     call printErrMsg(option)
  endif
#endif
 
end subroutine DenitrSetup

!****************************************************************************************!
subroutine DenitrReact(this,Residual,Jacobian,compute_derivative, &
                         rt_auxvar,global_auxvar,material_auxvar,reaction, &
                         option)
!
!---------------------------------------------------------------------------------------!

  use Option_module
  use Reaction_Aux_module
  use Material_Aux_class, only : material_auxvar_type
  use CLM_RspFuncs_module

#ifdef CLM_PFLOTRAN
#include "petsc/finclude/petscvec.h"
  use petscvec
  use clm_pflotran_interface_data
#endif

  implicit none

  class(reaction_sandbox_denitr_type) :: this
  type(option_type) :: option
  type(reaction_type) :: reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar

  PetscBool :: compute_derivative
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  PetscReal :: porosity
  PetscReal :: volume   ! m3 bulk
  PetscReal :: L_water  ! (kg)L/m3 bulk
  PetscInt :: ghosted_id
  PetscErrorCode :: ierr

  PetscReal :: temp_real

  PetscInt :: ires_no3, ires_n2, ires_n2o
  PetscInt :: ires_ngasdeni

  PetscScalar, pointer :: bsw(:)
  PetscScalar, pointer :: bulkdensity(:)

  PetscReal :: s_min
  PetscReal :: tc
  PetscReal :: f_t, f_w

  PetscReal :: c_no3         ! moles/L ==> moles/m3 bulk soil
  PetscReal :: fno3          ! c_no3 / (half_saturation + c_no3)
  PetscReal :: dfno3_dno3    ! d(fno3)/d(c_no3)
  PetscReal :: rate_deni, drate_deni_dno3
  PetscReal :: saturation
  PetscInt, parameter :: iphase = 1

! misc. local variables
  PetscInt :: ires
  PetscReal:: feps0, dfeps0_dx

!---------------------------------------------------------------------------------
  if(this%ispec_n2 < 0) return

  porosity = material_auxvar%porosity
  volume = material_auxvar%volume
  saturation = global_auxvar%sat(iphase)
  L_water = porosity * saturation * 1.d3

  tc = global_auxvar%temp

!---------------------------------------------------------------------------------
  ! indices for C and N species
  ires_no3 = this%ispec_no3 + reaction%offset_aqueous    ! as aq. species
  ires_n2o = this%ispec_n2o + reaction%offset_aqueous
  ires_n2  = this%ispec_n2 + reaction%offset_aqueous
  ires_ngasdeni = this%ispec_ngasdeni + reaction%offset_immobile

#ifdef CLM_PFLOTRAN
  ghosted_id = option%iflag
  call VecGetArrayReadF90(clm_pf_idata%bsw_pfs, bsw, ierr)
  CHKERRQ(ierr)
  temp_real = bsw(ghosted_id)
  call VecRestoreArrayReadF90(clm_pf_idata%bsw_pfs, bsw, ierr)
  CHKERRQ(ierr)

#else
  temp_real = 1.0d0
#endif

!---------------------------------------------------------------------------------
! denitrification (Dickinson et al. 2002)

  ! temperature response function
  f_t = exp(0.08d0 * (tc - 25.d0))

  ! moisture response function
  s_min = 0.6d0
  f_w = 0.d0
  if(saturation > s_min) then
     f_w = (saturation - s_min)/(1.0d0 - s_min)
     f_w = f_w ** temp_real
  endif

  c_no3 = rt_auxvar%total(ires_no3, iphase)*L_water         ! mol/Lw -> moles/m3 bulk
  if(this%x0eps>0.d0) then
    ! GP's cut-off approach (sort of Heaviside function)
    call HfunctionSmooth(c_no3, this%x0eps*10.d0, this%x0eps, feps0, dfeps0_dx)

  else
    feps0 = 1.d0
    dfeps0_dx = 0.d0
    if(c_no3 <= this%x0eps) return     ! this may bring in 'oscillation' around 'this%x0eps'
  endif

  ! rate dependence on substrate
  if (this%half_saturation > 0.0d0) then
    fno3      = funcMonod(c_no3, this%half_saturation, PETSC_FALSE)
    dfno3_dno3= funcMonod(c_no3, this%half_saturation, PETSC_TRUE)
  else
    fno3      = 1.0d0
    dfno3_dno3= 0.d0
  endif

  rate_deni = 0.d0
  drate_deni_dno3 = 0.d0
  if(f_t > 0.d0 .and. f_w > 0.d0) then
     ! unit: moles/second - 1/second * - * - * - * (moles/m3*m3*-)
     rate_deni = this%k_deni_max * f_t * f_w * fno3 * (c_no3*volume*feps0)

     Residual(ires_no3) = Residual(ires_no3) + rate_deni

     Residual(ires_n2) = Residual(ires_n2) - 0.5d0*rate_deni
     !(TODO) currently not separate denitrification gas into 'n2' and 'n2o', but needed soon.

     if(this%ispec_ngasdeni > 0) then
        Residual(ires_ngasdeni) = Residual(ires_ngasdeni) - rate_deni
     endif

    if (compute_derivative) then
      temp_real = dfno3_dno3 * (c_no3*volume*feps0) + &
                  fno3 * (c_no3*volume*dfeps0_dx + feps0)
      drate_deni_dno3 = this%k_deni_max * f_t * f_w * temp_real

      Jacobian(ires_no3,ires_no3) = Jacobian(ires_no3,ires_no3) +  &
        drate_deni_dno3 * &
        rt_auxvar%aqueous%dtotal(this%ispec_no3,this%ispec_no3,iphase)

      Jacobian(ires_n2,ires_no3) = Jacobian(ires_n2,ires_no3) - &
        0.5d0*drate_deni_dno3 * &
        rt_auxvar%aqueous%dtotal(this%ispec_n2,this%ispec_no3,iphase)

#ifndef nojacobian_track_vars
      ! for tracking
      if(this%ispec_ngasdeni > 0) then
        Jacobian(ires_ngasdeni,ires_no3) = &
          Jacobian(ires_ngasdeni,ires_no3) - drate_deni_dno3
      endif
#endif

    endif

  endif

#ifdef DEBUG
  if( (option%tran_dt<=option%dt_min .and. option%print_file_flag) .and. &
     rate_deni*option%dt_min >= c_no3) then

    write(option%fid_out, *) '----------------------------------------------'
    write(option%fid_out, *) 'Reaction Sandbox: DENITRIFICATION'
    write(option%fid_out, *) 'dt=',option%tran_dt, ' dt_min=',option%dt_min
    write(option%fid_out, *) 'ghosted_id=',ghosted_id, ' c_no3=',c_no3, &
    ' ratedt_denitri=',rate_deni*option%dt, ' drate_dno3=',drate_deni_dno3
  endif

  do ires=1, reaction%ncomp
    temp_real = Residual(ires)

    if (abs(temp_real) > huge(temp_real)) then
      write(option%fid_out, *) 'infinity of Residual matrix checking at ires=', ires
      write(option%fid_out, *) 'Reaction Sandbox: DENITRIFICATION'
      option%io_buffer = ' checking infinity of Residuals matrix @ DenitrReact '
      call printErrMsg(option)
    endif

    if (temp_real /= temp_real) then
      write(option%fid_out, *) 'NaN of Residual matrix checking at ires=', ires
      write(option%fid_out, *) 'Reaction Sandbox: DENITRIFICATION'
      option%io_buffer = ' checking NaN of Residuals matrix  @ DenitrReact '
      call printErrMsg(option)
    endif

  enddo
#endif

end subroutine DenitrReact

! ************************************************************************** !
!
! DenitrDestroy: Destroys allocatable or pointer objects created in this
!                  module
!
! ************************************************************************** !
subroutine DenitrDestroy(this)

  implicit none
  
  class(reaction_sandbox_denitr_type) :: this

end subroutine DenitrDestroy

end module Reaction_Sandbox_Denitr_class
