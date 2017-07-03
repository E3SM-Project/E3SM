module Reaction_Sandbox_PNNL_N_class

  use Reaction_Sandbox_Base_class
  
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none
  
  private
  
#include "finclude/petscsys.h"

  PetscReal, parameter :: d_to_s = 1.d0 / 3600.d0 / 24.d0
  PetscReal, parameter :: mM_to_M = 1.d-3
  PetscReal, parameter :: mM_g_d_to_M_g_s = mM_to_M * d_to_s
  PetscReal, parameter :: Q_NO3 = 88.d0 * mM_g_d_to_M_g_s   ! mM/g/d (13 - 88)
  PetscReal, parameter :: Q_NO2 = 125.d0 * mM_g_d_to_M_g_s  ! mM/g/d
  PetscReal, parameter :: Q_NO = 1200.d0 * mM_g_d_to_M_g_s  ! mM/g/d
  PetscReal, parameter :: Q_N2O = 480.d0 * mM_g_d_to_M_g_s  ! mM/g/d
  PetscReal, parameter :: Q_O2 = 500.d0 * mM_g_d_to_M_g_s   ! mM/g/d
  PetscReal, parameter :: K_S_OC_NO3 = 0.22d0 * mM_to_M     ! mM (0.03 - 0.22)
  PetscReal, parameter :: K_S_OC_NO2 = 0.1d0 * mM_to_M      ! mM (0.001 - 0.1)
  PetscReal, parameter :: K_S_OC_NO = 0.004d0 * mM_to_M     ! mM (0.0007 - 0.004)
  PetscReal, parameter :: K_S_OC_N2O = 0.02d0 * mM_to_M     ! mM (0.0002 - 0.02)
  PetscReal, parameter :: K_S_OC_O2 = 10.d0 * mM_to_M       ! mM (0.01 - 10)
  PetscReal, parameter :: K_S_NO3 = 1.1d0 * mM_to_M         ! mM (0.002 - 1.1)
  PetscReal, parameter :: K_S_NO2 = 4.1d-3 * mM_to_M        ! mM
  PetscReal, parameter :: K_S_NO = 1.1d-5 * mM_to_M         ! mM
  PetscReal, parameter :: K_S_N2O = 2.5d-2 * mM_to_M        ! mM
  PetscReal, parameter :: K_S_O2 = 0.06d0 * mM_to_M         ! mM (0.001 - 0.06)
  PetscReal, parameter :: K_I_O2 = 1.d-3 / 32.d0 * mM_to_M  ! mg/L -> M
  PetscReal, parameter :: K_I_NO3 = 20.d0 * 1.d-3           ! mM (0.039 - 20)
  PetscReal, parameter :: F_E_NO3_NO2 = 0.999d0             ! 0.5 - 0.999
  PetscReal, parameter :: F_E_NO2_NO = 0.94d0
  PetscReal, parameter :: F_E_NO_N2O = 0.99d0
  PetscReal, parameter :: F_E_N2O_N2 = 0.99d0
  PetscReal, parameter :: F_E_O2 = 0.4d0
  PetscReal, parameter :: B = 0.04d0 * d_to_s               ! 1/d (0.005 - 0.04)
  
  PetscReal, parameter :: STOICH_CH2O = -0.25d0

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_pnnl_n_type
    PetscInt :: o2_id
    PetscInt :: no3_id
    PetscInt :: no2_id
    PetscInt :: no_id
    PetscInt :: n2o_id
    PetscInt :: n2_id
    PetscInt :: h_id
    PetscInt :: oc_id
    PetscInt :: biomass_id
    PetscInt :: co2_id
    PetscReal :: stoich_1_no3
    PetscReal :: stoich_1_h
    PetscReal :: stoich_1_biomass
    PetscReal :: stoich_1_no2
    PetscReal :: stoich_1_co2
    PetscReal :: stoich_2_no2
    PetscReal :: stoich_2_h
    PetscReal :: stoich_2_biomass
    PetscReal :: stoich_2_no
    PetscReal :: stoich_2_co2
    PetscReal :: stoich_3_no
    PetscReal :: stoich_3_biomass
    PetscReal :: stoich_3_n2o
    PetscReal :: stoich_3_co2
    PetscReal :: stoich_4_n2o
    PetscReal :: stoich_4_biomass
    PetscReal :: stoich_4_n2
    PetscReal :: stoich_4_co2
    PetscReal :: stoich_5_no3
    PetscReal :: stoich_5_o2
    PetscReal :: stoich_5_h
    PetscReal :: stoich_5_biomass
    PetscReal :: stoich_5_co2
    PetscInt :: nrxn
    PetscInt :: nrow(5)
    PetscInt :: ncol(5)
    PetscInt :: irow(6,5)
    PetscInt :: icol(4,5)
    PetscReal :: stoich_row(6,5)
  contains
    procedure, public :: ReadInput => PNNL_NRead
    procedure, public :: Setup => PNNL_NSetup
    procedure, public :: Evaluate => PNNL_NReact
    procedure, public :: Destroy => PNNL_NDestroy
  end type reaction_sandbox_pnnl_n_type

  public :: PNNL_NCreate

contains

! ************************************************************************** !

function PNNL_NCreate()
  ! 
  ! Allocates PNNL N reaction object.
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/01/15
  ! 

  implicit none
  
  class(reaction_sandbox_pnnl_n_type), pointer :: PNNL_NCreate

  allocate(PNNL_NCreate)
  PNNL_NCreate%o2_id = UNINITIALIZED_INTEGER
  PNNL_NCreate%no3_id = UNINITIALIZED_INTEGER
  PNNL_NCreate%no2_id = UNINITIALIZED_INTEGER
  PNNL_NCreate%no_id = UNINITIALIZED_INTEGER
  PNNL_NCreate%n2o_id = UNINITIALIZED_INTEGER
  PNNL_NCreate%n2_id = UNINITIALIZED_INTEGER
  PNNL_NCreate%h_id = UNINITIALIZED_INTEGER
  PNNL_NCreate%oc_id = UNINITIALIZED_INTEGER
  PNNL_NCreate%biomass_id = UNINITIALIZED_INTEGER
  PNNL_NCreate%co2_id = UNINITIALIZED_INTEGER
  PNNL_NCreate%stoich_1_no3 = UNINITIALIZED_DOUBLE
  PNNL_NCreate%stoich_1_h = UNINITIALIZED_DOUBLE
  PNNL_NCreate%stoich_1_biomass = UNINITIALIZED_DOUBLE
  PNNL_NCreate%stoich_1_no2 = UNINITIALIZED_DOUBLE
  PNNL_NCreate%stoich_1_co2 = UNINITIALIZED_DOUBLE
  PNNL_NCreate%stoich_2_no2 = UNINITIALIZED_DOUBLE
  PNNL_NCreate%stoich_2_h = UNINITIALIZED_DOUBLE
  PNNL_NCreate%stoich_2_biomass = UNINITIALIZED_DOUBLE
  PNNL_NCreate%stoich_2_no = UNINITIALIZED_DOUBLE
  PNNL_NCreate%stoich_2_co2 = UNINITIALIZED_DOUBLE
  PNNL_NCreate%stoich_3_no = UNINITIALIZED_DOUBLE
  PNNL_NCreate%stoich_3_biomass = UNINITIALIZED_DOUBLE
  PNNL_NCreate%stoich_3_n2o = UNINITIALIZED_DOUBLE
  PNNL_NCreate%stoich_3_co2 = UNINITIALIZED_DOUBLE
  PNNL_NCreate%stoich_4_n2o = UNINITIALIZED_DOUBLE
  PNNL_NCreate%stoich_4_biomass = UNINITIALIZED_DOUBLE
  PNNL_NCreate%stoich_4_n2 = UNINITIALIZED_DOUBLE
  PNNL_NCreate%stoich_4_co2 = UNINITIALIZED_DOUBLE
  PNNL_NCreate%stoich_5_no3 = UNINITIALIZED_DOUBLE
  PNNL_NCreate%stoich_5_o2 = UNINITIALIZED_DOUBLE
  PNNL_NCreate%stoich_5_h = UNINITIALIZED_DOUBLE
  PNNL_NCreate%stoich_5_biomass = UNINITIALIZED_DOUBLE
  PNNL_NCreate%stoich_5_co2 = UNINITIALIZED_DOUBLE
  PNNL_NCreate%nrxn = UNINITIALIZED_INTEGER
  PNNL_NCreate%nrow = UNINITIALIZED_INTEGER
  PNNL_NCreate%ncol = UNINITIALIZED_INTEGER
  PNNL_NCreate%irow = UNINITIALIZED_INTEGER
  PNNL_NCreate%icol = UNINITIALIZED_INTEGER
  PNNL_NCreate%stoich_row = UNINITIALIZED_DOUBLE

  nullify(PNNL_NCreate%next)  
      
end function PNNL_NCreate

! ************************************************************************** !

subroutine PNNL_NRead(this,input,option)
  ! 
  ! Reads input deck for PNNL N reaction parameters (if any)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/01/15
  ! 

  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module, only : UnitsConvertToInternal
  
  implicit none
  
  class(reaction_sandbox_pnnl_n_type) :: this
  type(input_type) :: input
  type(option_type) :: option

  PetscInt :: i
  character(len=MAXWORDLENGTH) :: word, internal_units
  
  do 
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,REACTION_SANDBOX,PNNL_N')
    call StringToUpper(word)   

    select case(trim(word))
    end select
  enddo
  
end subroutine PNNL_NRead

! ************************************************************************** !

subroutine PNNL_NSetup(this,reaction,option)
  ! 
  ! Sets up the PNNL N reaction either with parameters either
  ! read from the input deck or hardwired.
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/01/15
  ! 

  use Reaction_Aux_module, only : reaction_type, GetPrimarySpeciesIDFromName
  use Reaction_Immobile_Aux_module, only : GetImmobileSpeciesIDFromName
  use Option_module

  implicit none
  
  class(reaction_sandbox_pnnl_n_type) :: this
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: irxn

  word = 'O2(aq)'
  this%o2_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'NO3-'
  this%no3_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'NO2-'
  this%no2_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'NO'
  this%no_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'N2O'
  this%n2o_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'N2'
  this%n2_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'H+'
  this%h_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'CH2O'
  this%oc_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'C5H7O2N'
  this%biomass_id = &
    GetImmobileSpeciesIDFromName(word,reaction%immobile,option) + reaction%offset_immobile
  word = 'CO2(aq)'
  this%co2_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
    
  this%stoich_1_no3 = -1.d0*(1.d0/28.d0 + 13.d0/28.d0*F_E_NO3_NO2)
  this%stoich_1_h = -1.d0*(1.d0 - F_E_NO3_NO2)/28.d0
  this%stoich_1_biomass = -1.d0*this%stoich_1_h
  this%stoich_1_no2 = F_E_NO3_NO2/2.d0
  this%stoich_1_co2 = 1.d0/14.d0 + 5.d0/28.d0*F_E_NO3_NO2
  this%stoich_2_no2 = -1.d0*(1.d0/26.d0 + 25.d0/26.d0*F_E_NO2_NO)
  this%stoich_2_h = this%stoich_2_no2
  this%stoich_2_biomass = (1.d0 - F_E_NO2_NO)/26.d0
  this%stoich_2_no = F_E_NO2_NO
  this%stoich_2_co2 = 3.d0/52.d0 + 5.d0/26.d0*F_E_NO2_NO
  this%stoich_3_no = -1.d0*(1.d0/25.d0 + 24.d0/25.d0*F_E_NO_N2O)
  this%stoich_3_biomass = (1.d0 - F_E_NO_N2O)/25.d0
  this%stoich_3_n2o = F_E_NO_N2O/2.d0
  this%stoich_3_co2 = 1.d0/20.d0 + 1.d0/5.d0*F_E_NO_N2O
  this%stoich_4_n2o = -1.d0*(1.d0/48.d0 + 23.d0/48.d0*F_E_N2O_N2)
  this%stoich_4_biomass = (1.d0 - F_E_N2O_N2)/24.d0
  this%stoich_4_n2 = F_E_N2O_N2/2.d0
  this%stoich_4_co2 = 1.d0/24.d0 + 5.d0/24.d0*F_E_N2O_N2 
  this%stoich_5_no3 = -1.d0*(1.d0-F_E_O2)/28.d0
  this%stoich_5_o2 = -0.25d0*F_E_O2
  this%stoich_5_h = this%stoich_5_no3
  this%stoich_5_biomass = -1.d0*this%stoich_5_no3
  this%stoich_5_co2 = 1.d0/14.d0 + 5.d0/28.d0*F_E_O2  
  
  this%nrxn = 5
  ! NO3- -> NO2-
  irxn = 1
  this%nrow(irxn) = 6
  this%irow(1,irxn) = this%oc_id
  this%irow(2,irxn) = this%no3_id
  this%irow(3,irxn) = this%h_id
  this%irow(4,irxn) = this%biomass_id
  this%irow(5,irxn) = this%no2_id
  this%irow(6,irxn) = this%co2_id
  this%stoich_row(1,irxn) = STOICH_CH2O
  this%stoich_row(2,irxn) = this%stoich_1_no3
  this%stoich_row(3,irxn) = this%stoich_1_h
  this%stoich_row(4,irxn) = this%stoich_1_biomass
  this%stoich_row(5,irxn) = this%stoich_1_no2
  this%stoich_row(6,irxn) = this%stoich_1_co2
  this%ncol(irxn) = 4
  this%icol(1,irxn) = this%oc_id
  this%icol(2,irxn) = this%no3_id
  this%icol(3,irxn) = this%biomass_id
  this%icol(4,irxn) = this%o2_id
  ! NO2- -> NO
  irxn = 2
  this%nrow(irxn) = 6
  this%irow(1,irxn) = this%oc_id
  this%irow(2,irxn) = this%no2_id
  this%irow(3,irxn) = this%h_id
  this%irow(4,irxn) = this%biomass_id
  this%irow(5,irxn) = this%no_id
  this%irow(6,irxn) = this%co2_id
  this%stoich_row(1,irxn) = STOICH_CH2O
  this%stoich_row(2,irxn) = this%stoich_2_no2
  this%stoich_row(3,irxn) = this%stoich_2_h
  this%stoich_row(4,irxn) = this%stoich_2_biomass
  this%stoich_row(5,irxn) = this%stoich_2_no
  this%stoich_row(6,irxn) = this%stoich_2_co2
  this%ncol(irxn) = 4
  this%icol(1,irxn) = this%oc_id
  this%icol(2,irxn) = this%no2_id
  this%icol(3,irxn) = this%biomass_id
  this%icol(4,irxn) = this%no3_id
  ! NO -> N2O
  irxn = 3
  this%nrow(irxn) = 5
  this%irow(1,irxn) = this%oc_id
  this%irow(2,irxn) = this%no_id
  this%irow(3,irxn) = this%biomass_id
  this%irow(4,irxn) = this%n2o_id
  this%irow(5,irxn) = this%co2_id
  this%stoich_row(1,irxn) = STOICH_CH2O
  this%stoich_row(2,irxn) = this%stoich_3_no
  this%stoich_row(3,irxn) = this%stoich_3_biomass
  this%stoich_row(4,irxn) = this%stoich_3_n2o
  this%stoich_row(5,irxn) = this%stoich_3_co2
  this%ncol(irxn) = 3
  this%icol(1,irxn) = this%oc_id
  this%icol(2,irxn) = this%no_id
  this%icol(3,irxn) = this%biomass_id
  ! N2O -> N2
  irxn = 4  
  this%nrow(irxn) = 5
  this%irow(1,irxn) = this%oc_id
  this%irow(2,irxn) = this%n2o_id
  this%irow(3,irxn) = this%biomass_id
  this%irow(4,irxn) = this%n2_id
  this%irow(5,irxn) = this%co2_id
  this%stoich_row(1,irxn) = STOICH_CH2O
  this%stoich_row(2,irxn) = this%stoich_4_n2o
  this%stoich_row(3,irxn) = this%stoich_4_biomass
  this%stoich_row(4,irxn) = this%stoich_4_n2
  this%stoich_row(5,irxn) = this%stoich_4_co2
  this%ncol(irxn) = 3
  this%icol(1,irxn) = this%oc_id
  this%icol(2,irxn) = this%n2o_id
  this%icol(3,irxn) = this%biomass_id
  ! O2 -> H2O + CO2
  irxn = 5  
  this%nrow(irxn) = 6
  this%irow(1,irxn) = this%oc_id
  this%irow(2,irxn) = this%no3_id
  this%irow(3,irxn) = this%o2_id
  this%irow(4,irxn) = this%h_id
  this%irow(5,irxn) = this%biomass_id
  this%irow(6,irxn) = this%co2_id
  this%stoich_row(1,irxn) = STOICH_CH2O
  this%stoich_row(2,irxn) = this%stoich_5_no3
  this%stoich_row(3,irxn) = this%stoich_5_o2
  this%stoich_row(4,irxn) = this%stoich_5_h
  this%stoich_row(5,irxn) = this%stoich_5_biomass
  this%stoich_row(6,irxn) = this%stoich_5_co2
  this%ncol(irxn) = 3
  this%icol(1,irxn) = this%oc_id
  this%icol(2,irxn) = this%o2_id
  this%icol(3,irxn) = this%biomass_id
      
end subroutine PNNL_NSetup

! ************************************************************************** !

subroutine PNNL_NReact(this,Residual,Jacobian,compute_derivative, &
                         rt_auxvar,global_auxvar,material_auxvar,reaction, &
                         option)
  ! 
  ! Evaluates reaction storing residual and/or Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/01/15
  ! 

  use Option_module
  use Reaction_Aux_module
  use Material_Aux_class
  
  implicit none
  
  class(reaction_sandbox_pnnl_n_type) :: this  
  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscBool :: compute_derivative
  ! the following arrays must be declared after reaction
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar

  PetscInt, parameter :: iphase = 1
  PetscReal :: L_water
  
  PetscInt :: i, j, irxn
  PetscReal :: cpk

  PetscReal :: Co2, Cno3, Cno2, Cno, Cn2o, Cn2, Coc, Cco2, X
  PetscReal :: rate_no3_no2, rate_no2_no, rate_no_n2o, rate_n2o_n2, rate_o2
  PetscReal :: drate_no3_no2_oc, drate_no3_no2_no3, drate_no3_no2_o2, drate_no3_no2_biomass
  PetscReal :: drate_no2_no_oc, drate_no2_no_no2, drate_no2_no_no3, drate_no2_no_biomass
  PetscReal :: drate_no_n2o_oc, drate_no_n2o_no, drate_no_n2o_biomass
  PetscReal :: drate_n2o_n2_oc, drate_n2o_n2_n2o, drate_n2o_n2_biomass
  PetscReal :: drate_o2_oc, drate_o2_o2, drate_o2_biomass
  
  PetscReal :: rate(5), derivative_col(6,5)
  
  L_water = material_auxvar%porosity*global_auxvar%sat(iphase)* &
            material_auxvar%volume*1.d3
    
  Co2 = rt_auxvar%pri_molal(this%o2_id)*rt_auxvar%pri_act_coef(this%o2_id)
  Cno3 = rt_auxvar%pri_molal(this%no3_id)*rt_auxvar%pri_act_coef(this%no3_id)
  Cno2 = rt_auxvar%pri_molal(this%no2_id)*rt_auxvar%pri_act_coef(this%no2_id)
  Cno = rt_auxvar%pri_molal(this%no_id)*rt_auxvar%pri_act_coef(this%no_id)
  Cn2o = rt_auxvar%pri_molal(this%n2o_id)*rt_auxvar%pri_act_coef(this%n2o_id)
  Cn2 = rt_auxvar%pri_molal(this%n2_id)*rt_auxvar%pri_act_coef(this%n2_id)
  Coc = rt_auxvar%pri_molal(this%oc_id)*rt_auxvar%pri_act_coef(this%oc_id)
  X = rt_auxvar%immobile(this%biomass_id-reaction%offset_immobile)
  
  ! NO3- -> NO2-
  rate_no3_no2 = Q_NO3 * &
                 X * &
                 (Coc/(Coc+K_S_OC_NO3)) * &
                 (Cno3/(Cno3+K_S_NO3)) * &
                 (K_I_O2/(Co2+K_I_O2))
  ! NO2- -> NO
  rate_no2_no =  Q_NO2 * &
                  X * &
                 (Coc/(Coc+K_S_OC_NO2)) * &
                 (Cno2/(Cno2+K_S_NO2)) * &
                 (K_I_NO3/(Cno3+K_I_NO3))
                 
  ! NO -> N2O
  rate_no_n2o =  Q_NO * &
                 X * &
                 (Coc/(Coc+K_S_OC_NO)) * &
                 (Cno/(Cno+K_S_NO))  
                 
  ! N2O -> N2
  rate_n2o_n2 =  Q_N2O * &
                 X * &
                 (Coc/(Coc+K_S_OC_N2O)) * &
                 (Cn2o/(Cn2o+K_S_N2O))                 
                 
  ! O2 -> 
  rate_o2 =      Q_O2 * &
                 X * &
                 (Coc/(Coc+K_S_OC_O2)) * &
                 (Co2/(Co2+K_S_O2))                 
                 
  rate(1) = rate_no3_no2
  rate(2) = rate_no2_no
  rate(3) = rate_no_n2o
  rate(4) = rate_n2o_n2
  rate(5) = rate_o2
  
  do irxn = 1, this%nrxn
    do i = 1, this%nrow(irxn)
      Residual(this%irow(i,irxn)) = Residual(this%irow(i,irxn)) + &
        this%stoich_row(i,irxn) * rate(irxn)
    enddo
  enddo
  
  ! decay of biomass
  Residual(this%biomass_id) = Residual(this%biomass_id) + &
        -1.d0 * B * X
                 
  if (compute_derivative) then
    ! NO3- -> NO2-
    cpk = Coc+K_S_OC_NO3
    drate_no3_no2_oc = rate_no3_no2 / (Coc/cpk) * (cpk - 1.d0) / (cpk*cpk) * rt_auxvar%pri_act_coef(this%oc_id)
    cpk = Cno3+K_S_NO3
    drate_no3_no2_no3 = rate_no3_no2 / (Cno3/cpk) * (cpk - 1.d0) / (cpk*cpk) * rt_auxvar%pri_act_coef(this%no3_id)
    drate_no3_no2_biomass = rate_no3_no2 / X
    cpk = Co2+K_I_O2
    drate_no3_no2_o2 = rate_no3_no2 / (K_I_O2/cpk) * -1.d0 * K_I_O2 / (cpk*cpk) * rt_auxvar%pri_act_coef(this%o2_id)
    
    irxn = 1
    derivative_col(1,irxn) = drate_no3_no2_oc
    derivative_col(2,irxn) = drate_no3_no2_no3
    derivative_col(3,irxn) = drate_no3_no2_biomass
    derivative_col(4,irxn) = drate_no3_no2_o2

    ! NO2- -> NO
    cpk = Coc+K_S_OC_NO2
    drate_no2_no_oc = rate_no2_no / (Coc/cpk) * (cpk - 1.d0) / (cpk*cpk) * rt_auxvar%pri_act_coef(this%oc_id)
    cpk = Cno2+K_S_NO2
    drate_no2_no_no2 = rate_no2_no / (Cno2/cpk) * (cpk - 1.d0) / (cpk*cpk) * rt_auxvar%pri_act_coef(this%no2_id)
    drate_no2_no_biomass = rate_no2_no / X
    cpk = Cno3+K_I_NO3
    drate_no2_no_no3 = rate_no2_no / (K_I_NO3/cpk) * -1.d0 * K_I_NO3 / (cpk*cpk) * rt_auxvar%pri_act_coef(this%no3_id)
    
    irxn = 2
    derivative_col(1,irxn) = drate_no2_no_oc
    derivative_col(2,irxn) = drate_no2_no_no2
    derivative_col(3,irxn) = drate_no2_no_biomass
    derivative_col(4,irxn) = drate_no2_no_no3
    
    ! NO -> N2O
    cpk = Coc+K_S_OC_NO
    drate_no_n2o_oc = rate_no_n2o / (Coc/cpk) * (cpk - 1.d0) / (cpk*cpk) * rt_auxvar%pri_act_coef(this%oc_id)
    cpk = Cno+K_S_NO
    drate_no_n2o_no = rate_no_n2o / (Cno/cpk) * (cpk - 1.d0) / (cpk*cpk) * rt_auxvar%pri_act_coef(this%no_id)
    drate_no_n2o_biomass = rate_no_n2o / X
    
    irxn = 3
    derivative_col(1,irxn) = drate_no_n2o_oc
    derivative_col(2,irxn) = drate_no_n2o_no
    derivative_col(3,irxn) = drate_no_n2o_biomass
    
    ! N2O -> N2
    cpk = Coc+K_S_OC_N2O
    drate_n2o_n2_oc = rate_n2o_n2 / (Coc/cpk) * (cpk - 1.d0) / (cpk*cpk) * rt_auxvar%pri_act_coef(this%oc_id)
    cpk = Cno+K_S_NO
    drate_n2o_n2_n2o = rate_n2o_n2 / (Cn2o/cpk) * (cpk - 1.d0) / (cpk*cpk) * rt_auxvar%pri_act_coef(this%n2o_id)
    drate_n2o_n2_biomass = rate_n2o_n2 / X
    
    irxn = 4
    derivative_col(1,irxn) = drate_n2o_n2_oc
    derivative_col(2,irxn) = drate_n2o_n2_n2o
    derivative_col(3,irxn) = drate_n2o_n2_biomass
    
    ! O2 ->
    cpk = Coc+K_S_OC_O2
    drate_o2_oc = rate_o2 / (Coc/cpk) * (cpk - 1.d0) / (cpk*cpk) * rt_auxvar%pri_act_coef(this%oc_id)
    cpk = Cno+K_S_O2
    drate_o2_o2 = rate_o2 / (Co2/cpk) * (cpk - 1.d0) / (cpk*cpk) * rt_auxvar%pri_act_coef(this%o2_id)
    drate_o2_biomass = rate_o2 / X
    
    irxn = 5
    derivative_col(1,irxn) = drate_o2_oc
    derivative_col(2,irxn) = drate_o2_o2
    derivative_col(3,irxn) = drate_o2_biomass
    
    ! fill the Jacobian
    do irxn = 1, this%nrxn
      do j = 1, this%ncol(irxn)
        do i = 1, this%nrow(irxn)
          Jacobian(this%irow(i,irxn),this%icol(j,irxn)) = &
            Jacobian(this%irow(i,irxn),this%icol(j,irxn)) + &
            this%stoich_row(i,irxn) * derivative_col(j,irxn)
        enddo
      enddo
    enddo
    
  ! decay of biomass
  Jacobian(this%biomass_id,this%biomass_id) = &
    Jacobian(this%biomass_id,this%biomass_id) + -1.d0 * B
    
    

! 11. If using an analytical Jacobian, add code for Jacobian evaluation

    ! always add contribution to Jacobian
    ! units = (mol/sec)*(kg water/mol) = kg water/sec
  !  Jacobian(this%species_id,this%species_id) = &
  !  Jacobian(this%species_id,this%species_id) + &
  !    this%rate_constant * & ! 1/sec
  !    L_water * & ! L water
                  ! kg water/L water
      ! rt_auxvar%aqueous%dtotal(this%species_id,this%species_id,iphase) = 
      !   derivative of total component concentration with respect to the
      !   free ion concentration of the same species.
   !   rt_auxvar%aqueous%dtotal(this%species_id,this%species_id,iphase) 

  endif
  
end subroutine PNNL_NReact

! ************************************************************************** !

subroutine PNNL_NDestroy(this)
  ! 
  ! Destroys allocatable or pointer objects created in this
  ! module
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/01/15
  ! 

  implicit none
  
  class(reaction_sandbox_pnnl_n_type) :: this  

! 12. Add code to deallocate contents of the PNNL N object

end subroutine PNNL_NDestroy

end module Reaction_Sandbox_PNNL_N_class
