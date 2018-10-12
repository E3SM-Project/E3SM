module CLM_Rxn_Base_class
  
  ! extended from reaction_sandbox_base to implement demand based 
  ! down regulation for use in CLM_Rxn t6g 10/06/2014 

  use PFLOTRAN_Constants_module

  implicit none
  
  private
  
#include "petsc/finclude/petscsys.h"

  type, abstract, public :: clm_rxn_base_type
    class(clm_rxn_base_type), pointer :: next
  contains
#if 0  
    procedure(Base_Read), public, deferred :: ReadInput
    procedure(Base_Setup), public, deferred :: Setup 
    procedure(Base_React), public, deferred :: Evaluate
    procedure(Base_Destroy), public, deferred :: Destroy
#else
    procedure, public :: ReadInput => Base_Read
    procedure, public :: Setup => Base_Setup
    procedure, public :: Evaluate => Base_React
    procedure, public :: Destroy => Base_Destroy    
#endif
  end type clm_rxn_base_type
  
! for some reason cannot use the interfaces when passing in "this"
! with Intel
#if 0 
  abstract interface
  
    subroutine Base_Setup(this,reaction,option)
    
      use Option_module
      use Reaction_Aux_module
  
      import clm_rxn_base_type
    
      implicit none
  
      class(clm_rxn_base_type) :: this
      type(reaction_type) :: reaction
      type(option_type) :: option
  
    end subroutine Base_Setup 

    subroutine Base_Read(this,input,option)
    
      use Option_module
      use Input_Aux_module
  
      import clm_rxn_base_type
    
      implicit none
  
      class(clm_rxn_base_type) :: this
      type(input_type), pointer :: input
      type(option_type) :: option
  
    end subroutine Base_Read 
    
    subroutine Base_SkipBlock(this,input,option)
    
      use Option_module
      use Input_Aux_module
  
      import clm_rxn_base_type
    
      implicit none
  
      class(clm_rxn_base_type) :: this
      type(input_type), pointer :: input
      type(option_type) :: option
  
    end subroutine Base_SkipBlock 
    
    subroutine Base_React(this,Res,Jac,compute_derivative,rt_auxvar, &
                          global_auxvar,material_auxvar,reaction,option, &
                          RateDemand_nh4,RateSupply_nh4, &
                          JacobianDemand_nh4,JacobianSupply_nh4, &
                          RateDemand_no3,RateSupply_no3, &
                          JacobianDemand_no3,JacobianSupply_no3, &
                          Rate_nh4_to_no3,Jacobian_nh4_to_no3)

      use Option_module
      use Reaction_Aux_module
      use Reactive_Transport_Aux_module
      use Global_Aux_module
      use Material_Aux_class
  
      import clm_rxn_base_type
    
      implicit none
  
      class(clm_rxn_base_type) :: this
      type(option_type) :: option
      type(reaction_type) :: reaction
      PetscBool :: compute_derivative
      PetscReal :: Res(reaction%ncomp)
      PetscReal :: Jac(reaction%ncomp,reaction%ncomp)
      PetscReal :: RateDemand_nh4(reaction%ncomp)
      PetscReal :: RateSupply_nh4(reaction%ncomp)
      PetscReal :: JacobianDemand_nh4(reaction%ncomp,reaction%ncomp)
      PetscReal :: JacobianSupply_nh4(reaction%ncomp,reaction%ncomp)
      PetscReal :: RateDemand_no3(reaction%ncomp)
      PetscReal :: RateSupply_no3(reaction%ncomp)
      PetscReal :: JacobianDemand_no3(reaction%ncomp,reaction%ncomp)
      PetscReal :: JacobianSupply_no3(reaction%ncomp,reaction%ncomp)
      PetscReal :: Rate_nh4_to_no3
      PetscReal :: Jacobian_nh4_to_no3(reaction%ncomp)
      type(reactive_transport_auxvar_type) :: rt_auxvar
      type(global_auxvar_type) :: global_auxvar
      class(material_auxvar_type) :: material_auxvar
      
    end subroutine
    
    subroutine Base_Destroy(this)

      import clm_rxn_base_type
    
      implicit none
  
      class(clm_rxn_base_type) :: this

    end subroutine Base_Destroy   
    
  end interface

#else

contains

! ************************************************************************** !

  subroutine Base_Setup(this,reaction,option)
    
    use Option_module
    use Reaction_Aux_module
  
    implicit none
  
    class(clm_rxn_base_type) :: this
    type(reaction_type) :: reaction
    type(option_type) :: option
  
  end subroutine Base_Setup 

! ************************************************************************** !

  subroutine Base_Read(this,input,option)
    
    use Option_module
    use Input_Aux_module
  
    implicit none
  
    class(clm_rxn_base_type) :: this
    type(input_type), pointer :: input
    type(option_type) :: option
  
  end subroutine Base_Read

! ************************************************************************** !

  subroutine Base_SkipBlock(this,input,option)
    
    use Option_module
    use Input_Aux_module
  
    implicit none
  
    class(clm_rxn_base_type) :: this
    type(input_type), pointer :: input
    type(option_type) :: option
  
  end subroutine Base_SkipBlock   

! ************************************************************************** !

  subroutine Base_React(this,Residual,Jacobian,compute_derivative,rt_auxvar, &
                        global_auxvar,material_auxvar,reaction,option, &
                        RateDemand_nh4,RateSupply_nh4, &
                        JacobianDemand_nh4,JacobianSupply_nh4, &
                        RateDemand_no3,RateSupply_no3, &
                        JacobianDemand_no3,JacobianSupply_no3, &
                        Rate_nh4_to_no3,Jacobian_nh4_to_no3)
    use Option_module
    use Reaction_Aux_module
    use Reactive_Transport_Aux_module
    use Global_Aux_module
    use Material_Aux_class
  
    implicit none
  
    class(clm_rxn_base_type) :: this
    type(option_type) :: option
    type(reaction_type) :: reaction
    PetscBool :: compute_derivative
    PetscReal :: Residual(reaction%ncomp)
    PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
    PetscReal :: RateDemand_nh4(reaction%ncomp)
    PetscReal :: RateSupply_nh4(reaction%ncomp)
    PetscReal :: RateDemand_no3(reaction%ncomp)
    PetscReal :: RateSupply_no3(reaction%ncomp)
    PetscReal :: JacobianDemand_nh4(reaction%ncomp,reaction%ncomp)
    PetscReal :: JacobianSupply_nh4(reaction%ncomp,reaction%ncomp)
    PetscReal :: JacobianDemand_no3(reaction%ncomp,reaction%ncomp)
    PetscReal :: JacobianSupply_no3(reaction%ncomp,reaction%ncomp)
    PetscReal :: Rate_nh4_to_no3
    PetscReal :: Jacobian_nh4_to_no3(reaction%ncomp)
    type(reactive_transport_auxvar_type) :: rt_auxvar
    type(global_auxvar_type) :: global_auxvar
    class(material_auxvar_type) :: material_auxvar
      
  end subroutine

! ************************************************************************** !

  subroutine Base_Destroy(this)

    implicit none
  
    class(clm_rxn_base_type) :: this

  end subroutine Base_Destroy  
#endif

end module CLM_Rxn_Base_class

module CLM_Rxn_Common_module

  implicit none
  
  private
  
#include "petsc/finclude/petscsys.h"

  public :: CalNLimitFunc

contains

subroutine CalNLimitFunc(c_n, ac_n, &
                         residual, half_saturation, &
                         cutoff_0, cutoff_1, f_n, d_n)

  PetscReal c_n, ac_n, residual, half_saturation, cutoff_0, cutoff_1, f_n, d_n
  PetscReal temp_real, regulator, dregulator, xxx, delta

  f_n = 1.0d0
  d_n = 0.0d0

  if (half_saturation >= 1.0d-20) then
    temp_real = (c_n - residual) * ac_n + half_saturation
    f_n       = (c_n - residual) * ac_n / temp_real 
    d_n       = ac_n * half_saturation / temp_real / temp_real
  endif    

  if (cutoff_0 > 0.0d0) then

    ! additional down regulation for N uptake / immobimization
    if (c_n <= cutoff_0) then
      regulator = 0.0d0
      dregulator = 0.0d0
    elseif (c_n >= cutoff_1 .or. cutoff_1 - cutoff_0 <= 1.0d-20) then
      regulator = 1.0d0
      dregulator = 0.0d0
    else
      xxx   = c_n - cutoff_0
      delta = cutoff_1 - cutoff_0
      regulator = 1.0d0 - (1.0d0 - xxx * xxx / delta / delta) ** 2
      dregulator = 4.0d0 * (1.0d0 - xxx * xxx / delta / delta) * xxx &
                 / delta / delta
    endif
    
    ! rate = rate_orginal * regulator
    ! drate = drate_original * regulator + rate_orginal * dregulator
    d_n = d_n * regulator + f_n * dregulator

    f_n = f_n * regulator

  endif

end subroutine CalNLimitFunc

end module CLM_Rxn_Common_module

module CLM_Rxn_Decomp_class

  use CLM_Rxn_Base_class
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  use PFLOTRAN_Constants_module

! ------------------------------------------------------------------------------
! Description
! extended from reaction_sandbox_clmdec to implement demand based down regulation
! for use in CLM_Rxn t6g 10/06/2014 

! following is a description for clm_dec
! to be used to implement CLM-CN, and CLM-Microbe decomposition reactions
! extended from clm_rxn_clm_cn  
! 1) pools can be either immobile or aqueous species (e.g., DOM, acetate-, )
! 2) separate N into NH3 (or NH4+) and NO3-; Must have NH3 or NH4+, NO3- is used 
!    if it is specified in the input file
! 3) include flexibilities to have multiple downstream pools, and variable 
!    respiration fraction as in CLM-Microbe; 
! 4) add residual concentrations for upstream pools, NH3, and NO3- to 
!    keep reactant concentrations above 0 (used if > 0); 
! 5) add shut off down regulation for NH3 and NO3- (used when the first > 0)
! 6) include NH3 oxidation in decomposition using Parton et al. 2001 (used when 
!    N2O(aq) is specified in the input file) 
! 7) add optional immobile species to track respiration, N mineralization, and 
!    immobilization
! Author: Guoping Tang
! Date:   07/08/14 
! -----------------------------------------------------------------------------

  implicit none
  
  private
  
#include "petsc/finclude/petscsys.h"

  PetscInt, parameter :: LITTER_DECOMP_CLMCN = 1 
  PetscInt, parameter :: LITTER_DECOMP_CLMMICROBE = 2 

                          ! 14.00674d0 / 12.011d0
  PetscReal, parameter :: CN_ratio_mass_to_mol = 1.16616d0

  ! Sinsabaugh et al. 2013 Ecology Letters, 16, 930-939
  PetscReal, parameter :: CN_ratio_microbe = 9.32928d0   ! 8.0d0 
  PetscReal, parameter :: CUE_max = 0.6d0

  type, public, &
    extends(clm_rxn_base_type) :: clm_rxn_clmdec_type

    PetscInt :: temperature_response_function
    PetscReal :: Q10
    PetscInt :: moisture_response_function

    PetscInt :: litter_decomp_type          ! CLM-CN or CLM-Microbe

    PetscReal :: half_saturation_nh4
    PetscReal :: half_saturation_no3
    PetscReal :: inhibition_nh4_no3
    PetscReal :: n2o_frac_mineralization     ! fraction of n2o from net N mineralization

    PetscReal :: residual_cpool
    PetscReal :: residual_nh4
    PetscReal :: residual_no3

    PetscReal :: cutoff_no3_0               ! shut off
    PetscReal :: cutoff_no3_1               ! start to decrease from 1
    PetscReal :: cutoff_nh4_0               ! shut off
    PetscReal :: cutoff_nh4_1               ! start to decrease from 1

    PetscReal :: net_n_min_rate_smooth_0     ! start from 0
    PetscReal :: net_n_min_rate_smooth_1     ! rise to 1

    PetscReal :: nc_bacteria
    PetscReal :: nc_fungi
    PetscReal :: fraction_bacteria

    PetscInt :: npool                        ! litter or variable CN ration pools
    PetscReal, pointer :: pool_nc_ratio(:)   ! NC ratio in mole  npool

    PetscInt :: nrxn
    PetscReal, pointer :: rate_constant(:)           ! nrxn

    PetscBool, pointer :: is_litter_decomp(:)        ! nrxn
    PetscInt,  pointer :: upstream_c_id(:)           ! nrxn
    PetscInt,  pointer :: upstream_n_id(:)           ! nrxn
    PetscReal, pointer :: upstream_nc(:)             ! nrxn
    PetscBool, pointer :: upstream_is_aqueous(:)     ! nrxn

    PetscInt,  pointer :: n_downstream_pools(:)      ! maximum # of downstream pools
    PetscInt,  pointer :: downstream_id(:,:)         ! nrxn by maximum # of downstream pools
    PetscBool, pointer :: downstream_is_aqueous(:,:) ! nrxn by maximum # of downstream pools
    PetscReal, pointer :: downstream_stoich(:,:)     ! nrxn by maximum # of downstream pools
    PetscReal, pointer :: downstream_nc(:,:)         ! nrxn by maximum # of downstream pools
    PetscReal, pointer :: mineral_c_stoich(:)        ! nrxn
    PetscReal, pointer :: mineral_n_stoich(:)        ! nrxn

    PetscInt :: species_id_co2
    PetscInt :: species_id_nh4
    PetscInt :: species_id_no3
    PetscInt :: species_id_n2o
    PetscInt :: species_id_dom
    PetscInt :: species_id_bacteria
    PetscInt :: species_id_fungi

    PetscInt :: species_id_hrimm
    PetscInt :: species_id_nmin
    PetscInt :: species_id_nimm
    PetscInt :: species_id_ngasmin
    PetscInt :: species_id_proton
    PetscBool :: bdebugoutput
    PetscBool :: bskipn2ojacobian
    PetscBool :: is_NH4_aqueous
    PetscBool :: is_NO3_aqueous

    type(pool_type), pointer :: pools
    type(clmdec_reaction_type), pointer :: reactions
  contains
    procedure, public :: ReadInput => CLMDec_Read
    procedure, public :: Setup => CLMDec_Setup
    procedure, public :: Evaluate => CLMDec_React
    procedure, public :: Destroy => CLMDec_Destroy
  end type clm_rxn_clmdec_type
  
  type :: pool_type
    character(len=MAXWORDLENGTH) :: name
    PetscReal :: stoich
    PetscReal :: nc_ratio
    type(pool_type), pointer :: next
  end type pool_type
  
  type :: clmdec_reaction_type
    character(len=MAXWORDLENGTH) :: upstream_pool_name
    type(pool_type), pointer :: downstream_pools
    PetscReal :: rate_constant
    type(clmdec_reaction_type), pointer :: next
  end type clmdec_reaction_type
  
  public :: CLMDec_Create

contains

! **************************************************************************** !

function CLMDec_Create()
  ! Allocates CLMDec reaction sandbox object.

  implicit none
  
  type(clm_rxn_clmdec_type), pointer :: CLMDec_Create
  
  allocate(CLMDec_Create)

  CLMDec_Create%Q10 = 1.5d0
  CLMDec_Create%litter_decomp_type=LITTER_DECOMP_CLMCN
  CLMDec_Create%half_saturation_nh4 =  1.0d-6
  CLMDec_Create%half_saturation_no3 =  1.0d-6
  CLMDec_Create%inhibition_nh4_no3 = -1.0d-15
  CLMDec_Create%n2o_frac_mineralization = 0.02d0  ! Parton et al. 2001

  CLMDec_Create%residual_cpool = 1.0d-20
  CLMDec_Create%residual_nh4 = 1.0d-10
  CLMDec_Create%residual_no3 = 1.0d-10

  CLMDec_Create%cutoff_no3_0 = -1.0d-9 
  CLMDec_Create%cutoff_no3_1 = 1.0d-7
  CLMDec_Create%cutoff_nh4_0 = -1.0d-9 
  CLMDec_Create%cutoff_nh4_1 = 1.0d-7

  CLMDec_Create%net_n_min_rate_smooth_0 = 0.0d0 
  CLMDec_Create%net_n_min_rate_smooth_1 = 1.0d-20

  CLMDec_Create%nc_bacteria = 0.17150d0

  ! CN_ratio_fungi = 17.4924d0     !15.0d0 ! or 10.0
  CLMDec_Create%nc_fungi = 0.05717d0

  CLMDec_Create%fraction_bacteria = 0.340927d0

  CLMDec_Create%npool = 0

  nullify(CLMDec_Create%pool_nc_ratio)

  CLMDec_Create%nrxn = 0
  nullify(CLMDec_Create%rate_constant)
  nullify(CLMDec_Create%is_litter_decomp)
  nullify(CLMDec_Create%upstream_c_id)
  nullify(CLMDec_Create%upstream_n_id)
  nullify(CLMDec_Create%upstream_nc)
  nullify(CLMDec_Create%upstream_is_aqueous)
  
  nullify(CLMDec_Create%n_downstream_pools)
  nullify(CLMDec_Create%downstream_id)
  nullify(CLMDec_Create%downstream_is_aqueous)
  nullify(CLMDec_Create%downstream_stoich)
  nullify(CLMDec_Create%mineral_c_stoich)
  nullify(CLMDec_Create%mineral_n_stoich)

  CLMDec_Create%species_id_co2 = 0
  CLMDec_Create%species_id_nh4 = 0
  CLMDec_Create%species_id_no3 = 0
  CLMDec_Create%species_id_n2o = 0
  CLMDec_Create%species_id_dom = 0
  CLMDec_Create%species_id_proton = 0
  CLMDec_Create%species_id_bacteria = 0
  CLMDec_Create%species_id_fungi = 0
  CLMDec_Create%species_id_hrimm = 0
  CLMDec_Create%species_id_nmin = 0
  CLMDec_Create%species_id_nimm = 0
  CLMDec_Create%species_id_ngasmin = 0

  CLMDec_Create%is_NH4_aqueous = PETSC_TRUE
  CLMDec_Create%is_NO3_aqueous = PETSC_TRUE
  CLMDec_Create%bdebugoutput = PETSC_FALSE

  nullify(CLMDec_Create%next)
  nullify(CLMDec_Create%pools)
  nullify(CLMDec_Create%reactions)

end function CLMDec_Create

! **************************************************************************** !

subroutine CLMDec_Read(this,input,option)
  ! 
  ! Reads input deck for reaction sandbox parameters
  ! 

  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module
  use Units_module, only : UnitsConvertToInternal
 
  implicit none
  
  class(clm_rxn_clmdec_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: word, internal_units
  
  type(pool_type), pointer :: new_pool, prev_pool
  type(pool_type), pointer :: new_pool_rxn, prev_pool_rxn
  type(clmdec_reaction_type), pointer :: new_reaction, prev_reaction
  
  PetscReal :: rate_constant, turnover_time
  PetscReal :: temp_real
  
  nullify(new_pool)
  nullify(prev_pool)

  nullify(new_pool_rxn)
  nullify(prev_pool_rxn)

  nullify(new_reaction)
  nullify(prev_reaction)
  
  do 
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword', &
      'CHEMISTRY,CLM_RXN,CLMDec')
    call StringToUpper(word)   
    select case(trim(word))

      case('CLM-MICROBE-LITTER-DECOMPOSITION')
        this%litter_decomp_type = LITTER_DECOMP_CLMMICROBE    

      case('RESIDUAL_CPOOL')
        call InputReadDouble(input,option,this%residual_cpool)
        call InputErrorMsg(input,option,'residual cpool', &
          'CHEMISTRY,CLM_RXN,CLMDec')

      case('RESIDUAL_NH4')
        call InputReadDouble(input,option,this%residual_nh4)
        call InputErrorMsg(input,option,'residual NH4+', &
          'CHEMISTRY,CLM_RXN,CLMDec')

      case('RESIDUAL_NO3')
        call InputReadDouble(input,option,this%residual_no3)
        call InputErrorMsg(input,option,'residual NO3-', &
          'CHEMISTRY,CLM_RXN,CLMDec')

      case('HALF_SATURATION_NH4')
        call InputReadDouble(input,option,this%half_saturation_nh4)
        call InputErrorMsg(input,option,'NH4 half saturation', &
          'CHEMISTRY,CLM_RXN,CLMDec')

      case('HALF_SATURATION_NO3')
        call InputReadDouble(input,option,this%half_saturation_no3)
        call InputErrorMsg(input,option,'NO3 half saturation', &
          'CHEMISTRY,CLM_RXN,CLMDec')

      case('CUTOFF_NH4')
        call InputReadDouble(input,option,this%cutoff_nh4_0)
        call InputErrorMsg(input,option,'cutoff_nh4_0', &
          'CHEMISTRY,CLM_RXN,CLMDec')
        call InputReadDouble(input,option,this%cutoff_nh4_1)
        call InputErrorMsg(input,option,'cutoff_nh4_1', &
          'CHEMISTRY,CLM_RXN,CLMDec')
        if (this%cutoff_nh4_0 > this%cutoff_nh4_1) then
          option%io_buffer = 'CHEMISTRY,CLM_RXN,CLMDec,' // &
            'NH4+ down regulation cut off concentration > concentration ' // &
            'where down regulation function = 1.'
          call printErrMsg(option)
        endif
      case('CUTOFF_NO3')
        call InputReadDouble(input,option,this%cutoff_no3_0)
        call InputErrorMsg(input,option,'cutoff_no3_0', &
          'CHEMISTRY,CLM_RXN,CLMDec,')
        call InputReadDouble(input,option,this%cutoff_no3_1)
        call InputErrorMsg(input,option,'cutoff_no3_1', &
          'CHEMISTRY,CLM_RXN,CLMDec,')
        if (this%cutoff_no3_0 > this%cutoff_no3_1) then
          option%io_buffer = 'CHEMISTRY,CLM_RXN,CLMDec' // &
            'NO3- down regulation cut off concentration > concentration ' // &
            'where down regulation function = 1.'
          call printErrMsg(option)
        endif

      case('SMOOTH_NET_N_MINERALIZATION')
        call InputReadDouble(input,option,this%net_n_min_rate_smooth_0)
        call InputErrorMsg(input,option,'net_n_min_rate_smooth_0', &
          'CHEMISTRY,CLM_RXN,CLMDec')
        call InputReadDouble(input,option,this%net_n_min_rate_smooth_1)
        call InputErrorMsg(input,option,'net_n_min_rate_smooth_1', &
          'CHEMISTRY,CLM_RXN,CLMDec')
        if (this%net_n_min_rate_smooth_0 > this%net_n_min_rate_smooth_1) then
          option%io_buffer = 'CHEMISTRY,CLM_RXN,CLMDec,' // &
            'Net N mineralization smooth 0 concentration > 1 concentration.'
          call printErrMsg(option)
        endif

     case('DEBUG_OUTPUT')
       this%bdebugoutput = PETSC_TRUE

     case('JACOBIAN_N2O_TRACKING_SKIP')
       this%bskipn2ojacobian = PETSC_TRUE

     case('NH4_INHIBITION_NO3')
       call InputReadDouble(input,option,this%inhibition_nh4_no3)
       call InputErrorMsg(input,option,'NH4 inhibition coefficient', &
         'CHEMISTRY,CLM_RXN,CLMDec')

     case('N2O_FRAC_MINERALIZATION')
       call InputReadDouble(input,option,this%n2o_frac_mineralization)
       call InputErrorMsg(input,option,'n2o fraction from mineralization', &
         'CHEMISTRY,CLM_RXN,CLMDec')

     case('POOLS')
       do
         call InputReadPflotranString(input,option)
         if (InputError(input)) exit
         if (InputCheckExit(input,option)) exit   

         allocate(new_pool)
         new_pool%name = ''
         new_pool%nc_ratio = -999.d0
         nullify(new_pool%next)

         call InputReadWord(input,option,new_pool%name,PETSC_TRUE)
         call InputErrorMsg(input,option,'pool name', &
           'CHEMISTRY,CLM_RXN,CLMDec,POOLS')
         call InputReadDouble(input,option,temp_real)
         if (InputError(input)) then
           new_pool%nc_ratio = -999.d0
         else
           ! convert CN ratio from mass C/mass N to mol C/mol N
           if (temp_real > 0.0d0 ) then
             new_pool%nc_ratio = 1.0d0/temp_real/CN_ratio_mass_to_mol
           endif
         endif

         if (associated(this%pools)) then
           prev_pool%next => new_pool
         else
           this%pools => new_pool
         endif
         prev_pool => new_pool
         nullify(new_pool)
       enddo

      case('REACTION')
      
        allocate(new_reaction)
        new_reaction%upstream_pool_name = ''
        new_reaction%rate_constant = -999.d0
        nullify(new_reaction%downstream_pools)
        nullify(new_reaction%next)
        
        ! need to set these temporarily in order to check that they
        ! are not both set.
        turnover_time = 0.d0
        rate_constant = 0.d0
        
        do 
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit

          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword', &
            'CHEMISTRY,CLM_RXN,CLMDec')
          call StringToUpper(word)   

          select case(trim(word))
            case('UPSTREAM_POOL')
              call InputReadWord(input,option, &
                new_reaction%upstream_pool_name,PETSC_TRUE)
              call InputErrorMsg(input,option,'upstream pool name', &
                'CHEMISTRY,CLM_RXN,CLMDec')
            case('DOWNSTREAM_POOL')
              allocate(new_pool_rxn)
              new_pool_rxn%name = ''
              new_pool_rxn%stoich = 0.d0
              nullify(new_pool_rxn%next)

              call InputReadWord(input,option, &
                 new_pool_rxn%name,PETSC_TRUE)
              call InputErrorMsg(input,option,'downstream pool name', &
                'CHEMISTRY,CLM_RXN,CLMDec')
              call InputReadDouble(input,option,new_pool_rxn%stoich)
              call InputErrorMsg(input,option,'Downstream pool stoich', &
                'CHEMISTRY,CLM_RXN,CLMDec' // &
                'TEMPERATURE RESPONSE FUNCTION')

              if (associated(new_reaction%downstream_pools)) then
                prev_pool_rxn%next => new_pool_rxn
              else
                new_reaction%downstream_pools => new_pool_rxn
              endif
              prev_pool_rxn => new_pool_rxn
              nullify(new_pool_rxn)

            case('RATE_CONSTANT')
              internal_units = 'mol/L-s|1/s|L/mol-s'
              call InputReadDouble(input,option,rate_constant)
              call InputErrorMsg(input,option,'rate constant', &
                     'CHEMISTRY,CLM_RXN,CLMDec,')
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (InputError(input)) then
                input%err_buf = 'CLMDec RATE CONSTANT UNITS'
                call InputDefaultMsg(input,option)
              else              
                rate_constant = rate_constant * &
                  UnitsConvertToInternal(word,internal_units,option)
              endif
            case('TURNOVER_TIME')
              internal_units = 'sec'
              call InputReadDouble(input,option,turnover_time)
              call InputErrorMsg(input,option,'turnover time', &
                     'CHEMISTRY,CLM_RXN,CLMDec')
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (InputError(input)) then
                input%err_buf = 'CLMDec TURNOVER TIME UNITS'
                call InputDefaultMsg(input,option)
              else              
                turnover_time = turnover_time * &
                  UnitsConvertToInternal(word,internal_units,option)
              endif
            case default
              call InputKeywordUnrecognized(word, &
                     'CHEMISTRY,CLM_RXN,CLMDec,REACTION',option)
          end select
        enddo
        
        ! check to ensure that one of turnover time or rate constant is set.
        if (turnover_time > 0.d0 .and. rate_constant > 0.d0) then
          option%io_buffer = 'Only TURNOVER_TIME or RATE_CONSTANT may ' // &
            'be included in a CLMDec reaction definition, but not both. ' // &
            'See reaction with upstream pool "' // &
            trim(new_reaction%upstream_pool_name) // '".'
          call printErrMsg(option)
        else if (turnover_time > 0.d0) then
          new_reaction%rate_constant = 1.d0 / turnover_time
        else
          new_reaction%rate_constant = rate_constant
        endif
        if (associated(this%reactions)) then
          prev_reaction%next => new_reaction
        else
          this%reactions => new_reaction
        endif
        prev_reaction => new_reaction
        nullify(new_reaction)        
      case default
        call InputKeywordUnrecognized(word,'CHEMISTRY,CLM_RXN,CLMDec',option)
    end select
  enddo
  
end subroutine CLMDec_Read

! **************************************************************************** !

subroutine CLMDec_Setup(this,reaction,option)
  ! 
  ! Sets up CLMDec reaction after it has been read from input
  ! 

  use Reaction_Aux_module
  use Option_module
  use String_module
  use Reaction_Immobile_Aux_module
  use Utility_module, only : DeallocateArray
  
  implicit none

  class(clm_rxn_clmdec_type) :: this
  type(option_type) :: option
  type(reaction_type) :: reaction
  
  character(len=MAXWORDLENGTH), allocatable :: pool_names(:)
  character(len=MAXWORDLENGTH) :: word
  
  PetscInt, pointer :: species_id_pool_c(:)
  PetscInt, pointer :: species_id_pool_n(:)
  PetscBool, pointer :: pool_is_aqueous(:)

  PetscInt :: icount, jcount, max_downstream_pools, ipool
  PetscReal :: stoich_c, stoich_n

  type(pool_type), pointer :: cur_pool
  type(clmdec_reaction_type), pointer :: cur_rxn
  
  ! count # pools
  icount = 0
  cur_pool => this%pools
  do
    if (.not.associated(cur_pool)) exit
    icount = icount + 1
    cur_pool => cur_pool%next
  enddo
  this%npool = icount
  
  ! count # reactions
  icount = 0
  cur_rxn => this%reactions
  do
    if (.not.associated(cur_rxn)) exit
    icount = icount + 1
    cur_rxn => cur_rxn%next
  enddo
  this%nrxn = icount
 
  allocate(this%n_downstream_pools(this%nrxn))
 
  ! count # downstream pools in each reaction
  max_downstream_pools = -1
  icount = 0
  cur_rxn => this%reactions
  do
    if (.not.associated(cur_rxn)) exit
    icount = icount + 1

    jcount = 0
    cur_pool => cur_rxn%downstream_pools

    do
      if (.not.associated(cur_pool)) exit
      jcount = jcount + 1
      cur_pool => cur_pool%next
    enddo

    this%n_downstream_pools(icount) = jcount

    if (max_downstream_pools < jcount) then
      max_downstream_pools = jcount
    endif 

    cur_rxn => cur_rxn%next
  enddo

  ! allocate and initialize arrays
  allocate(this%pool_nc_ratio(this%npool))

  allocate(this%rate_constant(this%nrxn))
  allocate(this%is_litter_decomp(this%nrxn))
  allocate(this%upstream_c_id(this%nrxn))
  allocate(this%upstream_n_id(this%nrxn))
  allocate(this%upstream_nc(this%nrxn))
  allocate(this%upstream_is_aqueous(this%nrxn))
  
  allocate(this%downstream_id(this%nrxn,max_downstream_pools))
  allocate(this%downstream_stoich(this%nrxn,max_downstream_pools))
  allocate(this%downstream_nc(this%nrxn,max_downstream_pools))
  allocate(this%downstream_is_aqueous(this%nrxn,max_downstream_pools))
  allocate(this%mineral_c_stoich(this%nrxn))
  allocate(this%mineral_n_stoich(this%nrxn))

  this%pool_nc_ratio = 0.d0
  this%rate_constant = 0.d0
  this%is_litter_decomp = PETSC_FALSE
  this%upstream_c_id = 0
  this%upstream_n_id = 0
  this%upstream_nc = -999.9
  this%upstream_is_aqueous = PETSC_FALSE

  this%downstream_id = 0
  this%downstream_is_aqueous = PETSC_FALSE
  this%downstream_stoich = 0.d0
  this%mineral_c_stoich = 0.d0
  this%mineral_n_stoich = 0.d0
  
  ! temporary array for mapping pools in reactions
  allocate(pool_names(this%npool))
  allocate(pool_is_aqueous(this%npool))
  allocate(species_id_pool_c(this%npool))
  allocate(species_id_pool_n(this%npool))

  pool_names = ''
  pool_is_aqueous = PETSC_FALSE
  species_id_pool_c = -999 
  species_id_pool_n = -999 

  ! pools
  icount = 0
  cur_pool => this%pools
  do
    if (.not.associated(cur_pool)) exit
    icount = icount + 1
    this%pool_nc_ratio(icount) = cur_pool%nc_ratio
    pool_names(icount) = cur_pool%name

    if (cur_pool%nc_ratio < 0.d0) then
      ! Since no CN ratio provided, must provide two species with the
      ! same name as the pool with C or N appended.
      word = trim(cur_pool%name) // 'C'
      species_id_pool_c(icount) = &
        GetImmobileSpeciesIDFromName(word,reaction%immobile, &
                                     PETSC_FALSE,option)
      word = trim(cur_pool%name) // 'N'
      species_id_pool_n(icount) = &
        GetImmobileSpeciesIDFromName(word,reaction%immobile, &
                                     PETSC_FALSE,option)
      if (species_id_pool_c(icount)<=0 .or. species_id_pool_n(icount)<=0) then
        option%io_buffer = 'For CLMDec pools with no CN ratio defined, ' // &
          'the user must define two immobile species with the same root ' // &
          'name as the pool with "C" or "N" appended, respectively.'
        call printErrMsg(option)
      endif
    else ! only one species (e.g. SOMX)
      species_id_pool_c(icount) = &
        GetImmobileSpeciesIDFromName(cur_pool%name,reaction%immobile, &
                                     PETSC_FALSE,option)
      if (species_id_pool_c(icount) <= 0) then
        species_id_pool_c(icount) = GetPrimarySpeciesIDFromName( &
          cur_pool%name, reaction, PETSC_FALSE,option)
        if (species_id_pool_c(icount) <= 0) then
          option%io_buffer = 'CLMDec pool: ' // cur_pool%name // 'is not ' // &
            'specified either in the IMMOBILE_SPECIES or PRIMARY_SPECIES!'
          call printErrMsg(option)
        else
          pool_is_aqueous(icount) = PETSC_TRUE
        endif
      endif
      
      if (StringCompare(cur_pool%name, 'Bacteria')) then
        this%nc_bacteria = cur_pool%nc_ratio
      endif 

      if (StringCompare(cur_pool%name, 'Fungi')) then
        this%nc_fungi = cur_pool%nc_ratio
      endif 
    endif
    cur_pool => cur_pool%next
  enddo
 
  ! reactions
  icount = 0
  cur_rxn => this%reactions
  do
    if (.not.associated(cur_rxn)) exit
    ! upstream pools
    icount = icount + 1
    ipool = StringFindEntryInList(cur_rxn%upstream_pool_name,pool_names)
    if (ipool == 0) then
      option%io_buffer = 'Upstream pool ' // &
        trim(cur_rxn%upstream_pool_name) // &
        'in reaction not found in list of pools.'
      call printErrMsg(option)
    else
      this%upstream_c_id(icount) = species_id_pool_c(ipool)
      this%upstream_n_id(icount) = species_id_pool_n(ipool)
      this%upstream_nc(icount) = this%pool_nc_ratio(ipool) 
      this%upstream_is_aqueous(icount) = pool_is_aqueous(ipool) 
      if (this%upstream_n_id(icount) > 0) then
        this%is_litter_decomp(icount) = PETSC_TRUE
      else
        if (this%upstream_nc(icount) < 0.0d0) then
          option%io_buffer = 'SOM decomp. reaction with upstream pool ' // &
            trim(cur_rxn%upstream_pool_name) // &
            'has negative C:N ratio in upstream pool.'
          call printErrMsg(option)
        endif
      endif
    endif

    ! downstream pools
    jcount = 0
    cur_pool => cur_rxn%downstream_pools

    do
      if (.not.associated(cur_pool)) exit
      jcount = jcount + 1

      if (len_trim(cur_pool%name) > 0) then
        ipool = StringFindEntryInList(cur_pool%name,pool_names)
        if (ipool == 0) then
          option%io_buffer = 'Downstream pool "' // trim(cur_pool%name) // &
            '" in reaction with upstream pool "' // &
            trim(cur_rxn%upstream_pool_name) // '" not found in list of pools.'
          call printErrMsg(option)
        else
          this%downstream_id(icount, jcount) = species_id_pool_c(ipool)
          this%downstream_stoich(icount, jcount) = cur_pool%stoich 
          this%downstream_nc(icount, jcount) = this%pool_nc_ratio(ipool) 
          this%downstream_is_aqueous(icount, jcount) = pool_is_aqueous(ipool) 

          if (this%downstream_nc(icount,jcount) < 0.d0) then
            option%io_buffer = 'For CLMDec reactions, downstream pools ' // &
              'must have a constant C:N ratio (i.e. C and N are not ' // &
              'tracked individually).  Therefore, pool "' // &
              trim(cur_pool%name) // &
             '" may not be used as a downstream pool.'
            call printErrMsg(option)
          endif
        endif
      endif

      cur_pool => cur_pool%next

    enddo

    this%rate_constant(icount) = cur_rxn%rate_constant
    cur_rxn => cur_rxn%next
  enddo 
  
  deallocate(pool_names)
  call DeallocateArray(pool_is_aqueous)
  call DeallocateArray(species_id_pool_c)
  call DeallocateArray(species_id_pool_n)

  ! set stoichiometric coefficients for som decomposition reactions  
  ! as they are constant due to fixed CN ratio
  do icount = 1, this%nrxn
    if (this%is_litter_decomp(icount)) then
      cycle
    else
      ! calculate respiration factor
      stoich_c = 1.0d0
      stoich_n = this%upstream_nc(icount)

      do jcount = 1, this%n_downstream_pools(icount)
        stoich_c = stoich_c - this%downstream_stoich(icount, jcount)
        stoich_n = stoich_n - this%downstream_stoich(icount, jcount) * &
                              this%downstream_nc(icount, jcount)
      enddo

      if (stoich_c < -1.0d-10) then
        option%io_buffer = 'CLMDec SOM decomposition reaction has negative' // &
          ' respiration fraction!'
        call printErrMsg(option)
      endif

      this%mineral_c_stoich(icount) = stoich_c
      this%mineral_n_stoich(icount) = stoich_n

     endif
  enddo

  word = 'HCO3-'
  this%species_id_co2 = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)

  if (this%species_id_co2 < 0) then
     word = 'CO2(aq)'
     this%species_id_co2 = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)
  endif

  if (this%species_id_co2 <= 0) then
    option%io_buffer = 'Neither HCO3- nor CO2(aq) is specified in the ' // &
      'input file for CLMDec!'
    call printErrMsg(option)
  endif

  word = 'NH4+'
  this%species_id_nh4 = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)

  if (this%species_id_nh4 < 0) then
    word = 'NH3(aq)'
    this%species_id_nh4 = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)
  endif

  if (this%species_id_nh4 < 0) then
    word = 'Ammonium'
    this%species_id_nh4 = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)
    if (this%species_id_nh4 > 0) then
      this%is_NH4_aqueous = PETSC_FALSE
    endif
  endif 

  if (this%species_id_nh4 <= 0) then
    option%io_buffer = 'NH4+, NH3(aq) or Ammonium is specified in the input' // &
      'file for CLMDec!'
    call printErrMsg(option)
  endif

  word = 'NO3-'
  this%species_id_no3 = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)
  
  if (this%species_id_no3 < 0) then
    word = 'Nitrate'
    this%species_id_no3 = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)
    if (this%species_id_no3 > 0) then
      this%is_NO3_aqueous = PETSC_FALSE
    endif
  endif 

  word = 'N2O(aq)'
  this%species_id_n2o = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)

  word = 'H+'
  this%species_id_proton = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)

  word = 'Bacteria'
  this%species_id_bacteria = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)

  word = 'Fungi'
  this%species_id_fungi = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)
 
  word = 'HRimm'
  this%species_id_hrimm = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)

  word = 'Nmin'
  this%species_id_nmin = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)
 
  word = 'Nimm'
  this%species_id_nimm = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)
 
  word = 'NGASmin'
  this%species_id_ngasmin = GetImmobileSpeciesIDFromName( &
    word,reaction%immobile,PETSC_FALSE,option)

  if (this%species_id_bacteria > 0 .and. this%species_id_fungi > 0 .and. & 
    this%nc_bacteria > 0.0d0 .and. this%nc_fungi > 0.0d0 ) then
    this%fraction_bacteria = (1.0d0/this%nc_bacteria) ** 0.6d0 / & 
      ((1.0d0/this%nc_bacteria) ** 0.6d0 + (1.0d0/this%nc_fungi) ** 0.6d0) 
  endif 

end subroutine CLMDec_Setup

! ************************************************************************** !
subroutine CLMDec_React(this,Residual,Jacobian,compute_derivative,rt_auxvar, &
                        global_auxvar,material_auxvar,reaction,option, &
                        RateDemand_nh4,RateSupply_nh4, &
                        JacobianDemand_nh4,JacobianSupply_nh4, &
                        RateDemand_no3,RateSupply_no3, &
                        JacobianDemand_no3,JacobianSupply_no3, &
                        Rate_nh4_to_no3,Jacobian_nh4_to_no3)
  ! 
  ! Evaluates reaction storing residual and/or Jacobian
  ! 

  use Option_module
  use Reaction_Aux_module
  use Material_Aux_class, only : material_auxvar_type
  use CLM_Rxn_Common_module, only: CalNLimitFunc

  implicit none

  class(clm_rxn_clmdec_type) :: this
  type(option_type) :: option
  type(reaction_type) :: reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar

  PetscBool :: compute_derivative
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  PetscReal :: RateDemand_nh4(reaction%ncomp)
  PetscReal :: RateSupply_nh4(reaction%ncomp)
  PetscReal :: JacobianDemand_nh4(reaction%ncomp,reaction%ncomp)
  PetscReal :: JacobianSupply_nh4(reaction%ncomp,reaction%ncomp)
  PetscReal :: RateDemand_no3(reaction%ncomp)
  PetscReal :: RateSupply_no3(reaction%ncomp)
  PetscReal :: JacobianDemand_no3(reaction%ncomp,reaction%ncomp)
  PetscReal :: JacobianSupply_no3(reaction%ncomp,reaction%ncomp)
  PetscReal :: Rate_nh4_to_no3
  PetscReal :: Jacobian_nh4_to_no3(reaction%ncomp)
  PetscReal :: porosity
  PetscReal :: volume
  PetscReal :: saturation
  PetscErrorCode :: ierr
  PetscInt :: local_id

  PetscReal :: theta
  PetscReal :: psi

  PetscReal :: c_nh4              ! concentration (mole/L)
  PetscReal :: ac_nh4             ! activity coefficient
  PetscReal :: f_nh4              ! nh4 / (half_saturation + nh4)
  PetscReal :: d_nh4              ! half_saturation / (half_saturation + nh4)^2
  PetscReal :: f_nh4_inhibit      ! inhibition_coef/(inhibition_coef + nh4)
  PetscReal :: d_nh4_inhibit_dnh4 ! -inhibition_coef/(inhibition_coef + nh4)^2

  PetscReal :: c_no3              ! concentration (mole/L)
  PetscReal :: ac_no3             ! activity coefficient 
  PetscReal :: f_no3              ! no3 / (half_saturation + no3)
  PetscReal :: d_no3              ! half_saturation/(no3 + half_saturation)^2 
  PetscReal :: temp_real

  PetscInt :: irxn
  PetscInt :: ipool_up, ipool_down
  PetscReal :: CN_ratio_up, CN_ratio_down
  PetscBool :: constant_CN_ratio_up
  PetscReal :: resp_frac

  PetscReal :: c_uc, c_un         ! upstream c pool, n pool concentration
  PetscReal :: ac_uc, ac_un       ! activity coefficient, if aqueous

  PetscInt :: ispec_uc, ispec_un, ispec_d   ! species id for upstream C, N, and downstream
  PetscInt :: ires_uc, ires_un, ires_d      ! id used for residual and Jacobian
  PetscInt :: ires_co2, ires_nh4, ires_n2o, ires_no3
  PetscInt :: ires_hrimm, ires_nmin, ires_nimm, ires_ngasmin
  PetscReal :: stoich_c, stoich_n

  PetscReal :: scaled_rate_const

  PetscReal :: rate_nh4        ! mole/s 
  PetscReal :: drate_nh4_duc   ! d Rate / d upstream c
  PetscReal :: drate_nh4_dnh4  ! d Rate / d nh4 ammonia limitation

  PetscReal :: Rdu_duc, Rdn_duc, Rdc_duc, Rdb_duc, Rdf_duc  ! u = Lit1N/Lit1C, c, b, f for CLM-Microbe
  PetscReal :: Rdu_dun, Rdn_dun, Rdc_dun, Rdb_dun, Rdf_dun
  PetscReal :: Rno3du_duc, Rno3dn_duc, Rno3dc_duc, Rno3db_duc, Rno3df_duc
  PetscReal :: Rno3du_dun, Rno3dn_dun, Rno3dc_dun, Rno3db_dun, Rno3df_dun

  ! for N immobilization reactions with NO3 as N source 
  PetscReal :: rate_no3       ! mole/s
  PetscReal :: drate_no3_dno3 ! d Rate_no3 / d no3 
  PetscReal :: drate_no3_duc  ! d Rate_no3 / d uc 
  PetscReal :: drate_no3_dnh4 ! d Rate_no3 / d nh4 

  PetscInt :: i, j
  PetscReal :: tc     ! temperature in C
  PetscReal :: f_t    ! temperature response function
  PetscReal :: f_w    ! moisture response function

  ! save mineral N fraction and decomposition rate for net N mineralization and N2O calculation 
  PetscReal :: net_n_mineralization_rate
  PetscReal :: dnet_n_mineralization_rate_dnh4
  PetscReal :: dnet_n_mineralization_rate_dno3
  PetscReal :: dnet_n_mineralization_rate_duc(this%nrxn)
  PetscReal :: ph, f_ph
  PetscReal :: rate_n2o, drate_n2o_dnh4, drate_n2o_dno3, drate_n2o_duc
  PetscReal :: f_rate_n2o, df_rate_n2o

  PetscInt :: ires_b, ires_f
  PetscReal :: xxx, delta, regulator, dregulator

  porosity = material_auxvar%porosity
  volume = material_auxvar%volume
  ires_nh4 = -999
  ires_no3 = -999
    
  if (this%is_NH4_aqueous) then   
    c_nh4    = rt_auxvar%pri_molal(this%species_id_nh4)
    ac_nh4   = rt_auxvar%pri_act_coef(this%species_id_nh4)
    ires_nh4 = this%species_id_nh4
  else
    c_nh4    = rt_auxvar%immobile(this%species_id_nh4)
    ac_nh4   = 1.0d0
    ires_nh4 = this%species_id_nh4 + reaction%offset_immobile
  endif

  call CalNLimitFunc(c_nh4, ac_nh4, this%residual_nh4, &
    this%half_saturation_nh4, this%cutoff_nh4_0, this%cutoff_nh4_1, &
    f_nh4, d_nh4)

  f_nh4_inhibit = 1.0d0
  d_nh4_inhibit_dnh4 = 0.0d0

  if (this%species_id_no3 > 0) then
    if (this%is_NO3_aqueous) then   
      c_no3     = rt_auxvar%pri_molal(this%species_id_no3)
      ac_no3    = rt_auxvar%pri_act_coef(this%species_id_no3)
      ires_no3 = this%species_id_no3
    else
      c_no3    = rt_auxvar%immobile(this%species_id_no3)
      ac_no3   = 1.0d0
      ires_no3 = this%species_id_no3 + reaction%offset_immobile
    endif

    call CalNLimitFunc(c_no3, ac_no3, this%residual_no3, &
      this%half_saturation_no3, this%cutoff_no3_0, this%cutoff_no3_1, &
      f_no3, d_no3)

    if (this%inhibition_nh4_no3 > this%residual_nh4) then 
      temp_real = this%inhibition_nh4_no3 + c_nh4 * ac_nh4
      f_nh4_inhibit = this%inhibition_nh4_no3/temp_real
      if (compute_derivative) then
        d_nh4_inhibit_dnh4 = -1.0d0 * this%inhibition_nh4_no3 * ac_nh4 &
                           / temp_real / temp_real
      endif
    endif
  endif 

  ires_co2 = this%species_id_co2
  ires_n2o = this%species_id_n2o

  ires_un = -999
  ires_hrimm = -999
  ires_nmin = -999
  ires_nimm = -999
  ires_b = -999
  ires_f = -999
  ires_ngasmin = -999

  if (this%species_id_hrimm > 0) then
    ires_hrimm = this%species_id_hrimm + reaction%offset_immobile
  endif

  if (this%species_id_nmin > 0) then
    ires_nmin = this%species_id_nmin + reaction%offset_immobile
  endif

  if (this%species_id_nimm > 0) then
    ires_nimm = this%species_id_nimm + reaction%offset_immobile
  endif

  if (this%species_id_ngasmin > 0) then
    ires_ngasmin = this%species_id_ngasmin + reaction%offset_immobile
  endif

  ! temperature response function
  tc = global_auxvar%temp

  f_t = 1.0d0

  saturation = global_auxvar%sat(1)
  theta = saturation * porosity 
  ! if positive, saturated soil's psi is nearly zero
  psi = min(global_auxvar%pres(1) - option%reference_pressure, -1.d-20)   

  ! moisture response function 
  f_w = 1.0d0

  if (f_t < 1.0d-20 .or. f_w < 1.0d-20) then
     return
  endif

  if (this%species_id_n2o > 0) then
    net_n_mineralization_rate = 0.0d0
    dnet_n_mineralization_rate_dnh4 = 0.0d0
    dnet_n_mineralization_rate_dno3 = 0.0d0
    do irxn = 1, this%nrxn
      dnet_n_mineralization_rate_duc(irxn) = 0.0d0
    enddo
  endif

  drate_nh4_duc = 0.0d0
  rdu_dun = 0.0d0
  rdu_duc = 0.0d0
  rdn_dun = 0.0d0
  rdn_duc = 0.0d0
  rdb_dun = 0.0d0
  rdb_duc = 0.0d0
  rdf_dun = 0.0d0
  rdf_duc = 0.0d0
  rdc_dun = 0.0d0
  rdc_duc = 0.0d0
 
  rate_no3 = 0.0d0
  drate_no3_duc = 0.0d0
  drate_no3_dnh4 = 0.0d0
  drate_no3_dno3 = 0.0d0
  rno3du_dun = 0.0d0
  rno3du_duc = 0.0d0
  rno3dn_dun = 0.0d0
  rno3dn_duc = 0.0d0
  rno3db_dun = 0.0d0
  rno3db_duc = 0.0d0
  rno3df_dun = 0.0d0
  rno3df_duc = 0.0d0
  rno3dc_dun = 0.0d0
  rno3dc_duc = 0.0d0

  resp_frac = 0.0d0

  do irxn = 1, this%nrxn
  
    ! upstream pool
    ispec_uc = this%upstream_c_id(irxn)

    if (this%upstream_is_aqueous(irxn)) then
      c_uc    = rt_auxvar%pri_molal(ispec_uc)
      ac_uc   = rt_auxvar%pri_act_coef(ispec_uc)
      ires_uc = ispec_uc
    else
      c_uc    = rt_auxvar%immobile(ispec_uc)
      ac_uc   = 1.0d0
      ires_uc = reaction%offset_immobile + ispec_uc
    endif

    ! for litter decomposition reactions, stoich needs to be calculated on the fly
    if (this%is_litter_decomp(irxn)) then

      ispec_un = this%upstream_n_id(irxn)
      if (this%upstream_is_aqueous(irxn)) then
        c_un    = rt_auxvar%pri_molal(ispec_un)
        ac_un   = rt_auxvar%pri_act_coef(ispec_un)
        ires_un = ispec_un
      else
        c_un    = rt_auxvar%immobile(ispec_un)
        ac_un   = 1.0d0
        ires_un = ispec_un + reaction%offset_immobile
      endif
      this%upstream_nc(irxn) = c_un / c_uc

      if (this%litter_decomp_type == LITTER_DECOMP_CLMCN) then

        ! calculate respiration factor (CO2 stoichiometry)
        stoich_c = 1.0d0

        do j = 1, this%n_downstream_pools(irxn)
          stoich_c = stoich_c - this%downstream_stoich(irxn, j)
        enddo

        if (stoich_c < 0.0d0) then
          option%io_buffer = 'CLMDec litter decomposition reaction has' // &
                             'negative respiration fraction!'
          call printErrMsg(option)
        endif

        this%mineral_c_stoich(irxn) = stoich_c

        ! calculate N stoichiometry
        stoich_n = this%upstream_nc(irxn)

        do j = 1, this%n_downstream_pools(irxn)
          stoich_n = stoich_n - this%downstream_stoich(irxn, j) * &
                              this%downstream_nc(irxn, j)
        enddo

        this%mineral_n_stoich(irxn) = stoich_n

      elseif (this%litter_decomp_type == LITTER_DECOMP_CLMMICROBE) then

        ! Sinsabaugh et al. 2013 Ecology Letters, 16, 930-939
        resp_frac = CN_ratio_microbe * this%upstream_nc(irxn) !c_un/c_uc

        if (resp_frac > CUE_max) then
          resp_frac = CUE_max
        endif

        ! c pools
        this%mineral_c_stoich(irxn) = resp_frac
       
        if (this%n_downstream_pools(irxn) .ne. 2) then
          option%io_buffer = 'CLM_Microbe litter decomposition reaction ' // &
                              'more than 2 (bacteria and fungi pools)!'
          call printErrMsg(option)
        endif

        do i = 1, this%n_downstream_pools(irxn)
          if (this%downstream_id(irxn, i) == this%species_id_bacteria) then
            this%downstream_stoich(irxn, i) = this%fraction_bacteria * &
              (1.0d0  - resp_frac)
          else
            this%downstream_stoich(irxn, i) = (1.0d0 - this%fraction_bacteria) &
                                            * (1.0d0  - resp_frac)
          endif
        enddo

        stoich_n = this%upstream_nc(irxn)

        do j = 1, this%n_downstream_pools(irxn)
          stoich_n = stoich_n - this%downstream_stoich(irxn, j) * &
                                this%downstream_nc(irxn, j)
        enddo

        this%mineral_n_stoich(irxn) = stoich_n

      endif

    endif

    if (this%upstream_is_aqueous(irxn)) then
      ! scaled_rate_const units: (kg water / s) = (1/s) * (kg water)
      scaled_rate_const = this%rate_constant(irxn) * volume * f_t * f_w &
                        * 1000.0d0 * theta
      ! will need to replace 1000.0d0 with water density
    else
      ! scaled_rate_const units: (m^3 bulk / s) = (1/s) * (m^3 bulk)
      scaled_rate_const = this%rate_constant(irxn)*volume*f_t*f_w
    endif

    ! residual units: (mol/sec) = (kg water/s) * (mol/kg water) or
    ! residual units: (mol/sec) = (m^3 bulk/s) * (mol/m^3 bulk)
    rate_nh4 = scaled_rate_const * (c_uc - this%residual_cpool) * ac_uc

    if (compute_derivative) then
      drate_nh4_duc = scaled_rate_const * ac_uc
    endif

    ! NH4 limiting
    if (this%mineral_n_stoich(irxn) < 0.0d0) then
      if (compute_derivative) then
        drate_nh4_dnh4 = rate_nh4 * d_nh4 
      endif
      rate_nh4       = rate_nh4 * f_nh4
      if (compute_derivative) then
        drate_nh4_duc  = drate_nh4_duc * f_nh4
      endif
    else
      drate_nh4_dnh4 = 0.d0
    endif 

    ! CO2
    Residual(ires_co2) = Residual(ires_co2) - &
                         this%mineral_c_stoich(irxn) * rate_nh4
    if (this%species_id_hrimm > 0) then
      Residual(ires_hrimm) = Residual(ires_hrimm) - &
                             this%mineral_c_stoich(irxn) * rate_nh4
    endif
    
    ! NH4
    Residual(ires_nh4) = Residual(ires_nh4) - &
                         this%mineral_n_stoich(irxn) * rate_nh4
    
    if (this%species_id_nimm > 0 .and. this%mineral_n_stoich(irxn) < 0.0d0) then
      Residual(ires_nimm) = Residual(ires_nimm) + &
                            this%mineral_n_stoich(irxn) * rate_nh4
    endif

    if (this%species_id_nmin > 0 .and. this%mineral_n_stoich(irxn) > 0.0d0) then
       Residual(ires_nmin) = Residual(ires_nmin) - &
                             this%mineral_n_stoich(irxn) * rate_nh4
    endif

    ! upstream c
    Residual(ires_uc) = Residual(ires_uc) - (-1.d0) * rate_nh4

    ! upstream n
    if (this%is_litter_decomp(irxn)) then
      Residual(ires_un) = Residual(ires_un) - &
                          (-1.d0) * this%upstream_nc(irxn) * rate_nh4
    endif
    
    ! downstream pools
    do j = 1, this%n_downstream_pools(irxn)
      ispec_d = this%downstream_id(irxn, j)
      if (this%downstream_is_aqueous(irxn, j)) then
        ires_d = ispec_d
      else
        ires_d = reaction%offset_immobile + ispec_d
      endif
      if (ispec_d > 0) then
        Residual(ires_d) = Residual(ires_d) - &
                           this%downstream_stoich(irxn, j) * rate_nh4
      endif
    enddo

    ! separate sink and source term for using NH4+ as nutrient
    if (this%mineral_n_stoich(irxn) >= 0.0d0 ) then
      RateSupply_nh4(ires_co2) = RateSupply_nh4(ires_co2) - &
                         this%mineral_c_stoich(irxn) * rate_nh4

      if (this%species_id_hrimm > 0) then
        RateSupply_nh4(ires_hrimm) = RateSupply_nh4(ires_hrimm) - &
                             this%mineral_c_stoich(irxn) * rate_nh4
      endif

      RateSupply_nh4(ires_nh4) = RateSupply_nh4(ires_nh4) - & 
                           this%mineral_n_stoich(irxn) * rate_nh4
    
      if (this%species_id_nmin > 0) then
         RateSupply_nh4(ires_nmin) = RateSupply_nh4(ires_nmin) - &
                               this%mineral_n_stoich(irxn) * rate_nh4
      endif

      ! upstream c
      RateSupply_nh4(ires_uc) = RateSupply_nh4(ires_uc) - (-1.d0) * rate_nh4

      ! upstream n
      if (this%is_litter_decomp(irxn)) then
        RateSupply_nh4(ires_un) = RateSupply_nh4(ires_un) - &
                            (-1.d0) * this%upstream_nc(irxn) * rate_nh4
      endif
    
      ! downstream pools
      do j = 1, this%n_downstream_pools(irxn)
        ispec_d = this%downstream_id(irxn, j)
        if (this%downstream_is_aqueous(irxn, j)) then
          ires_d = ispec_d
        else
          ires_d = reaction%offset_immobile + ispec_d
        endif
        if (ispec_d > 0) then
          RateSupply_nh4(ires_d) = RateSupply_nh4(ires_d) - &
                             this%downstream_stoich(irxn, j) * rate_nh4
        endif
      enddo

    else
      RateDemand_nh4(ires_co2) = RateDemand_nh4(ires_co2) - &
                         this%mineral_c_stoich(irxn) * rate_nh4

      if (this%species_id_hrimm > 0) then
        RateDemand_nh4(ires_hrimm) = RateDemand_nh4(ires_hrimm) - &
                             this%mineral_c_stoich(irxn) * rate_nh4
      endif

      RateDemand_nh4(ires_nh4) = RateDemand_nh4(ires_nh4) - & 
                           this%mineral_n_stoich(irxn) * rate_nh4
      
      if (this%species_id_nimm > 0) then
        RateDemand_nh4(ires_nimm) = RateDemand_nh4(ires_nimm) + &
                              this%mineral_n_stoich(irxn) * rate_nh4
      endif

      ! upstream c
      RateDemand_nh4(ires_uc) = RateDemand_nh4(ires_uc) - (-1.d0) * rate_nh4

      ! upstream n
      if (this%is_litter_decomp(irxn)) then
        RateDemand_nh4(ires_un) = RateDemand_nh4(ires_un) - &
                            (-1.d0) * this%upstream_nc(irxn) * rate_nh4
      endif
    
      ! downstream pools
      do j = 1, this%n_downstream_pools(irxn)
        ispec_d = this%downstream_id(irxn, j)
        if (this%downstream_is_aqueous(irxn, j)) then
          ires_d = ispec_d
        else
          ires_d = reaction%offset_immobile + ispec_d
        endif
        if (ispec_d > 0) then
          RateDemand_nh4(ires_d) = RateDemand_nh4(ires_d) - &
                             this%downstream_stoich(irxn, j) * rate_nh4
        endif
      enddo
    endif    

    if (this%species_id_n2o > 0) then
      net_n_mineralization_rate = net_n_mineralization_rate + &
        this%mineral_n_stoich(irxn) * rate_nh4

      if (compute_derivative) then
        dnet_n_mineralization_rate_dnh4 = dnet_n_mineralization_rate_dnh4 + &
          this%mineral_n_stoich(irxn) * drate_nh4_dnh4
        dnet_n_mineralization_rate_duc(irxn) = &
          dnet_n_mineralization_rate_duc(irxn) + &
          this%mineral_n_stoich(irxn) * drate_nh4_duc
      endif
    endif

    ! start residual calculation for N immobilization reaction with NO3 uptake
    ! if nitrate is available, N immobilization decomposition reactions occurs
    ! with rate depending on NH4, with reduced rate if NH4 is abundent  
    if (this%species_id_no3 > 0 .and. this%mineral_n_stoich(irxn) < 0.d0) then

      rate_no3       = scaled_rate_const * (c_uc - this%residual_cpool) &
                     * ac_uc * f_no3 * f_nh4_inhibit

      if (compute_derivative) then
        drate_no3_duc  = scaled_rate_const * ac_uc * f_no3 * f_nh4_inhibit
        drate_no3_dno3 = scaled_rate_const * (c_uc - this%residual_cpool) &
                       * ac_uc * d_no3 * f_nh4_inhibit
        drate_no3_dnh4 = scaled_rate_const * (c_uc - this%residual_cpool) &
                       * ac_uc * f_no3 * d_nh4_inhibit_dnh4
      endif

      ! carbon
      Residual(ires_co2) = Residual(ires_co2) - &
        this%mineral_c_stoich(irxn) * rate_no3
      if (this%species_id_hrimm > 0) then
        Residual(ires_hrimm) = Residual(ires_hrimm) - &
        this%mineral_c_stoich(irxn) * rate_no3
      endif
    
      ! NO3
      Residual(ires_no3) = Residual(ires_no3) - &
        this%mineral_n_stoich(irxn) * rate_no3

      if (this%species_id_nimm > 0) then
        Residual(ires_nimm) = Residual(ires_nimm) + &
          this%mineral_n_stoich(irxn) * rate_no3
      endif

      ! upstream c
      Residual(ires_uc) = Residual(ires_uc) - (-1.d0) * rate_no3
    
      ! upstream n
      if (this%is_litter_decomp(irxn)) then
        Residual(ires_un) = Residual(ires_un) + this%upstream_nc(irxn) *rate_no3
      endif
    
      ! downstream pools
      do j = 1, this%n_downstream_pools(irxn)
        ispec_d = this%downstream_id(irxn, j)
        if (this%downstream_is_aqueous(irxn, j)) then
          ires_d = ispec_d
        else
          ires_d = reaction%offset_immobile + ispec_d
        endif
        if (ispec_d > 0) then
          Residual(ires_d) = Residual(ires_d) - &
            this%downstream_stoich(irxn, j) * rate_no3
         endif
      enddo

      ! separate sink and source for using NO3- as nutrient
      RateDemand_no3(ires_co2) = RateDemand_no3(ires_co2) - &
        this%mineral_c_stoich(irxn) * rate_no3

      if (this%species_id_hrimm > 0) then
        RateDemand_no3(ires_hrimm) = RateDemand_no3(ires_hrimm) - &
        this%mineral_c_stoich(irxn) * rate_no3
      endif
    
      RateDemand_no3(ires_no3) = RateDemand_no3(ires_no3) - & 
                           this%mineral_n_stoich(irxn) * rate_no3

      if (this%species_id_nimm > 0) then
        RateDemand_no3(ires_nimm) = RateDemand_no3(ires_nimm) + &
          this%mineral_n_stoich(irxn) * rate_no3
      endif

      ! upstream c
      RateDemand_no3(ires_uc) = RateDemand_no3(ires_uc) - (-1.d0) * rate_no3
    
      ! upstream n
      if (this%is_litter_decomp(irxn)) then
        RateDemand_no3(ires_un) = RateDemand_no3(ires_un) + &
          this%upstream_nc(irxn) * rate_no3
      endif
    
      ! downstream pools
      do j = 1, this%n_downstream_pools(irxn)
        ispec_d = this%downstream_id(irxn, j)
        if (this%downstream_is_aqueous(irxn, j)) then
          ires_d = ispec_d
        else
          ires_d = reaction%offset_immobile + ispec_d
        endif
        if (ispec_d > 0) then
          RateDemand_no3(ires_d) = RateDemand_no3(ires_d) - &
            this%downstream_stoich(irxn, j) * rate_no3
         endif
      enddo

      if (this%species_id_n2o > 0) then
        net_n_mineralization_rate = net_n_mineralization_rate + &
          this%mineral_n_stoich(irxn)*rate_no3
        if (compute_derivative) then
          dnet_n_mineralization_rate_dnh4 = dnet_n_mineralization_rate_dnh4 + &
            this%mineral_n_stoich(irxn) * drate_no3_dnh4
          dnet_n_mineralization_rate_dno3 = dnet_n_mineralization_rate_dno3 + &
            this%mineral_n_stoich(irxn) * drate_no3_dno3
          dnet_n_mineralization_rate_duc(irxn) = &
            dnet_n_mineralization_rate_duc(irxn) + &
            this%mineral_n_stoich(irxn) * drate_no3_duc
        endif
      endif 
    endif

    if (compute_derivative) then

      if (this%is_litter_decomp(irxn)) then
        if (this%litter_decomp_type == LITTER_DECOMP_CLMCN) then
          ! LitC + u LitN -> di SOMi + (1 - di) CO2 + n N
          ! Rdu/duc = R (-1) LitN/LitC^2 = - u R / LitC 
          Rdu_duc = -1.0d0 * this%upstream_nc(irxn) * drate_nh4_duc
         
          ! n = u - (1 - di) ni
          ! dn/dLitC = du/dLitC
          Rdn_duc = Rdu_duc

          ! Rdu/dun = R /LitC 
          Rdu_dun = drate_nh4_duc 

          ! Rdn/dun = Rdu/dLitN = Rdu/dun
          Rdn_dun = Rdu_dun

        elseif (this%litter_decomp_type == LITTER_DECOMP_CLMMICROBE) then
          ! Lit1C + u Lit1N -> b Bacteria + f Fungi + c CO2 + n N
          ! c = min(CUEmax, Lit1N/Lit1C*CN_ratio_microbe)
          ! g = CNbacteria^0.6/(CNbacterial^0.6 + CNfungi^0.6)
          ! n = u - b nb - f nf
          ! b = g (1 - c)
          ! f = (1 - g) (1 - c)

          Rdu_duc = -1.0d0 * this%upstream_nc(irxn) * drate_nh4_duc
   
          if (resp_frac < CUE_max) then
            ! Rdc/dLit1C = -RLit1N/Lit1C^2*CN_ratio_microbe           
            Rdc_duc = -1.0d0 * this%upstream_nc(irxn) * drate_nh4_duc  &
                    * CN_ratio_microbe
          else
            Rdc_duc = 0.0d0 
          endif

          ! Rdb/dLitC = -g Rdc/dLitC
          Rdb_duc = -1.0d0 * this%fraction_bacteria * Rdc_duc
  
          ! Rdf/dLitC = -(1 - g) Rdc/dLitC
          Rdf_duc = -1.0d0 * (1.0d0 - this%fraction_bacteria) * Rdc_duc

          ! Rdn/dLitC = Rdu/dLitC - nb Rdb/dLitC - nf Rdf/dLitC
          Rdn_duc = Rdu_duc - this%nc_bacteria * Rdb_duc &
                  - this%nc_fungi * Rdf_duc

          ! Rdu/dun = R/LitC = dR/duc
          Rdu_dun = drate_nh4_duc 

          if (resp_frac < CUE_max) then
            ! Rdc/dLitN = R/LitC*CN_ratio_microbe
            Rdc_dun = drate_nh4_duc * CN_ratio_microbe
          else
            Rdc_dun = 0.0d0 
          endif

          ! Rdb/dLitN = -g Rdc/dLitN
          Rdb_dun = -1.0d0 * this%fraction_bacteria * Rdc_dun

          ! Rdf/dLitN = -(1 - g) Rdc/dLitN
          Rdf_dun = -1.0d0 * (1.0d0 - this%fraction_bacteria) * Rdc_dun

          ! Rdn/dLitN = Rdu/dLitN - nb Rdb/dLitN - nf Rdf/dLitN
          Rdn_dun = Rdu_dun - this%nc_bacteria * Rdb_dun &
                  - this%nc_fungi * Rdf_dun

          ires_b = reaction%offset_immobile + this%species_id_bacteria 
          ires_f = reaction%offset_immobile + this%species_id_fungi 
            
        endif

      endif   

      ! with respect to upstream C
      ! CO2
      Jacobian(ires_co2,ires_uc) = Jacobian(ires_co2,ires_uc) - &
        this%mineral_c_stoich(irxn) * drate_nh4_duc

      if (this%is_litter_decomp(irxn) .and. &
        this%litter_decomp_type == LITTER_DECOMP_CLMMICROBE) then
        ! dRco2/dLitC = dcR/dLitC = cdR/dLitC + R dc/dLitC
        Jacobian(ires_co2,ires_uc) = Jacobian(ires_co2,ires_uc) - Rdc_duc
      endif

      if (this%species_id_hrimm > 0) then
        Jacobian(ires_hrimm,ires_uc) = Jacobian(ires_hrimm,ires_uc) - &
          this%mineral_c_stoich(irxn) * drate_nh4_duc

        if (this%is_litter_decomp(irxn) .and. &
           this%litter_decomp_type == LITTER_DECOMP_CLMMICROBE) then
           Jacobian(ires_hrimm,ires_uc) = Jacobian(ires_hrimm,ires_uc) - Rdc_duc
        endif
      endif

      ! N
      ! dR_N/dC_u = d nR/dC_u = dn/dC_u R + n dR/dC_u
      ! first, n dR/dC_u
      Jacobian(ires_nh4,ires_uc) = Jacobian(ires_nh4,ires_uc) - &
        this%mineral_n_stoich(irxn) * drate_nh4_duc

      if (this%species_id_nimm > 0 .and. &
          this%mineral_n_stoich(irxn) < 0.0d0) then
        Jacobian(ires_nimm,ires_uc) = Jacobian(ires_nimm,ires_uc) + &
          this%mineral_n_stoich(irxn) * drate_nh4_duc
      endif

      if (this%species_id_nmin > 0 .and. &
          this%mineral_n_stoich(irxn) > 0.0d0) then
        Jacobian(ires_nmin,ires_uc) = Jacobian(ires_nmin,ires_uc) - &
          this%mineral_n_stoich(irxn) * drate_nh4_duc
      endif

      if (this%is_litter_decomp(irxn)) then
        ! litter pool is immobile
        ! second, Rdn/dC_u
        Jacobian(ires_nh4,ires_uc) = Jacobian(ires_nh4,ires_uc) - Rdn_duc

        if (this%species_id_nimm > 0 .and. &
            this%mineral_n_stoich(irxn) < 0.0d0) then
          Jacobian(ires_nimm,ires_uc) = Jacobian(ires_nimm,ires_uc) + Rdn_duc
        endif

        if (this%species_id_nmin > 0 .and. &
            this%mineral_n_stoich(irxn) > 0.0d0) then
          Jacobian(ires_nmin,ires_uc) = Jacobian(ires_nmin,ires_uc) - Rdn_duc
        endif

      endif

      ! upstream C pool
      Jacobian(ires_uc,ires_uc) = Jacobian(ires_uc,ires_uc) + drate_nh4_duc

      ! upstream N pool
      if (this%is_litter_decomp(irxn)) then
        ! litter pools are immobile
        ! R_Nu = Nu/Cu * R_Cu
        ! dR_Nu/dCu = Nu/Cu dR_Cu/dCu + R du/dCu
        ! = 0 only when residual_cpool = 0
        Jacobian(ires_un,ires_uc) = Jacobian(ires_un,ires_uc) + &
          this%upstream_nc(irxn) * drate_nh4_duc + Rdu_duc
        
      endif

      ! downstream pools
      do j = 1, this%n_downstream_pools(irxn)
        ispec_d = this%downstream_id(irxn, j)
        if (ispec_d < 0) then
          option%io_buffer = 'Downstream pool species not specified!'
          call printErrMsg(option)
        endif

        if (this%downstream_is_aqueous(irxn, j)) then
          ires_d = ispec_d
        else
          ires_d = reaction%offset_immobile + ispec_d
        endif
         
        Jacobian(ires_d,ires_uc) = Jacobian(ires_d,ires_uc) - &
          this%downstream_stoich(irxn, j) * drate_nh4_duc

        ! additional term if downstream stoich is variable
        if (this%is_litter_decomp(irxn) .and. &
            this%litter_decomp_type == LITTER_DECOMP_CLMMICROBE) then
          if (ispec_d == this%species_id_bacteria) then
            ! dRbacteria/dLit1C = b dR/dLit1C + R db/dLit1C
            Jacobian(ires_d,ires_uc) = Jacobian(ires_d,ires_uc) - Rdb_duc
          elseif (ispec_d == this%species_id_fungi) then
            Jacobian(ires_d,ires_uc) = Jacobian(ires_d,ires_uc) - Rdf_duc
          else
            option%io_buffer = 'Downstream pool for CLM-Microbe should be' // &
                               'either bacteria or fungi!'
            call printErrMsg(option)
          endif
        endif
      enddo

      ! with respect to upstream n (due to variable CN ratio)
      if (this%is_litter_decomp(irxn)) then
        ! upstream n, dR_Nu/dNu = d uR/dNu = R/Cu = dR/dCu
        Jacobian(ires_un,ires_un) = Jacobian(ires_un,ires_un) + Rdu_dun

        ! mineral N, dR_N/dNu = d nR/du = Rd (u - (1-f)d)/dNu = dR/dCu
        Jacobian(ires_nh4,ires_un) = Jacobian(ires_nh4,ires_un) - Rdn_dun

        if (this%species_id_nimm > 0 .and. &
            this%mineral_n_stoich(irxn) < 0.0d0) then
          Jacobian(ires_nimm,ires_un) = Jacobian(ires_nimm,ires_un) + Rdn_dun
        endif

        if (this%species_id_nmin > 0 .and. &
            this%mineral_n_stoich(irxn) > 0.0d0) then
          Jacobian(ires_nmin,ires_un) = Jacobian(ires_nmin,ires_un) - Rdn_dun
        endif

        if (this%litter_decomp_type == LITTER_DECOMP_CLMMICROBE) then
          ! CO2 dR_co2/dNu = d cR/dNu = Rdc/Nu
          Jacobian(ires_co2,ires_un) = Jacobian(ires_co2,ires_un) - Rdc_dun
       
          if (this%species_id_hrimm > 0) then
            Jacobian(ires_hrimm,ires_un) = Jacobian(ires_hrimm,ires_un) -Rdc_dun
          endif

          ! bacteria dRbacteria/dNu = Rdb/dNu
          Jacobian(ires_b,ires_un) = Jacobian(ires_b,ires_un) - Rdb_dun

          ! fungi dRfungi/dNu = Rdf/dNu
          Jacobian(ires_f,ires_un) = Jacobian(ires_f,ires_un) - Rdf_dun
        endif 
      endif

      ! with respect to nh4
      if (this%mineral_n_stoich(irxn) < 0.0d0) then
        ! CO2
        Jacobian(ires_co2,ires_nh4) = Jacobian(ires_co2,ires_nh4) - &
          this%mineral_c_stoich(irxn) * drate_nh4_dnh4

        if (this%species_id_hrimm > 0) then
          Jacobian(ires_hrimm,ires_nh4) = Jacobian(ires_hrimm,ires_nh4) - &
            this%mineral_c_stoich(irxn) * drate_nh4_dnh4
        endif

        ! N
        Jacobian(ires_nh4,ires_nh4) = Jacobian(ires_nh4,ires_nh4) - &
          this%mineral_n_stoich(irxn) * drate_nh4_dnh4
  
        if (this%species_id_nimm > 0) then
          Jacobian(ires_nimm,ires_nh4) = Jacobian(ires_nimm,ires_nh4) + &
            this%mineral_n_stoich(irxn) * drate_nh4_dnh4
        endif

        ! upstream C
        Jacobian(ires_uc,ires_nh4) = Jacobian(ires_uc,ires_nh4) - &
          (-1.d0) * drate_nh4_dnh4
 
        ! upstream N pool
        if (this%is_litter_decomp(irxn)) then
          Jacobian(ires_un,ires_nh4) = Jacobian(ires_un,ires_nh4) - &
            (-1.d0) * this%upstream_nc(irxn) * drate_nh4_dnh4 
        endif

        ! downstream pools
        do j = 1, this%n_downstream_pools(irxn)
          ispec_d = this%downstream_id(irxn, j)
          if (ispec_d < 0) then
            option%io_buffer = 'Downstream pool species not specified!'
            call printErrMsg(option)
          endif
          if (this%downstream_is_aqueous(irxn, j)) then
            ires_d = ispec_d
          else
            ires_d = reaction%offset_immobile + ispec_d
          endif
          Jacobian(ires_d,ires_nh4) = Jacobian(ires_d,ires_nh4) - &
            this%downstream_stoich(irxn, j) * drate_nh4_dnh4
        enddo

      endif

      ! separate sink and source for using NH4+ as nutrient
      if (this%mineral_n_stoich(irxn) >= 0.0d0) then
        ! with respect to upstream C
        JacobianSupply_nh4(ires_co2,ires_uc) = &
          JacobianSupply_nh4(ires_co2,ires_uc) - &
          this%mineral_c_stoich(irxn) * drate_nh4_duc

        if (this%is_litter_decomp(irxn) .and. &
          this%litter_decomp_type == LITTER_DECOMP_CLMMICROBE) then
          ! dRco2/dLitC = dcR/dLitC = cdR/dLitC + R dc/dLitC
          JacobianSupply_nh4(ires_co2,ires_uc) = &
            JacobianSupply_nh4(ires_co2,ires_uc) - Rdc_duc
        endif

        if (this%species_id_hrimm > 0) then
          JacobianSupply_nh4(ires_hrimm,ires_uc) = &
            JacobianSupply_nh4(ires_hrimm,ires_uc) - &
            this%mineral_c_stoich(irxn) * drate_nh4_duc

          if (this%is_litter_decomp(irxn) .and. &
            this%litter_decomp_type == LITTER_DECOMP_CLMMICROBE) then
            JacobianSupply_nh4(ires_hrimm,ires_uc) = &
              JacobianSupply_nh4(ires_hrimm,ires_uc) - Rdc_duc
          endif
        endif

        ! N
        ! dR_N/dC_u = d nR/dC_u = dn/dC_u R + n dR/dC_u
        ! first, n dR/dC_u
        JacobianSupply_nh4(ires_nh4,ires_uc) = &
          JacobianSupply_nh4(ires_nh4,ires_uc) - &
          this%mineral_n_stoich(irxn) * drate_nh4_duc

        if (this%species_id_nmin > 0) then
          JacobianSupply_nh4(ires_nmin,ires_uc) = &
            JacobianSupply_nh4(ires_nmin,ires_uc) - &
            this%mineral_n_stoich(irxn) * drate_nh4_duc
        endif

        if (this%is_litter_decomp(irxn)) then
          ! litter pool is immobile
          ! second, Rdn/dC_u
          JacobianSupply_nh4(ires_nh4,ires_uc) = &
            JacobianSupply_nh4(ires_nh4,ires_uc) - Rdn_duc

          if (this%species_id_nmin > 0) then
            JacobianSupply_nh4(ires_nmin,ires_uc) = &
              JacobianSupply_nh4(ires_nmin,ires_uc) - Rdn_duc
          endif

        endif

        ! upstream C pool
        JacobianSupply_nh4(ires_uc,ires_uc) = &
          JacobianSupply_nh4(ires_uc,ires_uc) + drate_nh4_duc

        ! upstream N pool
        if (this%is_litter_decomp(irxn)) then
          ! litter pools are immobile
          ! R_Nu = Nu/Cu * R_Cu
          ! dR_Nu/dCu = Nu/Cu dR_Cu/dCu + R du/dCu
          ! = 0 only when residual_cpool = 0
          JacobianSupply_nh4(ires_un,ires_uc) = &
            JacobianSupply_nh4(ires_un,ires_uc) + &
            this%upstream_nc(irxn) * drate_nh4_duc + Rdu_duc
        
        endif

        ! downstream pools
        do j = 1, this%n_downstream_pools(irxn)
          ispec_d = this%downstream_id(irxn, j)
          if (ispec_d < 0) then
            option%io_buffer = 'Downstream pool species not specified!'
            call printErrMsg(option)
          endif

          if (this%downstream_is_aqueous(irxn, j)) then
            ires_d = ispec_d
          else
            ires_d = reaction%offset_immobile + ispec_d
          endif
         
          JacobianSupply_nh4(ires_d,ires_uc) = &
            JacobianSupply_nh4(ires_d,ires_uc) - &
            this%downstream_stoich(irxn, j) * drate_nh4_duc

          ! additional term if downstream stoich is variable
          if (this%is_litter_decomp(irxn) .and. &
            this%litter_decomp_type == LITTER_DECOMP_CLMMICROBE) then
            if (ispec_d == this%species_id_bacteria) then
              ! dRbacteria/dLit1C = b dR/dLit1C + R db/dLit1C
              JacobianSupply_nh4(ires_d,ires_uc) = &
                JacobianSupply_nh4(ires_d,ires_uc) - Rdb_duc
            elseif (ispec_d == this%species_id_fungi) then
              JacobianSupply_nh4(ires_d,ires_uc) = &
                JacobianSupply_nh4(ires_d,ires_uc) - Rdf_duc
            else
              option%io_buffer ='Downstream pool for CLM-Microbe should be' // &
                               'either bacteria or fungi!'
              call printErrMsg(option)
            endif
          endif
        enddo

        ! with respect to upstream n (due to variable CN ratio)
        if (this%is_litter_decomp(irxn)) then
          ! upstream n, dR_Nu/dNu = d uR/dNu = R/Cu = dR/dCu
          JacobianSupply_nh4(ires_un,ires_un) = &
            JacobianSupply_nh4(ires_un,ires_un) + Rdu_dun

          ! mineral N, dR_N/dNu = d nR/du = Rd (u - (1-f)d)/dNu = dR/dCu
          JacobianSupply_nh4(ires_nh4,ires_un) = &
            JacobianSupply_nh4(ires_nh4,ires_un) - Rdn_dun

          if (this%species_id_nmin > 0) then
            JacobianSupply_nh4(ires_nmin,ires_un) = &
              JacobianSupply_nh4(ires_nmin,ires_un) - Rdn_dun
          endif

          if (this%litter_decomp_type == LITTER_DECOMP_CLMMICROBE) then
            ! CO2 dR_co2/dNu = d cR/dNu = Rdc/Nu
            JacobianSupply_nh4(ires_co2,ires_un) = &
              JacobianSupply_nh4(ires_co2,ires_un) - Rdc_dun
       
            if (this%species_id_hrimm > 0) then
              JacobianSupply_nh4(ires_hrimm,ires_un) = &
                JacobianSupply_nh4(ires_hrimm,ires_un) -Rdc_dun
            endif

            ! bacteria dRbacteria/dNu = Rdb/dNu
            JacobianSupply_nh4(ires_b,ires_un) = &
              JacobianSupply_nh4(ires_b,ires_un) - Rdb_dun

            ! fungi dRfungi/dNu = Rdf/dNu
            JacobianSupply_nh4(ires_f,ires_un) = &
              JacobianSupply_nh4(ires_f,ires_un) - Rdf_dun
          endif 
        endif

      else
        ! with respect to upstream C
        JacobianDemand_nh4(ires_co2,ires_uc) = &
          JacobianDemand_nh4(ires_co2,ires_uc) - &
          this%mineral_c_stoich(irxn) * drate_nh4_duc

        if (this%is_litter_decomp(irxn) .and. &
          this%litter_decomp_type == LITTER_DECOMP_CLMMICROBE) then
          ! dRco2/dLitC = dcR/dLitC = cdR/dLitC + R dc/dLitC
          JacobianDemand_nh4(ires_co2,ires_uc) = &
            JacobianDemand_nh4(ires_co2,ires_uc) - Rdc_duc
        endif

        if (this%species_id_hrimm > 0) then
          JacobianDemand_nh4(ires_hrimm,ires_uc) = &
            JacobianDemand_nh4(ires_hrimm,ires_uc) - &
            this%mineral_c_stoich(irxn) * drate_nh4_duc

          if (this%is_litter_decomp(irxn) .and. &
            this%litter_decomp_type == LITTER_DECOMP_CLMMICROBE) then
            JacobianDemand_nh4(ires_hrimm,ires_uc) = &
              JacobianDemand_nh4(ires_hrimm,ires_uc) - Rdc_duc
          endif
        endif

        ! N
        ! dR_N/dC_u = d nR/dC_u = dn/dC_u R + n dR/dC_u
        ! first, n dR/dC_u
        JacobianDemand_nh4(ires_nh4,ires_uc) = &
          JacobianDemand_nh4(ires_nh4,ires_uc) - &
          this%mineral_n_stoich(irxn) * drate_nh4_duc

        if (this%species_id_nimm > 0) then
          JacobianDemand_nh4(ires_nimm,ires_uc) = &
            JacobianDemand_nh4(ires_nimm,ires_uc) + &
            this%mineral_n_stoich(irxn) * drate_nh4_duc
        endif

        if (this%is_litter_decomp(irxn)) then
          ! litter pool is immobile
          ! second, Rdn/dC_u
          JacobianDemand_nh4(ires_nh4,ires_uc) = &
            JacobianDemand_nh4(ires_nh4,ires_uc) - Rdn_duc

          if (this%species_id_nimm > 0) then
            JacobianDemand_nh4(ires_nimm,ires_uc) = &
              JacobianDemand_nh4(ires_nimm,ires_uc) + Rdn_duc
          endif

        endif

        ! upstream C pool
        JacobianDemand_nh4(ires_uc,ires_uc) = &
          JacobianDemand_nh4(ires_uc,ires_uc) + drate_nh4_duc

        ! upstream N pool
        if (this%is_litter_decomp(irxn)) then
          ! litter pools are immobile
          ! R_Nu = Nu/Cu * R_Cu
          ! dR_Nu/dCu = Nu/Cu dR_Cu/dCu + R du/dCu
          ! = 0 only when residual_cpool = 0
          JacobianDemand_nh4(ires_un,ires_uc) = &
            JacobianDemand_nh4(ires_un,ires_uc) + &
            this%upstream_nc(irxn) * drate_nh4_duc + Rdu_duc
        
        endif

        ! downstream pools
        do j = 1, this%n_downstream_pools(irxn)
          ispec_d = this%downstream_id(irxn, j)
          if (ispec_d < 0) then
            option%io_buffer = 'Downstream pool species not specified!'
            call printErrMsg(option)
          endif

          if (this%downstream_is_aqueous(irxn, j)) then
            ires_d = ispec_d
          else
            ires_d = reaction%offset_immobile + ispec_d
          endif
         
          JacobianDemand_nh4(ires_d,ires_uc) = &
            JacobianDemand_nh4(ires_d,ires_uc) - &
            this%downstream_stoich(irxn, j) * drate_nh4_duc

          ! additional term if downstream stoich is variable
          if (this%is_litter_decomp(irxn) .and. &
            this%litter_decomp_type == LITTER_DECOMP_CLMMICROBE) then
            if (ispec_d == this%species_id_bacteria) then
              ! dRbacteria/dLit1C = b dR/dLit1C + R db/dLit1C
              JacobianDemand_nh4(ires_d,ires_uc) = &
                JacobianDemand_nh4(ires_d,ires_uc) - Rdb_duc
            elseif (ispec_d == this%species_id_fungi) then
              JacobianDemand_nh4(ires_d,ires_uc) = &
                JacobianDemand_nh4(ires_d,ires_uc) - Rdf_duc
            else
              option%io_buffer ='Downstream pool for CLM-Microbe should be' // &
                               'either bacteria or fungi!'
              call printErrMsg(option)
            endif
          endif
        enddo

        ! with respect to upstream n (due to variable CN ratio)
        if (this%is_litter_decomp(irxn)) then
          ! upstream n, dR_Nu/dNu = d uR/dNu = R/Cu = dR/dCu
          JacobianDemand_nh4(ires_un,ires_un) = &
            JacobianDemand_nh4(ires_un,ires_un) + Rdu_dun

          ! mineral N, dR_N/dNu = d nR/du = Rd (u - (1-f)d)/dNu = dR/dCu
          JacobianDemand_nh4(ires_nh4,ires_un) = &
            JacobianDemand_nh4(ires_nh4,ires_un) - Rdn_dun

          if (this%species_id_nimm > 0) then
            JacobianDemand_nh4(ires_nimm,ires_un) = &
              JacobianDemand_nh4(ires_nimm,ires_un) + Rdn_dun
          endif

          if (this%litter_decomp_type == LITTER_DECOMP_CLMMICROBE) then
            ! CO2 dR_co2/dNu = d cR/dNu = Rdc/Nu
            JacobianDemand_nh4(ires_co2,ires_un) = &
              JacobianDemand_nh4(ires_co2,ires_un) - Rdc_dun
       
            if (this%species_id_hrimm > 0) then
              JacobianDemand_nh4(ires_hrimm,ires_un) = &
                JacobianDemand_nh4(ires_hrimm,ires_un) - Rdc_dun
            endif

            ! bacteria dRbacteria/dNu = Rdb/dNu
            JacobianDemand_nh4(ires_b,ires_un) = &
              JacobianDemand_nh4(ires_b,ires_un) - Rdb_dun

            ! fungi dRfungi/dNu = Rdf/dNu
            JacobianDemand_nh4(ires_f,ires_un) = &
              JacobianDemand_nh4(ires_f,ires_un) - Rdf_dun
          endif 
        endif

        ! with respect to nh4
        ! CO2
        JacobianDemand_nh4(ires_co2,ires_nh4) = &
          JacobianDemand_nh4(ires_co2,ires_nh4) - &
          this%mineral_c_stoich(irxn) * drate_nh4_dnh4

        if (this%species_id_hrimm > 0) then
          JacobianDemand_nh4(ires_hrimm,ires_nh4) = &
            JacobianDemand_nh4(ires_hrimm,ires_nh4) - &
            this%mineral_c_stoich(irxn) * drate_nh4_dnh4
        endif

        ! N
        JacobianDemand_nh4(ires_nh4,ires_nh4) = &
          JacobianDemand_nh4(ires_nh4,ires_nh4) - &
          this%mineral_n_stoich(irxn) * drate_nh4_dnh4
  
        if (this%species_id_nimm > 0) then
          JacobianDemand_nh4(ires_nimm,ires_nh4) = &
            JacobianDemand_nh4(ires_nimm,ires_nh4) + &
            this%mineral_n_stoich(irxn) * drate_nh4_dnh4
        endif

        ! upstream C
        JacobianDemand_nh4(ires_uc,ires_nh4) = &
          JacobianDemand_nh4(ires_uc,ires_nh4) - (-1.d0) * drate_nh4_dnh4
 
        ! upstream N pool
        if (this%is_litter_decomp(irxn)) then
          JacobianDemand_nh4(ires_un,ires_nh4) = &
            JacobianDemand_nh4(ires_un,ires_nh4) - &
            (-1.d0) * this%upstream_nc(irxn) * drate_nh4_dnh4 
        endif

        ! downstream pools
        do j = 1, this%n_downstream_pools(irxn)
          ispec_d = this%downstream_id(irxn, j)
          if (ispec_d < 0) then
            option%io_buffer = 'Downstream pool species not specified!'
            call printErrMsg(option)
          endif
          if (this%downstream_is_aqueous(irxn, j)) then
            ires_d = ispec_d
          else
            ires_d = reaction%offset_immobile + ispec_d
          endif
          JacobianDemand_nh4(ires_d,ires_nh4) = &
            JacobianDemand_nh4(ires_d,ires_nh4) - &
            this%downstream_stoich(irxn, j) * drate_nh4_dnh4
        enddo

      endif 


      if (this%species_id_no3 > 0 .and. this%mineral_n_stoich(irxn) < 0.d0) then

        if (this%is_litter_decomp(irxn)) then
          if (this%litter_decomp_type == LITTER_DECOMP_CLMCN) then
            ! Lit1C + u Lit1N -> di SOMi + (1 - di) CO2 + n N
            ! Rdu/duc = R (-1) Lit1N/Lit1C^2 
            Rno3du_duc = -1.0d0 * this%upstream_nc(irxn) * drate_no3_duc

            ! n = u - (1 - di) ni
            ! dn/dLit1C = du/dLit1C
            Rno3dn_duc = Rno3du_duc

            ! Rdu/dun = R /Lit1C 
            Rno3du_dun = drate_no3_duc

            ! Rdn/dun = du/dLit1C = dR/duc 
            Rno3dn_dun = Rno3du_dun

          elseif (this%litter_decomp_type == LITTER_DECOMP_CLMMICROBE) then
            ! Lit1C + u Lit1N -> b Bacteria + f Fungi + c CO2 + n N
            ! c = min(CUEmax, Lit1N/Lit1C*CN_ratio_microbe)
            ! g = CNbacteria^0.6/(CNbacterial^0.6 + CNfungi^0.6)
            ! n = u - b nb - f nf
            ! b = g (1 - c)
            ! f = (1 - g) (1 - c)

            Rno3du_duc = -1.0d0 * this%upstream_nc(irxn) * drate_no3_duc

            if (resp_frac < CUE_max) then
              ! Rdc/dLit1C = -RLit1N/Lit1C^2*CN_ratio_microbe           
              Rno3dc_duc = -1.0d0 * this%upstream_nc(irxn) * drate_no3_duc &
                         * CN_ratio_microbe
            else
              Rno3dc_duc = 0.0d0
            endif

            ! Rdb/dLit1C = -g Rdc/dLit1C  
            Rno3db_duc = -1.0d0 * this%fraction_bacteria * Rno3dc_duc

            ! Rdf/dLit1C = -(1 - g) Rdc/dLit1C  
            Rno3df_duc = -1.0d0 * (1.0d0 - this%fraction_bacteria) * Rno3dc_duc

            ! Rdn/dLit1C = Rdu/dLit1C - nb Rdb/dLit1C - nf Rdf/dLit1C  
            Rno3dn_duc = Rno3du_duc - this%nc_bacteria * Rno3db_duc &
                       - this%nc_fungi * Rno3df_duc

            ! Rdu/dun = R /Lit1N 
            Rno3du_dun = drate_no3_duc

            if (resp_frac < CUE_max) then
              ! Rdc/dLit1N = R/Lit1C*CN_ratio_microbe = dR/dLit1C*CN_ratio_microbe           
              Rno3dc_dun = drate_no3_duc * CN_ratio_microbe
            else
              Rno3dc_dun = 0.0d0
            endif

            ! Rdb/dLit1N = -g Rdc/dLit1N 
            Rno3db_dun = -1.0d0 * this%fraction_bacteria * Rno3dc_dun

            ! Rdf/dLit1N = -(1 - g) Rdc/dLit1N  
            Rno3df_dun = -1.0d0 * (1.0d0 - this%fraction_bacteria) * Rno3dc_dun

            ! Rdn/dLit1N = Rdu/dLit1N - nb Rdb/dLit1N - nf Rdf/dLit1N  
            Rno3dn_dun = Rno3du_dun - this%nc_bacteria * Rno3db_dun &
                       - this%nc_fungi * Rno3df_dun

          endif
        endif

        ! with respect to upstream
        ! CO2

        Jacobian(ires_co2,ires_uc) = Jacobian(ires_co2,ires_uc) - &
          this%mineral_c_stoich(irxn) * drate_no3_duc

        if (this%is_litter_decomp(irxn) .and. &
          this%litter_decomp_type == LITTER_DECOMP_CLMMICROBE) then
          ! dRco2/dLitC = dcR/dLitC = cdR/dLitC + R dc/dLitC
          Jacobian(ires_co2,ires_uc) = Jacobian(ires_co2,ires_uc) - Rno3dc_duc
        endif

        if (this%species_id_hrimm > 0) then
          Jacobian(ires_hrimm,ires_uc) = Jacobian(ires_hrimm,ires_uc) - &
            this%mineral_c_stoich(irxn) * drate_no3_duc

          if (this%is_litter_decomp(irxn) .and. &
            this%litter_decomp_type == LITTER_DECOMP_CLMMICROBE ) then 
            Jacobian(ires_hrimm,ires_uc) = Jacobian(ires_hrimm,ires_uc) &
                                         - Rno3dc_duc
          endif
        endif

        ! NO3
        ! dR_N/dC_u = d nR/dC_u = dn/dC_u R + n dR/dC_u
        ! first term: n dR/dC_u
        Jacobian(ires_no3,ires_uc) = Jacobian(ires_no3,ires_uc) - &
          this%mineral_n_stoich(irxn) * drate_no3_duc

        if (this%species_id_nimm > 0) then
          Jacobian(ires_nimm,ires_uc) = Jacobian(ires_nimm,ires_uc) + &
            this%mineral_n_stoich(irxn) * drate_no3_duc
        endif

        ! second term: R dn/dC_u
        if (this%is_litter_decomp(irxn)) then
          Jacobian(ires_no3,ires_uc) = Jacobian(ires_no3,ires_uc) - Rno3dn_duc

          if (this%species_id_nimm > 0) then
            Jacobian(ires_nimm,ires_uc) = Jacobian(ires_nimm,ires_uc) &
                                        + Rno3dn_duc
          endif

        endif

        ! upstream C pool
        Jacobian(ires_uc,ires_uc) = Jacobian(ires_uc,ires_uc) + drate_no3_duc

        ! upstream N pool
        if (this%is_litter_decomp(irxn)) then
          Jacobian(ires_un,ires_uc) = Jacobian(ires_un,ires_uc) + &
            this%upstream_nc(irxn) * drate_no3_duc + Rno3du_duc
        endif

        ! downstream pools
        do j = 1, this%n_downstream_pools(irxn)
          ispec_d = this%downstream_id(irxn, j)
          if (ispec_d < 0) then
            option%io_buffer = 'Downstream pool species not specified!'
            call printErrMsg(option)
          endif

          if (this%downstream_is_aqueous(irxn, j)) then
            ires_d = ispec_d
          else
            ires_d = reaction%offset_immobile + ispec_d
          endif
         
          Jacobian(ires_d,ires_uc) = Jacobian(ires_d,ires_uc) - &
            this%downstream_stoich(irxn, j) * drate_no3_duc

          if (this%is_litter_decomp(irxn) .and. &
            this%litter_decomp_type == LITTER_DECOMP_CLMMICROBE) then
            if (ispec_d == this%species_id_bacteria) then
              Jacobian(ires_d,ires_uc) = Jacobian(ires_d,ires_uc) - Rno3db_duc
            elseif (ispec_d == this%species_id_fungi) then
              Jacobian(ires_d,ires_uc) = Jacobian(ires_d,ires_uc) - Rno3df_duc
            else
              option%io_buffer = 'Downstream pool for CLM-Microbe should ' // &
                                 'be either bacteria or fungi!'
            call printErrMsg(option)
            endif
          endif
        enddo

        ! with respect to upstream n (due to variable CN ratio)
        if (this%is_litter_decomp(irxn)) then

          Jacobian(ires_un,ires_un) = Jacobian(ires_un,ires_un) + Rno3du_dun

          Jacobian(ires_no3,ires_un) = Jacobian(ires_no3,ires_un) - Rno3dn_dun

          if (this%species_id_nimm > 0) then
            Jacobian(ires_nimm,ires_un) = Jacobian(ires_nimm,ires_un) &
                                        + Rno3dn_dun
          endif

          if (this%litter_decomp_type == LITTER_DECOMP_CLMMICROBE) then
            Jacobian(ires_co2,ires_un) = Jacobian(ires_co2,ires_un) - Rno3dc_dun

            if (this%species_id_hrimm > 0) then
              Jacobian(ires_hrimm,ires_un) = Jacobian(ires_hrimm,ires_un) &
                                           - Rno3dc_dun
            endif

            Jacobian(ires_b,ires_un) = Jacobian(ires_b,ires_un) - Rno3db_dun

            Jacobian(ires_f,ires_un) = Jacobian(ires_f,ires_un) - Rno3df_dun

          endif 
        endif

        ! with respect to no3
        ! CO2
        Jacobian(ires_co2,ires_no3) = Jacobian(ires_co2,ires_no3) - &
          this%mineral_c_stoich(irxn) * drate_no3_dno3
  
        if (this%species_id_hrimm > 0) then
          Jacobian(ires_hrimm,ires_no3) = Jacobian(ires_hrimm,ires_no3) - &
            this%mineral_c_stoich(irxn) * drate_no3_dno3
        endif

        ! N
        Jacobian(ires_no3,ires_no3) = Jacobian(ires_no3,ires_no3) - &
          this%mineral_n_stoich(irxn) * drate_no3_dno3
  
        if (this%species_id_nimm > 0) then
          Jacobian(ires_nimm,ires_no3) = Jacobian(ires_nimm,ires_no3) + &
            this%mineral_n_stoich(irxn) * drate_no3_dno3
        endif

        ! upstream C pool
        Jacobian(ires_uc,ires_no3) = Jacobian(ires_uc,ires_no3) - &
          (-1.d0) * drate_no3_dno3
 
        ! upstream N pool
        if (this%is_litter_decomp(irxn)) then
          Jacobian(ires_un,ires_no3) = Jacobian(ires_un,ires_no3) - &
            (-1.d0) * this%upstream_nc(irxn) * drate_no3_dno3
        endif
  
        ! downstream pools
        do j = 1, this%n_downstream_pools(irxn)

          ispec_d = this%downstream_id(irxn, j)

          if (ispec_d < 0) then
            option%io_buffer = 'Downstream pool species not specified!'
            call printErrMsg(option)
          endif

          if (this%downstream_is_aqueous(irxn, j)) then
            ires_d = ispec_d
          else
            ires_d = reaction%offset_immobile + ispec_d
          endif
          Jacobian(ires_d,ires_no3) = Jacobian(ires_d,ires_no3) - &
            this%downstream_stoich(irxn, j) * drate_no3_dno3
        enddo

        ! with respect to nh4 (due to nh4 inhibition on no3 immobilization)
        ! CO2
        Jacobian(ires_co2,ires_nh4) = Jacobian(ires_co2,ires_nh4) - &
          this%mineral_c_stoich(irxn) * drate_no3_dnh4

        if (this%species_id_hrimm > 0) then
          Jacobian(ires_hrimm,ires_nh4) = Jacobian(ires_hrimm,ires_nh4) - &
            this%mineral_c_stoich(irxn) * drate_no3_dnh4
        endif

        ! N
        Jacobian(ires_no3,ires_nh4) = Jacobian(ires_no3,ires_nh4) - &
          this%mineral_n_stoich(irxn) * drate_no3_dnh4
  
        if (this%species_id_nimm > 0) then
          Jacobian(ires_nimm,ires_nh4) = Jacobian(ires_nimm,ires_nh4) + &
            this%mineral_n_stoich(irxn) * drate_no3_dnh4
        endif

        ! upstream C pool
        Jacobian(ires_uc,ires_nh4) = Jacobian(ires_uc,ires_nh4) - &
          (-1.d0) * drate_no3_dnh4
  
        ! upstream N pool
        if (this%is_litter_decomp(irxn)) then
          Jacobian(ires_un,ires_nh4) = Jacobian(ires_un,ires_nh4) - &
            (-1.d0) * this%upstream_nc(irxn) * drate_no3_dnh4 
        endif
  
        ! downstream pools
        do j = 1, this%n_downstream_pools(irxn)

          ispec_d = this%downstream_id(irxn, j)

          if (ispec_d < 0) then
            option%io_buffer = 'Downstream pool species not specified!'
            call printErrMsg(option)
          endif

          if (this%downstream_is_aqueous(irxn, j)) then
            ires_d = ispec_d
          else
            ires_d = reaction%offset_immobile + ispec_d
          endif
          Jacobian(ires_d,ires_nh4) = Jacobian(ires_d,ires_nh4) - &
            this%downstream_stoich(irxn, j) * drate_no3_dnh4
        enddo

        ! separate sink and source (not exist in this case) term for using NO3-
        !  as nutrient
        ! with respect to upstream
        JacobianDemand_no3(ires_co2,ires_uc) = &
          JacobianDemand_no3(ires_co2,ires_uc) - &
          this%mineral_c_stoich(irxn) * drate_no3_duc

        if (this%is_litter_decomp(irxn) .and. &
          this%litter_decomp_type == LITTER_DECOMP_CLMMICROBE) then
          ! dRco2/dLitC = dcR/dLitC = cdR/dLitC + R dc/dLitC
          JacobianDemand_no3(ires_co2,ires_uc) = &
            JacobianDemand_no3(ires_co2,ires_uc) - Rno3dc_duc
        endif

        if (this%species_id_hrimm > 0) then
          JacobianDemand_no3(ires_hrimm,ires_uc) = &
            JacobianDemand_no3(ires_hrimm,ires_uc) - &
            this%mineral_c_stoich(irxn) * drate_no3_duc

          if (this%is_litter_decomp(irxn) .and. &
            this%litter_decomp_type == LITTER_DECOMP_CLMMICROBE ) then 
            JacobianDemand_no3(ires_hrimm,ires_uc) = &
              JacobianDemand_no3(ires_hrimm,ires_uc) - Rno3dc_duc
          endif
        endif

        ! dR_N/dC_u = d nR/dC_u = dn/dC_u R + n dR/dC_u
        ! first term: n dR/dC_u
        JacobianDemand_no3(ires_no3,ires_uc) = &
          JacobianDemand_no3(ires_no3,ires_uc) - &
          this%mineral_n_stoich(irxn) * drate_no3_duc

        if (this%species_id_nimm > 0) then
          JacobianDemand_no3(ires_nimm,ires_uc) = &
            JacobianDemand_no3(ires_nimm,ires_uc) + &
            this%mineral_n_stoich(irxn) * drate_no3_duc
        endif

        ! second term: R dn/dC_u
        if (this%is_litter_decomp(irxn)) then
          JacobianDemand_no3(ires_no3,ires_uc) = &
            JacobianDemand_no3(ires_no3,ires_uc) - Rno3dn_duc

          if (this%species_id_nimm > 0) then
            JacobianDemand_no3(ires_nimm,ires_uc) = &
              JacobianDemand_no3(ires_nimm,ires_uc) + Rno3dn_duc
          endif

        endif

        ! upstream C pool
        JacobianDemand_no3(ires_uc,ires_uc) = &
          JacobianDemand_no3(ires_uc,ires_uc) + drate_no3_duc

        ! upstream N pool
        if (this%is_litter_decomp(irxn)) then
          JacobianDemand_no3(ires_un,ires_uc) = &
            JacobianDemand_no3(ires_un,ires_uc) + &
            this%upstream_nc(irxn) * drate_no3_duc + Rno3du_duc
        endif

        ! downstream pools
        do j = 1, this%n_downstream_pools(irxn)
          ispec_d = this%downstream_id(irxn, j)
          if (ispec_d < 0) then
            option%io_buffer = 'Downstream pool species not specified!'
            call printErrMsg(option)
          endif

          if (this%downstream_is_aqueous(irxn, j)) then
            ires_d = ispec_d
          else
            ires_d = reaction%offset_immobile + ispec_d
          endif
         
          JacobianDemand_no3(ires_d,ires_uc) = &
            JacobianDemand_no3(ires_d,ires_uc) - &
            this%downstream_stoich(irxn, j) * drate_no3_duc

          if (this%is_litter_decomp(irxn) .and. &
            this%litter_decomp_type == LITTER_DECOMP_CLMMICROBE) then
            if (ispec_d == this%species_id_bacteria) then
              JacobianDemand_no3(ires_d,ires_uc) = &
                JacobianDemand_no3(ires_d,ires_uc) - Rno3db_duc
            elseif (ispec_d == this%species_id_fungi) then
              JacobianDemand_no3(ires_d,ires_uc) = &
                JacobianDemand_no3(ires_d,ires_uc) - Rno3df_duc
            else
              option%io_buffer = 'Downstream pool for CLM-Microbe should ' // &
                                 'be either bacteria or fungi!'
            call printErrMsg(option)
            endif
          endif
        enddo

        ! with respect to upstream n (due to variable CN ratio)
        if (this%is_litter_decomp(irxn)) then

          JacobianDemand_no3(ires_un,ires_un) = &
            JacobianDemand_no3(ires_un,ires_un) + Rno3du_dun

          JacobianDemand_no3(ires_no3,ires_un) = &
            JacobianDemand_no3(ires_no3,ires_un) - Rno3dn_dun

          if (this%species_id_nimm > 0) then
            JacobianDemand_no3(ires_nimm,ires_un) = &
              JacobianDemand_no3(ires_nimm,ires_un) + Rno3dn_dun
          endif

          if (this%litter_decomp_type == LITTER_DECOMP_CLMMICROBE) then
            JacobianDemand_no3(ires_co2,ires_un) = &
              JacobianDemand_no3(ires_co2,ires_un) - Rno3dc_dun

            if (this%species_id_hrimm > 0) then
              JacobianDemand_no3(ires_hrimm,ires_un) = &
                JacobianDemand_no3(ires_hrimm,ires_un) - Rno3dc_dun
            endif

            JacobianDemand_no3(ires_b,ires_un) = &
              JacobianDemand_no3(ires_b,ires_un) - Rno3db_dun

            JacobianDemand_no3(ires_f,ires_un) = &
              JacobianDemand_no3(ires_f,ires_un) - Rno3df_dun

          endif 
        endif

        ! with respect to no3
        ! CO2
        JacobianDemand_no3(ires_co2,ires_no3) = &
          JacobianDemand_no3(ires_co2,ires_no3) - &
          this%mineral_c_stoich(irxn) * drate_no3_dno3
  
        if (this%species_id_hrimm > 0) then
          JacobianDemand_no3(ires_hrimm,ires_no3) = &
            JacobianDemand_no3(ires_hrimm,ires_no3) - &
            this%mineral_c_stoich(irxn) * drate_no3_dno3
        endif

        ! N
        JacobianDemand_no3(ires_no3,ires_no3) = &
          JacobianDemand_no3(ires_no3,ires_no3) - &
          this%mineral_n_stoich(irxn) * drate_no3_dno3
  
        if (this%species_id_nimm > 0) then
          JacobianDemand_no3(ires_nimm,ires_no3) = &
            JacobianDemand_no3(ires_nimm,ires_no3) + &
            this%mineral_n_stoich(irxn) * drate_no3_dno3
        endif

        ! upstream C pool
        JacobianDemand_no3(ires_uc,ires_no3) = &
          JacobianDemand_no3(ires_uc,ires_no3) - (-1.d0) * drate_no3_dno3
 
        ! upstream N pool
        if (this%is_litter_decomp(irxn)) then
          JacobianDemand_no3(ires_un,ires_no3) = &
            JacobianDemand_no3(ires_un,ires_no3) - &
            (-1.d0) * this%upstream_nc(irxn) * drate_no3_dno3
        endif
  
        ! downstream pools
        do j = 1, this%n_downstream_pools(irxn)

          ispec_d = this%downstream_id(irxn, j)

          if (ispec_d < 0) then
            option%io_buffer = 'Downstream pool species not specified!'
            call printErrMsg(option)
          endif

          if (this%downstream_is_aqueous(irxn, j)) then
            ires_d = ispec_d
          else
            ires_d = reaction%offset_immobile + ispec_d
          endif
          JacobianDemand_no3(ires_d,ires_no3) = &
            JacobianDemand_no3(ires_d,ires_no3) - &
            this%downstream_stoich(irxn, j) * drate_no3_dno3
        enddo

        ! with respect to nh4 (due to nh4 inhibition on no3 immobilization)
        ! CO2
        JacobianDemand_no3(ires_co2,ires_nh4) = &
          JacobianDemand_no3(ires_co2,ires_nh4) - &
          this%mineral_c_stoich(irxn) * drate_no3_dnh4

        if (this%species_id_hrimm > 0) then
          JacobianDemand_no3(ires_hrimm,ires_nh4) = &
            JacobianDemand_no3(ires_hrimm,ires_nh4) - &
            this%mineral_c_stoich(irxn) * drate_no3_dnh4
        endif

        ! N
        JacobianDemand_no3(ires_no3,ires_nh4) = &
          JacobianDemand_no3(ires_no3,ires_nh4) - &
          this%mineral_n_stoich(irxn) * drate_no3_dnh4
  
        if (this%species_id_nimm > 0) then
          JacobianDemand_no3(ires_nimm,ires_nh4) = &
            JacobianDemand_no3(ires_nimm,ires_nh4) + &
            this%mineral_n_stoich(irxn) * drate_no3_dnh4
        endif

        ! upstream C pool
        JacobianDemand_no3(ires_uc,ires_nh4) = &
          JacobianDemand_no3(ires_uc,ires_nh4) - (-1.d0) * drate_no3_dnh4
  
        ! upstream N pool
        if (this%is_litter_decomp(irxn)) then
          JacobianDemand_no3(ires_un,ires_nh4) = &
            JacobianDemand_no3(ires_un,ires_nh4) - &
            (-1.d0) * this%upstream_nc(irxn) * drate_no3_dnh4 
        endif
  
        ! downstream pools
        do j = 1, this%n_downstream_pools(irxn)

          ispec_d = this%downstream_id(irxn, j)

          if (ispec_d < 0) then
            option%io_buffer = 'Downstream pool species not specified!'
            call printErrMsg(option)
          endif

          if (this%downstream_is_aqueous(irxn, j)) then
            ires_d = ispec_d
          else
            ires_d = reaction%offset_immobile + ispec_d
          endif
          JacobianDemand_no3(ires_d,ires_nh4) = &
            JacobianDemand_no3(ires_d,ires_nh4) - &
            this%downstream_stoich(irxn, j) * drate_no3_dnh4
        enddo

      endif

    endif
 
  enddo

  if (this%species_id_n2o > 0) then

    f_t = 1.0d0
    f_w = 1.0d0
    f_ph = 1.0d0

    if (f_t > 1.0d-20 .and. f_w > 1.0d-20 .and. f_ph > 1.0d-20) then
      temp_real = f_t * f_w * f_ph

      if (temp_real > 1.0d0) then
        temp_real = 1.0d0
      endif

      temp_real = temp_real * this%n2o_frac_mineralization 
      
      if (net_n_mineralization_rate <= this%net_n_min_rate_smooth_0) then
        f_rate_n2o = 0.0d0
        df_rate_n2o = 0.0d0
      elseif (net_n_mineralization_rate >= this%net_n_min_rate_smooth_1 .or. &
        this%net_n_min_rate_smooth_1-this%net_n_min_rate_smooth_1 > 1.d-20) then
        f_rate_n2o = 1.0d0
        df_rate_n2o = 0.0d0
      else
        xxx = net_n_mineralization_rate - this%net_n_min_rate_smooth_0
        delta = this%net_n_min_rate_smooth_1 - this%net_n_min_rate_smooth_0
        f_rate_n2o  = 1.0d0 - (1.0d0 - xxx * xxx / delta / delta) ** 2
        df_rate_n2o = 4.0d0 * (1.0d0 - xxx * xxx / delta / delta) * xxx &
                   / delta / delta
      endif

      ! residuals 
      rate_n2o = temp_real * net_n_mineralization_rate * f_nh4 * f_rate_n2o
 
      Residual(ires_nh4) = Residual(ires_nh4) + rate_n2o 

      Residual(ires_n2o) = Residual(ires_n2o) - 0.5d0 * rate_n2o

      if (this%species_id_ngasmin > 0) then
         Residual(ires_ngasmin) = Residual(ires_ngasmin) - 0.5d0 * rate_n2o
      endif

      RateDemand_nh4(ires_nh4) = RateDemand_nh4(ires_nh4) + rate_n2o 

      RateDemand_nh4(ires_n2o) = RateDemand_nh4(ires_n2o) - 0.5d0 * rate_n2o

      if (this%species_id_ngasmin > 0) then
         RateDemand_nh4(ires_ngasmin) = RateDemand_nh4(ires_ngasmin) &
                                      - 0.5d0 * rate_n2o
      endif

      if (compute_derivative) then
        drate_n2o_dnh4 = temp_real * dnet_n_mineralization_rate_dnh4 * f_nh4 &
                       + temp_real * net_n_mineralization_rate * d_nh4
   
        drate_n2o_dnh4 = drate_n2o_dnh4 * f_rate_n2o + rate_n2o * df_rate_n2o &
                       * dnet_n_mineralization_rate_dnh4

        Jacobian(ires_nh4,ires_nh4) = Jacobian(ires_nh4,ires_nh4)+drate_n2o_dnh4
        Jacobian(ires_n2o,ires_nh4) = Jacobian(ires_n2o,ires_nh4) &
                                    - 0.5d0 * drate_n2o_dnh4
        if (this%species_id_ngasmin > 0 .and. (.not.this%bskipn2ojacobian)) then
           Jacobian(ires_ngasmin,ires_nh4) = Jacobian(ires_ngasmin,ires_nh4) &
                                           - 0.5d0 * drate_n2o_dnh4
        endif

        JacobianDemand_nh4(ires_nh4,ires_nh4) = &
          JacobianDemand_nh4(ires_nh4,ires_nh4) + drate_n2o_dnh4
        JacobianDemand_nh4(ires_n2o,ires_nh4) = &
          JacobianDemand_nh4(ires_n2o,ires_nh4) - 0.5d0 * drate_n2o_dnh4
        if (this%species_id_ngasmin > 0 .and. (.not.this%bskipn2ojacobian)) then
           JacobianDemand_nh4(ires_ngasmin,ires_nh4) = &
             JacobianDemand_nh4(ires_ngasmin,ires_nh4) - 0.5d0 * drate_n2o_dnh4
        endif

        if (this%species_id_no3 > 0) then
          drate_n2o_dno3 = temp_real * dnet_n_mineralization_rate_dno3 * f_nh4
   
          drate_n2o_dno3 = drate_n2o_dno3 *f_rate_n2o + rate_n2o * df_rate_n2o &
                  * dnet_n_mineralization_rate_dno3

          Jacobian(ires_n2o,ires_no3) = Jacobian(ires_n2o,ires_no3) &
                                      - 0.5d0 * drate_n2o_dno3

          if (this%species_id_ngasmin > 0 .and. (.not.this%bskipn2ojacobian)) then
             Jacobian(ires_ngasmin,ires_no3) = Jacobian(ires_ngasmin,ires_no3) &
                                             - 0.5d0 * drate_n2o_dno3
          endif

          Jacobian(ires_n2o,ires_no3) = Jacobian(ires_n2o,ires_no3) &
                                      - 0.5d0 * drate_n2o_dno3

          if (this%species_id_ngasmin > 0 .and. (.not.this%bskipn2ojacobian)) then
             Jacobian(ires_ngasmin,ires_no3) = Jacobian(ires_ngasmin,ires_no3) &
                                             - 0.5d0 * drate_n2o_dno3
          endif

          JacobianDemand_nh4(ires_n2o,ires_no3) = &
            JacobianDemand_nh4(ires_n2o,ires_no3) - 0.5d0 * drate_n2o_dno3

          if (this%species_id_ngasmin > 0 .and. (.not.this%bskipn2ojacobian)) then
             JacobianDemand_nh4(ires_ngasmin,ires_no3) = &
               JacobianDemand_nh4(ires_ngasmin,ires_no3) - 0.5d0 *drate_n2o_dno3
          endif

          JacobianDemand_nh4(ires_n2o,ires_no3) = &
            JacobianDemand_nh4(ires_n2o,ires_no3) - 0.5d0 * drate_n2o_dno3

          if (this%species_id_ngasmin > 0 .and. (.not.this%bskipn2ojacobian)) then
             JacobianDemand_nh4(ires_ngasmin,ires_no3) = &
               JacobianDemand_nh4(ires_ngasmin,ires_no3) - 0.5d0 *drate_n2o_dno3
          endif
        endif       

        do irxn = 1, this%nrxn
          ispec_uc = this%upstream_c_id(irxn)

          if (this%upstream_is_aqueous(irxn)) then
            ires_uc = ispec_uc
          else
            ires_uc = reaction%offset_immobile + ispec_uc
          endif
      
          drate_n2o_duc = temp_real * dnet_n_mineralization_rate_duc(irxn)*f_nh4
   
          drate_n2o_duc = drate_n2o_duc * f_rate_n2o + rate_n2o * df_rate_n2o &
                        * dnet_n_mineralization_rate_duc(irxn)

          Jacobian(ires_n2o,ires_uc) = Jacobian(ires_n2o,ires_uc) &
                                     - 0.5d0 * drate_n2o_duc

          JacobianDemand_nh4(ires_n2o,ires_uc) = &
            JacobianDemand_nh4(ires_n2o,ires_uc) - 0.5d0 * drate_n2o_duc
          if (this%species_id_ngasmin > 0 .and. (.not.this%bskipn2ojacobian)) then
            Jacobian(ires_ngasmin,ires_uc) = Jacobian(ires_ngasmin,ires_uc) &
                                           - 0.5d0 * drate_n2o_duc
            JacobianDemand_nh4(ires_ngasmin,ires_uc) = &
              JacobianDemand_nh4(ires_ngasmin,ires_uc) - 0.5d0 * drate_n2o_duc
          endif
       
        enddo

        if (this%bdebugoutput) then
          write(*, *) 'CLMDEC N2O:', rate_n2o, drate_n2o_dnh4, drate_n2o_dno3
        endif
      endif

    endif

  endif

end subroutine CLMDec_React

! **************************************************************************** !
!
! CLMDecDestroy: Destroys allocatable or pointer objects created in this module
!
! **************************************************************************** !
subroutine CLMDec_Destroy(this)

  use Utility_module, only : DeallocateArray

  implicit none
  
  class(clm_rxn_clmdec_type) :: this
  
  type(pool_type), pointer :: cur_pool, prev_pool
  type(clmdec_reaction_type), pointer :: cur_reaction, prev_reaction
  
  cur_pool => this%pools
  do
    if (.not.associated(cur_pool)) exit
    prev_pool => cur_pool
    cur_pool => cur_pool%next
    deallocate(prev_pool)
    nullify(prev_pool)
  enddo
  
  cur_reaction => this%reactions
  do
    if (.not.associated(cur_reaction)) exit

    cur_pool => cur_reaction%downstream_pools  
    do
      if (.not.associated(cur_pool)) exit
      prev_pool => cur_pool
      cur_pool => cur_pool%next
      deallocate(prev_pool)
      nullify(prev_pool)
    enddo

    prev_reaction => cur_reaction
    cur_reaction => cur_reaction%next

    deallocate(prev_reaction)
    nullify(prev_reaction)
  enddo
  
  call DeallocateArray(this%pool_nc_ratio)
  call DeallocateArray(this%rate_constant)
  call DeallocateArray(this%is_litter_decomp)
  call DeallocateArray(this%upstream_c_id)
  call DeallocateArray(this%upstream_n_id)
  call DeallocateArray(this%upstream_nc)
  call DeallocateArray(this%upstream_is_aqueous)
  call DeallocateArray(this%downstream_id)
  call DeallocateArray(this%downstream_stoich)
  call DeallocateArray(this%downstream_is_aqueous)
  call DeallocateArray(this%mineral_c_stoich) 
  call DeallocateArray(this%mineral_n_stoich) 
 
end subroutine CLMDec_Destroy

end module CLM_Rxn_Decomp_class


module CLM_Rxn_PlantN_class

  use CLM_Rxn_Base_class
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  use PFLOTRAN_Constants_module
  
! ------------------------------------------------------------------------------
! Description
! extended from reaction_sandbox_plantn to implement demand based down regulation
! for use in CLM_Rxn t6g 10/06/2014 
! add NH4+ and NO3- deposition rates as supply

! to handle plant N uptake with
! 1) Monod type downregulation N/(hs + N)
! 2) Cut off downregulation 1 if N >= N1, 0 if N <= N0, 1 - [1 - (x/d)^2]^2
!     with x = (N - N0)/(N1 - N0)
! 3) inhibition of NH3 on NO3- uptake (assuming plant take NH3 preferentially)
! Author: Guoping Tang
! Date:   07/08/14 
! -----------------------------------------------------------------------------

  implicit none
  
  private
  
#include "petsc/finclude/petscsys.h"

  type, public, &
    extends(clm_rxn_base_type) :: clm_rxn_plantn_type
    PetscReal :: rate_plantntake
    PetscReal :: rate_plantntake_nh4
    PetscReal :: rate_plantntake_no3
    PetscReal :: rate_deposition_nh4
    PetscReal :: rate_deposition_no3
    PetscReal :: half_saturation_nh4
    PetscReal :: half_saturation_no3
    PetscReal :: inhibition_nh4_no3
    PetscReal :: residual_nh4
    PetscReal :: residual_no3
    PetscReal :: cutoff_no3_0 
    PetscReal :: cutoff_no3_1
    PetscReal :: cutoff_nh4_0
    PetscReal :: cutoff_nh4_1

    PetscInt :: ispec_nh4
    PetscInt :: ispec_no3
    PetscInt :: ispec_plantn
    PetscInt :: ispec_nh4in
    PetscInt :: ispec_no3in
    PetscInt :: ispec_plantndemand

    PetscBool :: bdebugoutput
    PetscBool :: enable_clm_n_in
    PetscBool :: bfixed_clm_n_in
    PetscBool :: disable_plantntake
    PetscBool :: is_NH4_aqueous
    PetscBool :: is_NO3_aqueous
    PetscBool :: bskippno3jacobian

  contains
    procedure, public :: ReadInput => PlantNRead
    procedure, public :: Setup => PlantNSetup
    procedure, public :: Evaluate => PlantNReact
    procedure, public :: Destroy => PlantNDestroy
  end type clm_rxn_plantn_type

  public :: PlantNCreate

contains

! **************************************************************************** !
!
! PlantNCreate: Allocates plantn reaction sandbox object.
!
! **************************************************************************** !
function PlantNCreate()

  implicit none
  
  class(clm_rxn_plantn_type), pointer :: PlantNCreate

  allocate(PlantNCreate)
  PlantNCreate%rate_plantntake = 1.d-10
  PlantNCreate%rate_plantntake_nh4 = 1.d-10
  PlantNCreate%rate_plantntake_no3 = 1.d-10
  PlantNCreate%rate_deposition_nh4 = 1.0d-11
  PlantNCreate%rate_deposition_no3 = 1.0d-11
  PlantNCreate%half_saturation_nh4 =  1.d-6
  PlantNCreate%half_saturation_no3 =  1.d-6
  PlantNCreate%inhibition_nh4_no3  = -1.d-15
  PlantNCreate%residual_nh4  = 1.d-10
  PlantNCreate%residual_no3  = 1.d-10
  PlantNCreate%cutoff_no3_0 = -1.0d-9 
  PlantNCreate%cutoff_no3_1 = 1.0d-7
  PlantNCreate%cutoff_nh4_0 = -1.0d-9 
  PlantNCreate%cutoff_nh4_1 = 1.0d-7
  PlantNCreate%ispec_nh4 = -1
  PlantNCreate%ispec_no3 = -1
  PlantNCreate%ispec_plantn = -1
  PlantNCreate%ispec_nh4in = -1
  PlantNCreate%ispec_no3in = -1
  PlantNCreate%ispec_plantndemand = -1

  PlantNCreate%bdebugoutput = PETSC_FALSE
  PlantNCreate%enable_clm_n_in = PETSC_TRUE
  PlantNCreate%bfixed_clm_n_in = PETSC_FALSE
  PlantNCreate%disable_plantntake = PETSC_FALSE

  PlantNCreate%is_NH4_aqueous = PETSC_TRUE
  PlantNCreate%is_NO3_aqueous = PETSC_TRUE
  PlantNCreate%bskippno3jacobian = PETSC_FALSE

  nullify(PlantNCreate%next)  
      
end function PlantNCreate

! **************************************************************************** !
!
! PlantNRead: Reads input deck for plantn reaction parameters
!
! **************************************************************************** !
subroutine PlantNRead(this,input,option)

  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module, only : UnitsConvertToInternal
  
  implicit none
  
  class(clm_rxn_plantn_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  PetscInt :: i
  character(len=MAXWORDLENGTH) :: word, internal_units
  
  do 
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,CLM_RXN,PLANTN')
    call StringToUpper(word)   

    select case(trim(word))
      case('RATE_PLANTNTAKE_NH4')
        call InputReadDouble(input,option,this%rate_plantntake_nh4)
        call InputErrorMsg(input,option,'rate plantntake nh4+', &
          'CHEMISTRY,CLM_RXN,PLANTN')
      case('RATE_PLANTNTAKE_NO3')
        call InputReadDouble(input,option,this%rate_plantntake_no3)
        call InputErrorMsg(input,option,'rate plantntake no3-', &
          'CHEMISTRY,CLM_RXN,PLANTN')
      case('RATE_DEPOSITION_NH4')
        call InputReadDouble(input,option,this%rate_deposition_nh4)
        call InputErrorMsg(input,option,'rate deposition NH4+', &
          'CHEMISTRY,CLM_RXN,PLANTN')
      case('RATE_DEPOSITION_NO3')
        call InputReadDouble(input,option,this%rate_deposition_no3)
        call InputErrorMsg(input,option,'rate deposition NO3-', &
          'CHEMISTRY,CLM_RXN,PLANTN')
      case('HALF_SATURATION_NH4')
        call InputReadDouble(input,option,this%half_saturation_nh4)
        call InputErrorMsg(input,option,'half saturation NH4', &
          'CHEMISTRY,CLM_RXN,PLANTN')
      case('HALF_SATURATION_NO3')
        call InputReadDouble(input,option,this%half_saturation_no3)
        call InputErrorMsg(input,option,'half saturation NO3-', &
          'CHEMISTRY,CLM_RXN,PLANTN')
      case('NH4_INHIBITION_NO3')
        call InputReadDouble(input,option,this%inhibition_nh4_no3)
        call InputErrorMsg(input,option,'NH4 inhibition on NO3-', &
          'CHEMISTRY,CLM_RXN,PLANTN')
      case('RESIDUAL_NH4')
        call InputReadDouble(input,option,this%residual_nh4)
        call InputErrorMsg(input,option,'residual concentration NH4+', &
          'CHEMISTRY,CLM_RXN,PLANTN')
      case('RESIDUAL_NO3')
        call InputReadDouble(input,option,this%residual_no3)
        call InputErrorMsg(input,option,'residual concentration NO3-', &
          'CHEMISTRY,CLM_RXN,PLANTN')
      case('CUTOFF_NH4')
        call InputReadDouble(input,option,this%cutoff_nh4_0)
        call InputErrorMsg(input,option,'cutoff_nh4_0', &
          'CHEMISTRY,CLM_RXN,PLANTN')
        call InputReadDouble(input,option,this%cutoff_nh4_1)
        call InputErrorMsg(input,option,'cutoff_nh4_1', &
          'CHEMISTRY,CLM_RXN,PLANTN')
        if (this%cutoff_nh4_0 > this%cutoff_nh4_1) then
          option%io_buffer = 'CHEMISTRY,CLM_RXN,PLANTN,' // &
            'NH4+ cut off concentration > concentration ' // &
            'where down regulation function = 1.'
          call printErrMsg(option)
        endif
      case('CUTOFF_NO3')
        call InputReadDouble(input,option,this%cutoff_no3_0)
        call InputErrorMsg(input,option,'cutoff_no3_0', &
          'CHEMISTRY,CLM_RXN,PLANTN')
        call InputReadDouble(input,option,this%cutoff_no3_1)
        call InputErrorMsg(input,option,'cutoff_no3_1', &
          'CHEMISTRY,CLM_RXN,PLANTN')
        if (this%cutoff_no3_0 > this%cutoff_no3_1) then
          option%io_buffer = 'CHEMISTRY,CLM_RXN,PLANTN,' // &
            'NO3- down regulation cut off concentration > concentration ' // &
            'where down regulation function = 1.'
          call printErrMsg(option)
        endif
      case('DEBUG_OUTPUT')
        this%bdebugoutput = PETSC_TRUE
      case('DISABLE_CLM_N_INPUT')
        this%enable_clm_n_in = PETSC_FALSE
      case('FIXED_CLM_N_INPUT')
        this%bfixed_clm_n_in = PETSC_TRUE
      case('DISABLE_PLANTNTAKE')
        this%disable_plantntake = PETSC_TRUE
      case('JACOBIAN_PLANT_NO3_SKIP')
        this%bskippno3jacobian = PETSC_TRUE
      case default
        call InputKeywordUnrecognized(word, &
               'CHEMISTRY,CLM_RXN,PLANTN,REACTION',option)
    end select
  enddo
  
end subroutine PlantNRead

! **************************************************************************** !
!
! PlantNSetup: Sets up the plantn reaction with parameters
!
! **************************************************************************** !
subroutine PlantNSetup(this,reaction,option)

  use Reaction_Aux_module
  use Option_module
  use Reaction_Immobile_Aux_module

  implicit none
  
  class(clm_rxn_plantn_type) :: this
  type(reaction_type) :: reaction
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word
 
  word = 'NH4+'
  this%ispec_nh4 = GetPrimarySpeciesIDFromName(word,reaction,PETSC_FALSE,option)

  if (this%ispec_nh4 < 0) then
    word = 'NH3(aq)'
    this%ispec_nh4 = GetPrimarySpeciesIDFromName(word,reaction,PETSC_FALSE, &
      option)
  endif
  
  if (this%ispec_nh4 < 0) then
    word = 'Ammonium'
    this%ispec_nh4 = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)
    if (this%ispec_nh4 > 0) then
      this%is_NH4_aqueous = PETSC_FALSE
    endif
  endif 

  if (this%ispec_nh4 < 0) then
    option%io_buffer = 'NH4+, NH3(aq) or Ammonium is specified in the input' // &
      'file for PlantN sandbox!'
    call printErrMsg(option)
  endif

  word = 'NO3-'
  this%ispec_no3 = GetPrimarySpeciesIDFromName(word,reaction,PETSC_FALSE,option)

  if (this%ispec_no3 < 0) then
    word = 'Nitrate'
    this%ispec_no3 = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)
    if (this%ispec_no3 > 0) then
      this%is_NO3_aqueous = PETSC_FALSE
    endif
  endif 

  word = 'PlantN'
  this%ispec_plantn = GetImmobileSpeciesIDFromName(word, reaction%immobile, &
    PETSC_FALSE,option)

  if (this%ispec_plantn < 0) then
    option%io_buffer = 'PlantN is specified in the input file!'
    call printErrMsg(option)
  endif

  word = 'Ain'
  this%ispec_nh4in = GetImmobileSpeciesIDFromName(word, reaction%immobile, &
    PETSC_FALSE,option)

  word = 'Tin'
  this%ispec_no3in = GetImmobileSpeciesIDFromName(word, reaction%immobile, &
    PETSC_FALSE,option)

  word = 'Plantndemand'
  this%ispec_plantndemand = GetImmobileSpeciesIDFromName(word, &
    reaction%immobile, PETSC_FALSE,option)

end subroutine PlantNSetup

! **************************************************************************** !
!
! PlantNReact: Evaluates reaction storing residual and/or Jacobian
!
! **************************************************************************** !
subroutine PlantNReact(this,Residual,Jacobian,compute_derivative,rt_auxvar, &
                       global_auxvar,material_auxvar,reaction,option, &
                       RateDemand_nh4,RateSupply_nh4, &
                       JacobianDemand_nh4,JacobianSupply_nh4, &
                       RateDemand_no3,RateSupply_no3, &
                       JacobianDemand_no3,JacobianSupply_no3, &
                       Rate_nh4_to_no3,Jacobian_nh4_to_no3)

  use Option_module
  use Reaction_Aux_module
  use Reaction_Immobile_Aux_module
  use Material_Aux_class, only : material_auxvar_type
  use CLM_Rxn_Common_module, only: CalNLimitFunc

  implicit none

  class(clm_rxn_plantn_type) :: this  
  type(option_type) :: option
  type(reaction_type) :: reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscBool :: compute_derivative

  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  PetscReal :: RateDemand_nh4(reaction%ncomp)
  PetscReal :: RateSupply_nh4(reaction%ncomp)
  PetscReal :: JacobianDemand_nh4(reaction%ncomp,reaction%ncomp)
  PetscReal :: JacobianSupply_nh4(reaction%ncomp,reaction%ncomp)
  PetscReal :: RateDemand_no3(reaction%ncomp)
  PetscReal :: RateSupply_no3(reaction%ncomp)
  PetscReal :: JacobianDemand_no3(reaction%ncomp,reaction%ncomp)
  PetscReal :: JacobianSupply_no3(reaction%ncomp,reaction%ncomp)
  PetscReal :: Rate_nh4_to_no3
  PetscReal :: Jacobian_nh4_to_no3(reaction%ncomp)
  PetscReal :: volume, porosity
  PetscErrorCode :: ierr
  PetscInt :: local_id

  character(len=MAXWORDLENGTH) :: word

  PetscInt, parameter :: iphase = 1
  PetscInt :: ires_nh4, ires_no3, ires_plantn
  PetscInt :: ires_nh4in, ires_no3in
  PetscInt :: ires_plantndemand

  PetscReal :: c_nh4         ! concentration (mole/L)
  PetscReal :: ac_nh4        ! activity coefficient
  PetscReal :: f_nh4         ! nh4 / (half_saturation + nh4)
  PetscReal :: d_nh4         ! half_saturation / (half_saturation + nh4)^2
  PetscReal :: f_nh4_inhibit ! inhibition_coef/(inhibition_coef + nh4)
  PetscReal :: d_nh4_inhibit ! d inhibition_coef/(inhibition_coef + nh4)
  PetscReal :: c_no3         ! concentration (mole/L)
  PetscReal :: ac_no3        ! activity coefficient
  PetscReal :: f_no3         ! no3 / (half_saturation + no3)
  PetscReal :: d_no3         ! half_saturation/(no3 + half_saturation)^2 
  PetscReal :: temp_real

  PetscReal :: rate_plantn
  PetscReal :: rate_nh4
  PetscReal :: rate_no3
  PetscReal :: drate_nh4_dnh4
  PetscReal :: drate_no3_dno3
  PetscReal :: drate_no3_dnh4
  PetscReal :: c_plantn, c_plantno3, c_plantnh4, c_plantndemand
  PetscReal :: xxx, delta, regulator, dregulator
  PetscReal :: rate_nh4_clm_input, rate_no3_clm_input

  porosity = material_auxvar%porosity
  volume = material_auxvar%volume

  ires_nh4 = this%ispec_nh4
  ires_no3 = this%ispec_no3
  ires_plantn = this%ispec_plantn + reaction%offset_immobile
  ires_nh4in = this%ispec_nh4in + reaction%offset_immobile
  ires_no3in = this%ispec_no3in + reaction%offset_immobile
  ires_plantndemand = this%ispec_plantndemand + reaction%offset_immobile

  if (this%ispec_plantn < 0) then
    option%io_buffer = 'PlantN is not specified in the input file!'
    call printErrMsg(option)
  endif

  ires_nh4 = -999
  ires_no3 = -999

  f_nh4 = 1.0d0
  d_nh4 = 0.0d0

  if (this%ispec_nh4 > 0) then
    if (this%is_NH4_aqueous) then   
      c_nh4    = rt_auxvar%pri_molal(this%ispec_nh4)
      ac_nh4   = rt_auxvar%pri_act_coef(this%ispec_nh4)
      ires_nh4 = this%ispec_nh4
    else
      c_nh4    = rt_auxvar%immobile(this%ispec_nh4)
      ac_nh4   = 1.0d0
      ires_nh4 = this%ispec_nh4 + reaction%offset_immobile
    endif

    call CalNLimitFunc(c_nh4, ac_nh4, this%residual_nh4, &
      this%half_saturation_nh4, this%cutoff_nh4_0, this%cutoff_nh4_1, &
      f_nh4, d_nh4)
  endif

  f_no3 = 1.0d0
  d_no3 = 0.0d0

  f_nh4_inhibit = 1.0d0
  d_nh4_inhibit = 0.0d0

  if (this%ispec_no3 > 0) then
    if (this%is_NO3_aqueous) then   
      c_no3     = rt_auxvar%pri_molal(this%ispec_no3)
      ac_no3    = rt_auxvar%pri_act_coef(this%ispec_no3)
      ires_no3 = this%ispec_no3
    else
      c_no3    = rt_auxvar%immobile(this%ispec_no3)
      ac_no3   = 1.0d0
      ires_no3 = this%ispec_no3 + reaction%offset_immobile
    endif

    call CalNLimitFunc(c_no3, ac_no3, this%residual_no3, &
      this%half_saturation_no3, this%cutoff_no3_0, this%cutoff_no3_1, &
      f_no3, d_no3)

    if (this%ispec_nh4 > 0 .and. &
      this%inhibition_nh4_no3 > this%residual_nh4) then
      temp_real = this%inhibition_nh4_no3 + c_nh4 * ac_nh4
      f_nh4_inhibit = this%inhibition_nh4_no3/temp_real
      d_nh4_inhibit = -1.0d0 * this%inhibition_nh4_no3 * ac_nh4 &
                    / temp_real / temp_real
    endif
  endif

  if (this%inhibition_nh4_no3 > this%residual_nh4) then
    rate_nh4 = this%rate_plantntake * volume
    rate_no3 = this%rate_plantntake * volume
  else
    rate_nh4 = this%rate_plantntake_nh4 * volume
    rate_no3 = this%rate_plantntake_no3 * volume
  endif

  rate_nh4_clm_input = this%rate_deposition_nh4 * volume
  rate_no3_clm_input = this%rate_deposition_no3 * volume

  if (this%ispec_plantndemand > 0) then
    Residual(ires_plantndemand) = Residual(ires_plantndemand) - rate_nh4
  endif

  if (this%ispec_nh4 > 0) then

    if (compute_derivative) then
      drate_nh4_dnh4 = rate_nh4 * d_nh4
    endif
    rate_nh4 = rate_nh4 * f_nh4

    Residual(ires_nh4) = Residual(ires_nh4) + rate_nh4
    Residual(ires_plantn) = Residual(ires_plantn) - rate_nh4

    if (this%ispec_nh4in > 0) then
      Residual(ires_nh4in) = Residual(ires_nh4in) - rate_nh4
    endif

    RateDemand_nh4(ires_nh4) = RateDemand_nh4(ires_nh4) + rate_nh4 
    RateDemand_nh4(ires_plantn) = RateDemand_nh4(ires_plantn) - rate_nh4

    if (this%ispec_nh4in > 0) then
      RateDemand_nh4(ires_nh4in) = RateDemand_nh4(ires_nh4in) - rate_nh4
    endif

    if (compute_derivative) then
      Jacobian(ires_nh4,ires_nh4) = Jacobian(ires_nh4,ires_nh4) + drate_nh4_dnh4

      Jacobian(ires_plantn,ires_nh4) = Jacobian(ires_plantn,ires_nh4) &
                                     - drate_nh4_dnh4

      if (this%ispec_nh4in > 0) then
        Jacobian(ires_nh4in,ires_nh4) = Jacobian(ires_nh4in,ires_nh4) &
                                      - drate_nh4_dnh4
      endif

      JacobianDemand_nh4(ires_nh4,ires_nh4) = &
        JacobianDemand_nh4(ires_nh4,ires_nh4) + drate_nh4_dnh4

      JacobianDemand_nh4(ires_plantn,ires_nh4) = &
        JacobianDemand_nh4(ires_plantn,ires_nh4) - drate_nh4_dnh4

      if (this%ispec_nh4in > 0) then
        JacobianDemand_nh4(ires_nh4in,ires_nh4) = &
          JacobianDemand_nh4(ires_nh4in,ires_nh4) - drate_nh4_dnh4
      endif

    endif
  endif

  if (this%ispec_no3 > 0) then
    if (compute_derivative) then
      drate_no3_dno3 = rate_no3 * f_nh4_inhibit * d_no3
      drate_no3_dnh4 = rate_no3 * d_nh4_inhibit * f_no3
    endif

    rate_no3 = rate_no3 * f_nh4_inhibit * f_no3

    Residual(ires_no3) = Residual(ires_no3) + rate_no3
    Residual(ires_plantn) = Residual(ires_plantn) - rate_no3

    if (this%ispec_no3in > 0) then
      Residual(ires_no3in) = Residual(ires_no3in) - rate_no3
    endif

    RateDemand_no3(ires_no3) = RateDemand_no3(ires_no3) + rate_no3 
    RateDemand_no3(ires_plantn) = RateDemand_no3(ires_plantn) - rate_no3

    if (this%ispec_no3in > 0) then
      RateDemand_no3(ires_no3in) = RateDemand_no3(ires_no3in) - rate_no3
    endif

    if (compute_derivative) then
      Jacobian(ires_no3,ires_no3) = Jacobian(ires_no3,ires_no3) + drate_no3_dno3

      Jacobian(ires_plantn,ires_no3) = Jacobian(ires_plantn,ires_no3) &
                                     - drate_no3_dno3

      if (this%ispec_no3in > 0 .and. (.not.this%bskippno3jacobian)) then
        Jacobian(ires_no3in,ires_no3) = Jacobian(ires_no3in,ires_no3) &
                                      - drate_no3_dno3
      endif

      Jacobian(ires_no3,ires_nh4) = Jacobian(ires_no3,ires_nh4) + drate_no3_dnh4

      Jacobian(ires_plantn,ires_nh4) = Jacobian(ires_plantn,ires_nh4) &
                                     - drate_no3_dnh4

      if (this%ispec_no3in > 0 .and. (.not.this%bskippno3jacobian)) then
        Jacobian(ires_no3in,ires_nh4) = Jacobian(ires_no3in,ires_nh4) &
                                      - drate_no3_dnh4
      endif

      JacobianDemand_no3(ires_no3,ires_no3) = &
        JacobianDemand_no3(ires_no3,ires_no3) + drate_no3_dno3

      JacobianDemand_no3(ires_plantn,ires_no3) = &
        JacobianDemand_no3(ires_plantn,ires_no3) - drate_no3_dno3

      if (this%ispec_no3in > 0 .and. (.not.this%bskippno3jacobian)) then
        JacobianDemand_no3(ires_no3in,ires_no3) = &
          JacobianDemand_no3(ires_no3in,ires_no3) - drate_no3_dno3
      endif

      JacobianDemand_no3(ires_no3,ires_nh4) = &
        JacobianDemand_no3(ires_no3,ires_nh4) + drate_no3_dnh4

      JacobianDemand_no3(ires_plantn,ires_nh4) = &
        JacobianDemand_no3(ires_plantn,ires_nh4) - drate_no3_dnh4

      if (this%ispec_no3in > 0 .and. (.not.this%bskippno3jacobian)) then
        JacobianDemand_no3(ires_no3in,ires_nh4) = &
          JacobianDemand_no3(ires_no3in,ires_nh4) - drate_no3_dnh4
      endif

      if (this%bdebugoutput) then
        c_plantn = rt_auxvar%immobile(this%ispec_plantn)
        write(*, *) c_nh4, c_no3, rate_plantn, rate_nh4_clm_input, &
          rate_no3_clm_input, rate_nh4, rate_no3, c_plantn
      endif
    endif
  endif

end subroutine PlantNReact

! **************************************************************************** !
!
! PlantNDestroy: Destroys allocatable or pointer objects created in this module
!
! **************************************************************************** !
subroutine PlantNDestroy(this)

  implicit none
  
  class(clm_rxn_plantn_type) :: this  

end subroutine PlantNDestroy

end module CLM_Rxn_PlantN_class


module CLM_Rxn_Nitr_class

! ------------------------------------------------------------------------------
! Description
! nitrification function following Dickinson et al. 2002
!   NH4+ -> NO3-
!   rate   = kmax ftheta fT NH4+
!   fT     = exp(0.08(T - 298))
!   ftheta = s (1 - s) / (0.25 + 1 / NH4+) 
! and Parton et al 1996
!   NH4+ -> 0.5 N2O
!   rate   = kmax ftheta fT fpH (1 - exp(-0.0104e6mN rhob/theta  NH4+)
! by t6g 10/06/2014 
!   1/(0.25 + 1 / NH4+) = 4 NH4+ /(NH4+ + 4)
!   simplifies to the general Monod function, add DICKINSON if not 
!   1 - exp(-x) = x + ... (remove high order terms)
!   simplify to first order rate, add PARTON if not
! by t6g 2/13/2015
! ------------------------------------------------------------------------------

  use CLM_Rxn_Base_class
  
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  use PFLOTRAN_Constants_module
  
  implicit none
  
  private
  
#include "petsc/finclude/petscsys.h"

  PetscInt, parameter :: TEMPERATURE_RESPONSE_FUNCTION_CLM4 = 1
  PetscInt, parameter :: TEMPERATURE_RESPONSE_FUNCTION_Q10 = 2

  type, public, &
    extends(clm_rxn_base_type) :: clm_rxn_nitr_type
    PetscInt :: ispec_proton
    PetscInt :: ispec_nh4
    PetscInt :: ispec_nh4sorb
    PetscInt :: ispec_no3
    PetscInt :: ispec_n2o
    PetscInt :: ispec_ngasnit
    PetscReal :: k_nitr_max
    PetscReal :: k_nitr_n2o
    PetscInt :: temperature_response_function
    PetscReal :: Q10
    PetscReal :: residual_conc
    PetscReal :: half_saturation
    PetscReal :: cutoff_nh4_0  ! shut off
    PetscReal :: cutoff_nh4_1  ! start to decrease from 1
    PetscReal :: c_nh4_ugg_0      
    PetscReal :: c_nh4_ugg_1    ! N2O production from nitr (Parton et al. 1996) 
    PetscBool :: disable_mrf    ! for testing purpose 
    PetscBool :: bdebugoutput
    ! to use 1/(0.25 + 1/NH4+) rather than the simple Monod substrate limiting function
    PetscBool :: bDickinson
    ! to use (1 - exp(-0.0104e6mN rhob/theta  NH4+) rather than first order
    PetscBool :: bParton
    PetscBool :: is_NH4_aqueous
    PetscBool :: is_NO3_aqueous
    PetscBool :: bskipnitrjacobian

  contains
    procedure, public :: ReadInput => NitrRead
    procedure, public :: Setup => NitrSetup
    procedure, public :: Evaluate => NitrReact
    procedure, public :: Destroy => NitrDestroy
  end type clm_rxn_nitr_type

  public :: NitrCreate

contains

! ************************************************************************** !
!
! NitrCreate: Allocates nitr reaction object.
!
! ************************************************************************** !
function NitrCreate()

  implicit none
  
  class(clm_rxn_nitr_type), pointer :: NitrCreate

  allocate(NitrCreate)
  NitrCreate%ispec_proton = 0
  NitrCreate%ispec_nh4 = 0
  NitrCreate%ispec_nh4sorb = 0
  NitrCreate%ispec_no3 = 0
  NitrCreate%ispec_ngasnit = 0
  NitrCreate%k_nitr_max = 1.d-6
  NitrCreate%k_nitr_n2o = 3.5d-8
  NitrCreate%temperature_response_function = TEMPERATURE_RESPONSE_FUNCTION_CLM4
  NitrCreate%Q10 = 1.5d0
  NitrCreate%residual_conc = 1.0d-10
  NitrCreate%half_saturation = -1.0d-6 
  NitrCreate%cutoff_nh4_0 =-1.0d-20 
  NitrCreate%cutoff_nh4_1 = 1.0d-20
  NitrCreate%c_nh4_ugg_0 = 2.9d0
  NitrCreate%c_nh4_ugg_1 = 3.0d0   ! N2O production from nitr (Parton et al. 1996) 
  NitrCreate%disable_mrf = PETSC_FALSE
  NitrCreate%bdebugoutput = PETSC_FALSE
  NitrCreate%bDickinson = PETSC_FALSE
  NitrCreate%bParton = PETSC_FALSE
  NitrCreate%is_NH4_aqueous = PETSC_TRUE
  NitrCreate%is_NO3_aqueous = PETSC_TRUE
  NitrCreate%bskipnitrjacobian = PETSC_FALSE
  nullify(NitrCreate%next)  
      
end function NitrCreate

! ************************************************************************** !
!
! NitrRead: Reads input deck for nitr reaction parameters (if any)
!
! ************************************************************************** !
subroutine NitrRead(this,input,option)

  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module, only : UnitsConvertToInternal
  
  implicit none
  
  class(clm_rxn_nitr_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  PetscInt :: i
  character(len=MAXWORDLENGTH) :: word, internal_units
  
  do 
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,CLM_RXN,NITRIFICATION')
    call StringToUpper(word)   

    select case(trim(word))
      case('TEMPERATURE_RESPONSE_FUNCTION')
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit

          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword', &
            'CHEMISTRY,CLM_RXN,NITRIFICATION,TEMPERATURE RESPONSE FUNCTION')
          call StringToUpper(word)   

          select case(trim(word))
            case('CLM4')
              this%temperature_response_function = &
                TEMPERATURE_RESPONSE_FUNCTION_CLM4 
            case('Q10')
              this%temperature_response_function = &
                TEMPERATURE_RESPONSE_FUNCTION_Q10    
              call InputReadDouble(input,option,this%Q10)  
              call InputErrorMsg(input,option,'Q10', &
                'CHEMISTRY,CLM_RXN_NITRIFICATION,TEMPERATURE RESPONSE FUNCTION')
            case default
              call InputKeywordUnrecognized(word, &
                'CHEMISTRY,CLM_RXN,NITRIFICATION,TEMPERATURE RESPONSE FUNCTION', &
                option)
          end select
        enddo 
      case('RATE_CONSTANT_NO3')
        call InputReadDouble(input,option,this%k_nitr_max)
        call InputErrorMsg(input,option,'nitr rate coefficient', &
          'CHEMISTRY,CLM_RXN,NITRIFICATION,REACTION')
      case('RATE_CONSTANT_N2O')
        call InputReadDouble(input,option,this%k_nitr_n2o)
        call InputErrorMsg(input,option,'N2O rate coefficient', &
                     'CHEMISTRY,CLM_RXN,NITRIFICATION,REACTION')
      case('RESIDUAL_NH4')
        call InputReadDouble(input,option,this%residual_conc)
        call InputErrorMsg(input,option,'residual NH4+', &
                  'CHEMISTRY,CLM_RXN,NITRIFICATION,REACTION')
      case('HALF_SATURATION_NH4')
        call InputReadDouble(input,option,this%half_saturation)
        call InputErrorMsg(input,option,'half saturation NH4+', &
                  'CHEMISTRY,CLM_RXN,NITRIFICATION,REACTION')
      case('CUTOFF_NH4')
        call InputReadDouble(input,option,this%cutoff_nh4_0)
        call InputErrorMsg(input,option,'cutoff_nh4_0', &
          'CHEMISTRY,CLM_RXN,NITRIFICATION,REACTION')
        call InputReadDouble(input,option,this%cutoff_nh4_1)
        call InputErrorMsg(input,option,'cutoff_nh4_1', &
          'CHEMISTRY,CLM_RXN,NITRIFICATION,REACTION')
        if (this%cutoff_nh4_0 > this%cutoff_nh4_1) then
          option%io_buffer = 'CHEMISTRY,CLM_RXN,NITRIFICATION,' // &
            'NH4+ cut off concentration > concentration ' // &
            'where down regulation function = 1.'
          call printErrMsg(option)
        endif
      case('SMOOTH_NH4_2_N2O')
        call InputReadDouble(input,option,this%c_nh4_ugg_0)
        call InputErrorMsg(input,option,'c_nh4_ugg_0', &
          'CHEMISTRY,CLM_RXN,NITRIFICATION,REACTION')
        call InputReadDouble(input,option,this%c_nh4_ugg_1)
        call InputErrorMsg(input,option,'c_nh4_ugg_1', &
          'CHEMISTRY,CLM_RXN,NITRIFICATION,REACTION')
      case('DICKINSON')
        this%bDickinson = PETSC_TRUE
      case('PARTON')
        this%bParton = PETSC_TRUE
      case('DISABLE_MRF')
        this%disable_mrf = PETSC_TRUE
      case('DEBUG_OUTPUT')
        this%bdebugoutput = PETSC_TRUE
      case('JACOBIAN_NITR_SKIP')
        this%bskipnitrjacobian = PETSC_TRUE
      case default
        call InputKeywordUnrecognized(word, &
                'CHEMISTRY,CLM_RXN,NITRIFICATION,REACTION',option)
    end select
  enddo
  
end subroutine NitrRead

! ************************************************************************** !
!
! NitrSetup: Sets up the nitr reaction either with parameters either
!                read from the input deck or hardwired.
!
! ************************************************************************** !
subroutine NitrSetup(this,reaction,option)

  use Reaction_Aux_module, only : reaction_type, GetPrimarySpeciesIDFromName
  use Option_module
  use Reaction_Immobile_Aux_module, only : GetImmobileSpeciesIDFromName 

  implicit none
  
  class(clm_rxn_nitr_type) :: this
  type(reaction_type) :: reaction
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word

  word = 'H+'
  this%ispec_proton = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)

  word = 'NH3(aq)'
  this%ispec_nh4 = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)

  if (this%ispec_nh4 < 0) then
    word = 'NH4+'
    this%ispec_nh4 = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)
  endif

  if (this%ispec_nh4 < 0) then
    word = 'Ammonium'
    this%ispec_nh4 = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)
    if (this%ispec_nh4 > 0) then
      this%is_NH4_aqueous = PETSC_FALSE
    endif
  endif 

  word = 'NO3-'
  this%ispec_no3 = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)

  if (this%ispec_no3 < 0) then
    word = 'Nitrate'
    this%ispec_no3 = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)
    if (this%ispec_no3 > 0) then
      this%is_NO3_aqueous = PETSC_FALSE
    endif
  endif 

  word = 'N2O(aq)'
  this%ispec_n2o = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)
  if (this%ispec_n2o < 0) then
    word = 'NO2-'
    this%ispec_n2o = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)
  endif

  if (this%ispec_nh4 < 0) then
     option%io_buffer = 'CHEMISTRY,CLM_RXN,NITRIFICATION: ' // &
       'NH3(aq), NH4+, or Ammonium is not specified in the input file.'
     call printErrMsg(option)
  endif

  if (this%ispec_no3 < 0) then
     option%io_buffer = 'CHEMISTRY,CLM_RXN,NITRIFICATION: ' // &
                        ' NO3- is not specified in the input file.'
     call printErrMsg(option)
  endif

!  if (this%ispec_n2o < 0) then
!     option%io_buffer = 'CHEMISTRY,CLM_RXN,NITRIFICATION: ' // &
!                        ' N2O(aq) is not specified in the input file.'
!     call printErrMsg(option)
!  endif

  word = 'NGASnitr'
  this%ispec_ngasnit = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)
 
end subroutine NitrSetup

! ************************************************************************** !
!
! NitrReact: Evaluates reaction storing residual and/or Jacobian
!
! ************************************************************************** !
subroutine NitrReact(this,Residual,Jacobian,compute_derivative, &
                     rt_auxvar,global_auxvar,material_auxvar,reaction,option, &
                     RateDemand_nh4,RateSupply_nh4, &
                     JacobianDemand_nh4,JacobianSupply_nh4, &
                     RateDemand_no3,RateSupply_no3, &
                     JacobianDemand_no3,JacobianSupply_no3, &
                     Rate_nh4_to_no3,Jacobian_nh4_to_no3)
  use Option_module
  use Reaction_Aux_module
  use Material_Aux_class, only : material_auxvar_type
  use CLM_Rxn_Common_module, only: CalNLimitFunc
  implicit none

  class(clm_rxn_nitr_type) :: this  
  type(option_type) :: option
  type(reaction_type) :: reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar

  PetscBool :: compute_derivative
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  PetscReal :: RateDemand_nh4(reaction%ncomp)
  PetscReal :: RateSupply_nh4(reaction%ncomp)
  PetscReal :: JacobianDemand_nh4(reaction%ncomp,reaction%ncomp)
  PetscReal :: JacobianSupply_nh4(reaction%ncomp,reaction%ncomp)
  PetscReal :: RateDemand_no3(reaction%ncomp)
  PetscReal :: RateSupply_no3(reaction%ncomp)
  PetscReal :: JacobianDemand_no3(reaction%ncomp,reaction%ncomp)
  PetscReal :: JacobianSupply_no3(reaction%ncomp,reaction%ncomp)
  PetscReal :: Rate_nh4_to_no3
  PetscReal :: Jacobian_nh4_to_no3(reaction%ncomp)
  PetscReal :: porosity
  PetscReal :: volume
  PetscInt :: local_id
  PetscErrorCode :: ierr

  PetscInt, parameter :: iphase = 1
  PetscReal, parameter :: rpi = 3.14159265358979323846
  PetscReal, parameter :: N_molecular_weight = 14.0067d0
  PetscReal :: M_2_ug_per_g
  PetscReal :: mol_m3_2_ug_per_g

  PetscInt :: ires_nh4, ires_no3, ires_n2o

  PetscScalar, pointer :: bulkdensity(:)
  PetscReal :: rho_b
  PetscReal :: theta
  PetscReal :: c_nh4      ! mole/L
  PetscReal :: ac_nh4      ! mole/L
  PetscReal :: s_nh4      ! mole/m3
  PetscReal :: c_nh4_ugg  ! ug ammonia N / g soil
  PetscReal :: ph
  PetscReal :: rate_n2o, drate_n2o
  PetscReal :: rate_nitri, drate_nitri
  PetscReal :: f_t, f_w, f_ph
  PetscReal :: dfw_dnh4
  PetscReal :: saturation
  PetscReal :: tc
  PetscReal :: kg_water
  PetscReal :: h2osoi
  PetscInt :: ires_ngasnit
  PetscReal :: xxx, delta, regulator, dregulator
  PetscReal :: f_nh4, d_nh4     ! for monod substrate limitation
  PetscReal :: f_n2o, d_n2o     ! for smoothing N2O production
  PetscReal :: c_nh4_0, c_nh4_1
  PetscReal :: temp_real
  PetscReal :: unitconv

  porosity = material_auxvar%porosity
  volume = material_auxvar%volume
  kg_water = material_auxvar%porosity*global_auxvar%sat(iphase)* &
            material_auxvar%volume*global_auxvar%den_kg(iphase)

  unitconv = 1.0d0
  M_2_ug_per_g = 1.0d0
  mol_m3_2_ug_per_g = 1.0d0
  h2osoi = 1.0d0

  local_id = option%iflag
  ! indices for C and N species
  ires_nh4 = this%ispec_nh4
  ires_no3 = this%ispec_no3
  ires_n2o = this%ispec_n2o
  ires_ngasnit = this%ispec_ngasnit + reaction%offset_immobile

  saturation = global_auxvar%sat(1)
  theta = saturation * porosity

  tc = global_auxvar%temp

  if (this%is_NH4_aqueous) then   
    c_nh4    = rt_auxvar%pri_molal(this%ispec_nh4)
    ac_nh4   = rt_auxvar%pri_act_coef(this%ispec_nh4)
    ires_nh4 = this%ispec_nh4
  else
    c_nh4    = rt_auxvar%immobile(this%ispec_nh4)
    ac_nh4   = 1.0d0
    ires_nh4 = this%ispec_nh4 + reaction%offset_immobile
  endif

  call CalNLimitFunc(c_nh4, ac_nh4, this%residual_conc, this%half_saturation, &
                     this%cutoff_nh4_0, this%cutoff_nh4_1, f_nh4, d_nh4)

  if (associated(rt_auxvar%total_sorb_eq)) then           ! original absorption-reactions in PF used
     s_nh4 = rt_auxvar%total_sorb_eq(this%ispec_nh4)
  else
     s_nh4 = 1.d-20
  endif

  c_nh4 = (c_nh4 - this%residual_conc) * ac_nh4

  ! nitrification (Dickinson et al. 2002)
  if (this%ispec_no3 > 0) then
    f_t = exp(0.08d0 * (tc - 298.0d0 + 273.15d0))

    if (tc < 0.0d0) f_t = 0.0d0   ! to be consistent in CLM  CNNDynamicsMod.F90 line 839

    if (this%disable_mrf) then
      f_w = 1.0d0
    else
      f_w = saturation * (1.0d0 - saturation)
    endif

    if (this%is_NH4_aqueous) then   
      temp_real = f_t * f_w * this%k_nitr_max * kg_water
    else
      temp_real = f_t * f_w * this%k_nitr_max * volume
    endif

    if (this%bDickinson) then
      ! to make is consistent with clm CNNDynamicsMod.F90 line 832
      if (this%is_NH4_aqueous) then
        unitconv = N_molecular_weight * h2osoi * 1000.0d0  ! from mol/L to g/m^3
      else
        unitconv = N_molecular_weight ! from mol/m^3 to g/m^3
      endif
      rate_nitri = temp_real * c_nh4 * c_nh4 / (0.25d0 * c_nh4 + 1.0d0/unitconv)
    else
      rate_nitri = temp_real * c_nh4
    endif

    Residual(ires_nh4) = Residual(ires_nh4) + rate_nitri * f_nh4
    Residual(ires_no3) = Residual(ires_no3) - rate_nitri * f_nh4

    RateDemand_nh4(ires_nh4) = RateDemand_nh4(ires_nh4) + rate_nitri * f_nh4
    RateDemand_nh4(ires_no3) = RateDemand_nh4(ires_no3) - rate_nitri * f_nh4

    Rate_nh4_to_no3 = Rate_nh4_to_no3 + rate_nitri * f_nh4

    if (compute_derivative) then
      if (this%bDickinson) then
        ! f = x^2/(x/4+1/u)
        ! f' = [2x(x/4+1/u) - x^2/4]/(x/4 + 1/u)^2 
        !    = (x^2/4 + 2x/u)/(x/4 + 1/u)^2
        drate_nitri = temp_real &
                  * (0.25d0 * c_nh4 * c_nh4 + 2.0d0 * c_nh4 / unitconv) &
                  / (0.25d0 * c_nh4 + 1.0d0/unitconv) &
                  / (0.25d0 * c_nh4 + 1.0d0/unitconv) * ac_nh4
      else
        drate_nitri = temp_real
      endif 
 
      drate_nitri = drate_nitri * f_nh4 + rate_nitri * d_nh4

      Jacobian(ires_nh4,ires_nh4) = Jacobian(ires_nh4,ires_nh4) + drate_nitri

      Jacobian(ires_no3,ires_nh4) = Jacobian(ires_no3,ires_nh4) - drate_nitri

      JacobianDemand_nh4(ires_nh4,ires_nh4) = &
        JacobianDemand_nh4(ires_nh4,ires_nh4) + drate_nitri

      JacobianDemand_nh4(ires_no3,ires_nh4) = &
        JacobianDemand_nh4(ires_no3,ires_nh4) - drate_nitri

      Jacobian_nh4_to_no3(ires_nh4) = Jacobian_nh4_to_no3(ires_nh4) + drate_nitri 
    endif
  endif

  ! N2O production from nitr (Parton et al. 1996)
  if (this%ispec_n2o > 0) then

    rho_b = 1.25d0

    if (this%is_NH4_aqueous) then   
      ! mole/L * 1000 L/m3 * g/mol / kg/m3 = g/kg = mg/g = 1000 ug/g  
      M_2_ug_per_g  = theta *1000.0d0 * N_molecular_weight / rho_b * 1000.0d0
      !c_nh4_ugg = (c_nh4 + s_nh4 / theta / 1000.0d0)* M_2_ug_per_g
      c_nh4_ugg = c_nh4 * M_2_ug_per_g

      c_nh4_0 = this%c_nh4_ugg_0 / M_2_ug_per_g
      c_nh4_1 = this%c_nh4_ugg_1 / M_2_ug_per_g
    else
      ! mole/m3 * g/mol / kg/m3 = g/kg = mg/g = 1000 ug/g  
      mol_m3_2_ug_per_g  = N_molecular_weight / rho_b * 1000.0d0
      !c_nh4_ugg = (c_nh4 + s_nh4 / theta / 1000.0d0)* M_2_ug_per_g
      c_nh4_ugg = c_nh4 * mol_m3_2_ug_per_g

      c_nh4_0 = this%c_nh4_ugg_0 / mol_m3_2_ug_per_g
      c_nh4_1 = this%c_nh4_ugg_1 / mol_m3_2_ug_per_g

    endif
  
    if (c_nh4 <= c_nh4_0) then
      f_n2o = 0.0d0
      d_n2o = 0.0d0
    elseif (c_nh4 >= c_nh4_1 .or. &
      this%c_nh4_ugg_0 - this%c_nh4_ugg_1 > 1.0d-20) then
      f_n2o = 1.0d0
      d_n2o = 0.0d0
    else
      xxx = c_nh4 - c_nh4_0
      delta = c_nh4_1 - c_nh4_0
      f_n2o = 1.0d0 - (1.0d0 - xxx * xxx / delta / delta) ** 2
      d_n2o = 4.0d0 * (1.0d0 - xxx * xxx / delta / delta) * xxx / delta / delta
    endif  

    ! temperature response function (Parton et al. 1996)
    f_t = -0.06d0 + 0.13d0 * exp( 0.07d0 * tc )

    f_w = ((1.27d0 - saturation)/0.67d0)**(3.1777d0) * &
        ((saturation - 0.0012d0)/0.5988d0)**2.84d0

    ph = 6.5d0       ! default
    if (this%ispec_proton > 0) then
      if (reaction%species_idx%h_ion_id > 0) then
        ph = &
          -log10(rt_auxvar%pri_molal(reaction%species_idx%h_ion_id)* &
                 rt_auxvar%pri_act_coef(reaction%species_idx%h_ion_id))
      else if (reaction%species_idx%h_ion_id < 0) then
        ph = &
          -log10(rt_auxvar%sec_molal(abs(reaction%species_idx%h_ion_id))* &
                 rt_auxvar%sec_act_coef(abs(reaction%species_idx%h_ion_id)))
      endif
    endif
    f_ph = 0.56 + atan(rpi * 0.45 * (-5.0 + ph))/rpi

    if (f_t > 0.0d0 .and. f_w > 0.0d0 .and. f_ph > 0.0d0) then
       if (f_w > 1.0d0) then
         f_w = 1.0d0
       endif

       if (f_ph > 1.0d0) then
         f_ph = 1.0d0
       endif
    
      if (this%is_NH4_aqueous) then   
        temp_real = f_t * f_w * f_ph * this%k_nitr_n2o * kg_water
      else
        temp_real = f_t * f_w * f_ph * this%k_nitr_n2o * volume
      endif


       rate_n2o = 1.0 - exp(-0.0105d0 * c_nh4_ugg)  ! need to change units 
       ! Parton et al. 1996 unit is g N ha^-1 d^-1
       rate_n2o = rate_n2o * temp_real 
       rate_n2o = rate_n2o * f_nh4 * f_n2o
    
       Residual(ires_nh4) = Residual(ires_nh4) + rate_n2o
       Residual(ires_n2o) = Residual(ires_n2o) - 0.5d0 * rate_n2o
       
       if (this%ispec_ngasnit > 0) then
         Residual(ires_ngasnit) = Residual(ires_ngasnit) - 0.5d0 * rate_n2o
       endif

       RateDemand_nh4(ires_nh4) = RateDemand_nh4(ires_nh4) + rate_n2o
       RateDemand_nh4(ires_n2o) = RateDemand_nh4(ires_n2o) - 0.5d0 * rate_n2o
       
       if (this%ispec_ngasnit > 0) then
         RateDemand_nh4(ires_ngasnit) = RateDemand_nh4(ires_ngasnit) &
                                      - 0.5d0 * rate_n2o
       endif

       if (compute_derivative) then
         if (this%is_NH4_aqueous) then   
           drate_n2o = 0.0105d0*exp(-0.0105d0*c_nh4_ugg) &
                     * M_2_ug_per_g
         else
           drate_n2o = 0.0105d0*exp(-0.0105d0*c_nh4_ugg) &
                     * mol_m3_2_ug_per_g
         endif

         drate_n2o = drate_n2o * temp_real
 
         drate_n2o = drate_n2o * f_nh4 + rate_n2o * d_nh4 
         drate_n2o = drate_n2o * f_n2o + rate_n2o * d_n2o 

         Jacobian(ires_nh4,ires_nh4) = Jacobian(ires_nh4,ires_nh4) + drate_n2o

         Jacobian(ires_n2o,ires_nh4) = Jacobian(ires_n2o,ires_nh4) - &
           0.5d0 * drate_n2o
      
         if (this%ispec_ngasnit > 0 .and. (.not.this%bskipnitrjacobian)) then
           Jacobian(ires_ngasnit,ires_nh4)=Jacobian(ires_ngasnit,ires_nh4) - &
             0.5d0 * drate_n2o
         endif

         JacobianDemand_nh4(ires_nh4,ires_nh4) = &
           JacobianDemand_nh4(ires_nh4,ires_nh4) + drate_n2o

         JacobianDemand_nh4(ires_n2o,ires_nh4) = &
           JacobianDemand_nh4(ires_n2o,ires_nh4) - 0.5d0 * drate_n2o
      
         if (this%ispec_ngasnit > 0 .and. (.not.this%bskipnitrjacobian)) then
           JacobianDemand_nh4(ires_ngasnit,ires_nh4) = &
             JacobianDemand_nh4(ires_ngasnit,ires_nh4) - 0.5d0 * drate_n2o
         endif

        if (this%bdebugoutput) then
          write(*, *) 'Nitri: N2O', rate_n2o, drate_n2o
        endif

       endif
     endif
  endif

  ! N2O production from nitr (Parton et al. 1996), simplify 1 - e^(-x) to x
  if (this%ispec_n2o > 0 .and. (.not.this%bParton)) then
    ! temperature response function (Parton et al. 1996)
    f_t = -0.06d0 + 0.13d0 * exp( 0.07d0 * tc )

    f_w = ((1.27d0 - saturation)/0.67d0)**(3.1777d0) * &
        ((saturation - 0.0012d0)/0.5988d0)**2.84d0

    f_ph = 1.0d0 ! not ready yet,  0.56 + atan(rpi * 0.45 * (-5.0 + ph))/rpi

    if (f_t > 0.0d0 .and. f_w > 0.0d0 .and. f_ph > 0.0d0) then
      if (f_w > 1.0d0) then
        f_w = 1.0d0
      endif

      if (this%is_NH4_aqueous) then   
        temp_real = f_t * f_w * f_ph * this%k_nitr_n2o * kg_water
      else
        temp_real = f_t * f_w * f_ph * this%k_nitr_n2o * volume
      endif

      rate_n2o = temp_real * c_nh4 * f_nh4
    
      Residual(ires_nh4) = Residual(ires_nh4) + rate_n2o
      Residual(ires_n2o) = Residual(ires_n2o) - 0.5d0 * rate_n2o
       
      if (this%ispec_ngasnit > 0) then
        Residual(ires_ngasnit) = Residual(ires_ngasnit) - 0.5d0 * rate_n2o
      endif

      RateDemand_nh4(ires_nh4) = RateDemand_nh4(ires_nh4) + rate_n2o
      RateDemand_nh4(ires_n2o) = RateDemand_nh4(ires_n2o) - 0.5d0 * rate_n2o
       
      if (this%ispec_ngasnit > 0) then
        RateDemand_nh4(ires_ngasnit) = RateDemand_nh4(ires_ngasnit) &
                                     - 0.5d0 * rate_n2o
      endif

      if (compute_derivative) then
        drate_n2o = temp_real * f_nh4 + temp_real * c_nh4 * d_nh4 

        Jacobian(ires_nh4,ires_nh4) = Jacobian(ires_nh4,ires_nh4) + drate_n2o

        Jacobian(ires_n2o,ires_nh4) = Jacobian(ires_n2o,ires_nh4) - &
          0.5d0 * drate_n2o
      
        if (this%ispec_ngasnit > 0 .and. (.not.this%bskipnitrjacobian)) then
          Jacobian(ires_ngasnit,ires_nh4)=Jacobian(ires_ngasnit,ires_nh4) - &
            0.5d0 * drate_n2o
        endif

        JacobianDemand_nh4(ires_nh4,ires_nh4) = &
          JacobianDemand_nh4(ires_nh4,ires_nh4) + drate_n2o

        JacobianDemand_nh4(ires_n2o,ires_nh4) = &
          JacobianDemand_nh4(ires_n2o,ires_nh4) - 0.5d0 * drate_n2o
      
        if (this%ispec_ngasnit > 0 .and. (.not.this%bskipnitrjacobian)) then
          JacobianDemand_nh4(ires_ngasnit,ires_nh4) = &
            JacobianDemand_nh4(ires_ngasnit,ires_nh4) - 0.5d0 * drate_n2o
        endif

        if (this%bdebugoutput) then
          write(*, *) 'Nitri: N2O', rate_n2o, drate_n2o
        endif

       endif
     endif
  endif

end subroutine NitrReact

! ************************************************************************** !
!
! NitrDestroy: Destroys allocatable or pointer objects created in this 
!                  module
!
! ************************************************************************** !
subroutine NitrDestroy(this)

  implicit none
  
  class(clm_rxn_nitr_type) :: this  

end subroutine NitrDestroy

end module CLM_Rxn_Nitr_class


module CLM_Rxn_Deni_class

! ------------------------------------------------------------------------------
! Description
! denitrification function following Dickinson et al. 2002
! NO3-   -> 0.5 N2
! rate   = kmax ftheta fT NO3-
! fT     = exp(0.08(T - 298))
! ftheta = [(s - smin)/(1 - smin)]^b  smin = 0.6
! kmax   = 2.5e-5
! by t6g 10/06/2014 
! ------------------------------------------------------------------------------

  use CLM_Rxn_Base_class
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  use PFLOTRAN_Constants_module
  
  implicit none
  
  private
  
#include "petsc/finclude/petscsys.h"

  PetscInt, parameter :: TEMPERATURE_RESPONSE_FUNCTION_CLM4 = 1
  PetscInt, parameter :: TEMPERATURE_RESPONSE_FUNCTION_Q10 = 2

  type, public, &
    extends(clm_rxn_base_type) :: clm_rxn_deni_type
    PetscInt :: ispec_no3
    PetscInt :: ispec_n2
    PetscInt :: ispec_ngasdeni
    PetscInt :: temperature_response_function
    PetscReal :: Q10
    PetscReal :: k_deni_max                 ! deni rate
    PetscReal :: half_saturation
    PetscReal :: cutoff_no3_0  ! shut off
    PetscReal :: cutoff_no3_1  ! start to decrease from 1
    PetscReal :: residual_conc
    PetscBool :: bdebugoutput
    PetscBool :: is_NH4_aqueous
    PetscBool :: is_NO3_aqueous
    PetscBool :: bskipdenijacobian

  contains
    procedure, public :: ReadInput => DeniRead
    procedure, public :: Setup => DeniSetup
    procedure, public :: Evaluate => DeniReact
    procedure, public :: Destroy => DeniDestroy
  end type clm_rxn_deni_type

  public :: DeniCreate

contains

! ************************************************************************** !
!
! DeniCreate: Allocates deni reaction object.
!
! ************************************************************************** !
function DeniCreate()

  implicit none
  
  class(clm_rxn_deni_type), pointer :: DeniCreate

  allocate(DeniCreate)
  DeniCreate%ispec_no3 = 0
  DeniCreate%ispec_n2 = 0
  DeniCreate%ispec_ngasdeni = 0
  DeniCreate%temperature_response_function = TEMPERATURE_RESPONSE_FUNCTION_CLM4
  DeniCreate%Q10 = 1.5d0
  DeniCreate%k_deni_max = 2.5d-6  ! deni rate
  DeniCreate%half_saturation =  -1.0d-6
  DeniCreate%cutoff_no3_0 =-1.0d-20 
  DeniCreate%cutoff_no3_1 = 1.0d-20
  DeniCreate%residual_conc = 1.0d-10
  DeniCreate%bdebugoutput = PETSC_FALSE
  DeniCreate%is_NH4_aqueous = PETSC_TRUE
  DeniCreate%is_NO3_aqueous = PETSC_TRUE
  DeniCreate%bskipdenijacobian = PETSC_FALSE

  nullify(DeniCreate%next)  
      
end function DeniCreate

! ************************************************************************** !
!
! DeniRead: Reads input deck for deni reaction parameters (if any)
!
! ************************************************************************** !
subroutine DeniRead(this,input,option)

  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module, only : UnitsConvertToInternal
  
  implicit none
  
  class(clm_rxn_deni_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  PetscInt :: i
  character(len=MAXWORDLENGTH) :: word, internal_units
  
  do 
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,CLM_RXN,DENITRIFICATION')
    call StringToUpper(word)   

    select case(trim(word))
      case('TEMPERATURE_RESPONSE_FUNCTION')
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit

          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword', &
            'CHEMISTRY,CLM_RXN,DENITRIFICATION,TEMPERATURE RESPONSE FUNCTION')
          call StringToUpper(word)   

          select case(trim(word))
            case('CLM4')
              this%temperature_response_function = &
                TEMPERATURE_RESPONSE_FUNCTION_CLM4
            case('Q10')
              this%temperature_response_function = &
                TEMPERATURE_RESPONSE_FUNCTION_Q10
              call InputReadDouble(input,option,this%Q10)  
              call InputErrorMsg(input,option,'Q10', &
                'CHEMISTRY,CLM_RXN,DENITRI,TEMPERATURE RESPONSE FUNCTION')
            case default
              call InputKeywordUnrecognized(word, &
                'CHEMISTRY,CLM_RXN,DENITRIFICATION,TEMPERATURE RESPONSE FUNCTION', &
                option)
          end select
        enddo 

      case('RATE_CONSTANT')
        call InputReadDouble(input,option,this%k_deni_max)
        call InputErrorMsg(input,option,'k_deni_max', &
                 'CHEMISTRY,CLM_RXN,DENITRIFICATION,REACTION')
      case('RESIDUAL_NO3')
        call InputReadDouble(input,option,this%residual_conc)
        call InputErrorMsg(input,option,'residual_NO3', &
                  'CHEMISTRY,CLM_RXN,NITRIFICATION,REACTION')
      case('HALF_SATURATION_NO3')
        call InputReadDouble(input,option,this%half_saturation)
        call InputErrorMsg(input,option,'half saturation no3-', &
                 'CHEMISTRY,CLM_RXN,DENITRIFICATION,REACTION')
      case('CUTOFF_NO3')
        call InputReadDouble(input,option,this%cutoff_no3_0)
        call InputErrorMsg(input,option,'cutoff_no3_0', &
          'CHEMISTRY,CLM_RXN,DENITRIFICATION,REACTION')
        call InputReadDouble(input,option,this%cutoff_no3_1)
        call InputErrorMsg(input,option,'cutoff_no3_1', &
          'CHEMISTRY,CLM_RXN,DENITRIFICATION,REACTION')
        if (this%cutoff_no3_0 > this%cutoff_no3_1) then
          option%io_buffer = 'CHEMISTRY,CLM_RXN,DENITRIFICATION,' // &
            'NO3- cut off concentration > concentration ' // &
            'where down regulation function = 1.'
          call printErrMsg(option)
        endif
      case('DEBUG_OUTPUT')
        this%bdebugoutput = PETSC_TRUE
      case('JACOBIAN_DENI_SKIP')
        this%bskipdenijacobian = PETSC_TRUE
      case default
        call InputKeywordUnrecognized(word, &
               'CHEMISTRY,CLM_RXN,DENITRIFICATION,REACTION',option)
    end select
  enddo
  
end subroutine DeniRead

! ************************************************************************** !
!
! DeniSetup: Sets up the deni reaction either with parameters either
!                read from the input deck or hardwired.
!
! ************************************************************************** !
subroutine DeniSetup(this,reaction,option)

  use Reaction_Aux_module, only : reaction_type, GetPrimarySpeciesIDFromName
  use Option_module
  use Reaction_Immobile_Aux_module, only : GetImmobileSpeciesIDFromName 

  implicit none
  
  class(clm_rxn_deni_type) :: this
  type(reaction_type) :: reaction
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word
 
  word = 'NO3-'
  this%ispec_no3 = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)

  if (this%ispec_no3 < 0) then
    word = 'Nitrate'
    this%ispec_no3 = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)
    if (this%ispec_no3 > 0) then
      this%is_NO3_aqueous = PETSC_FALSE
    endif
  endif 

  if (this%ispec_no3 < 0) then
     option%io_buffer = 'CHEMISTRY,CLM_RXN,DENITRIFICATION: ' // &
                        ' NO3- or nitrate is not specified in the input file.'
     call printErrMsg(option)
  endif

  word = 'N2(aq)'
  this%ispec_n2 = GetPrimarySpeciesIDFromName(word,reaction, &
                        PETSC_FALSE,option)

  if (this%ispec_n2 < 0) then
     option%io_buffer = 'CHEMISTRY,CLM_RXN,DENITRIFICATION: ' // &
                        ' N2(aq) is not specified in the input file.'
     call printErrMsg(option)
  endif

  word = 'NGASdeni'
  this%ispec_ngasdeni = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)
 
end subroutine DeniSetup

! ************************************************************************** !
subroutine DeniReact(this,Residual,Jacobian,compute_derivative, &
                     rt_auxvar,global_auxvar,material_auxvar,reaction,option, &
                     RateDemand_nh4,RateSupply_nh4, &
                     JacobianDemand_nh4,JacobianSupply_nh4, &
                     RateDemand_no3,RateSupply_no3, &
                     JacobianDemand_no3,JacobianSupply_no3, &
                     Rate_nh4_to_no3,Jacobian_nh4_to_no3)

  use Option_module
  use Reaction_Aux_module
  use Material_Aux_class, only : material_auxvar_type
  use CLM_Rxn_Common_module, only: CalNLimitFunc

  implicit none

  class(clm_rxn_deni_type) :: this
  type(option_type) :: option
  type(reaction_type) :: reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar

  PetscBool :: compute_derivative
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  PetscReal :: RateDemand_nh4(reaction%ncomp)
  PetscReal :: RateSupply_nh4(reaction%ncomp)
  PetscReal :: RateDemand_no3(reaction%ncomp)
  PetscReal :: RateSupply_no3(reaction%ncomp)
  PetscReal :: JacobianDemand_nh4(reaction%ncomp,reaction%ncomp)
  PetscReal :: JacobianSupply_nh4(reaction%ncomp,reaction%ncomp)
  PetscReal :: JacobianDemand_no3(reaction%ncomp,reaction%ncomp)
  PetscReal :: JacobianSupply_no3(reaction%ncomp,reaction%ncomp)
  PetscReal :: Rate_nh4_to_no3
  PetscReal :: Jacobian_nh4_to_no3(reaction%ncomp)
  PetscReal :: porosity
  PetscReal :: volume
  PetscReal :: kg_water
  PetscInt :: local_id
  PetscErrorCode :: ierr

  PetscReal :: temp_real

  PetscInt :: ires_no3, ires_n2o, ires_n2
  PetscInt :: ires_ngasdeni

  PetscScalar, pointer :: bsw(:)
  PetscScalar, pointer :: bulkdensity(:)

  PetscReal :: s_min
  PetscReal :: tc
  PetscReal :: f_t, f_w

  PetscReal :: c_no3         ! mole/kg
  PetscReal :: ac_no3        ! mole/kg
  PetscReal :: f_no3         ! no3 / (half_saturation + no3)
  PetscReal :: d_no3         ! half_saturation/(no3 + half_saturation)^2
  PetscReal :: rate_deni, drate_deni
  PetscReal :: saturation
  PetscInt, parameter :: iphase = 1

  porosity = material_auxvar%porosity
  volume = material_auxvar%volume

  kg_water = material_auxvar%porosity*global_auxvar%sat(iphase)* &
            material_auxvar%volume*global_auxvar%den_kg(iphase)

  ! indices for C and N species
  ires_no3 = this%ispec_no3
  ires_n2 = this%ispec_n2
  ires_ngasdeni = this%ispec_ngasdeni + reaction%offset_immobile

! denitrification (Dickinson et al. 2002)
  if (this%ispec_n2 < 0) return

  temp_real = 1.0d0

  tc = global_auxvar%temp
  ! f_t = exp(0.08d0 * (tc - 25.d0))
  ! make it consistent with CLM CNNDynamicsMod.F90 line 652
  f_t = exp(0.08d0 * (tc + 273.15d0 - 298.0d0))

  saturation = global_auxvar%sat(1)
  s_min = 0.6d0
  f_w = 0.d0
  if (saturation > s_min) then
     f_w = (saturation - s_min)/(1.0d0 - s_min)
     f_w = f_w ** temp_real
  endif

  if (this%is_NO3_aqueous) then   
    c_no3     = rt_auxvar%pri_molal(this%ispec_no3)
    ac_no3    = rt_auxvar%pri_act_coef(this%ispec_no3)
    ires_no3 = this%ispec_no3
  else
    c_no3    = rt_auxvar%immobile(this%ispec_no3)
    ac_no3   = 1.0d0
    ires_no3 = this%ispec_no3 + reaction%offset_immobile
  endif

  call CalNLimitFunc(c_no3, ac_no3, this%residual_conc, this%half_saturation, &
                     this%cutoff_no3_0, this%cutoff_no3_1, f_no3, d_no3)

  ! add first order rate
  d_no3 = (c_no3 - this%residual_conc) * ac_no3 * d_no3 + ac_no3 * f_no3
  f_no3 = (c_no3 - this%residual_conc) * ac_no3 * f_no3
 
  if (f_t > 0.d0 .and. f_w > 0.d0) then
    if (this%is_NO3_aqueous) then   
      rate_deni = this%k_deni_max * f_t * f_w * kg_water * f_no3
    else
      rate_deni = this%k_deni_max * f_t * f_w * volume * f_no3
    endif

    Residual(ires_no3) = Residual(ires_no3) + rate_deni
    Residual(ires_n2) = Residual(ires_n2) - 0.5d0 * rate_deni
    
    if (this%ispec_ngasdeni > 0) then
      Residual(ires_ngasdeni) = Residual(ires_ngasdeni) - 0.5d0 * rate_deni
    endif

    RateDemand_no3(ires_no3) = RateDemand_no3(ires_no3) + rate_deni
    RateDemand_no3(ires_n2) = RateDemand_no3(ires_n2) - 0.5d0 * rate_deni
    
    if (this%ispec_ngasdeni > 0) then
      RateDemand_no3(ires_ngasdeni) = RateDemand_no3(ires_ngasdeni) &
                                    - 0.5d0 * rate_deni
    endif

    if (compute_derivative) then

      if (this%is_NO3_aqueous) then   
        drate_deni = this%k_deni_max * f_t * f_w * kg_water * d_no3 
      else
        drate_deni = this%k_deni_max * f_t * f_w * volume * d_no3 
      endif

      Jacobian(ires_no3,ires_no3) = Jacobian(ires_no3,ires_no3) + drate_deni

      Jacobian(ires_n2,ires_no3)=Jacobian(ires_n2,ires_no3) - 0.5d0*drate_deni
    
      if (this%ispec_ngasdeni > 0 .and. (.not.this%bskipdenijacobian)) then
        Jacobian(ires_ngasdeni,ires_no3) = Jacobian(ires_ngasdeni,ires_no3) &
                                         - 0.5d0 * drate_deni
      endif

      JacobianDemand_no3(ires_no3,ires_no3) = &
        JacobianDemand_no3(ires_no3,ires_no3) + drate_deni

      JacobianDemand_no3(ires_n2,ires_no3) = &
        JacobianDemand_no3(ires_n2,ires_no3) - 0.5d0 * drate_deni
    
      if (this%ispec_ngasdeni > 0 .and. (.not.this%bskipdenijacobian)) then
        JacobianDemand_no3(ires_ngasdeni,ires_no3) = &
          JacobianDemand_no3(ires_ngasdeni,ires_no3) - 0.5d0 * drate_deni
      endif

      if (this%bdebugoutput) then
        write(*, *) 'Deni:', rate_deni, drate_deni
      endif

    endif
  endif

end subroutine DeniReact

! ************************************************************************** !
!
! DeniDestroy: Destroys allocatable or pointer objects created in this 
!                  module
!
! ************************************************************************** !
subroutine DeniDestroy(this)

  implicit none
  
  class(clm_rxn_deni_type) :: this  

end subroutine DeniDestroy

end module CLM_Rxn_Deni_class

module CLM_Rxn_module

  ! extended from reaction_sandbox to implement demand based down regulation
  ! in RCLMRxn t6g 10/06/2014 

  use CLM_Rxn_Base_class
  use CLM_Rxn_Decomp_class
  use CLM_Rxn_Deni_class
  use CLM_Rxn_Nitr_class
  use CLM_Rxn_PlantN_class
  
  use PFLOTRAN_Constants_module

  implicit none
  
  private
  
#include "petsc/finclude/petscsys.h"

  class(clm_rxn_base_type), pointer, public :: clmrxn_list

  PetscBool :: bdownreg
  PetscBool :: bdebugoutput
  PetscBool :: b_ignore_production
  PetscReal :: residual_nh4
  PetscReal :: residual_no3
  PetscReal :: accelerator
  PetscReal :: cutoff_nh4_0
  PetscReal :: cutoff_nh4_1
  PetscReal :: cutoff_no3_0
  PetscReal :: cutoff_no3_1

  interface RCLMRxnRead
    module procedure RCLMRxnRead1
    module procedure RCLMRxnRead2
  end interface
  
  interface RCLMRxnDestroy
    module procedure RCLMRxnDestroy1
    module procedure RCLMRxnDestroy2
  end interface
  
  public :: RCLMRxnInit, &
            RCLMRxnRead, &
            RCLMRxnSkipInput, &
            RCLMRxnSetup, &
            RCLMRxn, &
            RCLMRxnDestroy

contains

! ************************************************************************** !

subroutine RCLMRxnInit(option)
  ! 
  ! Initializes the clmrxn list
  ! 
  use Option_module
  implicit none
  type(option_type) :: option

  if (associated(clmrxn_list)) then
    call RCLMRxnDestroy()
  endif
  nullify(clmrxn_list)
  bdownreg = PETSC_FALSE
  bdebugoutput = PETSC_FALSE
  b_ignore_production = PETSC_FALSE
  
  residual_nh4 =  1.0d-20
  residual_no3 =  1.0d-20
  accelerator  =  1.0d0
  cutoff_nh4_0 = -1.0d-18
  cutoff_nh4_1 =  1.0d-18
  cutoff_no3_0 = -1.0d-15
  cutoff_no3_1 =  1.0d-15

end subroutine RCLMRxnInit

! ************************************************************************** !

subroutine RCLMRxnSetup(reaction,option)
  ! 
  ! Calls all the initialization routines for all reactions in
  ! the clmrxn list
  ! 

  use Option_module
  use Reaction_Aux_module, only : reaction_type 
  
  implicit none
  
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  class(clm_rxn_base_type), pointer :: cur_clmrxn  

  character(len=MAXWORDLENGTH) :: word

  ! clmrxn reactions
  cur_clmrxn => clmrxn_list
  do
    if (.not.associated(cur_clmrxn)) exit
    call cur_clmrxn%Setup(reaction,option)
    cur_clmrxn => cur_clmrxn%next
  enddo 


end subroutine RCLMRxnSetup

! ************************************************************************** !

subroutine RCLMRxnRead1(input,option)
  ! 
  ! Reads input deck for reaction clmrxn parameters
  ! 

  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module
  
  implicit none
  
  type(input_type), pointer :: input
  type(option_type) :: option

  call RCLMRxnRead(clmrxn_list,input,option)

end subroutine RCLMRxnRead1

! ************************************************************************** !

subroutine RCLMRxnRead2(local_clmrxn_list,input,option)
  ! 
  ! RCLMRxnRead: Reads input deck for reaction clmrxn parameters
  ! 

  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module
  
  implicit none
  
  class(clm_rxn_base_type), pointer :: local_clmrxn_list  
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  class(clm_rxn_base_type), pointer :: new_clmrxn, cur_clmrxn
  
  nullify(new_clmrxn)
  do 
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','CHEMISTRY,CLM_RXN')
    call StringToUpper(word)   

    select case(trim(word))
      case('DECOMPOSITION')
        new_clmrxn => CLMDec_Create()
      case('DENITRIFICATION')
        new_clmrxn => DeniCreate()
      case('NITRIFICATION')
        new_clmrxn => NitrCreate()
      case('PLANTNTAKE')
        new_clmrxn => PlantNCreate()
      case('ENABLE_DOWNREGULATION')
        bdownreg = PETSC_TRUE
      case('DEBUG_OUTPUT')
        bdebugoutput = PETSC_TRUE
      case('IGNORE_PRODUCTION')
        b_ignore_production = PETSC_TRUE
      case('RESIDUAL_NH4')
        call InputReadDouble(input,option,residual_nh4)
        call InputErrorMsg(input,option,'residual nh4','CHEMISTRY,CLMRXN')
      case('RESIDUAL_NO3')
        call InputReadDouble(input,option,residual_no3)
        call InputErrorMsg(input,option,'residual no3','CHEMISTRY,CLMRXN')
      case('ACCELERATOR')
        call InputReadDouble(input,option,accelerator)
        call InputErrorMsg(input,option,'accelerator','CHEMISTRY,CLMRXN')
      case('CUTOFF_NH4')
        call InputReadDouble(input,option,cutoff_nh4_0)
        call InputErrorMsg(input,option,'cutoff_nh4_0','CHEMISTRY,CLM_RXN')
        call InputReadDouble(input,option,cutoff_nh4_1)
        call InputErrorMsg(input,option,'cutoff_nh4_1','CHEMISTRY,CLM_RXN')
        if (cutoff_nh4_0 > cutoff_nh4_1) then
          option%io_buffer = 'CHEMISTRY,CLM_RXN,' // &
            'NH4+ cut off concentration > concentration ' // &
            'where down regulation function = 1.'
          call printErrMsg(option)
        endif
      case('CUTOFF_NO3')
        call InputReadDouble(input,option,cutoff_no3_0)
        call InputErrorMsg(input,option,'cutoff_no3_0','CHEMISTRY,CLM_RXN')
        call InputReadDouble(input,option,cutoff_no3_1)
        call InputErrorMsg(input,option,'cutoff_no3_1','CHEMISTRY,CLM_RXN')
        if (cutoff_no3_0 > cutoff_no3_1) then
          option%io_buffer = 'CHEMISTRY,CLM_RXN,' // &
            'NO3- down regulation cut off concentration > concentration ' // &
            'where down regulation function = 1.'
          call printErrMsg(option)
        endif

      case default
        call InputKeywordUnrecognized(word,'CHEMISTRY,CLM_RXN',option)
    end select
    
    call new_clmrxn%ReadInput(input,option)
    
    if (.not.associated(local_clmrxn_list)) then
      local_clmrxn_list => new_clmrxn
    else
      cur_clmrxn => local_clmrxn_list
      do
        if (.not.associated(cur_clmrxn%next)) exit
        cur_clmrxn => cur_clmrxn%next
      enddo
      cur_clmrxn%next => new_clmrxn
    endif
  enddo
  
end subroutine RCLMRxnRead2

! ************************************************************************** !

subroutine RCLMRxnSkipInput(input,option)
  ! 
  ! Intelligently skips over CLM_RXN block
  ! 

  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module
  
  implicit none
  
  type(input_type), pointer :: input
  type(option_type) :: option
  
  class(clm_rxn_base_type), pointer :: dummy_list
  
  nullify(dummy_list)
  call RCLMRxnRead(dummy_list,input,option)
  call RCLMRxnDestroy(dummy_list)
  
end subroutine RCLMRxnSkipInput

! ************************************************************************** !

subroutine RCLMRxn(Residual,Jacobian,compute_derivative,rt_auxvar, &
                    global_auxvar,material_auxvar,reaction,option)
  ! 
  ! Evaluates reaction storing residual and/or Jacobian
  ! 

  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Reaction_Immobile_Aux_module
  use Material_Aux_class, only: material_auxvar_type
  
  implicit none

  type(option_type) :: option
  type(reaction_type) :: reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  class(clm_rxn_base_type), pointer :: cur_reaction
  PetscBool :: compute_derivative
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)

  PetscReal :: ResidualPre(reaction%ncomp)
  PetscReal :: Jacobianpre(reaction%ncomp,reaction%ncomp)

  PetscReal :: RateDemand_nh4(reaction%ncomp)
  PetscReal :: RateSupply_nh4(reaction%ncomp)
  PetscReal :: RateDemand_no3(reaction%ncomp)
  PetscReal :: RateSupply_no3(reaction%ncomp)

  PetscReal :: JacobianDemand_nh4(reaction%ncomp,reaction%ncomp)
  PetscReal :: JacobianSupply_nh4(reaction%ncomp,reaction%ncomp)
  PetscReal :: JacobianDemand_no3(reaction%ncomp,reaction%ncomp)
  PetscReal :: JacobianSupply_no3(reaction%ncomp,reaction%ncomp)

  PetscReal :: ddownscale_nh4(reaction%ncomp)
  PetscReal :: ddownscale_no3(reaction%ncomp)

  PetscReal :: Rate_nh4_to_no3
  PetscReal :: Jacobian_nh4_to_no3(reaction%ncomp)
  PetscReal :: f_supply

  PetscReal :: dt
  PetscReal :: kg_water_or_volume
  PetscReal :: c_nh4,avail_nh4,davail_nh4
  PetscReal :: c_no3,avail_no3,davail_no3
  PetscReal :: demand_nh4,supply_nh4,downscale_nh4
  PetscReal :: demand_no3,supply_no3,downscale_no3
  PetscReal :: regulator,dregulator,xxx,delta

  PetscBool :: b_nh4_downscaled
  PetscBool :: b_no3_downscaled
  PetscBool :: is_nh4_aqueous, is_no3_aqueous

  PetscInt, parameter :: iphase = 1
  PetscInt :: i,j 
  PetscInt :: ispec_nh4
  PetscInt :: ispec_no3
  PetscInt :: ires_nh4
  PetscInt :: ires_no3

  character(len=MAXWORDLENGTH) :: word

  ResidualPre = Residual
  JacobianPre = Jacobian 

  RateDemand_nh4      = 0.0d0
  RateSupply_nh4      = 0.0d0
  JacobianDemand_nh4  = 0.0d0
  JacobianSupply_nh4  = 0.0d0

  RateDemand_no3      = 0.0d0
  RateSupply_no3      = 0.0d0
  JacobianDemand_no3  = 0.0d0
  JacobianSupply_no3  = 0.0d0

  Rate_nh4_to_no3     = 0.0d0
  Jacobian_nh4_to_no3 = 0.0d0

  ddownscale_no3      = 0.0d0
  ddownscale_nh4      = 0.0d0

  cur_reaction => clmrxn_list
  do
    if (.not.associated(cur_reaction)) exit
    call cur_reaction%Evaluate(Residual,Jacobian,compute_derivative, &
                               rt_auxvar,global_auxvar,material_auxvar, &
                               reaction,option, &
                               RateDemand_nh4,RateSupply_nh4, &
                               JacobianDemand_nh4,JacobianSupply_nh4, &
                               RateDemand_no3,RateSupply_no3, &
                               JacobianDemand_no3,JacobianSupply_no3, &
                               Rate_nh4_to_no3, Jacobian_nh4_to_no3)
    cur_reaction => cur_reaction%next
  enddo

  if (.not.bdownreg) return

  ! down regulate sink if sink * dt > source * dt + conc 

  is_nh4_aqueous = PETSC_TRUE 
  word = 'NH4+'
  ispec_nh4 = GetPrimarySpeciesIDFromName(word,reaction,PETSC_FALSE,option)

  ires_nh4 = -999 
  if (ispec_nh4 < 0) then
    word = 'NH3(aq)'
    ispec_nh4 = GetPrimarySpeciesIDFromName(word,reaction,PETSC_FALSE, option)
  endif
 
  if (ispec_nh4 > 0) ires_nh4 = ispec_nh4
 
  if (ispec_nh4 < 0) then
    word = 'Ammonium'
    ispec_nh4 = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)
    if (ispec_nh4 > 0) then
      is_nh4_aqueous = PETSC_FALSE
      ires_nh4 = ispec_nh4 + reaction%offset_immobile 
    endif
  endif 

  if (ispec_nh4 < 0) then
    option%io_buffer = 'NH4+, NH3(aq) or Ammonium is specified in the input' // &
      'file for clm_rxn!'
    call printErrMsg(option)
  endif

  word = 'NO3-'
  ispec_no3 = GetPrimarySpeciesIDFromName(word,reaction,PETSC_FALSE,option)
 
  ires_no3 = -999 
  if (ispec_no3 > 0) ires_no3 = ispec_no3

  if (ispec_no3 < 0) then
    word = 'Nitrate'
    ispec_no3 = GetImmobileSpeciesIDFromName( &
            word,reaction%immobile,PETSC_FALSE,option)
    if (ispec_no3 > 0) then
      is_no3_aqueous = PETSC_FALSE
      ires_no3 = ispec_no3 + reaction%offset_immobile 
    endif
  endif 

  if (ispec_nh4 > 0 .and. ispec_no3 > 0) then
    if ((is_nh4_aqueous .and. (.not.is_no3_aqueous)) .or. & 
        ((.not.is_nh4_aqueous) .and. is_no3_aqueous)) then
      option%io_buffer = 'ERROR: Ammonium and nitrate have different phases: one in aqueous, the other in immobile,' // &
        'please use the same in the input file!'
      call printErrMsg(option)
    endif
  endif

  if (is_nh4_aqueous) then 
    kg_water_or_volume = material_auxvar%porosity*global_auxvar%sat(iphase)* &
               material_auxvar%volume*global_auxvar%den_kg(iphase)
  else
    kg_water_or_volume = material_auxvar%volume
  endif

  b_nh4_downscaled = PETSC_FALSE
  b_no3_downscaled = PETSC_FALSE
  dt = option%tran_dt

  if (b_ignore_production) then
    f_supply = 0.0d0
  else
    f_supply = 1.0d0
  endif

  ! if there is NH4+ demand
  if (RateDemand_nh4(ires_nh4) > 0.0d0) then
    ! following residual calculation sign, sink/demand is positive, 
    !                                      source/production is negative
    if (is_nh4_aqueous) then
      c_nh4 = rt_auxvar%pri_molal(ispec_nh4)
    else
      c_nh4 = rt_auxvar%immobile(ispec_nh4)
    endif

    if (cutoff_nh4_0 > 0.0d0) then
      if (c_nh4 <= cutoff_nh4_0) then
        regulator = 0.0d0
        dregulator = 0.0d0
      elseif (c_nh4 >= cutoff_nh4_1 .or. &
              cutoff_nh4_1 - cutoff_nh4_0 <= 1.0d-20) then
        regulator = 1.0d0
        dregulator = 0.0d0
      else
        xxx = c_nh4 - cutoff_nh4_0
        delta = cutoff_nh4_1 - cutoff_nh4_0
        regulator = 1.0d0 - (1.0d0 - xxx * xxx / delta / delta) ** 2
        dregulator = 4.0d0 * (1.0d0 - xxx * xxx / delta / delta) * xxx &
                   / delta / delta
      endif
    else
      regulator = 1.0d0
      dregulator = 0.0d0
    endif
   
    avail_nh4 = (c_nh4 - residual_nh4) * regulator
    davail_nh4 = regulator + (c_nh4 - residual_nh4) * dregulator

    demand_nh4 = RateDemand_nh4(ires_nh4) * dt

    supply_nh4 = RateSupply_nh4(ires_nh4) * dt * f_supply &
                 - avail_nh4 * kg_water_or_volume 
 
    ! if no supply, demand reactions will not occur
    if (supply_nh4 >= 0.0d0) then
      downscale_nh4 = 0.0d0
      Residual = ResidualPre + RateSupply_nh4
      if (compute_derivative) then
        ddownscale_nh4 = 0.0d0
        Jacobian = JacobianPre + JacobianSupply_nh4
      endif
      b_nh4_downscaled = PETSC_TRUE 

    elseif (demand_nh4 + supply_nh4 > 0.0d0) then
      ! if demand < supply
      b_nh4_downscaled = PETSC_TRUE
      downscale_nh4 = -1.0d0 * supply_nh4 / demand_nh4
      downscale_nh4 = downscale_nh4 * accelerator
 
      Residual = ResidualPre + RateSupply_nh4 + downscale_nh4 * RateDemand_nh4


      if (compute_derivative) then
 
        Jacobian = JacobianPre + JacobianSupply_nh4 &
                               + downscale_nh4 * JacobianDemand_nh4 

        do i = 1, reaction%ncomp
          if (i == ires_nh4) then          
            ddownscale_nh4(i) =-1.0d0 * ( &
              (JacobianSupply_nh4(ires_nh4,i) * dt * f_supply - &
              davail_nh4 * kg_water_or_volume) * demand_nh4 - &
              supply_nh4 * JacobianDemand_nh4(ires_nh4,i) * dt) / &
              demand_nh4 / demand_nh4   
          else
            ddownscale_nh4(i) =-1.0d0 * ( &
              JacobianSupply_nh4(ires_nh4,i) * dt * f_supply * demand_nh4 - & 
              supply_nh4 * JacobianDemand_nh4(ires_nh4,i) * dt) / &
              demand_nh4 / demand_nh4   
          endif
        enddo
        ddownscale_nh4 = ddownscale_nh4 * accelerator

        do i = 1, reaction%ncomp
          do j = 1, reaction%ncomp
            Jacobian(i,j) = Jacobian(i,j) + RateDemand_nh4(i) *ddownscale_nh4(j)
          enddo
        enddo

      endif
    else
      ! no down regulation
      downscale_nh4 = 1.0d0
      if (compute_derivative) then
        ddownscale_nh4 = 0.0d0
      endif
    endif

    if (bdebugoutput) then
      write(*, *) 'Cell id = ', option%iflag, &
                  'downscale_nh4 = ', downscale_nh4, &
                  'NH4+ = ', rt_auxvar%pri_molal(ires_nh4), &
                  'supply = ', supply_nh4, &
                  'demand = ', demand_nh4, &
                  'residual = ', Residual(ires_nh4)
      !write(*, *) 'residual = '
      !write(*, *) (Residual(i), i = 1, reaction%ncomp)
      if (compute_derivative) then
        !write(*, *) 'jacobian = '
        !do i = 1, reaction%ncomp
        !  write(*, *) (Jacobian(i, j), j = 1, reaction%ncomp)
        !enddo 
      endif
    endif
  else
    ! no demand, no down regulation
    downscale_nh4 = 1.0d0 
    if (compute_derivative) then
      ddownscale_nh4 = 0.0d0
    endif
  endif

  if (ires_no3 > 0) then
    ! if there is NO3- demand
    if (RateDemand_no3(ires_no3) > 0.0d0) then

      if (is_no3_aqueous) then
        c_no3 = rt_auxvar%pri_molal(ispec_no3)
      else
        c_no3 = rt_auxvar%immobile(ispec_no3)
      endif 

      if (cutoff_no3_0 > 0.0d0) then
        if (c_no3 <= cutoff_no3_0) then
          regulator = 0.0d0
          dregulator = 0.0d0
        elseif (c_no3 >= cutoff_no3_1 .or. &
                cutoff_no3_1 - cutoff_no3_0 <= 1.0d-20) then
          regulator = 1.0d0
          dregulator = 0.0d0
        else
          xxx = c_no3 - cutoff_no3_0
          delta = cutoff_no3_1 - cutoff_no3_0
          regulator = 1.0d0 - (1.0d0 - xxx * xxx / delta / delta) ** 2
          dregulator = 4.0d0 * (1.0d0 - xxx * xxx / delta / delta) * xxx &
                     / delta / delta
        endif
      else
        regulator = 1.0d0
        dregulator = 0.0d0
      endif
   
      avail_no3 = (c_no3 - residual_no3) * regulator
      davail_no3 = regulator + (c_no3 - residual_no3) * dregulator

      demand_no3 = RateDemand_no3(ires_no3) * dt

      supply_no3 = RateSupply_no3(ires_no3) * dt * f_supply - &
        avail_no3 * kg_water_or_volume - &
        Rate_nh4_to_no3 * downscale_nh4 * dt * f_supply

      if (supply_no3 >= 0.0d0) then
        ! if no NO3- supply, demanding reactions won't occur
        if (.not. b_nh4_downscaled) then
          Residual = ResidualPre + RateSupply_nh4 + RateDemand_nh4
          if (compute_derivative) then
            Jacobian = JacobianPre + JacobianSupply_nh4 + JacobianDemand_nh4 
            ddownscale_nh4 =0.0d0
          endif
        endif

        downscale_no3 = 0.0d0
        Residual = Residual + RateSupply_no3

        if (compute_derivative) then
          ddownscale_no3 = 0.0d0
          Jacobian = Jacobian + JacobianSupply_no3   
        endif 

        b_no3_downscaled = PETSC_TRUE
      elseif (demand_no3 + supply_no3 > 0.0d0) then
        ! if NO3- demand > supply
        ! if NH4+ is not limiting, accumulate rate and jacobian for NH4+
        if (.not. b_nh4_downscaled) then
          Residual = ResidualPre + RateSupply_nh4 + RateDemand_nh4
          if (compute_derivative) then
            Jacobian = JacobianPre + JacobianSupply_nh4 + JacobianDemand_nh4 
            ddownscale_nh4 =0.0d0
          endif
        endif

        b_no3_downscaled = PETSC_TRUE

        downscale_no3 = -1.0d0 * supply_no3 / demand_no3
        downscale_no3 = downscale_no3 * accelerator

        Residual = Residual + RateSupply_no3 + downscale_no3 * RateDemand_no3

        if (compute_derivative) then
          Jacobian = Jacobian + JacobianSupply_no3 &
                              + downscale_no3 * JacobianDemand_no3 

          do i = 1, reaction%ncomp
            if (i == ires_no3) then          
              ddownscale_no3(i) =-1.0d0 * ( &
                (JacobianSupply_no3(ires_no3,i) * dt * f_supply - &
                davail_no3* kg_water_or_volume - &
                Jacobian_nh4_to_no3(i) * downscale_nh4 * dt * f_supply  - &
                Rate_nh4_to_no3 * ddownscale_nh4(i) * dt * f_supply) * demand_no3 - &
                supply_no3 * JacobianDemand_no3(ires_no3,i)* dt ) / &
                demand_no3 / demand_no3   
            else
              ddownscale_no3(i) =-1.0d0 * ( &
                (JacobianSupply_no3(ires_no3,i) * dt * f_supply - &
                Jacobian_nh4_to_no3(i) * downscale_nh4 * dt * f_supply  - &
                Rate_nh4_to_no3 * ddownscale_nh4(i) * dt * f_supply) * &
                demand_no3 - &
                supply_no3 * JacobianDemand_no3(ires_no3,i) * dt ) / &
                demand_no3 / demand_no3   
            endif
          enddo

          ddownscale_no3 = ddownscale_no3 * accelerator

          do i = 1, reaction%ncomp
            do j = 1, reaction%ncomp
              Jacobian(i,j) = Jacobian(i,j) + RateDemand_no3(i) * ddownscale_no3(j)
            enddo
          enddo
        endif
      else
        ! no down regulation
        downscale_no3 = 1.0d0
        if (compute_derivative) then
          ddownscale_no3 = 0.0d0
        endif
      endif

      if (bdebugoutput) then
        write(*, *) 'Cell id = ', option%iflag, &
                    'downscale_no3 = ', downscale_no3, &
                    'NO3- = ', rt_auxvar%pri_molal(ires_no3), &
                    'supply = ', supply_no3, &
                    'demand = ', demand_no3, &
                    'residual = ', Residual(ires_no3)
        !write(*, *) 'residual = '
        !write(*, *) (Residual(i), i = 1, reaction%ncomp)
        if (compute_derivative) then
          !write(*, *) 'jacobian = '
          !do i = 1, reaction%ncomp
          !  write(*, *) (Jacobian(i, j), j = 1, reaction%ncomp)
          !enddo 
        endif
      endif
      
    endif
  endif 

  if (b_nh4_downscaled .and. (.not.b_no3_downscaled)) then
    Residual = Residual + RateSupply_no3 + RateDemand_no3

    if (compute_derivative) then
      Jacobian = Jacobian + JacobianSupply_no3 + JacobianDemand_no3 
    endif
  endif

end subroutine RCLMRxn

! ************************************************************************** !

subroutine RCLMRxnDestroy1()
  ! 
  ! Destroys master clmrxn list
  ! 

  implicit none

  call RCLMRxnDestroy(clmrxn_list)
  
end subroutine RCLMRxnDestroy1

! ************************************************************************** !

subroutine RCLMRxnDestroy2(local_clmrxn_list)
  ! 
  ! Destroys arbitrary clmrxn list
  ! 

  implicit none

  class(clm_rxn_base_type), pointer :: local_clmrxn_list

  class(clm_rxn_base_type), pointer :: cur_clmrxn, prev_clmrxn
  
  ! clmrxn reactions
  cur_clmrxn => local_clmrxn_list
  do
    if (.not.associated(cur_clmrxn)) exit
    prev_clmrxn => cur_clmrxn%next
    call cur_clmrxn%Destroy()
    deallocate(cur_clmrxn)
    cur_clmrxn => prev_clmrxn
  enddo  

end subroutine RCLMRxnDestroy2

end module CLM_Rxn_module

