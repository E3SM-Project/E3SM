module WellSpec_Base_class
#if WELL_CLASS

  use PFLOTRAN_Constants_module

  implicit none

  private


#include "petsc/finclude/petscsys.h"

  PetscInt, parameter, public :: WATER_PROD_WELL_TYPE = 1
  PetscInt, parameter, public :: OIL_PROD_WELL_TYPE = 2
  PetscInt, parameter, public :: GAS_PROD_WELL_TYPE = 3
  PetscInt, parameter, public :: WATER_INJ_WELL_TYPE = 4
  PetscInt, parameter, public :: OIL_INJ_WELL_TYPE = 5
  PetscInt, parameter, public :: GAS_INJ_WELL_TYPE = 6

  ! well control variables
  PetscInt, parameter, public :: CNTRL_VAR_BHP = 1
  PetscInt, parameter, public :: CNTRL_VAR_MASS_RATE = 2
  PetscInt, parameter, public :: CNTRL_VAR_VOL_RATE = 3
  PetscInt, parameter, public :: CNTRL_VAR_WGR = 4
  PetscInt, parameter, public :: CNTRL_VAR_GWR = 5

  ! limiting vars not yet used
  !PetscInt, parameter, public :: LMT_VAR_NO_LMT = -1
  PetscInt, parameter, public :: num_well_limits = 5
  PetscInt, parameter, public :: LMT_BHP_MAX = 1
  PetscInt, parameter, public :: LMT_BHP_MIN = 2
  PetscInt, parameter, public :: LMT_MASS_RATE_MAX = 3
  PetscInt, parameter, public :: LMT_VOL_RATE_MAX = 4
  PetscInt, parameter, public :: LMT_WGR_MAX = 5
  PetscInt, parameter, public :: LMT_GWR_MAX = 5

  PetscInt, parameter, public :: WELL_FACTOR_CONST = 1
  PetscInt, parameter, public :: WELL_FACTOR_PEACEMAN = 2

  PetscInt, parameter, public :: WELL_STATUS_OPEN = 1
  PetscInt, parameter, public :: WELL_STATUS_CLOSE = 2
  !PetscInt, parameter, public :: WELL_STATUS_AUTO = 3  ! connection computed during the simulation

  PetscInt, parameter, public :: CONN_STATUS_OPEN = 1
  PetscInt, parameter, public :: CONN_STATUS_CLOSE = 2

  type, public :: well_spec_base_type
      PetscInt :: id
      character(len=MAXWORDLENGTH) :: name       ! well_spec name
      PetscInt  :: itype                         ! well integer type
      character(len=MAXWORDLENGTH) :: ctype      ! well char type
      PetscInt :: well_fact_itype                ! type of well factor     
      PetscReal :: const_well_fact               ! [m^3] constant well factor 
      PetscReal, pointer :: dxyz_const(:)        ! [m] (dx1,dx2,dh) extensions of perforated cells, required for expl. unstr.
      PetscReal :: radius                        ! [m] well radius
      PetscReal :: skin_factor                   ! [-] well skin factor
      PetscReal :: theta_frac                    ! [-] portion of the well exposure angle (betwen 0 and 1)
      PetscBool :: input_const_drill_dir         ! if(input_const_drill_dir) all conns have const_drill_dir 
      PetscInt :: const_drill_dir                ! constant drilling direction 
      PetscInt :: status                         ! well status (can be open, closed or auto)
      PetscInt :: cntrl_var                      ! controlling variable (e.g. vol/mass rate)
      PetscInt :: num_limits                     ! number of limiting parameter
      PetscBool, pointer :: lmt_var(:)           ! limiting variables (e.g. pressure )   
      PetscBool :: input_z_pw_ref                ! true if the z_pw_ref has been read from input
      PetscReal :: input_val_z_pw_ref            ! eevation input for z_pw_ref 
      class(well_spec_base_type), pointer :: next ! points to next link in the list
  contains  ! add here type-bound procedure 
    procedure, public :: Init => WellSpecBaseInit
    procedure, public :: Read => WellSpecBaseRead
    procedure, public :: Clear => WellSpecBaseClear
  end type well_spec_base_type

  type, public :: well_spec_list_type
    PetscInt :: num_well_specs
    class(well_spec_base_type), pointer :: first
    class(well_spec_base_type), pointer :: last
    class(well_spec_base_type), pointer :: array(:)
  end type well_spec_list_type

  public :: WellSpecBaseCreate, WellSpecInitList, WellSpecDestroyList, &
            WellSpecAddToList, WellSpecGetPtrFromList
            ! WellSpecBaseDestroy !, WellSpecBaseInit

contains

!*********************************************************************!
! well-spec constructors
function WellSpecBaseCreate()
  !  
  ! Allocate and initialise well_spec 
  ! Author: Paolo Orsini (OpenGoSim)  
  ! Date : 5/17/2016
  !  
  implicit none

  class(well_spec_base_type), pointer :: WellSpecBaseCreate
  class(well_spec_base_type), pointer :: well_spec

  allocate(well_spec);
  call well_spec%Init(); 
 
  WellSpecBaseCreate => well_spec

end function WellSpecBaseCreate

!*********************************************************************!
subroutine WellSpecBaseInit(this)
  !  
  ! Initialise well_spec  
  ! Author: Paolo Orsini (OpenGoSim)  
  ! Date : 5/17/2016
  !  
  implicit none

  class(well_spec_base_type) :: this

  this%id=0;
  this%name = "";
  this%itype = 0;
  this%ctype = "";
  this%well_fact_itype = WELL_FACTOR_PEACEMAN;
  this%const_well_fact = 0.0;               
  this%radius = 0.0;
  this%skin_factor = 0.0;
  this%theta_frac = 1.0; !as default the well spans 360 degrees
  this%input_const_drill_dir = PETSC_TRUE;
  this%const_drill_dir = Z_DIRECTION; 
  this%status = WELL_STATUS_OPEN;
  this%cntrl_var = CNTRL_VAR_BHP;
  this%num_limits = 0;
  allocate(this%lmt_var(num_well_limits))
  nullify(this%dxyz_const);
  this%lmt_var = PETSC_FALSE
  this%input_z_pw_ref = PETSC_FALSE 
  this%input_val_z_pw_ref = 0.0;
  nullify(this%next);

end subroutine WellSpecBaseInit

! ************************************************************************** !

subroutine WellSpecInitList(list)
  ! 
  ! Initializes a well_specs list
  ! 
  ! Author: Paolo Orsini - OpenGoSim
  ! Date: 5/17/2016
  ! 

  implicit none

  class(well_spec_list_type) :: list  

  nullify(list%first)
  nullify(list%last)
  nullify(list%array)
  list%num_well_specs = 0

end subroutine WellSpecInitList

! ************************************************************************** !
subroutine WellSpecBaseRead(this,input,option)
  ! 
  ! Reads a well:spec from the input file
  ! 
  ! Author: Paolo Orsini
  ! Date: 5/17/2016
  ! 

  use Input_Aux_module
  use String_module
  use Option_module
  use Units_module
  
  implicit none
  
  type(option_type) :: option
  class(well_spec_base_type) :: this
  type(input_type), pointer :: input
  
  character(len=MAXWORDLENGTH) :: keyword, word, units, sub_keyword
  character(len=MAXWORDLENGTH) :: internal_units
  PetscInt :: i_dim

  internal_units = 'not_assigned'
 
  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','WELL_SPEC')
    call StringToUpper(keyword)   

    select case(trim(keyword))
      case('RADIUS')  ! this has a conversion factor to account for
        call InputReadDouble(input,option,this%radius)
        call InputErrorMsg(input,option,'well-radius','WELL_SPEC')
        call InputReadWord(input,option,word,PETSC_TRUE) 
        internal_units = 'meter'
        if (InputError(input)) then
          word = trim(keyword) // ' UNITS'
          call InputDefaultMsg(input,option,word)
        else
          units = trim(word)
          this%radius = UnitsConvertToInternal(units,internal_units,option) * &
                                 this%radius
        endif
      case('SKIN_FACTOR')
        call InputReadDouble(input,option,this%skin_factor)
        call InputErrorMsg(input,option,'skin factor','WELL_SPEC')
        call InputReadWord(input,option,word,PETSC_TRUE) 
      case('DX_DY_DZ_CONST')
        allocate(this%dxyz_const(3))
        this%dxyz_const = 0.0d0
        !internal_units = 'meter,meter'
        do i_dim =1,3
          call InputReadDouble(input,option,this%dxyz_const(i_dim))
          call InputErrorMsg(input,option,'DX1_DX2_DZ','WELL_SPEC')
        end do
        call InputReadWord(input,option,word,PETSC_TRUE)
        if (InputError(input)) then
          word = trim(keyword) // ' UNITS'
          call InputDefaultMsg(input,option,word)
        else
          do i_dim =1,3
            internal_units = 'meter'
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,keyword,'WELL_SPEC')
              units = trim(word)
              this%dxyz_const(i_dim) = &
                      UnitsConvertToInternal(units,internal_units,option) * &
                      this%dxyz_const(i_dim)
          end do
        endif
      case('THETA_FRACTION')
        call InputReadDouble(input,option,this%theta_frac)
        call InputErrorMsg(input,option,'theta angle fraction','WELL_SPEC')
        call InputReadWord(input,option,word,PETSC_TRUE)      
      case('WELL_TYPE')
        call InputReadWord(input,option,keyword,PETSC_TRUE)
        call InputErrorMsg(input,option,'keyword','WELL_TYPE')
        select case(trim(keyword))
          case('WATER_PRODUCER','WAT_PROD')
            this%itype = WATER_PROD_WELL_TYPE;
            this%ctype="water_producer"
          case('GAS_PRODUCER','GAS_PROD')
            this%itype = GAS_PROD_WELL_TYPE;
            this%ctype="gas_producer"   
          case('WATER_INJECTOR','WAT_INJ')
            this%itype = WATER_INJ_WELL_TYPE;
            this%ctype="water_injector"
          case('GAS_INJECTOR','GAS_INJ')
            this%itype = GAS_INJ_WELL_TYPE;
            this%ctype="gas_injector"
          case('OIL_PRODUCER','OIL_PROD')
            this%itype = OIL_PROD_WELL_TYPE;
            this%ctype="oil_producer"
          case default
            option%io_buffer = 'WELL_SPEC keyword: '//trim(keyword)//' not recognized'
            call printErrMsg(option)
        end select
      case('VAR_CONTROL')
        call InputReadWord(input,option,keyword,PETSC_TRUE)
        call InputErrorMsg(input,option,'keyword','VAR_CONTROL')
        select case(trim(keyword))
          case('BHP')
            this%cntrl_var = CNTRL_VAR_BHP;
          case('MASS_RATE') 
            this%cntrl_var = CNTRL_VAR_MASS_RATE;
          case('VOL_RATE')
            this%cntrl_var = CNTRL_VAR_VOL_RATE;
          case default
            option%io_buffer = 'WELL_SPEC keyword: '//trim(keyword)//' not recognized'
            call printErrMsg(option)
        end select 
      case('VAR_LIMIT')
       do
         call InputReadPflotranString(input,option)
         if (InputError(input)) exit
         if (InputCheckExit(input,option)) exit
         call InputReadWord(input,option,sub_keyword,PETSC_TRUE)
         call InputErrorMsg(input,option,'sub_keyword','WELL VARLIMIT')
         call StringToUpper(sub_keyword) 
        ! each limiting variable can be prescribed fixed bounds, or be managed by VFP
        ! VFP (vertical flow performance) can be defined using tables
         select case(sub_keyword)
            case("BHP_MIN")
              this%lmt_var(LMT_BHP_MIN) = PETSC_TRUE
            case("BHP_MAX")
              this%lmt_var(LMT_BHP_MAX) = PETSC_TRUE
            case("MASS_RATE_MAX")
              this%lmt_var(LMT_MASS_RATE_MAX) = PETSC_TRUE
            case("VOL_RATE_MAX") 
              this%lmt_var(LMT_VOL_RATE_MAX) = PETSC_TRUE
            case("WGR_MAX")
              this%lmt_var(LMT_WGR_MAX) = PETSC_TRUE
            case("GWR_MAX")
              this%lmt_var(LMT_GWR_MAX) = PETSC_TRUE
         end select
       end do
      ! case('WELL_INIT_STATUS')
        ! for testing purpose let's assume initial status is open 
      case('CONSTANT_WELL_FACTOR')
        this%well_fact_itype = WELL_FACTOR_CONST
        call InputReadDouble(input,option,this%const_well_fact)
        call InputErrorMsg(input,option,'well factor','WELL_SPEC')
        call InputReadWord(input,option,word,PETSC_TRUE) 
        internal_units = 'm^3'
        if (InputError(input)) then
          word = trim(keyword) // ' UNITS'
          call InputDefaultMsg(input,option,word)
        else
          units = trim(word)
          this%const_well_fact = &
              UnitsConvertToInternal(units,internal_units,option) * &
                   this%const_well_fact
        end if
      case('CONST_DRILL_DIR')
        this%input_const_drill_dir = PETSC_TRUE
        call InputReadWord(input,option,keyword,PETSC_TRUE)
        call InputErrorMsg(input,option,'keyword','CONST_DRILL_DIR')
        select case(trim(keyword))
          case('DIR_X')
            this%const_drill_dir = X_DIRECTION;
          case('DIR_Y') 
            this%const_drill_dir = Y_DIRECTION;
          case('DIR_Z')
            this%const_drill_dir = Z_DIRECTION;
          case default
            option%io_buffer = 'WELL_SPEC keyword: '//trim(keyword)//' not recognized'
            call printErrMsg(option)
        end select    
      case('Z_REF')
        this%input_z_pw_ref = PETSC_TRUE
        call InputReadDouble(input,option,this%input_val_z_pw_ref)
        call InputErrorMsg(input,option,'z_ref','WELL_SPEC')
        call InputReadWord(input,option,word,PETSC_TRUE)
        internal_units = 'meter' 
        if (InputError(input)) then
          word = trim(keyword) // ' UNITS'
          call InputDefaultMsg(input,option,word)
        else
          units = trim(word)
          this%radius = UnitsConvertToInternal(units,internal_units,option) * &
                                 this%input_val_z_pw_ref
        endif
      case default
        option%io_buffer = 'WELL_SPEC keyword: '//trim(keyword)//' not recognized'
        call printErrMsg(option)
    end select
  enddo

 
end subroutine WellSpecBaseRead

! *************************************************************************** !

subroutine WellSpecAddToList(new_well_spec,list)
  ! 
  ! Adds a new well_spec to a well_spec list
  ! 
  ! Author: Paolo Orsini - OpenGoSim
  ! Date: 5/17/2016
  ! 

  implicit none
  
  class(well_spec_base_type), pointer :: new_well_spec
  type(well_spec_list_type), pointer :: list
  
  list%num_well_specs = list%num_well_specs + 1
  new_well_spec%id = list%num_well_specs
  if (.not.associated(list%first)) list%first => new_well_spec
  if (associated(list%last)) list%last%next => new_well_spec
  list%last => new_well_spec
  
end subroutine WellSpecAddToList

! ************************************************************************** !

function WellSpecGetPtrFromList(well_spec_name,well_spec_list)
  ! 
  ! Returns a pointer to the well_spec matching well_spec_name
  ! 
  ! Author: Paolo Orsini - OpenGoSim
  ! Date: 5/17/2016
  ! 

  use String_module

  implicit none
  
  class(well_spec_base_type), pointer :: WellSpecGetPtrFromList
  character(len=MAXWORDLENGTH) :: well_spec_name
  PetscInt :: length
  type(well_spec_list_type) :: well_spec_list

  class(well_spec_base_type), pointer :: well_spec
    
  nullify(WellSpecGetPtrFromList)
  well_spec => well_spec_list%first
  
  do 
    if (.not.associated(well_spec)) exit
    length = len_trim(well_spec_name)
    if (length == len_trim(well_spec%name) .and. &
        StringCompare(well_spec%name,well_spec_name,length)) then
      WellSpecGetPtrFromList => well_spec
      return
    endif
    well_spec => well_spec%next
  enddo
  
end function WellSpecGetPtrFromList

! ************************************************************************** !

subroutine WellSpecDestroyList(well_spec_list)
  ! 
  ! Deallocates a list of well_specs
  ! 
  ! Author: Paolo Orsini
  ! Date: 5/17/2016
  ! 

  implicit none
  
  type(well_spec_list_type), pointer :: well_spec_list
  
  class(well_spec_base_type), pointer :: well_spec, prev_well_spec
  
  if (.not.associated(well_spec_list)) return
  
  well_spec => well_spec_list%first
  do 
    if (.not.associated(well_spec)) exit
    prev_well_spec => well_spec
    well_spec => well_spec%next
    call WellSpecBaseDestroy(prev_well_spec)
  enddo
  
  well_spec_list%num_well_specs = 0
  nullify(well_spec_list%first)
  nullify(well_spec_list%last)
  if (associated(well_spec_list%array)) deallocate(well_spec_list%array)
  nullify(well_spec_list%array)
  
  deallocate(well_spec_list)
  nullify(well_spec_list)

end subroutine WellSpecDestroyList

! ************************************************************************** !

! WellSpecDestroy
subroutine WellSpecBaseDestroy(well_spec)
  ! 
  ! Deallocates a well_spec
  ! 
  ! Author: Paolo Orsini - OpenGoSim
  ! Date: 5/17/2016
  ! 
  
  implicit none
  
  class(well_spec_base_type), pointer :: well_spec
  
  if (.not.associated(well_spec)) return

  !call well_spec%WellSpecBaseClear();
  call well_spec%Clear();

  deallocate(well_spec)
  nullify(well_spec)

end subroutine WellSpecBaseDestroy

! ************************************************************************** !

subroutine WellSpecBaseClear(this)
  ! 
  ! Clear well data 
  ! 
  ! Author: Paolo Orsini
  ! Date: 5/17/2016
  ! 
  use Utility_module, only : DeallocateArray
  
  implicit none
  
  class(well_spec_base_type) :: this

  call DeallocateArray(this%lmt_var)   
  call DeallocateArray(this%dxyz_const)

  nullify(this%next);

end subroutine WellSpecBaseClear

! ************************************************************************** !

#endif  
end module WellSpec_Base_class

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!end of WELL_CLASS
