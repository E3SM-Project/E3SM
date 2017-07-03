module Well_module
#ifdef WELL_CLASS

  use PFLOTRAN_Constants_module
  use WellSpec_Base_class
  use Well_Base_class
  use Well_Flow_class
  use Well_FlowEnergy_class
  use Well_WaterInjector_class
  use Well_OilProducer_class
  use Well_TOilIms_class
  !add here other well classes, e.g. Wells_XXXX_class

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  public :: CreateWell, WellAuxVarSetUp, WellOutput, WellDestroy

contains

! ************************************************************************** !

function CreateWell(well_spec,option)
  ! 
  ! Create a toil ims well object based on type specified in the well_spec
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 05/18/16
  ! 

  use Connection_module
  use Option_module

  implicit none

  class(well_spec_base_type), pointer :: well_spec
  type(option_type) :: option
  class(well_base_type), pointer :: CreateWell

  select case(option%iflowmode)
    case(TOIL_IMS_MODE)
      CreateWell => CreateTOilImsWell(well_spec,option)
      !add here create well for other flow modes
    case default
      option%io_buffer = 'Well model supported supported for TOIL_IMS only'
      call printErrMsg(option)
  end select

  !Debug printing 
#ifdef WELL_DEBUG
  write(*,*) "well_factor type = ", CreateWell%spec%well_fact_itype 
  write(*,*) "well type = ", CreateWell%spec%ctype
  write(*,*) "radius = ", CreateWell%spec%radius

  select type(CreateWell)
    class is(well_toil_ims_wat_inj_type)
      write(*,*) "temp", CreateWell%tw_ref
      write(*,*) "well_press", CreateWell%pw_ref
  end select 

  call CreateWell%PrintMsg(); 
#endif

  !Create well outfile and write its header 
  !not here - otherwise will attempt to create a file for each process
  !at this stage a well is created in each process - even if empty   


end function CreateWell

! ************************************************************************** !
subroutine WellAuxVarSetUp(well,connection_set,flow_condition,aux, &
                           cpl_idx_start,ss_flow_vol_fluxes,option)
  ! 
  ! Create a toil ims well object based on type specified in the well_spec
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 06/03/16
  ! 

  use Option_module
  use Auxiliary_module
  use Condition_module
  use Connection_module

  implicit none

  class(well_base_type), pointer :: well
  type(connection_set_type), pointer :: connection_set 
  type(flow_condition_type), pointer :: flow_condition
  type(auxiliary_type) :: aux  
  PetscInt, intent(in) :: cpl_idx_start
  PetscReal, pointer, intent(in) :: ss_flow_vol_fluxes(:,:) 
  type(option_type) :: option

  PetscInt :: cpl_idx_end

#ifdef WELL_DEBUG
  !write(*,"('WS d11 before = ',e10.4)"), aux%TOil_ims%auxvars(0,1)%den(1)
  !write(*,"('WS d12 before = ',e10.4)"), aux%TOil_ims%auxvars(0,1)%den(2) 
  !write(*,"('WS p11 before = ',e10.4)"), aux%TOil_ims%auxvars(0,1)%pres(1) 
  !write(*,"('WS t1 before = ',e10.4)"), aux%TOil_ims%auxvars(0,1)%temp 
  write(*,"('WS p011 before = ',e10.4)"), aux%TOil_ims%auxvars(0,1)%pres(1)
  write(*,"('WS p111 before = ',e10.4)"), aux%TOil_ims%auxvars(1,1)%pres(1)
  write(*,"('WS p211 before = ',e10.4)"), aux%TOil_ims%auxvars(2,1)%pres(1)
  write(*,"('WS p311 before = ',e10.4)"), aux%TOil_ims%auxvars(3,1)%pres(1)
#endif

  select type(well)
    !if only auxvar_flow_energy needed can use class is(well_flow_energy_type)
    class is(well_toil_ims_wat_inj_type)
      well%flow_energy_auxvars => aux%TOil_ims%auxvars   
      well%flow_auxvars => aux%TOil_ims%auxvars
    class is(well_toil_ims_oil_prod_type)
      well%flow_energy_auxvars => aux%TOil_ims%auxvars   
      well%flow_auxvars => aux%TOil_ims%auxvars
    !when well implmented for other flow modes - add below
  end select 

  !allocate global (entire)  well arrays - flow components
  select type(well)
    class is(well_flow_type)
      well%flow_condition => flow_condition

      nullify(well%well_conn_den_kg)
      allocate(well%well_conn_den_kg(well%well_num_conns))
      well%well_conn_den_kg = 0.0d0
      nullify(well%well_conn_h_sorted)
      allocate(well%well_conn_h_sorted(well%well_num_conns))
      well%well_conn_h_sorted = 0.0d0

      cpl_idx_end = cpl_idx_start + connection_set%num_connections - 1  
      well%ss_flow_vol_fluxes => &
                ss_flow_vol_fluxes(1:option%nphase,cpl_idx_start:cpl_idx_end)

  end select

  !allocate global (entire)  well arrays - energy components
  select type(well)
    class is(well_flow_energy_type)
      nullify(well%well_conn_temp)
      allocate(well%well_conn_temp(well%well_num_conns))
      well%well_conn_temp = 0.0d0
  end select

  !for well base
  well%connection_set => connection_set


end subroutine WellAuxVarSetUp

! ************************************************************************** !

!subroutine WellOutput(well,output_option,src_name,option)
subroutine WellOutput(well,output_option,option)
  ! 
  ! Handle the well output part common to all type of wells
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 06/07/16
  ! 

  use Option_module
  use Output_Aux_module

  implicit none

  class(well_base_type), pointer :: well
  type(output_option_type), pointer :: output_option
  !character(len=MAXWORDLENGTH) :: src_name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: wfile_name
  PetscMPIInt :: cur_w_myrank
  PetscInt :: ios, ierr

  if( well%connection_set%num_connections > 0 ) then
    call MPI_Comm_rank(well%comm, cur_w_myrank, ierr )  
    if(well%cntr_rank == cur_w_myrank ) then
      wfile_name = trim(option%global_prefix) // "_" // &
                        trim(well%name) // ".tec"
                        !trim(src_name) // ".tec" 
      open(unit=IUNIT_TEMP,file=wfile_name,action="write", &
           position="append",status="old",iostat=ios)

      call well%output(IUNIT_TEMP,output_option,option)

      close(IUNIT_TEMP)
    end if
  end if 

end subroutine WellOutput

! ************************************************************************** !

subroutine WellDestroy(well)
  ! 
  ! Destroy well
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 06/07/16
  ! 

  class(well_base_type), pointer :: well 

  select type(well)
    class is(well_flow_energy_type)
      call FlowEnergyWellStrip(well)
    !for now data/data pointer present in well_flow_energy_type only
    !to be reviewed if data/data pointer also in doughter classes
    !in such case should include last extension, e.g.
    ! class is(well_toil_ims_wat_inj_type), etc 
  end select

  deallocate(well) 

end subroutine WellDestroy

! ************************************************************************** !

#endif  
end module Well_module
!end of WELL_CLASS

