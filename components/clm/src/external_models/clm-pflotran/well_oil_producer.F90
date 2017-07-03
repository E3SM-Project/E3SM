module Well_OilProducer_class
#ifdef WELL_CLASS

  use PFLOTRAN_Constants_module
  use WellSpec_Base_class
  use Well_Base_class
  use Well_Flow_class
  use Well_FlowEnergy_class

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  type, public, extends(well_flow_energy_type) :: well_oil_producer_type
    ! .................
  contains  ! add here type-bound procedure 
    procedure, public  :: PrintOutputHeader => WellOilProdPrintOutputHeader
    procedure, public :: VarsExplUpdate => WellOilProdVarsExplUpdate
    procedure, public :: LimitCheck => WellOilProdLimitCheck
    procedure, public :: InitDensity => OilProdInitDensity
  end type  well_oil_producer_type

  !public :: CreateTOilImsWell

contains

! ************************************************************************** !

subroutine WellOilProdPrintOutputHeader(this,output_option,file_unit)
  ! 
  ! Write header for well_TOilIms output file
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 05/18/16
  ! 
  use Output_Aux_module

  implicit none

  class(well_oil_producer_type) :: this
  type(output_option_type), intent(in) :: output_option
  PetscInt, intent(in) :: file_unit

  !character(len=MAXWORDLENGTH) :: tunit

  !tunit = trim(output_option%tunit)

  write(*,*) "Well WellOilProdPrintOutputHeader to be extended"
  

end subroutine WellOilProdPrintOutputHeader

! ************************************************************************** !

subroutine OilProdInitDensity(this,grid,option)
  !
  ! DOES NOTHING - eneable called to well%InitDensity 
  !
  ! Include commented implementation of init well fluid density
  ! computed as phase density avergage weighted 
  ! by phase saturations. Quantity taken at the control grid block.
  ! At this point the phase fluxes from the perforated grid blocks 
  ! into the well are not yet initialised.
  !
  ! Author: Paolo Orsini (OpenGoSim)  
  ! Date : 6/12/2016
  !

  use Grid_module
  use Option_module
  use EOS_Water_module

  implicit none

  class(well_oil_producer_type) :: this
  type(grid_type), pointer :: grid  
  type(option_type) :: option

  !PetscMPIInt :: cur_w_myrank
  !PetscInt :: ghost_cntrl_id,ierr

  !PetscInt :: iphase

  !THIS ROUTINE DOE NOTHING BECAUSE THE WELL FLUID DENSITY IS INITIALISE
  !IN WellFlowConnDenUpdate, called by FlowHydroCorrUpdate
  ! MUST BE HERE TO AVOID FAILING CALL TO well%InitDensity  

  !BELOW AN ALTERNATIVE TO INITALISE THE WELL FLUID DENSITY USING THE WELL
  ! CONTROL GRID BLOCK 

  !call MPI_Comm_rank(this%comm, cur_w_myrank, ierr )  

  !if(this%cntr_rank == cur_w_myrank ) then
  !  ghost_cntrl_id = grid%nL2G(this%cntrl_lcell_id); 
  !
  !  !Avergae on generic number of phases 
  !  this%dw_kg_ref(option%oil_phase) = 0.0d0
  !  do iphase = 1,option%nphase
  !    this%dw_kg_ref(option%oil_phase) = this%dw_kg_ref(option%oil_phase) +       
  !            this%flow_auxvars(ZERO_INTEGER,ghost_cntrl_id)%den_kg(iphase) * &
  !            this%flow_auxvars(ZERO_INTEGER,ghost_cntrl_id)%sat(iphase)
  !  end do
  !  !END of average loop
  !  !assuming perfect mixture both phase have the same densities
  !  this%dw_kg_ref(option%liquid_phase) = this%dw_kg_ref(option%oil_phase)
  !
  !end if
  !
  !call MPI_Bcast ( this%dw_kg_ref(option%oil_phase),1, &
  !                 MPI_DOUBLE_PRECISION, this%cntr_rank, this%comm, ierr )


end subroutine OilProdInitDensity

! ************************************************************************** !

subroutine WellOilProdVarsExplUpdate(this,grid,option)
  !
  ! Explicit update of well variable for a water injector
  !
  ! Author: Paolo Orsini (OpenGoSim)  
  ! Date : 6/04/2016
  !

  use Grid_module
  use Option_module

  class(well_oil_producer_type) :: this
  type(grid_type), pointer :: grid
  !PetscReal :: ss_fluxes(:,:)
  type(option_type) :: option

  !PetscReal :: enth_src_h2o, dw_h2o_kg,dw_h2o_mol

  PetscInt :: ierr
  PetscInt :: i_ph

  if(this%connection_set%num_connections == 0 ) return

  select case(this%spec%cntrl_var) 
    case(CNTRL_VAR_BHP)
      this%pw_ref = this%flow_condition%flow_well%pressure%dataset%rarray(1)
      call this%QPhase(grid,option%oil_phase,option)      
      call this%MRPhase(grid,option%oil_phase,option)

    case(CNTRL_VAR_MASS_RATE)      
      !call this%PressRef(grid,option%liquid_phase,option)
      call this%PressRefMRProd(grid,option%oil_phase,option)
      this%mr_fld(option%oil_phase) = &
                  this%flow_condition%flow_well%rate%dataset%rarray(1)
      call this%QPhase(grid,option%oil_phase,option)
 
    case(CNTRL_VAR_VOL_RATE)
      this%q_fld(option%oil_phase) = &
                this%flow_condition%flow_well%rate%dataset%rarray(1)
      !call this%PressRef(grid,option%liquid_phase,option)
      call this%PressRefQ(grid,option%oil_phase,option)
      call this%MRPhase(grid,option%oil_phase,option)
  end select

  !compute the other phase rates
  do i_ph = 1,option%nphase 
    if ( i_ph /= option%oil_phase ) then
      call this%QPhase(grid,i_ph,option)      
      call this%MRPhase(grid,i_ph,option)      
    end if
  end do

  !update well fluid temperature as perforated grid block average 
  !call this%AverageTemp(grid,ss_fluxes,option)
  call this%AverageTemp(grid,option)

  !The well fluid density is treated as full explicit and it is updated
  !when computing the well connection densities for hydro corrections


end subroutine WellOilProdVarsExplUpdate

! ************************************************************************** !

subroutine WellOilProdLimitCheck(this,pass,option)
  !
  ! Perform limit check for a water injector
  !
  ! Author: Paolo Orsini (OpenGoSim)  
  ! Date : 6/06/2016
  !

  use Option_module

  implicit none

  class(well_oil_producer_type) :: this
  PetscBool :: pass
  type(option_type) :: option

  PetscInt :: cntrl_var_tmp 
  PetscReal :: rate_max
  PetscReal :: bhp_min

  pass = PETSC_TRUE
  !if no changes cntrl_var mantains its initial value
  cntrl_var_tmp = this%spec%cntrl_var

  select case(this%spec%cntrl_var)
    case(CNTRL_VAR_BHP)
      if(this%spec%lmt_var(LMT_MASS_RATE_MAX)) then
        rate_max = this%flow_condition%flow_well%rate%dataset%rarray(1)
        if( this%mr_fld(option%oil_phase) > rate_max ) then
          print *, "oil_producer control switch: " // &
                   "BHP -> MASS_RATE, mass_rate - rate_max = ",&
                    & this%mr_fld(option%oil_phase), rate_max
          this%mr_fld(option%oil_phase) = rate_max
          cntrl_var_tmp = CNTRL_VAR_MASS_RATE   
          pass = PETSC_FALSE   
          !add limit on BHP min
          this%spec%lmt_var(LMT_BHP_MIN) = PETSC_TRUE
          !turn off control on max mass rate
          this%spec%lmt_var(LMT_MASS_RATE_MAX) = PETSC_FALSE
        end if 
      end if
      if(this%spec%lmt_var(LMT_VOL_RATE_MAX)) then
        rate_max = this%flow_condition%flow_well%rate%dataset%rarray(1)
        if( this%q_fld(option%oil_phase) > rate_max ) then
          print *, "oil_producer control switch: " // &
                   "BHP -> VOL_RATE, vol_rate - rate_max = ",&
                    & this%q_fld(option%oil_phase), rate_max
          this%q_fld(option%oil_phase) = rate_max
          cntrl_var_tmp = CNTRL_VAR_VOL_RATE   
          pass = PETSC_FALSE
          !add limit on BHP min
          this%spec%lmt_var(LMT_BHP_MIN) = PETSC_TRUE
          !turn off control on max_vol_rate
          this%spec%lmt_var(LMT_VOL_RATE_MAX) = PETSC_FALSE
        end if 
      end if
    case(CNTRL_VAR_MASS_RATE)
      if(this%spec%lmt_var(LMT_BHP_MIN)) then
        bhp_min = this%flow_condition%flow_well%pressure%dataset%rarray(1)
        if( this%pw_ref < bhp_min ) then
          print *, "oil_producer control switch: " // &
                   "LMT_MASS_RATE -> BHP, pw_well - bhp_min = ",&
                    & this%pw_ref, bhp_min
          this%pw_ref = bhp_min
          cntrl_var_tmp = CNTRL_VAR_BHP
          pass = PETSC_FALSE   
          !activate back the control on mass_rate
          this%spec%lmt_var(LMT_MASS_RATE_MAX) = PETSC_TRUE
          !turn off control on min pressure  
          this%spec%lmt_var(LMT_BHP_MIN) = PETSC_FALSE
        end if 
      end if
    case(CNTRL_VAR_VOL_RATE)
      if(this%spec%lmt_var(LMT_BHP_MIN)) then
        bhp_min = this%flow_condition%flow_well%pressure%dataset%rarray(1)
        if( this%pw_ref < bhp_min ) then
          print *, "oil_producer control switch: " // &
                   "LMT_VOL_RATE -> BHP, pw_well - bhp_min = ",&
                    & this%pw_ref, bhp_min
          this%pw_ref = bhp_min
          cntrl_var_tmp = CNTRL_VAR_BHP
          pass = PETSC_FALSE   
          !activate back the control on mass_rate
          this%spec%lmt_var(LMT_VOL_RATE_MAX) = PETSC_TRUE
          !turn off control on min pressure  
          this%spec%lmt_var(LMT_BHP_MIN) = PETSC_FALSE
        end if 
      end if

    !add here control on OWR (water cut) 

  end select


  ! update control variable
  this%spec%cntrl_var = cntrl_var_tmp

end subroutine WellOilProdLimitCheck

! ************************************************************************** !

#endif  
end module Well_OilProducer_class
!end of WELL_CLASS

