module Well_WaterInjector_class
#ifdef WELL_CLASS

  use PFLOTRAN_Constants_module
  use WellSpec_Base_class
  use Well_Base_class
  use Well_Flow_class
  use Well_FlowEnergy_class

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  type, public, extends(well_flow_energy_type) :: well_water_injector_type
    ! .................
  contains  ! add here type-bound procedure 
    procedure, public :: PrintMsg => PrintWatInj
    procedure, public  :: PrintOutputHeader => WellWatInjPrintOutputHeader
    procedure, public :: VarsExplUpdate => WellWatInjVarsExplUpdate
    procedure, public :: LimitCheck => WellWatInjLimitCheck
    procedure, public :: ConnDenUpdate => WellWatInjConnDenUpdate
    procedure, public :: HydrostaticUpdate => WatInjHydrostaticUpdate
    procedure, public :: ConnMob => WellWatInjConnMob
    procedure, public :: InitDensity => WatInjInitDensity
  end type  well_water_injector_type

  !public :: CreateTOilImsWell

contains

! ************************************************************************** !

subroutine PrintWatInj(this)

  implicit none

  class(well_water_injector_type) :: this

  write(*,*) "Well PrintWatInj Printing message"

end subroutine PrintWatInj

! ************************************************************************** !

subroutine WellWatInjPrintOutputHeader(this,output_option,file_unit)
  ! 
  ! Write header for well_TOilIms output file
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 05/18/16
  ! 
  use Output_Aux_module

  implicit none

  class(well_water_injector_type) :: this
  type(output_option_type), intent(in) :: output_option
  PetscInt, intent(in) :: file_unit

  !character(len=MAXWORDLENGTH) :: tunit

  !tunit = trim(output_option%tunit)

  write(*,*) "Well PrintOutputHeaderWellWatInj to be extended"
  

end subroutine WellWatInjPrintOutputHeader

! ************************************************************************** !

subroutine WatInjInitDensity(this,grid,option)
  !
  ! Init well water density using the pressure at the reference grid block
  !
  ! Author: Paolo Orsini (OpenGoSim)  
  ! Date : 6/08/2016
  !

  use Grid_module
  use Option_module
  use EOS_Water_module

  implicit none

  class(well_water_injector_type) :: this
  type(grid_type), pointer :: grid  
  type(option_type) :: option

  PetscMPIInt :: cur_w_myrank
  PetscInt :: ghost_cntrl_id,ierr
  PetscReal :: pw_ref_init,tw_ref_init 
  PetscReal :: dw_h2o_kg,dw_h2o_mol

  !TO DO - this should be generalised - leavign only the density calculation
  ! specific to the injector

  call MPI_Comm_rank(this%comm, cur_w_myrank, ierr )  

  if(this%cntr_rank == cur_w_myrank ) then
    ghost_cntrl_id = grid%nL2G(this%cntrl_lcell_id); 
    pw_ref_init = &
      this%flow_auxvars(ZERO_INTEGER,ghost_cntrl_id)%pres(option%liquid_phase)
    !for the temperature can also use injection temp 
    !the control grid block temp is used for consistency
    tw_ref_init = &
      this%flow_energy_auxvars(ZERO_INTEGER,ghost_cntrl_id)%temp

    call EOSWaterDensity(tw_ref_init,pw_ref_init, &
                          dw_h2o_kg,dw_h2o_mol,ierr) 
    this%dw_kg_ref(option%liquid_phase) = dw_h2o_kg
  end if

  call MPI_Bcast ( this%dw_kg_ref(option%liquid_phase),1, &
                   MPI_DOUBLE_PRECISION, this%cntr_rank, this%comm, ierr )


end subroutine WatInjInitDensity

! ************************************************************************** !

!subroutine WellWatInjVarsExplUpdate(this,grid,ss_fluxes,option)
subroutine WellWatInjVarsExplUpdate(this,grid,option)
  !
  ! Explicit update of well variable for a water injector
  !
  ! Author: Paolo Orsini (OpenGoSim)  
  ! Date : 6/04/2016
  !

  use Grid_module
  use Option_module
  use EOS_Water_module

  class(well_water_injector_type) :: this
  type(grid_type), pointer :: grid
  !PetscReal :: ss_fluxes(:,:)      !not curently used 
  type(option_type) :: option

  PetscReal :: enth_src_h2o, dw_h2o_kg,dw_h2o_mol
  PetscInt :: ierr

  !write(*,"('WellWatInjVarsExplUpdate d11 before = ',e10.4)"), this%flow_energy_auxvarsy(0,1)%den(1)
  !write(*,"('WellWatInjVarsExplUpdate d12 before = ',e10.4)"), this%flow_energy_auxvars(0,1)%den(2) 
  !write(*,"('WellWatInjVarsExplUpdate p11 before = ',e10.4)"), this%flow_energy_auxvars(0,1)%pres(1) 
  !write(*,"('WellWatInjVarsExplUpdate t1 before = ',e10.4)"), this%flow_energy_auxvars(0,1)%temp 

  !write(*,"('WellWatInjVarsExplUpdate rate = ',e10.4)"), &
  !        this%flow_condition%flow_well%rate%dataset%rarray(1)
  !write(*,"('WellWatInjVarsExplUpdate temp = ',e10.4)"), &
  !        this%flow_condition%flow_well%temperature%dataset%rarray(1) 

  !write(*,"('WellWatInjVarsExplUpdate press = ',e10.4)"), &
  !        this%flow_condition%flow_well%pressure%dataset%rarray(1) 

  ! note: cntrl_var is used instead of this%cntrl_var because the
  ! well control variable can change during the well press computation
  ! this%cntrl_var is the initial value 

  if(this%connection_set%num_connections == 0 ) return

  this%tw_ref = this%flow_condition%flow_well%temperature%dataset%rarray(1)

  select case(this%spec%cntrl_var) 
    case(CNTRL_VAR_BHP)
      this%pw_ref = this%flow_condition%flow_well%pressure%dataset%rarray(1)
      call this%QPhase(grid,option%liquid_phase,option)      
      call this%MRPhase(grid,option%liquid_phase,option)

    case(CNTRL_VAR_MASS_RATE)      
      !call this%PressRef(grid,option%liquid_phase,option)
      call this%PressRefMRInj(grid,option%liquid_phase,option)
      this%mr_fld(option%liquid_phase) = &
                  this%flow_condition%flow_well%rate%dataset%rarray(1)
      call this%QPhase(grid,option%liquid_phase,option)
 
    case(CNTRL_VAR_VOL_RATE)
      this%q_fld(option%liquid_phase) = &
                this%flow_condition%flow_well%rate%dataset%rarray(1)
      !call this%PressRef(grid,option%liquid_phase,option)
      call this%PressRefQ(grid,option%liquid_phase,option)
      call this%MRPhase(grid,option%liquid_phase,option)
  end select

  call EOSWaterDensity(this%tw_ref,this%pw_ref, &
                       dw_h2o_kg,dw_h2o_mol,ierr) 
  !call EOSWaterEnthalpy(toil_auxvar%temp,cell_pressure,toil_auxvar%H(lid),ierr)

  this%dw_kg_ref(option%liquid_phase) = dw_h2o_kg


end subroutine WellWatInjVarsExplUpdate

! ************************************************************************** !

subroutine WellWatInjLimitCheck(this,pass,option)
  !
  ! Perform limit check for a water injector
  !
  ! Author: Paolo Orsini (OpenGoSim)  
  ! Date : 6/06/2016
  !

  use Option_module

  implicit none

  class(well_water_injector_type) :: this
  PetscBool :: pass
  type(option_type) :: option

  PetscInt :: cntrl_var_tmp 
  PetscReal :: press_max

  pass = PETSC_TRUE
  !if no changes cntrl_var mantains its initial value
  cntrl_var_tmp = this%spec%cntrl_var

  select case(this%spec%cntrl_var)
 
    case(CNTRL_VAR_MASS_RATE)
      if(this%spec%lmt_var(LMT_BHP_MAX)) then
        press_max = this%flow_condition%flow_well%pressure%dataset%rarray(1) 
        if(this%pw_ref > press_max) then
          print *, "water_injector control switch: " // &
                   "MASS_RATE -> BHP, pw_ref/pw_max= ",&
                    & this%pw_ref, press_max
          ! updates after check pw_ref after check
          this%pw_ref = press_max 
          cntrl_var_tmp = CNTRL_VAR_BHP
          pass = PETSC_FALSE
        end if
      end if 

  end select

  ! update control variable
  this%spec%cntrl_var = cntrl_var_tmp

end subroutine WellWatInjLimitCheck

! ************************************************************************** !

subroutine WellWatInjConnDenUpdate(this,grid,option)
  !
  ! Compute connection densities for a water injector 
  !
  ! Author: Paolo Orsini (OpenGoSim)  
  ! Date : 6/06/2016
  !
  use Grid_module
  use Option_module

  use EOS_Water_module

  implicit none

  class(well_water_injector_type) :: this
  type(grid_type), pointer :: grid !not currently used
  type(option_type) :: option

  PetscReal :: dw_kg_injw, dw_mol_injw   
  PetscInt :: ierr

  if( this%connection_set%num_connections == 0 ) return

  call EOSWaterDensity(this%tw_ref,this%pw_ref, &
                       dw_kg_injw,dw_mol_injw,ierr) 

  ! all well conns are assigned the same density
  ! a more accurate iterative model could be implemented 
  ! (e.g. hydrostatic coupler)

  this%conn_den_kg = dw_kg_injw !conn_densities is an array

end subroutine WellWatInjConnDenUpdate

! ************************************************************************** !

!subroutine WatInjHydrostaticUpdate(this,grid,ss_fluxes,option)
subroutine WatInjHydrostaticUpdate(this,grid,option)
  !
  ! computes hydrostatic corrections for injectors computing first 
  ! the pressure/densities profiles as for the hydrostatic couplers 
  !
  ! Author: Paolo Orsini (OpenGoSim)  
  ! Date : 6/17/2016
  !
  use Hydrostatic_Common_module
  use Grid_module
  use Option_module

  implicit none

  class(well_water_injector_type) :: this
  type(grid_type), pointer :: grid
  !PetscReal :: ss_fluxes(:,:)
  type(option_type) :: option

  PetscReal :: xm_nacl
  PetscReal :: dist_x, dist_y, dist_z
  PetscReal :: dist_z_for_pressure
  !PetscReal :: pw_conn, po_conn
  PetscReal :: phase_pres(option%nphase)
  PetscReal :: dummy_pres_grad(3)
  PetscInt :: iconn, local_id, ghosted_id
  PetscInt :: ipressure
  !
  !PetscReal :: q_sum
  !PetscReal :: q_ph(option%nphase)
  PetscReal :: den_ph(option%nphase)
  !PetscReal :: den_ph_w(option%nphase)
  PetscReal :: ph_w(option%nphase)
  PetscReal :: conn_press
  PetscInt :: i_ph

  xm_nacl = option%m_nacl * FMWNACL
  xm_nacl = xm_nacl /(1.d3 + xm_nacl)

  !update well temperatures on local and global connections,
  !and on the vertical finer grid performing an interpolation
  !from flow_energy_auxvars 
  call this%TempUpdate(grid,option)

  do i_ph = 1,option%nphase

    call PhaseHydrostaticPressure(this%fine_grid,option%gravity, &
               option%phase_map(i_ph),this%pw_ref, &
               this%fine_grid%idatum,xm_nacl,this%well_fine_grid_temp, &
               this%well_fine_grid_pres(i_ph,:), &
               this%well_fine_grid_den_kg(i_ph,:) )
  end do

  dist_x = 0.0d0;
  dist_y = 0.0d0;
  dummy_pres_grad = 0.0d0;
  do iconn=1,this%connection_set%num_connections
    local_id = this%connection_set%id_dn(iconn)
    ghosted_id = grid%nL2G(local_id)
    dist_z = grid%z(ghosted_id) - this%z_pw_ref
    ipressure = this%fine_grid%idatum+int(dist_z/this%fine_grid%delta_z)
    dist_z_for_pressure = grid%z(ghosted_id) - this%fine_grid%z(ipressure) 
    do i_ph = 1,option%nphase
      phase_pres(i_ph) = PressInterp(ipressure,dist_x,dist_y, &
                          dist_z_for_pressure,option%gravity, &
                          this%well_fine_grid_pres(i_ph,:), &
                          this%well_fine_grid_den_kg(i_ph,:), &
                          dummy_pres_grad)
    end do
    !use saturation to weighting phase densities   
    do i_ph = 1,option%nphase
      ph_w(i_ph) = this%flow_auxvars(ZERO_INTEGER,ghosted_id)%sat(i_ph)
      den_ph(i_ph) = this%flow_auxvars(ZERO_INTEGER,ghosted_id)%den_kg(i_ph)
    end do

    !compute well fluid density - not used in computation (printing purposes)
    this%conn_den_kg(iconn) = den_ph(option%liquid_phase)

    !OPTION - 1
    !since ph_w(i_ph) are saturarions, when injecting water  
    !into a domain saturated by another phase only, the pressure profile
    !follows the dominating phase - this is not rigorous but help
    !stability and avoid initial reverse flow  during pressure build up 
    conn_press = 0.0d0
    do i_ph = 1,option%nphase
      conn_press = conn_press + phase_pres(i_ph) * ph_w(i_ph)
    end do
    !OPTION - 2 
    !commented below a more rigorous option where the well pressure follows
    !the pressure profile of the phase beinj injected 
    !conn_press = phase_pres(option%liquid_phase)

    this%conn_h(iconn) = conn_press - this%pw_ref

  end do


end subroutine WatInjHydrostaticUpdate

! ************************************************************************** !

function WellWatInjConnMob(this,mobility,iphase)
  !  
  ! Compute well connection mobility for mphase mode 
  ! For injection equal to total mobility (summ of all phase mobilities) of
  ! the perforated grid block
  !
  ! example of a function common to all injectors 
  ! (should create an injector class)
  ! otherwise use select case a move this to flow
  ! Other option: define different routines for injector and producers in
  !               well_flow 
  !
  ! for efficieny - ConnMob can be extended to well_xxx_mode. 
  ! Instead of performing a loop, a simple sum can be used 
  ! E.g. WellWatInjConnMob = mobility(liquid_phase) + mobility(oil_phase)
  !
  ! Author: Paolo Orsini (OpenGoSim)  
  ! Date : 1/07/2015
  !

  implicit none

  class(well_water_injector_type) :: this
  PetscInt :: iphase  !not required for this WellWatInjConnMob extension
  PetscReal :: mobility(:)

  PetscReal :: WellWatInjConnMob

  PetscInt :: iph

  WellWatInjConnMob = 0.0d0
  do iph=1,size(mobility)
    WellWatInjConnMob = WellWatInjConnMob + mobility(iph)
  end do


end function WellWatInjConnMob
!*****************************************************************************!

#endif 
end module Well_WaterInjector_class
!end of WELL_CLASS



