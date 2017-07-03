module Well_FlowEnergy_class
#ifdef WELL_CLASS

  use PFLOTRAN_Constants_module
  use Well_Base_class
  use Well_Flow_class
  use AuxVars_FlowEnergy_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  type, public, extends(well_flow_type) :: well_flow_energy_type
      PetscReal :: tw_ref                          ! [°C] well temperature at reference elevation
      PetscReal, pointer :: ent_ref(:)             ! ent_ref(iphase) MJ/kmol, well fluid enthalpy of iphase
      PetscReal, pointer :: conn_temp(:)           ! local well temperatures [°C]
      PetscReal, pointer :: well_conn_temp(:)      ! temperature on the global (entire) well [°C]
      PetscReal, pointer :: well_fine_grid_temp(:) ! temperature value on fine grid, well_fine_grid_temp(inode) [°C]
      class(auxvar_flow_energy_type), pointer :: flow_energy_auxvars(:,:)
  contains  ! add here type-bound procedure 
    procedure, public :: PrintMsg => PrintFlowEnergy
    procedure, public :: ConnInit => WellFlowEnergyConnInit
    procedure, public :: InitRun => WellFlowEnergyInitRun
    procedure, public :: InitTimeStep => WellFlowEnergyInitTimeStep
    procedure, public :: VarsExplUpdate => FlowEnergyVarsExplUpdate
    procedure, public :: ExplJDerivative => WellFlowEnergyExplJDerivative
    procedure, public :: AverageTemp => WellFlowEnergyAverageTemp
    procedure, public :: TempUpdate => FlowEnergyTempUpdate
    procedure, public :: OneDimGridVarsSetup => WellFlowEnergy1DGridVarsSetup
    procedure, public :: HydrostaticUpdate => FlowEnergyHydrostaticUpdate
  end type  well_flow_energy_type

  public :: WellFlowEnergyInit, FlowEnergyWellStrip

contains

! ************************************************************************** !

subroutine PrintFlowEnergy(this)

  implicit none

  class(well_flow_energy_type) :: this

  write(*,*) "Well FlowEnergy Printing message"

end subroutine PrintFlowEnergy

! ************************************************************************** !

subroutine WellFlowEnergyInit(this,option)
  ! 
  ! Initializes variables/objects in flow and energy well class
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 5/20/2016
  ! 

  use Option_module

  implicit none

  class(well_flow_energy_type) :: this
  type(option_type) :: option

  this%tw_ref = 0.0d0; 

  allocate( this%ent_ref(option%nphase) );
  this%ent_ref = 0.0d0;

  nullify(this%flow_energy_auxvars);
  nullify(this%well_conn_temp);

end subroutine WellFlowEnergyInit

! ************************************************************************** !

subroutine WellFlowEnergyConnInit(this,num_connections,option)
  ! 
  ! Allocate and initilize well_flow_energy connections arrays
  ! 
  ! Author: Paolo Orsini - OpenGoSim
  ! Date: 6/17/2016
  !

  use Option_module

  implicit none

  class(well_flow_energy_type) :: this
  PetscInt, intent(in) :: num_connections 
  type(option_type) :: option  

  call WellFlowConnInit(this,num_connections,option);

  nullify(this%conn_temp);
  allocate(this%conn_temp(num_connections));
  this%conn_temp = 0.0d0; 

end subroutine WellFlowEnergyConnInit

! ************************************************************************** !

subroutine WellFlowEnergyInitRun(this,grid,material_auxvars, &
                                 output_option,option)
  ! 
  ! Initialise well for a run
  ! 
  ! Author: Paolo Orsini
  ! Date: 4/08/16
  ! 

  use Grid_module
  use Material_Aux_class, only : material_auxvar_type
  use Output_Aux_module
  use Option_module

  implicit none

  class(well_flow_energy_type) :: this
  type(grid_type), pointer :: grid
  type(material_auxvar_type), intent(in) :: material_auxvars(:)
  type(output_option_type), intent(in) :: output_option
  type(option_type) :: option

  call WellBaseInitRun(this,grid,material_auxvars,output_option,option)

  !initialize pressure and well densities for injectectors
  call this%InitDensity(grid,option )

  call this%ExplUpdate(grid,option)
  !init well temperature 
  call this%TempUpdate(grid,option)
  !init well hydrostatic corrections
  call this%HydroCorrUpdates(grid,option)

  !update the pressure again after H correction, 
  ! only to print the right value at t=0
  ! move to InitRun/well_last_extension - as in 2. above
  call this%ExplUpdate(grid,option)

end subroutine WellFlowEnergyInitRun

! ************************************************************************** !

subroutine WellFlowEnergyInitTimeStep(this,grid,material_auxvars,option)
  ! 
  ! Initialise well time step
  ! 
  ! Author: Paolo Orsini
  ! Date: 4/08/16
  ! 

  use Grid_module
  use Material_Aux_class, only : material_auxvar_type
  use Option_module

  implicit none
  
  class(well_flow_energy_type) :: this
  type(grid_type), pointer :: grid
  type(material_auxvar_type), intent(in) :: material_auxvars(:)
  type(option_type) :: option

  ! update well connection factors if variable permeability 
  call WellBaseInitTimeStep(this,grid,material_auxvars,option)

  call this%TempUpdate(grid,option)

  call this%HydroCorrUpdates(grid,option)   

end subroutine WellFlowEnergyInitTimeStep

! ************************************************************************** !

subroutine WellFlowEnergy1DGridVarsSetup(this,option)

  use Option_module

  implicit none

  class(well_flow_energy_type) :: this
  type(option_type) :: option

  call WellFlow1DGridVarsSetup(this,option) 

  allocate(this%well_fine_grid_temp(size(this%fine_grid%z(:))))
  this%well_fine_grid_temp = 0.0d0 

end subroutine WellFlowEnergy1DGridVarsSetup

! ************************************************************************** !

subroutine WellFlowEnergyExplJDerivative(this,iconn,ghosted_id,isothermal, &
                                         energy_equation_index,option,Jac)
  ! 
  ! Computes the well derivatives terms for the jacobian
  ! 
  ! Author: Paolo Orsini
  ! Date: 6/06/16
  ! 

  use Option_module
  !use Condition_module

  implicit none

  class(well_flow_energy_type) :: this
  PetscInt :: iconn 
  PetscInt :: ghosted_id
  PetscBool :: isothermal
  PetscInt :: energy_equation_index
  type(option_type) :: option
  !PetscReal :: Jac(option%nflowdof,option%nflowdof)
  PetscReal :: Jac(:,:)

  !type(flow_toil_ims_condition_type), pointer :: src_sink_condition
  !type(toil_ims_auxvar_type) :: toil_auxvar(0:)
  !class(auxvar_toil_ims_type) :: toil_auxvar(0:)
  !type(auxvar_toil_ims_type) :: toil_auxvar(0:)
  !type(global_auxvar_type) :: global_auxvar
  !PetscReal :: scale
  
  
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscReal :: dummy_real(option%nphase)
  PetscInt :: idof, irow

  option%iflag = -3

  !call TOilImsSrcSink(option,src_sink_condition,toil_auxvar(ZERO_INTEGER), &
  !                        global_auxvar,dummy_real,scale,Res)

#ifdef WELL_DEBUG
  write(*,"('ExplJDerivative p011 = ',e16.10)") this%flow_energy_auxvars(0,1)%pres(1)
  write(*,"('ExplJDerivative p111 = ',e16.10)") this%flow_energy_auxvars(1,1)%pres(1)
  write(*,"('ExplJDerivative p211 = ',e16.10)") this%flow_energy_auxvars(2,1)%pres(1)
  write(*,"('ExplJDerivative p311 = ',e16.10)") this%flow_energy_auxvars(3,1)%pres(1) 
#endif


  call this%ExplRes(iconn,dummy_real,isothermal,ghosted_id,ZERO_INTEGER,&
                    option,res)

  ! downgradient derivatives
  do idof = 1, option%nflowdof

    !call TOilImsSrcSink(option,src_sink_condition,toil_auxvar(idof), &
    !                    global_auxvar,dummy_real,scale,res_pert)
    call this%ExplRes(iconn,dummy_real,isothermal,ghosted_id,idof, &
                      option,res_pert)

    do irow = 1, option%nflowdof
      !Jac(irow,idof) = (res_pert(irow)-res(irow))/toil_auxvar(idof)%pert
      Jac(irow,idof) = (res_pert(irow)-res(irow)) / &
                         this%flow_energy_auxvars(idof,ghosted_id)%pert
    enddo !irow
  enddo ! idof
  
  if (isothermal) then
    !Jac(TOIL_IMS_ENERGY_EQUATION_INDEX,:) = 0.d0
    !Jac(:,TOIL_IMS_ENERGY_EQUATION_INDEX) = 0.d0
    Jac(energy_equation_index,:) = 0.d0
    Jac(:,energy_equation_index) = 0.d0
  endif   

    !Jac =  0.0d0

end subroutine WellFlowEnergyExplJDerivative

! ************************************************************************** !

!subroutine WellFlowEnergyAverageTemp(this,grid,ss_fluxes,option)
subroutine WellFlowEnergyAverageTemp(this,grid,option)
  !
  ! Compute connection densities for producers
  ! to be overwritten for injectors 
  ! Compute well fluid average density from the well segment phase densities 
  ! using grid-blocks/well fluxes as weights 
  !
  ! Author: Paolo Orsini (OpenGoSim)  
  ! Date : 6/12/2016
  !
  use Grid_module
  use Option_module

  implicit none

  class(well_flow_energy_type) :: this
  type(grid_type), pointer :: grid 
  !PetscReal :: ss_fluxes(:,:)      
  type(option_type) :: option

  PetscReal, pointer :: ss_fluxes(:,:)
  PetscReal :: q_sum
  PetscReal :: q_ph(option%nphase)

  PetscInt :: iconn, local_id, ghosted_id, ierr
  PetscInt :: i_ph

  PetscReal :: q_sum_lc, q_sum_well 
  PetscReal :: temp_q_lc, temp_q_well
  PetscReal :: temp_lc, temp_well
 
  q_sum_lc = 0.0d0
  q_sum_well = 0.0d0
  temp_q_lc = 0.0d0
  temp_q_well = 0.0d0
  temp_lc = 0.0d0
  temp_well = 0.0d0

  ss_fluxes => this%ss_flow_vol_fluxes

  do iconn = 1,this%connection_set%num_connections
    local_id = this%connection_set%id_dn(iconn)
    ghosted_id = grid%nL2G(local_id)
    ! need to change signs because in the PFLOTRAN convention
    ! producing wells have ss_fluxes < 0
    q_sum = 0.0d0
    do i_ph = 1,option%nphase
      q_ph(i_ph) = 0.0d0
      if( ss_fluxes(i_ph,iconn) < 0.0d0) then
        !should not have ss_fluxes > 0.0d0, this reverse flow situation
        ! should be detected in Res computation 
        q_ph(i_ph) = -1.0d0 * ss_fluxes(i_ph,iconn)
      end if
      q_sum = q_sum + q_ph(i_ph)
    end do

    q_sum_lc = q_sum_lc + q_sum     
    temp_q_lc = temp_q_lc + &
       q_sum * this%flow_energy_auxvars(ZERO_INTEGER,ghosted_id)%temp 
    temp_lc = temp_lc + &
       this%flow_energy_auxvars(ZERO_INTEGER,ghosted_id)%temp
  end do

  call MPI_ALLREDUCE(q_sum_lc, q_sum_well, 1, MPI_DOUBLE_PRECISION, &
                     MPI_SUM,this%comm, ierr)

  call MPI_ALLREDUCE(temp_q_lc, temp_q_well, 1, MPI_DOUBLE_PRECISION, &
                     MPI_SUM,this%comm, ierr)

  call MPI_ALLREDUCE(temp_lc, temp_well, 1, MPI_DOUBLE_PRECISION, &
                     MPI_SUM,this%comm, ierr)
  
  if ( q_sum_well > wfloweps ) then
    this%tw_ref = temp_q_well / q_sum_well
  else
    this%tw_ref = temp_well / dble(this%well_num_conns)
  end if


end subroutine WellFlowEnergyAverageTemp

! ************************************************************************** !

subroutine FlowEnergyHydrostaticUpdate(this,grid,option)
  !
  ! computes hydrostatic corrections for producers computing first 
  ! the pressure/densities profiles as for the hydrostatic couplers 
  ! to be overwritten for injectors - where only the phase being injected 
  ! should be considered, at least in absence of cross flows  
  !
  ! Author: Paolo Orsini (OpenGoSim)  
  ! Date : 6/17/2016
  !
  use Hydrostatic_Common_module
  use Grid_module
  use Option_module

  implicit none

  class(well_flow_energy_type) :: this
  type(grid_type), pointer :: grid
  type(option_type) :: option

  PetscReal, pointer :: ss_fluxes(:,:)
  PetscReal :: xm_nacl
  PetscReal :: dist_x, dist_y, dist_z
  PetscReal :: dist_z_for_pressure
  !PetscReal :: pw_conn, po_conn
  PetscReal :: phase_pres(option%nphase)
  PetscReal :: dummy_pres_grad(3)
  PetscInt :: iconn, local_id, ghosted_id
  PetscInt :: ipressure
  !
  PetscReal :: q_sum
  PetscReal :: q_ph(option%nphase)
  PetscReal :: den_ph(option%nphase)
  !PetscReal :: den_ph_w(option%nphase)
  PetscReal :: ph_w(option%nphase)
  PetscReal :: conn_press
  PetscInt :: i_ph


  xm_nacl = option%m_nacl * FMWNACL
  xm_nacl = xm_nacl /(1.d3 + xm_nacl)

  !one dim grid already available

  ss_fluxes => this%ss_flow_vol_fluxes 

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

    !could move the fluxes weigthing factors computation to well_flow                         
    q_sum = 0.0d0
    do i_ph = 1,option%nphase
      q_ph(i_ph) = 0.0d0
      if( ss_fluxes(i_ph,iconn) < 0.0d0) then
        !should not have ss_fluxes > 0.0d0, this reverse flow situation
        ! should be detected in Res computation 
        q_ph(i_ph) = -1.0d0 * ss_fluxes(i_ph,iconn)
      end if
      q_sum = q_sum + q_ph(i_ph)
      den_ph(i_ph) = this%flow_auxvars(ZERO_INTEGER,ghosted_id)%den_kg(i_ph)
    end do
    if( q_sum > wfloweps ) then
      do i_ph = 1,option%nphase
        ph_w(i_ph) = q_ph(i_ph) / q_sum  
      end do
    else ! OK for nill fluxes 
         ! not accurate for back flow (should get warning in Res computation)
      do i_ph = 1,option%nphase
        ph_w(i_ph) = this%flow_auxvars(ZERO_INTEGER,ghosted_id)%sat(i_ph)
      end do
    end if   
    !end phase weigths phase computation  
    !compute fluid average density 
    this%conn_den_kg(iconn) = 0.0d0
    conn_press = 0.0d0
    do i_ph = 1,option%nphase
      this%conn_den_kg(iconn) = this%conn_den_kg(iconn) + &
                                den_ph(i_ph) * ph_w(i_ph)
      conn_press = conn_press + phase_pres(i_ph) * ph_w(i_ph)
    end do
    this%conn_h(iconn) = conn_press - this%pw_ref
    !this%conn_h(iconn) = phase_pres(2)
  end do


end subroutine FlowEnergyHydrostaticUpdate


! ************************************************************************** !

subroutine FlowEnergyTempUpdate(this,grid,option)
  !
  !update well flow temperature from flow_energy_auxvars
  !to be overwritten for injectors
  !
  ! Author: Paolo Orsini (OpenGoSim)  
  ! Date : 6/12/2016

  use Grid_module
  use Option_module

  implicit none

  class(well_flow_energy_type) :: this
  type(grid_type), pointer :: grid
  type(option_type) :: option

  PetscInt :: iconn, local_id, ghosted_id, ierr

  !update local well temperature connections from auxvars
  do iconn = 1,this%connection_set%num_connections
    local_id = this%connection_set%id_dn(iconn)
    ghosted_id = grid%nL2G(local_id)
    this%conn_temp(iconn) = &
      this%flow_energy_auxvars(ZERO_INTEGER,ghosted_id)%temp
  end do

  call MPI_Allgatherv(this%conn_temp,this%connection_set%num_connections, &
              MPI_DOUBLE_PRECISION, this%well_conn_temp, this%w_rank_conn, &
              this%disp_rank_conn,MPI_DOUBLE_PRECISION, this%comm,ierr)  
  
  if (this%hydrostatic_method == WELL_HYDROSTATIC_ITERATIVE) then
    !perform interpolation on the well_fine_grid
    call this%fine_grid%InterpFromWellConnTo1DGrid(this%w_conn_z, &
                         this%w_conn_order, &
                         this%well_conn_temp,this%well_fine_grid_temp)
  end if

end subroutine FlowEnergyTempUpdate

! ************************************************************************** !

!subroutine FlowEnergyVarsExplUpdate(this,grid,ss_fluxes,option)
subroutine FlowEnergyVarsExplUpdate(this,grid,option)

  use Grid_module
  use Option_module

  implicit none

  class(well_flow_energy_type) :: this
  type(grid_type), pointer :: grid
  !PetscReal :: ss_fluxes(:,:)
  type(option_type) :: option

  print *, "FlowEnergyVarsExplUpdate must be extended"
  stop

end subroutine FlowEnergyVarsExplUpdate


!*****************************************************************************!

!function WellFlowEnergyConnMob(this,mobility,iphase)
!
!  implicit none
!
!  class(well_flow_energy_type) :: this
!  PetscInt :: iphase  
!  PetscReal :: mobility(:)
!
!  PetscReal :: WellFlowEnergyConnMob
!
!  print *, "WellFlowEnergyConnMob must be extended"
!  stop
!
!end function WellFlowEnergyConnMob
!*****************************************************************************!

subroutine FlowEnergyWellStrip(well)
  !
  ! Strip well_flow_energy and all parent members
  !
  ! Author: Paolo Orsini (OpenGoSim)  
  ! Date : 6/14/2016
  !

  use Utility_module, only : DeallocateArray 

  implicit none

  class(well_flow_energy_type) :: well

  call DeallocateArray(well%ent_ref)
  !only pointers to auxvars 
  nullify(well%flow_energy_auxvars)

  call DeallocateArray(well%well_fine_grid_temp)
  call DeallocateArray(well%conn_temp)
  call DeallocateArray(well%well_conn_temp)

  ! this will strip all its parents too
  call FlowWellStrip(well) 


end subroutine FlowEnergyWellStrip

!*****************************************************************************!

#endif  
end module Well_FlowEnergy_class
!end of WELL_CLASS

