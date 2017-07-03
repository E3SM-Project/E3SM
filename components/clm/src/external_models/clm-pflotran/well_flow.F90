module Well_Flow_class
#ifdef WELL_CLASS

  use PFLOTRAN_Constants_module
  use WellSpec_Base_class
  use Well_Base_class
  use AuxVars_Flow_module
  use Condition_module
  implicit none

  private

#include "petsc/finclude/petscsys.h"

  PetscReal, parameter, public :: wfloweps = 1.D-24 
  PetscInt, parameter, public :: WELL_HYDROSTATIC_LINEAR = 1
  PetscInt, parameter, public :: WELL_HYDROSTATIC_ITERATIVE = 2

  type, public, extends(well_base_type) :: well_flow_type
    PetscReal :: pw_ref                          ! [Pa] well pressure at reference elevation
    PetscInt :: hydrostatic_method               ! it defines the mthod used in the well hydrostatic computation
    PetscReal, pointer :: dw_kg_ref(:)           ! dw_kg_ref(iphase) [kg/m3] well fluid density of iphase at reference elevation
    PetscReal, pointer :: q_fld(:)               ! q_fld(iphase)  [m3/s] well fluid flow rates of iphase
    PetscReal, pointer :: mr_fld(:)              ! mr_fld(iphase) [kg/s] well fluid mass rates of iphase
    PetscReal, pointer :: conn_h(:)              ! connection hydrostatic pressure corrections (local)
    PetscReal, pointer :: conn_den_kg(:)         ! connection densities (local)
    PetscReal, pointer :: well_conn_den_kg(:)    ! connection densities (global well) [kg/m3]
    PetscReal, pointer :: well_conn_h_sorted(:)  ! connection hydrostatic pressure corrections sorted by ascending elevation (entire well)
    PetscReal, pointer :: well_fine_grid_pres(:,:) !pressure values on fine grid, well_fine_grid_pres(iphase,inode)
    PetscReal, pointer :: well_fine_grid_den_kg(:,:) !density values on fine grid, well_fine_grid_den_kg(iphase,inode)
    PetscReal, pointer :: ss_flow_vol_fluxes(:,:) ! local volumetric fluxes (m^3/s), point to a section of patch%ss_flow_vol_fluxes 
    class(auxvar_flow_type), pointer :: flow_auxvars(:,:) !pointer to flow auxvars
    type(flow_condition_type), pointer :: flow_condition ! pointer to flow_condition associated with the well
    !PetscReal, pointer :: conn_mobs(:,:)       ! well connection mobilities ! TO REMOVE - computed when needed flight
  contains  ! add here type-bound procedure 
    procedure, public :: PrintMsg => PrintFlow
    procedure, public :: ConnInit => WellFlowConnInit
    procedure, public :: ExplUpdate => FlowExplUpdate
    procedure, public :: VarsExplUpdate => FlowVarsExplUpdate
    procedure, public :: ConnMob => WellFlowConnMob
    procedure, public :: PressRef => FlowPressRef
    procedure, public :: PressRefQ => FlowPressRefQ  
    procedure, public :: PressRefMRInj => FlowPressRefMRInj
    procedure, public :: PressRefMRProd => FlowPressRefMRProd 
    procedure, public :: QPhase => FlowQPhase
    procedure, public :: MRPhase => FlowMRPhase
    procedure, public :: LimitCheck => WellFlowLimitCheck
    procedure, public :: HydroCorrUpdates => FlowHydroCorrUpdate
    procedure, public :: ConnDenUpdate => WellFlowConnDenUpdate
    procedure, public :: InitDensity => WellFlowInitDensity
    procedure, public :: HydrostaticUpdate => FlowHydrostaticUpdate
    procedure, public :: OneDimGridVarsSetup => WellFlow1DGridVarsSetup
    procedure, public :: DataOutput => FlowDataOutput
  end type  well_flow_type

  public :: WellFlowInit, FlowWellStrip, WellFlow1DGridVarsSetup, &
            WellFlowConnInit

contains

! ************************************************************************** !

subroutine PrintFlow(this)

  implicit none

  class(well_flow_type) :: this

  write(*,*) "Well Flow Printing message"

end subroutine PrintFlow

! ************************************************************************** !

subroutine WellFlowInit(this,option)
  ! 
  ! Initializes variables/objects in flow well class
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 05/18/16
  ! 

  use Option_module

  implicit none

  class(well_flow_type) :: this

  type(option_type) :: option

  this%pw_ref = 0.0d0;
  !to be overwritten by input (for now default to linear)
  this%hydrostatic_method = WELL_HYDROSTATIC_LINEAR   

  allocate( this%dw_kg_ref(option%nphase) );
  this%dw_kg_ref = 0.0d0;
  allocate( this%q_fld(option%nphase) );
  this%q_fld = 0.0d0;
  allocate( this%mr_fld(option%nphase) );
  this%mr_fld = 0.0d0;

  nullify(this%ss_flow_vol_fluxes)
  nullify(this%flow_auxvars) 
  nullify(this%well_fine_grid_pres) 

end subroutine WellFlowInit

! ************************************************************************** !

subroutine WellFlowConnInit(this,num_connections,option)
  ! 
  ! Allocate and initilize well_flow connections arrays
  ! 
  ! Author: Paolo Orsini - OpenGoSim
  ! Date: 5/20/2016
  !

  use Option_module

  implicit none

  class(well_flow_type) :: this
  PetscInt, intent(in) :: num_connections 
  type(option_type) :: option  

  call WellBaseConnInit(this,num_connections,option);

  nullify(this%conn_h);
  allocate(this%conn_h(num_connections));
  this%conn_h = 0.0d0; 

  nullify(this%conn_den_kg);
  allocate( this%conn_den_kg(num_connections) )  
  this%conn_den_kg =0.0d0

end subroutine wellFlowConnInit

! ************************************************************************** !

subroutine WellFlow1DGridVarsSetup(this,option)

  use Option_module

  implicit none

  class(well_flow_type) :: this 
  type(option_type) :: option

  allocate(this%well_fine_grid_pres(option%nphase,size(this%fine_grid%z(:))))
  this%well_fine_grid_pres = 0.0d0
  allocate(this%well_fine_grid_den_kg(option%nphase,size(this%fine_grid%z(:))))
  this%well_fine_grid_den_kg = 0.0d0   

end subroutine WellFlow1DGridVarsSetup

! ************************************************************************** !

subroutine FlowExplUpdate(this,grid,option)
  ! 
  ! - Update FlowEnergy well vars
  ! - Perform a limit on well checks 
  ! - Update well control variable in case of switch when a limit is reached
  !
  ! Author: Paolo Orsini (OGS)
  ! Date: 6/03/2016
  ! 

  use Grid_module
  use Option_module

  implicit none

  class(well_flow_type) :: this
  type(grid_type), pointer :: grid
  type(option_type) :: option

  PetscBool :: pass

#ifdef WELL_DEBUG
  write(*,"('FlowExpl d11 before = ',e10.4)") this%flow_auxvars(0,1)%den(1)
  write(*,"('FlowExpl d12 before = ',e10.4)") this%flow_auxvars(0,1)%den(2) 
  write(*,"('FlowExpl p11 before = ',e10.4)") this%flow_auxvars(0,1)%pres(1) 
  !write(*,"('FlowExpl t1 before = ',e10.4)") this%flow_auxvars(0,1)%temp 
#endif

  if(this%connection_set%num_connections == 0 ) return

  pass = PETSC_FALSE

  !cntrl_var = this%cntrl_var ! initialise well control variable
  do
    if(pass) exit ! the well limits are satisfied

    call this%VarsExplUpdate(grid,option)

    call this%LimitCheck(pass,option)

    ! well check against VFPs - is pw admissible for volumtric rates given in VFPs?

  end do

end subroutine FlowExplUpdate

! ************************************************************************** !

!subroutine FlowVarsExplUpdate(this,grid,ss_fluxes,option)
subroutine FlowVarsExplUpdate(this,grid,option)

  use Grid_module
  use Option_module

  implicit none

  class(well_flow_type) :: this
  type(grid_type), pointer :: grid
  !PetscReal :: ss_fluxes(:,:)
  type(option_type) :: option

  print *, "FlowVarsExplUpdate must be extended"
  stop

end subroutine FlowVarsExplUpdate

!*****************************************************************************!
subroutine FlowPressRef(this,grid,phase,option)
  !  
  ! Compute well p_ref given the volumetric rate of one phase 
  ! IPR sign convention: rate > 0 for fluid being produced
  ! Tested for production well only but should work also for injectors
  !
  ! Author: Paolo Orsini (OpenGoSim)  
  ! Date : 25/06/2015
  !
  use Grid_module
  use Option_module

  implicit none

  class(well_flow_type) :: this
  type(grid_type), pointer :: grid
  PetscInt :: phase
  type(option_type) :: option
  
  !PetscInt :: phase, cntrl_var, ivar
  !type(connection_set_type), pointer :: connection_set
  !type(mphase_auxvar_type), pointer :: auxvars(:)

  PetscReal :: rate 

  PetscReal :: press_div_loc, press_div
  PetscReal :: conn_loc, conn_tot
  !PetscReal :: p_well_lc, p_well
  PetscInt :: iconn, local_id, ghosted_id, iph, ierr 
  PetscReal :: mob
 
  if(this%connection_set%num_connections > 0 ) then

#ifdef WELL_DEBUG
  write(*,"('FlowPressRef p011 = ',e10.4)") this%flow_auxvars(0,1)%pres(1)
  write(*,"('FlowPressRef p111 = ',e10.4)") this%flow_auxvars(1,1)%pres(1)
  write(*,"('FlowPressRef p211 = ',e10.4)") this%flow_auxvars(2,1)%pres(1)
  write(*,"('FlowPressRef p311 = ',e10.4)") this%flow_auxvars(3,1)%pres(1) 
#endif
    
    rate =  this%flow_condition%flow_well%rate%dataset%rarray(1)

    conn_tot = 0.0d0 
    ! divisor computation
    press_div_loc = 0.0d0
    do iconn = 1, this%connection_set%num_connections
      local_id = this%connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      mob = this%ConnMob(this%flow_auxvars(ZERO_INTEGER,ghosted_id)%mobility, &
                         phase)
      select case(this%spec%cntrl_var)
        case(CNTRL_VAR_VOL_RATE)
          press_div_loc = press_div_loc + this%conn_factors(iconn) * mob  
        case(CNTRL_VAR_MASS_RATE)
          press_div_loc = press_div_loc + this%conn_factors(iconn) * mob * &
                      this%flow_auxvars(ZERO_INTEGER,ghosted_id)%den_kg(phase)
      end select 
    end do

    call MPI_ALLREDUCE(press_div_loc, press_div, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM,this%comm, ierr)

    ! compute local well connection flux contribution
    conn_loc = 0.0d0 
    do iconn = 1, this%connection_set%num_connections      
      local_id = this%connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      mob = this%ConnMob(this%flow_auxvars(ZERO_INTEGER,ghosted_id)%mobility, &
                         phase)
      if(this%conn_status(iconn) == CONN_STATUS_OPEN ) then
        ! vol_rate > 0 for fluid entering the well  
#ifdef WELL_DEBUG
  write(*,"('gh = ',I5,'FlowExpl = ',e10.4)") ghosted_id, &
           this%flow_auxvars(ZERO_INTEGER,ghosted_id)%pres(phase) 
#endif
        select case(this%spec%cntrl_var)
          case(CNTRL_VAR_VOL_RATE)
            conn_loc = conn_loc + this%conn_factors(iconn) * mob * &
                 ( this%flow_auxvars(ZERO_INTEGER,ghosted_id)%pres(phase) - &
                   this%conn_h(iconn) )
          case(CNTRL_VAR_MASS_RATE)
            conn_loc = conn_loc + this%conn_factors(iconn) * mob * &
                 this%flow_auxvars(ZERO_INTEGER,ghosted_id)%den_kg(phase) * &
                 ( this%flow_auxvars(ZERO_INTEGER,ghosted_id)%pres(phase) - &
                   this%conn_h(iconn) )
        end select 
      end if
    end do ! end loop on well connections
    ! connection fluxes contibutions from well_segments in other ranks
    call MPI_ALLREDUCE(conn_loc, conn_tot, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM,this%comm, ierr)

  end if
 
  ! rate can be volume or mass rate depending on the control variable
  this%pw_ref = (conn_tot - rate) / press_div

#ifdef WELL_DEBUG
  write(*,"('FlowPressRef pw_ref = ',e10.4)") this%pw_ref 
#endif
  write(*,"('FlowPressRef pw_ref = ',e16.10)") this%pw_ref

end subroutine FlowPressRef


! ************************************************************************** !

!*****************************************************************************!
subroutine FlowPressRefQ(this,grid,phase,option)
  !  
  ! Compute well p_ref given the volumetric rate of one phase 
  ! IPR sign convention: rate > 0 for fluid being produced
  ! Tested for production well only but should work also for injectors
  !
  ! Author: Paolo Orsini (OpenGoSim)  
  ! Date : 25/06/2015
  !
  use Grid_module
  use Option_module

  implicit none

  class(well_flow_type) :: this
  type(grid_type), pointer :: grid
  PetscInt :: phase
  type(option_type) :: option
  
  !PetscInt :: phase, cntrl_var, ivar
  !type(connection_set_type), pointer :: connection_set
  !type(mphase_auxvar_type), pointer :: auxvars(:)

  PetscReal :: rate 

  PetscReal :: press_div_loc, press_div
  PetscReal :: conn_loc, conn_tot
  !PetscReal :: p_well_lc, p_well
  PetscInt :: iconn, local_id, ghosted_id, iph, ierr 
  PetscReal :: mob
 
  if(this%connection_set%num_connections > 0 ) then
  
    rate =  this%flow_condition%flow_well%rate%dataset%rarray(1)

    conn_tot = 0.0d0 
    ! divisor computation
    press_div_loc = 0.0d0
    do iconn = 1, this%connection_set%num_connections
      local_id = this%connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      mob = this%ConnMob(this%flow_auxvars(ZERO_INTEGER,ghosted_id)%mobility, &
                         phase)

      press_div_loc = press_div_loc + this%conn_factors(iconn) * mob  

    end do

    call MPI_ALLREDUCE(press_div_loc, press_div, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM,this%comm, ierr)

    ! compute local well connection flux contribution
    conn_loc = 0.0d0 
    do iconn = 1, this%connection_set%num_connections      
      local_id = this%connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      mob = this%ConnMob(this%flow_auxvars(ZERO_INTEGER,ghosted_id)%mobility, &
                         phase)
      if(this%conn_status(iconn) == CONN_STATUS_OPEN ) then
        ! vol_rate > 0 for fluid entering the well  
        conn_loc = conn_loc + this%conn_factors(iconn) * mob * &
                 ( this%flow_auxvars(ZERO_INTEGER,ghosted_id)%pres(phase) - &
                 this%conn_h(iconn) )
      end if
    end do ! end loop on well connections
    ! connection fluxes contibutions from well_segments in other ranks
    call MPI_ALLREDUCE(conn_loc, conn_tot, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM,this%comm, ierr)

  end if
 
  ! rate can be volume or mass rate depending on the control variable
  this%pw_ref = (conn_tot - rate) / press_div

#ifdef WELL_DEBUG
  write(*,"('FlowPressRefQ pw_ref = ',e10.4)"), this%pw_ref 
#endif

end subroutine FlowPressRefQ

! ************************************************************************** !

!*****************************************************************************!
subroutine FlowPressRefMRInj(this,grid,phase,option)
  !  
  ! Compute well p_ref given the mass rate of one phase 
  ! IPR sign convention: rate > 0 for fluid being produced
  ! Tested for production well only but should work also for injectors
  !
  ! Author: Paolo Orsini (OpenGoSim)  
  ! Date : 25/06/2015
  !
  use Grid_module
  use Option_module

  implicit none

  class(well_flow_type) :: this
  type(grid_type), pointer :: grid
  PetscInt :: phase
  type(option_type) :: option
  
  !PetscInt :: phase, cntrl_var, ivar
  !type(connection_set_type), pointer :: connection_set
  !type(mphase_auxvar_type), pointer :: auxvars(:)

  PetscReal :: rate 

  PetscReal :: press_div_loc, press_div
  PetscReal :: conn_loc, conn_tot
  !PetscReal :: p_well_lc, p_well
  PetscInt :: iconn, local_id, ghosted_id, iph, ierr 
  PetscReal :: mob
 
  if(this%connection_set%num_connections > 0 ) then
 
    rate =  this%flow_condition%flow_well%rate%dataset%rarray(1)

    conn_tot = 0.0d0 
    ! divisor computation
    press_div_loc = 0.0d0
    do iconn = 1, this%connection_set%num_connections
      local_id = this%connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      mob = this%ConnMob(this%flow_auxvars(ZERO_INTEGER,ghosted_id)%mobility, &
                         phase)

      press_div_loc = press_div_loc + this%conn_factors(iconn) * mob 

    end do
    press_div_loc = press_div_loc * this%dw_kg_ref(phase) 

    call MPI_ALLREDUCE(press_div_loc, press_div, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM,this%comm, ierr)

    ! compute local well connection flux contribution
    conn_loc = 0.0d0 
    do iconn = 1, this%connection_set%num_connections      
      local_id = this%connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      mob = this%ConnMob(this%flow_auxvars(ZERO_INTEGER,ghosted_id)%mobility, &
                         phase)
      if(this%conn_status(iconn) == CONN_STATUS_OPEN ) then
        ! vol_rate > 0 for fluid entering the well  
        conn_loc = conn_loc + this%conn_factors(iconn) * mob * &
                 ( this%flow_auxvars(ZERO_INTEGER,ghosted_id)%pres(phase) - &
                   this%conn_h(iconn) )
      end if
    end do ! end loop on well connections
  
    conn_loc = conn_loc * this%dw_kg_ref(phase)
    
    
    ! connection fluxes contibutions from well_segments in other ranks
    call MPI_ALLREDUCE(conn_loc, conn_tot, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM,this%comm, ierr)

  end if
 
  ! rate can be volume or mass rate depending on the control variable
  this%pw_ref = (conn_tot - rate) / press_div

#ifdef WELL_DEBUG
  write(*,"('FlowPressRefMRInj pw_ref = ',e10.4)") this%pw_ref 
#endif
  write(*,"('FlowPressRefMRInj pw_ref = ',e10.4)") this%pw_ref

end subroutine FlowPressRefMRInj


! ************************************************************************** !

!*****************************************************************************!
subroutine FlowPressRefMRProd(this,grid,phase,option)
  !  
  ! Compute well p_ref given the volumetric rate of one phase 
  ! IPR sign convention: rate > 0 for fluid being produced
  ! Tested for production well only but should work also for injectors
  !
  ! Author: Paolo Orsini (OpenGoSim)  
  ! Date : 25/06/2015
  !
  use Grid_module
  use Option_module

  implicit none

  class(well_flow_type) :: this
  type(grid_type), pointer :: grid
  PetscInt :: phase
  type(option_type) :: option
  
  !PetscInt :: phase, cntrl_var, ivar
  !type(connection_set_type), pointer :: connection_set
  !type(mphase_auxvar_type), pointer :: auxvars(:)

  PetscReal :: rate 

  PetscReal :: press_div_loc, press_div
  PetscReal :: conn_loc, conn_tot
  !PetscReal :: p_well_lc, p_well
  PetscInt :: iconn, local_id, ghosted_id, iph, ierr 
  PetscReal :: mob
 
  if(this%connection_set%num_connections > 0 ) then
   
    rate =  this%flow_condition%flow_well%rate%dataset%rarray(1)

    conn_tot = 0.0d0 
    ! divisor computation
    press_div_loc = 0.0d0
    do iconn = 1, this%connection_set%num_connections
      local_id = this%connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      mob = this%ConnMob(this%flow_auxvars(ZERO_INTEGER,ghosted_id)%mobility, &
                         phase)

      press_div_loc = press_div_loc + this%conn_factors(iconn) * mob * &
                  this%flow_auxvars(ZERO_INTEGER,ghosted_id)%den_kg(phase)

    end do

    call MPI_ALLREDUCE(press_div_loc, press_div, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM,this%comm, ierr)

    ! compute local well connection flux contribution
    conn_loc = 0.0d0 
    do iconn = 1, this%connection_set%num_connections      
      local_id = this%connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      mob = this%ConnMob(this%flow_auxvars(ZERO_INTEGER,ghosted_id)%mobility, &
                         phase)
      if(this%conn_status(iconn) == CONN_STATUS_OPEN ) then
        ! vol_rate > 0 for fluid entering the well  
         conn_loc = conn_loc + this%conn_factors(iconn) * mob * &
               this%flow_auxvars(ZERO_INTEGER,ghosted_id)%den_kg(phase) * &
               ( this%flow_auxvars(ZERO_INTEGER,ghosted_id)%pres(phase) - &
                this%conn_h(iconn) )
      end if
    end do ! end loop on well connections
    ! connection fluxes contibutions from well_segments in other ranks
    call MPI_ALLREDUCE(conn_loc, conn_tot, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM,this%comm, ierr)

  end if
 
  ! rate can be volume or mass rate depending on the control variable
  this%pw_ref = (conn_tot - rate) / press_div

#ifdef WELL_DEBUG
  write(*,"('FlowPressRefMRProd pw_ref = ',e10.4)") this%pw_ref 
#endif

end subroutine FlowPressRefMRProd

! ************************************************************************** !

subroutine FlowQPhase(this,grid,phase,option)
  !  
  ! Compute well volumetric rate for one phase given p_ref 
  ! IPR sign convention: vol_rate > 0 for fluid being produced 
  !
  ! Author: Paolo Orsini (OpenGoSim)  
  ! Date : 25/06/2015
  !

  use Grid_module
  use Option_module

  implicit none

  class(well_flow_type) :: this
  type(grid_type), pointer :: grid
  PetscInt :: phase
  type(option_type) :: option

  PetscReal :: vol_rate_lc, vol_rate 
  PetscInt :: iconn, local_id, ghosted_id, iph, ierr 
  PetscReal :: mob

  vol_rate = 0.0d0
  if(this%connection_set%num_connections > 0 ) then

    vol_rate_lc = 0.0d0 
    do iconn = 1, this%connection_set%num_connections      
      local_id = this%connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      mob = this%ConnMob(this%flow_auxvars(ZERO_INTEGER,ghosted_id)%mobility, &
                         phase)
      if(this%conn_status(iconn) == CONN_STATUS_OPEN ) then
        ! vol_rate > 0 for fluid entering the well  
        vol_rate_lc = vol_rate_lc + this%conn_factors(iconn) * mob * &
                  ( this%flow_auxvars(ZERO_INTEGER,ghosted_id)%pres(phase) - &                         
                    this%pw_ref - this%conn_h(iconn) )
      end if
    end do ! end loop on well connections
    ! vol_rate contibutions from well_segments in other ranks
    call MPI_ALLREDUCE(vol_rate_lc, vol_rate, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM,this%comm, ierr)

  end if

  this%q_fld(phase) = vol_rate

#ifdef WELL_DEBUG
  write(*,"('FlowQPhase pw_ref = ',e10.4)") this%pw_ref 
  write(*,"('FlowQPhase vol_rate = ',e10.4)") vol_rate
#endif

end subroutine FlowQPhase

!*****************************************************************************!

subroutine FlowMRPhase(this,grid,phase,option)
  !  
  ! Compute well mass rate for one phase given p_ref 
  ! IPR sign convention: mass_rate > 0 for fluid being produced 
  !
  ! Author: Paolo Orsini (OpenGoSim)  
  ! Date : 11/07/2015
  !
  use Grid_module
  use Option_module

  implicit none

  class(well_flow_type) :: this
  type(grid_type), pointer :: grid
  PetscInt :: phase
  type(option_type) :: option

  !type(connection_set_type), pointer :: connection_set

  PetscReal :: mass_rate_lc, mass_rate 
  PetscInt :: iconn, local_id, ghosted_id, iph, ierr 
  PetscReal :: mob

  mass_rate = 0.0d0
  if(this%connection_set%num_connections > 0 ) then

    mass_rate_lc = 0.0d0 
    do iconn = 1, this%connection_set%num_connections      
      local_id = this%connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      mob = this%ConnMob(this%flow_auxvars(ZERO_INTEGER,ghosted_id)%mobility, &
                         phase)
      if(this%conn_status(iconn) == CONN_STATUS_OPEN ) then
        ! mass_rate > 0 for fluid entering the well  
        mass_rate_lc = mass_rate_lc + this%conn_factors(iconn) * mob * &
                  this%flow_auxvars(ZERO_INTEGER,ghosted_id)%den_kg(phase) * &
                  (this%flow_auxvars(ZERO_INTEGER,ghosted_id)%pres(phase) - &
                   this%pw_ref - this%conn_h(iconn) )
      end if
    end do ! end loop on well connections
    ! mass_rate contibutions from well_segments in other ranks
    call MPI_ALLREDUCE(mass_rate_lc, mass_rate, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM,this%comm, ierr)

  end if

  this%mr_fld(phase) = mass_rate

end subroutine FlowMRPhase

!*****************************************************************************!

subroutine WellFlowLimitCheck(this,pass,option)
  ! 
  !
  ! Perform limit check for a water injector
  !
  ! Author: Paolo Orsini (OpenGoSim)  
  ! Date : 6/06/2016
  !

  use Option_module

  implicit none

  class(well_flow_type) :: this
  PetscBool :: pass
  type(option_type) :: option

  print *, "WellFlowLimitCheck must be extended"
  stop  

end subroutine WellFlowLimitCheck

!*****************************************************************************!

subroutine FlowHydroCorrUpdate(this,grid,option)
  !
  ! Updtae well hydrostatic correction for each well connection
  !
  ! Author: Paolo Orsini (OpenGoSim)  
  ! Date : 6/06/2016
  !

  use Grid_module
  use Option_module

  implicit none

  class(well_flow_type) :: this
  type(grid_type), pointer :: grid
  type(option_type) :: option

  !PetscReal, pointer :: ss_fluxes(:,:)
  PetscInt, pointer :: ord(:), l2w(:)
  PetscReal, pointer :: zcn(:) 
  PetscReal, pointer :: well_conn_h(:)
  PetscInt :: iconn, ierr
  PetscReal :: rho_up, rho_dn, p_up, p_dn, delta_z
  !for debugging
  PetscInt :: i_ph  

  PetscInt :: w_myrank

  if( this%connection_set%num_connections == 0 ) return

  ! add here call to compute all well connection densities and 
  ! and well connection (for z-ordered)
  ! note this operation is repeated in all processors sharing the well 
  ! it couldbe done by the control processors only and communicate results

  !create here a one
  if (this%hydrostatic_method == WELL_HYDROSTATIC_ITERATIVE) then

    call this%HydrostaticUpdate(grid,option)
  
    ! concatanate densities from different ranks
    call MPI_Allgatherv(this%conn_den_kg,this%connection_set%num_connections, &
                MPI_DOUBLE_PRECISION, this%well_conn_den_kg, this%w_rank_conn, &
                this%disp_rank_conn,MPI_DOUBLE_PRECISION, this%comm,ierr)  

    !for printing purposes only: assign dw_kg_ref avergare value 
    this%dw_kg_ref = 0.d0
    do iconn=1,this%well_num_conns
    
      do i_ph=1,option%nphase
        this%dw_kg_ref(i_ph) = this%dw_kg_ref(i_ph) + &
                               this%well_conn_den_kg(iconn)
      end do
    end do 
    do i_ph=1,option%nphase
      this%dw_kg_ref(i_ph) = this%dw_kg_ref(i_ph) / dble(this%well_num_conns)
    end do

    ord => this%w_conn_order
    l2w => this%conn_l2w      

    allocate(well_conn_h(this%well_num_conns))
    well_conn_h = 0.0d0
    ! concatanate hydrostatic corrections from different ranks for printing only
    call MPI_Allgatherv(this%conn_h,this%connection_set%num_connections, &
                MPI_DOUBLE_PRECISION, well_conn_h, this%w_rank_conn, &
                this%disp_rank_conn,MPI_DOUBLE_PRECISION, this%comm,ierr)  

#ifdef WELL_DEBUG
  print *,"After MPI_Allgatherv well_conn_h",well_conn_h(1:this%well_num_conns)
#endif

    !for printing purposes only: load ordered hydrostatic corrections
    do iconn=1,this%well_num_conns 
      this%well_conn_h_sorted(ord(iconn)) = well_conn_h(iconn)
    end do

    nullify(ord)
    nullify(l2w)
    deallocate(well_conn_h)
    nullify(well_conn_h)
  ! end of iterative method
  else if (this%hydrostatic_method == WELL_HYDROSTATIC_LINEAR) then 

!#ifdef WELL_DEBUG
!  print *,"After HUpdate MPI_Allgatherv" 
!#endif

    call this%ConnDenUpdate(grid,option) 

    ! concatanate densities from different ranks
    call MPI_Allgatherv(this%conn_den_kg,this%connection_set%num_connections, &
                MPI_DOUBLE_PRECISION, this%well_conn_den_kg, this%w_rank_conn, &
                this%disp_rank_conn,MPI_DOUBLE_PRECISION, this%comm,ierr)  

    ord => this%w_conn_order
    l2w => this%conn_l2w      
    zcn => this%w_conn_z !already sorted for acending z     

!#ifdef WELL_DEBUG
!  call MPI_Comm_rank(this%comm, w_myrank, ierr )
!  print *,"ConnHUpdate - Well_myrank", w_myrank
!#endif

!#ifdef WELL_DEBUG
!    print *,"ConnHUpdate - iwconn_ref", this%iwconn_ref
!#endif

    !Each process belonging to the well compute the hydrostatic correction
    !for all connections. Could use one rank only (e.g. control rank) 
    ! and communicate to the others

    !BEGINNING OF HYDROSTATIC CORRECTION COMPUTATION FOR THE ENTIRE WELL
    !hydro corrections
    !initialize cumulative well pressure to pw_ref
    p_up = this%pw_ref
    do iconn=1,this%well_num_conns
      if(iconn > this%iwconn_ref) then
      ! hydrostatic press corrections above z_ref
        delta_z = zcn(iconn) - zcn(iconn-1)
        rho_up = this%well_conn_den_kg(ord(iconn))
        rho_dn = this%well_conn_den_kg(ord(iconn - 1))
        p_up = p_up + 0.5d0*(rho_up+rho_dn) * option%gravity(Z_DIRECTION) * &
               delta_z
                                      ! this is a negative value 
        this%well_conn_h_sorted(iconn) = p_up - this%pw_ref   
      end if 
    end do  

    !initialize cumulative well pressure to pw_ref
    p_dn = this%pw_ref   
    do iconn=this%well_num_conns,1,-1
      if(iconn < this%iwconn_ref) then
      ! hydrostatic press corrections below z_ref
        delta_z = zcn(iconn + 1) - zcn(iconn)
        rho_up = this%well_conn_den_kg(ord(iconn + 1))
        rho_dn = this%well_conn_den_kg(ord(iconn))
        p_dn = p_dn - 0.5d0*(rho_up+rho_dn) * option%gravity(Z_DIRECTION) * &
               delta_z 
                                       ! this is a positive value
        this%well_conn_h_sorted(iconn) = p_dn - this%pw_ref 
      end if 
    end do
    !END OF HYDROSTATIC CORRECTION COMPUTATION FOR THE ENTIRE WELL

    ! load hydrostatic correction into the the well segment belonging 
    ! to the current well rank
    do iconn=1,this%connection_set%num_connections
       this%conn_h(iconn) = this%well_conn_h_sorted( ord(l2w(iconn)) )
    end do

#ifdef WELL_DEBUG
    do iconn=1,this%connection_set%num_connections
       print *, "conn i/h = ", iconn, this%conn_h(iconn)
    end do
#endif

    nullify(ord)
    nullify(l2w)
    nullify(zcn)
  
  end if !End WELL_HYDROSTATIC_LINEAR method

end subroutine FlowHydroCorrUpdate

!*****************************************************************************!

subroutine FlowHydrostaticUpdate(this,grid,option)
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

  class(well_flow_type) :: this
  type(grid_type), pointer :: grid
  type(option_type) :: option

  print *, "FlowHydrostaticUpdate must be extended"
  stop  


end subroutine FlowHydrostaticUpdate

!*****************************************************************************!

subroutine WellFlowConnDenUpdate(this,grid,option)
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

  class(well_flow_type) :: this
  type(grid_type), pointer :: grid 
  type(option_type) :: option

  PetscReal, pointer :: ss_fluxes(:,:)
  PetscReal :: q_sum
  PetscReal :: q_ph(option%nphase)
  PetscReal :: den_ph(option%nphase)
  PetscReal :: den_ph_w(option%nphase)
  PetscInt :: i_ph

  PetscReal :: q_sum_lc, q_sum_well 
  PetscReal :: den_q_lc, den_q_well
  PetscReal :: den_sum_well

  PetscInt :: iconn, local_id, ghosted_id, ierr

  q_sum_lc = 0.0d0
  q_sum_well = 0.0d0
  den_q_lc = 0.0d0
  den_q_well = 0.0d0
  den_sum_well = 0.0d0

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
      den_ph(i_ph) = this%flow_auxvars(ZERO_INTEGER,ghosted_id)%den_kg(i_ph)
    end do
    if( q_sum > wfloweps ) then
      do i_ph = 1,option%nphase
        den_ph_w(i_ph) = q_ph(i_ph) / q_sum  
      end do
    else ! OK for nill fluxes 
         ! not accurate for back flow (should get warning in Res computation)
      do i_ph = 1,option%nphase
        den_ph_w(i_ph) = this%flow_auxvars(ZERO_INTEGER,ghosted_id)%sat(i_ph)
      end do
    end if   
    !compute fluid averga density 
    this%conn_den_kg(iconn) = 0.0d0
    do i_ph = 1,option%nphase
      this%conn_den_kg(iconn) = this%conn_den_kg(iconn) + &
                                den_ph(i_ph) * den_ph_w(i_ph)
    end do
     
    q_sum_lc = q_sum_lc + q_sum
    den_q_lc = den_q_lc + this%conn_den_kg(iconn) * q_sum
    den_sum_well = den_sum_well + this%conn_den_kg(iconn)

  end do

  call MPI_ALLREDUCE(q_sum_lc, q_sum_well, 1, MPI_DOUBLE_PRECISION, &
                     MPI_SUM,this%comm, ierr)

  call MPI_ALLREDUCE(den_q_lc, den_q_well, 1, MPI_DOUBLE_PRECISION, &
                     MPI_SUM,this%comm, ierr)

  !write(*,*) "rank = ", option%myrank, "WellFlowConnDenUpdate, den_q_well, q_sum_well", den_q_well, q_sum_well
  !same well average density for all phase - assuming perfect mix
  do i_ph = 1,option%nphase
    if ( q_sum_well > wfloweps ) then
      this%dw_kg_ref(i_ph) = den_q_well / q_sum_well
    else
      this%dw_kg_ref(i_ph) = den_sum_well / &
                             dble(this%well_num_conns)
    end if 
  end do

end subroutine WellFlowConnDenUpdate

!*****************************************************************************!

subroutine WellFlowInitDensity(this,grid,option)
  !
  ! Author: Paolo Orsini (OpenGoSim)  
  ! Date : 8/04/2016
  !

  use Grid_module
  use Option_module

  implicit none

  class(well_flow_type) :: this
  type(grid_type), pointer :: grid  
  type(option_type) :: option

  print *, "Well WellFlowInitDensity must be extended"
  stop

end subroutine WellFlowInitDensity

!*****************************************************************************!

function WellFlowConnMob(this,mobility,iphase)

  ! mobilty computation for producers
  ! to be overwritten (or replaced) for injectors
  !
  ! Author: Paolo Orsini (OpenGoSim)  
  ! Date : 1/07/2015
  !

  implicit none

  class(well_flow_type) :: this 
  PetscInt :: iphase  !not required for this WellWatInjConnMob extension
  PetscReal :: mobility(:)

  PetscReal :: WellFlowConnMob

  WellFlowConnMob = mobility(iphase)

  !print *, "WellFlowConnMob must be extended"
  !stop

end function WellFlowConnMob

!*****************************************************************************!
subroutine FlowDataOutput(this,grid,src_name,option)
  !
  ! Write well pressure and perforated grid lock profile
  ! Overwrites previous file - currently for debugging
  ! TO DO - should add control at which time step to print the profiles 
  !
  ! Author: Paolo Orsini (OpenGoSim)  
  ! Date : 6/20/2016

  use Grid_module
  use Option_module
  
  implicit none

  class(well_flow_type) :: this
  type(grid_type), pointer :: grid
  character(len=MAXWORDLENGTH) :: src_name
  type(option_type) :: option

  PetscInt, pointer :: ord(:), l2w(:)
  PetscReal, pointer :: zcn(:) 
  PetscReal, pointer :: lc_perf_block_press(:,:)
  PetscReal, pointer :: perf_block_press(:,:)
  PetscInt :: iconn, local_id, ghosted_id, i_ph
  character(len=MAXSTRINGLENGTH) :: wfile_name
  PetscMPIInt :: cur_w_myrank
  PetscInt :: ierr

  ord => this%w_conn_order
  l2w => this%conn_l2w      
  zcn => this%w_conn_z !already sorted for acending z     

  allocate(perf_block_press(option%nphase,this%well_num_conns))
  perf_block_press = 0.0d0

  allocate(lc_perf_block_press(option%nphase,this%connection_set%num_connections))
  perf_block_press = 0.0d0

  !load local connections
  do iconn = 1, this%connection_set%num_connections      
    local_id = this%connection_set%id_dn(iconn)
    ghosted_id = grid%nL2G(local_id)
    do i_ph =1,option%nphase
      !perf_block_press(i_ph,ord(l2w(iconn))) = & 
        lc_perf_block_press(i_ph,iconn) = &
        this%flow_auxvars(ZERO_INTEGER,ghosted_id)%pres(i_ph)
    end do
  end do 

#ifdef WELL_DEBUG
  call MPI_Comm_rank(this%comm, cur_w_myrank, ierr ) 
  print *, "Before MPI_Allgatherv well myrank =",cur_w_myrank, &
           "myrank =", option%myrank, "perf_block_press", &
           perf_block_press(1,1:this%well_num_conns)
#endif
 
  do i_ph=1,option%nphase
    !MPI_Allgatherv because each well segment can have a different number of conns
    call MPI_Allgatherv(lc_perf_block_press(i_ph,:), &
               this%connection_set%num_connections, &
               MPI_DOUBLE_PRECISION, perf_block_press(i_ph,:), &
               this%w_rank_conn, this%disp_rank_conn,MPI_DOUBLE_PRECISION, &
               this%comm, ierr)
  end do

  !make sure all well ranks have filled in their part of perf_block_press
  !call MPI_Barrier( this%comm,ierr)

#ifdef WELL_DEBUG
  call MPI_Comm_rank(this%comm, cur_w_myrank, ierr ) 
  print *, "After MPI_Allgatherv well myrank =",cur_w_myrank, &
           "myrank =", option%myrank, "perf_block_press", &
           perf_block_press(1,1:this%well_num_conns)
#endif


  if( this%connection_set%num_connections > 0 ) then
    call MPI_Comm_rank(this%comm, cur_w_myrank, ierr )  
    if(this%cntr_rank == cur_w_myrank ) then
      wfile_name = trim(option%global_prefix) // "_" // &
                   trim(src_name) // ".dat" 
      open(unit=IUNIT_TEMP,file=wfile_name)
      write(IUNIT_TEMP,*) "z  press_well  press_block_ph1, press_block_ph2 "
      do iconn=1,this%well_num_conns
        write(IUNIT_TEMP,"(4(E16.10,1x))") zcn(iconn), &
                this%pw_ref + this%well_conn_h_sorted(iconn), &
                perf_block_press(1:option%nphase,ord(iconn)) 
      end do
      close(IUNIT_TEMP)
    end if
  end if 

  nullify(ord)
  nullify(l2w)
  nullify(zcn)

  deallocate(perf_block_press)
  nullify(perf_block_press)
  deallocate(lc_perf_block_press)
  nullify(lc_perf_block_press)

end subroutine FlowDataOutput
!*****************************************************************************!
subroutine FlowWellStrip(well)
  !
  ! Strip well_flow and all its parent members
  !
  ! Author: Paolo Orsini (OpenGoSim)  
  ! Date : 6/14/2016
  !

  use Utility_module, only : DeallocateArray 

  implicit none

  class(well_flow_type) :: well

  call DeallocateArray(well%dw_kg_ref)
  call DeallocateArray(well%q_fld)
  call DeallocateArray(well%mr_fld)
  call DeallocateArray(well%conn_h)
  call DeallocateArray(well%conn_den_kg)
  call DeallocateArray(well%well_conn_den_kg)
  call DeallocateArray(well%well_conn_h_sorted)

  call DeallocateArray(well%well_fine_grid_pres)
  call DeallocateArray(well%well_fine_grid_den_kg)

  !these are pointer only
  nullify(well%ss_flow_vol_fluxes) 
  nullify(well%flow_auxvars)
  nullify(well%flow_condition)

  !this will strip all its parents
  call BaseWellStrip(well)


end subroutine FlowWellStrip

!*****************************************************************************!

#endif  
end module Well_Flow_class
!end of WELL_CLASS
