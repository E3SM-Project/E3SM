module Surface_TH_module



#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use Surface_Global_Aux_module
  use Surface_TH_Aux_module
  
  use PFLOTRAN_Constants_module
  use Utility_module, only : Equal
  
  implicit none
  
  private

! Cutoff parameters
  PetscReal, parameter :: eps       = 1.D-12
  PetscReal, parameter :: perturbation_tolerance = 1.d-6

  public SurfaceTHSetup, &
         SurfaceTHRHSFunction, &
         SurfaceTHIFunction, &
         SurfaceTHComputeMaxDt, &
         SurfaceTHUpdateAuxVars, &
         SurfaceTHUpdateSolution, &
         SurfaceTHUpdateTemperature, &
         SurfaceTHUpdateSurfState, &
         SurfaceTHImplicitAtmForcing, &
         SurfaceTHDestroy

contains

! ************************************************************************** !

subroutine SurfaceTHSetup(surf_realization)
  ! 
  ! This routine sets up surface_TH_type
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 02/28/13
  ! 

  use Realization_Surface_class
  use Patch_module
  use Option_module
  use Grid_module
  use Region_module
  use Coupler_module
  use Connection_module
  use Fluid_module
  use Output_Aux_module
 
  implicit none
  
  class(realization_surface_type) :: surf_realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: boundary_condition
  type(Surface_TH_auxvar_type), pointer :: Surf_TH_auxvars(:)
  type(Surface_TH_auxvar_type), pointer :: Surf_TH_auxvars_bc(:)
  type(Surface_TH_auxvar_type), pointer :: Surf_TH_auxvars_ss(:)
  type(fluid_property_type), pointer :: cur_fluid_property
  type(coupler_type), pointer :: initial_condition
  type(output_variable_list_type), pointer :: list
  PetscReal :: area_per_vol

  PetscInt :: ghosted_id, iconn, sum_connection
  PetscInt :: i, iphase
  
  option => surf_realization%option
  patch => surf_realization%patch
  grid => patch%grid
    
  patch%surf_aux%SurfaceTH => SurfaceTHAuxCreate(option)

  ! allocate auxvar data structures for all grid cells
  allocate(Surf_TH_auxvars(grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    call SurfaceTHAuxVarInit(Surf_TH_auxvars(ghosted_id),option)
  enddo

  patch%surf_aux%SurfaceTH%auxvars => Surf_TH_auxvars
  patch%surf_aux%SurfaceTH%num_aux = grid%ngmax

  ! count the number of boundary connections and allocate
  ! auxvar data structures for them
  boundary_condition => patch%boundary_condition_list%first

  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    sum_connection = sum_connection + &
                     boundary_condition%connection_set%num_connections
    boundary_condition => boundary_condition%next
  enddo

  if (sum_connection > 0) then 
    allocate(Surf_TH_auxvars_bc(sum_connection))
    do iconn = 1, sum_connection
      call SurfaceTHAuxVarInit(Surf_TH_auxvars_bc(iconn),option)
    enddo
    patch%surf_aux%SurfaceTH%auxvars_bc => Surf_TH_auxvars_bc
  endif
  patch%surf_aux%SurfaceTH%num_aux_bc = sum_connection

  ! Create aux vars for source/sink
  sum_connection = CouplerGetNumConnectionsInList(patch%source_sink_list)
  if (sum_connection > 0) then
    allocate(Surf_TH_auxvars_ss(sum_connection))
    do iconn = 1, sum_connection
      call SurfaceTHAuxVarInit(Surf_TH_auxvars_ss(iconn),option)
    enddo
    patch%surf_aux%SurfaceTH%auxvars_ss => Surf_TH_auxvars_ss
  endif
  patch%surf_aux%SurfaceTH%num_aux_ss = sum_connection

  list => surf_realization%output_option%output_snap_variable_list
  call SurfaceTHSetPlotVariables(list)
  list => surf_realization%output_option%output_obs_variable_list
  call SurfaceTHSetPlotVariables(list)

end subroutine SurfaceTHSetup

! ************************************************************************** !

subroutine SurfaceTHSetPlotVariables(list)
  ! 
  ! This routine adds default variables to be printed to list
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 02/28/13
  ! 
  
  use Realization_Surface_class
  use Output_Aux_module
  use Variables_module
    
  implicit none
  
  type(output_variable_list_type), pointer :: list

  character(len=MAXWORDLENGTH) :: name, units
  
  if (associated(list%first)) then
    return
  endif

  name = 'H'
  units = 'm'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               SURFACE_LIQUID_HEAD)

  name = 'Temperature'
  units = 'C'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               SURFACE_LIQUID_TEMPERATURE)

  name = 'Material ID'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_DISCRETE,units, &
                               MATERIAL_ID)
  
end subroutine SurfaceTHSetPlotVariables

! ************************************************************************** !

subroutine SurfaceTHRHSFunction(ts,t,xx,ff,surf_realization,ierr)
  ! 
  ! This routine provides the function evaluation for PETSc TSSolve()
  ! Author: Gautam Bisht, LBNL
  ! 

#include <petsc/finclude/petscts.h>
  use petscts
  use EOS_Water_module
  use Connection_module
  use Realization_Surface_class
  use Discretization_module
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module  
  use Surface_Field_module
  use Debug_module
  use Surface_TH_Aux_module
  use Surface_Global_Aux_module

  implicit none
  
  TS :: ts
  PetscReal :: t
  Vec :: xx
  Vec :: ff
  class(realization_surface_type) :: surf_realization
  PetscErrorCode :: ierr

  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(surface_field_type), pointer :: surf_field
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set

  type(Surface_TH_auxvar_type), pointer :: surf_auxvars(:)
  type(Surface_TH_auxvar_type), pointer :: surf_auxvars_bc(:)
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars(:)
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars_bc(:)
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars_ss(:)

  PetscInt :: local_id_up, local_id_dn, local_id
  PetscInt :: ghosted_id_up, ghosted_id_dn, ghosted_id
  PetscInt :: iconn
  PetscInt :: sum_connection
  PetscInt :: istart, iend

  PetscReal :: dx, dy, dz
  PetscReal :: dist
  PetscReal :: vel
  PetscReal :: slope, slope_dn
  PetscReal :: rho          ! density      [kg/m^3]
  PetscReal :: hw_up, hw_dn ! water height [m]
  PetscReal :: Res(surf_realization%option%nflowdof), v_darcy
  PetscReal :: qsrc, qsrc_flow
  PetscReal :: esrc
  PetscReal :: den
  PetscReal :: dum1

  PetscViewer :: viewer
  character(len=MAXSTRINGLENGTH) :: string,string2

  PetscReal, pointer :: ff_p(:), mannings_loc_p(:),area_p(:)
  PetscReal, pointer :: xc(:),yc(:),zc(:)

  patch => surf_realization%patch
  grid => patch%grid
  option => surf_realization%option
  surf_field => surf_realization%surf_field

  surf_auxvars => patch%surf_aux%SurfaceTH%auxvars
  surf_auxvars_bc => patch%surf_aux%SurfaceTH%auxvars_bc
  surf_global_auxvars => patch%surf_aux%SurfaceGlobal%auxvars
  surf_global_auxvars_bc => patch%surf_aux%SurfaceGlobal%auxvars_bc
  surf_global_auxvars_ss => patch%surf_aux%SurfaceGlobal%auxvars_ss

  surf_realization%iter_count = surf_realization%iter_count+1
  if (surf_realization%iter_count < 10) then
    write(string2,'("00",i1)') surf_realization%iter_count
  else if (surf_realization%iter_count < 100) then
    write(string2,'("0",i2)') surf_realization%iter_count
  else if (surf_realization%iter_count < 1000) then
    write(string2,'(i3)') surf_realization%iter_count
  else if (surf_realization%iter_count < 10000) then
    write(string2,'(i4)') surf_realization%iter_count
  endif 

  ! First, update the solution vector
  call DiscretizationGlobalToLocal(surf_realization%discretization, &
                                   xx,surf_field%flow_xx_loc,NFLOWDOF)

  ! Then, update the aux vars
  ! RTM: This includes calculation of the accumulation terms, correct?
  call SurfaceTHUpdateTemperature(surf_realization)
  call SurfaceTHUpdateAuxVars(surf_realization)
  ! override flags since they will soon be out of date  
  patch%surf_aux%SurfaceTH%auxvars_up_to_date = PETSC_FALSE

  call VecGetArrayF90(ff,ff_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(surf_field%mannings_loc,mannings_loc_p,  &
                      ierr);CHKERRQ(ierr)
  call VecGetArrayF90(surf_field%area,area_p,ierr);CHKERRQ(ierr)

  ff_p = 0.d0
  Res  = 0.d0

  xc => surf_realization%discretization%grid%x
  yc => surf_realization%discretization%grid%y
  zc => surf_realization%discretization%grid%z

  ! Interior Flux Terms -----------------------------------
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0  
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1

      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)

      local_id_up = grid%nG2L(ghosted_id_up)
      local_id_dn = grid%nG2L(ghosted_id_dn)
      
      dx = xc(ghosted_id_dn) - xc(ghosted_id_up)
      dy = yc(ghosted_id_dn) - yc(ghosted_id_up)
      dz = zc(ghosted_id_dn) - zc(ghosted_id_up)
      dist = sqrt(dx*dx + dy*dy + dz*dz)
      slope = dz/dist
      
      call SurfaceTHFlux(surf_auxvars(ghosted_id_up), &
                         surf_global_auxvars(ghosted_id_up), &
                         zc(ghosted_id_up), &
                         mannings_loc_p(ghosted_id_up), &
                         surf_auxvars(ghosted_id_dn), &
                         surf_global_auxvars(ghosted_id_dn), &
                         zc(ghosted_id_dn), &
                         mannings_loc_p(ghosted_id_dn), &
                         dist, cur_connection_set%area(iconn), &
                         option,vel,dum1,Res)

      patch%internal_velocities(1,sum_connection) = vel
      patch%internal_flow_fluxes(:,sum_connection) = Res(:)

      if (local_id_up>0) then
        iend = local_id_up*option%nflowdof
        istart = iend-option%nflowdof+1
        ff_p(istart:iend) = ff_p(istart:iend) - Res(:)/area_p(local_id_up)
      endif
         
      if (local_id_dn>0) then
        iend = local_id_dn*option%nflowdof
        istart = iend-option%nflowdof+1
        ff_p(istart:iend) = ff_p(istart:iend) + Res(:)/area_p(local_id_dn)
      endif

    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  ! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection_set
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      local_id_dn = cur_connection_set%id_dn(iconn)
      ghosted_id_dn = grid%nL2G(local_id_dn)
  
      dx = xc(ghosted_id_dn) - cur_connection_set%intercp(1,iconn)
      dy = yc(ghosted_id_dn) - cur_connection_set%intercp(2,iconn)
      dz = zc(ghosted_id_dn) - cur_connection_set%intercp(3,iconn)
      dist = sqrt(dx*dx + dy*dy + dz*dz)
      slope_dn = dz/sqrt(dx*dx + dy*dy + dz*dz)

      call SurfaceTHBCFlux(boundary_condition%flow_condition%itype, &
                         boundary_condition%flow_aux_real_var(:,iconn), &
                         surf_auxvars_bc(sum_connection), &
                         surf_global_auxvars_bc(sum_connection), &
                         surf_auxvars(ghosted_id_dn), &
                         surf_global_auxvars(ghosted_id_dn), &
                         slope_dn, &
                         mannings_loc_p(ghosted_id_dn), &
                         dist, &
                         cur_connection_set%area(iconn), &
                         option,vel,dum1,Res)

      patch%boundary_velocities(1,sum_connection) = vel
      patch%boundary_flow_fluxes(:,sum_connection) = Res(:)
      
      iend = local_id_dn*option%nflowdof
      istart = iend-option%nflowdof+1
      ff_p(istart:iend) = ff_p(istart:iend) + Res(:)/area_p(local_id_dn)
    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sink_list%first
  sum_connection = 0
  do
    if (.not.associated(source_sink)) exit
    
    if (source_sink%flow_condition%rate%itype/=HET_VOL_RATE_SS.and. &
       source_sink%flow_condition%rate%itype/=HET_MASS_RATE_SS) &
    qsrc_flow = source_sink%flow_condition%rate%dataset%rarray(1)
      
    if (source_sink%flow_condition%rate%itype == ENERGY_RATE_SS) &
      esrc = source_sink%flow_condition%energy_rate%dataset%rarray(1)

    cur_connection_set => source_sink%connection_set
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle

      select case(source_sink%flow_condition%rate%itype)
        case(VOLUMETRIC_RATE_SS)  ! assume local density for now
          ! qsrc = m^3/sec
          qsrc = qsrc_flow*area_p(local_id)
        case(HET_VOL_RATE_SS)
          ! qsrc = m^3/sec
          qsrc = source_sink%flow_aux_real_var(ONE_INTEGER,iconn)*area_p(local_id)
        case default
          option%io_buffer = 'Source/Sink flow condition type not recognized'
          call printErrMsg(option)
      end select
      
      esrc = 0.d0
      select case(source_sink%flow_condition%itype(TH_TEMPERATURE_DOF))
        case (ENERGY_RATE_SS)
          esrc = source_sink%flow_condition%energy_rate%dataset%rarray(1)
        case (HET_ENERGY_RATE_SS)
          esrc = source_sink%flow_aux_real_var(TWO_INTEGER,iconn)
      end select

      iend = local_id*option%nflowdof
      istart = iend-option%nflowdof+1

      ff_p(istart) = ff_p(istart) + qsrc/area_p(local_id)
      ! RTM: TODO: What should the density term and specific heat capactiy be
      ! in the freezing case?
      ! I think using the weighted average of liquid and ice densities and Cwi 
      ! is correct here, but I should check.
      ff_p(iend) = ff_p(iend) + esrc + &
                    surf_global_auxvars_ss(sum_connection)%den_kg(1)* &
                    (surf_global_auxvars_ss(sum_connection)%temp + 273.15d0)* &
                    surf_auxvars(local_id)%Cwi* &
                    qsrc/area_p(local_id)
    enddo
    source_sink => source_sink%next
  enddo

  call VecRestoreArrayF90(ff,ff_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(surf_field%mannings_loc,mannings_loc_p, &
                          ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(surf_field%area,area_p,ierr);CHKERRQ(ierr)

  if (surf_realization%debug%vecview_solution) then
    string = 'Surf_xx_' // trim(adjustl(string2)) // '.bin'
    call PetscViewerBinaryOpen(surf_realization%option%mycomm,string, &
                              FILE_MODE_WRITE,viewer,ierr);CHKERRQ(ierr)
    call VecView(xx,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)

    string = 'Surf_ff_' // trim(adjustl(string2)) // '.bin'
    call PetscViewerBinaryOpen(surf_realization%option%mycomm,string, &
                              FILE_MODE_WRITE,viewer,ierr);CHKERRQ(ierr)
    call VecView(ff,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

end subroutine SurfaceTHRHSFunction

! ************************************************************************** !

subroutine SurfaceTHIFunction(ts,t,xx,xxdot,ff,surf_realization,ierr)
  ! 
  ! This routine provides the implicit function evaluation for PETSc TSSolve()
  ! Author: Nathan Collier, ORNL
  ! 

#include "petsc/finclude/petscts.h"
  use petscts
  use EOS_Water_module
  use Connection_module
  use Realization_Surface_class
  use Discretization_module
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module  
  use Surface_Field_module
  use Debug_module
  use Surface_TH_Aux_module
  use Surface_Global_Aux_module

  implicit none
  
  TS :: ts
  PetscReal :: t
  Vec :: xx,xxdot
  Vec :: ff
  class(realization_surface_type) :: surf_realization
  PetscErrorCode :: ierr

  ! Our equations are in the form: 
  !    xxdot = RHS(xx)
  ! or in residual form:
  !    ff = xxdot - RHS(xx)

  ! First we call RHS function: ff = RHS(xx)
  call SurfaceTHRHSFunction(ts,t,xx,ff,surf_realization,ierr);CHKERRQ(ierr)
  ! negate: RHS(xx) = -RHS(xx)
  call VecScale(ff,-1.d0,ierr);CHKERRQ(ierr)
  ! and finally: ff += xxdot
  call VecAYPX(ff,1.d0,xxdot,ierr);CHKERRQ(ierr)
  
end subroutine SurfaceTHIFunction

! ************************************************************************** !

subroutine SurfaceTHComputeMaxDt(surf_realization,max_allowable_dt)
  ! 
  ! This routine maximum allowable 'dt' for explicit time scheme.
  ! Author: Gautam Bisht, LBNL
  ! 

  use EOS_Water_module
  use Connection_module
  use Realization_Surface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module  
  use Surface_Field_module
  use Debug_module
  use Surface_TH_Aux_module
  use Surface_Global_Aux_module

  implicit none
  
  class(realization_surface_type) :: surf_realization
  PetscErrorCode :: ierr

  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(surface_field_type), pointer :: surf_field
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set

  type(Surface_TH_auxvar_type), pointer :: surf_auxvars(:)
  type(Surface_TH_auxvar_type), pointer :: surf_auxvars_bc(:)
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars(:)
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars_bc(:)

  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  PetscInt :: iconn
  PetscInt :: sum_connection
#ifdef SURFACE_TH_DEBUG
  PetscInt :: max_connection,max_iconn
#endif

  PetscReal :: dx, dy, dz
  PetscReal :: dist
  PetscReal :: vel
  PetscReal :: slope, slope_dn
  PetscReal :: hw_up, hw_dn ! water height [m]
  PetscReal :: Res(surf_realization%option%nflowdof), v_darcy
  PetscReal :: max_allowable_dt
  PetscReal :: dt

  PetscReal, pointer :: mannings_loc_p(:),area_p(:)
  PetscReal, pointer :: xc(:),yc(:),zc(:)

  patch => surf_realization%patch
  grid => patch%grid
  option => surf_realization%option
  surf_field => surf_realization%surf_field

  surf_auxvars => patch%surf_aux%SurfaceTH%auxvars
  surf_auxvars_bc => patch%surf_aux%SurfaceTH%auxvars_bc
  surf_global_auxvars => patch%surf_aux%SurfaceGlobal%auxvars
  surf_global_auxvars_bc => patch%surf_aux%SurfaceGlobal%auxvars_bc

  call VecGetArrayF90(surf_field%mannings_loc,mannings_loc_p,  &
                      ierr);CHKERRQ(ierr)
  call VecGetArrayF90(surf_field%area,area_p,ierr);CHKERRQ(ierr)

  Res  = 0.d0
  max_allowable_dt = 1.d10
  vel = 0.d0

  xc => surf_realization%discretization%grid%x
  yc => surf_realization%discretization%grid%y
  zc => surf_realization%discretization%grid%z

  ! Interior Flux Terms -----------------------------------
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0  
#ifdef SURFACE_TH_DEBUG
  max_connection = -1
  max_iconn      = -1
#endif
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1

      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)

      local_id_up = grid%nG2L(ghosted_id_up)
      local_id_dn = grid%nG2L(ghosted_id_dn)
      
      dx = xc(ghosted_id_dn) - xc(ghosted_id_up)
      dy = yc(ghosted_id_dn) - yc(ghosted_id_up)
      dz = zc(ghosted_id_dn) - zc(ghosted_id_up)
      dist = sqrt(dx*dx + dy*dy + dz*dz)
      slope = dz/dist
      
      call SurfaceTHFlux(surf_auxvars(ghosted_id_up), &
                         surf_global_auxvars(ghosted_id_up), &
                         zc(ghosted_id_up), &
                         mannings_loc_p(ghosted_id_up), &
                         surf_auxvars(ghosted_id_dn), &
                         surf_global_auxvars(ghosted_id_dn), &
                         zc(ghosted_id_dn), &
                         mannings_loc_p(ghosted_id_dn), &
                         dist, cur_connection_set%area(iconn), &
                         option,vel,dt,Res)

      patch%internal_velocities(1,sum_connection) = vel
      patch%internal_flow_fluxes(:,sum_connection) = Res(:)

#ifdef SURFACE_TH_DEBUG
      if (dt < max_allowable_dt) then
        max_connection = sum_connection
        max_iconn      = iconn
      endif
#endif
      max_allowable_dt = min(max_allowable_dt, dt)

    enddo
    cur_connection_set => cur_connection_set%next
  enddo

#ifdef SURFACE_TH_DEBUG
  if (max_allowable_dt < 1.d-1) then
    cur_connection_set => connection_set_list%first
    ghosted_id_up = cur_connection_set%id_up(max_iconn)
    ghosted_id_dn = cur_connection_set%id_dn(max_iconn)
    local_id_up = grid%nG2L(ghosted_id_up)
    local_id_dn = grid%nG2L(ghosted_id_dn)
    dx = xc(ghosted_id_dn) - xc(ghosted_id_up)
    dy = yc(ghosted_id_dn) - yc(ghosted_id_up)
    dz = zc(ghosted_id_dn) - zc(ghosted_id_up)
    dist = sqrt(dx*dx + dy*dy + dz*dz)
    slope = dz/dist
    print *,"--------------------------"
    print *,"max_allowable_dt:",max_allowable_dt
    print *,"connection:",max_iconn
    print *,"(dx,dy,dz):",dx,dy,dz
    print *,"dist:      ",dist
    print *,"slope:     ",slope
    print *,"flux:      ",patch%internal_velocities(1,max_connection)
    print *,"dt:        ",dist/abs(patch%internal_velocities(1,max_connection))/3.0d0
    print *,"up info:",ghosted_id_up
    print *,"  istate:",surf_global_auxvars(ghosted_id_up)%istate
    print *,"  head:  ",surf_global_auxvars(ghosted_id_up)%head(1)
    print *,"  zc:    ",zc(ghosted_id_up)
    print *,"  temp:  ",surf_global_auxvars(ghosted_id_up)%temp
    print *,"  is_dry:",surf_global_auxvars(ghosted_id_up)%is_dry
    print *,"dn info:",ghosted_id_dn
    print *,"  istate:",surf_global_auxvars(ghosted_id_dn)%istate
    print *,"  head:  ",surf_global_auxvars(ghosted_id_dn)%head(1)
    print *,"  zc:    ",zc(ghosted_id_dn)
    print *,"  temp:  ",surf_global_auxvars(ghosted_id_dn)%temp
    print *,"  is_dry:",surf_global_auxvars(ghosted_id_dn)%is_dry
  endif  
#endif

  ! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection_set

    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      local_id_dn = cur_connection_set%id_dn(iconn)
      ghosted_id_dn = grid%nL2G(local_id_dn)
  
      dx = xc(ghosted_id_dn) - cur_connection_set%intercp(1,iconn)
      dy = yc(ghosted_id_dn) - cur_connection_set%intercp(2,iconn)
      dz = zc(ghosted_id_dn) - cur_connection_set%intercp(3,iconn)
      dist = sqrt(dx*dx + dy*dy + dz*dz)
      slope_dn = dz/sqrt(dx*dx + dy*dy + dz*dz)

      call SurfaceTHBCFlux(boundary_condition%flow_condition%itype, &
                         boundary_condition%flow_aux_real_var(:,iconn), &
                         surf_auxvars_bc(sum_connection), &
                         surf_global_auxvars_bc(sum_connection), &
                         surf_auxvars(ghosted_id_dn), &
                         surf_global_auxvars(ghosted_id_dn), &
                         slope_dn, &
                         mannings_loc_p(ghosted_id_dn), &
                         dist, &
                         cur_connection_set%area(iconn), &
                         option,vel,dt,Res)

      patch%boundary_velocities(1,sum_connection) = vel
      patch%boundary_flow_fluxes(:,sum_connection) = Res(:)

      max_allowable_dt = min(max_allowable_dt, dt)
    enddo
    boundary_condition => boundary_condition%next
  enddo
  
  call VecRestoreArrayF90(surf_field%mannings_loc,mannings_loc_p, &
                          ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(surf_field%area,area_p,ierr);CHKERRQ(ierr)
  
  if (max_allowable_dt < 0.d0) then
    write(option%io_buffer, &
          '("surface_th.F90: SurfaceTHComputeMaxDt --> negative max_allowable_dt!",es15.7)') &
          max_allowable_dt
    call printErrMsg(option)     
  endif

end subroutine SurfaceTHComputeMaxDt

! ************************************************************************** !

subroutine SurfaceTHFlux(surf_auxvar_up, &
                         surf_global_auxvar_up, &
                         zc_up, &
                         mannings_up, &
                         surf_auxvar_dn, &
                         surf_global_auxvar_dn, &
                         zc_dn, &
                         mannings_dn, &
                         dist, &
                         length, &
                         option, &
                         vel, &
                         dt_max, &
                         Res)
  ! 
  ! This routine computes the internal flux term for under
  ! diffusion-wave assumption.
  ! 
  ! Author: Gautam Bisht, LBL
  ! Date: 08/03/12
  ! 

  use Surface_TH_Aux_module
  use Surface_Global_Aux_module
  use Option_module
  use PFLOTRAN_Constants_module, only : MIN_SURFACE_WATER_HEIGHT

  implicit none

  type(option_type) :: option
  type(Surface_TH_auxvar_type) :: surf_auxvar_up
  type(Surface_TH_auxvar_type) :: surf_auxvar_dn
  type(surface_global_auxvar_type) :: surf_global_auxvar_up
  type(surface_global_auxvar_type) :: surf_global_auxvar_dn
  PetscReal :: zc_up, zc_dn
  PetscReal :: mannings_up, mannings_dn

  PetscReal :: head_up, head_dn
  PetscReal :: dist, length
  PetscReal :: vel                      ! [m/s]
  PetscReal :: dt_max
  PetscReal :: Res(1:option%nflowdof)   ! [m^3/s]
  
  PetscReal :: hw_half
  PetscReal :: hw_liq_half
  PetscReal :: mannings_half
  PetscReal :: unfrozen_fraction_half
  PetscReal :: dhead
  PetscReal :: den_aveg
  PetscReal :: temp_half
  PetscReal :: dtemp
  PetscReal :: Cw
  PetscReal :: k_therm
  PetscReal :: dt

  ! Initialize
  dt_max  = PETSC_MAX_REAL

  ! We upwind Manning's coefficient, temperature, and the unfrozen head
  head_up = surf_global_auxvar_up%head(1) + zc_up
  head_dn = surf_global_auxvar_dn%head(1) + zc_dn
  if (head_up > head_dn) then
    mannings_half          = mannings_up
    temp_half              = surf_global_auxvar_up%temp + 273.15d0 ! [K]
    unfrozen_fraction_half = surf_auxvar_up%unfrozen_fraction
    hw_half                = surf_global_auxvar_up%head(1)
  else
    mannings_half          = mannings_dn
    temp_half              = surf_global_auxvar_dn%temp + 273.15d0 ! [K]
    unfrozen_fraction_half = surf_auxvar_dn%unfrozen_fraction
    hw_half                = surf_global_auxvar_dn%head(1)
  endif

  ! We clip to avoid problems later evaluating at negative water height
  hw_half     = max(hw_half,MIN_SURFACE_WATER_HEIGHT)
  if (Equal(hw_half,MIN_SURFACE_WATER_HEIGHT)) then
    temp_half = 0.d0
    hw_half   = 0.d0
  endif

  ! Frozen water doesn't contribute to the velocity
  hw_liq_half = unfrozen_fraction_half*hw_half

  ! Compute Manning's velocity
  dhead = head_up - head_dn
  vel   = sign(hw_liq_half**(2.d0/3.d0)/mannings_half*abs(dhead/dist)**0.5d0,dhead) ! [m/s]

  ! KLUDGE: To address high velocity oscillations of the surface water
  ! height, reduce this value to keep dt from shrinking too much. Add
  ! to options if we decide to keep it.
  vel = sign(min(option%max_manning_velocity,abs(vel)),vel)

  ! Load into residual
  Res(TH_PRESSURE_DOF) = vel*hw_liq_half*length ! [m^3/s]
  
  ! Temperature equation
  ! RTM: k_therm is the weighted average of the liquid and ice thermal 
  ! conductivities.  For the density and specific heat capacity in the 
  ! advection term, we want these for liquid water ONLY, as the ice portion 
  ! is immobile and thus should not make up part of the advection term. We 
  ! also multiply the ponded water depth (hw_half) by the unfrozen fraction 
  ! in the advection term but NOT the conduction term.
  ! We do the same in SurfaceTHBCFlux().

  ! Average density
  ! Here we only consider the LIQUID fraction.
  den_aveg = (surf_global_auxvar_up%den_kg(1) + &
              surf_global_auxvar_dn%den_kg(1))/2.d0
  ! Temperature difference
  if (surf_global_auxvar_up%is_dry .or. surf_global_auxvar_dn%is_dry) then
    dtemp = 0.d0
  else
    dtemp = surf_global_auxvar_up%temp - surf_global_auxvar_dn%temp
  endif

  ! We are not being careful with dry/wet conditions, so if the
  ! temperature change is greater than 100 [C] we will assuming that
  ! it was a wet/dry interface change that was missed.
  if (abs(dtemp) > 100.d0) then
    den_aveg = 0.d0
    dtemp    = 0.d0
  endif

  ! Note, Cw and k_therm are same for up and downwind
  Cw = surf_auxvar_up%Cw
  k_therm = surf_auxvar_up%k_therm
  
  ! Unfrozen fraction multiplies hw_half in advection term, but does NOT affect the 
  ! conduction therm.  
  ! RTM: Brookfield et al. 2009 also has dispersion term, which we are not using.
  Res(TH_TEMPERATURE_DOF) = (den_aveg*vel*temp_half*Cw*hw_liq_half + &
                             k_therm*dtemp/dist*hw_half)*length

  if (abs(vel)>eps) then
    ! 1) Restriction due to flow equation
    dt     = dist/abs(vel)/3.d0
    dt_max = min(dt_max, dt)
  endif

  if (abs(dtemp) > 1.0d-12) then
    ! 2) Restriction due to energy equation
    dt_max = min(dt_max,(dist**2.d0)*Cw*den_aveg/(2.d0*k_therm))
  endif

end subroutine SurfaceTHFlux

! ************************************************************************** !

subroutine SurfaceTHBCFlux(ibndtype, &
                           auxvars, &
                           surf_auxvar_up, &
                           surf_global_auxvar_up, &
                           surf_auxvar_dn, &
                           surf_global_auxvar_dn, &
                           slope, &
                           mannings, &
                           dist, &
                           length, &
                           option, &
                           vel, &
                           dt_max, &
                           Res)
  ! 
  ! This routine computes flux for boundary cells.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/07/13
  !

  use Option_module
  use PFLOTRAN_Constants_module, only : MIN_SURFACE_WATER_HEIGHT

  implicit none

  type(option_type) :: option
  type(Surface_TH_auxvar_type) :: surf_auxvar_up
  type(surface_global_auxvar_type) :: surf_global_auxvar_up
  type(Surface_TH_auxvar_type) :: surf_auxvar_dn
  type(surface_global_auxvar_type) :: surf_global_auxvar_dn
  PetscReal :: auxvars(:) ! from aux_real_var array
  PetscReal :: slope
  PetscReal :: mannings
  PetscReal :: length
  PetscReal :: flux
  PetscInt :: ibndtype(:)
  PetscReal :: vel
  PetscReal :: dt_max
  PetscReal :: Res(1:option%nflowdof)
  PetscReal :: dist

  PetscInt :: pressure_bc_type
  PetscReal :: head,dhead
  PetscReal :: head_liq
  PetscReal :: den
  PetscReal :: temp_half
  PetscReal :: Cw
  PetscReal :: dtemp
  PetscReal :: hw_half
  PetscReal :: k_therm
  PetscReal :: dt

  flux = 0.d0
  vel = 0.d0
  hw_half = 0.d0
  dtemp = 0.d0
  Cw = 0.d0
  dt_max = PETSC_MAX_REAL

  ! Flow  
  pressure_bc_type = ibndtype(TH_PRESSURE_DOF)
  head = surf_global_auxvar_dn%head(1)
  k_therm = surf_auxvar_dn%k_therm
  
  select case(pressure_bc_type)
    case (ZERO_GRADIENT_BC)
      if (slope<0.d0) then
        vel =  0.d0
        head_liq = 0.d0
      else
        head_liq = surf_auxvar_dn%unfrozen_fraction * head
        vel = -sqrt(dabs(slope))/mannings*(head_liq**(2.d0/3.d0))
        hw_half = head
      endif
      den = surf_global_auxvar_dn%den_kg(1)
      Cw = surf_auxvar_dn%Cw
    case (NEUMANN_BC)
      vel = auxvars(TH_PRESSURE_DOF)
      den = (surf_global_auxvar_up%den_kg(1) + &
             surf_global_auxvar_dn%den_kg(1))/2.d0
    case (SPILLOVER_BC)
      ! if liquid water height is above a user-defined value, then outflow can occur
      head_liq =  surf_auxvar_dn%unfrozen_fraction * head
      dhead    =  max(head_liq-auxvars(1),0.0d0)
      vel      = -dhead**(2.d0/3.d0)/mannings*abs(dhead/dist)**0.5d0
      hw_half  =  head
      Cw       =  surf_auxvar_dn%Cw 
      den      =  surf_global_auxvar_dn%den_kg(1)
    case default
      option%io_buffer = 'Unknown pressure_bc_type for surface flow '
      call printErrMsg(option)
  end select

  if (vel>0.d0) then
    temp_half = surf_global_auxvar_up%temp + 273.15d0
  else
    temp_half = surf_global_auxvar_dn%temp + 273.15d0
  endif

  if (pressure_bc_type /= ZERO_GRADIENT_BC) then
    select case (ibndtype(TH_TEMPERATURE_DOF))
      case (DIRICHLET_BC)
        dtemp = surf_global_auxvar_up%temp - surf_global_auxvar_dn%temp
      case default
        option%io_buffer = 'Unknown temperature_bc_type for surface flow '
        call printErrMsg(option)
    end select
  endif

  flux = head_liq*vel
  Res(TH_PRESSURE_DOF)    = flux*length
  Res(TH_TEMPERATURE_DOF) = den*temp_half*Cw*vel*head_liq*length + &
                            k_therm*dtemp/dist*hw_half*length

  ! Timestep restriction due to mass equation
  if (abs(vel)>eps) then
    dt     = dist/abs(vel)/3.d0
    dt_max = min(dt_max, dt)
  endif
  ! Timestep restriction due to energy equation
  if (head_liq > MIN_SURFACE_WATER_HEIGHT) then
    dt_max = min(dt_max,(dist**2.d0)*Cw*den/(2.d0*k_therm))
  endif

end subroutine SurfaceTHBCFlux

! ************************************************************************** !

subroutine SurfaceTHUpdateAuxVars(surf_realization)
  ! 
  ! This routine updates auxiliary variables
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/07/13
  ! 

  use Realization_Surface_class
  use Patch_module
  use Option_module
  use Surface_Field_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Surface_Material_module
  use PFLOTRAN_Constants_module, only : MIN_SURFACE_WATER_HEIGHT

  implicit none

  class(realization_surface_type) :: surf_realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(surface_field_type), pointer :: surf_field
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  type(Surface_TH_auxvar_type), pointer :: surf_th_auxvars(:)
  type(Surface_TH_auxvar_type), pointer :: surf_th_auxvars_bc(:)
  type(Surface_TH_auxvar_type), pointer :: surf_th_auxvars_ss(:)
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars(:)
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars_bc(:)
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars_ss(:)

  PetscInt :: ghosted_id, local_id, istart, iend, sum_connection, idof, iconn
  PetscInt :: iphasebc, iphase
  PetscReal, pointer :: xx_loc_p(:), icap_loc_p(:), iphase_loc_p(:)
  PetscReal, pointer :: perm_xx_loc_p(:), porosity_loc_p(:)
  PetscReal :: xxbc(surf_realization%option%nflowdof)
  PetscReal :: xxss(surf_realization%option%nflowdof)
  PetscReal :: tsrc1
  PetscErrorCode :: ierr
  PetscReal :: den,head

  option => surf_realization%option
  patch => surf_realization%patch
  grid => patch%grid
  surf_field => surf_realization%surf_field

  surf_th_auxvars => patch%surf_aux%SurfaceTH%auxvars
  surf_th_auxvars_bc => patch%surf_aux%SurfaceTH%auxvars_bc
  surf_th_auxvars_ss => patch%surf_aux%SurfaceTH%auxvars_ss
  surf_global_auxvars => patch%surf_aux%SurfaceGlobal%auxvars
  surf_global_auxvars_bc => patch%surf_aux%SurfaceGlobal%auxvars_bc
  surf_global_auxvars_ss => patch%surf_aux%SurfaceGlobal%auxvars_ss
  
  call VecGetArrayF90(surf_field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)

  ! Internal aux vars
  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells

    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    iend = ghosted_id*option%nflowdof
    istart = iend-option%nflowdof+1

    call SurfaceTHAuxVarCompute(xx_loc_p(istart:iend), &
                                surf_th_auxvars(ghosted_id), &
                                surf_global_auxvars(ghosted_id), &
                                option)
    ! [rho*h*T*Cwi]
    if (xx_loc_p(istart) >= MIN_SURFACE_WATER_HEIGHT) then
      xx_loc_p(istart+1) = surf_global_auxvars(ghosted_id)%den_kg(1)* &
                           xx_loc_p(istart)* &
                           (surf_global_auxvars(ghosted_id)%temp + 273.15d0)* &
                           surf_th_auxvars(ghosted_id)%Cwi
    else
      xx_loc_p(istart+1) = 0.d0
    endif
  enddo
   
  ! Boundary aux vars
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif

      do idof=1,option%nflowdof
        select case(boundary_condition%flow_condition%itype(idof))
          case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC,NEUMANN_BC, &
               HET_DIRICHLET_BC)
            xxbc(idof) = boundary_condition%flow_aux_real_var(idof,iconn)
          case(ZERO_GRADIENT_BC)
            xxbc(idof) = xx_loc_p((ghosted_id-1)*option%nflowdof+idof)
        end select
      enddo

      surf_global_auxvars_bc(sum_connection)%temp = xxbc(2)
      call SurfaceTHAuxVarCompute(xxbc, &
                                  surf_th_auxvars_bc(sum_connection), &
                                  surf_global_auxvars_bc(sum_connection), &
                                  option)

    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! Source/Sink aux vars
  ! source/sinks
  source_sink => patch%source_sink_list%first
  sum_connection = 0
  do
    if (.not.associated(source_sink)) exit
    cur_connection_set => source_sink%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle

      iend = ghosted_id*option%nflowdof
      istart = iend-option%nflowdof+1

      if (associated(source_sink%flow_condition%temperature)) then
        if (source_sink%flow_condition%temperature%itype /= &
            HET_DIRICHLET_BC) then
          tsrc1 = source_sink%flow_condition%temperature%dataset%rarray(1)
        else
          tsrc1 = source_sink%flow_aux_real_var(TWO_INTEGER,iconn)
        endif
      else
        tsrc1 = xx_loc_p((ghosted_id-1)*option%nflowdof+1)
        tsrc1 = surf_global_auxvars(ghosted_id)%temp
      endif

      xxss = xx_loc_p(istart:iend)
      head    = xxss(1)
      xxss(1) = 1.d0 ! set arbitrary amount of surface water so auxvar will evaluate
      xxss(2) = tsrc1

      surf_global_auxvars_ss(sum_connection)%temp = tsrc1
      call SurfaceTHAuxVarCompute(xxss, &
                                  surf_th_auxvars_ss(sum_connection), &
                                  surf_global_auxvars_ss(sum_connection), &
                                  option)
      surf_global_auxvars_ss(sum_connection)%head = head ! set head back just in case

    enddo
    source_sink => source_sink%next
  enddo

  patch%surf_aux%SurfaceTH%auxvars_up_to_date = PETSC_TRUE

  call VecRestoreArrayF90(surf_field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)

end subroutine SurfaceTHUpdateAuxVars

! ************************************************************************** !

subroutine EnergyToTemperatureBisection(T,TL,TR,h,energy,Cwi,Pr,option)
  ! 
  ! Solves the following nonlinear equation using the bisection method
  !
  ! R(T) = rho(T) Cwi hw T - energy = 0
  !
  ! Author: Nathan Collier, ORNL
  ! Date: 11/2014
  ! 
  use EOS_Water_module
  use Option_module

  implicit none

  PetscReal :: T,TL,TR,h,energy,Cwi,Pr
  type(option_type), pointer :: option

  PetscReal :: Tp,rho,rho_t,f,fR,fL,rtol
  PetscInt :: iter,niter
  PetscBool :: found
  PetscErrorCode :: ierr

  call EOSWaterdensity(TR,Pr,rho,rho_T,ierr)
  fR = rho*Cwi*h*(TR+273.15d0) - energy
  call EOSWaterdensity(TL,Pr,rho,rho_T,ierr)
  fL = rho*Cwi*h*(TL+273.15d0) - energy

  if (fL*fR > 0.d0) then
     print *,"[TL,TR] = ",TL,TR
     print *,"[fL,fR] = ",fL,fR
     write(option%io_buffer,'("surface_th.F90: EnergyToTemperatureBisection --> root is not bracketed")')
     call printErrMsg(option)
  endif

  T = 0.5d0*(TL+TR)
  call EOSWaterdensity(T,Pr,rho,rho_T,ierr)
  f = rho*Cwi*h*(T+273.15d0) - energy

  found = PETSC_FALSE
  niter = 200
  rtol  = 1.d-6
  do iter = 1,niter
     Tp = T
     if (fL*f < 0.d0) then
        TR = T
     else 
        TL = T
     endif

     T = 0.5d0*(TL+TR)

     call EOSWaterdensity(T,Pr,rho,rho_T,ierr)
     f = rho*Cwi*h*(T+273.15d0) - energy

     if (abs((T-Tp)/(T+273.15d0)) < rtol) then
        found = PETSC_TRUE
        exit
     endif
  enddo

  if (found .eqv. PETSC_FALSE) then
     print *,"[TL,T,TR] = ",TL,T,TR
     write(option%io_buffer,'("surface_th.F90: EnergyToTemperatureBisection --> root not found!")')
     call printErrMsg(option)
  endif

end subroutine EnergyToTemperatureBisection

! ************************************************************************** !

subroutine SurfaceTHUpdateTemperature(surf_realization)
  ! 
  ! This routine updates the temperature after TSSolve.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/25/13
  ! 

  use Realization_Surface_class
  use Patch_module
  use Option_module
  use Surface_Field_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Surface_Material_module
  use EOS_Water_module
  use PFLOTRAN_Constants_module, only : DUMMY_VALUE,MIN_SURFACE_WATER_HEIGHT

  implicit none

  class(realization_surface_type) :: surf_realization
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(surface_field_type), pointer :: surf_field
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  type(Surface_TH_auxvar_type), pointer :: surf_auxvars(:)
  type(Surface_TH_auxvar_type), pointer :: surf_auxvars_bc(:)
  type(Surface_TH_auxvar_type), pointer :: surf_auxvars_ss(:)
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars(:)
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars_bc(:)
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars_ss(:)

  PetscInt :: ghosted_id, local_id, istart, iend, sum_connection, idof, iconn
  PetscInt :: iphasebc, iphase
  PetscReal, pointer :: xx_loc_p(:), icap_loc_p(:), iphase_loc_p(:)
  PetscReal, pointer :: perm_xx_loc_p(:), porosity_loc_p(:)
  PetscReal :: xxbc(surf_realization%option%nflowdof)
  PetscReal :: xxss(surf_realization%option%nflowdof)
  PetscReal :: temp,TL,TR
  PetscReal :: den
  PetscReal :: dum1
  PetscErrorCode :: ierr

  option => surf_realization%option
  patch => surf_realization%patch
  grid => patch%grid
  surf_field => surf_realization%surf_field

  surf_global_auxvars => patch%surf_aux%SurfaceGlobal%auxvars
  surf_global_auxvars_bc => patch%surf_aux%SurfaceGlobal%auxvars_bc
  surf_global_auxvars_ss => patch%surf_aux%SurfaceGlobal%auxvars_ss
  surf_auxvars => patch%surf_aux%SurfaceTH%auxvars
  surf_auxvars_bc => patch%surf_aux%SurfaceTH%auxvars_bc

  !
  ! The unknown for the energy balance in the surface domain is
  ! energy. Thus we need to compute a temperature, which results in
  ! finding the root of the following nonlinear equation,
  !
  ! Residual(T) = rho(T) Cwi hw T - energy = 0
  !

  call VecGetArrayF90(surf_field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)

  do ghosted_id = 1,grid%ngmax
    istart = (ghosted_id-1)*option%nflowdof+1 ! surface water height dof
    iend   = istart+1                       ! surface energy dof
    if (xx_loc_p(istart) < MIN_SURFACE_WATER_HEIGHT) then
      ! If the cell is dry then we set temperature to a dummy value
      ! and then zero out the water height and energy.
      surf_global_auxvars(ghosted_id)%is_dry = PETSC_TRUE
      temp = DUMMY_VALUE
      xx_loc_p(istart) = 0.d0 ! no water 
      xx_loc_p(iend)   = 0.d0 ! no energy
    else
      TL = -100.d0
      TR =  100.d0
      call EnergyToTemperatureBisection(temp,TL,TR, &
                                        xx_loc_p(istart), &
                                        xx_loc_p(iend), &
                                        surf_auxvars(ghosted_id)%Cwi, &
                                        option%reference_pressure,option)
    endif
    surf_global_auxvars(ghosted_id)%temp = temp
    call EOSWaterdensity(temp,option%reference_pressure,den,dum1,ierr)
    surf_global_auxvars(ghosted_id)%den_kg(1) = den
  enddo

  call VecRestoreArrayF90(surf_field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)

end subroutine SurfaceTHUpdateTemperature

! ************************************************************************** !

subroutine SurfaceTHUpdateSurfState(surf_realization)
  ! 
  ! This routine updates the states for surface-model at the end of
  ! subsurface-model timestep.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/25/13
  ! 

#include "petsc/finclude/petscmat.h"
  use petscmat
  use Connection_module
  use Coupler_module
  use Discretization_module
  use DM_Kludge_module
  use Grid_module
  use Option_module
  use Patch_module
  use Realization_Subsurface_class
  use Realization_Base_class
  use String_module
  use Surface_Field_module
  use Realization_Surface_class
  use EOS_Water_module

  implicit none

  class(realization_surface_type) :: surf_realization

  type(coupler_list_type), pointer :: coupler_list
  type(coupler_type), pointer :: coupler
  type(connection_set_type), pointer :: cur_connection_set
  type(dm_ptr_type), pointer :: dm_ptr
  type(grid_type),pointer :: grid,surf_grid
  type(option_type), pointer :: option
  type(patch_type),pointer :: patch,surf_patch
  type(surface_field_type),pointer :: surf_field
  type(Surface_TH_auxvar_type), pointer :: surf_auxvars(:)

  PetscInt :: count
  PetscInt :: ghosted_id
  PetscInt :: local_id
  PetscInt :: ibeg
  PetscInt :: iend
  PetscInt :: iconn
  PetscInt :: sum_connection

  PetscReal :: den
  PetscReal :: dum1
  PetscReal, pointer :: avg_vdarcy_p(:)   ! avg darcy velocity [m/s]
  PetscReal, pointer :: xx_p(:)           ! head [m]
  PetscReal, pointer :: surfpress_p(:)
  PetscReal, pointer :: surftemp_p(:)
  PetscReal :: Cwi
  PetscReal :: temp_K
  PetscErrorCode :: ierr

  PetscBool :: coupler_found = PETSC_FALSE

  patch      => surf_realization%patch
  option     => surf_realization%option
  surf_field => surf_realization%surf_field
  surf_grid  => surf_realization%discretization%grid
  surf_auxvars => patch%surf_aux%SurfaceTH%auxvars

  call VecGetArrayF90(surf_field%flow_xx, xx_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(surf_field%press_subsurf, surfpress_p,  &
                      ierr);CHKERRQ(ierr)
  call VecGetArrayF90(surf_field%temp_subsurf, surftemp_p, ierr);CHKERRQ(ierr)

  count = 0
  do ghosted_id = 1,surf_grid%ngmax

    local_id = surf_grid%nG2L(ghosted_id)
    if (local_id <= 0) cycle

    iend = ghosted_id*option%nflowdof
    ibeg = iend - 1

    ! Compute density
    count = count + 1
    call EOSWaterdensity(surftemp_p(count),option%reference_pressure,den,dum1,ierr)
    xx_p(ibeg) = (surfpress_p(count)-option%reference_pressure)/ &
                        (abs(option%gravity(3)))/den
    if (surfpress_p(count)-option%reference_pressure < 1.0d-8) then
      xx_p(ibeg) = 0.d0
      xx_p(iend) = 0.d0
    else
      Cwi = surf_auxvars(ghosted_id)%Cwi
      temp_K = surftemp_p(count) + 273.15d0
      xx_p(iend) = den*Cwi*temp_K*xx_p(ibeg)
    endif

  enddo
  call VecRestoreArrayF90(surf_field%flow_xx, xx_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(surf_field%press_subsurf, surfpress_p,  &
                          ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(surf_field%temp_subsurf, surftemp_p,  &
                          ierr);CHKERRQ(ierr)

  call DiscretizationGlobalToLocal(surf_realization%discretization, &
                                   surf_field%flow_xx, &
                                   surf_field%flow_xx_loc, &
                                   NFLOWDOF)
  call SurfaceTHUpdateAuxVars(surf_realization)

end subroutine SurfaceTHUpdateSurfState

! ************************************************************************** !

subroutine AtmEnergyToTemperatureBisection(T,TL,TR,shift,RHS,Pr,option)
  ! 
  ! Solves the following nonlinear equation using the bisection method
  !
  ! R(T) = (rho(T)+shift)*T - RHS = 0
  !
  ! Author: Nathan Collier, ORNL
  ! Date: 11/2014
  ! 
  use EOS_Water_module
  use Option_module

  implicit none

  PetscReal :: T,TL,TR,shift,RHS,Pr
  type(option_type), pointer :: option

  PetscReal :: Tp,rho,rho_t,f,fR,fL,rtol
  PetscInt :: iter,niter
  PetscBool :: found
  PetscErrorCode :: ierr

  call EOSWaterdensity(TR,Pr,rho,rho_T,ierr)
  fR = (rho+shift)*(TR+273.15d0) - RHS
  call EOSWaterdensity(TL,Pr,rho,rho_T,ierr)
  fL = (rho+shift)*(TL+273.15d0) - RHS

  if (fL*fR > 0.d0) then
     print *,"[TL,TR] = ",TL,TR
     print *,"[fL,fR] = ",fL,fR
     write(option%io_buffer,'("surface_th.F90: AtmEnergyToTemperatureBisection --> root is not bracketed")')
     call printErrMsg(option)
  endif

  T = 0.5d0*(TL+TR)
  call EOSWaterdensity(T,Pr,rho,rho_T,ierr)
  f = (rho+shift)*(T+273.15d0) - RHS

  found = PETSC_FALSE
  niter = 200
  rtol  = 1.d-6
  do iter = 1,niter
     Tp = T
     if (fL*f < 0.d0) then
        TR = T
     else 
        TL = T
     endif

     T = 0.5d0*(TL+TR)

     call EOSWaterdensity(T,Pr,rho,rho_T,ierr)
     f = (rho+shift)*(T+273.15d0) - RHS

     if (abs((T-Tp)/(T+273.15d0)) < rtol) then
        found = PETSC_TRUE
        exit
     endif
  enddo

  if (found .eqv. PETSC_FALSE) then
     print *,"[TL,T,TR] = ",TL,T,TR
     write(option%io_buffer,'("surface_th.F90: AtmEnergyToTemperatureBisection --> root not found!")')
     call printErrMsg(option)
  endif

end subroutine AtmEnergyToTemperatureBisection

! ************************************************************************** !

subroutine SurfaceTHImplicitAtmForcing(surf_realization)
  !
  ! Updates the temperature of surface-water implicitly due to conduction.
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/24/2014
  !

  use Realization_Surface_class
  use Patch_module
  use Option_module
  use Surface_Field_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Surface_Material_module
  use EOS_Water_module
  use String_module
  use PFLOTRAN_Constants_module, only : MIN_SURFACE_WATER_HEIGHT
  implicit none

  class(realization_surface_type) :: surf_realization
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(surface_field_type), pointer :: surf_field
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  type(Surface_TH_auxvar_type), pointer :: surf_auxvars(:)
  type(Surface_TH_auxvar_type), pointer :: surf_auxvars_bc(:)
  type(Surface_TH_auxvar_type), pointer :: surf_auxvars_ss(:)
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars(:)
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars_bc(:)
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars_ss(:)

  PetscInt :: ghosted_id, local_id, istart, iend, sum_connection, idof, iconn
  PetscInt :: iphasebc, iphase
  PetscReal, pointer :: xx_loc_p(:), xx_p(:)
  PetscReal, pointer :: perm_xx_loc_p(:), porosity_loc_p(:)
  PetscReal :: xxbc(surf_realization%option%nflowdof)
  PetscReal :: xxss(surf_realization%option%nflowdof)
  PetscReal :: temp,ptemp,rtol
  PetscInt :: iter
  PetscInt :: niter
  PetscReal :: den
  PetscReal :: dum1
  PetscReal :: den_iter
  PetscReal :: den_old
  PetscReal :: k_therm
  PetscReal :: Cw
  PetscReal :: temp_old
  PetscReal :: head
  PetscReal :: beta,RHS,TL,TR
  PetscBool :: found
  PetscErrorCode :: ierr

  option => surf_realization%option
  patch => surf_realization%patch
  grid => patch%grid
  surf_field => surf_realization%surf_field

  surf_global_auxvars => patch%surf_aux%SurfaceGlobal%auxvars
  surf_global_auxvars_ss => patch%surf_aux%SurfaceGlobal%auxvars_ss
  surf_auxvars => patch%surf_aux%SurfaceTH%auxvars
  surf_auxvars_bc => patch%surf_aux%SurfaceTH%auxvars_bc

  ! niter = max(m)
  niter = 20
  rtol  = 1.d-12
  call VecGetArrayF90(surf_field%flow_xx,xx_p,ierr);CHKERRQ(ierr)

  ! Update source/sink aux vars
  source_sink => patch%source_sink_list%first
  sum_connection = 0
  do
    if (.not.associated(source_sink)) exit

    cur_connection_set => source_sink%connection_set

    if (StringCompare(source_sink%name,'atm_energy_ss')) then

      if (source_sink%flow_condition%itype(TH_TEMPERATURE_DOF) == &
          HET_DIRICHLET_BC) then

        do iconn = 1, cur_connection_set%num_connections

          sum_connection = sum_connection + 1

          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)

          head     = surf_global_auxvars(ghosted_id)%head(1)
          temp_old = surf_global_auxvars(ghosted_id)%temp
          k_therm  = surf_auxvars(ghosted_id)%k_therm
          Cw       = surf_auxvars(ghosted_id)%Cw

          if (head > MIN_SURFACE_WATER_HEIGHT) then

            call EOSWaterdensity(temp_old,option%reference_pressure,den_old,dum1,ierr)
            call EOSWaterdensity(temp_old,option%reference_pressure,den_iter,dum1,ierr)

            TL    = -100.d0
            TR    =  100.d0
            beta  = (2.d0*k_therm*option%surf_flow_dt)/(Cw*head**2.d0)
            RHS   =  den_old*(temp_old+273.15d0)+beta*(surf_global_auxvars_ss(sum_connection)%temp+273.15d0)
            call AtmEnergyToTemperatureBisection(temp,TL,TR,beta,RHS,option%reference_pressure,option)

            call EOSWaterdensity(temp,option%reference_pressure,den_iter,dum1,ierr)
            surf_global_auxvars(ghosted_id)%temp = temp

            iend = local_id*option%nflowdof
            istart = iend - option%nflowdof + 1
            xx_p(iend) = den_iter*Cw*(temp + 273.15d0)*xx_p(istart)
          endif

        enddo

      else
        sum_connection = sum_connection + cur_connection_set%num_connections
      endif

    else
      sum_connection = sum_connection + cur_connection_set%num_connections
    endif

    source_sink => source_sink%next
  enddo

  call VecRestoreArrayF90(surf_field%flow_xx,xx_p,ierr);CHKERRQ(ierr)

end subroutine SurfaceTHImplicitAtmForcing

! ************************************************************************** !

subroutine SurfaceTHUpdateSolution(surf_realization)
  ! 
  ! This routine updates solution after a successful time step
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/07/13
  ! 

  use Realization_Surface_class
  use Surface_Field_module

  implicit none

  class(realization_surface_type) :: surf_realization

  type(surface_field_type),pointer :: surf_field
  PetscErrorCode :: ierr

  surf_field => surf_realization%surf_field
  call VecCopy(surf_field%flow_xx,surf_field%flow_yy,ierr);CHKERRQ(ierr)

end subroutine SurfaceTHUpdateSolution


! ************************************************************************** !

subroutine SurfaceTHDestroy(surf_realization)
  ! 
  ! Deallocates variables associated with Richard
  ! 
  ! Author: ???
  ! Date: 02/14/08
  ! 

  use Realization_Surface_class

  implicit none
  
  class(realization_surface_type) :: surf_realization
  
  ! aux vars should be destroyed when surf_realization is destroyed.
  
end subroutine SurfaceTHDestroy

end module Surface_TH_module
