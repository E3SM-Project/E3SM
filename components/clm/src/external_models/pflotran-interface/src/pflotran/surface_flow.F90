module Surface_Flow_module

#include "petsc/finclude/petscmat.h"
  use petscmat
  use Global_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none
  
  private

! Cutoff parameters
  PetscReal, parameter :: eps       = 1.D-12
  PetscReal, parameter :: perturbation_tolerance = 1.d-6

  public SurfaceFlowSetup, &
         SurfaceFlowDiffusion, &
         SurfaceFlowUpdateSolution, &
         SurfaceFlowRHSFunction, &
         SurfaceFlowComputeMaxDt, &
         SurfaceFlowGetTecplotHeader, &
         SurfaceFlowUpdateAuxVars, &
         SurfaceFlowUpdateSurfState

contains

! ************************************************************************** !

subroutine SurfaceFlowSetup(surf_realization)
  ! 
  ! This routine sets up surfaceflow type
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 05/21/12
  ! 

  use Realization_Surface_class
  use Output_Aux_module
  
  class(realization_surface_type) :: surf_realization

  type(output_variable_list_type), pointer :: list

  list => surf_realization%output_option%output_snap_variable_list
  call SurfaceFlowSetPlotVariables(list)
  list => surf_realization%output_option%output_obs_variable_list
  call SurfaceFlowSetPlotVariables(list)
  
end subroutine SurfaceFlowSetup

! ************************************************************************** !

subroutine SurfaceFlowSetPlotVariables(list)
  ! 
  ! This routine adds default variables to be printed to list
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 10/30/12
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

  name = 'Material ID'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_DISCRETE,units, &
                               MATERIAL_ID)
  
end subroutine SurfaceFlowSetPlotVariables

! ************************************************************************** !

subroutine SurfaceFlowKinematic(hw_up, &
                                mannings_up, &
                                hw_dn, &
                                mannings_dn, &
                                slope, &
                                length, &
                                option, &
                                vel, &
                                Res)
  ! 
  ! This routine computes the internal flux term for the residual under
  ! kinematic-wave assumption.
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 05/21/12
  ! 

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  PetscReal :: hw_up, hw_dn
  PetscReal :: slope
  PetscReal :: mannings_up, mannings_dn
  PetscReal :: length
  PetscReal :: vel                      ! units: m/s
  PetscReal :: Res(1:option%nflowdof)   ! units: m^3/s
  
  PetscReal :: flux       ! units: m^2/s
  
  ! initialize
  flux = 0.d0
  vel  = 0.d0
  
  if (slope<0.d0) then
    vel =  sqrt(dabs(slope))/mannings_up*((hw_up)**(2.d0/3.d0))
    flux=  hw_up*vel
  else
    vel = -sqrt(dabs(slope))/mannings_dn*((hw_dn)**(2.d0/3.d0))
    flux=  hw_dn*vel
  endif

  Res(1) = flux*length

end subroutine SurfaceFlowKinematic

! ************************************************************************** !

subroutine SurfaceFlowDiffusion(hw_up, &
                                zc_up, &
                                mannings_up, &
                                hw_dn, &
                                zc_dn, &
                                mannings_dn, &
                                dist, &
                                length, &
                                option, &
                                vel, &
                                Res)
  ! 
  ! This routine computes the internal flux term for the residual under
  ! diffusion-wave assumption.
  ! 
  ! Author: Gautam Bisht, LBL
  ! Date: 08/03/12
  ! 

  use Option_module

  implicit none

  type(option_type) :: option
  PetscReal :: hw_up, hw_dn
  PetscReal :: zc_up, zc_dn
  PetscReal :: head_up, head_dn
  PetscReal :: mannings_up, mannings_dn
  PetscReal :: dist, length
  PetscReal :: vel                      ! units: m/s
  PetscReal :: Res(1:option%nflowdof)   ! units: m^3/s

  PetscReal :: flux       ! units: m^2/s
  PetscReal :: Cd
  PetscReal :: hw_half
  PetscReal :: mannings_half
  PetscReal :: dhead

  ! initialize
  flux = 0.d0
  Cd = 1.0d0

  head_up = hw_up + zc_up
  head_dn = hw_dn + zc_dn

  if (head_up>head_dn) then
    mannings_half = mannings_up
    if (hw_up>0.d0) then
      hw_half = hw_up
    else
      hw_half = 0.d0
    endif
  else
    mannings_half = mannings_dn
    if (hw_dn>0.d0) then
      hw_half = hw_dn
    else
      hw_half = 0.d0
    endif
  endif
  
  dhead=head_up-head_dn
  if (abs(dhead)<eps) then
    dhead=0.d0
    vel = 0.d0
  else
    vel = (hw_half**(2.d0/3.d0))/mannings_half* &
          dhead/(abs(dhead)**(1.d0/2.d0))* &
          1.d0/(dist**0.5d0)
  endif

  flux = hw_half*vel
  Res(1) = flux*length

end subroutine SurfaceFlowDiffusion

! ************************************************************************** !

subroutine SurfaceFlowUpdateSolution(surf_realization)
  ! 
  ! This routine updates data in module after a successful time step
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 05/22/12
  ! 

  use Realization_Surface_class
  use Surface_Field_module

  implicit none

  class(realization_surface_type) :: surf_realization

  type(surface_field_type),pointer :: surf_field
  PetscErrorCode :: ierr

  surf_field => surf_realization%surf_field
  call VecCopy(surf_field%flow_xx,surf_field%flow_yy,ierr);CHKERRQ(ierr)

end subroutine SurfaceFlowUpdateSolution

! ************************************************************************** !

subroutine SurfaceFlowRHSFunction(ts,t,xx,ff,surf_realization,ierr)
  ! 
  ! This routine provides the function evaluation for PETSc TSSolve()
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/07/13
  ! 
#include <petsc/finclude/petscts.h>
  use petscts
  use Realization_Surface_class
  use Surface_Field_module
  use Patch_module
  use Discretization_module
  use Option_module
  use Logging_module
  use Connection_module
  use Grid_module
  use Coupler_module
  use Surface_Field_module
  use Debug_module
  use Surface_Global_Aux_module

  implicit none
  
  TS :: ts
  PetscReal :: t
  Vec :: xx
  Vec :: ff
  class(realization_surface_type) :: surf_realization
  PetscErrorCode :: ierr

  PetscViewer :: viewer

  type(discretization_type), pointer :: discretization
  type(surface_field_type), pointer :: surf_field
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars(:)
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars_bc(:)

  PetscInt :: local_id_up, local_id_dn, local_id
  PetscInt :: ghosted_id_up, ghosted_id_dn, ghosted_id
  PetscInt :: iconn
  PetscInt :: sum_connection
  PetscReal :: dx, dy, dz
  PetscReal :: dist
  PetscReal :: vel
  PetscReal :: slope, slope_dn
  PetscReal :: hw_up, hw_dn ! water height [m]
  PetscReal :: Res(surf_realization%option%nflowdof), v_darcy
  PetscReal :: qsrc, qsrc_flow

  character(len=MAXSTRINGLENGTH) :: string,string2

  PetscReal, pointer :: ff_p(:), mannings_loc_p(:),area_p(:)
  PetscReal, pointer :: xc(:),yc(:),zc(:)

  patch => surf_realization%patch
  grid => patch%grid
  option => surf_realization%option
  surf_field => surf_realization%surf_field
  surf_global_auxvars => patch%surf_aux%SurfaceGlobal%auxvars
  surf_global_auxvars_bc => patch%surf_aux%SurfaceGlobal%auxvars_bc

  surf_field              => surf_realization%surf_field
  discretization          => surf_realization%discretization
  option                  => surf_realization%option
  patch                   => surf_realization%patch
  grid                    => patch%grid
  surf_global_auxvars    => patch%surf_aux%SurfaceGlobal%auxvars
  surf_global_auxvars_bc => patch%surf_aux%SurfaceGlobal%auxvars_bc
  
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

  call DiscretizationGlobalToLocal(discretization,xx,surf_field%flow_xx_loc,NFLOWDOF)
  ! Then, update the aux vars
  call SurfaceFlowUpdateAuxVars(surf_realization)

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
      
      if (surf_global_auxvars(ghosted_id_up)%head(1)<0.d0 .or. &
         surf_global_auxvars(ghosted_id_dn)%head(1)<0.d0) then
        write(*,*) 'In SurfaceFlowFlux: ', surf_global_auxvars(ghosted_id_up)%head(1), &
          surf_global_auxvars(ghosted_id_dn)%head(1),ghosted_id_up,ghosted_id_dn
          option%io_buffer='stopping: -ve head values '
          call printErrMsg(option)
      endif

      select case(option%surface_flow_formulation)
        case (KINEMATIC_WAVE)
          option%io_buffer='Explicit Surface flow not implemented for ' // &
            'Kinematic wave'
          call printErrMsg(option)
        case (DIFFUSION_WAVE)
          call SurfaceFlowFlux(surf_global_auxvars(ghosted_id_up), &
                               zc(ghosted_id_up), &
                               mannings_loc_p(ghosted_id_up), &
                               surf_global_auxvars(ghosted_id_dn), &
                               zc(ghosted_id_dn), &
                               mannings_loc_p(ghosted_id_dn), &
                               dist, cur_connection_set%area(iconn), &
                               option,vel,Res)
      end select

      patch%internal_velocities(1,sum_connection) = vel
      patch%internal_flow_fluxes(1,sum_connection) = Res(1)

      if (local_id_up>0) then
        ff_p(local_id_up) = ff_p(local_id_up) - Res(1)/area_p(local_id_up)
      endif
         
      if (local_id_dn>0) then
        ff_p(local_id_dn) = ff_p(local_id_dn) + Res(1)/area_p(local_id_dn)
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
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id_dn = grid%nL2G(local_id)
  
      dx = xc(ghosted_id_dn) - cur_connection_set%intercp(1,iconn)
      dy = yc(ghosted_id_dn) - cur_connection_set%intercp(2,iconn)
      dz = zc(ghosted_id_dn) - cur_connection_set%intercp(3,iconn)
      slope_dn = dz/sqrt(dx*dx + dy*dy + dz*dz)

      call SurfaceFlowBCFlux(boundary_condition%flow_condition%itype, &
                         surf_global_auxvars_bc(sum_connection), &
                         slope_dn, &
                         mannings_loc_p(ghosted_id_dn), &
                         cur_connection_set%area(iconn), &
                         option,vel,Res)

      patch%boundary_velocities(1,sum_connection) = vel
      patch%boundary_flow_fluxes(1,sum_connection) = Res(1)

      ff_p(local_id) = ff_p(local_id) + Res(1)/area_p(local_id)
    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sink_list%first
  sum_connection = 0
  do
    if (.not.associated(source_sink)) exit
    
    qsrc_flow = 0.d0
    if (source_sink%flow_condition%rate%itype/=HET_VOL_RATE_SS.and. &
       source_sink%flow_condition%rate%itype/=HET_MASS_RATE_SS) &
    qsrc_flow = source_sink%flow_condition%rate%dataset%rarray(1)
      
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
      
      ff_p(local_id) = ff_p(local_id) + qsrc/area_p(local_id)
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

end subroutine SurfaceFlowRHSFunction

! ************************************************************************** !

subroutine SurfaceFlowComputeMaxDt(surf_realization,max_allowable_dt)
  ! 
  ! This routine maximum allowable 'dt' for surface flow model.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/07/13
  ! 

  
  use Connection_module
  use Realization_Surface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module  
  use Surface_Field_module
  use Debug_module
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
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars(:)
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars_bc(:)

  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  PetscInt :: iconn
  PetscInt :: sum_connection

  PetscReal :: dx, dy, dz
  PetscReal :: dist
  PetscReal :: vel
  PetscReal :: slope, slope_dn
  PetscReal :: rho          ! density      [kg/m^3]
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
  surf_global_auxvars => patch%surf_aux%SurfaceGlobal%auxvars
  surf_global_auxvars_bc => patch%surf_aux%SurfaceGlobal%auxvars_bc

  call VecGetArrayF90(surf_field%mannings_loc,mannings_loc_p,  &
                      ierr);CHKERRQ(ierr)
  call VecGetArrayF90(surf_field%area,area_p,ierr);CHKERRQ(ierr)

  Res  = 0.d0
  max_allowable_dt = 1.d10

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
      
      select case(option%surface_flow_formulation)
        case (KINEMATIC_WAVE)
        case (DIFFUSION_WAVE)
          call SurfaceFlowFlux(surf_global_auxvars(ghosted_id_up), &
                               zc(ghosted_id_up), &
                               mannings_loc_p(ghosted_id_up), &
                               surf_global_auxvars(ghosted_id_dn), &
                               zc(ghosted_id_dn), &
                               mannings_loc_p(ghosted_id_dn), &
                               dist, cur_connection_set%area(iconn), &
                               option,vel,Res)
      end select

      patch%internal_velocities(1,sum_connection) = vel
      patch%internal_flow_fluxes(1,sum_connection) = Res(1)
      if (abs(vel)>eps) then
        dt = dist/abs(vel)/4.d0
        max_allowable_dt = min(max_allowable_dt,dt)
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
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id_dn = grid%nL2G(local_id)
  
      dx = xc(ghosted_id_dn) - cur_connection_set%intercp(1,iconn)
      dy = yc(ghosted_id_dn) - cur_connection_set%intercp(2,iconn)
      dz = zc(ghosted_id_dn) - cur_connection_set%intercp(3,iconn)
      dist = sqrt(dx*dx + dy*dy + dz*dz)
      slope_dn = dz/dist

      call SurfaceFlowBCFlux(boundary_condition%flow_condition%itype, &
                         surf_global_auxvars_bc(sum_connection), &
                         slope_dn, &
                         mannings_loc_p(ghosted_id_dn), &
                         cur_connection_set%area(iconn), &
                         option,vel,Res)
      patch%boundary_velocities(1,sum_connection) = vel
      patch%boundary_flow_fluxes(1,sum_connection) = Res(1)

      if (abs(vel)>eps) then
        dt = dist/abs(vel)/4.d0
        max_allowable_dt = min(max_allowable_dt,dt)
      endif
    enddo
    boundary_condition => boundary_condition%next
  enddo

  call VecRestoreArrayF90(surf_field%mannings_loc,mannings_loc_p, &
                          ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(surf_field%area,area_p,ierr);CHKERRQ(ierr)

end subroutine SurfaceFlowComputeMaxDt

! ************************************************************************** !

subroutine SurfaceFlowFlux(surf_global_auxvar_up, &
                         zc_up, &
                         mannings_up, &
                         surf_global_auxvar_dn, &
                         zc_dn, &
                         mannings_dn, &
                         dist, &
                         length, &
                         option, &
                         vel, &
                         Res)
  ! 
  ! This routine computes the internal flux term for under
  ! diffusion-wave assumption.
  ! 
  ! Author: Gautam Bisht, LBL
  ! Date: 08/03/12
  ! 

  use Surface_Global_Aux_module
  use Option_module

  implicit none

  type(option_type) :: option
  type(surface_global_auxvar_type) :: surf_global_auxvar_up
  type(surface_global_auxvar_type) :: surf_global_auxvar_dn
  PetscReal :: zc_up, zc_dn
  PetscReal :: mannings_up, mannings_dn
  PetscReal :: head_up, head_dn
  PetscReal :: dist, length
  PetscReal :: vel                      ! units: m/s
  PetscReal :: Res(1:option%nflowdof)   ! units: m^3/s

  PetscReal :: hw_half
  PetscReal :: mannings_half
  PetscReal :: dhead

  ! upwind the total head and Manning's coefficient
  head_up = surf_global_auxvar_up%head(1) + zc_up
  head_dn = surf_global_auxvar_dn%head(1) + zc_dn
  if (head_up > head_dn) then
    hw_half       = surf_global_auxvar_up%head(1)
    mannings_half = mannings_up
  else
    hw_half       = surf_global_auxvar_dn%head(1)
    mannings_half = mannings_dn
  endif
  
  ! compute Manning's velocity
  dhead = head_up - head_dn
  vel   = sign(hw_half**(2.d0/3.d0)/mannings_half*abs(dhead/dist)**0.5d0,dhead) ! [m/s]

  ! compute the volumetric flow rate
  Res(TH_PRESSURE_DOF) = vel*hw_half*length ! [m^3/s]

end subroutine SurfaceFlowFlux

! ************************************************************************** !

subroutine SurfaceFlowBCFlux(ibndtype, &
                           surf_global_auxvar, &
                           slope, &
                           mannings, &
                           length, &
                           option, &
                           vel, &
                           Res)
  ! 
  ! This routine computes the boundary term surface water equation
  ! 
  ! Author: Gautam Bisht, LBL
  ! Date: 08/03/12
  ! 

  use Option_module
  use Surface_Global_Aux_module
  
  implicit none

  type(option_type) :: option
  type(surface_global_auxvar_type) :: surf_global_auxvar
  PetscReal :: slope
  PetscReal :: mannings
  PetscReal :: length
  PetscReal :: flux
  PetscInt :: ibndtype(:)
  PetscReal :: vel
  PetscReal :: Res(1:option%nflowdof) 

  PetscInt :: pressure_bc_type
  PetscReal :: head

  flux = 0.d0
  vel = 0.d0
  
  ! Flow  
  pressure_bc_type = ibndtype(TH_PRESSURE_DOF)
  head = surf_global_auxvar%head(1)
  
  select case(pressure_bc_type)
    case (ZERO_GRADIENT_BC)
      if (slope<0.d0) then
        vel =  0.d0
      else
        vel = -sqrt(dabs(slope))/mannings*((head)**(2.d0/3.d0))
      endif
    case default
      option%io_buffer = 'Uknown pressure_bc_type for surface flow '
  end select
  
  flux = head*vel
  Res(TH_PRESSURE_DOF) = flux*length

end subroutine SurfaceFlowBCFlux

! ************************************************************************** !

subroutine SurfaceFlowUpdateAuxVars(surf_realization)
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
  use Surface_Global_Aux_module

  implicit none

  class(realization_surface_type) :: surf_realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(surface_field_type), pointer :: surf_field
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars(:)
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars_bc(:)
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars_ss(:)

  PetscInt :: ghosted_id, local_id, sum_connection, iconn
  PetscReal, pointer :: xx_loc_p(:), icap_loc_p(:), iphase_loc_p(:)
  PetscReal :: xxbc(1)
  PetscReal :: xxss(1)
  PetscReal :: tsrc1
  PetscErrorCode :: ierr

  option => surf_realization%option
  patch => surf_realization%patch
  grid => patch%grid
  surf_field => surf_realization%surf_field

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
    surf_global_auxvars(ghosted_id)%head(1) = xx_loc_p(ghosted_id)
  enddo
  call VecRestoreArrayF90(surf_field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)
   
  call VecGetArrayF90(surf_field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)
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

        select case(boundary_condition%flow_condition%itype(TH_PRESSURE_DOF))
          case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC,HET_DIRICHLET_BC)
            xxbc(1) = boundary_condition%flow_aux_real_var(TH_PRESSURE_DOF,iconn)
          case(NEUMANN_BC,ZERO_GRADIENT_BC)
            xxbc(1) = xx_loc_p(ghosted_id)
        end select
      
      surf_global_auxvars_bc(sum_connection)%head(1) = xxbc(1)
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

      xxss = xx_loc_p(ghosted_id)
      surf_global_auxvars_ss(sum_connection)%head(1) = xxss(1)
    enddo
    source_sink => source_sink%next
  enddo

  call VecRestoreArrayF90(surf_field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)

end subroutine SurfaceFlowUpdateAuxVars

! ************************************************************************** !

function SurfaceFlowGetTecplotHeader(surf_realization,icolumn)
  ! 
  ! This routine surface flow tecplot file header
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 05/29/12
  ! 

  use Realization_Surface_class
  use Option_module

  implicit none

  character(len=MAXSTRINGLENGTH) :: SurfaceFlowGetTecplotHeader
  class(realization_surface_type) :: surf_realization
  PetscInt :: icolumn

  character(len=MAXSTRINGLENGTH) :: string, string2

  string = ''

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-P [Pa]"'')') icolumn
  else
    write(string2,'('',"P [m]"'')')
  endif
  string = trim(string) // trim(string2)
#ifdef GLENN_NEW_IO
  !call OutputOptionAddPlotVariable(realization%output_option,PRESSURE, &
  !                           ZERO_INTEGER,ZERO_INTEGER)
#endif

  SurfaceFlowGetTecplotHeader = string

end function SurfaceFlowGetTecplotHeader

! ************************************************************************** !

subroutine SurfaceFlowUpdateSurfState(surf_realization)
  ! 
  ! This routine gets updated values of standing water at the end of
  ! subsurface-flow model timestep.
  ! 
  ! Author: Gautam Bisht, LBL
  ! Date: 07/30/13
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

  PetscInt :: iconn
  PetscInt :: local_id
  PetscInt :: sum_connection

  PetscReal :: den
  PetscReal :: dum1
  PetscReal, pointer :: avg_vdarcy_p(:)   ! avg darcy velocity [m/s]
  PetscReal, pointer :: hw_p(:)           ! head [m]
  PetscReal, pointer :: surfpress_p(:)
  PetscErrorCode :: ierr

  PetscBool :: coupler_found = PETSC_FALSE

  option     => surf_realization%option
  surf_field => surf_realization%surf_field
  surf_grid  => surf_realization%discretization%grid
  
  call EOSWaterdensity(option%reference_temperature, &
                       option%reference_pressure,den,dum1,ierr)

  call VecGetArrayF90(surf_field%flow_xx, hw_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(surf_field%press_subsurf, surfpress_p,  &
                      ierr);CHKERRQ(ierr)

  do local_id = 1,surf_grid%nlmax

    hw_p(local_id) = (surfpress_p(local_id)-option%reference_pressure)/ &
                        (abs(option%gravity(3)))/den
    if (hw_p(local_id)<1.d-15) hw_p(local_id) = 0.d0

  enddo
  call VecRestoreArrayF90(surf_field%flow_xx, hw_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(surf_field%press_subsurf, surfpress_p,  &
                          ierr);CHKERRQ(ierr)

  call DiscretizationGlobalToLocal(surf_realization%discretization, &
                                   surf_field%flow_xx, &
                                   surf_field%flow_xx_loc, &
                                   NFLOWDOF)
  call SurfaceFlowUpdateAuxVars(surf_realization)

end subroutine SurfaceFlowUpdateSurfState

! ************************************************************************** !

end module Surface_Flow_module
