module PMC_Surface_class

  use PMC_Base_class
  use Realization_Subsurface_class
  use Realization_Surface_class
  use Timestepper_Surface_class

  use PFLOTRAN_Constants_module

  implicit none

#include "petsc/finclude/petscsys.h"

  private

  type, public, extends(pmc_base_type) :: pmc_surface_type
    class(realization_subsurface_type), pointer :: subsurf_realization
    class(realization_surface_type), pointer :: surf_realization
  contains
    procedure, public :: Init => PMCSurfaceInit
    procedure, public :: RunToTime => PMCSurfaceRunToTime
    procedure, public :: GetAuxData => PMCSurfaceGetAuxData
    procedure, public :: SetAuxData => PMCSurfaceSetAuxData
    procedure, public :: PMCSurfaceGetAuxDataAfterRestart
    procedure, public :: Destroy => PMCSurfaceDestroy
  end type pmc_surface_type

  public :: PMCSurfaceCreate

contains

! ************************************************************************** !

function PMCSurfaceCreate()
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/27/13
  ! 

  implicit none
  
  class(pmc_surface_type), pointer :: PMCSurfaceCreate
  
  class(pmc_surface_type), pointer :: pmc

  print *, 'PMCSurfaceCreate%Create()'
  
  allocate(pmc)
  call pmc%Init()
  
  PMCSurfaceCreate => pmc  
  
end function PMCSurfaceCreate

! ************************************************************************** !

subroutine PMCSurfaceInit(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/27/13
  ! 

  implicit none
  
  class(pmc_surface_type) :: this
  
  print *, 'PMCSurfaceInit%Init()'
  
  call PMCBaseInit(this)
  nullify(this%surf_realization)
  nullify(this%subsurf_realization)
!  nullify(this%surf_timestepper)

end subroutine PMCSurfaceInit

! ************************************************************************** !

recursive subroutine PMCSurfaceRunToTime(this,sync_time,stop_flag)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/27/13
  ! 

#include "petsc/finclude/petscviewer.h"
  use petscsys
  use Timestepper_Base_class
  use Output_Aux_module
  use Output_module, only : Output
  use Realization_Subsurface_class, only : realization_subsurface_type
  use PM_Base_class
  use PM_Surface_Flow_class
  use Option_module
  use Surface_Flow_module
  use Surface_TH_module
  use Output_Surface_module
  use Checkpoint_module
  
  implicit none
  
  class(pmc_surface_type), target :: this
  PetscReal :: sync_time
  PetscInt :: stop_flag
  character(len=MAXSTRINGLENGTH) :: filename_append
  class(pmc_base_type), pointer :: pmc_base
  PetscInt :: local_stop_flag
  PetscBool :: failure
  PetscBool :: snapshot_plot_flag
  PetscBool :: observation_plot_flag
  PetscBool :: massbal_plot_flag
  PetscBool :: checkpoint_at_this_time_flag
  PetscBool :: checkpoint_at_this_timestep_flag
  PetscBool :: peer_already_run_to_time
  class(pm_base_type), pointer :: cur_pm
  PetscReal :: dt_max_loc
  PetscReal :: dt_max_glb
  PetscViewer :: viewer
  PetscErrorCode :: ierr

  if (stop_flag == TS_STOP_FAILURE) return
  
  this%option%io_buffer = trim(this%name) // ':' // trim(this%pm_list%name)
  call printVerboseMsg(this%option)
  
  ! Get data of other process-model
  if (this%option%restart_flag .and. this%option%first_step_after_restart) then
    this%option%first_step_after_restart = PETSC_FALSE
  else
    call this%GetAuxData()
  endif

  local_stop_flag = TS_CONTINUE
  do
    if (local_stop_flag /= TS_CONTINUE) exit ! end simulation
    if (this%timestepper%target_time >= sync_time) exit
    
    call SetOutputFlags(this)
    snapshot_plot_flag = PETSC_FALSE
    observation_plot_flag = PETSC_FALSE
    massbal_plot_flag = PETSC_FALSE
    checkpoint_at_this_time_flag = PETSC_FALSE
    checkpoint_at_this_timestep_flag = PETSC_FALSE
    
    cur_pm => this%pm_list

    select case(this%option%iflowmode)
      case (TH_MODE)
        call SurfaceTHComputeMaxDt(this%surf_realization,dt_max_loc)
    end select

    ! Find mininum allowable timestep across all processors
    call MPI_Allreduce(dt_max_loc,dt_max_glb,ONE_INTEGER_MPI, &
                       MPI_DOUBLE_PRECISION,MPI_MIN,this%option%mycomm,ierr)
    select type(timestepper => this%timestepper)
      class is(timestepper_surface_type)
        timestepper%dt_max_allowable = dt_max_glb
        timestepper%surf_subsurf_coupling_flow_dt = &
          this%option%surf_subsurf_coupling_flow_dt
    end select
    call this%timestepper%SetTargetTime(sync_time,this%option, &
                                        local_stop_flag,snapshot_plot_flag, &
                                        observation_plot_flag, &
                                        massbal_plot_flag, &
                                        checkpoint_at_this_time_flag)

    this%option%surf_flow_dt = this%timestepper%dt

    ! Accumulate data needed by process-model
    call this%AccumulateAuxData()

    call this%timestepper%StepDT(this%pm_list,local_stop_flag)

    if (local_stop_flag  == TS_STOP_FAILURE) exit ! failure
    ! Have to loop over all process models coupled in this object and update
    ! the time step size.  Still need code to force all process models to
    ! use the same time step size if tightly or iteratively coupled.
    cur_pm => this%pm_list
    do
      if (.not.associated(cur_pm)) exit
      ! have to update option%time for conditions
      this%option%time = this%timestepper%target_time
      call cur_pm%UpdateSolution()
      !! TODO(gb)
      !!!call this%timestepper%UpdateDT(cur_pm)
      cur_pm => cur_pm%next
    enddo

#if 0
    ! Run underlying process model couplers
    if (associated(this%child)) then
      call this%child%RunToTime(this%timestepper%target_time,local_stop_flag)
    endif
#endif
    
    ! only print output for process models of depth 0
    ! TODO(GB): Modify OutputSurface()
    !if (associated(this%Output)) then
      if (this%timestepper%time_step_cut_flag) then
        snapshot_plot_flag = PETSC_FALSE
      endif
      ! however, if we are using the modulus of the output_option%imod, we may
      ! still print
      if (mod(this%timestepper%steps,this%pm_list% &
              output_option%periodic_snap_output_ts_imod) == 0) then
        snapshot_plot_flag = PETSC_TRUE
      endif
      if (mod(this%timestepper%steps,this%pm_list%output_option% &
              periodic_obs_output_ts_imod) == 0) then
        observation_plot_flag = PETSC_TRUE
      endif
      if (mod(this%timestepper%steps,this%pm_list%output_option% &
              periodic_msbl_output_ts_imod) == 0) then
        massbal_plot_flag = PETSC_TRUE
      endif
      !call this%Output(this%pm_list%realization_base,snapshot_plot_flag, &
      !                 observation_plot_flag, massbal_plot_flag)
      call OutputSurface(this%surf_realization, this%subsurf_realization, &
                         snapshot_plot_flag, observation_plot_flag, &
                         massbal_plot_flag)
    !endif

    if (this%is_master .and. associated(this%checkpoint_option)) then
      if (this%checkpoint_option%periodic_ts_incr > 0 .and. &
          mod(this%timestepper%steps, &
              this%checkpoint_option%periodic_ts_incr) == 0) then
        checkpoint_at_this_timestep_flag = PETSC_TRUE
      endif
    endif

    peer_already_run_to_time = PETSC_FALSE
    if (checkpoint_at_this_time_flag .or. &
        checkpoint_at_this_timestep_flag) then
      ! if checkpointing, need to sync all other PMCs.  Those "below" are
      ! already in sync, but not those "next".
      ! Set data needed by process-model
      call this%SetAuxData()
      ! Run neighboring process model couplers
      if (associated(this%peer)) then
        call this%peer%RunToTime(this%timestepper%target_time,local_stop_flag)
        peer_already_run_to_time = PETSC_TRUE
      endif
      call this%GetAuxData()
      ! it is possible that two identical checkpoint files will be created,
      ! one at the time and another at the time step, but this is fine.
      if (checkpoint_at_this_time_flag) then
        filename_append = &
          CheckpointAppendNameAtTime(this%checkpoint_option, &
                                     this%option%time, &
                                     this%option)
        call this%Checkpoint(filename_append)
      endif
      if (checkpoint_at_this_timestep_flag) then
        filename_append = &
          CheckpointAppendNameAtTimestep(this%checkpoint_option, &
                                         this%timestepper%steps, &
                                         this%option)
        call this%Checkpoint(filename_append)
      endif
    endif                         
                         
  enddo
  
  this%option%surf_flow_time = this%timestepper%target_time

  ! Set data needed by process-model
  call this%SetAuxData()

  ! Run neighboring process model couplers
  if (associated(this%peer) .and. .not.peer_already_run_to_time) then
    call this%peer%RunToTime(sync_time,local_stop_flag)
  endif

  stop_flag = max(stop_flag,local_stop_flag)
  
end subroutine PMCSurfaceRunToTime

! ************************************************************************** !

subroutine PMCSurfaceGetAuxData(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 08/21/13
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Surface_Flow_module
  use Surface_TH_module
  use Surface_TH_module
  use Option_module

  implicit none

  class(pmc_surface_type) :: this

  PetscErrorCode :: ierr

#ifdef DEBUG
  print *, 'PMCSurfaceGetAuxData()'
#endif

  if (this%option%subsurf_surf_coupling == SEQ_COUPLED) then
    select type(pmc => this)
      class is(pmc_surface_type)
        select case(this%option%iflowmode)
          case (TH_MODE)
            call VecScatterBegin(pmc%sim_aux%subsurf_to_surf, &
                                 pmc%sim_aux%subsurf_pres_top_bc, &
                                 pmc%surf_realization%surf_field%press_subsurf, &
                                 INSERT_VALUES,SCATTER_FORWARD, &
                                 ierr);CHKERRQ(ierr)
            call VecScatterEnd(pmc%sim_aux%subsurf_to_surf, &
                               pmc%sim_aux%subsurf_pres_top_bc, &
                               pmc%surf_realization%surf_field%press_subsurf, &
                               INSERT_VALUES,SCATTER_FORWARD, &
                               ierr);CHKERRQ(ierr)
            call VecScatterBegin(pmc%sim_aux%subsurf_to_surf, &
                                 pmc%sim_aux%subsurf_temp_top_bc, &
                                 pmc%surf_realization%surf_field%temp_subsurf, &
                                 INSERT_VALUES,SCATTER_FORWARD, &
                                 ierr);CHKERRQ(ierr)
            call VecScatterEnd(pmc%sim_aux%subsurf_to_surf, &
                               pmc%sim_aux%subsurf_temp_top_bc, &
                               pmc%surf_realization%surf_field%temp_subsurf, &
                               INSERT_VALUES,SCATTER_FORWARD, &
                               ierr);CHKERRQ(ierr)
            call SurfaceTHUpdateSurfState(pmc%surf_realization)
        end select
    end select
  endif

end subroutine PMCSurfaceGetAuxData

! ************************************************************************** !

subroutine PMCSurfaceSetAuxData(this)
  ! 
  ! This routine extracts data from surface flow model and stores it sim-aux,
  ! which will be required by the subsurface flow model.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 08/21/13
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Connection_module
  use Coupler_module
  use Grid_module
  use Option_module
  use Patch_module
  use Surface_Global_Aux_module
  use Surface_Flow_module
  use Surface_TH_module
  use Surface_TH_Aux_module
  use Realization_Surface_class
  use String_module

  implicit none

  class(pmc_surface_type) :: this

  type(grid_type), pointer :: surf_grid
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars(:)
  type(Surface_TH_auxvar_type), pointer :: surf_auxvars(:)
  type(patch_type), pointer :: surf_patch
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  class(realization_surface_type), pointer :: surf_realization

  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: iend
  PetscInt :: istart
  PetscInt :: iconn

  PetscReal :: dt
  PetscReal, pointer :: xx_loc_p(:)
  PetscReal, pointer :: surf_head_p(:)
  PetscReal, pointer :: surf_temp_p(:)
  PetscReal, pointer :: surf_hflux_p(:)
  PetscBool :: found
  PetscReal :: esrc
  PetscReal :: atm_temp
  PetscErrorCode :: ierr

  dt = this%option%surf_subsurf_coupling_flow_dt

  if (this%option%subsurf_surf_coupling == SEQ_COUPLED) then
    select type(pmc => this)
      class is(pmc_surface_type)

        select case(this%option%iflowmode)

          case (TH_MODE)

            surf_realization => pmc%surf_realization
            surf_patch => surf_realization%patch
            surf_grid => surf_patch%grid
            surf_global_auxvars => surf_patch%surf_aux%SurfaceGlobal%auxvars
            surf_auxvars => surf_patch%surf_aux%SurfaceTH%auxvars

            call VecGetArrayF90(pmc%surf_realization%surf_field%flow_xx_loc, &
                                xx_loc_p,ierr);CHKERRQ(ierr)
            call VecGetArrayF90(pmc%sim_aux%surf_head, surf_head_p,  &
                                ierr);CHKERRQ(ierr)
            call VecGetArrayF90(pmc%sim_aux%surf_temp, surf_temp_p,  &
                                ierr);CHKERRQ(ierr)
            call VecGetArrayF90(pmc%sim_aux%surf_hflux_exchange_with_subsurf, &
                                surf_hflux_p, ierr);CHKERRQ(ierr)

            do ghosted_id = 1, surf_grid%ngmax
              local_id = surf_grid%nG2L(ghosted_id)
              if (local_id < 1) cycle
              iend = local_id*this%option%nflowdof
              istart = iend - this%option%nflowdof+1
              if (xx_loc_p(istart) < 1.d-8) then
                surf_head_p(local_id) = 0.d0
                surf_temp_p(local_id) = this%option%reference_temperature
              else
                surf_head_p(local_id) = xx_loc_p(istart)
                surf_temp_p(local_id) = surf_global_auxvars(ghosted_id)%temp
              endif
            enddo

            found = PETSC_FALSE
            source_sink => surf_patch%source_sink_list%first
            do
              if (.not.associated(source_sink)) exit

              if (associated(source_sink%flow_aux_real_var)) then
                cur_connection_set => source_sink%connection_set

                if (StringCompare(source_sink%name,'atm_energy_ss') .or. &
                    StringCompare(source_sink%name,'clm_energy_srf_ss')) then

                  do iconn = 1, cur_connection_set%num_connections

                    local_id = cur_connection_set%id_dn(iconn)
                    select case(source_sink%flow_condition% &
                                  itype(TH_TEMPERATURE_DOF))
                      case (ENERGY_RATE_SS)
                        esrc = source_sink%flow_condition%energy_rate% &
                                  dataset%rarray(1)
                      case (HET_ENERGY_RATE_SS)
                        esrc = source_sink%flow_aux_real_var(TWO_INTEGER,iconn)
                      case (DIRICHLET_BC)
                        esrc = source_sink%flow_condition%temperature% &
                                  dataset%rarray(1)
                      case (HET_DIRICHLET_BC)
                        esrc = source_sink%flow_aux_real_var(TWO_INTEGER,iconn)
                      case default
                        this%option%io_buffer = 'atm_energy_ss does not have &
                          &a temperature condition that is either a &
                          &ENERGY_RATE_SS/HET_ENERGY_RATE_SS/DIRICHLET_BC/ &
                          &HET_DIRICHLET_BC'
                        call printErrMsg(this%option)
                    end select

                    ! Only when no standing water is present, the atmospheric
                    ! energy flux is applied directly on subsurface domain.
                    if (surf_head_p(local_id) < 1.d-8) then
                      surf_hflux_p(local_id) = esrc
                    else
                      surf_hflux_p(local_id) = 0.d0
                    endif

                  enddo

                  found = PETSC_TRUE

                endif ! StringCompare()
              endif ! associate()

              source_sink => source_sink%next
            enddo

            call VecRestoreArrayF90(pmc%surf_realization%surf_field%flow_xx_loc, &
                                    xx_loc_p,ierr);CHKERRQ(ierr)
            call VecRestoreArrayF90(pmc%sim_aux%surf_head, surf_head_p, &
                                    ierr);CHKERRQ(ierr)
            call VecRestoreArrayF90(pmc%sim_aux%surf_temp, surf_temp_p, &
                                    ierr);CHKERRQ(ierr)
            call VecRestoreArrayF90(pmc%sim_aux%surf_hflux_exchange_with_subsurf, &
                                surf_hflux_p, ierr);CHKERRQ(ierr)

            if (.not.(found)) then
              this%option%io_buffer = 'atm_energy_ss/clm_energy_srf_ss not ' // &
                'found in surface-flow model'
              call printErrMsg(this%option)
            endif
        end select
    end select
  endif

end subroutine PMCSurfaceSetAuxData

! ************************************************************************** !

subroutine PMCSurfaceGetAuxDataAfterRestart(this)
  ! 
  ! This routine is called to set values in sim_aux PETSc vectors after restart
  ! checkpoint files is read.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 09/23/13
  ! 
#include "petsc/finclude/petscvec.h"
  use petscvec
  use Surface_Flow_module
  use Surface_TH_Aux_module
  use Surface_TH_module
  use Option_module
  use EOS_Water_module

  implicit none

  class(pmc_surface_type) :: this

  PetscInt :: ghosted_id
  PetscInt :: local_id
  PetscInt :: count
  PetscReal, pointer :: xx_p(:)
  PetscReal, pointer :: surfpress_p(:)
  PetscReal, pointer :: surftemp_p(:)
  PetscInt :: istart, iend
  PetscReal :: den
  PetscReal :: dum1
  PetscErrorCode :: ierr
  type(Surface_TH_auxvar_type), pointer :: surf_auxvars(:)

  print *, 'PMCSurfaceGetAuxDataAfterRestart()'

  if (this%option%subsurf_surf_coupling == SEQ_COUPLED) then
    select type(pmc => this)
      class is(pmc_surface_type)
        select case(this%option%iflowmode)

          case (TH_MODE)

            ! NOTE(GB:) This is strictly not correct since density should be
            ! computed based on surface-water temperature (not on
            ! reference-temperature). Presently, SurfaceCheckpointProcessModel()
            ! does not output surface-water temperature for TH-Mode and the
            ! subroutine needs to be modified in future.
            call EOSWaterdensity(this%option%reference_temperature, &
                                 this%option%reference_pressure,den,dum1,ierr)

            surf_auxvars => pmc%surf_realization%patch%surf_aux%SurfaceTH%auxvars

            call VecGetArrayF90(pmc%surf_realization%surf_field%flow_xx, xx_p,  &
                                ierr);CHKERRQ(ierr)
            call VecGetArrayF90(pmc%surf_realization%surf_field%press_subsurf, surfpress_p,  &
                                ierr);CHKERRQ(ierr)
            call VecGetArrayF90(pmc%surf_realization%surf_field%temp_subsurf, surftemp_p,  &
                                ierr);CHKERRQ(ierr)

            count = 0
            do ghosted_id = 1, pmc%surf_realization%discretization%grid%ngmax

              local_id = pmc%surf_realization%discretization%grid%nG2L(ghosted_id)
              if (local_id <= 0) cycle

              count = count + 1
              iend = ghosted_id*this%option%nflowdof
              istart = iend - this%option%nflowdof+1
              surfpress_p(count) = xx_p(istart)*den*abs(this%option%gravity(3)) + &
                                   this%option%reference_pressure
              surftemp_p = xx_p(iend)/xx_p(istart)/den/ &
                      surf_auxvars(ghosted_id)%Cwi - 273.15d0
            enddo
            call VecRestoreArrayF90(pmc%surf_realization%surf_field%flow_xx, xx_p,  &
                                    ierr);CHKERRQ(ierr)
            call VecRestoreArrayF90(pmc%surf_realization%surf_field%press_subsurf, surfpress_p,  &
                                    ierr);CHKERRQ(ierr)
            call VecRestoreArrayF90(pmc%surf_realization%surf_field%temp_subsurf, surftemp_p,  &
                                    ierr);CHKERRQ(ierr)

            call VecScatterBegin(pmc%sim_aux%subsurf_to_surf, &
                                 pmc%sim_aux%subsurf_pres_top_bc, &
                                 pmc%surf_realization%surf_field%press_subsurf, &
                                 INSERT_VALUES,SCATTER_FORWARD, &
                                 ierr);CHKERRQ(ierr)
            call VecScatterEnd(pmc%sim_aux%subsurf_to_surf, &
                               pmc%sim_aux%subsurf_pres_top_bc, &
                               pmc%surf_realization%surf_field%press_subsurf, &
                               INSERT_VALUES,SCATTER_FORWARD, &
                               ierr);CHKERRQ(ierr)
            call VecScatterBegin(pmc%sim_aux%subsurf_to_surf, &
                                 pmc%sim_aux%subsurf_temp_top_bc, &
                                 pmc%surf_realization%surf_field%temp_subsurf, &
                                 INSERT_VALUES,SCATTER_FORWARD, &
                                 ierr);CHKERRQ(ierr)
            call VecScatterEnd(pmc%sim_aux%subsurf_to_surf, &
                               pmc%sim_aux%subsurf_temp_top_bc, &
                               pmc%surf_realization%surf_field%temp_subsurf, &
                               INSERT_VALUES,SCATTER_FORWARD, &
                               ierr);CHKERRQ(ierr)
        end select
    end select
  endif

end subroutine PMCSurfaceGetAuxDataAfterRestart

! ************************************************************************** !

recursive subroutine PMCSurfaceFinalizeRun(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/27/13
  ! 

  use Option_module
  
  implicit none
  
  class(pmc_surface_type), target :: this
  
  call printMsg(this%option,'PMCSurface%FinalizeRun()')
  
  nullify(this%surf_realization)
!  nullify(this%surf_timestepper)
  
end subroutine PMCSurfaceFinalizeRun

! ************************************************************************** !

subroutine PMCSurfaceStrip(this)
  !
  ! Deallocates members of PMC Surface.
  !
  ! Author: Glenn Hammond
  ! Date: 12/02/14
  
  implicit none
  
  class(pmc_surface_type) :: this

  call PMCBaseStrip(this)
  ! realizations destroyed elsewhere
  nullify(this%subsurf_realization)
  nullify(this%surf_realization)

end subroutine PMCSurfaceStrip

! ************************************************************************** !

recursive subroutine PMCSurfaceDestroy(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/02/14
  ! 
  use Option_module

  implicit none
  
  class(pmc_surface_type) :: this
  
#ifdef DEBUG
  call printMsg(this%option,'PMCSurface%Destroy()')
#endif

  if (associated(this%child)) then
    call this%child%Destroy()
    ! destroy does not currently destroy; it strips
    deallocate(this%child)
    nullify(this%child)
  endif 
  
  if (associated(this%peer)) then
    call this%peer%Destroy()
    ! destroy does not currently destroy; it strips
    deallocate(this%peer)
    nullify(this%peer)
  endif
  
  call PMCSurfaceStrip(this)
  
end subroutine PMCSurfaceDestroy

end module PMC_Surface_class
