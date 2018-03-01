
module heat_transport_1D_problem

  implicit none

#include <petsc/finclude/petsc.h>
  PetscInt  , parameter :: nx       = 1
  PetscInt  , parameter :: ny       = 1
  PetscReal , parameter :: x_column = 1.d0
  PetscReal , parameter :: y_column = 1.d0
  PetscReal , parameter :: z_column = 1.d0
  PetscInt              :: nz
  PetscInt              :: ncells_local
  PetscInt              :: ncells_ghost

  public :: run_heat_transport_1D_problem
  public :: output_regression_heat_transport_1D_problem

contains

  !------------------------------------------------------------------------
  subroutine run_heat_transport_1D_problem()

    !
#include <petsc/finclude/petsc.h>
    use MultiPhysicsProbThermalEnthalpy , only : thermal_enthalpy_mpp
    use mpp_varpar                      , only : mpp_varpar_init
    use petscsys
    use petscvec
    use petscmat
    use petscts
    use petscdm
    use petscdmda
    !
    implicit none
    !
    !
    PetscBool          :: converged
    PetscInt           :: converged_reason
    PetscErrorCode     :: ierr
    PetscReal          :: dtime
    PetscInt           :: istep, nstep
    PetscBool          :: flg
    PetscBool          :: save_initial_soln, save_final_soln
    character(len=256) :: string
    character(len=256) :: output_suffix
    PetscViewer        :: viewer

    ! Set default settings
    nz                = 1
    nz                = 100
    dtime             = 3600.d0
    nstep             = 2
    save_initial_soln = PETSC_FALSE
    save_final_soln   = PETSC_FALSE
    output_suffix     = ''

    ! Get some command line options

    call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-nz',nz,flg,ierr)
    call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-dt',dtime,flg,ierr)
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-nstep',nstep,flg,ierr)
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-save_initial_soln',save_initial_soln,flg,ierr)
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-save_final_soln',save_final_soln,flg,ierr)
    call PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-output_suffix',output_suffix,flg,ierr)

    ! Initialize the problem
    call init()

    if (save_initial_soln) then
       if (len(trim(adjustl(output_suffix))) ==  0) then
          string = trim(output_suffix) // 'initial_soln.bin'
       else
          string = 'initial_soln_' // trim(output_suffix) // '.bin'
       endif

       call PetscViewerBinaryOpen(PETSC_COMM_SELF,trim(string),FILE_MODE_WRITE,viewer,ierr);CHKERRQ(ierr)
       call VecView(thermal_enthalpy_mpp%soe%solver%soln,viewer,ierr);CHKERRQ(ierr)
       call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
    endif

    do istep = 1, nstep
       ! Update BC
       call set_bondary_conditions()

       ! Run the model
       call thermal_enthalpy_mpp%soe%StepDT(dtime, istep, &
            converged, converged_reason, ierr); CHKERRQ(ierr)
    enddo

    if (save_final_soln) then
       if (len(trim(adjustl(output_suffix))) ==  0) then
          string = trim(output_suffix) // 'final_soln.bin'
       else
          string = 'final_soln_' // trim(output_suffix) // '.bin'
       endif
       call PetscViewerBinaryOpen(PETSC_COMM_SELF,trim(string),FILE_MODE_WRITE,viewer,ierr);CHKERRQ(ierr)
       call VecView(thermal_enthalpy_mpp%soe%solver%soln,viewer,ierr);CHKERRQ(ierr)
       call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
    endif

  end subroutine run_heat_transport_1D_problem

  !------------------------------------------------------------------------
  subroutine Init()
    !
    use MultiPhysicsProbThermalEnthalpy, only : thermal_enthalpy_mpp
    !
    implicit none
    !

    ! 1. Initialize the multi-physics-problem (MPP)
    call initialize_mpp()

    ! 2. Add all meshes needed for the MPP
    call add_meshes()

    ! 3. Add all governing equations
    call add_goveqns()

    ! 4. Add boundary and source-sink conditions to all governing equations
    call add_conditions_to_goveqns()

    ! 5. Allocate memory to hold auxvars
    call allocate_auxvars()

    ! 6. Setup the MPP
    call thermal_enthalpy_mpp%SetupProblem()

    ! 7. Add material properities associated with all governing equations
    call set_material_properties() 

    ! 8. Set initial conditions
    call set_initial_conditions()

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine initialize_mpp()
    !
    ! !DESCRIPTION:
    ! Initialization VSFM
    !
#include <petsc/finclude/petsc.h>
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : MPP_THERMAL_EBASED_SNES_CLM
    use MultiPhysicsProbThermalEnthalpy, only : thermal_enthalpy_mpp
    use petscsys
    !
    ! !ARGUMENTS
    implicit none
    !
    PetscInt       :: iam
    PetscErrorCode :: ierr

    call MPI_Comm_rank(MPI_COMM_WORLD, iam, ierr)

    !
    ! Set up the multi-physics problem
    !
    call thermal_enthalpy_mpp%Init       ()
    call thermal_enthalpy_mpp%SetName    ('1D heat conduction')
    call thermal_enthalpy_mpp%SetID      (MPP_THERMAL_EBASED_SNES_CLM)
    call thermal_enthalpy_mpp%SetMPIRank (iam)

  end subroutine initialize_mpp

  !------------------------------------------------------------------------
  subroutine add_meshes()
    !
#include <petsc/finclude/petsc.h>
    !
    use MultiPhysicsProbThermalEnthalpy , only : thermal_enthalpy_mpp
    use MultiPhysicsProbConstants       , only : MESH_ALONG_GRAVITY
    use MultiPhysicsProbConstants       , only : MESH_AGAINST_GRAVITY
    use MultiPhysicsProbConstants       , only : MESH_CLM_THERMAL_SOIL_COL
    use MultiPhysicsProbConstants       , only : VAR_XC
    use MultiPhysicsProbConstants       , only : VAR_YC
    use MultiPhysicsProbConstants       , only : VAR_ZC
    use MultiPhysicsProbConstants       , only : VAR_DX
    use MultiPhysicsProbConstants       , only : VAR_DY
    use MultiPhysicsProbConstants       , only : VAR_DZ
    use MultiPhysicsProbConstants       , only : VAR_AREA
    use MultiPhysicsProbConstants       , only : CONN_SET_INTERNAL
    use MultiPhysicsProbConstants       , only : CONN_SET_LATERAL
    use MultiPhysicsProbConstants       , only : CONN_VERTICAL
    use mpp_varpar                      , only : mpp_varpar_set_nlevsoi, mpp_varpar_set_nlevgrnd
    use petscsys
    !
    implicit none
    !
    PetscReal          :: dx, dy, dz
    PetscInt           :: imesh, kk
    PetscInt           :: nlev
    PetscInt           :: iconn, vert_nconn
    PetscReal, pointer :: soil_xc(:)           ! x-position of grid cell [m]
    PetscReal, pointer :: soil_yc(:)           ! y-position of grid cell [m]
    PetscReal, pointer :: soil_zc(:)           ! z-position of grid cell [m]
    PetscReal, pointer :: soil_dx(:)           ! layer thickness of grid cell [m]
    PetscReal, pointer :: soil_dy(:)           ! layer thickness of grid cell [m]
    PetscReal, pointer :: soil_dz(:)           ! layer thickness of grid cell [m]
    PetscReal, pointer :: soil_area(:)         ! area of grid cell [m^2]
    PetscInt , pointer :: soil_filter(:)       ! 

    PetscInt, pointer  :: vert_conn_id_up(:)   !
    PetscInt, pointer  :: vert_conn_id_dn(:)   !
    PetscReal, pointer :: vert_conn_dist_up(:) !
    PetscReal, pointer :: vert_conn_dist_dn(:) !
    PetscReal, pointer :: vert_conn_area(:)    !
    PetscInt , pointer :: vert_conn_type(:)    !

    PetscErrorCode     :: ierr

    call mpp_varpar_set_nlevsoi(nz)
    call mpp_varpar_set_nlevgrnd(nz)

    dx = x_column/nx
    dy = y_column/ny
    dz = z_column/nz

    imesh        = 1
    nlev         = nz
    ncells_local = nx*ny*nz
    ncells_ghost = 0

    allocate(soil_xc(nz))
    allocate(soil_yc(nz))
    allocate(soil_zc(nz))
    allocate(soil_dx(nz))
    allocate(soil_dy(nz))
    allocate(soil_dz(nz))
    allocate(soil_area(nz))
    allocate(soil_filter(nz))

    soil_filter (:) = 1
    soil_area   (:) = dx*dy
    soil_dx     (:) = dx
    soil_dy     (:) = dy
    soil_dz     (:) = dz
    soil_xc     (:) = dx/2.d0
    soil_xc     (:) = dy/2.d0

    do kk = 1,nz
       soil_zc(kk) = dz/2.d0 + dz * (kk - 1)
    enddo

    !
    ! Set up the meshes
    !    
    call thermal_enthalpy_mpp%SetNumMeshes(1)

    call thermal_enthalpy_mpp%MeshSetName        (imesh, 'Soil mesh')
    call thermal_enthalpy_mpp%MeshSetOrientation (imesh, MESH_AGAINST_GRAVITY)
    call thermal_enthalpy_mpp%MeshSetID          (imesh, MESH_CLM_THERMAL_SOIL_COL)
    call thermal_enthalpy_mpp%MeshSetDimensions  (imesh, ncells_local, ncells_ghost, nlev)

    call thermal_enthalpy_mpp%MeshSetGridCellFilter      (imesh, soil_filter)
    call thermal_enthalpy_mpp%MeshSetGeometricAttributes (imesh, VAR_XC   , soil_xc)
    call thermal_enthalpy_mpp%MeshSetGeometricAttributes (imesh, VAR_YC   , soil_yc)
    call thermal_enthalpy_mpp%MeshSetGeometricAttributes (imesh, VAR_ZC   , soil_zc)
    call thermal_enthalpy_mpp%MeshSetGeometricAttributes (imesh, VAR_DX   , soil_dx)
    call thermal_enthalpy_mpp%MeshSetGeometricAttributes (imesh, VAR_DY   , soil_dy)
    call thermal_enthalpy_mpp%MeshSetGeometricAttributes (imesh, VAR_DZ   , soil_dz)
    call thermal_enthalpy_mpp%MeshSetGeometricAttributes (imesh, VAR_AREA , soil_area)
    call thermal_enthalpy_mpp%MeshComputeVolume          (imesh)

    allocate (vert_conn_id_up   (nz-1))
    allocate (vert_conn_id_dn   (nz-1))
    allocate (vert_conn_dist_up (nz-1))
    allocate (vert_conn_dist_dn (nz-1))
    allocate (vert_conn_area    (nz-1))
    allocate (vert_conn_type    (nz-1))

    iconn = 0
    do kk = 1, nz-1

       iconn = iconn + 1
       vert_conn_id_up(iconn)   = kk
       vert_conn_id_dn(iconn)   = vert_conn_id_up(iconn) + 1
       vert_conn_dist_up(iconn) = 0.5d0*dz
       vert_conn_dist_dn(iconn) = 0.5d0*dz
       vert_conn_area(iconn)    = soil_area(kk)
       vert_conn_type(iconn)    = CONN_VERTICAL

    end do
    vert_nconn = iconn

    call thermal_enthalpy_mpp%MeshSetConnectionSet(imesh, CONN_SET_INTERNAL, &
         vert_nconn,  vert_conn_id_up, vert_conn_id_dn,                      &
         vert_conn_dist_up, vert_conn_dist_dn,  vert_conn_area,              &
         vert_conn_type)

    deallocate(soil_xc)
    deallocate(soil_yc)
    deallocate(soil_zc)
    deallocate(soil_dx)
    deallocate(soil_dy)
    deallocate(soil_dz)
    deallocate(soil_area)
    deallocate(soil_filter)  

  end subroutine add_meshes

  !------------------------------------------------------------------------

  subroutine add_goveqns()
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbThermalEnthalpy , only : thermal_enthalpy_mpp
    use MultiPhysicsProbConstants       , only : GE_THERM_SOIL_EBASED
    use MultiPhysicsProbConstants       , only : MESH_CLM_THERMAL_SOIL_COL
    !
    ! !ARGUMENTS
    implicit none

    call thermal_enthalpy_mpp%AddGovEqn(GE_THERM_SOIL_EBASED, 'Heat transport based on enthalpy ODE', &
         MESH_CLM_THERMAL_SOIL_COL)

    call thermal_enthalpy_mpp%SetMeshesOfGoveqns()

  end subroutine add_goveqns

  !------------------------------------------------------------------------
  subroutine add_conditions_to_goveqns()
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbThermalEnthalpy , only : thermal_enthalpy_mpp
    use MultiPhysicsProbConstants , only : SOIL_CELLS
    use MultiPhysicsProbConstants , only : SOIL_TOP_CELLS
    use MultiPhysicsProbConstants , only : SOIL_BOTTOM_CELLS
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_DIRICHLET
    !
    ! !ARGUMENTS
    implicit none
    !
    PetscInt :: ieqn

    ieqn = 1

    call thermal_enthalpy_mpp%soe%AddConditionInGovEqn(ieqn, COND_BC,   &
         'Constant temperature condition at top', 'K', COND_DIRICHLET, &
         SOIL_TOP_CELLS)

    call thermal_enthalpy_mpp%soe%AddConditionInGovEqn(ieqn, COND_BC,   &
         'Constant temperature condition at bottom', 'K', COND_DIRICHLET, &
         SOIL_BOTTOM_CELLS)

  end subroutine add_conditions_to_goveqns

  !------------------------------------------------------------------------
  subroutine allocate_auxvars()
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbThermalEnthalpy , only : thermal_enthalpy_mpp
    !
    implicit none

    !
    ! Allocate auxvars
    !
    call thermal_enthalpy_mpp%AllocateAuxVars()

  end subroutine allocate_auxvars

  !------------------------------------------------------------------------
  subroutine set_material_properties()
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbThermalEnthalpy , only : thermal_enthalpy_mpp
    use MultiPhysicsProbThermalEnthalpy , only : MPPThermalSetSoils
    use MultiPhysicsProbConstants       , only : GRAVITY_CONSTANT
    use EOSWaterMod                     , only : DENSITY_TGDPB01, DENSITY_CONSTANT
    use EOSWaterMod                     , only : INT_ENERGY_ENTHALPY_CONSTANT
    use mpp_varcon                      , only : denh2o
    use mpp_varcon                      , only : grav
    !
    implicit none
    !
    PetscReal        , parameter :: porosity            = 0.368d0
    PetscReal        , parameter :: lambda              = 0.5d0
    PetscReal        , parameter :: alpha               = 3.4257d-4
    PetscReal        , parameter :: perm                = 8.3913d-12
    !
    PetscReal        , pointer   :: watsat(:,:)
    PetscReal        , pointer   :: hksat(:,:)
    PetscReal        , pointer   :: bsw(:,:)
    PetscReal        , pointer   :: sucsat(:,:)
    PetscReal        , pointer   :: eff_porosity(:,:)
    PetscReal        , pointer   :: residual_sat(:,:)
    PetscInt         , pointer   :: lun_type(:)
    PetscReal        , pointer   :: csol(:,:)
    PetscReal        , pointer   :: tkmg(:,:)
    PetscReal        , pointer   :: tkdry(:,:)
    PetscReal        , parameter :: vish2o = 0.001002d0    ! [N s/m^2] @ 20 degC
    PetscInt                     :: begc , endc
    integer          , pointer   :: filter(:)
    character(len=32)            :: satfunc_type
    !-----------------------------------------------------------------------

    begc = 1
    endc = 1

    satfunc_type = 'van_genuchten'

    allocate(watsat       (1,nz))
    allocate(hksat        (1,nz))
    allocate(bsw          (1,nz))
    allocate(sucsat       (1,nz))
    allocate(eff_porosity (1,nz))
    allocate(residual_sat (1,nz))
    allocate(filter       (nz))
    allocate(lun_type     (nz))
    allocate(csol         (1,nz))
    allocate(tkmg         (1,nz))
    allocate(tkdry        (1,nz))

    ! Soil properties
    filter(:)         = 1
    watsat(:,:)       = porosity
    hksat(:,:)        = perm /vish2o * (denh2o * grav) / 0.001d0
    bsw(:,:)          = 1.d0/lambda
    sucsat(:,:)       = 1.d0/(alpha*GRAVITY_CONSTANT)
    eff_porosity(:,:) = 0.d0
    residual_sat(:,:) = 0.2772d0

    lun_type(:)       = 1
    csol(:,:)         = 837.d0
    tkmg(:,:)         = 0.5d0
    tkdry(:,:)        = 0.25d0

    call MPPThermalSetSoils(thermal_enthalpy_mpp, begc, endc, filter, &
         watsat, csol, tkdry,                           &
         hksat, bsw, sucsat, residual_sat,                &
         satfunc_type, DENSITY_CONSTANT, INT_ENERGY_ENTHALPY_CONSTANT)

    deallocate(watsat       )
    deallocate(hksat        )
    deallocate(bsw          )
    deallocate(sucsat       )
    deallocate(eff_porosity )
    deallocate(filter       )
    deallocate(residual_sat )
    deallocate(lun_type     )
    deallocate(csol         )
    deallocate(tkmg         )
    deallocate(tkdry        )

  end subroutine set_material_properties

  !------------------------------------------------------------------------
  subroutine set_initial_conditions()
    !
    ! !DESCRIPTION:
    !
#include <petsc/finclude/petsc.h>
    !
    use MultiPhysicsProbThermalEnthalpy , only : thermal_enthalpy_mpp
    use MultiPhysicsProbConstants       , only : AUXVAR_INTERNAL, VAR_PRESSURE
    use petscsys
    use petscvec
    use petscdm
    !
    implicit none
    !
    PetscErrorCode      :: ierr
    PetscInt            :: nDM
    DM        , pointer :: dms(:)
    Vec       , pointer :: soln_subvecs(:)
    PetscReal , pointer :: temp_p(:)
    PetscInt            :: ii
    PetscViewer         :: viewer
    PetscReal , pointer :: pressure(:)
    PetscInt            :: soe_auxvar_id

    ! Find number of GEs packed within the SoE
    call DMCompositeGetNumberDM(thermal_enthalpy_mpp%soe%solver%dm, nDM, ierr)

    ! Get DMs for each GE
    allocate (dms(nDM))
    call DMCompositeGetEntriesArray(thermal_enthalpy_mpp%soe%solver%dm, dms, ierr)

    ! Allocate vectors for individual GEs
    allocate(soln_subvecs(nDM))

    ! Get solution vectors for individual GEs
    call DMCompositeGetAccessArray(thermal_enthalpy_mpp%soe%solver%dm, thermal_enthalpy_mpp%soe%solver%soln, nDM, &
         PETSC_NULL_INTEGER, soln_subvecs, ierr)

    do ii = 1, nDM
       call VecGetArrayF90(soln_subvecs(ii), temp_p, ierr)
       temp_p(:) = 283.15d0
       call VecRestoreArrayF90(soln_subvecs(ii), temp_p, ierr)
    enddo

    ! Restore solution vectors for individual GEs
    call DMCompositeRestoreAccessArray(thermal_enthalpy_mpp%soe%solver%dm, thermal_enthalpy_mpp%soe%solver%soln, nDM, &
         PETSC_NULL_INTEGER, soln_subvecs, ierr)

    call VecCopy(thermal_enthalpy_mpp%soe%solver%soln, thermal_enthalpy_mpp%soe%solver%soln_prev, ierr); CHKERRQ(ierr)
    call VecCopy(thermal_enthalpy_mpp%soe%solver%soln, thermal_enthalpy_mpp%soe%solver%soln_prev_clm, ierr); CHKERRQ(ierr)

    allocate(pressure(ncells_local))
    pressure(:) = 091325.d0
    soe_auxvar_id = -1
    call thermal_enthalpy_mpp%soe%SetDataFromCLM(AUXVAR_INTERNAL,  &
         VAR_PRESSURE, soe_auxvar_id, pressure)
    deallocate(pressure)

  end subroutine set_initial_conditions

  !------------------------------------------------------------------------
  subroutine set_bondary_conditions()
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbThermalEnthalpy , only : thermal_enthalpy_mpp
    use MultiPhysicsProbConstants       , only : AUXVAR_BC, VAR_BC_SS_CONDITION
    use MultiPhysicsProbConstants       , only : AUXVAR_INTERNAL, VAR_PRESSURE
    !
    implicit none
    !
    PetscReal, pointer :: top_bc(:)
    PetscReal, pointer :: bot_bc(:)
    PetscReal, pointer :: pressure(:)
    PetscInt           :: soe_auxvar_id

    allocate(top_bc(1))
    allocate(bot_bc(1))
    allocate(pressure(ncells_local))

    top_bc(:) = 303.15d0
    bot_bc(:) = 293.15d0
    pressure(:) = 091325.d0

    soe_auxvar_id = 1
    call thermal_enthalpy_mpp%soe%SetDataFromCLM(AUXVAR_BC,  &
         VAR_BC_SS_CONDITION, soe_auxvar_id, top_bc)

    soe_auxvar_id = 2
    call thermal_enthalpy_mpp%soe%SetDataFromCLM(AUXVAR_BC,  &
         VAR_BC_SS_CONDITION, soe_auxvar_id, bot_bc)

    soe_auxvar_id = -1
    call thermal_enthalpy_mpp%soe%SetDataFromCLM(AUXVAR_INTERNAL,  &
         VAR_PRESSURE, soe_auxvar_id, pressure)

    deallocate(top_bc)
    deallocate(bot_bc)
    deallocate(pressure)

  end subroutine set_bondary_conditions

  !------------------------------------------------------------------------
  subroutine output_regression_heat_transport_1D_problem(filename_base, num_cells)
    !
    use MultiPhysicsProbThermalEnthalpy , only : thermal_enthalpy_mpp
    use MultiPhysicsProbConstants       , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants       , only : VAR_PRESSURE
    use MultiPhysicsProbConstants       , only : VAR_TEMPERATURE
    use regression_mod                  , only : regression_type
    !
    implicit none
    !
    character(len=256)    :: filename_base
    character(len=512)    :: filename
    PetscInt, intent(in)  :: num_cells
    !
    PetscInt              :: output
    character(len=64)     :: category
    character(len=64)     :: name
    PetscReal, pointer    :: data(:)
    type(regression_type) :: regression

    allocate(data(ncells_local))

    call regression%Init(filename_base, num_cells)
    call regression%OpenOutput()

    name = 'temperature'
    category = 'general'
    call thermal_enthalpy_mpp%soe%GetDataForCLM(AUXVAR_INTERNAL,  &
         VAR_TEMPERATURE, -1, data)
    call regression%WriteData(name, category, data)

    call regression%CloseOutput()
    
    deallocate(data)

  end subroutine output_regression_heat_transport_1D_problem

end module heat_transport_1D_problem
