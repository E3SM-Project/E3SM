module vsfm_spac_campbell_problem

  implicit none
  
#include <petsc/finclude/petsc.h>
  PetscInt  , parameter :: nx       = 1
  PetscInt  , parameter :: ny       = 1
  PetscReal , parameter :: x_column = 1.d0
  PetscReal , parameter :: y_column = 1.d0
  PetscReal , parameter :: z_column = 1.d0
  PetscInt              :: nz_xylem
  PetscInt              :: nz_root
  PetscInt              :: nz_soil
  PetscInt              :: ncells_local
  PetscInt              :: ncells_ghost
  PetscReal             :: Campbell_n
  PetscReal             :: Campbell_b
  PetscReal             :: Campbell_he
  PetscReal             :: theta_s
  PetscReal             :: VG_alpha
  PetscReal             :: VG_n
  PetscReal , parameter :: PI  = 4 * atan (1.0_8)
  PetscBool             :: multi_goveqns_formulation

  public :: run_vsfm_spac_campbell_problem
  public :: output_regression_vsfm_spac_campbell_problem

contains

  subroutine run_vsfm_spac_campbell_problem()
    !
#include <petsc/finclude/petsc.h>
    !
    use MultiPhysicsProbVSFM , only : vsfm_mpp
    use mpp_varpar           , only : mpp_varpar_init
    use petscsys
    use petscvec
    use petscmat
    use petscts
    use petscsnes
    use petscdm
    use petscdmda
    !    
    implicit none
    
    !
    PetscBool          :: converged
    PetscInt           :: converged_reason
    PetscErrorCode     :: ierr
    PetscReal          :: dtime
    PetscReal          :: time
    PetscInt           :: istep, nstep
    PetscBool          :: flg
    PetscBool          :: save_initial_soln, save_final_soln
    character(len=256) :: string
    character(len=256) :: output_suffix
    PetscViewer        :: viewer

    ! Set default settings
    nz_xylem               = 2
    nz_root                = 28
    nz_soil                = 50
    dtime                  = 3600.d0
    nstep                  = 24
    save_initial_soln      = PETSC_FALSE
    save_final_soln        = PETSC_FALSE
    output_suffix          = ''
    multi_goveqns_formulation = PETSC_FALSE

    ! Get some command line options

    call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-dt',dtime,flg,ierr)
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-nstep',nstep,flg,ierr)
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-save_initial_soln',save_initial_soln,flg,ierr)
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-save_final_soln',save_final_soln,flg,ierr)
    call PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-output_suffix',output_suffix,flg,ierr)
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-multi_goveqns_formulation',multi_goveqns_formulation,flg,ierr)

    ! Initialize the problem
    call Init()

    time = 0.d0

    do istep = 1, nstep

       call set_bondary_conditions(time)
       time = time + dtime

       ! Run the model
       call vsfm_mpp%soe%StepDT(dtime, istep, &
            converged, converged_reason, ierr); CHKERRQ(ierr)

    end do

  end subroutine run_vsfm_spac_campbell_problem

  !------------------------------------------------------------------------
  subroutine Init()
    !
    use MultiPhysicsProbVSFM , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : MPP_VSFM_SNES_CLM
    use SystemOfEquationsVSFMType , only : VSFMSOEUpdateConnections
    !
    implicit none
    !

    ! 1. Initialize the multi-physics-problem (MPP)
    call initialize_mpp()

    if (.not.multi_goveqns_formulation) then
       ! 2. Add all meshes needed for the MPP
       call add_single_mesh()

       ! 3. Add all governing equations
       call add_single_goveqn()

    else
       ! 2. Add all meshes needed for the MPP
       call add_multiple_meshes()

       ! 3. Add all governing equations
       call add_multiple_goveqns()

    end if

    ! 4. Add boundary and source-sink conditions to all governing equations
    call add_conditions_to_goveqns()

    ! 5. Allocate memory to hold auxvars
    call allocate_auxvars()

    ! 6. Setup the MPP
    call vsfm_mpp%SetupProblem()

    ! 7. Add material properities associated with all governing equations
    call set_material_properties() 

    ! 9. Set flux type
    call set_conn_flux_type()
    
    ! 8. Set initial conditions
    call set_initial_conditions()

    call VSFMSOEUpdateConnections(vsfm_mpp%soe, MPP_VSFM_SNES_CLM)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine initialize_mpp()
    !
    ! !DESCRIPTION:
    ! Initialization VSFM
    !
    !
#include <petsc/finclude/petsc.h>
    ! !USES:
    use MultiPhysicsProbConstants , only : MPP_VSFM_SNES_CLM
    use MultiPhysicsProbVSFM , only : vsfm_mpp
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
    call vsfm_mpp%Init       ()
    call vsfm_mpp%SetName    ('Variably-Saturated-Flow-Model')
    call vsfm_mpp%SetID      (MPP_VSFM_SNES_CLM)
    call vsfm_mpp%SetMPIRank (iam)

  end subroutine initialize_mpp

  !------------------------------------------------------------------------
  subroutine add_multiple_meshes()
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : MESH_ALONG_GRAVITY
    use MultiPhysicsProbConstants , only : MESH_AGAINST_GRAVITY
    use MultiPhysicsProbConstants , only : MESH_CLM_SOIL_COL
    use MultiPhysicsProbConstants , only : MESH_SPAC_ROOT_COL
    use MultiPhysicsProbConstants , only : MESH_SPAC_XYLEM_COL
    use MultiPhysicsProbConstants , only : VAR_XC
    use MultiPhysicsProbConstants , only : VAR_YC
    use MultiPhysicsProbConstants , only : VAR_ZC
    use MultiPhysicsProbConstants , only : VAR_DX
    use MultiPhysicsProbConstants , only : VAR_DY
    use MultiPhysicsProbConstants , only : VAR_DZ
    use MultiPhysicsProbConstants , only : VAR_AREA
    use MultiPhysicsProbConstants , only : VAR_VOLUME
    use MultiPhysicsProbConstants , only : CONN_SET_INTERNAL
    use MultiPhysicsProbConstants , only : CONN_SET_LATERAL
    use MultiPhysicsProbConstants , only : CONN_VERTICAL
    use mpp_varpar                , only : mpp_varpar_set_nlevsoi, mpp_varpar_set_nlevgrnd
    !
    implicit none
    !
#include <petsc/finclude/petsc.h>
    !
    PetscReal :: dx, dy, dz
    PetscInt :: imesh, kk
    PetscInt :: nlev
    PetscInt :: iconn, vert_nconn

    PetscInt :: ncells_xylem
    PetscInt :: ncells_root
    PetscInt :: ncells_soil

    PetscInt :: nconn_xylem
    PetscInt :: nconn_root
    PetscInt :: nconn_soil

    PetscReal, pointer :: xc_xylem(:)           ! x-position of grid cell [m]
    PetscReal, pointer :: yc_xylem(:)           ! y-position of grid cell [m]
    PetscReal, pointer :: zc_xylem(:)           ! z-position of grid cell [m]
    PetscReal, pointer :: dx_xylem(:)           ! layer thickness of grid cell [m]
    PetscReal, pointer :: dy_xylem(:)           ! layer thickness of grid cell [m]
    PetscReal, pointer :: dz_xylem(:)           ! layer thickness of grid cell [m]
    PetscReal, pointer :: area_xylem(:)         ! area of grid cell [m^2]
    PetscReal, pointer :: vol_xylem(:)          ! volume of grid cell [m^3]
    PetscInt , pointer :: filter_xylem(:)       ! 

    PetscReal, pointer :: xc_root(:)           ! x-position of grid cell [m]
    PetscReal, pointer :: yc_root(:)           ! y-position of grid cell [m]
    PetscReal, pointer :: zc_root(:)           ! z-position of grid cell [m]
    PetscReal, pointer :: dx_root(:)           ! layer thickness of grid cell [m]
    PetscReal, pointer :: dy_root(:)           ! layer thickness of grid cell [m]
    PetscReal, pointer :: dz_root(:)           ! layer thickness of grid cell [m]
    PetscReal, pointer :: area_root(:)         ! area of grid cell [m^2]
    PetscReal, pointer :: vol_root(:)          ! volume of grid cell [m^3]
    PetscInt , pointer :: filter_root(:)       ! 

    PetscReal, pointer :: xc_soil(:)           ! x-position of grid cell [m]
    PetscReal, pointer :: yc_soil(:)           ! y-position of grid cell [m]
    PetscReal, pointer :: zc_soil(:)           ! z-position of grid cell [m]
    PetscReal, pointer :: dx_soil(:)           ! layer thickness of grid cell [m]
    PetscReal, pointer :: dy_soil(:)           ! layer thickness of grid cell [m]
    PetscReal, pointer :: dz_soil(:)           ! layer thickness of grid cell [m]
    PetscReal, pointer :: area_soil(:)         ! area of grid cell [m^2]
    PetscReal, pointer :: vol_soil(:)          ! volume of grid cell [m^3]
    PetscInt , pointer :: filter_soil(:)       ! 
    
    PetscInt, pointer  :: vert_conn_id_up_xylem(:)   !
    PetscInt, pointer  :: vert_conn_id_dn_xylem(:)   !
    PetscReal, pointer :: vert_conn_dist_up_xylem(:) !
    PetscReal, pointer :: vert_conn_dist_dn_xylem(:) !
    PetscReal, pointer :: vert_conn_area_xylem(:)    !
    PetscInt , pointer :: vert_conn_type_xylem(:)    !

    PetscInt, pointer  :: vert_conn_id_up_soil(:)   !
    PetscInt, pointer  :: vert_conn_id_dn_soil(:)   !
    PetscReal, pointer :: vert_conn_dist_up_soil(:) !
    PetscReal, pointer :: vert_conn_dist_dn_soil(:) !
    PetscReal, pointer :: vert_conn_area_soil(:)    !
    PetscInt , pointer :: vert_conn_type_soil(:)    !

    PetscErrorCode :: ierr

    dx = x_column/nx
    dy = y_column/ny
    dz = z_column/nz_soil

    ncells_ghost = 0
    ncells_xylem = nx * ny * nz_xylem
    ncells_root  = nx * ny * nz_root
    ncells_soil  = nx * ny * nz_soil

    allocate(xc_xylem     (ncells_xylem ))
    allocate(yc_xylem     (ncells_xylem ))
    allocate(zc_xylem     (ncells_xylem ))
    allocate(dx_xylem     (ncells_xylem ))
    allocate(dy_xylem     (ncells_xylem ))
    allocate(dz_xylem     (ncells_xylem ))
    allocate(area_xylem   (ncells_xylem ))
    allocate(filter_xylem (ncells_xylem ))
    allocate(vol_xylem    (ncells_xylem ))

    allocate(xc_root      (ncells_root  ))
    allocate(yc_root      (ncells_root  ))
    allocate(zc_root      (ncells_root  ))
    allocate(dx_root      (ncells_root  ))
    allocate(dy_root      (ncells_root  ))
    allocate(dz_root      (ncells_root  ))
    allocate(area_root    (ncells_root  ))
    allocate(filter_root  (ncells_root  ))
    allocate(vol_root     (ncells_root  ))

    allocate(xc_soil      (ncells_soil  ))
    allocate(yc_soil      (ncells_soil  ))
    allocate(zc_soil      (ncells_soil  ))
    allocate(dx_soil      (ncells_soil  ))
    allocate(dy_soil      (ncells_soil  ))
    allocate(dz_soil      (ncells_soil  ))
    allocate(area_soil    (ncells_soil  ))
    allocate(filter_soil  (ncells_soil  ))
    allocate(vol_soil     (ncells_soil  ))
    
    filter_xylem (:) = 1
    area_xylem   (:) = 1.d0
    dx_xylem     (:) = 1.d0
    dy_xylem     (:) = 1.d0
    dz_xylem     (:) = 1.d0/50.d0
    xc_xylem     (:) = 1.d0/2.d0
    yc_xylem     (:) = 1.d0/2.d0
    vol_xylem    (:) = 1.d0/50.d0

    filter_root (:) = 1
    area_root   (:) = 1.d0
    dx_root     (:) = 1.d0
    dy_root     (:) = 1.d0
    dz_root     (:) = 1.d0/50.d0
    xc_root     (:) = 1.d0/2.d0
    yc_root     (:) = 1.d0/2.d0
    vol_root    (:) = 1.d0/50.d0

    filter_soil (:) = 1
    area_soil   (:) = 1.d0
    dx_soil     (:) = 1.d0
    dy_soil     (:) = 1.d0
    dz_soil     (:) = 1.d0/50.d0
    xc_soil     (:) = 1.d0/2.d0
    yc_soil     (:) = 1.d0/2.d0
    vol_soil    (:) = 1.d0/50.d0

    zc_xylem(1) = 0.d0
    zc_xylem(2) = 0.d0

    vol_soil(1) = vol_soil(1) / 2.d0
    
    do kk = 1,nz_root
       zc_root(kk) = -(dz/2.d0 + dz * (kk - 1 + 2))
    enddo

    do kk = 1, nz_soil
       zc_soil(kk) = -(dz/2.d0 + dz * (kk - 1))
    end do

    call mpp_varpar_set_nlevsoi(nz_soil)
    call mpp_varpar_set_nlevgrnd(nz_soil)

    !
    ! Set up the meshes
    !    
    call vsfm_mpp%SetNumMeshes(3)

    ! Xylem Mesh
    imesh        = 1
    call vsfm_mpp%MeshSetName                (imesh, 'Xylem mesh')
    call vsfm_mpp%MeshSetOrientation         (imesh, MESH_AGAINST_GRAVITY)
    call vsfm_mpp%MeshSetID                  (imesh, MESH_SPAC_XYLEM_COL)
    call vsfm_mpp%MeshSetDimensions          (imesh, ncells_xylem, ncells_ghost, nz_xylem)

    call vsfm_mpp%MeshSetGridCellFilter      (imesh, filter_xylem)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_XC     , xc_xylem)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_YC     , yc_xylem)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_ZC     , zc_xylem)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_DX     , dx_xylem)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_DY     , dy_xylem)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_DZ     , dz_xylem)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_AREA   , area_xylem)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_VOLUME , vol_xylem)

    ! Root Mesh
    
    imesh        = 2
    call vsfm_mpp%MeshSetName                (imesh, 'Root mesh')
    call vsfm_mpp%MeshSetOrientation         (imesh, MESH_AGAINST_GRAVITY)
    call vsfm_mpp%MeshSetID                  (imesh, MESH_SPAC_ROOT_COL)
    call vsfm_mpp%MeshSetDimensions          (imesh, ncells_root, ncells_ghost, nz_root)

    call vsfm_mpp%MeshSetGridCellFilter      (imesh, filter_root)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_XC     , xc_root)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_YC     , yc_root)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_ZC     , zc_root)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_DX     , dx_root)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_DY     , dy_root)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_DZ     , dz_root)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_AREA   , area_root)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_VOLUME , vol_root)
    
    ! Soil Mesh
    imesh        = 3
    call vsfm_mpp%MeshSetName                (imesh, 'Soil mesh')
    call vsfm_mpp%MeshSetOrientation         (imesh, MESH_AGAINST_GRAVITY)
    call vsfm_mpp%MeshSetID                  (imesh, MESH_CLM_SOIL_COL)
    call vsfm_mpp%MeshSetDimensions          (imesh, ncells_soil, ncells_ghost, nz_soil)

    call vsfm_mpp%MeshSetGridCellFilter      (imesh, filter_soil)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_XC     , xc_soil)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_YC     , yc_soil)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_ZC     , zc_soil)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_DX     , dx_soil)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_DY     , dy_soil)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_DZ     , dz_soil)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_AREA   , area_soil)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_VOLUME , vol_soil)


    nconn_xylem = 1
    nconn_soil  = nz_soil - 1

    allocate (vert_conn_id_up_xylem   (nconn_xylem))
    allocate (vert_conn_id_dn_xylem   (nconn_xylem))
    allocate (vert_conn_dist_up_xylem (nconn_xylem))
    allocate (vert_conn_dist_dn_xylem (nconn_xylem))
    allocate (vert_conn_area_xylem    (nconn_xylem))
    allocate (vert_conn_type_xylem    (nconn_xylem))

    allocate (vert_conn_id_up_soil   (nconn_soil))
    allocate (vert_conn_id_dn_soil   (nconn_soil))
    allocate (vert_conn_dist_up_soil (nconn_soil))
    allocate (vert_conn_dist_dn_soil (nconn_soil))
    allocate (vert_conn_area_soil    (nconn_soil))
    allocate (vert_conn_type_soil    (nconn_soil))

    iconn = 0
    iconn = iconn + 1
    vert_conn_id_up_xylem(iconn)   = 1
    vert_conn_id_dn_xylem(iconn)   = 2
    vert_conn_dist_up_xylem(iconn) = 0.5d0*dz
    vert_conn_dist_dn_xylem(iconn) = 0.5d0*dz
    vert_conn_area_xylem(iconn)    = area_soil(1)
    vert_conn_type_xylem(iconn)    = CONN_VERTICAL

    iconn = 0
    do kk = 1, nz_soil - 1
       iconn = iconn + 1
       vert_conn_id_up_soil(iconn)   = kk 
       vert_conn_id_dn_soil(iconn)   = kk + 1
       vert_conn_dist_up_soil(iconn) = 0.5d0*dz
       vert_conn_dist_dn_soil(iconn) = 0.5d0*dz
       vert_conn_area_soil(iconn)    = area_soil(kk)
       vert_conn_type_soil(iconn)    = CONN_VERTICAL
    end do

    imesh = 1
    call vsfm_mpp%MeshSetConnectionSet(imesh, CONN_SET_INTERNAL, &
         nconn_xylem,  vert_conn_id_up_xylem, vert_conn_id_dn_xylem,          &
         vert_conn_dist_up_xylem, vert_conn_dist_dn_xylem,  vert_conn_area_xylem,  &
         vert_conn_type_xylem)

    imesh = 3
    vert_nconn = nconn_soil
    call vsfm_mpp%MeshSetConnectionSet(imesh, CONN_SET_INTERNAL, &
         nconn_soil,  vert_conn_id_up_soil, vert_conn_id_dn_soil,          &
         vert_conn_dist_up_soil, vert_conn_dist_dn_soil,  vert_conn_area_soil,  &
         vert_conn_type_soil)

    deallocate (xc_xylem                )
    deallocate (yc_xylem                )
    deallocate (zc_xylem                )
    deallocate (dx_xylem                )
    deallocate (dy_xylem                )
    deallocate (dz_xylem                )
    deallocate (area_xylem              )
    deallocate (filter_xylem            )
    deallocate (vol_xylem               )

    deallocate (xc_root                 )
    deallocate (yc_root                 )
    deallocate (zc_root                 )
    deallocate (dx_root                 )
    deallocate (dy_root                 )
    deallocate (dz_root                 )
    deallocate (area_root               )
    deallocate (filter_root             )
    deallocate (vol_root                )

    deallocate (xc_soil                 )
    deallocate (yc_soil                 )
    deallocate (zc_soil                 )
    deallocate (dx_soil                 )
    deallocate (dy_soil                 )
    deallocate (dz_soil                 )
    deallocate (area_soil               )
    deallocate (filter_soil             )
    deallocate (vol_soil                )

    deallocate (vert_conn_id_up_xylem   )
    deallocate (vert_conn_id_dn_xylem   )
    deallocate (vert_conn_dist_up_xylem )
    deallocate (vert_conn_dist_dn_xylem )
    deallocate (vert_conn_area_xylem    )
    deallocate (vert_conn_type_xylem    )

    deallocate (vert_conn_id_up_soil    )
    deallocate (vert_conn_id_dn_soil    )
    deallocate (vert_conn_dist_up_soil  )
    deallocate (vert_conn_dist_dn_soil  )
    deallocate (vert_conn_area_soil     )
    deallocate (vert_conn_type_soil     )

  end subroutine add_multiple_meshes

  !------------------------------------------------------------------------
  subroutine add_single_mesh()
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : MESH_ALONG_GRAVITY
    use MultiPhysicsProbConstants , only : MESH_AGAINST_GRAVITY
    use MultiPhysicsProbConstants , only : MESH_CLM_SOIL_COL
    use MultiPhysicsProbConstants , only : VAR_XC
    use MultiPhysicsProbConstants , only : VAR_YC
    use MultiPhysicsProbConstants , only : VAR_ZC
    use MultiPhysicsProbConstants , only : VAR_DX
    use MultiPhysicsProbConstants , only : VAR_DY
    use MultiPhysicsProbConstants , only : VAR_DZ
    use MultiPhysicsProbConstants , only : VAR_AREA
    use MultiPhysicsProbConstants , only : VAR_VOLUME
    use MultiPhysicsProbConstants , only : CONN_SET_INTERNAL
    use MultiPhysicsProbConstants , only : CONN_SET_LATERAL
    use MultiPhysicsProbConstants , only : CONN_VERTICAL
    use mpp_varpar                , only : mpp_varpar_set_nlevsoi, mpp_varpar_set_nlevgrnd
    !
    implicit none
    !
#include <petsc/finclude/petsc.h>
    !
    PetscReal :: dx, dy, dz
    PetscInt :: imesh, kk
    PetscInt :: nlev, nz, nconn
    PetscInt :: iconn, vert_nconn
    PetscReal, pointer :: soil_xc(:)           ! x-position of grid cell [m]
    PetscReal, pointer :: soil_yc(:)           ! y-position of grid cell [m]
    PetscReal, pointer :: soil_zc(:)           ! z-position of grid cell [m]
    PetscReal, pointer :: soil_dx(:)           ! layer thickness of grid cell [m]
    PetscReal, pointer :: soil_dy(:)           ! layer thickness of grid cell [m]
    PetscReal, pointer :: soil_dz(:)           ! layer thickness of grid cell [m]
    PetscReal, pointer :: soil_area(:)         ! area of grid cell [m^2]
    PetscReal, pointer :: soil_vol(:)          ! volume of grid cell [m^3]
    PetscInt , pointer :: soil_filter(:)       ! 

    PetscInt, pointer  :: vert_conn_id_up(:)   !
    PetscInt, pointer  :: vert_conn_id_dn(:)   !
    PetscReal, pointer :: vert_conn_dist_up(:) !
    PetscReal, pointer :: vert_conn_dist_dn(:) !
    PetscReal, pointer :: vert_conn_area(:)    !
    PetscInt , pointer :: vert_conn_type(:)    !

    PetscErrorCode :: ierr

    dx = x_column/nx
    dy = y_column/ny
    dz = z_column/nz_soil

    imesh        = 1
    nlev         = nz_xylem + nz_root + nz_soil
    nz           = nlev
    ncells_local = nx*ny*nlev
    ncells_ghost = 0

    allocate(soil_xc     (nz))
    allocate(soil_yc     (nz))
    allocate(soil_zc     (nz))
    allocate(soil_dx     (nz))
    allocate(soil_dy     (nz))
    allocate(soil_dz     (nz))
    allocate(soil_area   (nz))
    allocate(soil_filter (nz))
    allocate(soil_vol    (nz))

    soil_filter (:) = 1
    soil_area   (:) = 1.d0
    soil_dx     (:) = 1.d0
    soil_dy     (:) = 1.d0
    soil_dz     (:) = 1.d0/50.d0
    soil_xc     (:) = 1.d0/2.d0
    soil_yc     (:) = 1.d0/2.d0
    soil_vol    (:) = 1.d0/50.d0

    soil_zc(1) = 0.d0
    soil_zc(2) = 0.d0
    soil_vol(31) = soil_vol(1) / 2.d0
    
    do kk = 3,nz_xylem + nz_root
       soil_zc(kk) = -(dz/2.d0 + dz * (kk - 1))
    enddo

    do kk = nz_xylem + nz_root + 1, nz_xylem + nz_root + nz_soil
       soil_zc(kk) = -(dz/2.d0 + dz * (kk - nz_xylem - nz_root - 1))
    end do

    call mpp_varpar_set_nlevsoi(nz)
    call mpp_varpar_set_nlevgrnd(nz)

    !
    ! Set up the meshes
    !    
    call vsfm_mpp%SetNumMeshes(1)

    call vsfm_mpp%MeshSetName        (imesh, 'Soil mesh')
    call vsfm_mpp%MeshSetOrientation (imesh, MESH_AGAINST_GRAVITY)
    call vsfm_mpp%MeshSetID          (imesh, MESH_CLM_SOIL_COL)
    call vsfm_mpp%MeshSetDimensions  (imesh, ncells_local, ncells_ghost, nlev)

    call vsfm_mpp%MeshSetGridCellFilter      (imesh, soil_filter)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_XC   , soil_xc)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_YC   , soil_yc)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_ZC   , soil_zc)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_DX   , soil_dx)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_DY   , soil_dy)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_DZ   , soil_dz)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_AREA , soil_area)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_VOLUME , soil_vol)
    !call vsfm_mpp%MeshComputeVolume          (imesh)

    nconn = 1 + nz_root * 2 + nz_soil - 1

    allocate (vert_conn_id_up   (nconn))
    allocate (vert_conn_id_dn   (nconn))
    allocate (vert_conn_dist_up (nconn))
    allocate (vert_conn_dist_dn (nconn))
    allocate (vert_conn_area    (nconn))
    allocate (vert_conn_type    (nconn))

    iconn = 0
    iconn = iconn + 1
    vert_conn_id_up(iconn)   = 1
    vert_conn_id_dn(iconn)   = 2
    vert_conn_dist_up(iconn) = 0.5d0*dz
    vert_conn_dist_dn(iconn) = 0.5d0*dz
    vert_conn_area(iconn)    = soil_area(1)
    vert_conn_type(iconn)    = CONN_VERTICAL

    do kk = 2, nz_xylem + nz_root -1

       iconn = iconn + 1
       vert_conn_id_up(iconn)   = 2
       vert_conn_id_dn(iconn)   = kk + 1
       vert_conn_dist_up(iconn) = 0.5d0*dz
       vert_conn_dist_dn(iconn) = 0.5d0*dz
       vert_conn_area(iconn)    = soil_area(kk)
       vert_conn_type(iconn)    = CONN_VERTICAL

    end do

    do kk = 2, nz_xylem + nz_root -1

       iconn = iconn + 1
       vert_conn_id_up(iconn)   = kk + 1
       vert_conn_id_dn(iconn)   = kk + 1 + nz_xylem + nz_root
       vert_conn_dist_up(iconn) = 0.5d0*dz
       vert_conn_dist_dn(iconn) = 0.5d0*dz
       vert_conn_area(iconn)    = soil_area(kk)
       vert_conn_type(iconn)    = CONN_VERTICAL

    end do

    do kk = 1, nz_soil - 1
       iconn = iconn + 1
       vert_conn_id_up(iconn)   = kk + nz_xylem + nz_root
       vert_conn_id_dn(iconn)   = kk + nz_xylem + nz_root + 1
       vert_conn_dist_up(iconn) = 0.5d0*dz
       vert_conn_dist_dn(iconn) = 0.5d0*dz
       vert_conn_area(iconn)    = soil_area(kk)
       vert_conn_type(iconn)    = CONN_VERTICAL
    end do
    
    vert_nconn = iconn

    call vsfm_mpp%MeshSetConnectionSet(imesh, CONN_SET_INTERNAL, &
         vert_nconn,  vert_conn_id_up, vert_conn_id_dn,          &
         vert_conn_dist_up, vert_conn_dist_dn,  vert_conn_area,  &
         vert_conn_type)

    deallocate(soil_xc)
    deallocate(soil_yc)
    deallocate(soil_zc)
    deallocate(soil_dx)
    deallocate(soil_dy)
    deallocate(soil_dz)
    deallocate(soil_area)
    deallocate(soil_filter)  

  end subroutine add_single_mesh

  !------------------------------------------------------------------------

  subroutine add_single_goveqn()
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbVSFM , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : GE_RE
    use MultiPhysicsProbConstants , only : MESH_CLM_SOIL_COL
    !
    ! !ARGUMENTS
    implicit none

    call vsfm_mpp%AddGovEqn(GE_RE, 'Richards Equation ODE', MESH_CLM_SOIL_COL)

    call vsfm_mpp%SetMeshesOfGoveqns()

  end subroutine add_single_goveqn

  !------------------------------------------------------------------------

  subroutine add_multiple_goveqns()
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbVSFM     , only : vsfm_mpp
    use MultiPhysicsProbConstants, only : GE_RE
    use MultiPhysicsProbConstants, only : MESH_CLM_SOIL_COL
    use MultiPhysicsProbConstants, only : MESH_SPAC_ROOT_COL
    use MultiPhysicsProbConstants, only : MESH_SPAC_XYLEM_COL
    !
    ! !ARGUMENTS
    implicit none

    call vsfm_mpp%AddGovEqn(GE_RE, 'Richards Equation ODE for Xylem', MESH_SPAC_XYLEM_COL)
    call vsfm_mpp%AddGovEqn(GE_RE, 'Richards Equation ODE for Root' , MESH_SPAC_ROOT_COL )
    call vsfm_mpp%AddGovEqn(GE_RE, 'Richards Equation ODE for Soil' , MESH_CLM_SOIL_COL  )

    call vsfm_mpp%SetMeshesOfGoveqns()

  end subroutine add_multiple_goveqns

  !------------------------------------------------------------------------
  subroutine add_conditions_to_goveqns()
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : ALL_CELLS
    use MultiPhysicsProbConstants , only : SOIL_TOP_CELLS
    use MultiPhysicsProbConstants , only : SOIL_BOTTOM_CELLS
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_SS
    use MultiPhysicsProbConstants , only : COND_DIRICHLET
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    use MultiPhysicsProbConstants , only : COND_DOWNREG_MASS_RATE_CAMPBELL
    use MultiPhysicsProbConstants , only : CONN_VERTICAL
    use ConnectionSetType         , only : connection_set_type
    use ConnectionSetType         , only : ConnectionSetDestroy
    use MeshType                  , only : MeshCreateConnectionSet
    use petscsys
    !
    ! !ARGUMENTS
    implicit none
    !
    PetscInt                            :: kk
    PetscInt                            :: ieqn
    PetscInt                            :: nconn

    PetscInt                  , pointer :: id_up(:)
    PetscInt                  , pointer :: id_dn(:)
    PetscReal                 , pointer :: dist_up(:)
    PetscReal                 , pointer :: dist_dn(:)
    PetscReal                 , pointer :: area(:)
    PetscInt                  , pointer :: itype(:)
    PetscReal                 , pointer :: unit_vec(:,:)
    PetscInt                            :: num_other_goveqs
    PetscInt                  , pointer :: ieqn_others(:)

    type(connection_set_type) , pointer :: conn_set

    ieqn       = 1
    call vsfm_mpp%soe%AddConditionInGovEqn(ieqn, COND_SS,   &
         'Potential Mass_Flux', 'kg/s', COND_DOWNREG_MASS_RATE_CAMPBELL, &
         SOIL_BOTTOM_CELLS)

    if (.not.multi_goveqns_formulation) return

    nconn = nz_root
    allocate(id_up    (nconn   ))
    allocate(id_dn    (nconn   ))
    allocate(dist_up  (nconn   ))
    allocate(dist_dn  (nconn   ))
    allocate(area     (nconn   ))
    allocate(itype    (nconn   ))
    allocate(unit_vec (nconn,3 ))

    do kk = 1, nz_root
       id_up(kk)      = 0
       dist_up(kk)    = z_column/nz_soil/2.d0
       dist_dn(kk)    = z_column/nz_soil/2.d0
       area(kk)       = 1.d0
       unit_vec(kk,1) = 1.d0
       unit_vec(kk,2) = 0.d0
       unit_vec(kk,3) = 0.d0
       itype(kk)      = CONN_VERTICAL
    enddo
    
    allocate(conn_set)
    
    num_other_goveqs = 1 
    allocate(ieqn_others(num_other_goveqs))

    ! BoundaryCondition:
    ! Xylem mass ---> Root mass
    ieqn           = 1
    ieqn_others(1) = 2

    do kk = 1, nz_root
       id_up(kk)      = 0
       id_dn(kk)      = 2
    end do

    call MeshCreateConnectionSet(vsfm_mpp%meshes(2), &
         nconn, id_up, id_dn, &
         dist_up, dist_dn, area, itype, unit_vec, conn_set)

    call vsfm_mpp%soe%AddCouplingBCsInGovEqn( ieqn   , &
         name               = 'Root BC in xylem equation'  , &
         unit               = 'Pa'                         , &
         region_type        = SOIL_TOP_CELLS               , &
         num_other_goveqs   = num_other_goveqs             , &
         id_of_other_goveqs = ieqn_others                  , &
         conn_set           = conn_set)

    ! BoundaryCondition:
    ! Root mass ---> Xylem mass
    !           ---> Soil mass
    ieqn           = 2
    ieqn_others(1) = 1
    
    do kk = 1, nz_root
       id_up(kk)      = 0
       id_dn(kk)      = kk
    end do

    call MeshCreateConnectionSet(vsfm_mpp%meshes(2), &
         nconn, id_up, id_dn, &
         dist_up, dist_dn, area, itype, unit_vec, conn_set)

    call vsfm_mpp%soe%AddCouplingBCsInGovEqn( ieqn   , &
         name               = 'Xylem BC in root equation'  , &
         unit               = 'Pa'                         , &
         region_type        = SOIL_TOP_CELLS               , &
         num_other_goveqs   = num_other_goveqs             , &
         id_of_other_goveqs = ieqn_others                  , &
         conn_set           = conn_set)

    ieqn_others(1) = 3
    
    call MeshCreateConnectionSet(vsfm_mpp%meshes(2), &
         nconn, id_up, id_dn, &
         dist_up, dist_dn, area, itype, unit_vec, conn_set)

    call vsfm_mpp%soe%AddCouplingBCsInGovEqn( ieqn   , &
         name               = 'Soil BC in root equation'   , &
         unit               = 'Pa'                         , &
         region_type        = SOIL_TOP_CELLS               , &
         num_other_goveqs   = num_other_goveqs             , &
         id_of_other_goveqs = ieqn_others                  , &
         conn_set           = conn_set)

    ! BoundaryCondition:
    ! Soil mass ---> Root mass

    ieqn           = 3
    ieqn_others(1) = 2
    
    do kk = 1, nz_root
       id_up(kk)      = 0
       id_dn(kk)      = kk + 2
    end do

    call MeshCreateConnectionSet(vsfm_mpp%meshes(3), &
         nconn, id_up, id_dn, &
         dist_up, dist_dn, area, itype, unit_vec, conn_set)

    call vsfm_mpp%soe%AddCouplingBCsInGovEqn( ieqn   , &
         name               = 'Root BC in soil equation'   , &
         unit               = 'Pa'                         , &
         region_type        = SOIL_TOP_CELLS               , &
         num_other_goveqs   = num_other_goveqs             , &
         id_of_other_goveqs = ieqn_others                  , &
         conn_set           = conn_set)


    call ConnectionSetDestroy(conn_set)

    deallocate(id_up    )
    deallocate(id_dn    )
    deallocate(dist_up  )
    deallocate(dist_dn  )
    deallocate(area     )
    deallocate(itype    )
    deallocate(unit_vec )
    deallocate(ieqn_others)

  end subroutine add_conditions_to_goveqns

  !------------------------------------------------------------------------
  subroutine allocate_auxvars()
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbVSFM , only : vsfm_mpp
    use MultiPhysicsProbConstants, only : VAR_PRESSURE
    !
    implicit none
    !
    integer           :: ieqn
    integer           :: nvars_for_coupling
    integer, pointer  :: var_ids_for_coupling(:)
    integer, pointer  :: goveqn_ids_for_coupling(:)
    integer, pointer  :: is_bc(:)

    !
    ! Allocate auxvars
    !

    if (multi_goveqns_formulation) then
       !
       ! XYLEM <---> ROOT
       !
       ieqn               = 1
       nvars_for_coupling = 1

       allocate (var_ids_for_coupling    (nvars_for_coupling))
       allocate (goveqn_ids_for_coupling (nvars_for_coupling))

       var_ids_for_coupling    (1) = VAR_PRESSURE
       goveqn_ids_for_coupling (1) = 2

       call vsfm_mpp%GovEqnSetCouplingVars(ieqn, nvars_for_coupling, &
            var_ids_for_coupling, goveqn_ids_for_coupling)
       
       deallocate(var_ids_for_coupling   )
       deallocate(goveqn_ids_for_coupling)

       !
       ! ROOT <---> XYLEM
       ! ROOT <---> SOIL
       !
       ieqn               = 2
       nvars_for_coupling = 2

       allocate (var_ids_for_coupling    (nvars_for_coupling))
       allocate (goveqn_ids_for_coupling (nvars_for_coupling))

       var_ids_for_coupling    (1) = VAR_PRESSURE
       var_ids_for_coupling    (2) = VAR_PRESSURE
       goveqn_ids_for_coupling (1) = 1
       goveqn_ids_for_coupling (2) = 3

       call vsfm_mpp%GovEqnSetCouplingVars(ieqn, nvars_for_coupling, &
            var_ids_for_coupling, goveqn_ids_for_coupling)

       deallocate(var_ids_for_coupling   )
       deallocate(goveqn_ids_for_coupling)

       !
       ! SOIL <---> ROOT
       !
       ieqn               = 3
       nvars_for_coupling = 1

       allocate (var_ids_for_coupling    (nvars_for_coupling))
       allocate (goveqn_ids_for_coupling (nvars_for_coupling))

       var_ids_for_coupling    (1) = VAR_PRESSURE
       goveqn_ids_for_coupling (1) = 2

       call vsfm_mpp%GovEqnSetCouplingVars(ieqn, nvars_for_coupling, &
            var_ids_for_coupling, goveqn_ids_for_coupling)

       deallocate(var_ids_for_coupling   )
       deallocate(goveqn_ids_for_coupling)

    end if

    call vsfm_mpp%AllocateAuxVars()
    

  end subroutine allocate_auxvars

  !------------------------------------------------------------------------
  subroutine set_material_properties()
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSourceSinkAuxVarRealValue
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSoilPorosity
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSaturationFunction
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSoilPermeability
    use MultiPhysicsProbConstants , only : GRAVITY_CONSTANT
    use MultiPhysicsProbConstants , only : VAR_POT_MASS_SINK_PRESSURE
    use MultiPhysicsProbConstants , only : VAR_POT_MASS_SINK_EXPONENT
    use EOSWaterMod               , only : DENSITY_TGDPB01
    use mpp_varcon                , only : denh2o
    use mpp_varcon                , only : grav
    use SaturationFunction        , only : SAT_FUNC_VAN_GENUCHTEN, SAT_FUNC_BROOKS_COREY
    !
    implicit none
    !
    PetscReal , pointer   :: vsfm_watsat(:,:)
    PetscReal , pointer   :: vsfm_hksat(:,:)
    PetscReal , pointer   :: vsfm_bsw(:,:)
    PetscReal , pointer   :: vsfm_sucsat(:,:)
    PetscReal , pointer   :: vsfm_eff_porosity(:,:)
    PetscReal , pointer   :: vsfm_residual_sat(:,:)
    PetscReal , pointer   :: por(:)
    PetscReal , pointer   :: alpha(:)
    PetscReal , pointer   :: lambda(:)
    PetscReal , pointer   :: sat_res(:)
    PetscReal , pointer   :: perm(:)
    PetscReal , parameter :: vish2o = 0.001002d0    ! [N s/m^2] @ 20 degC
    PetscReal             :: Ks
    PetscInt              :: begc , endc
    integer   , pointer   :: vsfm_filter(:)
    PetscInt  , pointer   :: satfunc_type(:)
    PetscReal , pointer   :: ss_auxvar_value(:)
    PetscInt              :: nz
    PetscInt              :: ieqn
    !-----------------------------------------------------------------------

    Ks          = 0.001d0 ! [kg s m^{-3}]
    theta_s     = 0.46d0  !
    Campbell_b  = 4.58d0  ! [-]
    Campbell_he = -4.2d0  ! [J kg^{-1}]

    VG_n        = 1.35d0  ! [-]
    VG_alpha    = 0.15d0  ! [kg J^{-1}]

    Campbell_n  = 2.d0 + 3.d0/Campbell_b

    begc = 1
    endc = 1

    if (.not.multi_goveqns_formulation) then
       nz = nz_xylem + nz_root + nz_soil

       allocate (por          (nz))
       allocate (alpha        (nz))
       allocate (lambda       (nz))
       allocate (sat_res      (nz))
       allocate (satfunc_type (nz))
       allocate (perm         (nz))

       por          (1                      : nz_xylem + nz_root           ) = 0.d0
       sat_res      (1                      : nz_xylem + nz_root           ) = 0.d0
       lambda       (1                      : nz_xylem + nz_root           ) = 1.d0/Campbell_b
       alpha        (1                      : nz_xylem + nz_root           ) = 1.d-3/(-Campbell_he)
       perm         (1                      : nz_xylem + nz_root           ) = Ks/1.d6*8.904156d-4
       satfunc_type (1                      : nz_xylem + nz_root           ) = SAT_FUNC_BROOKS_COREY

       por          (nz_xylem + nz_root + 1 : nz_xylem + nz_root + nz_soil ) = theta_s
       sat_res      (nz_xylem + nz_root + 1 : nz_xylem + nz_root + nz_soil ) = 0.01d0
       lambda       (nz_xylem + nz_root + 1 : nz_xylem + nz_root + nz_soil ) = 1.d0 - 1.d0/VG_n
       alpha        (nz_xylem + nz_root + 1 : nz_xylem + nz_root + nz_soil ) = VG_alpha * 1.d-3
       perm         (nz_xylem + nz_root + 1 : nz_xylem + nz_root + nz_soil ) = Ks/1.d6*8.904156d-4
       satfunc_type (nz_xylem + nz_root + 1 : nz_xylem + nz_root + nz_soil ) = SAT_FUNC_VAN_GENUCHTEN


       call VSFMMPPSetSoilPorosity(vsfm_mpp, 1, por)

       call VSFMMPPSetSaturationFunction(vsfm_mpp, 1, satfunc_type, &
            alpha, lambda, sat_res)

       call VSFMMPPSetSoilPermeability(vsfm_mpp, 1, perm, perm, perm)

       deallocate(por          )
       deallocate(sat_res      )
       deallocate(lambda       )
       deallocate(alpha        )
       deallocate(perm         )
       deallocate(satfunc_type )

       allocate(ss_auxvar_value(1))

       ss_auxvar_value(:) = 10.d0
       call VSFMMPPSetSourceSinkAuxVarRealValue(vsfm_mpp, 1, &
            VAR_POT_MASS_SINK_EXPONENT, ss_auxvar_value)

       ss_auxvar_value(:) = -1500000.d0
       call VSFMMPPSetSourceSinkAuxVarRealValue(vsfm_mpp, 1, &
            VAR_POT_MASS_SINK_PRESSURE, ss_auxvar_value)

       deallocate(ss_auxvar_value)

    else

       !
       ! Xylem equation
       !
       ieqn = 1
       nz   = nz_xylem

       allocate (por          (nz))
       allocate (alpha        (nz))
       allocate (lambda       (nz))
       allocate (sat_res      (nz))
       allocate (satfunc_type (nz))
       allocate (perm         (nz))

       por          (:) = 0.d0
       sat_res      (:) = 0.d0
       lambda       (:) = 1.d0/Campbell_b
       alpha        (:) = 1.d-3/(-Campbell_he)
       perm         (:) = Ks/1.d6*8.904156d-4
       satfunc_type (:) = SAT_FUNC_BROOKS_COREY
       
       call VSFMMPPSetSoilPorosity(vsfm_mpp, ieqn, por)
       call VSFMMPPSetSaturationFunction(vsfm_mpp, ieqn, satfunc_type, alpha, lambda, sat_res)
       call VSFMMPPSetSoilPermeability(vsfm_mpp, ieqn, perm, perm, perm)

       deallocate(por          )
       deallocate(sat_res      )
       deallocate(lambda       )
       deallocate(alpha        )
       deallocate(perm         )
       deallocate(satfunc_type )
       
       allocate(ss_auxvar_value(1))

       ss_auxvar_value(:) = 10.d0
       call VSFMMPPSetSourceSinkAuxVarRealValue(vsfm_mpp, ieqn, &
            VAR_POT_MASS_SINK_EXPONENT, ss_auxvar_value)

       ss_auxvar_value(:) = -1500000.d0
       call VSFMMPPSetSourceSinkAuxVarRealValue(vsfm_mpp, ieqn, &
            VAR_POT_MASS_SINK_PRESSURE, ss_auxvar_value)

       deallocate(ss_auxvar_value)

       !
       ! Root equation
       !
       ieqn = 2
       nz   = nz_root

       allocate (por          (nz))
       allocate (alpha        (nz))
       allocate (lambda       (nz))
       allocate (sat_res      (nz))
       allocate (satfunc_type (nz))
       allocate (perm         (nz))

       por          (:) = 0.d0
       sat_res      (:) = 0.d0
       lambda       (:) = 1.d0/Campbell_b
       alpha        (:) = 1.d-3/(-Campbell_he)
       perm         (:) = Ks/1.d6*8.904156d-4
       satfunc_type (:) = SAT_FUNC_BROOKS_COREY
       
       call VSFMMPPSetSoilPorosity(vsfm_mpp, ieqn, por)
       call VSFMMPPSetSaturationFunction(vsfm_mpp, ieqn, satfunc_type, alpha, lambda, sat_res)
       call VSFMMPPSetSoilPermeability(vsfm_mpp, ieqn, perm, perm, perm)

       deallocate(por          )
       deallocate(sat_res      )
       deallocate(lambda       )
       deallocate(alpha        )
       deallocate(perm         )
       deallocate(satfunc_type )
       
       !
       ! Soil equation
       !
       ieqn = 3
       nz   = nz_soil

       allocate (por          (nz))
       allocate (alpha        (nz))
       allocate (lambda       (nz))
       allocate (sat_res      (nz))
       allocate (satfunc_type (nz))
       allocate (perm         (nz))

       por          (:) = theta_s
       sat_res      (:) = 0.01d0
       lambda       (:) = 1.d0 - 1.d0/VG_n
       alpha        (:) = VG_alpha * 1.d-3
       perm         (:) = Ks/1.d6*8.904156d-4
       satfunc_type (:) = SAT_FUNC_VAN_GENUCHTEN
       
       call VSFMMPPSetSoilPorosity(vsfm_mpp, ieqn, por)
       call VSFMMPPSetSaturationFunction(vsfm_mpp, ieqn, satfunc_type, alpha, lambda, sat_res)
       call VSFMMPPSetSoilPermeability(vsfm_mpp, ieqn, perm, perm, perm)

       deallocate(por          )
       deallocate(sat_res      )
       deallocate(lambda       )
       deallocate(alpha        )
       deallocate(perm         )
       deallocate(satfunc_type )
       
    endif

  end subroutine set_material_properties

  !------------------------------------------------------------------------
  subroutine set_initial_conditions()
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use petscsys
    use petscvec
    use petscmat
    use petscts
    use petscsnes
    use petscdm
    use petscdmda
    !
    implicit none
    !
    PetscReal          :: theta
    PetscReal          :: Se
    PetscInt           :: ii
    PetscReal, pointer :: press_ic(:)
    PetscErrorCode     :: ierr

    call VecGetArrayF90(vsfm_mpp%soe%solver%soln, press_ic, ierr); CHKERRQ(ierr)

    theta = 0.20d0

    Se = theta/theta_s
    do ii = 1, nz_xylem + nz_root + nz_soil
       press_ic(ii) = (Campbell_he * Se**(-Campbell_b))* 1.d3 + 101325.d0
    enddo

    call VecRestoreArrayF90(vsfm_mpp%soe%solver%soln, press_ic, ierr); CHKERRQ(ierr)
    call VecCopy(vsfm_mpp%soe%solver%soln, vsfm_mpp%soe%solver%soln_prev, ierr); CHKERRQ(ierr)
    call VecCopy(vsfm_mpp%soe%solver%soln, vsfm_mpp%soe%solver%soln_prev_clm, ierr); CHKERRQ(ierr)
    
  end subroutine set_initial_conditions

  !------------------------------------------------------------------------
  subroutine set_bondary_conditions(time)
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : AUXVAR_BC, VAR_BC_SS_CONDITION
    use MultiPhysicsProbConstants , only : AUXVAR_SS
    !
    implicit none
    !
    PetscReal          :: time
    !
    PetscReal, pointer :: ss_value   (:)
    PetscInt           :: soe_auxvar_id
    PetscReal          :: TimeOfDay
    PetscReal          :: fi
    PetscReal          :: ETp
    PetscReal          :: tp

    TimeOfDay = mod(time,(3600.d0*24.d0))/3600.d0
    fi        = 0.9d0
    ETp       = 5.55555555556d-05
    
    tp = fi * ETp * 2.3d0 * (0.05d0 + sin(0.0175d0 * 7.5d0 * TimeOfDay))** 4.d0
    
    allocate(ss_value(1))
    ss_value(:) = -tp 

    soe_auxvar_id = 1
    call vsfm_mpp%soe%SetDataFromCLM(AUXVAR_SS,  &
         VAR_BC_SS_CONDITION, soe_auxvar_id, ss_value)

    deallocate(ss_value   )

  end subroutine set_bondary_conditions

  !------------------------------------------------------------------------
  subroutine set_conn_flux_type()
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetAuxVarConnRealValue 
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetAuxVarConnIntValue
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSaturationFunctionAuxVarConn
    use MultiPhysicsProbConstants , only : AUXVAR_CONN_INTERNAL
    use MultiPhysicsProbConstants , only : AUXVAR_CONN_BC
    use MultiPhysicsProbConstants , only : AUXVAR_CONN_BC
    use MultiPhysicsProbConstants , only : VAR_FLUX_TYPE
    use MultiPhysicsProbConstants , only : VAR_CONDUCTANCE
    use MultiPhysicsProbConstants , only : CONDUCTANCE_FLUX_TYPE
    use MultiPhysicsProbConstants , only : DARCY_FLUX_TYPE
    use MultiPhysicsProbConstants , only : VAR_CAMPBELL_HE
    use MultiPhysicsProbConstants , only : VAR_CAMPBELL_N
    use SaturationFunction        , only : RELPERM_FUNC_CAMPBELL
    use MultiPhysicsProbConstants , only : AUXVAR_CONN_BC
    use petscsys
    !
    implicit none
    !
    PetscReal , pointer   :: cond_conn_in(:)
    PetscInt  , pointer   :: flux_type_conn_in(:)
    PetscReal , pointer   :: campbell_he_conn_in(:)
    PetscReal , pointer   :: campbell_n_conn_in(:)
    PetscInt  , pointer   :: satfunc_itype_conn_in(:)
    PetscBool , pointer   :: upwind_auxvar_conn_in(:)

    PetscReal , pointer   :: cond_conn_in_x(:)
    PetscInt  , pointer   :: flux_type_conn_in_x(:)
    PetscReal , pointer   :: campbell_he_conn_in_x(:)
    PetscReal , pointer   :: campbell_n_conn_in_x(:)
    PetscInt  , pointer   :: satfunc_itype_conn_in_x(:)
    PetscBool , pointer   :: upwind_auxvar_conn_in_x(:)

    PetscReal , pointer   :: cond_conn_bc_x(:)
    PetscInt  , pointer   :: flux_type_conn_bc_x(:)
    PetscReal , pointer   :: campbell_he_conn_bc_x(:)
    PetscReal , pointer   :: campbell_n_conn_bc_x(:)
    PetscInt  , pointer   :: satfunc_itype_conn_bc_x(:)
    PetscBool , pointer   :: upwind_auxvar_conn_bc_x(:)

    PetscReal , pointer   :: cond_conn_in_r(:)
    PetscInt  , pointer   :: flux_type_conn_in_r(:)
    PetscReal , pointer   :: campbell_he_conn_in_r(:)
    PetscReal , pointer   :: campbell_n_conn_in_r(:)
    PetscInt  , pointer   :: satfunc_itype_conn_in_r(:)
    PetscBool , pointer   :: upwind_auxvar_conn_in_r(:)

    PetscReal , pointer   :: cond_conn_bc_r(:)
    PetscInt  , pointer   :: flux_type_conn_bc_r(:)
    PetscReal , pointer   :: campbell_he_conn_bc_r(:)
    PetscReal , pointer   :: campbell_n_conn_bc_r(:)
    PetscInt  , pointer   :: satfunc_itype_conn_bc_r(:)
    PetscBool , pointer   :: upwind_auxvar_conn_bc_r(:)

    PetscReal , pointer   :: cond_conn_in_s(:)
    PetscInt  , pointer   :: flux_type_conn_in_s(:)
    PetscReal , pointer   :: campbell_he_conn_in_s(:)
    PetscReal , pointer   :: campbell_n_conn_in_s(:)
    PetscInt  , pointer   :: satfunc_itype_conn_in_s(:)
    PetscBool , pointer   :: upwind_auxvar_conn_in_s(:)
    PetscReal , pointer   :: cond_conn_bc_s(:)
    PetscInt  , pointer   :: flux_type_conn_bc_s(:)
    PetscReal , pointer   :: campbell_he_conn_bc_s(:)
    PetscReal , pointer   :: campbell_n_conn_bc_s(:)
    PetscInt  , pointer   :: satfunc_itype_conn_bc_s(:)
    PetscBool , pointer   :: upwind_auxvar_conn_bc_s(:)

    PetscInt              :: ieqn
    PetscInt              :: nconn_in
    PetscInt              :: nconn_in_x
    PetscInt              :: nconn_bc_x
    PetscInt              :: nconn_in_r
    PetscInt              :: nconn_bc_r
    PetscInt              :: nconn_in_s
    PetscInt              :: nconn_bc_s
    PetscInt              :: kk
    PetscReal             :: rootDepth
    PetscReal             :: rootMin
    PetscReal             :: rw
    PetscReal             :: r1
    PetscReal             :: RL
    PetscReal             :: L
    PetscInt              :: nz_loc
    PetscReal             :: dz_loc
    PetscReal             :: Ks
    PetscReal             :: theta,psi_soil,Se,K
    PetscReal , pointer   :: z_int(:)
    PetscReal , pointer   :: Rr(:)
    PetscReal , pointer   :: bz(:)


    if (.not.multi_goveqns_formulation) then
       nconn_in = nz_xylem - 1 + nz_root * 2 + nz_soil - 1

       allocate (cond_conn_in          (nconn_in))
       allocate (flux_type_conn_in     (nconn_in))
       allocate (campbell_he_conn_in   (nconn_in))
       allocate (campbell_n_conn_in    (nconn_in))
       allocate (satfunc_itype_conn_in (nconn_in))
       allocate (upwind_auxvar_conn_in (nconn_in))

       satfunc_itype_conn_in(:) = 0
       upwind_auxvar_conn_in(:) = PETSC_FALSE

       flux_type_conn_in(1                             :nz_xylem - 1 + nz_root * 2) = CONDUCTANCE_FLUX_TYPE
       flux_type_conn_in(nz_xylem - 1 + nz_root * 2 + 1:nconn_in                  ) = DARCY_FLUX_TYPE

       call VSFMMPPSetAuxVarConnIntValue( vsfm_mpp, 1, AUXVAR_CONN_INTERNAL, VAR_FLUX_TYPE, flux_type_conn_in)

       nz_loc = 50
       dz_loc = 1.d0/nz_loc
       allocate(z_int(nz_loc+1))
       allocate(Rr(nz_loc))
       allocate(bz(nz_loc))

       Ks          = 0.001d0 !
       rootDepth   = 0.6d0
       rootMin     = 0.02d0
       rw          = 25000000000.d0
       r1          = 0.001d0

       !pc          = -1500.d0
       RL          = 1/(3.d6 * 1.d6)

       do kk = 1, nz_loc + 1
          z_int(kk) = (kk-1)*dz_loc
       end do

       cond_conn_in(1) = RL

       do kk = 1, nz_loc
          if (z_int(kk) > rootMin .and. z_int(kk) < rootDepth) then
             L = 40000.d0 * (rootDepth - z_int(kk)) / rootDepth
             Rr(kk) = 2.d0 * rw / (L * (z_int(kk+1) - z_int(kk-1)))
             bz(kk) = ((1.d0 - Campbell_n) * log(PI * r1 * r1 * L) / (2 * PI * L * (z_int(kk+1) - z_int(kk-1))))
          else
             Rr(kk) = 0.d0
             bz(kk) = 0.d0
          endif


          if (kk >= 3 .and. kk <= 30) then
             theta    = 0.1d0

             Se       = theta/theta_s
             psi_soil = Campbell_he * Se**(-Campbell_b)
             K        = Ks !* (Campbell_he / psi_soil)**Campbell_n;

             cond_conn_in          (kk-1     ) = 1.d-6/Rr(kk)
             cond_conn_in          (kk-2 + 29) = 1.d-6/(bz(kk)/K)
             campbell_he_conn_in   (kk-2 + 29) = -campbell_he * 1.d3
             campbell_n_conn_in    (kk-2 + 29) = campbell_n
             satfunc_itype_conn_in (kk-2 + 29) = RELPERM_FUNC_CAMPBELL

          endif
       enddo

       call VSFMMPPSetAuxVarConnRealValue(vsfm_mpp, 1, AUXVAR_CONN_INTERNAL, VAR_CONDUCTANCE, cond_conn_in)

       call VSFMMPPSetSaturationFunctionAuxVarConn(vsfm_mpp, 1, AUXVAR_CONN_INTERNAL, upwind_auxvar_conn_in, &
            satfunc_itype_conn_in, campbell_he_conn_in, campbell_n_conn_in)

    else

       Ks          = 0.001d0 !
       rootDepth   = 0.6d0
       rootMin     = 0.02d0
       rw          = 25000000000.d0
       r1          = 0.001d0
       !pc          = -1500.d0
       RL          = 1/(3.d6 * 1.d6)

       !
       ! Xylem: Internal connection auxvars
       !
       nconn_in_x = 1
       nconn_bc_x = nz_root
       nconn_bc_r = nz_root * 2
       nconn_bc_s = nz_root

       allocate (cond_conn_in_x          (nconn_in_x))
       allocate (flux_type_conn_in_x     (nconn_in_x))
       allocate (campbell_he_conn_in_x   (nconn_in_x))
       allocate (campbell_n_conn_in_x    (nconn_in_x))
       allocate (satfunc_itype_conn_in_x (nconn_in_x))
       allocate (upwind_auxvar_conn_in_x (nconn_in_x))

       allocate (cond_conn_bc_x          (nconn_bc_x))
       allocate (flux_type_conn_bc_x     (nconn_bc_x))
       allocate (campbell_he_conn_bc_x   (nconn_bc_x))
       allocate (campbell_n_conn_bc_x    (nconn_bc_x))
       allocate (satfunc_itype_conn_bc_x (nconn_bc_x))
       allocate (upwind_auxvar_conn_bc_x (nconn_bc_x))

       allocate (cond_conn_bc_r          (nconn_bc_r))
       allocate (flux_type_conn_bc_r     (nconn_bc_r))
       allocate (campbell_he_conn_bc_r   (nconn_bc_r))
       allocate (campbell_n_conn_bc_r    (nconn_bc_r))
       allocate (satfunc_itype_conn_bc_r (nconn_bc_r))
       allocate (upwind_auxvar_conn_bc_r (nconn_bc_r))

       allocate (cond_conn_bc_s          (nconn_bc_s))
       allocate (flux_type_conn_bc_s     (nconn_bc_s))
       allocate (campbell_he_conn_bc_s   (nconn_bc_s))
       allocate (campbell_n_conn_bc_s    (nconn_bc_s))
       allocate (satfunc_itype_conn_bc_s (nconn_bc_s))
       allocate (upwind_auxvar_conn_bc_s (nconn_bc_s))

       flux_type_conn_in_x(:)     = CONDUCTANCE_FLUX_TYPE
       cond_conn_in_x(:)          = RL
       satfunc_itype_conn_in_x(:) = 0
       campbell_he_conn_in_x(:)   = 0.d0
       campbell_n_conn_in_x(:)    = 0.d0
       upwind_auxvar_conn_in_x(:) = PETSC_FALSE

       flux_type_conn_bc_x(:)     = CONDUCTANCE_FLUX_TYPE
       satfunc_itype_conn_bc_x(:) = 0
       campbell_he_conn_bc_x(:)   = 0.d0
       campbell_n_conn_bc_x(:)    = 0.d0       
       upwind_auxvar_conn_bc_x(:) = PETSC_FALSE

       flux_type_conn_bc_r(:)     = CONDUCTANCE_FLUX_TYPE
       satfunc_itype_conn_bc_r(:) = 0
       campbell_he_conn_bc_r(:)   = 0.d0
       campbell_n_conn_bc_r(:)    = 0.d0       

       flux_type_conn_bc_s(:)     = CONDUCTANCE_FLUX_TYPE
       satfunc_itype_conn_bc_s(:) = 0
       campbell_he_conn_bc_s(:)   = 0.d0
       campbell_n_conn_bc_s(:)    = 0.d0       


       nz_loc = 50
       dz_loc = 1.d0/nz_loc
       allocate(z_int(nz_loc+1))
       allocate(Rr(nz_loc))
       allocate(bz(nz_loc))

       do kk = 1, nz_loc + 1
          z_int(kk) = (kk-1)*dz_loc
       end do

       do kk = 1, nz_loc
          if (z_int(kk) > rootMin .and. z_int(kk) < rootDepth) then
             L = 40000.d0 * (rootDepth - z_int(kk)) / rootDepth
             Rr(kk) = 2.d0 * rw / (L * (z_int(kk+1) - z_int(kk-1)))
             bz(kk) = ((1.d0 - Campbell_n) * log(PI * r1 * r1 * L) / (2 * PI * L * (z_int(kk+1) - z_int(kk-1))))
          else
             Rr(kk) = 0.d0
             bz(kk) = 0.d0
          endif


          if (kk >= 3 .and. kk <= 30) then
             theta    = 0.1d0

             Se       = theta/theta_s
             psi_soil = Campbell_he * Se**(-Campbell_b)
             K        = Ks

             cond_conn_bc_x          (kk-2         ) = 1.d-6/Rr(kk)
             cond_conn_bc_r          (kk-2         ) = 1.d-6/Rr(kk)
             upwind_auxvar_conn_bc_r (kk-2         ) = PETSC_FALSE
             
             cond_conn_bc_r          (kk-2+nz_root ) = 1.d-6/(bz(kk)/K)
             campbell_he_conn_bc_r   (kk-2+nz_root ) = -campbell_he * 1.d3
             campbell_n_conn_bc_r    (kk-2+nz_root ) = campbell_n
             satfunc_itype_conn_bc_r (kk-2+nz_root ) = RELPERM_FUNC_CAMPBELL
             upwind_auxvar_conn_bc_r (kk-2+nz_root ) = PETSC_TRUE
             
             cond_conn_bc_s          (kk-2         ) = 1.d-6/(bz(kk)/K)
             campbell_he_conn_bc_s   (kk-2         ) = -campbell_he * 1.d3
             campbell_n_conn_bc_s    (kk-2         ) = campbell_n
             satfunc_itype_conn_bc_s (kk-2         ) = RELPERM_FUNC_CAMPBELL
             upwind_auxvar_conn_bc_s (kk-2         ) = PETSC_FALSE

          endif
       enddo

       ieqn     = 1
       call VSFMMPPSetAuxVarConnIntValue( vsfm_mpp, ieqn, AUXVAR_CONN_INTERNAL, VAR_FLUX_TYPE, flux_type_conn_in_x)
       call VSFMMPPSetAuxVarConnRealValue(vsfm_mpp, ieqn, AUXVAR_CONN_INTERNAL, VAR_CONDUCTANCE, cond_conn_in_x)
       call VSFMMPPSetSaturationFunctionAuxVarConn(vsfm_mpp, ieqn, AUXVAR_CONN_INTERNAL, upwind_auxvar_conn_in_x, &
            satfunc_itype_conn_in_x, campbell_he_conn_in_x, campbell_n_conn_in_x)

       call VSFMMPPSetAuxVarConnIntValue( vsfm_mpp, ieqn, AUXVAR_CONN_BC, VAR_FLUX_TYPE, flux_type_conn_bc_x)
       call VSFMMPPSetAuxVarConnRealValue(vsfm_mpp, ieqn, AUXVAR_CONN_BC, VAR_CONDUCTANCE, cond_conn_bc_x)
       call VSFMMPPSetSaturationFunctionAuxVarConn(vsfm_mpp, ieqn, AUXVAR_CONN_BC, upwind_auxvar_conn_bc_x, &
            satfunc_itype_conn_bc_x, campbell_he_conn_bc_x, campbell_n_conn_bc_x)

       ieqn     = 2
       call VSFMMPPSetAuxVarConnIntValue( vsfm_mpp, ieqn, AUXVAR_CONN_BC, VAR_FLUX_TYPE, flux_type_conn_bc_r)
       call VSFMMPPSetAuxVarConnRealValue(vsfm_mpp, ieqn, AUXVAR_CONN_BC, VAR_CONDUCTANCE, cond_conn_bc_r)
       call VSFMMPPSetSaturationFunctionAuxVarConn(vsfm_mpp, ieqn, AUXVAR_CONN_BC, upwind_auxvar_conn_bc_r, &
            satfunc_itype_conn_bc_r, campbell_he_conn_bc_r, campbell_n_conn_bc_r)

       ieqn     = 3
       call VSFMMPPSetAuxVarConnIntValue( vsfm_mpp, ieqn, AUXVAR_CONN_BC, VAR_FLUX_TYPE, flux_type_conn_bc_s)
       call VSFMMPPSetAuxVarConnRealValue(vsfm_mpp, ieqn, AUXVAR_CONN_BC, VAR_CONDUCTANCE, cond_conn_bc_s)
       call VSFMMPPSetSaturationFunctionAuxVarConn(vsfm_mpp, ieqn, AUXVAR_CONN_BC, upwind_auxvar_conn_bc_s, &
            satfunc_itype_conn_bc_s, campbell_he_conn_bc_s, campbell_n_conn_bc_s)

    endif

  end subroutine set_conn_flux_type


  !------------------------------------------------------------------------
  subroutine output_regression_vsfm_spac_campbell_problem(filename_base, num_cells)
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants , only : VAR_PRESSURE
    use MultiPhysicsProbConstants , only : VAR_LIQ_SAT
    use regression_mod            , only : regression_type
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

    name = 'liquid_pressure'
    category = 'pressure'
    call vsfm_mpp%soe%GetDataForCLM(AUXVAR_INTERNAL,  &
         VAR_PRESSURE, -1, data)
    call regression%WriteData(name, category, data)

    name = 'liquid_saturation'
    category = 'general'
    call vsfm_mpp%soe%GetDataForCLM(AUXVAR_INTERNAL,  &
         VAR_LIQ_SAT, -1, data)
    call regression%WriteData(name, category, data)

    call regression%CloseOutput()
    
    deallocate(data)

  end subroutine output_regression_vsfm_spac_campbell_problem

end module vsfm_spac_campbell_problem
