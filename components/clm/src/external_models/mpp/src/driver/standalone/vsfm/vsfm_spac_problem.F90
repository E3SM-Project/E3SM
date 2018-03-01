module vsfm_spac_problem

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
  PetscReal             :: Campbell_n
  PetscReal             :: Campbell_b
  PetscReal             :: Campbell_he
  PetscReal             :: theta_s

  public :: run_vsfm_spac_problem
  public :: output_regression_vsfm_spac_problem
  
contains

  subroutine run_vsfm_spac_problem()
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
    PetscInt           :: istep, nstep
    PetscBool          :: flg
    PetscBool          :: save_initial_soln, save_final_soln
    character(len=256) :: string
    character(len=256) :: output_suffix
    PetscViewer        :: viewer

    ! Set default settings
    nz                = 30
    dtime             = 3600.d0
    nstep             = 1
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
    call Init()  

    call set_bondary_conditions()
    
    ! Run the model
    call vsfm_mpp%soe%StepDT(dtime, istep, &
         converged, converged_reason, ierr); CHKERRQ(ierr)

     end subroutine run_vsfm_spac_problem

  !------------------------------------------------------------------------
  subroutine Init()
    !
    use MultiPhysicsProbVSFM , only : vsfm_mpp
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
    call vsfm_mpp%SetupProblem()

    ! 7. Add material properities associated with all governing equations
    call set_material_properties() 

    ! 9. Set flux type
    call set_conn_flux_type()
    
    ! 8. Set initial conditions
    call set_initial_conditions()

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
  subroutine add_meshes()
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

    PetscErrorCode :: ierr

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
    soil_area   (:) = 1.d0
    soil_dx     (:) = 1.d0
    soil_dy     (:) = 1.d0
    soil_dz     (:) = 1.d0/50.d0
    soil_xc     (:) = 1.d0/2.d0
    soil_yc     (:) = 1.d0/2.d0

    do kk = 1,nz
       soil_zc(kk) = dz/2.d0 + dz * (kk - 1)
    enddo

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
    call vsfm_mpp%MeshComputeVolume          (imesh)

    allocate (vert_conn_id_up   (nz-1))
    allocate (vert_conn_id_dn   (nz-1))
    allocate (vert_conn_dist_up (nz-1))
    allocate (vert_conn_dist_dn (nz-1))
    allocate (vert_conn_area    (nz-1))
    allocate (vert_conn_type    (nz-1))

    iconn = 0
    iconn = iconn + 1
    vert_conn_id_up(iconn)   = 1
    vert_conn_id_dn(iconn)   = 2
    vert_conn_dist_up(iconn) = 0.5d0*dz
    vert_conn_dist_dn(iconn) = 0.5d0*dz
    vert_conn_area(iconn)    = soil_area(1)
    vert_conn_type(iconn)    = CONN_VERTICAL

    do kk = 2, nz-1

       iconn = iconn + 1
       vert_conn_id_up(iconn)   = 2
       vert_conn_id_dn(iconn)   = kk + 1
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

  end subroutine add_meshes

  !------------------------------------------------------------------------

  subroutine add_goveqns()
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

  end subroutine add_goveqns

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
    use MultiPhysicsProbConstants , only : COND_DOWNREG_MASS_RATE_CAMPBELL
    use MultiPhysicsProbConstants , only : CONN_VERTICAL
    use ConnectionSetType         , only : connection_set_type
    use ConnectionSetType         , only : ConnectionSetDestroy
    use MeshType                  , only : MeshCreateConnectionSet
    !
    ! !ARGUMENTS
    implicit none
    !
    PetscInt                            :: ieqn, ieqn_other
    PetscInt                            :: kk
    PetscInt                            :: nconn
    PetscInt                            :: ncells_root
    PetscReal                           :: root_len
    PetscReal                 , pointer :: root_len_den(:)
    PetscReal                 , pointer :: root_vol(:)
    PetscReal                 , pointer :: root_surf_area(:)
    PetscInt                  , pointer :: id_up(:)
    PetscInt                  , pointer :: id_dn(:)
    PetscReal                 , pointer :: dist_up(:)
    PetscReal                 , pointer :: dist_dn(:)
    PetscReal                 , pointer :: root_zc(:)           ! z-position of grid cell [m]
    PetscReal                 , pointer :: soil_dz(:)           ! layer thickness of grid cell [m]
    PetscReal                 , pointer :: area(:)
    PetscInt                  , pointer :: itype(:)
    PetscReal                 , pointer :: unit_vec(:,:)
    type(connection_set_type) , pointer :: conn_set


    nconn         = 28

    allocate(id_up    (nconn   ))
    allocate(id_dn    (nconn   ))
    allocate(dist_up  (nconn   ))
    allocate(dist_dn  (nconn   ))
    allocate(area     (nconn   ))
    allocate(itype    (nconn   ))
    allocate(unit_vec (nconn,3 ))

    do kk = 1, nconn
       id_up(kk)      = 0
       id_dn(kk)      = kk + 2
       dist_up(kk)    = 0.d0
       dist_dn(kk)    = 1.d0
       area(kk)       = 1.d0
       unit_vec(kk,1) = -1.d0
       unit_vec(kk,2) = 0.d0
       unit_vec(kk,3) = 0.d0
       itype(kk)      = CONN_VERTICAL
    enddo

    allocate(conn_set)

    ieqn       = 1

    call MeshCreateConnectionSet(vsfm_mpp%meshes(1), &
         nconn, id_up, id_dn, &
         dist_up, dist_dn, area, itype, unit_vec, conn_set)
    
    call vsfm_mpp%soe%AddConditionInGovEqn(ieqn, ss_or_bc_type=COND_BC,   &
         name='Root BC in soil equation', unit='Pa', cond_type=COND_DIRICHLET, &
         region_type=SOIL_TOP_CELLS, &
         conn_set=conn_set)

    call vsfm_mpp%soe%AddConditionInGovEqn(ieqn, COND_SS,   &
         'Potential Mass_Flux', 'kg/s', COND_DOWNREG_MASS_RATE_CAMPBELL, &
         SOIL_BOTTOM_CELLS)


    call ConnectionSetDestroy(conn_set)


    deallocate(id_up          )
    deallocate(id_dn          )
    deallocate(dist_up        )
    deallocate(dist_dn        )
    deallocate(area           )
    deallocate(itype          )
    deallocate(unit_vec       )

  end subroutine add_conditions_to_goveqns

  !------------------------------------------------------------------------
  subroutine allocate_auxvars()
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbVSFM , only : vsfm_mpp
    !
    implicit none

    !
    ! Allocate auxvars
    !
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
    PetscInt  , pointer   :: satfunc_type(:)
    PetscReal , parameter :: vish2o = 0.001002d0    ! [N s/m^2] @ 20 degC
    PetscReal             :: Ks
    PetscInt              :: begc , endc
    integer   , pointer   :: vsfm_filter(:)
    PetscReal, pointer    :: ss_auxvar_value(:)
    !-----------------------------------------------------------------------

    Ks          = 0.001d0 !
    theta_s     = 0.46d0  !
    Campbell_b  = 4.58d0  ! [-]
    Campbell_he = -4.2d0  ! [J kg^{-1}]

    Campbell_n  = 2.d0 + 3.d0/Campbell_b

    begc = 1
    endc = 1

    allocate(por(nz))
    allocate(alpha(nz))
    allocate(lambda(nz))
    allocate(sat_res(nz))
    allocate(satfunc_type(nz))

    por     (:) = 0.d0
    sat_res (:) = 0.d0
    lambda  (:) = 1.d0/Campbell_b
    alpha   (:) = 1.d-3/(-Campbell_he)
    satfunc_type(:) = SAT_FUNC_BROOKS_COREY
    
    call VSFMMPPSetSoilPorosity(vsfm_mpp, 1, por)

    call VSFMMPPSetSaturationFunction(vsfm_mpp, 1, satfunc_type, &
       alpha, lambda, sat_res)

    allocate(ss_auxvar_value(1))
    
    ss_auxvar_value(:) = 10.d0
    call VSFMMPPSetSourceSinkAuxVarRealValue(vsfm_mpp, 1, &
         VAR_POT_MASS_SINK_EXPONENT, ss_auxvar_value)

    ss_auxvar_value(:) = -1500000.d0
    call VSFMMPPSetSourceSinkAuxVarRealValue(vsfm_mpp, 1, &
         VAR_POT_MASS_SINK_PRESSURE, ss_auxvar_value)

    deallocate(ss_auxvar_value )
    deallocate(por             )
    deallocate(alpha           )
    deallocate(sat_res         )
    deallocate(satfunc_type    )

  end subroutine set_material_properties

  !------------------------------------------------------------------------
  subroutine set_initial_conditions()
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    !
    implicit none
    !
    PetscReal :: theta
    PetscReal :: Se
    PetscInt  :: ii
    PetscReal, pointer :: press_ic(:)

    allocate(press_ic(nz))

    theta = 0.30d0

    Se = theta/theta_s
    do ii = 1, nz
       press_ic(ii) = (Campbell_he * Se**(-Campbell_b))* 1.d3 + 101325.d0
    enddo

    call vsfm_mpp%Restart(press_ic)

    deallocate(press_ic)

  end subroutine set_initial_conditions

  !------------------------------------------------------------------------
  subroutine set_bondary_conditions()
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : AUXVAR_BC, VAR_BC_SS_CONDITION
    use MultiPhysicsProbConstants , only : AUXVAR_SS
    !
    implicit none
    !
    PetscReal, pointer :: pressure_bc(:)
    PetscReal, pointer :: ss_value   (:)
    PetscInt           :: soe_auxvar_id
    PetscReal :: theta
    PetscReal :: Se
    PetscInt  :: ii

    allocate(pressure_bc(28))

    theta = 0.30d0

    Se = theta/theta_s
    do ii = 1, 28
       pressure_bc(ii) = Campbell_he * Se**(-Campbell_b) * 1.d3 + 101325.d0
    end do

    soe_auxvar_id = 1
    call vsfm_mpp%soe%SetDataFromCLM(AUXVAR_BC,  &
         VAR_BC_SS_CONDITION, soe_auxvar_id, pressure_bc)

    allocate(ss_value(1))
    ss_value(:) = 7.1875e-10 * 1e3

    soe_auxvar_id = 1
    call vsfm_mpp%soe%SetDataFromCLM(AUXVAR_SS,  &
         VAR_BC_SS_CONDITION, soe_auxvar_id, ss_value)

    deallocate(pressure_bc)
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
    use MultiPhysicsProbConstants , only : VAR_FLUX_TYPE
    use MultiPhysicsProbConstants , only : VAR_CONDUCTANCE
    use MultiPhysicsProbConstants , only : CONDUCTANCE_FLUX_TYPE
    use MultiPhysicsProbConstants , only : VAR_CAMPBELL_HE
    use MultiPhysicsProbConstants , only : VAR_CAMPBELL_N
    use SaturationFunction        , only : RELPERM_FUNC_CAMPBELL
    use MultiPhysicsProbConstants , only : AUXVAR_CONN_BC
    use petscsys
    !
    implicit none
    !
    PetscReal , pointer   :: cond_conn_in(:)
    PetscReal , pointer   :: cond_conn_bc(:)
    PetscInt  , pointer   :: flux_type_conn_in(:)
    PetscInt  , pointer   :: flux_type_conn_bc(:)
    PetscReal , pointer   :: campbell_he_conn_bc(:)
    PetscReal , pointer   :: campbell_n_conn_bc(:)
    PetscInt  , pointer   :: satfunc_itype_conn_bc(:)
    PetscBool , pointer   :: upwind_auxvar_conn_bc(:)
    PetscInt              :: nconn_in
    PetscInt              :: nconn_bc
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
    PetscReal , parameter :: PI  = 4 * atan (1.0_8)


    nconn_in = 29
    nconn_bc = 28

    allocate (cond_conn_in          (nconn_in))
    allocate (flux_type_conn_in     (nconn_in))
    allocate (cond_conn_bc          (nconn_bc))
    allocate (flux_type_conn_bc     (nconn_bc))
    allocate (campbell_he_conn_bc   (nconn_bc))
    allocate (campbell_n_conn_bc    (nconn_bc))
    allocate (satfunc_itype_conn_bc (nconn_bc))
    allocate (upwind_auxvar_conn_bc (nconn_bc))

    flux_type_conn_in(:) = CONDUCTANCE_FLUX_TYPE
    flux_type_conn_bc(:) = CONDUCTANCE_FLUX_TYPE
    
    call VSFMMPPSetAuxVarConnIntValue( vsfm_mpp, 1, AUXVAR_CONN_INTERNAL, VAR_FLUX_TYPE, flux_type_conn_in)
    call VSFMMPPSetAuxVarConnIntValue( vsfm_mpp, 1, AUXVAR_CONN_BC      , VAR_FLUX_TYPE, flux_type_conn_bc)

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
    RL          = 1/(3.d6 * 1.d0)

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
          theta    = 0.30d0

          Se       = theta/theta_s
          psi_soil = Campbell_he * Se**(-Campbell_b)
          K        = Ks !* (Campbell_he / psi_soil)**Campbell_n;

          cond_conn_in          (kk-1) = 1.d-3/Rr(kk)
          cond_conn_bc          (kk-2) = 1.d-3/(bz(kk)/K)
          campbell_he_conn_bc   (kk-2) = -campbell_he * 1.d3
          campbell_n_conn_bc    (kk-2) = campbell_n
          satfunc_itype_conn_bc (kk-2) = RELPERM_FUNC_CAMPBELL
          upwind_auxvar_conn_bc (kk-2) = PETSC_FALSE

       endif
    enddo

    call VSFMMPPSetAuxVarConnRealValue(vsfm_mpp, 1, AUXVAR_CONN_INTERNAL, VAR_CONDUCTANCE, cond_conn_in        )
    call VSFMMPPSetAuxVarConnRealValue(vsfm_mpp, 1, AUXVAR_CONN_BC      , VAR_CONDUCTANCE, cond_conn_bc        )

    call VSFMMPPSetSaturationFunctionAuxVarConn(vsfm_mpp, 1, AUXVAR_CONN_BC, upwind_auxvar_conn_bc, &
         satfunc_itype_conn_bc, campbell_he_conn_bc, campbell_n_conn_bc)
        
    deallocate(cond_conn_in      )
    deallocate(flux_type_conn_in )
    deallocate(cond_conn_bc      )
    deallocate(flux_type_conn_bc )

  end subroutine set_conn_flux_type


  !------------------------------------------------------------------------
  subroutine output_regression_vsfm_spac_problem(filename_base, num_cells)
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

  end subroutine output_regression_vsfm_spac_problem
  
end module vsfm_spac_problem
