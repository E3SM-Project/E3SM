module vsfm_spac_fetch2_problem

  implicit none

#include <petsc/finclude/petsc.h>
  PetscInt  , parameter :: nx       = 1           ! -
  PetscInt  , parameter :: ny       = 1           ! -
  PetscInt  , parameter :: nz       = 60           ! -
  PetscReal , parameter :: Az       = 9.8715e-3   ! m^2
  PetscReal , parameter :: dx       = 9.8715e-3   ! m
  PetscReal , parameter :: porosity = 0.55d0      ! -
  PetscReal , parameter :: dy       = 1.d0        ! m
  PetscReal , parameter :: dz       = 0.2d0       ! m
  PetscReal , parameter :: Acrown   = 28.8415d0   ! m^2
  PetscReal , parameter :: phis50   = -1.0d6      ! Pa
  PetscReal , parameter :: phi50    = -2.1048d6   ! Pa
  PetscReal , parameter :: phi88    = -0.5876d6   ! Pa
  PetscReal , parameter :: c1       = 1.278d6     ! Pa
  PetscReal , parameter :: c2       = 2.8299d0    ! -
  PetscReal , parameter :: c3       = 5.d0        ! -
  PetscReal , parameter :: kmax     = 9.8513d-7   ! s^{-1}
  PetscReal , parameter :: vis      = 8.904156d-4 ! Pa s
  PetscReal , parameter :: rho      = 1000.d0     ! kg m^{-3}
  PetscReal , parameter :: grav     = 9.81        ! m s^{-2}
  PetscInt              :: ncells_local
  PetscInt              :: ncells_ghost
  PetscReal             :: Campbell_n
  PetscReal             :: Campbell_b
  PetscReal             :: Campbell_he
  PetscReal             :: theta_s
  PetscReal             :: VG_alpha
  PetscReal             :: VG_n
  PetscReal , parameter :: PI  = 4 * atan (1.0_8)

  public :: run_vsfm_spac_fetch2_problem
  
contains

  subroutine run_vsfm_spac_fetch2_problem()
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
    PetscBool          :: print_actual_ET
    character(len=256) :: string
    character(len=256) :: output_suffix
    PetscViewer        :: viewer

    ! Set default settings
    dtime                  = 180.d0
    nstep                  = 24
    save_initial_soln      = PETSC_FALSE
    save_final_soln        = PETSC_FALSE
    print_actual_ET        = PETSC_FALSE
    output_suffix          = ''

    ! Get some command line options

    call PetscOptionsGetInt    (PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-nstep             ',nstep             ,flg,ierr)
    call PetscOptionsGetBool   (PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-save_initial_soln ',save_initial_soln ,flg,ierr)
    call PetscOptionsGetBool   (PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-save_final_soln   ',save_final_soln   ,flg,ierr)
    call PetscOptionsGetBool   (PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-print_actual_ET   ',print_actual_ET   ,flg,ierr)
    call PetscOptionsGetString (PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-output_suffix     ',output_suffix     ,flg,ierr)
    
    ! Initialize the problem
    call Init()

    time = 0.d0

    do istep = 1, nstep

       call set_bondary_conditions(istep)
       time = time + dtime

       ! Run the model
       call vsfm_mpp%soe%StepDT(dtime, istep, &
            converged, converged_reason, ierr); CHKERRQ(ierr)

       if (print_actual_ET) call diagnose_actual_sink()

    end do

  end subroutine run_vsfm_spac_fetch2_problem
  
  !------------------------------------------------------------------------
  subroutine Init()
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : MPP_VSFM_SNES_CLM
    use SystemOfEquationsVSFMType , only : VSFMSOEUpdateConnections
    !
    implicit none
    !

    ! 1. Initialize the multi-physics-problem (MPP)
    call initialize_mpp()

    ! 2. Add all meshes needed for the MPP
    call add_mesh()

    ! 3. Add all governing equations
    call add_goveqn()

    ! 4. Add boundary and source-sink conditions to all governing equations
    call add_conditions_to_goveqns()

    ! 5. Allocate memory to hold auxvars
    call allocate_auxvars()

    ! 6. Setup the MPP
    call vsfm_mpp%SetupProblem()

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
    !
#include <petsc/finclude/petsc.h>
    ! !USES:
    use MultiPhysicsProbConstants , only : MPP_VSFM_SNES_CLM
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
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
    call vsfm_mpp%SetName    ('Tree level hydrodynamics')
    call vsfm_mpp%SetID      (MPP_VSFM_SNES_CLM)
    call vsfm_mpp%SetMPIRank (iam)

  end subroutine initialize_mpp

  !------------------------------------------------------------------------
  subroutine add_mesh()
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
    PetscInt            :: imesh, kk
    PetscInt            :: nlev
    PetscInt            :: iconn, vert_nconn

    PetscInt            :: ncells_xylem

    PetscInt            :: nconn_xylem

    PetscReal , pointer :: xc_xylem(:)                ! x-position of grid cell [m]
    PetscReal , pointer :: yc_xylem(:)                ! y-position of grid cell [m]
    PetscReal , pointer :: zc_xylem(:)                ! z-position of grid cell [m]
    PetscReal , pointer :: dx_xylem(:)                ! layer thickness of grid cell [m]
    PetscReal , pointer :: dy_xylem(:)                ! layer thickness of grid cell [m]
    PetscReal , pointer :: dz_xylem(:)                ! layer thickness of grid cell [m]
    PetscReal , pointer :: area_xylem(:)              ! area of grid cell [m^2]
    PetscReal , pointer :: vol_xylem(:)               ! volume of grid cell [m^3]
    PetscInt  , pointer :: filter_xylem(:)            ! 
    
    PetscInt  , pointer :: vert_conn_id_up_xylem(:)   !
    PetscInt  , pointer :: vert_conn_id_dn_xylem(:)   !
    PetscReal , pointer :: vert_conn_dist_up_xylem(:) !
    PetscReal , pointer :: vert_conn_dist_dn_xylem(:) !
    PetscReal , pointer :: vert_conn_area_xylem(:)    !
    PetscInt  , pointer :: vert_conn_type_xylem(:)    !

    PetscErrorCode      :: ierr

    ncells_ghost = 0
    ncells_xylem = nz

    allocate(xc_xylem     (ncells_xylem ))
    allocate(yc_xylem     (ncells_xylem ))
    allocate(zc_xylem     (ncells_xylem ))
    allocate(dx_xylem     (ncells_xylem ))
    allocate(dy_xylem     (ncells_xylem ))
    allocate(dz_xylem     (ncells_xylem ))
    allocate(area_xylem   (ncells_xylem ))
    allocate(filter_xylem (ncells_xylem ))
    allocate(vol_xylem    (ncells_xylem ))
    
    filter_xylem (:) = 1
    area_xylem   (:) = dx * dy
    dx_xylem     (:) = dx
    dy_xylem     (:) = dy
    dz_xylem     (:) = dz
    xc_xylem     (:) = dx/2.d0
    yc_xylem     (:) = dy/2.d0
    vol_xylem    (:) = 1.d0/50.d0

    zc_xylem(1)  = 0.17d0
    vol_xylem(1) = area_xylem(1) * dz * porosity
    do kk = 2, nz
       zc_xylem(kk)  = -(dz/2.d0 + dz * (kk - 1))
       vol_xylem(kk) = area_xylem(kk) * dz * porosity
    enddo

    call mpp_varpar_set_nlevsoi(nz)
    call mpp_varpar_set_nlevgrnd(nz)

    !
    ! Set up the meshes
    !    
    call vsfm_mpp%SetNumMeshes(1)

    ! Xylem Mesh
    imesh        = 1
    call vsfm_mpp%MeshSetName                (imesh, 'Xylem mesh')
    call vsfm_mpp%MeshSetOrientation         (imesh, MESH_ALONG_GRAVITY)
    call vsfm_mpp%MeshSetID                  (imesh, MESH_CLM_SOIL_COL)
    call vsfm_mpp%MeshSetDimensions          (imesh, ncells_xylem, ncells_ghost, nz)

    call vsfm_mpp%MeshSetGridCellFilter      (imesh, filter_xylem)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_XC     , xc_xylem)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_YC     , yc_xylem)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_ZC     , zc_xylem)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_DX     , dx_xylem)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_DY     , dy_xylem)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_DZ     , dz_xylem)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_AREA   , area_xylem)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_VOLUME , vol_xylem)

    nconn_xylem = nz - 1

    allocate (vert_conn_id_up_xylem   (nconn_xylem))
    allocate (vert_conn_id_dn_xylem   (nconn_xylem))
    allocate (vert_conn_dist_up_xylem (nconn_xylem))
    allocate (vert_conn_dist_dn_xylem (nconn_xylem))
    allocate (vert_conn_area_xylem    (nconn_xylem))
    allocate (vert_conn_type_xylem    (nconn_xylem))

    iconn = 0
    do kk = 1, nz - 1
       iconn = iconn + 1
       vert_conn_id_up_xylem(iconn)   = kk 
       vert_conn_id_dn_xylem(iconn)   = kk + 1
       vert_conn_dist_up_xylem(iconn) = 0.5d0*dz
       vert_conn_dist_dn_xylem(iconn) = 0.5d0*dz
       vert_conn_area_xylem(iconn)    = area_xylem(kk)
       vert_conn_type_xylem(iconn)    = CONN_VERTICAL
    end do

    imesh = 1
    call vsfm_mpp%MeshSetConnectionSet(imesh, CONN_SET_INTERNAL, &
         nconn_xylem,  vert_conn_id_up_xylem, vert_conn_id_dn_xylem,          &
         vert_conn_dist_up_xylem, vert_conn_dist_dn_xylem,  vert_conn_area_xylem,  &
         vert_conn_type_xylem)

    deallocate (xc_xylem                )
    deallocate (yc_xylem                )
    deallocate (zc_xylem                )
    deallocate (dx_xylem                )
    deallocate (dy_xylem                )
    deallocate (dz_xylem                )
    deallocate (area_xylem              )
    deallocate (filter_xylem            )
    deallocate (vol_xylem               )

    deallocate (vert_conn_id_up_xylem   )
    deallocate (vert_conn_id_dn_xylem   )
    deallocate (vert_conn_dist_up_xylem )
    deallocate (vert_conn_dist_dn_xylem )
    deallocate (vert_conn_area_xylem    )
    deallocate (vert_conn_type_xylem    )

  end subroutine add_mesh

  !------------------------------------------------------------------------

  subroutine add_goveqn()
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

  end subroutine add_goveqn

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
    use MultiPhysicsProbConstants , only : COND_DOWNREG_MASS_RATE_FETCH2
    use MultiPhysicsProbConstants , only : CONN_VERTICAL
    use ConnectionSetType         , only : connection_set_type
    use ConnectionSetType         , only : ConnectionSetDestroy
    use MeshType                  , only : MeshCreateConnectionSet
    use petscsys
    !
    ! !ARGUMENTS
    implicit none
    !
    PetscInt                            :: ieqn

    ieqn       = 1
    call vsfm_mpp%soe%AddConditionInGovEqn(ieqn, COND_SS,   &
         'Potential Mass_Flux', 'kg/s', COND_DOWNREG_MASS_RATE_FETCH2, &
         ALL_CELLS)

    ieqn       = 1
    call vsfm_mpp%soe%AddConditionInGovEqn(ieqn, COND_BC,   &
         'Bottom BC', 'Pa', COND_DIRICHLET, SOIL_BOTTOM_CELLS)

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
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSoilPermeability
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetRelativePermeability
    use MultiPhysicsProbConstants , only : VAR_POT_MASS_SINK_PRESSURE
    use MultiPhysicsProbConstants , only : VAR_POT_MASS_SINK_EXPONENT
    use EOSWaterMod               , only : DENSITY_TGDPB01
    use SaturationFunction        , only : SAT_FUNC_FETCH2
    use SaturationFunction        , only : RELPERM_FUNC_WEIBULL
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
    PetscInt  , pointer   :: vsfm_filter(:)
    PetscInt  , pointer   :: satfunc_type(:)
    PetscInt  , pointer   :: relperm_type(:)
    PetscReal , pointer   :: ss_auxvar_value(:)
    PetscReal , pointer   :: weibull_d(:)
    PetscReal , pointer   :: weibull_c(:)
    PetscInt              :: ieqn
    !-----------------------------------------------------------------------

    allocate (por            (nz))
    allocate (alpha          (nz))
    allocate (lambda         (nz))
    allocate (sat_res        (nz))
    allocate (satfunc_type   (nz))
    allocate (perm           (nz))
    allocate (relperm_type   (nz))
    allocate (weibull_d      (nz))
    allocate (weibull_c      (nz))
    allocate(ss_auxvar_value (nz))

    por          (1 : nz ) = 1.d0
    call VSFMMPPSetSoilPorosity(vsfm_mpp, 1, por)

    satfunc_type (1 : nz ) = SAT_FUNC_FETCH2
    alpha        (1 : nz ) = phi88
    lambda       (1 : nz ) = phi50
    sat_res      (1 : nz ) = 0.d0
    call VSFMMPPSetSaturationFunction(vsfm_mpp, 1, satfunc_type, &
         alpha, lambda, sat_res)


    perm         (1 : nz ) = kmax  * vis / rho
    call VSFMMPPSetSoilPermeability(vsfm_mpp, 1, perm, perm, perm)

    relperm_type (:) = RELPERM_FUNC_WEIBULL
    weibull_d    (:) = c1
    weibull_c    (:) = c2
    
    call VSFMMPPSetRelativePermeability(vsfm_mpp, 1, relperm_type, weibull_d, weibull_c)

    ss_auxvar_value(:) = c3
    call VSFMMPPSetSourceSinkAuxVarRealValue(vsfm_mpp, 1, &
         VAR_POT_MASS_SINK_EXPONENT, ss_auxvar_value)

    ss_auxvar_value(:) = phis50
    call VSFMMPPSetSourceSinkAuxVarRealValue(vsfm_mpp, 1, &
         VAR_POT_MASS_SINK_PRESSURE, ss_auxvar_value)

    deallocate(por          )
    deallocate(sat_res      )
    deallocate(lambda       )
    deallocate(alpha        )
    deallocate(perm         )
    deallocate(satfunc_type )
    deallocate(relperm_type )
    deallocate(weibull_c    )
    deallocate(weibull_d    )
    deallocate(ss_auxvar_value)

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
    PetscReal          :: phi_root_mean_times_beta_s
    PetscInt           :: ii
    PetscReal, pointer :: press_ic(:)
    PetscErrorCode     :: ierr

    phi_root_mean_times_beta_s = 5.831916333333334e+03
    
    call VecGetArrayF90(vsfm_mpp%soe%solver%soln, press_ic, ierr); CHKERRQ(ierr)
    do ii = 1, nz
       press_ic(nz-ii+1) = -phi_root_mean_times_beta_s - rho * grav * (0.17d0 + (ii-1)*dz) + 101325.d0
    enddo
    call VecRestoreArrayF90(vsfm_mpp%soe%solver%soln, press_ic, ierr); CHKERRQ(ierr)

    call VecCopy(vsfm_mpp%soe%solver%soln, vsfm_mpp%soe%solver%soln_prev, ierr); CHKERRQ(ierr)
    call VecCopy(vsfm_mpp%soe%solver%soln, vsfm_mpp%soe%solver%soln_prev_clm, ierr); CHKERRQ(ierr)
    
  end subroutine set_initial_conditions

  !------------------------------------------------------------------------
  subroutine set_bondary_conditions(nstep)
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : AUXVAR_BC, VAR_BC_SS_CONDITION
    use MultiPhysicsProbConstants , only : AUXVAR_SS
    use petscsys
    use petscvec
    !
    implicit none
    !
    PetscInt           :: nstep
    !
    PetscReal, pointer :: ss_value(:)
    PetscReal, pointer :: bc_value(:)
    PetscViewer        :: viewer
    Vec                :: ET
    PetscReal, pointer :: et_p(:)
    PetscErrorCode     :: ierr
    PetscInt           :: soe_auxvar_id
    PetscInt           :: kk

    allocate(ss_value(nz))

    call PetscViewerBinaryOpen(PETSC_COMM_WORLD, '../../src/driver/standalone/vsfm/data/PotentialET.bin', FILE_MODE_READ, viewer, ierr)
    call VecCreate(PETSC_COMM_WORLD, ET, ierr)
    call VecLoad(ET, viewer, ierr)
    call PetscViewerDestroy(viewer, ierr)

    call VecGetArrayF90(ET, et_p, ierr)
    do kk = 1, nz
       ss_value(nz-kk+1) = -et_p(nz*(nstep-1) + kk) * dz
    end do
    call VecRestoreArrayF90(ET, et_p, ierr)
    call VecDestroy(ET, ierr)

    soe_auxvar_id = 1
    call vsfm_mpp%soe%SetDataFromCLM(AUXVAR_SS,  &
         VAR_BC_SS_CONDITION, soe_auxvar_id, ss_value)

    deallocate(ss_value)

    allocate(bc_value(1))
    bc_value(1) = 9.382538342475891e+04
    soe_auxvar_id = 1
    call vsfm_mpp%soe%SetDataFromCLM(AUXVAR_BC,  &
         VAR_BC_SS_CONDITION, soe_auxvar_id, bc_value)
    deallocate(bc_value)

  end subroutine set_bondary_conditions

  !------------------------------------------------------------------------
  subroutine diagnose_actual_sink()
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : AUXVAR_SS
    use MultiPhysicsProbConstants , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants , only : VAR_MASS_FLUX
    use MultiPhysicsProbConstants , only : VAR_PRESSURE
    use petscsys
    use petscvec
    !
    implicit none
    !
    PetscInt           :: nstep
    !
    PetscReal, pointer :: ss_value(:)
    PetscReal, pointer :: press(:)
    PetscInt           :: kk

    allocate(ss_value(nz))
    allocate(press(nz))

    call vsfm_mpp%soe%GetDataForCLM(AUXVAR_SS, VAR_MASS_FLUX, 1, ss_value)
    call vsfm_mpp%soe%GetDataForCLM(AUXVAR_INTERNAL,VAR_PRESSURE, -1, press)
    do kk = 1, nz
       write(*,*)kk,ss_value(kk),press(kk)
    end do
    
    deallocate(ss_value)
    deallocate(press)

  end subroutine diagnose_actual_sink

end module vsfm_spac_fetch2_problem
