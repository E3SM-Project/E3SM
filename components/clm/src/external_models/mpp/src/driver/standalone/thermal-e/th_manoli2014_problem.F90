module th_manoli2014_problem

#include <petsc/finclude/petsc.h>
  use mpp_varctl      , only : iulog
  use mpp_abortutils  , only : endrun
  use mpp_shr_log_mod , only : errMsg => shr_log_errMsg
  use petscsys
  use petscvec
  use petscmat
  use petscts
  use petscdm
  use petscdmda

  implicit none

  PetscInt  , parameter :: nx       = 1
  PetscInt  , parameter :: ny       = 1
  PetscReal , parameter :: x_column = 1.d0
  PetscReal , parameter :: y_column = 1.d0
  PetscReal , parameter :: z_column = 1.d0
  PetscInt              :: nz
  PetscInt              :: ncells_ghost


  PetscInt, parameter  :: nx_soil        =   1
  PetscInt, parameter  :: ny_soil        =   1
  PetscInt, parameter  :: nz_soil        =  50
  PetscInt, parameter  :: nz_soil_top    =   3
  PetscInt, parameter  :: nz_soil_mid    =   3
  PetscInt, parameter  :: nz_soil_bot    =  44

  PetscInt, parameter  :: nx_root        =   1
  PetscInt, parameter  :: ny_root        =   1
  PetscInt, parameter  :: nz_root        =  30

  PetscInt, parameter  :: nx_xylem       =   1
  PetscInt, parameter  :: ny_xylem       =   1
  PetscInt, parameter  :: nz_xylem       = 170

  PetscReal, parameter  :: dx            = 1.d0         ! [m]
  PetscReal, parameter  :: dy            = 1.d0         ! [m]
  PetscReal, parameter  :: dz            = 0.1d0        ! [m]

  PetscReal, parameter  :: dx_xylem      = 0.25d0       ! [m]
  PetscReal, parameter  :: dy_xylem      = 0.25d0       ! [m]
  PetscReal, parameter  :: xylem_height  = 17.d0        ! [m]
  PetscReal, parameter  :: xylem_LAI     = 5.d0         ! [leaf area m^2/ soil area m^2]

  PetscReal, parameter  :: scalez        = 0.025d0

  ! Parameters for root length density = length-of-root/volume-of-soil  [m_root/m^3_soil]
  PetscReal, parameter :: qz             = 9.d0          ! [-]
  PetscReal, parameter :: d              = 3.2d0         ! [m]
  PetscReal, parameter :: rld0           = 4.d4          ! [m/m^3]
  PetscReal, parameter :: root_radius    = 2.d-3         ! [m]

  PetscReal, parameter :: pi              = 3.1416d0
  PetscBool :: single_pde_formulation

  PetscReal, parameter :: perm_xy_top    = 6.83d-11      ! [m^2]
  PetscReal, parameter :: perm_z_top     = 6.83d-11      ! [m^2]
  PetscReal, parameter :: sat_res_top    = 0.06d0        ! [-]
  PetscReal, parameter :: alpha_top      = 0.00005d0     ! [Pa^{-1}]
  PetscReal, parameter :: vg_m_top       = 0.33d0        ! [-]
  PetscReal, parameter :: por_top        = 0.5d0         ! [-]

  PetscReal, parameter :: perm_xy_mid    = 6.83d-11      ! [m^2]
  PetscReal, parameter :: perm_z_mid     = 6.83d-11      ! [m^2]
  PetscReal, parameter :: sat_res_mid    = 0.06d0        ! [-]
  PetscReal, parameter :: alpha_mid      = 0.00005d0     ! [Pa^{-1}]
  PetscReal, parameter :: vg_m_mid       = 0.33d0        ! [-]
  PetscReal, parameter :: por_mid        = 0.5d0         ! [-]

  PetscReal, parameter :: perm_xy_bot    = 6.83d-11      ! [m^2]
  PetscReal, parameter :: perm_z_bot     = 6.83d-11      ! [m^2]
  PetscReal, parameter :: sat_res_bot    = 0.06d0        ! [-]
  PetscReal, parameter :: alpha_bot      = 0.00005d0     ! [Pa^{-1}]
  PetscReal, parameter :: vg_m_bot       = 0.33d0        ! [-]
  PetscReal, parameter :: por_bot        = 0.5d0         ! [-]
  
  PetscReal, parameter :: press_initial  = 3.5355d3      ! [Pa]

  public :: run_th_manoli2014_problem

contains
!------------------------------------------------------------------------

  subroutine run_th_manoli2014_problem()
    !
    use MultiPhysicsProbTH , only : th_mpp
    use mpp_varpar         , only : mpp_varpar_init
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
    dtime                  = 8640.d0
    nstep                  = 3
    save_initial_soln      = PETSC_FALSE
    save_final_soln        = PETSC_FALSE
    single_pde_formulation = PETSC_TRUE
    single_pde_formulation = PETSC_FALSE
    output_suffix          = ''

    call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-dt',dtime,flg,ierr)
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-nstep',nstep,flg,ierr)
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-save_initial_soln',save_initial_soln,flg,ierr)
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-save_final_soln',save_final_soln,flg,ierr)
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-single_pde_formulation',single_pde_formulation,flg,ierr)
    call PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-output_suffix',output_suffix,flg,ierr)

    ! Initialize the problem
    call Init()

    if (save_initial_soln) then
       if (len(trim(adjustl(output_suffix))) ==  0) then
          string = trim(output_suffix) // 'initial_soln.bin'
       else
          string = 'initial_soln_' // trim(output_suffix) // '.bin'
       endif

       call PetscViewerBinaryOpen(PETSC_COMM_SELF,trim(string),FILE_MODE_WRITE,viewer,ierr);CHKERRQ(ierr)
       call VecView(th_mpp%soe%solver%soln,viewer,ierr);CHKERRQ(ierr)
       call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
    endif

    do istep = 1, nstep
       ! Update BC

       ! Run the model
       call th_mpp%soe%StepDT(dtime, istep, &
            converged, converged_reason, ierr); CHKERRQ(ierr)
    enddo

    if (save_final_soln) then
       if (len(trim(adjustl(output_suffix))) ==  0) then
          string = trim(output_suffix) // 'final_soln.bin'
       else
          string = 'final_soln_' // trim(output_suffix) // '.bin'
       endif
       call PetscViewerBinaryOpen(PETSC_COMM_SELF,trim(string),FILE_MODE_WRITE,viewer,ierr);CHKERRQ(ierr)
       call VecView(th_mpp%soe%solver%soln,viewer,ierr);CHKERRQ(ierr)
       call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
    endif

  end subroutine run_th_manoli2014_problem

  !------------------------------------------------------------------------
  subroutine Init()
    !
    use MultiPhysicsProbTH        , only : th_mpp
    use SystemOfEquationsTHType   , only : SOETHUpdateConnections
    use SystemOfEquationsBaseType , only : sysofeqns_base_type
    use SystemOfEquationsTHType   , only : sysofeqns_th_type
    !
    implicit none
    !
    class(sysofeqns_base_type) , pointer :: base_soe

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
    call th_mpp%SetupProblem()
 
    ! 7. Add material properities associated with all governing equations
    call set_material_properties() 

    ! 8. Set initial conditions
    call set_initial_conditions()

    base_soe => th_mpp%soe
    select type(base_soe)
    class is(sysofeqns_th_type)
       call SOETHUpdateConnections(base_soe, th_mpp%id)
    class default
       write(iulog,*) 'Unsupported class type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine initialize_mpp()
    !
    ! !DESCRIPTION:
    ! Initialization VSFM
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : MPP_TH_SNES_CLM
    use MultiPhysicsProbTH        , only : th_mpp
    !
    ! !ARGUMENTS
    implicit none
    !
#include <petsc/finclude/petsc.h>
    !
    PetscInt       :: iam
    PetscErrorCode :: ierr

    call MPI_Comm_rank(MPI_COMM_WORLD, iam, ierr)

    !
    ! Set up the multi-physics problem
    !
    call th_mpp%Init       ()
    call th_mpp%SetName    ('Thermal-Hydrology For SPAC')
    call th_mpp%SetID      (MPP_TH_SNES_CLM)
    call th_mpp%SetMPIRank (iam)

  end subroutine initialize_mpp

  !------------------------------------------------------------------------
  subroutine add_meshes()
    !
    use MultiPhysicsProbTH        , only : th_mpp
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
    use MultiPhysicsProbConstants , only : CONN_HORIZONTAL
    use mpp_varpar                , only : mpp_varpar_set_nlevsoi, mpp_varpar_set_nlevgrnd
    !
    implicit none
    !
#include <petsc/finclude/petsc.h>
    !
    PetscInt            :: imesh, ii, jj, kk
    PetscInt            :: ncells_soil, ncells_root, ncells_xylem
    PetscInt            :: ncells_srx
    PetscInt            :: count
    PetscInt            :: iconn, nconn

    PetscReal , pointer :: soil_xc(:)           ! x-position of grid cell [m]
    PetscReal , pointer :: soil_yc(:)           ! y-position of grid cell [m]
    PetscReal , pointer :: soil_zc(:)           ! z-position of grid cell [m]
    PetscReal , pointer :: soil_dx(:)           ! layer thickness of grid cell [m]
    PetscReal , pointer :: soil_dy(:)           ! layer thickness of grid cell [m]
    PetscReal , pointer :: soil_dz(:)           ! layer thickness of grid cell [m]
    PetscReal , pointer :: soil_area(:)         ! area of grid cell [m^2]
    PetscInt  , pointer :: soil_filter(:)       ! 
    PetscReal , pointer :: soil_vol(:)

    PetscReal , pointer :: root_xc(:)           ! x-position of grid cell [m]
    PetscReal , pointer :: root_yc(:)           ! y-position of grid cell [m]
    PetscReal , pointer :: root_zc(:)           ! z-position of grid cell [m]
    PetscReal , pointer :: root_dx(:)           ! layer thickness of grid cell [m]
    PetscReal , pointer :: root_dy(:)           ! layer thickness of grid cell [m]
    PetscReal , pointer :: root_dz(:)           ! layer thickness of grid cell [m]
    PetscReal , pointer :: root_area(:)         ! area of grid cell [m^2]
    PetscInt  , pointer :: root_filter(:)       ! 
    PetscReal , pointer :: root_vol(:)

    PetscReal , pointer :: xylem_xc(:)           ! x-position of grid cell [m]
    PetscReal , pointer :: xylem_yc(:)           ! y-position of grid cell [m]
    PetscReal , pointer :: xylem_zc(:)           ! z-position of grid cell [m]
    PetscReal , pointer :: xylem_dx(:)           ! layer thickness of grid cell [m]
    PetscReal , pointer :: xylem_dy(:)           ! layer thickness of grid cell [m]
    PetscReal , pointer :: xylem_dz(:)           ! layer thickness of grid cell [m]
    PetscReal , pointer :: xylem_area(:)         ! area of grid cell [m^2]
    PetscInt  , pointer :: xylem_filter(:)       ! 
    PetscReal , pointer :: xylem_vol(:)

    PetscReal , pointer :: srx_xc(:)           ! x-position of grid cell [m]
    PetscReal , pointer :: srx_yc(:)           ! y-position of grid cell [m]
    PetscReal , pointer :: srx_zc(:)           ! z-position of grid cell [m]
    PetscReal , pointer :: srx_dx(:)           ! layer thickness of grid cell [m]
    PetscReal , pointer :: srx_dy(:)           ! layer thickness of grid cell [m]
    PetscReal , pointer :: srx_dz(:)           ! layer thickness of grid cell [m]
    PetscReal , pointer :: srx_area(:)         ! area of grid cell [m^2]
    PetscInt  , pointer :: srx_filter(:)       ! 
    PetscReal , pointer :: srx_vol(:)

    PetscReal , pointer :: root_len_den(:)
    PetscReal , pointer :: root_surf_area(:)
    PetscReal           :: root_len

    PetscInt  , pointer :: soil_conn_id_up(:)   !
    PetscInt  , pointer :: soil_conn_id_dn(:)   !
    PetscReal , pointer :: soil_conn_dist_up(:) !
    PetscReal , pointer :: soil_conn_dist_dn(:) !
    PetscReal , pointer :: soil_conn_area(:)    !
    PetscInt  , pointer :: soil_conn_type(:)    !

    PetscInt  , pointer :: root_conn_id_up(:)   !
    PetscInt  , pointer :: root_conn_id_dn(:)   !
    PetscReal , pointer :: root_conn_dist_up(:) !
    PetscReal , pointer :: root_conn_dist_dn(:) !
    PetscReal , pointer :: root_conn_area(:)    !
    PetscInt  , pointer :: root_conn_type(:)    !

    PetscInt  , pointer :: xylem_conn_id_up(:)   !
    PetscInt  , pointer :: xylem_conn_id_dn(:)   !
    PetscReal , pointer :: xylem_conn_dist_up(:) !
    PetscReal , pointer :: xylem_conn_dist_dn(:) !
    PetscReal , pointer :: xylem_conn_area(:)    !
    PetscInt  , pointer :: xylem_conn_type(:)    !

    PetscInt  , pointer :: srx_conn_id_up(:)   !
    PetscInt  , pointer :: srx_conn_id_dn(:)   !
    PetscReal , pointer :: srx_conn_dist_up(:) !
    PetscReal , pointer :: srx_conn_dist_dn(:) !
    PetscReal , pointer :: srx_conn_area(:)    !
    PetscInt  , pointer :: srx_conn_type(:)    !

    PetscReal , pointer :: zv_x(:)
    PetscReal , pointer :: zv_y(:)
    PetscReal , pointer :: xv_2d(:,:)
    PetscReal , pointer :: yv_2d(:,:)
    PetscReal , pointer :: zv_2d(:,:)
    PetscReal , pointer :: xc_3d(:,:,:)
    PetscReal , pointer :: yc_3d(:,:,:)
    PetscReal , pointer :: zc_3d(:,:,:)
    PetscInt  , pointer :: id_3d(:,:,:)
    PetscReal           :: dist_x, dist_y, dist_z, dist

    PetscErrorCode :: ierr

    call mpp_varpar_set_nlevsoi(nz_soil)
    call mpp_varpar_set_nlevgrnd(nz_soil)

    !
    ! Set up the meshes
    !    
    call th_mpp%SetNumMeshes(3)

    ncells_ghost = 0
    ncells_soil  = nx_soil *ny_soil *nz_soil
    ncells_root  = nx_root *ny_root *nz_root
    ncells_xylem = nx_xylem*ny_xylem*nz_xylem
    ncells_srx   = ncells_soil + ncells_root + ncells_xylem

    allocate(soil_xc            (ncells_soil))
    allocate(soil_yc            (ncells_soil))
    allocate(soil_zc            (ncells_soil))
    allocate(soil_dx            (ncells_soil))
    allocate(soil_dy            (ncells_soil))
    allocate(soil_dz            (ncells_soil))
    allocate(soil_area          (ncells_soil))
    allocate(soil_filter        (ncells_soil))
    allocate(soil_vol           (ncells_soil))

    allocate(root_xc            (ncells_root))
    allocate(root_yc            (ncells_root))
    allocate(root_zc            (ncells_root))
    allocate(root_dx            (ncells_root))
    allocate(root_dy            (ncells_root))
    allocate(root_dz            (ncells_root))
    allocate(root_area          (ncells_root))
    allocate(root_filter        (ncells_root))
    allocate(root_vol           (ncells_root))
    allocate(root_len_den       (ncells_root))
    allocate(root_surf_area     (ncells_root))

    allocate(xylem_xc           (ncells_xylem))
    allocate(xylem_yc           (ncells_xylem))
    allocate(xylem_zc           (ncells_xylem))
    allocate(xylem_dx           (ncells_xylem))
    allocate(xylem_dy           (ncells_xylem))
    allocate(xylem_dz           (ncells_xylem))
    allocate(xylem_area         (ncells_xylem))
    allocate(xylem_filter       (ncells_xylem))
    allocate(xylem_vol          (ncells_xylem))

    allocate(srx_xc             (ncells_soil + ncells_root + ncells_xylem))
    allocate(srx_yc             (ncells_soil + ncells_root + ncells_xylem))
    allocate(srx_zc             (ncells_soil + ncells_root + ncells_xylem))
    allocate(srx_dx             (ncells_soil + ncells_root + ncells_xylem))
    allocate(srx_dy             (ncells_soil + ncells_root + ncells_xylem))
    allocate(srx_dz             (ncells_soil + ncells_root + ncells_xylem))
    allocate(srx_area           (ncells_soil + ncells_root + ncells_xylem))
    allocate(srx_filter         (ncells_soil + ncells_root + ncells_xylem))
    allocate(srx_vol            (ncells_soil + ncells_root + ncells_xylem))

    allocate(soil_conn_id_up    (ncells_soil-1))
    allocate(soil_conn_id_dn    (ncells_soil-1))
    allocate(soil_conn_dist_up  (ncells_soil-1))
    allocate(soil_conn_dist_dn  (ncells_soil-1))
    allocate(soil_conn_area     (ncells_soil-1))
    allocate(soil_conn_type     (ncells_soil-1))

    allocate(root_conn_id_up    (ncells_root-1))
    allocate(root_conn_id_dn    (ncells_root-1))
    allocate(root_conn_dist_up  (ncells_root-1))
    allocate(root_conn_dist_dn  (ncells_root-1))
    allocate(root_conn_area     (ncells_root-1))
    allocate(root_conn_type     (ncells_root-1))

    allocate(xylem_conn_id_up   (ncells_xylem-1))
    allocate(xylem_conn_id_dn   (ncells_xylem-1))
    allocate(xylem_conn_dist_up (ncells_xylem-1))
    allocate(xylem_conn_dist_dn (ncells_xylem-1))
    allocate(xylem_conn_area    (ncells_xylem-1))
    allocate(xylem_conn_type    (ncells_xylem-1))


    allocate(srx_conn_id_up     (ncells_soil-1 + ncells_root-1 + ncells_xylem-1 + ncells_soil + 1))
    allocate(srx_conn_id_dn     (ncells_soil-1 + ncells_root-1 + ncells_xylem-1 + ncells_soil + 1))
    allocate(srx_conn_dist_up   (ncells_soil-1 + ncells_root-1 + ncells_xylem-1 + ncells_soil + 1))
    allocate(srx_conn_dist_dn   (ncells_soil-1 + ncells_root-1 + ncells_xylem-1 + ncells_soil + 1))
    allocate(srx_conn_area      (ncells_soil-1 + ncells_root-1 + ncells_xylem-1 + ncells_soil + 1))
    allocate(srx_conn_type      (ncells_soil-1 + ncells_root-1 + ncells_xylem-1 + ncells_soil + 1))

    ! Soil mesh
    soil_filter (:) = 1
    soil_area   (:) = dx*dy
    soil_dx     (:) = dx
    soil_dy     (:) = dy
    soil_dz     (:) = dz
    soil_xc     (:) = dx/2.d0
    soil_yc     (:) = dy/2.d0
    soil_zc     (:) = 0.d0
    soil_zc     (1) = -dz/2.d0
    do kk = 2,nz_soil
       soil_zc(kk) = soil_zc(kk-1) - dz
    enddo

    root_filter (:) = 1
    root_area   (:) = pi*(root_radius**2.d0)
    root_dx     (:) = dx
    root_dy     (:) = dy
    root_dz     (:) = dz
    root_xc     (:) = dx/2.d0 + root_radius
    root_yc     (:) = dy/2.d0
    root_zc     (:) = 0.d0
    do kk = 1,nz_root
       root_zc(kk) = soil_zc(kk)
    enddo

    xylem_filter (:) = 1
    xylem_area   (:) = dx_xylem*dy_xylem
    xylem_dx     (:) = dx_xylem
    xylem_dy     (:) = dy_xylem
    xylem_dz     (:) = dz
    xylem_xc     (:) = dx/2.d0 + root_radius
    xylem_yc     (:) = dy/2.d0
    xylem_zc     (1) = nz_xylem*dz - dz/2.d0
    kk = 1;
    do kk = 2, nz_xylem
       xylem_zc(kk) = xylem_zc(kk-1) - dz
    enddo

    do kk = 1, nz_soil

       if (kk <= nz_root) then
          root_len_den(kk) = rld0 * (1.d0 - abs(root_zc(kk))/d)*exp(-qz*abs(root_zc(kk))/d)

          root_len = root_len_den(kk)*     &       ! [m/m^3]
                     (dx*dy*soil_dz(kk))           ! [m^3]

          root_surf_area(kk) = 2.d0*pi*root_radius*root_len   ! [m^2]

          root_vol(kk) = pi*(root_radius**2.d0)*root_len      ! [m^3]

          soil_vol(kk) = dx * dy * soil_dz(kk) - root_vol(kk) ! [m^3]
       else
          soil_vol(kk) = dx * dy * soil_dz(kk)
       endif

    end do
    xylem_vol(:) = dx*dy*dz
    soil_vol(:)  = dx*dy*dz
    root_vol(:)  = dx*dy*dz
    root_area(:) = dx*dy
    xylem_area(:) = dx*dy

    srx_filter (:) = 1
    srx_area   (:) = dx*dy
    srx_dx     (:) = dx
    srx_dy     (:) = dy
    srx_dz     (:) = dz
    srx_vol    (:) = dx*dy*dz

    srx_xc     (1                        :ncells_soil                         ) = soil_xc(1:ncells_soil )
    srx_xc     (ncells_soil+1            :ncells_soil+ncells_root             ) = root_xc(1:ncells_root )
    srx_xc     (ncells_soil+ncells_root+1:ncells_soil+ncells_root+ncells_xylem) = xylem_xc(1:ncells_xylem)

    srx_yc     (1                        :ncells_soil                         ) = soil_yc(1:ncells_soil )
    srx_yc     (ncells_soil+1            :ncells_soil+ncells_root             ) = root_yc(1:ncells_root )
    srx_yc     (ncells_soil+ncells_root+1:ncells_soil+ncells_root+ncells_xylem) = xylem_yc(1:ncells_xylem)

    srx_zc     (1                        :ncells_soil                         ) = soil_zc(1:ncells_soil )
    srx_zc     (ncells_soil+1            :ncells_soil+ncells_root             ) = root_zc(1:ncells_root )
    srx_zc     (ncells_soil+ncells_root+1:ncells_soil+ncells_root+ncells_xylem) = xylem_zc(1:ncells_xylem)
    
    ! *********** Soil mesh ***********
    
    if (.not. single_pde_formulation) then

       imesh        = 1

       call th_mpp%MeshSetName                (imesh, 'Soil mesh')
       call th_mpp%MeshSetOrientation         (imesh, MESH_ALONG_GRAVITY)
       call th_mpp%MeshSetID                  (imesh, MESH_CLM_SOIL_COL)
       call th_mpp%MeshSetDimensions          (imesh, ncells_soil, ncells_ghost, nz_soil)

       call th_mpp%MeshSetGridCellFilter      (imesh, soil_filter)
       call th_mpp%MeshSetGeometricAttributes (imesh, VAR_XC   , soil_xc)
       call th_mpp%MeshSetGeometricAttributes (imesh, VAR_YC   , soil_yc)
       call th_mpp%MeshSetGeometricAttributes (imesh, VAR_ZC   , soil_zc)
       call th_mpp%MeshSetGeometricAttributes (imesh, VAR_DX   , soil_dx)
       call th_mpp%MeshSetGeometricAttributes (imesh, VAR_DY   , soil_dy)
       call th_mpp%MeshSetGeometricAttributes (imesh, VAR_DZ   , soil_dz)
       call th_mpp%MeshSetGeometricAttributes (imesh, VAR_AREA , soil_area)
       call th_mpp%MeshSetGeometricAttributes (imesh, VAR_VOLUME  , soil_vol)

       !
       ! Vertical connections
       !
       iconn = 0
       do kk = 1, nz_soil-1

          iconn               = iconn + 1
          soil_conn_id_up(iconn)   = kk
          soil_conn_id_dn(iconn)   = kk + 1
          soil_conn_dist_up(iconn) = 0.5d0*soil_dz(kk  )
          soil_conn_dist_dn(iconn) = 0.5d0*soil_dz(kk+1)
          soil_conn_area(iconn)    = soil_area(kk)
          soil_conn_type(iconn)    = CONN_VERTICAL
       end do
       nconn = iconn

       call th_mpp%MeshSetConnectionSet(imesh, CONN_SET_INTERNAL, &
            nconn,  soil_conn_id_up, soil_conn_id_dn,             &
            soil_conn_dist_up, soil_conn_dist_dn, soil_conn_area, &
            soil_conn_type)

       ! *********** Root mesh *********** 

       imesh        = 2

       ncells_ghost = 0

       call th_mpp%MeshSetName                (imesh, 'Root mesh')
       call th_mpp%MeshSetOrientation         (imesh, MESH_ALONG_GRAVITY)
       call th_mpp%MeshSetID                  (imesh, MESH_SPAC_ROOT_COL)
       call th_mpp%MeshSetDimensions          (imesh, ncells_root, ncells_ghost, nz_root)
       call th_mpp%MeshSetGridCellFilter      (imesh, root_filter)
       call th_mpp%MeshSetGeometricAttributes (imesh, VAR_XC   , root_xc)
       call th_mpp%MeshSetGeometricAttributes (imesh, VAR_YC   , root_yc)
       call th_mpp%MeshSetGeometricAttributes (imesh, VAR_ZC   , root_zc)
       call th_mpp%MeshSetGeometricAttributes (imesh, VAR_DX   , root_dx)
       call th_mpp%MeshSetGeometricAttributes (imesh, VAR_DY   , root_dy)
       call th_mpp%MeshSetGeometricAttributes (imesh, VAR_DZ   , root_dz)
       call th_mpp%MeshSetGeometricAttributes (imesh, VAR_AREA , root_area)
       call th_mpp%MeshSetGeometricAttributes (imesh, VAR_VOLUME  , root_vol)

       !
       ! Vertical connections
       !
       iconn = 0
       do kk = 1, nz_root-1
          iconn                    = iconn + 1
          root_conn_id_up(iconn)   = kk
          root_conn_id_dn(iconn)   = kk + 1
          root_conn_dist_up(iconn) = 0.5d0*root_dz(kk  )
          root_conn_dist_dn(iconn) = 0.5d0*root_dz(kk+1)
          root_conn_area(iconn)    = root_area(kk)
          root_conn_type(iconn)    = CONN_VERTICAL
       end do
       nconn = iconn

       call th_mpp%MeshSetConnectionSet(imesh, CONN_SET_INTERNAL, &
            nconn,  root_conn_id_up, root_conn_id_dn,             &
            root_conn_dist_up, root_conn_dist_dn, root_conn_area, &
            root_conn_type)

       ! *********** Xylem mesh *********** 
       imesh        = 3

       ncells_ghost = 0

       call th_mpp%MeshSetName                (imesh , 'Xylem mesh'                         )
       call th_mpp%MeshSetOrientation         (imesh , MESH_ALONG_GRAVITY                   )
       call th_mpp%MeshSetID                  (imesh , MESH_SPAC_XYLEM_COL)
       call th_mpp%MeshSetDimensions          (imesh , ncells_xylem, ncells_ghost, nz_xylem )
       call th_mpp%MeshSetGridCellFilter      (imesh , xylem_filter                         )
       call th_mpp%MeshSetGeometricAttributes (imesh , VAR_XC     , xylem_xc                )
       call th_mpp%MeshSetGeometricAttributes (imesh , VAR_YC     , xylem_yc                )
       call th_mpp%MeshSetGeometricAttributes (imesh , VAR_ZC     , xylem_zc                )
       call th_mpp%MeshSetGeometricAttributes (imesh , VAR_DX     , xylem_dx                )
       call th_mpp%MeshSetGeometricAttributes (imesh , VAR_DY     , xylem_dy                )
       call th_mpp%MeshSetGeometricAttributes (imesh , VAR_DZ     , xylem_dz                )
       call th_mpp%MeshSetGeometricAttributes (imesh , VAR_AREA   , xylem_area              )
       call th_mpp%MeshSetGeometricAttributes (imesh , VAR_VOLUME , xylem_vol               )

       !
       ! Vertical connections
       !
       iconn = 0
       do kk = 1, nz_xylem-1

          iconn                     = iconn + 1
          xylem_conn_id_up(iconn)   = kk
          xylem_conn_id_dn(iconn)   = kk + 1
          xylem_conn_dist_up(iconn) = 0.5d0*xylem_dz(kk  )
          xylem_conn_dist_dn(iconn) = 0.5d0*xylem_dz(kk+1)
          xylem_conn_area(iconn)    = xylem_area(kk)
          xylem_conn_type(iconn)    = CONN_VERTICAL
       end do
       nconn = iconn

       call th_mpp%MeshSetConnectionSet(imesh, CONN_SET_INTERNAL,     &
            nconn,  xylem_conn_id_up, xylem_conn_id_dn,               &
            xylem_conn_dist_up, xylem_conn_dist_dn,  xylem_conn_area, &
            xylem_conn_type)

    else

       imesh        = 1

       call th_mpp%MeshSetName                (imesh, 'Soil+Root+Xylem mesh')
       call th_mpp%MeshSetOrientation         (imesh, MESH_ALONG_GRAVITY)
       call th_mpp%MeshSetID                  (imesh, MESH_CLM_SOIL_COL)
       call th_mpp%MeshSetDimensions          (imesh, ncells_srx, ncells_ghost, ncells_soil)

       call th_mpp%MeshSetGridCellFilter      (imesh, srx_filter)
       call th_mpp%MeshSetGeometricAttributes (imesh, VAR_XC   , srx_xc)
       call th_mpp%MeshSetGeometricAttributes (imesh, VAR_YC   , srx_yc)
       call th_mpp%MeshSetGeometricAttributes (imesh, VAR_ZC   , srx_zc)
       call th_mpp%MeshSetGeometricAttributes (imesh, VAR_DX   , srx_dx)
       call th_mpp%MeshSetGeometricAttributes (imesh, VAR_DY   , srx_dy)
       call th_mpp%MeshSetGeometricAttributes (imesh, VAR_DZ   , srx_dz)
       call th_mpp%MeshSetGeometricAttributes (imesh, VAR_AREA , srx_area)
       call th_mpp%MeshSetGeometricAttributes (imesh, VAR_VOLUME  , srx_vol)

       !
       ! Vertical connections
       !
       iconn = 0
       do kk = 1, nz_soil-1

          iconn                   = iconn + 1
          srx_conn_id_up(iconn)   = kk
          srx_conn_id_dn(iconn)   = kk + 1
          srx_conn_dist_up(iconn) = 0.5d0*soil_dz(kk  )
          srx_conn_dist_dn(iconn) = 0.5d0*soil_dz(kk+1)
          srx_conn_area(iconn)    = soil_area(kk)
          srx_conn_type(iconn)    = CONN_VERTICAL
       end do

       do kk = 1, nz_root-1
          iconn                    = iconn + 1
          srx_conn_id_up(iconn)   = kk + ncells_soil
          srx_conn_id_dn(iconn)   = kk + 1 + ncells_soil
          srx_conn_dist_up(iconn) = 0.5d0*root_dz(kk  )
          srx_conn_dist_dn(iconn) = 0.5d0*root_dz(kk+1)
          srx_conn_area(iconn)    = root_area(kk)
          srx_conn_type(iconn)    = CONN_VERTICAL
       end do

       do kk = 1, nz_xylem-1

          iconn                   = iconn + 1
          srx_conn_id_up(iconn)   = kk + ncells_soil + ncells_root
          srx_conn_id_dn(iconn)   = kk + 1 + ncells_soil + ncells_root
          srx_conn_dist_up(iconn) = 0.5d0*xylem_dz(kk  )
          srx_conn_dist_dn(iconn) = 0.5d0*xylem_dz(kk+1)
          srx_conn_area(iconn)    = xylem_area(kk)
          srx_conn_type(iconn)    = CONN_VERTICAL
       end do
       
       do kk = 1, nz_root
          iconn                   = iconn + 1
          srx_conn_id_up(iconn)   = kk
          srx_conn_id_dn(iconn)   = kk + nz_soil
          srx_conn_dist_up(iconn) = 0.5d0*root_radius
          srx_conn_dist_dn(iconn) = 0.5d0*root_radius
          srx_conn_area(iconn)    = soil_area(kk)
          srx_conn_type(iconn)    = CONN_VERTICAL
       end do

       iconn                   = iconn + 1
       srx_conn_id_up(iconn)   = nz_xylem + nz_root + nz_soil
       srx_conn_id_dn(iconn)   = nz_soil + 1
       srx_conn_dist_up(iconn) = 0.5d0*xylem_dz(nz_xylem)
       srx_conn_dist_dn(iconn) = 0.5d0*root_dz(1)
       srx_conn_area(iconn)    = xylem_area(kk)
       srx_conn_type(iconn)    = CONN_VERTICAL

       nconn = iconn

       call th_mpp%MeshSetConnectionSet(imesh, CONN_SET_INTERNAL, &
            nconn,  srx_conn_id_up, srx_conn_id_dn,               &
            srx_conn_dist_up, srx_conn_dist_dn, srx_conn_area,    &
            srx_conn_type)
    endif

    deallocate(soil_xc            )
    deallocate(soil_yc            )
    deallocate(soil_zc            )
    deallocate(soil_dx            )
    deallocate(soil_dy            )
    deallocate(soil_dz            )
    deallocate(soil_area          )
    deallocate(soil_filter        )
    deallocate(soil_vol           )

    deallocate(root_xc            )
    deallocate(root_yc            )
    deallocate(root_zc            )
    deallocate(root_dx            )
    deallocate(root_dy            )
    deallocate(root_dz            )
    deallocate(root_area          )
    deallocate(root_filter        )
    deallocate(root_vol           )
    deallocate(root_len_den       )
    deallocate(root_surf_area     )

    deallocate(xylem_xc            )
    deallocate(xylem_yc            )
    deallocate(xylem_zc            )
    deallocate(xylem_dx            )
    deallocate(xylem_dy            )
    deallocate(xylem_dz            )
    deallocate(xylem_area          )
    deallocate(xylem_filter        )
    deallocate(xylem_vol           )

    deallocate(soil_conn_id_up    )
    deallocate(soil_conn_id_dn    )
    deallocate(soil_conn_dist_up  )
    deallocate(soil_conn_dist_dn  )
    deallocate(soil_conn_area     )
    deallocate(soil_conn_type     )
    deallocate(root_conn_id_up    )
    deallocate(root_conn_id_dn    )
    deallocate(root_conn_dist_up  )
    deallocate(root_conn_dist_dn  )
    deallocate(root_conn_area     )
    deallocate(root_conn_type     )
    deallocate(xylem_conn_id_up   )
    deallocate(xylem_conn_id_dn   )
    deallocate(xylem_conn_dist_up )
    deallocate(xylem_conn_dist_dn )
    deallocate(xylem_conn_area    )
    deallocate(xylem_conn_type    )

  end subroutine add_meshes

  !------------------------------------------------------------------------

  subroutine add_goveqns()
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbTH       , only : th_mpp
    use MultiPhysicsProbConstants, only : GE_RE
    use MultiPhysicsProbConstants, only : MESH_CLM_SOIL_COL
    use MultiPhysicsProbConstants, only : MESH_SPAC_ROOT_COL
    use MultiPhysicsProbConstants, only : MESH_SPAC_XYLEM_COL
    use MultiPhysicsProbConstants, only : GE_THERM_SOIL_EBASED
    use MultiPhysicsProbConstants, only : GE_RE
    use MultiPhysicsProbConstants, only : MESH_CLM_THERMAL_SOIL_COL
    !
    !
    ! !ARGUMENTS
    implicit none

    if (.not. single_pde_formulation) then
       call th_mpp%AddGovEqn(GE_RE, 'Mass Equation ODE for Soil', MESH_CLM_SOIL_COL)
       call th_mpp%AddGovEqn(GE_RE, 'Mass Equation ODE for Root', MESH_SPAC_ROOT_COL)
       call th_mpp%AddGovEqn(GE_RE, 'Mass Equation ODE for Xylem', MESH_SPAC_XYLEM_COL)
       call th_mpp%AddGovEqn(GE_THERM_SOIL_EBASED, 'Enthalpy-based ODE for heat transport based for Soil', MESH_CLM_SOIL_COL)
       call th_mpp%AddGovEqn(GE_THERM_SOIL_EBASED, 'Enthalpy-based ODE for heat transport based for Root', MESH_SPAC_ROOT_COL)
       call th_mpp%AddGovEqn(GE_THERM_SOIL_EBASED, 'Enthalpy-based ODE for heat transport based for Xylem', MESH_SPAC_XYLEM_COL)
    else
       call th_mpp%AddGovEqn(GE_RE, 'Richards Equation ODE for Soil+Root+Xylem', MESH_CLM_SOIL_COL)
       call th_mpp%AddGovEqn(GE_THERM_SOIL_EBASED, 'Enthalpy-based ODE for Soil+Root+Xylem', MESH_CLM_SOIL_COL)
    endif

    call th_mpp%SetMeshesOfGoveqns()

  end subroutine add_goveqns

  !------------------------------------------------------------------------
  subroutine add_conditions_to_goveqns()
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbTH        , only : th_mpp
    use MultiPhysicsProbConstants , only : ALL_CELLS
    use MultiPhysicsProbConstants , only : SOIL_TOP_CELLS
    use MultiPhysicsProbConstants , only : SOIL_BOTTOM_CELLS
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_DIRICHLET
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    use MultiPhysicsProbConstants , only : CONN_VERTICAL
    use ConnectionSetType         , only : connection_set_type
    use ConnectionSetType         , only : ConnectionSetDestroy
    use MeshType                  , only : MeshCreateConnectionSet
    !
    ! !ARGUMENTS
    implicit none
    !
    PetscInt                            :: ieqn, num_ieqn_others
    PetscInt                  , pointer :: ieqn_others(:)
    PetscBool                 , pointer :: icoupling_others(:)
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

    if (single_pde_formulation) return

    ncells_root  = nx_root*ny_root*nz_root

    allocate(soil_dz            (ncells_root))
    allocate(root_zc            (ncells_root))
    allocate(root_len_den       (ncells_root))
    allocate(root_surf_area     (ncells_root))

    soil_dz     (:) = dz
    root_zc     (:) = 0.d0
    do kk = 2,nz_root
       root_zc(kk) = root_zc(kk-1) - dz
    enddo

    do kk = 1, nz_root

       root_len_den(kk) = rld0 * (1.d0 - &
            abs(root_zc(kk))/d)*exp(-qz*abs(root_zc(kk))/d)

       root_len = root_len_den(kk)*     &       ! [m/m^3]
                  (dx*dy*soil_dz(kk))           ! [m^3]
       root_surf_area(kk) = 2.d0*pi*root_radius*root_len   ! [m^2]

    end do
    root_surf_area(:) = dx*dy

    nconn         = nz_root

    allocate(id_up    (nconn   ))
    allocate(id_dn    (nconn   ))
    allocate(dist_up  (nconn   ))
    allocate(dist_dn  (nconn   ))
    allocate(area     (nconn   ))
    allocate(itype    (nconn   ))
    allocate(unit_vec (nconn,3 ))

    do kk = 1, nz_root
       id_up(kk)      = 0
       id_dn(kk)      = kk
       dist_up(kk)    = 0.d0
       dist_dn(kk)    = root_radius/2.d0
       area(kk)       = root_surf_area(kk)
       unit_vec(kk,1) = -1.d0
       unit_vec(kk,2) = 0.d0
       unit_vec(kk,3) = 0.d0
       itype(kk)      = CONN_VERTICAL
    enddo

    allocate(conn_set)

    ! BoundaryCondition:
    ! SOIL mass ---> ROOT mass
    !           ---> ROOT temperature
    num_ieqn_others = 2

    allocate(ieqn_others(     num_ieqn_others))
    allocate(icoupling_others(num_ieqn_others))

    ieqn            = 1
    ieqn_others(1)  = 2
    ieqn_others(2)  = 5
    icoupling_others(1) = PETSC_FALSE
    icoupling_others(2) = PETSC_FALSE

    call MeshCreateConnectionSet(th_mpp%meshes(1), &
         nconn, id_up, id_dn, &
         dist_up, dist_dn, area, itype, unit_vec, conn_set)
    
    call th_mpp%soe%AddCouplingBCsInGovEqn(ieqn, &
         name='Root BC in soil mass equation',      &
         unit='Pa',                                 &
         region_type=SOIL_TOP_CELLS,                &
         num_other_goveqs=num_ieqn_others,         &
         id_of_other_goveqs=ieqn_others,            &
         icoupling_of_other_goveqns = icoupling_others, &
         conn_set=conn_set)

    deallocate(ieqn_others     )
    deallocate(icoupling_others)

    ! BoundaryCondition:
    ! SOIL temperature ---> ROOT temperature
    !                  ---> ROOT mass
    num_ieqn_others = 3

    allocate(ieqn_others(     num_ieqn_others))
    allocate(icoupling_others(num_ieqn_others))

    ieqn            = 4
    ieqn_others(1)  = 5
    ieqn_others(2)  = 2
    ieqn_others(3)  = 1
    icoupling_others(1) = PETSC_FALSE
    icoupling_others(2) = PETSC_FALSE
    icoupling_others(3) = PETSC_TRUE

    call MeshCreateConnectionSet(th_mpp%meshes(1), &
         nconn, id_up, id_dn, &
         dist_up, dist_dn, area, itype, unit_vec, conn_set)
    
    call th_mpp%soe%AddCouplingBCsInGovEqn(ieqn, &
         name='Root BC in soil energy equation',    &
         unit='Pa',                                 &
         region_type=SOIL_TOP_CELLS,                &
         num_other_goveqs=num_ieqn_others,         &
         id_of_other_goveqs=ieqn_others,            &
         icoupling_of_other_goveqns = icoupling_others, &
         conn_set=conn_set)

    call ConnectionSetDestroy(conn_set)
    deallocate(ieqn_others     )
    deallocate(icoupling_others)

    ! BoundaryCondition:
    ! Root mass ---> Soil mass
    !           ---> Soil temperature
    num_ieqn_others = 2

    allocate(ieqn_others(     num_ieqn_others))
    allocate(icoupling_others(num_ieqn_others))

    ieqn            = 2
    ieqn_others(1)  = 1
    ieqn_others(2)  = 4
    icoupling_others(1) = PETSC_FALSE
    icoupling_others(2) = PETSC_FALSE

    unit_vec(:,1) = 1.d0

    call MeshCreateConnectionSet(th_mpp%meshes(2), &
         nconn, id_up, id_dn, &
         dist_up, dist_dn, area, itype, unit_vec, conn_set)
    
    call th_mpp%soe%AddCouplingBCsInGovEqn(ieqn, &
         name='Soil BC in root mass equation',      &
         unit='Pa',                                 &
         region_type=SOIL_TOP_CELLS,                &
         num_other_goveqs=num_ieqn_others,         &
         id_of_other_goveqs=ieqn_others,            &
         icoupling_of_other_goveqns = icoupling_others, &
         conn_set=conn_set)

    deallocate(ieqn_others     )
    deallocate(icoupling_others)

    ! BoundaryCondition:
    ! Root temperature ---> Soil temperature
    !                  ---> Soil mass
    num_ieqn_others = 3

    allocate(ieqn_others(     num_ieqn_others))
    allocate(icoupling_others(num_ieqn_others))

    ieqn            = 5
    ieqn_others(1)  = 4
    ieqn_others(2)  = 1
    ieqn_others(3)  = 2
    icoupling_others(1) = PETSC_FALSE
    icoupling_others(2) = PETSC_FALSE
    icoupling_others(3) = PETSC_TRUE

    unit_vec(:,1) = 1.d0

    call MeshCreateConnectionSet(th_mpp%meshes(2), &
         nconn, id_up, id_dn, &
         dist_up, dist_dn, area, itype, unit_vec, conn_set)
    
    call th_mpp%soe%AddCouplingBCsInGovEqn(ieqn, &
         name='Soil BC in root energy equation',    &
         unit='Pa',                                 &
         region_type=SOIL_TOP_CELLS,                &
         num_other_goveqs=num_ieqn_others,         &
         id_of_other_goveqs=ieqn_others,            &
         icoupling_of_other_goveqns = icoupling_others, &
         conn_set=conn_set)

    call ConnectionSetDestroy(conn_set)
    deallocate(ieqn_others     )
    deallocate(icoupling_others)

    ! BoundaryCondition:
    ! Root mass ---> Xylem mass
    !           ---> Xylem temperature
    num_ieqn_others = 2

    allocate(ieqn_others(     num_ieqn_others))
    allocate(icoupling_others(num_ieqn_others))

    ieqn            = 2
    ieqn_others(1)  = 3
    ieqn_others(2)  = 6
    icoupling_others(1) = PETSC_FALSE
    icoupling_others(2) = PETSC_FALSE

    call th_mpp%soe%AddCouplingBCsInGovEqn(ieqn, &
         name='Xylem BC in root mass equation',     &
         unit='Pa',                                 &
         region_type=SOIL_TOP_CELLS,                &
         num_other_goveqs=num_ieqn_others,         &
         id_of_other_goveqs=ieqn_others,            &
         icoupling_of_other_goveqns = icoupling_others)

    deallocate(ieqn_others     )
    deallocate(icoupling_others)

    ! BoundaryCondition:
    ! Root temperature ---> Xylem temperature
    !                  ---> Xylem mass
    num_ieqn_others = 3

    allocate(ieqn_others(     num_ieqn_others))
    allocate(icoupling_others(num_ieqn_others))

    ieqn            = 5
    ieqn_others(1)  = 6
    ieqn_others(2)  = 3
    ieqn_others(3)  = 2
    icoupling_others(1) = PETSC_FALSE
    icoupling_others(2) = PETSC_FALSE
    icoupling_others(3) = PETSC_TRUE

    call th_mpp%soe%AddCouplingBCsInGovEqn(ieqn, &
         name='Xylem BC in root energy equation',   &
         unit='K',                                 &
         region_type=SOIL_TOP_CELLS,                &
         num_other_goveqs=num_ieqn_others,         &
         id_of_other_goveqs=ieqn_others,            &
         icoupling_of_other_goveqns = icoupling_others)

    deallocate(ieqn_others     )
    deallocate(icoupling_others)

    ! BoundaryCondition:
    ! Xylem mass ---> Root mass
    !            ---> Root temperature
    num_ieqn_others = 2

    allocate(ieqn_others(     num_ieqn_others))
    allocate(icoupling_others(num_ieqn_others))

    ieqn            = 3
    ieqn_others(1)  = 2
    ieqn_others(2)  = 5
    icoupling_others(1) = PETSC_FALSE
    icoupling_others(2) = PETSC_FALSE

    call th_mpp%soe%AddCouplingBCsInGovEqn(ieqn, &
         name='Root BC in xylem mass equation',     &
         unit='K',                                  &
         region_type=SOIL_BOTTOM_CELLS,             &
         num_other_goveqs=num_ieqn_others,         &
         id_of_other_goveqs=ieqn_others,            &
         icoupling_of_other_goveqns = icoupling_others)

    deallocate(ieqn_others     )
    deallocate(icoupling_others)

    ! BoundaryCondition:
    ! Xylem temperature ---> Root temperature
    !                   ---> Root mass
    num_ieqn_others = 3

    allocate(ieqn_others(     num_ieqn_others))
    allocate(icoupling_others(num_ieqn_others))

    ieqn            = 6
    ieqn_others(1)  = 5
    ieqn_others(2)  = 2
    ieqn_others(3)  = 3
    icoupling_others(1) = PETSC_FALSE
    icoupling_others(2) = PETSC_FALSE
    icoupling_others(3) = PETSC_TRUE

    call th_mpp%soe%AddCouplingBCsInGovEqn(ieqn, &
         name='Root BC in xylem energy equation',   &
         unit='K',                                  &
         region_type=SOIL_BOTTOM_CELLS,             &
         num_other_goveqs=num_ieqn_others,         &
         id_of_other_goveqs=ieqn_others,            &
         icoupling_of_other_goveqns = icoupling_others)

    deallocate(ieqn_others     )
    deallocate(icoupling_others)

    deallocate(soil_dz        )
    deallocate(root_zc        )
    deallocate(root_len_den   )
    deallocate(root_surf_area )
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
    use MultiPhysicsProbTH       , only : th_mpp
    use MultiPhysicsProbConstants, only : VAR_PRESSURE
    use MultiPhysicsProbConstants, only : VAR_TEMPERATURE
    !
    implicit none
    !
    integer           :: ieqn
    integer           :: nvars_for_coupling
    integer, pointer  :: var_ids_for_coupling(:)
    integer, pointer  :: goveqn_ids_for_coupling(:)
    integer, pointer  :: is_bc(:)

    if (.not. single_pde_formulation) then
       !
       ! SOIL mass ---> ROOT mass
       !           ---> ROOT temperature
       !           ---> SOIL temperuature
       ieqn               = 1
       nvars_for_coupling = 3

       allocate (var_ids_for_coupling    (nvars_for_coupling))
       allocate (goveqn_ids_for_coupling (nvars_for_coupling))
       allocate (is_bc                   (nvars_for_coupling))

       var_ids_for_coupling    (1) = VAR_PRESSURE
       var_ids_for_coupling    (2) = VAR_TEMPERATURE
       var_ids_for_coupling    (3) = VAR_TEMPERATURE
       goveqn_ids_for_coupling (1) = 2
       goveqn_ids_for_coupling (2) = 5
       goveqn_ids_for_coupling (3) = 4
       is_bc                   (1) = 1
       is_bc                   (2) = 1
       is_bc                   (3) = 0

       call th_mpp%GovEqnSetBothCouplingVars(ieqn, nvars_for_coupling, &
            var_ids_for_coupling, goveqn_ids_for_coupling, is_bc)

       deallocate(var_ids_for_coupling   )
       deallocate(goveqn_ids_for_coupling)
       deallocate (is_bc                 )

       !
       ! SOIL temperature ---> ROOT temperature
       !                  ---> ROOT mass
       !                  ---> SOIL mass
       ieqn               = 4
       nvars_for_coupling = 3

       allocate (var_ids_for_coupling    (nvars_for_coupling))
       allocate (goveqn_ids_for_coupling (nvars_for_coupling))
       allocate (is_bc                   (nvars_for_coupling))

       var_ids_for_coupling    (1) = VAR_TEMPERATURE
       var_ids_for_coupling    (2) = VAR_PRESSURE
       var_ids_for_coupling    (3) = VAR_PRESSURE
       goveqn_ids_for_coupling (1) = 5
       goveqn_ids_for_coupling (2) = 2
       goveqn_ids_for_coupling (3) = 1
       is_bc                   (1) = 1
       is_bc                   (2) = 1
       is_bc                   (3) = 0

       call th_mpp%GovEqnSetBothCouplingVars(ieqn, nvars_for_coupling, &
            var_ids_for_coupling, goveqn_ids_for_coupling, is_bc)

       deallocate(var_ids_for_coupling   )
       deallocate(goveqn_ids_for_coupling)
       deallocate (is_bc                 )

       !
       ! ROOT mass ---> SOIL  mass
       !           ---> XYLEM mass
       !           ---> SOIL  temperature
       !           ---> XYLEM temperature
       !           ---> ROOT  temperature
       !
       !
       ieqn               = 2
       nvars_for_coupling = 5

       allocate (var_ids_for_coupling    (nvars_for_coupling))
       allocate (goveqn_ids_for_coupling (nvars_for_coupling))
       allocate (is_bc                   (nvars_for_coupling))

       var_ids_for_coupling    (1) = VAR_PRESSURE
       var_ids_for_coupling    (2) = VAR_PRESSURE
       var_ids_for_coupling    (3) = VAR_TEMPERATURE
       var_ids_for_coupling    (4) = VAR_TEMPERATURE
       var_ids_for_coupling    (5) = VAR_TEMPERATURE
       goveqn_ids_for_coupling (1) = 1
       goveqn_ids_for_coupling (2) = 3
       goveqn_ids_for_coupling (3) = 4
       goveqn_ids_for_coupling (4) = 6
       goveqn_ids_for_coupling (5) = 5
       is_bc                   (1) = 1
       is_bc                   (2) = 1
       is_bc                   (3) = 1
       is_bc                   (4) = 1
       is_bc                   (5) = 0

       call th_mpp%GovEqnSetBothCouplingVars(ieqn, nvars_for_coupling, &
            var_ids_for_coupling, goveqn_ids_for_coupling, is_bc)

       deallocate(var_ids_for_coupling   )
       deallocate(goveqn_ids_for_coupling)
       deallocate (is_bc                 )

       !
       ! ROOT temperature ---> SOIL  temperature
       !                  ---> XYLEM temperature
       !                  ---> SOIL  mass
       !                  ---> XYLEM mass
       !                  ---> ROOT  mass
       !
       !
       ieqn               = 5
       nvars_for_coupling = 5

       allocate (var_ids_for_coupling    (nvars_for_coupling))
       allocate (goveqn_ids_for_coupling (nvars_for_coupling))
       allocate (is_bc                   (nvars_for_coupling))

       var_ids_for_coupling    (1) = VAR_TEMPERATURE
       var_ids_for_coupling    (2) = VAR_TEMPERATURE
       var_ids_for_coupling    (3) = VAR_PRESSURE
       var_ids_for_coupling    (4) = VAR_PRESSURE
       var_ids_for_coupling    (5) = VAR_PRESSURE
       goveqn_ids_for_coupling (1) = 4
       goveqn_ids_for_coupling (2) = 6
       goveqn_ids_for_coupling (3) = 1
       goveqn_ids_for_coupling (4) = 3
       goveqn_ids_for_coupling (5) = 2
       is_bc                   (1) = 1
       is_bc                   (2) = 1
       is_bc                   (3) = 1
       is_bc                   (4) = 1
       is_bc                   (5) = 0

       call th_mpp%GovEqnSetBothCouplingVars(ieqn, nvars_for_coupling, &
            var_ids_for_coupling, goveqn_ids_for_coupling, is_bc)

       deallocate(var_ids_for_coupling   )
       deallocate(goveqn_ids_for_coupling)
       deallocate (is_bc                 )

       !
       ! XYLEM mass ---> ROOT  mass
       !            ---> ROOT  temperature
       !            ---> XYLEM temperature
       !
       ieqn               = 3
       nvars_for_coupling = 3

       allocate (var_ids_for_coupling    (nvars_for_coupling))
       allocate (goveqn_ids_for_coupling (nvars_for_coupling))
       allocate (is_bc                   (nvars_for_coupling))

       var_ids_for_coupling    (1) = VAR_PRESSURE
       var_ids_for_coupling    (2) = VAR_TEMPERATURE
       var_ids_for_coupling    (3) = VAR_TEMPERATURE
       goveqn_ids_for_coupling (1) = 2
       goveqn_ids_for_coupling (2) = 5
       goveqn_ids_for_coupling (3) = 6
       is_bc                   (1) = 1
       is_bc                   (2) = 1
       is_bc                   (3) = 0

       call th_mpp%GovEqnSetBothCouplingVars(ieqn, nvars_for_coupling, &
            var_ids_for_coupling, goveqn_ids_for_coupling, is_bc)

       deallocate(var_ids_for_coupling   )
       deallocate(goveqn_ids_for_coupling)
       deallocate (is_bc                 )

       !
       ! XYLEM temperature ---> ROOT  temperature
       !                   ---> ROOT  mass
       !                   ---> XYLEM mass
       !
       ieqn               = 6
       nvars_for_coupling = 3

       allocate (var_ids_for_coupling    (nvars_for_coupling))
       allocate (goveqn_ids_for_coupling (nvars_for_coupling))
       allocate (is_bc                   (nvars_for_coupling))

       var_ids_for_coupling    (1) = VAR_TEMPERATURE
       var_ids_for_coupling    (2) = VAR_PRESSURE
       var_ids_for_coupling    (3) = VAR_PRESSURE
       goveqn_ids_for_coupling (1) = 5
       goveqn_ids_for_coupling (2) = 2
       goveqn_ids_for_coupling (3) = 3
       is_bc                   (1) = 1
       is_bc                   (2) = 1
       is_bc                   (3) = 0

       call th_mpp%GovEqnSetBothCouplingVars(ieqn, nvars_for_coupling, &
            var_ids_for_coupling, goveqn_ids_for_coupling, is_bc)

       deallocate(var_ids_for_coupling   )
       deallocate(goveqn_ids_for_coupling)
       deallocate (is_bc                 )

    else

       ! mass equation (id=1) needs temperature from energy equation (id=2)
       ieqn = 1
       nvars_for_coupling = 1
       allocate (var_ids_for_coupling    (nvars_for_coupling))
       allocate (goveqn_ids_for_coupling (nvars_for_coupling))
       allocate (is_bc                   (nvars_for_coupling))

       var_ids_for_coupling    (1) = VAR_TEMPERATURE
       goveqn_ids_for_coupling (1) = 2
       is_bc                   (1) = 0

       call th_mpp%GovEqnSetBothCouplingVars(ieqn, nvars_for_coupling, &
            var_ids_for_coupling, goveqn_ids_for_coupling, is_bc)

       deallocate(var_ids_for_coupling   )
       deallocate(goveqn_ids_for_coupling)

       ! energy equation (id=2) needs pressure from energy equation (id=1)
       ieqn = 2
       nvars_for_coupling = 1
       allocate (var_ids_for_coupling    (nvars_for_coupling))
       allocate (goveqn_ids_for_coupling (nvars_for_coupling))

       var_ids_for_coupling    (1) = VAR_PRESSURE
       goveqn_ids_for_coupling (1) = 1
       is_bc                   (1) = 0

       call th_mpp%GovEqnSetBothCouplingVars(ieqn, nvars_for_coupling, &
            var_ids_for_coupling, goveqn_ids_for_coupling, is_bc)

       deallocate(var_ids_for_coupling   )
       deallocate(goveqn_ids_for_coupling)

    endif

    !
    ! Allocate auxvars
    !
    call th_mpp%AllocateAuxVars()

  end subroutine allocate_auxvars

  !------------------------------------------------------------------------
  subroutine set_material_properties()
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbTH               , only : th_mpp
    use MultiPhysicsProbTH               , only : MPPTHSetSoils
    use MultiPhysicsProbConstants        , only : GRAVITY_CONSTANT
    use EOSWaterMod                      , only : DENSITY_TGDPB01
    use mpp_varcon                       , only : denh2o
    use mpp_varcon                       , only : grav
    use GoveqnThermalEnthalpySoilType    , only : goveqn_thermal_enthalpy_soil_type
    use SaturationFunction               , only : SatFunc_Set_VG
    use SaturationFunction               , only : SatFunc_Set_Weibull_RelPerm
    use PorosityFunctionMod              , only : PorosityFunctionSetConstantModel
    use ConditionType                    , only : condition_type
    use ConnectionSetType                , only : connection_set_type
    use MultiPhysicsProbConstants        , only : MESH_CLM_SOIL_COL
    use MultiPhysicsProbConstants        , only : MESH_SPAC_ROOT_COL
    use MultiPhysicsProbConstants        , only : MESH_SPAC_XYLEM_COL
    use GoverningEquationBaseType
    use GoveqnRichardsODEPressureType
    use RichardsODEPressureAuxType
    !
    implicit none
    !
    PetscReal        , parameter :: porosity            = 0.368d0
    PetscReal        , parameter :: lambda              = 0.5d0
    PetscReal        , parameter :: alpha               = 3.4257d-4
    PetscReal        , parameter :: perm                = 8.3913d-12
    !
    PetscReal        , pointer   :: vsfm_watsat(:,:)
    PetscReal        , pointer   :: vsfm_hksat(:,:)
    PetscReal        , pointer   :: vsfm_bsw(:,:)
    PetscReal        , pointer   :: vsfm_sucsat(:,:)
    PetscReal        , pointer   :: vsfm_eff_porosity(:,:)
    PetscReal        , pointer   :: vsfm_residual_sat(:,:)
    PetscReal        , parameter :: vish2o = 0.001002d0    ! [N s/m^2] @ 20 degC
    PetscInt                     :: begc , endc
    integer          , pointer   :: vsfm_filter(:)
    character(len=32)            :: satfunc_type
    !
    class(goveqn_base_type),pointer                     :: cur_goveq
    type (rich_ode_pres_auxvar_type), pointer           :: aux_vars_in(:)  !!< Internal state.
    type (rich_ode_pres_auxvar_type), pointer           :: aux_vars_bc(:)  !!< Boundary conditions.
    type (rich_ode_pres_auxvar_type), pointer           :: aux_vars_ss(:)  !!< Source-sink.
    type(condition_type),pointer                        :: cur_cond
    type(connection_set_type), pointer                  :: cur_conn_set
    PetscInt                                            :: ghosted_id
    PetscInt                                            :: sum_conn
    PetscInt                                            :: iconn
    PetscInt                                            :: ncells
    PetscInt                                            :: neqn
    !-----------------------------------------------------------------------

    begc = 1
    endc = nx*ny

    satfunc_type = 'van_genuchten'


    cur_goveq => th_mpp%soe%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       select type(cur_goveq)
       class is (goveqn_richards_ode_pressure_type)
          call set_material_properties_for_mass_equation(cur_goveq)

       class is (goveqn_thermal_enthalpy_soil_type)
          call set_material_properties_for_energy_equation(cur_goveq)

       class default
          write(*,*)'SPACMPPSetupPetscSNESSetup: Unknown goveq type.'
          stop
       end select

       cur_goveq => cur_goveq%next
    enddo

  end subroutine set_material_properties

  !------------------------------------------------------------------------
  subroutine set_material_properties_for_mass_equation(cur_goveq)
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbTH               , only : th_mpp
    use MultiPhysicsProbTH               , only : MPPTHSetSoils
    use MultiPhysicsProbConstants        , only : GRAVITY_CONSTANT
    use EOSWaterMod                      , only : DENSITY_TGDPB01
    use mpp_varcon                       , only : denh2o
    use mpp_varcon                       , only : grav
    use GoveqnThermalEnthalpySoilType    , only : goveqn_thermal_enthalpy_soil_type
    use SaturationFunction               , only : SatFunc_Set_VG
    use SaturationFunction               , only : SatFunc_Set_Weibull_RelPerm
    use PorosityFunctionMod              , only : PorosityFunctionSetConstantModel
    use ConditionType                    , only : condition_type
    use ConnectionSetType                , only : connection_set_type
    use MultiPhysicsProbConstants        , only : MESH_CLM_SOIL_COL
    use MultiPhysicsProbConstants        , only : MESH_SPAC_ROOT_COL
    use MultiPhysicsProbConstants        , only : MESH_SPAC_XYLEM_COL
    use GoverningEquationBaseType
    use GoveqnRichardsODEPressureType
    use RichardsODEPressureAuxType
    !
    implicit none
    !
    PetscReal        , parameter :: porosity            = 0.368d0
    PetscReal        , parameter :: lambda              = 0.5d0
    PetscReal        , parameter :: alpha               = 3.4257d-4
    PetscReal        , parameter :: perm                = 8.3913d-12
    !
    PetscReal        , pointer   :: vsfm_watsat(:,:)
    PetscReal        , pointer   :: vsfm_hksat(:,:)
    PetscReal        , pointer   :: vsfm_bsw(:,:)
    PetscReal        , pointer   :: vsfm_sucsat(:,:)
    PetscReal        , pointer   :: vsfm_eff_porosity(:,:)
    PetscReal        , pointer   :: vsfm_residual_sat(:,:)
    PetscReal        , parameter :: vish2o = 0.001002d0    ! [N s/m^2] @ 20 degC
    PetscInt                     :: begc , endc
    integer          , pointer   :: vsfm_filter(:)
    character(len=32)            :: satfunc_type
    !
    class(goveqn_richards_ode_pressure_type),pointer    :: cur_goveq
    type (rich_ode_pres_auxvar_type), pointer           :: aux_vars_in(:)  !!< Internal state.
    type (rich_ode_pres_auxvar_type), pointer           :: aux_vars_bc(:)  !!< Boundary conditions.
    type (rich_ode_pres_auxvar_type), pointer           :: aux_vars_ss(:)  !!< Source-sink.
    type(condition_type),pointer                        :: cur_cond
    type(connection_set_type), pointer                  :: cur_conn_set
    PetscInt                                            :: ghosted_id
    PetscInt                                            :: sum_conn
    PetscInt                                            :: iconn
    PetscInt                                            :: ncells
    PetscInt                                            :: neqn

    ! Set properties for internal auxvars
    aux_vars_in => cur_goveq%aux_vars_in

    select case(cur_goveq%GetMeshIType())
    case (MESH_CLM_SOIL_COL)
       do ghosted_id = 1,cur_goveq%mesh%ncells_local

          if (ghosted_id <= nx_soil*ny_soil*nz_soil_bot) then

             aux_vars_in(ghosted_id)%perm(1:2)            = perm_xy_bot
             aux_vars_in(ghosted_id)%perm(3)              = perm_z_bot
             aux_vars_in(ghosted_id)%por                  = por_bot

             call PorosityFunctionSetConstantModel(    &
                  aux_vars_in(ghosted_id)%porParams, &
                  por_bot)

             call SatFunc_Set_VG(aux_vars_in(ghosted_id)%satParams,  &
                  sat_res_bot,                        &
                  alpha_bot,                          &
                  vg_m_bot)

          else if (ghosted_id <= nx_soil*ny_soil*(nz_soil_bot+nz_soil_mid)) then

             aux_vars_in(ghosted_id)%perm(1:2)            = perm_xy_mid
             aux_vars_in(ghosted_id)%perm(3)              = perm_z_mid
             aux_vars_in(ghosted_id)%por                  = por_mid

             call PorosityFunctionSetConstantModel(    &
                  aux_vars_in(ghosted_id)%porParams, &
                  por_mid)

             call SatFunc_Set_VG(aux_vars_in(ghosted_id)%satParams,  &
                  sat_res_mid,                        &
                  alpha_mid,                          &
                  vg_m_mid)
          else

             aux_vars_in(ghosted_id)%perm(1:2)            = perm_xy_top
             aux_vars_in(ghosted_id)%perm(3)              = perm_z_top
             aux_vars_in(ghosted_id)%por                  = por_top

             call PorosityFunctionSetConstantModel(    &
                  aux_vars_in(ghosted_id)%porParams, &
                  por_top)

             call SatFunc_Set_VG(aux_vars_in(ghosted_id)%satParams,  &
                  sat_res_top,                        &
                  alpha_top,                          &
                  vg_m_top)
          end if

          aux_vars_in(ghosted_id)%pressure_prev        = press_initial
       enddo
    case (MESH_SPAC_ROOT_COL , &
         MESH_SPAC_XYLEM_COL)

       do ghosted_id = 1,cur_goveq%mesh%ncells_local
          !                                                [1/s] *[m]*[Ns/m^2]/[kg/m^3]/[m/s^2]
          !                aux_vars_in(ghosted_id)%perm(1)              = cond_plant*dx*8.9d-4/denh2o/9.8068d0
          !                aux_vars_in(ghosted_id)%perm(2)              = cond_plant*dy*8.9d-4/denh2o/9.8068d0
          !                aux_vars_in(ghosted_id)%perm(3)              = cond_plant*dz*8.9d-4/denh2o/9.8068d0
          !                aux_vars_in(ghosted_id)%por                  = por_plant

          aux_vars_in(ghosted_id)%perm(1:2)            = perm_xy_top
          aux_vars_in(ghosted_id)%perm(3)              = perm_z_top
          aux_vars_in(ghosted_id)%por                  = por_top

          call PorosityFunctionSetConstantModel(    &
               aux_vars_in(ghosted_id)%porParams, &
               por_top)

          call SatFunc_Set_VG(aux_vars_in(ghosted_id)%satParams,  &
               sat_res_top,                        &
               alpha_top,                          &
               vg_m_top)

          ! Use Weibull relative perm model
          !call SatFunc_Set_Weibull_RelPerm(aux_vars_in(ghosted_id)%satParams, &
          !       weibull_d, weibull_c)

          aux_vars_in(ghosted_id)%pressure_prev        = press_initial
       enddo
    case default
       write(*,*)'SPACMPPSetupPetscSNESSetup: Unknown mesh_itype.'
       stop
    end select

    ! Set hydraulic properties for boundary-condition auxvars
    aux_vars_bc => cur_goveq%aux_vars_bc
    sum_conn = 0
    cur_cond => cur_goveq%boundary_conditions%first
    do
       if (.not.associated(cur_cond)) exit
       cur_conn_set => cur_cond%conn_set

       do iconn = 1, cur_conn_set%num_connections
          sum_conn = sum_conn + 1
          ghosted_id = cur_conn_set%conn(iconn)%GetIDDn()

          aux_vars_bc(sum_conn)%perm(:)             = aux_vars_in(ghosted_id)%perm(:)
          aux_vars_bc(sum_conn)%por                 = aux_vars_in(ghosted_id)%por
          aux_vars_bc(sum_conn)%satParams           = aux_vars_in(ghosted_id)%satParams
          aux_vars_bc(sum_conn)%porParams           = aux_vars_in(ghosted_id)%porParams

       enddo
       cur_cond => cur_cond%next
    enddo

    ! Set hydraulic properties for source-sink auxvars
    aux_vars_ss => cur_goveq%aux_vars_ss
    sum_conn = 0
    cur_cond => cur_goveq%source_sinks%first
    do
       if (.not.associated(cur_cond)) exit
       cur_conn_set => cur_cond%conn_set

       do iconn = 1, cur_conn_set%num_connections
          sum_conn = sum_conn + 1
          ghosted_id = cur_conn_set%conn(iconn)%GetIDDn()

          aux_vars_ss(sum_conn)%perm(:)             = aux_vars_in(ghosted_id)%perm(:)
          aux_vars_ss(sum_conn)%por                 = aux_vars_in(ghosted_id)%por
          aux_vars_ss(sum_conn)%satParams           = aux_vars_in(ghosted_id)%satParams
          aux_vars_ss(sum_conn)%porParams           = aux_vars_in(ghosted_id)%porParams

       enddo
       cur_cond => cur_cond%next
    enddo

  end subroutine set_material_properties_for_mass_equation

  !------------------------------------------------------------------------
  subroutine set_material_properties_for_energy_equation(cur_goveq)
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbTH               , only : th_mpp
    use MultiPhysicsProbTH               , only : MPPTHSetSoils
    use MultiPhysicsProbConstants        , only : GRAVITY_CONSTANT
    use EOSWaterMod                      , only : DENSITY_TGDPB01
    use mpp_varcon                       , only : denh2o
    use mpp_varcon                       , only : grav
    use GoveqnThermalEnthalpySoilType    , only : goveqn_thermal_enthalpy_soil_type
    use SaturationFunction               , only : SatFunc_Set_VG
    use SaturationFunction               , only : SatFunc_Set_Weibull_RelPerm
    use PorosityFunctionMod              , only : PorosityFunctionSetConstantModel
    use ConditionType                    , only : condition_type
    use ConnectionSetType                , only : connection_set_type
    use MultiPhysicsProbConstants        , only : MESH_CLM_SOIL_COL
    use MultiPhysicsProbConstants        , only : MESH_SPAC_ROOT_COL
    use MultiPhysicsProbConstants        , only : MESH_SPAC_XYLEM_COL
    use ThermalEnthalpySoilAuxType       , only : therm_enthalpy_soil_auxvar_type
    use GoverningEquationBaseType
    !
    implicit none
    !
    PetscReal                                , parameter :: porosity = 0.368d0
    PetscReal                                , parameter :: lambda   = 0.5d0
    PetscReal                                , parameter :: alpha    = 3.4257d-4
    PetscReal                                , parameter :: perm     = 8.3913d-12
    PetscReal                                , parameter :: vish2o   = 0.001002d0    ! [N s/m^2] @ 20 degC
    !
    PetscReal                                , pointer   :: vsfm_watsat(:,:)
    PetscReal                                , pointer   :: vsfm_hksat(:,:)
    PetscReal                                , pointer   :: vsfm_bsw(:,:)
    PetscReal                                , pointer   :: vsfm_sucsat(:,:)
    PetscReal                                , pointer   :: vsfm_eff_porosity(:,:)
    PetscReal                                , pointer   :: vsfm_residual_sat(:,:)
    PetscInt                                             :: begc     , endc
    integer                                  , pointer   :: vsfm_filter(:)
    character(len=32)                                    :: satfunc_type
    !
    class(goveqn_thermal_enthalpy_soil_type) , pointer   :: cur_goveq
    type (therm_enthalpy_soil_auxvar_type)   , pointer   :: aux_vars_in(:)
    type (therm_enthalpy_soil_auxvar_type)   , pointer   :: aux_vars_bc(:)
    type (therm_enthalpy_soil_auxvar_type)   , pointer   :: aux_vars_ss(:)
    type(condition_type)                     , pointer   :: cur_cond
    type(connection_set_type)                , pointer   :: cur_conn_set
    PetscInt                                             :: ghosted_id
    PetscInt                                             :: sum_conn
    PetscInt                                             :: iconn
    PetscInt                                             :: ncells
    PetscInt                                             :: neqn

    ! Set properties for internal auxvars
    aux_vars_in => cur_goveq%aux_vars_in

    select case(cur_goveq%GetMeshIType())
    case (MESH_CLM_SOIL_COL)
       do ghosted_id = 1,cur_goveq%mesh%ncells_local

          if (ghosted_id <= nx_soil*ny_soil*nz_soil_bot) then

             aux_vars_in(ghosted_id)%perm(1:2)            = perm_xy_bot
             aux_vars_in(ghosted_id)%perm(3)              = perm_z_bot
             aux_vars_in(ghosted_id)%por                  = por_bot

             call PorosityFunctionSetConstantModel(    &
                  aux_vars_in(ghosted_id)%porParams, &
                  por_bot)

             call SatFunc_Set_VG(aux_vars_in(ghosted_id)%satParams,  &
                  sat_res_bot,                        &
                  alpha_bot,                          &
                  vg_m_bot)

          else if (ghosted_id <= nx_soil*ny_soil*(nz_soil_bot+nz_soil_mid)) then

             aux_vars_in(ghosted_id)%perm(1:2)            = perm_xy_mid
             aux_vars_in(ghosted_id)%perm(3)              = perm_z_mid
             aux_vars_in(ghosted_id)%por                  = por_mid

             call PorosityFunctionSetConstantModel(    &
                  aux_vars_in(ghosted_id)%porParams, &
                  por_mid)

             call SatFunc_Set_VG(aux_vars_in(ghosted_id)%satParams,  &
                  sat_res_mid,                        &
                  alpha_mid,                          &
                  vg_m_mid)
          else

             aux_vars_in(ghosted_id)%perm(1:2)            = perm_xy_top
             aux_vars_in(ghosted_id)%perm(3)              = perm_z_top
             aux_vars_in(ghosted_id)%por                  = por_top

             call PorosityFunctionSetConstantModel(    &
                  aux_vars_in(ghosted_id)%porParams, &
                  por_top)

             call SatFunc_Set_VG(aux_vars_in(ghosted_id)%satParams,  &
                  sat_res_top,                        &
                  alpha_top,                          &
                  vg_m_top)
          end if

          !          aux_vars_in(ghosted_id)%pressure_prev        = press_initial
          aux_vars_in(ghosted_id)%therm_alpha     = 0.45d0
          aux_vars_in(ghosted_id)%therm_cond_wet  = 1.3d0
          aux_vars_in(ghosted_id)%therm_cond_dry  = 0.25d0
          aux_vars_in(ghosted_id)%heat_cap_soil   = 837.d0
          aux_vars_in(ghosted_id)%den_soil        = 2700.d0
       enddo

    case (MESH_SPAC_ROOT_COL , &
         MESH_SPAC_XYLEM_COL)

       do ghosted_id = 1,cur_goveq%mesh%ncells_local
          !                                                [1/s] *[m]*[Ns/m^2]/[kg/m^3]/[m/s^2]
          !                aux_vars_in(ghosted_id)%perm(1)              = cond_plant*dx*8.9d-4/denh2o/9.8068d0
          !                aux_vars_in(ghosted_id)%perm(2)              = cond_plant*dy*8.9d-4/denh2o/9.8068d0
          !                aux_vars_in(ghosted_id)%perm(3)              = cond_plant*dz*8.9d-4/denh2o/9.8068d0
          !                aux_vars_in(ghosted_id)%por                  = por_plant

          aux_vars_in(ghosted_id)%perm(1:2)            = perm_xy_top
          aux_vars_in(ghosted_id)%perm(3)              = perm_z_top
          aux_vars_in(ghosted_id)%por                  = por_top

          call PorosityFunctionSetConstantModel(    &
               aux_vars_in(ghosted_id)%porParams, &
               por_top)

          call SatFunc_Set_VG(aux_vars_in(ghosted_id)%satParams,  &
               sat_res_top,                        &
               alpha_top,                          &
               vg_m_top)

          ! Use Weibull relative perm model
          !call SatFunc_Set_Weibull_RelPerm(aux_vars_in(ghosted_id)%satParams, &
          !       weibull_d, weibull_c)

!          aux_vars_in(ghosted_id)%pressure_prev        = press_initial
          aux_vars_in(ghosted_id)%therm_alpha     = 0.45d0
          aux_vars_in(ghosted_id)%therm_cond_wet  = 1.3d0
          aux_vars_in(ghosted_id)%therm_cond_dry  = 0.25d0
          aux_vars_in(ghosted_id)%heat_cap_soil   = 837.d0
          aux_vars_in(ghosted_id)%den_soil        = 2700.d0
       enddo
    case default
       write(*,*)'SPACMPPSetupPetscSNESSetup: Unknown mesh_itype.'
       stop
    end select

    ! Set hydraulic properties for boundary-condition auxvars
    aux_vars_bc => cur_goveq%aux_vars_bc
    sum_conn = 0
    cur_cond => cur_goveq%boundary_conditions%first
    do
       if (.not.associated(cur_cond)) exit
       cur_conn_set => cur_cond%conn_set

       do iconn = 1, cur_conn_set%num_connections
          sum_conn = sum_conn + 1
          ghosted_id = cur_conn_set%conn(iconn)%GetIDDn()

          aux_vars_bc(sum_conn)%perm(:)       = aux_vars_in(ghosted_id)%perm(:)
          aux_vars_bc(sum_conn)%por           = aux_vars_in(ghosted_id)%por
          aux_vars_bc(sum_conn)%satParams     = aux_vars_in(ghosted_id)%satParams
          aux_vars_bc(sum_conn)%porParams     = aux_vars_in(ghosted_id)%porParams

         aux_vars_bc(sum_conn)%therm_alpha    = aux_vars_in(ghosted_id)%therm_alpha
         aux_vars_bc(sum_conn)%therm_cond_wet = aux_vars_in(ghosted_id)%therm_cond_wet
         aux_vars_bc(sum_conn)%therm_cond_dry = aux_vars_in(ghosted_id)%therm_cond_dry
         aux_vars_bc(sum_conn)%heat_cap_soil  = aux_vars_in(ghosted_id)%heat_cap_soil
         aux_vars_bc(sum_conn)%den_soil       = aux_vars_in(ghosted_id)%den_soil
       enddo
       cur_cond => cur_cond%next
    enddo

    ! Set hydraulic properties for source-sink auxvars
    aux_vars_ss => cur_goveq%aux_vars_ss
    sum_conn = 0
    cur_cond => cur_goveq%source_sinks%first
    do
       if (.not.associated(cur_cond)) exit
       cur_conn_set => cur_cond%conn_set

       do iconn = 1, cur_conn_set%num_connections
          sum_conn = sum_conn + 1
          ghosted_id = cur_conn_set%conn(iconn)%GetIDDn()

          aux_vars_ss(sum_conn)%perm(:)       = aux_vars_in(ghosted_id)%perm(:)
          aux_vars_ss(sum_conn)%por           = aux_vars_in(ghosted_id)%por
          aux_vars_ss(sum_conn)%satParams     = aux_vars_in(ghosted_id)%satParams
          aux_vars_ss(sum_conn)%porParams     = aux_vars_in(ghosted_id)%porParams

         aux_vars_bc(sum_conn)%therm_alpha    = aux_vars_in(ghosted_id)%therm_alpha
         aux_vars_bc(sum_conn)%therm_cond_wet = aux_vars_in(ghosted_id)%therm_cond_wet
         aux_vars_bc(sum_conn)%therm_cond_dry = aux_vars_in(ghosted_id)%therm_cond_dry
         aux_vars_bc(sum_conn)%heat_cap_soil  = aux_vars_in(ghosted_id)%heat_cap_soil
         aux_vars_bc(sum_conn)%den_soil       = aux_vars_in(ghosted_id)%den_soil

      enddo
       cur_cond => cur_cond%next
    enddo

  end subroutine set_material_properties_for_energy_equation

  !------------------------------------------------------------------------
  subroutine set_initial_conditions()
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbTH, only : th_mpp
    !
    implicit none
    !
    PetscErrorCode                                    :: ierr

    !
    PetscInt            :: nDM
    DM        , pointer :: dms(:)
    Vec       , pointer :: soln_subvecs(:)
    PetscReal , pointer :: v_p(:)
    PetscInt            :: ii
    PetscViewer         :: viewer
    PetscInt            :: soe_auxvar_id

    ! Find number of GEs packed within the SoE
    call DMCompositeGetNumberDM(th_mpp%soe%solver%dm, nDM, ierr)

    ! Get DMs for each GE
    allocate (dms(nDM))
    call DMCompositeGetEntriesArray(th_mpp%soe%solver%dm, dms, ierr)

    ! Allocate vectors for individual GEs
    allocate(soln_subvecs(nDM))

    ! Get solution vectors for individual GEs
    call DMCompositeGetAccessArray(th_mpp%soe%solver%dm, &
         th_mpp%soe%solver%soln, nDM, &
         PETSC_NULL_INTEGER, soln_subvecs, ierr)

    do ii = 1, nDM
       call VecGetArrayF90(soln_subvecs(ii), v_p, ierr)
       if (.not.single_pde_formulation) then
       if (ii < 4) then
          v_p(:) = press_initial
       else
          v_p(:) = 283.15d0
       endif
       else
       if (ii ==1) then
          v_p(:) = press_initial
       else
          v_p(:) = 283.15d0
       endif
       endif
       call VecRestoreArrayF90(soln_subvecs(ii), v_p, ierr)
    enddo

    ! Restore solution vectors for individual GEs
    call DMCompositeRestoreAccessArray(th_mpp%soe%solver%dm, &
         th_mpp%soe%solver%soln, nDM, &
         PETSC_NULL_INTEGER, soln_subvecs, ierr)

    deallocate(dms)
    
    call VecCopy(th_mpp%soe%solver%soln, th_mpp%soe%solver%soln_prev, ierr); CHKERRQ(ierr)
    call VecCopy(th_mpp%soe%solver%soln, th_mpp%soe%solver%soln_prev_clm, ierr); CHKERRQ(ierr)

  end subroutine set_initial_conditions

end module th_manoli2014_problem
