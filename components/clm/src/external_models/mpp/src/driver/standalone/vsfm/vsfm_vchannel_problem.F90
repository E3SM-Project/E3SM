module vsfm_vchannel_problem

  implicit none

#include <petsc/finclude/petsc.h>

  PetscInt  , parameter :: nx       = 20
  PetscInt  , parameter :: ny       = 10
  PetscInt              :: nz
  PetscInt  , parameter :: nz_hex   = 87
  PetscReal , parameter :: dx       = 10.0d0
  PetscReal , parameter :: dy       = 10.0d0
  PetscReal , parameter :: dz       =  0.5d0
  PetscReal , parameter :: slope_x  =  0.1d0
  PetscReal , parameter :: slope_y  =  0.2d0
  PetscBool             :: with_seepage_bc
  PetscBool             :: orthogonal_hex
  PetscInt              :: ncells_local
  PetscInt              :: ncells_ghost
  PetscInt, pointer     :: soil_filter_2d(:,:)
  PetscInt, pointer     :: kk_lower_idx(:)

  public :: run_vsfm_vchannel_problem
  public :: output_regression_vsfm_vchannel_problem

contains

!------------------------------------------------------------------------
  subroutine run_vsfm_vchannel_problem()

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
    dtime             = 8640.d0
    nstep             = 3
    save_initial_soln = PETSC_FALSE
    save_final_soln   = PETSC_FALSE
    with_seepage_bc   = PETSC_FALSE
    orthogonal_hex    = PETSC_FALSE
    output_suffix     = ''
    nz                = 30

    ! Get some command line options

    call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-dt',dtime,flg,ierr)
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-nstep',nstep,flg,ierr)
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-save_initial_soln',save_initial_soln,flg,ierr)
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-save_final_soln',save_final_soln,flg,ierr)
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-with_seepage_bc',with_seepage_bc,flg,ierr)
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-orthogonal_hex',orthogonal_hex,flg,ierr)
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
       call VecView(vsfm_mpp%soe%solver%soln,viewer,ierr);CHKERRQ(ierr)
       call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
    endif

    do istep = 1, nstep
       ! Update BC
       if (with_seepage_bc) call set_bondary_conditions()

       ! Run the model
       call vsfm_mpp%soe%StepDT(dtime, istep, &
            converged, converged_reason, ierr); CHKERRQ(ierr)
    enddo

    if (save_final_soln) then
       if (len(trim(adjustl(output_suffix))) ==  0) then
          string = trim(output_suffix) // 'final_soln.bin'
       else
          string = 'final_soln_' // trim(output_suffix) // '.bin'
       endif
       call PetscViewerBinaryOpen(PETSC_COMM_SELF,trim(string),FILE_MODE_WRITE,viewer,ierr);CHKERRQ(ierr)
       call VecView(vsfm_mpp%soe%solver%soln,viewer,ierr);CHKERRQ(ierr)
       call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
    endif

  end subroutine run_vsfm_vchannel_problem


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
    if (.not.orthogonal_hex) then
       call add_meshes()
    else
       call add_orthogonal_hex_mesh()
    endif

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
    call vsfm_mpp%SetName    ('Variably-Saturated-Flow-Model For V Channel')
    call vsfm_mpp%SetID      (MPP_VSFM_SNES_CLM)
    call vsfm_mpp%SetMPIRank (iam)

  end subroutine initialize_mpp

  !------------------------------------------------------------------------
  subroutine add_meshes()
    !
#include <petsc/finclude/petsc.h>
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
    use MultiPhysicsProbConstants , only : CONN_HORIZONTAL
    use mpp_varpar                , only : mpp_varpar_set_nlevsoi, mpp_varpar_set_nlevgrnd
    use petscsys
    !
    implicit none
    !
    PetscInt           :: imesh, ii, jj, kk
    PetscInt           :: nlev
    PetscInt           :: count
    PetscInt           :: iconn, nconn
    PetscReal, pointer :: soil_xc(:)           ! x-position of grid cell [m]
    PetscReal, pointer :: soil_yc(:)           ! y-position of grid cell [m]
    PetscReal, pointer :: soil_zc(:)           ! z-position of grid cell [m]
    PetscReal, pointer :: soil_dx(:)           ! layer thickness of grid cell [m]
    PetscReal, pointer :: soil_dy(:)           ! layer thickness of grid cell [m]
    PetscReal, pointer :: soil_dz(:)           ! layer thickness of grid cell [m]
    PetscReal, pointer :: soil_area(:)         ! area of grid cell [m^2]
    PetscInt , pointer :: soil_filter(:)       ! 

    PetscInt, pointer  :: conn_id_up(:)   !
    PetscInt, pointer  :: conn_id_dn(:)   !
    PetscReal, pointer :: conn_dist_up(:) !
    PetscReal, pointer :: conn_dist_dn(:) !
    PetscReal, pointer :: conn_area(:)    !
    PetscInt , pointer :: conn_type(:)    !
    PetscReal, pointer :: zv_x(:)
    PetscReal, pointer :: zv_y(:)
    PetscReal, pointer :: xv_2d(:,:)
    PetscReal, pointer :: yv_2d(:,:)
    PetscReal, pointer :: zv_2d(:,:)
    PetscReal, pointer :: xc_3d(:,:,:)
    PetscReal, pointer :: yc_3d(:,:,:)
    PetscReal, pointer :: zc_3d(:,:,:)
    PetscInt, pointer  :: id_3d(:,:,:)
    PetscReal :: dist_x, dist_y, dist_z, dist

    PetscErrorCode :: ierr

    call mpp_varpar_set_nlevsoi(nz)
    call mpp_varpar_set_nlevgrnd(nz)

    imesh        = 1
    nlev         = nz
    ncells_local = nx*ny*nz
    ncells_ghost = 0

    allocate(soil_xc     (ncells_local))
    allocate(soil_yc     (ncells_local))
    allocate(soil_zc     (ncells_local))
    allocate(soil_dx     (ncells_local))
    allocate(soil_dy     (ncells_local))
    allocate(soil_dz     (ncells_local))
    allocate(soil_area   (ncells_local))
    allocate(soil_filter (ncells_local))

    allocate(zv_x(nx+1))
    allocate(zv_y(ny+1))

    allocate(xv_2d(nx+1,ny+1))
    allocate(yv_2d(nx+1,ny+1))
    allocate(zv_2d(nx+1,ny+1))

    allocate(xc_3d(nx,ny,nz))
    allocate(yc_3d(nx,ny,nz))
    allocate(zc_3d(nx,ny,nz))
    allocate(id_3d(nx,ny,nz))

    allocate (conn_id_up   ((nz-1)*nx*ny + nx*(ny-1)*nz + (nx-1)*ny*nz))
    allocate (conn_id_dn   ((nz-1)*nx*ny + nx*(ny-1)*nz + (nx-1)*ny*nz))
    allocate (conn_dist_up ((nz-1)*nx*ny + nx*(ny-1)*nz + (nx-1)*ny*nz))
    allocate (conn_dist_dn ((nz-1)*nx*ny + nx*(ny-1)*nz + (nx-1)*ny*nz))
    allocate (conn_area    ((nz-1)*nx*ny + nx*(ny-1)*nz + (nx-1)*ny*nz))
    allocate (conn_type    ((nz-1)*nx*ny + nx*(ny-1)*nz + (nx-1)*ny*nz))
    allocate(soil_filter_2d(nx*ny,nz_hex))
    allocate(kk_lower_idx(nx*ny))

    soil_filter (:) = 1
    soil_area   (:) = dx*dy
    soil_dx     (:) = dx
    soil_dy     (:) = dy
    soil_dz     (:) = dz
    soil_xc     (:) = dx/2.d0
    soil_xc     (:) = dy/2.d0
    soil_filter_2d(:,:) = 1
    kk_lower_idx(:)     = 1

    ! Compute elevation at vertices in x-direction
    do ii = 1, nx/2+1
       zv_x(ii) = slope_x*dx*(nx/2) - (ii-1)*slope_x*dx
    enddo
    do ii = nx/2+2, nx+1
       zv_x(ii) = (ii-nx/2-1)*slope_x*dx
    enddo

    ! Compute elevation at vertices in y-direction
    do jj = 1, ny+1
       zv_y(jj) = (jj-1)*slope_y*dy
    enddo

    do jj = 1,ny+1
       do ii = 1,nx+1
          zv_2d(ii,jj) = zv_x(ii) + zv_y(jj);
          xv_2d(ii,jj) = (ii-1)*dx;
          yv_2d(ii,jj) = (jj-1)*dy;
       enddo
    enddo

    do kk = 1,nz
       do jj = 1,ny
          do ii = 1,nx
             xc_3d(ii,jj,kk) = (xv_2d(ii,jj) + xv_2d(ii+1,jj) + xv_2d(ii,jj+1) + xv_2d(ii+1,jj+1))/4.d0
             yc_3d(ii,jj,kk) = (yv_2d(ii,jj) + yv_2d(ii+1,jj) + yv_2d(ii,jj+1) + yv_2d(ii+1,jj+1))/4.d0
             zc_3d(ii,jj,kk) = (zv_2d(ii,jj) + zv_2d(ii+1,jj) + zv_2d(ii,jj+1) + zv_2d(ii+1,jj+1))/4.d0 -(dz/2.0d0 + (nz-kk)*dz)
          enddo
       enddo
    enddo

    count = 0
    do kk = 1,nz
       do jj = 1,ny
          do ii = 1,nx
             count = count + 1
             id_3d(ii,jj,kk) = count
             soil_xc(count) = xc_3d(ii,jj,kk)
             soil_yc(count) = yc_3d(ii,jj,kk)
             soil_zc(count) = zc_3d(ii,jj,kk)
          enddo
       enddo
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

    iconn = 0

    !
    ! Vertical connections
    !
    do kk = 1, nz-1
       do ii = 1,nx
          do jj = 1,ny

             iconn = iconn + 1
             conn_id_up(iconn)   = id_3d(ii,jj,kk  )
             conn_id_dn(iconn)   = id_3d(ii,jj,kk+1)

             dist_x = soil_xc(conn_id_up(iconn)) - soil_xc(conn_id_dn(iconn))
             dist_y = soil_yc(conn_id_up(iconn)) - soil_yc(conn_id_dn(iconn))
             dist_z = soil_zc(conn_id_up(iconn)) - soil_zc(conn_id_dn(iconn))

             dist   = (dist_x**2.d0 + dist_y**2.d0 + dist_z**2.d0)**0.5d0

             conn_dist_up(iconn) = 0.5d0*dz !dist
             conn_dist_dn(iconn) = 0.5d0*dz !dist

             conn_area(iconn)    = soil_area(kk)
             conn_type(iconn)    = CONN_VERTICAL

          end do
       end do
    end do

    !
    ! Horizontal connections
    !
    do ii = 1,nx-1
       do kk = 1, nz
          do jj = 1,ny

             iconn = iconn + 1
             conn_id_up(iconn)   = id_3d(ii  ,jj,kk)
             conn_id_dn(iconn)   = id_3d(ii+1,jj,kk)

             dist_x = soil_xc(conn_id_up(iconn)) - soil_xc(conn_id_dn(iconn))
             dist_y = soil_yc(conn_id_up(iconn)) - soil_yc(conn_id_dn(iconn))
             dist_z = soil_zc(conn_id_up(iconn)) - soil_zc(conn_id_dn(iconn))

             dist   = (dist_x**2.d0 + dist_y**2.d0 + dist_z**2.d0)**0.5d0

             conn_dist_up(iconn) = 0.5d0*dist
             conn_dist_dn(iconn) = 0.5d0*dist

             conn_area(iconn)    = dz*dy
             conn_type(iconn)    = CONN_HORIZONTAL

          end do
       end do
    end do

    do jj = 1,ny-1
       do kk = 1, nz
          do ii = 1,nx

             iconn = iconn + 1
             conn_id_up(iconn)   = id_3d(ii,jj,kk)
             conn_id_dn(iconn)   = id_3d(ii,jj+1,kk)

             dist_x = soil_xc(conn_id_up(iconn)) - soil_xc(conn_id_dn(iconn))
             dist_y = soil_yc(conn_id_up(iconn)) - soil_yc(conn_id_dn(iconn))
             dist_z = soil_zc(conn_id_up(iconn)) - soil_zc(conn_id_dn(iconn))

             dist   = (dist_x**2.d0 + dist_y**2.d0 + dist_z**2.d0)**0.5d0

             conn_dist_up(iconn) = 0.5d0*dist
             conn_dist_dn(iconn) = 0.5d0*dist

             conn_area(iconn)    = dz*dx
             conn_type(iconn)    = CONN_HORIZONTAL

          end do
       end do
    end do

    nconn = iconn

    call vsfm_mpp%MeshSetConnectionSet(imesh, CONN_SET_INTERNAL, &
         nconn,  conn_id_up, conn_id_dn,                         &
         conn_dist_up, conn_dist_dn,  conn_area, conn_type)

    deallocate(soil_xc)
    deallocate(soil_yc)
    deallocate(soil_zc)
    deallocate(soil_dx)
    deallocate(soil_dy)
    deallocate(soil_dz)
    deallocate(soil_area)
    deallocate(soil_filter)

    deallocate(zv_x)
    deallocate(zv_y)

    deallocate(xv_2d)
    deallocate(yv_2d)
    deallocate(zv_2d)

    deallocate(xc_3d)
    deallocate(yc_3d)
    deallocate(zc_3d)
    deallocate(id_3d)

    deallocate (conn_id_up)
    deallocate (conn_id_dn)
    deallocate (conn_dist_up)
    deallocate (conn_dist_dn)
    deallocate (conn_area    )
    deallocate (conn_type    )

  end subroutine add_meshes

  !------------------------------------------------------------------------
  subroutine add_orthogonal_hex_mesh()
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
    use MultiPhysicsProbConstants , only : CONN_HORIZONTAL
    use mpp_varpar                , only : mpp_varpar_set_nlevsoi, mpp_varpar_set_nlevgrnd
    !
    implicit none
    !
#include <petsc/finclude/petsc.h>
    !
    PetscInt           :: imesh, ii, jj, kk
    PetscInt           :: nlev
    PetscInt           :: count
    PetscInt           :: iconn, nconn
    PetscReal, pointer :: soil_xc(:)           ! x-position of grid cell [m]
    PetscReal, pointer :: soil_yc(:)           ! y-position of grid cell [m]
    PetscReal, pointer :: soil_zc(:)           ! z-position of grid cell [m]
    PetscReal, pointer :: soil_dx(:)           ! layer thickness of grid cell [m]
    PetscReal, pointer :: soil_dy(:)           ! layer thickness of grid cell [m]
    PetscReal, pointer :: soil_dz(:)           ! layer thickness of grid cell [m]
    PetscReal, pointer :: soil_area(:)         ! area of grid cell [m^2]
    PetscInt , pointer :: soil_filter(:)       !

    PetscInt, pointer  :: conn_id_up(:)   !
    PetscInt, pointer  :: conn_id_dn(:)   !
    PetscReal, pointer :: conn_dist_up(:) !
    PetscReal, pointer :: conn_dist_dn(:) !
    PetscReal, pointer :: conn_area(:)    !
    PetscInt , pointer :: conn_type(:)    !
    PetscReal, pointer :: zv_x(:)
    PetscReal, pointer :: zv_y(:)
    PetscReal, pointer :: xv_2d(:,:)
    PetscReal, pointer :: yv_2d(:,:)
    PetscReal, pointer :: zv_2d(:,:)
    PetscReal, pointer :: xc_3d(:,:,:)
    PetscReal, pointer :: yc_3d(:,:,:)
    PetscReal, pointer :: zc_3d(:,:,:)
    PetscInt, pointer  :: id_3d(:,:,:)
    PetscReal :: dist_x, dist_y, dist_z, dist
    PetscReal, pointer :: xc_hex_3d(:,:,:)
    PetscReal, pointer :: yc_hex_3d(:,:,:)
    PetscReal, pointer :: zc_hex_3d(:,:,:)
    PetscReal, pointer :: zc_hex(:)
    PetscInt, pointer  :: id_hex_3d(:,:,:)
    PetscInt, pointer  :: act_hex_3d(:,:,:)
    PetscInt :: kk_upper, idx, count_2d

    PetscErrorCode :: ierr

    call mpp_varpar_set_nlevsoi(nz_hex)
    call mpp_varpar_set_nlevgrnd(nz_hex)

    nz = 30

    allocate(zv_x(nx+1))
    allocate(zv_y(ny+1))

    allocate(xv_2d(nx+1,ny+1))
    allocate(yv_2d(nx+1,ny+1))
    allocate(zv_2d(nx+1,ny+1))

    allocate(xc_3d(nx,ny,nz))
    allocate(yc_3d(nx,ny,nz))
    allocate(zc_3d(nx,ny,nz))
    allocate(id_3d(nx,ny,nz))

    ! Compute elevation at vertices in x-direction
    do ii = 1, nx/2+1
       zv_x(ii) = slope_x*dx*(nx/2) - (ii-1)*slope_x*dx
    enddo
    do ii = nx/2+2, nx+1
       zv_x(ii) = (ii-nx/2-1)*slope_x*dx
    enddo

    ! Compute elevation at vertices in y-direction
    do jj = 1, ny+1
       zv_y(jj) = (jj-1)*slope_y*dy
    enddo

    do jj = 1,ny+1
       do ii = 1,nx+1
          zv_2d(ii,jj) = zv_x(ii) + zv_y(jj);
          xv_2d(ii,jj) = (ii-1)*dx;
          yv_2d(ii,jj) = (jj-1)*dy;
       enddo
    enddo

    do kk = 1,nz
       do jj = 1,ny
          do ii = 1,nx
             xc_3d(ii,jj,kk) = (xv_2d(ii,jj) + xv_2d(ii+1,jj) + xv_2d(ii,jj+1) + xv_2d(ii+1,jj+1))/4.d0
             yc_3d(ii,jj,kk) = (yv_2d(ii,jj) + yv_2d(ii+1,jj) + yv_2d(ii,jj+1) + yv_2d(ii+1,jj+1))/4.d0
             zc_3d(ii,jj,kk) = (zv_2d(ii,jj) + zv_2d(ii+1,jj) + zv_2d(ii,jj+1) + zv_2d(ii+1,jj+1))/4.d0 -(dz/2.0d0 + (nz-kk)*dz)
          enddo
       enddo
    enddo


       imesh        = 1
       nlev         = nz_hex
       ncells_local = nx*ny*nz_hex
       ncells_ghost = 0

       allocate(soil_xc     (ncells_local))
       allocate(soil_yc     (ncells_local))
       allocate(soil_zc     (ncells_local))
       allocate(soil_dx     (ncells_local))
       allocate(soil_dy     (ncells_local))
       allocate(soil_dz     (ncells_local))
       allocate(soil_area   (ncells_local))
       allocate(soil_filter (ncells_local))

       soil_filter (:) = 0
       soil_area   (:) = dx*dy
       soil_dx     (:) = dx
       soil_dy     (:) = dy
       soil_dz     (:) = dz
       soil_xc     (:) = dx/2.d0
       soil_yc     (:) = dy/2.d0

       allocate(zc_hex(nz_hex))

       do kk = 1, nz_hex
          zc_hex(kk) = -13.5d0 + dz/2.d0 + (kk-1)*dz
       enddo

       allocate(xc_hex_3d(nx,ny,nz_hex))
       allocate(yc_hex_3d(nx,ny,nz_hex))
       allocate(zc_hex_3d(nx,ny,nz_hex))
       allocate(id_hex_3d(nx,ny,nz_hex))
       allocate(act_hex_3d(nx,ny,nz_hex))

       act_hex_3d(:,:,:) = 0
       id_hex_3d (:,:,:) = 0

       count_2d = 0
       allocate(kk_lower_idx(nx*ny))
       do jj = 1,ny
          do ii = 1,nx
             kk_upper = -1
             do kk = 1,nz_hex
                xc_hex_3d(ii,jj,kk) = xc_3d(ii,jj,nz)
                yc_hex_3d(ii,jj,kk) = yc_3d(ii,jj,nz)
                zc_hex_3d(ii,jj,kk) = zc_hex(kk)

                if (zc_hex(kk) <= zc_3d(ii,jj,nz)) kk_upper = kk
             enddo
             count_2d = count_2d + 1
             kk_lower_idx(count_2d) = kk_upper-30+1
             do kk = kk_upper-30+1,kk_upper
                act_hex_3d(ii,jj,kk) = 1
             enddo
          enddo
       enddo

       allocate(soil_filter_2d(nx*ny,nz_hex))
       soil_filter_2d(:,:) = 0

       idx = 0
       count = 0
       do kk = 1,nz_hex
          count_2d = 0
          do jj = 1,ny
             do ii = 1,nx
                count    = count + 1
                count_2d = count_2d + 1

                if (act_hex_3d(ii,jj,kk) == 1) then
                   idx                         = idx + 1
                   id_hex_3d(ii,jj,kk)         = idx
                   soil_filter(count)          = 1
                   soil_filter_2d(count_2d,kk) = 1
                endif

                soil_xc(count) = xc_hex_3d(ii,jj,kk)
                soil_yc(count) = yc_hex_3d(ii,jj,kk)
                soil_zc(count) = zc_hex_3d(ii,jj,kk)
             enddo
          enddo
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

    !
    ! Vertical connections
    !
    nz          = nz_hex
    allocate (conn_id_up   ((nz-1)*nx*ny + nx*(ny-1)*nz + (nx-1)*ny*nz))
    allocate (conn_id_dn   ((nz-1)*nx*ny + nx*(ny-1)*nz + (nx-1)*ny*nz))
    allocate (conn_dist_up ((nz-1)*nx*ny + nx*(ny-1)*nz + (nx-1)*ny*nz))
    allocate (conn_dist_dn ((nz-1)*nx*ny + nx*(ny-1)*nz + (nx-1)*ny*nz))
    allocate (conn_area    ((nz-1)*nx*ny + nx*(ny-1)*nz + (nx-1)*ny*nz))
    allocate (conn_type    ((nz-1)*nx*ny + nx*(ny-1)*nz + (nx-1)*ny*nz))


    iconn = 0
    do ii = 1,nx
       do jj = 1,ny
          do kk = 1, nz_hex-1
             if (act_hex_3d(ii,jj,kk) == 0 .or. act_hex_3d(ii,jj,kk+1) == 0 ) cycle

             iconn = iconn + 1
             conn_id_up(iconn)   = id_hex_3d(ii,jj,kk  )
             conn_id_dn(iconn)   = id_hex_3d(ii,jj,kk+1)

             dist_x = soil_xc(conn_id_up(iconn)) - soil_xc(conn_id_dn(iconn))
             dist_y = soil_yc(conn_id_up(iconn)) - soil_yc(conn_id_dn(iconn))
             dist_z = soil_zc(conn_id_up(iconn)) - soil_zc(conn_id_dn(iconn))

             !write(*,*)'conn: ',iconn,conn_id_up(iconn),conn_id_dn(iconn),dist_x,dist_y,dist_z
             dist   = (dist_x**2.d0 + dist_y**2.d0 + dist_z**2.d0)**0.5d0

             conn_dist_up(iconn) = 0.5d0*dz !dist
             conn_dist_dn(iconn) = 0.5d0*dz !dist

             conn_area(iconn)    = soil_area(kk)
             conn_type(iconn)    = CONN_VERTICAL

          end do
       end do
    end do

    !
    ! Horizontal connections
    !
    do ii = 1,nx-1
       do kk = 1, nz_hex
          do jj = 1,ny

             if (act_hex_3d(ii,jj,kk) == 0 .or. act_hex_3d(ii+1,jj,kk) == 0 ) cycle

             iconn = iconn + 1
             conn_id_up(iconn)   = id_hex_3d(ii  ,jj,kk)
             conn_id_dn(iconn)   = id_hex_3d(ii+1,jj,kk)

             dist_x = soil_xc(conn_id_up(iconn)) - soil_xc(conn_id_dn(iconn))
             dist_y = soil_yc(conn_id_up(iconn)) - soil_yc(conn_id_dn(iconn))
             dist_z = soil_zc(conn_id_up(iconn)) - soil_zc(conn_id_dn(iconn))

             dist   = (dist_x**2.d0 + dist_y**2.d0 + dist_z**2.d0)**0.5d0

             conn_dist_up(iconn) = 0.5d0*dist
             conn_dist_dn(iconn) = 0.5d0*dist

             conn_area(iconn)    = dz*dy
             conn_type(iconn)    = CONN_HORIZONTAL

          end do
       end do
    end do

    do jj = 1,ny-1
       do kk = 1, nz_hex
          do ii = 1,nx

             if (act_hex_3d(ii,jj,kk) == 0 .or. act_hex_3d(ii,jj+1,kk) == 0 ) cycle

             iconn = iconn + 1
             conn_id_up(iconn)   = id_hex_3d(ii,jj,kk)
             conn_id_dn(iconn)   = id_hex_3d(ii,jj+1,kk)

             dist_x = soil_xc(conn_id_up(iconn)) - soil_xc(conn_id_dn(iconn))
             dist_y = soil_yc(conn_id_up(iconn)) - soil_yc(conn_id_dn(iconn))
             dist_z = soil_zc(conn_id_up(iconn)) - soil_zc(conn_id_dn(iconn))

             dist   = (dist_x**2.d0 + dist_y**2.d0 + dist_z**2.d0)**0.5d0

             conn_dist_up(iconn) = 0.5d0*dist
             conn_dist_dn(iconn) = 0.5d0*dist

             conn_area(iconn)    = dz*dx
             conn_type(iconn)    = CONN_HORIZONTAL

          end do
       end do
    end do
    nconn = iconn

    call vsfm_mpp%MeshSetConnectionSet(imesh, CONN_SET_INTERNAL, &
         nconn,  conn_id_up, conn_id_dn,                         &
         conn_dist_up, conn_dist_dn,  conn_area, conn_type)

    deallocate(soil_xc)
    deallocate(soil_yc)
    deallocate(soil_zc)
    deallocate(soil_dx)
    deallocate(soil_dy)
    deallocate(soil_dz)
    deallocate(soil_area)
    deallocate(soil_filter)

    deallocate(zv_x)
    deallocate(zv_y)

    deallocate(xv_2d)
    deallocate(yv_2d)
    deallocate(zv_2d)

    deallocate(xc_3d)
    deallocate(yc_3d)
    deallocate(zc_3d)
    deallocate(id_3d)

    deallocate (conn_id_up)
    deallocate (conn_id_dn)
    deallocate (conn_dist_up)
    deallocate (conn_dist_dn)
    deallocate (conn_area    )
    deallocate (conn_type    )

  end subroutine add_orthogonal_hex_mesh

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
    use MultiPhysicsProbConstants , only : SOIL_CELLS
    use MultiPhysicsProbConstants , only : SOIL_TOP_CELLS
    use MultiPhysicsProbConstants , only : SOIL_BOTTOM_CELLS
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_SEEPAGE_BC
    use MultiPhysicsProbConstants , only : CONN_VERTICAL
    use ConnectionSetType         , only : connection_set_type
    use ConnectionSetType         , only : ConnectionSetDestroy
    use MeshType                  , only : MeshCreateConnectionSet
    !
    ! !ARGUMENTS
    implicit none
    !
    PetscInt                            :: ieqn
    PetscInt                            :: icond
    PetscInt                            :: iconn, nconn
    PetscInt                            :: ii,jj
    PetscInt                  , pointer :: conn_id_up(:)   !
    PetscInt                  , pointer :: conn_id_dn(:)   !
    PetscReal                 , pointer :: conn_dist_up(:) !
    PetscReal                 , pointer :: conn_dist_dn(:) !
    PetscReal                 , pointer :: conn_area(:)    !
    PetscInt                  , pointer :: conn_type(:)    !
    PetscReal                 , pointer :: conn_unitvec(:,:)
    type(connection_set_type) , pointer :: conn_set

    if (.not. with_seepage_bc) return

    allocate (conn_id_up   (nx*ny))
    allocate (conn_id_dn   (nx*ny))
    allocate (conn_dist_up (nx*ny))
    allocate (conn_dist_dn (nx*ny))
    allocate (conn_area    (nx*ny))
    allocate (conn_type    (nx*ny))
    allocate (conn_unitvec (nx*ny, 3))

    nconn = nx*ny
    iconn = 0
    do jj = 1,ny
       do ii = 1,nx
          iconn = iconn + 1
          conn_id_dn   (iconn)   = nx*ny*(nz-1) + iconn
          conn_id_up   (iconn)   = -1
          conn_dist_dn (iconn)   = 0.5d0*dz
          conn_dist_up (iconn)   = 0.d0
          conn_area    (iconn)   = dx*dy
          conn_type    (iconn)   = CONN_VERTICAL
          conn_unitvec (iconn,1) = 0.d0
          conn_unitvec (iconn,2) = 0.d0
          conn_unitvec (iconn,3) = -1.d0
       enddo
    enddo

    allocate(conn_set)

    call MeshCreateConnectionSet(vsfm_mpp%meshes(1), &
         nconn, conn_id_up, conn_id_dn, &
         conn_dist_up, conn_dist_dn, conn_area, &
         conn_type, conn_unitvec, conn_set)

    ieqn = 1
    call vsfm_mpp%soe%AddConditionInGovEqn(ieqn, COND_BC,   &
         'Constant head condition at top', 'Pa', COND_SEEPAGE_BC, &
         SOIL_TOP_CELLS, conn_set)

    call ConnectionSetDestroy(conn_set)

    deallocate (conn_id_up   )
    deallocate (conn_id_dn   )
    deallocate (conn_dist_up )
    deallocate (conn_dist_dn )
    deallocate (conn_area    )
    deallocate (conn_type    )
    deallocate (conn_unitvec )

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
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSoils
    use MultiPhysicsProbConstants , only : GRAVITY_CONSTANT
    use EOSWaterMod               , only : DENSITY_TGDPB01
    use mpp_varcon                , only : denh2o
    use mpp_varcon                , only : grav
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
    !-----------------------------------------------------------------------

    begc = 1
    endc = nx*ny

    satfunc_type = 'van_genuchten'

    allocate(vsfm_watsat       (1,ncells_local))
    allocate(vsfm_hksat        (1,ncells_local))
    allocate(vsfm_bsw          (1,ncells_local))
    allocate(vsfm_sucsat       (1,ncells_local))
    allocate(vsfm_eff_porosity (1,ncells_local))
    allocate(vsfm_residual_sat (1,ncells_local))
    allocate(vsfm_filter       (ncells_local))

    ! Soil properties
    vsfm_filter(:)         = 1
    vsfm_watsat(:,:)       = porosity
    vsfm_hksat(:,:)        = perm /vish2o * (denh2o * grav) / 0.001d0
    vsfm_bsw(:,:)          = 1.d0/lambda
    vsfm_sucsat(:,:)       = 1.d0/(alpha*grav)
    vsfm_eff_porosity(:,:) = 0.d0
    vsfm_residual_sat(:,:) = 0.2772d0

    call VSFMMPPSetSoils(vsfm_mpp, begc, endc, &
         ncells_ghost, vsfm_filter, &
         vsfm_watsat, vsfm_hksat, vsfm_bsw, &
         vsfm_sucsat, vsfm_residual_sat, &
         satfunc_type, DENSITY_TGDPB01)

    deallocate(vsfm_watsat)
    deallocate(vsfm_hksat   )
    deallocate(vsfm_bsw   )
    deallocate(vsfm_sucsat)
    deallocate(vsfm_eff_porosity)
    deallocate(vsfm_filter      )
    deallocate(vsfm_residual_sat)

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
    PetscReal, pointer :: press_ic(:)
    PetscInt :: jj, c, icell

    allocate(press_ic(ncells_local))

    icell = 0
    do jj = 1,nz
       do c = 1,nx*ny
          icell = icell + 1
          if (soil_filter_2d(c,jj) == 1) then
             press_ic(icell) = (18.75d0 -0.5d0*(jj-kk_lower_idx(c)) - 2.d0)*997.18d0*9.8d0 + 101325.d0
          else
             press_ic(icell) = -1.0d10
          endif
       enddo
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
    !
    implicit none
    !
    PetscReal, pointer :: top_pressure_bc(:)
    PetscInt           :: soe_auxvar_id

    allocate(top_pressure_bc(nx*ny))

    top_pressure_bc(:) = 101325.d0

    soe_auxvar_id = 1
    call vsfm_mpp%soe%SetDataFromCLM(AUXVAR_BC,  &
         VAR_BC_SS_CONDITION, soe_auxvar_id, top_pressure_bc)

    deallocate(top_pressure_bc)

  end subroutine set_bondary_conditions

  !------------------------------------------------------------------------
  subroutine output_regression_vsfm_vchannel_problem(filename_base, num_cells)
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

  end subroutine output_regression_vsfm_vchannel_problem

end module vsfm_vchannel_problem
