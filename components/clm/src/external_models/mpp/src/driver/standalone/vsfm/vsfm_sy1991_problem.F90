module vsfm_sy1991_problem

  implicit none

#include <petsc/finclude/petsc.h>
  PetscInt  , parameter :: nx       = 1
  PetscInt  , parameter :: ny       = 1
  PetscReal , parameter :: x_column = 1.d0
  PetscReal , parameter :: y_column = 1.d0
  PetscReal , parameter :: z_column = 2.d0
  PetscInt              :: nz
  PetscInt              :: problem_number
  PetscInt  , parameter :: WETTING_PROBLEM = 1
  PetscInt  , parameter :: DRYING_PROBLEM  = 2
  PetscInt              :: ncells_local
  PetscInt              :: ncells_ghost
  PetscReal , parameter :: press_ic_wetting(200) = (/ &
       101281.100000d0, 101193.700000d0, 101106.800000d0, 101020.600000d0, 100935.100000d0,  &
       100850.500000d0, 100766.800000d0, 100684.100000d0, 100602.600000d0, 100522.300000d0,  &
       100443.400000d0, 100366.100000d0, 100290.400000d0, 100216.500000d0, 100144.500000d0,  &
       100074.600000d0, 100006.800000d0, 99941.280000d0, 99878.170000d0, 99817.550000d0,  &
       99759.500000d0, 99704.080000d0, 99651.350000d0, 99601.330000d0, 99554.040000d0,  &
       99509.460000d0, 99467.580000d0, 99428.350000d0, 99391.710000d0, 99357.590000d0,  &
       99325.900000d0, 99296.550000d0, 99269.430000d0, 99244.430000d0, 99221.440000d0,  &
       99200.350000d0, 99181.020000d0, 99163.360000d0, 99147.230000d0, 99132.530000d0,  &
       99119.160000d0, 99107.010000d0, 99095.980000d0, 99085.980000d0, 99076.920000d0,  &
       99068.730000d0, 99061.320000d0, 99054.630000d0, 99048.590000d0, 99043.140000d0,  &
       99038.220000d0, 99033.790000d0, 99029.800000d0, 99026.210000d0, 99022.980000d0,  &
       99020.070000d0, 99017.450000d0, 99015.100000d0, 99012.990000d0, 99011.090000d0,  &
       99009.380000d0, 99007.850000d0, 99006.470000d0, 99005.230000d0, 99004.120000d0,  &
       99003.130000d0, 99002.230000d0, 99001.430000d0, 99000.700000d0, 99000.060000d0,  &
       98999.480000d0, 98998.950000d0, 98998.490000d0, 98998.070000d0, 98997.690000d0,  &
       98997.350000d0, 98997.050000d0, 98996.780000d0, 98996.530000d0, 98996.310000d0,  &
       98996.120000d0, 98995.940000d0, 98995.780000d0, 98995.640000d0, 98995.510000d0,  &
       98995.400000d0, 98995.300000d0, 98995.210000d0, 98995.120000d0, 98995.050000d0,  &
       98994.980000d0, 98994.920000d0, 98994.870000d0, 98994.820000d0, 98994.780000d0,  &
       98994.740000d0, 98994.710000d0, 98994.680000d0, 98994.650000d0, 98994.620000d0,  &
       98953.280000d0, 98866.850000d0, 98781.600000d0, 98697.630000d0, 98615.040000d0,  &
       98533.930000d0, 98454.380000d0, 98376.510000d0, 98300.410000d0, 98226.160000d0,  &
       98153.860000d0, 98083.590000d0, 98015.410000d0, 97949.410000d0, 97885.640000d0,  &
       97824.140000d0, 97764.970000d0, 97708.140000d0, 97653.680000d0, 97601.580000d0,  &
       97551.860000d0, 97504.500000d0, 97459.470000d0, 97416.730000d0, 97376.250000d0,  &
       97337.980000d0, 97301.860000d0, 97267.820000d0, 97235.810000d0, 97205.730000d0,  &
       97177.520000d0, 97151.110000d0, 97126.400000d0, 97103.320000d0, 97081.780000d0,  &
       97061.720000d0, 97043.030000d0, 97025.660000d0, 97009.520000d0, 96994.540000d0,  &
       96980.640000d0, 96967.760000d0, 96955.830000d0, 96944.790000d0, 96934.580000d0,  &
       96925.140000d0, 96916.430000d0, 96908.370000d0, 96900.940000d0, 96894.090000d0,  &
       96887.770000d0, 96881.940000d0, 96876.570000d0, 96871.620000d0, 96867.070000d0,  &
       96862.870000d0, 96859.010000d0, 96855.460000d0, 96852.190000d0, 96849.180000d0,  &
       96846.420000d0, 96843.870000d0, 96841.530000d0, 96839.380000d0, 96837.410000d0,  &
       96835.590000d0, 96833.920000d0, 96832.390000d0, 96830.980000d0, 96829.690000d0,  &
       96828.500000d0, 96827.400000d0, 96826.400000d0, 96825.480000d0, 96824.630000d0,  &
       96823.860000d0, 96823.140000d0, 96822.490000d0, 96821.880000d0, 96821.330000d0,  &
       96820.830000d0, 96820.360000d0, 96819.930000d0, 96819.540000d0, 96819.180000d0,  &
       96818.850000d0, 96818.540000d0, 96818.260000d0, 96818.010000d0, 96817.770000d0,  &
       96817.560000d0, 96817.360000d0, 96817.180000d0, 96817.010000d0, 96816.860000d0,  &
       96816.720000d0, 96816.590000d0, 96816.470000d0, 96816.360000d0, 96816.260000d0  &
       /)

  PetscReal , parameter :: press_ic_drying(200) = (/ &
       101320.200000d0, 101310.800000d0, 101301.700000d0, 101292.900000d0, 101284.500000d0,  &
       101276.300000d0, 101268.500000d0, 101261.000000d0, 101253.700000d0, 101246.800000d0,  &
       101240.200000d0, 101233.900000d0, 101227.900000d0, 101222.100000d0, 101216.600000d0,  &
       101211.400000d0, 101206.400000d0, 101201.600000d0, 101197.100000d0, 101192.800000d0,  &
       101188.700000d0, 101184.900000d0, 101181.200000d0, 101177.700000d0, 101174.400000d0,  &
       101171.200000d0, 101168.200000d0, 101165.400000d0, 101162.700000d0, 101160.200000d0,  &
       101157.800000d0, 101155.500000d0, 101153.400000d0, 101151.300000d0, 101149.400000d0,  &
       101147.600000d0, 101145.800000d0, 101144.200000d0, 101142.700000d0, 101141.200000d0,  &
       101139.800000d0, 101138.500000d0, 101137.300000d0, 101136.100000d0, 101135.000000d0,  &
       101134.000000d0, 101133.000000d0, 101132.100000d0, 101131.200000d0, 101130.400000d0,  &
       101129.600000d0, 101128.900000d0, 101128.200000d0, 101127.600000d0, 101126.900000d0,  &
       101126.400000d0, 101125.800000d0, 101125.300000d0, 101124.800000d0, 101124.300000d0,  &
       101123.900000d0, 101123.500000d0, 101123.100000d0, 101122.700000d0, 101122.400000d0,  &
       101122.100000d0, 101121.800000d0, 101121.500000d0, 101121.200000d0, 101120.900000d0,  &
       101120.700000d0, 101120.500000d0, 101120.300000d0, 101120.100000d0, 101119.900000d0,  &
       101119.700000d0, 101119.500000d0, 101119.400000d0, 101119.200000d0, 101119.100000d0,  &
       101118.900000d0, 101118.800000d0, 101118.700000d0, 101118.600000d0, 101118.500000d0,  &
       101118.400000d0, 101118.300000d0, 101118.200000d0, 101118.100000d0, 101118.000000d0,  &
       101117.900000d0, 101117.900000d0, 101117.800000d0, 101117.700000d0, 101117.700000d0,  &
       101117.600000d0, 101117.600000d0, 101117.500000d0, 101117.500000d0, 101117.400000d0,  &
       101074.900000d0, 100987.800000d0, 100901.300000d0, 100815.700000d0, 100731.000000d0,  &
       100647.200000d0, 100564.600000d0, 100483.200000d0, 100403.200000d0, 100324.600000d0,  &
       100247.600000d0, 100172.300000d0, 100098.900000d0, 100027.500000d0, 99958.250000d0,  &
       99891.210000d0, 99826.530000d0, 99764.310000d0, 99704.630000d0, 99647.570000d0,  &
       99593.180000d0, 99541.510000d0, 99492.580000d0, 99446.380000d0, 99402.900000d0,  &
       99362.110000d0, 99323.950000d0, 99288.350000d0, 99255.240000d0, 99224.530000d0,  &
       99196.110000d0, 99169.870000d0, 99145.710000d0, 99123.510000d0, 99103.150000d0,  &
       99084.500000d0, 99067.470000d0, 99051.930000d0, 99037.770000d0, 99024.890000d0,  &
       99013.190000d0, 99002.570000d0, 98992.940000d0, 98984.220000d0, 98976.330000d0,  &
       98969.200000d0, 98962.750000d0, 98956.930000d0, 98951.680000d0, 98946.950000d0,  &
       98942.680000d0, 98938.840000d0, 98935.370000d0, 98932.260000d0, 98929.450000d0,  &
       98926.930000d0, 98924.660000d0, 98922.610000d0, 98920.780000d0, 98919.130000d0,  &
       98917.650000d0, 98916.310000d0, 98915.120000d0, 98914.040000d0, 98913.070000d0,  &
       98912.210000d0, 98911.430000d0, 98910.730000d0, 98910.100000d0, 98909.540000d0,  &
       98909.030000d0, 98908.580000d0, 98908.170000d0, 98907.800000d0, 98907.470000d0,  &
       98907.180000d0, 98906.910000d0, 98906.670000d0, 98906.460000d0, 98906.270000d0,  &
       98906.100000d0, 98905.940000d0, 98905.810000d0, 98905.680000d0, 98905.570000d0,  &
       98905.470000d0, 98905.380000d0, 98905.300000d0, 98905.230000d0, 98905.160000d0,  &
       98905.100000d0, 98905.050000d0, 98905.010000d0, 98904.960000d0, 98904.930000d0,  &
       98904.890000d0, 98904.860000d0, 98904.830000d0, 98904.810000d0, 98904.790000d0   &
       /)
  
  public :: run_vsfm_sy1991_problem

contains
  
!------------------------------------------------------------------------

  subroutine run_vsfm_sy1991_problem()
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
    !
    !
    PetscBool          :: converged
    PetscInt           :: converged_reason
    PetscErrorCode     :: ierr
    PetscReal          :: dtime
    PetscInt           :: istep, nstep
    PetscBool           :: flg
    PetscBool          :: save_initial_soln, save_final_soln
    character(len=256) :: string
    character(len=256) :: output_suffix
    PetscViewer        :: viewer

    ! Set default settings
    nz                = 200
    dtime             = 3600.d0
    nstep             = 24
    save_initial_soln = PETSC_FALSE
    save_final_soln   = PETSC_FALSE
    output_suffix     = ''
    problem_number    = DRYING_PROBLEM

    ! Get some command line options

    call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-problem_number',problem_number,flg,ierr)
    call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-dt',dtime,flg,ierr)
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-nstep',nstep,flg,ierr)
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-save_initial_soln',save_initial_soln,flg,ierr)
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-save_final_soln',save_final_soln,flg,ierr)
    call PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-output_suffix',output_suffix,flg,ierr)

    select case(problem_number)
    case (WETTING_PROBLEM)
    case (DRYING_PROBLEM)
    case default
       write(*,*)'Invalid problem no. Valid problem number are: 1 or 2'
       stop
    end select

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
       call set_bondary_conditions()

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

  end subroutine run_vsfm_sy1991_problem

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
    soil_area   (:) = dx*dy
    soil_dx     (:) = dx
    soil_dy     (:) = dy
    soil_dz     (:) = dz
    soil_xc     (:) = dx/2.d0
    soil_yc     (:) = dy/2.d0

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
    use MultiPhysicsProbVSFM , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : SOIL_CELLS
    use MultiPhysicsProbConstants , only : SOIL_TOP_CELLS
    use MultiPhysicsProbConstants , only : SOIL_BOTTOM_CELLS
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_SS
    use MultiPhysicsProbConstants , only : COND_DIRICHLET
    use MultiPhysicsProbConstants , only : COND_MASS_RATE
    !
    ! !ARGUMENTS
    implicit none
    !
    PetscInt :: ieqn

    ieqn = 1

    call vsfm_mpp%soe%AddConditionInGovEqn(ieqn, COND_SS,   &
         'Constant flux condition at top', 'kg/s/m^2', COND_MASS_RATE, &
         SOIL_TOP_CELLS)

    call vsfm_mpp%soe%AddConditionInGovEqn(ieqn, COND_BC,   &
         'Constant head condition at bottom', 'Pa', COND_DIRICHLET, &
         SOIL_BOTTOM_CELLS)

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
    PetscReal        , parameter :: porosity            = 0.4d0
    PetscReal        , parameter :: lambda              = 0.5455d0
    PetscReal        , parameter :: alpha               = 4.d-4
    PetscReal        , parameter :: perm_high           = 2.5281d-12
    PetscReal        , parameter :: perm_low            = 2.5281d-13
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
    endc = 1

    satfunc_type = 'van_genuchten'

    allocate(vsfm_watsat       (1,nz))
    allocate(vsfm_hksat        (1,nz))
    allocate(vsfm_bsw          (1,nz))
    allocate(vsfm_sucsat       (1,nz))
    allocate(vsfm_eff_porosity (1,nz))
    allocate(vsfm_residual_sat (1,nz))
    allocate(vsfm_filter       (nz))

    ! Soil properties
    vsfm_filter(:)         = 1
    vsfm_watsat(:,:)       = porosity
    vsfm_hksat(:,001:100)  = perm_low  /vish2o * (denh2o * grav) / 0.001d0
    vsfm_hksat(:,101:200)  = perm_high /vish2o * (denh2o * grav) / 0.001d0
    vsfm_bsw(:,:)          = 1.d0/lambda
    vsfm_sucsat(:,:)       = 1.d0/(alpha*GRAVITY_CONSTANT)
    vsfm_eff_porosity(:,:) = 0.d0
    vsfm_residual_sat(:,:) = 0.15d0

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

    select case(problem_number)
    case (DRYING_PROBLEM)
       call vsfm_mpp%Restart(press_ic_drying)
    case (WETTING_PROBLEM)
       call vsfm_mpp%Restart(press_ic_wetting)
    end select


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
    PetscReal, pointer :: top_recharge(:)
    PetscReal, pointer :: bot_pressure_bc(:)
    PetscInt           :: soe_auxvar_id
    PetscReal, parameter :: recharge_wetting = 2.5d-6
    PetscReal, parameter :: recharge_drying  = 2.7778d-7

    allocate(top_recharge(1))
    allocate(bot_pressure_bc(1))

    select case(problem_number)
    case (DRYING_PROBLEM)
       top_recharge(:)    = recharge_drying * 997.16d0
       bot_pressure_bc(:) = press_ic_drying(1)

    case (WETTING_PROBLEM)
       top_recharge(:)    = recharge_wetting * 997.16d0
       bot_pressure_bc(:) = press_ic_wetting(1)
    end select

    soe_auxvar_id = 1
    call vsfm_mpp%soe%SetDataFromCLM(AUXVAR_SS,  &
         VAR_BC_SS_CONDITION, soe_auxvar_id, top_recharge)

    soe_auxvar_id = 1
    call vsfm_mpp%soe%SetDataFromCLM(AUXVAR_BC,  &
         VAR_BC_SS_CONDITION, soe_auxvar_id, bot_pressure_bc)

    deallocate(top_recharge)
    deallocate(bot_pressure_bc)

  end subroutine set_bondary_conditions

  !------------------------------------------------------------------------
  subroutine output_regression_vsfm_celia1990_problem(filename_base, num_cells)
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

  end subroutine output_regression_vsfm_celia1990_problem

end module vsfm_sy1991_problem
