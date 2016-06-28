
module MPPVSFMALM_Initialize

#ifdef USE_PETSC_LIB

  !
  ! !USES:
  use MultiPhysicsProbVSFM , only : vsfm_mpp
  use shr_kind_mod         , only : r8 => shr_kind_r8
  use shr_log_mod          , only : errMsg => shr_log_errMsg
  use abortutils           , only : endrun
  use mpp_varctl           , only : iulog
  use GridcellType         , only : grc
  use LandunitType         , only : lun
  use ColumnType           , only : col
  !
  implicit none
  !
#include "finclude/petscsys.h"
  !
  ! IDs to indentify the conditions for VSFM
  integer :: vsfm_cond_id_for_infil
  integer :: vsfm_cond_id_for_et
  integer :: vsfm_cond_id_for_dew
  integer :: vsfm_cond_id_for_drainage
  integer :: vsfm_cond_id_for_snow
  integer :: vsfm_cond_id_for_sublimation
  integer :: vsfm_cond_id_for_lateral_flux
  !
  public :: MPPVSFMALM_Init

  !------------------------------------------------------------------------
contains

  subroutine MPPVSFMALM_Init
    !
    ! !DESCRIPTION:
    ! Initialization PETSc-based thermal model
    !
    ! !USES:
    use clm_varctl                , only : use_vsfm
    !
    ! !ARGUMENTS
    implicit none

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

    ! 9. Determine IDs for various source-sink condition
    call determine_condition_ids()

  end subroutine MPPVSFMALM_Init

  !------------------------------------------------------------------------
  subroutine initialize_mpp()
    !
    ! !DESCRIPTION:
    ! Initialization VSFM
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : MPP_VSFM_SNES_CLM
    use spmdMod                   , only : iam
    !
    ! !ARGUMENTS
    implicit none

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
    ! !DESCRIPTION:
    ! Add meshes used in the VSFM MPP
    !
    ! !USES:
    use decompMod                 , only : get_proc_clumps
    use decompMod                 , only : bounds_type, get_proc_bounds 
    use landunit_varcon           , only : istcrop, istsoil
    use column_varcon             , only : icol_road_perv
    use clm_varpar                , only : nlevgrnd
    use clm_varctl                , only : lateral_connectivity
    use spmdMod                   , only : iam
    use landunit_varcon           , only : max_lunit
    use clm_varctl                , only : vsfm_lateral_model_type

    use MultiPhysicsProbConstants , only : MESH_ALONG_GRAVITY
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
    use MultiPhysicsProbConstants , only : GE_RE
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_SS
    use MultiPhysicsProbConstants , only : COND_HEAT_FLUX
    use MultiPhysicsProbConstants , only : COND_HEAT_RATE
    use MultiPhysicsProbConstants , only : SNOW_TOP_CELLS
    use MultiPhysicsProbConstants , only : SNOW_BOTTOM_CELLS
    use MultiPhysicsProbConstants , only : SSW_TOP_CELLS
    use MultiPhysicsProbConstants , only : SOIL_TOP_CELLS
    use MultiPhysicsProbConstants , only : ALL_CELLS
    use MultiPhysicsProbConstants , only : CONN_VERTICAL
    use MultiPhysicsProbConstants , only : VAR_FRAC
    use MultiPhysicsProbConstants , only : VAR_ACTIVE
    use MultiPhysicsProbConstants , only : VAR_DIST_UP
    use MultiPhysicsProbConstants , only : VAR_DIST_DN
    use MultiPhysicsProbConstants , only : DISCRETIZATION_VERTICAL_ONLY
    use MultiPhysicsProbConstants , only : DISCRETIZATION_VERTICAL_WITH_SS
    use MultiPhysicsProbConstants , only : DISCRETIZATION_THREE_DIM
    !
    implicit none
    !
    !
    ! !LOCAL VARIABLES:
    integer            :: c,g,fc,j,l           ! do loop indices
    integer            :: nclumps              ! number of clumps on this processor

    integer            :: imesh
    integer            :: nlev
    integer            :: first_active_soil_col_id
    integer            :: col_id
    integer            :: icell
    integer            :: iconn
    integer            :: horz_nconn
    integer            :: vert_nconn
    integer            :: comb_nconn
    integer            :: ieqn
    integer            :: ieqn_1
    integer            :: ieqn_2
    integer            :: icond
    integer            :: ncells_local

    type(bounds_type)  :: bounds_proc

    real(r8), pointer  :: z(:,:)               ! centroid at "z" level [m]
    real(r8), pointer  :: zi(:,:)              ! interface level below a "z" level (m)
    real(r8), pointer  :: dz(:,:)              ! layer thickness at "z" level (m)

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
    PetscReal, pointer :: vert_conn_type(:)    !

    PetscInt, pointer  :: horz_conn_id_up(:)   !
    PetscInt, pointer  :: horz_conn_id_dn(:)   !
    PetscReal, pointer :: horz_conn_dist_up(:) !
    PetscReal, pointer :: horz_conn_dist_dn(:) !
    PetscReal, pointer :: horz_conn_area(:)    !
    PetscInt, pointer  :: horz_conn_type(:)    !

    PetscInt, pointer  :: comb_conn_id_up(:)   !
    PetscInt, pointer  :: comb_conn_id_dn(:)   !
    PetscReal, pointer :: comb_conn_dist_up(:) !
    PetscReal, pointer :: comb_conn_dist_dn(:) !
    PetscReal, pointer :: comb_conn_area(:)    !
    PetscInt, pointer  :: comb_conn_type(:)    !

    real(r8), pointer  :: xc_col(:)            ! x-position of grid cell [m]
    real(r8), pointer  :: yc_col(:)            ! y-position of grid cell [m]
    real(r8), pointer  :: zc_col(:)            ! z-position of grid cell [m]
    real(r8), pointer  :: area_col(:)          ! area of grid cell [m^2]
    integer, pointer   :: grid_owner(:)        ! MPI rank owner of grid cell

    integer            :: ncells_ghost         ! total number of ghost gridcells on the processor
    integer            :: ncols_ghost

    integer, pointer   :: vsfm_landunit_ind(:,:)
    integer, pointer   :: vsfm_coli(:)
    integer, pointer   :: vsfm_colf(:)

    PetscInt           :: discretization_type  !

    !----------------------------------------------------------------------

    z          =>    col%z
    zi         =>    col%zi
    dz         =>    col%dz

    call get_proc_bounds(bounds_proc)

    nclumps = get_proc_clumps()

    if (nclumps /= 1) then
       call endrun(msg='ERROR VSFM only supported for clumps = 1')
    endif
    
    allocate(xc_col            (bounds_proc%begc_all:bounds_proc%endc_all                    ))
    allocate(yc_col            (bounds_proc%begc_all:bounds_proc%endc_all                    ))
    allocate(zc_col            (bounds_proc%begc_all:bounds_proc%endc_all                    ))
    allocate(area_col          (bounds_proc%begc_all:bounds_proc%endc_all                    ))
    allocate(grid_owner        (bounds_proc%begg_all:bounds_proc%endg_all                    ))

    allocate (soil_xc           ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*nlevgrnd     ))
    allocate (soil_yc           ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*nlevgrnd     ))
    allocate (soil_zc           ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*nlevgrnd     ))
    allocate (soil_dx           ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*nlevgrnd     ))
    allocate (soil_dy           ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*nlevgrnd     ))
    allocate (soil_dz           ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*nlevgrnd     ))
    allocate (soil_area         ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*nlevgrnd     ))
    allocate (soil_filter       ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*nlevgrnd     ))

    allocate (vert_conn_id_up   ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*(nlevgrnd-1) ))
    allocate (vert_conn_id_dn   ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*(nlevgrnd-1) ))
    allocate (vert_conn_dist_up ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*(nlevgrnd-1) ))
    allocate (vert_conn_dist_dn ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*(nlevgrnd-1) ))
    allocate (vert_conn_area    ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*(nlevgrnd-1) ))
    allocate (vert_conn_type    ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*(nlevgrnd-1) ))

    xc_col(:)     = 0.d0
    yc_col(:)     = 0.d0
    zc_col(:)     = 0.d0
    area_col(:)   = 0.d0
    grid_owner(:) = 0

    first_active_soil_col_id = -1
    do c = bounds_proc%begc, bounds_proc%endc
       l = col%landunit(c)

       if (col%active(c) .and. .not.lun%lakpoi(l) .and. .not.lun%urbpoi(l)) then
          first_active_soil_col_id = c
          exit
       endif
    end do

    if (first_active_soil_col_id == -1) then
       write(iulog,*)'No active soil column found'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif    

    xc_col(:)     = 0._r8
    yc_col(:)     = 0._r8
    zc_col(:)     = 0._r8
    grid_owner(:) = 0
    area_col(:)   = 1._r8
    ncols_ghost   = 0
    horz_nconn    = 0
    
    if (lateral_connectivity) then

       call update_mesh_information (ncells_ghost, &
            ncols_ghost, xc_col, yc_col, zc_col,   &
            area_col, grid_owner )

       call setup_lateral_connections (grid_owner, &
            horz_nconn,                            &
            horz_conn_id_up, horz_conn_id_dn,      &
            horz_conn_dist_up, horz_conn_dist_dn,  &
            horz_conn_area, horz_conn_type)

    endif

    ! Save geometric attributes for soil mesh
    icell = 0
    do c = bounds_proc%begc, bounds_proc%endc
       do j = 1, nlevgrnd
          icell = icell + 1

          if (col%active(c) .and. &
               (lun%itype(l) == istsoil .or. col%itype(c) == icol_road_perv .or. &
                lun%itype(l) == istcrop)) then
             col_id = c
             soil_filter(icell) = 1
          else
             col_id = first_active_soil_col_id
             soil_filter(icell) = 0
          end if

          soil_xc(icell)   = xc_col(c)
          soil_yc(icell)   = yc_col(c)
          soil_zc(icell)   = -0.5d0*(zi(col_id,j-1) + zi(col_id,j)) + zc_col(c)
          soil_dx(icell)   = 1.d0
          soil_dy(icell)   = 1.d0
          soil_dz(icell)   = dz(col_id,j)
          soil_area(icell) = area_col(c)

       end do
    end do

    ! Save information about internal connections for soil mesh
    iconn = 0
    do c = bounds_proc%begc, bounds_proc%endc

       do j = 1, nlevgrnd-1

          iconn = iconn + 1
          vert_conn_id_up(iconn)   = (c-bounds_proc%begc)*nlevgrnd + j
          vert_conn_id_dn(iconn)   = vert_conn_id_up(iconn) + 1
          vert_conn_dist_up(iconn) = 0.5d0*dz(c,j  ) !zi(c,j)   - z(c,j)
          vert_conn_dist_dn(iconn) = 0.5d0*dz(c,j+1) !z( c,j+1) - zi(c,j)
          vert_conn_area(iconn)    = area_col(c)
          vert_conn_type(iconn)    = CONN_VERTICAL

       end do
    end do
    vert_nconn = iconn

    !
    ! Set up the meshes
    !    
    call vsfm_mpp%SetNumMeshes(1)

    !
    ! Set mesh value
    !
    imesh        = 1
    nlev         = nlevgrnd
    ncells_local = (bounds_proc%endc - bounds_proc%begc + 1)*nlev
    ncells_ghost = ncols_ghost*nlev

    call vsfm_mpp%MeshSetName        (imesh, 'Soil mesh')
    call vsfm_mpp%MeshSetOrientation (imesh, MESH_ALONG_GRAVITY)
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
    ! Set connection
    !
    if (vsfm_lateral_model_type == 'none') then
       discretization_type = DISCRETIZATION_VERTICAL_ONLY

       call vsfm_mpp%MeshSetConnectionSet(imesh, CONN_SET_INTERNAL, &
            vert_nconn,  vert_conn_id_up, vert_conn_id_dn, &
            vert_conn_dist_up, vert_conn_dist_dn,  vert_conn_area)

    else if (vsfm_lateral_model_type == 'source_sink') then
       discretization_type = DISCRETIZATION_VERTICAL_WITH_SS

       call vsfm_mpp%MeshSetConnectionSet(imesh, CONN_SET_INTERNAL, &
            vert_nconn,  vert_conn_id_up, vert_conn_id_dn, &
            vert_conn_dist_up, vert_conn_dist_dn, vert_conn_area)

       call vsfm_mpp%MeshSetConnectionSet(imesh, CONN_SET_LATERAL, &
            horz_nconn,  horz_conn_id_up, horz_conn_id_dn, &
            horz_conn_dist_up, horz_conn_dist_dn, horz_conn_area)

    else if (vsfm_lateral_model_type == 'three_dimensional') then
       discretization_type = DISCRETIZATION_THREE_DIM

       comb_nconn = vert_nconn + horz_nconn

       allocate (comb_conn_id_up   (comb_nconn))
       allocate (comb_conn_id_dn   (comb_nconn))
       allocate (comb_conn_dist_up (comb_nconn))
       allocate (comb_conn_dist_dn (comb_nconn))
       allocate (comb_conn_area    (comb_nconn))
       allocate (comb_conn_type    (comb_nconn))

       comb_conn_id_up   (1:vert_nconn) = vert_conn_id_up   (1:vert_nconn)
       comb_conn_id_dn   (1:vert_nconn) = vert_conn_id_dn   (1:vert_nconn)
       comb_conn_dist_up (1:vert_nconn) = vert_conn_dist_up (1:vert_nconn)
       comb_conn_dist_dn (1:vert_nconn) = vert_conn_dist_dn (1:vert_nconn)
       comb_conn_area    (1:vert_nconn) = vert_conn_area    (1:vert_nconn)
       comb_conn_type    (1:vert_nconn) = vert_conn_type    (1:vert_nconn)

       comb_conn_id_up   (vert_nconn+1:comb_nconn) = horz_conn_id_up   (1:horz_nconn)
       comb_conn_id_dn   (vert_nconn+1:comb_nconn) = horz_conn_id_dn   (1:horz_nconn)
       comb_conn_dist_up (vert_nconn+1:comb_nconn) = horz_conn_dist_up (1:horz_nconn)
       comb_conn_dist_dn (vert_nconn+1:comb_nconn) = horz_conn_dist_dn (1:horz_nconn)
       comb_conn_area    (vert_nconn+1:comb_nconn) = horz_conn_area    (1:horz_nconn)
       comb_conn_type    (vert_nconn+1:comb_nconn) = horz_conn_type    (1:horz_nconn)

       call vsfm_mpp%MeshSetConnectionSet(imesh, CONN_SET_INTERNAL, &
            comb_nconn,  comb_conn_id_up, comb_conn_id_dn, &
            comb_conn_dist_up, comb_conn_dist_dn, comb_conn_area)

       deallocate (comb_conn_id_up  )
       deallocate (comb_conn_id_dn  )
       deallocate (comb_conn_dist_up)
       deallocate (comb_conn_dist_dn)
       deallocate (comb_conn_area   )
       deallocate (comb_conn_type   )

    else
       call endrun(msg='ERROR: ' // &
            'Unknown vsfm_lateral_model_type = ' // trim(vsfm_lateral_model_type) // &
            errMsg(__FILE__, __LINE__))
    endif

    ! Free up memory
    deallocate(xc_col            )
    deallocate(yc_col            )
    deallocate(zc_col            )
    deallocate(area_col          )
    deallocate(grid_owner        )

    deallocate (soil_xc          )
    deallocate (soil_yc          )
    deallocate (soil_zc          )
    deallocate (soil_dx          )
    deallocate (soil_dy          )
    deallocate (soil_dz          )
    deallocate (soil_area        )
    deallocate (soil_filter      )

    deallocate (vert_conn_id_up  )
    deallocate (vert_conn_id_dn  )
    deallocate (vert_conn_dist_up)
    deallocate (vert_conn_dist_dn)
    deallocate (vert_conn_area   )
    deallocate (vert_conn_type   )

    if (lateral_connectivity) then
       deallocate (horz_conn_id_up   )
       deallocate (horz_conn_id_dn   )
       deallocate (horz_conn_dist_up )
       deallocate (horz_conn_dist_dn )
       deallocate (horz_conn_area    )
       deallocate (horz_conn_type    )
    endif

  end subroutine add_meshes

  !------------------------------------------------------------------------

  subroutine update_mesh_information (ncells_ghost, ncols_ghost,  &
       xc_col, yc_col, zc_col, area_col, &
       grid_owner &
    )
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use domainMod                 , only : ldomain
    use domainLateralMod          , only : ldomain_lateral
    use clm_varpar                , only : nlevgrnd
    use landunit_varcon           , only : istsoil
    use initGridCellsMod          , only : initGhostGridCells
    use decompMod                 , only : get_proc_total_ghosts
    use decompMod                 , only : bounds_type, get_proc_bounds 
    use spmdMod                   , only : iam
    use UnstructuredGridType      , only : ScatterDataG2L
    use clm_instMod               , only : soilstate_vars
    !
    implicit none
    !
    integer            :: ncols_ghost    ! number of ghost columns
    real(r8), pointer  :: xc_col(:)      ! x-position of grid cell [m]
    real(r8), pointer  :: yc_col(:)      ! y-position of grid cell [m]
    real(r8), pointer  :: zc_col(:)      ! z-position of grid cell [m]
    real(r8), pointer  :: area_col(:)    ! area of grid cell [m^2]
    integer, pointer   :: grid_owner(:)  ! MPI rank owner of grid cell
                                         !
    integer            :: c,g,fc,j,l     ! do loop indices
    type(bounds_type)  :: bounds_proc
    integer            :: nblocks
    integer            :: beg_idx
    integer            :: ndata_send     ! number of data sent by local mpi rank
    integer            :: ndata_recv     ! number of data received by local mpi rank
    real(r8), pointer  :: data_send(:)   ! data sent by local mpi rank
    real(r8), pointer  :: data_recv(:)   ! data received by local mpi rank

    integer            :: ncells_ghost   ! total number of ghost gridcells on the processor
    integer            :: nlunits_ghost  ! total number of ghost landunits on the processor
    integer            :: npfts_ghost    ! total number of ghost pfts on the processor
    integer            :: nCohorts_ghost ! total number of ghost cohorts on the processor
    !----------------------------------------------------------------------

    call get_proc_bounds(bounds_proc)

    call initGhostGridCells()

    call soilstate_vars%InitColdGhost(bounds_proc)

    nblocks    = 4
    ndata_send = nblocks*ldomain_lateral%ugrid%ngrid_local
    ndata_recv = nblocks*ldomain_lateral%ugrid%ngrid_ghosted

    allocate(data_send(ndata_send))
    allocate(data_recv(ndata_recv))

    ! Aggregate the data to send
    do g = bounds_proc%begg, bounds_proc%endg

       beg_idx = (g-bounds_proc%begg)*nblocks

       if (isnan(ldomain%xCell(g))) then
          call endrun(msg='ERROR initialize3: xCell = NaN')
       endif
       beg_idx = beg_idx + 1;
       data_send(beg_idx) = ldomain%xCell(g)

       if (isnan(ldomain%yCell(g))) then
          call endrun(msg='ERROR initialize3: yCell = NaN')
       endif
       beg_idx = beg_idx + 1;
       data_send(beg_idx) = ldomain%yCell(g)

       if (isnan(ldomain%topo(g))) then
          call endrun(msg='ERROR initialize3: topo = NaN')
       endif
       beg_idx = beg_idx + 1;
       data_send(beg_idx) = ldomain%topo(g)

       beg_idx = beg_idx + 1;
       data_send(beg_idx) = real(iam)

    enddo

    ! Scatter: Global-to-Local
    call ScatterDataG2L(ldomain_lateral%ugrid, nblocks, &
         ndata_send, data_send, ndata_recv, data_recv)

    ! Save data for ghost subgrid category
    do c = bounds_proc%begc_all, bounds_proc%endc_all

       g       = col%gridcell(c)
       beg_idx = (g-bounds_proc%begg)*nblocks

       beg_idx = beg_idx + 1; xc_col(c) = data_recv(beg_idx)
       beg_idx = beg_idx + 1; yc_col(c) = data_recv(beg_idx)
       beg_idx = beg_idx + 1; zc_col(c) = data_recv(beg_idx)
       beg_idx = beg_idx + 1; grid_owner(g) = data_recv(beg_idx)

       area_col(c) = ldomain_lateral%ugrid%areaGrid_ghosted(g-bounds_proc%begg + 1)

    enddo

    deallocate(data_send)
    deallocate(data_recv)

    call get_proc_total_ghosts(ncells_ghost, nlunits_ghost, &
         ncols_ghost, npfts_ghost, nCohorts_ghost)

  end subroutine update_mesh_information

  !------------------------------------------------------------------------

  subroutine setup_lateral_connections (grid_owner, &
       nconn_horz,                                  &
       horz_conn_id_up, horz_conn_id_dn,            &
       horz_conn_dist_up, horz_conn_dist_dn,        &
       horz_conn_area, horz_conn_type               &
       )
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    !
    use domainLateralMod          , only : ldomain_lateral
    use clm_varpar                , only : nlevgrnd
    use MultiPhysicsProbConstants , only : CONN_HORIZONTAL
    use landunit_varcon           , only : istsoil
    !
    implicit none
    !
    integer, pointer   :: grid_owner(:)        ! MPI rank owner of grid cell
    PetscInt           :: nconn_horz
    PetscInt, pointer  :: horz_conn_type(:)    !
    PetscInt, pointer  :: horz_conn_id_up(:)   !
    PetscInt, pointer  :: horz_conn_id_dn(:)   !
    PetscReal, pointer :: horz_conn_dist_up(:) !
    PetscReal, pointer :: horz_conn_dist_dn(:) !
    PetscReal, pointer :: horz_conn_area(:)    !
    !
    ! !LOCAL VARIABLES:
    PetscInt           :: c,j,l                !indices
    integer            :: begg,endg
    integer            :: begc,endc
    PetscInt           :: icell
    PetscInt           :: iconn
    PetscInt           :: nconn
    PetscInt           :: id_up, id_dn
    PetscInt           :: first_active_hydro_col_id
    PetscInt           :: col_id
    PetscInt           :: g_up, g_dn, iedge
    PetscInt           :: c_idx_up, c_idx_dn
    PetscInt           :: l_idx_up, l_idx_dn
    PetscInt           :: ltype, ctype
    PetscInt           :: tmp
    PetscReal          :: dist_x, dist_y, dist_z, dist
    PetscReal          :: dc, dv
    PetscErrorCode     :: ierr
    real(r8), pointer  :: dz(:,:)              ! layer thickness at "z" level (m)
    !----------------------------------------------------------------------

    dz         =>    col%dz

    !
    ! Sets up lateral connection between columns of type 'istsoil'
    !
    ! Assumptions:
    ! - There is only ONE 'istsoil' column per landunit per grid cell.
    ! - Grid cells that are laterally connected via cellsOnCell
    !   field defined in domain netcdf file have at least ONE
    !   column of 'istsoil' type
    
    nconn_horz = 0

    ltype = istsoil
    ctype = istsoil

    ! Determine number of lateral connections

    do icell = 1, ldomain_lateral%ugrid%ngrid_local
       do iedge = 1, ldomain_lateral%ugrid%maxEdges

          if (ldomain_lateral%ugrid%gridsOnGrid_local(iedge,icell) > icell) then
             g_up = icell + begg - 1
             g_dn = ldomain_lateral%ugrid%gridsOnGrid_local(iedge,icell) + begg - 1

             l_idx_up = grc%landunit_indices(ltype, g_up)
             l_idx_up = grc%landunit_indices(ltype, g_dn)

             c_idx_up = -1
             c_idx_dn = -1

             do c = lun%coli(l_idx_up), lun%colf(l_idx_up)
                if (col%itype(c) == ctype) then
                   if (c_idx_up /= -1) then
                      write(iulog,*)'CreateFromCLMCols: More than one column found for ' // &
                           'ctype = ', ctype, ' for ltype = ', ltype, ' in grid cell ', g_up
                      call endrun(msg=errMsg(__FILE__, __LINE__))
                   endif
                   c_idx_up = c
                endif
             enddo

             do c = lun%coli(l_idx_dn), lun%colf(l_idx_dn)
                if (col%itype(c) == ctype) then
                   if (c_idx_dn /= -1) then
                      write(iulog,*)'CreateFromCLMCols: More than one column found for ' // &
                           'ctype = ', ctype, ' for ltype = ', ltype, ' in grid cell ', g_dn
                      call endrun(msg=errMsg(__FILE__, __LINE__))
                   endif
                   c_idx_dn = c
                endif
             enddo

             if (c_idx_up > -1 .and. c_idx_dn > -1) then
                nconn_horz = nconn_horz + 1
             else
                write(iulog,*)'CreateFromCLMCols: No column of ctype = ', ctype, &
                     ' found between following grid cells: ',g_up,g_dn
                call endrun(msg=errMsg(__FILE__, __LINE__))
             endif
          endif
       enddo
    enddo

    nconn_horz = nconn_horz * nlevgrnd

    allocate (horz_conn_id_up   (nconn_horz))
    allocate (horz_conn_id_dn   (nconn_horz))
    allocate (horz_conn_dist_up (nconn_horz))
    allocate (horz_conn_dist_dn (nconn_horz))
    allocate (horz_conn_area    (nconn_horz))
    allocate (horz_conn_type    (nconn_horz))

    do icell = 1, ldomain_lateral%ugrid%ngrid_local

       do iedge = 1, ldomain_lateral%ugrid%maxEdges
          if (ldomain_lateral%ugrid%gridsOnGrid_local(iedge,icell) > icell) then
             g_up = icell + begg - 1
             g_dn = ldomain_lateral%ugrid%gridsOnGrid_local(iedge,icell) + begg - 1

             l_idx_up = grc%landunit_indices(ltype, g_up)
             l_idx_up = grc%landunit_indices(ltype, g_dn)

             c_idx_up = -1
             c_idx_dn = -1

             do c = lun%coli(l_idx_up), lun%colf(l_idx_up)
                if (col%itype(c) == ctype) then
                   if (c_idx_up /= -1) then
                      write(iulog,*)'CreateFromCLMCols: More than one column found for ' // &
                           'ctype = ', ctype, ' for ltype = ', ltype, ' in grid cell ', g_up
                      call endrun(msg=errMsg(__FILE__, __LINE__))
                   endif
                   c_idx_up = c
                endif
             enddo

             do c = lun%coli(l_idx_dn), lun%colf(l_idx_dn)
                if (col%itype(c) == ctype) then
                   if (c_idx_dn /= -1) then
                      write(iulog,*)'CreateFromCLMCols: More than one column found for ' // &
                           'ctype = ', ctype, ' for ltype = ', ltype, ' in grid cell ', g_dn
                      call endrun(msg=errMsg(__FILE__, __LINE__))
                   endif
                   c_idx_dn = c
                endif
             enddo

             if (c_idx_up > -1 .and. c_idx_dn > -1) then

                if (grid_owner(g_up) > grid_owner(g_dn)) then
                   tmp      = g_up;
                   g_up     = g_dn
                   g_dn     = tmp

                   tmp      = l_idx_up;
                   l_idx_up = l_idx_dn
                   l_idx_dn = tmp

                   tmp      = c_idx_up;
                   c_idx_up = c_idx_dn
                   c_idx_dn = tmp
                endif

                dc = ldomain_lateral%ugrid%dcOnGrid_local(iedge, icell)
                dv = ldomain_lateral%ugrid%dvOnGrid_local(iedge, icell)

                do j = 1, nlevgrnd
                   iconn = iconn + 1

                   id_up = (c_idx_up - begc)*nlevgrnd + j
                   id_dn = (c_idx_dn - begc)*nlevgrnd + j

                   horz_conn_type         = CONN_HORIZONTAL
                   horz_conn_id_up(iconn) = id_up
                   horz_conn_id_dn(iconn) = id_dn

                   !this%is_active(id_up) = PETSC_TRUE
                   !this%is_active(id_dn) = PETSC_TRUE

                   horz_conn_area(iconn) = dz(c_idx_up,j)*dv

                   horz_conn_dist_up(iconn) = 0.5d0*dist
                   horz_conn_dist_dn(iconn) = 0.5d0*dist

                enddo
             endif ! if (c_idx_up > -1 .and. c_idx_dn > -1)
          endif ! if (ugrid%gridsOnGrid_local(iedge,icell) > icell)
       enddo
    enddo

  end subroutine setup_lateral_connections

  !------------------------------------------------------------------------

  subroutine add_goveqns()
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : GE_RE
    !
    ! !ARGUMENTS
    implicit none

    call vsfm_mpp%AddGovEqn(GE_RE, 'Richards Equation ODE')

    call vsfm_mpp%SetMeshesOfGoveqns()
    
  end subroutine add_goveqns

  !------------------------------------------------------------------------
  subroutine add_conditions_to_goveqns()
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use clm_varctl                , only : vsfm_lateral_model_type
    use MultiPhysicsProbConstants , only : SOIL_CELLS
    use MultiPhysicsProbConstants , only : SOIL_TOP_CELLS
    use MultiPhysicsProbConstants , only : COND_SS
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_MASS_RATE
    use MultiPhysicsProbConstants , only : COND_SEEPAGE_BC
    !
    ! !ARGUMENTS
    implicit none
    !
    PetscInt :: ieqn

    ieqn = 1

    call vsfm_mpp%GovEqnAddCondition(ieqn, COND_SS,   &
         'Infiltration_Flux', 'kg/s', COND_MASS_RATE, &
         SOIL_TOP_CELLS)

    call vsfm_mpp%GovEqnAddCondition(ieqn, COND_SS,   &
         'Evapotranspiration_Flux', 'kg/s', COND_MASS_RATE, &
         SOIL_CELLS)

    call vsfm_mpp%GovEqnAddCondition(ieqn, COND_SS,   &
         'Dew_Flux', 'kg/s', COND_MASS_RATE, &
         SOIL_TOP_CELLS)

    call vsfm_mpp%GovEqnAddCondition(ieqn, COND_SS,   &
         'Drainage_Flux', 'kg/s', COND_MASS_RATE, &
         SOIL_CELLS)

    call vsfm_mpp%GovEqnAddCondition(ieqn, COND_SS,   &
         'Snow_Disappearance_Flux', 'kg/s', COND_MASS_RATE, &
         SOIL_TOP_CELLS)

    call vsfm_mpp%GovEqnAddCondition(ieqn, COND_SS,   &
         'Sublimation_Flux', 'kg/s', COND_MASS_RATE, &
         SOIL_TOP_CELLS)

    if (vsfm_lateral_model_type == 'source_sink' .or. &
        vsfm_lateral_model_type == 'three_dimensional') then

       call vsfm_mpp%GovEqnAddCondition(ieqn, COND_SS,   &
            'Lateral_Flux', 'kg/s', COND_MASS_RATE, &
            SOIL_CELLS)

       call vsfm_mpp%GovEqnAddCondition(ieqn, COND_BC,   &
            'Seepage_bc', 'kg/s', COND_SEEPAGE_BC, &
            SOIL_TOP_CELLS)

    endif

  end subroutine add_conditions_to_goveqns

  !------------------------------------------------------------------------
  subroutine allocate_auxvars()
    !
    ! !DESCRIPTION:
    !
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
    use decompMod                 , only : bounds_type
    use decompMod                 , only : get_proc_bounds 
    use clm_instMod               , only : soilstate_vars
    use clm_instMod               , only : soilhydrology_vars
    use decompMod                 , only : get_proc_total_ghosts
    use landunit_varcon           , only : istcrop
    use landunit_varcon           , only : istsoil
    use column_varcon             , only : icol_road_perv
    use clm_varpar                , only : nlevgrnd
    use clm_varctl                , only : vsfm_satfunc_type
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSoils
    use EOSWaterMod               , only : DENSITY_TGDPB01
    !
    implicit none
    !
    real(r8), pointer    :: clm_watsat(:,:)
    real(r8), pointer    :: clm_hksat(:,:)
    real(r8), pointer    :: clm_bsw(:,:)
    real(r8), pointer    :: clm_sucsat(:,:)
    real(r8), pointer    :: clm_eff_porosity(:,:)
    real(r8), pointer    :: clm_zwt(:)
    real(r8), pointer    :: vsfm_watsat(:,:)
    real(r8), pointer    :: vsfm_hksat(:,:)
    real(r8), pointer    :: vsfm_bsw(:,:)
    real(r8), pointer    :: vsfm_sucsat(:,:)
    real(r8), pointer    :: vsfm_eff_porosity(:,:)
    real(r8), pointer    :: vsfm_residual_sat(:,:)
    integer, pointer     :: vsfm_filter(:)
    !
    type(bounds_type)    :: bounds_proc
    !
    integer              :: c,g,fc,j,l            ! do loop indices
    integer              :: ncells_ghost          ! total number of ghost gridcells on the processor
    integer              :: nlunits_ghost         ! total number of ghost landunits on the processor
    integer              :: ncols_ghost           ! total number of ghost columns on the processor
    integer              :: npfts_ghost           ! total number of ghost pfts on the processor
    integer              :: nCohorts_ghost        ! total number of ghost cohorts on the processor
    !-----------------------------------------------------------------------

    clm_watsat       => soilstate_vars%watsat_col       ! volumetric soil water at saturation (porosity)
    clm_hksat        => soilstate_vars%hksat_col        ! hydraulic conductivity at saturation (mm H2O /s)
    clm_bsw          => soilstate_vars%bsw_col          ! Clapp and Hornberger "b"
    clm_sucsat       => soilstate_vars%sucsat_col       ! minimum soil suction (mm)
    clm_eff_porosity => soilstate_vars%eff_porosity_col ! effective porosity = porosity - vol_ice
    clm_zwt          => soilhydrology_vars%zwt_col      ! water table depth

    call get_proc_bounds(bounds_proc)

    call get_proc_total_ghosts(ncells_ghost, nlunits_ghost, &
         ncols_ghost, npfts_ghost, nCohorts_ghost)

    ! Allocate memory
    allocate(vsfm_filter       (bounds_proc%begc_all:bounds_proc%endc_all           ))
    allocate(vsfm_watsat       (bounds_proc%begc_all:bounds_proc%endc_all, nlevgrnd ))
    allocate(vsfm_hksat        (bounds_proc%begc_all:bounds_proc%endc_all, nlevgrnd ))
    allocate(vsfm_bsw          (bounds_proc%begc_all:bounds_proc%endc_all, nlevgrnd ))
    allocate(vsfm_sucsat       (bounds_proc%begc_all:bounds_proc%endc_all, nlevgrnd ))
    allocate(vsfm_eff_porosity (bounds_proc%begc_all:bounds_proc%endc_all, nlevgrnd ))
    allocate(vsfm_residual_sat (bounds_proc%begc_all:bounds_proc%endc_all, nlevgrnd ))

    ! Initialize
    vsfm_filter       (:)   = 0
    vsfm_watsat       (:,:) = 0._r8
    vsfm_hksat        (:,:) = 0._r8
    vsfm_bsw          (:,:) = 0._r8
    vsfm_sucsat       (:,:) = 0._r8
    vsfm_residual_sat (:,:) = 0._r8

    ! Save data to initialize VSFM
    do c = bounds_proc%begc, bounds_proc%endc
       l = col%landunit(c)

       if (col%active(c) .and. &
           (lun%itype(l) == istsoil .or. col%itype(c) == icol_road_perv .or. &
           lun%itype(l) == istcrop)) then

          vsfm_filter    (c) = 1

          do j = 1 ,nlevgrnd
             vsfm_watsat(c,j)       = clm_watsat(c,j)
             vsfm_hksat(c,j)        = clm_hksat(c,j)
             vsfm_bsw(c,j)          = clm_bsw(c,j)
             vsfm_sucsat(c,j)       = clm_sucsat(c,j)
             vsfm_eff_porosity(c,j) = clm_eff_porosity(c,j)
          enddo

       endif
    enddo

    call VSFMMPPSetSoils(vsfm_mpp, bounds_proc%begc, bounds_proc%endc, &
         ncols_ghost, vsfm_filter, &
         vsfm_watsat, vsfm_hksat, vsfm_bsw, vsfm_sucsat, vsfm_eff_porosity, &
         vsfm_residual_sat, vsfm_satfunc_type, DENSITY_TGDPB01)

    ! Free up memory
    deallocate(vsfm_filter       )
    deallocate(vsfm_watsat       )
    deallocate(vsfm_hksat        )
    deallocate(vsfm_bsw          )
    deallocate(vsfm_sucsat       )
    deallocate(vsfm_eff_porosity )
    deallocate(vsfm_residual_sat )

  end subroutine set_material_properties

  !------------------------------------------------------------------------
  subroutine set_initial_conditions()
    !
    ! !DESCRIPTION:
    !
    use decompMod                 , only : bounds_type
    use decompMod                 , only : get_proc_bounds 
    use clm_instMod               , only : soilstate_vars
    use clm_instMod               , only : soilhydrology_vars
    use decompMod                 , only : get_proc_total_ghosts
    use landunit_varcon           , only : istcrop
    use landunit_varcon           , only : istsoil
    use column_varcon             , only : icol_road_perv
    use clm_varpar                , only : nlevgrnd
    use clm_varctl                , only : vsfm_satfunc_type
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSoils
    use MultiPhysicsProbConstants , only : GRAVITY_CONSTANT
    use MultiPhysicsProbConstants , only : PRESSURE_REF
    !
    implicit none
    !
    real(r8), pointer    :: clm_zi(:,:)           ! interface level below a "z" level (m)
    real(r8), pointer    :: clm_zwt(:)            ! 
    !
    type(bounds_type)    :: bounds_proc
    !
    integer              :: c,g,fc,j,l            ! do loop indices
    integer              :: ncells_ghost          ! total number of ghost gridcells on the processor
    integer              :: nlunits_ghost         ! total number of ghost landunits on the processor
    integer              :: ncols_ghost           ! total number of ghost columns on the processor
    integer              :: npfts_ghost           ! total number of ghost pfts on the processor
    integer              :: nCohorts_ghost        ! total number of ghost cohorts on the processor
    real(r8), pointer    :: press_ic_1d(:)        ! pressure initial condition (m)
    integer :: icell
    !-----------------------------------------------------------------------

    clm_zwt          => soilhydrology_vars%zwt_col      ! water table depth
    clm_zi           => col%zi

    call get_proc_bounds(bounds_proc)

    ! Allocate memory
    allocate(press_ic_1d ((bounds_proc%endc_all - bounds_proc%begc_all + 1)*nlevgrnd))

    ! Initialize
    press_ic_1d(:) = 101325.d0

    ! Save data to initialize VSFM
    do c = bounds_proc%begc, bounds_proc%endc
       l = col%landunit(c)
       
       do j = 1, nlevgrnd
          icell = (c - bounds_proc%begc)*nlevgrnd + j
          if (col%active(c) .and. &
               (lun%itype(l) == istsoil .or. col%itype(c) == icol_road_perv .or. &
               lun%itype(l) == istcrop)) then

             press_ic_1d(icell) = PRESSURE_REF + &
                  997.16d0*GRAVITY_CONSTANT * &
                  (-clm_zwt(c) - (-0.5d0*(clm_zi(c,j-1) + clm_zi(c,j))))
          endif
       enddo
    enddo

    !call VSFMMPPSetICs(vsfm_mpp, bounds_proc%begc, bounds_proc%endc, &
    !zc_col, vsfm_filter, vsfm_zwt)
    call vsfm_mpp%Restart(press_ic_1d)

    ! Free up memory
    deallocate(press_ic_1d)

  end subroutine set_initial_conditions

  !-----------------------------------------------------------------------
  subroutine determine_condition_ids()
    !
    !DESCRIPTION
    !  Determines the IDs of various source-sink conditions in VSFM
    !
    use clm_varctl                       , only : use_vsfm, lateral_connectivity
    use MultiPhysicsProbVSFM             , only : vsfm_mpp
    use MultiPhysicsProbConstants        , only : COND_SS
    use MultiPhysicsProbConstants        , only : COND_NULL
    use clm_varctl                       , only : iulog
    use abortutils                       , only : endrun
    use shr_log_mod                      , only : errMsg => shr_log_errMsg
    ! !ARGUMENTS:
    implicit none
    integer :: ier ! error status
    character (len=256), pointer :: cond_names(:)
    integer                      :: num_conds
    integer                      :: num_conds_expected
    integer                      :: nn
    integer                      :: kk
    character (len=256)          :: cond_name
    !------------------------------------------------------------------------------


    num_conds_expected = 6

    if (lateral_connectivity) num_conds_expected = num_conds_expected + 1

    ! Get the number of conditions
    call vsfm_mpp%sysofeqns%GetConditionNames(COND_SS, COND_NULL, num_conds, cond_names)

    if (num_conds /= num_conds_expected) then
      write(iulog,*)'In init_vsfm_condition_ids: Source-sink conditions /= ', num_conds_expected
      call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    do nn = 1, num_conds
       select case(trim(cond_names(nn)))
       case ("Infiltration_Flux")
         vsfm_cond_id_for_infil = nn
       case ("Evapotranspiration_Flux")
         vsfm_cond_id_for_et = nn
       case ("Dew_Flux")
         vsfm_cond_id_for_dew = nn
       case ("Drainage_Flux")
         vsfm_cond_id_for_drainage = nn
       case ("Snow_Disappearance_Flux")
         vsfm_cond_id_for_snow = nn
       case ("Sublimation_Flux")
         vsfm_cond_id_for_sublimation = nn
       case ("Lateral_flux")
         vsfm_cond_id_for_lateral_flux = nn
       case default
         write(iulog,*)'In init_vsfm_condition_ids: Unknown flux.'
         call endrun(msg=errMsg(__FILE__, __LINE__))
       end select
    enddo

    if (vsfm_cond_id_for_infil == -1) then
      write(iulog,*)'In init_vsfm_condition_ids: vsfm_cond_id_for_infil not defined.'
      call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    if (vsfm_cond_id_for_et == -1) then
      write(iulog,*)'In init_vsfm_condition_ids: vsfm_cond_id_for_et not defined.'
      call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    if (vsfm_cond_id_for_dew == -1) then
      write(iulog,*)'In init_vsfm_condition_ids: vsfm_cond_id_for_dew not defined.'
      call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    if (vsfm_cond_id_for_drainage == -1) then
      write(iulog,*)'In init_vsfm_condition_ids: vsfm_cond_id_for_drainage not defined.'
      call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    if (vsfm_cond_id_for_snow == -1) then
      write(iulog,*)'In init_vsfm_condition_ids: vsfm_cond_id_for_snow not defined.'
      call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    if (vsfm_cond_id_for_sublimation == -1) then
      write(iulog,*)'In init_vsfm_condition_ids: vsfm_cond_id_for_sublimation not defined.'
      call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    if (lateral_connectivity .and. vsfm_cond_id_for_lateral_flux == -1) then
      write(iulog,*)'In init_vsfm_condition_ids: vsfm_cond_id_for_lateral_flux not defined.'
      call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    deallocate(cond_names)

  end subroutine determine_condition_ids

#endif

end module MPPVSFMALM_Initialize
