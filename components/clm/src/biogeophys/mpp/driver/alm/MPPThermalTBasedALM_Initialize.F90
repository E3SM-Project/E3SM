
module MPPThermalTBasedALM_Initialize

#ifdef USE_PETSC_LIB

  !-----------------------------------------------------------------------
  ! Performs land model initialization
  !
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use spmdMod          , only : masterproc, iam
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use decompMod        , only : bounds_type, get_proc_bounds 
  use abortutils       , only : endrun
  use mpp_varctl       , only : iulog
  ! 
  !-----------------------------------------
  ! Definition of component types
  !-----------------------------------------
  use LandunitType           , only : lun                
  use ColumnType             , only : col                
  use clm_instMod            , only : soilstate_vars

  implicit none
  
  public :: MPPThermalTBasedALM_Init

  !------------------------------------------------------------------------
contains

  subroutine MPPThermalTBasedALM_Init()
    !
    ! !DESCRIPTION:
    ! Initialization PETSc-based thermal model
    !
    ! !USES:
    use clm_varctl                , only : use_petsc_thermal_model
    use MultiPhysicsProbThermal   , only : thermal_mpp
    !
    ! !ARGUMENTS
    implicit none

    if (.not.use_petsc_thermal_model) return

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
    call thermal_mpp%SetupProblem()

    ! 7. Add material properities associated with all governing equations
    call add_material_properties()

  end subroutine MPPThermalTBasedALM_Init
    
  !------------------------------------------------------------------------
  subroutine initialize_mpp()
    !
    ! !DESCRIPTION:
    ! Initialization PETSc-based thermal model
    !
    ! !USES:
    use MultiPhysicsProbThermal   , only : thermal_mpp
    use MultiPhysicsProbConstants , only : MPP_THERMAL_TBASED_KSP_CLM
    use spmdMod                   , only : iam
    !
    ! !ARGUMENTS
    implicit none

    !
    ! Set up the multi-physics problem
    !
    call thermal_mpp%Init       ()
    call thermal_mpp%SetName    ('Snow + Standing water + Soil thermal model using temperature')
    call thermal_mpp%SetID      (MPP_THERMAL_TBASED_KSP_CLM)
    call thermal_mpp%SetMPIRank (iam)

  end subroutine initialize_mpp

  !------------------------------------------------------------------------
  subroutine add_meshes()
    !
    ! !DESCRIPTION:
    ! Add meshes to the thermal MPP problem
    !
    ! !USES:
    use spmdMod                   , only : mpicom, npes
    use filterMod                 , only : filter
    use decompMod                 , only : get_proc_clumps
    use clm_varpar                , only : nlevgrnd, nlevsno
    use clm_varctl                , only : finidat
    use shr_infnan_mod            , only : shr_infnan_isnan
    use abortutils                , only : endrun
    use shr_infnan_mod            , only : isnan => shr_infnan_isnan
    use clm_varcon                , only : spval
    use MultiPhysicsProbThermal   , only : thermal_mpp, MPPThermalSetSoils
    use MultiPhysicsProbConstants , only : MPP_THERMAL_TBASED_KSP_CLM
    use MultiPhysicsProbConstants , only : MESH_ALONG_GRAVITY
    use MultiPhysicsProbConstants , only : MESH_CLM_SNOW_COL
    use MultiPhysicsProbConstants , only : MESH_CLM_SSW_COL
    use MultiPhysicsProbConstants , only : MESH_CLM_THERMAL_SOIL_COL
    use MultiPhysicsProbConstants , only : VAR_XC
    use MultiPhysicsProbConstants , only : VAR_YC
    use MultiPhysicsProbConstants , only : VAR_ZC
    use MultiPhysicsProbConstants , only : VAR_DX
    use MultiPhysicsProbConstants , only : VAR_DY
    use MultiPhysicsProbConstants , only : VAR_DZ
    use MultiPhysicsProbConstants , only : VAR_AREA
    use MultiPhysicsProbConstants , only : CONN_SET_INTERNAL
    use MultiPhysicsProbConstants , only : GE_THERM_SNOW_TBASED
    use MultiPhysicsProbConstants , only : GE_THERM_SSW_TBASED
    use MultiPhysicsProbConstants , only : GE_THERM_SOIL_TBASED
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_SS
    use MultiPhysicsProbConstants , only : COND_HEAT_FLUX
    use MultiPhysicsProbConstants , only : COND_HEAT_RATE
    use MultiPhysicsProbConstants , only : SNOW_TOP_CELLS
    use MultiPhysicsProbConstants , only : SNOW_BOTTOM_CELLS
    use MultiPhysicsProbConstants , only : SSW_TOP_CELLS
    use MultiPhysicsProbConstants , only : SOIL_TOP_CELLS
    use MultiPhysicsProbConstants , only : ALL_CELLS
    use MultiPhysicsProbConstants , only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants , only : VAR_THERMAL_COND
    use MultiPhysicsProbConstants , only : VAR_FRAC
    use MultiPhysicsProbConstants , only : VAR_ACTIVE
    use MultiPhysicsProbConstants , only : VAR_DIST_UP
    use MultiPhysicsProbConstants , only : VAR_DIST_DN
    !
    ! !ARGUMENTS
    implicit none
    !
#include "finclude/petscsys.h"
    !
    ! !LOCAL VARIABLES:
    integer           :: nclumps                  ! number of clumps on this processor
    integer           :: nc                       ! clump index
    integer           :: c,g,fc,j,l               ! do loop indices
    type(bounds_type) :: bounds_proc
    real(r8), pointer :: z(:,:)                   ! centroid at "z" level [m]
    real(r8), pointer :: zi(:,:)                  ! interface level below a "z" level (m)
    real(r8), pointer :: dz(:,:)                  ! layer thickness at "z" level (m)

    integer           :: imesh
    integer           :: ncells_local
    integer           :: ncells_ghost
    integer           :: nlev
    integer           :: first_active_soil_col_id
    integer           :: col_id
    integer           :: icell
    integer           :: iconn
    integer           :: snow_nconn
    integer           :: soil_nconn
    integer           :: ieqn
    integer           :: ieqn_1
    integer           :: ieqn_2
    integer           :: icond
    integer           :: nvars_for_coupling
    integer, pointer  :: var_ids_for_coupling(:)
    integer, pointer  :: goveqn_ids_for_coupling(:)

    real(r8), pointer :: snow_xc(:)               ! x-position of grid cell [m]
    real(r8), pointer :: snow_yc(:)               ! y-position of grid cell [m]
    real(r8), pointer :: snow_zc(:)               ! z-position of grid cell [m]
    real(r8), pointer :: snow_dx(:)               ! layer thickness of grid cell [m]
    real(r8), pointer :: snow_dy(:)               ! layer thickness of grid cell [m]
    real(r8), pointer :: snow_dz(:)               ! layer thickness of grid cell [m]
    real(r8), pointer :: snow_area(:)             ! area of grid cell [m^2]
    integer , pointer :: snow_filter(:)           ! 
    integer , pointer :: snow_conn_id_up(:)       !
    integer , pointer :: snow_conn_id_dn(:)       !
    real(r8), pointer :: snow_conn_dist_up(:)     !
    real(r8), pointer :: snow_conn_dist_dn(:)     !
    real(r8), pointer :: snow_conn_area(:)        !

    real(r8), pointer :: ssw_xc(:)                ! x-position of grid cell [m]
    real(r8), pointer :: ssw_yc(:)                ! y-position of grid cell [m]
    real(r8), pointer :: ssw_zc(:)                ! z-position of grid cell [m]
    real(r8), pointer :: ssw_dx(:)                ! layer thickness of grid cell [m]
    real(r8), pointer :: ssw_dy(:)                ! layer thickness of grid cell [m]
    real(r8), pointer :: ssw_dz(:)                ! layer thickness of grid cell [m]
    real(r8), pointer :: ssw_area(:)              ! area of grid cell [m^2]
    integer , pointer :: ssw_filter(:)            ! 

    real(r8), pointer :: soil_xc(:)               ! x-position of grid cell [m]
    real(r8), pointer :: soil_yc(:)               ! y-position of grid cell [m]
    real(r8), pointer :: soil_zc(:)               ! z-position of grid cell [m]
    real(r8), pointer :: soil_dx(:)               ! layer thickness of grid cell [m]
    real(r8), pointer :: soil_dy(:)               ! layer thickness of grid cell [m]
    real(r8), pointer :: soil_dz(:)               ! layer thickness of grid cell [m]
    real(r8), pointer :: soil_area(:)             ! area of grid cell [m^2]
    integer , pointer :: soil_filter(:)           ! 
    integer, pointer  :: soil_conn_id_up(:)       !
    integer, pointer  :: soil_conn_id_dn(:)       !
    real(r8), pointer :: soil_conn_dist_up(:)     !
    real(r8), pointer :: soil_conn_dist_dn(:)     !
    real(r8), pointer :: soil_conn_area(:)        !

    !----------------------------------------------------------------------

    z          =>    col%z
    zi         =>    col%zi
    dz         =>    col%dz

    call get_proc_bounds(bounds_proc)
    nclumps = get_proc_clumps()

    if (nclumps /= 1) then
       call endrun(msg='ERROR clm_initializeMod: '//&
           'PETSc-based thermal model only supported for clumps = 1')
    endif

    nc = 1

    ! Allocate memory and setup data structure for VSFM-MPP

    allocate (snow_xc           ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*nlevsno      ))
    allocate (snow_yc           ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*nlevsno      ))
    allocate (snow_zc           ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*nlevsno      ))
    allocate (snow_dx           ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*nlevsno      ))
    allocate (snow_dy           ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*nlevsno      ))
    allocate (snow_dz           ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*nlevsno      ))
    allocate (snow_area         ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*nlevsno      ))
    allocate (snow_filter       ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*nlevsno      ))
    allocate (snow_conn_id_up   ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*(nlevsno-1)  ))
    allocate (snow_conn_id_dn   ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*(nlevsno-1)  ))
    allocate (snow_conn_dist_up ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*(nlevsno-1)  ))
    allocate (snow_conn_dist_dn ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*(nlevsno-1)  ))
    allocate (snow_conn_area    ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*(nlevsno-1)  ))

    allocate (ssw_xc            ((bounds_proc%endc_all-bounds_proc%begc_all+1                )))
    allocate (ssw_yc            ((bounds_proc%endc_all-bounds_proc%begc_all+1                )))
    allocate (ssw_zc            ((bounds_proc%endc_all-bounds_proc%begc_all+1                )))
    allocate (ssw_dx            ((bounds_proc%endc_all-bounds_proc%begc_all+1                )))
    allocate (ssw_dy            ((bounds_proc%endc_all-bounds_proc%begc_all+1                )))
    allocate (ssw_dz            ((bounds_proc%endc_all-bounds_proc%begc_all+1                )))
    allocate (ssw_area          ((bounds_proc%endc_all-bounds_proc%begc_all+1                )))
    allocate (ssw_filter        ((bounds_proc%endc_all-bounds_proc%begc_all+1                )))

    allocate (soil_xc           ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*nlevgrnd     ))
    allocate (soil_yc           ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*nlevgrnd     ))
    allocate (soil_zc           ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*nlevgrnd     ))
    allocate (soil_dx           ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*nlevgrnd     ))
    allocate (soil_dy           ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*nlevgrnd     ))
    allocate (soil_dz           ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*nlevgrnd     ))
    allocate (soil_area         ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*nlevgrnd     ))
    allocate (soil_filter       ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*nlevgrnd     ))
    allocate (soil_conn_id_up   ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*(nlevgrnd-1) ))
    allocate (soil_conn_id_dn   ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*(nlevgrnd-1) ))
    allocate (soil_conn_dist_up ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*(nlevgrnd-1) ))
    allocate (soil_conn_dist_dn ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*(nlevgrnd-1) ))
    allocate (soil_conn_area    ((bounds_proc%endc_all-bounds_proc%begc_all+1 )*(nlevgrnd-1) ))

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

    ! Save geometric attributes for snow mesh
    icell = 0
    do c = bounds_proc%begc, bounds_proc%endc
       do j = -nlevsno+1, 0
          icell = icell + 1

          if (col%active(c) .and. .not.lun%lakpoi(l) .and. .not.lun%urbpoi(l)) then
             col_id = c
             snow_filter(icell) = 1
          else
             col_id = first_active_soil_col_id
             snow_filter(icell) = 0
          end if

          snow_xc(icell)   = 0.d0
          snow_yc(icell)   = 0.d0
          snow_zc(icell)   = -0.5d0*(zi(col_id,j-1) + zi(col_id,j))
          snow_dx(icell)   = 1.d0
          snow_dy(icell)   = 1.d0
          snow_dz(icell)   = dz(col_id,j)
          snow_area(icell) = 1.d0

       end do
    end do

    ! Save geometric attributes for standing surface water mesh
    icell = 0
    j     = 0
    do c = bounds_proc%begc, bounds_proc%endc

       icell = icell + 1

       if (col%active(c) .and. .not.lun%lakpoi(l) .and. .not.lun%urbpoi(l)) then
          col_id = c
          ssw_filter(icell) = 1
       else
          col_id = first_active_soil_col_id
          ssw_filter(icell) = 0
       end if

       ssw_xc(icell)   = 0.d0
       ssw_yc(icell)   = 0.d0
       ssw_zc(icell)   = -0.5d0*(1.d-6 + zi(col_id,j))
       ssw_dx(icell)   = 1.d0
       ssw_dy(icell)   = 1.d0
       ssw_dz(icell)   = 1.d-6
       ssw_area(icell) = 1.d0

    end do

    ! Save geometric attributes for soil mesh
    icell = 0
    do c = bounds_proc%begc, bounds_proc%endc
       do j = 1, nlevgrnd
          icell = icell + 1

          if (col%active(c) .and. .not.lun%lakpoi(l) .and. .not.lun%urbpoi(l)) then
             col_id = c
             soil_filter(icell) = 1
          else
             col_id = first_active_soil_col_id
             soil_filter(icell) = 0
          end if

          soil_xc(icell)   = 0.d0
          soil_yc(icell)   = 0.d0
          soil_zc(icell)   = -0.5d0*(zi(col_id,j-1) + zi(col_id,j))
          soil_dx(icell)   = 1.d0
          soil_dy(icell)   = 1.d0
          soil_dz(icell)   = dz(col_id,j)
          soil_area(icell) = 1.d0

       end do
    end do

    ! Save information about internal connections for snow mesh
    iconn = 0
    do c = bounds_proc%begc, bounds_proc%endc
       do j = -nlevsno+1, -1

          iconn = iconn + 1
          snow_conn_id_up(iconn)   = (c-bounds_proc%begc)*nlevsno + j  + nlevsno
          snow_conn_id_dn(iconn)   = snow_conn_id_up(iconn) + 1
          snow_conn_dist_up(iconn) = zi(c,j)   - z(c,j)
          snow_conn_dist_dn(iconn) = z( c,j+1) - zi(c,j)
          snow_conn_area(iconn)    = 1.d0

       end do
    end do
    snow_nconn = iconn
    
    ! Save information about internal connections for soil mesh
    iconn = 0
    do c = bounds_proc%begc, bounds_proc%endc

       do j = 1, nlevgrnd-1

          iconn = iconn + 1
          soil_conn_id_up(iconn)   = (c-bounds_proc%begc)*nlevgrnd + j
          soil_conn_id_dn(iconn)   = soil_conn_id_up(iconn) + 1
          soil_conn_dist_up(iconn) = zi(c,j)   - z(c,j)
          soil_conn_dist_dn(iconn) = z( c,j+1) - zi(c,j)
          soil_conn_area(iconn)    = 1.d0

       end do
    end do
    soil_nconn = iconn

    !
    ! Set up the meshes
    !
    
    call thermal_mpp%SetNumMeshes (3)

    !
    ! Set mesh for snow
    !
    imesh        = 1
    nlev         = nlevsno
    ncells_local = (bounds_proc%endc - bounds_proc%begc + 1)*nlev
    ncells_ghost = 0

    call thermal_mpp%MeshSetName        (imesh, 'CLM snow thermal mesh')
    call thermal_mpp%MeshSetOrientation (imesh, MESH_ALONG_GRAVITY)
    call thermal_mpp%MeshSetID          (imesh, MESH_CLM_SNOW_COL)
    call thermal_mpp%MeshSetDimensions  (imesh, ncells_local, ncells_ghost, nlev)

    call thermal_mpp%MeshSetGridCellFilter      (imesh, snow_filter)
    call thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_XC   , snow_xc)
    call thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_YC   , snow_yc)
    call thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_ZC   , snow_zc)
    call thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_DX   , snow_dx)
    call thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_DY   , snow_dy)
    call thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_DZ   , snow_dz)
    call thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_AREA , snow_area)
    call thermal_mpp%MeshComputeVolume          (imesh)

    call thermal_mpp%MeshSetConnectionSet(imesh, CONN_SET_INTERNAL, snow_nconn,  &
         snow_conn_id_up, snow_conn_id_dn, snow_conn_dist_up, snow_conn_dist_dn, &
         snow_conn_area)

    !
    ! Set mesh for standing water
    !
    imesh        = 2
    nlev         = 1
    ncells_local = (bounds_proc%endc - bounds_proc%begc + 1)*nlev
    ncells_ghost = 0    
    call thermal_mpp%MeshSetName(imesh, 'CLM standing water thermal snow mesh')
    call thermal_mpp%MeshSetOrientation (imesh, MESH_ALONG_GRAVITY)
    call thermal_mpp%MeshSetID          (imesh, MESH_CLM_SSW_COL)
    call thermal_mpp%MeshSetDimensions  (imesh, ncells_local, ncells_ghost, nlev)

    call thermal_mpp%MeshSetGridCellFilter      (imesh, ssw_filter)
    call thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_XC   , ssw_xc)
    call thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_YC   , ssw_yc)
    call thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_ZC   , ssw_zc)
    call thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_DX   , ssw_dx)
    call thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_DY   , ssw_dy)
    call thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_DZ   , ssw_dz)
    call thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_AREA , ssw_area)
    call thermal_mpp%MeshComputeVolume          (imesh)

    !
    ! Set mesh for standing water
    !
    imesh        = 3
    nlev         = nlevgrnd
    ncells_local = (bounds_proc%endc - bounds_proc%begc + 1)*nlev
    ncells_ghost = 0    

    call thermal_mpp%MeshSetName        (imesh, 'CLM soil thermal mesh')
    call thermal_mpp%MeshSetOrientation (imesh, MESH_ALONG_GRAVITY)
    call thermal_mpp%MeshSetID          (imesh, MESH_CLM_THERMAL_SOIL_COL)
    call thermal_mpp%MeshSetDimensions  (imesh, ncells_local, ncells_ghost, nlev)

    call thermal_mpp%MeshSetGridCellFilter      (imesh, soil_filter)
    call thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_XC   , soil_xc)
    call thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_YC   , soil_yc)
    call thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_ZC   , soil_zc)
    call thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_DX   , soil_dx)
    call thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_DY   , soil_dy)
    call thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_DZ   , soil_dz)
    call thermal_mpp%MeshSetGeometricAttributes (imesh, VAR_AREA , soil_area)
    call thermal_mpp%MeshComputeVolume          (imesh)

    call thermal_mpp%MeshSetConnectionSet(imesh, CONN_SET_INTERNAL, soil_nconn,  &
         soil_conn_id_up, soil_conn_id_dn, soil_conn_dist_up, soil_conn_dist_dn, &
         soil_conn_area)

    ! Free up memory

    deallocate (snow_xc           )
    deallocate (snow_yc           )
    deallocate (snow_zc           )
    deallocate (snow_dz           )
    deallocate (snow_area         )
    deallocate (snow_filter       )
    deallocate (snow_conn_id_up   )
    deallocate (snow_conn_id_dn   )
    deallocate (snow_conn_dist_up )
    deallocate (snow_conn_dist_dn )
    deallocate (snow_conn_area    )

    deallocate (ssw_xc            )
    deallocate (ssw_yc            )
    deallocate (ssw_zc            )
    deallocate (ssw_dz            )
    deallocate (ssw_area          )
    deallocate (ssw_filter        )

    deallocate (soil_xc           )
    deallocate (soil_yc           )
    deallocate (soil_zc           )
    deallocate (soil_dz           )
    deallocate (soil_area         )
    deallocate (soil_filter       )
    deallocate (soil_conn_id_up   )
    deallocate (soil_conn_id_dn   )
    deallocate (soil_conn_dist_up )
    deallocate (soil_conn_dist_dn )
    deallocate (soil_conn_area    )

  end subroutine add_meshes

  !------------------------------------------------------------------------
  subroutine add_goveqns()
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbThermal   , only : thermal_mpp, MPPThermalSetSoils
    use MultiPhysicsProbConstants , only : GE_THERM_SNOW_TBASED
    use MultiPhysicsProbConstants , only : GE_THERM_SSW_TBASED
    use MultiPhysicsProbConstants , only : GE_THERM_SOIL_TBASED
    !
    implicit none

    call thermal_mpp%AddGovEqn(GE_THERM_SNOW_TBASED, &
         'Thermal equation using temprature formulation in snow (KSP formulation)')

    call thermal_mpp%AddGovEqn(GE_THERM_SSW_TBASED, &
         'Thermal equation using temprature formulation in standing surface water (KSP formulation)')

    call thermal_mpp%AddGovEqn(GE_THERM_SOIL_TBASED, &
         'Thermal equation using temprature formulation in soil (KSP formulation)')

    call thermal_mpp%SetMeshesOfGoveqns()
    
  end subroutine add_goveqns

  !------------------------------------------------------------------------
  subroutine add_conditions_to_goveqns()
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbThermal   , only : thermal_mpp, MPPThermalSetSoils
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_SS
    use MultiPhysicsProbConstants , only : COND_HEAT_FLUX
    use MultiPhysicsProbConstants , only : COND_HEAT_RATE
    use MultiPhysicsProbConstants , only : SNOW_TOP_CELLS
    use MultiPhysicsProbConstants , only : SNOW_BOTTOM_CELLS
    use MultiPhysicsProbConstants , only : SSW_TOP_CELLS
    use MultiPhysicsProbConstants , only : SOIL_TOP_CELLS
    use MultiPhysicsProbConstants , only : ALL_CELLS
    use MultiPhysicsProbConstants , only : VAR_DIST_DN
    !
    implicit none
    !
    integer           :: c                        ! do loop indices
    integer           :: ieqn, ieqn_1, ieqn_2
    integer           :: icond
    type(bounds_type) :: bounds_proc
    real(r8), pointer :: z(:,:)                   ! centroid at "z" level [m]
    real(r8), pointer :: zi(:,:)                  ! interface level below a "z" level (m)
    real(r8), pointer :: soil_top_conn_dist_dn(:) !

    z          =>    col%z
    zi         =>    col%zi

    call get_proc_bounds(bounds_proc)

    allocate (soil_top_conn_dist_dn(bounds_proc%endc_all-bounds_proc%begc_all+1))

    do c = bounds_proc%begc, bounds_proc%endc
       soil_top_conn_dist_dn(c - bounds_proc%begc + 1) = z(c,1) - zi(c,0)
    end do

    !
    ! Add BC/SS
    !
    ieqn = 1
    call thermal_mpp%GovEqnAddCondition(ieqn, COND_BC,           &
         'Heat_flux_BC_at_top_of_snow', 'W/m^2', COND_HEAT_FLUX, &
         SNOW_TOP_CELLS)
    call thermal_mpp%GovEqnAddCondition(ieqn, COND_SS,           &
         'Absorbed_solar_radiation', 'W/m^2', COND_HEAT_RATE, &
         ALL_CELLS)

    ieqn = 2
    call thermal_mpp%GovEqnAddCondition(ieqn, COND_BC,           &
         'Heat_flux_BC_at_top_of_standing_surface_water', 'W/m^2', COND_HEAT_FLUX, &
         SSW_TOP_CELLS)

    ieqn = 3
    call thermal_mpp%GovEqnAddCondition(ieqn, COND_BC,           &
         'Heat_flux_BC_at_top_of_soil', 'W/m^2', COND_HEAT_FLUX, &
         SOIL_TOP_CELLS)
    call thermal_mpp%GovEqnAddCondition(ieqn, COND_SS,           &
         'Absorbed_solar_radiation', 'W/m^2', COND_HEAT_RATE, &
         ALL_CELLS)

    !
    ! Add coupling boundary condition
    !
    ieqn_1 = 1
    ieqn_2 = 3
    call thermal_mpp%GovEqnAddCouplingCondition(ieqn_1, ieqn_2, &
         SNOW_BOTTOM_CELLS, SOIL_TOP_CELLS)

    ieqn_1 = 2
    ieqn_2 = 3
    call thermal_mpp%GovEqnAddCouplingCondition(ieqn_1, ieqn_2, &
         SSW_TOP_CELLS, SOIL_TOP_CELLS)

    !
    ! Hack to update the interface distance for soil governing equation because
    ! the centroid of a soil layer is not located at a depth which is
    ! average of the bounding interface depths
    !
    !    Location of       Centroid
    !  CLM's centroid          that is equidistant
    !  ---------             ---------
    !      |                     |
    !      o                     |
    !      |                     o
    !      |                     |
    !      |                     |
    !  ---------             ---------
    ! 
    ieqn_2 = 3
    icond  = 2
    call thermal_mpp%GovEqnUpdateBCConnectionSet(ieqn_2, icond, VAR_DIST_DN, &
         bounds_proc%endc - bounds_proc%begc + 1, &
         soil_top_conn_dist_dn)
    icond  = 3
    call thermal_mpp%GovEqnUpdateBCConnectionSet(ieqn_2, icond, VAR_DIST_DN, &
         bounds_proc%endc - bounds_proc%begc + 1, &
         soil_top_conn_dist_dn)

  end subroutine add_conditions_to_goveqns

  !------------------------------------------------------------------------
  subroutine allocate_auxvars()
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbThermal   , only : thermal_mpp
    use MultiPhysicsProbConstants , only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants , only : VAR_THERMAL_COND
    use MultiPhysicsProbConstants , only : VAR_FRAC
    use MultiPhysicsProbConstants , only : VAR_ACTIVE
    use MultiPhysicsProbConstants , only : VAR_DIST_UP
    use MultiPhysicsProbConstants , only : VAR_DIST_DN
    use MultiPhysicsProbConstants , only : VAR_DZ
    !
    implicit none
    !
    integer           :: ieqn
    integer           :: nvars_for_coupling
    integer, pointer  :: var_ids_for_coupling(:)
    integer, pointer  :: goveqn_ids_for_coupling(:)

    !
    ! Allocate auxvars
    !
    call thermal_mpp%AllocateAuxVars()

    !
    ! Set variables to be exchanged among the coupled equations
    !

    ieqn = 1
    nvars_for_coupling = 2
    allocate (var_ids_for_coupling    (nvars_for_coupling))
    allocate (goveqn_ids_for_coupling (nvars_for_coupling))

    var_ids_for_coupling    (1) = VAR_TEMPERATURE
    var_ids_for_coupling    (2) = VAR_THERMAL_COND
    goveqn_ids_for_coupling (1) = 3
    goveqn_ids_for_coupling (2) = 3

    call thermal_mpp%GovEqnSetCouplingVars(ieqn, nvars_for_coupling, &
         var_ids_for_coupling, goveqn_ids_for_coupling)
    

    ieqn = 2
    call thermal_mpp%GovEqnSetCouplingVars(ieqn, nvars_for_coupling, &
         var_ids_for_coupling, goveqn_ids_for_coupling)


    deallocate(var_ids_for_coupling   )
    deallocate(goveqn_ids_for_coupling)

    ieqn = 3
    nvars_for_coupling = 10
    allocate (var_ids_for_coupling    (nvars_for_coupling))
    allocate (goveqn_ids_for_coupling (nvars_for_coupling))

    var_ids_for_coupling    ( 1) = VAR_TEMPERATURE
    var_ids_for_coupling    ( 2) = VAR_THERMAL_COND
    var_ids_for_coupling    ( 3) = VAR_FRAC
    var_ids_for_coupling    ( 4) = VAR_ACTIVE
    var_ids_for_coupling    ( 5) = VAR_DIST_UP
    var_ids_for_coupling    ( 6) = VAR_TEMPERATURE
    var_ids_for_coupling    ( 7) = VAR_THERMAL_COND
    var_ids_for_coupling    ( 8) = VAR_FRAC
    var_ids_for_coupling    ( 9) = VAR_ACTIVE
    var_ids_for_coupling    (10) = VAR_DZ

    goveqn_ids_for_coupling ( 1) = 1
    goveqn_ids_for_coupling ( 2) = 1
    goveqn_ids_for_coupling ( 3) = 1
    goveqn_ids_for_coupling ( 4) = 1
    goveqn_ids_for_coupling ( 5) = 1
    goveqn_ids_for_coupling ( 6) = 2
    goveqn_ids_for_coupling ( 7) = 2
    goveqn_ids_for_coupling ( 8) = 2
    goveqn_ids_for_coupling ( 9) = 2
    goveqn_ids_for_coupling (10) = 2
    
    call thermal_mpp%GovEqnSetCouplingVars(ieqn, nvars_for_coupling, &
         var_ids_for_coupling, goveqn_ids_for_coupling)

    deallocate(var_ids_for_coupling   )
    deallocate(goveqn_ids_for_coupling)

  end subroutine allocate_auxvars
    
  !------------------------------------------------------------------------
  subroutine add_material_properties()
    !
    ! !DESCRIPTION:
    ! Initialization PETSc-based thermal model
    !
    ! !USES:
    use MultiPhysicsProbThermal   , only : thermal_mpp
    use MultiPhysicsProbThermal   , only : MPPThermalSetSoils
    use clm_varcon                , only : spval
    use clm_varpar                , only : nlevgrnd, nlevsno
    !
    ! !ARGUMENTS
    implicit none
    !
    integer           :: c,g,fc,j,l               ! do loop indices
    type(bounds_type) :: bounds_proc
    integer, pointer  :: thermal_filter(:)
    integer, pointer  :: thermal_lun_type(:)
    real(r8), pointer :: thermal_watsat(:,:)
    real(r8), pointer :: thermal_csol(:,:)
    real(r8), pointer :: thermal_tkmg(:,:)
    real(r8), pointer :: thermal_tkdry(:,:)
    real(r8), pointer :: clm_watsat(:,:)
    real(r8), pointer :: clm_csol(:,:)
    real(r8), pointer :: clm_tkmg(:,:)
    real(r8), pointer :: clm_tkdry(:,:)

    clm_watsat =>    soilstate_vars%watsat_col
    clm_csol   =>    soilstate_vars%csol_col
    clm_tkmg   =>    soilstate_vars%tkmg_col
    clm_tkdry  =>    soilstate_vars%tkdry_col

    call get_proc_bounds(bounds_proc)

    ! Allocate memory and setup data structure for VSFM-MPP

    allocate (thermal_filter   (bounds_proc%begc:bounds_proc%endc_all      ))
    allocate (thermal_lun_type (bounds_proc%begc:bounds_proc%endc_all      ))
    allocate (thermal_watsat   (bounds_proc%begc:bounds_proc%endc,nlevgrnd ))
    allocate (thermal_csol     (bounds_proc%begc:bounds_proc%endc,nlevgrnd ))
    allocate (thermal_tkmg     (bounds_proc%begc:bounds_proc%endc,nlevgrnd ))
    allocate (thermal_tkdry    (bounds_proc%begc:bounds_proc%endc,nlevgrnd ))

    thermal_filter(:)   = 0
    thermal_lun_type(:) = 0
    thermal_watsat(:,:) = 0._r8
    thermal_csol(:,:)   = 0._r8
    thermal_tkmg(:,:)   = 0._r8
    thermal_tkdry(:,:)  = 0._r8

    do c = bounds_proc%begc, bounds_proc%endc
       l = col%landunit(c)

       thermal_lun_type(c)  = lun%itype(l)

       if (col%active(c) .and. .not.lun%lakpoi(l) .and. .not.lun%urbpoi(l)) then
          thermal_filter(c) = 1
          do j = 1,nlevgrnd
             if (clm_watsat(c,j) /= spval) thermal_watsat(c,j) = clm_watsat(c,j)
             if (clm_csol(  c,j) /= spval) thermal_csol  (c,j) = clm_csol  (c,j)
             if (clm_tkmg(  c,j) /= spval) thermal_tkmg  (c,j) = clm_tkmg  (c,j)
             if (clm_tkdry( c,j) /= spval) thermal_tkdry (c,j) = clm_tkdry (c,j)
          enddo

       endif
    enddo

    call MPPThermalSetSoils(thermal_mpp, bounds_proc%begc, &
                           bounds_proc%endc,               &
                           thermal_filter,                 &
                           thermal_lun_type,               &
                           thermal_watsat,                 &
                           thermal_csol,                   &
                           thermal_tkmg,                   &
                           thermal_tkdry                   &
                           )

  end subroutine add_material_properties

#endif

end module MPPThermalTBasedALM_Initialize
