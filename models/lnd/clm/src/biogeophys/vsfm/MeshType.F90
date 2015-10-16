#ifdef USE_PETSC_LIB


module MeshType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Mesh data type allocation
  ! --------------------------------------------------------

  ! !USES:
  use clm_varctl         , only : iulog
  use abortutils         , only : endrun
  use shr_log_mod        , only : errMsg => shr_log_errMsg
  use shr_kind_mod       , only : r8 => shr_kind_r8
  use WaterstateType     , only : waterstate_type
  use SoilHydrologyType  , only : soilhydrology_type
  use ConnectionSetType  , only : connection_set_list_type
  !
  ! !PUBLIC TYPES:
  implicit none
  private
#include "finclude/petscsys.h"

  type, public :: mesh_type
     character (len=256) :: name         ! indentifer for the mesh

     PetscInt            :: ncells       ! Number of cells within the mesh
     PetscInt            :: nlev         ! Number of cells in z
     PetscInt            :: itype        ! identifier

                                         ! centroid
     PetscReal, pointer  :: x(:)         ! [m]
     PetscReal, pointer  :: y(:)         ! [m]
     PetscReal, pointer  :: z(:)         ! [m]

     PetscReal, pointer  :: z_m(:)       ! [m]
     PetscReal, pointer  :: z_p(:)       ! [m]

     PetscReal, pointer  :: dx(:)        ! [m]
     PetscReal, pointer  :: dy(:)        ! [m]
     PetscReal, pointer  :: dz(:)        ! [m]
     PetscReal, pointer  :: vol(:)       ! [m^3]

     PetscBool, pointer  :: is_active(:) ! [true/false]

     PetscInt            :: orientation

     type(connection_set_list_type) :: intrn_conn_set_list

   contains
     procedure, public :: Init
     procedure, public :: Clean
     procedure, public :: Create
     procedure, public :: CreateFromCLMCols
  end type mesh_type

  public :: MeshCreateConnectionSet
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this)
    !
    ! !DESCRIPTION:
    ! Initializes a mesh object
    !
    ! !USES:
    use ConnectionSetType  , only : ConnectionSetListInit
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mesh_type) :: this

    this%name        = ""
    this%ncells      = 0
    this%nlev        = 0
    this%itype       = -1
    this%orientation = -1

    nullify(this%x         )
    nullify(this%y         )
    nullify(this%z         )

    nullify(this%z_m       )
    nullify(this%z_p       )

    nullify(this%dx        )
    nullify(this%dy        )
    nullify(this%dz        )
    nullify(this%vol       )

    nullify(this%is_active )

    call ConnectionSetListInit(this%intrn_conn_set_list)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine Create(this, mesh_itype, begc, endc, waterstate_vars, soilhydrology_vars)
    !
    ! !DESCRIPTION:
    ! Creates a mesh from CLM column level data structure
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : MESH_CLM_SOIL_COL
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mesh_type)                      :: this
    integer                  , intent(in) :: mesh_itype
    integer                  , intent(in) :: begc,endc
    type(waterstate_type)    , intent(in) :: waterstate_vars
    type(soilhydrology_type) , intent(in) :: soilhydrology_vars

    select case(mesh_itype)
    case (MESH_CLM_SOIL_COL)
       call this%CreateFromCLMCols(begc, endc)
    case default
       write(iulog,*)'MeshType: Create() Unknown mesh_type = ',mesh_itype
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine Create

  !------------------------------------------------------------------------
  subroutine CreateFromCLMCols(this, begc, endc)
    !
    ! !DESCRIPTION:
    ! Creates a mesh from CLM column level data structure
    !
    ! !USES:
    use clm_varpar                  , only : nlevgrnd
    use ColumnType                  , only : col
    use LandunitType                , only : lun
    use landunit_varcon             , only : istcrop, istsoil
    use column_varcon               , only : icol_road_perv
    use ConnectionSetType           , only : ConnectionSetNew
    use MultiPhysicsProbConstants   , only : MESH_ALONG_GRAVITY
    use MultiPhysicsProbConstants   , only : MESH_CLM_SOIL_COL
    use ConnectionSetType           , only : connection_set_type
    use ConnectionSetType           , only : ConnectionSetListAddSet
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mesh_type)    :: this
    integer, intent(in) :: begc,endc
    !
    ! !LOCAL VARIABLES:
    PetscInt                          :: c,j,l                              !indices
    PetscInt                          :: icell
    PetscInt                          :: iconn
    PetscInt                          :: nconn
    PetscInt                          :: id_up, id_dn
    PetscInt                          :: first_active_hydro_col_id
    PetscInt                          :: col_id
    type(connection_set_type),pointer :: conn_set

    call this%Init()

    this%name        = "Soil Mesh"
    this%itype       = MESH_CLM_SOIL_COL
    this%ncells      = (endc - begc + 1)*nlevgrnd
    this%nlev        = nlevgrnd
    this%orientation = MESH_ALONG_GRAVITY

    allocate(this%x(this%ncells         ))
    allocate(this%y(this%ncells         ))
    allocate(this%z(this%ncells         ))

    allocate(this%z_m(this%ncells       ))
    allocate(this%z_p(this%ncells       ))

    allocate(this%dx(this%ncells        ))
    allocate(this%dy(this%ncells        ))
    allocate(this%dz(this%ncells        ))
    allocate(this%vol(this%ncells       ))
    allocate(this%is_active(this%ncells ))

    !
    ! Populate location of cell centroids
    !

    first_active_hydro_col_id = -1
    do c = begc, endc
       l = col%landunit(c)

       if (col%active(c) .and. &
            (lun%itype(l) == istsoil .or. col%itype(c) == icol_road_perv .or. &
            lun%itype(l) == istcrop)) then
          if (first_active_hydro_col_id == -1) then
             first_active_hydro_col_id = c
             exit
          endif
       endif
    enddo

    if (first_active_hydro_col_id == -1) then
       write(iulog,*)'No active soil hydrology column found'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    icell = 0
    do c = begc, endc
       do j = 1, this%nlev

          icell = icell + 1
          l = col%landunit(c)

          this%x (icell) = 0.0_r8
          this%y (icell) = 0.0_r8

          this%dx(icell) = 1.0_r8
          this%dy(icell) = 1.0_r8


          if (col%active(c) .and. &
               (lun%itype(l) == istsoil        .or. &
                col%itype(c) == icol_road_perv .or. &
                lun%itype(l) == istcrop)) then
             col_id = c
             this%is_active(icell) = PETSC_TRUE
          else
             col_id = first_active_hydro_col_id
             this%is_active(icell) = PETSC_FALSE
          endif

          this%z(icell) = -0.5_r8*(col%zi(col_id,j-1) + col%zi(col_id,j))
          this%z_m(icell) = -col%zi(col_id,j-1)
          this%z_p(icell) = -col%zi(col_id,j  )

          this%dz(icell) = col%dz(col_id,j)
          this%vol(icell) = this%dx(icell)*this%dy(icell)*this%dz(icell)
       enddo
    enddo

    nconn = (this%nlev - 1)*(endc - begc + 1)
    conn_set => ConnectionSetNew(nconn)

    iconn = 0
    do c = begc, endc
       do j = 1, this%nlev-1
          iconn = iconn + 1

          id_up = (c-begc)*this%nlev + j
          id_dn = id_up + 1

          conn_set%id_up(iconn) = id_up
          conn_set%id_dn(iconn) = id_dn

          conn_set%area(iconn) = 1.0_r8

          conn_set%dist_up(iconn) = 0.5_r8*this%dz(id_up)
          conn_set%dist_dn(iconn) = 0.5_r8*this%dz(id_dn)

          conn_set%dist_unitvec(iconn)%arr(1) = 0._r8
          conn_set%dist_unitvec(iconn)%arr(2) = 0._r8
          conn_set%dist_unitvec(iconn)%arr(3) = -1._r8

       enddo
    enddo

    call ConnectionSetListAddSet(this%intrn_conn_set_list, conn_set)

  end subroutine CreateFromCLMCols

  !------------------------------------------------------------------------
  subroutine MeshCreateConnectionSet(mesh, region_itype, conn_set, ncells, &
       soil_top_cell_offset)
    !
    ! !DESCRIPTION:
    ! Creates a connection set for a mesh that stores:
    ! - ID of control volumes (CVs) upwind and downwind of a connection,
    ! - Area of the face shared between upwind and downwind CVs,
    ! - Unit normal vector from up to down control volume,
    ! - Distance from cell centroid to face centroid for upwind and downwind CV.
    !
    ! !USES:
    use ConnectionSetType           , only : ConnectionSetNew
    use MultiPhysicsProbConstants   , only : SOIL_TOP_CELLS, SOIL_BOTTOM_CELLS, SOIL_CELLS
    use MultiPhysicsProbConstants   , only : SNOW_TOP_CELLS, SSW_TOP_CELLS, ALL_CELLS
    use MultiPhysicsProbConstants   , only : MESH_ALONG_GRAVITY, MESH_AGAINST_GRAVITY
    use ConnectionSetType           , only : connection_set_type
    !
    implicit none
    !
    ! !ARGUMENTS
    type(mesh_type)                   :: mesh
    PetscInt                          :: region_itype
    type(connection_set_type),pointer :: conn_set
    PetscInt, intent(out)             :: ncells
    PetscInt,optional                 :: soil_top_cell_offset
    !
    ! !LOCAL VARIABLES:
    PetscInt                          :: c,j
    PetscInt                          :: icell
    PetscInt                          :: iconn
    PetscInt                          :: nconn
    PetscInt                          :: id_up, id_dn
    PetscInt                          :: ncols
    PetscInt                          :: offset

    ncols = mesh%ncells/mesh%nlev

    offset = 0
    if (present(soil_top_cell_offset)) offset = soil_top_cell_offset

    select case (region_itype)
    case (SOIL_BOTTOM_CELLS, SOIL_TOP_CELLS, SNOW_TOP_CELLS, SSW_TOP_CELLS)

       nconn = ncols
       conn_set => ConnectionSetNew(nconn)

       iconn = 0
       do c = 1, ncols
          iconn = iconn + 1

          select case (mesh%orientation)
          case (MESH_AGAINST_GRAVITY)

             select case(region_itype)
             case (SOIL_TOP_CELLS, SNOW_TOP_CELLS, SSW_TOP_CELLS)
                id_up = -1
                id_dn = mesh%nlev*c + offset

                conn_set%dist_unitvec(iconn)%arr(1) =  0._r8
                conn_set%dist_unitvec(iconn)%arr(2) =  0._r8
                conn_set%dist_unitvec(iconn)%arr(3) =  -1._r8

             case (SOIL_BOTTOM_CELLS)

                id_up = -1
                id_dn = mesh%nlev*(c-1) + 1 + offset

                conn_set%dist_unitvec(iconn)%arr(1) =  0._r8
                conn_set%dist_unitvec(iconn)%arr(2) =  0._r8
                conn_set%dist_unitvec(iconn)%arr(3) =  1._r8
             end select

          case (MESH_ALONG_GRAVITY)
             select case(region_itype)
             case (SOIL_BOTTOM_CELLS)
                id_up = -1
                id_dn = mesh%nlev*c + offset

                conn_set%dist_unitvec(iconn)%arr(1) =  0._r8
                conn_set%dist_unitvec(iconn)%arr(2) =  0._r8
                conn_set%dist_unitvec(iconn)%arr(3) =  1._r8

             case (SOIL_TOP_CELLS, SNOW_TOP_CELLS, SSW_TOP_CELLS)
                id_up = -1
                id_dn = mesh%nlev*(c-1) + 1 + offset

                conn_set%dist_unitvec(iconn)%arr(1) =  0._r8
                conn_set%dist_unitvec(iconn)%arr(2) =  0._r8
                conn_set%dist_unitvec(iconn)%arr(3) = -1._r8
             end select

          case default
             write(iulog,*)'MeshCreateConnectionSet: Unknown mesh%orientation'
             call endrun(msg=errMsg(__FILE__, __LINE__))
          end select

          conn_set%id_up(iconn) = id_up
          conn_set%id_dn(iconn) = id_dn
          conn_set%area(iconn)  = mesh%dx(id_dn)*mesh%dy(id_dn)
          conn_set%dist_up(iconn) = 0.0_r8
          conn_set%dist_dn(iconn) = 0.5_r8*mesh%dz(id_dn)
       enddo

    case (SOIL_CELLS, ALL_CELLS)

       nconn = mesh%ncells
       conn_set => ConnectionSetNew(nconn)

       iconn = 0
       do c = 1,ncols
          do j = 1,mesh%nlev
             iconn = iconn + 1

             id_up = -1
             id_dn = iconn

             conn_set%id_up(iconn) = id_up
             conn_set%id_dn(iconn) = id_dn

             conn_set%dist_unitvec(iconn)%arr(1) =  0._r8
             conn_set%dist_unitvec(iconn)%arr(2) =  0._r8
             conn_set%dist_unitvec(iconn)%arr(3) =  0._r8

             conn_set%area(iconn)  = mesh%dx(id_dn)*mesh%dy(id_dn)

             conn_set%dist_up(iconn) = 0.0_r8
             conn_set%dist_dn(iconn) = 0.0_r8

          enddo
       enddo

    case default

    end select

    ncells = nconn

  end subroutine MeshCreateConnectionSet

  !------------------------------------------------------------------------
  subroutine Clean(this)
    !
    ! !DESCRIPTION:
    ! Release memory allocated to the mesh
    !
    ! !USES:
    use ConnectionSetType, only : ConnectionSetListClean
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mesh_type) :: this

    deallocate(this%x  )
    deallocate(this%y  )
    deallocate(this%z  )

    deallocate(this%z_m)
    deallocate(this%z_p)

    deallocate(this%dx )
    deallocate(this%dy )
    deallocate(this%dz )

    call ConnectionSetListClean(this%intrn_conn_set_list)

  end subroutine Clean

end module MeshType

#endif
