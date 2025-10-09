module ColumnConnectionSetType


  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : isnan => shr_infnan_isnan, nan => shr_infnan_nan, assignment(=)
  use decompMod      , only : bounds_type
  use abortutils     , only : endrun
  use ColumnType     , only : col_pp
  implicit none
  save
  public

  type, public :: col_connection_set_type
     Integer           :: nconn                          ! number of connections
     Integer, pointer  :: col_id_up(:)         => null() ! list of ids of upwind cells
     Integer, pointer  :: col_id_dn(:)         => null() ! list of ids of downwind cells
     Integer, pointer  :: grid_id_up(:)        => null() ! list of ids of upwind cells
     Integer, pointer  :: grid_id_dn(:)        => null() ! list of ids of downwind cells
     integer, pointer  :: grid_id_up_norder(:) => null() ! list of ids of upwind cells in natural order
     integer, pointer  :: grid_id_dn_norder(:) => null() ! list of ids of downwind cells in natural order
     integer, pointer  :: col_up_forder(:)     => null() ! the order in which the lateral flux should be added for upwind cells
     integer, pointer  :: col_dn_forder(:)     => null() ! the order in which the lateral flux should be added for downwind cells
     Real(r8), pointer :: dist(:)              => null() ! list of distance vectors
     Real(r8), pointer :: face_length(:)       => null() ! list of edge of faces normal to distance vectors
     Real(r8), pointer :: uparea(:)            => null() ! list of up cell areas of horizaontal faces 
     Real(r8), pointer :: downarea(:)          => null() ! list of down cell areas of horizaontal faces 
     Real(r8), pointer :: dzg(:)               => null() ! list of areas of dz between downwind and upwind cells
     Real(r8), pointer :: facecos(:)           => null() ! dot product of the cell face normal vector and cell centroid vector
     Real(r8), pointer :: vertcos(:)           => null() ! dot product of the cell face normal vector and cell centroid vector for vertical flux, the rank for vertcos
     ! is from 1 to column size which is different from rank of lateral faces
   contains
#ifdef HAVE_MOAB
   procedure, public :: Init => InitViaMOAB
#endif
  end type col_connection_set_type

  type (col_connection_set_type), public, target :: c2c_connections   ! connection type

contains

#ifdef HAVE_MOAB
  !------------------------------------------------------------------------
  subroutine InitViaMOAB(this, bounds_proc)
    !
    use MOABGridType, only : moab_edge_internal, moab_gcell
    use decompMod   , only : bounds_type
    !
    implicit none
    !
    class (col_connection_set_type) :: this
    type(bounds_type), intent(in)   :: bounds_proc       ! bound information at processor level
    !
    integer            :: g, g_up_moab, g_dn_moab, g_up_elm, g_dn_elm
    integer            :: c, c_up, c_dn
    integer            :: iconn, nconn
    integer, parameter :: nat_veg_col_itype = 1
    integer, pointer   :: nat_col_id(:)


    ! allocate memory and initialize
    allocate(nat_col_id(bounds_proc%begg_all:bounds_proc%endg_all))
    nat_col_id(:) = -1

    ! loop over columns to determine the naturally-vegetated column for each grid cell.
    do c = bounds_proc%begc_all, bounds_proc%endc_all
       if (col_pp%itype(c) == nat_veg_col_itype) then
          g = col_pp%gridcell(c)

          if (nat_col_id(g) /= -1) then
             call endrun('ERROR: More than one naturally vegetated column found.')
          end if

          nat_col_id(g) = c
       end if
    end do

    ! loop over grid level connections and determine number of column level connections
    nconn = 0
    do iconn = 1, moab_edge_internal%num
       g_up_moab = moab_edge_internal%cell_ids(iconn, 1)
       g_dn_moab = moab_edge_internal%cell_ids(iconn, 2)

       g_up_elm = moab_gcell%moab2elm(g_up_moab)
       g_dn_elm = moab_gcell%moab2elm(g_dn_moab)

       if (nat_col_id(g_up_elm) /= -1 .and. nat_col_id(g_dn_elm) /= -1) then
          nconn = nconn + 1
       end if
    end do

    ! allocate and initialize data structure
    this%nconn = nconn
    allocate(this%col_id_up(nconn))         ;  this%col_id_up(:)         = 0
    allocate(this%col_id_dn(nconn))         ;  this%col_id_dn(:)         = 0
    allocate(this%grid_id_up(nconn))        ;  this%grid_id_up(:)        = 0
    allocate(this%grid_id_dn(nconn))        ;  this%grid_id_dn(:)        = 0
    allocate(this%grid_id_up_norder(nconn)) ;  this%grid_id_up_norder(:) = 0
    allocate(this%grid_id_dn_norder(nconn)) ;  this%grid_id_dn_norder(:) = 0
    allocate(this%col_up_forder(nconn))     ;  this%col_up_forder(:)     = 0
    allocate(this%col_dn_forder(nconn))     ;  this%col_dn_forder(:)     = 0
    allocate(this%face_length(nconn))       ;  this%face_length(:)       = 0
    allocate(this%uparea(nconn))            ;  this%uparea(:)            = 0
    allocate(this%downarea(nconn))          ;  this%downarea(:)          = 0
    allocate(this%dist(nconn))              ;  this%dist(:)              = 0
    allocate(this%dzg(nconn))               ;  this%dzg(:)               = 0
    allocate(this%facecos(nconn))           ;  this%facecos(:)           = 0

    nconn = 0
    do iconn = 1, moab_edge_internal%num
       g_up_moab = moab_edge_internal%cell_ids(iconn, 1)
       g_dn_moab = moab_edge_internal%cell_ids(iconn, 2)

       g_up_elm = moab_gcell%moab2elm(g_up_moab)
       g_dn_elm = moab_gcell%moab2elm(g_dn_moab)

       if (nat_col_id(g_up_elm) /= -1 .and. nat_col_id(g_dn_elm) /= -1) then
          nconn = nconn + 1

          this%col_id_up(nconn) = nat_col_id(g_up_elm)
          this%col_id_dn(nconn) = nat_col_id(g_dn_elm)

          this%grid_id_up(nconn) = g_up_elm
          this%grid_id_dn(nconn) = g_dn_elm

          this%grid_id_up_norder(nconn) = moab_gcell%natural_id(g_up_moab)
          this%grid_id_dn_norder(nconn) = moab_gcell%natural_id(g_dn_moab)
       end if
    end do

    ! free up memory
    deallocate(nat_col_id)

  end subroutine InitViaMOAB
#endif

end module ColumnConnectionSetType

