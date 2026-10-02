module ColumnConnectionSetType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Column-to-column connections for lateral flow between grid cells.
  !
  ! A connection joins the naturally vegetated columns of two grid cells that
  ! share an edge. Its "up" column belongs to the cell with the smaller natural
  ! (global) id, so the orientation does not depend on the domain decomposition.
  ! Connections between two ghost columns are not stored.
  !
  ! For each owned column, col_conn_id(col_conn_beg(c):col_conn_beg(c)+col_nconn(c)-1)
  ! lists its connections sorted by the natural id of the neighboring cell. Summing
  ! per-connection fluxes in this order gives results that are independent of the
  ! number of MPI ranks.
  !
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use decompMod      , only : bounds_type
  use abortutils     , only : endrun
  use ColumnType     , only : col_pp
  implicit none
  save
  public

  type, public :: col_connection_set_type
     integer           :: nconn                          ! number of local connections
     integer           :: nconn_global                   ! number of unique connections across all MPI ranks
     integer , pointer :: col_id_up(:)         => null() ! [nconn] up column
     integer , pointer :: col_id_dn(:)         => null() ! [nconn] down column
     integer , pointer :: grid_id_up(:)        => null() ! [nconn] up grid cell
     integer , pointer :: grid_id_dn(:)        => null() ! [nconn] down grid cell
     integer , pointer :: grid_id_up_norder(:) => null() ! [nconn] natural id of the up grid cell
     integer , pointer :: grid_id_dn_norder(:) => null() ! [nconn] natural id of the down grid cell
     real(r8), pointer :: dist(:)              => null() ! [nconn] distance between cell centers along the land surface [m]
     real(r8), pointer :: face_length(:)       => null() ! [nconn] horizontal length of the shared edge [m]
     real(r8), pointer :: uparea(:)            => null() ! [nconn] horizontal area of the up cell [m^2]
     real(r8), pointer :: downarea(:)          => null() ! [nconn] horizontal area of the down cell [m^2]
     real(r8), pointer :: dzg(:)               => null() ! [nconn] elevation of the down cell minus the up cell [m]
     real(r8), pointer :: facecos(:)           => null() ! [nconn] cosine of the slope between cell centers [-]

     integer , pointer :: col_nconn(:)         => null() ! [begc:endc] number of connections of each owned column
     integer , pointer :: col_conn_beg(:)      => null() ! [begc:endc] position of the column's first entry in col_conn_id
     integer , pointer :: col_conn_id(:)       => null() ! connection ids, sorted by neighbor natural id within each column
     real(r8), pointer :: col_conn_sign(:)     => null() ! +1 if the column is the down column of the connection, -1 if up
   contains
#ifdef MOAB_LATERAL
     procedure, public :: Init => InitViaMOAB
#endif
  end type col_connection_set_type

  type (col_connection_set_type), public, target :: c2c_connections   ! connection type

contains

#ifdef MOAB_LATERAL
  !------------------------------------------------------------------------
  subroutine InitViaMOAB(this, bounds_proc)
    !
    use MOABGridType    , only : moab_edge_internal, moab_gcell
    use landunit_varcon , only : istsoil
    use elm_varctl      , only : iulog, use_lateral_subsurface_flow
    use spmdMod         , only : masterproc, mpicom, MPI_INTEGER, MPI_SUM
    !
    implicit none
    !
    class (col_connection_set_type) :: this
    type(bounds_type), intent(in)   :: bounds_proc       ! bound information at processor level
    !
    integer            :: g, e, g_up_moab, g_dn_moab, g_up_elm, g_dn_elm
    integer            :: c, c_up, c_dn
    integer            :: iconn, nconn, nconn_owned_up, n, k, i, beg, ierr
    integer            :: id_tmp, norder_tmp
    real(r8)           :: dc, sign_tmp
    integer, pointer   :: nat_col_id(:)
    integer, pointer   :: nbr_norder(:)
    integer, pointer   :: fill_count(:)

    if (use_lateral_subsurface_flow .and. .not. moab_gcell%has_elevation) then
       call endrun('ERROR: use_lateral_subsurface_flow requires an "elevation" variable '// &
            'in the land domain file.')
    end if

    ! determine the naturally vegetated column of each (owned or ghost) grid cell
    allocate(nat_col_id(bounds_proc%begg_all:bounds_proc%endg_all))
    nat_col_id(:) = -1

    do c = bounds_proc%begc_all, bounds_proc%endc_all
       if (col_pp%itype(c) == istsoil) then
          g = col_pp%gridcell(c)
          if (nat_col_id(g) /= -1) then
             call endrun('ERROR: More than one naturally vegetated column found.')
          end if
          nat_col_id(g) = c
       end if
    end do

    ! count connections, skipping those between two ghost columns
    nconn = 0
    do e = 1, moab_edge_internal%num
       if (.not. IsConnected(e, c_up, c_dn)) cycle
       nconn = nconn + 1
    end do

    this%nconn = nconn
    allocate(this%col_id_up(nconn))         ;  this%col_id_up(:)         = 0
    allocate(this%col_id_dn(nconn))         ;  this%col_id_dn(:)         = 0
    allocate(this%grid_id_up(nconn))        ;  this%grid_id_up(:)        = 0
    allocate(this%grid_id_dn(nconn))        ;  this%grid_id_dn(:)        = 0
    allocate(this%grid_id_up_norder(nconn)) ;  this%grid_id_up_norder(:) = 0
    allocate(this%grid_id_dn_norder(nconn)) ;  this%grid_id_dn_norder(:) = 0
    allocate(this%face_length(nconn))       ;  this%face_length(:)       = 0._r8
    allocate(this%uparea(nconn))            ;  this%uparea(:)            = 0._r8
    allocate(this%downarea(nconn))          ;  this%downarea(:)          = 0._r8
    allocate(this%dist(nconn))              ;  this%dist(:)              = 0._r8
    allocate(this%dzg(nconn))               ;  this%dzg(:)               = 0._r8
    allocate(this%facecos(nconn))           ;  this%facecos(:)           = 0._r8

    ! fill connection data
    iconn = 0
    nconn_owned_up = 0
    do e = 1, moab_edge_internal%num
       if (.not. IsConnected(e, c_up, c_dn)) cycle
       iconn = iconn + 1

       g_up_moab = moab_edge_internal%cell_ids(e, 1)
       g_dn_moab = moab_edge_internal%cell_ids(e, 2)

       this%col_id_up(iconn)         = c_up
       this%col_id_dn(iconn)         = c_dn
       this%grid_id_up(iconn)        = moab_gcell%moab2elm(g_up_moab)
       this%grid_id_dn(iconn)        = moab_gcell%moab2elm(g_dn_moab)
       this%grid_id_up_norder(iconn) = moab_gcell%natural_id(g_up_moab)
       this%grid_id_dn_norder(iconn) = moab_gcell%natural_id(g_dn_moab)

       dc                      = moab_edge_internal%dc(e)
       this%face_length(iconn) = moab_edge_internal%lv(e)
       this%uparea(iconn)      = moab_gcell%area(g_up_moab)
       this%downarea(iconn)    = moab_gcell%area(g_dn_moab)
       this%dzg(iconn)         = moab_gcell%elevation(g_dn_moab) - moab_gcell%elevation(g_up_moab)
       this%dist(iconn)        = sqrt(dc**2._r8 + this%dzg(iconn)**2._r8)
       this%facecos(iconn)     = dc / this%dist(iconn)

       ! each connection is counted once globally: by the rank that owns its up column
       if (c_up <= bounds_proc%endc) nconn_owned_up = nconn_owned_up + 1
    end do

    ! build per-column connection lists for owned columns
    allocate(this%col_nconn   (bounds_proc%begc:bounds_proc%endc)) ; this%col_nconn(:)    = 0
    allocate(this%col_conn_beg(bounds_proc%begc:bounds_proc%endc)) ; this%col_conn_beg(:) = 0

    do iconn = 1, nconn
       c_up = this%col_id_up(iconn)
       c_dn = this%col_id_dn(iconn)
       if (c_up <= bounds_proc%endc) this%col_nconn(c_up) = this%col_nconn(c_up) + 1
       if (c_dn <= bounds_proc%endc) this%col_nconn(c_dn) = this%col_nconn(c_dn) + 1
    end do

    n = 0
    do c = bounds_proc%begc, bounds_proc%endc
       this%col_conn_beg(c) = n + 1
       n = n + this%col_nconn(c)
    end do

    allocate(this%col_conn_id(n))   ; this%col_conn_id(:)   = 0
    allocate(this%col_conn_sign(n)) ; this%col_conn_sign(:) = 0._r8
    allocate(nbr_norder(n))         ; nbr_norder(:)         = 0
    allocate(fill_count(bounds_proc%begc:bounds_proc%endc)) ; fill_count(:) = 0

    do iconn = 1, nconn
       c_up = this%col_id_up(iconn)
       c_dn = this%col_id_dn(iconn)
       if (c_up <= bounds_proc%endc) then
          k = this%col_conn_beg(c_up) + fill_count(c_up)
          fill_count(c_up)    = fill_count(c_up) + 1
          this%col_conn_id(k)   = iconn
          this%col_conn_sign(k) = -1._r8
          nbr_norder(k)         = this%grid_id_dn_norder(iconn)
       end if
       if (c_dn <= bounds_proc%endc) then
          k = this%col_conn_beg(c_dn) + fill_count(c_dn)
          fill_count(c_dn)    = fill_count(c_dn) + 1
          this%col_conn_id(k)   = iconn
          this%col_conn_sign(k) = +1._r8
          nbr_norder(k)         = this%grid_id_up_norder(iconn)
       end if
    end do

    ! sort each column's list by the natural id of the neighboring cell (insertion sort)
    do c = bounds_proc%begc, bounds_proc%endc
       beg = this%col_conn_beg(c)
       do i = beg + 1, beg + this%col_nconn(c) - 1
          norder_tmp = nbr_norder(i)
          id_tmp     = this%col_conn_id(i)
          sign_tmp   = this%col_conn_sign(i)
          k = i - 1
          do while (k >= beg)
             if (nbr_norder(k) <= norder_tmp) exit
             nbr_norder(k+1)         = nbr_norder(k)
             this%col_conn_id(k+1)   = this%col_conn_id(k)
             this%col_conn_sign(k+1) = this%col_conn_sign(k)
             k = k - 1
          end do
          nbr_norder(k+1)         = norder_tmp
          this%col_conn_id(k+1)   = id_tmp
          this%col_conn_sign(k+1) = sign_tmp
       end do
    end do

    call MPI_Allreduce(nconn_owned_up, this%nconn_global, 1, MPI_INTEGER, MPI_SUM, mpicom, ierr)
    if (masterproc) then
       write(iulog,*) 'c2c_connections: number of column-to-column connections (global) = ', &
            this%nconn_global
    end if

    deallocate(nat_col_id)
    deallocate(nbr_norder)
    deallocate(fill_count)

  contains

    logical function IsConnected(e, c_up, c_dn)
      !
      ! True if both cells of internal edge 'e' have a naturally vegetated column
      ! and at least one of the two columns is owned. Returns the two columns.
      !
      integer, intent(in)  :: e
      integer, intent(out) :: c_up, c_dn

      c_up = nat_col_id(moab_gcell%moab2elm(moab_edge_internal%cell_ids(e, 1)))
      c_dn = nat_col_id(moab_gcell%moab2elm(moab_edge_internal%cell_ids(e, 2)))

      IsConnected = (c_up /= -1 .and. c_dn /= -1)
      if (IsConnected) then
         IsConnected = (c_up <= bounds_proc%endc .or. c_dn <= bounds_proc%endc)
      end if

    end function IsConnected

  end subroutine InitViaMOAB
#endif

end module ColumnConnectionSetType
