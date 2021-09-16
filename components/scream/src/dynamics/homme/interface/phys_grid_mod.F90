module phys_grid_mod

  use iso_c_binding, only: c_int, c_double
  use parallel_mod,  only: abortmp, MPIinteger_t, MPIreal_t
  use kinds,         only: iulog

  implicit none
  private

  public :: phys_grid_init, cleanup_grid_init_data
  public :: finalize_phys_grid
  public :: get_my_phys_data
  public :: get_num_local_columns, get_num_global_columns

  ! Available options for pg balancing
  integer, public, parameter :: pg_default   = 0
  integer, public, parameter :: pg_twin_cols = 1  ! Not yet supported in SCREAM!

  ! Arrays of size nprocs
  integer (kind=c_int), pointer :: g_dofs_per_rank(:)     ! Number of cols on each rank
  integer (kind=c_int), pointer :: g_elem_per_rank(:)     ! Number of elems on each rank
  integer (kind=c_int), pointer :: g_elem_offsets(:)      ! Offset of each rank in global nelem-sized arrays
  integer (kind=c_int), pointer :: g_dofs_offsets(:)      ! Offset of each rank in global ngcols-sized arrays

  ! Global arrays are spliced by rank, meaning that all entries on rank N come *before*
  ! all entries on rank N+1. That's because the most common use is to grab the values
  ! on the a given rank (usually, my rank), for which one can use the g_dofs_offsets
  ! to establish the start/end points.

  ! Arrays of size ngcols
  real (kind=c_double), pointer :: g_lat(:)     ! Latitude coordinates
  real (kind=c_double), pointer :: g_lon(:)     ! Longitude coordinates
  real (kind=c_double), pointer :: g_area(:)    ! Column area
  integer (kind=c_int), pointer :: g_dofs(:)    ! Column global index

  ! Arrays of size nelem
  integer (kind=c_int), pointer :: g_elem_gids(:)         ! List of elem gids, spliced by rank
  integer (kind=c_int), pointer :: g_dofs_per_elem(:)     ! List of num dofs on each elem, spliced by rank

  integer, public :: fv_nphys = 0

  logical :: is_phys_grid_inited = .false.

! To get MPI_IN_PLACE and MPI_DATATYPE_NULL
#include <mpif.h>

contains

  subroutine check_phys_grid_inited ()
    if (.not. is_phys_grid_inited) then
      call abortmp ("Error! Physics grid was not inited.\n")
    endif
  end subroutine check_phys_grid_inited

  subroutine phys_grid_init (pgN)
    use gllfvremap_mod,    only: gfr_init
    use homme_context_mod, only: elem, par
    use dimensions_mod,    only: nelem, nelemd
    !
    ! Input(s)
    !
    integer (kind=c_int), intent(in) :: pgN
    !
    ! Local(s)
    !
    integer :: ldofs, ierr, proc, ngcols, ie, gid, offset
    integer, allocatable :: dofs_per_elem (:)

    ! Note: in the following, we often use MPI_IN_PLACE,0,MPI_DATATYPE_NULL
    !       for the src array specs in Allgatherv calls. These special values 
    !       inform MPI that src array is aliasing the dst one, so MPI will
    !       grab the src data from the dst array, using the offsets info

    if (pgN>0) then
      fv_nphys = pgN
      call gfr_init(par, elem, fv_nphys)
    endif

    ! Set this right away, so calls to get_num_[local|global]_columns doesn't crap out
    is_phys_grid_inited = .true.

    ! Gather num elems on each rank
    allocate(g_elem_per_rank(par%nprocs))
    call MPI_Allgather(nelemd, 1, MPIinteger_t, g_elem_per_rank, 1, MPIinteger_t, par%comm, ierr)

    ! Compute offset of each rank in the g_elem_per_rank array
    ! This allows each rank to know where to fill arrays of size num_global_elements
    allocate(g_elem_offsets(par%nprocs))
    g_elem_offsets = 0
    do proc=2,par%nprocs
      g_elem_offsets(proc) = g_elem_offsets(proc-1) + g_elem_per_rank(proc-1)
    enddo

    ! Gather elem GIDs on each rank
    allocate(g_elem_gids(nelem))
    do ie=1,nelemd
      g_elem_gids(g_elem_offsets(par%rank+1)+ie) = elem(ie)%globalID
    enddo
    call MPI_Allgatherv( MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
                         g_elem_gids, g_elem_per_rank, g_elem_offsets, MPIinteger_t, par%comm, ierr)

    ! Gather the number of unique cols on each rank
    allocate(g_dofs_per_rank(par%nprocs))
    call MPI_Allgather(get_num_local_columns(), 1, MPIinteger_t, g_dofs_per_rank, 1, MPIinteger_t, par%comm, ierr)

    ! Gather the number of unique cols on each element
    ! NOTE: we use a temp (dofs_per_elem) rather than g_dofs, since in g_dofs we order by
    !       elem GID. But elem GIDs may not be distributed linearly across ranks
    !       (e.g., proc0 may own [1-20,41-60],and proc1 may owns [21-40]), while in order
    !       to use Allgatherv, we need the send buffer to be contiguous. So in dofs_per_elem
    !       we order the data *by rank*. Once the gather is complete, we copy the data from
    !       dofs_per_elem into g_dofs_per_elem, ordering by elem GID instead.
    !
    allocate(dofs_per_elem(nelem))
    do ie=1,nelemd
      dofs_per_elem(g_elem_offsets(par%rank+1)+ie) = elem(ie)%idxP%NumUniquePts
    enddo
    call MPI_Allgatherv( MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
                         dofs_per_elem, g_elem_per_rank, g_elem_offsets, MPIinteger_t, par%comm, ierr)

    allocate(g_dofs_per_elem(nelem))
    do proc=1,par%nprocs
      offset = g_elem_offsets(proc)
      do ie=1,g_elem_per_rank(proc)
        g_dofs_per_elem(g_elem_gids(offset+ie)) = dofs_per_elem(offset+ie)
      enddo
    enddo

    ! Compute global dofs
    call compute_global_dofs ()

    ! Compute area and lat/lon coords
    call compute_global_coords ()
    call compute_global_area ()

  end subroutine phys_grid_init

  subroutine cleanup_grid_init_data ()
    if (is_phys_grid_inited) then
      ! nprocs-sized arrays
      deallocate(g_dofs_per_rank)
      deallocate(g_elem_per_rank)
      deallocate(g_elem_offsets)

      ! nelem-sized arrays
      deallocate(g_dofs_per_elem)
      deallocate(g_elem_gids)
    endif
  end subroutine cleanup_grid_init_data

  subroutine finalize_phys_grid ()
    use gllfvremap_mod, only: gfr_finish

    if (fv_nphys>0) then
      call gfr_finish()
    endif

    deallocate(g_area)
    deallocate(g_lat)
    deallocate(g_lon)
    deallocate(g_dofs)

    deallocate(g_elem_gids)

    is_phys_grid_inited = .false.
  end subroutine finalize_phys_grid

  function get_num_local_columns () result (ncols)
    use dimensions_mod,    only: nelemd
    use homme_context_mod, only: elem
    !
    ! Local(s)
    !
    integer :: ncols, ie

    ! Sanity check
    call check_phys_grid_inited()

    if (fv_nphys>0) then
      ncols = nelemd*fv_nphys*fv_nphys
    else
      ncols = 0
      do ie=1,nelemd
        ncols = ncols + elem(ie)%idxP%NumUniquePts
      enddo
    endif
  end function get_num_local_columns

  function get_num_global_columns () result (num_cols) bind(c)
    use homme_context_mod, only: is_geometry_inited, elem
    use dimensions_mod,    only: nelem, np
    !
    ! Local(s)
    !
    integer (kind=c_int) :: num_cols
    integer :: ie

    ! Sanity check
    call check_phys_grid_inited()

    if (fv_nphys>0) then
      num_cols = nelem*fv_nphys*fv_nphys
    else
      num_cols = nelem*(np-1)*(np-1)+2
    endif
  end function get_num_global_columns

  subroutine get_my_phys_data (gids, lat, lon, area, pg_type)
    use homme_context_mod, only: elem, iam
    use dimensions_mod,    only: nelemd
    use shr_const_mod,     only: pi=>SHR_CONST_PI
    !
    ! Input(s)
    !
    integer(kind=c_int), pointer, intent(in) :: gids(:)
    real(kind=c_double), pointer, intent(in) :: lat(:), lon(:), area(:)
    integer(kind=c_int),          intent(in) :: pg_type
    !
    ! Local(s)
    !
    integer :: pgN, balancing_opt, idof, ndofs

    call check_phys_grid_inited ()

    ! Possible values for pg_type:
    !   0 : physics grid on GLL nodes
    !   2 : physics grid on FV points in 2x2 subcell (PG2)
    !   3 : physics grid on FV points in 3x3 subcell (PG3)
    !   4 : physics grid on FV points in 4x4 subcell (PG4)
    !  10 : GLL points, use twin columns
    !  12 : PG2 points, use twin columns
    !  13 : PG3 points, use twin columns
    !  14 : PG4 points, use twin columns
    ! In general:
    !   *0 => GLL nodes, *N,N>0 => FV points
    !   1* => redistribute using twin columns

    pgN = mod(pg_type,10)
    if (pgN .ne. fv_nphys) then
      call abortmp ("Error! Requested N for phys grid different from what was used at init time.")
    endif

    balancing_opt = pg_type / 10

    if (balancing_opt .ne. pg_default) then
      call abortmp ("Error! Twin columns not yet implemented for physics grid.")
    else
      ! TODO: when you enable twin columns, you'll have to manually
      !       do the search, since you can't just grab the offset-ed entries
      ndofs = get_num_local_columns ()
      do idof=1,ndofs
        gids(idof) = g_dofs(g_dofs_offsets(iam+1) + idof)
        lat(idof)  = g_lat(g_dofs_offsets(iam+1) + idof) * 180.0_c_double / pi
        lon(idof)  = g_lon(g_dofs_offsets(iam+1) + idof) * 180.0_c_double / pi
        area(idof) = g_area(g_dofs_offsets(iam+1) + idof)
      enddo

    endif
  end subroutine get_my_phys_data

  subroutine compute_global_dofs ()
    use dimensions_mod,    only: nelemd,nelem
    use homme_context_mod, only: elem, par, iam
    !
    ! Local(s)
    !
    integer (kind=c_int), pointer :: dofs_l(:)
    integer :: ie, icol, idof, proc, elem_gid, elem_offset
    integer, allocatable :: exclusive_scan_dofs_per_elem(:)

    if (associated(g_dofs)) then
      call abortmp ("Error! compute_global_dofs was already called.")
    endif

    allocate(exclusive_scan_dofs_per_elem(nelem))
    allocate(g_dofs(get_num_global_columns()))

    exclusive_scan_dofs_per_elem(1) = 0
    do ie=2,nelem
      exclusive_scan_dofs_per_elem(ie) = exclusive_scan_dofs_per_elem(ie-1) + g_dofs_per_elem(ie-1)
    enddo

    if (fv_nphys .gt. 0) then
      ! TODO
      call abortmp ("Error! Phys grid N>0 not yet implemented in scream.")
    else
      allocate (g_dofs_offsets(par%nprocs))
      g_dofs_offsets(1)=0
      do proc=2,par%nprocs
        g_dofs_offsets(proc) = g_dofs_offsets(proc-1) + g_dofs_per_rank(proc-1)
      enddo

      idof = 1
      do proc=1,par%nprocs
        elem_offset = g_elem_offsets(proc)
        do ie=1,g_elem_per_rank(proc)
          elem_gid = g_elem_gids(elem_offset+ie)
          do icol=1,g_dofs_per_elem(elem_gid)
            g_dofs(idof) = exclusive_scan_dofs_per_elem(elem_gid)+icol
            idof = idof+1
          enddo
        enddo
      enddo
    endif

  end subroutine compute_global_dofs

  subroutine compute_global_area()
    use dof_mod,           only: UniquePoints
    use dimensions_mod,    only: np, nelemd
    use homme_context_mod, only: elem, par, iam, masterproc
    !
    ! Local(s)
    !
    real(kind=c_double), pointer :: area_l(:)
    real(kind=c_double), dimension(np,np)  :: areaw
    integer  :: ie, offset, start, ierr, ncols

    if (masterproc) then
      write(iulog,*) 'INFO: Non-scalable action: Computing global area in SE dycore.'
    endif

    if (associated(g_area)) then
      call abortmp ("Error! compute_global_area was already called.")
    endif

    allocate(g_area(get_num_global_columns()))

    if (fv_nphys > 0) then
      ! TODO
      call abortmp("Error! Not yet implemented. Why didn't you get an error before?")
    else
      ! physics is on GLL grid
      offset = g_dofs_offsets(iam+1)

      ncols = get_num_local_columns()
      area_l => g_area(offset+1 : offset+ncols)
      start = 1
      do ie=1,nelemd
        areaw = 1.0_c_double / elem(ie)%rspheremp(:,:)
        ncols = elem(ie)%idxP%NumUniquePts
        call UniquePoints(elem(ie)%idxP, areaw, area_l(start:start+ncols-1))
        start = start + ncols
      enddo

      ncols = get_num_local_columns()
      ! Note: using MPI_IN_PLACE,0,MPI_DATATYPE_NULL for the src array
      !       informs MPI that src array is aliasing the dst one, so
      !       MPI will grab the src part from dst, using the offsets info
      call MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
                          g_area, g_dofs_per_rank, g_dofs_offsets, MPIreal_t, &
                          par%comm, ierr)
    endif

  end subroutine compute_global_area

  subroutine compute_global_coords()
    use dimensions_mod,    only: nelemd
    use dof_mod,           only: UniqueCoords
    use gllfvremap_mod,    only: gfr_f_get_latlon
    use homme_context_mod, only: elem, par, iam, masterproc
    use parallel_mod,      only: abortmp
    !
    ! Local(s)
    !
    real(kind=c_double), pointer :: lat_l(:), lon_l(:)
    integer  :: ncols, ie, offset, start, ierr

    if (associated(g_lat) .or. associated(g_lon)) then
      call abortmp ("Error! compute_global_coords was already called.")
    endif

    allocate(g_lat(get_num_global_columns()))
    allocate(g_lon(get_num_global_columns()))

    if (masterproc) then
      write(iulog,*) 'INFO: Non-scalable action: Computing global coords in SE dycore.'
    end if

    if (fv_nphys > 0) then
      ! TODO
      call abortmp("Error! Not yet implemented. Why didn't you get an error before?")
    else
      ! physics is on GLL grid

      offset = g_dofs_offsets(iam+1)

      ncols = get_num_local_columns()
      lat_l => g_lat(offset+1 : offset+ncols)
      lon_l => g_lon(offset+1 : offset+ncols)
      start = 1
      do ie=1,nelemd
        ncols = elem(ie)%idxP%NumUniquePts
        call UniqueCoords(elem(ie)%idxP, elem(ie)%spherep, &
                          lat_l(start:start+ncols-1),      &
                          lon_l(start:start+ncols-1))
        start = start + ncols
      enddo

      ncols = get_num_local_columns()

      ! Note: using MPI_IN_PLACE,0,MPI_DATATYPE_NULL for the src array
      !       informs MPI that src array is aliasing the dst one, so
      !       MPI will grab the src part from dst, using the offsets info
      call MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
                          g_lat, g_dofs_per_rank, g_dofs_offsets, MPIreal_t, &
                          par%comm, ierr)
      call MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
                          g_lon, g_dofs_per_rank, g_dofs_offsets, MPIreal_t, &
                          par%comm, ierr)
    end if

  end subroutine compute_global_coords

end module phys_grid_mod
