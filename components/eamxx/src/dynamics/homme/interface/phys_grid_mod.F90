module phys_grid_mod

  use iso_c_binding, only: c_int, c_double
  use parallel_mod,  only: abortmp, MPIinteger_t, MPIreal_t
  use kinds,         only: iulog
  use mpi

  implicit none
  private

  public :: phys_grids_init, cleanup_grid_init_data
  public :: finalize_phys_grid
  public :: get_my_phys_data
  public :: get_num_local_columns, get_num_global_columns

  ! Available options for pg balancing
  integer, public, parameter :: load_bal_none      = 0
  integer, public, parameter :: load_bal_twin_cols = 1  ! Not yet supported in SCREAM!

  ! Max phys grid N. Convention: pgN=0 -> GLL grid, pgN=1,2,3,4 -> FV grid
  integer, public, parameter :: pgN_min = 0
  integer, public, parameter :: pgN_max = 4
  integer, public, parameter :: num_pgN = pgN_max - pgN_min + 1

  ! Global grid info
  ! Arrays of size nprocs
  integer (kind=c_int), pointer :: g_elem_per_rank(:)     ! Number of elems on each rank
  integer (kind=c_int), pointer :: g_elem_offsets(:)      ! Offset of each rank in global nelem-sized arrays
  ! Arrays of size nelem
  integer (kind=c_int), pointer :: g_elem_gids(:)         ! List of elem gids, spliced by rank

  ! This struct holds the global info of a physics grid.
  ! NOTE: there is one struct for each pg type,
  type :: pg_specs_t
    ! Arrays of size nprocs
    integer (kind=c_int), pointer :: g_dofs_per_rank(:)     ! Number of cols on each rank
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
    integer (kind=c_int), pointer :: g_dofs_per_elem(:)     ! List of num dofs on each elem, spliced by rank

    integer :: pgN
    logical :: inited = .false.
  end type pg_specs_t

  type(pg_specs_t), target :: pg_specs (pgN_min:pgN_max)

  ! Note: in this module, we often use MPI_IN_PLACE,0,MPI_DATATYPE_NULL
  !       for the src array specs in Allgatherv calls. These special values 
  !       inform MPI that src array is aliasing the dst one, so MPI will
  !       grab the src data from the dst array, using the offsets info

contains

  subroutine check_phys_grid_inited (pgN)
    !
    ! Inputs(s)
    !
    integer, intent(in) :: pgN  ! The phys grid type to check
    !
    ! Local(s)
    !
    character(1) :: str

    write(str, '(I1)') pgN
    if (pgN .lt. pgN_min .or. pgN .gt. pgN_max) then
      call abortmp ("Error! Physics grid '"//str//"' is not supported.\n")
    endif

    if (.not. pg_specs(pgN)%inited) then
      call abortmp ("Error! Physics grid '"//str//"' was not inited.\n")
    endif
  end subroutine check_phys_grid_inited

  subroutine phys_grids_init (pg_types)
    use gllfvremap_mod,    only: gfr_init
    use homme_context_mod, only: elem, par
    use dimensions_mod,    only: nelem, nelemd
    !
    ! Input(s)
    !
    integer (kind=c_int), intent(in) :: pg_types(:)
    !
    ! Local(s)
    !
    integer :: ie, ierr, proc, i, load_bal, pgN, fvN
    logical :: do_gll
    character(2) :: str

    fvN = -1
    do_gll = .false.
    do i=1, size(pg_types)
      pgN = mod(pg_types(i),10)
      load_bal = pg_types(i) / 10

      if (load_bal .ne. load_bal_none) then
        write(str,'(I2)') load_bal
        call abortmp("Error! Load balancing option '"//str//"' not supported.")
      endif

      if (pgN .lt. pgN_min .or. pgN .gt. pgN_max) then
        write(str,'(I2)') pgN
        call abortmp("Error! Physics grid type '"//str//"' not supported.")
      endif

      if (pgN > 0) then
        if (fvN .ne. -1 .and. fvN .ne. pgN) then
          call abortmp("Error! Only ONE finite volume phys grid can be enabled.")
        endif
        fvN = pgN
      else
        do_gll = .true.
      endif
    enddo
    
    ! Compute elem-related quantities, which are common to all phys grids

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

    ! Now init grid-specific quantities for the requested pg types
    if (do_gll) then
      call phys_grid_init (0)
    endif
    if (fvN .gt. 0) then
      call phys_grid_init (fvN)
    endif
  end subroutine phys_grids_init

  !   integer :: ldofs, ierr, proc, ngcols, ie, gid, offset
  !   integer, allocatable :: dofs_per_elem (:)

  !   ! Note: in the following, we often use MPI_IN_PLACE,0,MPI_DATATYPE_NULL
  !   !       for the src array specs in Allgatherv calls. These special values 
  !   !       inform MPI that src array is aliasing the dst one, so MPI will
  !   !       grab the src data from the dst array, using the offsets info

  !   if (pgN>0) then
  !     fv_nphys = pgN
  !     call gfr_init(par, elem, fv_nphys)
  !   endif

  !   ! Set this right away, so calls to get_num_[local|global]_columns doesn't crap out
  !   is_phys_grid_inited = .true.

  !   ! Gather num elems on each rank
  !   allocate(g_elem_per_rank(par%nprocs))
  !   call MPI_Allgather(nelemd, 1, MPIinteger_t, g_elem_per_rank, 1, MPIinteger_t, par%comm, ierr)

  !   ! Compute offset of each rank in the g_elem_per_rank array
  !   ! This allows each rank to know where to fill arrays of size num_global_elements
  !   allocate(g_elem_offsets(par%nprocs))
  !   g_elem_offsets = 0
  !   do proc=2,par%nprocs
  !     g_elem_offsets(proc) = g_elem_offsets(proc-1) + g_elem_per_rank(proc-1)
  !   enddo

  !   ! Gather elem GIDs on each rank
  !   allocate(g_elem_gids(nelem))
  !   do ie=1,nelemd
  !     g_elem_gids(g_elem_offsets(par%rank+1)+ie) = elem(ie)%globalID
  !   enddo
  !   call MPI_Allgatherv( MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
  !                        g_elem_gids, g_elem_per_rank, g_elem_offsets, MPIinteger_t, par%comm, ierr)

  !   ! Gather the number of unique cols on each rank
  !   allocate(g_dofs_per_rank(par%nprocs))
  !   call MPI_Allgather(get_num_local_columns(), 1, MPIinteger_t, g_dofs_per_rank, 1, MPIinteger_t, par%comm, ierr)

  !   ! Gather the number of unique cols on each element
  !   ! NOTE: we use a temp (dofs_per_elem) rather than g_dofs, since in g_dofs we order by
  !   !       elem GID. But elem GIDs may not be distributed linearly across ranks
  !   !       (e.g., proc0 may own [1-20,41-60],and proc1 may owns [21-40]), while in order
  !   !       to use Allgatherv, we need the send buffer to be contiguous. So in dofs_per_elem
  !   !       we order the data *by rank*. Once the gather is complete, we copy the data from
  !   !       dofs_per_elem into g_dofs_per_elem, ordering by elem GID instead.
  !   !
  !   allocate(dofs_per_elem(nelem))
  !   do ie=1,nelemd
  !     dofs_per_elem(g_elem_offsets(par%rank+1)+ie) = elem(ie)%idxP%NumUniquePts
  !   enddo
  !   call MPI_Allgatherv( MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
  !                        dofs_per_elem, g_elem_per_rank, g_elem_offsets, MPIinteger_t, par%comm, ierr)

  !   allocate(g_dofs_per_elem(nelem))
  !   do proc=1,par%nprocs
  !     offset = g_elem_offsets(proc)
  !     do ie=1,g_elem_per_rank(proc)
  !       g_dofs_per_elem(g_elem_gids(offset+ie)) = dofs_per_elem(offset+ie)
  !     enddo
  !   enddo

  !   ! Compute global dofs
  !   call compute_global_dofs ()

  !   ! Compute area and lat/lon coords
  !   call compute_global_coords ()
  !   call compute_global_area ()

  ! end subroutine phys_grid_init

  subroutine cleanup_grid_init_data ()
    !
    ! Local(s)
    !
    integer :: pgN

    ! Cleanup data no longer needed in each pg_specs_t struct
    do pgN = pgN_min,pgN_max
      if (pg_specs(pgN)%inited) then
        ! nprocs-sized arrays
        deallocate(pg_specs(pgN)%g_dofs_per_rank)

        ! nelem-sized arrays
        deallocate(pg_specs(pgN)%g_dofs_per_elem)
      endif
    enddo

    ! Cleanup global elem-related arrays
    deallocate(g_elem_per_rank)
    deallocate(g_elem_offsets)
    deallocate(g_elem_gids)
  end subroutine cleanup_grid_init_data

  subroutine finalize_phys_grid ()
    use gllfvremap_mod, only: gfr_finish
    !
    ! Local(s)
    !
    integer :: pgN

    do pgN = pgN_min,pgN_max
      if (pg_specs(pgN)%inited) then
        deallocate(pg_specs(pgN)%g_area)
        deallocate(pg_specs(pgN)%g_lat)
        deallocate(pg_specs(pgN)%g_lon)
        deallocate(pg_specs(pgN)%g_dofs)

        deallocate(pg_specs(pgN)%g_dofs_offsets)

        if (pgN>0) then
          call gfr_finish()
        endif

        pg_specs(pgN)%inited = .false.
      endif
    enddo
  end subroutine finalize_phys_grid

  function get_num_local_columns (pgN) result (ncols)
    use dimensions_mod,    only: nelemd
    use homme_context_mod, only: elem
    !
    ! Input(s)
    !
    integer, intent(in) :: pgN
    !
    ! Local(s)
    !
    integer :: ncols, ie

    ! Sanity check
    call check_phys_grid_inited(pgN)

    if (pgN>0) then
      ncols = nelemd*pgN*pgN
    else
      ncols = 0
      do ie=1,nelemd
        ncols = ncols + elem(ie)%idxP%NumUniquePts
      enddo
    endif
  end function get_num_local_columns

  function get_num_global_columns (pgN) result (num_cols) bind(c)
    use dimensions_mod,    only: nelem, np
    !
    ! Input(s)
    !
    integer (kind=c_int), intent(in) :: pgN
    !
    ! Local(s)
    !
    integer (kind=c_int) :: num_cols

    ! Sanity check
    call check_phys_grid_inited(pgN)

    if (pgN>0) then
      num_cols = nelem*pgN*pgN
    else
      num_cols = nelem*(np-1)*(np-1)+2
    endif
  end function get_num_global_columns

  subroutine get_my_phys_data (gids, lat, lon, area, pg_type)
    use homme_context_mod, only: iam
    use shr_const_mod,     only: pi=>SHR_CONST_PI
    use control_mod,       only: geometry
    !
    ! Input(s)
    !
    integer(kind=c_int), intent(out) :: gids(:)
    real(kind=c_double), intent(out) :: lat(:), lon(:), area(:)
    integer(kind=c_int), intent(in) :: pg_type
    !
    ! Local(s)
    !
    integer :: pgN, load_bal, idof, ndofs
    type(pg_specs_t), pointer :: pgs
    logical :: is_sphere

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

    call check_phys_grid_inited (pgN)

    pgs => pg_specs(pgN)

    load_bal = pg_type / 10

    if (load_bal .ne. load_bal_none) then
      call abortmp ("Error! No load balancing option implemented in scream yet.")
    else
      ! TODO: when you enable twin columns, you'll have to manually
      !       do the search, since you can't just grab the offset-ed entries
      ndofs = get_num_local_columns (pgN)
      is_sphere = trim(geometry) /= 'plane'
      do idof=1,ndofs
        gids(idof) = pgs%g_dofs(pgs%g_dofs_offsets(iam+1) + idof)
        lat(idof)  = pgs%g_lat (pgs%g_dofs_offsets(iam+1) + idof)
        lon(idof)  = pgs%g_lon (pgs%g_dofs_offsets(iam+1) + idof)
        if (is_sphere) then
          lat(idof) = lat(idof) * 180.0_c_double / pi
          lon(idof) = lon(idof) * 180.0_c_double / pi
        end if
        area(idof) = pgs%g_area(pgs%g_dofs_offsets(iam+1) + idof)
      enddo

    endif
  end subroutine get_my_phys_data

  subroutine compute_global_dofs (pg)
    use dimensions_mod,    only: nelem
    use homme_context_mod, only: par
    !
    ! Input(s)
    !
    type(pg_specs_t), intent(inout) :: pg
    !
    ! Local(s)
    !
    integer :: ie, icol, idof, proc, elem_gid, elem_offset
    integer, allocatable :: exclusive_scan_dofs_per_elem(:)

    ! This routine's calculations are independent of physics grid type.
    allocate(exclusive_scan_dofs_per_elem(nelem))
    allocate(pg%g_dofs(get_num_global_columns(pg%pgN)))

    exclusive_scan_dofs_per_elem(1) = 0
    do ie=2,nelem
      exclusive_scan_dofs_per_elem(ie) = exclusive_scan_dofs_per_elem(ie-1) + pg%g_dofs_per_elem(ie-1)
    enddo

    allocate (pg%g_dofs_offsets(par%nprocs))

    pg%g_dofs_offsets(1)=0
    do proc=2,par%nprocs
      pg%g_dofs_offsets(proc) = pg%g_dofs_offsets(proc-1) + pg%g_dofs_per_rank(proc-1)
    enddo

    idof = 1
    do proc=1,par%nprocs
      elem_offset = g_elem_offsets(proc)
      do ie=1,g_elem_per_rank(proc)
        elem_gid = g_elem_gids(elem_offset+ie)
        do icol=1,pg%g_dofs_per_elem(elem_gid)
          pg%g_dofs(idof) = exclusive_scan_dofs_per_elem(elem_gid)+icol
          idof = idof+1
        enddo
      enddo
    enddo
  end subroutine compute_global_dofs

  subroutine compute_global_area(pg)
    use dof_mod,           only: UniquePoints
    use dimensions_mod,    only: np, nelemd
    use homme_context_mod, only: elem, par, iam, masterproc
    use gllfvremap_mod,    only: gfr_f_get_area
    !
    ! Input(s)
    !
    type(pg_specs_t), intent(inout) :: pg
    !
    ! Local(s)
    !
    real(kind=c_double), pointer :: area_l(:)
    real(kind=c_double), dimension(np,np)  :: areaw
    integer :: ie, offset, start, ierr, ncols, i, j, k

    if (masterproc) then
      write(iulog,*) 'INFO: Non-scalable action: Computing global area in SE dycore.'
    endif

    allocate(pg%g_area(get_num_global_columns(pg%pgN)))

    offset = pg%g_dofs_offsets(iam+1)

    ncols = get_num_local_columns(pg%pgN)
    area_l => pg%g_area(offset+1 : offset+ncols)

    if (pg%pgN > 0) then
       start = 0
       do ie=1,nelemd
          do j = 1,pg%pgN
             do i = 1,pg%pgN
                k = start + (j-1)*pg%pgN + i
                area_l(k) = gfr_f_get_area(ie, i, j)
             end do
          end do
          start = start + pg%pgN*pg%pgN
       end do
    else
      ! physics is on GLL grid
      start = 1
      do ie=1,nelemd
        areaw = 1.0_c_double / elem(ie)%rspheremp(:,:)
        ncols = elem(ie)%idxP%NumUniquePts
        call UniquePoints(elem(ie)%idxP, areaw, area_l(start:start+ncols-1))
        start = start + ncols
      enddo
    endif

    ncols = get_num_local_columns(pg%pgN)
    ! Note: using MPI_IN_PLACE,0,MPI_DATATYPE_NULL for the src array
    !       informs MPI that src array is aliasing the dst one, so
    !       MPI will grab the src part from dst, using the offsets info
    call MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
                        pg%g_area, pg%g_dofs_per_rank, pg%g_dofs_offsets, MPIreal_t, &
                        par%comm, ierr)
  end subroutine compute_global_area

  subroutine compute_global_coords(pg)
    use dimensions_mod,    only: nelemd
    use dof_mod,           only: UniqueCoords
    use homme_context_mod, only: elem, par, iam, masterproc
    use gllfvremap_mod,    only: gfr_f_get_latlon
    !
    ! Input(s)
    !
    type(pg_specs_t), intent(inout) :: pg
    !
    ! Local(s)
    !
    real(kind=c_double), pointer :: lat_l(:), lon_l(:)
    integer  :: ncols, ie, offset, start, ierr, i, j, k

    if (masterproc) then
      write(iulog,*) 'INFO: Non-scalable action: Computing global coords in SE dycore.'
    end if

    allocate(pg%g_lat(get_num_global_columns(pg%pgN)))
    allocate(pg%g_lon(get_num_global_columns(pg%pgN)))

    offset = pg%g_dofs_offsets(iam+1)

    ncols = get_num_local_columns(pg%pgN)
    lat_l => pg%g_lat(offset+1 : offset+ncols)
    lon_l => pg%g_lon(offset+1 : offset+ncols)

    if (pg%pgN > 0) then
       start = 0
       do ie=1,nelemd
          do j = 1,pg%pgN
             do i = 1,pg%pgN
                k = start + (j-1)*pg%pgN + i
                call gfr_f_get_latlon(ie, i, j, lat_l(k), lon_l(k))
             end do
          end do
          start = start + pg%pgN*pg%pgN
       end do
    else
      ! physics is on GLL grid
      start = 1
      do ie=1,nelemd
        ncols = elem(ie)%idxP%NumUniquePts
        call UniqueCoords(elem(ie)%idxP, elem(ie)%spherep, &
                          lat_l(start:start+ncols-1),      &
                          lon_l(start:start+ncols-1))
        start = start + ncols
      enddo
    endif

    ncols = get_num_local_columns(pg%pgN)

    ! Note: using MPI_IN_PLACE,0,MPI_DATATYPE_NULL for the src array
    !       informs MPI that src array is aliasing the dst one, so
    !       MPI will grab the src part from dst, using the offsets info
    call MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
                        pg%g_lat, pg%g_dofs_per_rank, pg%g_dofs_offsets, MPIreal_t, &
                        par%comm, ierr)
    call MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
                        pg%g_lon, pg%g_dofs_per_rank, pg%g_dofs_offsets, MPIreal_t, &
                        par%comm, ierr)
  end subroutine compute_global_coords

  subroutine phys_grid_init (pgN)
    use gllfvremap_mod,    only: gfr_init
    use homme_context_mod, only: elem, par
    use dimensions_mod,    only: nelem, nelemd
#ifdef HAVE_MOAB
    use seq_comm_mct,      only: MHID, MHFID  ! id of homme moab coarse and fine applications
    use seq_comm_mct,      only: ATMID
    use seq_comm_mct,      only: mhpgid       ! id of pgx moab application
    use semoab_mod,        only: create_moab_meshes
    use iMOAB, only : iMOAB_RegisterApplication
    use iso_c_binding
#endif
    !
    ! Input(s)
    !
    integer, intent(in), target :: pgN
    !
    ! Local(s)
    !
    integer, allocatable :: dofs_per_elem (:)
    integer :: ierr, proc, ie, offset
    character(2) :: str
    type(pg_specs_t), pointer :: pg

#ifdef HAVE_MOAB
    integer :: ATM_ID1
    character*32  appname
#endif
    pg => pg_specs(pgN)

    if (pg%inited) then
      write (str, '(I1)') pgN
      call abortmp("Error! Physics grid with pgN="//str//" was already inited.")
    endif

    ! Set this to true immediately, so we can call get_num_xyz_columns
    pg%inited = .true.
    pg%pgN = pgN

    if (pgN>0) call gfr_init(par, elem, pgN)

    ! Gather the number of unique cols on each rank
    allocate(pg%g_dofs_per_rank(par%nprocs))
    call MPI_Allgather(get_num_local_columns(pg%pgN), 1, MPIinteger_t, &
                       pg%g_dofs_per_rank, 1, MPIinteger_t, par%comm, ierr)

    ! Gather the number of unique cols on each element
    ! NOTE: we use a temp (dofs_per_elem) rather than g_dofs, since in g_dofs we order by
    !       elem GID. But elem GIDs may not be distributed linearly across ranks
    !       (e.g., proc0 may own [1-20,41-60],and proc1 may owns [21-40]), while in order
    !       to use Allgatherv, we need the send buffer to be contiguous. So in dofs_per_elem
    !       we order the data *by rank*. Once the gather is complete, we copy the data from
    !       dofs_per_elem into g_dofs_per_elem, ordering by elem GID instead.
    !
    allocate(dofs_per_elem(nelem))
    if (pgN>0) then
       do ie=1,nelem
          dofs_per_elem(ie) = pgN*pgN
       end do
    else
       do ie=1,nelemd
          dofs_per_elem(g_elem_offsets(par%rank+1)+ie) = elem(ie)%idxP%NumUniquePts
       enddo
       call MPI_Allgatherv( MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
            dofs_per_elem, g_elem_per_rank, g_elem_offsets, MPIinteger_t, par%comm, ierr)
    end if

    allocate(pg%g_dofs_per_elem(nelem))
    do proc=1,par%nprocs
      offset = g_elem_offsets(proc)
      do ie=1,g_elem_per_rank(proc)
        pg%g_dofs_per_elem(g_elem_gids(offset+ie)) = dofs_per_elem(offset+ie)
      enddo
    enddo

    ! Compute arrays of global coords, area, dofs
    call compute_global_dofs (pg)
    call compute_global_coords (pg)
    call compute_global_area (pg)
#ifdef HAVE_MOAB
    if (pgN > 0) then
       appname="HM_COARSE"//C_NULL_CHAR
       ATM_ID1 = 120 !
       ierr = iMOAB_RegisterApplication(appname, par%comm, ATM_ID1, MHID)
       if (ierr > 0 )  &
       call abortmp('Error: cannot register moab app')
       if(par%masterproc) then
           write(iulog,*) " "
           write(iulog,*) "register MOAB app:", trim(appname), "  MHID=", MHID
           write(iulog,*) " "
       endif
       appname="HM_FINE"//C_NULL_CHAR
       ATM_ID1 = 119 ! this number should not conflict with other components IDs; how do we know?
       ierr = iMOAB_RegisterApplication(appname, par%comm, ATM_ID1, MHFID)
       if (ierr > 0 )  &
       call abortmp('Error: cannot register moab app for fine mesh')
       if(par%masterproc) then
           write(iulog,*) " "
           write(iulog,*) "register MOAB app:", trim(appname), "  MHFID=", MHFID
           write(iulog,*) " "
       endif
         appname="HM_PGX"//C_NULL_CHAR
         ATM_ID1 =  ATMID(1) ! this number should not conflict with other components IDs; how do we know?
         !
         ! in this case, we reuse the main atm id, mhid will not be used for intersection anymore
         ! still, need to be careful
         ierr = iMOAB_RegisterApplication(appname, par%comm, ATM_ID1, mhpgid)
         if (ierr > 0 )  &
             call abortmp('Error: cannot register moab app for fine mesh')
         if(par%masterproc) then
             write(iulog,*) " "
             write(iulog,*) "register MOAB app:", trim(appname), "  MHPGID=", mhpgid
             write(iulog,*) " "
         endif
! instance distributed moab meshes from elem structures
!    1 ) spectral coarse mesh
!    2 ) GLL fine quad mesh (used mostly for visualization)
!    3 ) pgN FV type mesh, (most of the time pg2 mesh), used for coupling with other components;
       call create_moab_meshes(par, elem, pgN)
     endif
#endif
  end subroutine phys_grid_init


end module phys_grid_mod
